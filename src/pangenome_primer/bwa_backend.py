"""Genome-wide primer binding search on a full haplotype assembly.

`bwa aln`/`samse` locates candidate binding positions across a multi-Gb assembly quickly;
we then extract each candidate window and recompute the exact mismatches and 3'-end offsets
with the naive comparator from `binding.py` (which the tests verify). bwa prunes the search
space; the trusted code decides the details. Requires `bwa` + `samtools` + `pysam`.
"""
from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
from pathlib import Path

from .binding import find_binding_sites_naive
from .model import BindingSite


def ensure_index(fasta: str) -> None:
    if not Path(fasta + ".bwt").exists():
        subprocess.run(["bwa", "index", fasta], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    if not Path(fasta + ".fai").exists():
        subprocess.run(["samtools", "faidx", fasta], check=True)


def _candidate_positions(primer_seq: str, fasta: str, max_mm: int) -> set[tuple[str, int]]:
    """(chrom, 0-based pos) candidates from bwa primary + XA alternative hits."""
    import pysam

    with tempfile.TemporaryDirectory() as td:
        q = Path(td) / "q.fa"
        q.write_text(f">p\n{primer_seq}\n")
        sai = Path(td) / "q.sai"
        L = len(primer_seq)
        with open(sai, "wb") as fh:
            subprocess.run(
                ["bwa", "aln", "-n", str(max_mm), "-o", "0", "-l", str(L),
                 "-k", str(max_mm), "-N", fasta, str(q)],
                check=True, stdout=fh, stderr=subprocess.DEVNULL,
            )
        sam = subprocess.run(
            ["bwa", "samse", "-n", "1000", fasta, str(sai), str(q)],
            check=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True,
        ).stdout
    cands: set[tuple[str, int]] = set()
    for line in sam.splitlines():
        if line.startswith("@"):
            continue
        f = line.split("\t")
        flag = int(f[1])
        if flag & 4:  # unmapped
            continue
        chrom, pos = f[2], int(f[3]) - 1
        cands.add((chrom, pos))
        for tag in f[11:]:
            if tag.startswith("XA:Z:"):
                for alt in tag[5:].split(";"):
                    if not alt:
                        continue
                    a = alt.split(",")
                    cands.add((a[0], abs(int(a[1])) - 1))
    return cands


def find_binding_sites_bwa(
    primer_name: str,
    primer_seq: str,
    fasta: str,
    haplotype_id: str,
    max_mismatches: int,
    *,
    slop: int = 3,
) -> list[BindingSite]:
    ensure_index(fasta)
    import pysam

    L = len(primer_seq)
    cands = _candidate_positions(primer_seq, fasta, max_mismatches)
    fa = pysam.FastaFile(fasta)
    seen: set[tuple[str, int, str]] = set()
    out: list[BindingSite] = []
    for chrom, pos in cands:
        clen = fa.get_reference_length(chrom)
        ws = max(0, pos - slop)
        we = min(clen, pos + L + slop)
        window = fa.fetch(chrom, ws, we)
        for s in find_binding_sites_naive(
            primer_name, primer_seq, window, haplotype_id, chrom, max_mismatches
        ):
            gstart, gend = ws + s.start, ws + s.end
            key = (chrom, gstart, s.strand.value)
            if key in seen:
                continue
            seen.add(key)
            s.start, s.end = gstart, gend
            out.append(s)
    fa.close()
    return out
