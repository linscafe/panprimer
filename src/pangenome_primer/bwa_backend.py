"""Genome-wide primer binding search on a full haplotype assembly.

`bwa aln`/`samse` locates candidate binding positions across a multi-Gb assembly quickly;
we then extract each candidate window and recompute the exact mismatches and 3'-end offsets
with the naive comparator from `binding.py` (which the tests verify). bwa prunes the search
space; the trusted code decides the details. Requires `bwa` + `samtools` + `pysam`.
"""
from __future__ import annotations

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


def _candidate_positions_batch(
    seqs: list[str], fasta: str, max_mm: int
) -> dict[int, set[tuple[str, int]]]:
    """Run ONE bwa aln/samse over all primer sequences (index loaded once). Returns
    {query_index -> {(chrom, 0-based pos)}} from primary + XA alternative hits."""
    with tempfile.TemporaryDirectory() as td:
        q = Path(td) / "q.fa"
        # bwa seeding needs -l/-k sized to the shortest primer
        minlen = min(len(s) for s in seqs)
        q.write_text("".join(f">{i}\n{s}\n" for i, s in enumerate(seqs)))
        sai = Path(td) / "q.sai"
        with open(sai, "wb") as fh:
            subprocess.run(
                ["bwa", "aln", "-n", str(max_mm), "-o", "0", "-l", str(minlen),
                 "-k", str(max_mm), "-N", fasta, str(q)],
                check=True, stdout=fh, stderr=subprocess.DEVNULL,
            )
        sam = subprocess.run(
            ["bwa", "samse", "-n", "1000", fasta, str(sai), str(q)],
            check=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True,
        ).stdout
    cands: dict[int, set[tuple[str, int]]] = {i: set() for i in range(len(seqs))}
    for line in sam.splitlines():
        if line.startswith("@"):
            continue
        f = line.split("\t")
        if int(f[1]) & 4:  # unmapped
            continue
        qi = int(f[0])
        cands[qi].add((f[2], int(f[3]) - 1))
        for tag in f[11:]:
            if tag.startswith("XA:Z:"):
                for alt in tag[5:].split(";"):
                    if alt:
                        a = alt.split(",")
                        cands[qi].add((a[0], abs(int(a[1])) - 1))
    return cands


def find_binding_sites_batch(
    seqs: list[str],
    fasta: str,
    haplotype_id: str,
    max_mismatches: int,
    *,
    slop: int = 3,
) -> dict[str, list[BindingSite]]:
    """Genome-wide binding sites for many primers in one index load. Deduplicates identical
    sequences (candidate pairs share primers). Returns {sequence -> [BindingSite]}. Each
    candidate position's exact mismatches/3' offsets are recomputed with the naive
    comparator, so bwa only prunes the search space."""
    ensure_index(fasta)
    import pysam

    uniq = sorted(set(seqs))
    cands = _candidate_positions_batch(uniq, fasta, max_mismatches)
    fa = pysam.FastaFile(fasta)
    out: dict[str, list[BindingSite]] = {}
    for i, primer_seq in enumerate(uniq):
        L = len(primer_seq)
        seen: set[tuple[str, int, str]] = set()
        sites: list[BindingSite] = []
        for chrom, pos in cands[i]:
            clen = fa.get_reference_length(chrom)
            ws, we = max(0, pos - slop), min(clen, pos + L + slop)
            window = fa.fetch(chrom, ws, we)
            for s in find_binding_sites_naive(
                primer_seq, primer_seq, window, haplotype_id, chrom, max_mismatches
            ):
                gstart, gend = ws + s.start, ws + s.end
                key = (chrom, gstart, s.strand.value)
                if key in seen:
                    continue
                seen.add(key)
                s.start, s.end = gstart, gend
                sites.append(s)
        out[primer_seq] = sites
    fa.close()
    return out
