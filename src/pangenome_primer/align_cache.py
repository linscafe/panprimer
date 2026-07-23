"""Whole-genome CHM13<->haplotype alignment cache for fast locus projection (item 1).

Building a per-haplotype PAF once (in prep) turns per-query projection from a multi-GB
minimap2 index build into a cheap coordinate lift: find the alignment block(s) covering the
CHM13 target interval and map to haplotype coordinates. No per-query alignment, no index load.

PAF orientation: we align `minimap2 chm13 hap` so the *target* columns are CHM13 and the
*query* columns are the haplotype. Given a CHM13 interval we lift target->query.
"""
from __future__ import annotations

import subprocess
from dataclasses import dataclass
from pathlib import Path

from .binding import revcomp
from .model import Locus
from .project import Projection


def paf_path(hap_fasta: str) -> str:
    return hap_fasta + ".chm13.paf"


def build_alignment(chm13_fasta: str, hap_fasta: str, out_paf: str | None = None) -> str:
    """minimap2 asm5 alignment of a haplotype to CHM13, written as PAF (cached once)."""
    out = out_paf or paf_path(hap_fasta)
    tmp = out + ".part"
    with open(tmp, "wb") as fh:
        subprocess.run(
            ["minimap2", "-cx", "asm5", "--secondary=no", chm13_fasta, hap_fasta],
            check=True, stdout=fh, stderr=subprocess.DEVNULL,
        )
    Path(tmp).replace(out)
    return out


@dataclass
class _Block:
    q_name: str
    q_start: int
    q_end: int
    strand: str
    t_name: str
    t_start: int
    t_end: int


def _overlapping_blocks(paf: str, chrom: str, start: int, end: int) -> list[_Block]:
    blocks: list[_Block] = []
    with open(paf) as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) < 12 or f[5] != chrom:
                continue
            t_start, t_end = int(f[7]), int(f[8])
            if t_start < end and start < t_end:  # overlaps the CHM13 target
                blocks.append(
                    _Block(f[0], int(f[2]), int(f[3]), f[4], f[5], t_start, t_end)
                )
    return blocks


def _lift(b: _Block, t: int) -> int:
    """Linearly map a CHM13 (target) coordinate into haplotype (query) coordinates within a
    collinear asm5 block. Exact to within local indel sizes, which the flank absorbs."""
    span_t = max(1, b.t_end - b.t_start)
    frac = (min(max(t, b.t_start), b.t_end) - b.t_start) / span_t
    if b.strand == "+":
        return round(b.q_start + frac * (b.q_end - b.q_start))
    return round(b.q_end - frac * (b.q_end - b.q_start))


def project_from_paf(
    paf: str, chrom: str, start: int, end: int, hap_fasta: str
) -> Projection:
    """Lift a CHM13 target interval to a haplotype window using the cached PAF."""
    import pysam

    blocks = _overlapping_blocks(paf, chrom, start, end)
    if not blocks:
        return Projection(None, reason="no cached alignment block covers the target")
    # best = block with the largest overlap of the target interval
    best = max(blocks, key=lambda b: min(end, b.t_end) - max(start, b.t_start))
    a, c = _lift(best, start), _lift(best, end)
    q_lo, q_hi = min(a, c), max(a, c)
    fa = pysam.FastaFile(hap_fasta)
    clen = fa.get_reference_length(best.q_name)
    q_lo, q_hi = max(0, q_lo), min(clen, q_hi)
    seq = fa.fetch(best.q_name, q_lo, q_hi)
    fa.close()
    if best.strand == "-":
        # extracted haplotype window is antiparallel to the CHM13 template; orient to match
        seq = revcomp(seq)
    locus = Locus(hap_fasta, best.q_name, q_lo, q_hi)
    return Projection(locus, haplotype_seq=seq, reason="paf-cache")
