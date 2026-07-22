"""In-silico PCR: pair competent binding sites into amplicons, then assign a per-haplotype
status. This is where "does it amplify, on target, uniquely?" is decided.

Inputs are already-competent binding sites (dropout at the site level was decided in
classify.py). An amplicon forms between a PLUS-strand site (3' pointing right) and a
downstream MINUS-strand site (3' pointing left) whose product size falls in range. Either
primer of the pair may act as the plus-strand primer at a given locus.
"""
from __future__ import annotations

from .model import Amplicon, BindingSite, HaplotypeResult, HaplotypeStatus, Locus, Strand


def pair_amplicons(
    sites: list[BindingSite],
    expected_locus: Locus | None,
    min_size: int,
    max_size: int,
) -> list[Amplicon]:
    plus = [s for s in sites if s.strand is Strand.PLUS]
    minus = [s for s in sites if s.strand is Strand.MINUS]
    amplicons: list[Amplicon] = []
    for p in plus:
        for m in minus:
            if m.chrom != p.chrom:
                continue
            # convergent: plus primer to the left, minus primer to the right,
            # 3' ends facing each other.
            if p.start > m.start:
                continue
            size = m.end - p.start
            if size < min_size or size > max_size:
                continue
            on_target = expected_locus is not None and _contained(
                p.chrom, p.start, m.end, expected_locus
            )
            amplicons.append(
                Amplicon(
                    haplotype_id=p.haplotype_id,
                    chrom=p.chrom,
                    start=p.start,
                    end=m.end,
                    size=size,
                    forward=p,
                    reverse=m,
                    on_target=on_target,
                )
            )
    return amplicons


def _contained(chrom: str, start: int, end: int, locus: Locus) -> bool:
    return chrom == locus.chrom and start >= locus.start and end <= locus.end


def classify_haplotype(
    haplotype_id: str,
    amplicons: list[Amplicon],
    *,
    projection_ok: bool,
    dropout_reason: str = "",
) -> HaplotypeResult:
    """Assign one status. Precedence follows CONTEXT.md:

    uncertain (no reliable projection) > multi_product (>1 product) >
    off_target (1 product, off target) > pass (1 product, on target) > dropout (0).
    """
    if not projection_ok:
        return HaplotypeResult(
            haplotype_id,
            HaplotypeStatus.UNCERTAIN,
            amplicons,
            reason="locus not reliably projected onto this haplotype",
        )
    if len(amplicons) == 0:
        return HaplotypeResult(
            haplotype_id,
            HaplotypeStatus.DROPOUT,
            [],
            reason=dropout_reason or "no amplification-competent product",
        )
    if len(amplicons) > 1:
        n_on = sum(1 for a in amplicons if a.on_target)
        return HaplotypeResult(
            haplotype_id,
            HaplotypeStatus.MULTI_PRODUCT,
            amplicons,
            reason=f"{len(amplicons)} products ({n_on} on-target)",
        )
    (only,) = amplicons
    if only.on_target:
        return HaplotypeResult(
            haplotype_id,
            HaplotypeStatus.PASS,
            amplicons,
            reason=f"single on-target product {only.size} bp",
        )
    return HaplotypeResult(
        haplotype_id,
        HaplotypeStatus.OFF_TARGET,
        amplicons,
        reason=f"single off-target product {only.size} bp",
    )
