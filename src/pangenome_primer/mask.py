"""Build a per-position variability track over the CHM13 target from the haplotype
projections, and turn it into Primer3 excluded regions (decision 5).

Each projected haplotype sequence is aligned back to the CHM13 template; substitutions and
deletions relative to the template are tallied per template position. Positions whose
alt-allele frequency across aligned haplotypes reaches `min_allele_freq` are masked so
Primer3 will not place a primer over them. This is a steering heuristic — the downstream
thermodynamic classifier still catches real dropouts on every haplotype regardless.
"""
from __future__ import annotations

import re

_CS_TOKEN = re.compile(r"(:\d+|\*[a-z][a-z]|[+\-][a-z]+|=[A-Za-z]+)")


def _accumulate_variants(cs: str, r_start: int, counts: list[int]) -> None:
    """Walk a minimap2 short `cs` string, incrementing `counts` at variable template pos."""
    rpos = r_start
    for tok in _CS_TOKEN.findall(cs):
        op = tok[0]
        if op == ":":
            rpos += int(tok[1:])
        elif op == "=":
            rpos += len(tok) - 1
        elif op == "*":  # substitution: ref base -> query base, 1 ref position
            if 0 <= rpos < len(counts):
                counts[rpos] += 1
            rpos += 1
        elif op == "-":  # deletion in query: these ref positions are missing here
            for _ in range(len(tok) - 1):
                if 0 <= rpos < len(counts):
                    counts[rpos] += 1
                rpos += 1
        elif op == "+":  # insertion in query: mark the flanking ref position
            if 0 <= rpos < len(counts):
                counts[rpos] += 1
            # ref position does not advance


def variability_counts(
    template_seq: str, haplotype_seqs: list[str], *, min_frac: float = 0.8
) -> tuple[list[int], int]:
    """Return (per-position alt counts, n_aligned haplotypes)."""
    import mappy

    counts = [0] * len(template_seq)
    aligner = mappy.Aligner(seq=template_seq.upper(), preset="asm5")
    if not aligner:
        return counts, 0
    n_aligned = 0
    for hseq in haplotype_seqs:
        if not hseq:
            continue
        best = None
        for h in aligner.map(hseq.upper(), cs=True):
            if (h.q_en - h.q_st) >= min_frac * len(hseq):
                best = h
                break
        if best is None:
            continue
        n_aligned += 1
        _accumulate_variants(best.cs, best.r_st, counts)
    return counts, n_aligned


def _merge(positions: list[int]) -> list[tuple[int, int]]:
    regions: list[tuple[int, int]] = []
    for p in sorted(positions):
        if regions and p == regions[-1][0] + regions[-1][1]:
            s, ln = regions[-1]
            regions[-1] = (s, ln + 1)
        else:
            regions.append((p, 1))
    return regions


def build_excluded_regions(
    template_seq: str,
    haplotype_seqs: list[str],
    *,
    min_allele_freq: float = 0.05,
    min_frac: float = 0.8,
) -> list[tuple[int, int]]:
    """Masked (start, length) intervals for Primer3 `SEQUENCE_EXCLUDED_REGION`."""
    counts, n = variability_counts(template_seq, haplotype_seqs, min_frac=min_frac)
    if n == 0:
        return []
    masked = [i for i, c in enumerate(counts) if c / n >= min_allele_freq]
    return _merge(masked)
