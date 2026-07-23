"""Anchor inputs on CHM13 and project the target locus onto each haplotype.

Both operations are minimap2 alignments (via `mappy`, lazy-imported so the package imports
without it):

* `anchor_sequence` — a GRCh38-coord slice or a raw FASTA locus is aligned to CHM13 to get
  its CHM13 anchor. Multiple strong hits (segmental dup / gene family) => ambiguous => we
  fail loud and ask the user to disambiguate (decision 2), never silently pick one.
* `project_locus` — the CHM13 target ± flank is aligned to a haplotype assembly to find the
  expected homologous locus. No confident hit => projection fails => the haplotype is
  scored `uncertain`, not `dropout`.
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

from .model import Locus

_COORD_RE = re.compile(r"^([^:\s]+):([\d,]+)-([\d,]+)$")


def parse_coords(target: str) -> tuple[str, int, int] | None:
    """Parse 'chr:start-end' (1-based inclusive-looking, treated as 0-based half-open here)
    into (chrom, start, end), or None if it is not a coordinate string."""
    m = _COORD_RE.match(target.strip())
    if not m:
        return None
    chrom, s, e = m.group(1), int(m.group(2).replace(",", "")), int(m.group(3).replace(",", ""))
    return chrom, s, e


@dataclass
class Projection:
    """Result of projecting the CHM13 target onto one haplotype."""

    locus: Locus | None       # None => projection failed (uncertain)
    haplotype_seq: str = ""   # aligned haplotype subsequence, for the variability mask
    reason: str = ""


class AmbiguousAnchor(Exception):
    """Input maps to more than one place on CHM13; the user must disambiguate."""


def _aligner(fasta: str, preset: str):
    """Build a mappy aligner, preferring a prebuilt `<fasta>.mmi` if present. Building a
    full-genome minimap2 index in-process costs minutes and several GB; a cached .mmi (see
    scripts/prepare_haplotypes.sh) loads in seconds at the same memory footprint."""
    import mappy

    from pathlib import Path

    mmi = fasta + ".mmi"
    if Path(mmi).exists():
        return mappy.Aligner(fn_idx_in=mmi, preset=preset)
    return mappy.Aligner(fasta, preset=preset)


def _confident_hits(aligner, seq: str, min_frac: float, min_mapq: int):
    hits = [
        h
        for h in aligner.map(seq)
        if h.mapq >= min_mapq and (h.q_en - h.q_st) >= min_frac * len(seq)
    ]
    hits.sort(key=lambda h: (h.q_en - h.q_st), reverse=True)
    return hits


def anchor_sequence(
    query_seq: str,
    chm13_fasta: str,
    *,
    min_frac: float = 0.9,
    min_mapq: int = 30,
) -> Locus:
    """Locate `query_seq` on CHM13 -> anchor Locus. Raises AmbiguousAnchor on ties."""
    aligner = _aligner(chm13_fasta, preset="asm5")
    if not aligner:
        raise RuntimeError(f"failed to index {chm13_fasta}")
    hits = _confident_hits(aligner, query_seq, min_frac, min_mapq)
    if not hits:
        raise ValueError("input did not map confidently to CHM13")
    best = hits[0]
    # a near-equal second hit => ambiguous locus
    strong = [h for h in hits if (h.q_en - h.q_st) >= 0.95 * (best.q_en - best.q_st)]
    if len(strong) > 1:
        locs = ", ".join(f"{h.ctg}:{h.r_st}-{h.r_en}" for h in strong[:5])
        raise AmbiguousAnchor(
            f"input maps to {len(strong)} CHM13 loci ({locs}); "
            "specify CHM13 coordinates to disambiguate"
        )
    return Locus("CHM13v2.0", best.ctg, best.r_st, best.r_en)


def resolve_target(
    target: str,
    chm13_fasta: str,
    *,
    assembly: str = "chm13",
    grch38_fasta: str | None = None,
) -> tuple[str, int, int]:
    """Normalize any accepted target to a CHM13 (chrom, start, end):

    * `chr:start-end` with assembly=chm13  -> used directly (already anchored).
    * `chr:start-end` with assembly=grch38 -> slice extracted from `grch38_fasta`, then
      aligned to CHM13 (same path as FASTA input).
    * a FASTA path                          -> its sequence aligned to CHM13.

    Raises AmbiguousAnchor when a sequence maps to multiple CHM13 loci (decision 2)."""
    coords = parse_coords(target)
    if coords is not None:
        chrom, start, end = coords
        if assembly == "chm13":
            return chrom, start, end
        if assembly == "grch38":
            if not grch38_fasta:
                raise ValueError("GRCh38 coordinate input requires --grch38 <GRCh38 FASTA>")
            import pysam

            fa = pysam.FastaFile(grch38_fasta)
            if chrom not in fa.references:
                fa.close()
                raise ValueError(
                    f"{chrom} not in {grch38_fasta} (naming mismatch? e.g. 'chr1' vs '1')"
                )
            seq = fa.fetch(chrom, start, end)
            fa.close()
            loc = anchor_sequence(seq, chm13_fasta)
            return loc.chrom, loc.start, loc.end
        raise ValueError(f"unknown target assembly: {assembly!r}")
    # FASTA path
    if not Path(target).exists():
        raise FileNotFoundError(f"target is neither chr:start-end coords nor a file: {target}")
    seq = "".join(l for l in Path(target).read_text().splitlines() if not l.startswith(">"))
    loc = anchor_sequence(seq, chm13_fasta)
    return loc.chrom, loc.start, loc.end


def project_target(
    chrom: str,
    tstart: int,
    tend: int,
    template_seq: str,
    haplotype_fasta: str,
) -> Projection:
    """Project a CHM13 target interval onto a haplotype. Prefers the whole-genome PAF cache
    (a fast coordinate lift; see align_cache) and falls back to on-the-fly alignment of the
    template sequence when no cache exists."""
    from pathlib import Path

    from . import align_cache

    paf = align_cache.paf_path(haplotype_fasta)
    if Path(paf).exists():
        return align_cache.project_from_paf(paf, chrom, tstart, tend, haplotype_fasta)
    return project_locus(template_seq, haplotype_fasta)


def project_locus(
    target_seq: str,
    haplotype_fasta: str,
    *,
    min_frac: float = 0.8,
    min_mapq: int = 20,
) -> Projection:
    """Project the CHM13 target (± flank) sequence onto one haplotype assembly."""
    aligner = _aligner(haplotype_fasta, preset="asm5")
    if not aligner:
        return Projection(None, reason=f"failed to index {haplotype_fasta}")
    hits = _confident_hits(aligner, target_seq, min_frac, min_mapq)
    if not hits:
        return Projection(None, reason="target did not map confidently to this haplotype")
    best = hits[0]
    locus = Locus(haplotype_fasta, best.ctg, best.r_st, best.r_en)
    # pull the haplotype subsequence for the variability mask
    seq = ""
    try:
        seq = aligner.seq(best.ctg, best.r_st, best.r_en) or ""
    except Exception:
        seq = ""
    return Projection(locus, haplotype_seq=seq, reason="ok")
