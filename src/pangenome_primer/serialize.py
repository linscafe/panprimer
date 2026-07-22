"""JSON (de)serialization for the pipeline's inter-stage artifacts. Each Nextflow process
reads and writes these so stages compose without sharing Python state."""
from __future__ import annotations

from .design import Candidate
from .model import (
    Amplicon,
    BindingSite,
    HaplotypeResult,
    HaplotypeStatus,
    Locus,
    Primer,
    PrimerPair,
    Strand,
)


def pair_to_dict(p: PrimerPair) -> dict:
    return {
        "name": p.name,
        "forward": {"name": p.forward.name, "sequence": p.forward.sequence},
        "reverse": {"name": p.reverse.name, "sequence": p.reverse.sequence},
        "product_size_chm13": p.product_size_chm13,
    }


def pair_from_dict(d: dict) -> PrimerPair:
    return PrimerPair(
        d["name"],
        Primer(d["forward"]["name"], d["forward"]["sequence"]),
        Primer(d["reverse"]["name"], d["reverse"]["sequence"]),
        product_size_chm13=d["product_size_chm13"],
    )


def candidate_to_dict(c: Candidate) -> dict:
    return {
        "pair": pair_to_dict(c.pair),
        "left_start": c.left_start,
        "left_len": c.left_len,
        "right_start": c.right_start,
        "right_len": c.right_len,
        "penalty": c.penalty,
    }


def candidate_from_dict(d: dict) -> Candidate:
    return Candidate(
        pair_from_dict(d["pair"]),
        d["left_start"], d["left_len"], d["right_start"], d["right_len"], d["penalty"],
    )


def locus_to_dict(loc: Locus | None) -> dict | None:
    if loc is None:
        return None
    return {"assembly": loc.assembly, "chrom": loc.chrom, "start": loc.start, "end": loc.end}


def locus_from_dict(d: dict | None) -> Locus | None:
    if d is None:
        return None
    return Locus(d["assembly"], d["chrom"], d["start"], d["end"])


def _amp_to_dict(a: Amplicon) -> dict:
    return {
        "haplotype_id": a.haplotype_id, "chrom": a.chrom, "start": a.start,
        "end": a.end, "size": a.size, "on_target": a.on_target,
    }


def result_to_dict(r: HaplotypeResult) -> dict:
    return {
        "haplotype_id": r.haplotype_id,
        "status": r.status.value,
        "reason": r.reason,
        "amplicons": [_amp_to_dict(a) for a in r.amplicons],
    }


def result_from_dict(d: dict) -> HaplotypeResult:
    # amplicons are re-hydrated shallowly (binding sites are not needed post-aggregation).
    amps = [
        Amplicon(a["haplotype_id"], a["chrom"], a["start"], a["end"], a["size"],
                 _stub_site(a), _stub_site(a), a["on_target"])
        for a in d["amplicons"]
    ]
    return HaplotypeResult(
        d["haplotype_id"], HaplotypeStatus(d["status"]), amps, reason=d["reason"]
    )


def _stub_site(a: dict) -> BindingSite:
    return BindingSite("", a["haplotype_id"], a["chrom"], a["start"], a["end"], Strand.PLUS, 0)
