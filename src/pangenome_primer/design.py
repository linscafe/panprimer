"""Candidate primer generation with Primer3, steered away from variable sites.

The template is the CHM13 target sequence. Positions that are polymorphic across the
haplotype subset (see mask.py) are passed as Primer3 excluded regions so no primer is
placed over them — that is the variability-aware design that minimizes dropout by
construction (decision 5).
"""
from __future__ import annotations

from dataclasses import dataclass

from .model import Primer, PrimerPair


@dataclass
class DesignConfig:
    product_size_min: int = 100
    product_size_max: int = 300
    primer_len_min: int = 18
    primer_len_opt: int = 20
    primer_len_max: int = 27
    tm_min: float = 57.0
    tm_opt: float = 60.0
    tm_max: float = 63.0
    n_candidates: int = 20


@dataclass
class Candidate:
    pair: PrimerPair
    left_start: int   # 0-based on the template
    left_len: int
    right_start: int
    right_len: int
    penalty: float


def design_candidates(
    template_seq: str,
    excluded_regions: list[tuple[int, int]],
    cfg: DesignConfig | None = None,
    seq_id: str = "target",
) -> list[Candidate]:
    """Return up to `cfg.n_candidates` primer pairs. `excluded_regions` are (start, length)
    intervals (0-based) that primers must avoid."""
    import primer3

    cfg = cfg or DesignConfig()
    seq_args = {
        "SEQUENCE_ID": seq_id,
        "SEQUENCE_TEMPLATE": template_seq.upper(),
    }
    if excluded_regions:
        seq_args["SEQUENCE_EXCLUDED_REGION"] = [list(r) for r in excluded_regions]
    global_args = {
        "PRIMER_TASK": "generic",
        "PRIMER_PICK_LEFT_PRIMER": 1,
        "PRIMER_PICK_RIGHT_PRIMER": 1,
        "PRIMER_PICK_INTERNAL_OLIGO": 0,
        "PRIMER_OPT_SIZE": cfg.primer_len_opt,
        "PRIMER_MIN_SIZE": cfg.primer_len_min,
        "PRIMER_MAX_SIZE": cfg.primer_len_max,
        "PRIMER_OPT_TM": cfg.tm_opt,
        "PRIMER_MIN_TM": cfg.tm_min,
        "PRIMER_MAX_TM": cfg.tm_max,
        "PRIMER_PRODUCT_SIZE_RANGE": [[cfg.product_size_min, cfg.product_size_max]],
        "PRIMER_NUM_RETURN": cfg.n_candidates,
    }
    res = primer3.design_primers(seq_args, global_args)
    n = res.get("PRIMER_PAIR_NUM_RETURNED", 0)
    candidates: list[Candidate] = []
    for i in range(n):
        left_seq = res[f"PRIMER_LEFT_{i}_SEQUENCE"]
        right_seq = res[f"PRIMER_RIGHT_{i}_SEQUENCE"]
        left_start, left_len = res[f"PRIMER_LEFT_{i}"]
        right_start, right_len = res[f"PRIMER_RIGHT_{i}"]
        penalty = float(res.get(f"PRIMER_PAIR_{i}_PENALTY", 0.0))
        size = int(res[f"PRIMER_PAIR_{i}_PRODUCT_SIZE"])
        pair = PrimerPair(
            name=f"{seq_id}_pair{i}",
            forward=Primer(f"{seq_id}_pair{i}_F", left_seq),
            reverse=Primer(f"{seq_id}_pair{i}_R", right_seq),
            product_size_chm13=size,
        )
        candidates.append(
            Candidate(pair, left_start, left_len, right_start, right_len, penalty)
        )
    return candidates
