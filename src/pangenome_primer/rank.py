"""Rank primer pairs: hard filters first (every rejection has a stated reason), then a
transparent tie-break over the survivors. No opaque composite score (decision 9)."""
from __future__ import annotations

from dataclasses import dataclass

from .model import PairResult


@dataclass
class RankConfig:
    min_coverage: float = 0.95  # doc suggests >=0.99 for general assays
    max_off_target: int = 0  # off-target-product haplotypes allowed to still pass
    max_multi_product: int = 0  # multi-product haplotypes allowed to still pass
    min_evaluable: int = 1  # need at least this many non-uncertain haplotypes


@dataclass
class RankedPair:
    result: PairResult
    passed: bool
    reject_reasons: list[str]


def _reasons(r: PairResult, cfg: RankConfig) -> list[str]:
    reasons: list[str] = []
    if len(r.evaluable) < cfg.min_evaluable:
        reasons.append(
            f"only {len(r.evaluable)} evaluable haplotypes "
            f"(< {cfg.min_evaluable}); {r.n_uncertain} uncertain"
        )
    if r.on_target_coverage < cfg.min_coverage:
        reasons.append(
            f"coverage {r.on_target_coverage:.2%} < {cfg.min_coverage:.2%} "
            f"({r.n_dropout} dropout)"
        )
    if r.n_off_target > cfg.max_off_target:
        reasons.append(f"{r.n_off_target} off-target haplotypes > {cfg.max_off_target}")
    if r.n_multi_product > cfg.max_multi_product:
        reasons.append(
            f"{r.n_multi_product} multi-product haplotypes > {cfg.max_multi_product}"
        )
    return reasons


def _tie_break_key(r: PairResult):
    # higher coverage, higher uniqueness, fewer off-target, fewer multi, lower penalty,
    # then name for determinism. Negate the "higher is better" terms for ascending sort.
    return (
        -r.on_target_coverage,
        -r.unique_product_rate,
        r.n_off_target,
        r.n_multi_product,
        r.primer3_penalty,
        r.pair.name,
    )


def rank_pairs(results: list[PairResult], cfg: RankConfig | None = None) -> list[RankedPair]:
    cfg = cfg or RankConfig()
    ranked = [RankedPair(r, not (rs := _reasons(r, cfg)), rs) for r in results]
    # survivors first (sorted by tie-break), then the rejected (also stably ordered).
    survivors = sorted((rp for rp in ranked if rp.passed), key=lambda rp: _tie_break_key(rp.result))
    rejected = sorted((rp for rp in ranked if not rp.passed), key=lambda rp: _tie_break_key(rp.result))
    return survivors + rejected
