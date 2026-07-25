"""Evaluate a primer pair across haplotypes.

This module holds the in-memory (naive-backend) evaluation used for the subset path and
the correctness tests. At full scale the binding search is done by `bwa` inside a Nextflow
process and the resulting BindingSites are fed to `pcr.pair_amplicons` directly; the
classification/pairing logic is identical either way.
"""
from __future__ import annotations

from dataclasses import dataclass

from . import binding, classify, pcr
from .classify import RuleConfig, ThermoConfig
from .model import BindingSite, HaplotypeResult, HaplotypeStatus, Locus, PairResult, PrimerPair


@dataclass
class HaplotypeContext:
    """Everything needed to evaluate one haplotype in memory."""

    haplotype_id: str
    chrom: str
    sequence: str
    expected_locus: Locus | None  # None => projection failed => uncertain
    projection_ok: bool = True


@dataclass
class EvalConfig:
    mode: str = "thermo"  # or "rule"
    max_mismatches: int = 3  # search budget (competence is decided separately)
    min_product: int = 100
    max_product: int = 300
    #: Cap on amplification-competent binding sites per primer. Above it, scanning stops
    #: and the haplotype is reported as ">N binding sites" rather than enumerating
    #: products. Guards against repeat-derived primers (an Alu 20-mer has ~106k competent
    #: sites). `None` or 0 disables the cap. See `_filter_competent`.
    #: Default 100: measured competent-site counts for well-behaved primers top out at 28,
    #: so 100 clears them comfortably while sitting ~3 orders of magnitude below a repeat.
    max_binding_sites: int | None = 100
    rule_cfg: RuleConfig = None  # type: ignore[assignment]
    thermo_cfg: ThermoConfig = None  # type: ignore[assignment]

    def __post_init__(self) -> None:
        if self.rule_cfg is None:
            self.rule_cfg = RuleConfig()
        if self.thermo_cfg is None:
            self.thermo_cfg = ThermoConfig()


def _filter_competent(
    primer_seq: str,
    sites: list[BindingSite],
    window_fn,
    cfg: EvalConfig,
) -> tuple[list[BindingSite], list[str], bool]:
    """Keep amplification-competent sites; collect reasons for the rest. `window_fn(chrom,
    start, end) -> str` yields the top-strand reference the primer matched (for thermo).

    Returns `(kept, reasons, capped)`. `capped` is True when the primer exceeded
    `cfg.max_binding_sites` competent sites, at which point scanning STOPS: a primer that
    binds that many places is promiscuous regardless of what the remaining sites say, and
    scoring them is pure cost. The demo's Alu primer has ~330k raw sites (~106k competent);
    thermodynamically scoring every one of them is what made a verify run take 25 min.
    Once capped, `kept` holds exactly `max_binding_sites` sites and the true total is
    reported only as "more than the cap" -- deliberately, since we stopped counting.
    """
    kept: list[BindingSite] = []
    reasons: list[str] = []
    limit = cfg.max_binding_sites
    for s in sites:
        template_window = window_fn(s.chrom, s.start, s.end).upper()
        ok, why = classify.is_competent(
            s,
            mode=cfg.mode,
            primer_seq=primer_seq,
            template_window=template_window,
            rule_cfg=cfg.rule_cfg,
            thermo_cfg=cfg.thermo_cfg,
        )
        if ok:
            kept.append(s)
            if limit is not None and limit > 0 and len(kept) > limit:
                return kept[:limit], reasons, True
        else:
            reasons.append(f"{s.primer_name}@{s.chrom}:{s.start} {why}")
    return kept, reasons, False


def evaluate_with_sites(
    pair: PrimerPair,
    haplotype_id: str,
    f_sites: list[BindingSite],
    r_sites: list[BindingSite],
    expected_locus,
    window_fn,
    cfg: EvalConfig,
) -> HaplotypeResult:
    """Shared core: classify competence, pair into amplicons, assign status. Backend-
    agnostic — naive (in-memory) and bwa (genome-wide) paths both funnel through here."""
    if expected_locus is None:
        return HaplotypeResult(
            haplotype_id,
            HaplotypeStatus.UNCERTAIN,
            [],
            reason="locus not reliably projected onto this haplotype",
        )
    f_kept, f_reasons, f_capped = _filter_competent(
        pair.forward.sequence, f_sites, window_fn, cfg
    )
    r_kept, r_reasons, r_capped = _filter_competent(
        pair.reverse.sequence, r_sites, window_fn, cfg
    )
    if f_capped or r_capped:
        # A primer binding this many places amplifies indiscriminately; enumerating the
        # resulting product sizes would be a wall of noise, not a result. Report the fact
        # instead. `which` names the offending primer(s) so the cause is actionable.
        which = " and ".join(
            n for n, c in (("forward", f_capped), ("reverse", r_capped)) if c
        )
        return HaplotypeResult(
            haplotype_id,
            HaplotypeStatus.MULTI_PRODUCT,
            [],
            reason=(
                f">{cfg.max_binding_sites} binding sites "
                f"({which} primer); likely repeat-derived"
            ),
        )
    amplicons = pcr.pair_amplicons(
        f_kept + r_kept, expected_locus, cfg.min_product, cfg.max_product
    )
    dropout_reason = ""
    if not amplicons:
        if f_reasons or r_reasons:
            dropout_reason = "; ".join((f_reasons + r_reasons)[:3])
        elif not f_sites or not r_sites:
            dropout_reason = "a primer has no binding site in the expected locus"
        else:
            dropout_reason = "no convergent product within the product-size window"
    return pcr.classify_haplotype(
        haplotype_id, amplicons, projection_ok=True, dropout_reason=dropout_reason
    )


def evaluate_haplotype(
    pair: PrimerPair, ctx: HaplotypeContext, cfg: EvalConfig
) -> HaplotypeResult:
    """Naive in-memory path: search the whole (small) sequence, then the shared core."""
    if not ctx.projection_ok or ctx.expected_locus is None:
        return HaplotypeResult(
            ctx.haplotype_id,
            HaplotypeStatus.UNCERTAIN,
            [],
            reason="locus not reliably projected onto this haplotype",
        )
    f_sites = binding.find_binding_sites(
        pair.forward.name, pair.forward.sequence, ctx.sequence,
        ctx.haplotype_id, ctx.chrom, cfg.max_mismatches, backend="naive",
    )
    r_sites = binding.find_binding_sites(
        pair.reverse.name, pair.reverse.sequence, ctx.sequence,
        ctx.haplotype_id, ctx.chrom, cfg.max_mismatches, backend="naive",
    )
    return evaluate_with_sites(
        pair, ctx.haplotype_id, f_sites, r_sites, ctx.expected_locus,
        lambda c, s, e: ctx.sequence[s:e], cfg,
    )


def evaluate_pair(
    pair: PrimerPair,
    contexts: list[HaplotypeContext],
    cfg: EvalConfig,
    primer3_penalty: float = 0.0,
) -> PairResult:
    results = [evaluate_haplotype(pair, ctx, cfg) for ctx in contexts]
    return PairResult(pair=pair, per_haplotype=results, primer3_penalty=primer3_penalty)
