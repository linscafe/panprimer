"""Decide whether a primer binding site is amplification-competent.

Two modes (decision 10):

* **thermodynamic** (primary) — score the primer/template duplex, including mismatches
  and 3'-end stability, with primer3's thermodynamic alignment. A site is competent when
  its predicted annealing Tm is within `tm_drop_max` of the design Tm and its 3'-end
  duplex ΔG is stable enough. Requires `primer3` (conda env).
* **rule** (option) — pure-Python gates on mismatch count / 3'-window / indels. Needs no
  deps, so it is the backbone the correctness tests exercise.

Both return a bool (competent) plus a short reason string used in the report.
"""
from __future__ import annotations

from dataclasses import dataclass

from .binding import revcomp
from .model import BindingSite, Strand


@dataclass
class RuleConfig:
    max_total_mm: int = 2
    max_3prime_mm: int = 0
    suffix: int = 5  # size of the 3' window that must be (near) clean


@dataclass
class ThermoConfig:
    design_tm: float = 60.0
    tm_drop_max: float = 5.0  # competent if site Tm >= design_tm - tm_drop_max
    three_prime_hard_nt: int = 2  # a mismatch within this terminal window is a hard fail
    mv_conc: float = 50.0
    dv_conc: float = 1.5
    dntp_conc: float = 0.6
    dna_conc: float = 250.0


def competent_rule(site: BindingSite, cfg: RuleConfig) -> tuple[bool, str]:
    if site.has_indel:
        return False, "indel in binding site"
    if site.mismatches > cfg.max_total_mm:
        return False, f"{site.mismatches} mismatches > max {cfg.max_total_mm}"
    n_3prime = sum(1 for o in site.mismatch_offsets_3p if o < cfg.suffix)
    if n_3prime > cfg.max_3prime_mm:
        return False, (
            f"{n_3prime} mismatch(es) within 3' {cfg.suffix}nt "
            f"> max {cfg.max_3prime_mm}"
        )
    return True, "ok"


def competent_thermo(
    primer_seq: str,
    template_window: str,
    cfg: ThermoConfig,
    strand: Strand = Strand.PLUS,
) -> tuple[bool, str]:
    """`template_window` is the top-strand reference the primer aligned to (same length as
    the primer). The primer is identical to `window` on a PLUS site and to `revcomp(window)`
    on a MINUS site; either way it hybridizes to the reverse complement of the strand it is
    identical to. Orienting by `strand` is essential — a perfect minus-strand primer scored
    against the wrong strand yields a nonsense ~4 C Tm."""
    import primer3  # lazy: only when thermo mode is actually used

    primer_aligned_ref = template_window if strand is Strand.PLUS else revcomp(template_window)
    target_strand = revcomp(primer_aligned_ref)
    res = primer3.calc_heterodimer(
        primer_seq,
        target_strand,
        mv_conc=cfg.mv_conc,
        dv_conc=cfg.dv_conc,
        dntp_conc=cfg.dntp_conc,
        dna_conc=cfg.dna_conc,
        output_structure=False,
    )
    site_tm = res.tm if res.structure_found else -273.0
    if site_tm < cfg.design_tm - cfg.tm_drop_max:
        return False, (
            f"annealing Tm {site_tm:.1f}C < {cfg.design_tm - cfg.tm_drop_max:.1f}C"
        )
    # 3'-end stability: an unstable 3' end fails to extend even at adequate overall Tm.
    end = primer3.calc_end_stability(
        primer_seq,
        target_strand,
        mv_conc=cfg.mv_conc,
        dv_conc=cfg.dv_conc,
        dntp_conc=cfg.dntp_conc,
        dna_conc=cfg.dna_conc,
    )
    if end.tm < cfg.design_tm - cfg.tm_drop_max - 10.0:
        return False, f"3' end unstable (Tm {end.tm:.1f}C)"
    return True, f"annealing Tm {site_tm:.1f}C"


def is_competent(
    site: BindingSite,
    *,
    mode: str = "thermo",
    primer_seq: str | None = None,
    template_window: str | None = None,
    rule_cfg: RuleConfig | None = None,
    thermo_cfg: ThermoConfig | None = None,
) -> tuple[bool, str]:
    if mode == "rule":
        return competent_rule(site, rule_cfg or RuleConfig())
    if mode == "thermo":
        if primer_seq is None or template_window is None:
            raise ValueError("thermo mode requires primer_seq and template_window")
        tcfg = thermo_cfg or ThermoConfig()
        # Non-negotiable 3'-terminal gate: whole-duplex Tm is nearly blind to a terminal
        # mismatch, but polymerase extension is not. Reuse the per-strand 3' offsets the
        # search already computed (verified logic) before trusting the Tm model.
        if site.has_indel:
            return False, "indel in binding site"
        if any(o < tcfg.three_prime_hard_nt for o in site.mismatch_offsets_3p):
            return False, f"mismatch within 3' terminal {tcfg.three_prime_hard_nt}nt"
        return competent_thermo(primer_seq, template_window, tcfg, strand=site.strand)
    raise ValueError(f"unknown classify mode: {mode!r}")
