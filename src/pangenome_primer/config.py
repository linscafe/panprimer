"""Load run parameters from config/defaults.yaml (+ optional overrides) into the typed
config objects each stage uses."""
from __future__ import annotations

from pathlib import Path

import yaml

from .classify import RuleConfig, ThermoConfig
from .design import DesignConfig
from .rank import RankConfig

_DEFAULTS = Path(__file__).resolve().parent.parent.parent / "config" / "defaults.yaml"


def load_raw(path: str | None = None) -> dict:
    return yaml.safe_load(Path(path or _DEFAULTS).read_text())


def design_config(raw: dict) -> DesignConfig:
    d = raw["design"]
    return DesignConfig(
        product_size_min=d["product_size_min"],
        product_size_max=d["product_size_max"],
        primer_len_min=d["primer_len_min"],
        primer_len_opt=d["primer_len_opt"],
        primer_len_max=d["primer_len_max"],
        tm_min=d["primer_tm_min"],
        tm_opt=d["primer_tm_opt"],
        tm_max=d["primer_tm_max"],
        n_candidates=d["n_candidates"],
    )


def thermo_config(raw: dict, design_tm: float) -> ThermoConfig:
    t = raw["dropout"]["thermo"]
    return ThermoConfig(
        design_tm=design_tm,
        tm_drop_max=t["tm_drop_max"],
        three_prime_hard_nt=t["three_prime_hard_nt"],
        mv_conc=t["mv_conc"],
        dv_conc=t["dv_conc"],
        dntp_conc=t["dntp_conc"],
        dna_conc=t["dna_conc"],
    )


def rule_config(raw: dict) -> RuleConfig:
    r = raw["dropout"]["rule"]
    return RuleConfig(
        max_total_mm=r["max_total_mm"],
        max_3prime_mm=r["max_3prime_mm"],
        suffix=r["suffix"],
    )


def rank_config(raw: dict) -> RankConfig:
    r = raw["rank"]
    return RankConfig(
        min_coverage=r["min_coverage"],
        max_off_target=r["max_off_target"],
        max_multi_product=r["max_multi_product"],
        min_evaluable=r["min_evaluable"],
    )
