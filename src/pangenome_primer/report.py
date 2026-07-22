"""Emit results: ranked TSV + JSON (source of truth) and a self-contained HTML report
(decision 11). The JSON is authoritative and scriptable; the HTML renders the ranked pairs
and the per-pair x per-haplotype status heatmap for choosing primers by eye.
"""
from __future__ import annotations

import json
from dataclasses import asdict, is_dataclass
from pathlib import Path

from .model import Amplicon, HaplotypeStatus
from .rank import RankedPair


def _amplicon_dict(a: Amplicon) -> dict:
    return {
        "chrom": a.chrom,
        "start": a.start,
        "end": a.end,
        "size": a.size,
        "on_target": a.on_target,
    }


def results_to_dict(ranked: list[RankedPair], provenance: dict | None = None) -> dict:
    pairs = []
    for rp in ranked:
        r = rp.result
        pairs.append(
            {
                "name": r.pair.name,
                "forward": r.pair.forward.sequence,
                "reverse": r.pair.reverse.sequence,
                "product_size_chm13": r.pair.product_size_chm13,
                "passed": rp.passed,
                "reject_reasons": rp.reject_reasons,
                "on_target_coverage": round(r.on_target_coverage, 4),
                "unique_product_rate": round(r.unique_product_rate, 4),
                "n_evaluable": len(r.evaluable),
                "n_uncertain": r.n_uncertain,
                "n_dropout": r.n_dropout,
                "n_off_target": r.n_off_target,
                "n_multi_product": r.n_multi_product,
                "primer3_penalty": r.primer3_penalty,
                "per_haplotype": [
                    {
                        "haplotype_id": hr.haplotype_id,
                        "status": hr.status.value,
                        "reason": hr.reason,
                        "amplicons": [_amplicon_dict(a) for a in hr.amplicons],
                    }
                    for hr in r.per_haplotype
                ],
            }
        )
    return {"provenance": provenance or {}, "pairs": pairs}


def write_json(ranked: list[RankedPair], path: str, provenance: dict | None = None) -> None:
    Path(path).write_text(json.dumps(results_to_dict(ranked, provenance), indent=2))


_TSV_COLS = [
    "name", "passed", "on_target_coverage", "unique_product_rate", "n_evaluable",
    "n_uncertain", "n_dropout", "n_off_target", "n_multi_product",
    "primer3_penalty", "forward", "reverse", "product_size_chm13", "reject_reasons",
]


def write_tsv(ranked: list[RankedPair], path: str, provenance: dict | None = None) -> None:
    d = results_to_dict(ranked, provenance)
    lines = ["\t".join(_TSV_COLS)]
    for p in d["pairs"]:
        row = dict(p)
        row["reject_reasons"] = "; ".join(p["reject_reasons"])
        lines.append("\t".join(str(row[c]) for c in _TSV_COLS))
    Path(path).write_text("\n".join(lines) + "\n")


def write_html(
    ranked: list[RankedPair], path: str, provenance: dict | None = None
) -> None:
    from jinja2 import Environment, FileSystemLoader, select_autoescape

    tpl_dir = Path(__file__).resolve().parent.parent.parent / "report"
    env = Environment(
        loader=FileSystemLoader(str(tpl_dir)),
        autoescape=select_autoescape(["html", "j2"]),
    )
    data = results_to_dict(ranked, provenance)
    # stable haplotype column order from the first pair
    hap_order = (
        [h["haplotype_id"] for h in data["pairs"][0]["per_haplotype"]]
        if data["pairs"]
        else []
    )
    status_map = {
        p["name"]: {h["haplotype_id"]: h for h in p["per_haplotype"]}
        for p in data["pairs"]
    }
    html = env.get_template("template.html.j2").render(
        data=data,
        hap_order=hap_order,
        status_map=status_map,
        statuses=[s.value for s in HaplotypeStatus],
    )
    Path(path).write_text(html)


def write_all(
    ranked: list[RankedPair], outdir: str, provenance: dict | None = None
) -> dict[str, str]:
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    paths = {
        "json": str(out / "results.json"),
        "tsv": str(out / "results.tsv"),
        "html": str(out / "report.html"),
    }
    write_json(ranked, paths["json"], provenance)
    write_tsv(ranked, paths["tsv"], provenance)
    write_html(ranked, paths["html"], provenance)
    return paths


def _default(o):  # pragma: no cover - json fallback for stray dataclasses
    if is_dataclass(o):
        return asdict(o)
    raise TypeError(o)
