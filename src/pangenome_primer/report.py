"""Emit results: ranked TSV + JSON (source of truth) and a self-contained HTML report
(decision 11). The JSON is authoritative and scriptable; the HTML renders the ranked pairs
and the per-pair x per-haplotype status heatmap for choosing primers by eye.
"""
from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path

from .model import Amplicon, HaplotypeStatus
from .rank import RankedPair

_STATUSES = [s.value for s in HaplotypeStatus]


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


# --- Markdown intermediate + optional Quarto render -------------------------
#
# The Markdown is a readable/diffable source that Quarto can render to a self-contained HTML
# (opt-in). Cell colours survive as Pandoc/Quarto span classes ([text]{.class}); a raw-HTML
# <style> block (passed through to HTML output) defines them. Fixed light theme — white
# background, dark text — to match the Jinja reports.

_MD_STYLE = """```{=html}
<style>
:root{--pass:#1a7f37;--dropout:#cf222e;--off_target:#bc4c00;--multi_product:#8250df;
      --uncertain:#6e7781;--ok:#1a7f37;--off:#cf222e;--dev:#bc4c00;}
body{background:#fff;color:#1f2328;}
.legend span{display:inline-block;padding:2px 8px;border-radius:3px;color:#fff;margin-right:6px;font-size:.8rem;}
span.cell{display:inline-block;min-width:1.1rem;text-align:center;color:#fff;font-weight:600;border-radius:3px;padding:1px 6px;}
.pass{background:var(--pass)}.dropout{background:var(--dropout)}.off_target{background:var(--off_target)}
.multi_product{background:var(--multi_product)}.uncertain{background:var(--uncertain)}
.passed{color:var(--pass);font-weight:700}.failed{color:var(--dropout);font-weight:700}
.ok{color:var(--ok);font-weight:700}.off{color:var(--off);font-weight:700}
.drop{color:var(--uncertain);font-style:italic}.unc{color:var(--uncertain)}
.dev{text-decoration:underline dotted var(--dev)}
</style>
```"""

_GLOSSARY_MD = """### Column & status definitions

coverage
:   **On-target coverage** — fraction of *evaluable* haplotypes that produce the intended
    amplicon (≥1 on-target band). The anti-dropout metric: does the target amplify across
    diversity?

unique
:   **Unique product rate** — fraction of evaluable haplotypes giving exactly one product and
    it is the on-target one (pass). The specificity metric: is the reaction clean?

eval
:   **Evaluable** haplotypes = total − uncertain; the denominator for all fractions.

uncertain
:   Locus broken/unmappable in the assembly; excluded from the denominator, never a failure.

dropout
:   Target should amplify but a binding-site SNP/indel prevents it (no on-target product).

off-target
:   Single product forms outside the expected homologous locus.

multi
:   More than one product (extra bands), on- and/or off-target.

penalty
:   Primer3 pair penalty (lower is better); a tie-break among passing pairs only.

product bp
:   Designed amplicon size on the CHM13 template.
"""


def _md_table(header: list[str], rows: list[list[str]]) -> str:
    def line(cells):
        return "| " + " | ".join(cells) + " |"

    out = [line(header), line(["---"] * len(header))]
    out += [line(r) for r in rows]
    return "\n".join(out)


def _provenance_md(prov: dict) -> str:
    return "**Provenance:** " + " ".join(f"`{k}={v}`" for k, v in prov.items())


def _front_matter(title: str) -> str:
    return (
        "---\n"
        f'title: "{title}"\n'
        "format:\n"
        "  html:\n"
        "    embed-resources: true\n"
        "    toc: false\n"
        "---"
    )


def results_to_markdown(data: dict) -> str:
    """Render the design results dict (see results_to_dict) as a Quarto-ready Markdown report."""
    pairs = data["pairs"]
    hap_order = (
        [h["haplotype_id"] for h in pairs[0]["per_haplotype"]] if pairs else []
    )
    status_map = {
        p["name"]: {h["haplotype_id"]: h for h in p["per_haplotype"]} for p in pairs
    }
    parts = [
        _front_matter("Pangenome PCR primer design — results"),
        _MD_STYLE,
        "",
        "# Pangenome PCR primer design — results",
        "",
        "::: {.legend}",
        " ".join(f"[{s}]{{.cell .{s}}}" for s in _STATUSES),
        ":::",
        "",
    ]
    if data.get("provenance"):
        parts += [_provenance_md(data["provenance"]), ""]
    parts += [_GLOSSARY_MD, "## Ranked primer pairs", ""]

    header = ["pair", "verdict", "coverage", "unique", "eval", "uncertain",
              "dropout", "off-target", "multi", "penalty", "product bp"]
    rows, rejects = [], []
    for p in pairs:
        rows.append([
            f"`{p['name']}`",
            "[PASS]{.passed}" if p["passed"] else "[reject]{.failed}",
            f"{p['on_target_coverage'] * 100:.1f}%",
            f"{p['unique_product_rate'] * 100:.1f}%",
            str(p["n_evaluable"]), str(p["n_uncertain"]), str(p["n_dropout"]),
            str(p["n_off_target"]), str(p["n_multi_product"]),
            f"{p['primer3_penalty']:.2f}", str(p["product_size_chm13"]),
        ])
        if not p["passed"] and p["reject_reasons"]:
            rejects.append(f"- `{p['name']}` — " + "; ".join(p["reject_reasons"]))
    parts.append(_md_table(header, rows))
    if rejects:
        parts += ["", "**Rejected pairs:**", "", *rejects]

    parts += ["", "## Per-haplotype status matrix", ""]
    mrows = []
    for p in pairs:
        cells = [f"`{p['name']}`"]
        for h in hap_order:
            st = status_map[p["name"]][h]["status"]
            cells.append(f"[{st[0].upper()}]{{.cell .{st}}}")
        mrows.append(cells)
    parts.append(_md_table(["pair \\ haplotype", *hap_order], mrows))
    parts += ["", "*Cell letter: P=pass, D=dropout, O=off_target, M=multi_product, "
              "U=uncertain.*", ""]
    return "\n".join(parts) + "\n"


def write_markdown(
    ranked: list[RankedPair], path: str, provenance: dict | None = None
) -> None:
    Path(path).write_text(results_to_markdown(results_to_dict(ranked, provenance)))


def quarto_available() -> bool:
    return shutil.which("quarto") is not None


def render_quarto(md_path: str) -> str | None:
    """Render a Markdown report to a self-contained HTML via Quarto, next to the input
    (report.md -> report.html). Returns the HTML path, or None if quarto is not installed or
    the render fails (callers fall back to the Jinja HTML)."""
    if not quarto_available():
        return None
    p = Path(md_path)
    proc = subprocess.run(
        ["quarto", "render", p.name, "--to", "html"],
        cwd=str(p.parent), capture_output=True, text=True,
    )
    html = p.with_suffix(".html")
    if proc.returncode == 0 and html.exists():
        return str(html)
    return None


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
    ranked: list[RankedPair], outdir: str, provenance: dict | None = None,
    *, quarto: bool = False, warn=None,
) -> dict[str, str]:
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    paths = {
        "json": str(out / "results.json"),
        "tsv": str(out / "results.tsv"),
        "md": str(out / "report.md"),
        "html": str(out / "report.html"),
    }
    write_json(ranked, paths["json"], provenance)
    write_tsv(ranked, paths["tsv"], provenance)
    write_markdown(ranked, paths["md"], provenance)
    rendered = render_quarto(paths["md"]) if quarto else None
    if quarto and rendered is None and warn:
        warn("quarto not found or render failed; wrote the Jinja HTML report instead")
    if rendered is None:  # default path, or Quarto unavailable/failed
        write_html(ranked, paths["html"], provenance)
    return paths


def rerender_from_json(
    json_path: str, outdir: str, rank_cfg=None, *, quarto: bool = False, warn=None
) -> dict[str, str]:
    """Rebuild PairResults from a saved results.json and re-emit all outputs with the
    current code (metrics, template). Lets a metric/template fix be applied without
    re-running the pipeline — the per-haplotype amplicons in the JSON carry everything the
    recomputed metrics need."""
    from .model import PairResult, Primer, PrimerPair
    from .rank import RankConfig, rank_pairs
    from .serialize import result_from_dict

    d = json.loads(Path(json_path).read_text())
    results = []
    for p in d["pairs"]:
        pair = PrimerPair(
            p["name"],
            Primer(f"{p['name']}_F", p["forward"]),
            Primer(f"{p['name']}_R", p["reverse"]),
            product_size_chm13=p.get("product_size_chm13", 0),
        )
        hrs = [result_from_dict(h) for h in p["per_haplotype"]]
        results.append(PairResult(pair, hrs, primer3_penalty=p.get("primer3_penalty", 0.0)))
    ranked = rank_pairs(results, rank_cfg or RankConfig())
    return write_all(ranked, outdir, provenance=d.get("provenance"), quarto=quarto, warn=warn)


# --- verify mode (screen user-supplied primers) ------------------------------

def verify_to_dict(rows, provenance: dict | None = None) -> dict:
    haplotypes = [c.haplotype_id for c in rows[0].cells] if rows else []
    return {
        "provenance": provenance or {},
        "haplotypes": haplotypes,
        "rows": [
            {
                "primer_id": r.primer_id,
                "forward": r.forward,
                "reverse": r.reverse,
                "target_input": r.target_input,
                "target_chm13": r.target_chm13,
                "expected_size": r.expected_size,
                "cells": [
                    {
                        "haplotype_id": c.haplotype_id, "status": c.status,
                        "on_target": c.on_target, "off_target": c.off_target,
                        "size_flag": c.size_flag, "reason": c.reason,
                        # Must be serialized: every renderer below (and the Jinja template)
                        # reads this dict, not the dataclass. Omitting it made capped cells
                        # fall through to the "dropout" branch -- the exact mislabel the cap
                        # exists to prevent.
                        "site_cap": c.site_cap,
                    }
                    for c in r.cells
                ],
            }
            for r in rows
        ],
    }


def _cell_text(c) -> str:
    if c["status"] == "uncertain":
        return "?"
    # Before the empty-product test: a capped cell enumerates no products by design.
    if c.get("site_cap"):
        return f">{c['site_cap']} binding sites"
    if not c["on_target"] and not c["off_target"]:
        return "dropout"
    parts = [",".join(map(str, c["on_target"]))] if c["on_target"] else []
    if c["off_target"]:
        parts.append("off:" + ",".join(map(str, c["off_target"])))
    return " | ".join(parts)


def _verify_cell_md(c: dict) -> str:
    """One verify matrix cell as Quarto Markdown: coloured product sizes via span classes."""
    if c["status"] == "uncertain":
        return "[?]{.unc}"
    if c.get("site_cap"):  # before the empty-product test; see _cell_text
        return f"[>{c['site_cap']} binding sites]{{.cap}}"
    if not c["on_target"] and not c["off_target"]:
        return "[dropout]{.drop}"
    cls = "{.ok .dev}" if c["size_flag"] else "{.ok}"
    parts = []
    if c["on_target"]:
        parts.append(", ".join(f"[{s}]{cls}" for s in c["on_target"]))
    if c["off_target"]:
        parts.append("off: " + ", ".join(f"[{s}]{{.off}}" for s in c["off_target"]))
    return " · ".join(parts)  # not '|', which would break the Markdown table


# haplotype header columns start at the 4th <th> (after primer_id, expected bp, target);
# rotate them to vertical only when there are many, to keep the table from getting too wide.
_MD_VERTICAL_HEADERS = (
    "```{=html}\n"
    "<style>thead th:nth-child(n+4){writing-mode:vertical-rl;transform:rotate(180deg);"
    "white-space:nowrap;}</style>\n"
    "```"
)


def verify_to_markdown(data: dict) -> str:
    """Render the verify dict (see verify_to_dict) as a Quarto-ready Markdown matrix.

    Haplotype headers stay horizontal when there are fewer than 6; at 6+ they rotate to
    vertical so the matrix stays narrow (mirrors the Jinja HTML report)."""
    haps = data["haplotypes"]
    parts = [
        _front_matter("Primer verification across the pangenome"),
        _MD_STYLE,
        *( [_MD_VERTICAL_HEADERS] if len(haps) >= 6 else [] ),
        "",
        "# Primer verification across the pangenome",
        "",
        "::: {.legend}",
        "[correct size]{.ok} [off-target]{.off} [dropout]{.drop} [? = not projectable]{.unc}",
        ":::",
        "",
    ]
    if data.get("provenance"):
        parts += [_provenance_md(data["provenance"]), ""]
    header = ["primer_id", "expected bp", "target (CHM13)", *haps]
    rows = []
    for r in data["rows"]:
        cellmap = {c["haplotype_id"]: c for c in r["cells"]}
        row = [f"`{r['primer_id']}`", str(r["expected_size"]), r["target_chm13"]]
        row += [_verify_cell_md(cellmap[h]) for h in haps]
        rows.append(row)
    parts.append(_md_table(header, rows))
    parts += ["", "*Cells show predicted PCR product sizes (bp). Dotted underline = size "
              "differs from expected.*", ""]
    return "\n".join(parts) + "\n"


def write_verify(
    rows, outdir: str, provenance: dict | None = None, *, quarto: bool = False, warn=None
) -> dict[str, str]:
    from jinja2 import Environment, FileSystemLoader, select_autoescape

    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    data = verify_to_dict(rows, provenance)
    paths = {
        "json": str(out / "verify.json"),
        "tsv": str(out / "verify.tsv"),
        "md": str(out / "verify_matrix.md"),
        "html": str(out / "verify_matrix.html"),
    }
    Path(paths["json"]).write_text(json.dumps(data, indent=2))

    cols = ["primer_id", "expected_size", "target_chm13", *data["haplotypes"]]
    lines = ["\t".join(cols)]
    for r in data["rows"]:
        cellmap = {c["haplotype_id"]: c for c in r["cells"]}
        row = [r["primer_id"], str(r["expected_size"]), r["target_chm13"]]
        row += [_cell_text(cellmap[h]) for h in data["haplotypes"]]
        lines.append("\t".join(row))
    Path(paths["tsv"]).write_text("\n".join(lines) + "\n")

    Path(paths["md"]).write_text(verify_to_markdown(data))
    rendered = render_quarto(paths["md"]) if quarto else None
    if quarto and rendered is None and warn:
        warn("quarto not found or render failed; wrote the Jinja HTML matrix instead")
    if rendered is None:
        tpl_dir = Path(__file__).resolve().parent.parent.parent / "report"
        env = Environment(loader=FileSystemLoader(str(tpl_dir)),
                          autoescape=select_autoescape(["html", "j2"]))
        html = env.get_template("verify_matrix.html.j2").render(data=data)
        Path(paths["html"]).write_text(html)
    return paths
