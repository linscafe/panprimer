"""Tests for the Markdown report intermediate and the optional Quarto render.

The Markdown builders operate on the plain result dicts (results_to_dict / verify_to_dict
shapes), so they need no engine objects, minimap2, or network."""
from __future__ import annotations

from pathlib import Path

import pytest

from pangenome_primer import report

_DESIGN = {
    "provenance": {"reference_build": "CHM13v2.0", "target": "chr12:1-201"},
    "pairs": [
        {
            "name": "p1", "forward": "ACGT", "reverse": "TTGC",
            "product_size_chm13": 201, "passed": True, "reject_reasons": [],
            "on_target_coverage": 1.0, "unique_product_rate": 1.0,
            "n_evaluable": 2, "n_uncertain": 0, "n_dropout": 0,
            "n_off_target": 0, "n_multi_product": 0, "primer3_penalty": 0.5,
            "per_haplotype": [
                {"haplotype_id": "h1", "status": "pass", "reason": "", "amplicons": []},
                {"haplotype_id": "h2", "status": "dropout", "reason": "3' SNP", "amplicons": []},
            ],
        },
        {
            "name": "p2", "forward": "AAAA", "reverse": "CCCC",
            "product_size_chm13": 180, "passed": False, "reject_reasons": ["off-target > 0"],
            "on_target_coverage": 1.0, "unique_product_rate": 0.5,
            "n_evaluable": 2, "n_uncertain": 0, "n_dropout": 0,
            "n_off_target": 1, "n_multi_product": 0, "primer3_penalty": 1.2,
            "per_haplotype": [
                {"haplotype_id": "h1", "status": "pass", "reason": "", "amplicons": []},
                {"haplotype_id": "h2", "status": "off_target", "reason": "paralog", "amplicons": []},
            ],
        },
    ],
}

def _attach_summaries(data: dict) -> dict:
    """Fill each row's `summary` using the production summarizer instead of hand-written
    literals, so these fixtures track verify_to_dict's shape as it changes."""
    from dataclasses import asdict

    from pangenome_primer.verify import VerifyCell, summarize

    for r in data["rows"]:
        r["summary"] = asdict(summarize([
            VerifyCell(c["haplotype_id"], c["status"], c["on_target"], c["off_target"],
                       c["size_flag"], c["reason"], c.get("site_cap"))
            for c in r["cells"]
        ]))
    return data


_VERIFY = _attach_summaries({
    "provenance": {"reference_build": "CHM13v2.0"},
    "haplotypes": ["h1", "h2"],
    "rows": [
        {
            "primer_id": "GAPDH", "forward": "A", "reverse": "C",
            "target_input": "chr12:1-201", "target_chm13": "chr12:1-201", "expected_size": 201,
            "cells": [
                {"haplotype_id": "h1", "status": "pass", "on_target": [201],
                 "off_target": [], "size_flag": False, "reason": ""},
                {"haplotype_id": "h2", "status": "off_target", "on_target": [201],
                 "off_target": [596], "size_flag": False, "reason": "paralog"},
            ],
        },
        {
            "primer_id": "DROP", "forward": "A", "reverse": "C",
            "target_input": "chr22:1-122", "target_chm13": "chr22:1-122", "expected_size": 122,
            "cells": [
                {"haplotype_id": "h1", "status": "dropout", "on_target": [],
                 "off_target": [], "size_flag": False, "reason": "3' indel"},
                {"haplotype_id": "h2", "status": "uncertain", "on_target": [],
                 "off_target": [], "size_flag": False, "reason": "unmappable"},
            ],
        },
    ],
})


def test_design_markdown_structure():
    md = report.results_to_markdown(_DESIGN)
    assert md.startswith("---\n")                       # Quarto YAML front matter
    assert "embed-resources: true" in md
    assert "# Pangenome PCR primer design" in md
    assert "[PASS]{.passed}" in md and "[reject]{.failed}" in md
    assert "[P]{.cell .pass}" in md and "[D]{.cell .dropout}" in md
    assert "**Rejected pairs:**" in md and "off-target > 0" in md
    assert "background:#fff" in md                       # forced light theme, dark text
    # provenance surfaced
    assert "`target=chr12:1-201`" in md


def test_verify_markdown_cells():
    md = report.verify_to_markdown(_VERIFY)
    assert "[201]{.ok}" in md                            # on-target size, green
    assert "off: [596]{.off}" in md                      # off-target size, red
    assert "[dropout]{.drop}" in md and "[?]{.unc}" in md
    assert " | " in md                                   # a Markdown table is present
    assert " · " in md or "off:" in md                   # cells use '·', never a bare '|'


def test_verify_cell_never_contains_bare_pipe():
    """A '|' inside a cell would corrupt the Markdown table; the separator must be '·'."""
    c = {"status": "off_target", "on_target": [201], "off_target": [596],
         "size_flag": True, "reason": ""}
    cell = report._verify_cell_md(c)
    assert "|" not in cell
    assert "{.ok .dev}" in cell                          # size_flag -> dotted underline


def _verify_data(n_haps: int, n_primers: int = 1) -> dict:
    haps = [f"h{i}" for i in range(n_haps)]
    return _attach_summaries({
        "provenance": {}, "haplotypes": haps,
        "rows": [{
            "primer_id": f"P{p}", "forward": "A", "reverse": "C",
            "target_input": "chr1:1-100", "target_chm13": "chr1:1-100", "expected_size": 100,
            "cells": [{"haplotype_id": h, "status": "pass", "on_target": [100],
                       "off_target": [], "size_flag": False, "reason": ""} for h in haps],
        } for p in range(n_primers)],
    })


def test_verify_markdown_headers_horizontal_when_few():
    # The table is transposed, so rotation keys on the PRIMER count, not the haplotype
    # count: 30 haplotypes are rows now and never widen the table.
    assert "writing-mode" not in report.verify_to_markdown(_verify_data(30, n_primers=5))


def test_verify_markdown_headers_vertical_when_many():
    # >= 6 primer pairs -> rotate the primer header columns (2nd <th> onward; column 1 is
    # the haplotype label)
    md = report.verify_to_markdown(_verify_data(3, n_primers=6))
    assert "writing-mode:vertical-rl" in md
    assert "th:nth-child(n+2)" in md


def test_verify_markdown_is_transposed():
    """Haplotypes are rows and primers are columns. This is what keeps the table readable
    as the haplotype count grows -- it is the axis that scales (30, 60, 464)."""
    md = report.verify_to_markdown(_verify_data(30, n_primers=2))
    header = next(ln for ln in md.splitlines() if ln.startswith("| haplotype"))
    assert "`P0`" in header and "`P1`" in header      # primers across the top
    assert "| h29 " in md                              # ...and every haplotype down the side
    assert "**expected bp**" in md and "**summary**" in md


def _render_verify_html(data: dict) -> str:
    from jinja2 import Environment, FileSystemLoader, select_autoescape
    tpl_dir = Path(report.__file__).resolve().parent.parent.parent / "report"
    env = Environment(loader=FileSystemLoader(str(tpl_dir)),
                      autoescape=select_autoescape(["html", "j2"]))
    return env.get_template("verify_matrix.html.j2").render(data=data)


def test_verify_html_headers_horizontal_when_few():
    html = _render_verify_html(_verify_data(30, n_primers=5))
    assert 'class="pid"' in html                 # primer headers present, not rotated
    assert 'class="pid vertical"' not in html    # the rotation modifier is absent


def test_verify_html_headers_vertical_when_many():
    # Rotation follows the primer count now: haplotypes are rows and never widen the table.
    html = _render_verify_html(_verify_data(3, n_primers=6))
    assert 'class="pid vertical"' in html


def test_verify_html_is_transposed():
    html = _render_verify_html(_verify_data(30, n_primers=2))
    # one haplotype row-label per haplotype, one primer column header per pair
    assert html.count('<th class="hap">') == 30
    assert html.count('class="pid') == 2
    assert '<th class="rowlab">expected bp</th>' in html
    assert '<th class="rowlab">summary</th>' in html


def test_quarto_render_when_available(tmp_path):
    if not report.quarto_available():
        pytest.skip("quarto not installed")
    md = tmp_path / "report.md"
    md.write_text(report.results_to_markdown(_DESIGN))
    html = report.render_quarto(str(md))
    assert html is not None and Path(html).exists()
    body = Path(html).read_text()
    assert "<table" in body                              # tables rendered
    assert "Pangenome PCR primer design" in body         # title carried through
    # embed-resources -> a single self-contained file, no sidecar _files dir
    assert not (tmp_path / "report_files").exists()
