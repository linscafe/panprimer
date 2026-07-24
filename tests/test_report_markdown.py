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

_VERIFY = {
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
}


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
