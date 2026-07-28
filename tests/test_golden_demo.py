"""Golden-output harness (Phase 0 safety net).

`tests/golden/verify.json` and `tests/golden/results.json` are frozen copies of
`demo-verify-pipeline/verify.json` and `demo-design-pipeline/gapdh/results.json` -- the
existing, human-reviewed demo output. `golden_compare.py` diffs a fresh run's output
against them cell for cell. This is the gate every backend and storage change must clear:
the regenerated matrix must reproduce the frozen one exactly, including:

* `CYP2D6_paralog` in verify.json -- the CYP2D7-off-target row (status `multi_product`,
  an off-target product alongside the correct one on every haplotype).
* `CYP2D6_dropout` in verify.json -- the allele-dropout row. It is HAPLOTYPE-DIFFERENTIAL:
  `dropout` in HG00097#hap1 and `pass` in HG01884#hap1 / HG00408#hap1, driven by rs1058164
  (CYP2D6 c.1661G>C, exon 3) sitting under the forward primer's 3' terminal base. Do not
  weaken this to "dropout everywhere" -- see the note on that assertion for why a
  fails-everywhere row is indistinguishable from a broken primer.

  These three carry AFR/EUR/EAS labels, but the split is NOT population-structured: rs1058164
  is globally common, and across 30 haplotypes the same pair drops out in 11/24 evaluable in
  every superpopulation (AFR 2/5, AMR 3/6, EAS 1/3, EUR 3/5, SAS 2/5). The fixture picks up a
  genotype difference between three individuals, not a population signal -- so do not add
  assertions tying a superpopulation label to an expected status.

Two speeds of test live here:

* Fast (default `pytest -q`): the comparator's own logic against synthetic dicts, plus a
  sanity check that the frozen golden fixtures themselves still contain the off-target and
  dropout rows in the shape the comparator expects, plus (if present) a diff of the frozen
  golden against the live, gitignored `demo-*/*.json` outputs in this checkout.
* Slow (`pytest -m slow --runslow`, ~9 min, needs real HPRC genome data): actually invokes
  the CLI to regenerate verify.json/results.json from scratch and diffs the result against
  the frozen golden.
"""
from __future__ import annotations

import json
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent))
from golden_compare import compare_results_dicts, compare_verify_dicts  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parent.parent
GOLDEN_VERIFY = Path(__file__).parent / "golden" / "verify.json"
GOLDEN_RESULTS = Path(__file__).parent / "golden" / "results.json"


def _load(path: Path) -> dict:
    return json.loads(path.read_text())


# =========================================================================================
# Fast: comparator logic against synthetic dicts (no genome, no I/O beyond these dicts)
# =========================================================================================


def _verify_dict(*, status="pass", on_target=None, off_target=None, size_flag=False):
    return {
        "provenance": {},
        "haplotypes": ["hapA"],
        "rows": [
            {
                "primer_id": "p1",
                "forward": "ACGT",
                "reverse": "TGCA",
                "target_input": "chr1:1-100",
                "target_chm13": "chr1:1-100",
                "expected_size": 100,
                "cells": [
                    {
                        "haplotype_id": "hapA",
                        "status": status,
                        "on_target": on_target or [],
                        "off_target": off_target or [],
                        "size_flag": size_flag,
                        "reason": "x",
                    }
                ],
            }
        ],
    }


def test_compare_verify_dicts_identical_is_no_diff():
    d = _verify_dict(status="pass", on_target=[100])
    assert compare_verify_dicts(d, json.loads(json.dumps(d))) == []


def test_compare_verify_dicts_catches_status_mismatch_and_names_the_cell():
    golden = _verify_dict(status="pass", on_target=[100])
    candidate = _verify_dict(status="dropout", on_target=[])
    diffs = compare_verify_dicts(golden, candidate)
    assert diffs
    assert any("row 'p1' haplotype 'hapA'" in d and "status" in d for d in diffs)


def test_compare_verify_dicts_catches_off_target_size_change():
    """Mirrors the CYP2D7-paralog failure mode: an extra off-target product size."""
    golden = _verify_dict(status="multi_product", on_target=[282], off_target=[279, 285])
    candidate = _verify_dict(status="multi_product", on_target=[282], off_target=[279])
    diffs = compare_verify_dicts(golden, candidate)
    assert any("off_target" in d and "row 'p1' haplotype 'hapA'" in d for d in diffs)


def test_compare_verify_dicts_off_target_order_insensitive():
    golden = _verify_dict(status="multi_product", on_target=[282], off_target=[279, 285])
    candidate = _verify_dict(status="multi_product", on_target=[282], off_target=[285, 279])
    assert compare_verify_dicts(golden, candidate) == []


def test_compare_verify_dicts_catches_dropout_regression():
    """Mirrors the engineered-dropout row: a candidate that now (wrongly) amplifies."""
    golden = _verify_dict(status="dropout", on_target=[], off_target=[])
    candidate = _verify_dict(status="pass", on_target=[122], off_target=[])
    diffs = compare_verify_dicts(golden, candidate)
    assert any("status" in d for d in diffs)
    assert any("on_target" in d for d in diffs)


def test_compare_verify_dicts_missing_row_and_haplotype_reported():
    golden = _verify_dict()
    candidate = {"provenance": {}, "haplotypes": [], "rows": []}
    diffs = compare_verify_dicts(golden, candidate)
    assert any("missing primer_id" in d for d in diffs)


def _results_dict(*, status="pass", amplicons=None, passed=True):
    return {
        "provenance": {},
        "pairs": [
            {
                "name": "pair0",
                "forward": "ACGT",
                "reverse": "TGCA",
                "product_size_chm13": 201,
                "passed": passed,
                "reject_reasons": [],
                "on_target_coverage": 1.0,
                "unique_product_rate": 1.0,
                "n_evaluable": 1,
                "n_uncertain": 0,
                "n_dropout": 0,
                "n_off_target": 0,
                "n_multi_product": 0,
                "primer3_penalty": 0.1,
                "per_haplotype": [
                    {
                        "haplotype_id": "hapA",
                        "status": status,
                        "reason": "x",
                        "amplicons": amplicons or [],
                    }
                ],
            }
        ],
    }


def test_compare_results_dicts_identical_is_no_diff():
    d = _results_dict(amplicons=[{"chrom": "c1", "start": 0, "end": 201, "size": 201, "on_target": True}])
    assert compare_results_dicts(d, json.loads(json.dumps(d))) == []


def test_compare_results_dicts_catches_amplicon_mismatch():
    golden = _results_dict(
        amplicons=[{"chrom": "c1", "start": 0, "end": 201, "size": 201, "on_target": True}]
    )
    candidate = _results_dict(
        amplicons=[{"chrom": "c1", "start": 5, "end": 206, "size": 201, "on_target": True}]
    )
    diffs = compare_results_dicts(golden, candidate)
    assert any("amplicons" in d and "pair 'pair0' haplotype 'hapA'" in d for d in diffs)


def test_compare_results_dicts_amplicons_order_insensitive():
    a = {"chrom": "c1", "start": 0, "end": 100, "size": 100, "on_target": True}
    b = {"chrom": "c1", "start": 500, "end": 600, "size": 100, "on_target": False}
    golden = _results_dict(amplicons=[a, b])
    candidate = _results_dict(amplicons=[b, a])
    assert compare_results_dicts(golden, candidate) == []


# =========================================================================================
# Fast: the frozen golden fixtures themselves still have the shape/rows we depend on
# =========================================================================================


def test_golden_verify_fixture_exists_and_has_off_target_and_dropout_rows():
    assert GOLDEN_VERIFY.exists(), (
        f"{GOLDEN_VERIFY} missing -- re-freeze from demo-verify-pipeline/verify.json"
    )
    data = _load(GOLDEN_VERIFY)
    rows = {r["primer_id"]: r for r in data["rows"]}

    assert "CYP2D6_paralog" in rows, "CYP2D7-off-target row missing from golden verify.json"
    paralog_cells = rows["CYP2D6_paralog"]["cells"]
    assert paralog_cells, "CYP2D6_paralog has no haplotype cells"
    assert all(c["status"] == "multi_product" for c in paralog_cells)
    assert all(c["off_target"] for c in paralog_cells), (
        "CYP2D6_paralog cells must carry off-target product sizes"
    )

    assert "CYP2D6_dropout" in rows, "dropout row missing from golden verify.json"
    dropout_cells = rows["CYP2D6_dropout"]["cells"]
    assert dropout_cells, "CYP2D6_dropout has no haplotype cells"

    # This row is DIFFERENTIAL BETWEEN HAPLOTYPES, not a universal failure, and asserting the
    # latter is how the previous row's defect hid. It exploits rs1058164 (CYP2D6 c.1661G>C,
    # exon 3) under the forward primer's 3' terminal base: HG00097#hap1 carries C and drops
    # out; HG01884#hap1 and HG00408#hap1 carry G and amplify. The superpopulation labels are
    # incidental -- the variant is globally common (11/24 dropout across 30 haplotypes, all
    # five superpopulations), so this is a genotype difference, not a population one.
    #
    # The pair it replaced was an Alu consensus (~330k binding sites) that the old bwa
    # backend's XA cap
    # truncated to a single reported hit, so the pipeline saw no product and called it
    # "dropout" everywhere. `all(status == "dropout")` passed happily on that artefact.
    # Requiring BOTH outcomes is what makes this row evidence of a real 3'-end variant.
    statuses = {c["status"] for c in dropout_cells}
    assert "dropout" in statuses, (
        f"CYP2D6_dropout must drop out in at least one haplotype; got {statuses}"
    )
    assert "pass" in statuses, (
        "CYP2D6_dropout must PASS in at least one haplotype -- a row that fails everywhere "
        f"cannot distinguish a real variant from a broken primer; got {statuses}"
    )
    for c in dropout_cells:
        if c["status"] == "dropout":
            assert not c["on_target"], "a dropout cell must have no on-target product"
        elif c["status"] == "pass":
            assert c["on_target"], "a passing cell must carry an on-target product size"


def test_golden_results_fixture_exists_and_is_well_formed():
    assert GOLDEN_RESULTS.exists(), (
        f"{GOLDEN_RESULTS} missing -- re-freeze from demo-design-pipeline/gapdh/results.json"
    )
    data = _load(GOLDEN_RESULTS)
    assert data["pairs"], "golden results.json has no pairs"
    for p in data["pairs"]:
        assert p["per_haplotype"], f"pair {p['name']} has no per-haplotype results"


# =========================================================================================
# Fast: frozen golden vs. the live (gitignored) demo output in THIS checkout, if present
# =========================================================================================


def test_frozen_golden_matches_live_demo_output_if_present():
    """demo-verify-pipeline/verify.json and demo-design-pipeline/gapdh/results.json are
    gitignored (large, environment-specific), so a fresh clone won't have them without
    running the demo. If they ARE present here, they must still match the committed
    tests/golden/*.json cell for cell -- a mismatch means either the golden fixture is
    stale or something changed the pipeline's output. Skips with a clear message if the
    live files are absent."""
    live_verify = REPO_ROOT / "demo-verify-pipeline" / "verify.json"
    live_results = REPO_ROOT / "demo-design-pipeline" / "gapdh" / "results.json"
    if not live_verify.exists() or not live_results.exists():
        pytest.skip(
            f"live demo outputs not present in this checkout (gitignored): "
            f"{live_verify} exists={live_verify.exists()}, "
            f"{live_results} exists={live_results.exists()}"
        )

    diffs = compare_verify_dicts(_load(GOLDEN_VERIFY), _load(live_verify))
    assert not diffs, "verify.json drifted from tests/golden/verify.json:\n" + "\n".join(diffs)

    diffs = compare_results_dicts(_load(GOLDEN_RESULTS), _load(live_results))
    assert not diffs, "results.json drifted from tests/golden/results.json:\n" + "\n".join(diffs)


# =========================================================================================
# Slow: actually run the pipeline against real HPRC data and diff against golden
# =========================================================================================


def _hprc_data_present(manifest: str) -> bool:
    """True when `manifest`'s haplotypes and CHM13 are all on disk.

    Resolved through `load_haplotypes` against the manifest the test actually passes to the
    CLI, rather than probing a hardcoded filename. A previous version checked for
    `<sample>.fa`; when Phase 3 repointed the manifests at the BGZF `.fa.gz` and deleted the
    uncompressed copies, that check went false and **both slow gates silently skipped** —
    the whole suite reported green in 0.03 s. A guard that hardcodes a path the pipeline no
    longer uses does not protect the test, it disables it.
    """
    from pangenome_primer.samples import load_haplotypes

    if not (REPO_ROOT / "hprc-r2" / "references" / "chm13v2.0.fa").exists():
        return False
    try:
        return bool(load_haplotypes(str(REPO_ROOT / manifest)))
    except (FileNotFoundError, OSError):
        return False


@pytest.mark.parametrize(
    "manifest", ["demo-verify-pipeline/samples.tsv", "demo-design-pipeline/samples.tsv"]
)
def test_slow_gate_is_not_silently_disabled(manifest):
    """If the genome data is on disk, the slow gates must actually run.

    This is the watchdog for the skip guard itself. `_hprc_data_present` once probed a
    hardcoded `<sample>.fa`; when the manifests were repointed at `.fa.gz` and the
    uncompressed copies deleted, it started returning False and both slow gates skipped
    while the suite still reported green. Nothing else in the suite could notice, because a
    skipped test looks exactly like a passing one in the summary line.

    Deliberately keyed on the presence of *any* haplotype file rather than on the guard's
    own logic, so it fails when the two disagree instead of restating the guard.
    """
    have_data = any((REPO_ROOT / "hprc-r2" / "assemblies").glob("*.fa*")) and (
        REPO_ROOT / "hprc-r2" / "references" / "chm13v2.0.fa"
    ).exists()
    if not have_data:
        pytest.skip("no hprc-r2/ genome data on this machine at all")
    assert _hprc_data_present(manifest), (
        f"genome data is present but the slow-gate guard reports it missing for {manifest}. "
        f"The slow gates are being skipped, so the golden matrices are NOT being verified."
    )


@pytest.mark.slow
def test_verify_pipeline_reproduces_golden_verify_matrix(tmp_path):
    """Regenerates verify.json from scratch (real genome-wide search + thermo classification
    across 3 real HPRC haplotypes) and diffs it against tests/golden/verify.json cell for
    cell, including the CYP2D7-off-target and engineered-dropout rows. ~minutes; needs
    real genome data under hprc-r2/."""
    if not _hprc_data_present("demo-verify-pipeline/samples.tsv"):
        pytest.skip("hprc-r2/ genome data not present; see README 'Get the data'")
    if shutil.which("pangenome-primer") is None:
        pytest.skip("'pangenome-primer' console script not on PATH (pip install -e .)")

    outdir = tmp_path / "verify_out"
    subprocess.run(
        [
            "pangenome-primer", "verify",
            "--primers", "demo-verify-pipeline/primers.csv",
            "--chm13", "hprc-r2/references/chm13v2.0.fa",
            "--samples", "demo-verify-pipeline/samples.tsv",
            "--target-assembly", "chm13",
            "--outdir", str(outdir),
        ],
        check=True, cwd=str(REPO_ROOT),
    )
    candidate = _load(outdir / "verify.json")
    diffs = compare_verify_dicts(_load(GOLDEN_VERIFY), candidate)
    assert not diffs, "regenerated verify.json differs from golden:\n" + "\n".join(diffs)


@pytest.mark.slow
def test_design_pipeline_reproduces_golden_results(tmp_path):
    """Regenerates results.json from scratch (design + evaluate + rank for GAPDH across 3
    real HPRC haplotypes) and diffs it against tests/golden/results.json cell for cell.
    ~minutes; needs real genome data under hprc-r2/."""
    if not _hprc_data_present("demo-design-pipeline/samples.tsv"):
        pytest.skip("hprc-r2/ genome data not present; see README 'Get the data'")
    if shutil.which("pangenome-primer") is None:
        pytest.skip("'pangenome-primer' console script not on PATH (pip install -e .)")

    outdir = tmp_path / "gapdh_out"
    subprocess.run(
        [
            "pangenome-primer", "design",
            "--target", "chr12:6544868-6548730",
            "--chm13", "hprc-r2/references/chm13v2.0.fa",
            "--samples", "demo-design-pipeline/samples.tsv",
            "--outdir", str(outdir),
            "--mode", "thermo", "--top-k", "5",
        ],
        check=True, cwd=str(REPO_ROOT),
    )
    candidate = _load(outdir / "results.json")
    diffs = compare_results_dicts(_load(GOLDEN_RESULTS), candidate)
    assert not diffs, "regenerated results.json differs from golden:\n" + "\n".join(diffs)
