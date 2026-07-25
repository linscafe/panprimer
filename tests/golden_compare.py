"""Cell-for-cell comparators for the two golden demo outputs (Phase 0 safety net).

`compare_verify_dicts` diffs `demo-verify-pipeline/verify.json`-shaped dicts (see
`report.verify_to_dict`); `compare_results_dicts` diffs `demo-design-pipeline/*/results.json`-
shaped dicts (see `report.results_to_dict`). Both match rows/haplotypes by KEY
(`primer_id`/`name` and `haplotype_id`), not list position, so a harmless reordering never
produces a spurious diff -- only a real content difference does. Every mismatch is reported
as one human-readable string naming the exact row and cell, so a failing `assert not diffs`
prints something actionable rather than a giant dict diff.

Pure functions, no I/O, no pytest dependency -- importable from both the fast synthetic-dict
unit tests and the slow real-pipeline test.
"""
from __future__ import annotations


def _fmt(value) -> str:
    return repr(value)


def _compare_scalar(diffs: list[str], where: str, field: str, golden, candidate) -> None:
    if golden != candidate:
        diffs.append(f"{where}: {field} golden={_fmt(golden)} candidate={_fmt(candidate)}")


def _compare_size_list(diffs: list[str], where: str, field: str, golden, candidate) -> None:
    """Compare size lists (on_target/off_target) order-insensitively -- these are sorted by
    the producer today, but the comparator shouldn't fail on a harmless reordering."""
    g, c = sorted(golden), sorted(candidate)
    if g != c:
        diffs.append(f"{where}: {field} golden={_fmt(g)} candidate={_fmt(c)}")


# --- verify.json (demo-verify-pipeline) --------------------------------------------------


def compare_verify_dicts(golden: dict, candidate: dict) -> list[str]:
    """Diff two `verify_to_dict`-shaped dicts. Returns a list of human-readable mismatch
    descriptions, empty if they match cell for cell."""
    diffs: list[str] = []

    g_haps, c_haps = set(golden.get("haplotypes", [])), set(candidate.get("haplotypes", []))
    if g_haps != c_haps:
        diffs.append(
            f"haplotypes: golden={sorted(g_haps)} candidate={sorted(c_haps)} "
            f"(missing={sorted(g_haps - c_haps)}, extra={sorted(c_haps - g_haps)})"
        )

    g_rows = {r["primer_id"]: r for r in golden.get("rows", [])}
    c_rows = {r["primer_id"]: r for r in candidate.get("rows", [])}
    missing = g_rows.keys() - c_rows.keys()
    extra = c_rows.keys() - g_rows.keys()
    if missing:
        diffs.append(f"rows: missing primer_id(s) in candidate: {sorted(missing)}")
    if extra:
        diffs.append(f"rows: unexpected extra primer_id(s) in candidate: {sorted(extra)}")

    for primer_id in sorted(g_rows.keys() & c_rows.keys()):
        gr, cr = g_rows[primer_id], c_rows[primer_id]
        where = f"row '{primer_id}'"
        for field in ("forward", "reverse", "target_input", "target_chm13", "expected_size"):
            _compare_scalar(diffs, where, field, gr.get(field), cr.get(field))

        g_cells = {c["haplotype_id"]: c for c in gr.get("cells", [])}
        c_cells = {c["haplotype_id"]: c for c in cr.get("cells", [])}
        missing_h = g_cells.keys() - c_cells.keys()
        extra_h = c_cells.keys() - g_cells.keys()
        if missing_h:
            diffs.append(f"{where}: missing haplotype cell(s) {sorted(missing_h)}")
        if extra_h:
            diffs.append(f"{where}: unexpected extra haplotype cell(s) {sorted(extra_h)}")

        for hid in sorted(g_cells.keys() & c_cells.keys()):
            gc, cc = g_cells[hid], c_cells[hid]
            cell_where = f"row '{primer_id}' haplotype '{hid}'"
            _compare_scalar(diffs, cell_where, "status", gc.get("status"), cc.get("status"))
            _compare_size_list(
                diffs, cell_where, "on_target", gc.get("on_target", []), cc.get("on_target", [])
            )
            _compare_size_list(
                diffs, cell_where, "off_target", gc.get("off_target", []), cc.get("off_target", [])
            )
            _compare_scalar(
                diffs, cell_where, "size_flag", gc.get("size_flag"), cc.get("size_flag")
            )
            # `reason` is a free-text explanation (may include e.g. a Tm value formatted
            # to 1 decimal place); status + sizes are the load-bearing fields, so a reason
            # mismatch is reported but does not need separate coverage here.

    return diffs


# --- results.json (demo-design-pipeline) -------------------------------------------------


def _amplicon_key(a: dict) -> tuple:
    return (a.get("chrom"), a.get("start"), a.get("end"), a.get("size"), a.get("on_target"))


def _compare_amplicons(diffs: list[str], where: str, golden: list[dict], candidate: list[dict]) -> None:
    g = sorted((_amplicon_key(a) for a in golden))
    c = sorted((_amplicon_key(a) for a in candidate))
    if g != c:
        diffs.append(f"{where}: amplicons golden={g} candidate={c}")


def compare_results_dicts(golden: dict, candidate: dict) -> list[str]:
    """Diff two `results_to_dict`-shaped dicts (the design pipeline's results.json).
    Returns a list of human-readable mismatch descriptions, empty if they match cell for
    cell."""
    diffs: list[str] = []

    g_pairs = {p["name"]: p for p in golden.get("pairs", [])}
    c_pairs = {p["name"]: p for p in candidate.get("pairs", [])}
    missing = g_pairs.keys() - c_pairs.keys()
    extra = c_pairs.keys() - g_pairs.keys()
    if missing:
        diffs.append(f"pairs: missing name(s) in candidate: {sorted(missing)}")
    if extra:
        diffs.append(f"pairs: unexpected extra name(s) in candidate: {sorted(extra)}")

    for name in sorted(g_pairs.keys() & c_pairs.keys()):
        gp, cp = g_pairs[name], c_pairs[name]
        where = f"pair '{name}'"
        for field in (
            "forward", "reverse", "product_size_chm13", "passed",
            "n_evaluable", "n_uncertain", "n_dropout", "n_off_target", "n_multi_product",
        ):
            _compare_scalar(diffs, where, field, gp.get(field), cp.get(field))

        g_ph = {h["haplotype_id"]: h for h in gp.get("per_haplotype", [])}
        c_ph = {h["haplotype_id"]: h for h in cp.get("per_haplotype", [])}
        missing_h = g_ph.keys() - c_ph.keys()
        extra_h = c_ph.keys() - g_ph.keys()
        if missing_h:
            diffs.append(f"{where}: missing haplotype(s) {sorted(missing_h)}")
        if extra_h:
            diffs.append(f"{where}: unexpected extra haplotype(s) {sorted(extra_h)}")

        for hid in sorted(g_ph.keys() & c_ph.keys()):
            gh, ch = g_ph[hid], c_ph[hid]
            cell_where = f"pair '{name}' haplotype '{hid}'"
            _compare_scalar(diffs, cell_where, "status", gh.get("status"), ch.get("status"))
            _compare_amplicons(
                diffs, cell_where, gh.get("amplicons", []), ch.get("amplicons", [])
            )

    return diffs
