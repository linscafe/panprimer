"""Unit tests for verify-mode CSV parsing and matrix rendering (no external tools)."""
from __future__ import annotations

from pangenome_primer import report
from pangenome_primer.verify import VerifyCell, VerifyRow, parse_primer_csv, summarize


def test_parse_primer_csv(tmp_path):
    p = tmp_path / "primers.csv"
    p.write_text(
        "primer_id,target,forward,reverse\n"
        "P1,chr22:100-380,acgtACGTacgtACGTacgt,TTGCTTGCTTGCTTGCTTGC\n"
        "P2,chr7:5-405,AAAACCCCGGGGTTTTAAAA,CCCCGGGGTTTTAAAACCCC\n"
    )
    specs = parse_primer_csv(str(p))
    assert [s.primer_id for s in specs] == ["P1", "P2"]
    assert specs[0].target == "chr22:100-380"
    assert specs[0].forward == "ACGTACGTACGTACGTACGT"  # upper-cased


def test_parse_primer_csv_header_aliases(tmp_path):
    p = tmp_path / "primers.csv"
    p.write_text("id,region,fwd,rev\nX,chr1:1-100,ACGT,TGCA\n")
    specs = parse_primer_csv(str(p))
    assert specs[0].primer_id == "X" and specs[0].reverse == "TGCA"


def test_parse_primer_csv_missing_column(tmp_path):
    p = tmp_path / "bad.csv"
    p.write_text("primer_id,target,forward\nX,chr1:1-100,ACGT\n")
    import pytest

    with pytest.raises(ValueError, match="reverse"):
        parse_primer_csv(str(p))


def _rows():
    return [
        VerifyRow("P1", "F", "R", "chr22:100-380", "chr22:90-370", 280, cells=[
            VerifyCell("HGA#hap1", "pass", on_target=[280]),
            VerifyCell("HGB#hap1", "multi_product", on_target=[279], off_target=[282, 282]),
            VerifyCell("HGC#hap1", "dropout"),
            VerifyCell("HGD#hap1", "uncertain", reason="not projected"),
            VerifyCell("HGE#hap1", "pass", on_target=[305], size_flag=True),
        ]),
    ]


def test_write_verify_outputs(tmp_path):
    paths = report.write_verify(_rows(), str(tmp_path), provenance={"target_assembly": "grch38"})
    html = open(paths["html"]).read()
    tsv = open(paths["tsv"]).read()
    import json
    data = json.load(open(paths["json"]))

    # matrix has all haplotype columns and the primer row
    assert data["haplotypes"] == ["HGA#hap1", "HGB#hap1", "HGC#hap1", "HGD#hap1", "HGE#hap1"]
    # on-target green, off-target red, dropout, uncertain, size-deviation flag all present
    assert 'class="ok' in html and 'class="off"' in html
    assert "dropout" in html and ">?<" in html
    assert "dev" in html          # size_flag on P1/HGE
    # TSV cell text: on-target size, off:list, dropout, ?
    assert "280" in tsv and "off:282,282" in tsv and "dropout" in tsv and "?" in tsv


class TestRowSummary:
    """The per-row summary is what makes a 30-haplotype matrix readable: at that width the
    verdict is not recoverable by scanning cells. These lock the definitions it shares with
    the design pipeline's `model.PairResult`."""

    def test_counts_and_denominator_exclude_uncertain(self):
        s = summarize(_rows()[0].cells)
        assert s.total == 5
        # 'uncertain' means the locus could not be projected -- no evidence either way, so
        # it leaves the denominator rather than counting as a failure.
        assert s.n_uncertain == 1 and s.evaluable == 4
        assert (s.n_pass, s.n_dropout, s.n_multi_product, s.n_off_target) == (2, 1, 1, 0)

    def test_coverage_credits_multi_product_not_just_pass(self):
        """The mislabel this guards: a pair whose target always amplifies but always throws
        an extra band is a real, usable-with-caveats result, not a 0%-coverage dead pair.
        `model.on_target_coverage` says the same thing; verify must not disagree with it."""
        cells = [VerifyCell(f"h{i}", "multi_product", on_target=[280], off_target=[400])
                 for i in range(4)]
        s = summarize(cells)
        assert s.coverage == 1.0        # every haplotype makes the intended product
        assert s.unique_rate == 0.0     # ...and none of them makes *only* it
        assert s.n_covered == 4

    def test_coverage_and_unique_rate_on_the_mixed_row(self):
        s = summarize(_rows()[0].cells)
        # covered = pass(280) + multi_product(279) + pass(305); dropout has no on-target.
        assert s.n_covered == 3
        assert s.coverage == 0.75       # 3 of 4 evaluable
        assert s.unique_rate == 0.5     # 2 of 4 evaluable are 'pass'

    def test_capped_cells_are_reported_separately_from_dropouts(self):
        """A capped cell enumerates no products by design. Folding it into n_dropout would
        make a promiscuous primer indistinguishable from one that fails to amplify."""
        cells = [VerifyCell("h0", "multi_product", site_cap=100),
                 VerifyCell("h1", "dropout")]
        s = summarize(cells)
        assert s.n_capped == 1 and s.n_dropout == 1
        assert s.n_covered == 0

    def test_all_uncertain_row_does_not_divide_by_zero(self):
        s = summarize([VerifyCell("h0", "uncertain"), VerifyCell("h1", "uncertain")])
        assert s.evaluable == 0 and s.coverage == 0.0 and s.unique_rate == 0.0

    def test_summary_is_derived_from_cells_not_stored(self):
        """Mutating the cells changes the summary, because it is recomputed. A stored copy
        could disagree with the matrix printed beside it."""
        row = _rows()[0]
        before = row.summary.n_dropout
        row.cells.append(VerifyCell("HGF#hap1", "dropout"))
        assert row.summary.n_dropout == before + 1


def test_summary_reaches_every_renderer(tmp_path):
    """JSON, TSV, Markdown and HTML must all carry the verdict. The TSV gets sortable
    numeric columns rather than one packed string, because that file is meant to be loaded
    and ranked by coverage -- which '75% (3/4) · 1 dropout' would defeat."""
    import json

    paths = report.write_verify(_rows(), str(tmp_path))
    data = json.load(open(paths["json"]))
    s = data["rows"][0]["summary"]
    assert s["coverage"] == 0.75 and s["n_covered"] == 3 and s["evaluable"] == 4

    tsv = open(paths["tsv"]).read()
    header, first = tsv.splitlines()[0].split("\t"), tsv.splitlines()[1].split("\t")
    assert "coverage" in header and "n_dropout" in header
    assert first[header.index("coverage")] == "0.7500"
    assert first[header.index("n_dropout")] == "1"

    md = report.verify_to_markdown(data)
    assert "**summary**" in md and "75% (3/4)" in md and "1 dropout" in md

    # Transposed: 'summary' is a row label, not a column header.
    html = open(paths["html"]).read()
    assert '<th class="rowlab">summary</th>' in html and "75%" in html and "1 dropout" in html
    # Below the 0.95 default rank.min_coverage gate, so the figure is flagged, not neutral.
    assert "covlow" in html
