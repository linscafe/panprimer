"""Unit tests for verify-mode CSV parsing and matrix rendering (no external tools)."""
from __future__ import annotations

from pangenome_primer import report
from pangenome_primer.verify import VerifyCell, VerifyRow, parse_primer_csv


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
