"""The binding-site cap: repeat-derived primers report ">N binding sites", not products.

Motivation is concrete. The demo's original `CYP2D6_dropout` forward primer is an Alu
consensus with ~330k binding sites (~106k of them amplification-competent). Thermodynamically
scoring every one pushed a 3-haplotype verify run to 25 min, and the resulting cell was a
wall of ~2000 product sizes spanning 85-2000 bp -- noise, not a result.

The cap stops scanning once a primer exceeds `EvalConfig.max_binding_sites` competent sites
and reports that fact instead. Two properties matter and are both tested here:

1. it must SHORT-CIRCUIT (stop scoring), or it saves display noise but not the runtime; and
2. a capped cell must NOT render as "dropout" -- a promiscuous primer read as a failure to
   amplify is precisely the mislabel that put the wrong pair in the demo to begin with.
"""
from __future__ import annotations

import pytest

from pangenome_primer.engine import EvalConfig, _filter_competent, evaluate_with_sites
from pangenome_primer.model import (
    BindingSite,
    HaplotypeStatus,
    Locus,
    Primer,
    PrimerPair,
    Strand,
)

PRIMER_F = "TTGCACAGTCCAGATTGCAA"
PRIMER_R = "CATTGCGATTGACACTTGCG"
HAP = "CAP#hap1"


def _site(
    start: int, primer: str, strand: Strand = Strand.PLUS, chrom: str = "ctg"
) -> BindingSite:
    """A perfect (0-mismatch) site -- competent under any sane config.

    `chrom` doubles as the key `_pair_window` uses to hand each site the right template:
    both demo primers are 20 nt, so window length cannot distinguish them.
    """
    return BindingSite(
        primer_name=primer,
        haplotype_id=HAP,
        chrom=chrom,
        start=start,
        end=start + len(primer),
        strand=strand,
        mismatches=0,
        mismatch_offsets_3p=[],
    )


def _revcomp(s: str) -> str:
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _pair_window(chrom, start, end):
    """Perfect top-strand template for whichever primer the site belongs to.

    Keyed by contig because both primers are 20 nt. `ctgR` holds MINUS-strand reverse-primer
    sites, and `classify.competent_thermo` revcomps the window for MINUS -- so the top-strand
    template there must be `revcomp(PRIMER_R)` for the primer to score as a perfect match.
    """
    return _revcomp(PRIMER_R) if chrom == "ctgR" else PRIMER_F


class TestFilterCompetentCap:
    def test_returns_capped_false_below_the_limit(self):
        cfg = EvalConfig(mode="thermo", max_binding_sites=10)
        sites = [_site(i * 100, PRIMER_F) for i in range(5)]
        kept, _, capped = _filter_competent(PRIMER_F, sites, _pair_window, cfg)
        assert capped is False
        assert len(kept) == 5

    def test_exactly_at_the_limit_is_not_capped(self):
        """Boundary: 10 competent sites with a limit of 10 is fine; the 11th trips it."""
        cfg = EvalConfig(mode="thermo", max_binding_sites=10)
        sites = [_site(i * 100, PRIMER_F) for i in range(10)]
        kept, _, capped = _filter_competent(PRIMER_F, sites, _pair_window, cfg)
        assert capped is False
        assert len(kept) == 10

    def test_one_over_the_limit_trips_the_cap(self):
        cfg = EvalConfig(mode="thermo", max_binding_sites=10)
        sites = [_site(i * 100, PRIMER_F) for i in range(11)]
        kept, _, capped = _filter_competent(PRIMER_F, sites, _pair_window, cfg)
        assert capped is True
        assert len(kept) == 10, "kept should be truncated to exactly the cap"

    def test_short_circuits_rather_than_scoring_everything(self):
        """The whole point: scanning must STOP at the cap. A cap that still scored all
        330k sites would fix the display and none of the 25-minute runtime."""
        cfg = EvalConfig(mode="thermo", max_binding_sites=10)
        calls = {"n": 0}

        def counting_window(chrom, start, end):
            calls["n"] += 1
            return PRIMER_F

        sites = [_site(i * 100, PRIMER_F) for i in range(100_000)]
        _, _, capped = _filter_competent(PRIMER_F, sites, counting_window, cfg)
        assert capped is True
        assert calls["n"] <= 12, (
            f"scored {calls['n']} sites; must stop at the cap, not scan all 100k"
        )

    @pytest.mark.parametrize("disabled", [None, 0])
    def test_cap_can_be_disabled(self, disabled):
        cfg = EvalConfig(mode="thermo", max_binding_sites=disabled)
        sites = [_site(i * 100, PRIMER_F) for i in range(50)]
        kept, _, capped = _filter_competent(PRIMER_F, sites, _pair_window, cfg)
        assert capped is False
        assert len(kept) == 50


class TestEvaluateWithSitesReportsTheCap:
    def _pair(self):
        return PrimerPair("capped", Primer("F", PRIMER_F), Primer("R", PRIMER_R))

    def test_capped_primer_yields_multi_product_not_dropout(self):
        """A promiscuous primer must never be reported as a dropout -- that mislabel is
        exactly how an Alu primer ended up in the demo as the dropout example."""
        cfg = EvalConfig(mode="thermo", max_binding_sites=10)
        f_sites = [_site(i * 500, PRIMER_F) for i in range(200)]
        r_sites = [_site(300, PRIMER_R, Strand.MINUS, chrom="ctgR")]

        res = evaluate_with_sites(
            self._pair(), HAP, f_sites, r_sites,
            Locus("hap", "ctg", 0, 100_000), _pair_window, cfg,
        )
        assert res.status is HaplotypeStatus.MULTI_PRODUCT
        assert res.status is not HaplotypeStatus.DROPOUT
        assert ">10 binding sites" in res.reason
        assert "forward" in res.reason
        assert res.amplicons == [], "capped cells enumerate no products by design"

    def test_reason_names_which_primer_is_promiscuous(self):
        cfg = EvalConfig(mode="thermo", max_binding_sites=5)
        r_sites = [
            _site(i * 500, PRIMER_R, Strand.MINUS, chrom="ctgR") for i in range(50)
        ]
        f_sites = [_site(10, PRIMER_F)]

        res = evaluate_with_sites(
            self._pair(), HAP, f_sites, r_sites,
            Locus("hap", "ctg", 0, 100_000), _pair_window, cfg,
        )
        assert "reverse" in res.reason and "forward" not in res.reason
        assert ">5 binding sites" in res.reason


class TestCappedCellSurvivesSerialization:
    """The engine-level tests above pass even when the cap never reaches the report.

    That gap was real: `VerifyCell.site_cap` was added to the dataclass and the Jinja
    template, but `report.verify_to_dict` builds its cell dicts from an explicit field list
    and did not include it. Every renderer reads that dict, not the dataclass, so capped
    cells silently rendered as "dropout" in JSON, TSV, Markdown and HTML -- the precise
    mislabel the cap exists to prevent. These tests walk the full path.
    """

    def _rows(self, site_cap):
        from pangenome_primer.verify import VerifyCell, VerifyRow

        return [
            VerifyRow(
                primer_id="repeaty", forward=PRIMER_F, reverse=PRIMER_R,
                target_input="chr1:100-300", target_chm13="chr1:100-300",
                expected_size=200,
                cells=[
                    VerifyCell(
                        "H#hap1", "multi_product", [], [], False,
                        f">{site_cap} binding sites (forward primer)", site_cap=site_cap,
                    )
                ],
            )
        ]

    def test_site_cap_is_present_in_the_serialized_dict(self):
        from pangenome_primer.report import verify_to_dict

        d = verify_to_dict(self._rows(100))
        cell = d["rows"][0]["cells"][0]
        assert "site_cap" in cell, "site_cap must be serialized or every renderer loses it"
        assert cell["site_cap"] == 100

    def test_tsv_text_reports_the_cap_not_dropout(self):
        from pangenome_primer.report import _cell_text, verify_to_dict

        cell = verify_to_dict(self._rows(100))["rows"][0]["cells"][0]
        assert _cell_text(cell) == ">100 binding sites"
        assert _cell_text(cell) != "dropout"

    def test_markdown_reports_the_cap_not_dropout(self):
        from pangenome_primer.report import _verify_cell_md, verify_to_dict

        cell = verify_to_dict(self._rows(100))["rows"][0]["cells"][0]
        out = _verify_cell_md(cell)
        assert ">100 binding sites" in out
        assert "dropout" not in out

    def test_uncapped_cell_still_reports_dropout(self):
        """The guard must not swallow genuine dropouts."""
        from pangenome_primer.report import _cell_text, _verify_cell_md, verify_to_dict

        cell = verify_to_dict(self._rows(None))["rows"][0]["cells"][0]
        assert _cell_text(cell) == "dropout"
        assert "dropout" in _verify_cell_md(cell)
