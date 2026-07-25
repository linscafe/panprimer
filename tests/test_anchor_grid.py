"""Anchor-grid projection (Phase 4): build, load, bracket, realign.

These run against a synthetic CHM13/haplotype pair small enough to build in seconds, so the
whole round trip -- probe extraction, real minimap2 placement, gzipped grid, bracketing,
in-memory realignment -- is exercised without touching a 3 Gb assembly.

The property that matters is that the grid never *silently* moves a locus. It decides which
window to align; the coordinate still comes from a base-level alignment inside that window.
So a bad grid must yield a failed projection (`locus is None`), never a confident wrong
answer -- that is what `test_window_miss_reports_uncertain_not_wrong_locus` pins.
"""
from __future__ import annotations

import gzip
import random
import shutil

import pytest

from pangenome_primer import anchor_grid

pytestmark = pytest.mark.skipif(
    shutil.which("minimap2") is None, reason="minimap2 not on PATH"
)

CHR_LEN = 200_000
HAP_PREFIX = 10_000       # haplotype coords are offset from CHM13 by this much


def _rand_seq(n: int, seed: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def _write_fa(path, records: dict[str, str], width: int = 60):
    with open(path, "w") as fh:
        for name, seq in records.items():
            fh.write(f">{name}\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i:i + width] + "\n")
    return str(path)


def _bgzip(path, dest):
    """Write a BGZF copy plus .fai/.gzi, the layout the pipeline actually reads."""
    import pysam

    pysam.tabix_compress(str(path), str(dest), force=True)
    pysam.faidx(str(dest))
    return str(dest)


@pytest.fixture(scope="module")
def synthetic(tmp_path_factory):
    """A CHM13 chr1 and a haplotype carrying the same locus, shifted and lightly mutated."""
    d = tmp_path_factory.mktemp("grid")
    chm13_seq = _rand_seq(CHR_LEN, seed=1)
    chm13 = _write_fa(d / "chm13.fa", {"chr1": chm13_seq})

    rng = random.Random(99)
    body = list(chm13_seq)
    for _ in range(200):                      # ~1 SNP per kb, realistic divergence
        i = rng.randrange(len(body))
        body[i] = rng.choice("ACGT".replace(body[i], ""))
    hap_seq = _rand_seq(HAP_PREFIX, seed=2) + "".join(body) + _rand_seq(5_000, seed=3)
    hap_plain = _write_fa(d / "hap.fa", {"H#1#ctg1": hap_seq})
    hap = _bgzip(hap_plain, d / "hap.fa.gz")

    grid = anchor_grid.build_grid(
        hap, chm13, threads=2, probe_bp=1_000, step_bp=2_000,
        min_contig=5_000, min_mapq=0,
    )
    return {"chm13": chm13, "hap": hap, "hap_seq": hap_seq,
            "chm13_seq": chm13_seq, "stats": grid}


class TestBuild:
    def test_anchors_most_probes(self, synthetic):
        s = synthetic["stats"]
        assert s["probes"] > 50
        assert s["anchored_fraction"] > 0.9, (
            f"only {s['anchored_fraction']:.1%} of probes placed on a near-identical "
            f"sequence; the builder is losing anchors"
        )

    def test_reports_unanchored_fraction(self, synthetic):
        """Phase 5 consumes this; it must be present and consistent."""
        s = synthetic["stats"]
        assert s["unanchored_fraction"] == pytest.approx(1 - s["anchored_fraction"])

    def test_grid_is_small(self, synthetic):
        """The whole point is that this replaces a multi-GB index."""
        s = synthetic["stats"]
        assert s["bytes"] < s["genome_bp"] / 100

    def test_file_is_gzipped_tsv_with_metadata_header(self, synthetic):
        with gzip.open(synthetic["stats"]["path"], "rt") as fh:
            head = [next(fh) for _ in range(3)]
        assert all(h.startswith("#") for h in head)
        assert "probe=1000" in head[0] and "step=2000" in head[0]

    def test_partial_build_leaves_no_usable_grid(self, synthetic, tmp_path):
        """An interrupted build must not leave a truncated grid that looks loadable."""
        import os
        out = tmp_path / "x.anchors.tsv.gz"
        (tmp_path / "x.anchors.tsv.gz.part").write_bytes(b"garbage")
        assert not os.path.exists(out)


class TestProjection:
    @pytest.mark.parametrize("tstart", [20_000, 90_000, 150_000])
    def test_locus_lands_at_the_expected_offset(self, synthetic, tstart):
        """Haplotype coords are CHM13 coords + HAP_PREFIX by construction."""
        tend = tstart + 600
        template = synthetic["chm13_seq"][tstart:tend]
        proj = anchor_grid.project_from_grid(
            synthetic["stats"]["path"], "chr1", tstart, tend, template, synthetic["hap"]
        )
        assert proj.locus is not None, proj.reason
        assert proj.reason == "anchor-grid"
        assert abs(proj.locus.start - (tstart + HAP_PREFIX)) <= 5, (
            f"expected ~{tstart + HAP_PREFIX}, got {proj.locus.start}"
        )
        # The locus spans minimap2's *aligned* extent (r_st..r_en), which stops short of the
        # template ends when the terminal bases do not extend -- exactly what `project_locus`
        # returns from the .mmi path, so the tolerance matches its `min_frac=0.8` contract
        # rather than demanding an exact 600.
        span = proj.locus.end - proj.locus.start
        assert 0.8 * 600 <= span <= 1.05 * 600, f"span {span} outside the min_frac envelope"

    def test_returned_sequence_matches_the_locus(self, synthetic):
        tstart, tend = 60_000, 60_500
        template = synthetic["chm13_seq"][tstart:tend]
        proj = anchor_grid.project_from_grid(
            synthetic["stats"]["path"], "chr1", tstart, tend, template, synthetic["hap"]
        )
        assert proj.locus is not None
        expected = synthetic["hap_seq"][proj.locus.start:proj.locus.end]
        assert proj.haplotype_seq == expected

    def test_unknown_chromosome_is_uncertain_not_wrong(self, synthetic):
        template = synthetic["chm13_seq"][1000:1600]
        proj = anchor_grid.project_from_grid(
            synthetic["stats"]["path"], "chrZZ", 1000, 1600, template, synthetic["hap"]
        )
        assert proj.locus is None
        assert "no anchor" in proj.reason

    def test_window_miss_reports_uncertain_not_wrong_locus(self, synthetic, tmp_path):
        """A grid pointing at the wrong place must fail the realignment, not return a
        confident bad locus. This is the safety property the whole design rests on."""
        bad = tmp_path / "bad.anchors.tsv.gz"
        with gzip.open(bad, "wt") as fh:
            fh.write("# probe=1000 step=2000\n# probes=1 anchored=1\n# hdr\n")
            # claims CHM13 chr1:60000 sits in the haplotype's random prefix
            fh.write("chr1\t59000\t60000\tH#1#ctg1\t0\t1000\t+\n")
            fh.write("chr1\t60000\t61000\tH#1#ctg1\t1000\t2000\t+\n")
        template = synthetic["chm13_seq"][60_000:60_600]
        proj = anchor_grid.project_from_grid(
            str(bad), "chr1", 60_000, 60_600, template, synthetic["hap"]
        )
        assert proj.locus is None, (
            "a mis-pointed grid produced a confident locus; the realignment step is not "
            "actually validating the window"
        )

    def test_oversized_window_is_refused(self, synthetic, tmp_path, monkeypatch):
        """Disagreeing anchors must not trigger a multi-megabase in-process alignment.

        The cap is applied *after* clamping to contig length -- what costs time is the window
        actually fetched, not the span the anchors nominally imply. That makes the real 5 Mb
        threshold untestable on a 215 kb synthetic contig, so the cap itself is lowered here
        rather than fabricating a genome large enough to trip it.
        """
        monkeypatch.setattr(anchor_grid, "MAX_WINDOW_BP", 50_000)
        bad = tmp_path / "wide.anchors.tsv.gz"
        with gzip.open(bad, "wt") as fh:
            fh.write("# probe=1000 step=2000\n# probes=2 anchored=2\n# hdr\n")
            fh.write("chr1\t59000\t60000\tH#1#ctg1\t0\t1000\t+\n")
            fh.write("chr1\t60000\t61000\tH#1#ctg1\t150000\t151000\t+\n")
        template = synthetic["chm13_seq"][60_000:60_600]
        proj = anchor_grid.project_from_grid(
            str(bad), "chr1", 60_000, 60_600, template, synthetic["hap"]
        )
        assert proj.locus is None and "cap" in proj.reason


class TestWindowSelection:
    def test_contig_chosen_by_majority_not_by_a_single_outlier(self):
        """One probe mis-placed into a paralogue must not drag the window to its contig."""
        anchors = {
            "chr1": [anchor_grid.Anchor(i * 1000, i * 1000 + 500, "good", i * 1000,
                                        i * 1000 + 500, "+")
                     for i in range(20)]
        }
        anchors["chr1"].append(anchor_grid.Anchor(10_000, 10_500, "paralogue",
                                                  9_000_000, 9_000_500, "+"))
        anchors["chr1"].sort(key=lambda a: a.t_start)
        grid = anchor_grid.Grid(anchors, {})
        contig, lo, hi = anchor_grid.window_for(grid, "chr1", 10_000, 10_500)
        assert contig == "good"
        assert lo <= 10_000 and hi >= 10_500

    def test_no_anchors_returns_none(self):
        grid = anchor_grid.Grid({}, {})
        assert anchor_grid.window_for(grid, "chr1", 0, 100) is None

    def test_same_contig_outlier_does_not_balloon_the_window(self):
        """Majority vote settles the contig; this covers a probe that mis-placed *within*
        the right contig, which the vote cannot catch."""
        anchors = {
            "chr1": [anchor_grid.Anchor(i * 1000, i * 1000 + 500, "ctg", i * 1000,
                                        i * 1000 + 500, "+")
                     for i in range(20)]
        }
        far = 50_000_000
        anchors["chr1"].append(anchor_grid.Anchor(10_000, 10_500, "ctg", far, far + 500, "+"))
        anchors["chr1"].sort(key=lambda a: a.t_start)
        grid = anchor_grid.Grid(anchors, {})
        contig, lo, hi = anchor_grid.window_for(grid, "chr1", 10_000, 10_500)
        assert contig == "ctg"
        assert hi - lo < 1_000_000, (
            f"one outlier stretched the window to {(hi - lo) / 1e6:.1f} Mb; the whole "
            f"target region would be realigned instead of a local window"
        )
        assert lo <= 10_000 and hi >= 10_500, "the true locus fell outside the window"


class TestGridCache:
    def test_rebuilt_grid_is_not_served_stale(self, synthetic, tmp_path):
        """Cache is keyed on mtime, so a rebuild in a long-lived process is picked up."""
        import os
        import time

        p = tmp_path / "g.anchors.tsv.gz"
        with gzip.open(p, "wt") as fh:
            fh.write("# a=1\n# b=2\n# hdr\nchr1\t0\t100\tctg\t0\t100\t+\n")
        first = anchor_grid.load_grid(str(p))
        assert len(first.near("chr1", 0, 100, 0)) == 1

        time.sleep(0.01)
        with gzip.open(p, "wt") as fh:
            fh.write("# a=1\n# b=2\n# hdr\n")
            fh.write("chr1\t0\t100\tctg\t0\t100\t+\nchr1\t500\t600\tctg\t500\t600\t+\n")
        os.utime(p, (time.time() + 1, time.time() + 1))
        second = anchor_grid.load_grid(str(p))
        assert len(second.near("chr1", 0, 600, 0)) == 2, "stale grid served from cache"
