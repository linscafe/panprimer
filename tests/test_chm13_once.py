"""CHM13-once candidate discovery (Phase 5).

The mode scans CHM13 once, lifts the hits onto each haplotype through the anchor grid, and
runs the exact comparator over only those windows. Its correctness claim is narrow and
worth stating precisely:

    Within a window the search is exhaustive. A site is missed only by not being in any
    window -- never by an approximation inside one. So `chm13-once` may return FEWER sites
    than `exhaustive`, and must never return a site `exhaustive` does not, nor a site with
    different coordinates.

That asymmetry is what these tests pin: `EXTRA == 0` is a hard invariant, `MISSED > 0` is
the documented envelope.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent))

from fixtures.mini_genome import build_mini_genome  # noqa: E402

from pangenome_primer import chm13_once  # noqa: E402
from pangenome_primer.model import BindingSite, Strand  # noqa: E402


def _site(chrom, start, length=20):
    return BindingSite("p", "h", chrom, start, start + length, Strand.PLUS, 0, [])


class TestMergeRegions:
    def test_overlapping_sites_collapse(self):
        sites = {"p": [_site("chr1", 1000), _site("chr1", 1500)]}
        got = chm13_once.merge_regions(sites, pad=2000)
        assert len(got) == 1, "sites 500 bp apart with 2 kb padding must merge"
        assert got[0].start == 0 and got[0].end == 3520

    def test_distant_sites_stay_separate(self):
        sites = {"p": [_site("chr1", 1000), _site("chr1", 500_000)]}
        assert len(chm13_once.merge_regions(sites, pad=2000)) == 2

    def test_regions_are_split_per_chromosome(self):
        sites = {"p": [_site("chr1", 1000), _site("chr2", 1000)]}
        got = chm13_once.merge_regions(sites, pad=2000)
        assert {r.chrom for r in got} == {"chr1", "chr2"}

    def test_padding_never_goes_negative(self):
        got = chm13_once.merge_regions({"p": [_site("chr1", 10)]}, pad=2000)
        assert got[0].start == 0

    def test_no_sites_gives_no_regions(self):
        assert chm13_once.merge_regions({"p": []}) == []


class TestMergeWindows:
    def test_overlaps_collapse_and_contigs_stay_separate(self):
        got = chm13_once._merge_windows(
            [("a", 0, 100), ("a", 50, 200), ("a", 500, 600), ("b", 0, 100)]
        )
        assert sorted(got) == [("a", 0, 200), ("a", 500, 600), ("b", 0, 100)]


@pytest.fixture(scope="module")
def mini(tmp_path_factory):
    """One BGZF mini-genome per module: building it costs a bgzip + faidx round trip."""
    return build_mini_genome(tmp_path_factory.mktemp("chm13once"))


class TestWindowScan:
    """The concatenate-and-scan path, against a real BGZF fixture with known sites."""

    def test_finds_planted_sites_at_correct_global_coordinates(self, mini):
        """Coordinates must come back in CONTIG space, not window or concatenated space."""
        planted = [p for p in mini.planted if p.chrom == "ctgB"][:3]
        assert planted, "fixture should plant sites on ctgB"
        windows = [(p.chrom, max(0, p.start - 500), p.start + 500) for p in planted]
        seqs = [p.primer_seq for p in planted]
        got = chm13_once.find_binding_sites_in_windows(
            seqs, mini.fasta_gz, mini.haplotype_id, 3, windows
        )
        for p in planted:
            hits = {(s.chrom, s.start) for s in got[p.primer_seq]}
            assert (p.chrom, p.start) in hits, (
                f"planted site {p.chrom}:{p.start} not recovered at its global coordinate; "
                f"got {sorted(hits)[:5]}"
            )

    def test_site_length_is_preserved(self, mini):
        """Regression: `end` was once computed from an already-reassigned `start`, leaving
        it in concatenated coordinates. That surfaced only much later, as a pysam
        'start > stop' on an unrelated fetch -- so assert the invariant here at the source."""
        p = mini.planted[0]
        windows = [(p.chrom, max(0, p.start - 500), p.start + 500)]
        got = chm13_once.find_binding_sites_in_windows(
            [p.primer_seq], mini.fasta_gz, mini.haplotype_id, 3, windows
        )
        for s in got[p.primer_seq]:
            assert s.end > s.start, f"end {s.end} <= start {s.start}"
            assert s.end - s.start == len(p.primer_seq), (
                f"site spans {s.end - s.start} bp for a {len(p.primer_seq)} bp primer"
            )

    def test_no_hit_spans_a_window_boundary(self, mini):
        """Concatenation must not manufacture sites at the junctions between windows."""
        p = mini.planted[0]
        # Two widely separated windows; anything reported between them is fabricated.
        windows = [("ctgA", 1000, 3000), ("ctgC", 1000, 3000)]
        got = chm13_once.find_binding_sites_in_windows(
            [p.primer_seq], mini.fasta_gz, mini.haplotype_id, 3, windows
        )
        for s in got[p.primer_seq]:
            in_a = s.chrom == "ctgA" and 1000 <= s.start and s.end <= 3000
            in_c = s.chrom == "ctgC" and 1000 <= s.start and s.end <= 3000
            assert in_a or in_c, (
                f"site {s.chrom}:{s.start}-{s.end} lies outside every requested window -- "
                f"a separator-spanning match leaked through"
            )

    def test_overlapping_windows_do_not_duplicate_sites(self, mini):
        p = mini.planted[0]
        w = (p.chrom, max(0, p.start - 500), p.start + 500)
        overlapping = [w, (p.chrom, max(0, p.start - 400), p.start + 600)]
        got = chm13_once.find_binding_sites_in_windows(
            [p.primer_seq], mini.fasta_gz, mini.haplotype_id, 3, overlapping
        )
        keys = [(s.chrom, s.start, s.strand.value) for s in got[p.primer_seq]]
        assert len(keys) == len(set(keys)), "overlapping windows produced duplicate sites"

    def test_empty_window_list_is_not_an_error(self, mini):
        got = chm13_once.find_binding_sites_in_windows(
            ["ACGTACGTACGTACGTACGT"], mini.fasta_gz, mini.haplotype_id, 3, []
        )
        assert got == {"ACGTACGTACGTACGTACGT": []}

    def test_subset_of_exhaustive_over_the_same_span(self, mini):
        """The core invariant: scanning a window must agree exactly with scanning the whole
        contig, restricted to that window. Never extra, never shifted."""
        from pangenome_primer.search import find_binding_sites_batch

        p = mini.planted[0]
        lo, hi = max(0, p.start - 2000), p.start + 2000
        windowed = chm13_once.find_binding_sites_in_windows(
            [p.primer_seq], mini.fasta_gz, mini.haplotype_id, 3, [(p.chrom, lo, hi)]
        )
        full = find_binding_sites_batch(
            [p.primer_seq], mini.fasta_gz, mini.haplotype_id, 3, backend="naive"
        )
        def key(s):
            return (s.chrom, s.start, s.end, s.strand.value, s.mismatches)
        got = {key(s) for s in windowed[p.primer_seq]}
        expected = {
            key(s) for s in full[p.primer_seq]
            if s.chrom == p.chrom and s.start >= lo and s.end <= hi
        }
        assert got == expected, (
            f"windowed scan disagrees with the exhaustive scan over the same span; "
            f"extra={got - expected} missing={expected - got}"
        )


class TestGridRequirement:
    def test_missing_grid_fails_loudly(self, tmp_path):
        """Silently falling back to a whole-genome scan would hide the exact cost the
        caller opted in to avoid."""
        searcher = chm13_once.Chm13OnceSearcher.__new__(chm13_once.Chm13OnceSearcher)
        searcher.seqs = ["ACGTACGTACGTACGTACGT"]
        searcher.max_mismatches = 3
        searcher.regions = [chm13_once.Region("chr1", 0, 100)]
        fake = tmp_path / "nogrid.fa.gz"
        fake.write_bytes(b"\x1f\x8b")
        with pytest.raises(FileNotFoundError, match="anchor grid"):
            searcher.sites_for(str(fake), "H#hap1")
