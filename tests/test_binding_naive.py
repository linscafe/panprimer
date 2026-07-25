"""Characterization tests for `binding.find_binding_sites_naive` (Phase 0 safety net).

These pin down the contract `binding.py`'s module docstring declares and that
the backend contract names explicitly, since this backend is the reference every
future backend (rust, anchor-lift, pangenome index) is compared against:

  1. `chrom`/`start`/`end` are haplotype genome coordinates.
  2. `end - start == len(primer)` for every returned site.
  3. `mismatch_offsets_3p` is measured from the primer 3' end (0 == 3' terminal base).
  4. `strand` is correct: PLUS matches T[start:end]; MINUS matches revcomp(T[start:end]).
  5. (not directly testable here) a superset return is acceptable for backends that
     delegate rescoring; naive itself is exact, which trivially satisfies "superset".
  6. `N` is never special-cased (counts as a mismatch via plain `!=`), no indels
     (`has_indel` is always False), emission order is all-PLUS-then-all-MINUS, and there
     is NO dedupe (a palindromic primer yields duplicate entries on both strands).

No fixtures/dependencies: every test builds its reference/primer strings inline so the
offset arithmetic is easy to verify by eye.
"""
from __future__ import annotations

from pangenome_primer.binding import find_binding_sites_naive, revcomp
from pangenome_primer.model import Strand


def _plus(sites):
    return [s for s in sites if s.strand is Strand.PLUS]


def _minus(sites):
    return [s for s in sites if s.strand is Strand.MINUS]


# --- clause 1: chrom/start/end are haplotype coordinates; basic PLUS/MINUS fields ------


def test_plus_site_basic_fields():
    primer = "ACGTACGTAC"  # 10 bp
    ref = "TTTT" + primer + "GGGG"  # site at [4, 14)
    sites = find_binding_sites_naive("p1", primer, ref, "hapX", "chr1", max_mismatches=0)
    plus = _plus(sites)
    assert len(plus) == 1
    s = plus[0]
    assert s.primer_name == "p1"
    assert s.haplotype_id == "hapX"
    assert s.chrom == "chr1"
    assert (s.start, s.end) == (4, 14)
    assert s.mismatches == 0
    assert s.mismatch_offsets_3p == []
    assert s.has_indel is False


def test_minus_site_basic_fields():
    primer = "ACGTACGTAC"
    ref = "TTTT" + revcomp(primer) + "GGGG"  # site at [4, 14)
    sites = find_binding_sites_naive("p1", primer, ref, "hapX", "chr1", max_mismatches=0)
    minus = _minus(sites)
    assert len(minus) == 1
    s = minus[0]
    assert (s.start, s.end) == (4, 14)
    assert s.mismatches == 0
    assert s.mismatch_offsets_3p == []


# --- clause 2: end - start == len(primer) for every returned site ---------------------


def test_end_minus_start_equals_primer_length_for_every_site():
    primer = "ACGTACGTAC"
    ref = primer + "NNNN" + revcomp(primer)
    sites = find_binding_sites_naive("p", primer, ref, "h", "c", max_mismatches=0)
    assert len(sites) == 2  # sanity: found both
    assert all(s.end - s.start == len(primer) for s in sites)


# --- clause 3: 3'-offset arithmetic, verified independently on both strands -----------


def test_plus_offset_from_mismatch_at_known_primer_index():
    # primer index k -> offset (L-1-k); mismatch at index 2 of a 10bp primer -> offset 7.
    primer = "AAAAAAAAAA"
    window = list(primer)
    window[2] = "G"
    ref = "".join(window)
    sites = find_binding_sites_naive("p", primer, ref, "h", "c", max_mismatches=1)
    plus = _plus(sites)
    assert len(plus) == 1
    assert plus[0].mismatches == 1
    assert plus[0].mismatch_offsets_3p == [7]


def test_minus_offset_zero_is_genomic_position_start():
    """For a MINUS site, offset 0 (the primer's 3' base) corresponds to genomic position
    `start` -- the FIRST base of the top-strand window, not the last."""
    primer = "AAAAAAAAAA"
    window = list(revcomp(primer))  # all T's
    window[0] = "C"  # mutate genomic position `start`
    ref = "".join(window)
    sites = find_binding_sites_naive("p", primer, ref, "h", "c", max_mismatches=1)
    minus = _minus(sites)
    assert len(minus) == 1
    assert minus[0].start == 0
    assert minus[0].mismatch_offsets_3p == [0]


def test_minus_offset_L_minus_1_is_genomic_position_end_minus_1():
    """The last base of the top-strand window (genomic position end-1) is the primer's
    5'-most base -> offset L-1."""
    primer = "AAAAAAAAAA"
    window = list(revcomp(primer))
    window[-1] = "C"  # mutate genomic position end-1
    ref = "".join(window)
    sites = find_binding_sites_naive("p", primer, ref, "h", "c", max_mismatches=1)
    minus = _minus(sites)
    assert len(minus) == 1
    assert minus[0].mismatch_offsets_3p == [len(primer) - 1]


def test_offsets_sorted_ascending_multi_mismatch():
    primer = "AAAAAAAAAA"
    window = list(primer)
    window[7] = "C"  # offset (9-7) = 2
    window[1] = "C"  # offset (9-1) = 8
    ref = "".join(window)
    sites = find_binding_sites_naive("p", primer, ref, "h", "c", max_mismatches=2)
    plus = _plus(sites)
    assert len(plus) == 1
    assert plus[0].mismatch_offsets_3p == [2, 8]


# --- mismatch counting exactly at the max_mismatches cap boundary ---------------------


def _three_mismatch_ref(primer: str) -> str:
    window = list(primer)
    for i in (0, 1, 2):
        window[i] = "C" if window[i] != "C" else "G"
    return "".join(window)


def test_mismatch_cap_boundary_accepted_at_exact_max():
    primer = "AAAAAAAAAA"
    ref = _three_mismatch_ref(primer)
    sites = find_binding_sites_naive("p", primer, ref, "h", "c", max_mismatches=3)
    plus = _plus(sites)
    assert len(plus) == 1
    assert plus[0].mismatches == 3


def test_mismatch_cap_boundary_rejected_above_max():
    primer = "AAAAAAAAAA"
    ref = _three_mismatch_ref(primer)  # 3 actual mismatches
    sites = find_binding_sites_naive("p", primer, ref, "h", "c", max_mismatches=2)
    plus = _plus(sites)
    assert plus == []  # cap is 2; 3 mismatches must be rejected


# --- clause 6: N handling, never special-cased -----------------------------------------


def test_n_in_reference_counts_as_mismatch():
    primer = "ACGTACGTAC"
    window = list(primer)
    window[5] = "N"
    ref = "".join(window)
    sites = find_binding_sites_naive("p", primer, ref, "h", "c", max_mismatches=1)
    plus = _plus(sites)
    assert len(plus) == 1
    assert plus[0].mismatches == 1
    assert plus[0].mismatch_offsets_3p == [4]  # (10-1)-5


def test_n_in_primer_counts_as_mismatch():
    primer = "ACGTNCGTAC"  # N at index 4
    ref = "ACGTACGTAC"  # a real base at the same position
    sites = find_binding_sites_naive("p", primer, ref, "h", "c", max_mismatches=1)
    plus = _plus(sites)
    assert len(plus) == 1
    assert plus[0].mismatches == 1


def test_n_vs_n_is_literal_equality_not_a_mismatch():
    """Documents current behavior rather than asserting it is correct: mismatch counting
    is plain string `!=`, so an N in the reference aligned against an N in the primer at
    the same position is NOT flagged as a mismatch -- there is no biological ambiguity
    handling. This only matters if a primer sequence itself ever contains 'N'."""
    primer = "ACGTNCGTAC"
    ref = "ACGTNCGTAC"  # identical, including the N at the same position
    sites = find_binding_sites_naive("p", primer, ref, "h", "c", max_mismatches=0)
    plus = _plus(sites)
    assert len(plus) == 1
    assert plus[0].mismatches == 0


# --- clause 6: no indels ----------------------------------------------------------------


def test_has_indel_always_false():
    primer = "ACGTACGTAC"
    ref = primer + "NNNN" + revcomp(primer)
    sites = find_binding_sites_naive("p", primer, ref, "h", "c", max_mismatches=0)
    assert sites  # sanity
    assert all(s.has_indel is False for s in sites)


# --- clause 6: emission order is all-PLUS-then-all-MINUS --------------------------------


def test_emission_order_all_plus_then_all_minus():
    primer = "ACGTACGTAC"
    # two PLUS sites bracketing one MINUS site in reference order, to prove the *output*
    # order is strand-grouped rather than position-ordered.
    ref = primer + "NNNN" + revcomp(primer) + "NNNN" + primer
    sites = find_binding_sites_naive("p", primer, ref, "h", "c", max_mismatches=0)
    strands = [s.strand for s in sites]
    assert strands.count(Strand.PLUS) == 2
    assert strands.count(Strand.MINUS) == 1
    first_minus = strands.index(Strand.MINUS)
    assert all(st is Strand.PLUS for st in strands[:first_minus])
    assert all(st is Strand.MINUS for st in strands[first_minus:])


# --- clause 6: no dedupe -- a palindromic primer duplicates on both strands ------------


def test_palindromic_primer_duplicates_on_both_strands():
    primer = "ACGT"
    assert revcomp(primer) == primer  # self-complementary => palindromic
    ref = "TT" + primer + "TT"
    sites = find_binding_sites_naive("p", primer, ref, "h", "c", max_mismatches=0)
    assert len(sites) == 2  # current no-dedupe behavior: one entry per strand
    assert {s.strand for s in sites} == {Strand.PLUS, Strand.MINUS}
    assert all((s.start, s.end) == (2, 6) for s in sites)


# --- boundary conditions: site at position 0 and at the very end of the reference -----


def test_site_at_start_of_reference():
    primer = "ACGTACGTAC"
    ref = primer + "NNNNNNNNNN"
    sites = find_binding_sites_naive("p", primer, ref, "h", "c", max_mismatches=0)
    plus = _plus(sites)
    assert len(plus) == 1
    assert plus[0].start == 0
    assert plus[0].end == len(primer)


def test_site_at_end_of_reference():
    primer = "ACGTACGTAC"
    ref = "NNNNNNNNNN" + primer
    sites = find_binding_sites_naive("p", primer, ref, "h", "c", max_mismatches=0)
    plus = _plus(sites)
    assert len(plus) == 1
    assert plus[0].end == len(ref)
    assert plus[0].start == len(ref) - len(primer)


def test_minus_site_at_start_and_end_of_reference():
    primer = "ACGTACGTAC"
    rc = revcomp(primer)
    ref = rc + "NNNNNNNNNN" + rc
    sites = find_binding_sites_naive("p", primer, ref, "h", "c", max_mismatches=0)
    minus = _minus(sites)
    assert len(minus) == 2
    starts = sorted(s.start for s in minus)
    assert starts == [0, len(ref) - len(primer)]
