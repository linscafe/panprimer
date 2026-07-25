"""Backend-conformance suite (Phase 0 safety net).

Every binding-search backend must satisfy the contract in `binding.py`'s module docstring
(the backend contract): haplotype-genome coordinates,
`end - start == len(primer)`, 3'-end-relative `mismatch_offsets_3p`, correct `strand`, and
(for backends that only prune, like `bwa_backend`) it is fine to return a SUPERSET of the
true positions since `binding.find_binding_sites_naive` re-scores every candidate window.

This file runs one battery of assertions against every backend in `BACKENDS` below, using
the synthetic mini-genome from `fixtures/mini_genome.py` as ground truth. Adding a new
backend (rust `pgp_scan`, an anchor-lift backend, a pangenome index) is meant to be one
entry in that list; a backend whose dependency isn't available (or isn't wired up yet)
`pytest.skip`s cleanly instead of failing, so this file does not need to change when a
platform lacks a wheel or a phase hasn't landed.

`naive` is the only backend live today; it is also the REFERENCE implementation every
other backend's output must be validated against once more backends land.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent))
from fixtures.mini_genome import MiniGenome, build_mini_genome  # noqa: E402

from pangenome_primer.model import BindingSite  # noqa: E402

MAX_MISMATCHES = 3


# --- backend registry -------------------------------------------------------------------
# Each entry is (name, loader). `loader()` returns a callable
#   run(seqs, fasta_path, haplotype_id, max_mismatches) -> dict[str, list[BindingSite]]
# matching `bwa_backend.find_binding_sites_batch`'s signature/return shape, or raises
# `pytest.skip.Exception` (via `pytest.skip(...)`) if the backend isn't available.


def _load_naive():
    """Wrap `find_binding_sites_naive` (in-memory, per-string) to the batch/whole-genome
    signature every backend in this registry is expected to expose, by scanning every
    contig in the FASTA for every query sequence."""
    import pysam

    from pangenome_primer.binding import find_binding_sites_naive

    def run(seqs, fasta_path, haplotype_id, max_mismatches):
        fa = pysam.FastaFile(fasta_path)
        try:
            out: dict[str, list[BindingSite]] = {s: [] for s in seqs}
            for chrom in fa.references:
                ref = fa.fetch(chrom)
                for s in seqs:
                    out[s].extend(
                        find_binding_sites_naive(s, s, ref, haplotype_id, chrom, max_mismatches)
                    )
            return out
        finally:
            fa.close()

    return run


def _load_rust():
    """The compiled scanner, via its Python wrapper rather than the raw extension: the
    wrapper is what `verify.py`/`cli.py` actually call, so it is what needs conforming.

    Both failure modes must SKIP, not error: the wrapper module may not exist yet (a commit
    before the backend landed) and the compiled extension may not be built (a checkout
    without a wheel). Erroring on either would make this file's "adding a backend is one
    list entry" promise false for anyone who hasn't built the Rust side.
    """
    try:
        from pangenome_primer import rust_backend
    except ImportError:
        pytest.skip("pangenome_primer.rust_backend is not present in this checkout")

    if not rust_backend.available():
        pytest.skip("compiled pangenome_primer._scan extension is not installed")

    def run(seqs, fasta_path, haplotype_id, max_mismatches):
        return rust_backend.find_binding_sites_batch(
            seqs, fasta_path, haplotype_id, max_mismatches
        )

    return run


BACKENDS: list[tuple[str, "object"]] = [
    ("naive", _load_naive),
    ("rust", _load_rust),
]


@pytest.fixture(scope="session")
def mini_genome(tmp_path_factory) -> MiniGenome:
    tmp = tmp_path_factory.mktemp("mini_genome_conformance")
    return build_mini_genome(tmp)


def _sites_at(sites: list[BindingSite], planted) -> list[BindingSite]:
    return [s for s in sites if s.strand is planted.strand and s.start == planted.start]


@pytest.fixture(
    scope="module",
    params=[name for name, _ in BACKENDS],
    ids=[name for name, _ in BACKENDS],
)
def backend(request):
    """Yields (name, run_callable) for a live backend; skips cleanly if unavailable."""
    name = request.param
    loader = dict(BACKENDS)[name]
    run = loader()  # may call pytest.skip() internally
    return name, run


@pytest.fixture(scope="module")
def result(backend, mini_genome) -> tuple[str, dict[str, list[BindingSite]]]:
    """Runs the backend exactly ONCE per module (all planted-site primers, whole genome)
    and shares the output across every assertion test below -- naive's O(N*L) pure-Python
    scan makes re-running per test function prohibitively slow (see mini_genome.py's
    sizing note)."""
    name, run = backend
    seqs = sorted({p.primer_seq for p in mini_genome.planted})
    out = run(seqs, mini_genome.fasta_gz, mini_genome.haplotype_id, MAX_MISMATCHES)
    return name, out


# --- the shared battery ------------------------------------------------------------------


def test_finds_every_planted_site(result, mini_genome):
    """Every planted site must appear in the result (superset semantics: extra candidate
    hits are fine, missing the true positive is not)."""
    name, out = result
    for planted in mini_genome.planted:
        found = _sites_at(out[planted.primer_seq], planted)
        assert found, f"[{name}] missing planted site {planted.label}"


def test_site_shape_and_coordinates(result, mini_genome):
    """chrom/start/end are haplotype coordinates and end-start == primer length for every
    returned site touching a planted locus (clauses 1 and 2 of the contract)."""
    name, out = result
    for planted in mini_genome.planted:
        for s in _sites_at(out[planted.primer_seq], planted):
            assert s.chrom == planted.chrom
            assert s.end - s.start == len(planted.primer_seq), (
                f"[{name}] {planted.label}: end-start={s.end - s.start} != "
                f"{len(planted.primer_seq)}"
            )


def test_strand_and_mismatch_offsets(result, mini_genome):
    """strand and mismatch_offsets_3p (3'-end convention) match the planted ground truth
    (clauses 3 and 4). Superset backends may return extra low-confidence hits elsewhere,
    but the planted-site record itself must be exactly right, since classify.py trusts
    this field directly."""
    name, out = result
    for planted in mini_genome.planted:
        found = _sites_at(out[planted.primer_seq], planted)
        assert found, f"[{name}] missing planted site {planted.label}"
        # Multiple superset hits at the exact planted coordinate/strand shouldn't happen,
        # but if a backend ever does return >1, every one of them must agree.
        for s in found:
            assert s.strand == planted.strand
            assert s.mismatches == planted.mismatches, (
                f"[{name}] {planted.label}: mismatches={s.mismatches} != "
                f"{planted.mismatches}"
            )
            assert s.mismatch_offsets_3p == planted.expected_offsets_3p, (
                f"[{name}] {planted.label}: offsets={s.mismatch_offsets_3p} != "
                f"{planted.expected_offsets_3p}"
            )


def test_no_indels_reported(result):
    """has_indel is always False (clause 6) -- no backend in this registry does indel
    detection."""
    name, out = result
    for sites in out.values():
        assert all(s.has_indel is False for s in sites), f"[{name}] has_indel set"


def test_close_pair_sites_both_recovered_independently(result, mini_genome):
    """Two independent binding sites separated by a 2 bp gap (< slop=3) must both be
    reported at their own coordinates -- a windowed/slop-based backend must not merge or
    drop one of them just because their candidate windows overlap."""
    name, out = result
    close = [p for p in mini_genome.planted if p.label.startswith("ctgB:120")]
    assert len(close) == 2  # sanity on the fixture itself
    for planted in close:
        found = _sites_at(out[planted.primer_seq], planted)
        assert found, f"[{name}] missing close-pair site {planted.label}"
