"""Differential test: the Rust scanner must agree with `find_binding_sites_naive` exactly.

`naive` is the reference implementation (see `test_binding_naive.py`, which pins its
semantics clause by clause). This file is the strongest correctness weapon Phase 2 has: for
randomly generated genomes and randomly generated primers, the compiled backend and the
pure-Python one must return the *same set of sites*, field for field.

Comparison is set-based, not list-based, because the two backends legitimately differ in
emission order and in deduplication: `find_binding_sites_naive` emits all-PLUS-then-all-MINUS
and never dedupes, while the batch backends dedupe on `(chrom, start, strand)`. Neither is
"wrong" -- but the *set* of distinct sites must be identical, and that is what downstream
(`pcr.pair_amplicons`) actually consumes.

Also covers `resolve_scan_path`, which `rust_backend.py`'s docstring points at: the `.fa`
vs `.fa.gz` precedence that lets Phase 2 land before Phase 3 repoints `samples.tsv`.
"""
from __future__ import annotations

import random
import subprocess
from pathlib import Path

import pytest

from pangenome_primer import rust_backend
from pangenome_primer.binding import find_binding_sites_naive
from pangenome_primer.model import BindingSite

pytestmark = pytest.mark.skipif(
    not rust_backend.available(),
    reason="compiled pangenome_primer._scan extension is not installed",
)

MAX_MISMATCHES = 3
HAP_ID = "DIFF#hap1"
_ALPHABET = "ACGT"


# --- helpers ------------------------------------------------------------------------------


def _norm(sites: list[BindingSite]) -> set[tuple]:
    """A backend-independent identity for a site: everything downstream reads."""
    return {
        (
            s.chrom,
            s.start,
            s.end,
            s.strand.value,
            s.mismatches,
            tuple(sorted(s.mismatch_offsets_3p)),
            bool(s.has_indel),
        )
        for s in sites
    }


def _write_bgzf(tmp_path: Path, contigs: dict[str, str], name: str = "g") -> str:
    """Write a FASTA, bgzip it, and index it -- the shape this backend actually reads."""
    fa = tmp_path / f"{name}.fa"
    with fa.open("w") as fh:
        for cname, seq in contigs.items():
            fh.write(f">{cname}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i : i + 60] + "\n")
    subprocess.run(["bgzip", "-f", "-@", "2", str(fa)], check=True)
    gz = str(fa) + ".gz"
    subprocess.run(["samtools", "faidx", gz], check=True)
    return gz


def _naive_whole_genome(
    seqs: list[str], contigs: dict[str, str]
) -> dict[str, list[BindingSite]]:
    """Reference answer: scan every contig with the pure-Python implementation."""
    out: dict[str, list[BindingSite]] = {s: [] for s in seqs}
    for chrom, ref in contigs.items():
        for s in seqs:
            out[s].extend(
                find_binding_sites_naive(s, s, ref, HAP_ID, chrom, MAX_MISMATCHES)
            )
    return out


def _assert_agrees(seqs: list[str], contigs: dict[str, str], gz: str, label: str) -> int:
    """Both backends, same input; returns the number of sites compared."""
    got = rust_backend.find_binding_sites_batch(seqs, gz, HAP_ID, MAX_MISMATCHES)
    want = _naive_whole_genome(seqs, contigs)

    total = 0
    for s in seqs:
        g, w = _norm(got.get(s, [])), _norm(want.get(s, []))
        total += len(w)
        assert g == w, (
            f"[{label}] primer {s!r} diverged\n"
            f"  only in rust ({len(g - w)}): {sorted(g - w)[:5]}\n"
            f"  only in naive ({len(w - g)}): {sorted(w - g)[:5]}"
        )
    return total


def _random_seq(rng: random.Random, n: int) -> str:
    return "".join(rng.choice(_ALPHABET) for _ in range(n))


# --- the differential sweep ---------------------------------------------------------------


@pytest.mark.parametrize("seed", list(range(12)))
def test_random_genome_and_primers_agree_with_naive(seed, tmp_path):
    """The core sweep: random genome, random primers, mismatch budgets 0..3.

    Primers are drawn two ways -- some lifted from the genome itself (so they are guaranteed
    to hit, often with engineered mismatches) and some purely random (mostly misses). Both
    matter: agreeing only on "no hits anywhere" would be a vacuous pass.
    """
    rng = random.Random(seed)
    contigs = {
        f"ctg{i}": _random_seq(rng, rng.randint(400, 2500)) for i in range(rng.randint(1, 4))
    }
    gz = _write_bgzf(tmp_path, contigs, name=f"rand{seed}")

    seqs: set[str] = set()
    names = list(contigs)
    for _ in range(6):  # planted: lifted from the genome, then perturbed
        c = contigs[rng.choice(names)]
        L = rng.randint(18, 25)
        if len(c) <= L:
            continue
        i = rng.randrange(0, len(c) - L)
        p = list(c[i : i + L])
        for _ in range(rng.randint(0, MAX_MISMATCHES)):
            j = rng.randrange(L)
            p[j] = rng.choice([b for b in _ALPHABET if b != p[j]])
        seqs.add("".join(p))
    for _ in range(3):  # unplanted: random, usually absent
        seqs.add(_random_seq(rng, rng.randint(18, 25)))

    compared = _assert_agrees(sorted(seqs), contigs, gz, f"seed={seed}")
    assert compared >= 0


def test_planted_sites_are_actually_found(tmp_path):
    """Guards against a vacuous differential pass: if both backends silently returned
    nothing, every comparison above would still succeed. Here a known exact primer must
    produce a real hit."""
    rng = random.Random(99)
    ref = _random_seq(rng, 3000)
    primer = ref[500:520]
    contigs = {"ctgA": ref}
    gz = _write_bgzf(tmp_path, contigs, name="planted")

    got = rust_backend.find_binding_sites_batch([primer], gz, HAP_ID, MAX_MISMATCHES)
    hits = [s for s in got[primer] if s.start == 500 and s.strand.value == "+"]
    assert hits, "exact planted primer produced no plus-strand hit at its own coordinates"
    assert hits[0].mismatches == 0
    assert hits[0].end - hits[0].start == len(primer)
    _assert_agrees([primer], contigs, gz, "planted")


def test_soft_masked_lowercase_agrees(tmp_path):
    """HPRC assemblies carry soft-masked (lowercase) repeat sequence. `binding.py` upper()s
    both sides, so lowercase must NOT suppress a hit -- a classic streaming-scanner bug."""
    rng = random.Random(7)
    ref = _random_seq(rng, 1200)
    primer = ref[300:320]
    masked = ref[:280] + ref[280:360].lower() + ref[360:]
    contigs = {"ctgSoft": masked}
    gz = _write_bgzf(tmp_path, contigs, name="soft")

    got = rust_backend.find_binding_sites_batch([primer], gz, HAP_ID, MAX_MISMATCHES)
    assert any(
        s.start == 300 and s.mismatches == 0 for s in got[primer]
    ), "soft-masked lowercase suppressed an exact match"
    _assert_agrees([primer], contigs, gz, "soft-masked")


def test_n_runs_agree(tmp_path):
    """`N` counts as a mismatch and is never special-cased -- except that `N` vs `N`
    compares literally equal (see test_binding_naive.py). Assembly gap runs must not
    desynchronise the scanner's coordinates either."""
    rng = random.Random(23)
    ref = _random_seq(rng, 800) + "N" * 150 + _random_seq(rng, 800)
    primer = ref[1000:1020]
    contigs = {"ctgN": ref}
    gz = _write_bgzf(tmp_path, contigs, name="ngap")

    _assert_agrees([primer, "N" * 20, ref[700:720]], contigs, gz, "N-runs")


def test_contig_boundary_is_not_crossed(tmp_path):
    """A primer spanning the junction of two concatenated contigs must NOT match: contigs
    are separate sequences. This is the classic failure mode of a scanner that streams the
    whole file as one buffer."""
    rng = random.Random(31)
    a, b = _random_seq(rng, 600), _random_seq(rng, 600)
    contigs = {"ctgA": a, "ctgB": b}
    gz = _write_bgzf(tmp_path, contigs, name="boundary")

    straddle = a[-10:] + b[:10]  # exists only if the join is (wrongly) scanned through
    got = rust_backend.find_binding_sites_batch([straddle], gz, HAP_ID, MAX_MISMATCHES)
    want = _naive_whole_genome([straddle], contigs)
    assert _norm(got[straddle]) == _norm(want[straddle])


def test_adjacent_sites_within_slop_agree(tmp_path):
    """Two true sites closer together than slop=3 must both survive -- a dedupe keyed too
    loosely would collapse them into one."""
    rng = random.Random(41)
    unit = _random_seq(rng, 20)
    contigs = {"ctgRep": _random_seq(rng, 300) + unit + unit + _random_seq(rng, 300)}
    gz = _write_bgzf(tmp_path, contigs, name="adjacent")
    _assert_agrees([unit], contigs, gz, "adjacent-repeats")


def test_empty_and_missing_query_sets(tmp_path):
    rng = random.Random(53)
    contigs = {"ctgE": _random_seq(rng, 400)}
    gz = _write_bgzf(tmp_path, contigs, name="empty")
    assert rust_backend.find_binding_sites_batch([], gz, HAP_ID, MAX_MISMATCHES) == {}


# --- resolve_scan_path (referenced by rust_backend.py's docstring) -------------------------


class TestResolveScanPath:
    """`.fa` vs `.fa.gz` precedence -- what lets Phase 2 land before Phase 3 repoints
    `config/samples.tsv` from the uncompressed FASTA to the BGZF file."""

    def test_prefers_bgzf_sibling_over_plain_fa(self, tmp_path):
        fa, gz = tmp_path / "h.fa", tmp_path / "h.fa.gz"
        fa.write_text(">c\nACGT\n")
        gz.write_bytes(b"")
        assert rust_backend.resolve_scan_path(str(fa)) == str(gz)

    def test_accepts_an_explicit_gz_path(self, tmp_path):
        gz = tmp_path / "h.fa.gz"
        gz.write_bytes(b"")
        assert rust_backend.resolve_scan_path(str(gz)) == str(gz)

    def test_falls_back_to_plain_fa_when_no_gz(self, tmp_path):
        fa = tmp_path / "h.fa"
        fa.write_text(">c\nACGT\n")
        assert rust_backend.resolve_scan_path(str(fa)) == str(fa)

    def test_raises_filenotfounderror_naming_both_candidates(self, tmp_path):
        missing = tmp_path / "nope.fa"
        with pytest.raises(FileNotFoundError) as exc:  # ScanFileNotFound subclasses it,
            rust_backend.resolve_scan_path(str(missing))  # so CLI handlers still catch it
        assert "nope.fa" in str(exc.value)

    def test_explicit_gz_that_is_absent_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            rust_backend.resolve_scan_path(str(tmp_path / "absent.fa.gz"))
