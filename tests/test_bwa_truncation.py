"""bwa's XA cap silently truncates repeat primers; this pins the detection.

The bug: `bwa samse -n N` does not emit a truncated XA list for a read with more than N
hits -- it omits the XA tag entirely. Parsing primary+XA then yields exactly ONE position,
with nothing to distinguish "this primer binds in one place" from "this primer binds in
330,000 places and we were told about one of them".

Downstream that lone site has no partner primer nearby, so the pair is reported as a
**dropout** -- a confident wrong answer. That is not hypothetical: the demo's original
CYP2D6_dropout forward primer is an Alu consensus that bwa reports as `X0=12760 X1=317208`
with no XA tag, and it sat in the demo as the dropout exemplar.

The fix reads X0 (best hits) and X1 (suboptimal), which bwa still reports when it suppresses
XA. Their sum matches the exhaustive scanner exactly (12760+317208 = 329,968 for the Alu
primer; 3+113 = 116 for a well-behaved one), so it is a sound signal.

These tests use synthetic SAM rather than a real 3 Gb index so they run in milliseconds.
"""
from __future__ import annotations

import warnings

import pytest

from pangenome_primer import bwa_backend


def _sam(qi: int, rname: str, pos: int, tags: str) -> str:
    """One minimal mapped SAM record for query index `qi`."""
    return f"{qi}\t0\t{rname}\t{pos}\t0\t20M\t*\t0\t0\t*\t*\t{tags}"


def _parse(sam_lines, n_queries):
    """Drive the parser over canned SAM by faking the subprocess boundary."""
    return sam_lines, n_queries


class TestTruncationDetection:
    """Exercises the X0/X1 arithmetic in `_candidate_positions_batch`'s parse loop."""

    def _run(self, monkeypatch, sam_body: str, seqs: list[str]):
        import subprocess

        class _R:
            def __init__(self, stdout=""):
                self.stdout = stdout

        def fake_run(cmd, **kw):
            if cmd[1] == "samse":
                return _R(sam_body)
            return _R("")

        monkeypatch.setattr(subprocess, "run", fake_run)
        monkeypatch.setattr(bwa_backend, "subprocess", subprocess)
        return bwa_backend._candidate_positions_batch(seqs, "/fake.fa", 3)

    def test_missing_xa_with_huge_x0_is_flagged(self, monkeypatch):
        """The exact Alu case: X0=12760, X1=317208, no XA tag, one parsed position."""
        sam = _sam(0, "ctg", 100, "X0:i:12760\tX1:i:317208")
        cands, trunc = self._run(monkeypatch, sam, ["A" * 20])
        assert len(cands[0]) == 1, "only the primary hit is recoverable without XA"
        assert trunc == {0: 329968}, (
            "a primer with 329,968 hits reported as one site MUST be flagged; "
            "silently returning it is the dropout bug"
        )

    def test_complete_xa_list_is_not_flagged(self, monkeypatch):
        """X0+X1 == parsed positions => nothing was dropped."""
        xa = "XA:Z:" + "".join(f"ctg,+{200 + i},20M,0;" for i in range(2))
        sam = _sam(0, "ctg", 100, f"X0:i:1\tX1:i:2\t{xa}")
        cands, trunc = self._run(monkeypatch, sam, ["A" * 20])
        assert len(cands[0]) == 3
        assert trunc == {}, "a fully-reported hit list must not be flagged"

    def test_flags_only_the_offending_query_in_a_batch(self, monkeypatch):
        """One repeat primer must not taint its well-behaved batch-mates."""
        good_xa = "XA:Z:ctg,+400,20M,0;"
        sam = "\n".join([
            _sam(0, "ctg", 100, "X0:i:12760\tX1:i:317208"),          # truncated
            _sam(1, "ctg", 300, f"X0:i:1\tX1:i:1\t{good_xa}"),        # complete
        ])
        cands, trunc = self._run(monkeypatch, sam, ["A" * 20, "C" * 20])
        assert set(trunc) == {0}
        assert len(cands[1]) == 2

    def test_unmapped_read_is_skipped_not_flagged(self, monkeypatch):
        sam = "0\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\tX0:i:0"
        cands, trunc = self._run(monkeypatch, sam, ["A" * 20])
        assert cands[0] == set()
        assert trunc == {}

    def test_absent_x0_x1_tags_do_not_false_positive(self, monkeypatch):
        """Older/other bwa output without X0/X1 must not be reported as truncated."""
        sam = _sam(0, "ctg", 100, "NM:i:0")
        cands, trunc = self._run(monkeypatch, sam, ["A" * 20])
        assert len(cands[0]) == 1
        assert trunc == {}, "no X0/X1 means no evidence of truncation, not evidence of it"


class TestTruncationIsNeverSilent:
    """The defining property of the bug was silence. A caller that ignores the new
    `truncated` parameter must still be told."""

    def test_warns_even_when_caller_does_not_opt_in(self, monkeypatch):
        monkeypatch.setattr(bwa_backend, "ensure_index", lambda fasta: None)
        monkeypatch.setattr(
            bwa_backend, "_candidate_positions_batch",
            lambda seqs, fasta, mm: ({0: {("ctg", 100)}}, {0: 329968}),
        )

        class _FA:
            def get_reference_length(self, c):
                return 1000

            def fetch(self, c, s, e):
                return "A" * (e - s)

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            bwa_backend.find_binding_sites_batch(
                ["A" * 20], "/fake.fa", "H#hap1", 3, fa=_FA()
            )
        assert any(issubclass(x.category, RuntimeWarning) for x in w), (
            "truncation must warn even without the `truncated` dict -- silence is the bug"
        )
        msg = str(w[0].message)
        assert "329968" in msg and "INCOMPLETE" in msg

    def test_populates_the_truncated_dict_when_supplied(self, monkeypatch):
        monkeypatch.setattr(bwa_backend, "ensure_index", lambda fasta: None)
        monkeypatch.setattr(
            bwa_backend, "_candidate_positions_batch",
            lambda seqs, fasta, mm: ({0: {("ctg", 100)}}, {0: 329968}),
        )

        class _FA:
            def get_reference_length(self, c):
                return 1000

            def fetch(self, c, s, e):
                return "A" * (e - s)

        got: dict[str, int] = {}
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            bwa_backend.find_binding_sites_batch(
                ["A" * 20], "/fake.fa", "H#hap1", 3, fa=_FA(), truncated=got
            )
        assert got == {"A" * 20: 329968}


def test_rust_backend_accepts_truncated_and_never_populates_it():
    """Contract parity: the exhaustive backend takes the same kwarg and leaves it empty."""
    rust_backend = pytest.importorskip("pangenome_primer.rust_backend")
    if not rust_backend.available():
        pytest.skip("compiled extension not installed")
    import inspect

    sig = inspect.signature(rust_backend.find_binding_sites_batch)
    assert "truncated" in sig.parameters, (
        "both backends must accept `truncated` or the dispatcher cannot pass it through"
    )
