"""Sidecar resolution across the `.fa` -> `.fa.gz` repoint (Phase 3).

The failure this guards against is not an exception -- it is a silent slowdown. When
`samples.tsv` `local_path` becomes `X.fa.gz`, `project._aligner` used to look for
`X.fa.gz.mmi`, miss, and fall through to `mappy.Aligner(X.fa.gz)`, building a whole-genome
minimap2 index in-process for every haplotype. The run still produces correct answers, just
minutes per haplotype slower, so no test that only checks *output* would ever catch it.
Hence these assert on the resolved path itself.
"""
from __future__ import annotations

from pangenome_primer import align_cache
from pangenome_primer.samples import sidecar_path


class TestSidecarPath:
    def test_plain_fa_is_unchanged(self, tmp_path):
        """The pre-Phase-3 layout must keep resolving exactly as it did."""
        fa = tmp_path / "h.fa"
        fa.write_text(">c\nACGT\n")
        (tmp_path / "h.fa.mmi").write_bytes(b"idx")
        assert sidecar_path(str(fa), ".mmi") == str(tmp_path / "h.fa.mmi")

    def test_gz_path_finds_the_fa_named_sidecar(self, tmp_path):
        """The Phase 3 case: local_path is the BGZF, the .mmi is named for the .fa."""
        gz = tmp_path / "h.fa.gz"
        gz.write_bytes(b"\x1f\x8b")
        (tmp_path / "h.fa.mmi").write_bytes(b"idx")
        assert sidecar_path(str(gz), ".mmi") == str(tmp_path / "h.fa.mmi")

    def test_exact_named_sidecar_wins_over_the_stem(self, tmp_path):
        """A genuinely .fa.gz-derived index must never be shadowed by a stale .fa one."""
        gz = tmp_path / "h.fa.gz"
        gz.write_bytes(b"\x1f\x8b")
        (tmp_path / "h.fa.mmi").write_bytes(b"stale")
        (tmp_path / "h.fa.gz.mmi").write_bytes(b"fresh")
        assert sidecar_path(str(gz), ".mmi") == str(tmp_path / "h.fa.gz.mmi")

    def test_missing_sidecar_returns_the_writer_path(self, tmp_path):
        """Nothing exists => hand back the path a writer should create, next to the input."""
        gz = tmp_path / "h.fa.gz"
        gz.write_bytes(b"\x1f\x8b")
        assert sidecar_path(str(gz), ".mmi") == str(tmp_path / "h.fa.gz.mmi")

    def test_does_not_strip_gz_from_a_non_gz_path(self, tmp_path):
        """`.gz` stripping applies only to the suffix, never mid-name."""
        fa = tmp_path / "h.gz.fa"
        fa.write_text(">c\nACGT\n")
        assert sidecar_path(str(fa), ".mmi") == str(tmp_path / "h.gz.fa.mmi")


class TestPafPathUsesIt:
    def test_paf_cache_resolves_across_the_repoint(self, tmp_path):
        """`verify.run_verify` gates on `Path(paf_path(fasta)).exists()`; if that misses, a
        built PAF cache is silently ignored and every projection reloads the .mmi."""
        gz = tmp_path / "h.fa.gz"
        gz.write_bytes(b"\x1f\x8b")
        paf = tmp_path / "h.fa.chm13.paf"
        paf.write_text("")
        assert align_cache.paf_path(str(gz)) == str(paf)
