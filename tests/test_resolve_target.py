"""Tests for target resolution, including the new GRCh38-coordinate input path."""
from __future__ import annotations

import random

import pytest

from pangenome_primer import project


def test_parse_coords():
    assert project.parse_coords("chr12:6,544,868-6,548,730") == ("chr12", 6544868, 6548730)
    assert project.parse_coords("HG002#1#ctg:10-20") == ("HG002#1#ctg", 10, 20)
    assert project.parse_coords("not-coords.fa") is None
    assert project.parse_coords("chr1:100") is None


def test_chm13_coords_passthrough():
    # CHM13 coords need no alignment (no mappy required)
    assert project.resolve_target("chr7:100-200", "unused.fa", assembly="chm13") == ("chr7", 100, 200)


def test_grch38_requires_reference():
    with pytest.raises(ValueError, match="requires --grch38"):
        project.resolve_target("chr7:100-200", "unused.fa", assembly="grch38")


def test_missing_file_target():
    with pytest.raises(FileNotFoundError):
        project.resolve_target("/no/such/locus.fa", "unused.fa")


def _write_fa(path, name, seq):
    import pysam

    path.write_text(f">{name}\n{seq}\n")
    pysam.faidx(str(path))
    return str(path)


def test_grch38_coords_anchor_to_chm13(tmp_path):
    """GRCh38 chr:start-end -> slice extracted -> aligned to CHM13 -> CHM13 coords."""
    pytest.importorskip("mappy")
    pytest.importorskip("pysam")
    random.seed(7)
    chm = "".join(random.choice("ACGT") for _ in range(3000))
    chm13 = _write_fa(tmp_path / "chm13.fa", "c1", chm)
    # GRCh38 "locus" is a 600 bp window that corresponds to chm[1000:1600]
    grch38 = _write_fa(tmp_path / "grch38.fa", "g1", chm[1000:1600])

    chrom, start, end = project.resolve_target(
        "g1:0-600", chm13, assembly="grch38", grch38_fasta=grch38)
    assert chrom == "c1"
    assert abs(start - 1000) <= 5 and abs(end - 1600) <= 5


def test_empty_mmi_is_ignored(tmp_path):
    """A 0-byte/truncated .mmi (interrupted `minimap2 -d`) must be skipped, not loaded as an
    empty aligner that maps nothing (would silently turn every projection uncertain)."""
    pytest.importorskip("mappy")
    pytest.importorskip("pysam")
    random.seed(11)
    seq = "".join(random.choice("ACGT") for _ in range(3000))
    fa = _write_fa(tmp_path / "hap.fa", "h1", seq)
    (tmp_path / "hap.fa.mmi").write_bytes(b"")  # corrupt/empty index next to the FASTA
    aln = project._aligner(fa, preset="asm5")  # must fall back to building from the FASTA
    hits = aln.map(seq[1000:1200])
    assert any(True for _ in hits), "aligner built from FASTA should map the query"
