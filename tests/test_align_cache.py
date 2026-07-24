"""Unit tests for the PAF-based projection lift (item 1). No minimap2/network: a synthetic
PAF and haplotype FASTA exercise the coordinate lift on both strands."""
from __future__ import annotations

import random

import pytest

from pangenome_primer import align_cache
from pangenome_primer.binding import revcomp

pytest.importorskip("pysam")


def _hap(tmp_path) -> tuple[str, str]:
    random.seed(1)
    seq = "".join(random.choice("ACGT") for _ in range(2000))
    fa = tmp_path / "hap.fa"
    fa.write_text(f">hapctg\n{seq}\n")
    import pysam

    pysam.faidx(str(fa))
    return str(fa), seq


def test_lift_plus_strand(tmp_path):
    fa, seq = _hap(tmp_path)
    # CHM13 chr1:[200,400] <-> hapctg:[500,700] on +
    paf = tmp_path / "hap.chm13.paf"
    paf.write_text("\t".join([
        "hapctg", "2000", "500", "700", "+", "chr1", "300000000", "200", "400",
        "200", "200", "60"]) + "\n")
    proj = align_cache.project_from_paf(str(paf), "chr1", 250, 300, fa)
    assert proj.locus is not None
    assert proj.locus.chrom == "hapctg"
    assert (proj.locus.start, proj.locus.end) == (550, 600)
    assert proj.haplotype_seq == seq[550:600]


def test_lift_minus_strand_revcomps(tmp_path):
    fa, seq = _hap(tmp_path)
    paf = tmp_path / "hap.chm13.paf"
    paf.write_text("\t".join([
        "hapctg", "2000", "500", "700", "-", "chr1", "300000000", "200", "400",
        "200", "200", "60"]) + "\n")
    proj = align_cache.project_from_paf(str(paf), "chr1", 250, 300, fa)
    # minus strand: 250 -> 700-(50/200*200)=650 ; 300 -> 700-100=600  => window [600,650]
    assert (proj.locus.start, proj.locus.end) == (600, 650)
    assert proj.haplotype_seq == revcomp(seq[600:650])


def test_cigar_precise_lift_beats_linear(tmp_path):
    """A 50 bp query insertion inside the block shifts coordinates: CIGAR lifting places it
    exactly, linear interpolation is off by the insertion size."""
    fa, seq = _hap(tmp_path)
    # CHM13 chr1:[200,400] <-> hapctg:[500,750] on +, CIGAR 100M 50I 100M (50 bp query insert)
    paf = tmp_path / "hap.chm13.paf"
    paf.write_text("\t".join([
        "hapctg", "2000", "500", "750", "+", "chr1", "300000000", "200", "400",
        "200", "250", "60", "cg:Z:100M50I100M"]) + "\n")
    b = align_cache._overlapping_blocks(str(paf), "chr1", 290, 310)[0]
    assert b.cigar == [(100, "M"), (50, "I"), (100, "M")]
    # target 300 sits right after the insertion
    assert align_cache._lift_cigar(b, 300) == 650   # precise: 100 M + 50 I
    assert align_cache._lift(b, 300) == 625          # linear: off by the 25 bp half-insertion
    proj = align_cache.project_from_paf(str(paf), "chr1", 290, 310, fa)
    assert proj.locus.start == 590 and proj.locus.end == 660
    assert proj.haplotype_seq == seq[590:660]


def test_no_block_is_uncertain(tmp_path):
    fa, _ = _hap(tmp_path)
    paf = tmp_path / "hap.chm13.paf"
    paf.write_text("\t".join([
        "hapctg", "2000", "500", "700", "+", "chr2", "300000000", "200", "400",
        "200", "200", "60"]) + "\n")  # different chrom
    proj = align_cache.project_from_paf(str(paf), "chr1", 250, 300, fa)
    assert proj.locus is None
