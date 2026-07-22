"""A synthetic mini-pangenome for smoke-testing the whole pipeline with no downloads.

Builds a handful of haplotype contexts engineered to land in each status, plus a primer
pair. Used by `pangenome-primer selftest` and available to tests.
"""
from __future__ import annotations

from .binding import revcomp
from .engine import HaplotypeContext
from .model import Locus, Primer, PrimerPair

F = "TTGCACAGTCCAGATTGCAA"
R = "CATTGCGATTGACACTTGCG"
MID = "GATACCATGCTGACGTTAAC" * 5
LEFT = "AACCGGTTAACCGGTTAACC"
RIGHT = "TTGGCCAATTGGCCAATTGG"
SPACER = "ACGACTACGT" * 50

REGION = LEFT + F + MID + revcomp(R) + RIGHT
F_START = len(LEFT)
PRODUCT = 20 + len(MID) + 20


def _mut_3p(seq: str, at: int) -> str:
    pos = at + len(F) - 1
    b = seq[pos]
    return seq[:pos] + ("C" if b != "C" else "G") + seq[pos + 1 :]


def demo_pair() -> PrimerPair:
    return PrimerPair(
        "demo", Primer("demo_F", F), Primer("demo_R", R), product_size_chm13=PRODUCT
    )


def demo_contexts() -> list[HaplotypeContext]:
    dead = _mut_3p(REGION, F_START)
    return [
        HaplotypeContext("HG_pass#hap1", "chr1", REGION, Locus("HG_pass#hap1", "chr1", 0, len(REGION))),
        HaplotypeContext("HG_drop#hap1", "chr1", dead, Locus("HG_drop#hap1", "chr1", 0, len(REGION))),
        HaplotypeContext("HG_off#hap1", "chr1", dead + SPACER + REGION, Locus("HG_off#hap1", "chr1", 0, len(dead))),
        HaplotypeContext("HG_multi#hap1", "chr1", REGION + SPACER + REGION, Locus("HG_multi#hap1", "chr1", 0, len(REGION))),
        HaplotypeContext("HG_unc#hap1", "chr1", REGION, None, projection_ok=False),
    ]
