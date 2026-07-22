"""Pangenome PCR primer design — design universal primers against the HPRC R2 pangenome.

See CONTEXT.md for the ubiquitous language and the plan for architecture/decisions.
"""
from __future__ import annotations

__version__ = "0.1.0"

from .model import (
    Amplicon,
    BindingSite,
    Haplotype,
    HaplotypeResult,
    HaplotypeStatus,
    Locus,
    PairResult,
    Primer,
    PrimerPair,
    Strand,
)

__all__ = [
    "Amplicon",
    "BindingSite",
    "Haplotype",
    "HaplotypeResult",
    "HaplotypeStatus",
    "Locus",
    "PairResult",
    "Primer",
    "PrimerPair",
    "Strand",
    "__version__",
]
