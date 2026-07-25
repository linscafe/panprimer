"""Genome-wide primer binding search backed by the compiled Rust scanner.

This is the Phase 2 replacement for `bwa_backend.find_binding_sites_batch`. It exposes the
same signature and the same return shape (`{primer_sequence: [BindingSite, ...]}`, with
`primer_name` set to the sequence), so the three lazy import sites in `verify.py`/`cli.py`
swap between them without any other change. `search.py` is the dispatcher that picks.

Differences from the bwa path, all deliberate:

* **No index.** It streams the BGZF `.fa.gz` HPRC already distributes, so the 3.08 GB `.fa`
  and the 5.31 GB bwa index per haplotype become deletable (that is Phase 3's job, not
  this module's).
* **Exhaustive, not heuristic.** `bwa aln` prunes and `bwa samse -n 1000` caps alternative
  hits; the extension proves <= `max_mismatches` or > `max_mismatches` at every position for
  every primer in both orientations. `bwa_backend` re-scores its candidates with
  `find_binding_sites_naive` for exactly this reason -- here the exact comparator IS the
  search, so there is no re-scoring pass.
* **`slop` is inert.** It exists in the signature only because the seam has it; there are no
  candidate windows to widen. See `rust/pgp-scan/src/lib.rs`.

The extension is optional. `available()` is false when it could not be imported (no wheel
for this platform, no Rust toolchain at install time), and `search.py` degrades to `bwa`.
Nothing here raises at import.
"""
from __future__ import annotations

from pathlib import Path

from .model import BindingSite, Strand

try:  # pragma: no cover - the "absent" branch is exercised by monkeypatching in tests
    from . import _scan as _ext
except ImportError:  # pragma: no cover
    _ext = None


def available() -> bool:
    """True when the compiled `pangenome_primer._scan` extension imported successfully."""
    return _ext is not None


def version() -> str | None:
    return getattr(_ext, "__version__", None) if _ext is not None else None


class ScanFileNotFound(FileNotFoundError):
    """No readable sequence file for a haplotype path. Subclasses FileNotFoundError so the
    CLI's existing `except FileNotFoundError -> ClickException` handlers still catch it."""


def resolve_scan_path(fasta: str) -> str:
    """Map a haplotype path from `samples.tsv` to the file this backend actually reads.

    `config/samples.tsv` still points `local_path` at the uncompressed `.fa` (Phase 3
    repoints it at the `.fa.gz`), so both must work today. Precedence, in order:

    1. the path itself, if it already ends in `.gz` (it is expected to be BGZF);
    2. `<path>.gz` if it exists -- the HPRC-distributed BGZF file sitting beside the `.fa`.
       Preferred over the `.fa` even when both exist: 0.90 GB of I/O instead of 3.08 GB;
    3. the path itself, if it exists (a plain uncompressed FASTA; supported, just slower).

    Raises `ScanFileNotFound` naming both candidates if neither is present. This is
    deliberately a plain, testable function rather than a fallback buried in the caller --
    see `tests/test_rust_backend_differential.py::TestResolveScanPath`.
    """
    p = Path(fasta)
    if p.suffix == ".gz":
        if not p.exists():
            raise ScanFileNotFound(f"missing sequence file {p}")
        return str(p)
    gz = Path(str(p) + ".gz")
    if gz.exists():
        return str(gz)
    if p.exists():
        return str(p)
    raise ScanFileNotFound(
        f"no sequence file for {p}: neither {gz} (preferred; BGZF, as HPRC distributes it) "
        f"nor {p} exists. Fetch it with scripts/prepare_haplotypes.sh."
    )


def _to_sites(
    primer_seq: str, haplotype_id: str, tuples: list[tuple]
) -> list[BindingSite]:
    return [
        BindingSite(
            primer_name=primer_seq,
            haplotype_id=haplotype_id,
            chrom=chrom,
            start=start,
            end=end,
            strand=Strand(strand),
            mismatches=mm,
            mismatch_offsets_3p=list(offsets),
        )
        for (chrom, start, end, strand, mm, offsets) in tuples
    ]


def find_binding_sites_batch(
    seqs: list[str],
    fasta: str,
    haplotype_id: str,
    max_mismatches: int,
    *,
    slop: int = 3,
    fa=None,
    truncated: dict[str, int] | None = None,
) -> dict[str, list[BindingSite]]:
    """Genome-wide binding sites for many primers in one pass over the assembly.

    Signature-compatible with `bwa_backend.find_binding_sites_batch`. `truncated` is accepted
    and never populated: this backend enumerates every site exhaustively, so its result is
    never incomplete. Only bwa can truncate (it drops the XA tag above `samse -n`).

    `fa` (an already-open `pysam.FastaFile` the caller owns) is accepted and ignored: this
    backend reads the BGZF
    file itself and never needs random access.
    """
    if _ext is None:  # pragma: no cover - guarded by search.resolve_backend
        raise RuntimeError(
            "the compiled pangenome_primer._scan extension is not available; "
            "build it with `maturin develop --release` or select search.backend: bwa"
        )
    del fa
    path = resolve_scan_path(fasta)
    raw = _ext.scan(list(seqs), path, int(max_mismatches), int(slop))
    return {seq: _to_sites(seq, haplotype_id, tups) for seq, tups in raw.items()}


def find_binding_sites_in_seq(
    primer_name: str,
    primer_seq: str,
    ref: str,
    haplotype_id: str,
    chrom: str,
    max_mismatches: int,
) -> list[BindingSite]:
    """In-memory single-primer search, matching `binding.find_binding_sites_naive`'s
    signature. This is what `binding.find_binding_sites(backend="rust")` dispatches to."""
    if _ext is None:  # pragma: no cover
        raise RuntimeError("pangenome_primer._scan is not available")
    raw = _ext.scan_seq([primer_seq], ref, int(max_mismatches), chrom)
    sites = _to_sites(primer_seq, haplotype_id, raw.get(primer_seq, []))
    for s in sites:
        s.primer_name = primer_name
    return sites
