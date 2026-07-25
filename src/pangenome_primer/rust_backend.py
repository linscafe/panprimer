"""Genome-wide primer binding search backed by the compiled Rust scanner.

The genome-wide search engine. `search.py` is the dispatcher that selects it.

Two properties worth stating plainly:

* **No index.** It streams the BGZF `.fa.gz` HPRC already distributes, so no per-haplotype
  index exists on disk at all.
* **Exhaustive, not heuristic.** It proves <= `max_mismatches` or > `max_mismatches` at
  every position for every primer in both orientations. The exact comparator IS the search,
  so there is no candidate-generation or re-scoring pass, and no cap on reported hits.
* **`slop` is inert.** It exists in the signature only because the seam has it; there are no
  candidate windows to widen. See `rust/pgp-scan/src/lib.rs`.

`available()` is false when the extension could not be imported (no wheel for this platform,
no Rust toolchain at install time); `search.resolve_backend` then raises with a build hint.
Nothing here raises at import.
"""
from __future__ import annotations

import warnings
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


def _threads(n: int | None) -> int | None:
    """Normalise a `search.threads` value. 0, null and absent all mean "let the scanner
    choose" -- one knob, one meaning, rather than 0 being a distinct third state."""
    if n is None:
        return None
    n = int(n)
    return n if n > 0 else None


_pool_warned = False


def _warn_if_pool_already_sized(requested: int | None) -> None:
    """Warn once if the pool was already built with a different thread count.

    The pool is a process-wide `OnceLock`: the first scan fixes the worker count, and every
    later `threads=` is quietly ignored. That is fine for the pipelines here, which pass one
    configured value for a whole run -- but "quietly ignored" is how a `search.threads: 4`
    meant to keep a 64-way fan-out from thrashing becomes a no-op that nobody notices. Say so
    instead.
    """
    global _pool_warned
    if requested is None or _pool_warned or _ext is None:
        return
    live = int(_ext.pool_threads(requested))
    if live != requested:
        _pool_warned = True
        warnings.warn(
            f"search.threads={requested} was ignored: the scanner's thread pool was already "
            f"built with {live} workers earlier in this process and is not resized. The "
            f"first scan in a process wins.",
            RuntimeWarning,
            stacklevel=3,
        )


def pool_threads(requested: int | None = None) -> int | None:
    """Workers in the scanner's rayon pool -- the live count once a scan has run, otherwise
    the count that *would* be chosen.

    The pool is built on first use and **not rebuilt**, so in a process that scans several
    haplotypes the first call fixes the thread count for all of them. Every caller here
    passes the same configured value for a whole run, but this function exists so that
    assumption can be checked rather than trusted.

    Returns None when the extension is absent, so a caller can report "unknown" instead of a
    fabricated number.
    """
    if _ext is None:
        return None
    return int(_ext.pool_threads(_threads(requested)))


class ScanFileNotFound(FileNotFoundError):
    """No readable sequence file for a haplotype path. Subclasses FileNotFoundError so the
    CLI's existing `except FileNotFoundError -> ClickException` handlers still catch it."""


def resolve_scan_path(fasta: str) -> str:
    """Map a haplotype path from `samples.tsv` to the file this backend actually reads.

`samples.tsv` normally points `local_path` at the `.fa.gz`, but an uncompressed `.fa` is
    still supported. Precedence, in order:

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
    threads: int | None = None,
) -> dict[str, list[BindingSite]]:
    """Genome-wide binding sites for many primers in one pass over the assembly.

    Exhaustive: every position is proved to be within or outside the mismatch budget, so the
    result is never a truncation or a sample.

    `fa` (an already-open `pysam.FastaFile` the caller owns) is accepted and ignored: this
    backend reads the BGZF
    file itself and never needs random access.

    `threads` sizes the scanner's rayon pool; see `pool_threads` for what `None` means and
    for the caveat that the pool is built once per process.
    """
    if _ext is None:  # pragma: no cover - guarded by search.resolve_backend
        raise RuntimeError(
            "the compiled pangenome_primer._scan extension is not available; "
            "build it with `maturin develop --release`"
        )
    del fa
    path = resolve_scan_path(fasta)
    want = _threads(threads)
    _warn_if_pool_already_sized(want)
    raw = _ext.scan(list(seqs), path, int(max_mismatches), int(slop), want)
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
    raw = _ext.scan_seq([primer_seq], ref, int(max_mismatches), chrom, None)
    sites = _to_sites(primer_seq, haplotype_id, raw.get(primer_seq, []))
    for s in sites:
        s.primer_name = primer_name
    return sites
