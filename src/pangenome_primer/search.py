"""The genome-wide search seam: which backend answers "where does this primer bind?".

`config/defaults.yaml`'s `search.backend` used to be marked RESERVED with no reader. This
module is that reader, and the single place `verify.py` / `cli.py` go through instead of
importing a concrete backend.

Backends
--------
``rust``   the compiled `pangenome_primer._scan` extension (`rust_backend.py`). Streams the
           BGZF `.fa.gz`, no index on disk, exhaustive. The default, and the only
           genome-scale option.
``naive``  pure Python (`binding.find_binding_sites_naive`), for tiny references only. It is
           the correctness reference, not a genome-scale option: an O(N*L) Python scan of a
           3.1 Gb assembly does not finish in useful time, so `find_binding_sites_batch`
           refuses it above a small size rather than appearing to hang.

Requires the compiled extension
-------------------------------
There is no genome-scale fallback. The `bwa` backend used to serve that role, and removing
it removed the "degrade, never break" path with it: a platform with no wheel and no Rust
toolchain now gets a clear error naming `maturin develop --release`, rather than a slower
but working search. That is deliberate -- `bwa` could not be trusted on repeat-derived
primers (it silently dropped the XA tag above `samse -n`, reporting a 330k-site Alu primer
as a single unique hit), so keeping it as a fallback meant keeping a path that could answer
*wrongly*. Failing loudly is better than degrading to that.
"""
from __future__ import annotations

import warnings

from .model import BindingSite

#: Fallback chain per requested backend. `None` means "no fallback: use it or fail".
_CHAIN: dict[str, tuple[str, ...]] = {
    "auto": ("rust",),
    "rust": ("rust",),
    "naive": ("naive",),
}

#: Above this many reference bases the naive backend is refused rather than run. A 3.1 Gb
#: assembly at naive's pure-Python O(N*L) rate is a multi-hour "hang"; the mini-genome
#: fixture (450 kb) and projected windows stay well under it.
NAIVE_MAX_BASES = 50_000_000


def _rust_ok() -> bool:
    from . import rust_backend

    return rust_backend.available()


_PROBES = {"rust": _rust_ok, "naive": lambda: True}


def resolve_backend(name: str | None = None, *, warn=warnings.warn) -> str:
    """Turn a configured `search.backend` value into a backend that is actually usable.

    `warn` is injected so the CLI can route the message through `click.echo`; tests use it
    to assert the fallback happened rather than parsing stderr.
    """
    requested = (name or "auto").strip().lower()
    chain = _CHAIN.get(requested)
    if chain is None:
        raise ValueError(
            f"unknown search.backend {name!r}; expected one of "
            f"{sorted(_CHAIN)} (see config/defaults.yaml)"
        )
    for candidate in chain:
        if _PROBES[candidate]():
            return candidate
    raise RuntimeError(
        f"no usable genome-wide search backend for search.backend={requested!r}: tried "
        f"{list(chain)}. Build the compiled extension with `maturin develop --release`."
    )


def backend_from_config(raw: dict | None, *, warn=warnings.warn) -> str:
    """Read `search.backend` out of a loaded `config/defaults.yaml` dict and resolve it."""
    name = None
    if raw:
        name = (raw.get("search") or {}).get("backend")
    return resolve_backend(name, warn=warn)


def threads_from_config(raw: dict | None) -> int | None:
    """Read `search.threads` out of a loaded config dict.

    None (absent, null, or 0) means "let the scanner choose": all cores up to its own cap.
    That is the right default for the CLI as it stands, which loops haplotypes one at a time
    -- a solo scan should take the machine. It is the wrong default the moment haplotypes run
    concurrently, which is why the key exists. See ISSUE-002 in `docs/scanner_notes.md` for the
    measured curves; the short version is `threads x concurrent-haplotypes ~= core count`.
    """
    if not raw:
        return None
    n = (raw.get("search") or {}).get("threads")
    if n is None:
        return None
    n = int(n)
    return n if n > 0 else None


def find_binding_sites_batch(
    seqs: list[str],
    fasta: str,
    haplotype_id: str,
    max_mismatches: int,
    *,
    slop: int = 3,
    fa=None,
    backend: str | None = None,
    threads: int | None = None,
) -> dict[str, list[BindingSite]]:
    """Genome-wide binding sites for many primers against one haplotype.

    Returns
    `{primer_sequence: [BindingSite, ...]}` keyed by SEQUENCE (not primer name), with
    `primer_name` set to the sequence, deduplicated on `(chrom, start, strand)`.

    Every backend here is **exhaustive**: it proves <= `max_mismatches` or > `max_mismatches`
    at every position, so the returned list is never a sample or a truncation. Callers may
    rely on that. (The removed `bwa` backend could not promise it -- above `samse -n` it
    dropped the XA tag wholesale and reported a ~330k-site Alu primer as a single unique
    hit, which the pipeline then called an allele dropout.)

    `threads` sizes the compiled scanner's worker pool; None lets it choose. It is ignored by
    the naive backend, which is single-threaded.
    """
    chosen = resolve_backend(backend)
    if chosen == "rust":
        from .rust_backend import find_binding_sites_batch as run

        return run(seqs, fasta, haplotype_id, max_mismatches, slop=slop, fa=fa,
                   threads=threads)
    return _naive_batch(seqs, fasta, haplotype_id, max_mismatches, fa=fa)


def _naive_batch(
    seqs: list[str],
    fasta: str,
    haplotype_id: str,
    max_mismatches: int,
    *,
    fa=None,
) -> dict[str, list[BindingSite]]:
    """Whole-FASTA scan with the pure-Python comparator. Correct, and far too slow for a
    real assembly -- hence the explicit size refusal instead of a silent multi-hour run."""
    import pysam

    from .binding import find_binding_sites_naive

    owns = fa is None
    if owns:
        fa = pysam.FastaFile(fasta)
    try:
        total = sum(fa.get_reference_length(c) for c in fa.references)
        if total > NAIVE_MAX_BASES:
            raise ValueError(
                f"search.backend='naive' refused on {fasta}: {total:,} bases exceeds the "
                f"{NAIVE_MAX_BASES:,}-base guard. The naive backend is the correctness "
                f"reference for small references, not a genome-scale search; use 'rust'."
            )
        uniq = sorted(set(seqs))
        out: dict[str, list[BindingSite]] = {s: [] for s in uniq}
        for chrom in fa.references:
            ref = fa.fetch(chrom)
            for s in uniq:
                out[s].extend(
                    find_binding_sites_naive(s, s, ref, haplotype_id, chrom, max_mismatches)
                )
        return out
    finally:
        if owns:
            fa.close()
