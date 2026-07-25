"""The genome-wide search seam: which backend answers "where does this primer bind?".

`config/defaults.yaml`'s `search.backend` used to be marked RESERVED with no reader. This
module is that reader, and the single place `verify.py` / `cli.py` go through instead of
importing a concrete backend.

Backends
--------
``rust``   the compiled `pangenome_primer._scan` extension (`rust_backend.py`). Streams the
           BGZF `.fa.gz`, no index on disk, exhaustive. The default.
``bwa``    `bwa aln`/`samse` over the per-haplotype bwa index (`bwa_backend.py`). Needs the
           uncompressed `.fa` plus a ~5.3 GB index. Kept working, and kept selectable, as
           the reference the rust path was validated against.
``naive``  pure Python (`binding.find_binding_sites_naive`), for tiny references only. It is
           the correctness reference, not a genome-scale option: an O(N*L) Python scan of a
           3.1 Gb assembly does not finish in useful time, so `find_binding_sites_batch`
           refuses it above a small size rather than appearing to hang.

Degrade, never break
--------------------
`resolve_backend` never fails just because the extension is missing. ``rust`` (and
``auto``) fall through to ``bwa`` with a warning; only an explicitly-named backend that is
unavailable, or an unknown name, is an error. A platform without a wheel and without a Rust
toolchain keeps working exactly as it did before Phase 2.
"""
from __future__ import annotations

import warnings

from .model import BindingSite

#: Fallback chain per requested backend. `None` means "no fallback: use it or fail".
_CHAIN: dict[str, tuple[str, ...]] = {
    "auto": ("rust", "bwa"),
    "rust": ("rust", "bwa"),
    "bwa": ("bwa",),
    "naive": ("naive",),
}

#: Above this many reference bases the naive backend is refused rather than run. A 3.1 Gb
#: assembly at naive's pure-Python O(N*L) rate is a multi-hour "hang"; the mini-genome
#: fixture (450 kb) and projected windows stay well under it.
NAIVE_MAX_BASES = 50_000_000


def _rust_ok() -> bool:
    from . import rust_backend

    return rust_backend.available()


def _bwa_ok() -> bool:
    import shutil

    return shutil.which("bwa") is not None


_PROBES = {"rust": _rust_ok, "bwa": _bwa_ok, "naive": lambda: True}


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
    for i, candidate in enumerate(chain):
        if _PROBES[candidate]():
            if i:
                why = (
                    "the compiled pangenome_primer._scan extension is not installed"
                    if chain[0] == "rust"
                    else f"{chain[0]!r} is not available"
                )
                warn(
                    f"search.backend={requested!r}: {why}; falling back to {candidate!r}."
                )
            return candidate
    raise RuntimeError(
        f"no usable genome-wide search backend for search.backend={requested!r}: tried "
        f"{list(chain)}. Install the compiled extension (`maturin develop --release`) or "
        f"put `bwa` on PATH."
    )


def backend_from_config(raw: dict | None, *, warn=warnings.warn) -> str:
    """Read `search.backend` out of a loaded `config/defaults.yaml` dict and resolve it."""
    name = None
    if raw:
        name = (raw.get("search") or {}).get("backend")
    return resolve_backend(name, warn=warn)


def find_binding_sites_batch(
    seqs: list[str],
    fasta: str,
    haplotype_id: str,
    max_mismatches: int,
    *,
    slop: int = 3,
    fa=None,
    backend: str | None = None,
    truncated: dict[str, int] | None = None,
) -> dict[str, list[BindingSite]]:
    """Genome-wide binding sites for many primers against one haplotype.

    Identical contract to `bwa_backend.find_binding_sites_batch`, which is what this
    replaces at `verify.py`, `cli.py:run` and `cli.py:evaluate`: returns
    `{primer_sequence: [BindingSite, ...]}` keyed by SEQUENCE (not primer name), with
    `primer_name` set to the sequence, deduplicated on `(chrom, start, strand)`.

    `truncated`, if given, collects {primer_sequence -> true hit count} for primers whose
    site list is INCOMPLETE. Only the bwa backend can populate it: above `samse -n` bwa drops
    the XA tag wholesale, leaving one recoverable position for a primer that may bind
    hundreds of thousands of places. Callers must not treat those lists as exhaustive -- that
    is precisely how a ~330k-site Alu primer was reported as a dropout. `rust` and `naive`
    enumerate everything and never truncate.
    """
    chosen = resolve_backend(backend)
    if chosen == "rust":
        from .rust_backend import find_binding_sites_batch as run

        return run(seqs, fasta, haplotype_id, max_mismatches,
                   slop=slop, fa=fa, truncated=truncated)
    if chosen == "bwa":
        from .bwa_backend import find_binding_sites_batch as run

        return run(seqs, fasta, haplotype_id, max_mismatches,
                   slop=slop, fa=fa, truncated=truncated)
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
                f"reference for small references, not a genome-scale search; use 'rust' "
                f"or 'bwa'."
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
