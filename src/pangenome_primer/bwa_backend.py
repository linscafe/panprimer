"""Genome-wide primer binding search on a full haplotype assembly.

`bwa aln`/`samse` locates candidate binding positions across a multi-Gb assembly quickly;
we then extract each candidate window and recompute the exact mismatches and 3'-end offsets
with the naive comparator from `binding.py` (which the tests verify). bwa prunes the search
space; the trusted code decides the details. Requires `bwa` + `samtools` + `pysam`.
"""
from __future__ import annotations

import subprocess
import tempfile
import warnings
from pathlib import Path

from .binding import find_binding_sites_naive
from .model import BindingSite


def ensure_index(fasta: str) -> None:
    """Verify the bwa index for `fasta` exists. Building one in-process would be a silent
    ~55 min surprise on a multi-Gb assembly, so this fails loudly instead — indexes are meant
    to be built ahead of time by scripts/prepare_haplotypes.sh. `samtools faidx` IS still
    auto-built here: it takes seconds, not tens of minutes."""
    if not Path(fasta + ".bwt").exists():
        raise FileNotFoundError(
            f"missing bwa index for {fasta} (no {fasta}.bwt).\n"
            "Build it with: bash scripts/prepare_haplotypes.sh "
            "(bwa index takes ~55 min on a ~3 Gb assembly; not built automatically)."
        )
    if not Path(fasta + ".fai").exists():
        subprocess.run(["samtools", "faidx", fasta], check=True)


#: `bwa samse -n` cap. BWA omits the XA tag ENTIRELY for a read with more hits than this --
#: it does not emit a truncated list -- so above the cap only the primary hit is recoverable.
XA_CAP = 1000


def _candidate_positions_batch(
    seqs: list[str], fasta: str, max_mm: int
) -> tuple[dict[int, set[tuple[str, int]]], dict[int, int]]:
    """Run ONE bwa aln/samse over all primer sequences (index loaded once).

    Returns `(candidates, truncated)`:
      * `candidates` = {query_index -> {(chrom, 0-based pos)}} from primary + XA hits.
      * `truncated`  = {query_index -> true total hit count} for queries whose hit list bwa
        did NOT fully report.

    The second value exists because of a silent-false-negative bug. When a read has more
    than `-n` hits, `bwa samse` drops the XA tag rather than truncating it, so parsing
    primary+XA yields exactly ONE position with no indication that thousands were discarded.
    Downstream that looks like "this primer binds in one place", and if no partner primer
    sits nearby the pair is reported as a DROPOUT -- a confident wrong answer.

    That is not hypothetical: the demo's original CYP2D6_dropout forward primer is an Alu
    consensus reported as `X0=12760 X1=317208` with **no XA tag**. It was recorded as a
    dropout on every haplotype for years' worth of runs; it actually binds ~330k places.

    The true count is still recoverable, because BWA reports it in the X0 (best hits) and
    X1 (suboptimal hits) tags even when it suppresses XA. X0+X1 matches the exhaustive
    scanner exactly (12760+317208 = 329,968; 3+113 = 116), so it is a reliable signal.
    """
    with tempfile.TemporaryDirectory() as td:
        q = Path(td) / "q.fa"
        # bwa seeding needs -l/-k sized to the shortest primer
        minlen = min(len(s) for s in seqs)
        q.write_text("".join(f">{i}\n{s}\n" for i, s in enumerate(seqs)))
        sai = Path(td) / "q.sai"
        with open(sai, "wb") as fh:
            subprocess.run(
                ["bwa", "aln", "-n", str(max_mm), "-o", "0", "-l", str(minlen),
                 "-k", str(max_mm), "-N", fasta, str(q)],
                check=True, stdout=fh, stderr=subprocess.DEVNULL,
            )
        sam = subprocess.run(
            ["bwa", "samse", "-n", str(XA_CAP), fasta, str(sai), str(q)],
            check=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True,
        ).stdout
    cands: dict[int, set[tuple[str, int]]] = {i: set() for i in range(len(seqs))}
    truncated: dict[int, int] = {}
    for line in sam.splitlines():
        if line.startswith("@"):
            continue
        f = line.split("\t")
        if int(f[1]) & 4:  # unmapped
            continue
        qi = int(f[0])
        cands[qi].add((f[2], int(f[3]) - 1))
        x0 = x1 = 0
        for tag in f[11:]:
            if tag.startswith("XA:Z:"):
                for alt in tag[5:].split(";"):
                    if alt:
                        a = alt.split(",")
                        cands[qi].add((a[0], abs(int(a[1])) - 1))
            elif tag.startswith("X0:i:"):
                x0 = int(tag[5:])
            elif tag.startswith("X1:i:"):
                x1 = int(tag[5:])
        # X0+X1 is bwa's own count of hits. If it exceeds what we could actually parse, the
        # hit list is incomplete and any verdict built on it is unsafe -- record the real
        # total so callers can say so instead of silently under-reporting.
        total = x0 + x1
        if total > len(cands[qi]):
            truncated[qi] = total
    return cands, truncated


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
    """Genome-wide binding sites for many primers in one index load. Deduplicates identical
    sequences (candidate pairs share primers). Returns {sequence -> [BindingSite]}. Each
    candidate position's exact mismatches/3' offsets are recomputed with the naive
    comparator, so bwa only prunes the search space. Pass an already-open `fa` (pysam.
    FastaFile on `fasta`) to reuse a handle a caller holds open across other calls on the
    same haplotype; the caller retains ownership and closes it. Left as None (default), a
    handle is opened and closed here, as before.

    Pass a dict as `truncated` to receive {primer_sequence -> true total hit count} for any
    primer whose hit list bwa could not fully report (see `_candidate_positions_batch`). For
    those primers the returned site list is INCOMPLETE and must not be read as "these are all
    the binding sites" -- doing so is what made a ~330k-site Alu primer look like a dropout.
    A warning is emitted regardless of whether the dict is supplied, so the condition is
    never silent even for callers that do not opt in."""
    ensure_index(fasta)
    import pysam

    uniq = sorted(set(seqs))
    cands, trunc_idx = _candidate_positions_batch(uniq, fasta, max_mismatches)
    for qi, total in sorted(trunc_idx.items()):
        seq = uniq[qi]
        if truncated is not None:
            truncated[seq] = total
        warnings.warn(
            f"bwa reported {total} hits for primer {seq} on {haplotype_id} but omitted the "
            f"XA tag (samse -n {XA_CAP} cap), so only {len(cands[qi])} position(s) are "
            f"recoverable. This primer is likely repeat-derived; its binding-site list is "
            f"INCOMPLETE and any pass/dropout verdict from it would be unreliable. Use "
            f"search.backend=rust for an exhaustive scan.",
            RuntimeWarning,
            stacklevel=2,
        )
    owns_fa = fa is None
    if owns_fa:
        fa = pysam.FastaFile(fasta)
    out: dict[str, list[BindingSite]] = {}
    for i, primer_seq in enumerate(uniq):
        L = len(primer_seq)
        seen: set[tuple[str, int, str]] = set()
        sites: list[BindingSite] = []
        for chrom, pos in cands[i]:
            clen = fa.get_reference_length(chrom)
            ws, we = max(0, pos - slop), min(clen, pos + L + slop)
            window = fa.fetch(chrom, ws, we)
            for s in find_binding_sites_naive(
                primer_seq, primer_seq, window, haplotype_id, chrom, max_mismatches
            ):
                gstart, gend = ws + s.start, ws + s.end
                key = (chrom, gstart, s.strand.value)
                if key in seen:
                    continue
                seen.add(key)
                s.start, s.end = gstart, gend
                sites.append(s)
        out[primer_seq] = sites
    if owns_fa:
        fa.close()
    return out
