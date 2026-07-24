# Verify mode plan — screen user-supplied primer pairs

## Context

This is **Decision 1 option D** from the original design grilling ("screen existing primers: user supplies primers; tool only evaluates specificity/dropout across the pangenome; no candidate generation") — consciously deferred at the time, now requested.

The user brings primer pairs and asks: across the HPRC haplotypes, does each pair (1) create off-target products, and (2) what amplicon size does it give (which may differ by a few bp per haplotype due to indels)? Output is a matrix: one row per pair, one column per haplotype, cells showing amplicon sizes — green for correct (on-target), red for off-target.

## Input

CSV (`--primers`) with header `primer_id,target,forward,reverse`:
- `primer_id` — label
- `target` — **GRCh38 `chr:start-end`**, the *intended amplicon span* (decision)
- `forward`, `reverse` — primer sequences, 5′→3′

## Decisions (from review)

- **Off-target size bound**: F and R convergent within **≤ 2 kb** (default `--max-amplicon`, configurable), outside the target. Standard-PCR reach.
- **Coordinate meaning**: `target` = intended amplicon span; expected size = `end − start`. On-target = amplicon within the projected target window; on-target sizes deviating from the expected span by > `--size-tolerance` (default 20 bp) are flagged.
- **Dropout**: a haplotype with no correct amplicon (binding-site variant) is shown as a distinct grey `dropout` cell, even when off-targets also form.
- **Uncertain**: target not projectable in an assembly → `?`.

## Pipeline (reuses ~80% of the engine)

1. Parse CSV → `PrimerSpec` list.
2. **One genome-wide `bwa` search per haplotype** over *all* pairs' primers (batch `2N` sequences) — `find_binding_sites_batch`. Cost ~flat in the number of pairs.
3. Per pair: `resolve_target` (GRCh38→CHM13) → expected amplicon region + size.
4. Per pair × haplotype: `project_target` (cache-backed) → expected window; `evaluate_with_ sites` (competence + `pair_amplicons`) → amplicons with `size` + `on_target`.
5. Build the matrix and render.

Per-pair the pairing product-size cap is `max(--max-amplicon, expected_size + pad)`, so a correct amplicon larger than the off-target cap is never filtered out.

## New code (small)

- `samples.py` — shared `load_haplotypes` (also de-dupes the loader used by `run`).
- `verify.py` — `parse_primer_csv`, `run_verify` → `list[VerifyRow]` (rows of `VerifyCell`).
- `report.py` — `write_verify` (verify.json / verify.tsv / verify_matrix.md / verify_matrix.html) + `report/verify_matrix.html.j2` (green on-target / red off-target sizes, grey dropout, `?`).
- `cli.py` — `pangenome-primer verify` subcommand.
- `classify.pair_tm` — promote the tiny Tm helper for reuse.

Reused as-is: `resolve_target`, `project_target`, `find_binding_sites_batch`, `is_competent`, `pair_amplicons`, `evaluate_with_sites`.

## Output

- **verify_matrix.html** — the colored matrix (primary deliverable).
- **verify_matrix.md** — a readable/diffable Markdown intermediate (same matrix). Pass `--quarto` to render the HTML from it via Quarto (needs `quarto` on PATH); otherwise the built-in template renders the HTML directly. Colors carry as span classes (`[…]{.ok}`).
- **verify.json / verify.tsv** — scriptable source of truth (per-cell on/off sizes, status).

## Verification

- Unit: `parse_primer_csv`; `write_verify` renders green/red/dropout/`?` and a size-deviation flag from hand-built `VerifyRow`s (no external tools).
- Integration: `pangenome-primer verify` on the synthetic mini-pangenome (`--target-assembly chm13`) — a pair in a unique region shows a green on-target size; a pair whose region is duplicated shows a red off-target size; a variant that kills a primer shows `dropout`.
