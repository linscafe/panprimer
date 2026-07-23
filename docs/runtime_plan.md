# Runtime reduction plan (per-query)

## Context

The first real run (GAPDH, 10 haplotypes) took **~83 min**, dominated by *locus projection*
rebuilding a full-genome minimap2 index per haplotype **on every query**. The one-time costs
(`bwa index`, ~1 hr/haplotype) are already amortized into `scripts/prepare_haplotypes.sh` and
cached on disk. This plan removes the *repeated per-query* cost so each new target locus is
fast.

Measured per-query costs (before this plan):
1. Locus projection — ~40 min (10 × full-genome minimap2 index build). **Biggest cost.**
2. Genome-wide `bwa` off-target search — one BWT load per haplotype (batched already).
3. Template extraction, mask, Primer3, ranking — negligible.

## Item 1 — Whole-genome alignment cache (projection → lookup)

Precompute the **CHM13 ↔ haplotype alignment once per haplotype** (in prep), cached as PAF.
Each query then projects the target locus by a **PAF coordinate lift** (find alignment blocks
overlapping the CHM13 target ± flank, map to haplotype coordinates) instead of aligning at
query time. No per-query minimap2, no multi-GB index load.

- `align_cache.py`:
  - `build_alignment(chm13, hap_fasta, out_paf)` → `minimap2 -cx asm5 chm13 hap` (target=CHM13).
  - `project_from_paf(paf, chrom, start, end, hap_fasta, flank)` → `Projection` (homologous
    window + haplotype sequence for the mask), handling strand.
- `project.project_locus` prefers the PAF cache (`<hap>.chm13.paf`) and falls back to the
  on-the-fly `.mmi` path when no cache exists.
- Prep builds the PAF cache per haplotype (replaces the `.mmi` step as the projection index).

## Item 2 — Two-stage search (cheap on-target, then expensive off-target)

Split evaluation so the expensive genome-wide search only runs on a shortlist:

- **Stage A (cheap, no BWT):** evaluate each candidate against the *projected homologous
  window* only (naive search on the ~40 kb window we already extracted for the mask), with the
  thermodynamic dropout model. Yields **on-target coverage / dropout** per candidate. This is
  the anti-dropout ranking, computed in seconds with no genome-wide index load.
- Rank by Stage A coverage (+ Primer3 penalty); keep the top-K (`--top-k`, default 5).
- **Stage B (expensive, BWT):** run the genome-wide `bwa` search only for the top-K to detect
  `off_target` / `multi_product`, finalize per-haplotype status, `unique_product_rate`, and the
  final ranking.

Non-shortlisted candidates are reported with Stage-A metrics (coverage known, specificity "not
assessed"), so nothing is hidden.

## Items 4 & 5 — RAM fit (MVP)

MVP uses **4 haplotypes** (one per AFR/EUR/EAS/AMR) so peak memory fits comfortably in 15 GB
(projection is now a PAF lookup; Stage B loads one BWT at a time). This is a capability demo,
not the full subset. Persistent-worker / heavier parallelism stay a cloud-profile concern.

## Demo layout

```
demo/
├── samples.tsv        # 4 diverse haplotypes (subset of config/samples.tsv)
├── run_demo.sh        # builds PAF caches (once) then runs the two-stage pipeline
└── results/<locus>/   # results.tsv / results.json / report.html
```

## Verification

- Unit: PAF lift maps a known CHM13 interval to the expected haplotype window; two-stage
  produces the same top pair as the full search on the synthetic fixture.
- End-to-end: `demo/run_demo.sh` on GAPDH (`chr12:6544868-6548730`) — expect the same best
  universal primer as the 10-haplotype run, in a fraction of the time; then CYP2D6 hard case.
```
