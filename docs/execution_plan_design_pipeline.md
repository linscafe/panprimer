# Runtime reduction plan (per-query)

## Context

The first real run (GAPDH, 10 haplotypes) took **~83 min**, dominated by *locus projection* rebuilding a full-genome minimap2 index per haplotype **on every query**. The one-time costs (`bwa index`, ~1 hr/haplotype) are already amortized into `scripts/prepare_haplotypes.sh` and cached on disk. This plan removes the *repeated per-query* cost so each new target locus is fast.

Measured per-query costs (before this plan):
1. Locus projection — ~40 min (10 × full-genome minimap2 index build). **Biggest cost.**
2. Genome-wide `bwa` off-target search — one BWT load per haplotype (batched already).
3. Template extraction, mask, Primer3, ranking — negligible.

## Item 1 — Whole-genome alignment cache (projection → lookup)

Precompute the **CHM13 ↔ haplotype alignment once per haplotype** (in prep), cached as PAF. Each query then projects the target locus by a **PAF coordinate lift** (find alignment blocks overlapping the CHM13 target ± flank, map to haplotype coordinates) instead of aligning at query time. No per-query minimap2, no multi-GB index load.

- `align_cache.py`:
  - `build_alignment(chm13, hap_fasta, out_paf)` → `minimap2 -cx asm5 chm13 hap` (target=CHM13).
  - `project_from_paf(paf, chrom, start, end, hap_fasta)` → `Projection` (homologous window + haplotype sequence for the mask), handling strand.
- `project.project_target` prefers the PAF cache (`<hap>.chm13.paf`) and falls back to the `.mmi`-cached on-the-fly path when no PAF exists.

### Measured tradeoff (important)

On this 15 GB box the whole-genome PAF is **too expensive as a one-time cost**: building it took **114 min** for one haplotype and **OOM-killed** two others (a 3 Gb-vs-3 Gb `minimap2` alignment with 16 threads exceeds 15 GB). The instant per-query lift is real, but the prep is worse than the problem here.

**Resolution:** the memory-safe `.mmi` projection index (build ~12 min, loads in seconds per query) is the **local default**; the whole-genome **PAF is opt-in** (`BUILD_PAF=1`, capped threads) and only worth it on high-RAM / cloud where the one-time alignment is affordable. Both code paths ship; `project_target` uses whichever cache is present.

## Item 2 — Two-stage search (cheap on-target, then expensive off-target)

Split evaluation so the expensive genome-wide search only runs on a shortlist:

- **Stage A (cheap, no BWT):** evaluate each candidate against the *projected homologous window* only (naive search on the ~40 kb window we already extracted for the mask), with the thermodynamic dropout model. Yields **on-target coverage / dropout** per candidate. This is the anti-dropout ranking, computed in seconds with no genome-wide index load.
- Rank by Stage A coverage (+ Primer3 penalty); keep the top-K (`--top-k`, default 5).
- **Stage B (expensive, BWT):** run the genome-wide `bwa` search only for the top-K to detect `off_target` / `multi_product`, finalize per-haplotype status, `unique_product_rate`, and the final ranking.

Non-shortlisted candidates are reported with Stage-A metrics (coverage known, specificity "not assessed"), so nothing is hidden.

## Items 4 & 5 — RAM fit (MVP): avoiding OOM

*(Item 3 — the `.mmi` projection cache — is folded into Item 1's resolution above. Items 4 & 5 from the original review, persistent worker / parallelism, are RAM-bound and reframed here as the fit-in-15 GB strategy.)*

MVP uses **3 haplotypes** (AFR/EUR/EAS: HG01884, HG00097, HG00408), all with prebuilt caches, so a demo run does **no index building** and touches one genome-scale index at a time. This is a capability demo, not the full subset.

**Why OOM happened and the plan to avoid it.** The kills came from *building* the whole-genome PAF: a 3 Gb-vs-3 Gb `minimap2` alignment with 16 threads holds both genomes' data plus per-thread alignment buffers, exceeding 15 GB. The rules that keep the pipeline inside 15 GB:

1. **One genome-scale index in RAM at a time.** Projection processes haplotypes strictly sequentially (load index → project → release → next); Stage B loads one BWT per haplotype in turn. **Peak ≈ a single ~6 GB index — independent of haplotype count**, so adding haplotypes costs time, never peak memory.
2. **Prefer the light cache.** `.mmi` *loading* (~6 GB, one at a time) instead of the whole-genome PAF *alignment* (both genomes + thread buffers — the thing that OOM'd).
3. **Cap threads on any alignment.** Memory scales with thread count; the opt-in PAF build uses `-t 4`, and heavy steps should never run at full 16-thread width on this box.
4. **Two-stage isolates the heavy step.** Stage A works only on the small projected windows (KB in RAM, no genome index); only Stage B loads a genome index, one haplotype at a time, and only for the top-K shortlist.
5. **Scale horizontally, not by RAM.** More haplotypes → more sequential time locally, or one haplotype per worker on the cloud profile — never a larger local working set.

## Demo layout

```
demo-design-pipeline/
├── samples.tsv        # 3 diverse haplotypes (AFR/EUR/EAS; subset of config/samples.tsv)
├── run_demo.sh        # ensures projection caches, then runs the two-stage pipeline
└── <locus>/           # results.tsv / results.json / report.md / report.html
```

`report.md` is a readable/diffable Markdown intermediate; the HTML is rendered from the built-in template by default, or from the Markdown via Quarto with `--quarto` (needs `quarto` on PATH — an optional dependency, so nothing regresses without it). Both HTML reports use a fixed light theme (white background, dark text).

## Verification

- Unit: PAF lift maps a known CHM13 interval to the expected haplotype window; two-stage produces the same top pair as the full search on the synthetic fixture.
- End-to-end: `demo-design-pipeline/run_demo.sh` on GAPDH (`chr12:6544868-6548730`) — expect the same best universal primer as the 10-haplotype run, in a fraction of the time; then CYP2D6 hard case.
