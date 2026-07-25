# Scaling the verify pipeline: storage and runtime optimization

**Status:** proposal / analysis. Nothing here is implemented yet.
**Scope:** the [verify pipeline](intro_verify_pipeline.md) only. The design pipeline shares the same engine and would benefit similarly, but is not analysed here.

The demo runs against 3 haplotypes. HPRC Release 2 has ~464. This page measures what stands in the way of that jump and proposes a path. All figures marked **measured** were taken on this machine (16 CPU, 15 GB RAM, single SATA/NVMe volume) against `hprc-r2/assemblies/HG00097_hap1.fa`; the rest are estimates and are labelled as such.

---

## 1. The problem, stated once

Both the storage problem and the runtime problem have the same root cause:

> The pipeline treats every haplotype as an independent genome to be indexed, when 464 HPRC haplotypes are ~99.9% the same sequence.

Everything expensive follows from that. Storage scales linearly with haplotype count because each genome gets its own full-text index. Runtime scales linearly with haplotype count because each of those indexes has to be loaded from disk for a query that returns a handful of hits.

### 1.1 Storage — 15.1 GB per haplotype

Measured, for `HG00097_hap1`:

| File | Size | Required by |
|---|---:|---|
| `.fa` (uncompressed) | 3.08 GB | **only `bwa aln`** |
| `.fa.gz` (as downloaded from HPRC) | 0.90 GB | every sequence read in the pipeline |
| `.bwt` | 3.03 GB | `bwa aln` |
| `.sa` | 1.52 GB | `bwa aln` |
| `.pac` | 0.76 GB | `bwa aln` |
| `.mmi` | 5.80 GB | coordinate projection only |
| `.fai` / `.fa.gz.fai` | ~8 KB | random access |
| **Total** | **15.10 GB** | |

At ~464 haplotypes this is **~7.0 TB**.

The decisive finding is in the second row. **The `.fa.gz` that HPRC distributes is already BGZF-compressed**, so it supports random access directly. Verified:

```
magic 1f8b   FEXTRA flag True   BC subfield b'BC'   → BGZF
pysam.FastaFile('HG00097_hap1.fa.gz') → 75 refs
500 random 26 bp fetches: 0.0013 s
```

Every sequence read the verify pipeline performs is a small random slice — ~26 bp binding windows (`bwa_backend.py:89-91`), ~700 bp projection templates (`align_cache.py:135`), the CHM13 target template (`verify.py:125`). All of them work against the 0.90 GB BGZF file with no code change.

So the `gunzip -kf` in `scripts/prepare_haplotypes.sh:47` produces 3.08 GB that exists **solely** to feed `bwa aln`, and the 5.80 GB `.mmi` exists solely for coordinate projection. Between them, `.fa` + bwa index + `.mmi` account for **14.2 GB of the 15.1 GB — 94% of per-haplotype storage serves just those two functions.**

### 1.2 Runtime — 47 s per haplotype, nearly all of it index loading

Measured, 8 primer sequences against one haplotype using the exact flags from `bwa_backend.py:38-44`:

| Step | Wall time | Peak RSS |
|---|---:|---:|
| `bwa aln -n 3 -o 0 -l 20 -k 3 -N` | 12.2 s | 2.97 GB |
| `bwa samse -n 1000` | 35.1 s | 4.45 GB |
| **Total** | **47.3 s** | **4.45 GB** |

The search returned 8 hits. Forty-seven seconds of work to produce eight records is not a search cost — it is the cost of reading a 5.3 GB index off disk and into memory. This is the same point `intro_verify_pipeline.md` §5 already makes ("the expensive step is loading a whole genome's index into memory"), but the mitigation there — batch all primers into one pass per haplotype — only removes the *per-primer* factor. The *per-haplotype* factor remains, and it is the one that matters at 464.

Extrapolating: **~6.1 hours of serial wall time per verify request** at full pangenome scale, before any of the actual analysis.

> **Correction (measured 2026-07-24).** The 47.3 s above is the *search component only*. A full 3-haplotype demo verify run takes **17 min 05 s — ~5.7 min per haplotype end to end**, because coordinate projection loads a 5.80 GB `.mmi` per haplotype on top of the search. So the full-pangenome figure is **~44 h**, not ~6.1 h, and **projection is currently the larger cost, not search.** This does not change the diagnosis or the proposed steps — it raises the value of Step 3 (anchor grid, which eliminates the `.mmi`) relative to Step 2.

### 1.3 A live bug found while measuring

`hprc-r2/prepare.log` shows the whole-genome projection PAF (`prepare_haplotypes.sh:63`) took **114 minutes** for the one haplotype that completed, at **14.5 GB peak RSS**, and two others logged `PAF FAILED`. On a 15 GB machine those failures are almost certainly OOM. Passing `-K 100M` to cap minimap2's query batch size should fix it. This matters independently of everything below, because the PAF route is currently the only alternative to the 5.8 GB `.mmi`.

---

## 2. Does `bwa-mem4` help?

**No — it makes problem (1) worse and does not address problem (2).**

[`bwa-mem4`](https://lib.rs/crates/bwa-mem4) (v4.1.1, July 2026) is a Rust reimplementation of **bwa-mem2**, aiming for byte-identical SAM output. Three reasons it does not fit:

1. **Index size.** It requires a prebuilt index in bwa-mem2's format, which is *larger* than bwa's 5.3 GB. Storage gets worse.
2. **Wrong algorithm for 20-mers.** It is BWA-MEM, validated on 150 bp paired-end reads. MEM-based seeding is unreliable for 20 bp primers — which is precisely why `bwa_backend.py:38` uses `bwa aln` and not `bwa mem`.
3. **Speeds up the wrong thing.** Its claimed 2.62× is on the alignment kernel. Our 47 s is I/O-bound index loading, which a faster kernel does not touch.

If an in-process Rust aligner is wanted later for other reasons, [`bwa-mem3-rs`](https://lib.rs/crates/bwa-mem3-rs) or [10X Genomics' `rust-bwa`](https://github.com/10XGenomics/rust-bwa) (a wrapper over the actual BWA C API, including `aln`) are the relevant options — avoiding subprocess spawn and re-reading the index per call. Neither reduces index size.

---

## 3. Proposal

Four steps. Steps 1–3 are coupled and should land together; step 4 is the scale-out speed layer built on top.

### Step 1 — Keep only the BGZF file

**Storage: 15.10 GB → 0.90 GB per haplotype (~16.7×). At 464 haplotypes, ~7.0 TB → ~430 GB.**

- Stop `gunzip`-ing. Point `local_path` in `config/samples.tsv` at the `.fa.gz`, and drop the `.fa.gz` → `.fa` rewrite at `prepare_haplotypes.sh:81`.
- Generate the `.gzi` companion once (`samtools faidx` on the `.gz`, ~1 MB).
- Stop building the bwa index and the `.mmi` entirely.
- Add one shared CHM13 bundle: `.fa` + `.fai` + bwa index + `.mmi` ≈ 14.4 GB, amortized across all haplotypes. **Note the CHM13 bwa index does not currently exist** (`hprc-r2/references/` has only `.fa`, `.fai`, and `.asm5.mmi`) and must be built once.

Prep time per haplotype drops from ~60 min (bwa index alone is ~55 min, measured from `prepare.log` timestamps) to a download plus ~1 min of indexing.

`samples.py:18` hard-codes `Path(path).exists()`, which still works for a `.gz` path — but `hprc.py:249` already downloads the published `.fa.gz.fai` and nothing reads it. That becomes live.

### Step 2 — Replace the per-haplotype `bwa` search with a streaming scan

This is what makes step 1 possible: it removes the only consumer of the `.fa` and the bwa index.

Decompressing 0.90 GB of BGZF is ~10–20 s multithreaded. Finding ≤3-mismatch 20-mers in that stream is a 2-bit encode plus a rolling-hash pass, using pigeonhole seeding on exact 10-mers (~2,861 expected genomic hits per 10-mer seed in 3.1 Gb; for a typical panel this yields ~200k candidate positions). Those candidate windows then go to the **existing** `binding.find_binding_sites_naive`.

**Estimated ~30–45 s per haplotype** — roughly today's 47 s — with **zero index on disk** and ~200 MB RAM instead of 4.45 GB.

The integration point is clean. `bwa_backend.find_binding_sites_batch(seqs, fasta, haplotype_id, max_mismatches) -> dict[str, list[BindingSite]]` is imported lazily at exactly three call sites (`verify.py:106`, `cli.py:108`, `cli.py:387`), so swapping the backend is an import redirect. Critically, **bwa is only used to prune** — `bwa_backend.py:92-94` re-scores every candidate window in Python anyway — so a replacement backend need only return a *superset* of candidate positions, not exact alignments. `binding.py:108-125` already declares a `backend=` dispatcher whose `"bwa"` branch raises `NotImplementedError`; that is the intended registry point.

**Caveat:** the matcher must be numpy-vectorized or a small compiled helper. A pure-Python scan of 3.1 Gb will not hit these numbers.

### Step 3 — Replace the `.mmi` with a sparse anchor grid

The other consumer of discarded files is coordinate projection (`project.make_aligner`, `project.py:63`). Neither a 5.80 GB `.mmi` nor a 114-minute whole-genome PAF is necessary.

HPRC R2 haplotypes are near-chromosome-level. Measured:

| Assembly | Contigs | ≥10 Mb | Share of sequence |
|---|---:|---:|---:|
| `HG00097_hap1` | 75 | 26 | 99.3% |
| `HG01884_hap1` | 89 | 26 | 98.9% |

Instead of aligning 3 Gb of haplotype to CHM13 base-by-base, sample ~1 kb probes every 10 kb along the haplotype and map those against the **shared** CHM13 `.mmi`. That is ~300k anchors per haplotype — a few MB gzipped, buildable in an estimated 10–15 min at low RAM, versus 114 min and 14.5 GB today.

To project a locus: look up the bracketing anchors, fetch that window from the BGZF file, and realign locally with an in-memory `mappy.Aligner(seq=window, preset="asm5")` — an idiom `mask.py:48` already uses. `project.project_target` (`project.py:152-172`) already dispatches between two projection strategies, so this is a third branch rather than a rewrite.

### Step 4 — Search CHM13 once per run, not once per haplotype

With steps 1–3 in place, the per-haplotype `bwa` load is gone but each haplotype is still scanned. Step 4 removes that too.

Off-target amplification requires **both** primers binding within ~2 kb and pointing at each other. Candidate loci are therefore, by construction, homologous to the target — paralogs, pseudogenes, repeats. So:

1. Find candidate loci **once per run** on CHM13, using the shared bwa index (~47 s, amortized across every haplotype in the run instead of paid per haplotype).
2. Lift each candidate CHM13 locus onto each haplotype through the step-3 anchor grid.
3. Fetch only those windows from the BGZF file (measured ~2.6 µs per fetch) and run the existing naive scanner.

**Estimated well under 1 s per haplotype**, turning a 464-haplotype run from ~6 hours into minutes.

The many-to-one nature of the lift is a bonus rather than a compromise: a haplotype carrying three copies of a CYP2D7-like segment, all aligning to one CHM13 locus, is recovered naturally — which is exactly the CYP2D6 case the README leads with.

**Honest limitation.** Haplotype-private sequence absent from CHM13 cannot be reached this way. This fits the existing semantics rather than breaking them: the anchor grid identifies precisely which regions are unanchored, and `intro_verify_pipeline.md` §5 already argues for keeping **`uncertain` distinct from `dropout`** for exactly this class of "we could not confidently examine this locus" statement. Step 2 also remains available as an exact genome-wide path, so step 4's approximation can be validated against it on a handful of haplotypes.

---

## 4. Summary of options

Per-haplotype figures. **M** = measured, **E** = estimated.

| Approach | Storage/hap | Query/hap | 464-hap run | Precompute/hap |
|---|---:|---:|---:|---:|
| **Current** (`bwa aln` on per-hap index) | 15.10 GB **M** | 47 s **M** | ~6.1 h **E** | ~60 min **M** |
| On-the-fly `bwa index` per request | 0.90 GB | ~55 min **E** | infeasible | none |
| **Step 1+2** (BGZF + streaming scan) | 0.90 GB **E** | ~30–45 s **E** | ~4–6 h **E** | ~1 min **E** |
| **Step 1–4** (+ CHM13-once + anchor lift) | ~0.91 GB **E** | <1 s **E** | ~minutes **E** | ~15 min **E** |
| Pangenome index (r-index / GBZ) | ~0.1 GB amortized **E** | <1 s **E** | ~minutes **E** | one shared build |

### The long-term destination

If true genome-wide rigor at full scale is wanted, the principled answer is a **run-length-compressed pangenome index** — an r-index / [Movi](https://github.com/mohsenzakeri/Movi)-style structure, or `vg`'s GBZ. These compress in proportion to *distinct* sequence rather than total sequence, so 464 near-identical haplotypes fit in tens of GB as **one** index, loaded once, serving every haplotype. That is the textbook fix for this exact situation.

It is not the next step: per-haplotype coordinate reporting through a graph is a substantial rewrite, and graph mappers like `vg giraffe` are tuned for reads rather than 20-mers. Treat it as the destination, with steps 1–4 as the path that makes the current architecture viable in the meantime.

---

## 5. Sequencing and risk

**Do not delete any `.fa` or `.bwt` until the replacement search path is green.** Rebuilding a bwa index costs ~55 min per haplotype. Gate deletion on:

- `pytest -q` (in particular `tests/test_verify.py`, `tests/test_engine.py`, `tests/test_align_cache.py`)
- `pangenome-primer selftest` — all 5 statuses on the synthetic mini-pangenome
- The real-locus checks (GAPDH, CYP2D6) reproducing the current `demo-verify-pipeline/` matrix **cell for cell**, including the CYP2D7 off-target and the engineered dropout row

Other risks worth tracking:

| Risk | Mitigation |
|---|---|
| Streaming matcher is slower than estimated | Prototype the numpy rolling hash against the measured 47 s baseline *before* deleting indexes |
| Anchor-grid projection less precise than `.mmi` | Local realignment on the fetched window recovers base-level precision; diff against current output on the demo loci |
| Step 4 misses haplotype-private loci | Report as `uncertain` (existing status); validate against step 2's exact path |
| `main.nf:121` stages `['fai','mmi','amb','ann','bwt','pac','sa']` as sibling files | Must be updated in the same change, or Nextflow runs will fail to stage |

## 6. Housekeeping found along the way

Small, independent of the above:

- `prepare_haplotypes.sh:63` needs `-K 100M` — two PAF builds OOM'd at 14.5 GB (§1.3).
- Zero-byte artifacts are left behind on failure: `HG00097_hap1.fa.chm13.paf.part`, `HG01884_hap1.fa.chm13.paf.part`, `HG00639_hap1.fa.mmi`. The size>0 guard at `project.py:55-59` exists to work around exactly this — better to clean up on failure.
- `bwa_backend.ensure_index` (`bwa_backend.py:18-23`) silently launches a ~55-minute in-process `bwa index` if a FASTA arrives unindexed. It should fail loudly with a pointer to `prepare_haplotypes.sh` instead.
- The haplotype FASTA is opened twice per haplotype (`verify.py:144`, `bwa_backend.py:82`), and `align_cache.py:135` reopens it once per *primer* inside the inner loop. Hoisting that handle is free.
