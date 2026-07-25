# Sizing the pipeline for your machine

Every default in this repo is tuned for the machine it was developed on: **16 logical cores,
16 GB RAM**. Most of them are wrong for a bigger host, and one of them — the scanner's thread
count — used to be wrong in a way that makes a 256-core server *slower than the laptop*
(ISSUE-002 in [`issues.md`](issues.md)).

Nothing here is auto-detected on purpose. Auto-tuning would guess, and a wrong guess on a
shared cluster node is expensive and silent. The knobs are explicit, they have defaults that
are safe on a small machine, and this page says what to change.

> **Measure, don't copy.** Every number below came from one WSL2 laptop that is probably
> 8 physical cores plus hyperthreading. `scripts/sizing_sweep.sh <hap.fa.gz>` reproduces the
> whole measurement on your host in about two minutes. Run it before committing a value.

---

## The one rule

> **threads × concurrent-haplotypes ≈ core count**, with threads in the **2–4** range.

The scanner parallelises *within* one haplotype, and the pipeline parallelises *across*
haplotypes. Both draw on the same cores, so the two numbers multiply. Getting this wrong in
the generous direction is what produces 16,384 threads on 256 cores.

Why threads stays low, measured on 16 cores with a 3.03 Gbp haplotype:

| threads | wall | speedup | efficiency |
|---:|---:|---:|---:|
| 1 | 20.05 s | 1.00× | 100% |
| 2 | 12.23 s | 1.64× | 82% |
| 4 | 8.10 s | 2.48× | 62% |
| 8 | 7.95 s | 2.52× | 32% |
| 16 | 5.84 s | 3.43× | **21%** |

At 16 threads a scan burns **2.2× the CPU to go 3.4× faster**. Fine when it is the only thing
running; across 464 haplotypes it wastes more than half the machine. End-to-end on 16 cores,
**4 threads × 4 concurrent beat 16 threads × 2 concurrent by 19%**.

**Peak RSS is ~400 MB at every thread count** — the scanner streams the BGZF and builds no
index — so memory never argues for a high thread count. RAM matters for a different job
entirely: see *Anchor-grid builds* below.

---

## Recommended settings

| | **16-core laptop**<br>16 GB | **~64-core server**<br>256 GB | **~256-core server**<br>700 GB |
|---|---|---|---|
| Nextflow profile | `-profile local,laptop16` | `-profile local,server64` | `-profile local,server256` |
| `scan_threads` (= EVALUATE `cpus`) | 4 | 4 | 4 |
| concurrent EVALUATE tasks | 4 | 16 | 64 |
| `scan_memory` | 2 GB | 2 GB | 2 GB |
| `grid_threads` (minimap2 `-t`) | 4 | 8 | 8 |
| concurrent anchor-grid builds | **1** | 16 | 64 |
| 464-haplotype grid build | ~62 h | ~4 h | ~1 h |

Single-haplotype runs are the exception: there is nothing to run concurrently, so give the
scan more threads — 8–16, and **never 256**, because efficiency is already 21% at 16.

The scanner's own default, when nothing sets a value, is `min(cores, 16)`. That leaves a
16-core laptop exactly as it was while stopping a 256-core host from spawning 256 workers.

---

## Where each knob lives

### Nextflow (the fan-out path — this is the one that matters at scale)

Sizing profiles compose with execution profiles:

```bash
nextflow run main.nf -profile local,server256 --target chr12:… --chm13 …
```

Override individually without a profile:

```bash
nextflow run main.nf -profile local --scan_threads 8 --scan_memory '4 GB'
```

`EVALUATE` declares `cpus = params.scan_threads` and exports `PGP_SCAN_THREADS=${task.cpus}`.
That is deliberate: with the local executor Nextflow will not start more tasks than the
declared `cpus` fit, so **threads × concurrency is bounded by construction** — one knob, not
two that must be kept consistent. On a cluster executor the same `cpus` value is what the
scheduler allocates, so the guarantee carries over.

Profiles available: `laptop16`, `server64`, `server256`, `bigmem` (fewer, fatter tasks —
8 threads, 4 GB). They set only params, so a `--scan_threads` on the command line still wins.

### The CLI (sequential — the default is already right)

`config/defaults.yaml`:

```yaml
search:
  threads: null    # null = all cores, capped at 16
```

Leave it `null` while `pangenome-primer run`/`verify` loop haplotypes one at a time: a solo
scan *should* take the machine. Set it only if you run several CLI processes at once.

`PGP_SCAN_THREADS` and `RAYON_NUM_THREADS` both override it, in that order.

> The scanner's pool is built **once per process**. The first scan fixes the count; a later
> different value is ignored, with a `RuntimeWarning` rather than in silence.
> `rust_backend.pool_threads()` reports what is actually live.

### Anchor-grid builds (RAM-bound, not CPU-bound)

`scripts/prepare_haplotypes.sh` honours `MM_THREADS` (default 4), and
`pangenome-primer build-anchor-grid` takes `--threads`.

The binding constraint here is **~8.3 GB peak per concurrent build**, not cores:

| host | concurrent builds | why |
|---|---|---|
| 16 GB laptop | **1** | a second build OOMs |
| 256 GB server | ~16 | 16 × 8.3 GB ≈ 133 GB, with headroom |
| 700 GB server | ~64 | ~530 GB |

The script builds haplotypes sequentially. To parallelise, split the manifest and run several
instances — each haplotype is independent and every step is skipped when its output already
exists, so this is safe and resumable:

```bash
split -n l/8 --numeric-suffixes=1 --additional-suffix=.tsv config/samples.tsv part
for f in part*.tsv; do MM_THREADS=8 scripts/prepare_haplotypes.sh "$f" & done; wait
```

Keep `instances × MM_THREADS` under the core count **and** `instances × 8.3 GB` under RAM.
RAM is almost always the one that binds first.

---

## Re-measuring

```bash
scripts/sizing_sweep.sh hprc-r2/assemblies/HG00097_hap1.fa.gz
```

Part 1 sweeps threads for a single scan (the number for a one-haplotype run). Part 2 sweeps
`threads × concurrency` shapes at constant total work (the number for a real multi-haplotype
job) — these usually differ, and Part 2 is the one to apply. The script ends by printing the
exact flag for each entry point.

Record what you measure here, so the next person inherits a number rather than a method.

| host | date | best single-scan threads | best fan-out shape | source |
|---|---|---|---|---|
| 16-core / 16 GB WSL2 laptop | 2026-07-25 | 16 (5.84 s) | 4 threads × 4 concurrent | ISSUE-002 |
| 16-core / 16 GB WSL2 laptop | 2026-07-25 | 16 (5.89 s) | 4 threads × 4 concurrent | `sizing_sweep.sh` |

The second row is the script re-deriving the first independently — 1 thread 19.95 s vs 20.05 s,
4 threads 7.94 s vs 8.10 s, efficiency 100/63/21% vs 100/62/21%. Reassuring for the method,
and a reminder that ±3% run-to-run is the noise floor: do not chase differences that small.
