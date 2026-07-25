# Known issues

Open defects with enough evidence recorded that someone else — or a later session — can pick
them up cold. Each entry states what is **established** separately from what is
**hypothesised**, because the difference decides how much of the diagnosis to trust.

**Currently open: none.** ISSUE-001 and ISSUE-002 are resolved; their records are kept below
in full, because the reasoning is the only thing that makes the fix reviewable and the
measurements are the basis for how the scanner should be sized on a bigger host.

---

## ISSUE-001 — Intermittent segfault in `pangenome_primer._scan` (stack overflow)

| | |
|---|---|
| **Status** | **Resolved** 2026-07-25. Root cause established, not inferred. |
| **Severity** | High — crashed the process; no wrong answers observed. |
| **Found** | 2026-07-25, while profiling scan cost for the Phase 5/6 review. |
| **Component** | `rust/pgp-scan` (`pangenome_primer._scan`), v0.1.0 |
| **Environment** | Linux 6.18.35.2-microsoft-standard-WSL2, 16 cores / 16 GB; Python 3.12.13; rustc 1.89.0; pyo3 0.29, rayon 1.10, flate2 1.1.9 (miniz_oxide 0.8.9) |
| **Fixed by** | `bgzf.rs`: `#[inline(never)] new_decompressor` + `map_init` reuse. `pool.rs`: explicit 16 MiB worker stacks. |
| **Regression test** | `tests/test_scan_stack_depth.py` |

### Symptom

A Python process calling `rust_backend.find_binding_sites_batch` repeatedly died with
`SIGSEGV`. The script performed a warm-up scan and then five more scans of
`hprc-r2/assemblies/HG00097_hap1.fa.gz` with 1, 2, 8, 32 and 128 primers.

```
/bin/bash: line 32: 1131932 Segmentation fault (core dumped) python - <<'EOF'
```

### Reproduction: was NOT reproducible on demand

This was the most important fact about the issue while it was open. After the initial crash:

| attempt | result |
|---|---|
| same script, unbuffered (`python -u`) | completed, 5.81–6.11 s per scan |
| same script, 3 further runs with different seeds | all completed |
| 8 consecutive scans in one process, 4 primers each | all completed |
| each primer count (1/2/8/32/128) in its own process | all completed |
| every golden gate run to date (16 passed, several runs) | no crash |

So it was **schedule-dependent, not input-dependent**, and "I ran it and it was fine" was
never evidence of absence. That is also why the fix is not certified by a clean run count —
see *Verification* below.

### Established: a stack overflow, not memory corruption

From `dmesg`:

```
python[4483]: segfault at 7dc099809c68 ip 00007dc09c78923b sp 00007dc099809c68 error 6
  in _scan.cpython-312-x86_64-linux-gnu.so[2b23b,7dc09c767000+69000]
Code: ... 49 81 eb 00 d0 00 00     sub  r11, 0xd000      ; probe target, 52 KB down
          48 81 ec 00 10 00 00     sub  rsp, 0x1000      ; step down one page
          <48> c7 04 24 00 00 00 00  mov  qword [rsp], 0  ; <-- FAULTS
          4c 39 dc                 cmp  rsp, r11
          75 ec                    jne  (loop)
RSP: 002b:00007dc099809c68
FS:  00007dc099a096c0  GS: 0000000000000000
```

1. **The faulting address equals `sp`.** A write to the stack pointer's own page failing
   means the guard page was reached.
2. **The instruction bytes are LLVM's stack-probe loop** — emitted when a function reserves a
   frame larger than one page, touching each page on the way down so the guard page is hit
   deterministically. The probe walks a frame of `0xd000` (52 KB), with a further
   `sub rsp, 0x298` after the loop.
3. **`FS:` is a thread-local base distinct from the main thread**, and the fault is inside
   `_scan...so` — a rayon pool thread.

### Established: which frame

The earlier draft of this record listed the frame as *hypothesised* and named `inflate_raw`
as a guess. It has since been attributed properly, on the exact binary that crashed.

The kernel's `[2b23b,...]` is the **file offset**, so `addr2line -e _scan...so 0x2b23b`:

```
rayon::iter::plumbing::bridge_producer_consumer::helper
```

and `objdump` at that symbol reproduces the dmesg bytes exactly, `0xd000` probe and all:

```
2b22d:  49 81 eb 00 d0 00 00   sub  $0xd000,%r11
2b234:  48 81 ec 00 10 00 00   sub  $0x1000,%rsp
2b23b:  48 c7 04 24 00 00 00   movq $0x0,(%rsp)      <-- the faulting ip
2b248:  48 81 ec 98 02 00 00   sub  $0x298,%rsp
```

Total frame **0xd298 = 53,912 bytes**. Two facts make this fatal:

* **`bridge_producer_consumer::helper` is recursive.** It is rayon's splitter: it halves the
  work and re-enters itself through `join`. A worker blocked in `join` also *steals* other
  tasks and runs them on the same stack, so the depth reached is a property of the
  work-stealing schedule — which is precisely why the same input crashed once and then
  succeeded a dozen times.
* **Rust's default spawned-thread stack is 2 MiB** and rayon does not raise it. 2 MiB ÷ 54 KB
  ≈ 38 frames. With `READ_CHUNK` = 32 MB yielding ~500 BGZF members, split depth alone is
  ~9 (≈ 486 KB); a couple of nested steals reach the guard page.

Which of the two parallel regions it was, confirmed from the calls inside that frame
(`drop_in_place<ListVecFolder<Vec<u8>>>`, `io::Error::new`): the **inflate** region,
`bgzf.rs`'s `members.par_iter()`, not `scanner.rs`'s `scan_contig`.

**Why the frame was 54 KB.** `Decompress::new(false)` → flate2 `Inflate::make` →
`InflateState::new_boxed`, which is `Box::default()`:

```rust
pub fn new_boxed(data_format: DataFormat) -> Box<InflateState> {
    let mut b: Box<InflateState> = Box::default();   // built on the STACK, then memcpy'd
```

`InflateState` is a 32 KB LZ dictionary plus three Huffman tables ≈ 44 KB, and LLVM does not
elide the temporary. Inlined into the `par_iter().map(...)` closure, that temporary became
part of `helper`'s frame.

Neither crate is at fault on its own. flate2 constructing a large value and rayon recursing
are both reasonable; the combination was ours.

### The fix

1. **`bgzf.rs` — take the frame off the recursive stack.** `new_decompressor()` is
   `#[inline(never)]`, so the 44 KB temporary lives in a leaf frame that pops immediately and
   only the boxed state crosses back. `#[inline(never)]` here is load-bearing, not a hint.
2. **`bgzf.rs` — reuse the decompressor.** `map_init(new_decompressor, ...)` gives each rayon
   job one `Decompress`, reset per member, so a ~500-member batch costs one 44 KB allocation
   per worker rather than 500.
3. **`pool.rs` — stop inheriting rayon's defaults.** A dedicated `ThreadPoolBuilder` with an
   explicit `stack_size` (16 MiB; address space, not RSS) and an explicit thread count. Kept
   even though (1) removes the pressure, because a stack limit chosen by accident is what let
   a foreseeable overflow become an intermittent SIGSEGV. `PGP_SCAN_STACK` overrides it.

After the fix, both `helper` monomorphisations have **440- and 424-byte** frames, and the
53 KB probe appears only in `pgp_scan::bgzf::new_decompressor`, a non-recursive leaf.

### Verification

A clean run count cannot certify a schedule-dependent bug — the bug already survived several
clean runs. So the fix was verified by **shrinking the stack until the difference is
deterministic** (`PGP_SCAN_STACK` exists for this):

| build | worker stack | result |
|---|---|---|
| inlined `Decompress::new` (the bug) | 256 KB | **SIGSEGV, 3 of 3 runs** |
| inlined `Decompress::new` | 2 MiB (what shipped) | completes — the intermittency |
| `#[inline(never)]` + `map_init` | 256 KB | completes |
| `#[inline(never)]` + `map_init` | 128 KB | completes |

That is a causal demonstration: same binary, same input, same thread count, one line of
attribute changed. Against the shipped 16 MiB stack the fixed build has ≥128× headroom.

`tests/test_scan_stack_depth.py` pins it at a 256 KB stack, and was itself checked to **fail
on the buggy build** (returncode −11) rather than merely pass on the fixed one.

Also confirmed: 50 consecutive runs of the original triggering script, `pytest -q`
(150 passed), and the golden gate unchanged.

### Why it mattered more than it looked

`verify.run_verify` issues **one scan call per haplotype in a single process**. The exposure
actually exercised was 3 calls. A 464-haplotype run is 150× that, and the crash would land
after hours of work with no partial output. A defect that is rare at demo scale is not rare
at production scale — it is merely untested there.

---

## ISSUE-002 — `rayon` thread pool is uncapped; **blocks deployment to many-core hosts**

| | |
|---|---|
| **Status** | **Resolved** 2026-07-25 (same fix site as ISSUE-001). |
| **Severity** | Low on the laptop, **High on any server/cloud host.** Correctness unaffected; throughput was not. |
| **Found** | 2026-07-25, while sizing the pipeline for a 256-core / 700 GB server. |
| **Component** | `rust/pgp-scan` (`pangenome_primer._scan`), v0.1.0 |
| **Fixed by** | `rust/pgp-scan/src/pool.rs`; `search.threads` in `config/defaults.yaml` |

### The defect

The crate never constructed a `rayon::ThreadPoolBuilder`. Both parallel regions
(`bgzf.rs` inflate, `scanner.rs` scan) ran on rayon's **global pool, which defaults to one
worker per logical CPU**. Nothing in the crate, the CLI, or `config/defaults.yaml` capped or
exposed it.

**This was benign at the time, which is exactly why it was easy to miss.**
`verify.run_verify` loops haplotypes sequentially, so on the 16-core laptop exactly one scan
runs at a time and taking all 16 cores is the *right* behaviour. The default was correct for
the only configuration it had ever been run in.

It became a defect the moment either changed:

1. **haplotypes are processed in parallel** (the Nextflow `EVALUATE` fan-out already does
   this); or
2. **the host has many cores.**

On a 256-core server running 64 concurrent haplotypes, each worker process spawns 256 rayon
threads: **64 × 256 = 16,384 threads competing for 256 cores.** Not a suboptimal setting, a
thrashing one — the large server can finish behind the laptop.

### Measured: why a low cap is right for throughput

Single scan of `HG00097_hap1.fa.gz` (3.03 Gbp), varying `RAYON_NUM_THREADS`:

| threads | wall | speedup | core-seconds | efficiency | peak RSS |
|---:|---:|---:|---:|---:|---:|
| 1 | 20.05 s | 1.00× | 19.8 | 100% | 397 MB |
| 2 | 12.23 s | 1.64× | 21.0 | 82% | 408 MB |
| 4 | 8.10 s | 2.48× | 23.1 | 62% | 399 MB |
| 6 | 6.88 s | 2.91× | 25.9 | 49% | 401 MB |
| 8 | 7.95 s | 2.52× | 36.0 | 32% | 404 MB |
| 12 | 6.83 s | 2.94× | 41.3 | 24% | 406 MB |
| 16 | 5.84 s | 3.43× | 44.4 | **21%** | 410 MB |

At 16 threads the scan burns **2.2× the CPU to go 3.4× faster**. Fine for a single scan;
across 464 haplotypes it wastes more than half the machine.

**Peak RSS is flat (~400 MB) at every thread count**, so memory never justifies a high cap —
the constraint is purely CPU efficiency.

End-to-end, 8 real haplotypes on 16 cores:

| configuration | wall | vs default |
|---|---:|---:|
| default (16 threads) × 2 concurrent | 43.93 s | — |
| **4 threads × 4 concurrent** | **36.82 s** | **19% faster** |
| 2 threads × 8 concurrent | 38.45 s | 14% faster |
| 1 thread × 8 concurrent | 71.60 s | worse (only 8 of 16 cores occupied) |

### Recommended settings

**These now live in [`sizing.md`](sizing.md)**, together with the Nextflow sizing profiles and `scripts/sizing_sweep.sh` for re-measuring on a target host. Kept here for the record.

Rule of thumb: **threads × concurrent-haplotypes ≈ core count**, threads in the 2–4 range.

| scenario | `search.threads` | concurrent haplotypes | peak RAM |
|---|---|---|---|
| 16-core laptop, CLI as it is today (sequential loop) | `null` (→16) | 1 | 0.4 GB |
| 16-core laptop, once haplotypes run in parallel | **4** | 4 | 1.6 GB |
| **256-core server, many haplotypes** | **4** | **64** | ~26 GB |
| 256-core server, single haplotype | 8–16 — **never 256** | 1 | 0.4 GB |

For the anchor-grid build the equivalent knob is minimap2's `-t` (`MM_THREADS` in
`scripts/prepare_haplotypes.sh`), and there the binding constraint is RAM, not cores: ~8.3 GB
per concurrent build means **1 at a time on a 16 GB laptop** but ~64 on a 700 GB server — the
difference between ~62 h and ~1 h for 464 haplotypes.

### The fix

`pool.rs` resolves the worker count in this order, and builds one pool for the process:

1. the `threads` argument — `search.threads` in `config/defaults.yaml`;
2. `PGP_SCAN_THREADS`;
3. `RAYON_NUM_THREADS` (kept working: it is the knob anyone tuning a rayon program reaches
   for, and silently ignoring it would be its own defect);
4. `min(available_parallelism, 16)`.

Rule 4 is the part that closes the issue: **the default no longer scales with host core
count.** It leaves the 16-core laptop exactly where it was — a solo scan still gets all 16,
and capping lower would be a real regression (5.84 s at 16 vs 7.95 s at 8) — while making a
256-core host behave like a 16-core one.

`search.threads` is `null` by default, which is correct only while the CLI scans haplotypes
one at a time. **Set it when scans run concurrently.**

The pool is built once per process and is not resized; a later, different `threads` is
ignored. `rust_backend` emits a `RuntimeWarning` when that happens rather than dropping it
silently, and `rust_backend.pool_threads()` reports the live count.

### Measurement caveat — re-measure on the target host

All figures above come from **16 logical cores on a WSL2 laptop, most likely 8 physical cores
plus hyperthreading**. The anomalous dip at 8 threads (7.95 s, slower than 6.88 s at 6) has
the shape of an HT/physical boundary. A many-core server will have different NUMA topology
and a different curve; the optimum could reasonably be 8 rather than 4.

**Re-run this sweep on the target machine before fixing a value.** It takes about two minutes
and the numbers here should not be copied to a server unexamined.

In Nextflow, `EVALUATE` exports `PGP_SCAN_THREADS=${task.cpus}`, so the `cpus` directive in
`nextflow.config` *is* the scan-parallelism knob and a fan-out cannot oversubscribe by
construction. It is set to 4, and `EVALUATE`'s memory reservation dropped 8 GB → 2 GB
(measured peak RSS is ~400 MB; the scanner streams the BGZF and holds no index), as did
`PROJECT`'s, which was sized for the 5.80 GB per-haplotype `.mmi` that the ~4 MB anchor grid
replaced.

### Still open as a follow-up (not a defect)

* **De-nest the two parallel regions.** The scan uses **738% of a possible 1600% CPU**
  because inflate and scan alternate rather than overlap. Pipelining them is the largest
  remaining performance item.
