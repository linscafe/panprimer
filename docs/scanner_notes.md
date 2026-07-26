# Scanner notes — why `rust/pgp-scan` looks the way it does

Two defects (ISSUE-001, ISSUE-002) left behind constructs that read as optional and are not.
This file is the evidence. **Issue status and discussion live in the GitHub tracker; both are
closed.** What is kept here is what the tracker is a bad home for: reasoning a future reader
needs *in the clone*, reachable by `grep`, when deciding whether a line can be deleted.

| construct | where | deleting it costs |
|---|---|---|
| `#[inline(never)]` on `new_decompressor` | `bgzf.rs` | reinstates ISSUE-001 (intermittent SIGSEGV) |
| `map_init` decompressor reuse | `bgzf.rs` | 500 × 44 KB allocations per batch |
| explicit `stack_size` / thread count | `pool.rs` | reinstates both defects' preconditions |
| `PGP_SCAN_STACK` | `pool.rs` | removes the only falsifiability handle on ISSUE-001 |
| `d.reset(false)` | `bgzf.rs` | mis-inflates the next member in a job |

Measured on Linux 6.18.35 (WSL2), 16 cores / 16 GB; Python 3.12.13; rustc 1.89.0; pyo3 0.29,
rayon 1.10, flate2 1.1.9 (miniz_oxide 0.8.9). Regression tests: `tests/test_scan_stack_depth.py`.

---

## ISSUE-001 — Intermittent segfault in `_scan` (stack overflow)

### It was not reproducible on demand

This was the most important fact while the issue was open. After the initial `SIGSEGV`:

| attempt | result |
|---|---|
| same script, unbuffered | completed, 5.81–6.11 s per scan |
| 3 further runs, different seeds | all completed |
| 8 consecutive scans in one process | all completed |
| each primer count in its own process | all completed |
| every golden gate run | no crash |

So it was **schedule-dependent, not input-dependent**, and "I ran it and it was fine" was
never evidence of absence — which is also why the fix is not certified by a clean run count.

### Established: a stack overflow

From `dmesg`:

```
python[4483]: segfault at 7dc099809c68 ip 00007dc09c78923b sp 00007dc099809c68 error 6
  in _scan.cpython-312-x86_64-linux-gnu.so[2b23b,7dc09c767000+69000]
Code: ... 49 81 eb 00 d0 00 00     sub  r11, 0xd000      ; probe target, 52 KB down
          48 81 ec 00 10 00 00     sub  rsp, 0x1000      ; step down one page
          <48> c7 04 24 00 00 00 00  mov  qword [rsp], 0  ; <-- FAULTS
```

The faulting address **equals `sp`**, so the guard page was reached. The bytes are LLVM's
**stack-probe loop**, emitted for frames larger than a page. `FS:` shows a thread-local base
distinct from the main thread — a rayon worker.

### Established: which frame

The kernel's `[2b23b,...]` is a **file offset**, not an ip delta. `addr2line` on it gives:

```
rayon::iter::plumbing::bridge_producer_consumer::helper
```

and `objdump` there reproduces the dmesg bytes exactly, `0xd000` probe and all — total frame
**0xd298 = 53,912 bytes**. Two facts make that fatal:

- **`helper` is recursive.** It is rayon's splitter, re-entering itself through `join`. A
  worker blocked in `join` also *steals* tasks and runs them on the same stack, so depth is a
  property of the work-stealing schedule — precisely why the same input crashed once and then
  succeeded a dozen times.
- **Rust's default spawned-thread stack is 2 MiB** and rayon does not raise it. 2 MiB ÷ 54 KB
  ≈ 38 frames. Split depth alone is ~9; a couple of nested steals reach the guard page.

Calls inside the frame (`drop_in_place<ListVecFolder<Vec<u8>>>`, `io::Error::new`) identify it
as the **inflate** region, `bgzf.rs`'s `members.par_iter()`, not `scanner.rs`.

**Why 54 KB.** `Decompress::new` → `InflateState::new_boxed`, which is `Box::default()`:

```rust
let mut b: Box<InflateState> = Box::default();   // built on the STACK, then memcpy'd
```

`InflateState` is a 32 KB LZ dictionary plus Huffman tables ≈ 44 KB, and LLVM does not elide
the temporary. Inlined into the `par_iter().map(...)` closure, it became part of `helper`'s
frame. Neither crate is wrong alone — flate2 constructing a large value and rayon recursing
are both reasonable. The combination was ours.

### The fix

1. **Take the frame off the recursive stack.** `new_decompressor()` is `#[inline(never)]`, so
   the 44 KB temporary lives in a leaf frame that pops immediately and only the boxed state
   crosses back. Load-bearing, not a hint.
2. **Reuse the decompressor.** `map_init(new_decompressor, ...)` gives each rayon job one
   `Decompress`, reset per member — one 44 KB allocation per worker, not 500.
3. **Stop inheriting rayon's defaults.** A dedicated `ThreadPoolBuilder` with explicit
   `stack_size` (16 MiB of address space, not RSS) and thread count. Kept even though (1)
   removes the pressure, because a stack limit chosen by accident is what turned a foreseeable
   overflow into an intermittent SIGSEGV. `PGP_SCAN_STACK` overrides it.

After the fix both `helper` monomorphisations have **440- and 424-byte** frames, and the 53 KB
probe survives only in `new_decompressor`, a non-recursive leaf.

### Verification

A clean run count cannot certify a schedule-dependent bug — this one already survived several.
So the fix was verified by **shrinking the stack until the difference is deterministic**:

| build | worker stack | result |
|---|---|---|
| inlined `Decompress::new` (the bug) | 256 KB | **SIGSEGV, 3 of 3 runs** |
| inlined `Decompress::new` | 2 MiB (shipped) | completes — the intermittency |
| `#[inline(never)]` + `map_init` | 256 KB | completes |
| `#[inline(never)]` + `map_init` | 128 KB | completes |

Same binary, same input, same thread count, one attribute changed. Against the shipped 16 MiB
stack the fixed build has ≥128× headroom. The regression test pins a 256 KB stack and was
itself checked to **fail on the buggy build** (returncode −11), not merely pass on the fixed
one. Also confirmed: 50 consecutive runs of the triggering script, and the golden gate.

### Why it mattered more than it looked

`verify.run_verify` issues one scan call per haplotype in a single process. The exposure
actually exercised was 3 calls; a 464-haplotype run is 150× that, and the crash would land
after hours with no partial output. A defect that is rare at demo scale is not rare at
production scale — it is merely untested there.

---

## ISSUE-002 — uncapped `rayon` pool blocks many-core hosts

Correctness was never affected; throughput on a large host was. Same fix site (`pool.rs`),
plus `search.threads` in `config/defaults.yaml`.

### The defect

The crate never constructed a `ThreadPoolBuilder`, so both parallel regions ran on rayon's
**global pool — one worker per logical CPU**, with nothing in the crate, CLI or config
capping or exposing it.

**This was benign at the time, which is why it was easy to miss.** `verify.run_verify` loops
haplotypes sequentially, so on a 16-core laptop exactly one scan runs and taking all 16 cores
is *correct*. The default was right for the only configuration it had ever run in.

It became a defect the moment haplotypes ran in parallel (the Nextflow `EVALUATE` fan-out
does) or the host got big. On a 256-core server running 64 concurrent haplotypes, each process
spawns 256 rayon threads: **64 × 256 = 16,384 threads for 256 cores.** Not suboptimal —
thrashing. The large server can finish behind the laptop.

### Measured: why a low cap is right

Single scan of a 3.03 Gbp haplotype, varying `RAYON_NUM_THREADS`:

| threads | wall | speedup | core-seconds | efficiency | peak RSS |
|---:|---:|---:|---:|---:|---:|
| 1 | 20.05 s | 1.00× | 19.8 | 100% | 397 MB |
| 2 | 12.23 s | 1.64× | 21.0 | 82% | 408 MB |
| 4 | 8.10 s | 2.48× | 23.1 | 62% | 399 MB |
| 6 | 6.88 s | 2.91× | 25.9 | 49% | 401 MB |
| 8 | 7.95 s | 2.52× | 36.0 | 32% | 404 MB |
| 12 | 6.83 s | 2.94× | 41.3 | 24% | 406 MB |
| 16 | 5.84 s | 3.43× | 44.4 | **21%** | 410 MB |

At 16 threads the scan burns **2.2× the CPU to go 3.4× faster** — fine alone, wasteful across
464 haplotypes. **Peak RSS is flat (~400 MB)**, so memory never justifies a high cap.

End-to-end, 8 real haplotypes on 16 cores:

| configuration | wall | vs default |
|---|---:|---:|
| default (16 threads) × 2 concurrent | 43.93 s | — |
| **4 threads × 4 concurrent** | **36.82 s** | **19% faster** |
| 2 threads × 8 concurrent | 38.45 s | 14% faster |
| 1 thread × 8 concurrent | 71.60 s | worse (only 8 of 16 cores busy) |

> **threads × concurrent-haplotypes ≈ core count**, threads in the **2–4** range.

Everything derived from these numbers — per-host recommendations, Nextflow profiles, the
anchor-grid RAM constraint, `scripts/sizing_sweep.sh` — lives in [`sizing.md`](sizing.md) and
is **deliberately not restated here**. An earlier revision duplicated it and the copies
drifted. This file keeps the measurements; `sizing.md` cites them as their source.

### The fix

`pool.rs` resolves the worker count in order, building one pool per process:

1. the `threads` argument (`search.threads`);
2. `PGP_SCAN_THREADS`;
3. `RAYON_NUM_THREADS` — kept working, because it is the knob anyone tuning a rayon program
   reaches for, and silently ignoring it would be its own defect;
4. `min(available_parallelism, 16)`.

Rule 4 closes the issue: **the default no longer scales with host core count.** The laptop is
unchanged — a solo scan still gets all 16, and capping lower would be a real regression
(5.84 s at 16 vs 7.95 s at 8) — while a 256-core host now behaves like a 16-core one.

The pool is built once per process and never resized; a later, different `threads` is ignored
with a `RuntimeWarning` rather than in silence, and `rust_backend.pool_threads()` reports the
live count.

In Nextflow, `EVALUATE` exports `PGP_SCAN_THREADS=${task.cpus}`, so the `cpus` directive *is*
the scan-parallelism knob and a fan-out cannot oversubscribe by construction. `EVALUATE`'s
memory reservation dropped 8 GB → 2 GB (measured peak RSS ~400 MB), as did `PROJECT`'s, which
had been sized for the 5.80 GB per-haplotype `.mmi` the ~4 MB anchor grid replaced.

### Measurement caveat

All figures come from **16 logical cores on a WSL2 laptop, probably 8 physical plus
hyperthreading**. The dip at 8 threads (7.95 s, slower than 6.88 s at 6) has the shape of an
HT/physical boundary. A server will have different NUMA topology and a different curve — the
optimum could reasonably be 8 rather than 4. **Re-run the sweep before fixing a value**; it
takes two minutes.

### Open follow-up (not a defect)

**De-nest the two parallel regions.** The scan uses **738% of a possible 1600% CPU** because
inflate and scan alternate rather than overlap. Pipelining them is the largest remaining
performance item.
