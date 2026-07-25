# Known issues

Open defects with enough evidence recorded that someone else — or a later session — can pick
them up cold. Each entry states what is **established** separately from what is
**hypothesised**, because the difference decides how much of the diagnosis to trust.

---

## ISSUE-001 — Intermittent segfault in `pangenome_primer._scan` (stack overflow)

| | |
|---|---|
| **Status** | Open. Not yet fixed. |
| **Severity** | High — crashes the process; no wrong answers observed. |
| **Found** | 2026-07-25, while profiling scan cost for the Phase 5/6 review. |
| **Component** | `rust/pgp-scan` (`pangenome_primer._scan`), v0.1.0 |
| **Environment** | Linux 6.18.35.2-microsoft-standard-WSL2, 16 cores / 16 GB; Python 3.12.13; rustc 1.89.0; pyo3 0.29, rayon 1.10, flate2 1.0 (miniz_oxide) |

### Symptom

A Python process calling `rust_backend.find_binding_sites_batch` repeatedly died with
`SIGSEGV`. The script performed a warm-up scan and then five more scans of
`hprc-r2/assemblies/HG00097_hap1.fa.gz` with 1, 2, 8, 32 and 128 primers.

```
/bin/bash: line 32: 1131932 Segmentation fault (core dumped) python - <<'EOF'
```

### Reproduction: currently NOT reproducible on demand

This is the most important fact about the issue. After the initial crash:

| attempt | result |
|---|---|
| same script, unbuffered (`python -u`) | completed, 5.81–6.11 s per scan |
| same script, 3 further runs with different seeds | all completed |
| 8 consecutive scans in one process, 4 primers each | all completed |
| each primer count (1/2/8/32/128) in its own process | all completed |
| every golden gate run to date (16 passed, several runs) | no crash |

So it is **schedule-dependent, not input-dependent**. Do not treat "I ran it and it was
fine" as evidence the bug is absent.

### Established: this is a stack overflow, not memory corruption

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

Three things follow, and all three are solid:

1. **The faulting address equals `sp`** (`segfault at 7dc099809c68`, `sp 00007dc099809c68`).
   A write to the stack pointer's own page failing means the guard page was reached.
2. **The instruction bytes are LLVM's stack-probe loop** — the sequence emitted when a
   function reserves a frame larger than one page, touching each page on the way down so the
   guard page is hit deterministically rather than skipped. The probe is walking a frame of
   at least `0xd000` (52 KB), with a further `sub rsp, 0x298` after the loop.
3. **`FS:` is set to a thread-local base distinct from the main thread**, and the fault is
   inside `_scan...so`, so this is a worker thread in the extension — a rayon pool thread.

Together: a rayon worker exhausted its stack. This is *not* a data race, use-after-free, or
heap corruption, and it should not be investigated as one.

### Hypothesised: which frame, and why it only sometimes overflows

**This part is inference and has not been confirmed.**

Two rayon parallel regions exist and they are adjacent in the call graph:

- `bgzf.rs:175` — `members.par_iter().map(|..| inflate_raw(..))`, inflating BGZF blocks.
- `scanner.rs:292` — `(0..n_chunks).into_par_iter()` inside `scan_contig`.

They are chained through the streaming sink: `lib.rs:112` calls `bgzf::stream(path, |chunk|
parser.feed(chunk))`, and the parser calls `index.scan_contig(&c.codes)` at `lib.rs:109`
when a contig completes.

`inflate_raw` (`bgzf.rs:60`) constructs `flate2::Decompress::new(false)` per block. flate2's
pure-Rust backend holds a large inflate state (Huffman tables plus a 32 KB LZ window), which
matches the ~52 KB probe. Rust's default spawned-thread stack is 2 MiB, and rayon does not
raise it. Under recursive splitting and work-stealing a single worker can accumulate several
such frames, and *how many* depends on the steal schedule — which would explain why the same
input crashes once and then succeeds a dozen times.

**Not established:** that `inflate_raw` is the frame in question. Nobody has attributed
`ip 0x7dc09c78923b` to a symbol. Do that first — `addr2line`/`objdump` against the built
`.so`, or reproduce under `rust-gdb` — rather than trusting the paragraph above.

### Why this matters more than it looks

`verify.run_verify` issues **one scan call per haplotype in a single process**. The exposure
that has actually been exercised is 3 calls. A 464-haplotype run is 150× that, and the crash
would land after hours of work with no partial output. A defect that is rare at demo scale is
not rare at production scale — it is merely untested there.

No incorrect results have been attributed to this. The failure mode observed is a hard crash,
which is the benign direction.

### Fix directions, cheapest first

1. **Give the pool a bigger stack.** Build a dedicated
   `rayon::ThreadPoolBuilder::new().stack_size(N)` and run both parallel regions inside it.
   Addresses stack exhaustion regardless of which frame is responsible; does not require
   root-causing first. Lowest risk, and the right move if a fix is needed before the
   attribution work is done.
2. **Heap-allocate the decompressor / hoist it out of the hot path.** Reusing one
   `Decompress` per worker via `reset()` removes both the large frame and the per-block
   allocation. Likely a small throughput win as well.
3. **De-nest the two parallel regions.** Worth doing independently: the scan currently uses
   **738% of a possible 1600% CPU** because inflate and scan alternate rather than overlap
   (see Phase 6 notes in `implementation-plan.md`). Pipelining them addresses the stack
   depth *and* is the largest remaining performance item.

### Verification before closing

- Attribute the faulting `ip` to a symbol; record it here.
- Run the triggering profile script **50×** without a crash. One clean run proves nothing —
  the bug already survived several.
- `pytest -q` and the backend conformance + differential suites unchanged.
- Re-measure scan wall time; fixes 2 and 3 should not regress it, and may improve it.

---

## ISSUE-002 — `rayon` thread pool is uncapped; **blocks deployment to many-core hosts**

| | |
|---|---|
| **Status** | Open. Not yet fixed. |
| **Severity** | Low today, **High on any server/cloud host.** Correctness is unaffected; throughput is not. |
| **Found** | 2026-07-25, while sizing the pipeline for a 256-core / 700 GB server. |
| **Component** | `rust/pgp-scan` (`pangenome_primer._scan`), v0.1.0 |
| **Trigger** | Host with many cores **and** more than one scan running concurrently. |

### The defect

The crate never constructs a `rayon::ThreadPoolBuilder`. Both parallel regions
(`bgzf.rs:175` inflate, `scanner.rs:292` scan) therefore run on rayon's **global pool, which
defaults to one worker per logical CPU**. Nothing in the crate, the CLI, or
`config/defaults.yaml` caps or exposes it.

**This is currently benign, and that is exactly why it is easy to miss.** `verify.run_verify`
loops haplotypes sequentially (`for hid, fasta in haplos:`), so on the 16-core development
laptop exactly one scan runs at a time and taking all 16 cores is the *right* behaviour. The
default is correct for the only configuration it has ever been run in.

It becomes a defect the moment either of these changes:

1. **haplotypes are processed in parallel** (the Nextflow `EVALUATE` fan-out already does
   this, and parallelising the CLI loop is the obvious next optimisation); or
2. **the host has many cores.**

On a 256-core server running 64 concurrent haplotypes, each worker process spawns 256 rayon
threads: **64 × 256 = 16,384 threads competing for 256 cores.** That is not a suboptimal
setting, it is a thrashing one — the large server can end up *slower* than the laptop.

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

At 16 threads the scan burns **2.2× the CPU to go 3.4× faster**. For a single scan that is a
fine trade; across 464 haplotypes it wastes more than half the machine.

**Peak RSS is flat (~400 MB) regardless of thread count**, so memory never justifies a high
cap — the constraint is purely CPU efficiency.

End-to-end, 8 real haplotypes on 16 cores:

| configuration | wall | vs default |
|---|---:|---:|
| default (16 threads) × 2 concurrent | 43.93 s | — |
| **4 threads × 4 concurrent** | **36.82 s** | **19% faster** |
| 2 threads × 8 concurrent | 38.45 s | 14% faster |
| 1 thread × 8 concurrent | 71.60 s | worse (only 8 of 16 cores occupied) |

### Recommended settings

Rule of thumb: **threads × concurrent-haplotypes ≈ core count**, threads in the 2–4 range.

| scenario | `RAYON_NUM_THREADS` | concurrent haplotypes | peak RAM |
|---|---|---|---|
| 16-core laptop, CLI as it is today (sequential loop) | unset (16) — currently correct | 1 | 0.4 GB |
| 16-core laptop, once haplotypes run in parallel | **4** | 4 | 1.6 GB |
| **256-core server, many haplotypes** | **4** | **64** | ~26 GB |
| 256-core server, single haplotype | 8–16 — **never 256** | 1 | 0.4 GB |

For the anchor-grid build, the equivalent knob is minimap2's `-t` (`MM_THREADS` in
`scripts/prepare_haplotypes.sh`), and there the binding constraint is RAM, not cores: ~8.3 GB
per concurrent build means **1 at a time on a 16 GB laptop** but ~64 on a 700 GB server —
the difference between ~62 h and ~1 h for 464 haplotypes.

### Measurement caveat — re-measure on the target host

All figures above come from **16 logical cores on a WSL2 laptop, most likely 8 physical cores
plus hyperthreading**. The anomalous dip at 8 threads (7.95 s, slower than 6.88 s at 6) has
the shape of an HT/physical boundary. A many-core server will have different NUMA topology
and a different curve; the optimum could reasonably be 8 rather than 4.

**Re-run this sweep on the target machine before fixing a value.** It takes about two
minutes and the numbers here should not be copied to a server unexamined.

### Fix directions

1. **Set an explicit default inside the crate** rather than inheriting rayon's. A dedicated
   `ThreadPoolBuilder` with a sane cap makes behaviour independent of host core count, and is
   the same change ISSUE-001 needs for `stack_size` — do both at once.
2. **Expose it as configuration**, e.g. `search.threads` in `config/defaults.yaml`, so it is
   not reachable only through an environment variable that no documentation mentions.
   `RAYON_NUM_THREADS` should remain an override.
3. **Set it per task in `main.nf`**, derived from the Nextflow `cpus` directive, so a cluster
   run cannot oversubscribe by construction.
4. Document the threads × concurrency rule wherever the scale-out path is described.

### Verification before closing

- On a many-core host, confirm total thread count stays near the core count with N concurrent
  workers (`ps -eLf | wc -l`), not N × cores.
- Reproduce the 8-haplotype configuration sweep on that host and record the curve here.
- Confirm the laptop single-haplotype path does not regress: a solo scan should still use all
  available cores unless explicitly capped.
