//! The crate's own rayon pool: an explicit thread count and an explicit stack size.
//!
//! Both parallel regions (`bgzf::stream`'s inflate, `scanner::scan_contig`) used to run on
//! rayon's *global* pool, which the crate never configured. That inherited two defaults, and
//! each was a filed defect:
//!
//! * **one worker per logical CPU** (ISSUE-002). Correct on the 16-core laptop, where
//!   `verify.run_verify` loops haplotypes sequentially and one scan at a time should take the
//!   machine. Wrong the moment haplotypes run concurrently -- the Nextflow `EVALUATE` fan-out
//!   already does -- or the host is large: 64 workers on a 256-core server each spawning 256
//!   rayon threads is 16,384 threads over 256 cores, i.e. thrashing, and the big server can
//!   finish behind the laptop.
//! * **a 2 MiB worker stack** (ISSUE-001), which a recursive `bridge_producer_consumer::helper`
//!   carrying a ~54 KB frame could exhaust depending on the work-stealing schedule.
//!
//! `bgzf.rs` fixes the frame that made the stack the binding constraint; `STACK_SIZE` here is
//! deliberately kept as well, because a stack limit chosen by accident is what let a
//! foreseeable overflow become an intermittent SIGSEGV in the first place.
//!
//! # The default is a cap, not a target
//!
//! Resolution order, highest priority first:
//!
//! 1. the `threads` argument (from `search.threads` in `config/defaults.yaml`);
//! 2. `PGP_SCAN_THREADS`;
//! 3. `RAYON_NUM_THREADS`, which is what anyone reaching for rayon's usual knob will set;
//! 4. `min(available_parallelism, DEFAULT_MAX_THREADS)`.
//!
//! Rule 4 leaves the 16-core laptop exactly where it was (a solo scan still gets all 16 --
//! measured 5.84 s at 16 threads against 7.95 s at 8, so capping lower would be a real
//! regression) while making a 256-core host behave like a 16-core one instead of spawning 256
//! workers for a job whose parallel efficiency is already 21% at 16.
//!
//! It is a ceiling for the *single-scan* case only. Running many haplotypes at once wants a
//! much lower number -- measured optimum on 16 cores was 4 threads x 4 concurrent, 19% faster
//! than the default -- and that is what `search.threads` is for. See ISSUE-002 in
//! `docs/scanner_notes.md` for the measured curves and the threads x concurrency rule.

use std::sync::OnceLock;

use rayon::{ThreadPool, ThreadPoolBuilder};

/// Ceiling applied when nothing else specifies a thread count. See the module docs: this
/// exists so behaviour does not scale with the host's core count by default.
const DEFAULT_MAX_THREADS: usize = 16;

/// 16 MiB per worker, against rayon's 2 MiB default.
///
/// Address space, not memory: thread stacks are mapped lazily, so the resident cost is the
/// pages actually touched (peak RSS was flat at ~400 MB across every thread count measured).
/// With the `bgzf.rs` frame fix the deep-recursion headroom should no longer be needed; it is
/// held in reserve because ISSUE-001 was schedule-dependent and therefore cheap to believe
/// fixed on insufficient evidence.
const STACK_SIZE: usize = 16 << 20;

/// Overrides `STACK_SIZE`, in bytes, via `PGP_SCAN_STACK`.
///
/// This exists to make ISSUE-001 *falsifiable*. "It did not crash in 50 runs" is weak
/// evidence for a schedule-dependent bug that already survived several clean runs; shrinking
/// the stack to a value the old frame cannot fit and the new one easily can turns the same
/// question into a deterministic experiment. `tests/test_scan_stack_depth.py` uses it that
/// way. It is also a genuine escape hatch if some future host needs more headroom.
fn stack_size() -> usize {
    env_usize("PGP_SCAN_STACK").unwrap_or(STACK_SIZE)
}

static POOL: OnceLock<ThreadPool> = OnceLock::new();

fn env_usize(var: &str) -> Option<usize> {
    std::env::var(var).ok()?.trim().parse::<usize>().ok().filter(|n| *n > 0)
}

fn resolve(threads: Option<usize>) -> usize {
    if let Some(n) = threads.filter(|n| *n > 0) {
        return n;
    }
    if let Some(n) = env_usize("PGP_SCAN_THREADS").or_else(|| env_usize("RAYON_NUM_THREADS")) {
        return n;
    }
    let cores = std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1);
    cores.min(DEFAULT_MAX_THREADS)
}

/// The process-wide scan pool, built on first use.
///
/// **Built once.** A later call asking for a different `threads` gets the existing pool; the
/// count is not renegotiated mid-process. Every caller in this codebase passes the same
/// configured value for the whole run, so this is not a live concern -- but it is silent if
/// it ever becomes one, which is why `threads()` exposes what was actually built so a caller
/// can compare rather than assume.
pub fn pool(threads: Option<usize>) -> &'static ThreadPool {
    POOL.get_or_init(|| {
        ThreadPoolBuilder::new()
            .num_threads(resolve(threads))
            .stack_size(stack_size())
            .thread_name(|i| format!("pgp-scan-{i}"))
            .build()
            .expect("failed to build the pgp-scan rayon pool")
    })
}

/// Worker count of the live pool, or the count that *would* be chosen if it is not built yet.
///
/// Exposed to Python so a run can report the number it is actually using, and so the
/// many-core verification in ISSUE-002 (total threads should track the core count, not
/// N x cores) can be checked without `ps`.
pub fn threads(requested: Option<usize>) -> usize {
    match POOL.get() {
        Some(p) => p.current_num_threads(),
        None => resolve(requested),
    }
}
