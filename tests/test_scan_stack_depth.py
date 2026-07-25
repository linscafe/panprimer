"""ISSUE-001 regression: the inflate path must not put a large frame on a recursive stack.

The bug was an intermittent `SIGSEGV` in `pangenome_primer._scan` -- a rayon worker
overflowing its 2 MiB stack. It was never reproducible on demand (the same input crashed
once and then completed a dozen times), because the depth reached depends on rayon's
work-stealing schedule rather than on the data.

That is exactly the kind of bug a "ran it 50 times, no crash" test certifies falsely. So this
file does not try to reproduce the crash at the shipped stack size. It **shrinks the stack**
until the difference between the old and new code is deterministic, and asserts on that:

    build                                     stack     result
    inlined `Decompress::new` (the bug)        256 KB    SIGSEGV, 3/3 runs
    inlined `Decompress::new`                  2 MiB     completes -- the intermittency
    `#[inline(never)] new_decompressor`        256 KB    completes
    `#[inline(never)] new_decompressor`        128 KB    completes

The mechanism, established by disassembling the `.so` that actually crashed:
`Decompress::new` bottoms out in `Box::default()`, which materialises a ~44 KB `InflateState`
as a stack temporary. Inlined into the closure, it landed in the frame of rayon's
*recursive* `bridge_producer_consumer::helper` -- measured at 0xd298 = 53,912 bytes, against
2 MiB of stack. Out of line, that helper's frame is 440 bytes.

`STACK_FLOOR` below is therefore a proxy for frame size, which is the thing that regressed.
Anything that reintroduces a multi-KB frame into the parallel inflate closure fails here.
"""
from __future__ import annotations

import os
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

from pangenome_primer import rust_backend

# Small enough that the old, buggy build died on it every time; large enough that the fixed
# build clears it with room to spare. Not a tuning knob -- see the module docstring.
STACK_FLOOR = 256 * 1024

# Enough BGZF members in one read batch for rayon to split several levels deep. A handful of
# members would make the test pass on the buggy build too, since the recursion would never
# get deep enough to matter.
FIXTURE_MB = 6

pytestmark = pytest.mark.skipif(
    not rust_backend.available(),
    reason="needs the compiled pangenome_primer._scan extension",
)


@pytest.fixture(scope="module")
def big_bgzf(tmp_path_factory) -> str:
    """A multi-MB BGZF FASTA -- ~100 members, so the inflate iterator really recurses."""
    import random

    d = tmp_path_factory.mktemp("stack")
    plain = d / "big.fa"
    rng = random.Random(11)
    n_lines = FIXTURE_MB * 1_000_000 // 60
    # A 20 nt marker on every 2000th line, so markers are spread across every BGZF member
    # rather than clustered in the first one. `test_every_bgzf_member_inflates_correctly`
    # depends on that spread: it is what makes a stale-decompressor bug visible.
    marker = "GGATCCTTAAGGCCTTAAGC"
    with open(plain, "w") as fh:
        fh.write(">ctg1\n")
        for i in range(n_lines):
            line = "".join(rng.choice("ACGT") for _ in range(60))
            if i % 2000 == 1000:
                line = marker + line[len(marker):]
            fh.write(line + "\n")
    gz = d / "big.fa.gz"
    with open(gz, "wb") as out:
        subprocess.run(["bgzip", "-f", "-c", str(plain)], stdout=out, check=True)
    subprocess.run(["samtools", "faidx", str(gz)], check=True)
    plain.unlink()
    return str(gz)


def _scan_in_subprocess(gz: str, stack: int) -> subprocess.CompletedProcess:
    """Run one scan in a fresh process. Subprocess, not in-process: a stack overflow is a
    SIGSEGV, so the failure has to be observed as an exit status rather than an exception.
    A fresh process is also required because the pool is built once and `PGP_SCAN_STACK` is
    only read at that moment."""
    code = textwrap.dedent(
        f"""
        from pangenome_primer import rust_backend
        got = rust_backend.find_binding_sites_batch(
            ["ACGTACGTACGTACGTACGT", "TTGCACAGTCCAGATTGCAA"], {gz!r}, "H#hap1", 2,
            threads=4)
        print("OK", sum(len(v) for v in got.values()))
        """
    )
    env = {**os.environ, "PGP_SCAN_STACK": str(stack)}
    return subprocess.run([sys.executable, "-c", code], capture_output=True, text=True,
                          env=env, cwd=Path(__file__).resolve().parents[1], timeout=300)


class TestInflateFrameStaysOffTheRecursiveStack:
    def test_scan_survives_a_small_worker_stack(self, big_bgzf):
        """A 256 KB worker stack killed the buggy build on every run. If this segfaults,
        something has put a large frame back inside the parallel inflate closure -- check for
        a value-returning constructor inlined into the `map_init` closure in `bgzf.rs`."""
        r = _scan_in_subprocess(big_bgzf, STACK_FLOOR)
        assert r.returncode == 0, (
            f"scan died with returncode {r.returncode} "
            f"({'SIGSEGV -- stack overflow, see ISSUE-001' if r.returncode == -11 else 'error'})"
            f" at a {STACK_FLOOR // 1024} KB worker stack.\nstderr: {r.stderr[-2000:]}"
        )
        assert r.stdout.startswith("OK"), r.stdout + r.stderr

    def test_every_bgzf_member_inflates_correctly(self, big_bgzf):
        """The frame fix reuses one `Decompress` per rayon job via `reset()` instead of
        constructing one per member, so every member after the first in a job now depends on
        that reset being right.

        Deleting the reset was tried, to check the risk is real: the next member in each job
        inflates wrong (65,276 bytes against an ISIZE of 65,270), and `inflate_raw`'s
        per-member ISIZE check turns that into a hard error. So the loud failure mode is
        already covered. What ISIZE cannot see is wrong *content* of the right length, which
        is what this test pins: a marker planted at known offsets spanning the whole file
        (~100 BGZF members), with `str.find` over the same sequence as the reference.
        """
        import pysam

        fa = pysam.FastaFile(big_bgzf)
        seq = fa.fetch(fa.references[0]).upper()
        marker = "GGATCCTTAAGGCCTTAAGC"

        expect, at = [], seq.find(marker)
        while at != -1:
            expect.append(at)
            at = seq.find(marker, at + 1)
        assert len(expect) >= 20, (
            f"fixture planted only {len(expect)} markers; the test cannot show that later "
            "members inflate correctly"
        )

        got = rust_backend.find_binding_sites_batch([marker], big_bgzf, "H#hap1", 0)
        starts = sorted(s.start for s in got.get(marker, []) if s.strand.value == "+")
        assert starts == expect, (
            f"scanner found {len(starts)} exact sites, sequence contains {len(expect)}. "
            "A shortfall concentrated in later offsets means BGZF members are being "
            "decompressed with a stale decompressor state."
        )

    def test_results_are_unchanged_by_the_stack_size(self, big_bgzf):
        """Stack size must be a performance/robustness knob only. If the two runs disagree,
        the scan's output depends on how rayon split the work -- which would mean a data race
        rather than the stack-depth problem this file is about, and a much worse bug."""
        small = _scan_in_subprocess(big_bgzf, STACK_FLOOR)
        default = _scan_in_subprocess(big_bgzf, 16 << 20)
        assert small.returncode == 0 and default.returncode == 0
        assert small.stdout == default.stdout, (
            "same scan, different worker stack size, different hits -- the inflate output "
            f"is not deterministic.\nsmall: {small.stdout!r}\ndefault: {default.stdout!r}"
        )


def _in_fresh_process(body: str, **env_extra) -> str:
    """Run `body` in a new interpreter and return its stdout.

    Every assertion here has to be made in a fresh process. The pool is a process-wide
    `OnceLock`, so `pool_threads(3)` answers 3 only while no scan has run yet -- asserting it
    in-process makes the result depend on which tests ran first, which is how the first draft
    of this file passed alone and failed in the full suite.
    """
    env = {**os.environ, **{k: str(v) for k, v in env_extra.items()}}
    for k in ("PGP_SCAN_THREADS", "RAYON_NUM_THREADS"):
        if k not in env_extra:
            env.pop(k, None)
    r = subprocess.run([sys.executable, "-c", textwrap.dedent(body)], capture_output=True,
                       text=True, env=env, timeout=120,
                       cwd=Path(__file__).resolve().parents[1])
    assert r.returncode == 0, r.stdout + r.stderr
    return r.stdout.strip()


class TestPoolIsExplicitlySized:
    """ISSUE-002: the crate must not inherit rayon's one-worker-per-CPU global default."""

    def test_configured_thread_count_is_what_gets_built(self):
        out = _in_fresh_process(
            "from pangenome_primer import rust_backend\n"
            "print(rust_backend.pool_threads(3))\n"
        )
        assert out == "3"

    def test_zero_and_none_mean_auto_not_zero_threads(self):
        """A pool of 0 workers would deadlock, so 0 must fold into "auto" rather than reach
        rayon as a literal count."""
        out = _in_fresh_process(
            "from pangenome_primer import rust_backend\n"
            "print(rust_backend.pool_threads(None), rust_backend.pool_threads(0))\n"
        )
        auto, zero = out.split()
        assert int(auto) >= 1 and zero == auto

    def test_default_does_not_scale_with_host_core_count(self):
        """The point of the cap: a 256-core host must not spawn 256 workers for a job whose
        parallel efficiency is already 21% at 16 threads."""
        out = _in_fresh_process(
            "from pangenome_primer import rust_backend\n"
            "print(rust_backend.pool_threads(None))\n"
        )
        assert int(out) <= 16

    def test_rayon_num_threads_is_still_honoured(self):
        """Anyone tuning a rayon program reaches for this variable; it must not be silently
        ignored now that the crate builds its own pool instead of using the global one."""
        out = _in_fresh_process(
            "from pangenome_primer import rust_backend\n"
            "print(rust_backend.pool_threads())\n",
            RAYON_NUM_THREADS=2,
        )
        assert out == "2"

    def test_pgp_scan_threads_beats_rayon_num_threads(self):
        out = _in_fresh_process(
            "from pangenome_primer import rust_backend\n"
            "print(rust_backend.pool_threads())\n",
            PGP_SCAN_THREADS=5, RAYON_NUM_THREADS=2,
        )
        assert out == "5"

    def test_single_threaded_pool_does_not_deadlock(self, big_bgzf):
        """`threads: 1` is a legal setting a user reaching for "use one core" will try, and
        the scan now runs *inside* `pool.install()` -- so the serial streaming work and both
        parallel regions share that one worker. Rayon handles it (a `join` runs both halves
        on the calling worker), but a one-worker pool is the shape where a nested-parallelism
        mistake turns into a hang rather than a slowdown. The subprocess timeout is the
        assertion."""
        out = _in_fresh_process(
            f"""
            from pangenome_primer import rust_backend
            got = rust_backend.find_binding_sites_batch(
                ["ACGTACGTACGTACGTACGT"], {big_bgzf!r}, "H", 1, threads=1)
            print(rust_backend.pool_threads(), sum(len(v) for v in got.values()) >= 0)
            """
        )
        assert out.split()[0] == "1", out

    def test_a_second_different_thread_count_warns_instead_of_silently_losing(self, big_bgzf):
        """The pool is built once. A later `search.threads` is ignored -- which is tolerable,
        but only if it is audible: on a 64-way fan-out, a silently ignored cap is the
        thrashing ISSUE-002 is about."""
        out = _in_fresh_process(
            f"""
            import warnings
            from pangenome_primer import rust_backend
            gz = {big_bgzf!r}
            rust_backend.find_binding_sites_batch(["ACGTACGTACGTACGTACGT"], gz, "H", 1,
                                                  threads=2)
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                rust_backend.find_binding_sites_batch(["ACGTACGTACGTACGTACGT"], gz, "H", 1,
                                                      threads=8)
            print([str(x.message) for x in w if x.category is RuntimeWarning][:1])
            """
        )
        assert "threads=8 was ignored" in out and "built with 2 workers" in out, out
