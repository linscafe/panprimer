# 30-haplotype **VERIFY** run — measured process times

Host: 16 logical cores, 15 GB RAM, WSL2. Disk /dev/sdd, 695 GB free at start.
Manifest: config/samples_30.tsv (30 haplotypes, 6 per superpopulation x AFR/AMR/EAS/EUR/SAS).
Starting state: 10 of 30 assemblies already downloaded, 3 of 30 anchor grids already built.

| step | scope | start | end | wall |
|---|---|---|---|---|
| download | 20 new assemblies (~18 GB); 10 already present, md5-reverified | 2026-07-26T00:00:14 | 2026-07-26T00:30:29 | **30m15s** (1815 s, ~10 MB/s) |
| anchor grids | 27 builds, sequential, MM_THREADS=4 (3 already present, skipped) | 2026-07-26T00:31:00 | 2026-07-26T02:39:29 | **2h08m29s** (7709 s) |
| verify run | 4 primer pairs x 30 haplotypes, rust backend, grid projection | 2026-07-26T02:40:34 | 2026-07-26T02:47:13 | **6m39s** (399 s) |

**Total of the three phases: 2h45m23s** (wall clock 00:00:14 -> 02:47:13 = 2h46m59s, the
difference being the gaps where I inspected results between phases). Anchor-grid builds are
2h08m of it — **78%**. Everything else is comparatively free.

Note this run started with 10 of 30 assemblies already downloaded and 3 of 30 grids already
built. A true cold start would be ~45 min of download and ~2h15m of grid builds, ~3h05m total.

## Per-anchor-grid build (27 builds, derived from log timestamps, not the script's
## truncated integer-minute figures)

| min | median | mean | max | sum |
|---|---|---|---|---|
| 3.93 min | 4.48 min | 4.52 min | 5.40 min | 2.04 h |

Grid quality was uniform: 90.6-92.8% of probes anchored, ~4.3 MB per grid, 0 failures.

## **VERIFY** run resources

- wall 399 s for 4 pairs x 30 haplotypes = **13.3 s per haplotype** (all 4 pairs)
- 381% CPU (the scan parallelises; projection and pairing do not)
- peak RSS **3.2 GB** -- higher than the ~1 GB README quotes for search+projection,
  because verify holds the CHM13 handle and the per-haplotype anchor grid alongside
  the scanner. Worth re-checking before sizing a many-concurrent-verify job.
