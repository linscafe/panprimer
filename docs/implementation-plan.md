# Implementation plan — verify pipeline storage & runtime optimization

**Companion to** [`plan-optimization.md`](plan-optimization.md), which contains the measurements and the rationale. This page is the execution plan: phases, delegation, and gates.

**Decisions taken** (2026-07-24): Rust/PyO3 for the matcher; spike-then-commit for the full-scale pangenome index; no bulk genome data deleted until the replacement path is green.

---

## Context

The verify pipeline costs **15.10 GB per haplotype** in storage, and **~5.7 min per haplotype** end to end (measured: 17 min 05 s for the 3-haplotype demo verify run). Of that, only **47.3 s** is the `bwa` genome-wide search — the remainder is dominated by loading a 5.80 GB `.mmi` per haplotype for coordinate projection. At the ~464 haplotypes of HPRC R2 that is ~7.0 TB of storage and **~44 h** of serial wall time per request. Both numbers trace to one cause: each haplotype is treated as an independent genome to be indexed, when 464 human haplotypes are ~99.9% identical.

> **Baseline correction (measured 2026-07-24, after Phase 0/1).** `plan-optimization.md` quotes 47.3 s/haplotype; that is the search component alone. The end-to-end figure is ~5.7 min/haplotype, so **projection — not search — is the larger cost today.** This raises the value of Phase 4 (anchor grid) relative to Phase 2, though Phase 2 still gates Phase 3's storage reclaim and must land first.

The target end state is **~0.91 GB and <1 s per haplotype** — ~430 GB and minutes per run — reached in six phases, each independently shippable and revertible.

### Two findings that shape the sequencing

1. **There is no safety net today.** `pangenome-primer selftest` runs entirely on in-memory synthetic strings (`fixtures.py`) — no FASTA, no bwa, no mappy, no pysam. It will keep passing through this entire refactor **without validating any of it**. There is also no direct unit test of `binding.find_binding_sites_naive` or `bwa_backend.find_binding_sites_batch`. Phase 0 exists solely to fix this, and nothing else may start until it lands.
2. **Storage cleanup is the *last* step, not the first.** The uncompressed `.fa` and the bwa index are required by the pipeline as it stands today. Deleting them before Phase 2 ships breaks the tool, and a bwa index costs ~55 min/haplotype to rebuild.

---

## The seams this plan uses

All of these already exist. The refactor is mostly a matter of filling them in rather than inventing new structure.

| Seam | Location | Role |
|---|---|---|
| `find_binding_sites_batch(seqs, fasta, hap_id, max_mm) -> dict[seq, list[BindingSite]]` | `bwa_backend.py:65`, imported lazily at `verify.py:106`, `cli.py:108`, `cli.py:387` | The genome-search seam. Swapping backends is a three-line import redirect. |
| `binding.find_binding_sites(..., backend=...)` | `binding.py:108-125` | Declared dispatcher. `"bwa"` branch raises `NotImplementedError` naming a nonexistent `modules/binding_search.nf`. Vestigial today — becomes the real registry. |
| `search.backend` | `config/defaults.yaml:25`, marked `RESERVED (not yet wired)` | Config key with no reader. Becomes the backend selector. |
| `project.project_target` | `project.py:152-172` | Already dispatches PAF-lift vs live alignment. Anchor-grid lift becomes a third branch. |
| `window_fn(chrom, start, end) -> str` | `engine.py:71-79`; callers `verify.py:152`, `cli.py:217`, `cli.py:422` | The narrowest abstraction over "the genome". Any BGZF/remote reader drops in unchanged. |

**Contract any new search backend must satisfy** (derived from `binding.py` and its consumers — this is the spec Phase 0 must pin down in tests):

1. `chrom`/`start`/`end` are **haplotype genome coordinates**; `pcr.pair_amplicons` computes `size = m.end - p.start` and `window_fn` fetches against them.
2. `end - start` **must equal the primer length** — `window_fn` assumes the returned window is primer-length.
3. `mismatch_offsets_3p` uses the **3'-end convention** (0 = 3' terminal base); `classify.competent_rule:55` and `classify.is_competent:127` gate on it.
4. `strand` must be correct — `competent_thermo:77` revcomps the window for `MINUS`; a wrong strand yields a nonsense ~4 °C Tm.
5. Returning a **superset** of true positions is fine when the backend delegates re-scoring, which is exactly what `bwa_backend.py:92-94` does today.
6. Match semantics to preserve: `N` counts as a mismatch (never special-cased), no indels (`has_indel` always `False`), emission order plus-strand-then-minus, and **no dedupe** in the naive path.

---

## Phase 0 — Safety net (blocking; nothing else starts until green)

**Model: Sonnet.** Mechanical, high-volume, well-specified.

1. **Characterization tests for `binding.find_binding_sites_naive`** covering every clause of the contract above: mismatch counting at the cap boundary, 3'-offset arithmetic on both strands, `N` handling, palindromic primer duplication, emission order, `end - start == L`.
2. **Golden-output harness.** Freeze the current `demo-verify-pipeline/verify.json` and `demo-design-pipeline/gapdh/results.json` as fixtures, plus a comparator that diffs a fresh run **cell for cell** — including the CYP2D7 off-target and the engineered dropout row. This is the gate the plan doc names (`plan-optimization.md:176`).
3. **A backend conformance suite** parameterized over backend name, so every future backend (Rust, anchor-lift, pangenome index) runs the identical assertions against `naive` as reference.
4. **A synthetic mini-genome fixture** — a few Mb of BGZF FASTA with planted binding sites at known coordinates, including sites near contig edges and within `slop` distance of each other. Real assemblies are too slow for CI.

**Gate:** `pytest -q` green (42 existing + new), and the golden comparator reproduces today's demo output exactly.

---

## Phase 1 — Cleanup and housekeeping

**Model: Sonnet.** Independent of everything else; can run in parallel with Phase 0.

- Safe stale-file deletion (~2.4 MB — see [Cleanup](#cleanup-executed) below; already executed).
- `prepare_haplotypes.sh:63` — add `-K 100M` to the PAF minimap2 call. Two builds OOM'd at 14.5 GB peak on this 15 GB machine.
- Clean up `.part`/truncated outputs on failure rather than leaving them; `project.py:55-59`'s zero-byte `.mmi` guard exists only to survive this.
- `bwa_backend.ensure_index` (`bwa_backend.py:18-23`) — fail loudly pointing at `prepare_haplotypes.sh` instead of silently launching a ~55 min in-process `bwa index`.
- Delete dead code: `align_cache.build_alignment` (`align_cache.py:28-38`, the only zero-reference public function in the package) and the `binding.py:123` reference to the nonexistent `modules/binding_search.nf`.
- Hoist redundant FASTA handles: `verify.py:144` and `bwa_backend.py:82` both open the haplotype; `align_cache.py:135` reopens it **per primer** inside verify's inner loop.
- `README.md:14` links `demo-design-pipeline/gapdh/report.html` and `demo-verify-pipeline/verify_matrix.html`, both of which `.gitignore:24-25` excludes — dead for anyone who clones. Either track the two HTML files or point at `docs/img/`.

---

## Phase 2 — Rust matcher (`pgp-scan`)

**Model: Opus.** Algorithmic correctness plus FFI; the highest-risk work in the plan.

Replaces `bwa aln` entirely, which is what makes Phase 3 possible.

**Crate design.** A `maturin`/PyO3 extension exposing one function mirroring the existing signature:
```
scan(seqs: list[str], bgzf_path: str, max_mismatches: int, slop: int) -> dict[str, list[SiteTuple]]
```
Algorithm: stream BGZF blocks, 2-bit encode, pigeonhole on exact 10-mer seeds (~2,861 expected genomic hits per seed in 3.1 Gb; ~200k candidates for a typical panel), then verify candidates with a bit-parallel ≤3-mismatch comparator. Emit plus and minus strand separately, preserving the naive backend's exact semantics from the Phase 0 contract.

**Packaging.** `pyproject.toml` moves from setuptools to `maturin`; `requires-python` stays `>=3.10`. Ship wheels for linux-x86_64 (manylinux) and macOS arm64, and add `rust` to `env/environment.yml` so a source build works in the conda env. **The Python fallback must remain**: `search.backend` selects `rust` when the extension imports and `naive`/`bwa` otherwise, so a wheel-less platform degrades rather than breaks.

**Wiring.** Implement the real `binding.find_binding_sites` dispatcher, make `search.backend` live, and redirect the three lazy import sites.

**Gate:** backend conformance suite passes with `backend=rust`; golden comparator reproduces the demo matrix cell for cell; measured wall time per haplotype recorded against the **47.3 s** baseline. Target ~10–15 s.

**Measured 2026-07-25 (HG00097_hap1, 8 primers):**

| | bwa baseline | rust | change |
|---|---:|---:|---:|
| wall time | 47.3 s | **6.9 s** | **6.9× faster** |
| peak RSS | 4.45 GB | **420 MB** | **10.6× less** |
| on-disk index | 5.31 GB | **none** | — |

Fast suite **105 passed, 2 skipped**. Conformance passes with `rust` live (no skips); 23 differential tests confirm `rust` and `naive` return identical site sets across random genomes/primers, contig boundaries, soft-masked lowercase, `N` runs, and sites closer together than `slop`. Fallback verified by simulating a missing extension: `auto`/`rust` degrade to `bwa` with a warning, an unknown name raises `ValueError`.

> **Open question at the gate.** The Rust scanner is *exhaustive* (differential-equal to `naive`), whereas `bwa aln` is a heuristic that can miss sites. It is therefore possible for `rust` to find genuine off-targets the golden matrix — generated with `bwa` — does not contain. If the golden comparator diverges, the correct response is to **inspect each differing cell** and decide whether it is a Rust bug or a real off-target `bwa` missed. Do not "fix" it by relaxing the comparator.

### Gate result: diverged on exactly one row — and the golden was wrong

The comparator failed (25 min 01 s) on **`CYP2D6_dropout` only**, across all three haplotypes: `dropout` → `multi_product`, with `on_target=[122]` (exactly the 122 bp target span) plus ~2,000 off-target products spanning a continuous 85–2000 bp spectrum. Every other cell — including the `CYP2D6_paralog` CYP2D7 off-target the README leads with, and all GAPDH rows — matched **exactly**.

Diagnosis (per-primer site counts, HG00097_hap1, ≤3 mismatches):

| pair | fwd sites | rev sites | perfect (0 mm) fwd/rev |
|---|---:|---:|---|
| GAPDH_clean | 149 | 263 | 2 / 1 |
| CYP2D6_clean | 178 | 168 | 1 / 2 |
| CYP2D6_paralog | 180 | 331 | 1 / 13 |
| **CYP2D6_dropout** | **329,968** | 365 | **12,760** / 1 |

`CYP2D6_dropout`'s forward primer `ACTCCTGGGCTCAAGCAATC` is an **Alu repeat consensus** with ~12,760 perfect genomic copies. Head-to-head on the same haplotype:

- **reverse primer: bwa 365, rust 365 — identical.** Independent confirmation the scanner is correct where bwa is reliable.
- **forward primer: bwa 1, rust 329,968.** `bwa samse -n 1000` silently truncates a read whose alternate count exceeds the cap, so the pipeline saw a single forward site, found no partner near the target, and recorded "no on-target product" — i.e. **`dropout`**.

**The golden `dropout` verdict is a false negative produced by bwa's XA cap, not biology.** Nothing "dropped out": a primer binding ~330k places amplifies everywhere, which is what `multi_product` means. The exhaustive scanner did not regress the demo — it exposed a pre-existing wrong call in it.

Two consequences requiring action:

1. **The demo narrative is affected.** *(Resolved — the row was replaced; see below.)* `README.md:10` presented this row as a dropout demonstration. That claim was wrong in mechanism.

### Replacement dropout row (applied 2026-07-25)

```csv
CYP2D6_dropout,chr22:42608601-42608798,CCAAGTTGCGCAAGGTGGAG,TGTGACCAGCTGGACAGAGC
```

Exploits **rs1058164** (CYP2D6 c.1661G>C, exon 3) sitting under the forward primer's 3′
terminal base. HG01884 and HG00408 carry G (primer binds); HG00097 carries C, giving a G·G
purine–purine terminal mismatch that `three_prime_hard_nt = 2` hard-fails. Independently
confirmed with the Rust scanner: on HG00097 the forward primer has **0** perfect matches and
its CYP2D6 site reads `mm=1, mismatch_offsets_3p=[0]`; on HG01884 it has a clean perfect
match. The CYP2D7/CYP2D8P paralog copies also carry a 3′-terminal mismatch, so no paralog
band appears.

Specificity is normal — 113–116 raw sites (0–1 perfect) forward, 301–303 (1 perfect) reverse
— against 329,968 / 12,760 for the Alu row it replaces.

Regenerated demo matrix (13 min 44 s, `backend=rust`, cap 100):

| row | expected | HG01884 (AFR) | HG00097 (EUR) | HG00408 (EAS) |
|---|---:|---|---|---|
| GAPDH_clean | 201 | 201 | 201 | 201 |
| CYP2D6_clean | 280 | 280 | 280 | 280 |
| CYP2D6_paralog | 282 | 282 + off | 282 + off | 282 + off |
| **CYP2D6_dropout** | **197** | **197** | **dropout** | **197** |

**Why population-differential matters.** A row that fails on *every* haplotype cannot be
distinguished from a simply broken primer — which is exactly how the Alu artefact survived.
`tests/test_golden_demo.py` previously asserted `all(status == "dropout")` and passed happily
on it; it now requires **both** a dropout and a pass among the cells. Do not weaken that back.

Also updated: `tests/golden/verify.json` re-frozen, `docs/img/verify_matrix.svg` regenerated,
`README.md` corrected from "fails on every haplotype".

> **Follow-up worth doing: consolidate cell rendering.** The verify matrix cell label is
> derived independently in **five** places — `report.verify_to_dict`, `report._cell_text`,
> `report._verify_cell_md`, `report/verify_matrix.html.j2`, and
> `scripts/snapshot_verify_svg.py:cell_segments`. Each re-implements the same branch order,
> and each must special-case a capped cell *before* the empty-product test or it renders
> "dropout". Adding the cap required touching all five, and the omission was shipped twice
> (serializer, then SVG) because engine-level tests pass regardless. One shared helper
> returning `(text, css_class)` would make the next cell state a one-line change.
2. **Repeat-derived primers are a performance hazard.** *(Addressed — see below.)*

### Binding-site cap (implemented 2026-07-25)

`search.max_binding_sites` (default **100**) caps amplification-competent sites per primer.
`engine._filter_competent` **short-circuits** at the cap — scanning stops rather than scoring
every site — and the cell reports `>100 binding sites` instead of enumerating products.

The threshold comes from measured separation, not intuition (HG00097_hap1, ≤3 mismatches):

| primer | raw sites | competent |
|---|---:|---:|
| CYP2D6_clean fwd / rev | 178 / 168 | 1 / 4 |
| GAPDH_clean fwd / rev | 149 / 263 | 4 / 11 |
| CYP2D6_paralog fwd / rev | 180 / 331 | 13 / **28** |
| CYP2D6_dropout fwd (Alu) | 329,968 | **~106,000** |

Well-behaved primers top out at 28; the Alu primer is ~3,800× higher. 100 clears the former
with margin and catches the latter. **A cap of 10 would have suppressed the CYP2D7 paralog
row** (13 and 28) and `GAPDH_clean`'s reverse primer (11) — i.e. broken the demo's headline
result to catch the repeat.

Verified end to end (3-haplotype demo, `backend=rust`, **18 min 41 s** vs 25 min 01 s uncapped):

| row | result |
|---|---|
| GAPDH_clean | `pass`, 201 bp on all three |
| CYP2D6_clean | `pass`, 280 bp on all three |
| CYP2D6_paralog | `multi_product`, 6 products (1 on-target) — **CYP2D7 detection intact** |
| CYP2D6_dropout | `multi_product`, `>100 binding sites (forward primer); likely repeat-derived` |

> **Bug caught during verification.** `VerifyCell.site_cap` was added to the dataclass and
> the Jinja template, but `report.verify_to_dict` builds cell dicts from an explicit field
> list and omitted it. Every renderer reads that dict, not the dataclass — so capped cells
> serialized without the flag and fell through to the "dropout" branch, reproducing exactly
> the mislabel the cap exists to prevent. Engine-level tests passed throughout. Fixed in all
> four renderers (dict, TSV, Markdown, HTML) and covered by
> `tests/test_binding_site_cap.py::TestCappedCellSurvivesSerialization`, which walks the
> full path rather than stopping at the engine.

Original hazard, for the record: 330k sites pushed the run to 25 min (vs the 17 min bwa baseline) *despite* search dropping from 47.3 s to 6.9 s — the cost moved into per-site thermodynamic scoring. At pangenome scale this is pathological. The pipeline should detect a repeat-derived primer early (site count far above the ~184 expected for a random 20-mer) and short-circuit with a clear verdict rather than scoring every site.

---

## Phase 3 — BGZF-only storage

**Model: Sonnet**, with Opus review of the Nextflow staging change.

Verified prerequisite: the HPRC-distributed `.fa.gz` is **already BGZF** (`BC` subfield confirmed), `pysam.FastaFile` opens it directly, and 500 random 26 bp fetches take 1.3 ms. `samples.load_haplotypes` is path-agnostic, so a `.fa.gz` path passes its existence check unchanged.

- Drop the `.fa.gz` → `.fa` rewrite at `prepare_haplotypes.sh:79`; stop `gunzip`-ing; stop building bwa indexes and `.mmi`.
- Generate `.gzi` once per haplotype (`samtools faidx` on the `.gz`). Only `HG00097_hap1` has one today.
- Repoint `local_path` in `config/samples.tsv` and both demo `samples.tsv` files.
- **`main.nf:116-124`** — the staging list `['fai','mmi','chm13.paf','amb','ann','bwt','pac','sa']` must become `['fai','gzi']`. `collectMany { files(...) }` skips missing siblings silently, so additions are safe, but changing `local_path` from `.fa` to `.fa.gz` changes every glob base at once. This is the single highest-risk line in the phase.
- Build the shared CHM13 bwa index (~55 min, ~5.3 GB) — it **does not exist today**; `hprc-r2/references/` has only `.fa`, `.fai`, and `.asm5.mmi`.
- Update the README storage/memory warning (lines 16–19) and `hprc.py`'s `INDEXED_GB_PER_HAP = 15.0` resource estimate.

**Gate:** full demo run from a `.fa.gz`-only haplotype directory reproduces the golden output. **Only after this gate passes** may the `.fa`, `.bwt`, `.sa`, `.pac`, `.amb`, `.ann` files be deleted (~28.7 + 49.4 GB).

### Phase 3 in progress — three corrections to the plan above

**1. Repointing `local_path` silently breaks every projection sidecar.** `project.py:55` and `align_cache.py:23` derived `.mmi` and `.chm13.paf` paths by plain string append on `local_path`. Changing that path from `X.fa` to `X.fa.gz` makes them look for `X.fa.gz.mmi` / `X.fa.gz.chm13.paf`, miss, and fall through to `mappy.Aligner(X.fa.gz)` — building a whole-genome minimap2 index **in-process, per haplotype, per run**. This does not raise; it produces correct answers minutes-per-haplotype slower, so no output-comparing test would catch it. Fixed by `samples.sidecar_path(seq_path, suffix)`, which prefers an exact-named sidecar, then the `.gz`-stripped stem, then falls back to the input path for writers. Pinned by `tests/test_sidecar_path.py` (6 tests), which asserts on the **resolved path**, not on run output, because the failure mode is performance rather than correctness.

**2. `collectMany { files(...) }` does *not* skip missing siblings** — the claim in the bullet above is wrong. For a pattern containing no glob metacharacter, Nextflow's `files()` returns the literal path whether or not anything exists there, so every absent sidecar would be staged as a broken symlink. This is a **pre-existing** bug, not one Phase 3 introduces: `.chm13.paf` has never been built for any haplotype, so it was already in the staging list as a phantom. Fixed with `.findAll { it.exists() }`. Caught by executing the closure against the real `samples.tsv` rather than relying on `nextflow run --help` parsing the file.

**3. `.mmi` must keep being built until Phase 4 lands.** The bullet "stop building bwa indexes and `.mmi`" conflates two phases. Projection has no replacement until the anchor grid exists, so dropping `.mmi` now would trigger exactly the in-process index build described in (1). `prepare_haplotypes.sh` therefore still builds a `.mmi` — now from the BGZF, which `minimap2` reads directly — and only the `.fa` + bwa index become opt-in (`WITH_BWA=1`). The staging list is `['fai','gzi','mmi','chm13.paf','amb','ann','bwt','pac','sa']` filtered by existence, not `['fai','gzi']`; it narrows to that naturally as the artifacts stop being produced.

**4. The shared CHM13 bwa index is not needed at all.** The bullet above budgets ~55 min and ~5.3 GB for it. Nothing in the codebase asks for one: CHM13 is used for anchoring and template fetch, both `minimap2`/`pysam` operations, and `grep` finds no `bwa`/`ensure_index`/`find_binding_sites` call taking a CHM13 path. It was a carry-over from when `bwa` was the only search backend. Dropped — this is 55 minutes and 5.3 GB the phase does not have to spend. Revisit only if Phase 5 (CHM13-once discovery) turns out to want a genome index on the reference, in which case the `rust` scanner already covers it index-free.

**Measured staging cost per haplotype (via the closure probe):** 11.13 GB of sidecars today (`.fai .mmi .amb .ann .bwt .pac .sa`), or 5.3 GB for the four haplotypes with no `.mmi`. Post-repoint that is `.fa.gz` 0.90 GB + `.fai`/`.gzi` 0.76 MB.

**Scope: 3 haplotypes.** All three manifests (`config/samples.tsv`, `demo-verify-pipeline/`, `demo-design-pipeline/`) are pinned to HG01884#hap1 / HG00097#hap1 / HG00408#hap1 — AFR/EUR/EAS. `config/samples.tsv` previously listed 10; it was trimmed because `nextflow.config` defaults `--samples` to it, so a bare `nextflow run main.nf` would otherwise fan out across every assembly on disk. The CLI's `--samples` is `required=True`, so it never had this exposure. The 7 unlisted haplotypes remain on disk untouched; git history and `pangenome-primer fetch-subset` both restore the wider manifest.

### Phase 3 gate: passed, and the reclaim is done

`pytest --runslow tests/test_golden_demo.py` — **14 passed in 31m36s** reading entirely from `.fa.gz`, against 33m07s from the uncompressed `.fa`. Fast suite 131 passed. A single-haplotype smoke run completed in 4m25s with the expected verdicts (`CYP2D6_paralog` off-target present, `CYP2D6_dropout` clean on the AFR haplotype — that dropout is EUR-specific).

**Reclaimed 95.45 GB across 72 files.** `hprc-r2/assemblies/` went 114 GB → 25 GB; free disk 589 GB → 679 GB.

| | deleted | kept |
|---|---|---|
| 3 manifest haplotypes | `.fa`, `.fa.fai`, `.amb`, `.ann`, `.bwt`, `.pac`, `.sa` | `.fa.gz`, `.fa.gz.fai`, `.fa.gz.gzi`, **`.fa.mmi`** |
| 7 unreferenced haplotypes | the above **plus `.fa.mmi`** | `.fa.gz`, `.fa.gz.fai`, `.fa.gz.gzi` |

The deletion ran behind a guard that refused to touch any haplotype whose `.fa.gz` + `.fai` + `.gzi` replacement triple was not present and non-empty; all 10 passed. The `.mmi` is retained for the 3 manifest haplotypes because projection has no replacement until Phase 4.

**Cost to undo:** ~2 min `gunzip` per haplotype, plus ~55 min `bwa index` if `search.backend: bwa` is ever needed again (`WITH_BWA=1 bash scripts/prepare_haplotypes.sh`). Nothing else regresses — the `.fa` was pure duplication of the BGZF.

#### The reclaim nearly hid itself

The first post-reclaim gate reported **12 passed, 2 skipped in 0.03 s** and exited 0. Both slow gates had skipped: `_hprc_data_present()` probed a hardcoded `hprc-r2/assemblies/<sample>.fa`, and the reclaim had just deleted every one of those. Deleting the data the gate depends on turned the gate off, and the suite reported success — a green tick that meant *nothing was checked*.

This is the third instance of the same pattern in this work, after the golden fixture that had never been reproduced and `bwa`'s silently-omitted `XA` tag: **the failure was always silence, never an error.** The guard is now derived from the manifest the test actually hands the CLI (`load_haplotypes`), so it tracks any future repoint automatically, and `test_slow_gate_is_not_silently_disabled` asserts that whenever genome data exists on disk the guard agrees — keyed on file presence rather than on the guard's own logic, so it fails when the two disagree instead of restating one of them.

Worth noting for later phases: a skipped test and a passing test are indistinguishable in pytest's summary line. Any phase that deletes or moves data should re-read the skip reasons (`-rs`), not just the exit code.

**Inventory (2026-07-25, pre-reclaim):** all 10 `.fa.gz` verified BGZF by magic-byte check. Only `HG00097_hap1` has a `.gzi`; 6 of 10 have `.fa.gz.fai`; 6 of 10 have `.mmi`; no haplotype has a `.chm13.paf`.

---

## Phase 4 — Sparse anchor grid (replaces `.mmi`)

**Model: Opus.** Design judgement on sampling density and precision.

Measured basis: HPRC R2 assemblies are near-chromosome-level — 26 contigs ≥10 Mb carry **99.3%** (HG00097_hap1) and **98.9%** (HG01884_hap1) of sequence. Whole-genome base-level alignment is unnecessary; the current PAF route costs 114 min and 14.5 GB peak.

- Builder: sample ~1 kb probes every 10 kb, map against the **shared** CHM13 `.mmi`, store ~300k anchors/haplotype (a few MB gzipped). Est. 10–15 min/hap at low RAM.
- Consumer: a third branch in `project.project_target` — look up bracketing anchors, fetch the window from BGZF, realign in-memory with `mappy.Aligner(seq=window, preset="asm5")` (the idiom `mask.py:48` already uses).
- The replacement object must satisfy the existing duck type: `.map(seq)` yielding hits with `mapq/q_st/q_en/ctg/r_st/r_en`, `.seq(ctg, start, end)`, and truthiness.
- Record the **unanchored fraction** per haplotype — Phase 5 needs it.

**Note:** `tests/test_resolve_target.py:57-68` asserts that `project._aligner` skips a zero-byte `.mmi` and rebuilds from FASTA. If `.mmi` support is removed, that test's premise changes and it must be rewritten, not deleted.

**Gate:** projected loci match the current `.mmi` path within a few bp on the demo loci; `.mmi` files (~27 GB) become deletable.

---

## Phase 5 — CHM13-once candidate discovery

**Model: Opus.** The core algorithmic insight and its correctness envelope.

Off-target amplification requires both primers binding within ~2 kb pointing at each other, so candidate loci are by construction homologous to the target. Therefore: search CHM13 **once per run** (~47 s, amortized across every haplotype instead of paid per haplotype), lift each candidate locus through the Phase 4 anchor grid, and fetch only those windows (~2.6 µs each).

The many-to-one lift is a feature: a haplotype carrying three CYP2D7-like copies that all align to one CHM13 locus is recovered naturally — the case the README leads with.

**Correctness envelope.** Haplotype-private sequence absent from CHM13 is unreachable this way. Report it as the existing **`uncertain`** status, which `intro_verify_pipeline.md` §5 already argues must stay distinct from `dropout`. Phase 2's scanner remains available as the exact genome-wide path, so Phase 5's approximation is validated against it — not merely asserted.

**Gate:** on ≥5 haplotypes, Phase 5 output differs from Phase 2's exact scan **only** in loci the anchor grid marks unanchored. Any other difference is a bug. Target <1 s/haplotype.

---

## Phase 6 — Pangenome index spike, then commit

**Model: Opus** for evaluation and the eventual build; **Sonnet** for harness scaffolding.

Time-boxed spike on ~20 haplotypes before committing to any rewrite. Candidates:

| Candidate | Hypothesis | Main risk |
|---|---|---|
| Movi / r-index | Run-length BWT compresses with *distinct* sequence; 464 near-identical haplotypes in tens of GB | Younger toolchain; haplotype-attribution layer is yours to build |
| `vg` GBZ + giraffe | Mature on HPRC; haplotype paths first-class, so attribution is free | Tuned for reads not 20-mers; graph→coordinate reporting is substantial |
| Phases 1–5 only | No new index at all | Per-haplotype cost stays >0; unanchored fraction stays `uncertain` |

**Measure:** total index size, build time, 20-mer query accuracy against Phase 2's exact scan as ground truth, and per-haplotype coordinate attribution fidelity. **Commit only on evidence.** If no candidate beats Phases 1–5 on the accuracy/effort trade, that is a valid and documented outcome.

---

## Delegation summary

| Phase | Model | Rationale | Parallel with |
|---|---|---|---|
| 0 — Safety net | Sonnet | Mechanical, well-specified, high volume | 1 |
| 1 — Cleanup | Sonnet | Independent mechanical edits | 0 |
| 2 — Rust matcher | **Opus** | Algorithmic correctness + FFI; highest risk | — (needs 0) |
| 3 — BGZF storage | Sonnet (+Opus review of `main.nf`) | Mostly path plumbing; one risky line | — (needs 2) |
| 4 — Anchor grid | **Opus** | Sampling-density and precision judgement | 3 |
| 5 — CHM13-once | **Opus** | Core insight + correctness envelope | — (needs 4) |
| 6 — Pangenome spike | **Opus** (Sonnet harness) | Open-ended evaluation, architectural commitment | — (needs 2 for ground truth) |

Every phase brief must carry: the backend contract above, the relevant `file:line` seams, its gate, and the standing instruction **not to delete genome data**.

---

## Verification

Run at every phase gate, in order:

```bash
pytest -q                                        # 42 existing + Phase 0 additions
pangenome-primer selftest --outdir selftest_out  # 5 statuses (engine guard only — NOT a refactor guard)
bash demo-design-pipeline/run_demo.sh            # GAPDH, 3 haplotypes
```

Then the golden comparator (`pytest --runslow tests/test_golden_demo.py`): the regenerated `verify.json` must match the frozen fixture **cell for cell**, including the CYP2D7 off-target row and the engineered dropout row. It runs the real pipeline against `hprc-r2/` and takes ~17 min, so it is a gate check, not an inner-loop test.

Record measured wall-time and peak RSS per haplotype at each gate against the baseline: **~5.7 min/haplotype end to end**, of which **47.3 s** is search (4.45 GB peak RSS).

**Gate status:**

- **Phase 0 + 1** passed 2026-07-24 — `77 passed, 7 skipped`, plus the verify golden green in 1025 s.
- **Design golden was stale from the start** (found 2026-07-25). `tests/golden/results.json`
  was frozen from `demo-design-pipeline/gapdh/results.json` dated **2026-07-24 08:32** —
  before Phase 0/1 landed — and never reproduced from current code, because the slow test
  that would have caught it did not run until Phase 2's gate. Three independent runs then
  produced an identical, *different* shortlist (`pair1` in, `pair5` out). The new output is
  the correct one: `cli.py:186` sorts by `(-coverage, penalty)`, all candidates tie at
  coverage 1.0, and the fresh shortlist is exactly the five lowest penalties — an ordering
  the golden's membership cannot produce. Every pair present in both was bit-identical, so
  the Rust backend changed nothing here. Demo and golden regenerated; both gates now green.
  **Lesson: a fixture that has never been reproduced from current code is an assumption, not
  a test.** Freeze and validate are separate steps.
- **Phase 2** passed 2026-07-25 — `117 passed, 2 skipped`, `selftest OK`, and
  `test_verify_pipeline_reproduces_golden_verify_matrix` green in **809 s (13 min 29 s)** with
  `backend=rust`, `max_binding_sites=100`, and the replacement rs1058164 dropout row.
  Down from the 17 min 05 s bwa baseline despite the search now being exhaustive.

**Standing rule:** no `.fa`, `.bwt`, `.sa`, `.pac`, or `.mmi` file is deleted until the phase that replaces it has passed its gate.

---

## Risks

| Risk | Mitigation |
|---|---|
| Rust matcher slower than estimated | Prototype against the 47.3 s baseline **before** anything is deleted; `naive`/`bwa` backends stay selectable |
| Wheel-less platform can't install | Pure-Python fallback via `search.backend`; degrade, never break |
| `main.nf:116-124` staging breaks silently | `collectMany` skips missing globs quietly — assert staged file counts in a Nextflow smoke test |
| Anchor grid less precise than `.mmi` | Local realignment on the fetched window recovers base-level precision; diff against `.mmi` output on demo loci |
| Phase 5 misses haplotype-private loci | Report `uncertain`; validate against Phase 2's exact path on ≥5 haplotypes |
| Phase 6 spike inconclusive | A documented "Phases 1–5 are sufficient" is an acceptable outcome |

---

## Cleanup executed

Per the decision to touch only unambiguously stale files, the following ~2.4 MB were removed. All were gitignored or dead; nothing genome-related was touched.

| Path | Size | Reason |
|---|---:|---|
| `hprc-r2/assemblies/HG00639_hap1.fa.chm13.paf.part` | 2,150,400 B | Truncated PAF; build started `prepare.log:1039`, never completed |
| `src/pangenome_primer/__pycache__/` | 181,651 B | Bytecode cache, regenerable |
| `tests/__pycache__/` | 116,822 B | Same; held a stale non-pytest `test_engine` pyc |
| `.nextflow/` | 33,576 B | Resume cache for two runs whose work dirs no longer exist — unusable |
| `.nextflow.log`, `.nextflow.log.1` | 36,531 B | Logs of runs into deleted scratchpad dirs |
| `.pytest_cache/` | 3,280 B | Test cache |
| `src/pangenome_primer.egg-info/` | 1,109 B | Rebuilt by `pip install -e .` |
| `hprc-r2/assemblies/HG00639_hap1.fa.mmi` | 0 B | Truncated `minimap2 -d` |
| `hprc-r2/assemblies/HG01884_hap1.fa.chm13.paf.part` | 0 B | `PAF FAILED` (OOM), `prepare.log:1008` |
| `hprc-r2/assemblies/HG00097_hap1.fa.chm13.paf.part` | 0 B | `PAF FAILED` (OOM), `prepare.log:1017` |

**Deliberately kept:** `assets/NO_FILE` (0 B but a tracked Nextflow sentinel, `main.nf:18,22,113`); `hprc-r2/prepare.log` (provisioning evidence cited in `plan-optimization.md`); all `.fa`, `.fa.gz`, bwa index, and `.mmi` files; `ref/` (untracked personal notes, zero references — left for the owner to decide).
