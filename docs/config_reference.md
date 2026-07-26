# Configuration reference (`config/defaults.yaml`)

Every value is overridable on the CLI or via `--config <file.yaml>`. Terminology is in [`CONTEXT.md`](../CONTEXT.md).

Two thresholds are easy to confuse:

- **`dropout.thermo.tm_drop_max`** — a per-**site** test: "does *this* primer stick *here*?"
- **`rank.min_coverage`** — a per-**pair** gate: "does the pair work across enough people?"

## `design:` — Primer3 candidate generation

| key | meaning | unit |
|---|---|---|
| `product_size_min` / `_max` | Allowed amplicon size range. | bp |
| `primer_len_min` / `_opt` / `_max` | Primer length: min / optimum / max. | nt |
| `primer_tm_min` / `_opt` / `_max` | Melting temperature (**Tm**): min / optimum / max. | °C |
| `n_candidates` | Pairs Primer3 returns per locus. | count |
| `flank` | CHM13 sequence projected onto each haplotype either side of the target; defines the design template and the on-target window. | bp |

## `mask:` — steer primers off variable sites

| key | meaning | unit |
|---|---|---|
| `min_allele_freq` | Mask a template position when its variant frequency across the subset reaches this. `0.05` = mask sites varying in ≥5%. | fraction |
| `three_prime_weight_nt` | *(Reserved — not wired.)* Intended width of an extra-penalty zone near the 3′ end. v1 masking uses hard `SEQUENCE_EXCLUDED_REGION`s, not a graded penalty. | nt |
| `merge_gap` | Variant clusters within this distance coalesce into one masked zone — a cluster is one bad zone for a ~20 nt primer. | bp |
| `max_regions` | Cap on masked intervals; Primer3 allows at most 200. Over the cap, the widest zones win. | count |

## `search:` — genome-wide binding-site search

| key | meaning |
|---|---|
| `backend` | `rust` (default) is the compiled scanner: streams the BGZF, needs no index, and is **exhaustive** — every position is proven to match or not, so the site list is never a sample or a truncation. `naive` is pure Python and refuses references above 50 Mb (a correctness reference for fixtures, not a genome-scale option). `auto` behaves like `rust`. **No fallback:** without the compiled extension this raises, naming `maturin develop --release`. |
| `threads` | Worker threads for the scanner. `null` (default) = every core, capped at 16. See below. |
| `max_mismatches` | Search budget — report sites with up to this many mismatches. Widens the *net* only; whether a site amplifies is decided by `dropout`, not by this number. |
| `max_binding_sites` | Cap on amplification-**competent** sites per primer. Above it the pair is not scored and the cell reads `>N binding sites`, because a primer binding that many places amplifies indiscriminately and enumerating its products is meaningless. `0`/null disables. The default `100` comes from measured separation: well-behaved demo primers top out at 28 competent sites, an *Alu* 20-mer at ~106,000. Values ≤28 would suppress the genuine CYP2D7 paralog row. |

### Setting `search.threads`

Leave it `null` while haplotypes are scanned **one at a time** — the CLI loops them sequentially, and a solo scan should take the machine.

**Set it low when scans run concurrently** (the Nextflow `EVALUATE` fan-out does). Aim for `threads × concurrent-haplotypes ≈ core count`, threads in the 2–4 range. Measured on 16 cores, 4 threads × 4 concurrent beat the 16-thread default by 19%. Peak RSS is ~400 MB regardless of thread count, so memory never argues for a high value.

`PGP_SCAN_THREADS` then `RAYON_NUM_THREADS` override it. The pool is built once per process, so the first scan fixes the count.

**Running on anything bigger than a laptop?** See [`sizing.md`](sizing.md) — per-host recommendations, Nextflow sizing profiles, and `scripts/sizing_sweep.sh` to re-measure on your own machine.

## `dropout:` — does a primer actually bind at a site?

`mode` — `thermo` (physics-based, primary) or `rule` (mismatch counts).

### `dropout.thermo` — nearest-neighbor Tm model (via primer3)

The first two are **decision thresholds**; the rest are **wet-lab conditions** fed to primer3's salt-corrected (SantaLucia) Tm calculation. Set them to match your actual protocol.

| key | meaning | unit |
|---|---|---|
| `tm_drop_max` | How far a site's annealing Tm may fall **below the design Tm** before it is a dropout. `5.0` = a site >5 °C under the intended Tm won't amplify. | °C |
| `three_prime_hard_nt` | Size of the 3′ window where *any* mismatch is an automatic dropout — the polymerase cannot extend off a mismatched 3′ end. | nt |
| `mv_conc` | Monovalent cations (Na⁺/K⁺). Stabilises the duplex, raising Tm. | mM |
| `dv_conc` | Divalent cations (**Mg²⁺**). Strongly stabilising, required by the polymerase; a major Tm driver. | mM |
| `dntp_conc` | Free dNTPs. They chelate Mg²⁺, lowering effective Mg²⁺, so Tm is corrected for it. | mM |
| `dna_conc` | Primer concentration. Annealing Tm depends on strand concentration; 250 nM (0.25 µM) is standard. | nM |

### `dropout.rule` — the simple alternative

Reads as: *binds iff ≤ `max_total_mm` mismatches total **and** ≤ `max_3prime_mm` within the 3′-terminal `suffix` bases.*

| key | meaning |
|---|---|
| `max_total_mm` | Max mismatches over the whole primer. |
| `max_3prime_mm` | Max mismatches inside the 3′ window. |
| `suffix` | Size of that window (nt). |

Defaults `2 / 0 / 5` — "≤2 mismatches total, with a perfectly matched 3′ pentamer", the conservative policy from the reference notes.

## `rank:` — hard filters a pair must clear to *pass*

Every gate applies over **evaluable** haplotypes (total − uncertain).

| key | meaning |
|---|---|
| `min_coverage` | Minimum on-target coverage. `0.95` = amplifies in ≥95%; ≥0.99 is suggested for general assays. |
| `max_off_target` | Max haplotypes showing an off-target product. `0` = zero tolerance. |
| `max_multi_product` | Max haplotypes showing extra bands. |
| `min_evaluable` | Minimum haplotypes that must project for the pair to be judged at all. |

Survivors are ordered by a transparent tie-break (coverage → uniqueness → off-target count → Primer3 penalty); see `src/pangenome_primer/rank.py`.

## Verify-mode knobs (CLI flags, not in this file)

`pangenome-primer verify` reuses the `dropout` model but takes its size bounds on the command line: `--max-amplicon` (max off-target product size, default 2000 bp) and `--size-tolerance` (flag on-target sizes deviating from the expected span by more than this, default 20 bp).
