# Configuration reference (`config/defaults.yaml`)

Every value is overridable on the CLI or via `--config <file.yaml>`. Terminology is defined
in [`CONTEXT.md`](../CONTEXT.md). This page explains each parameter — especially the
abbreviations in the `dropout` and `rank` sections.

Two thresholds are easy to confuse:
- **`dropout.thermo.tm_drop_max`** is a per-**site** binding test — "does *this* primer stick
  *here*?"
- **`rank.min_coverage`** is a per-**pair** summary gate — "does the pair work across enough
  people?"

---

## `design:` — Primer3 candidate generation

| key | meaning | unit |
|---|---|---|
| `product_size_min` / `product_size_max` | Allowed amplicon (PCR product) size range. | bp |
| `primer_len_min` / `primer_len_opt` / `primer_len_max` | Primer length: minimum / optimum / maximum. | nt |
| `primer_tm_min` / `primer_tm_opt` / `primer_tm_max` | Primer melting temperature (**Tm**): minimum / optimum / maximum. | °C |
| `n_candidates` | How many primer pairs Primer3 returns per locus. | count |
| `flank` | How much CHM13 sequence on each side of the target is projected onto each haplotype (defines the design template and the on-target window). | bp |

## `mask:` — variability-aware design (steer primers off variable sites)

| key | meaning | unit |
|---|---|---|
| `min_allele_freq` | A template position is masked (excluded from primer placement) when its variant frequency across the haplotype subset reaches this. `0.05` = mask sites varying in ≥5%. | fraction |
| `three_prime_weight_nt` | *(Reserved — not yet wired.)* Intended width of an extra-penalty zone near a candidate primer's 3′ end. v1 masking uses hard `SEQUENCE_EXCLUDED_REGION`s, not a graded 3′ penalty. | nt |
| `merge_gap` | Variant clusters within this many bp are coalesced into one masked zone (a cluster is one bad zone for a ~20 nt primer). | bp |
| `max_regions` | Cap on masked intervals — Primer3 allows at most 200 `SEQUENCE_EXCLUDED_REGION`s. Over the cap, the widest zones are kept. | count |

## `search:` — genome-wide binding-site search

| key | meaning |
|---|---|
| `backend` | *(Reserved — not yet wired.)* Intended to select the binding-site search engine. Today the engine picks automatically: `bwa` genome-wide on real assemblies, `naive` (pure-Python) for the in-memory/tiny-input path. |
| `max_mismatches` | Search budget — bwa reports candidate sites with up to this many mismatches. This only widens the *net*; whether a site actually amplifies is decided by the `dropout` model, not by this number. |

---

## `dropout:` — does a primer actually bind at a site?

`mode` — which model decides bind-vs-dropout: `thermo` (physics-based, primary) or `rule`
(simple mismatch counts).

### `dropout.thermo` — nearest-neighbor melting-temperature model (via primer3)

The first two are **decision thresholds**; the last four are **wet-lab reaction conditions**
fed to primer3's salt-corrected (SantaLucia) Tm calculation. Set the four concentrations to
match your actual PCR protocol.

| key | meaning | unit |
|---|---|---|
| `tm_drop_max` | How far the site's predicted annealing **Tm** may fall **below the design Tm** before it's a dropout. `5.0` = a site whose Tm is >5 °C under the intended Tm won't amplify. | °C |
| `three_prime_hard_nt` | Size of the **3′-terminal window** where *any* mismatch is an automatic dropout regardless of Tm — the polymerase cannot extend off a mismatched 3′ end. `2` = the last 2 bases must match. | nt |
| `mv_conc` | **M**ono**v**alent cation **conc**entration — Na⁺/K⁺ (e.g. KCl buffer). Salt stabilizes the duplex → raises Tm. | mM |
| `dv_conc` | **D**i**v**alent cation **conc**entration — **Mg²⁺** (magnesium). Strongly stabilizes DNA and is required by the polymerase; a major Tm driver. | mM |
| `dntp_conc` | Free **dNTP** (deoxynucleotide triphosphate) **conc**entration. dNTPs chelate Mg²⁺, lowering the *effective* Mg²⁺, so Tm is corrected for it. | mM |
| `dna_conc` | Primer (oligo/**DNA**) **conc**entration in the reaction. Annealing Tm depends on strand concentration; `250` nM (0.25 µM) is a standard primer amount. | nM |

### `dropout.rule` — the simple alternative (no thermodynamics)

`mm` = **mismatches**. Reads as: *a primer binds iff it has ≤ `max_total_mm` mismatches total
**and** ≤ `max_3prime_mm` mismatches within its 3′-terminal `suffix` bases.*

| key | meaning |
|---|---|
| `max_total_mm` | Max total mismatches over the whole primer. `2` = up to 2 anywhere. |
| `max_3prime_mm` | Max mismatches allowed inside the 3′ window. `0` = none tolerated near the 3′ end. |
| `suffix` | Size of that 3′ window (nt). `5` = the last 5 bases. |

Defaults `2 / 0 / 5` = "≤2 mismatches total, with a perfectly matched 3′ pentamer" — the
conservative policy from the reference notes (`max_total_mm=2; max_3prime_mm=0; suffix=5`).

---

## `rank:` — hard filters a primer pair must clear to *pass*

Every gate is applied over the **evaluable** haplotypes (total − `uncertain`).

| key | meaning |
|---|---|
| `min_coverage` | Minimum **on-target coverage** — the fraction of haplotypes that must produce the target amplicon. `0.95` = amplifies in ≥95%; the reference notes suggest ≥0.99 for general assays. |
| `max_off_target` | Max haplotypes allowed to show an **off-target** product while still passing. `0` = zero tolerance. |
| `max_multi_product` | Max haplotypes allowed to show **multiple products** (extra bands). `0` = zero tolerance. |
| `min_evaluable` | Minimum haplotypes that must project successfully to judge the pair at all. `1` = need ≥1 non-uncertain haplotype. |

Survivors of these gates are then ordered by a transparent tie-break (coverage → uniqueness →
off-target count → Primer3 penalty); see the ranking logic in `src/pangenome_primer/rank.py`.

## Verify-mode knobs (CLI flags, not in this file)

`pangenome-primer verify` reuses the `dropout` model but takes its size bounds on the command
line: `--max-amplicon` (max off-target product size, default 2000 bp) and `--size-tolerance`
(flag on-target sizes deviating from the expected span by more than this, default 20 bp).
