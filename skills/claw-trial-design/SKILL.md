---
name: claw-trial-design
description: Randomise and lay out field trials with agricolae and FielDHub - RCBD, alpha-lattice (incomplete blocks), row-column, Latin square, CRD, augmented RCBD with checks, partially replicated (p-rep), diagonal-check, split-plot, strip-plot and factorial designs. Writes a field book CSV (plot, row, column, rep, block, entry, check flag, factor levels), a field map PNG and a design summary with efficiency and seed. Use when the user asks to design or randomise a trial, make a field layout, plot plan, planting map or field book, choose between RCBD and alpha-lattice, place checks in an early-generation yield trial, or set up a split-plot or strip-plot experiment.
license: MIT
compatibility: R >= 4.2 with agricolae and FielDHub (openssl optional, for SHA-256 checksums). No internet needed.
metadata:
  version: "1.0.0"
  domain: trials
  language: R
  status: tested
  repository: https://github.com/bgtamang/AgriClaw
---

# Field trial design and layout

## What it does

Randomises a trial with the established design functions in agricolae and FielDHub, then places every plot on the field grid so the physical layout respects the design: blocks and replicates are compact rectangles, split-plot whole plots are contiguous, strip-plot strips stay intact, and every augmented block holds every check. Plot numbers (serpentine or cartesian) are assigned after the layout, so numbering never moves an entry. Every design is checked for its defining properties before any file is written.

## When to use / when not to use

- Use to design a new trial, re-create a layout from a recorded seed, or compare designs (for example RCBD against alpha) by efficiency and error degrees of freedom.
- Use `claw-fieldbook` next, to turn `field_book.csv` into a Field Book app import file and later to merge and QC the collected data.
- Don't use it to analyse data: use `claw-trial-analysis` (one trial) or `claw-met-gxe` (several environments).
- Don't use it for designs it does not build (balanced incomplete blocks, Youden, split-split-plot, multi-location sparse testing). Say so rather than forcing another design.

## Inputs

Either `--entries N` (named G01..GN; add `--checks N` for CHK1..CHKN) or `--entry-list FILE.csv`:

| Column | Meaning |
|---|---|
| `entry` | required; unique entry names (any characters) |
| `is_check` | optional; TRUE/FALSE, yes/no or 1/0 |
| `reps` | optional, p-rep only: copies of each entry |
| anything else | carried into the field book (pedigree, source, seed lot) |

Design parameters: `--reps` (replicates; squares for `latin`), `--blocks` (augmented), `--block-size` (alpha), `--rep-rows` (rowcol), `--factor-a/-b/-c 'name:level1,level2'` (split, strip, factorial), `--prep-fraction`, `--check-reps`. Field: `--field-rows`, `--field-cols`, `--block-rows`, `--numbering`, `--plot-start`. Always `--seed`.

**Ask the user, never guess:**

- the design, or the facts that decide it: number of entries, replicates the seed supply allows, whether this is an unreplicated early-generation trial, known field gradients;
- which entries are checks, and how many blocks (augmented) or copies (p-rep);
- the field grid (plots per row and number of rows) if the trial must fit a field; without it the script picks the most compact layout;
- for split-plot, which factor is hard to change on small plots (it becomes the whole-plot factor); for strip-plot, which factor runs along rows;
- the seed, if they want to reproduce an earlier layout.

## How to run

Commands use `${CLAUDE_SKILL_DIR}`, this skill's folder: Claude Code fills it in; in other tools, replace it with the folder that contains this SKILL.md.

```bash
Rscript "${CLAUDE_SKILL_DIR}/scripts/trial_design.R" --demo
Rscript "${CLAUDE_SKILL_DIR}/scripts/trial_design.R" --design rcbd --entries 20 --reps 4 --seed 42
Rscript "${CLAUDE_SKILL_DIR}/scripts/trial_design.R" --design alpha --entry-list entries.csv --reps 3 --block-size 6 --field-rows 12 --field-cols 30 --seed 7 --out alpha_2027/
Rscript "${CLAUDE_SKILL_DIR}/scripts/trial_design.R" --design augmented --entries 200 --checks 4 --blocks 10 --seed 11
Rscript "${CLAUDE_SKILL_DIR}/scripts/trial_design.R" --design prep --entry-list lines.csv --prep-fraction 0.25 --check-reps 4 --field-rows 20 --field-cols 18 --seed 3
Rscript "${CLAUDE_SKILL_DIR}/scripts/trial_design.R" --design diagonal --entries 300 --checks 5 --field-rows 22 --field-cols 15 --seed 5
Rscript "${CLAUDE_SKILL_DIR}/scripts/trial_design.R" --design split --factor-a "irrigation:dry,wet" --entry-list varieties.csv --reps 4 --seed 9
Rscript "${CLAUDE_SKILL_DIR}/scripts/trial_design.R" --design strip --factor-a "tillage:conv,notill" --factor-b "N:0,60,120" --reps 4 --seed 10
```

The demo lays out the 30 entries in `demo/entries.csv` as an alpha design (3 replicates, blocks of 5, seed 2026). `--help` lists every option.

## Outputs

- `tables/field_book.csv`: one row per plot. `trial`, `plot_id` (`<trial>-<plot>`, unique), `plot`, `row`, `col`, `design` (the design name, read by `claw-trial-analysis`), `rep` (replicate or complete block; the square for `latin`; NA when the design has none), `block` (incomplete or augmented block, numbered across the trial), `whole_plot` (split), `entry_no`, `entry`, `is_check`, `is_filler`, factor columns, and extra entry-list columns. For split/strip with entries as factor B, `entry` is the combination and `entry_b` holds the entry.
- `tables/design_summary.csv`: design, **package function used**, package version, seed, RNG kind, layout, efficiency factor and bound (alpha), A-efficiencies (rowcol), error df.
- `tables/design_checks.csv`: the design-property checks; the run stops if any fails.
- `figures/field_map.png`: row 1 at the top; thick outlines = replicates (or augmented blocks), thin = incomplete blocks, grey = whole plots; checks orange, fillers grey.
- `report.md`: warnings first, then design, efficiency, checks, map, references. `reproducibility/commands.sh` always contains the seed.

## Expert decisions

Choosing a design (confirm with the user):

| Situation | Design | Notes |
|---|---|---|
| Few entries (rule of thumb: up to about 20-30) on a uniform area | `rcbd` | Blocks grow with the entry count and stop being uniform; the script warns above 30 entries. |
| Many entries, replicated | `alpha` | Resolvable incomplete blocks (Patterson & Williams 1976). Block size near sqrt(entries) is a rule of thumb: small blocks with many of them lower the efficiency factor, which the summary reports. |
| Gradients in two directions | `rowcol` | Rows and columns within each replicate are both blocking factors (Williams, John & Whitaker 2006). |
| Two gradients, 3-12 treatments | `latin` | Each treatment once per row and column; one square of 4 has only 6 error df, so use `--reps 2` or more squares. |
| Unreplicated early-generation lines with checks | `augmented` | Checks in every block (Federer 1956, 1961). |
| Unreplicated lines, spatial analysis planned | `prep` or `diagonal` | p-rep replicates a share of the test lines (Cullis et al. 2006); diagonal designs spread checks systematically (Clarke & Stefanova 2011). Analyse both with a spatial model in `claw-trial-analysis`. |
| Factor that needs large plots (irrigation, planting date) x a second factor | `split` | The hard-to-change factor goes on whole plots; it gets fewer error df. |
| Both factors need strips (tillage x N with machinery) | `strip` | Three error strata; the interaction is estimated most precisely. |
| All combinations fit small plots | `factorial` | Factorial treatments in RCBD. |
| Uniform conditions (growth chamber) | `crd` | No blocking. |

Defaults and thresholds:

- **Package functions.** agricolae: `design.crd`, `design.rcbd`, `design.lsd`, `design.alpha`, `design.dau`, `design.split`, `design.strip`, `design.ab`. FielDHub: `alpha_lattice`, `row_column`, `partially_replicated`, `diagonal_arrangement`. The function used is recorded in the summary. agricolae's `zigzag()` renumbers agricolae's own one-block-per-row sketch; because this skill lays blocks out as rectangles that can span several rows, it numbers plots itself after the layout (same serpentine rule, applied to the final grid).
- **Alpha engine (`--alpha-engine auto`).** agricolae builds alpha(0,1) designs (no pair of entries shares more than one block) only for r = 2-4 under its series conditions. FielDHub's `alpha_lattice` covers more cases but gives no concurrence guarantee. In auto mode both are generated with the same seed, and the one with the higher efficiency factor is kept; both values are reported. The efficiency factor is the harmonic mean of the canonical efficiency factors, computed from the incidence matrix; agricolae's printed "Efficiency" is the upper bound (v-1)(r-1)/[(v-1)(r-1)+r(s-1)], which square lattices attain (tested).
- **Entries are never dropped.** If the block size does not divide the entry count, the run stops and lists the block sizes that do, or how many entries (for example checks) to add.
- **Augmented error df** = (blocks - 1)(checks - 1). The script warns below 10. "About 10-12 or more" is a rule of thumb: the 5% two-sided t value is 2.23 at 10 df and 2.18 at 12 df, against 1.96 with unlimited df, so extra df beyond that buy little. Uneven blocks are padded with `FILLER` plots, which are flagged and never analysed.
- **p-rep share.** The default is 25% of test entries planted twice (`--prep-fraction 0.25`). FielDHub's vignette describes a typical 1:4 ratio, p = 25%, for Cullis et al. (2006); treat 20-30% as a rule of thumb and warn outside it. When the field size is given, the number of replicated entries is set to fill it, and the realised share is reported. Checks get `--check-reps 4` copies, an adjustable rule of thumb.
- **Layout search.** Each nesting level (replicate, block, plot) is tiled as a grid; all factorisations are searched, and the layout with the most compact blocks is chosen (then replicates, then the whole field), subject to `--field-rows/--field-cols/--block-rows`. Compactness is measured in plot units; if plots are long and narrow, give the field grid or `--block-rows` explicitly.
- **Seeds.** agricolae sets its own RNG kind (Super-Duper) through `set.seed(seed, kinds)`; the script restores the session RNG afterwards and records the kind. FielDHub uses R's default RNG. Without `--seed`, a seed is drawn and written to the summary and `commands.sh`.
- **Known FielDHub 1.5.0 bug.** `row_column()` fails ("object 'best_efficiencies' not found") for some seeds. The script retries with seed + 1, + 2, ... (deterministic), records the seed FielDHub used, and warns.

## Checks before reporting

- Read the Warnings in `report.md` and pass every one on.
- `design_checks.csv` must be all `TRUE`; the script stops otherwise. Never hand-edit a field book after randomisation (for example moving checks); re-run with new parameters instead.
- Alpha: efficiency factor below 0.9 x the upper bound means the blocks are too small. Suggest a block size nearer sqrt(entries).
- Augmented: error df below 10, or a single check, means weak tests of entry differences. Suggest more checks or blocks.
- FILLER plots must be planted (border seed) but excluded from analysis; tell the user how many there are.
- The map draws row 1 at the top. Confirm which field corner is row 1, column 1 before planting, and which way the planter moves.
- Stop and ask if the user's field cannot hold the design: the error lists field shapes that work.

## References

- de Mendiburu F (2023). agricolae: Statistical Procedures for Agricultural Research. R package version 1.3-7. https://doi.org/10.32614/CRAN.package.agricolae
- Murillo D, Gezan S, Heilman A, Walk T, Aparicio J, Horsley R (2021). FielDHub: A Shiny App for Design of Experiments in Life Sciences. Journal of Open Source Software 6(61):3122. https://doi.org/10.21105/joss.03122
- Patterson HD, Williams ER (1976). A new class of resolvable incomplete block designs. Biometrika 63(1):83-92. https://doi.org/10.1093/biomet/63.1.83
- Federer WT (1956). Augmented (or hoonuiaku) designs. Hawaiian Planters' Record 55(2):191-208. Cornell Biometrics Unit report BU-74-M: https://hdl.handle.net/1813/32841
- Federer WT (1961). Augmented designs with one-way elimination of heterogeneity. Biometrics 17(3):447-473. https://doi.org/10.2307/2527837
- Cullis BR, Smith AB, Coombes NE (2006). On the design of early generation variety trials with correlated data. Journal of Agricultural, Biological, and Environmental Statistics 11(4):381-393. https://doi.org/10.1198/108571106X154443
- Clarke GPY, Stefanova KT (2011). Optimal design for early-generation plant-breeding trials with unreplicated or partially replicated test lines. Australian & New Zealand Journal of Statistics 53(4):461-480. https://doi.org/10.1111/j.1467-842X.2011.00642.x
- Williams ER, John JA, Whitaker D (2006). Construction of resolvable spatial row-column designs. Biometrics 62(1):103-108. https://doi.org/10.1111/j.1541-0420.2005.00393.x
