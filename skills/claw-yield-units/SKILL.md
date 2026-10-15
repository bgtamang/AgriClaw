---
name: claw-yield-units
description: Convert crop grain yields between bushels per acre (bu/ac), lb/ac, kg/ha, t/ha (Mg/ha) and g/m2 using USDA standard bushel weights (corn 56 lb, soybean 60 lb, wheat 60 lb, barley 48 lb, oats 32 lb, ...), adjust yields to a standard moisture content (corn 15.5%, soybean 13%, wheat 13.5%) with the dry-matter formula, and scale plot harvest weights to per-area yield from harvested rows, row spacing and length (ft, in, m), for one value or a CSV of plots. Use when the user asks to convert yield units, convert bushels to kilograms per hectare, adjust yield for grain moisture or shrink, calculate yield from plot weights, or compare yields reported in different units or moisture bases.
license: MIT
compatibility: Python >= 3.9 with matplotlib (figure for plot and batch runs); standard library otherwise. No internet needed.
metadata:
  version: "1.0.0"
  domain: field
  language: Python
  status: tested
  repository: https://github.com/bgtamang/AgriClaw
---

# Yield units, moisture and plot-to-area conversion

## What it does

Converts grain yields between bu/ac, lb/ac, cwt/ac, kg/ha, t/ha and g/m² with exact unit definitions and crop-specific bushel weights. It adjusts yields to a standard (or chosen) moisture content, keeping dry matter constant, and turns plot harvest weights into per-area yield from the harvested area. It works on a single value or a CSV of plots, and every output carries its bushel weight and moisture basis.

## When to use / when not to use

- Use for unit conversions, moisture adjustment ("yield at 15.5%"), and yield from combine plot weights in breeding or agronomy trials.
- Use **claw-trial-analysis** (and its design-aware models) to compare entries. This skill only converts numbers.
- Don't use it for forage or biomass dry-matter yields with a different moisture convention, or for test-weight discounts. Elevator shrink tables also add handling loss and differ by buyer.

## Inputs

| Mode | Flags | Notes |
|---|---|---|
| One value | `--value 60 --from bu/ac` | Units: `kg/ha`, `t/ha`, `g/m2`, `lb/ac`, `cwt/ac`, `bu/ac`. |
| One plot | `--plot-weight 9.8 --weight-unit lb` plus `--rows 2 --row-spacing 30 --row-spacing-unit in --length 17.5 --length-unit ft`, or `--plot-area 7.6 --area-unit m2` | Length units `m`, `cm`, `in`, `ft`. Area units `m2`, `ft2`, `ha`, `ac`. |
| Many plots | `--input plots.csv` | Column `plot_weight` (or `yield` with `--from`). Optional per-row `harvest_moisture_pct`, `rows`, `row_spacing`, `length`, `plot_area`; otherwise the CLI values apply. Other columns are kept. |
| Crop basis | `--crop corn\|soybean\|wheat\|sorghum\|barley\|oats\|rye\|flaxseed\|rice\|canola\|sunflower` | Sets the bushel weight and standard moisture. `--lb-per-bu` and `--to-moisture` override them. |
| Moisture | `--moisture 21.3` (harvest, % wet basis) and `--to-moisture` | Without a harvest moisture nothing is adjusted, and the basis is reported as unknown. |

Always ask:

- the crop, and whether the user's moisture is wet-basis percent;
- the harvested area: rows actually harvested, row spacing, and harvested length after end trimming (not the planted plot length), with units;
- which moisture basis the results should be on (market standard, a trial network's standard, or dry matter, `--to-moisture 0`);
- for sorghum, barley, rye, flax and rice, which standard moisture applies. There is no default, because the sources disagree or are silent.

## How to run

Commands use `${CLAUDE_SKILL_DIR}`, this skill's folder: Claude Code fills it in; in other tools, replace it with the folder that contains this SKILL.md.

```bash
python "${CLAUDE_SKILL_DIR}/scripts/yield_units.py" --demo
python "${CLAUDE_SKILL_DIR}/scripts/yield_units.py" --value 60 --from bu/ac --crop soybean
python "${CLAUDE_SKILL_DIR}/scripts/yield_units.py" --value 11.2 --from t/ha --crop corn --moisture 21 --to-moisture 15.5
python "${CLAUDE_SKILL_DIR}/scripts/yield_units.py" --plot-weight 9.8 --weight-unit lb --rows 2 --row-spacing 30 --row-spacing-unit in --length 17.5 --length-unit ft --moisture 14.1 --crop soybean
python "${CLAUDE_SKILL_DIR}/scripts/yield_units.py" --input combine_plots.csv --crop corn --weight-unit lb --rows 2 --row-spacing 30 --row-spacing-unit in --length 17.5 --length-unit ft --out results/yield
```

The demo uses 12 synthetic soybean plots (weights in lb at harvest moisture), harvested as 2 rows × 30 in × 17.5 ft and reported at 13% moisture.

## Outputs

- `report.md`: warnings, then a **Basis** section (bushel weight and source, moisture basis and source, the lb/ac to kg/ha factor), a results table, the method and references.
- `tables/yields.csv`: input columns plus `area_m2`, `harvest_moisture_pct`, `moisture_basis`, `lb_per_bu`, `yield_kg_ha`, `yield_t_ha`, `yield_g_m2`, `yield_lb_ac`, `yield_bu_ac`, `dry_matter_kg_ha` and `flags` (plausibility checks).
- `figures/yields.png` (plot and batch runs only; flagged records in red). A single-value conversion has no figure.
- `reproducibility/`.

**When you report any number from this skill, always state the bushel weight (lb/bu) and the moisture basis it is on**, for example "4,035 kg/ha (60 bu/ac at 13% moisture, 60 lb/bu)". Conversions without a stated basis are a common source of 10–20% errors.

## Expert decisions

- **Exact units.** 1 lb = 0.45359237 kg (NIST SP 811). 1 ft = 0.3048 m, the international foot, which since 2023 replaces the US survey foot in all applications (85 FR 62698). So 1 acre = 43,560 ft² = 4046.8564224 m² and 1 lb/ac = 1.12085116 kg/ha. AH-697 Table 5 prints "1 pound = 0.455 924 kilogram", a typo; its 453.5924 g line and its bushel conversions match the exact value.
- **Bushel weights** (USDA AH-697, Tables 5–6): corn (shelled), sorghum grain, rye and flaxseed 56 lb; wheat and soybeans 60 lb; barley 48 lb; oats 32 lb; rough rice 45 lb.
  - AH-697 gives only ranges for rapeseed (50–60 lb) and sunflower (24–32 lb), so canola 50 lb and sunflower 25 lb come from UF/IFAS AG442.
  - Arkansas FSA-1078 lists grain sorghum at 55 lb/bu.
- **Standard moistures** (UF/IFAS AG442 Table 1): corn 15.5%, soybean 13%, wheat 13.5%, oats 13.5%, canola 8.5%, sunflower 10%.
  - Corn, soybean and wheat agree with Arkansas FSA-1078.
  - Sorghum has no default (13% in AG442 versus 14% in FSA-1078), and nor do barley, rye, flax or rice. The script stops and asks for `--to-moisture`.
- **Moisture adjustment.** W_std = W × (100 − M_harvest)/(100 − M_std), on a wet basis with dry matter conserved (AG442; FSA-1078).
  - Tested on AG442's example (10,000 lb/ac at 20% = 9,467 lb/ac at 15.5%) and FSA-1078's (261,000 lb dry matter = 300,000 lb at 13%, 310,714 lb at 16%, 290,000 lb at 10%).
  - Buyers' shrink factors add handling loss; this skill reports water shrink only.
- **Bushels and moisture.** A bushel weight applies at the standard moisture. Adjust to the standard first, then convert to bu/ac. A non-standard `--to-moisture` triggers a warning.
- **Plot area.** Area = harvested rows × row spacing × harvested length, or a given area. Yield = weight / area. Plausibility flags (area 0.5 m²–1 ha, row spacing 0.1–2 m, length 1–500 m, harvest moisture 5–40%) are adjustable rules of thumb meant to catch unit mistakes. They never change values.

## Checks before reporting

- Quote the bushel weight and moisture basis with every number (see Outputs).
- Read the warnings, and stop and ask the user if:
  - the harvest moisture is missing but results will be compared with other trials or market yields;
  - a crop without a default standard moisture is used;
  - any `flags` appear (a wrong unit is the usual cause, e.g. row spacing in inches given as metres);
  - the user's "plot length" may be the planted length rather than the harvested length.
- Plausibility: US soybean trials are usually 2–6 t/ha and corn 8–16 t/ha (rule of thumb). Values far outside suggest a unit or area error.

## References

- USDA Economic Research Service (1992). Weights, Measures, and Conversion Factors for Agricultural Commodities and Their Products. Agricultural Handbook No. 697. https://www.ers.usda.gov/publications/pub-details?pubid=41881
- Mulvaney, M.J. & Devkota, P. (2020). Adjusting Crop Yield to a Standard Moisture Content. UF/IFAS Extension SS-AGR-443/AG442. EDIS 2020(3). https://doi.org/10.32473/edis-ag442-2020
- Sadaka, S., Atungulu, G.G. & Olatunde, G. Understanding Grain Shrinkage and Expansion. FSA1078, University of Arkansas Division of Agriculture Research & Extension. https://www.uaex.uada.edu/publications/pdf/FSA-1078.pdf
- Thompson, A. & Taylor, B.N. (2008). Guide for the Use of the International System of Units (SI). NIST Special Publication 811. https://doi.org/10.6028/NIST.SP.811e2008
- NIST & NOAA (2020). Deprecation of the United States (U.S.) Survey Foot. Federal Register 85 FR 62698. https://www.federalregister.gov/documents/2020/10/05/2020-21902/deprecation-of-the-united-states-us-survey-foot
