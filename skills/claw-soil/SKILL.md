---
name: claw-soil
description: Look up mapped soil properties for a US field location (latitude/longitude or map unit key) from the USDA-NRCS SSURGO soil survey through the Soil Data Access API - soil series/components, texture class by horizon, sand/silt/clay, organic matter, pH, bulk density, CEC, available water capacity and available water storage (AWS) to 25/50/100/150 cm, drainage class, water table depth, hydrologic group - with component-weighted map-unit values. Also classifies lab sand/silt/clay results into USDA texture classes. Use when the user asks what soil is at a field or trial site, its texture, drainage, pH or water-holding capacity, needs total available water for a water balance, or wants the USDA texture class of particle-size data.
license: MIT
compatibility: Python >= 3.9 with matplotlib; standard library otherwise. Internet access to sdmdataaccess.sc.egov.usda.gov for --lat/--lon and --mukey; --demo and --texture-input run offline. SSURGO covers the USA and some territories only.
metadata:
  version: "1.0.0"
  domain: field
  language: Python
  status: tested
  repository: https://github.com/bgtamang/AgriClaw
---

# Soil properties from SSURGO

## What it does

Finds the SSURGO map unit at a point and reports what the soil survey records for it. That covers its components (soil series and their share of the map unit), each horizon's texture (SSURGO's own class), sand/silt/clay, organic matter, pH, water retention and AWC. It also reports component-weighted map-unit averages for 0–30 and 0–100 cm, available water storage by depth, drainage class and water-table depth. A separate mode classifies the user's sand/silt/clay data with the USDA texture triangle.

## When to use / when not to use

- Use for site characterisation of US fields and trials, total available water (TAW) for **claw-water-balance**, and USDA texture classes of lab particle-size results (`--texture-input`).
- Don't use SSURGO values as a field soil test. They are mapped estimates. Fertiliser and lime decisions need a current soil test.
- Outside the USA, SSURGO has no data. Ask for local survey data or measured θFC/θWP instead (for claw-water-balance).
- For weather or thermal time, use **claw-weather** or **claw-gdd**.

## Inputs

| Input | Flag | Notes |
|---|---|---|
| Location | `--lat`, `--lon` (WGS84 decimal degrees) | Longitude is negative in the western hemisphere. A swapped pair usually lands outside SSURGO. |
| Map unit | `--mukey 242963` | Use this instead of a point if the user knows the map unit. |
| Lab texture data | `--texture-input FILE.csv` | Columns `sand_pct`, `silt_pct`, `clay_pct` (% of the < 2 mm fine-earth fraction). Other columns are kept. |

Ask the user:

- For the field's coordinates, not a town name. Several fields in one town can sit on very different soils.
- Whether the field is tile-drained. SSURGO's natural drainage class and water table describe undrained conditions.
- For lab data, whether the percentages are of the fine-earth fraction (USDA) and which particle-size limits were used. USDA uses sand 2–0.05 mm and silt 0.05–0.002 mm. Data with the international 0.02 mm sand/silt limit can't be classified with the USDA triangle.

## How to run

Commands use `${CLAUDE_SKILL_DIR}`, this skill's folder: Claude Code fills it in; in other tools, replace it with the folder that contains this SKILL.md.

```bash
python "${CLAUDE_SKILL_DIR}/scripts/soil.py" --demo
python "${CLAUDE_SKILL_DIR}/scripts/soil.py" --lat 40.08 --lon -88.20 --out results/soil
python "${CLAUDE_SKILL_DIR}/scripts/soil.py" --mukey 426520
python "${CLAUDE_SKILL_DIR}/scripts/soil.py" --texture-input lab_particle_size.csv --out results/texture
```

The demo is Drummer silty clay loam (mukey 242963) near Champaign, Illinois, from bundled Soil Data Access responses.

## Outputs

- `report.md`: warnings first, then the map unit, components, the dominant component's horizons, water storage and drainage, general guidance, figures and references.
- `tables/horizons.csv`: every horizon of every component. Key columns:
  - `texdesc` (SSURGO's representative texture), `texture_class` (the texture reported), `texture_source`;
  - `texture_computed` (triangle cross-check) and `texture_agrees`;
  - sand/silt/clay %, `om_pct`, `ph_h2o` and `ph_class`;
  - `awc_cm_cm`, `w33kpa_pct_vol`, `w1500kpa_pct_vol`, `ksat_um_s`, `db_g_cm3`, `cec7_meq_100g`.
- `tables/components.csv`: component name, percent, major flag, drainage class, hydrologic group, taxonomy.
- `tables/mapunit_weighted_properties.csv`: depth-weighted then component-weighted means for 0–30 and 0–100 cm, with the % of the map unit that has data.
- `tables/aws_profile.csv`: cumulative available water storage (mm) every 10 cm to 200 cm, component-weighted and for the dominant component. **This is the soil input of claw-water-balance** (`--soil-dir`).
- `tables/mapunit_summary.csv`: SSURGO muaggatt AWS 0–25/50/100/150 cm (mm), the same values recomputed, drainage class, hydrologic group, minimum water-table depth (cm), bedrock depth, slope, capability class and survey version date.
- `figures/soil_profile.png`, `raw/sda_responses.json` (online runs) and `reproducibility/`.
- With `--texture-input`: `tables/texture_classes.csv`, with the class and a note for rows that can't be classified.

## Expert decisions

- **Texture.** SSURGO's own representative texture description (`chtexturegrp.texdesc` where `rvindicator = 'Yes'`, joined on `chkey`) is the reported texture. It includes modifiers (gravelly, mucky), stratified layers and organic materials that sand/silt/clay alone can't express.
- **Fallback classifier.** The 12 USDA texture classes are coded word for word from the Soil Survey Manual (Soil Science Division Staff 2017, pp. 122–124). Examples: sands, > 85% sand and silt + 1.5 × clay < 15; sandy loams, 7 to < 20% clay with > 52% sand and silt + 2 × clay ≥ 30, or < 7% clay and < 50% silt with silt + 2 × clay ≥ 30.
  - Tests cover one point in each class, 21 boundary points, and a 0.5% grid (20,301 points) on which exactly one definition holds.
  - The old classifier's errors are fixed: 10/70/20 is silt loam and 60/35/5 is sandy loam.
  - Sums within 100 ± 2 points are rescaled to 100 (adjustable rule of thumb, `SUM_TOLERANCE`). Larger mismatches are refused.
- **Map-unit weighting.** Each property is thickness-weighted over the depth range within each component, then weighted by component percent (`comppct_r`). Components without data are excluded, and the share that has data is reported.
  - pH is averaged arithmetically, as Web Soil Survey's weighted average does. Treat averaged pH as approximate.
- **Available water storage.** AWS to depth d = Σ AWC × horizon thickness above d, component-weighted. This is how SSURGO's `muaggatt.aws0100wta` and the other `aws*wta` columns are defined ("weighted average of all components", expressed in cm).
  - The recomputation reproduces muaggatt to within 0.1 mm at 25, 50, 100 and 150 cm for both test map units. The report flags any map unit where they disagree by more than 5%.
  - AWC is SSURGO's estimate of θFC − θWP (1/10 or 1/3 bar minus 15 bar), adjusted for rock fragments and salinity, so TAW(Zr) = AWS(Zr) (FAO-56 eq. 82).
  - Storage stops at the deepest described horizon.
- **Interpretations** (general guidance, not recommendations):
  - Reaction classes follow Soil Survey Manual Table 3-15 (e.g. 6.6–7.3 neutral).
  - Drainage-class text paraphrases the Manual's definitions (pp. 213–215).
  - Organic-matter words (<1 very low, 1–2 low, 2–4 moderate, 4–6 high, ≥6 very high %) are an adjustable rule of thumb.
- **Red flags written to the report.**
  - Water table shallower than 200 cm: FAO-56 assumes capillary rise ≈ 0 only when the water table is more than about 1 m below the root zone.
  - Bedrock shallower than 150 cm.
  - Components without data, component percentages not summing to 100, and texture disagreements.

## Checks before reporting

- Read the Warnings section of `report.md` and pass every warning on.
- Always say that SSURGO is mapped survey data for a map unit, not a sample of the user's field, and give the survey version date.
- Stop and tell the user if:
  - the point falls outside SSURGO, or the coordinates look swapped;
  - the map unit is a complex with several major components of contrasting properties (one value hides real variation);
  - the water table is shallow or the soil is poorly drained, and the user plans a water balance (ask about tile drainage);
  - lab particle-size data use non-USDA size limits.
- For the water balance, pass `tables/aws_profile.csv` (or the folder) to **claw-water-balance** and check that the chosen root depth doesn't exceed the described profile depth.

## References

- Soil Science Division Staff (2017). Soil survey manual. C. Ditzler, K. Scheffe & H.C. Monger (eds.). USDA Handbook 18. Government Printing Office, Washington, D.C. Texture classes pp. 122–124; reaction classes Table 3-15; drainage classes pp. 213–215. https://www.nrcs.usda.gov/resources/guides-and-instructions/soil-survey-manual
- Soil Survey Staff, Natural Resources Conservation Service, USDA. Soil Survey Geographic (SSURGO) Database, accessed through Soil Data Access. https://sdmdataaccess.sc.egov.usda.gov/
- USDA-NRCS (2025). Soil Data Mart 4.3 Table Column Descriptions (muaggatt, chorizon, chtexturegrp). https://sdmdataaccess.nrcs.usda.gov/documents/TableColumnDescriptionsReport.pdf
- Allen, R.G., Pereira, L.S., Raes, D. & Smith, M. (1998). Crop evapotranspiration: guidelines for computing crop water requirements. FAO Irrigation and Drainage Paper 56. FAO, Rome. https://www.fao.org/4/x0490e/x0490e00.htm
