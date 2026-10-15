---
name: claw-fluorescence
description: PAM chlorophyll-fluorescence parameters from a CSV or LI-COR LI-6800 fluorometer files - Fv/Fm (maximum quantum efficiency of PSII), PhiPSII (PhiPS2, Y(II), operating efficiency), Fv'/Fm', qP, qL, NPQ, estimated Fo' (Oxborough-Baker) and electron transport rate (ETR) - with quality-control flags for dark adaptation, implausible Fv/Fm, Fm' greater than Fm, impossible signal orders and stale dark-adapted values reused by the instrument. Use when the user has Fo, Fm, Fs and Fm' (or an LI-6800 file with fluorometer data) and asks for Fv/Fm, photoinhibition or stress screening, PSII efficiency, NPQ, qL, qP, ETR, photoprotection, or wants to check fluorescence data quality.
license: MIT
compatibility: R >= 4.2. Reading LI-6800 files needs the PhotoGEA package (>= 1.4); CSV input needs only base R; openssl optional for SHA-256 checksums. No internet needed.
metadata:
  version: "1.0.0"
  domain: physiology
  language: R
  status: tested
  summary: "PAM chlorophyll fluorescence: Fv/Fm, PhiPSII, qP, qL, NPQ, ETR"
  repository: https://github.com/bgtamang/AgriClaw
---

# PAM chlorophyll-fluorescence parameters

## What it does

Computes the standard saturating-pulse fluorescence parameters row by row from the four basic signals (dark-adapted Fo and Fm; light-adapted Fs and Fm', plus Fo' if measured), summarises them by group, and flags rows that break the definitions or look implausible. LI-6800 files are read with PhotoGEA, and the results are checked against the instrument's own values.

## When to use / when not to use

- Use for PAM measurements from any fluorometer that exports Fo, Fm, Fs and Fm' (LI-6800, Walz PAM instruments, other modulated fluorometers), or to audit an LI-6800 file's fluorescence columns.
- Not for gas exchange or A-Ci curves: use `claw-licor-6800` and `claw-aci-fit`. For fitting light-response curves of ETR or PhiPSII use `claw-light-response` on this skill's output.
- Not for fast OJIP / JIP-test induction curves, per-pixel fluorescence images (use exported leaf means), or remotely sensed solar-induced fluorescence.
- Instruments that export only derived values (for example Phi2 without Fm') cannot be recomputed or checked here.

## Inputs

**CSV** (one row per light-adapted flash): `fs` and `fm_prime` (required); `fo` and `fm` (dark-adapted, needed for Fv/Fm, NPQ and an estimated Fo'); optional `fo_prime` (measured with far-red), `ppfd_umol_m2_s` (incident PPFD, needed for ETR), `absorptance` (per-row measured leaf absorptance), and any grouping columns (`--group genotype,treatment`). Repeat the leaf's Fo and Fm on each of its rows.

**LI-6800** (`.xlsx` or plaintext; `--format licor`, detected automatically). Column names verified against the example files shipped with PhotoGEA:

| Meaning | LI-6800 column |
|---|---|
| dark-adapted Fo, Fm | `Fo`, `Fm` |
| light-adapted Fs, Fm', measured Fo' | `Fs`, `Fm'`, `Fo'` |
| incident PPFD; absorbed PPFD of the fluorometer light; PSII fraction | `Qin`; `Qabs_fs`; `PS2/1` |
| instrument results compared | `Fv/Fm`, `PhiPS2`, `Fv'/Fm'`, `qP`, `qL`, `NPQ`, `ETR`, `alt. Fo'` |
| flash IDs with timestamps | `DarkAdaptedID`, `LightAdaptedID`, `DarkPulseID` |

**Ask the user before running** (never guess):

1. Were Fo and Fm measured on dark-adapted leaves, for how long, and on the same leaf as the light-adapted flashes? Pass `--dark-adapted yes` or `no`; without it the report warns that Fv/Fm and NPQ are unconfirmed. Overnight or >= 20-30 min dark adaptation is common practice (rule of thumb; species and prior light history matter).
2. Was Fo' measured (far-red dark pulse) or should it be estimated?
3. Leaf absorptance: measured (per row or one value) or the default 0.84? Is the PSII share 0.5 reasonable (it can differ in C4 leaves and under far-red-rich light)?
4. Which columns identify genotype, treatment or plant for the summary, and should PPFD levels be summarised separately (`--by-ppfd`)?

## How to run

Commands use `${CLAUDE_SKILL_DIR}`, this skill's folder: Claude Code fills it in; in other tools, replace it with the folder that contains this SKILL.md.

```bash
Rscript "${CLAUDE_SKILL_DIR}/scripts/fluorescence.R" --demo --out demo_fluorescence
Rscript "${CLAUDE_SKILL_DIR}/scripts/fluorescence.R" --input fluor.csv --group genotype --dark-adapted yes --out results/
Rscript "${CLAUDE_SKILL_DIR}/scripts/fluorescence.R" --input leaf1.xlsx --format licor --dark-adapted yes --out results/
Rscript "${CLAUDE_SKILL_DIR}/scripts/fluorescence.R" --input fluor.csv --group genotype --by-ppfd \
  --absorptance 0.86 --exclude-flagged --out results/
```

Options: `--absorptance`, `--f-psii` (override the defaults or the LI-6800's own values), `--fvfm-high` (0.86), `--fvfm-low` (0.75), `--exclude-flagged` (leave flagged rows out of the summary).

## Outputs

- `report.md`: warnings first (dark adaptation, flags, stale dark-adapted values), formulas and sources, group summary, agreement with the LI-6800's values (LI-6800 input), figures, limitations.
- `tables/fluorescence_parameters.csv`: one row per input row: inputs, `fv_fm`, `phi_psii`, `fo_prime_used` and `fo_prime_source` (`measured` or `estimated`), `fo_prime_est`, `fvp_fmp`, `qp`, `ql`, `npq`, `qabs_umol_m2_s`, `etr_umol_m2_s`; for LI-6800 files the instrument's values as `licor_*` and `dark_to_light_h`; `qc_invalid` (row excluded from summaries) and `qc_check` (row kept, must be checked).
- `tables/group_summary.csv`: per group, `n_rows`, `n_used`, `n_invalid`, `n_flagged`, and mean, SE and n of each parameter.
- `figures/fv_fm.png` (one point per dark-adapted measurement; dashed line 0.83, dotted lines the flag limits) and `figures/phipsii_etr_npq.png` (when PPFD varies).
- `reproducibility/`: `commands.sh`, `sessionInfo.txt`, checksums.

## Expert decisions

- **Formulas**: Fv/Fm = (Fm - Fo)/Fm (Kitajima & Butler 1975); PhiPSII = (Fm' - Fs)/Fm' (Genty et al. 1989); Fv'/Fm' = (Fm' - Fo')/Fm'; qP = (Fm' - Fs)/(Fm' - Fo'); qL = qP x Fo'/Fs (Kramer et al. 2004); NPQ = (Fm - Fm')/Fm' (Bilger & Bjorkman 1990). Reviews: Maxwell & Johnson (2000), Baker (2008), Murchie & Lawson (2013). These are exactly the LI-6800's definitions (test 3).
- **qL rather than qP** as the fraction of open PSII centres for land plants: qP assumes isolated antennae, qL connected ones (Kramer et al. 2004). Both are reported.
- **qN is not reported**: published definitions differ (the LI-6800's qN is (Fm - Fm')/(Fm - Fo')); use NPQ.
- **Fo'**: a measured value is used when present (> 0); otherwise Fo' = Fo / (Fv/Fm + Fo/Fm') (Oxborough & Baker 1997), which needs dark-adapted Fo and Fm from the same leaf. The source is recorded per row.
- **ETR** = PhiPSII x absorbed PPFD x f_PSII. Defaults: absorptance 0.84, a conventional standard value (Stemke & Santiago 2011, who measured 0.55-0.92 across 58 species), and f_PSII = 0.5, the usual equal-partitioning assumption (also the LI-6800 default `PS2/1`). Both are adjustable. For LI-6800 files the instrument's absorbed PPFD (`Qabs_fs`, from its own absorptance for each LED colour) and `PS2/1` are used unless `--absorptance` or `--f-psii` is given.
- **QC thresholds** (rules of thumb, adjustable): Fv/Fm > 0.86 is flagged as implausibly high (unstressed C3 leaves average about 0.83: 0.832 +/- 0.004 over 37 species in Bjorkman & Demmig 1987, measured at 77 K); Fv/Fm < 0.75 is flagged (photoinhibition, stress, or incomplete dark adaptation). Fm' > Fm (negative NPQ) and Fs < Fo' (qP > 1) are flagged. Fo >= Fm, Fs >= Fm', Fo' >= Fm' and non-positive signals are invalid and excluded from summaries.
- **Stale dark-adapted values** (LI-6800): the instrument keeps the last Fo/Fm and writes them on every later row. A `DarkAdaptedID` timestamp more than 12 h from the `LightAdaptedID` timestamp (rule of thumb) is flagged: Fv/Fm, NPQ and the estimated Fo' of that row describe another leaf or session. The PhotoGEA example files show this (dark-adapted values from 24 July reused on 4 and 23 August, with Fv/Fm 0.925).
- **Fm' on the LI-6800** comes from a multiphase flash that extrapolates to infinite intensity (Loriaux et al. 2013); dark-adapted Fm measured with a different flash type can then exceed or fall short of Fm' for instrumental reasons. Ask whether both used the same flash settings when Fm' > Fm.
- Summaries include flagged rows unless `--exclude-flagged`; they always exclude invalid rows.

## Checks before reporting

Stop and tell the user when:

- dark adaptation is unconfirmed, or any row is flagged "dark-adapted Fo/Fm from another session": do not report Fv/Fm or NPQ for those rows;
- Fv/Fm is above 0.86 (Fo underestimated: measuring-light or background problems) or many leaves are below 0.75 without a known stress;
- Fm' > Fm: NPQ is invalid for that leaf; ask whether Fm was measured after full dark adaptation, on the same leaf, with the same flash settings;
- Fo' was estimated for many rows and qL or qP drive the conclusion: estimated Fo' assumes the same quenching of Fo and Fm;
- ETR is compared between species, leaf types or light sources: absorptance and f_PSII differ; ask for measured absorptance;
- for LI-6800 input, the report shows any difference from the instrument's values above 1e-6: constants were changed or the file was recomputed.

## Validation

`tests/test_fluorescence.R` (124 checks): a hand-worked example (Fo 300, Fm 1500, Fs 600, Fm' 1000, Fo' 250, PPFD 1000 gives Fv/Fm 0.8, PhiPSII 0.4, Fv'/Fm' 0.75, qP 0.533, qL 0.222, NPQ 0.5, ETR 168; with Fo' estimated 272.73, qP 0.55, qL 0.25); identities on 2000 random valid rows (PhiPSII = qP x Fv'/Fm'; qL <= qP when Fs >= Fo'; the estimated Fo' satisfies 1/Fo' - 1/Fm' = 1/Fo - 1/Fm; Fm' = Fm gives Fo' = Fo and NPQ = 0); all 424 rows of the eight LI-6800 example files in PhotoGEA, where Fv/Fm, PhiPSII, Fv'/Fm', qP, qL, NPQ, ETR and the Oxborough-Baker Fo' match the instrument's values (largest relative difference 1e-15); QC flags on constructed rows and on the planted problems in the demo; the `--dark-adapted no` and light-adapted-only paths; and end-to-end runs on the demo and on an LI-6800 file. Not validated: instruments other than the LI-6800 (only the formulas, which are instrument-independent), and the plausibility thresholds, which are rules of thumb.

## References

- Baker NR (2008) Chlorophyll fluorescence: a probe of photosynthesis in vivo. Annual Review of Plant Biology 59:89-113. https://doi.org/10.1146/annurev.arplant.59.032607.092759
- Bilger W, Bjorkman O (1990) Role of the xanthophyll cycle in photoprotection elucidated by measurements of light-induced absorbance changes, fluorescence and photosynthesis in leaves of Hedera canariensis. Photosynthesis Research 25:173-185. https://doi.org/10.1007/BF00033159
- Bjorkman O, Demmig B (1987) Photon yield of O2 evolution and chlorophyll fluorescence characteristics at 77 K among vascular plants of diverse origins. Planta 170:489-504. https://doi.org/10.1007/BF00402983
- Genty B, Briantais J-M, Baker NR (1989) The relationship between the quantum yield of photosynthetic electron transport and quenching of chlorophyll fluorescence. Biochimica et Biophysica Acta 990:87-92. https://doi.org/10.1016/S0304-4165(89)80016-9
- Kitajima M, Butler WL (1975) Quenching of chlorophyll fluorescence and primary photochemistry in chloroplasts by dibromothymoquinone. Biochimica et Biophysica Acta 376:105-115. https://doi.org/10.1016/0005-2728(75)90209-1
- Kramer DM, Johnson G, Kiirats O, Edwards GE (2004) New fluorescence parameters for the determination of QA redox state and excitation energy fluxes. Photosynthesis Research 79:209-218. https://doi.org/10.1023/B:PRES.0000015391.99477.0d
- Lochocki EB, Salesse-Smith CE, McGrath JM (2025) PhotoGEA: an R package for closer fitting of photosynthetic gas exchange data with non-Gaussian confidence interval estimation. Plant, Cell & Environment 48:5104-5119. https://doi.org/10.1111/pce.15501
- Loriaux SD, Avenson TJ, Welles JM, McDermitt DK, Eckles RD, Riensche B, Genty B (2013) Closing in on maximum yield of chlorophyll fluorescence using a single multiphase flash of sub-saturating intensity. Plant, Cell & Environment 36:1755-1770. https://doi.org/10.1111/pce.12115
- Maxwell K, Johnson GN (2000) Chlorophyll fluorescence: a practical guide. Journal of Experimental Botany 51:659-668. https://doi.org/10.1093/jexbot/51.345.659
- Murchie EH, Lawson T (2013) Chlorophyll fluorescence analysis: a guide to good practice and understanding some new applications. Journal of Experimental Botany 64:3983-3998. https://doi.org/10.1093/jxb/ert208
- Oxborough K, Baker NR (1997) Resolving chlorophyll a fluorescence images of photosynthetic efficiency into photochemical and non-photochemical components: calculation of qP and Fv'/Fm' without measuring Fo'. Photosynthesis Research 54:135-142. https://doi.org/10.1023/A:1005936823310
- Stemke JA, Santiago LS (2011) Consequences of light absorptance in calculating electron transport rate of desert and succulent plants. Photosynthetica 49:195-200. https://doi.org/10.1007/s11099-011-0026-y
