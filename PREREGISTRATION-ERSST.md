# Pre-registration: the AMOC forward test on a second SST record (NOAA ERSSTv5)

Study: amoc-dd23-forward. Arm: ERSSTv5. Frozen 2026-11-24 (UTC), before this study read any ERSSTv5
data. Licence: Apache-2.0 (code), CC BY 4.0 (text). Made with AI agents working in the name of MV2
of Munim, Inc.; no human expert has checked it yet.

## Why

The v2 grading (GRADING-v2.md) found, on HadISST1 with the authors' own penalised estimator, that
the Ditlevsen & Ditlevsen (2023, corrected 2025) forecast held out of sample (T1 WITHIN_FORECAST,
q = 0.70; T2 CALIBRATED) and that 69 new months moved the tipping year later (T3 LATER,
2064.77 -> 2069.99). One SST product is one choice. Ben-Yami et al. (2024, Science Advances,
Table S2) and Morr et al. (2026, arXiv:2604.20341) showed that the same fingerprint built from
ERSSTv5 gives a very different tipping year on data to 2020 (Ben-Yami et al.: 3047.7). Nobody has
run the out-of-sample test on ERSSTv5. This arm asks: **is the v2 out-of-sample result specific to
HadISST1?**

## Data and fingerprint (fixed)

- NOAA ERSSTv5 yearly ASCII files, https://www.ncei.noaa.gov/pub/data/cmb/ersst/v5/ascii/
  (2-degree grid, 0.01 degC, -9999 missing; layout from the NCEI Readme and read.ersst.f).
- SG region: Caesar et al. (2018) sg_grid.txt 1-degree cells; each 2-degree ERSST cell is weighted
  by the number of SG 1-degree centres inside it times cos(latitude). GM: cos(latitude)-weighted
  mean of all ERSST cells with data. Monthly anomalies against 1870-2020 calendar-month means.
- Harmonised exactly as v2: SG and GM anomalies are each mapped onto the DD23 file's AMOC0 and GM
  columns by a least-squares affine fit over the overlap (1870-2020); fingerprint = SG - 2 x GM.
- Primary window: January 1870 to September 2026 (the v2 HadISST1 window). "2020 fit" = data to
  December 2020.
- Estimator: the v2 modules ddfit.py, ddrep.py, grade.py, hadisst.py, fetched unchanged at run time
  from commit 74c83b30f81c4b32ee7ad03d1a8460ca05550aa9 (fetch_deps.py). Primary = P, the authors'
  corrected penalised estimator (pen 0.004).

## Gate

- **E0 read:** READ_OK if the record is contiguous from 1870-01 to at least 2026-06 and at least 99%
  of months pass the fixed layout check (four open-ocean cells hold data, four continental cells
  are missing, equatorial Pacific 15-35 degC). Otherwise READ_FAIL and nothing is graded.

## Grades (same rules as v2, in grade.py and ersst.py)

- **E1 = T1** forecast envelope: from the 2020 fit's state at December 2020, 2,000 simulated paths;
  q = share of paths whose 2021-01..2026-09 mean is below the observed mean.
  BELOW_FORECAST q < 0.05, ABOVE_FORECAST q > 0.95, else WITHIN_FORECAST.
- **E2 = T2** one-step residual flags (BIAS, VARIANCE, KS at 2.576 / p 0.01) or CALIBRATED.
- **E3 = T3** tipping-year shift, latest fit minus 2020 fit: <= -5 y EARLIER, >= +5 y LATER, else
  STABLE; UNBOUNDED if either fit sits at tau's upper bound.
- **Headline (ersst.grade_consistency):** CONSISTENT if E1 and E3 equal the HadISST1 v2 grades
  (WITHIN_FORECAST, LATER); T1_ONLY if only E1 matches; DIFFERENT if E1 differs.
- **Convergence (ersst.grade_convergence):** with gap = |tc(ERSST) - tc(HadISST v2)| for the 2020
  fits (HadISST 2064.77) and for the latest fits (HadISST 2069.99): CONVERGING if the gap shrinks
  by more than 10%, DIVERGING if it grows by more than 10%, else UNCHANGED.

Expectation written down before running: because ERSSTv5's 2020 tipping year is far from
HadISST1's (Ben-Yami et al. 3047.7), T3 on ERSSTv5 measures little about the mid-century forecast.
E1/E2 are the informative grades, and the convergence grade says whether six new years pull the
two records together.

## Report only (not graded)

Unpenalised estimator U (forward test and fits); the 2023 ("original") code on ERSST data to 2020
next to Ben-Yami et al.'s 3047.7; P on all available months; unharmonised (degC) fingerprint;
beta = 1 and 3; tipping year by data end 2020..2026; P parametric bootstraps (500 each, skipped if
tau is at its bound); annual raw SG and GM anomalies (degC), which give an independent check of the
post-hoc v2 note that the 2025-2026 fall of the fingerprint comes mostly from the global-mean term.

## Outputs

results_ersst.json (point results in results_ersst_point.json), fingerprint_ersst_monthly.csv,
fingerprint_ersst_annual.csv, forecast_fan_ersst_P.csv, tc_by_data_end_ersst.csv, deps/ (the exact
estimator modules used), inputs/.

## Credits

Ditlevsen & Ditlevsen 2023 (Nat Commun 14:4254; 2025 correction); Ben-Yami, Morr, Bathiany & Boers
2024 (Sci Adv); Morr et al. 2026 (arXiv:2604.20341); Caesar et al. 2018 (SG grid); Huang et al. 2017
(ERSSTv5).
