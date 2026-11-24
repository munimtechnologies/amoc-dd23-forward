# Grading of v2 (job amoc-dd23-forward-v2, ran 2026-11-24 20:26-20:44 UTC)

Graded 2026-11-24 against [PREREGISTRATION.md](PREREGISTRATION.md) and
[DEVIATION-2026-11-24.md](DEVIATION-2026-11-24.md) (written after v1 ran, before v2 ran), unchanged.
Files: `results/v2/`. Made by AI agents working in the name of MV2 (Munim, Inc.); no human expert has
checked this yet.

HadISST1 covered January 1870 to **September 2026** (1,881 months; **69 months after 2020**).

## Gates (both passed, so P leads)

| check | grade | numbers |
|---|---|---|
| Replica R: `original`+penalty on the authors' AMOC2 file | within 1.0 y of 2056.52 | **2056.80** (diff 0.28 y) |
| Replica R: `corrected`+penalty on the authors' AMOC2 file | within 3 y of 2065 | **2064.77** (diff 0.23 y) |
| Replica overall | **REPLICA_OK** | report-only: 1xGM 2123.55 (their Table 1: 2123); 3xGM 2067.96 (Table 1: 2068) |
| Fingerprint H vs their AMOC2, 1870-2020 | **HARMONISE_OK** | RMS 0.00077 K (limit 0.01 K), r = 0.9999977; affine maps SG k = 0.5843, GM k = 0.6870 |

Report-only on the replica: without the penalty the same code gives 2086.30 (corrected) and 2078.08
(original), the numbers the penalty critique of Morr et al. (arXiv:2604.20341) is about. The `mu_next`
variant of the corrected step gives 2064.71.

Our P fit to the harmonised HadISST rebuild through December 2020 gives **2064.77**, the same as on the
authors' own file. Its 1,000-draw parametric bootstrap gives a 95% range of **2036.2-2107.0** (authors'
published 2037-2109).

## Forward test, P = the authors' penalised, corrected estimator (leads)

| test | pre-registered grade | numbers |
|---|---|---|
| T1 forecast check, 2021-01..2026-09 | **WITHIN_FORECAST** | q = 0.70; observed mean -0.535 vs model 5-95%: -0.864 to -0.403 (2,000 paths; 1 path tipped) |
| T2 one-step residuals after 2020 | **CALIBRATED** | variance 1.35 (z = 2.02, limit 2.576); mean +0.04 (t = 0.31); KS p = 0.91 |
| T3 tipping-year update | **LATER** | 2064.77 (data to 2020-12) -> 2069.99 (data to 2026-09), shift **+5.22 y** (threshold +5) |

Bootstrap 95% range: 2036.2-2107.0 (data to 2020) -> 2040.9-2111.1 (data to 2026-09); medians 2062.0 -> 2065.3.

Per calendar year, the observed mean fingerprint sat at these percentiles of the model's forecast:
2021 0.87, 2022 0.81, 2023 0.63, 2024 0.63, 2025 0.39, 2026 (9 months) 0.22; the last 12 months 0.29.

Tipping year refitted with data ending each December (report-only):
2020 2064.77, 2021 2066.71, 2022 2067.07, 2023 2066.95, 2024 2071.21, 2025 2072.40, 2026-09 2069.99.

## Forward test, U = unpenalised estimator (reported alongside)

| test | grade | numbers |
|---|---|---|
| T1 | **WITHIN_FORECAST** | q = 0.67 |
| T2 | **CALIBRATED** | variance 1.36 (z = 2.09); mean +0.03 (t = 0.22); KS p = 0.93 |
| T3 | **LATER** | 2086.50 -> 2093.99, shift +7.49 y |

Bootstrap 95% (500 draws): 2043.7 to beyond 59,000 (data to 2020) -> 2047.6 to beyond 100,000 (data to
2026-09): without the penalty the upper end is unbounded in practice, as Morr et al. found.
By data end: 2086.50, 2089.10, 2089.49, 2089.68, 2096.78, 2099.61, 2093.99.

## Report-only variants (P)

| series | t_c to 2020 | t_c to 2026-09 | shift |
|---|---|---|---|
| SG - 1xGM | 2123.52 | 2136.29 | +12.77 |
| SG - 3xGM | 2068.05 | 2066.30 | -1.75 |
| SG - 2xGM, unharmonised scale | 2059.78 | 2065.25 | +5.46 |

## What it means

1. **The published forecast survived its first out-of-sample test.** With the authors' exact estimator
   (reproduced to 0.3 years), 69 months of new data stayed inside the model's forecast range, and the
   one-step residuals are consistent with the model's noise (v1's VARIANCE_HIGH came from v1's slightly
   different series; on the series that matches the authors' file the residual variance is within limits).
2. **The new data moved the date later, but only just past the threshold.** +5.2 years against a
   pre-registered +5. The path is not monotone: the estimate rose to 2072.4 with data to 2025 and fell back
   to 2070.0 when the cold-fingerprint months of 2026 were added. The 95% range moved about 4-5 years later
   and kept its width, so six years of data neither confirm nor sharpen the forecast.
3. **2025-2026 were below the model's median path** (q 0.39 and 0.22), the opposite direction from
   2021-2022 (0.87, 0.81). Because the fingerprint is SG minus twice the global mean, global warm spikes
   (such as 2023-24) push it down whatever the AMOC does; the forward test itself does not separate them.
   **Post-hoc (not pre-registered), by hand from v1's plain-mean terms in `results/v1/fingerprint_monthly.csv`:**
   the subpolar-gyre anomaly averaged 0.378 K in 2021, 0.305 K in 2025 and 0.332 K in January-September 2026,
   while the global-mean anomaly averaged 0.380 K, 0.581 K and 0.673 K. So most of the 2025-26 fall of the
   fingerprint comes from the -2xGM term (global warmth), not from a colder gyre. The v1 series is the
   plain-mean variant, not the harmonised v2 series, so treat these as approximate.
4. **The penalty matters more than the new data.** P and U differ by 22-27 years at every data end; the
   new data move each by 5-7 years. This agrees with Morr et al.'s point; it is theirs, not ours.

v1's pre-registered grades ([GRADING-v1.md](GRADING-v1.md)) stand as recorded; v2 is the test of the
published forecast because only v2 passes both gates.
