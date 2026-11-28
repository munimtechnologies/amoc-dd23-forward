# Grading of the ERSSTv5 arm (PREREGISTRATION-ERSST.md), 2026-11-28

Made with AI agents working in MV2's name (MV2, Munim, Inc.); no human expert has checked it yet.

The frozen job (`job_ersst.py`; estimator modules fetched unchanged from commit 74c83b30, copies in
`results/ersst-v0.1.0/deps/`) ran once, on 2026-11-28 05:01-05:17 UTC (compete job `amoc-dd23-ersst-v1`; unit
tests passed first). Its output is committed unchanged in `results/ersst-v0.1.0/`.

## Gate

**E0: READ_OK.** NOAA ERSSTv5 read contiguous 1870-01 to 2026-10 (1,882 months), no layout failures. The SG
region covers 208 ERSST 2-degree cells (744 SG 1-degree centres). Harmonised SG - 2 GM fingerprint against the
HadISST-based DD23 file: r = 0.918, rms 0.137 over 1870-2020.

## Grades (primary estimator P, the authors' corrected penalised fit)

| grade | result | HadISST1 v2 grade |
|---|---|---|
| E1 = T1 forecast envelope, 2021-01..2026-09 (69 months) | **WITHIN_FORECAST**, q = 0.51 (observed mean -0.566; 90% band -0.762 to -0.394) | WITHIN_FORECAST, q = 0.70 |
| E2 = T2 one-step residuals | **CALIBRATED** (mean -0.03, variance 0.83, KS p = 0.42) | CALIBRATED |
| E3 = T3 tipping-year shift | **LATER**: 2098.3 (data to 2020) -> 2105.9 (data to 2026-09), +7.7 y | LATER, 2064.77 -> 2069.99 |
| Headline | **CONSISTENT** (E1 and E3 match HadISST1) | |
| Convergence | **UNCHANGED**: gap to HadISST 33.5 y (2020 fits) -> 36.0 y (latest), +7% | |

**Reading.** On a second sea-surface-temperature record the Ditlevsen & Ditlevsen forecast also held out of
sample: the 69 months since the 2020 fit sit near the middle of its own forecast fan, its one-step residuals are
calibrated, and six more years moved the estimated tipping year later, not earlier, as on HadISST1. So the v2
out-of-sample result is not specific to HadISST1. What does depend on the record is the tipping year itself:
about 2065-2070 on HadISST1 and about 2098-2106 on ERSSTv5 with the same estimator, a gap that six new years did
not close.

Within the window the last 12 months are low (q = 0.25) and 2026 is the lowest year (9 months, mean -0.80,
q = 0.17), still inside the fan.

## Report-only items (not graded)

- P bootstraps (500 each, none tipped in simulation, none failed): tipping year 2070 (2041-2114, 95%) for the
  2020 fit and 2074 (2047-2117) for the latest fit. The point estimates (2098, 2106) lie in the upper part of
  their own bootstrap distributions.
- Unpenalised estimator U: T1 WITHIN_FORECAST (q = 0.36), T2 CALIBRATED, tipping year 2306 -> 2384.
- Original 2023 code on ERSSTv5 to 2020: 2088.0. Ben-Yami et al. (2024) report 3047.7 for ERSSTv5. **We do not
  reproduce their number**; our fingerprint is harmonised onto the DD23 file's columns (their choice may differ),
  and this difference is not resolved here.
- beta = 1: 2171 -> 2191; beta = 3: 2067 -> 2064. Unharmonised degC fingerprint: 2083 -> 2092. P on all months to
  2026-10: 2099.3.
- Tipping year (P / U) by data end, Dec 2020 to Sep 2026: `tc_by_data_end_ersst.csv`.

## Limits

One estimator family (the authors'), one fingerprint definition (SG - 2 GM), harmonised onto the DD23 columns.
69 out-of-sample months are a short test of a forecast whose distinctive part lies decades ahead; WITHIN_FORECAST
says the model is not contradicted by 2021-2026, not that it is right. ERSSTv5 is a coarser (2-degree) and more
smoothed record than HadISST1.

## Credit

Ditlevsen & Ditlevsen 2023 (Nat Commun 14:4254; 2025 correction); Ben-Yami, Morr, Bathiany & Boers 2024 (Sci Adv);
Morr et al. 2026 (arXiv:2604.20341); Caesar et al. 2018 (SG grid); Huang et al. 2017 (ERSSTv5, NOAA NCEI).
