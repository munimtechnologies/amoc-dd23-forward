# Deviation 2026-11-24: v2 uses the authors' own (penalised) estimator and the matching SST variant

Written 2026-11-24 after v1 had run and been graded ([GRADING-v1.md](GRADING-v1.md)), before v2 has run.
The v1 pre-registered grades stand and will be reported next to v2. Made by AI agents working in the name of
MV2 (Munim, Inc.); no human expert has reviewed this design.

## Why

v1's control failed (2086.5 vs the authors' 2065) because v1's likelihood has no penalty, while the authors'
code adds `pen*n*(1/a - 1)` for a < 1 with pen = 0.004 (analysed by Morr et al., arXiv:2604.20341). v1's gate
was weak because its variant rule chose by RMS, though the area-weighted, ice-at-−1.8 C variant matches the
authors' file to r = 0.999999. The question "did the DD23 forecast hold up after 2020?" needs the authors'
estimator on the authors' series. v1 answered it for a different estimator and a slightly different series.

## What v2 does (frozen now)

- **Replica R** (`ddrep.py`): a line-by-line Python copy of the authors' R `loglik.fun`. That means the
  -2 x log-likelihood objective, the penalty, the floor a >= 0.1, baseline data time <= 1924 and ramp data
  time > 1924, Nelder-Mead from (tau, a) = (100, 1), plus a 15 x 15 grid restart (the better optimum is kept and
  both are reported). Variants: `original` (the 2023 final step, with the bug) and `corrected` (the 2025 fix:
  the final half-flow uses the observation at t_i).
  - **REPLICA_OK** if `original`+penalty on their AMOC2 file is within 1.0 year of 2056.52 (printed by their
    2023 code) AND `corrected`+penalty is within 3 years of 2065. **REPLICA_NEAR** if within 3 and 6 years.
    Otherwise **REPLICA_FAIL**.
- **Fingerprint H:** area-weighted SG and GM means with 100%-ice cells at -1.8 C, anomalies against the
  1870-2020 monthly means. SG and GM are each mapped onto the authors' scale by a least-squares affine fit to
  their AMOC0 and GM columns over 1870-2020. The map is fixed by 1870-2020 only and then applied to every
  month. Fingerprint = SG_h - 2 GM_h.
  - **HARMONISE_OK** if RMS vs their AMOC2 <= 0.01 K, **HARMONISE_NEAR** if <= 0.03 K, otherwise
    **HARMONISE_FAIL**.
- **Forward test F:** T1 (forecast check, 2000 paths from the observed December 2020 value), T2 (one-step
  residual flags) and T3 (shift of t_c, thresholds +-5 years), all exactly as pre-registered. They run for two
  estimators on the harmonised series, through 2020-12 and through the last month:
  - **P** = the authors' penalised, corrected estimator (bootstrap 1000 draws);
  - **U** = v1's unpenalised estimator (bootstrap 500 draws).
- **Which result leads:** P leads if REPLICA_OK (or NEAR) and HARMONISE_OK (or NEAR), because only then is it a
  test of the published forecast. Otherwise the write-up leads with the failed check, as in v1. U is reported
  alongside in every case, because Morr et al. argue the unpenalised estimate is the more defensible one.
- **Report only:** P with SG-1xGM and SG-3xGM; P on the unharmonised series (the penalty depends on the units
  of the series); the 1xGM and 3xGM replicas on the authors' file (their Table 1: 2123, 2068); the P and U
  tipping year refitted with data ending each December 2020-2025 and at the last month.

## Unchanged

Data sources, t0 = 1924, the monthly time step, the T1/T2/T3 thresholds, and the caveats in the
pre-registration (six years is short; the 2023-24 global SST spike enters through the 2xGM term; HadISST1's
operational input stream after 2020).
