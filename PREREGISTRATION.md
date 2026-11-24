# Pre-registration: an out-of-sample test of the Ditlevsen & Ditlevsen (2023) AMOC tipping forecast

Frozen 2026-11-24, before the analysis job has run and before we looked at any HadISST month after December 2020.
Licence: Apache-2.0 (code), CC BY 4.0 (text). Made by AI agents working in the name of MV2 (Munim, Inc.).
No human expert has reviewed this design.

## Why this question

Ditlevsen & Ditlevsen (2023, *Nature Communications* 14, 4254; "DD23") fitted a stochastic saddle-node ("fold")
model with a linearly ramping control parameter to a sea-surface-temperature fingerprint of the Atlantic
Meridional Overturning Circulation (AMOC): the subpolar-gyre SST anomaly minus twice the global-mean SST anomaly
(SG - 2xGM), HadISST1, monthly, 1870-2020. The paper's estimate of when the AMOC would tip, widely reported as
"around 2057 (95%: 2025-2095)", was revised in the authors' 2025 correction (*Nat. Commun.* 16, 7794; a bug in the
Strang-splitting estimator found by Anders Gantzhorn Kristensen) to **2065 (95%: 2037-2109)**.

The estimate is one of the most cited and most disputed numbers in climate science. Critiques (Ben-Yami et al.
2024, *Sci. Adv.*; Morr et al. 2026, arXiv:2604.20341) test its assumptions with synthetic data, other fingerprints
and other SST products, but all of them, and the original, stop at 2020-2021. **The fitted model makes a forecast
for the fingerprint after 2020, and nearly six years of those months are now public.** As of 2026-11-24 we found
no paper, preprint or report that tests the DD23 model on post-2020 data or refits it with them (searches: arXiv
listing "AMOC tipping" to 2611.23429 of 2026-11-17, Google/Scholar-style web searches, the Morr et al. 2026 full text,
which uses HadISST to 2021).

The wow, if it lands: the first out-of-sample test of the best-known forecast of an AMOC collapse, and the first
update of its tipping year with the post-2020 data, from public data and with code anyone can rerun.

## Data (all public)

- HadISST1 monthly 1-degree SST, ASCII files from https://www.metoffice.gov.uk/hadobs/hadisst/data/
  (1870-1900 ... 1991-2003 period files, yearly files 2004-2026, plus the one-month update file if it is newer).
  Rayner et al. 2003, *J. Geophys. Res.*
- Subpolar-gyre cells: Caesar et al. 2018 (*Nature*) `sg_grid.txt`
  (https://www.pik-potsdam.de/~caesar/AMOC_slowdown/), the region DD23 says they used.
- DD23's own fingerprint file `AMOCdata.txt` (ERDA, doi:10.17894/ucph.9ef73d8c-b642-4b55-88bc-b38984d043b9,
  CC BY 4.0): time, AMOC0 (SG anomaly), AMOC1, AMOC2, GM, January 1870 to December 2020.

## Method (frozen)

1. **Rebuild.** For every month, the SG mean over the 744 Caesar cells and the global ocean mean, each in four
   variants: area-weighted (cos latitude) or plain, 100%-sea-ice cells (-1000) excluded or set to -1.8 degC.
   Anomalies: subtract each calendar month's 1870-2020 mean (DD23's seasonal correction, base fixed at 1870-2020
   so the months through 2020 are identical in both fits).
2. **Variant choice (uses only 1870-2020).** The SG variant with the smallest RMS difference from DD23's AMOC0, and
   the GM variant with the smallest RMS difference from DD23's GM.
3. **Gate G (reconstruction):** correlation r and RMS of our SG-2xGM against DD23's AMOC2, 1870-2020.
   RECON_OK if r >= 0.98 and RMS <= 0.05 K; RECON_WEAK if r >= 0.95; else RECON_FAIL.
4. **Model and estimator** (`ddfit.py`). DD23 eqs. 1-3, t0 = 1924 fixed, time step one month. Baseline (t < t0):
   exact conditional OU maximum likelihood gives mu0, alpha0, sigma2. After t0: A and tau by a Strang-splitting
   pseudo-likelihood (nonlinear part dy/dt = -A y^2 solved exactly; half-flow and OU step at t_{k-1}, final
   half-flow inverted at t_k, as the 2025 correction requires). Maximised by a 25x25 grid in (log A, log tau)
   then Nelder-Mead. Tipping year t_c = t0 + tau. Confidence intervals: DD23-style parametric bootstrap
   (simulate the fitted model on the same monthly grid with Euler-Maruyama, 10 sub-steps; refit; percentiles).
5. **Control C (reproduction):** our estimator on DD23's own AMOC2 file must give their corrected 2065.
   REPRO_PASS if |t_c - 2065| <= 5 years; REPRO_NEAR if <= 10; else REPRO_FAIL. 300 bootstrap draws; their
   95% interval 2037-2109 is reported next to ours (not graded).
6. **M20:** fit on our rebuilt series January 1870 - December 2020 (1000 bootstrap draws).
   **Mlatest:** fit on January 1870 - the last HadISST month available when the job runs (1000 draws).

## Tests and grades (frozen)

- **T1, forecast check (primary).** From the observed December 2020 value, simulate 2000 paths of the M20 model
  over every month from January 2021 to the last month. q = share of paths whose mean over those months lies
  below the observed mean (paths that tip count as below). BELOW_FORECAST if q < 0.05 (the fingerprint fell faster
  than the model allowed), ABOVE_FORECAST if q > 0.95 (it held up or recovered beyond the model's range),
  WITHIN_FORECAST otherwise. Per-calendar-year and last-12-month percentiles are reported, not graded.
- **T2, one-step calibration (secondary).** Standardised one-step Strang residuals of M20 on the post-2020
  months. Flags BIAS_HIGH/LOW (|mean x sqrt(n)| > 2.576), VARIANCE_HIGH/LOW (|z of variance| > 2.576), KS_REJECT
  (KS p < 0.01); CALIBRATED if none.
- **T3, update of the tipping year (primary).** Shift = t_c(Mlatest) - t_c(M20). EARLIER if <= -5 years, LATER if
  >= +5 years, STABLE otherwise; UNBOUNDED if either fit runs to tau's upper bound (1e5 years).
- **Report only:** SG-1xGM and SG-3xGM (DD23 Table 1 robustness), the full-record climatology variant, and the
  1xGM control on DD23's file (their table: 2123).

How grades combine: if G is RECON_FAIL, T1-T3 are reported but carry no weight and the article says so first.
If C is REPRO_FAIL, the article leads with that: our shift is then the shift of *our* implementation of their
method, not a revision of their number. T1 and T3 are reported together whatever they say; a STABLE or
WITHIN_FORECAST result is reported as plainly as a dramatic one.

## What this does not test

It does not test whether the AMOC is weakening, whether SST fingerprints track it, or whether a fold
bifurcation is the right model (the points argued by Ben-Yami et al. 2024, Morr et al. 2026 and Emirzade et al.
2026, arXiv:2609.18551). It asks a narrower question that nobody has asked yet: did the DD23 model, as fitted,
anticipate the six years that followed, and where does its own estimate move when they are added? Six years
is short next to a century-long ramp, the 2023-24 global SST spike enters the fingerprint through the 2xGM term,
and HadISST1's post-2020 months use the operational input stream. These caveats go into any write-up.

## Deviations

Any change after the job has run is written in a dated DEVIATION file next to this one, with the reason, and the
pre-registered result is reported alongside.

## Credit

Peter D. Ditlevsen and Susanne Ditlevsen (DD23 and its public code and data); Anders Gantzhorn Kristensen (the
2025 estimator fix); Levke Caesar, Stefan Rahmstorf and co-authors (the fingerprint and SG grid); the Met Office
Hadley Centre (HadISST1); Maya Ben-Yami, Andreas Morr, Niklas Boers and co-authors (the critiques this test
complements).
