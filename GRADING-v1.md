# Grading of v1 (job amoc-dd23-forward-v1, ran 2026-11-24, finished 18:32 UTC)

Graded 2026-11-24 against [PREREGISTRATION.md](PREREGISTRATION.md), unchanged. Files: `results/v1/`.
Made by AI agents working in the name of MV2 (Munim, Inc.); no human expert has checked this.

HadISST1 covered January 1870 to **September 2026** (1,881 months; **69 months after 2020**).

| test | pre-registered grade | numbers |
|---|---|---|
| Gate G (rebuild vs the authors' AMOC2, 1870-2020) | **RECON_WEAK** | r = 0.982, RMS = 0.167 K (variant chosen by the frozen rule: plain mean, ice at -1.8 C, for SG and GM) |
| Control C (our estimator on the authors' own file) | **REPRO_FAIL** | t_c = 2086.5 (authors: 2065); bootstrap 95% 2040.5-5239 (300 draws) |
| T1 forecast check, 2021-01..2026-09 | **WITHIN_FORECAST** | q = 0.802; observed mean -0.587 vs model 5-95%: -1.119 to -0.447 |
| T2 one-step residuals after 2020 | **VARIANCE_HIGH** | variance 1.49 (z = 2.84); mean +0.08 (t = 0.67); KS p = 0.51. In-sample 1924-2020: variance 0.99 |
| T3 tipping-year update | **LATER** | 2073.8 (data to 2020) -> 2081.1 (data to 2026-09), shift +7.2 years |

Report-only: SG-1xGM shift +27.2 years (2187.1 -> 2214.3); SG-3xGM shift -0.3 years (2080.7 -> 2080.5);
full-record climatology 2081.7; the 1xGM control on the authors' file gives 2123.9 (their Table 1: 2123).

Per calendar year, the observed mean fingerprint sat at these percentiles of the model's forecast: 2021 0.86,
2022 0.79, 2023 0.77, 2024 0.79, 2025 0.40, 2026 (9 months) 0.29; the last 12 months 0.35.

## What the two failed checks mean (diagnosed after grading)

1. **Control C failed because the authors' estimator is penalised and ours (v1) is not.** Their public R code
   (`AMOCestimation.html`, doi:10.17894/ucph.b8f99b67-d4e6-4a2e-b518-00bddeed323b) adds `pen*n*(1/a - 1)` for
   a < 1 with pen = 0.004 to -2 x the log-likelihood. Our v1 fit on their file has a = 0.713, where the penalty
   is active; on their 1xGM series (where a > 1 and the penalty is inactive) we reproduce their Table 1 value
   (2123.9 vs 2123). Morr et al. (arXiv:2604.20341, 2026) had already analysed this penalty in detail and found
   that the unpenalised bootstrap gives 2042-5108; our unpenalised control bootstrap gives 2040.5-5239. So v1
   matches the unpenalised version of the method that Morr et al. examined, not the authors' published number.
   This is not a new finding about the original paper; it tells us which estimator v1 actually tested.
2. **Gate G was weak because the frozen variant rule chose by RMS.** The area-weighted, ice-at-−1.8 C
   variant matches the authors' file almost exactly in shape (SG r = 0.999999, GM r = 0.999997) but differs in
   scale or offset, so the RMS rule picked the plain-mean variant instead. Without the penalty, the tipping year
   does not change under an affine rescaling of the series, but the variant does change the series' shape.
   That is why v1's 2020 fit (2073.8) differs from its control (2086.5).

Both were fixed in v2 ([DEVIATION-2026-11-24.md](DEVIATION-2026-11-24.md)), written after v1 ran and before v2
ran. The v1 grades above stand as the pre-registered result and are reported next to v2.
