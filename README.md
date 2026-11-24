# amoc-dd23-forward

An open, pre-registered **out-of-sample test of the Ditlevsen & Ditlevsen (2023) forecast of an Atlantic
overturning (AMOC) collapse**. The forecast puts the tipping year at 2065 (95%: 2037-2109) after the authors'
2025 correction; it was first reported as about 2057. It was fitted to sea-surface temperatures from 1870 to
2020. This project asks two questions nobody has published answers to:

1. Did the fitted model anticipate what the fingerprint did from January 2021 to the latest month? The
   fingerprint is subpolar-gyre SST minus twice global-mean SST, from HadISST1. (T1, T2)
2. Where does the model's own tipping-year estimate move when those months are added? (T3)

The design, the thresholds and every grade label were frozen in [PREREGISTRATION.md](PREREGISTRATION.md) on
2026-11-24, before the job ran.

## Result (v2, graded 2026-11-24): [GRADING-v2.md](GRADING-v2.md), files in `results/v2/`

- **Gates passed.** A line-by-line replica of the authors' penalised R estimator (`ddrep.py`) reproduces
  their published numbers from their own file: 2056.80 (published 2056.52) and 2064.77 (published 2065),
  plus Table 1's 2123 and 2068. Our HadISST rebuild matches their series to 0.0008 K RMS.
- **T1 WITHIN_FORECAST:** with data to September 2026 (69 new months) the fingerprint stayed inside the
  model's forecast range (q = 0.70). 2021-22 ran above the model's median path; 2025-26 ran below it.
- **T2 CALIBRATED:** the post-2020 one-step residuals are consistent with the model's noise.
- **T3 LATER:** the refitted tipping year moved from 2064.8 to 2070.0 (+5.2 years; threshold +5). The 95%
  range moved from 2036-2107 to 2041-2111.
- The unpenalised estimator (reported alongside) moves from 2086.5 to 2094.0; the penalty changes the date
  far more than six years of data do, as Morr et al. (2026) found for the penalty.

v1 (unpenalised, a slightly different SST variant) is graded in [GRADING-v1.md](GRADING-v1.md); its
reproduction control failed, which is why v2 was declared in
[DEVIATION-2026-11-24.md](DEVIATION-2026-11-24.md) before v2 ran.

## Files

| file | what it does |
|---|---|
| `ddfit.py` | the DD23 model: OU baseline MLE, Strang-splitting likelihood (unpenalised), Nelder-Mead fit, simulation, bootstrap, forecast checks (standard library only) |
| `ddrep.py` | v2: replica of the authors' R objective (penalty, a floor, both Strang variants), its fit and bootstrap, affine harmonisation, v2 grades |
| `hadisst.py` | HadISST1 ASCII reader, Caesar et al. subpolar-gyre cells, monthly means, anomalies, DD23 file parser |
| `grade.py` | the frozen v1 grading rules |
| `job.py` / `job_v2.py` | the v1 and v2 jobs (download public data, fit, test, write results) |
| `test_*.py` | unit tests with hand-checked values (`python -m unittest`) |
| `results/v1/`, `results/v2/` | outputs of the two runs |

Run: `python -m unittest -v && python job_v2.py --data-dir data --out-dir out` (Python 3.10+, no packages
needed; about 1 GB of downloads).

## Credit and data

- P. D. Ditlevsen & S. Ditlevsen, *Nature Communications* 14, 4254 (2023), and Author Correction (2025). Code
  and data: doi:10.17894/ucph.9ef73d8c-b642-4b55-88bc-b38984d043b9 (CC BY 4.0). The estimator fix is credited
  there to Anders Gantzhorn Kristensen.
- A. Morr, M. Ben-Yami, B. Groenke, C. Schötz, A. Cotronei, E. Myrvoll-Nilsen, S. Bathiany, M. Rypdal &
  N. Boers, arXiv:2604.20341 (2026). Their analysis of the penalty and of the method's coverage is the context
  for v2.
- L. Caesar et al., *Nature* (2018): the SST fingerprint and the subpolar-gyre grid (downloaded at run time from
  PIK, not redistributed here).
- HadISST1, Met Office Hadley Centre (Rayner et al. 2003).
- Ben-Yami et al. 2024 (*Science Advances*).

This is research software, not a climate projection service. Made with AI agents working in the name of MV2
(Munim, Inc., mv2@munimtech.com); no human expert has checked it yet. Licence: Apache-2.0.
