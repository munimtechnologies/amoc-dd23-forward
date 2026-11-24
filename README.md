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

## Status

- **v1 (pre-registered) has run and is graded:** [GRADING-v1.md](GRADING-v1.md), files in `results/v1/`.
  - With data to September 2026 (69 new months), the fingerprint stayed inside the model's forecast range
    (T1 WITHIN_FORECAST, q = 0.80).
  - The post-2020 months were noisier than the model allows (T2 VARIANCE_HIGH).
  - The refitted tipping year moved 7.2 years later, from 2073.8 to 2081.1 (T3 LATER).
  - The reproduction control failed: v1 gives 2086.5 on the authors' own file, not 2065. v1's likelihood
    has no penalty, while the authors' code penalises a < 1 (see Morr et al. 2026).
- **v2** ([DEVIATION-2026-11-24.md](DEVIATION-2026-11-24.md), written before v2 ran) repeats the forward test
  with a line-by-line replica of the authors' penalised estimator (`ddrep.py`). The replica must first
  reproduce both numbers they published from their own file (2056.52 and 2065). v2 also uses the SST variant
  that matches their file, put on their file's scale. The unpenalised estimate is reported alongside.

## Files

| file | what it does |
|---|---|
| `ddfit.py` | the DD23 model: OU baseline MLE, Strang-splitting likelihood (unpenalised), Nelder-Mead fit, simulation, bootstrap, forecast checks (standard library only) |
| `ddrep.py` | v2: replica of the authors' R objective (penalty, a floor, both Strang variants), its fit and bootstrap, affine harmonisation, v2 grades |
| `hadisst.py` | HadISST1 ASCII reader, Caesar et al. subpolar-gyre cells, monthly means, anomalies, DD23 file parser |
| `grade.py` | the frozen v1 grading rules |
| `job.py` / `job_v2.py` | the v1 and v2 jobs (download public data, fit, test, write results) |
| `test_*.py` | unit tests with hand-checked values (`python -m unittest`) |

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
