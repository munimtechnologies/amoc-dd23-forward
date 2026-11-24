# amoc-dd23-forward

An open, pre-registered **out-of-sample test of the Ditlevsen & Ditlevsen (2023) forecast of an Atlantic
overturning (AMOC) collapse**. That forecast is a tipping year of 2065 (95%: 2037-2109) after the authors' 2025
correction, first reported as about 2057. It was fitted to sea-surface temperatures from 1870 to 2020. This project
asks two questions nobody has published answers to:

1. Did the fitted model anticipate what the fingerprint (subpolar-gyre SST minus twice global-mean SST, HadISST1)
   did from January 2021 to the latest month? (T1, T2)
2. Where does the model's own tipping-year estimate move when those months are added? (T3)

The design, the thresholds and every grade label were frozen in [PREREGISTRATION.md](PREREGISTRATION.md) on
2026-11-24, before the job ran.

## Files

| file | what it does |
|---|---|
| `ddfit.py` | the DD23 model: OU baseline MLE, Strang-splitting likelihood, Nelder-Mead fit, simulation, bootstrap, forecast checks (standard library only) |
| `hadisst.py` | HadISST1 ASCII reader, Caesar et al. subpolar-gyre cells, monthly means, anomalies, DD23 file parser |
| `grade.py` | the frozen grading rules |
| `job.py` | downloads the public data, runs the gate, the reproduction control, the fits and the tests; writes `results.json` and CSVs |
| `test_ddfit.py`, `test_hadisst.py` | unit tests with hand-checked values (`python -m unittest`) |

Run: `python -m unittest -v && python job.py --data-dir data --out-dir out` (Python 3.10+, no packages needed;
about 1 GB of downloads).

## Credit and data

- P. D. Ditlevsen & S. Ditlevsen, *Nature Communications* 14, 4254 (2023), and Author Correction (2025); code and
  data at doi:10.17894/ucph.9ef73d8c-b642-4b55-88bc-b38984d043b9 (CC BY 4.0). The estimator fix is credited there
  to Anders Gantzhorn Kristensen.
- L. Caesar et al., *Nature* (2018): the SST fingerprint and the subpolar-gyre grid (downloaded at run time from
  PIK, not redistributed here).
- HadISST1, Met Office Hadley Centre (Rayner et al. 2003).
- Critiques this test complements: Ben-Yami et al. 2024 (*Science Advances*), Morr et al. 2026 (arXiv:2604.20341).

This is research software, not a climate projection service. Made with AI agents working in the name of MV2
(Munim, Inc., mv2@munimtech.com); no human expert has checked it yet. Licence: Apache-2.0.
