# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Munim, Inc. (MV2). Made with AI agents working in MV2's name.
import math
import random
import unittest

import ddfit
import grade


class TestOU(unittest.TestCase):
    def test_pure_decay(self):
        # x_k = 0.5 x_{k-1} exactly: rho 0.5, mu 0, no noise
        r = ddfit.ou_mle([1.0, 0.5, 0.25, 0.125, 0.0625], h=1.0)
        self.assertAlmostEqual(r["rho"], 0.5, places=12)
        self.assertAlmostEqual(r["mu0"], 0.0, places=12)
        self.assertAlmostEqual(r["alpha0"], math.log(2.0), places=12)
        self.assertAlmostEqual(r["sigma2"], 0.0, places=12)

    def test_with_intercept(self):
        # x_k = 1 + 0.5 x_{k-1}: mu = 1/(1-0.5) = 2
        r = ddfit.ou_mle([0.0, 1.0, 1.5, 1.75, 1.875], h=1.0)
        self.assertAlmostEqual(r["rho"], 0.5, places=12)
        self.assertAlmostEqual(r["mu0"], 2.0, places=12)

    def test_rejects_negative_rho(self):
        with self.assertRaises(ValueError):
            ddfit.ou_mle([0.0, 1.0, 0.0, 1.0, 0.0, 1.0])


class TestRamp(unittest.TestCase):
    def setUp(self):
        self.r = ddfit.Ramp(A=1.0, tau=100.0, mu0=0.5, alpha0=2.0, sigma2=0.3, t0=1924.0)

    def test_hand_values(self):
        r = self.r
        self.assertAlmostEqual(r.lam0, -1.0)
        self.assertAlmostEqual(r.m, -0.5)
        self.assertAlmostEqual(r.tc, 2024.0)
        self.assertAlmostEqual(r.mu(1900.0), 0.5)
        self.assertAlmostEqual(r.alpha(1900.0), 2.0)
        self.assertAlmostEqual(r.lam(1974.0), -0.5)
        self.assertAlmostEqual(r.alpha(1974.0), math.sqrt(2.0), places=12)
        self.assertAlmostEqual(r.mu(1974.0), -0.5 + math.sqrt(0.5), places=12)
        # DD23 identities: mu = m + sqrt(-lambda/A), alpha = 2 sqrt(-A lambda)
        for t in (1950.0, 2000.0, 2020.0):
            self.assertAlmostEqual(r.mu(t), r.m + math.sqrt(-r.lam(t) / r.A), places=12)
            self.assertAlmostEqual(r.alpha(t), 2.0 * math.sqrt(-r.A * r.lam(t)), places=12)

    def test_flow_inverse(self):
        y = ddfit._flow(0.3, 0.05)
        self.assertAlmostEqual(y, 0.3 / 1.015, places=12)
        self.assertAlmostEqual(ddfit._flow(y, -0.05), 0.3, places=12)
        self.assertIsNone(ddfit._flow(-30.0, 0.05))

    def test_strang_at_fixed_point(self):
        r = ddfit.Ramp(A=1.5, tau=120.0, mu0=0.3, alpha0=3.0, sigma2=0.3)
        z, mean, var, lj = ddfit.strang_terms(r, 1900.0, 0.3, 1900.0 + ddfit.H, 0.3)
        self.assertAlmostEqual(z, 0.3, places=12)
        self.assertAlmostEqual(mean, 0.3, places=12)
        self.assertAlmostEqual(lj, 0.0, places=12)
        v = 0.3 * (1.0 - math.exp(-0.5)) / 6.0      # sigma2 (1 - e^{-2 alpha h}) / (2 alpha)
        self.assertAlmostEqual(var, v, places=12)
        lp = ddfit.strang_logpdf(r, 1900.0, 0.3, 1900.0 + ddfit.H, 0.3)
        self.assertAlmostEqual(lp, -0.5 * math.log(2.0 * math.pi * v), places=10)
        self.assertAlmostEqual(lp, 1.0453, delta=2e-3)

    def test_strang_density_integrates_to_one(self):
        r = ddfit.Ramp(A=1.5, tau=120.0, mu0=0.3, alpha0=3.0, sigma2=0.3)
        tp, xp = 2000.0, 0.1
        t = tp + ddfit.H
        _, mean, var, _ = ddfit.strang_terms(r, tp, xp, t, 0.0)
        dx = 0.0005
        total = 0.0
        x = mean - 3.0
        while x < mean + 3.0:
            total += math.exp(ddfit.strang_logpdf(r, tp, xp, t, x)) * dx
            x += dx
        self.assertAlmostEqual(total, 1.0, delta=1e-3)


class TestOptimiser(unittest.TestCase):
    def test_quadratic(self):
        f = lambda u: (u[0] - 1.0) ** 2 + 10.0 * (u[1] + 2.0) ** 2
        u, fv, it = ddfit.nelder_mead(f, [0.0, 0.0])
        self.assertAlmostEqual(u[0], 1.0, delta=1e-3)
        self.assertAlmostEqual(u[1], -2.0, delta=1e-3)
        self.assertLess(fv, 1e-6)

    def test_quantile(self):
        self.assertEqual(ddfit.quantile([1.0, 2.0, 3.0, 4.0], 0.5), 2.5)
        self.assertEqual(ddfit.quantile([1.0, 2.0, 3.0, 4.0], 0.0), 1.0)
        self.assertEqual(ddfit.quantile([1.0, 2.0, 3.0, 4.0], 1.0), 4.0)
        self.assertIsNone(ddfit.quantile([], 0.5))


def _grid():
    return [1870.0 + (k + 1) / 12.0 for k in range(151 * 12)]   # Jan 1870 .. Dec 2020


class TestSimulationAndFit(unittest.TestCase):
    def test_deterministic_path(self):
        r = ddfit.Ramp(A=1.5, tau=120.0, mu0=0.3, alpha0=3.0, sigma2=0.0)
        ts = _grid()
        p = ddfit.simulate(r, ts, 0.3, random.Random(0), nsub=10)
        i_pre = max(i for i, t in enumerate(ts) if t < 1924.0)
        self.assertAlmostEqual(p[i_pre], 0.3, delta=1e-9)
        lag = p[-1] - r.mu(ts[-1])
        self.assertGreater(lag, 0.0)        # the path trails the falling stable point
        self.assertLess(lag, 0.02)

    def test_recovers_tipping_year_with_tiny_noise(self):
        r = ddfit.Ramp(A=1.5, tau=120.0, mu0=0.3, alpha0=3.0, sigma2=1e-6)
        ts = _grid()
        p = ddfit.simulate(r, ts, 0.3, random.Random(12345), nsub=10)
        f = ddfit.fit(ts, p)
        self.assertAlmostEqual(f["tc"], 2044.0, delta=8.0)
        self.assertFalse(f["tau_at_upper_bound"])

    def test_forecast_check_percentiles(self):
        r = ddfit.Ramp(A=1.5, tau=120.0, mu0=0.3, alpha0=3.0, sigma2=1e-8)
        ts = [2020.0 + k / 12.0 for k in range(25)]
        det = ddfit.simulate(ddfit.Ramp(1.5, 120.0, 0.3, 3.0, 0.0), ts, r.mu(2020.0), random.Random(0))
        fore_t, fore_x = ts[1:], det[1:]
        mid = ddfit.forecast_check(r, ts[0], det[0], fore_t, fore_x, 200, 7)
        q = mid["groups"]["all"]["q"]
        self.assertTrue(0.0 < q < 1.0)
        hi = ddfit.forecast_check(r, ts[0], det[0], fore_t, [x + 0.01 for x in fore_x], 200, 7)
        lo = ddfit.forecast_check(r, ts[0], det[0], fore_t, [x - 0.01 for x in fore_x], 200, 7)
        self.assertEqual(hi["groups"]["all"]["q"], 1.0)
        self.assertEqual(lo["groups"]["all"]["q"], 0.0)


class TestStats(unittest.TestCase):
    def test_phi_and_ks(self):
        self.assertAlmostEqual(ddfit.phi(0.0), 0.5)
        self.assertAlmostEqual(ddfit.phi(1.959964), 0.975, places=5)
        d, p = ddfit.ks_uniform([0.5])
        self.assertAlmostEqual(d, 0.5)
        self.assertTrue(0.0 <= p <= 1.0)

    def test_residual_summary(self):
        s = ddfit.residual_summary([1.0, -1.0, 1.0, -1.0, None])
        self.assertEqual(s["n"], 4)
        self.assertEqual(s["n_outside_domain"], 1)
        self.assertAlmostEqual(s["mean"], 0.0)
        self.assertAlmostEqual(s["var"], 4.0 / 3.0)


class TestGrades(unittest.TestCase):
    def test_recon(self):
        self.assertEqual(grade.grade_recon(0.995, 0.02), "RECON_OK")
        self.assertEqual(grade.grade_recon(0.995, 0.08), "RECON_WEAK")
        self.assertEqual(grade.grade_recon(0.90, 0.02), "RECON_FAIL")
        self.assertEqual(grade.grade_recon(None, None), "RECON_FAIL")

    def test_control(self):
        self.assertEqual(grade.grade_control(2068.0), "REPRO_PASS")
        self.assertEqual(grade.grade_control(2057.0), "REPRO_NEAR")
        self.assertEqual(grade.grade_control(2080.0), "REPRO_FAIL")

    def test_forecast_and_shift(self):
        self.assertEqual(grade.grade_forecast(0.01), "BELOW_FORECAST")
        self.assertEqual(grade.grade_forecast(0.5), "WITHIN_FORECAST")
        self.assertEqual(grade.grade_forecast(0.99), "ABOVE_FORECAST")
        self.assertEqual(grade.grade_shift(2065.0, 2058.0), "EARLIER")
        self.assertEqual(grade.grade_shift(2065.0, 2067.0), "STABLE")
        self.assertEqual(grade.grade_shift(2065.0, 2075.0), "LATER")
        self.assertEqual(grade.grade_shift(2065.0, 2075.0, unbounded=True), "UNBOUNDED")

    def test_onestep(self):
        ok = {"n": 60, "t_mean": 0.5, "z_var": -0.3, "ks_p": 0.4}
        self.assertEqual(grade.grade_onestep(ok), ["CALIBRATED"])
        bad = {"n": 60, "t_mean": -3.0, "z_var": 3.0, "ks_p": 0.001}
        self.assertEqual(grade.grade_onestep(bad), ["BIAS_LOW", "VARIANCE_HIGH", "KS_REJECT"])


if __name__ == "__main__":
    unittest.main()
