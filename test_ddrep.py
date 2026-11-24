# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Munim, Inc. (MV2). Made with AI agents working in MV2's name.
import math
import random
import unittest

import ddfit
import ddrep

T0 = 1924.0
H = 1.0 / 12.0


def _ramp_series(n=240):
    """Calendar times t0 + h k (k = 1..n) and values wiggling around the stable point."""
    r = ddfit.Ramp(A=1.2, tau=150.0, mu0=0.3, alpha0=3.0, sigma2=0.3, t0=T0)
    ts = [T0 + H * k for k in range(1, n + 1)]
    xs = [r.mu(t) + 0.05 * math.sin(0.7 * k) for k, t in enumerate(ts)]
    return r, ts, xs


class TestObjective(unittest.TestCase):
    def test_matches_ddfit_likelihood(self):
        # with the final half-flow at t_i (as in ddfit), the authors' objective is
        # -2 x ddfit's log-likelihood minus the dropped constants (n-1) log(2 pi)
        r, ts, xs = _ramp_series()
        ll = ddfit.ramp_loglik(r, ts, xs)
        obj = ddrep.minus2_loglik(r.tau, r.A, xs, r.mu0, r.alpha0, r.sigma2, pen=0.0,
                                  variant="corrected_mu_next")
        expect = -2.0 * ll - (len(xs) - 1) * math.log(2.0 * math.pi)
        self.assertAlmostEqual(obj, expect, delta=1e-6 * (1.0 + abs(expect)))

    def test_variants_differ_and_are_close(self):
        r, ts, xs = _ramp_series()
        args = (r.tau, r.A, xs, r.mu0, r.alpha0, r.sigma2)
        c = ddrep.minus2_loglik(*args, variant="corrected")
        o = ddrep.minus2_loglik(*args, variant="original")
        n = ddrep.minus2_loglik(*args, variant="corrected_mu_next")
        self.assertNotAlmostEqual(c, o, places=6)
        self.assertLess(abs(c - n), 0.05 * (1.0 + abs(c)))
        with self.assertRaises(ValueError):
            ddrep.minus2_loglik(*args, variant="nope")

    def test_penalty(self):
        r, ts, xs = _ramp_series()
        n = len(xs)
        base = ddrep.minus2_loglik(150.0, 0.5, xs, r.mu0, r.alpha0, r.sigma2, pen=0.0)
        pen = ddrep.minus2_loglik(150.0, 0.5, xs, r.mu0, r.alpha0, r.sigma2, pen=0.004)
        self.assertAlmostEqual(pen - base, 0.004 * n * (1.0 / 0.5 - 1.0), places=6)
        b2 = ddrep.minus2_loglik(150.0, 1.5, xs, r.mu0, r.alpha0, r.sigma2, pen=0.0)
        p2 = ddrep.minus2_loglik(150.0, 1.5, xs, r.mu0, r.alpha0, r.sigma2, pen=0.004)
        self.assertEqual(b2, p2)

    def test_floor_and_invalid(self):
        r, ts, xs = _ramp_series()
        f05 = ddrep.minus2_loglik(150.0, 0.05, xs, r.mu0, r.alpha0, r.sigma2, pen=0.004)
        f10 = ddrep.minus2_loglik(150.0, 0.1, xs, r.mu0, r.alpha0, r.sigma2, pen=0.004)
        self.assertEqual(f05, f10)
        # the record after t0 spans (240-1)/12 = 19.9 years: tau = 10 crosses lambda = 0
        self.assertEqual(ddrep.minus2_loglik(10.0, 1.0, xs, r.mu0, r.alpha0, r.sigma2), ddrep.BAD)
        self.assertEqual(ddrep.minus2_loglik(-5.0, 1.0, xs, r.mu0, r.alpha0, r.sigma2), ddrep.BAD)


class TestSplitAndAffine(unittest.TestCase):
    def test_split(self):
        ts = [1923.75, 1923.8333, 1923.9167, 1924.0, 1924.0833, 1924.1667]
        base, ramp = ddrep.split_dd23(ts, [1, 2, 3, 4, 5, 6])
        self.assertEqual(base, [1, 2, 3, 4])
        self.assertEqual(ramp, [5, 6])

    def test_affine(self):
        k, c, rms, n = ddrep.affine_fit([0.0, 1.0, 2.0, None], [1.0, 3.0, 5.0, 7.0])
        self.assertAlmostEqual(k, 2.0)
        self.assertAlmostEqual(c, 1.0)
        self.assertAlmostEqual(rms, 0.0)
        self.assertEqual(n, 3)


def _grid():
    return [1870.0 + (k + 1) / 12.0 for k in range(151 * 12)]   # Jan 1870 .. Dec 2020


class TestFit(unittest.TestCase):
    def test_recovers_and_is_local_optimum(self):
        r = ddfit.Ramp(A=1.5, tau=120.0, mu0=0.3, alpha0=3.0, sigma2=1e-4)
        ts = _grid()
        p = ddfit.simulate(r, ts, 0.3, random.Random(2026), nsub=10)
        f = ddrep.fit_dd23(ts, p, pen=0.0)
        self.assertAlmostEqual(f["tc"], 2044.0, delta=10.0)
        self.assertFalse(f["tau_at_upper_bound"])
        base, ramp = ddrep.split_dd23(ts, p)
        b = ddfit.ou_mle(base)

        def obj(tau, a):
            return ddrep.minus2_loglik(tau, a, ramp, b["mu0"], b["alpha0"], b["sigma2"])
        best = obj(f["tau"], f["a"])
        self.assertAlmostEqual(best, f["objective"], places=6)
        for dt, da in ((0.5, 0.0), (-0.5, 0.0), (0.0, 0.005), (0.0, -0.005)):
            self.assertGreaterEqual(obj(f["tau"] + dt, f["a"] + da), best - 1e-6)

    def test_ramp_of(self):
        f = {"a": 0.9, "tau": 140.0, "mu0": 0.25, "alpha0": 3.06, "sigma2": 0.3, "t0": T0}
        r = ddrep.ramp_of(f)
        self.assertAlmostEqual(r.tc, 2064.0)
        self.assertAlmostEqual(r.lam0, -3.06 ** 2 / 3.6)


class TestGrades(unittest.TestCase):
    def test_replica(self):
        self.assertEqual(ddrep.grade_replica(2056.6, 2064.0), "REPLICA_OK")
        self.assertEqual(ddrep.grade_replica(2058.5, 2069.0), "REPLICA_NEAR")
        self.assertEqual(ddrep.grade_replica(2070.0, 2086.0), "REPLICA_FAIL")
        self.assertEqual(ddrep.grade_replica(None, 2065.0), "REPLICA_FAIL")

    def test_harmonise(self):
        self.assertEqual(ddrep.grade_harmonise(0.004), "HARMONISE_OK")
        self.assertEqual(ddrep.grade_harmonise(0.02), "HARMONISE_NEAR")
        self.assertEqual(ddrep.grade_harmonise(0.2), "HARMONISE_FAIL")
        self.assertEqual(ddrep.grade_harmonise(None), "HARMONISE_FAIL")


if __name__ == "__main__":
    unittest.main()
