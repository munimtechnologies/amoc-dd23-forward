# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Munim, Inc. (MV2). Made with AI agents working in MV2's name.
"""A line-by-line Python replica of Ditlevsen & Ditlevsen's own R estimator (v2 of this study).

Source read on 2026-11-24: AMOCestimation.html in the authors' public archive
(doi:10.17894/ucph.b8f99b67-d4e6-4a2e-b518-00bddeed323b, original 2023 code; the 2025 replacement
archive doi:10.17894/ucph.9ef73d8c-b642-4b55-88bc-b38984d043b9 changes only the final Strang step).
Their R function loglik.fun(pars = (tau, a), data, delta, alpha0, mu0, sigma20, pen):

    a = max(0.1, a);  m = mu0 - alpha0/(2a);  lambda0 = -alpha0^2/(4a)
    Xupp = data[2:n];  Xlow = data[1:(n-1)];  time = delta*(1:(n-1))
    lam = lambda0*(1 - time/tau);  alpha = 2 sqrt(-a lam);  mu = m + sqrt(-lam/a)
    gamma2 = sigma2/(2 alpha);  rho = exp(-alpha delta)
    tmp = a delta (Xlow - mu)/2;  fh.half = (mu tmp + Xlow)/(tmp + 1)
    fh.half.inv = (mu tmp' - Xupp)/(tmp' - 1)      (tmp' = tmp in the 2023 code: the bug;
                                                    tmp' = a delta (Xupp - mu)/2 once corrected)
    mu.h = fh.half rho + mu (1 - rho);  m.part = fh.half.inv - mu.h;  var.part = gamma2 (1 - rho^2)
    det.Dfh.half.inv = 1/(a delta (Xupp - mu)/2 - 1)^2
    loglik = -sum(log var.part) - sum(m.part^2/var.part) + 2 sum(log det.Dfh.half.inv)
             - pen n (1/a - 1) (a < 1)
    return -loglik                                   (that is, -2 x log-likelihood + penalty)

with data = the series after t0 = 1924 (time > t0), baseline OU estimates from time <= t0,
pen = 0.004 and optim(c(tau = 100, a = 1), Nelder-Mead).  The penalty is the one Morr et al.
(arXiv:2604.20341) analysed; v1 of this study used the unpenalised likelihood.

Pure standard library.
"""
import math
import random

import ddfit

H = 1.0 / 12.0
T0 = 1924.0
PEN_DD23 = 0.004
A_FLOOR = 0.1
TAU_MAX = 1e5
BAD = 1e300                    # R returns 10000 for NA; the optimum lies far inside the valid region
DD23_ORIGINAL_TC = 2056.52     # printed by the authors' 2023 code (tau 132.52 + 1924)
DD23_CORRECTED_TC = 2065.0     # corrected article (2025): "the optimal fit is tc = 2065"
VARIANTS = ("corrected", "original", "corrected_mu_next")


def split_dd23(ts, xs, t0=T0):
    """Baseline = time <= t0, ramp = time > t0 (as in the authors' R code)."""
    base = [x for t, x in zip(ts, xs) if t <= t0 + 1e-6]
    ramp = [x for t, x in zip(ts, xs) if t > t0 + 1e-6]
    return base, ramp


def minus2_loglik(tau, a, ramp, mu0, alpha0, sigma2, h=H, pen=0.0, variant="corrected"):
    """The authors' objective (-2 x Strang log-likelihood without constants, plus penalty)."""
    if variant not in VARIANTS:
        raise ValueError("unknown variant %r" % variant)
    if not (tau > 0.0) or tau > TAU_MAX:
        return BAD
    a = max(A_FLOOR, a)
    n = len(ramp)
    m = mu0 - alpha0 / (2.0 * a)
    lam0 = -alpha0 * alpha0 / (4.0 * a)
    s_var = 0.0
    s_m2 = 0.0
    s_det = 0.0
    for i in range(1, n):
        time = h * i
        lam = lam0 * (1.0 - time / tau)
        if not lam < 0.0:
            return BAD
        alpha = 2.0 * math.sqrt(-a * lam)
        mu = m + math.sqrt(-lam / a)
        if variant == "corrected_mu_next":
            lam_n = lam0 * (1.0 - (time + h) / tau)
            if not lam_n < 0.0:
                return BAD
            mu_inv = m + math.sqrt(-lam_n / a)
        else:
            mu_inv = mu
        xl = ramp[i - 1]
        xu = ramp[i]
        gamma2 = sigma2 / (2.0 * alpha)
        rho = math.exp(-alpha * h)
        tmp = a * h * (xl - mu) / 2.0
        if tmp + 1.0 == 0.0:
            return BAD
        fh = (mu * tmp + xl) / (tmp + 1.0)
        tmp_inv = tmp if variant == "original" else a * h * (xu - mu_inv) / 2.0
        if tmp_inv - 1.0 == 0.0:
            return BAD
        fh_inv = (mu_inv * tmp_inv - xu) / (tmp_inv - 1.0)
        mu_h = fh * rho + mu * (1.0 - rho)
        var = gamma2 * (1.0 - rho * rho)
        dd = a * h * (xu - mu_inv) / 2.0 - 1.0
        if not var > 0.0 or dd == 0.0:
            return BAD
        s_var += math.log(var)
        s_m2 += (fh_inv - mu_h) ** 2 / var
        s_det += math.log(1.0 / (dd * dd))
    ll = -s_var - s_m2 + 2.0 * s_det
    if a < 1.0:
        ll -= pen * n * (1.0 / a - 1.0)
    return -ll


def _geom(a, b, n):
    la, lb = math.log(a), math.log(b)
    return [math.exp(la + (lb - la) * i / (n - 1)) for i in range(n)]


def fit_dd23(ts, xs, pen=PEN_DD23, variant="corrected", t0=T0, h=H, start=(100.0, 1.0),
             grid=True, grid_n=15, restarts=2, xtol=1e-9, ftol=1e-13):
    """Fit tau and a with the authors' objective. Returns a dict with tc = t0 + tau.

    Optimisation: Nelder-Mead in (tau/100, a) from the authors' start (100, 1) (moved just past
    the record length if the record is longer than 100 years after t0), then, if grid=True,
    again from the best point of a 15 x 15 grid; the better optimum is kept and both are reported.
    """
    base, ramp = split_dd23(ts, xs, t0)
    if len(base) < 4 or len(ramp) < 4:
        raise ValueError("t0 not inside the record")
    b = ddfit.ou_mle(base, h)
    span = h * (len(ramp) - 1)

    def obj(tau, a):
        return minus2_loglik(tau, a, ramp, b["mu0"], b["alpha0"], b["sigma2"], h, pen, variant)

    def f(u):
        return obj(100.0 * u[0], u[1])

    def run(u0, step):
        u, v, it = ddfit.nelder_mead(f, u0, step=step, max_iter=3000, ftol=ftol, xtol=xtol)
        for _ in range(restarts):
            u, v, k = ddfit.nelder_mead(f, u, step=0.02, max_iter=3000, ftol=ftol, xtol=xtol)
            it += k
        return u, v, it

    tau_s = max(float(start[0]), span + 5.0)
    u_r, v_r, it_r = run([tau_s / 100.0, float(start[1])], 0.1)
    best_u, best_v, how = u_r, v_r, "authors_start"
    v_g = None
    if grid:
        gbest = (BAD, None)
        for a in _geom(A_FLOOR, 10.0, grid_n):
            for tau in _geom(span + 0.5, 3000.0, grid_n):
                v = obj(tau, a)
                if v < gbest[0]:
                    gbest = (v, [tau / 100.0, a])
        if gbest[1] is not None:
            u_g, v_g, _ = run(gbest[1], 0.1)
            if v_g < best_v - 1e-9:
                best_u, best_v, how = u_g, v_g, "grid_start"
    if best_v >= BAD:
        raise ValueError("no finite objective")
    tau = 100.0 * best_u[0]
    a = max(A_FLOOR, best_u[1])
    m = b["mu0"] - b["alpha0"] / (2.0 * a)
    return {
        "tau": tau, "a": a, "A": a, "tc": t0 + tau, "t0": t0, "m": m,
        "lambda0": -b["alpha0"] ** 2 / (4.0 * a), "mu0": b["mu0"], "alpha0": b["alpha0"],
        "sigma2": b["sigma2"], "objective": best_v, "objective_authors_start": v_r,
        "objective_grid_start": v_g, "chosen": how, "pen": pen, "variant": variant,
        "n_base": len(base), "n_ramp": len(ramp), "penalty_active": a < 1.0 and pen > 0.0,
        "tau_at_upper_bound": tau > 0.99 * TAU_MAX, "iterations_authors_start": it_r,
    }


def ramp_of(fitres):
    return ddfit.Ramp(fitres["a"], fitres["tau"], fitres["mu0"], fitres["alpha0"], fitres["sigma2"],
                      fitres["t0"])


def bootstrap_dd23(fitres, ts, B, seed, nsub=10, log=None):
    """Parametric bootstrap as in the authors' code: simulate the fitted model, refit (same pen)."""
    r = ramp_of(fitres)
    rng = random.Random(seed)
    tcs = []
    tipped = 0
    failed = 0
    gamma = math.sqrt(r.sigma2 / (2.0 * r.alpha0))
    for k in range(B):
        x0 = r.mu0 + gamma * rng.gauss(0.0, 1.0)
        path = ddfit.simulate(r, ts, x0, rng, nsub)
        if path is None:
            tipped += 1
            continue
        try:
            fb = fit_dd23(ts, path, pen=fitres["pen"], variant=fitres["variant"], t0=fitres["t0"],
                          start=(fitres["tau"], fitres["a"]), grid=False, restarts=1,
                          xtol=1e-6, ftol=1e-10)
        except ValueError:
            failed += 1
            continue
        tcs.append(fb["tc"])
        if log and (k + 1) % 100 == 0:
            log("penalised bootstrap %d/%d" % (k + 1, B))
    tcs.sort()
    q = ddfit.quantile
    return {
        "B": B, "n_ok": len(tcs), "n_tipped_in_simulation": tipped, "n_fit_failed": failed,
        "tc_mean": (math.fsum(tcs) / len(tcs)) if tcs else None,
        "tc_q025": q(tcs, 0.025), "tc_q17": q(tcs, 1.0 / 6.0), "tc_q50": q(tcs, 0.5),
        "tc_q83": q(tcs, 5.0 / 6.0), "tc_q975": q(tcs, 0.975), "tcs": tcs,
    }


def affine_fit(x, y):
    """Least squares y = c + k x over pairs where both exist. Returns (k, c, rms_resid, n)."""
    pairs = [(a, b) for a, b in zip(x, y) if a is not None and b is not None]
    n = len(pairs)
    if n < 3:
        raise ValueError("too few pairs")
    mx = math.fsum(a for a, _ in pairs) / n
    my = math.fsum(b for _, b in pairs) / n
    sxx = math.fsum((a - mx) ** 2 for a, _ in pairs)
    sxy = math.fsum((a - mx) * (b - my) for a, b in pairs)
    if sxx <= 0.0:
        raise ValueError("constant x")
    k = sxy / sxx
    c = my - k * mx
    rms = math.sqrt(math.fsum((b - c - k * a) ** 2 for a, b in pairs) / n)
    return k, c, rms, n


# ----------------------------------------------------------------------------- v2 grades (frozen)
def grade_replica(tc_original, tc_corrected):
    """Does the replica reproduce both numbers the authors published from their own file?"""
    if tc_original is None or tc_corrected is None:
        return "REPLICA_FAIL"
    d1 = abs(tc_original - DD23_ORIGINAL_TC)
    d2 = abs(tc_corrected - DD23_CORRECTED_TC)
    if d1 <= 1.0 and d2 <= 3.0:
        return "REPLICA_OK"
    if d1 <= 3.0 and d2 <= 6.0:
        return "REPLICA_NEAR"
    return "REPLICA_FAIL"


def grade_harmonise(rms):
    """Our harmonised HadISST rebuild vs the authors' AMOC2 column, 1870-2020 (K)."""
    if rms is None:
        return "HARMONISE_FAIL"
    if rms <= 0.01:
        return "HARMONISE_OK"
    if rms <= 0.03:
        return "HARMONISE_NEAR"
    return "HARMONISE_FAIL"
