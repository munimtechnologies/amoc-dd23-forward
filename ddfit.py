# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Munim, Inc. (MV2). Made with AI agents working in MV2's name.
"""The Ditlevsen & Ditlevsen (2023) saddle-node ramp model, re-implemented.

Model (DD23, Nat. Commun. 14, 4254, eqs. 1-3, as corrected 2025):

    dX = -(A (X - m)^2 + lambda(t)) dt + sigma dB
    lambda(t) = lambda0 * (1 - Theta(t - t0) (t - t0) / tau),   t_c = t0 + tau

Stable point mu(t) = m + sqrt(-lambda(t)/A), linear rate alpha(t) = 2 sqrt(-A lambda(t)).
Baseline (t < t0): exact conditional Ornstein-Uhlenbeck MLE gives mu0, alpha0, sigma2;
then lambda0 = -alpha0^2/(4A) and m = mu0 - alpha0/(2A).  After t0, A and tau are
estimated by a Strang-splitting pseudo-likelihood:

    X_k = f_{h/2}( OU_h( f_{h/2}(X_{k-1}) ) ),   f_s(x) = mu + (x-mu)/(1 + A s (x-mu))

with the first half-flow and the OU step at t_{k-1} and the final half-flow inverted at
t_k (the 2025 correction says the final step must be evaluated at t_i).

Pure standard library on purpose: the tests run in a sandbox without numpy.
"""
import math
import random

H = 1.0 / 12.0          # one month, in years (DD23: Delta_t = 1 month)
T0 = 1924.0             # DD23 ramp onset, fixed
NEG_INF = float("-inf")
BIG = 1e300
A_MIN, A_MAX = 1e-4, 1e4
TAU_MAX = 1e5


# ----------------------------------------------------------------------------- baseline
def ou_mle(xs, h=H):
    """Exact conditional MLE of a stationary OU process sampled every h.

    AR(1) form x_k = c + rho x_{k-1} + e, rho = exp(-alpha h), mu = c/(1-rho),
    stationary variance gamma2 = sigma2/(2 alpha) = v/(1-rho^2).
    """
    n = len(xs) - 1
    if n < 3:
        raise ValueError("need at least 4 points")
    x0 = xs[:-1]
    x1 = xs[1:]
    sx = math.fsum(x0)
    sy = math.fsum(x1)
    sxx = math.fsum(a * a for a in x0)
    sxy = math.fsum(a * b for a, b in zip(x0, x1))
    den = n * sxx - sx * sx
    if den <= 0:
        raise ValueError("degenerate baseline")
    rho = (n * sxy - sx * sy) / den
    if not (0.0 < rho < 1.0):
        raise ValueError("baseline rho outside (0,1): %r" % rho)
    c = (sy - rho * sx) / n
    mu = c / (1.0 - rho)
    v = math.fsum((b - c - rho * a) ** 2 for a, b in zip(x0, x1)) / n
    alpha = -math.log(rho) / h
    sigma2 = 2.0 * alpha * v / (1.0 - rho * rho)
    return {"mu0": mu, "alpha0": alpha, "sigma2": sigma2, "rho": rho, "n": n}


# ----------------------------------------------------------------------------- model
class Ramp(object):
    def __init__(self, A, tau, mu0, alpha0, sigma2, t0=T0):
        self.A = float(A)
        self.tau = float(tau)
        self.mu0 = float(mu0)
        self.alpha0 = float(alpha0)
        self.sigma2 = float(sigma2)
        self.t0 = float(t0)
        self.lam0 = -self.alpha0 * self.alpha0 / (4.0 * self.A)
        self.K = self.alpha0 / (2.0 * self.A)      # sqrt(-lambda0/A): drop of mu at tipping
        self.m = self.mu0 - self.K

    def frac(self, t):
        if t <= self.t0:
            return 1.0
        return 1.0 - (t - self.t0) / self.tau

    def lam(self, t):
        return self.lam0 * self.frac(t)

    def alpha(self, t):
        f = self.frac(t)
        return self.alpha0 * math.sqrt(f) if f > 0 else 0.0

    def mu(self, t):
        f = self.frac(t)
        return self.m + self.K * math.sqrt(f) if f > 0 else self.m

    @property
    def tc(self):
        return self.t0 + self.tau

    def as_dict(self):
        return {"A": self.A, "tau": self.tau, "tc": self.tc, "t0": self.t0,
                "lambda0": self.lam0, "m": self.m, "mu0": self.mu0,
                "alpha0": self.alpha0, "sigma2": self.sigma2,
                "minus_A_lambda0": -self.A * self.lam0}


def _flow(y, a_s):
    """Exact flow of dy/dt = -A y^2 over time s, with a_s = A*s (s may be negative)."""
    d = 1.0 + a_s * y
    if d <= 0.0:
        return None
    return y / d


def strang_terms(r, tp, xp, t, x, h=H):
    """(z, mean, var, logjac) of the Strang transition (tp, xp) -> (t, x), or None."""
    a_p = r.alpha(tp)
    if a_p <= 0.0:
        return None
    mu_p = r.mu(tp)
    yf = _flow(xp - mu_p, r.A * h / 2.0)
    if yf is None:
        return None
    e = math.exp(-a_p * h)
    mean = mu_p + e * yf
    var = r.sigma2 * (1.0 - e * e) / (2.0 * a_p)
    if not var > 0.0:
        return None
    mu_c = r.mu(t)
    yc = x - mu_c
    dc = 1.0 - r.A * h / 2.0 * yc
    if dc <= 0.0:
        return None
    z = mu_c + yc / dc
    return z, mean, var, -2.0 * math.log(dc)


_LOG2PI = math.log(2.0 * math.pi)


def strang_logpdf(r, tp, xp, t, x, h=H):
    s = strang_terms(r, tp, xp, t, x, h)
    if s is None:
        return NEG_INF
    z, mean, var, lj = s
    return -0.5 * (_LOG2PI + math.log(var)) - (z - mean) ** 2 / (2.0 * var) + lj


def ramp_loglik(r, ts, xs, h=H):
    ll = 0.0
    for i in range(1, len(xs)):
        v = strang_logpdf(r, ts[i - 1], xs[i - 1], ts[i], xs[i], h)
        if v == NEG_INF:
            return NEG_INF
        ll += v
    return ll


# ----------------------------------------------------------------------------- optimiser
def nelder_mead(f, x0, step=0.3, max_iter=600, ftol=1e-10, xtol=1e-7):
    n = len(x0)
    pts = [list(x0)]
    for i in range(n):
        p = list(x0)
        p[i] += step
        pts.append(p)
    vals = [f(p) for p in pts]
    it = 0
    while it < max_iter:
        it += 1
        order = sorted(range(n + 1), key=lambda k: vals[k])
        pts = [pts[k] for k in order]
        vals = [vals[k] for k in order]
        spread = max(abs(pts[k][j] - pts[0][j]) for k in range(1, n + 1) for j in range(n))
        if abs(vals[-1] - vals[0]) <= ftol * (1.0 + abs(vals[0])) and spread <= xtol:
            break
        cen = [sum(pts[k][j] for k in range(n)) / n for j in range(n)]
        worst = pts[-1]
        xr = [cen[j] + (cen[j] - worst[j]) for j in range(n)]
        fr = f(xr)
        if fr < vals[0]:
            xe = [cen[j] + 2.0 * (cen[j] - worst[j]) for j in range(n)]
            fe = f(xe)
            if fe < fr:
                pts[-1], vals[-1] = xe, fe
            else:
                pts[-1], vals[-1] = xr, fr
        elif fr < vals[-2]:
            pts[-1], vals[-1] = xr, fr
        else:
            if fr < vals[-1]:
                xc = [cen[j] + 0.5 * (xr[j] - cen[j]) for j in range(n)]
            else:
                xc = [cen[j] + 0.5 * (worst[j] - cen[j]) for j in range(n)]
            fc = f(xc)
            if fc < min(fr, vals[-1]):
                pts[-1], vals[-1] = xc, fc
            else:
                best = pts[0]
                for k in range(1, n + 1):
                    pts[k] = [best[j] + 0.5 * (pts[k][j] - best[j]) for j in range(n)]
                    vals[k] = f(pts[k])
    k = min(range(n + 1), key=lambda q: vals[q])
    return pts[k], vals[k], it


def _linspace(a, b, n):
    return [a + (b - a) * i / (n - 1) for i in range(n)]


# ----------------------------------------------------------------------------- fit
def fit(ts, xs, t0=T0, h=H, start=None, grid_n=25, step=0.3, xtol=1e-7):
    """Fit the DD23 model to (ts, xs). Returns a dict (see Ramp.as_dict plus diagnostics)."""
    if len(ts) != len(xs):
        raise ValueError("length mismatch")
    i0 = None
    for i, t in enumerate(ts):
        if t >= t0:
            i0 = i
            break
    if i0 is None or i0 < 4 or i0 >= len(ts) - 3:
        raise ValueError("t0 not inside the record")
    base = ou_mle(xs[:i0], h)
    seg_t = ts[i0 - 1:]
    seg_x = xs[i0 - 1:]
    tspan = ts[-1] - t0

    def nll(u):
        A = math.exp(u[0])
        tau = math.exp(u[1])
        if not (A_MIN <= A <= A_MAX) or tau <= tspan * (1.0 + 1e-9) or tau > TAU_MAX:
            return BIG
        r = Ramp(A, tau, base["mu0"], base["alpha0"], base["sigma2"], t0)
        ll = ramp_loglik(r, seg_t, seg_x, h)
        return BIG if ll == NEG_INF else -ll

    if start is None:
        best = (BIG, None)
        for la in _linspace(math.log(1e-3), math.log(1e3), grid_n):
            for lt in _linspace(math.log(tspan + 0.5), math.log(5000.0), grid_n):
                v = nll([la, lt])
                if v < best[0]:
                    best = (v, [la, lt])
        if best[1] is None:
            raise ValueError("no finite likelihood on the grid")
        u0 = best[1]
    else:
        u0 = [math.log(start[0]), math.log(start[1])]
    u, fval, nit = nelder_mead(nll, u0, step=step, xtol=xtol)
    if fval >= BIG:
        raise ValueError("optimiser found no finite likelihood")
    r = Ramp(math.exp(u[0]), math.exp(u[1]), base["mu0"], base["alpha0"], base["sigma2"], t0)
    out = r.as_dict()
    out.update({
        "nll": fval, "iterations": nit, "n_base": base["n"] + 1, "n_ramp_transitions": len(seg_x) - 1,
        "t_first": ts[0], "t_last": ts[-1],
        "tau_at_upper_bound": r.tau > 0.99 * TAU_MAX,
        "A_at_bound": r.A < 1.01 * A_MIN or r.A > 0.99 * A_MAX,
    })
    return out


def ramp_from(fitres):
    return Ramp(fitres["A"], fitres["tau"], fitres["mu0"], fitres["alpha0"], fitres["sigma2"], fitres["t0"])


# ----------------------------------------------------------------------------- simulation
def simulate(r, ts, x_start, rng, nsub=10):
    """Euler-Maruyama of the full nonlinear SDE, recorded at ts (ts[0] is the start time).

    Returns the list of values at ts (first = x_start), or None if the path tipped
    (ran far below the unstable branch).
    """
    out = [x_start]
    x = x_start
    floor = r.m - 10.0 * abs(r.K) - 5.0
    for i in range(1, len(ts)):
        dt = (ts[i] - ts[i - 1]) / nsub
        sq = math.sqrt(max(r.sigma2, 0.0) * dt)
        t = ts[i - 1]
        for j in range(nsub):
            tj = t + j * dt
            x += -(r.A * (x - r.m) ** 2 + r.lam(tj)) * dt + sq * rng.gauss(0.0, 1.0)
            if not (x > floor):
                return None
        out.append(x)
    return out


def quantile(sorted_vals, q):
    """Type-7 (linear) quantile of an already sorted list."""
    n = len(sorted_vals)
    if n == 0:
        return None
    if n == 1:
        return sorted_vals[0]
    pos = q * (n - 1)
    lo = int(math.floor(pos))
    hi = min(lo + 1, n - 1)
    w = pos - lo
    return sorted_vals[lo] * (1.0 - w) + sorted_vals[hi] * w


def bootstrap(fitres, ts, B, seed, t0=T0, nsub=10, log=None):
    """DD23-style parametric bootstrap: simulate B records on the same time grid and refit."""
    r = ramp_from(fitres)
    rng = random.Random(seed)
    tcs = []
    tipped = 0
    failed = 0
    gamma = math.sqrt(r.sigma2 / (2.0 * r.alpha0))
    for b in range(B):
        x0 = r.mu0 + gamma * rng.gauss(0.0, 1.0)
        path = simulate(r, ts, x0, rng, nsub)
        if path is None:
            tipped += 1
            continue
        try:
            fb = fit(ts, path, t0, start=(fitres["A"], fitres["tau"]), step=0.1, xtol=1e-5)
        except ValueError:
            failed += 1
            continue
        tcs.append(fb["tc"])
        if log and (b + 1) % 100 == 0:
            log("bootstrap %d/%d" % (b + 1, B))
    tcs.sort()
    return {
        "B": B, "n_ok": len(tcs), "n_tipped_in_simulation": tipped, "n_fit_failed": failed,
        "tc_mean": (math.fsum(tcs) / len(tcs)) if tcs else None,
        "tc_q025": quantile(tcs, 0.025), "tc_q17": quantile(tcs, 1.0 / 6.0),
        "tc_q50": quantile(tcs, 0.5), "tc_q83": quantile(tcs, 5.0 / 6.0),
        "tc_q975": quantile(tcs, 0.975), "tcs": tcs,
    }


# ----------------------------------------------------------------------------- forecast checks
def phi(z):
    return 0.5 * math.erfc(-z / math.sqrt(2.0))


def forecast_check(r, t_start, x_start, ts_fore, xs_obs, N, seed, nsub=10, groups=None):
    """Simulate N paths from (t_start, x_start) over ts_fore; compare observed means.

    groups: dict name -> list of indices into ts_fore; each gets its own percentile.
    Returns percentiles q = P(sim mean < obs mean) (ties count half); tipped paths count
    as below every observation.
    """
    rng = random.Random(seed)
    full = [t_start] + list(ts_fore)
    if groups is None:
        groups = {}
    groups = dict(groups)
    groups["all"] = list(range(len(ts_fore)))
    sims = {g: [] for g in groups}
    tipped = 0
    for _ in range(N):
        p = simulate(r, full, x_start, rng, nsub)
        if p is None:
            tipped += 1
            for g in groups:
                sims[g].append(NEG_INF)
            continue
        p = p[1:]
        for g, idx in groups.items():
            sims[g].append(math.fsum(p[i] for i in idx) / len(idx))
    out = {"N": N, "n_tipped": tipped, "groups": {}}
    for g, idx in groups.items():
        obs = math.fsum(xs_obs[i] for i in idx) / len(idx)
        s = sorted(sims[g])
        below = sum(1 for v in s if v < obs)
        ties = sum(1 for v in s if v == obs)
        finite = [v for v in s if v != NEG_INF]
        out["groups"][g] = {
            "n_months": len(idx), "observed_mean": obs,
            "q": (below + 0.5 * ties) / N,
            "sim_q05": quantile(finite, 0.05), "sim_q50": quantile(finite, 0.5),
            "sim_q95": quantile(finite, 0.95),
        }
    return out


def one_step_residuals(r, ts, xs, h=H):
    """Standardised one-step Strang residuals for transitions ts[i-1] -> ts[i]."""
    res = []
    for i in range(1, len(xs)):
        s = strang_terms(r, ts[i - 1], xs[i - 1], ts[i], xs[i], h)
        if s is None:
            res.append(None)
            continue
        z, mean, var, _ = s
        res.append((z - mean) / math.sqrt(var))
    return res


def ks_uniform(us):
    """Kolmogorov-Smirnov D against U(0,1) and its asymptotic p-value."""
    u = sorted(us)
    n = len(u)
    if n == 0:
        return None, None
    d = 0.0
    for i, v in enumerate(u):
        d = max(d, (i + 1) / n - v, v - i / n)
    lam = (math.sqrt(n) + 0.12 + 0.11 / math.sqrt(n)) * d
    p = 0.0
    for k in range(1, 101):
        p += 2.0 * (-1) ** (k - 1) * math.exp(-2.0 * k * k * lam * lam)
    return d, max(0.0, min(1.0, p))


def residual_summary(res):
    vals = [v for v in res if v is not None]
    n = len(vals)
    if n < 3:
        return {"n": n}
    mean = math.fsum(vals) / n
    var = math.fsum((v - mean) ** 2 for v in vals) / (n - 1)
    d, p = ks_uniform([phi(v) for v in vals])
    return {
        "n": n, "n_outside_domain": len(res) - n, "mean": mean, "var": var,
        "t_mean": mean * math.sqrt(n), "z_var": (var - 1.0) / math.sqrt(2.0 / (n - 1)),
        "ks_d": d, "ks_p": p,
    }
