# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Munim, Inc. (MV2). Made with AI agents working in MV2's name.
"""amoc-dd23-forward v2 (DEVIATION-2026-11-24.md): the forward test with the authors' own estimator.

v1 (pre-registered) used an unpenalised likelihood and picked the plain-mean HadISST variant; its
control gave 2086 on the authors' file instead of their 2065. v2, declared after v1 ran:
  R. replica check: a line-by-line copy of the authors' R estimator (penalty pen = 0.004 on a < 1)
     must give their two published numbers from their own file: 2056.52 (2023 code) and 2065 (2025).
  H. fingerprint: area-weighted, ice at -1.8 C (the variant with r = 0.99999 to their file in v1),
     mapped onto their file's scale by an affine fit over 1870-2020 (SG and GM separately).
  F. forward test as pre-registered (T1, T2, T3), for P = their penalised estimator and
     U = the unpenalised MLE (v1's estimator), on the harmonised series through 2020 and through
     the last HadISST month.
Pure standard library.
"""
import argparse
import json
import math
import os
import sys
import time
import urllib.request

import ddfit
import ddrep
import grade
import hadisst as hs

HADISST_BASE = "https://www.metoffice.gov.uk/hadobs/hadisst/data/"
SG_GRID_URL = "https://www.pik-potsdam.de/~caesar/AMOC_slowdown/sg_grid.txt"
DD23_URL = "https://erda.ku.dk/archives/afce4e1c3ac6f0db61c27cb45c2e9b14/AMOCdata.txt"
UA = "mv2-research-job/1.0 (open science; HadISST reader)"
VARIANT = "w_ice"
T_START = time.time()


def log(msg):
    print("[%8.1fs] %s" % (time.time() - T_START, msg), flush=True)


def download(url, path, refresh=False, tries=5):
    if os.path.exists(path) and os.path.getsize(path) > 0 and not refresh:
        return path
    err = None
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=300) as resp, open(path + ".part", "wb") as f:
                while True:
                    chunk = resp.read(1 << 20)
                    if not chunk:
                        break
                    f.write(chunk)
            os.replace(path + ".part", path)
            return path
        except Exception as e:
            err = e
            log("download failed (%d/%d) %s: %r" % (k + 1, tries, url, e))
            time.sleep(10 * (k + 1))
    raise RuntimeError("could not download %s: %r" % (url, err))


def write_csv(path, header, rows):
    with open(path, "w", encoding="utf-8") as f:
        f.write(",".join(header) + "\n")
        for r in rows:
            f.write(",".join("" if v is None else (("%.6f" % v) if isinstance(v, float) else str(v))
                             for v in r) + "\n")


def strip_tcs(b):
    d = dict(b)
    d.pop("tcs", None)
    return d


def load_hadisst(data_dir, last_year, sg_cells):
    files = ["HadISST1_SST_1870-1900.txt.gz", "HadISST1_SST_1901-1930.txt.gz",
             "HadISST1_SST_1931-1960.txt.gz", "HadISST1_SST_1961-1990.txt.gz",
             "HadISST1_SST_1991-2003.txt.gz"]
    files += ["HadISST1_SST_%d.txt.gz" % y for y in range(2004, last_year + 1)]
    weights = hs.row_weights()
    stats = {}
    for name in files:
        p = download(HADISST_BASE + name, os.path.join(data_dir, name),
                     refresh=name.endswith("%d.txt.gz" % last_year))
        with hs.open_text(p) as f:
            for yr, mon, vals in hs.iter_months(f):
                if (yr, mon) not in stats:
                    stats[(yr, mon)] = hs.month_stats(vals, sg_cells, weights)
        log("%s read; %d months so far" % (name, len(stats)))
    try:
        p = download(HADISST_BASE + "HadISST1_SST_update.txt.gz",
                     os.path.join(data_dir, "HadISST1_SST_update.txt.gz"), refresh=True, tries=2)
        with hs.open_text(p) as f:
            for yr, mon, vals in hs.iter_months(f):
                if (yr, mon) not in stats:
                    stats[(yr, mon)] = hs.month_stats(vals, sg_cells, weights)
                    log("update file added %d-%02d" % (yr, mon))
    except Exception as e:
        log("update file skipped: %r" % (e,))
    keys = sorted(stats)
    if keys[0] != (1870, 1):
        raise RuntimeError("record does not start in 1870-01")
    for k1, k2 in zip(keys, keys[1:]):
        if hs.next_month(*k1) != k2:
            raise RuntimeError("gap in HadISST record after %r" % (k1,))
    return keys, stats


def forward(name, fit20, fitlt, ramp20, ts20, x20, fore_t, fore_x, groups, nfore, seed, nsub):
    t1 = ddfit.forecast_check(ramp20, ts20[-1], x20[-1], fore_t, fore_x, nfore, seed, nsub=nsub,
                              groups=groups)
    t1_grade = grade.grade_forecast(t1["groups"]["all"]["q"])
    res = ddfit.one_step_residuals(ramp20, [ts20[-1]] + fore_t, [x20[-1]] + fore_x)
    t2 = ddfit.residual_summary(res)
    unb = fit20["tau_at_upper_bound"] or fitlt["tau_at_upper_bound"]
    t3 = grade.grade_shift(fit20["tc"], fitlt["tc"], unbounded=unb)
    log("%s: tc %.2f -> %.2f (%s); T1 q=%.3f %s; T2 %s" % (
        name, fit20["tc"], fitlt["tc"], t3, t1["groups"]["all"]["q"], t1_grade, grade.grade_onestep(t2)))
    return {"T1_forecast": dict(t1, grade=t1_grade),
            "T2_one_step": {"flags": grade.grade_onestep(t2), "out_of_sample": t2},
            "T3_shift": {"grade": t3, "tc_2020": fit20["tc"], "tc_latest": fitlt["tc"],
                         "shift_years": fitlt["tc"] - fit20["tc"]}}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="/data")
    ap.add_argument("--out-dir", default="/out")
    ap.add_argument("--boot-p", type=int, default=1000)
    ap.add_argument("--boot-u", type=int, default=500)
    ap.add_argument("--nfore", type=int, default=2000)
    ap.add_argument("--nsub", type=int, default=10)
    ap.add_argument("--seed", type=int, default=20261125)
    ap.add_argument("--last-year", type=int, default=2026)
    a = ap.parse_args(argv)
    os.makedirs(a.data_dir, exist_ok=True)
    os.makedirs(a.out_dir, exist_ok=True)
    inputs = os.path.join(a.out_dir, "inputs")
    os.makedirs(inputs, exist_ok=True)

    sg_path = download(SG_GRID_URL, os.path.join(inputs, "sg_grid.txt"))
    dd_path = download(DD23_URL, os.path.join(inputs, "dd23_AMOCdata.txt"))
    with open(sg_path, "r", encoding="latin-1") as f:
        sg_cells = hs.sg_cells_from_grid(f.read())
    with open(dd_path, "r", encoding="latin-1") as f:
        dd = hs.parse_dd23(f.read())
    if len(sg_cells) < 500 or len(dd) < 1500:
        raise RuntimeError("unexpected SG grid or DD23 file")
    dd_t = [r["t"] for r in dd]
    out = {"study": "amoc-dd23-forward", "version": "v2",
           "deviation": "DEVIATION-2026-11-24.md (written after v1 ran, before v2 ran)"}

    # ------------------------------------------------------------------ R. replica check
    rep = {}
    for key, col, variant, pen in (("original_pen", "AMOC2", "original", ddrep.PEN_DD23),
                                   ("corrected_pen", "AMOC2", "corrected", ddrep.PEN_DD23),
                                   ("corrected_mu_next_pen", "AMOC2", "corrected_mu_next", ddrep.PEN_DD23),
                                   ("corrected_nopen", "AMOC2", "corrected", 0.0),
                                   ("original_nopen", "AMOC2", "original", 0.0),
                                   ("corrected_pen_beta1", "AMOC1", "corrected", ddrep.PEN_DD23)):
        xs = [r[col] for r in dd]
        rep[key] = ddrep.fit_dd23(dd_t, xs, pen=pen, variant=variant)
        log("replica %s: tc=%.2f a=%.4f (%s)" % (key, rep[key]["tc"], rep[key]["a"], rep[key]["chosen"]))
    x3 = [r["AMOC0"] - 3.0 * r["GM"] for r in dd]
    rep["corrected_pen_beta3"] = ddrep.fit_dd23(dd_t, x3, pen=ddrep.PEN_DD23, variant="corrected")
    u_dd = ddfit.fit(dd_t, [r["AMOC2"] for r in dd])
    rep["v1_estimator"] = u_dd
    replica_grade = ddrep.grade_replica(rep["original_pen"]["tc"], rep["corrected_pen"]["tc"])
    out["replica"] = {"grade": replica_grade, "fits": rep,
                      "published": {"original_2023_code": ddrep.DD23_ORIGINAL_TC,
                                    "corrected_2025": ddrep.DD23_CORRECTED_TC,
                                    "table1_beta1": 2123.0, "table1_beta3": 2068.0}}
    log("REPLICA %s" % replica_grade)

    # ------------------------------------------------------------------ H. fingerprint
    keys, stats = load_hadisst(a.data_dir, a.last_year, sg_cells)
    last = keys[-1]
    sg_raw = [hs.variant_means(stats[k]["sg"])[VARIANT] for k in keys]
    gm_raw = [hs.variant_means(stats[k]["gm"])[VARIANT] for k in keys]
    sg_an = hs.monthly_anomalies(keys, sg_raw)
    gm_an = hs.monthly_anomalies(keys, gm_raw)
    pos = {k: i for i, k in enumerate(keys)}
    dd_idx = [pos.get((r["yr"], r["mon"])) for r in dd]

    def aligned(series):
        return [series[i] if i is not None else None for i in dd_idx]

    k1, c1, rms1, n1 = ddrep.affine_fit(aligned(sg_an), [r["AMOC0"] for r in dd])
    k2, c2, rms2, n2 = ddrep.affine_fit(aligned(gm_an), [r["GM"] for r in dd])
    sg_h = [c1 + k1 * v for v in sg_an]
    gm_h = [c2 + k2 * v for v in gm_an]

    def fp(beta):
        return [s - beta * g for s, g in zip(sg_h, gm_h)]

    xh = fp(2.0)
    r_h, rms_h, n_h = hs.corr_rms(aligned(xh), [r["AMOC2"] for r in dd])
    harm_grade = ddrep.grade_harmonise(rms_h)
    xraw = [s - 2.0 * g for s, g in zip(sg_an, gm_an)]
    r_raw, rms_raw, _ = hs.corr_rms(aligned(xraw), [r["AMOC2"] for r in dd])
    out["harmonise"] = {"grade": harm_grade, "variant": VARIANT,
                        "sg": {"k": k1, "c": c1, "rms_resid": rms1, "n": n1},
                        "gm": {"k": k2, "c": c2, "rms_resid": rms2, "n": n2},
                        "amoc2_r": r_h, "amoc2_rms": rms_h, "n": n_h,
                        "unharmonised_amoc2_r": r_raw, "unharmonised_amoc2_rms": rms_raw}
    log("HARMONISE %s: SG k=%.4f c=%.4f; GM k=%.4f c=%.4f; AMOC2 r=%.6f rms=%.5f" % (
        harm_grade, k1, c1, k2, c2, r_h, rms_h))

    ts = [hs.t_of(y, m) for (y, m) in keys]
    i20 = max(i for i, k in enumerate(keys) if k[0] <= 2020)
    ts20, x20 = ts[:i20 + 1], xh[:i20 + 1]
    fore_t, fore_x = ts[i20 + 1:], xh[i20 + 1:]
    nf = len(fore_t)
    if nf < 12:
        raise RuntimeError("fewer than 12 out-of-sample months")
    out.update({"hadisst_last_month": "%d-%02d" % last, "n_months": len(keys), "n_out_of_sample_months": nf})
    years = sorted(set(k[0] for k in keys[i20 + 1:]))
    groups = {"year_%d" % y: [j for j, k in enumerate(keys[i20 + 1:]) if k[0] == y] for y in years}
    groups["last_12_months"] = list(range(max(0, nf - 12), nf))

    # ------------------------------------------------------------------ F. forward test, P and U
    p20 = ddrep.fit_dd23(ts20, x20, pen=ddrep.PEN_DD23, variant="corrected")
    plt = ddrep.fit_dd23(ts, xh, pen=ddrep.PEN_DD23, variant="corrected")
    u20 = ddfit.fit(ts20, x20)
    ult = ddfit.fit(ts, xh)
    fP = forward("P", p20, plt, ddrep.ramp_of(p20), ts20, x20, fore_t, fore_x, groups,
                 a.nfore, a.seed + 1, a.nsub)
    fU = forward("U", u20, ult, ddfit.ramp_from(u20), ts20, x20, fore_t, fore_x, groups,
                 a.nfore, a.seed + 2, a.nsub)
    out["P_penalised"] = dict(fP, fit_2020=p20, fit_latest=plt)
    out["U_unpenalised"] = dict(fU, fit_2020=u20, fit_latest=ult)

    # ------------------------------------------------------------------ secondary (report only)
    sec = {}
    for beta in (1.0, 3.0):
        xb = fp(beta)
        try:
            b20 = ddrep.fit_dd23(ts20, xb[:i20 + 1], pen=ddrep.PEN_DD23)
            blt = ddrep.fit_dd23(ts, xb, pen=ddrep.PEN_DD23)
            sec["P_beta_%d" % int(beta)] = {"tc_2020": b20["tc"], "tc_latest": blt["tc"],
                                            "shift": blt["tc"] - b20["tc"]}
        except ValueError as e:
            sec["P_beta_%d" % int(beta)] = {"error": str(e)}
    try:
        r20 = ddrep.fit_dd23(ts20, xraw[:i20 + 1], pen=ddrep.PEN_DD23)
        rlt = ddrep.fit_dd23(ts, xraw, pen=ddrep.PEN_DD23)
        sec["P_unharmonised_scale"] = {"tc_2020": r20["tc"], "tc_latest": rlt["tc"],
                                       "shift": rlt["tc"] - r20["tc"]}
    except ValueError as e:
        sec["P_unharmonised_scale"] = {"error": str(e)}
    # month-by-month path of the penalised tipping year as data are added (yearly steps)
    path = []
    for y in range(2020, last[0] + 1):
        j = max(i for i, k in enumerate(keys) if k[0] <= y)
        fy = ddrep.fit_dd23(ts[:j + 1], xh[:j + 1], pen=ddrep.PEN_DD23, start=(plt["tau"], plt["a"]))
        uy = ddfit.fit(ts[:j + 1], xh[:j + 1])
        path.append(["%d-%02d" % keys[j], fy["tc"], uy["tc"]])
        log("data to %d-%02d: P tc=%.2f, U tc=%.2f" % (keys[j][0], keys[j][1], fy["tc"], uy["tc"]))
    sec["tc_by_data_end"] = path
    out["secondary"] = sec

    with open(os.path.join(a.out_dir, "results_v2_point.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, indent=1, sort_keys=True)

    rows = []
    ddmap = {(r["yr"], r["mon"]): r["AMOC2"] for r in dd}
    for i, k in enumerate(keys):
        rows.append(["%d-%02d" % k, ts[i], xh[i], ddmap.get(k)])
    write_csv(os.path.join(a.out_dir, "fingerprint_v2_monthly.csv"),
              ["month", "t", "amoc2_harmonised", "amoc2_dd23"], rows)
    ann = []
    rP20, rPlt = ddrep.ramp_of(p20), ddrep.ramp_of(plt)
    for y in sorted(set(k[0] for k in keys)):
        idx = [i for i, k in enumerate(keys) if k[0] == y]
        ann.append([y, len(idx), math.fsum(xh[i] for i in idx) / len(idx),
                    math.fsum(rP20.mu(ts[i]) for i in idx) / len(idx),
                    math.fsum(rPlt.mu(ts[i]) for i in idx) / len(idx)])
    write_csv(os.path.join(a.out_dir, "fingerprint_v2_annual.csv"),
              ["year", "n_months", "amoc2_harmonised", "mu_P2020", "mu_Platest"], ann)
    for tag, fw in (("P", fP), ("U", fU)):
        fan = []
        for y in years:
            g = fw["T1_forecast"]["groups"]["year_%d" % y]
            fan.append([y, g["n_months"], g["observed_mean"], g["sim_q05"], g["sim_q50"], g["sim_q95"], g["q"]])
        write_csv(os.path.join(a.out_dir, "forecast_fan_v2_%s.csv" % tag),
                  ["year", "n_months", "observed", "model_q05", "model_q50", "model_q95", "q"], fan)
    write_csv(os.path.join(a.out_dir, "tc_by_data_end_v2.csv"), ["data_end", "tc_penalised", "tc_unpenalised"], path)
    log("POINT replica=%s harmonise=%s P: %s %s | U: %s %s" % (
        replica_grade, harm_grade, fP["T1_forecast"]["grade"], fP["T3_shift"]["grade"],
        fU["T1_forecast"]["grade"], fU["T3_shift"]["grade"]))

    # ------------------------------------------------------------------ bootstraps (slow, last)
    bP20 = ddrep.bootstrap_dd23(p20, ts20, a.boot_p, a.seed + 11, nsub=a.nsub, log=log)
    log("P 2020 bootstrap: %s" % strip_tcs(bP20))
    bPlt = ddrep.bootstrap_dd23(plt, ts, a.boot_p, a.seed + 12, nsub=a.nsub, log=log)
    log("P latest bootstrap: %s" % strip_tcs(bPlt))
    bU20 = ddfit.bootstrap(u20, ts20, a.boot_u, a.seed + 13, nsub=a.nsub, log=log)
    log("U 2020 bootstrap: %s" % strip_tcs(bU20))
    bUlt = ddfit.bootstrap(ult, ts, a.boot_u, a.seed + 14, nsub=a.nsub, log=log)
    log("U latest bootstrap: %s" % strip_tcs(bUlt))
    out["P_penalised"]["bootstrap_2020"] = strip_tcs(bP20)
    out["P_penalised"]["bootstrap_latest"] = strip_tcs(bPlt)
    out["U_unpenalised"]["bootstrap_2020"] = strip_tcs(bU20)
    out["U_unpenalised"]["bootstrap_latest"] = strip_tcs(bUlt)
    out["complete"] = True
    out["runtime_s"] = time.time() - T_START
    with open(os.path.join(a.out_dir, "results_v2.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, indent=1, sort_keys=True)
    with open(os.path.join(a.out_dir, "bootstrap_tc_v2.json"), "w", encoding="utf-8") as f:
        json.dump({"P_2020": bP20["tcs"], "P_latest": bPlt["tcs"], "U_2020": bU20["tcs"],
                   "U_latest": bUlt["tcs"]}, f)
    q = ("tc_q025", "tc_q17", "tc_q50", "tc_q83", "tc_q975")
    write_csv(os.path.join(a.out_dir, "tc_summary_v2.csv"), ["fit", "tc"] + [k[3:] for k in q],
              [["P_2020", p20["tc"]] + [bP20[k] for k in q], ["P_latest", plt["tc"]] + [bPlt[k] for k in q],
               ["U_2020", u20["tc"]] + [bU20[k] for k in q], ["U_latest", ult["tc"]] + [bUlt[k] for k in q]])
    log("OUTCOME replica=%s harmonise=%s P: T1=%s T3=%s (%.1f -> %.1f) | U: T1=%s T3=%s (%.1f -> %.1f)" % (
        replica_grade, harm_grade, fP["T1_forecast"]["grade"], fP["T3_shift"]["grade"], p20["tc"], plt["tc"],
        fU["T1_forecast"]["grade"], fU["T3_shift"]["grade"], u20["tc"], ult["tc"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
