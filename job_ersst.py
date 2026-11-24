# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Munim, Inc. (MV2). Made with AI agents working in MV2's name.
"""amoc-dd23-forward, ERSSTv5 arm (PREREGISTRATION-ERSST.md, frozen 2026-11-24 before any ERSST
data was read by this study).

The same forward test as v2 (T1 forecast envelope, T2 one-step residuals, T3 tipping-year shift),
with the same estimator code (fetched unchanged from the v2 commit), on the DD23 fingerprint
SG - 2 x GM rebuilt from NOAA ERSSTv5 instead of HadISST1.
Primary window: January 1870 to September 2026 (the HadISST1 v2 window).  Pure standard library.
"""
import argparse
import json
import math
import os
import sys
import time
import urllib.request

import ersst
import fetch_deps

ERSST_BASE = "https://www.ncei.noaa.gov/pub/data/cmb/ersst/v5/ascii/"
SG_GRID_URL = "https://www.pik-potsdam.de/~caesar/AMOC_slowdown/sg_grid.txt"
DD23_URL = "https://erda.ku.dk/archives/afce4e1c3ac6f0db61c27cb45c2e9b14/AMOCdata.txt"
UA = "mv2-research-job/1.0 (open science; ERSSTv5 reader)"
PRIMARY_END = (2026, 9)
BENYAMI_ERSST_TC = 3047.7     # Ben-Yami et al. 2024, Table S2, SPG-2xGMT ERSSTv5, data to 2020
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


def load_ersst(data_dir, last_year, sg_cells):
    weights = ersst.row_weights()
    stats = {}
    n_checked = 0
    n_fail = 0
    fail_months = []
    for y in range(1870, last_year + 1):
        name = "ersst.v5.%d.asc" % y
        try:
            p = download(ERSST_BASE + name, os.path.join(data_dir, name), refresh=(y >= last_year - 1))
        except RuntimeError:
            if y == last_year:
                log("no file for %d yet" % y)
                break
            raise
        with open(p, "r", encoding="ascii", errors="replace") as f:
            lines = f.readlines()
        for yr, mon, vals in ersst.iter_year(lines, y):
            n_checked += 1
            if not ersst.layout_ok(vals):
                n_fail += 1
                fail_months.append("%d-%02d" % (yr, mon))
            stats[(yr, mon)] = ersst.month_stats(vals, sg_cells, weights)
        if y % 20 == 0 or y >= last_year - 1:
            log("%s read; %d months so far" % (name, len(stats)))
    return stats, n_checked, n_fail, fail_months


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="/data")
    ap.add_argument("--out-dir", default="/out")
    ap.add_argument("--boot-p", type=int, default=500)
    ap.add_argument("--nfore", type=int, default=2000)
    ap.add_argument("--nsub", type=int, default=10)
    ap.add_argument("--seed", type=int, default=20261126)
    ap.add_argument("--last-year", type=int, default=2026)
    a = ap.parse_args(argv)
    os.makedirs(a.data_dir, exist_ok=True)
    os.makedirs(a.out_dir, exist_ok=True)
    inputs = os.path.join(a.out_dir, "inputs")
    os.makedirs(inputs, exist_ok=True)

    deps = os.path.join(a.out_dir, "deps")
    fetch_deps.main(deps)
    sys.path.insert(0, deps)
    import ddfit
    import ddrep
    import grade
    import hadisst as hs

    sg_path = download(SG_GRID_URL, os.path.join(inputs, "sg_grid.txt"))
    dd_path = download(DD23_URL, os.path.join(inputs, "dd23_AMOCdata.txt"))
    with open(sg_path, "r", encoding="latin-1") as f:
        sg_cells = ersst.sg_cells_2deg(f.read())
    with open(dd_path, "r", encoding="latin-1") as f:
        dd = hs.parse_dd23(f.read())
    if sum(sg_cells.values()) < 500 or len(dd) < 1500:
        raise RuntimeError("unexpected SG grid or DD23 file")
    out = {"study": "amoc-dd23-forward", "arm": "ersstv5", "prereg": "PREREGISTRATION-ERSST.md",
           "estimator_commit": fetch_deps.COMMIT,
           "sg_2deg_cells": len(sg_cells), "sg_1deg_cells": sum(sg_cells.values())}

    # ------------------------------------------------------------------ E0. read gate
    stats, n_checked, n_fail, fail_months = load_ersst(a.data_dir, a.last_year, sg_cells)
    keys_all = sorted(stats)
    contiguous = all(hs.next_month(*k1) == k2 for k1, k2 in zip(keys_all, keys_all[1:]))
    read_grade = ersst.grade_read(n_checked, n_fail, keys_all[0], keys_all[-1], contiguous)
    sg_frac = [stats[k]["sg_frac"] for k in keys_all]
    out["read"] = {"grade": read_grade, "first": "%d-%02d" % keys_all[0], "last": "%d-%02d" % keys_all[-1],
                   "n_months": len(keys_all), "layout_fail": n_fail, "layout_fail_months": fail_months[:50],
                   "contiguous": contiguous, "sg_frac_min": min(sg_frac), "sg_frac_max": max(sg_frac)}
    log("READ %s: %s to %s, %d layout failures, SG data share %.3f-%.3f" % (
        read_grade, out["read"]["first"], out["read"]["last"], n_fail, min(sg_frac), max(sg_frac)))
    if read_grade != "READ_OK":
        out["headline"] = "NOT_RUN"
        with open(os.path.join(a.out_dir, "results_ersst.json"), "w", encoding="utf-8") as f:
            json.dump(out, f, indent=1, sort_keys=True)
        log("OUTCOME read=%s headline=NOT_RUN" % read_grade)
        return 0

    # ------------------------------------------------------------------ fingerprint, harmonised as v2
    sg_raw = [stats[k]["sg"] for k in keys_all]
    gm_raw = [stats[k]["gm"] for k in keys_all]
    sg_an = hs.monthly_anomalies(keys_all, sg_raw)
    gm_an = hs.monthly_anomalies(keys_all, gm_raw)
    pos = {k: i for i, k in enumerate(keys_all)}
    dd_idx = [pos.get((r["yr"], r["mon"])) for r in dd]

    def aligned(series):
        return [series[i] if i is not None else None for i in dd_idx]

    k1, c1, rms1, n1 = ddrep.affine_fit(aligned(sg_an), [r["AMOC0"] for r in dd])
    k2, c2, rms2, n2 = ddrep.affine_fit(aligned(gm_an), [r["GM"] for r in dd])
    sg_h = [c1 + k1 * v for v in sg_an]
    gm_h = [c2 + k2 * v for v in gm_an]

    def fp(beta, s=sg_h, g=gm_h):
        return [x - beta * y for x, y in zip(s, g)]

    xh_all = fp(2.0)
    xraw_all = fp(2.0, sg_an, gm_an)
    r_h, rms_h, n_h = hs.corr_rms(aligned(xh_all), [r["AMOC2"] for r in dd])
    out["harmonise"] = {"sg": {"k": k1, "c": c1, "rms_resid": rms1, "n": n1},
                        "gm": {"k": k2, "c": c2, "rms_resid": rms2, "n": n2},
                        "amoc2_r_vs_hadisst_dd23": r_h, "amoc2_rms_vs_hadisst_dd23": rms_h, "n": n_h}
    log("HARMONISE SG k=%.4f c=%.4f rms=%.4f; GM k=%.4f c=%.4f rms=%.4f; AMOC2 r=%.4f rms=%.4f" % (
        k1, c1, rms1, k2, c2, rms2, r_h, rms_h))

    ts_all = [hs.t_of(y, m) for (y, m) in keys_all]
    n_prim = max(i for i, k in enumerate(keys_all) if k <= PRIMARY_END) + 1
    keys = keys_all[:n_prim]
    ts, xh = ts_all[:n_prim], xh_all[:n_prim]
    i20 = max(i for i, k in enumerate(keys) if k[0] <= 2020)
    ts20, x20 = ts[:i20 + 1], xh[:i20 + 1]
    fore_t, fore_x = ts[i20 + 1:], xh[i20 + 1:]
    nf = len(fore_t)
    years = sorted(set(k[0] for k in keys[i20 + 1:]))
    groups = {"year_%d" % y: [j for j, k in enumerate(keys[i20 + 1:]) if k[0] == y] for y in years}
    groups["last_12_months"] = list(range(max(0, nf - 12), nf))
    out.update({"primary_last_month": "%d-%02d" % keys[-1], "n_out_of_sample_months": nf})

    # ------------------------------------------------------------------ E1-E3. forward test, P (primary)
    def forward(name, fit20, fitlt, ramp20, seed):
        t1 = ddfit.forecast_check(ramp20, ts20[-1], x20[-1], fore_t, fore_x, a.nfore, seed, nsub=a.nsub,
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

    p20 = ddrep.fit_dd23(ts20, x20, pen=ddrep.PEN_DD23, variant="corrected")
    plt = ddrep.fit_dd23(ts, xh, pen=ddrep.PEN_DD23, variant="corrected")
    fP = forward("P", p20, plt, ddrep.ramp_of(p20), a.seed + 1)
    out["P_penalised"] = dict(fP, fit_2020=p20, fit_latest=plt)
    head = ersst.grade_consistency(fP["T1_forecast"]["grade"], fP["T3_shift"]["grade"])
    conv = ersst.grade_convergence(p20["tc"], plt["tc"])
    out["headline"] = head
    out["convergence"] = {"grade": conv, "hadisst_v2_tc_2020": 2064.77, "hadisst_v2_tc_latest": 2069.99}
    log("HEADLINE %s; CONVERGENCE %s" % (head, conv))

    # ------------------------------------------------------------------ secondary (report only)
    sec = {}

    def safe(fn):
        try:
            return fn()
        except (ValueError, ZeroDivisionError, OverflowError) as e:
            return {"error": repr(e)}

    try:
        u20 = ddfit.fit(ts20, x20)
        ult = ddfit.fit(ts, xh)
        fU = forward("U", u20, ult, ddfit.ramp_from(u20), a.seed + 2)
        sec["U_unpenalised"] = dict(fU, fit_2020=u20, fit_latest=ult)
    except ValueError as e:
        u20 = ult = None
        sec["U_unpenalised"] = {"error": repr(e)}
    sec["P_original_2023_code_2020"] = safe(lambda: ddrep.fit_dd23(ts20, x20, pen=ddrep.PEN_DD23,
                                                                   variant="original"))
    sec["benyami_2024_ersst_tc"] = BENYAMI_ERSST_TC
    sec["P_all_months"] = safe(lambda: ddrep.fit_dd23(ts_all, xh_all, pen=ddrep.PEN_DD23))
    if isinstance(sec["P_all_months"], dict) and "tc" in sec["P_all_months"]:
        sec["P_all_months"]["last_month"] = "%d-%02d" % keys_all[-1]
    xraw = xraw_all[:n_prim]
    sec["P_unharmonised_degC"] = safe(lambda: {
        "tc_2020": ddrep.fit_dd23(ts20, xraw[:i20 + 1], pen=ddrep.PEN_DD23)["tc"],
        "tc_latest": ddrep.fit_dd23(ts, xraw, pen=ddrep.PEN_DD23)["tc"]})
    for beta in (1.0, 3.0):
        xb = fp(beta)[:n_prim]
        sec["P_beta_%d" % int(beta)] = safe(lambda xb=xb: {
            "tc_2020": ddrep.fit_dd23(ts20, xb[:i20 + 1], pen=ddrep.PEN_DD23)["tc"],
            "tc_latest": ddrep.fit_dd23(ts, xb, pen=ddrep.PEN_DD23)["tc"]})
    path = []
    for y in range(2020, keys[-1][0] + 1):
        j = max(i for i, k in enumerate(keys) if k[0] <= y)
        fy = safe(lambda j=j: ddrep.fit_dd23(ts[:j + 1], xh[:j + 1], pen=ddrep.PEN_DD23,
                                             start=(plt["tau"], plt["a"])))
        uy = safe(lambda j=j: ddfit.fit(ts[:j + 1], xh[:j + 1]))
        path.append(["%d-%02d" % keys[j], fy.get("tc"), uy.get("tc")])
        log("data to %d-%02d: P tc=%s, U tc=%s" % (keys[j][0], keys[j][1], fy.get("tc"), uy.get("tc")))
    sec["tc_by_data_end"] = path
    out["secondary"] = sec

    with open(os.path.join(a.out_dir, "results_ersst_point.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, indent=1, sort_keys=True)

    # ------------------------------------------------------------------ files
    ddmap = {(r["yr"], r["mon"]): r["AMOC2"] for r in dd}
    write_csv(os.path.join(a.out_dir, "fingerprint_ersst_monthly.csv"),
              ["month", "t", "amoc2_ersst_harmonised", "sg_anom_degC", "gm_anom_degC", "amoc2_dd23_hadisst"],
              [["%d-%02d" % k, ts_all[i], xh_all[i], sg_an[i], gm_an[i], ddmap.get(k)]
               for i, k in enumerate(keys_all)])
    ann = []
    for y in sorted(set(k[0] for k in keys_all)):
        idx = [i for i, k in enumerate(keys_all) if k[0] == y]
        m = lambda s: math.fsum(s[i] for i in idx) / len(idx)
        ann.append([y, len(idx), m(xh_all), m(sg_an), m(gm_an), m(sg_h), m(gm_h)])
    write_csv(os.path.join(a.out_dir, "fingerprint_ersst_annual.csv"),
              ["year", "n_months", "amoc2_ersst_harmonised", "sg_anom_degC", "gm_anom_degC",
               "sg_harmonised", "gm_harmonised"], ann)
    fan = []
    for y in years:
        g = fP["T1_forecast"]["groups"]["year_%d" % y]
        fan.append([y, g["n_months"], g["observed_mean"], g["sim_q05"], g["sim_q50"], g["sim_q95"], g["q"]])
    write_csv(os.path.join(a.out_dir, "forecast_fan_ersst_P.csv"),
              ["year", "n_months", "observed", "model_q05", "model_q50", "model_q95", "q"], fan)
    write_csv(os.path.join(a.out_dir, "tc_by_data_end_ersst.csv"),
              ["data_end", "tc_penalised", "tc_unpenalised"], path)
    log("POINT read=%s P: T1=%s T2=%s T3=%s (%.1f -> %.1f) headline=%s convergence=%s" % (
        read_grade, fP["T1_forecast"]["grade"], fP["T2_one_step"]["flags"], fP["T3_shift"]["grade"],
        p20["tc"], plt["tc"], head, conv))

    # ------------------------------------------------------------------ bootstraps (slow, last)
    if not p20["tau_at_upper_bound"]:
        bP20 = ddrep.bootstrap_dd23(p20, ts20, a.boot_p, a.seed + 11, nsub=a.nsub, log=log)
        out["P_penalised"]["bootstrap_2020"] = strip_tcs(bP20)
        log("P 2020 bootstrap: %s" % strip_tcs(bP20))
    if not plt["tau_at_upper_bound"]:
        bPlt = ddrep.bootstrap_dd23(plt, ts, a.boot_p, a.seed + 12, nsub=a.nsub, log=log)
        out["P_penalised"]["bootstrap_latest"] = strip_tcs(bPlt)
        log("P latest bootstrap: %s" % strip_tcs(bPlt))
    out["complete"] = True
    out["runtime_s"] = time.time() - T_START
    with open(os.path.join(a.out_dir, "results_ersst.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, indent=1, sort_keys=True)
    log("OUTCOME read=%s P: T1=%s T3=%s (%.1f -> %.1f) headline=%s convergence=%s" % (
        read_grade, fP["T1_forecast"]["grade"], fP["T3_shift"]["grade"], p20["tc"], plt["tc"], head, conv))
    return 0


if __name__ == "__main__":
    sys.exit(main())
