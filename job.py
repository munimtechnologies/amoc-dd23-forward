# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Munim, Inc. (MV2). Made with AI agents working in MV2's name.
"""amoc-dd23-forward: out-of-sample test and refit of the Ditlevsen & Ditlevsen (2023)
AMOC tipping forecast with HadISST1 months that were not in it (2021-01 onward).

Steps (all frozen in PREREGISTRATION.md):
  1. download HadISST1 ASCII (1870-latest), Caesar et al. 2018 SG grid, DD23 AMOCdata.txt
  2. rebuild SG and GM monthly means (4 variants each), anomalies vs 1870-2020 months
  3. pick the SG and GM variants closest to DD23's file on 1870-2020 (gate G)
  4. control C: our estimator on DD23's own series -> must give their corrected 2065
  5. M20 = fit on our rebuild through 2020-12; Mlatest = fit through the last month
  6. T1 forecast check of 2021-latest against M20; T2 one-step residuals; T3 shift of t_c
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
import grade
import hadisst as hs

HADISST_BASE = "https://www.metoffice.gov.uk/hadobs/hadisst/data/"
SG_GRID_URL = "https://www.pik-potsdam.de/~caesar/AMOC_slowdown/sg_grid.txt"
DD23_URL = "https://erda.ku.dk/archives/afce4e1c3ac6f0db61c27cb45c2e9b14/AMOCdata.txt"
UA = "mv2-research-job/1.0 (open science; HadISST reader)"

T_START = time.time()


def log(msg):
    print("[%8.1fs] %s" % (time.time() - T_START, msg), flush=True)


def hadisst_files(last_year):
    files = ["HadISST1_SST_1870-1900.txt.gz", "HadISST1_SST_1901-1930.txt.gz",
             "HadISST1_SST_1931-1960.txt.gz", "HadISST1_SST_1961-1990.txt.gz",
             "HadISST1_SST_1991-2003.txt.gz"]
    files += ["HadISST1_SST_%d.txt.gz" % y for y in range(2004, last_year + 1)]
    return files


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
        except Exception as e:  # network errors are retried, then raised
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


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="/data")
    ap.add_argument("--out-dir", default="/out")
    ap.add_argument("--boot", type=int, default=1000)
    ap.add_argument("--boot-control", type=int, default=300)
    ap.add_argument("--nfore", type=int, default=2000)
    ap.add_argument("--nsub", type=int, default=10)
    ap.add_argument("--seed", type=int, default=20261124)
    ap.add_argument("--last-year", type=int, default=2026)
    a = ap.parse_args(argv)
    os.makedirs(a.data_dir, exist_ok=True)
    os.makedirs(a.out_dir, exist_ok=True)
    inputs = os.path.join(a.out_dir, "inputs")
    os.makedirs(inputs, exist_ok=True)

    # ------------------------------------------------------------------ 1. data
    sg_path = download(SG_GRID_URL, os.path.join(inputs, "sg_grid.txt"))
    dd_path = download(DD23_URL, os.path.join(inputs, "dd23_AMOCdata.txt"))
    with open(sg_path, "r", encoding="latin-1") as f:
        sg_cells = hs.sg_cells_from_grid(f.read())
    with open(dd_path, "r", encoding="latin-1") as f:
        dd = hs.parse_dd23(f.read())
    log("SG cells: %d; DD23 rows: %d (%s .. %s)" % (len(sg_cells), len(dd), dd[0]["t"], dd[-1]["t"]))
    if len(sg_cells) < 500 or len(dd) < 1500:
        raise RuntimeError("unexpected SG grid or DD23 file")

    weights = hs.row_weights()
    stats = {}
    dup = 0
    sources = {}
    names = hadisst_files(a.last_year)
    for name in names:
        refresh = name.endswith("%d.txt.gz" % a.last_year)
        p = download(HADISST_BASE + name, os.path.join(a.data_dir, name), refresh=refresh)
        n0 = len(stats)
        with hs.open_text(p) as f:
            for yr, mon, vals in hs.iter_months(f):
                if (yr, mon) in stats:
                    dup += 1
                    continue
                stats[(yr, mon)] = hs.month_stats(vals, sg_cells, weights)
                sources[(yr, mon)] = name
        log("%s: %d months" % (name, len(stats) - n0))
    # newest month, if the update file is ahead of the yearly file (optional)
    try:
        p = download(HADISST_BASE + "HadISST1_SST_update.txt.gz",
                     os.path.join(a.data_dir, "HadISST1_SST_update.txt.gz"), refresh=True, tries=2)
        with hs.open_text(p) as f:
            for yr, mon, vals in hs.iter_months(f):
                if (yr, mon) not in stats:
                    stats[(yr, mon)] = hs.month_stats(vals, sg_cells, weights)
                    sources[(yr, mon)] = "HadISST1_SST_update.txt.gz"
                    log("update file added %d-%02d" % (yr, mon))
    except Exception as e:
        log("update file skipped: %r" % (e,))

    keys = sorted(stats)
    if keys[0] != (1870, 1):
        raise RuntimeError("record does not start in 1870-01: %r" % (keys[0],))
    for k1, k2 in zip(keys, keys[1:]):
        if hs.next_month(*k1) != k2:
            raise RuntimeError("gap in HadISST record after %r" % (k1,))
    last = keys[-1]
    log("HadISST months: %d, 1870-01 .. %d-%02d, duplicates skipped %d" % (len(keys), last[0], last[1], dup))

    # ------------------------------------------------------------------ 2. variants and anomalies
    sg_raw = {v: [hs.variant_means(stats[k]["sg"])[v] for k in keys] for v in hs.VARIANTS}
    gm_raw = {v: [hs.variant_means(stats[k]["gm"])[v] for k in keys] for v in hs.VARIANTS}
    sg_an = {v: hs.monthly_anomalies(keys, sg_raw[v]) for v in hs.VARIANTS}
    gm_an = {v: hs.monthly_anomalies(keys, gm_raw[v]) for v in hs.VARIANTS}
    sg_full = {v: hs.monthly_anomalies(keys, sg_raw[v], base=(1870, last[0])) for v in hs.VARIANTS}
    gm_full = {v: hs.monthly_anomalies(keys, gm_raw[v], base=(1870, last[0])) for v in hs.VARIANTS}

    # ------------------------------------------------------------------ 3. gate G
    pos = {k: i for i, k in enumerate(keys)}
    dd_idx = [pos.get((r["yr"], r["mon"])) for r in dd]
    dd_t = [r["t"] for r in dd]

    def aligned(series):
        return [series[i] if i is not None else None for i in dd_idx]

    cmp_sg = {}
    for v in hs.VARIANTS:
        r, rms, n = hs.corr_rms(aligned(sg_an[v]), [x["AMOC0"] for x in dd])
        cmp_sg[v] = {"r": r, "rms": rms, "n": n}
    cmp_gm = {}
    for v in hs.VARIANTS:
        r, rms, n = hs.corr_rms(aligned(gm_an[v]), [x["GM"] for x in dd])
        cmp_gm[v] = {"r": r, "rms": rms, "n": n}
    sel_sg = min(hs.VARIANTS, key=lambda v: cmp_sg[v]["rms"] if cmp_sg[v]["rms"] is not None else 9e9)
    sel_gm = min(hs.VARIANTS, key=lambda v: cmp_gm[v]["rms"] if cmp_gm[v]["rms"] is not None else 9e9)

    def fingerprint(beta, sg=sg_an, gm=gm_an):
        return [None if (s is None or g is None) else s - beta * g
                for s, g in zip(sg[sel_sg], gm[sel_gm])]

    x2 = fingerprint(2.0)
    r2, rms2, n2 = hs.corr_rms(aligned(x2), [x["AMOC2"] for x in dd])
    gate = grade.grade_recon(r2, rms2)
    log("gate: SG=%s GM=%s r=%s rms=%s -> %s" % (sel_sg, sel_gm, r2, rms2, gate))

    ts = [hs.t_of(y, m) for (y, m) in keys]
    if any(v is None for v in x2):
        raise RuntimeError("missing values in the rebuilt fingerprint")
    i20 = max(i for i, k in enumerate(keys) if k[0] <= 2020)
    ts20, x20 = ts[:i20 + 1], x2[:i20 + 1]
    fore_t, fore_x = ts[i20 + 1:], x2[i20 + 1:]
    log("through-2020 months: %d; out-of-sample months: %d" % (len(ts20), len(fore_t)))
    if len(fore_t) < 12:
        raise RuntimeError("fewer than 12 out-of-sample months")

    # ------------------------------------------------------------------ 4. control C
    dd_x = [x["AMOC2"] for x in dd]
    fc = ddfit.fit(dd_t, dd_x)
    log("control fit on DD23 file: tc=%.2f tau=%.2f A=%.4f" % (fc["tc"], fc["tau"], fc["A"]))
    control = grade.grade_control(fc["tc"])

    # ------------------------------------------------------------------ 5. M20 and Mlatest (points)
    m20 = ddfit.fit(ts20, x20)
    log("M20: tc=%.2f tau=%.2f A=%.4f" % (m20["tc"], m20["tau"], m20["A"]))
    mlt = ddfit.fit(ts, x2)
    log("Mlatest: tc=%.2f tau=%.2f A=%.4f" % (mlt["tc"], mlt["tau"], mlt["A"]))

    # ------------------------------------------------------------------ 6. T1, T2, T3
    r20 = ddfit.ramp_from(m20)
    years = sorted(set(k[0] for k in keys[i20 + 1:]))
    groups = {}
    for y in years:
        groups["year_%d" % y] = [j for j, k in enumerate(keys[i20 + 1:]) if k[0] == y]
    nf = len(fore_t)
    groups["last_12_months"] = list(range(max(0, nf - 12), nf))
    t1 = ddfit.forecast_check(r20, ts20[-1], x20[-1], fore_t, fore_x, a.nfore, a.seed + 4,
                              nsub=a.nsub, groups=groups)
    t1_grade = grade.grade_forecast(t1["groups"]["all"]["q"])
    res = ddfit.one_step_residuals(r20, [ts20[-1]] + fore_t, [x20[-1]] + fore_x)
    t2 = ddfit.residual_summary(res)
    t2_flags = grade.grade_onestep(t2)
    # in-sample residuals of M20 for comparison
    i0 = next(i for i, t in enumerate(ts20) if t >= ddfit.T0)
    t2_in = ddfit.residual_summary(ddfit.one_step_residuals(r20, ts20[i0 - 1:], x20[i0 - 1:]))
    unb = m20["tau_at_upper_bound"] or mlt["tau_at_upper_bound"]
    t3_grade = grade.grade_shift(m20["tc"], mlt["tc"], unbounded=unb)
    log("T1 q=%.3f -> %s; T2 %s; T3 dtc=%.2f -> %s" % (
        t1["groups"]["all"]["q"], t1_grade, t2_flags, mlt["tc"] - m20["tc"], t3_grade))

    # ------------------------------------------------------------------ secondary (report only)
    secondary = {}
    for beta in (1.0, 3.0):
        xb = fingerprint(beta)
        try:
            f20 = ddfit.fit(ts20, xb[:i20 + 1])
            flt = ddfit.fit(ts, xb)
            secondary["beta_%d" % int(beta)] = {"tc_2020": f20["tc"], "tc_latest": flt["tc"],
                                                "shift": flt["tc"] - f20["tc"],
                                                "tau_bound": f20["tau_at_upper_bound"] or flt["tau_at_upper_bound"]}
        except ValueError as e:
            secondary["beta_%d" % int(beta)] = {"error": str(e)}
    xfull = fingerprint(2.0, sg_full, gm_full)
    try:
        ff = ddfit.fit(ts, xfull)
        secondary["full_record_climatology"] = {"tc_latest": ff["tc"], "shift_vs_M20": ff["tc"] - m20["tc"]}
    except ValueError as e:
        secondary["full_record_climatology"] = {"error": str(e)}
    try:
        fc1 = ddfit.fit(dd_t, [x["AMOC1"] for x in dd])
        secondary["control_beta_1_dd23_file"] = {"tc": fc1["tc"], "dd23_table1": 2123}
    except ValueError as e:
        secondary["control_beta_1_dd23_file"] = {"error": str(e)}

    # ------------------------------------------------------------------ outputs
    def build(bc, b20, blt):
        def q(b, key):
            return b[key] if b else None
        return {
            "study": "amoc-dd23-forward",
            "preregistration": "PREREGISTRATION.md (committed 2026-11-24, before this job ran)",
            "complete": bool(bc and b20 and blt),
            "hadisst_last_month": "%d-%02d" % last,
            "n_months": len(keys), "n_out_of_sample_months": nf,
            "gate": {"grade": gate, "sg_variant": sel_sg, "gm_variant": sel_gm,
                     "amoc2_r": r2, "amoc2_rms": rms2, "n": n2, "sg_variants": cmp_sg, "gm_variants": cmp_gm},
            "control": {"grade": control, "fit": fc, "bootstrap": strip_tcs(bc) if bc else None,
                        "dd23_tc": grade.DD23_TC, "dd23_ci95": list(grade.DD23_CI95)},
            "M20": {"fit": m20, "bootstrap": strip_tcs(b20) if b20 else None},
            "Mlatest": {"fit": mlt, "bootstrap": strip_tcs(blt) if blt else None},
            "T1_forecast": dict(t1, grade=t1_grade),
            "T2_one_step": {"flags": t2_flags, "out_of_sample": t2, "in_sample_1924_2020": t2_in},
            "T3_shift": {"grade": t3_grade, "tc_2020": m20["tc"], "tc_latest": mlt["tc"],
                         "shift_years": mlt["tc"] - m20["tc"],
                         "ci95_2020": [q(b20, "tc_q025"), q(b20, "tc_q975")],
                         "ci95_latest": [q(blt, "tc_q025"), q(blt, "tc_q975")]},
            "secondary": secondary,
            "runtime_s": time.time() - T_START,
        }

    with open(os.path.join(a.out_dir, "results_point.json"), "w", encoding="utf-8") as f:
        json.dump(build(None, None, None), f, indent=1, sort_keys=True)
    log("POINT gate=%s control=%s T1=%s T2=%s T3=%s (tc %.1f -> %.1f)" % (
        gate, control, t1_grade, ",".join(t2_flags), t3_grade, m20["tc"], mlt["tc"]))

    # monthly series
    rows = []
    ddmap = {(r["yr"], r["mon"]): r["AMOC2"] for r in dd}
    for i, k in enumerate(keys):
        rows.append(["%d-%02d" % k, ts[i], sg_an[sel_sg][i], gm_an[sel_gm][i], x2[i], ddmap.get(k),
                     sources[k]])
    write_csv(os.path.join(a.out_dir, "fingerprint_monthly.csv"),
              ["month", "t", "sg_anom", "gm_anom", "amoc2_ours", "amoc2_dd23", "source"], rows)
    # annual means
    ann = []
    rlt = ddfit.ramp_from(mlt)
    for y in sorted(set(k[0] for k in keys)):
        idx = [i for i, k in enumerate(keys) if k[0] == y]
        ours = math.fsum(x2[i] for i in idx) / len(idx)
        ddv = [ddmap[keys[i]] for i in idx if keys[i] in ddmap]
        ann.append([y, len(idx), ours, (math.fsum(ddv) / len(ddv)) if len(ddv) == len(idx) else None,
                    math.fsum(r20.mu(ts[i]) for i in idx) / len(idx),
                    math.fsum(rlt.mu(ts[i]) for i in idx) / len(idx)])
    write_csv(os.path.join(a.out_dir, "fingerprint_annual.csv"),
              ["year", "n_months", "amoc2_ours", "amoc2_dd23", "mu_M20", "mu_Mlatest"], ann)
    # forecast fan per year
    fan = []
    for y in years:
        g = t1["groups"]["year_%d" % y]
        fan.append([y, g["n_months"], g["observed_mean"], g["sim_q05"], g["sim_q50"], g["sim_q95"], g["q"]])
    write_csv(os.path.join(a.out_dir, "forecast_fan.csv"),
              ["year", "n_months", "observed", "model_q05", "model_q50", "model_q95", "q"], fan)

    # ------------------------------------------------------------------ bootstraps (slow, last)
    log("starting bootstraps")
    bc = ddfit.bootstrap(fc, dd_t, a.boot_control, a.seed + 1, nsub=a.nsub, log=log)
    log("control bootstrap: %s" % strip_tcs(bc))
    b20 = ddfit.bootstrap(m20, ts20, a.boot, a.seed + 2, nsub=a.nsub, log=log)
    log("M20 bootstrap: %s" % strip_tcs(b20))
    blt = ddfit.bootstrap(mlt, ts, a.boot, a.seed + 3, nsub=a.nsub, log=log)
    log("Mlatest bootstrap: %s" % strip_tcs(blt))
    with open(os.path.join(a.out_dir, "results.json"), "w", encoding="utf-8") as f:
        json.dump(build(bc, b20, blt), f, indent=1, sort_keys=True)
    # bootstrap distributions
    with open(os.path.join(a.out_dir, "bootstrap_tc.json"), "w", encoding="utf-8") as f:
        json.dump({"control": bc["tcs"], "M20": b20["tcs"], "Mlatest": blt["tcs"]}, f)
    write_csv(os.path.join(a.out_dir, "tc_summary.csv"),
              ["fit", "tc", "q025", "q17", "q50", "q83", "q975"],
              [["DD23_file_control", fc["tc"], bc["tc_q025"], bc["tc_q17"], bc["tc_q50"], bc["tc_q83"], bc["tc_q975"]],
               ["rebuild_to_2020", m20["tc"], b20["tc_q025"], b20["tc_q17"], b20["tc_q50"], b20["tc_q83"], b20["tc_q975"]],
               ["rebuild_to_latest", mlt["tc"], blt["tc_q025"], blt["tc_q17"], blt["tc_q50"], blt["tc_q83"], blt["tc_q975"]]])
    log("OUTCOME gate=%s control=%s T1=%s T2=%s T3=%s (tc %.1f -> %.1f)" % (
        gate, control, t1_grade, ",".join(t2_flags), t3_grade, m20["tc"], mlt["tc"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
