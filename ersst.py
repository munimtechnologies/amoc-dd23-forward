# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Munim, Inc. (MV2). Made with AI agents working in MV2's name.
"""Read NOAA ERSSTv5 yearly ASCII files and build the DD23 fingerprint terms (SG, GM).

ERSSTv5 ASCII (NCEI /pub/data/cmb/ersst/v5/ascii/Readme and read.ersst.f, read 2026-11-24):
one file per year, ersst.v5.YYYY.asc; per month 180 records, one per longitude (0, 2, ..., 358 E),
each holding 89 integers in i6 format for latitudes -88, -86, ..., 88; units 0.01 degC;
missing (land) = -9999.  So value k of a month is longitude k // 89, latitude k % 89.
The subpolar-gyre (SG) region is Caesar et al. (2018)'s 1-degree sg_grid.txt; each 1-degree cell
centre lies inside exactly one 2-degree ERSST cell, and a 2-degree cell is weighted by how many
SG 1-degree cells it holds (times cos latitude).
Pure standard library.
"""
import math

NLON, NLAT = 180, 89
NCELL = NLON * NLAT
MISSING_BELOW = -900          # read.ersst.f: valid if > -900


def idx(lat, lon):
    """Index of the 2-degree cell centred at (lat, lon), lat in -88..88, lon in 0..358 (even)."""
    return (int(lon) // 2) * NLAT + (int(lat) + 88) // 2


def lat_of_index(k):
    return -88.0 + 2.0 * (k % NLAT)


def iter_year(lines, year):
    """Yield (year, month, values[16020]) from the lines of one yearly ERSSTv5 ASCII file."""
    rows = [ln for ln in lines if ln.strip()]
    if len(rows) % NLON != 0 or not rows:
        raise ValueError("year %d: %d records is not a whole number of months" % (year, len(rows)))
    nmon = len(rows) // NLON
    if nmon > 12:
        raise ValueError("year %d: %d months" % (year, nmon))
    for m in range(nmon):
        vals = []
        for ln in rows[m * NLON:(m + 1) * NLON]:
            v = [int(q) for q in ln.split()] if "-" not in ln.replace(" -", "") else []
            if len(v) != NLAT and len(ln.rstrip("\r\n")) >= 6 * NLAT:
                v = [int(ln[6 * q:6 * q + 6]) for q in range(NLAT)]      # fixed-width i6 fallback
            if len(v) != NLAT:
                raise ValueError("year %d month %d: record with %d values" % (year, m + 1, len(v)))
            vals.extend(v)
        yield year, m + 1, vals


def sg_cells_2deg(text):
    """Caesar et al. sg_grid.txt (lat lon of 1-degree centres) -> {2-degree cell index: count}."""
    cells = {}
    for line in text.splitlines():
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        p = s.split()
        lat, lon = float(p[0]), float(p[1])
        i = int(round((lon % 360.0) / 2.0)) % NLON
        j = int(round((lat + 88.0) / 2.0))
        j = min(max(j, 0), NLAT - 1)
        k = i * NLAT + j
        cells[k] = cells.get(k, 0) + 1
    return cells


def row_weights():
    return [math.cos(math.radians(-88.0 + 2.0 * j)) for j in range(NLAT)]


def month_stats(vals, sg_cells, weights):
    """Area-weighted global mean (GM) and SG mean in degC, plus the SG weight share with data."""
    sw = 0.0
    w = 0.0
    for k, v in enumerate(vals):
        if v > MISSING_BELOW:
            wk = weights[k % NLAT]
            sw += wk * v
            w += wk
    gm = sw / w / 100.0 if w > 0 else None
    ssw = 0.0
    sww = 0.0
    tot = 0.0
    for k, c in sg_cells.items():
        wk = c * weights[k % NLAT]
        tot += wk
        v = vals[k]
        if v > MISSING_BELOW:
            ssw += wk * v
            sww += wk
    sg = ssw / sww / 100.0 if sww > 0 else None
    return {"gm": gm, "sg": sg, "sg_frac": (sww / tot) if tot > 0 else 0.0}


# Fixed layout check (frozen in PREREGISTRATION-ERSST.md): open-ocean cells must hold data and
# continental cells must be missing.  A transposed or shifted reading fails it.
OCEAN_POINTS = ((0, 330), (0, 180), (-40, 80), (50, 330))
LAND_POINTS = ((30, 80), (-80, 0), (40, 260), (0, 20))
WARM_POINT = (0, 180)


def layout_ok(vals):
    for lat, lon in OCEAN_POINTS:
        if not vals[idx(lat, lon)] > MISSING_BELOW:
            return False
    for lat, lon in LAND_POINTS:
        if vals[idx(lat, lon)] > MISSING_BELOW:
            return False
    v = vals[idx(*WARM_POINT)] / 100.0
    return 15.0 <= v <= 35.0


def grade_read(n_months_checked, n_layout_fail, first_key, last_key, contiguous):
    """READ_OK: contiguous 1870-01 to at least 2026-06 and >= 99% of months pass the layout check."""
    if not contiguous or first_key != (1870, 1) or last_key < (2026, 6) or n_months_checked == 0:
        return "READ_FAIL"
    if n_layout_fail > 0.01 * n_months_checked:
        return "READ_FAIL"
    return "READ_OK"


def grade_consistency(t1_ersst, t3_ersst, t1_had="WITHIN_FORECAST", t3_had="LATER"):
    """Headline: does ERSSTv5 give the same pre-registered T1 and T3 grades as HadISST1 (v2)?"""
    if t1_ersst in (None, "NOT_RUN"):
        return "NOT_RUN"
    if t1_ersst != t1_had:
        return "DIFFERENT"
    if t3_ersst == t3_had:
        return "CONSISTENT"
    return "T1_ONLY"


def grade_convergence(tc_e2020, tc_elatest, tc_h2020=2064.77, tc_hlatest=2069.99):
    """Does adding 2021-2026 bring the ERSST tipping year closer to the HadISST one?"""
    if tc_e2020 is None or tc_elatest is None:
        return "NOT_RUN"
    g0 = abs(tc_e2020 - tc_h2020)
    g1 = abs(tc_elatest - tc_hlatest)
    if g1 < 0.9 * g0:
        return "CONVERGING"
    if g1 > 1.1 * g0:
        return "DIVERGING"
    return "UNCHANGED"
