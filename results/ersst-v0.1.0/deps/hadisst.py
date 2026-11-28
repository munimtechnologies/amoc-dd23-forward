# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Munim, Inc. (MV2). Made with AI agents working in MV2's name.
"""Read HadISST1 ASCII files and build the DD23 AMOC fingerprint (SG - beta x GM).

HadISST1 ASCII (Met Office Read_instructions_sst.txt): integers, degC*100; 100% sea-ice
cells = -1000; land = -32768; a header with day, month, year starts each month; the
360x180 array runs west to east from 179.5W and north to south from 89.5N.
The subpolar-gyre (SG) cells are Caesar et al. (2018)'s sg_grid.txt (lat, lon centres).
Pure standard library.
"""
import gzip
import math
import re

INT_RE = re.compile(r"-?\d+")
NLON, NLAT = 360, 180
NCELL = NLON * NLAT
ICE = -1000
VARIANTS = ("w_excl", "w_ice", "u_excl", "u_ice")


def iter_months(lines):
    """Yield (year, month, values[64800]) from HadISST1 ASCII lines.

    The first three integers of a header are day, month, year. Data integers are read
    with a regex, so fixed-width fields that touch ("-32768-32768") and wrapped rows
    both parse.
    """
    header = None
    buf = []
    for line in lines:
        s = line.strip()
        if not s:
            continue
        if header is None:
            toks = INT_RE.findall(s)
            if len(toks) < 3:
                raise ValueError("bad HadISST header: %r" % s[:100])
            mon, yr = int(toks[1]), int(toks[2])
            if not (1 <= mon <= 12 and 1800 <= yr <= 2200):
                raise ValueError("bad HadISST header values: %r" % s[:100])
            header = (yr, mon)
            buf = []
            continue
        buf.extend(int(v) for v in INT_RE.findall(s))
        if len(buf) >= NCELL:
            if len(buf) > NCELL:
                raise ValueError("month %r has %d values" % (header, len(buf)))
            yield header[0], header[1], buf
            header = None
    if header is not None:
        raise ValueError("truncated month %r (%d values)" % (header, len(buf)))


def open_text(path):
    """Open a .gz text file; fall back to plain text if it is not really gzipped."""
    try:
        f = gzip.open(path, "rt", encoding="ascii", errors="replace")
        f.read(1)
        f.close()
        return gzip.open(path, "rt", encoding="ascii", errors="replace")
    except (OSError, EOFError):
        return open(path, "rt", encoding="ascii", errors="replace")


def sg_cells_from_grid(text):
    """Caesar et al. 2018 sg_grid.txt (lat lon per line, '#' header) -> sorted cell indices."""
    cells = set()
    for line in text.splitlines():
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        p = s.split()
        lat, lon = float(p[0]), float(p[1])
        r = int(round(89.5 - lat))
        c = int(round(lon + 179.5))
        if not (0 <= r < NLAT and 0 <= c < NLON):
            raise ValueError("cell outside grid: %r" % s)
        cells.add(r * NLON + c)
    return sorted(cells)


def row_weights():
    return [math.cos(math.radians(89.5 - r)) for r in range(NLAT)]


def _new_acc():
    return {"sw": 0.0, "w": 0.0, "sx": 0.0, "n": 0, "wi": 0.0, "ni": 0, "nbad": 0}


def month_stats(vals, sg_cells, weights):
    """Sums for the global ocean (GM) and the SG region, in degC*100 units."""
    g = _new_acc()
    for r in range(NLAT):
        w = weights[r]
        row = vals[r * NLON:(r + 1) * NLON]
        ocean = [v for v in row if v > -999]
        nice = row.count(ICE)
        nland = sum(1 for v in row if v <= -30000)
        so = sum(ocean)
        g["sw"] += w * so
        g["w"] += w * len(ocean)
        g["sx"] += so
        g["n"] += len(ocean)
        g["wi"] += w * nice
        g["ni"] += nice
        g["nbad"] += NLON - len(ocean) - nice - nland
    s = _new_acc()
    for idx in sg_cells:
        v = vals[idx]
        w = weights[idx // NLON]
        if v > -999:
            s["sw"] += w * v
            s["w"] += w
            s["sx"] += v
            s["n"] += 1
        elif v == ICE:
            s["wi"] += w
            s["ni"] += 1
    return {"gm": g, "sg": s}


def variant_means(acc):
    """Area-weighted or plain means, with 100%-ice cells excluded or set to -1.8 degC."""
    out = {}
    out["w_excl"] = acc["sw"] / acc["w"] / 100.0 if acc["w"] > 0 else None
    out["w_ice"] = ((acc["sw"] - 180.0 * acc["wi"]) / (acc["w"] + acc["wi"]) / 100.0
                    if acc["w"] + acc["wi"] > 0 else None)
    out["u_excl"] = acc["sx"] / acc["n"] / 100.0 if acc["n"] > 0 else None
    out["u_ice"] = ((acc["sx"] - 180.0 * acc["ni"]) / (acc["n"] + acc["ni"]) / 100.0
                    if acc["n"] + acc["ni"] > 0 else None)
    return out


def monthly_anomalies(keys, values, base=(1870, 2020)):
    """Subtract the calendar-month mean over base years (DD23 seasonal correction)."""
    sums = [0.0] * 13
    cnt = [0] * 13
    for (yr, mon), v in zip(keys, values):
        if base[0] <= yr <= base[1] and v is not None:
            sums[mon] += v
            cnt[mon] += 1
    clim = [None] + [(sums[m] / cnt[m]) if cnt[m] else None for m in range(1, 13)]
    out = []
    for (yr, mon), v in zip(keys, values):
        out.append(None if (v is None or clim[mon] is None) else v - clim[mon])
    return out


def t_of(yr, mon):
    """DD23 time convention: January 1870 -> 1870.0833, December 1870 -> 1871.0."""
    return yr + mon / 12.0


def ym_of_t(t):
    k = int(round(12.0 * t)) - 12 * 1870 - 1
    return 1870 + k // 12, k % 12 + 1


def parse_dd23(text):
    """DD23 AMOCdata.txt: time AMOC0 AMOC1 AMOC2 GM (R write.table; row names tolerated)."""
    rows = []
    for line in text.splitlines():
        s = line.strip()
        if not s:
            continue
        p = [q.strip('"') for q in s.split()]
        if len(p) == 6:
            p = p[1:]
        if len(p) != 5:
            continue
        try:
            vals = [float(q) for q in p]
        except ValueError:
            continue          # header
        yr, mon = ym_of_t(vals[0])
        rows.append({"t": vals[0], "yr": yr, "mon": mon, "AMOC0": vals[1],
                     "AMOC1": vals[2], "AMOC2": vals[3], "GM": vals[4]})
    return rows


def corr_rms(a, b):
    pairs = [(x, y) for x, y in zip(a, b) if x is not None and y is not None]
    n = len(pairs)
    if n < 3:
        return None, None, n
    ma = math.fsum(x for x, _ in pairs) / n
    mb = math.fsum(y for _, y in pairs) / n
    sab = math.fsum((x - ma) * (y - mb) for x, y in pairs)
    saa = math.fsum((x - ma) ** 2 for x, _ in pairs)
    sbb = math.fsum((y - mb) ** 2 for _, y in pairs)
    r = sab / math.sqrt(saa * sbb) if saa > 0 and sbb > 0 else None
    rms = math.sqrt(math.fsum((x - y) ** 2 for x, y in pairs) / n)
    return r, rms, n


def next_month(yr, mon):
    return (yr + 1, 1) if mon == 12 else (yr, mon + 1)
