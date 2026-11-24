# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Munim, Inc. (MV2). Made with AI agents working in MV2's name.
import math
import unittest

import hadisst as hs


def month_text(mon, yr, vals, per_line=360, header=None):
    lines = [header or ("%12d%12d%12d" % (1, mon, yr))]
    for i in range(0, len(vals), per_line):
        lines.append("".join("%6d" % v for v in vals[i:i + per_line]))
    return lines


class TestReader(unittest.TestCase):
    def setUp(self):
        v = [-32768] * hs.NCELL
        v[1] = -1000
        v[2] = 123
        v[3] = -180
        v[hs.NCELL - 1] = -5
        self.v = v

    def test_two_months_fixed_width(self):
        lines = month_text(1, 1870, self.v) + month_text(2, 1870, self.v)
        # fields touch: "-32768-32768" must still split
        self.assertIn("-32768-32768", lines[1])
        got = list(hs.iter_months(lines))
        self.assertEqual([(y, m) for y, m, _ in got], [(1870, 1), (1870, 2)])
        vals = got[0][2]
        self.assertEqual(len(vals), hs.NCELL)
        self.assertEqual(vals[:4], [-32768, -1000, 123, -180])
        self.assertEqual(vals[-1], -5)

    def test_wrapped_rows_and_text_header(self):
        lines = month_text(12, 2025, self.v, per_line=10,
                           header="     1    12  2025   180 rows   360 columns")
        got = list(hs.iter_months(lines))
        self.assertEqual(len(got), 1)
        self.assertEqual(got[0][:2], (2025, 12))
        self.assertEqual(got[0][2], self.v)

    def test_truncated(self):
        lines = month_text(1, 1870, self.v)[:50]
        with self.assertRaises(ValueError):
            list(hs.iter_months(lines))


class TestGrid(unittest.TestCase):
    def test_sg_cells(self):
        text = "# On a grid from -90S to 90N and -180W and 180E\n  65.50  -59.50\n  42.50  -28.50\n"
        cells = hs.sg_cells_from_grid(text)
        self.assertEqual(cells, [24 * 360 + 120, 47 * 360 + 151])

    def test_row_weights(self):
        w = hs.row_weights()
        self.assertAlmostEqual(w[0], math.cos(math.radians(89.5)))
        self.assertAlmostEqual(w[90], math.cos(math.radians(-0.5)))


class TestStats(unittest.TestCase):
    def test_month_stats_and_variants(self):
        v = [-32768] * hs.NCELL
        v[0] = 1000        # row 0 (89.5N): 10 degC
        v[360] = -1000     # row 1 (88.5N): sea ice
        v[720] = 2000      # row 2 (87.5N): 20 degC
        w = hs.row_weights()
        st = hs.month_stats(v, [720], w)
        g = hs.variant_means(st["gm"])
        self.assertAlmostEqual(g["u_excl"], 15.0)
        self.assertAlmostEqual(g["u_ice"], 9.4)
        w0, w1, w2 = w[0], w[1], w[2]
        self.assertAlmostEqual(g["w_excl"], (w0 * 1000 + w2 * 2000) / (w0 + w2) / 100.0)
        self.assertAlmostEqual(g["w_ice"], (w0 * 1000 + w2 * 2000 - 180 * w1) / (w0 + w1 + w2) / 100.0)
        self.assertEqual(st["gm"]["nbad"], 0)
        s = hs.variant_means(st["sg"])
        for k in hs.VARIANTS:
            self.assertAlmostEqual(s[k], 20.0)

    def test_anomalies(self):
        keys = [(1870, 1), (1870, 2), (1871, 1), (1871, 2), (2021, 1)]
        an = hs.monthly_anomalies(keys, [1.0, 2.0, 3.0, 4.0, 10.0])
        self.assertEqual(an, [-1.0, -1.0, 1.0, 1.0, 8.0])

    def test_corr_rms(self):
        r, rms, n = hs.corr_rms([1.0, 2.0, 3.0], [2.0, 3.0, 4.0])
        self.assertAlmostEqual(r, 1.0)
        self.assertAlmostEqual(rms, 1.0)
        self.assertEqual(n, 3)


class TestDD23(unittest.TestCase):
    def test_time_convention(self):
        self.assertEqual(hs.ym_of_t(1870.0833), (1870, 1))
        self.assertEqual(hs.ym_of_t(1871.0), (1870, 12))
        self.assertEqual(hs.ym_of_t(2018.0), (2017, 12))
        self.assertEqual(hs.ym_of_t(2021.0), (2020, 12))
        self.assertAlmostEqual(hs.t_of(2020, 12), 2021.0)
        self.assertEqual(hs.next_month(2020, 12), (2021, 1))

    def test_parse(self):
        text = ('"time" "AMOC0" "AMOC1" "AMOC2" "GM"\n'
                '1870.0833 0.16395064 0.25571298 0.34747531 -0.09176234\n'
                '"2" 2021 -0.091480104 -0.31887462 -0.54626913 0.227394516\n')
        rows = hs.parse_dd23(text)
        self.assertEqual(len(rows), 2)
        self.assertEqual((rows[0]["yr"], rows[0]["mon"]), (1870, 1))
        self.assertEqual((rows[1]["yr"], rows[1]["mon"]), (2020, 12))
        self.assertAlmostEqual(rows[1]["AMOC2"], rows[1]["AMOC0"] - 2 * rows[1]["GM"], places=6)


if __name__ == "__main__":
    unittest.main()
