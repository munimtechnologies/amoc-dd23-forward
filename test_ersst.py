# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Munim, Inc. (MV2). Made with AI agents working in MV2's name.
import unittest

import ersst


def month_lines(fn):
    """180 records of 89 i6 integers; fn(i_lon, j_lat) -> int."""
    return ["".join("%6d" % fn(i, j) for j in range(ersst.NLAT)) for i in range(ersst.NLON)]


def blank():
    return [-9999] * ersst.NCELL


def good_field():
    v = blank()
    for lat, lon in ersst.OCEAN_POINTS:
        v[ersst.idx(lat, lon)] = 1500
    v[ersst.idx(*ersst.WARM_POINT)] = 2800
    return v


class TestParse(unittest.TestCase):
    def test_two_months_layout(self):
        lines = month_lines(lambda i, j: i * 100 + j) + month_lines(lambda i, j: -(i + j))
        months = list(ersst.iter_year(lines, 2001))
        self.assertEqual([(y, m) for y, m, _ in months], [(2001, 1), (2001, 2)])
        v1 = months[0][2]
        self.assertEqual(len(v1), ersst.NCELL)
        self.assertEqual(v1[5 * 89 + 7], 507)
        self.assertEqual(v1[ersst.idx(-88, 0)], 0)
        self.assertEqual(v1[ersst.idx(88, 358)], 179 * 100 + 88)
        self.assertEqual(months[1][2][ersst.idx(0, 2)], -(1 + 44))

    def test_missing_values_parse(self):
        lines = month_lines(lambda i, j: -9999 if j == 0 else -180)
        (_, _, v), = list(ersst.iter_year(lines, 1870))
        self.assertEqual(v[0], -9999)
        self.assertEqual(v[1], -180)

    def test_partial_month_rejected(self):
        lines = month_lines(lambda i, j: 0)[:100]
        with self.assertRaises(ValueError):
            list(ersst.iter_year(lines, 1999))

    def test_short_record_rejected(self):
        lines = month_lines(lambda i, j: 0)
        lines[3] = lines[3][:-6]
        with self.assertRaises(ValueError):
            list(ersst.iter_year(lines, 1999))


class TestRegion(unittest.TestCase):
    def test_sg_mapping(self):
        text = "# lat lon\n55.5 -40.5\n56.5 -39.5\n57.5 -41.5\n\n"
        cells = ersst.sg_cells_2deg(text)
        self.assertEqual(cells, {160 * 89 + 72: 2, 159 * 89 + 73: 1})
        self.assertEqual(ersst.lat_of_index(160 * 89 + 72), 56.0)

    def test_dateline_wraps(self):
        cells = ersst.sg_cells_2deg("10.5 -0.5\n")
        self.assertEqual(list(cells), [ersst.idx(10, 0)])

    def test_month_stats(self):
        w = ersst.row_weights()
        v = blank()
        k_eq = ersst.idx(0, 100)
        k_60 = ersst.idx(60, 100)
        k_60m = ersst.idx(60, 104)
        v[k_eq] = 1000
        v[k_60] = 400
        s = ersst.month_stats(v, {k_60: 2, k_eq: 1, k_60m: 3}, w)
        self.assertAlmostEqual(s["gm"], 8.0, places=9)
        self.assertAlmostEqual(s["sg"], 7.0, places=9)
        self.assertAlmostEqual(s["sg_frac"], 2.0 / 3.5, places=9)


class TestChecksAndGrades(unittest.TestCase):
    def test_layout(self):
        self.assertTrue(ersst.layout_ok(good_field()))
        self.assertFalse(ersst.layout_ok(blank()))
        v = good_field()
        v[ersst.idx(40, 260)] = 1200          # data on Kansas: wrong layout
        self.assertFalse(ersst.layout_ok(v))
        v = good_field()
        v[ersst.idx(*ersst.WARM_POINT)] = 500   # 5 degC on the equator
        self.assertFalse(ersst.layout_ok(v))

    def test_grade_read(self):
        self.assertEqual(ersst.grade_read(1880, 0, (1870, 1), (2026, 10), True), "READ_OK")
        self.assertEqual(ersst.grade_read(1880, 18, (1870, 1), (2026, 10), True), "READ_OK")
        self.assertEqual(ersst.grade_read(1880, 19, (1870, 1), (2026, 10), True), "READ_FAIL")
        self.assertEqual(ersst.grade_read(1880, 0, (1870, 1), (2026, 5), True), "READ_FAIL")
        self.assertEqual(ersst.grade_read(1880, 0, (1870, 2), (2026, 10), True), "READ_FAIL")
        self.assertEqual(ersst.grade_read(1880, 0, (1870, 1), (2026, 10), False), "READ_FAIL")

    def test_grade_consistency(self):
        self.assertEqual(ersst.grade_consistency("WITHIN_FORECAST", "LATER"), "CONSISTENT")
        self.assertEqual(ersst.grade_consistency("WITHIN_FORECAST", "STABLE"), "T1_ONLY")
        self.assertEqual(ersst.grade_consistency("WITHIN_FORECAST", "UNBOUNDED"), "T1_ONLY")
        self.assertEqual(ersst.grade_consistency("BELOW_FORECAST", "LATER"), "DIFFERENT")
        self.assertEqual(ersst.grade_consistency("NOT_RUN", "LATER"), "NOT_RUN")

    def test_grade_convergence(self):
        self.assertEqual(ersst.grade_convergence(3047.0, 2500.0), "CONVERGING")
        self.assertEqual(ersst.grade_convergence(3047.0, 3500.0), "DIVERGING")
        self.assertEqual(ersst.grade_convergence(3047.0, 3050.0), "UNCHANGED")
        self.assertEqual(ersst.grade_convergence(None, 3050.0), "NOT_RUN")


if __name__ == "__main__":
    unittest.main()
