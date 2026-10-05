"""Maitra's rotating dust, stationary and cylindrically symmetric, with shear and no closed
timelike curve: .venv.noindex/bin/python -m unittest discover -s _tools

What its texts state, held on the published files, in units of its length a: the twist k is less
than r at every radius, so every circle about the axis is spacelike and g^tt is negative; the
density is greatest on the axis, c^2/8 pi G a^2, and falls with distance; on a cylinder of t and
phi the null curves are c dt = (k +- r) dphi, one rising and one falling toward +phi, so no cone
has tipped past the horizontal; the rays of no angular momentum run at
c dt/dr = +-e^(gamma/2) sqrt(1 - k^2/r^2); every cone of the figure in three dimensions is null
for the metric and rises in t all the way round; and the embedded plane z = 0 has the circles
sqrt(r^2 - k^2), which grow at every radius, at heights that make every distance the metric's.
The tests read the published files alone.
"""
import json
import math
import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "MFS" / "assets" / "data"


def root(r):
    return math.sqrt(1 + 4 * r * r)


def twist(r):
    """Maitra's k at a = 1."""
    s = root(r)
    return (s - 1 - math.log((s + 1) / 2)) / 2


def gamma(r):
    s = root(r)
    return 0.25 - 1 / (2 * (s + 1)) - math.log((s + 1) / 2) / 2


def density(r):
    """8 pi G rho a^2/c^2, the Ricci scalar the file publishes, at a = 1."""
    s = root(r)
    return 4 * math.exp(-gamma(r)) / (s * (s + 1) ** 2)


class Chart(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.metric = json.loads((DATA / "metrics" / "maitra_dust.json").read_text(encoding="utf-8"))
        cls.chart = cls.metric["coordinates"][0]

    def test_the_file_publishes_the_functions_these_tests_use(self):
        symbols = [p["symbol"] for p in self.chart["parameters"]]
        self.assertEqual(symbols[0], "a")
        self.assertIn("\\sqrt{1 + \\dfrac{4r^2}{a^2}}", symbols[1])
        self.assertIn("s - 1 - \\ln\\left(\\dfrac{s + 1}{2}\\right)", symbols[2])
        self.assertIn("\\dfrac{1}{4} - \\dfrac{1}{2\\left(s + 1\\right)} - \\dfrac{1}{2}\\ln", symbols[3])
        self.assertEqual(self.chart["ricci_scalar"], "R = \\dfrac{4e^{-\\gamma}}{a^2\\,s\\left(s + 1\\right)^2}")
        published = {tuple(e["indices"]): e["value"] for e in self.chart["inverse_metric_components"]}
        self.assertEqual(published["t", "t"], "-\\dfrac{r^2 - k^2}{r^2}")
        self.assertEqual({tuple(e["indices"]): e["value"] for e in self.chart["metric_components"]}["\\phi", "\\phi"],
                         "r^2 - k^2")

    def test_the_twist_is_less_than_the_radius_so_no_circle_is_null_and_t_is_a_time(self):
        for exponent in range(-30, 61):
            r = 10 ** (exponent / 10)
            k = twist(r)
            self.assertGreater(k, 0, r)
            self.assertLess(k, r, r)
            # g_phiphi = r^2 - k^2 and g^tt = -(r^2 - k^2)/r^2.
            self.assertGreater(r * r - k * k, 0, r)
            self.assertLess(-(r * r - k * k) / (r * r), 0, r)

    def test_the_slope_of_the_twist_is_the_stated_one_and_less_than_one(self):
        for r in (0.01, 0.3, 1.0, 5.0, 80.0, 3000.0):
            h = 1e-5 * r
            slope = (twist(r + h) - twist(r - h)) / (2 * h)
            stated = 2 * r / (root(r) + 1)
            self.assertAlmostEqual(slope, stated, delta=1e-7)
            self.assertLess(stated, 1)
            slope = (gamma(r + h) - gamma(r - h)) / (2 * h)
            self.assertAlmostEqual(slope, -2 * r / (root(r) + 1) ** 2, delta=1e-7)

    def test_the_density_is_greatest_on_the_axis_and_falls_with_distance(self):
        self.assertAlmostEqual(density(1e-6), 1.0, places=9)
        radii = [10 ** (e / 10) for e in range(-30, 41)]
        values = [density(r) for r in radii]
        self.assertTrue(all(b < a for a, b in zip(values, values[1:])))
        self.assertLess(values[-1], 1e-9)

    def test_on_the_axis_the_functions_begin_as_van_stockums_at_R_equal_to_2a(self):
        for r in (1e-3, 1e-2):
            self.assertAlmostEqual(twist(r) / (r * r / 2), 1, delta=10 * r * r)
            self.assertAlmostEqual(gamma(r) / (-r * r / 4), 1, delta=10 * r * r)

    def test_the_history_keeps_to_three_questions_and_five_paragraphs(self):
        history = self.metric["history"]
        self.assertLessEqual(history.count("?"), 3)
        self.assertEqual(len(history.split("¶")), 5)
        cited = set(k for group in re.findall(r"\[([a-z0-9, ]+)\]", history) for k in group.split(", "))
        self.assertLessEqual(cited, set(self.metric["references"]))


class Drawings(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        data = json.loads((DATA / "diagrams" / "maitra_dust.json").read_text(encoding="utf-8"))
        cls.views = {v["id"]: v for v in data["systems"]["cylindrical"]}
        cls.figure = data["projections"]["cylindrical"][0]
        cls.embedding = json.loads((DATA / "embedding" / "maitra_dust.json").read_text(encoding="utf-8"))

    def slopes(self, view):
        """d(ct)/d(r phi) of every ray from end to end, by family; a ray is straight, and its
        points are rounded to a ten thousandth of the drawing, so a short ray is passed over."""
        X0, X1, Y0, Y1 = view["box"]
        found = {}
        for family, lines in view["rays"].items():
            for line in lines:
                (ua, wa), (ub, wb) = line[0], line[-1]
                if abs(ub - ua) > 0.1:
                    found.setdefault(family, []).append((wb - wa) * (Y1 - Y0) / ((ub - ua) * (X1 - X0)))
        return found

    def test_the_rays_of_each_cylinder_are_the_null_lines_of_the_metric(self):
        for view, r in (("near", 1.0), ("far", 5.0)):
            k = twist(r)
            found = self.slopes(self.views[view])
            self.assertEqual(len(found), 2, view)
            want = {k / r + 1, k / r - 1}
            for family, slopes in found.items():
                nearest = min(want, key=lambda x: abs(x - slopes[0]))
                want.discard(nearest)
                for slope in slopes:
                    self.assertAlmostEqual(slope, nearest, delta=2e-3, msg=f"{view} {family}")
                # -(c dt - k dphi)^2 + r^2 dphi^2 = 0 with c dt = slope r dphi.
                self.assertAlmostEqual(-(nearest * r - k) ** 2 + r * r, 0.0, places=10)
            self.assertEqual(want, set(), view)

    def test_on_each_cylinder_one_family_rises_and_one_falls_so_no_cone_has_tipped_over(self):
        for view in ("near", "far"):
            found = self.slopes(self.views[view])
            self.assertEqual(sorted(slopes[0] > 0 for slopes in found.values()), [False, True], view)

    def test_the_numbers_the_captions_state(self):
        self.assertAlmostEqual(twist(1.0), 0.377, places=3)
        self.assertAlmostEqual(twist(5.0), 3.67, places=2)
        self.assertAlmostEqual((1 + twist(1.0)) / (1 - twist(1.0)), 2.2, delta=0.05)
        self.assertAlmostEqual((5 + twist(5.0)) / (5 - twist(5.0)), 6.5, delta=0.05)
        self.assertAlmostEqual(twist(2.0) / 2, 0.55, delta=0.005)
        self.assertAlmostEqual(twist(5.0) / 5, 0.73, delta=0.005)
        self.assertIn("0.377", " ".join(self.views["near"]["caption"]))
        self.assertIn("3.67", " ".join(self.views["far"]["caption"]))

    def test_the_rays_of_no_angular_momentum_run_at_the_stated_speed(self):
        view = self.views["radial"]
        X0, X1, Y0, Y1 = view["box"]
        checked = 0
        for lines in view["rays"].values():
            for line in lines:
                for (ua, wa), (ub, wb) in zip(line, line[1:]):
                    ra, rb = X0 + ua * (X1 - X0), X0 + ub * (X1 - X0)
                    if abs(rb - ra) < 0.05 or min(ra, rb) < 0.1:
                        continue
                    r = 0.5 * (ra + rb)
                    want = math.exp(gamma(r) / 2) * math.sqrt(1 - (twist(r) / r) ** 2)
                    got = abs((wb - wa) * (Y1 - Y0) / (rb - ra))
                    self.assertAlmostEqual(got, want, delta=0.02)
                    checked += 1
        self.assertGreater(checked, 50)

    def test_the_moment_embedded_is_the_line_t_equal_to_zero_of_each_view(self):
        for name, view in self.views.items():
            self.assertEqual(len(view["slices"]), 1, name)
            for line in view["slices"][0]["lines"]:
                for _, w in line:
                    self.assertAlmostEqual(w, 0.5, places=9, msg=name)

    def proper(self, r, n=2000):
        """The proper distance from the axis to r, the radius the figure is drawn with."""
        h = r / n
        return sum(math.exp(gamma((i + 0.5) * h) / 2) for i in range(n)) * h

    def test_every_cone_of_the_figure_is_null_and_rises_all_the_way_round(self):
        turn = self.figure["turn"]
        self.assertEqual({line["class"] for line in turn["lines"]}, {"axis", "floor"})
        self.assertEqual(len(turn["cones"]), 9)
        radius = {r: self.proper(r) for r in (2.0, 5.0)}
        for cone in turn["cones"]:
            X, Y, T = cone["apex"]
            drawn = math.hypot(X, Y)
            if drawn < 1e-3:
                r, k, stretch = 0.0, 0.0, 1.0
            else:
                r = min(radius, key=lambda x: abs(radius[x] - drawn))
                self.assertAlmostEqual(drawn, radius[r], delta=2e-3)
                k = twist(r)
                # The drawn radius is the proper distance, so a step along it is e^(gamma/2) dr.
                stretch = drawn / r
            for x, y, t in cone["rim"]:
                dX, dY, dT = x - X, y - Y, t - T
                self.assertGreater(dT, 0)
                size = dX * dX + dY * dY + dT * dT
                if r == 0.0:
                    null = -dT * dT + dX * dX + dY * dY
                else:
                    radial = (X * dX + Y * dY) / drawn
                    round_ = (X * dY - Y * dX) / drawn          # the drawn radius times dphi
                    dphi = round_ / drawn
                    null = -(dT - k * dphi) ** 2 + r * r * dphi * dphi + radial * radial
                self.assertAlmostEqual(null / size, 0.0, delta=3e-3)
            self.assertGreater(stretch, 0)

    def test_the_embedded_plane_has_the_circles_of_the_metric_which_grow_at_every_radius(self):
        view = self.embedding["views"][0]
        points = view["surfaces"][0]["pieces"][0]["points"]
        self.assertEqual(points[0], [0.0, 0.0, 0.0])
        self.assertAlmostEqual(points[-1][0], 6.0, places=12)
        for r, rho, _ in points:
            self.assertAlmostEqual(rho, math.sqrt(max(r * r - twist(r) ** 2, 0.0)), delta=2e-6)
        self.assertTrue(all(b[1] > a[1] and b[2] > a[2] for a, b in zip(points, points[1:])))
        self.assertAlmostEqual(next(rho for r, rho, _ in points if r == 1.0), 0.93, delta=0.005)
        self.assertAlmostEqual(points[-1][1], 3.87, delta=0.005)

    def test_every_distance_along_the_embedded_surface_is_the_metric_distance(self):
        points = self.embedding["views"][0]["surfaces"][0]["pieces"][0]["points"]
        along = sum(math.hypot(b[1] - a[1], b[2] - a[2]) for a, b in zip(points, points[1:]))
        self.assertAlmostEqual(along, self.proper(6.0, 20000), delta=2e-3)


if __name__ == "__main__":
    unittest.main()
