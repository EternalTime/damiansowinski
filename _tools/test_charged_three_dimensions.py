"""The point charge in three dimensions:
.venv.noindex/bin/python -m unittest discover -s _tools

What its charts, texts and drawings state, held on the published files. The static chart is held to
the line element of Gott, Simon and Alpert's abstract and its published Einstein tensor to kT for the
charge's field; the two null charts to the static chart pulled back through c du = c dt - dr/f and
c dv = c dt + dr/f; each chart's own check to failing on a metric that is not the solution; the History
to the static chart's line element; the embedding diagram's sheets to the proper distances of the slice
t = 0; the conformal diagram to its hexagon; the curvature figure to K = 3R^2; and the spacetime
diagrams to the horizon r = r_c. The tests that read the mathematics need sympy and are skipped where it is
absent.
"""
import importlib.util
import json
import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "MFS" / "assets" / "data"
HAS_SYMPY = importlib.util.find_spec("sympy") is not None
sys.path.insert(0, str(Path(__file__).resolve().parent / "derivations"))
ID = "charged_three_dimensions"
METRIC = json.loads((DATA / "metrics" / f"{ID}.json").read_text(encoding="utf-8"))
KQ2 = math.pi          # every drawing is made at kQ^2 = pi with r_c = 1
LEVEL = math.exp(-2)   # g_rr = 1 on the slice t = 0 at r = e^{-2pi/kQ^2} r_c


def chart(system):
    return next(c for c in METRIC["coordinates"] if c["id"] == system)


def drawing(kind):
    return json.loads((DATA / kind / f"{ID}.json").read_text(encoding="utf-8"))


def f(r, kq2=KQ2, rc=1.0):
    return kq2 / (2 * math.pi) * math.log(rc / r)


@unittest.skipUnless(HAS_SYMPY, "needs sympy")
class Charts(unittest.TestCase):
    """Each chart read through the checker's own Reader at points."""

    @classmethod
    def setUpClass(cls):
        import sympy as sp
        import verify_metrics as vm
        cls.sp, cls.vm = sp, vm
        cls.readers = {}
        for c in METRIC["coordinates"]:
            R = vm.Reader(c["coords"], [p["symbol"] for p in c["parameters"]], ())
            cls.readers[c["id"]] = (R, c)

    def components(self, system, field="metric_components"):
        sp = self.sp
        R, c = self.readers[system]
        n = len(c["coords"])
        g = sp.zeros(n)
        entries = c[field] if isinstance(c[field], list) else c[field]["variants"]["ll"]["nonzero"]
        for m in entries:
            i, j = (c["coords"].index(x) for x in m["indices"])
            g[i, j] = R(m["value"])
        return R, c, g

    def at(self, R, c, value, r, k=0.5, Q=None, rc=1.0):
        sp = self.sp
        Q = math.sqrt(2.0 / k) if Q is None else Q        # kQ^2 = 2
        point = {R.symbol[c["coords"][0]]: 0.3, R.symbol["r"]: r, R.symbol["\\phi"]: 1.0,
                 R.parameters["k"]: k, R.parameters["Q"]: Q, R.parameters["r_c"]: rc}
        return float(sp.N(sp.sympify(value).subs(point), 30))

    def test_the_static_chart_is_the_abstracts_line_element(self):
        R, c, g = self.components("static")
        for r in (0.2, 0.7, 1.6, 4.0):
            kq2 = 2.0
            self.assertAlmostEqual(self.at(R, c, g[0, 0], r), -f(r, kq2), places=12)
            self.assertAlmostEqual(self.at(R, c, g[1, 1], r), 1 / f(r, kq2), places=9)
            self.assertAlmostEqual(self.at(R, c, g[2, 2], r), r * r, places=12)

    def test_the_einstein_tensor_is_k_times_the_stress_of_the_field(self):
        # G_ab = k (F_ac F_b^c - g_ab F^2/4)/2 pi with F_rt = Q/r: in the static chart, with f the lapse
        # squared, G_tt = kQ^2 f/(4 pi r^2), G_rr = -kQ^2/(4 pi r^2 f) and G_phiphi = kQ^2/(4 pi).
        R, c, G = self.components("static", "einstein_tensor")
        k, Q = 0.5, 2.0
        for r in (0.3, 0.8, 2.5):
            a = k * Q * Q
            lapse = f(r, a)
            self.assertAlmostEqual(self.at(R, c, G[0, 0], r, k, Q), a * lapse / (4 * math.pi * r * r), places=12)
            self.assertAlmostEqual(self.at(R, c, G[1, 1], r, k, Q), -a / (4 * math.pi * r * r * lapse), places=9)
            self.assertAlmostEqual(self.at(R, c, G[2, 2], r, k, Q), a / (4 * math.pi), places=12)

    def test_every_chart_publishes_the_same_scalars(self):
        for system in ("static", "eddington_finkelstein_outgoing", "eddington_finkelstein_ingoing"):
            R, c = self.readers[system]
            ricci = R(c["ricci_scalar"].split("=", 1)[1])
            K = R(c["kretschmann"].split("=", 1)[1])
            for r in (0.25, 1.0, 3.0):
                value = self.at(R, c, ricci, r, 0.5, 2.0)
                self.assertAlmostEqual(value, 0.5 * 4 / (2 * math.pi * r * r), places=12)
                self.assertAlmostEqual(self.at(R, c, K, r, 0.5, 2.0), 3 * value ** 2, places=10)

    def test_the_null_charts_are_the_static_chart_pulled_back(self):
        sp = self.sp
        _, _, gs = self.components("static")
        Rs, cs = self.readers["static"]
        for system, sign in (("eddington_finkelstein_outgoing", -1), ("eddington_finkelstein_ingoing", 1)):
            R, c, g = self.components(system)
            for r in (0.3, 0.9, 1.4, 6.0):
                static = sp.Matrix(3, 3, lambda i, j: self.at(Rs, cs, gs[i, j], r))
                lapse = -static[0, 0]
                jacobian = sp.Matrix([[1, -sign / lapse, 0], [0, 1, 0], [0, 0, 1]])
                pulled = jacobian.T * static * jacobian
                for i in range(3):
                    for j in range(3):
                        self.assertAlmostEqual(float(pulled[i, j]), self.at(R, c, g[i, j], r), places=9,
                                               msg=f"{system} at r = {r}, slot {i}{j}")

    def test_each_charts_check_fails_on_a_metric_that_is_not_the_solution(self):
        import chart_printer as cp
        import charged_three_dimensions as module
        for system in module.SYSTEMS:
            line = module.line_element(system, c="").replace("{2\\pi}", "{3\\pi}", 1)
            chart_ = cp.Chart([module.TIME[system], "r", "\\phi"], module.PARAMETERS, line)
            with self.assertRaises(AssertionError, msg=system):
                module.check(chart_, system)


class Texts(unittest.TestCase):

    def test_the_history_prints_the_static_chart_line_element(self):
        line = chart("static")["line_element"]
        self.assertIn(f"${line}$", METRIC["history"])

    def test_the_three_charts_name_one_horizon(self):
        for c in METRIC["coordinates"]:
            self.assertIn("r = r_c \\;\\text{(the horizon)}", c["domains"])


class Embedding(unittest.TestCase):
    """Each sheet against the proper distance of the slice t = 0, g_rr = 2/ln(r_c/r) at kQ^2 = pi."""

    def test_the_sheets_keep_proper_distance_and_meet_at_the_horizon(self):
        view = drawing("embedding")["views"][0]
        pieces = {p["id"]: p for p in view["surfaces"][0]["pieces"]}
        for pid, sign in (("static", -1), ("other_static", 1)):
            pts = pieces[pid]["points"]
            self.assertAlmostEqual(pts[0][0], LEVEL, places=12)
            self.assertAlmostEqual(pts[-1][0], 1.0, places=12)
            self.assertAlmostEqual(pts[-1][2], 0.0, places=6)
            for (x0, rho0, z0), (x1, rho1, z1) in zip(pts, pts[1:]):
                self.assertAlmostEqual(rho0, x0, places=6)
                if x1 < 0.99:
                    mid = (x0 + x1) / 2
                    proper = (x1 - x0) * math.sqrt(1 / f(mid))
                    self.assertLess(abs(math.hypot(rho1 - rho0, z1 - z0) / proper - 1), 2e-3)
            self.assertTrue(all(sign * p[2] >= -1e-9 for p in pts))

    def test_the_surface_lies_level_where_g_rr_is_one(self):
        pts = {p["id"]: p for p in drawing("embedding")["views"][0]["surfaces"][0]["pieces"]}["static"]["points"]
        (x0, _, z0), (x1, _, z1) = pts[0], pts[1]
        self.assertLess(abs((z1 - z0) / (x1 - x0)), 0.05)


class Conformal(unittest.TestCase):

    def test_the_spacetime_is_a_hexagon_with_the_charges_on_its_vertical_sides(self):
        for view in drawing("conformal")["views"]:
            region = next(l for l in view["layers"] if l["kind"] == "fill" and l["class"] == "region")
            corners = sorted(tuple(round(c, 4) for c in p) for p in region["points"])
            h = round(math.pi / 2, 4)
            self.assertEqual(corners, sorted([(-h, -h), (-h, h), (0.0, round(math.pi, 4)), (h, h), (h, -h),
                                              (0.0, round(-math.pi, 4))]))
            zigs = [l for l in view["layers"] if l["kind"] == "zig" and l["class"] == "singular"]
            self.assertEqual(len(zigs), 2)
            for z in zigs:
                xs = {round(p[0], 4) for p in z["points"]}
                self.assertEqual(len(xs), 1)
                self.assertAlmostEqual(abs(xs.pop()), math.pi / 2, places=4)

    def test_the_tortoise_coordinate_has_the_slope_one_over_f(self):
        import charged_three_dimensions_drawings as d
        for r in (0.05, 0.4, 0.9, 1.2, 3.0, 20.0):
            h = 1e-6 * r
            slope = (d.rstar(r + h) - d.rstar(r - h)) / (2 * h)
            self.assertLess(abs(float(slope) * f(r) - 1), 1e-6)
        self.assertEqual(float(d.rstar(0.0)), 0.0)


class Curvature(unittest.TestCase):

    def test_the_curves_are_the_published_scalars(self):
        view = drawing("curvature")["views"][0]
        for curve in view["curves"]:
            for r, value in curve["points"][::50]:
                want = 3 / (4 * r ** 4) if curve["class"] == "K" else 1 / (2 * r * r)
                self.assertLess(abs(value / want - 1), 1e-9)
        self.assertEqual(view["marks"], [{"class": "horizon", "x": 1.0}])


class Diagrams(unittest.TestCase):

    def test_every_view_marks_the_horizon_at_r_c(self):
        for system, views in drawing("diagrams")["systems"].items():
            for view in views:
                grr = next(m for m in view["markers"] if m["kind"] == "grr")
                x0, x1, *_ = view["box"]
                for line in grr["lines"]:
                    for X, _ in line:
                        self.assertAlmostEqual(x0 + X * (x1 - x0), 1.0, places=3, msg=system)
                self.assertTrue(any(m["kind"] == "singular" and "left" in m.get("edges", [])
                                    for m in view["markers"]), system)


if __name__ == "__main__":
    unittest.main()
