"""The black holes of string theory with three and four charges, held on the published files:
.venv.noindex/bin/python -m unittest discover -s _tools

The published line elements, with f and the harmonic functions written out, are held to what their
sources state: the area of the horizon, 2 pi^2 prod sqrt(r_0^2 + r_i^2) in five dimensions and
4 pi prod sqrt(r_0 + r_i) in four, which at r_0 = 0 is 2 pi^2 r_1 r_2 r_3 and 4 pi sqrt(r_1 r_2 r_3 r_4)
and vanishes when a charge is missing; the mass each chart's convention states, read off g_tt far
away; the extreme chart as the chart off extremality at r_0 = 0; the areal chart as the chart of
three equal charges in rho^2 = r^2 + r_q^2; and the throat of the extreme hole of four charges as
Bertotti and Robinson's metric of radius (r_1 r_2 r_3 r_4)^(1/4). Those need sympy and are skipped
where it is absent. The drawings are held from their numbers alone: the embedded surfaces to their
circumference radii and their throats, the tortoise coordinate to its slope, and each conformal
view to its chart.
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
HAS_NUMPY = importlib.util.find_spec("numpy") is not None and importlib.util.find_spec("scipy") is not None
sys.path.insert(0, str(Path(__file__).resolve().parent / "derivations"))
ID = "string_bh_three_four_charges"
CHARTS = ["five_charges", "five_extreme", "five_areal", "four_charges", "four_extreme"]
DRAWN = {"five_charges": (1.0, [0.5, 1.0, 1.5]), "five_extreme": (0.0, [0.5, 1.0, 2.0]),
         "four_charges": (1.0, [0.5, 1.0, 1.5, 2.0]), "four_extreme": (0.0, [0.5, 1.0, 1.5, 2.0])}


def read(folder):
    return json.loads((DATA / folder / f"{ID}.json").read_text(encoding="utf-8"))


def circumference(chart, r):
    """The circumference radius of the sphere r of a chart at the parameters it is drawn at."""
    _, q = DRAWN[chart]
    if chart.startswith("five"):
        return r * math.prod(1 + a * a / (r * r) for a in q) ** (1 / 6)
    return r * math.prod(1 + a / r for a in q) ** (1 / 4)


class TheFile(unittest.TestCase):
    def test_the_five_charts_in_order(self):
        self.assertEqual([c["id"] for c in read("metrics")["coordinates"]], CHARTS)

    def test_each_charge_has_a_harmonic_function(self):
        for chart in read("metrics")["coordinates"]:
            symbols = [p["symbol"] for p in chart["parameters"]]
            radii = [s for s in symbols if s.startswith("r_") and s not in ("r_0", "r_q")]
            functions = [s for s in symbols if s.startswith("H_")]
            if chart["id"] == "five_areal":
                self.assertEqual(symbols, ["r_0", "r_q"])
                continue
            self.assertEqual(len(radii), 3 if chart["id"].startswith("five") else 4, chart["id"])
            self.assertEqual(len(functions), len(radii), chart["id"])
            self.assertEqual(any(s.startswith("f =") for s in symbols), not chart["id"].endswith("extreme"))


@unittest.skipUnless(HAS_SYMPY, "sympy is not installed")
class TheMetrics(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import sympy
        import verify_metrics
        cls.sp, cls.vm = sympy, verify_metrics
        cls.charts = {}
        for entry in read("metrics")["coordinates"]:
            held = verify_metrics.HELD.get((ID, entry["id"]), ())
            reader = verify_metrics.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]], (), held=held)
            g = verify_metrics.metric_from_line_element(reader, entry["line_element"], entry["coords"])
            cls.charts[entry["id"]] = (reader, g.xreplace(reader.held).subs(reader.c, 1))

    def numbers(self, chart, r0=None):
        """The chart's metric at the parameters it is drawn at, a function of r alone on the equator."""
        sp = self.sp
        reader, g = self.charts[chart]
        r0_drawn, q = DRAWN[chart]
        at = {reader.parameters[f"r_{i + 1}"]: sp.Rational(str(a)) for i, a in enumerate(q)}
        if "r_0" in reader.parameters:
            at[reader.parameters["r_0"]] = sp.Rational(str(r0_drawn if r0 is None else r0))
        at.update({reader.symbol[c]: sp.pi / 2 for c in reader.coords[2:]})
        return reader, g.subs(at), reader.symbol["r"]

    def test_the_area_of_the_horizon(self):
        sp = self.sp
        for chart, (r0, q) in DRAWN.items():
            reader, g, r = self.numbers(chart)
            five = chart.startswith("five")
            radius = sp.sqrt(g[2, 2])
            at_horizon = sp.limit(radius, r, 0, "+") if r0 == 0 else radius.subs(r, 1)
            k = 2 if five else 1
            expected = math.prod(r0 ** k + a ** k for a in q) ** (1 / (6 if five else 4))
            self.assertAlmostEqual(float(at_horizon), expected, places=12, msg=chart)

    def test_the_extreme_area_vanishes_when_a_charge_is_missing(self):
        sp = self.sp
        for chart in ("five_extreme", "four_extreme"):
            reader, g = self.charts[chart]
            r = reader.symbol["r"]
            at = {reader.parameters["r_1"]: 0, **{reader.parameters[f"r_{i}"]: 1 for i in range(2, len(DRAWN[chart][1]) + 1)},
                  **{reader.symbol[c]: sp.pi / 2 for c in reader.coords[2:]}}
            self.assertEqual(sp.limit(g[2, 2].subs(at), r, 0, "+"), 0, chart)

    def test_the_mass_of_each_convention(self):
        sp = self.sp
        for chart, (r0, q) in DRAWN.items():
            reader, g, r = self.numbers(chart)
            five = chart.startswith("five")
            power = 2 if five else 1
            fall = sp.limit((1 + g[0, 0]) * r ** power, r, sp.oo)
            # 8GM/(3 pi c^2) = r_0^2 + (2/3) sum r_i^2 in five dimensions, 2GM/c^2 = r_0 + (1/2) sum r_i in four.
            expected = r0 ** 2 + 2 * sum(a * a for a in q) / 3 if five else r0 + sum(q) / 2
            self.assertAlmostEqual(float(fall), expected, places=12, msg=chart)

    def test_the_extreme_chart_is_r_0_equal_to_zero(self):
        for extreme, general in (("five_extreme", "five_charges"), ("four_extreme", "four_charges")):
            _, g, r = self.numbers(extreme)
            DRAWN_BACKUP = DRAWN[general]
            DRAWN[general] = DRAWN[extreme]
            try:
                _, full, rr = self.numbers(general, r0=0)
            finally:
                DRAWN[general] = DRAWN_BACKUP
            for i in range(g.shape[0]):
                for x in (0.3, 1.7):
                    self.assertAlmostEqual(float(g[i, i].subs(r, x)), float(full[i, i].subs(rr, x)), places=12)

    def test_the_areal_chart_is_three_equal_charges(self):
        sp = self.sp
        reader, g = self.charts["five_areal"]
        other, full = self.charts["five_charges"]
        r0, rq, rho = sp.Rational(3, 4), sp.Integer(1), sp.Rational(9, 5)
        mine = g.subs({reader.parameters["r_0"]: r0, reader.parameters["r_q"]: rq,
                       **{reader.symbol[c]: sp.pi / 2 for c in reader.coords[2:]}})
        r = other.symbol["r"]
        theirs = full.subs({other.parameters["r_0"]: r0, **{other.parameters[f"r_{i}"]: rq for i in (1, 2, 3)},
                            **{other.symbol[c]: sp.pi / 2 for c in other.coords[2:]}})
        x = sp.sqrt(rho ** 2 - rq ** 2)
        jacobian = rho / x                                     # dr/d rho
        for i in range(5):
            scale = jacobian ** 2 if i == 1 else 1
            self.assertAlmostEqual(float(mine[i, i].subs(reader.symbol["\\rho"], rho)), float(theirs[i, i].subs(r, x) * scale),
                                   places=12)

    def test_the_throat_of_four_extreme_charges_is_bertotti_robinson(self):
        sp = self.sp
        _, g, r = self.numbers("four_extreme")
        b2 = math.sqrt(math.prod(DRAWN["four_extreme"][1]))
        self.assertAlmostEqual(float(sp.limit(-g[0, 0] / r ** 2, r, 0, "+")), 1 / b2, places=12)
        self.assertAlmostEqual(float(sp.limit(g[1, 1] * r ** 2, r, 0, "+")), b2, places=12)
        self.assertAlmostEqual(float(sp.limit(g[2, 2], r, 0, "+")), b2, places=12)


class TheEmbeddings(unittest.TestCase):
    def test_one_view_for_each_chart(self):
        self.assertEqual([v["id"] for v in read("embedding")["views"]], CHARTS)

    def test_every_circle_has_its_circumference_radius(self):
        for view in read("embedding")["views"]:
            for surface in view["surfaces"]:
                for piece in surface["pieces"]:
                    for r, rho, _ in piece["points"]:
                        expected = r if view["id"] == "five_areal" else circumference(view["id"], r)
                        self.assertAlmostEqual(rho, expected, places=5, msg=view["id"])

    def test_the_throats(self):
        views = {v["id"]: v for v in read("embedding")["views"]}
        for chart, radius in (("five_charges", (1.25 * 2 * 3.25) ** (1 / 6)), ("four_charges", (1.5 * 2 * 2.5 * 3) ** 0.25),
                              ("five_areal", 1.25)):
            pieces = views[chart]["surfaces"][0]["pieces"]
            self.assertEqual(len(pieces), 2, chart)
            for piece in pieces:
                self.assertAlmostEqual(piece["points"][0][1], radius, places=5, msg=chart)
                self.assertEqual(piece["points"][0][2], 0.0)
        for chart, radius in (("five_extreme", 1.0), ("four_extreme", 1.5 ** 0.25)):
            points = views[chart]["surfaces"][0]["pieces"][0]["points"]
            self.assertAlmostEqual(points[0][0], 1 / 50, places=9)
            # Down the throat the circles close on the geometric mean of the charge radii from above.
            self.assertGreater(points[0][1], radius)
            self.assertLess(points[0][1], 1.03 * radius)
            self.assertLess(points[0][2], points[-1][2] - 3.0)


class TheOtherDrawings(unittest.TestCase):
    def test_one_spacetime_diagram_for_each_chart(self):
        systems = read("diagrams")["systems"]
        ids = [s["id"] if isinstance(s, dict) else s for s in systems]
        self.assertEqual(sorted(ids), sorted(CHARTS))

    def test_one_conformal_view_for_each_chart(self):
        views = read("conformal")["views"]
        self.assertEqual([v["id"] for v in views], CHARTS)
        self.assertEqual([v["system"] for v in views], CHARTS)
        for v in views:
            self.assertTrue(v["caption"] and v["settings"], v["id"])


@unittest.skipUnless(HAS_NUMPY, "numpy and scipy are not installed")
class TheTortoiseCoordinate(unittest.TestCase):
    def test_its_slope_is_the_ratio_of_the_metric_components(self):
        import slices
        for chart, (r0, q) in DRAWN.items():
            for r in (0.3, 0.8, 2.5):
                h = 1e-5
                slope = (float(slices.sbc_rstar(chart, r + h)) - float(slices.sbc_rstar(chart, r - h))) / (2 * h)
                if chart.startswith("five"):
                    product = math.sqrt(math.prod(1 + a * a / (r * r) for a in q))
                    f = 1 - r0 * r0 / (r * r)
                else:
                    product = math.sqrt(math.prod(1 + a / r for a in q))
                    f = 1 - r0 / r
                self.assertAlmostEqual(slope / (product / f), 1.0, places=7, msg=chart)

    def test_the_surface_gravity(self):
        import slices
        self.assertAlmostEqual(slices.sbc_kappa("five_charges"), 1 / math.sqrt(1.25 * 2 * 3.25), places=12)
        self.assertAlmostEqual(slices.sbc_kappa("four_charges"), 1 / (2 * math.sqrt(1.5 * 2 * 2.5 * 3)), places=12)
        self.assertEqual(slices.sbc_kappa("five_extreme"), 0.0)
        self.assertEqual(slices.sbc_kappa("four_extreme"), 0.0)

    def test_the_areal_chart(self):
        import slices
        for rho in (0.5, 1.1, 3.0):
            h = 1e-6
            slope = (float(slices.sbc_areal_rstar(rho + h)) - float(slices.sbc_areal_rstar(rho - h))) / (2 * h)
            self.assertAlmostEqual(slope * (1 - 1.5625 / rho ** 2) * (1 - 1 / rho ** 2), 1.0, places=7)


if __name__ == "__main__":
    unittest.main()
