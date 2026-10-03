"""Bah and Heidmann's topological star:
python3 -m unittest discover -s _tools

What its texts and drawings state, held on the published files. The numbers its captions quote,
the times light takes, the clocks on the bubble, the curvature there, the periods of the swinging
particles and the length of the cigar, are computed again here by quadrature from the line
element alone, and the embedding diagram's surfaces are held to the slopes they are built from.
Where sympy is installed, every published chart is held to Einstein's equations with Maxwell's
field, to the published black string and Schwarzschild's metric at the two ends of the family,
and to Bah and Heidmann's chart under its map; those tests are skipped where it is absent.
"""
import importlib.util
import json
import math
import sys
import types
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "MFS" / "assets" / "data"
HAS_SYMPY = importlib.util.find_spec("sympy") is not None
sys.path.insert(0, str(Path(__file__).resolve().parent / "derivations"))

STAR, STRING = (0.75, 1.0), (1.0, 0.75)


def simpson(of, a, b, n=4000):
    h = (b - a) / n
    return h / 3 * (of(a) + of(b) + sum((4 if k % 2 else 2) * of(a + k * h) for k in range(1, n)))


def tortoise(r, rS, rB):
    """The integral of dr/((1 - r_S/r) sqrt(1 - r_B/r)) from r_B, with r = r_B + u^2, whose
    integrand 2 (u^2 + r_B)^(3/2)/(u^2 + r_B - r_S) is finite on the bubble."""
    return simpson(lambda u: 2 * (u * u + rB) ** 1.5 / (u * u + rB - rS), 0.0, math.sqrt(r - rB))


def kretschmann(r, rS, rB):
    return (48 * r * r * (rS * rS + rS * rB + rB * rB) - 120 * r * rS * rB * (rS + rB)
            + 115 * rS * rS * rB * rB) / (4 * r ** 8)


def captions(kind):
    return {(v["id"] if kind != "diagrams" else f"{system}/{v['id']}"): " ".join(v["caption"])
            for system, views in _views(kind) for v in views}


def _views(kind):
    data = json.loads((DATA / kind / "topological_star.json").read_text(encoding="utf-8"))
    if kind == "diagrams":
        return list(data["systems"].items())
    return [(None, data["views"])]


class Numbers(unittest.TestCase):
    def test_the_period_of_the_circle_makes_the_bubble_smooth(self):
        # A circle a small proper distance d from the bubble has the circumference 2 pi R_y sqrt(f_B),
        # and with R_y = 2 sqrt(r_B^3/(r_B - r_S)) that is 2 pi d.
        rS, rB = STAR
        Ry = 2 * math.sqrt(rB ** 3 / (rB - rS))
        self.assertAlmostEqual(Ry, 4.0, places=12)
        for eps in (1e-5, 1e-6):
            r = rB + eps
            d = simpson(lambda u: 2 * (u * u + rB) / math.sqrt(u * u + rB - rS), 0.0, math.sqrt(eps), 200)
            self.assertAlmostEqual(Ry * math.sqrt(1 - rB / r) / d, 1.0, places=4)

    def test_light_from_the_bubble_and_across_it(self):
        rS, rB = STAR
        out = tortoise(2.0, rS, rB)
        self.assertEqual(f"{out:.2f}", "5.92")
        self.assertEqual(f"{2 * out:.2f}", "11.84")
        text = captions("diagrams")
        self.assertIn("$ct = 5.92\\,r_B$", text["bah_heidmann/radial"])
        self.assertIn("$ct = 11.84\\,r_B$", text["bubble/through"])
        # rho = 4 is r = r_B + (r_B - r_S) rho^2/4 = 2 r_B.
        self.assertEqual(rB + (rB - rS) * 16 / 4, 2.0)

    def test_the_closed_form_of_the_tortoise_coordinate(self):
        import slices
        for (rS, rB), points in ((STAR, (1.3, 2.0, 7.0)), (STRING, (0.8, 0.95, 1.6, 7.0))):
            for r in points:
                if rS > rB and r > rS:
                    # Across the horizon the integral is taken from 1.6 r_S outward.
                    continue
                self.assertAlmostEqual(float(slices.topological_star_rstar(r, rS, rB)), tortoise(r, rS, rB), places=8)
            self.assertEqual(float(slices.topological_star_rstar(rB, rS, rB)), 0.0)
        a, b = 1.6, 7.0
        rS, rB = STRING
        quad = simpson(lambda r: 1 / ((1 - rS / r) * math.sqrt(1 - rB / r)), a, b)
        self.assertAlmostEqual(float(slices.topological_star_rstar(b, rS, rB) - slices.topological_star_rstar(a, rS, rB)),
                               quad, places=8)
        quad = simpson(lambda x: (1 + 1 / x) ** 1.5, 0.3, 2.5)
        self.assertAlmostEqual(float(slices.extremal_string_rstar(2.5) - slices.extremal_string_rstar(0.3)), quad, places=8)

    def test_clocks_and_curvature_on_the_bubble(self):
        rS, rB = STAR
        self.assertEqual(math.sqrt(1 - rS / rB), 0.5)
        # c dt/drho on the bubble is r_B^(3/2)/sqrt(r_B - r_S).
        self.assertAlmostEqual(rB ** 1.5 / math.sqrt(rB - rS), 2.0, places=14)
        self.assertEqual(f"{kretschmann(rB, rS, rB):.2f}", "4.55")
        self.assertEqual(kretschmann(rB, rS, rB), (48 * rB ** 2 - 72 * rB * rS + 43 * rS ** 2) / (4 * rB ** 6))
        # Greatest on the bubble.
        for r in (1.01, 1.5, 3.0, 10.0):
            self.assertLess(kretschmann(r, rS, rB), kretschmann(rB, rS, rB))
        rS, rB = STRING
        self.assertEqual(f"{kretschmann(rB, rS, rB):.1f}", "22.5")
        self.assertIn("$22.5/r_S^4$", captions("diagrams")["eddington_finkelstein_ingoing/finkelstein"])
        # The extremal string's scalar m^2 (144 rho^2 + 48 rho m + 19 m^2)/(4 (rho + m)^8) on its horizon.
        self.assertEqual(19 / 4, 4.75)
        self.assertIn("$19/(4m^4)$", captions("conformal")["extremal"])

    def test_the_surface_gravity_of_the_black_string(self):
        rS, rB = STRING
        kappa = math.sqrt(rS - rB) / (2 * rS ** 1.5)
        self.assertEqual(kappa, 0.25)
        # 1/(f_S sqrt f_B) near r_S is 1/(2 kappa (r - r_S)).
        r = rS + 1e-7
        self.assertAlmostEqual((r - rS) / ((1 - rS / r) * math.sqrt(1 - rB / r)), 1 / (2 * kappa), places=5)

    def test_the_swings_through_the_bubble(self):
        # Released from rest at rho_0, with N^2 = -g_tt and h = g_rhorho at r_B = 1, r_S = 3/4, the
        # period in ct is 4 times the integral of E sqrt(h)/(N^2 sqrt(E^2/N^2 - 1)), rho = rho_0 sin(a).
        N2 = lambda x: (4 + x * x) / (16 + x * x)
        h = lambda x: (16 + x * x) ** 2 / (64 * (4 + x * x))

        def period(x0):
            E2 = N2(x0)

            def of(a):
                x = x0 * math.sin(a)
                return math.sqrt(E2 * h(x)) / N2(x) / math.sqrt(E2 / N2(x) - 1) * x0 * math.cos(a)
            # The integrand is finite at the turning point, a = pi/2, where it is 0/0, so the
            # midpoint rule, which never evaluates an end.
            n = 4000
            return 4 * (math.pi / 2) / n * sum(of((k + 0.5) * (math.pi / 2) / n) for k in range(n))
        self.assertAlmostEqual(period(1e-3), 16 * math.pi / math.sqrt(3), places=3)
        self.assertEqual(f"{period(2.0):.1f}", "35.6")
        self.assertEqual(f"{period(4.0):.1f}", "52.8")
        self.assertIn("$ct = 35.6\\,r_B$ and $52.8\\,r_B$", captions("conformal")["bubble"])

    def test_mass_and_double_wick_rotation(self):
        # M = 2 pi (2 r_S + r_B)/kappa_4^2 with kappa_4^2 = 8 pi G/c^4 is c^2 (2 r_S + r_B)/(4G), and at
        # r_B = 0 it is Schwarzschild's c^2 r_S/(2G).
        for rS, rB in ((1.0, 0.0), (0.75, 1.0), (1.0, 0.75)):
            self.assertAlmostEqual(2 * math.pi * (2 * rS + rB) / (8 * math.pi), (2 * rS + rB) / 4, places=15)
        self.assertEqual((2 * 1.0 + 0.0) / 4, 1.0 / 2)


class Surface(unittest.TestCase):
    def setUp(self):
        views = {v["id"]: v for v in json.loads((DATA / "embedding" / "topological_star.json").read_text(
            encoding="utf-8"))["views"]}
        self.cigar = views["cigar"]["surfaces"][0]["pieces"][0]
        self.equator = {p["id"]: p for p in views["equator"]["surfaces"][0]["pieces"]}

    def test_the_cigar_has_the_circles_of_the_fifth_dimension(self):
        points = self.cigar["points"]
        self.assertEqual(points[0][:2], [0, 0])
        self.assertEqual(points[-1][0], 8.0)
        for rho, R, Z in points[1::7]:
            self.assertLess(abs(R - 4 * rho / math.sqrt(16 + rho * rho)), 2e-6)

    def test_the_cigar_rises_at_its_slope(self):
        def slope(u):
            r = 1 + u * u / 16
            dR = 2 * math.sqrt(1 - 0.75 / r) / r ** 2
            return 2 * r / math.sqrt(4 + u * u) * math.sqrt(max(1 - dR * dR, 0.0))
        for rho, R, Z in self.cigar["points"][1::9] + [self.cigar["points"][-1]]:
            self.assertLess(abs(Z - simpson(slope, 0.0, rho, 600)), 2e-5, f"at {rho}")

    def test_the_cigar_is_7_79_long(self):
        r = 5.0
        self.assertEqual(f"{math.sqrt((r - 0.75) * (r - 1)) + 1.75 * math.asinh(4.0):.2f}", "7.79")
        quad = simpson(lambda u: 2 * (1 + u * u / 16) / math.sqrt(4 + u * u), 0.0, 8.0)
        self.assertAlmostEqual(quad, math.sqrt(4.25 * 4) + 1.75 * math.asinh(4.0), places=8)

    def test_the_equator_stands_vertical_on_the_bubble_and_meets_its_mirror_there(self):
        near, far = self.equator["near"]["points"], self.equator["far"]["points"]
        self.assertEqual(near[0][0], 1.0)
        self.assertEqual(near[0][1:], far[0][1:])
        for (r, R, Z), (_, R2, Z2) in zip(near, far):
            self.assertLess(abs(R - r), 2e-6)
            self.assertAlmostEqual(Z, -Z2, places=7)
        # The first chord out of the bubble climbs far faster than it widens.
        (r0, R0, Z0), (r1, R1, Z1) = near[0], near[1]
        self.assertGreater(abs(Z1 - Z0), 3 * abs(R1 - R0))


@unittest.skipUnless(HAS_SYMPY, "sympy is not installed")
class FieldEquations(unittest.TestCase):
    def test_every_published_chart_is_the_star(self):
        import print_charts
        import verify_metrics as vm
        metric = json.loads((DATA / "metrics" / "topological_star.json").read_text(encoding="utf-8"))
        self.assertEqual([c["id"] for c in metric["coordinates"]], print_charts.TOPOLOGICAL_STAR_CHARTS)
        for entry in metric["coordinates"]:
            reader = vm.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]], ())
            symbols = [reader.symbol[c] for c in entry["coords"]]
            n = len(symbols)
            published = {tuple(e["indices"]): reader(e["value"]) for e in entry["metric_components"]}
            g = vm.sp.Matrix(n, n, lambda i, j: published.get((entry["coords"][i], entry["coords"][j]), 0))
            chart = types.SimpleNamespace(geo=vm.Geometry(g, symbols, 600), reader=reader, symbols=symbols,
                                          coords_tex=entry["coords"])
            # Raises where the published metric misses Einstein's equations with Maxwell's field, the
            # published black string or Schwarzschild's metric at the two ends of the family, the
            # exchange of t with y, Bah and Heidmann's chart under its map, the plane at the bubble, or
            # the equations of four dimensions.
            print_charts.topological_star_check(chart, entry["id"])


if __name__ == "__main__":
    unittest.main()
