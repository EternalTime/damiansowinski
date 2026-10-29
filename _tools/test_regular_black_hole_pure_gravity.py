"""Bueno, Cano and Hennigar's regular black hole of five dimensions:
.venv.noindex/bin/python -m unittest discover -s _tools

What its texts and drawings state, held on the published files: the metric function solves the
algebraic field equation of the quasi-topological tower summed as a geometric series, the
horizons and surface gravities of the hole drawn, the centre's de Sitter curvature, the
tortoise coordinate the drawings use, and the embedding surfaces against the closed form of
their height. Where sympy is installed, the published metric is held to the general-D function
at D = 5 and to Hayward's published metric at D = 4.
"""
import importlib.util
import json
import math
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "MFS" / "assets" / "data"
HAS_SYMPY = importlib.util.find_spec("sympy") is not None
ID = "regular_black_hole_pure_gravity"

# The hole as every drawing takes it, in units of r_0.
R0, ELL = 1.0, 12 / 25
R_PLUS, R_MINUS = 4 / 5, 3 / 5


def f(r, r0=R0, ell=ELL):
    return 1 - r0 ** 2 * r ** 2 / (r ** 4 + ell ** 2 * r0 ** 2)


def rstar(r):
    return r + 10 / 7 * math.log(abs((r - 0.8) / (r + 0.8))) - 15 / 14 * math.log(abs((r - 0.6) / (r + 0.6)))


def load(kind):
    return json.loads((DATA / kind / f"{ID}.json").read_text(encoding="utf-8"))


class FieldEquation(unittest.TestCase):
    def test_the_metric_function_solves_the_summed_tower(self):
        # h(psi) = psi/(1 - alpha psi) = m/r^4 at D = 5, with m = r_0^2 and alpha = l^2.
        for r0, ell in ((1.0, 0.48), (1.3, 0.2), (2.0, 1.1)):
            for r in (0.05, 0.5, 1.0, 3.0):
                psi = (1 - f(r, r0, ell)) / r ** 2
                # 1 - f loses digits to cancellation at small r, so the agreement is relative.
                self.assertAlmostEqual(psi / (1 - ell ** 2 * psi) * r ** 4 / r0 ** 2, 1.0, places=8)

    def test_every_finite_truncation_falls_short_of_the_series(self):
        # The partial sums of the geometric series approach h from below for 0 < alpha psi < 1.
        psi, a = 1.7, 0.48 ** 2
        partial = [sum(a ** (n - 1) * psi ** n for n in range(1, N + 1)) for N in (1, 2, 5, 40)]
        self.assertTrue(all(x < y for x, y in zip(partial, partial[1:])))
        self.assertAlmostEqual(partial[-1], psi / (1 - a * psi), places=10)


class Hole(unittest.TestCase):
    def test_the_horizons_and_surface_gravities_drawn(self):
        self.assertAlmostEqual(f(R_PLUS), 0.0, places=14)
        self.assertAlmostEqual(f(R_MINUS), 0.0, places=14)
        self.assertAlmostEqual(R_PLUS ** 2 * R_MINUS ** 2, ELL ** 2, places=14)
        self.assertAlmostEqual(R_PLUS ** 2 + R_MINUS ** 2, R0 ** 2, places=14)
        h = 1e-7
        for r, kappa in ((R_PLUS, 7 / 20), (R_MINUS, -7 / 15)):
            self.assertAlmostEqual((f(r + h) - f(r - h)) / (4 * h), kappa, places=7)
            self.assertAlmostEqual((r ** 2 - 2 * ELL ** 2) / r ** 3, kappa, places=12)

    def test_the_extremal_hole_is_at_r0_equal_2l(self):
        ell = 0.5
        x = 2 * ell ** 2
        self.assertAlmostEqual(x * x - x + ell ** 2, 0.0, places=14)
        self.assertAlmostEqual(2 * x - 1, 0.0, places=14)

    def test_the_tortoise_coordinate_is_the_integral_of_1_over_f(self):
        self.assertEqual(rstar(0.0), 0.0)
        for a, b in ((1.0, 2.0), (0.1, 0.5), (0.65, 0.75)):
            n = 20000
            h = (b - a) / n
            simpson = h / 3 * (1 / f(a) + 1 / f(b) + sum((4 if k % 2 else 2) / f(a + k * h) for k in range(1, n)))
            self.assertAlmostEqual(rstar(b) - rstar(a), simpson, places=8)

    def test_the_centre_is_de_sitter_in_five_dimensions(self):
        metric = json.loads((DATA / "metrics" / f"{ID}.json").read_text(encoding="utf-8"))
        chart = next(c for c in metric["coordinates"] if c["id"] == "spherical")
        self.assertIn("5\\ell^8\\,r_0^8", chart["kretschmann"])
        # K and R as published, at r = 0: 8 r_0^4 5 l^8 r_0^8/(l^2 r_0^2)^6 and 4 l^2 r_0^4 5 l^2 r_0^2/(l^2 r_0^2)^3.
        K0 = 8 * R0 ** 4 * 5 * ELL ** 8 * R0 ** 8 / (ELL ** 2 * R0 ** 2) ** 6
        R00 = 4 * ELL ** 2 * R0 ** 4 * 5 * ELL ** 2 * R0 ** 2 / (ELL ** 2 * R0 ** 2) ** 3
        self.assertAlmostEqual(K0 * ELL ** 4, 40.0, places=10)
        self.assertAlmostEqual(R00 * ELL ** 2, 20.0, places=10)


class Drawings(unittest.TestCase):
    def test_the_embedding_heights_are_the_closed_form(self):
        views = {v["id"]: v for v in load("embedding")["views"]}
        a2, b2 = R_PLUS ** 2, R_MINUS ** 2
        outside = views["outside"]["surfaces"][0]["pieces"][0]["points"]
        z0 = outside[0][2]
        for r, _, z in outside[5::40]:
            if r > R_PLUS + 0.05:
                self.assertAlmostEqual(abs(z - z0), 0.5 * math.acosh((2 * r * r - a2 - b2) / (a2 - b2)), delta=2e-3)
        inside = views["inside"]["surfaces"][0]["pieces"][0]["points"]
        top = inside[-1][2]
        for r, _, z in inside[5::40]:
            if 0.05 < r < R_MINUS - 0.05:
                self.assertAlmostEqual(abs(z - top), 0.5 * math.acosh((a2 + b2 - 2 * r * r) / (a2 - b2)), delta=2e-3)

    def test_the_curvature_text_fits_in_200_characters(self):
        for v in load("curvature")["views"]:
            text = list(v["caption"]) + [v["clip"], v["settings"]] + [e[2] for e in v["legend"]]
            self.assertLessEqual(sum(len(t) for t in text), 200, v["id"])


@unittest.skipUnless(HAS_SYMPY, "sympy is not installed")
class Published(unittest.TestCase):
    def test_the_chart_is_their_function_at_d5_and_haywards_at_d4(self):
        import sympy as sp
        D, m, a, x = sp.symbols("D m alpha x", positive=True)
        general = 1 - m * x ** 2 / (x ** (D - 1) + a * m)
        r0, ell = sp.symbols("r_0 ell", positive=True)
        five = 1 - r0 ** 2 * x ** 2 / (x ** 4 + ell ** 2 * r0 ** 2)
        self.assertEqual(sp.simplify(general.subs({D: 5, m: r0 ** 2, a: ell ** 2}) - five), 0)
        mh, l = sp.symbols("m_H l", positive=True)
        hayward = 1 - 2 * mh * x ** 2 / (x ** 3 + 2 * mh * l ** 2)
        self.assertEqual(sp.simplify(general.subs({D: 4, m: 2 * mh, a: l ** 2}) - hayward), 0)


if __name__ == "__main__":
    unittest.main()
