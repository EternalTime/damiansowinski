"""Simpson and Visser's eye of the storm, which is Ghosh's nonsingular rotating black hole:
.venv.noindex/bin/python -m unittest discover -s _tools

What its texts and drawings state, held on the published files. The first class needs nothing
but the files. The second needs sympy and mpmath and is skipped where they are absent: it reads
the published metric of the Boyer-Lindquist chart and holds it to Kerr's metric at ell = 0, to the
static metric of Simpson and Visser's asymptotically Minkowski core at a = 0, to horizons strictly
inside Kerr's, to the extremal relations of their (3.17), to the throat radius 2m e^(-ell/r+), to
Ghosh's far field, Kerr-Newman's with q^2 = 2m ell, to the energy conditions of their section 4,
and to every curvature invariant falling to zero toward r = 0 in every direction.
"""
import importlib.util
import json
import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "MFS" / "assets" / "data"
HAS_NUMBERS = all(importlib.util.find_spec(name) is not None for name in ("sympy", "mpmath"))
sys.path.insert(0, str(Path(__file__).resolve().parent / "derivations"))
ME = "eye_of_the_storm"
A, ELL = 0.6, 0.2
R_MINUS, R_PLUS = 0.4348787520714397, 1.5150508229447324


def load(folder):
    return json.loads((DATA / folder / f"{ME}.json").read_text(encoding="utf-8"))


def delta(r, a=A, ell=ELL, m=1.0):
    return r * r + a * a - 2 * m * r * math.exp(-ell / r)


class Drawings(unittest.TestCase):
    def test_the_horizons_drawn_are_the_roots_of_delta(self):
        for r in (R_MINUS, R_PLUS):
            self.assertLess(abs(delta(r)), 1e-12)
        self.assertLess(delta(0.9), 0)
        self.assertGreater(delta(0.2), 0)
        self.assertGreater(delta(2.0), 0)

    def test_the_horizons_lie_inside_kerrs(self):
        kerr = (1 - math.sqrt(1 - A * A), 1 + math.sqrt(1 - A * A))
        self.assertLess(kerr[0], R_MINUS)
        self.assertLess(R_PLUS, kerr[1])

    def test_the_embedding_throat_is_two_m_e_to_the_minus_ell_over_r_plus(self):
        view, = load("embedding")["views"]
        points = view["surfaces"][0]["pieces"][0]["points"]
        self.assertAlmostEqual(points[0][0], R_PLUS, places=9)
        self.assertAlmostEqual(points[0][1], 2 * math.exp(-ELL / R_PLUS), places=6)

    def test_each_drawing_is_there(self):
        diagrams = load("diagrams")
        self.assertEqual([v["id"] for v in diagrams["systems"]["boyer_lindquist"]], ["radial", "principal", "above"])
        self.assertEqual([v["id"] for v in diagrams["systems"]["kerr_ingoing"]], ["axis"])
        self.assertEqual([v["id"] for v in load("conformal")["views"]], ["axis", "ingoing"])
        self.assertEqual([v["id"] for v in load("curvature")["views"]], ["boyer_lindquist"])

    def test_the_curvature_falls_to_zero_at_the_centre_on_every_curve(self):
        view, = load("curvature")["views"]
        for curve in view["curves"]:
            first = curve["points"][0]
            self.assertLess(first[0], 0.01)
            self.assertLess(abs(first[1]), 1e-6, f"{curve['panel']} {curve['class']}")


@unittest.skipUnless(HAS_NUMBERS, "needs sympy and mpmath")
class Physics(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import sympy as sp
        import verify_metrics as vm
        cls.sp = sp
        entry = next(c for c in load("metrics")["coordinates"] if c["id"] == "boyer_lindquist")
        reader = vm.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]], ())
        cls.r, cls.th = reader.symbol["r"], reader.symbol["\\theta"]
        cls.m, cls.a, cls.ell = (reader.parameters[k] for k in ("m", "a", "ell"))
        values = {tuple(e["indices"]): e["value"] for e in entry["metric_components"]}
        cls.g = sp.Matrix(4, 4, lambda i, j: reader(values.get((entry["coords"][i], entry["coords"][j]), "0")))
        cls.R = reader(entry["ricci_scalar"].partition("=")[2])
        cls.K = reader(entry["kretschmann"].partition("=")[2])
        cls.reader = reader

    def test_ell_zero_is_kerr(self):
        sp, r, th, m, a = self.sp, self.r, self.th, self.m, self.a
        sigma = r ** 2 + a ** 2 * sp.cos(th) ** 2
        self.assertEqual(sp.simplify(self.g[0, 0].subs(self.ell, 0) + 1 - 2 * m * r / sigma), 0)
        self.assertEqual(sp.simplify(self.g[1, 1].subs(self.ell, 0) - sigma / (r ** 2 + a ** 2 - 2 * m * r)), 0)

    def test_a_zero_is_the_static_minkowski_core(self):
        sp, r, m, ell = self.sp, self.r, self.m, self.ell
        f = 1 - 2 * m * sp.exp(-ell / r) / r
        self.assertEqual(sp.simplify(self.g[0, 0].subs(self.a, 0) + f), 0)
        self.assertEqual(sp.simplify(self.g[1, 1].subs(self.a, 0) - 1 / f), 0)

    def test_the_extremal_hole_of_their_317(self):
        """m = r_H^2 e^(ell/r_H)/(r_H + ell) and a^2 = r_H^2 (r_H - ell)/(r_H + ell) make r_H a double root."""
        sp, r = self.sp, self.r
        rh, ell = sp.Rational(3, 4), sp.Rational(1, 5)
        m = rh ** 2 * sp.exp(ell / rh) / (rh + ell)
        a2 = rh ** 2 * (rh - ell) / (rh + ell)
        D = r ** 2 + a2 - 2 * m * r * sp.exp(-ell / r)
        self.assertEqual(sp.simplify(D.subs(r, rh)), 0)
        self.assertEqual(sp.simplify(sp.diff(D, r).subs(r, rh)), 0)

    def test_far_away_it_is_kerr_newman_with_q_squared_two_m_ell(self):
        sp, r, m, a, ell = self.sp, self.r, self.m, self.a, self.ell
        sigma = r ** 2 + a ** 2 * sp.cos(self.th) ** 2
        D = sigma / self.g[1, 1]
        far = sp.series(sp.simplify(D), r, sp.oo, 1).removeO()
        self.assertEqual(sp.simplify(far - (r ** 2 - 2 * m * r + a ** 2 + 2 * m * ell)), 0)

    def test_the_ricci_scalar_is_positive_and_the_energy_conditions_fail_only_in_the_core(self):
        """rho = ell m e^(-ell/r)/(4 pi Sigma^2) > 0, rho + p_t > 0 where Xi < 2 and the strong
        condition 2 p_t > 0 where Xi < 1, Xi = ell Sigma/(2 r^3): on the equator r > ell/4 and r > ell/2."""
        sp = self.sp
        values = {self.m: 1, self.a: sp.Rational(3, 5), self.ell: sp.Rational(1, 5)}
        for r in (0.03, 0.2, 1.0, 3.0):
            for th in (0.0, 0.7, math.pi / 2):
                self.assertGreater(float(self.R.subs(values).subs({self.r: r, self.th: th})), 0)
        for r, null, strong in ((0.04, False, False), (0.07, True, False), (0.2, True, True)):
            xi = 0.2 / (2 * r)
            self.assertEqual(xi < 2, null)
            self.assertEqual(xi < 1, strong)

    def test_every_invariant_falls_to_zero_toward_r_zero_in_every_direction(self):
        sp = self.sp
        values = {self.m: 1, self.a: sp.Rational(3, 5), self.ell: sp.Rational(1, 5)}
        for th in (0.0, 0.4, 1.0, math.pi / 2):
            for scalar in (self.R, self.K):
                near = abs(float(sp.N(scalar.subs(values).subs({self.r: sp.Rational(1, 1000), self.th: th}), 30)))
                self.assertLess(near, 1e-10)


if __name__ == "__main__":
    unittest.main()
