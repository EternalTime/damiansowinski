"""Israel's disc source of the Kerr metric:
.venv.noindex/bin/python -m unittest discover -s _tools

What the page and its drawings state of the physics, held on the published files. The first class
needs nothing but the files. The second needs sympy and is skipped where it is absent: it reads the
published Boyer-Lindquist metric, takes the jump in its extrinsic curvature across the disc r = 0
whose two faces are glued, and holds the layer that jump carries to what the history and
kerr_disc_source.md state: its form lambda (w w + e e), the superluminal stream at speed a/rho,
the rigid observer who sees no energy, the surface mass density -ma/2 pi (a^2 - rho^2)^(3/2) and
angular momentum, and the layer vanishing with m.
"""
import importlib.util
import json
import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "MFS" / "assets" / "data"
HAS_NUMBERS = importlib.util.find_spec("sympy") is not None
sys.path.insert(0, str(Path(__file__).resolve().parent / "derivations"))


def load(folder):
    return json.loads((DATA / folder / "kerr_disc_source.json").read_text(encoding="utf-8"))


class Files(unittest.TestCase):
    def test_both_charts_keep_r_at_least_zero_and_glue_the_disc(self):
        charts = {c["id"]: c for c in load("metrics")["coordinates"]}
        self.assertEqual(set(charts), {"boyer_lindquist", "kerr_schild"})
        for chart in charts.values():
            self.assertIn("r \\in [0, \\infty)", chart["domains"])
            self.assertTrue(any("\\sim (0, \\pi - \\theta)" in d for d in chart["domains"]))
            names = [p["symbol"].split(" = ")[0] for p in chart["parameters"]]
            self.assertEqual(names, ["m", "a", "\\Sigma", "\\Delta"])

    def test_every_drawing_takes_a_past_m(self):
        for text in ([load("embedding")["views"][0]["settings"], load("conformal")["views"][0]["settings"]]):
            self.assertIn("a = 1.5\\,m", text)
        self.assertEqual(load("curvature")["views"][0]["map"]["params"], {"m": 1.0, "a": 1.5})


@unittest.skipUnless(HAS_NUMBERS, "needs sympy")
class Layer(unittest.TestCase):
    """The Lanczos tensor of the glued disc, from the published metric at points of the disc."""

    @classmethod
    def setUpClass(cls):
        import sympy as sp
        import verify_metrics as vm
        cls.sp = sp
        chart = next(c for c in load("metrics")["coordinates"] if c["id"] == "boyer_lindquist")
        names = [p["symbol"] for p in chart["parameters"]]
        reader = vm.Reader(chart["coords"], names, ())
        cls.r, cls.theta = reader.symbol["r"], reader.symbol["\\theta"]
        cls.m, cls.a = reader.parameters["m"], reader.parameters["a"]
        index = {name: k for k, name in enumerate(chart["coords"])}
        g = sp.zeros(4)
        for entry in chart["metric_components"]:
            i, j = (index[x] for x in entry["indices"])
            g[i, j] = reader(entry["value"])
        cls.g = g

    def layer(self, m, a, theta):
        """S^i_j in (t, theta, phi) on the upper face at r = 0, with G = c = 1:
        S_ij = -(1/4 pi)(K_ij - h_ij K), K_ij = d_r g_ij/(2 sqrt(g_rr)), each face's normal into r > 0."""
        sp = self.sp
        at = {self.m: m, self.a: a, self.theta: theta}
        idx = [0, 2, 3]
        h = sp.Matrix(3, 3, lambda i, j: self.g[idx[i], idx[j]].subs(self.r, 0).subs(at))
        grr = self.g[1, 1].subs(self.r, 0).subs(at)
        K = sp.Matrix(3, 3, lambda i, j: sp.diff(self.g[idx[i], idx[j]], self.r).subs(self.r, 0).subs(at)
                      / (2 * sp.sqrt(grr)))
        mixed = h.inv() * K
        S = -(mixed - sp.eye(3) * mixed.trace()) / (4 * sp.pi)
        return h, S

    def points(self):
        sp = self.sp
        for m, a, theta in ((1, sp.Rational(3, 2), sp.Rational(1, 3)), (sp.Rational(2, 3), 2, sp.Rational(9, 10)),
                            (3, sp.Rational(5, 4), sp.Rational(13, 10))):
            yield m, a, theta

    def test_the_layer_is_lambda_w_w_plus_e_e(self):
        sp = self.sp
        for m, a, theta in self.points():
            h, S = self.layer(m, a, theta)
            rho = a * sp.sin(theta)
            lam = -m / (4 * sp.pi * a * sp.sqrt(a ** 2 - rho ** 2))
            W = sp.Matrix([1, 0, a / rho ** 2])
            W = W / sp.sqrt((W.T * h * W)[0])
            E = sp.Matrix([0, 1 / sp.sqrt(h[1, 1]), 0])
            expected = lam * (W * (h * W).T + E * (h * E).T)
            self.assertLess(max(abs(sp.N(v, 40)) for v in S - expected), 1e-30)

    def test_the_stream_outruns_light_at_a_over_rho(self):
        sp = self.sp
        for m, a, theta in self.points():
            h, S = self.layer(m, a, theta)
            rho = a * sp.sin(theta)
            # The stream turns at Omega = a/rho^2, so its speed past the static frame is rho Omega.
            self.assertGreater(sp.N(rho * a / rho ** 2), 1)
            W = sp.Matrix([1, 0, a / rho ** 2])
            self.assertGreater(sp.N((W.T * h * W)[0]), 0)

    def test_a_rigid_observer_at_one_over_a_sees_no_energy(self):
        sp = self.sp
        for m, a, theta in self.points():
            h, S = self.layer(m, a, theta)
            U = sp.Matrix([1, 0, 1 / sp.sympify(a)])
            self.assertLess(sp.N((U.T * h * U)[0]), 0)
            self.assertLess(abs(sp.N((U.T * h * S * U)[0], 40)), 1e-30)

    def test_the_mass_and_angular_momentum_on_the_disc(self):
        sp = self.sp
        for m, a, theta in self.points():
            h, S = self.layer(m, a, theta)
            root = sp.sqrt(-h.det())
            rho = a * sp.sin(theta)
            drho = a * sp.cos(theta)
            sigma = -m * a / (2 * sp.pi * (a ** 2 - rho ** 2) ** sp.Rational(3, 2))
            tolman = (-S[0, 0] + S[1, 1] + S[2, 2]) * root * 2 * sp.pi
            self.assertLess(abs(sp.N(tolman - 2 * sp.pi * rho * drho * sigma, 40)), 1e-30)
            spin = S[0, 2] * root * 2 * sp.pi
            self.assertLess(abs(sp.N(spin + m * rho ** 3 / (2 * (a ** 2 - rho ** 2) ** sp.Rational(3, 2)) * drho, 40)),
                            1e-30)

    def test_the_layer_vanishes_with_m(self):
        sp = self.sp
        h, S = self.layer(0, sp.Rational(3, 2), sp.Rational(1, 2))
        self.assertTrue(all(sp.simplify(v) == 0 for v in S))

    def test_the_mass_density_diverges_at_the_rim(self):
        sp = self.sp
        values = []
        for theta in (sp.Rational(1, 2), sp.Rational(3, 2), sp.Rational(157, 100)):
            h, S = self.layer(1, sp.Rational(3, 2), theta)
            values.append(sp.N((-S[0, 0] + S[1, 1] + S[2, 2]) * sp.sqrt(-h.det())))
        self.assertTrue(values[0] > values[1] > values[2])
        self.assertLess(values[2], -1000)


if __name__ == "__main__":
    unittest.main()
