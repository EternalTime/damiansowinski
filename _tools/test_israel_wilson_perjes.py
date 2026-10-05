"""Israel, Wilson and Perjes's stationary fields, Majumdar and Papapetrou's charges set spinning:
.venv.noindex/bin/python -m unittest discover -s _tools

What its texts state, held on the published files and on the two sources its drawings declare.
The published spheroidal chart has no horizon and a g_tt that vanishes only on the ring r = 0 of
its equator; its circles of phi turn null on the equator where the captions say; the published
spherical chart has its degenerate horizon at r = m, of area 4 pi (m^2 + l^2). The declared two
sources solve Laplace's equation and the equation for omega, omega vanishes on the axis beyond
them and is 4l(1 + m/d) between them, Hartle and Hawking's (4.24), and far away it is the field
of the angular momentum 2ld, their Na. The tests need sympy and are skipped where it is absent.
"""
import importlib.util
import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "MFS" / "assets" / "data"
HAS_SYMPY = importlib.util.find_spec("sympy") is not None and importlib.util.find_spec("scipy") is not None
sys.path.insert(0, str(Path(__file__).resolve().parent / "derivations"))


@unittest.skipUnless(HAS_SYMPY, "sympy and scipy are not installed")
class Published(unittest.TestCase):
    def chart(self, system, values):
        """The published metric of a chart as a sympy matrix, G = c = 1, at the parameter values given."""
        import sympy as sp
        import verify_metrics as vm
        metric = json.loads((DATA / "metrics" / "israel_wilson_perjes.json").read_text(encoding="utf-8"))
        entry = next(c for c in metric["coordinates"] if c["id"] == system)
        reader = vm.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]], ())
        there = {tuple(e["indices"]): e["value"] for e in entry["metric_components"]}
        g = sp.Matrix(4, 4, lambda i, j: reader(there.get((entry["coords"][i], entry["coords"][j]), "0")))
        at = {reader.parameters[k]: v for k, v in values.items()}
        at[reader.c] = 1
        return g.subs(reader.held).subs(at) if getattr(reader, "held", None) else g.subs(at), reader

    def test_the_spinning_source_has_no_horizon_and_no_ergoregion(self):
        import sympy as sp
        g, reader = self.chart("spheroidal", {"m": 1, "a": 1})
        r, theta = reader.symbol["r"], reader.symbol["\\theta"]
        grr = sp.simplify(g[1, 1])
        gtt = sp.simplify(g[0, 0])
        for x in (-3, -1, sp.Rational(-1, 2), 0, 1, 5):
            for angle in (0, sp.pi / 5, sp.pi / 3):
                self.assertGreater(1 / grr.subs({r: x, theta: angle}), 0)
                self.assertLess(gtt.subs({r: x, theta: angle}), 0)
        # On the equator g_tt = -r^2/(r + m)^2, which vanishes on the ring r = 0 and is negative beside it.
        self.assertEqual(sp.simplify(gtt.subs(theta, sp.pi / 2) + r ** 2 / (r + 1) ** 2), 0)

    def test_the_circles_of_the_equator_turn_null_where_the_captions_say(self):
        import sympy as sp
        g, reader = self.chart("spheroidal", {"m": 1, "a": 1})
        r, theta = reader.symbol["r"], reader.symbol["\\theta"]
        gpp = sp.simplify(g[3, 3].subs(theta, sp.pi / 2))
        R = sp.Symbol("R")
        root = sp.nsolve(R ** 4 + R ** 2 + 2 * R - 1, R, 0.4) - 1
        self.assertAlmostEqual(float(root), -0.595, places=3)
        self.assertAlmostEqual(float(gpp.subs(r, root)), 0.0, places=10)
        self.assertLess(float(gpp.subs(r, root - 0.01)), 0)
        self.assertGreater(float(gpp.subs(r, root + 0.01)), 0)

    def test_the_source_of_complex_mass_has_a_degenerate_horizon_of_area_four_pi_m2_plus_l2(self):
        import sympy as sp
        g, reader = self.chart("spherical", {})
        r, theta = reader.symbol["r"], reader.symbol["\\theta"]
        m, l = reader.parameters["m"], reader.parameters["l"]
        inverse_grr = sp.simplify(1 / g[1, 1])
        self.assertEqual(sp.simplify(inverse_grr.subs(r, m)), 0)
        self.assertEqual(sp.simplify(sp.diff(inverse_grr, r).subs(r, m)), 0)
        # At the horizon g_tt = g_tphi = 0, so the sphere's area is the integral of sqrt(g_thth g_phph).
        element = sp.sqrt(sp.simplify((g[2, 2] * g[3, 3]).subs(r, m) / sp.sin(theta) ** 2))
        self.assertEqual(sp.simplify(4 * sp.pi * element - 4 * sp.pi * (m ** 2 + l ** 2)), 0)


@unittest.skipUnless(HAS_SYMPY, "sympy and scipy are not installed")
class TwoSources(unittest.TestCase):
    """The W and omega the drawings declare, m = 1, l = 1/2 and d = 2, as null_rays.IWP_TWO writes them."""

    def setUp(self):
        import sympy as sp
        import null_rays as nr
        self.sp = sp
        self.rho, self.z = sp.symbols("rho z", positive=True)
        names = {"rho": self.rho, "z": self.z}
        self.W = sp.sympify(nr.IWP_TWO["W"], locals=names)
        self.omega = sp.sympify(nr.IWP_TWO["omega"], locals=names)
        r1, r2 = sp.sqrt(self.rho ** 2 + (self.z - 2) ** 2), sp.sqrt(self.rho ** 2 + (self.z + 2) ** 2)
        self.P, self.Q = 1 + 1 / r1 + 1 / r2, sp.Rational(1, 2) / r2 - sp.Rational(1, 2) / r1

    def at(self, expression, rho, z):
        return float(expression.subs({self.rho: self.sp.Rational(rho), self.z: self.sp.Rational(z)}))

    def test_w_is_the_modulus_of_a_harmonic_potential(self):
        sp, rho, z = self.sp, self.rho, self.z
        for where in (("7/5", "3/4"), ("1/3", "5/2"), ("4", "1/10")):
            self.assertAlmostEqual(self.at(self.W ** 2 - self.P ** 2 - self.Q ** 2, *where), 0.0, places=12)
            for f in (self.P, self.Q):
                laplacian = sp.diff(f, rho, 2) + sp.diff(f, rho) / rho + sp.diff(f, z, 2)
                self.assertAlmostEqual(self.at(laplacian, *where), 0.0, places=12)

    def test_omega_solves_its_equation(self):
        sp, rho, z, P, Q = self.sp, self.rho, self.z, self.P, self.Q
        along_z = sp.diff(self.omega, z) + 2 * rho * (P * sp.diff(Q, rho) - Q * sp.diff(P, rho))
        along_rho = sp.diff(self.omega, rho) - 2 * rho * (P * sp.diff(Q, z) - Q * sp.diff(P, z))
        for where in (("7/5", "3/4"), ("1/3", "5/2"), ("4", "1/10")):
            self.assertAlmostEqual(self.at(along_z, *where), 0.0, places=12)
            self.assertAlmostEqual(self.at(along_rho, *where), 0.0, places=12)

    def test_omega_on_the_axis_and_far_away(self):
        # Beyond the sources omega vanishes on the axis; between them it is 4l(1 + m/d) = 3.
        self.assertAlmostEqual(self.at(self.omega, "1/1000000", 3), 0.0, places=9)
        self.assertAlmostEqual(self.at(self.omega, "1/1000000", -5), 0.0, places=9)
        self.assertAlmostEqual(self.at(self.omega, "1/1000000", 0), 3.0, places=9)
        # Far out on the midplane omega rho tends to twice the angular momentum 2ld = 2.
        self.assertAlmostEqual(self.at(self.omega * self.rho, 100000, 0), 4.0, places=3)

    def test_the_midplane_as_the_caption_writes_it(self):
        sp, rho = self.sp, self.rho
        s = sp.sqrt(rho ** 2 + 4)
        self.assertEqual(sp.simplify(self.W.subs(self.z, 0) - (1 + 2 / s)), 0)
        self.assertEqual(sp.simplify(self.omega.subs(self.z, 0) - (4 / s + 4 / s ** 2)), 0)
        null = sp.nsolve((1 + 2 / s) ** 4 * rho ** 2 - (4 / s + 4 / s ** 2) ** 2, rho, 0.7)
        self.assertAlmostEqual(float(null), 0.734, places=3)


if __name__ == "__main__":
    unittest.main()
