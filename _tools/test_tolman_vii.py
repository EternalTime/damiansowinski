"""Tolman's solution VII, held on the published file with the checker's own reader:
.venv.noindex/bin/python -m unittest discover -s _tools

The reader reads an arctangent and holds a listed name as a function. The star's chart is a
perfect fluid whose density is 15 beta (1 - r^2/R^2)/R^2, it meets Schwarzschild's exterior at
its surface, its central pressure is infinite at beta = 0.3862 and sound at its centre reaches
the speed of light at beta = 0.2698, Lattimer and Prakash's numbers, and Tolman's chart with his
constants for the same star is the same metric. Needs sympy, and is skipped where it is absent.
"""
import importlib.util
import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "MFS" / "assets" / "data"
HAS_SYMPY = importlib.util.find_spec("sympy") is not None
sys.path.insert(0, str(Path(__file__).resolve().parent / "derivations"))


@unittest.skipUnless(HAS_SYMPY, "sympy is not installed")
class TheReader(unittest.TestCase):
    def setUp(self):
        import sympy
        import verify_metrics
        self.sp, self.vm = sympy, verify_metrics

    def test_an_arctangent_is_read(self):
        reader = self.vm.Reader(["t", "x"], ["a"], ())
        x, a = reader.symbol["x"], reader.parameters["a"]
        self.assertEqual(reader("\\arctan\\left(\\sqrt{x/a}\\right)"), self.sp.atan(self.sp.sqrt(x / a)))
        self.assertEqual(reader("a\\tan x"), a * self.sp.tan(x))

    def test_a_held_name_is_a_function_whose_definition_surface_writes_out(self):
        reader = self.vm.Reader(["t", "x"], ["a", "f = \\ln\\left(a + x^2\\right)"], (), held=("f",))
        x, a, f = reader.symbol["x"], reader.parameters["a"], reader.parameters["f"]
        self.assertEqual(f.func.__name__, "f")
        self.assertEqual(f.args, (x,))
        rate = self.sp.Derivative(f, x)
        self.assertEqual(self.sp.simplify(reader.surface(rate) - 2 * x / (a + x ** 2)), 0)

    def test_a_held_name_has_to_be_defined(self):
        with self.assertRaises(self.vm.LatexError):
            self.vm.Reader(["t", "x"], ["a"], (), held=("a",))


@unittest.skipUnless(HAS_SYMPY, "sympy is not installed")
class TheStar(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import sympy
        import verify_metrics
        cls.sp, cls.vm = sympy, verify_metrics
        metric = json.loads((DATA / "metrics" / "tolman_vii.json").read_text(encoding="utf-8"))
        cls.charts = {}
        for entry in metric["coordinates"]:
            reader = cls.vm.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]], ())
            cls.charts[entry["id"]] = (reader, cls.vm.metric_from_line_element(reader, entry["line_element"],
                                                                               entry["coords"]))

    def star(self, compactness):
        """g_tt/c^2 and g_rr of the star's chart as functions of r, at R = 1."""
        reader, g = self.charts["spherical"]
        at = {reader.parameters["R"]: 1, reader.parameters["r_s"]: 2 * self.sp.Rational(compactness), reader.c: 1}
        return reader.symbol["r"], g[0, 0].subs(at), g[1, 1].subs(at)

    def fluid(self, compactness, x):
        """(8 pi G rho/c^2, 8 pi G p/c^4, the same from G^theta_theta) at r = x, R = 1, from the
        static sphere's Einstein tensor in the metric functions: G^t_t = -(r(1 - Z))'/r^2,
        G^r_r = Z(nu'/r + 1/r^2) - 1/r^2 and G^theta_theta = Z(nu''/2 + nu'^2/4 + nu'/2r) + Z'(nu'/4 + 1/2r)."""
        sp = self.sp
        r, gtt, grr = self.star(compactness)
        Z, nu = 1 / grr, sp.log(-gtt)
        d1, d2 = sp.diff(nu, r), sp.diff(nu, r, 2)
        density = sp.diff(r * (1 - Z), r) / r ** 2
        radial = Z * (d1 / r + 1 / r ** 2) - 1 / r ** 2
        across = Z * (d2 / 2 + d1 ** 2 / 4 + d1 / (2 * r)) + sp.diff(Z, r) * (d1 / 4 + 1 / (2 * r))
        return [sp.N(e.subs(r, x), 30) for e in (density, radial, across)]

    def test_the_fluid_is_perfect_and_its_density_falls_as_the_square_of_the_radius(self):
        for compactness, x in (("1/4", "1/3"), ("1/10", "4/5"), ("3/10", "1/2")):
            beta, x = self.sp.Rational(compactness), self.sp.Rational(x)
            density, radial, across = self.fluid(compactness, x)
            self.assertAlmostEqual(float(radial - across), 0, places=20)
            self.assertAlmostEqual(float(density - 15 * beta * (1 - x ** 2)), 0, places=20)
            self.assertGreater(radial, 0)

    def test_the_surface_is_schwarzschilds(self):
        r, gtt, grr = self.star("1/4")
        self.assertAlmostEqual(float(self.sp.N(gtt.subs(r, 1), 30)), -0.5, places=20)
        self.assertAlmostEqual(float(self.sp.N(grr.subs(r, 1), 30)), 2.0, places=20)
        self.assertAlmostEqual(float(self.fluid("1/4", 1)[1]), 0, places=20)

    def central(self, beta):
        """tan(psi) at the centre, from g_tt = -(1 - 5 beta/3) cos^2(psi)."""
        sp = self.sp
        r, gtt, _ = self.star(beta)
        cos2 = -gtt.subs(r, 0) / (1 - 5 * sp.Rational(beta) / 3)
        return sp.sqrt(1 / cos2 - 1)

    def test_the_central_pressure_is_infinite_at_lattimer_and_prakashs_compactness(self):
        # psi = pi/2 at the centre, where g_tt vanishes: between beta = 0.3861 and 0.3863.
        sp = self.sp
        before, after = (float(sp.N(self.central(b), 30)) for b in ("3861/10000", "38619/100000"))
        self.assertGreater(after, 10 * before)
        self.assertGreater(after, 1e4)

    def test_sound_at_the_centre_reaches_light_at_lattimer_and_prakashs_compactness(self):
        # Their (18): c_s^2/c^2 = tan(psi_c)(tan(psi_c)/5 + sqrt(beta/3)).
        sp = self.sp

        def speed(beta):
            tangent = self.central(beta)
            return float(sp.N(tangent * (tangent / 5 + sp.sqrt(sp.Rational(beta) / 3)), 30))
        self.assertLess(speed("2697/10000"), 1)
        self.assertGreater(speed("2699/10000"), 1)

    def test_tolmans_constants_for_the_same_star_give_the_same_metric(self):
        sp = self.sp
        reader, g = self.charts["tolman"]
        p = reader.parameters
        C = sp.sqrt(3) / 2 * (sp.Rational(1, 6) + sp.sqrt(6) / 3) * sp.exp(2 * sp.atan(1 / sp.sqrt(6)) - sp.pi)
        at = {p["R"]: 4 / sp.sqrt(5), p["A"]: 4 / sp.Rational(3) ** sp.Rational(1, 4), p["B"]: sp.sqrt(sp.Rational(7, 12)),
              p["C"]: C, reader.c: 1}
        own, radius = self.charts["spherical"]
        star = {own.parameters["R"]: 2, own.parameters["r_s"]: 1, own.c: 1}
        for x in ("0", "1/2", "9/7", "2"):
            for k in (0, 1):
                mine = sp.N(radius[k, k].subs(star).subs(own.symbol["r"], sp.Rational(x)), 30)
                his = sp.N(g[k, k].subs(at).subs(reader.symbol["r"], sp.Rational(x)), 30)
                self.assertAlmostEqual(float(mine - his), 0, places=20, msg=f"g_{k}{k} at r = {x}")


if __name__ == "__main__":
    unittest.main()
