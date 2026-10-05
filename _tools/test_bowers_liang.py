"""Bowers and Liang's anisotropic star, held on the published file with the checker's own reader:
.venv.noindex/bin/python -m unittest discover -s _tools

The density is uniform, the pressure along the radius is their (3.4) and the pressure across it
exceeds that by their (3.1) with (3.2). The surface is Schwarzschild's, Q = 1/2 is Schwarzschild's
published interior and Q -> 0 Florides's published cluster. The central pressure is infinite at
their (3.6), the heaviest star holds 19 percent more than Schwarzschild's, and a redshift of
2.358 takes Q = 0.45. The star the diagrams draw, Q = 1/4 at R = 9 r_s/8, has a central pressure
of rho c^2/sqrt(3) and a pressure across the radius equal to its energy density at the surface.
Needs sympy, and is skipped where it is absent.
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


def published(vm, metric_id, system_id):
    """(reader, g) of a published chart in x^0 = ct, from its metric components as they are
    written, every defined name written out. The components are read as they stand, with no
    canonical form taken: the canonical form of a root or of a power Q of 1 - r_s r^2/R^3 is
    sound in a comparison and is not the number the power has."""
    import sympy as sp
    metric = json.loads((DATA / "metrics" / f"{metric_id}.json").read_text(encoding="utf-8"))
    entry = next(c for c in metric["coordinates"] if c["id"] == system_id)
    reader = vm.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]], (),
                       held=vm.HELD.get((metric_id, system_id), ()))
    g = sp.zeros(len(entry["coords"]), len(entry["coords"]))
    for component in entry["metric_components"]:
        i, j = (entry["coords"].index(k) for k in component["indices"])
        g[i, j] = reader(component["value"])
    return reader, g.subs(reader.held).doit()


@unittest.skipUnless(HAS_SYMPY, "sympy is not installed")
class TheStar(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import sympy
        import verify_metrics
        cls.sp, cls.vm = sympy, verify_metrics
        cls.reader, cls.g = published(verify_metrics, "bowers_liang", "areal")

    def star(self, exponent, compactness):
        """r, -g_tt and g_rr in x^0 = ct as functions of r, at R = 1, for Q and r_s/R given as fractions."""
        sp, p = self.sp, self.reader.parameters
        at = {p["R"]: 1, p["r_s"]: sp.Rational(compactness), p["Q"]: sp.Rational(exponent)}
        return self.reader.symbol["r"], -self.g[0, 0].subs(at), self.g[1, 1].subs(at)

    def stresses(self, exponent, compactness, x):
        """(8 pi G rho/c^2, 8 pi G p_r/c^4, 8 pi G p_perp/c^4) at r = x, R = 1, from the static
        sphere's Einstein tensor in the metric functions: G^t_t = -(r(1 - Z))'/r^2, G^r_r =
        Z(nu'/r + 1/r^2) - 1/r^2 and G^theta_theta = Z(nu''/2 + nu'^2/4 + nu'/2r) + Z'(nu'/4 + 1/2r)."""
        sp = self.sp
        r, lapse, grr = self.star(exponent, compactness)
        Z, nu = 1 / grr, sp.log(lapse)
        d1, d2 = sp.diff(nu, r), sp.diff(nu, r, 2)
        density = sp.diff(r * (1 - Z), r) / r ** 2
        radial = Z * (d1 / r + 1 / r ** 2) - 1 / r ** 2
        across = Z * (d2 / 2 + d1 ** 2 / 4 + d1 / (2 * r)) + sp.diff(Z, r) * (d1 / 4 + 1 / (2 * r))
        return [sp.N(e.subs(r, sp.Rational(x)), 30) for e in (density, radial, across)]

    def test_the_density_is_uniform_and_the_radial_pressure_is_bowers_and_liangs(self):
        sp = self.sp
        for exponent, compactness, x in (("1/4", "8/9", "1/3"), ("1/2", "2/3", "4/5"), ("3/4", "1/2", "1/2"),
                                         ("1/10", "19/20", "9/10")):
            Q, k, radius = sp.Rational(exponent), sp.Rational(compactness), sp.Rational(x)
            density, radial, _ = self.stresses(exponent, compactness, x)
            self.assertAlmostEqual(float(density - 3 * k), 0, places=20)
            y, y1 = (1 - k * radius ** 2) ** Q, (1 - k) ** Q
            # Their (3.4): p = rho_0 (y - y1)/(3 y1 - y).
            self.assertAlmostEqual(float(radial - 3 * k * (y - y1) / (3 * y1 - y)), 0, places=20)

    def test_the_pressure_across_the_radius_is_their_anisotropy(self):
        # Their (3.1) and (3.2) with n = 2 and C = (2 pi/3)(1 - 2Q), in 8 pi G/c^4 times the stresses.
        sp = self.sp
        for exponent, compactness, x in (("1/4", "8/9", "1/3"), ("3/4", "1/2", "1/2"), ("1/10", "19/20", "9/10")):
            Q, k, radius = sp.Rational(exponent), sp.Rational(compactness), sp.Rational(x)
            density, radial, across = self.stresses(exponent, compactness, x)
            wanted = (1 - 2 * Q) * radius ** 2 * (density + radial) * (density + 3 * radial) / (12 * (1 - k * radius ** 2))
            self.assertAlmostEqual(float(across - radial - wanted), 0, places=20)
            self.assertEqual(across > radial, Q < sp.Rational(1, 2))

    def test_equal_pressures_at_one_half(self):
        _, radial, across = self.stresses("1/2", "2/3", "4/5")
        self.assertAlmostEqual(float(across - radial), 0, places=20)

    def test_the_surface_is_schwarzschilds(self):
        sp = self.sp
        for exponent, compactness in (("1/4", "8/9"), ("3/4", "1/2")):
            r, lapse, grr = self.star(exponent, compactness)
            k = sp.Rational(compactness)
            self.assertAlmostEqual(float(sp.N(lapse.subs(r, 1), 30) - (1 - k)), 0, places=20)
            self.assertAlmostEqual(float(sp.N(grr.subs(r, 1), 30) - 1 / (1 - k)), 0, places=20)
            self.assertAlmostEqual(float(self.stresses(exponent, compactness, 1)[1]), 0, places=20)

    def test_one_half_is_schwarzschilds_published_interior(self):
        sp = self.sp
        reader, g = published(self.vm, "interior_schwarzschild", "spherical")
        r, lapse, _ = self.star("1/2", "2/3")
        his = {reader.parameters["R"]: 1, reader.parameters["r_s"]: sp.Rational(2, 3)}
        for x in ("0", "1/3", "9/10", "1"):
            mine = sp.N(lapse.subs(r, sp.Rational(x)), 30)
            theirs = sp.N(-g[0, 0].subs(his).subs(reader.symbol["r"], sp.Rational(x)), 30)
            self.assertAlmostEqual(float(mine - theirs), 0, places=20, msg=f"g_tt at r = {x}")

    def test_the_limit_of_no_radial_pressure_is_florides_published_cluster(self):
        sp = self.sp
        reader, g = published(self.vm, "einstein_cluster", "uniform")
        p = self.reader.parameters
        r, Q = self.reader.symbol["r"], p["Q"]
        k = sp.Rational(4, 5)
        lapse = -self.g[0, 0].subs({p["R"]: 1, p["r_s"]: k})
        his = {reader.parameters["R"]: 1, reader.parameters["r_s"]: k}
        for x in ("1/3", "9/10"):
            mine = sp.limit(lapse.subs(r, sp.Rational(x)), Q, 0, "+")
            theirs = -g[0, 0].subs(his).subs(reader.symbol["r"], sp.Rational(x))
            self.assertAlmostEqual(float(sp.N(mine - theirs, 30)), 0, places=20, msg=f"g_tt at r = {x}")
        # The pressure along the radius goes with Q.
        small = self.stresses("1/1000", "4/5", "1/2")
        self.assertLess(abs(float(small[1])), 1e-2)
        self.assertGreater(float(small[2]), 0.1)

    def test_the_central_pressure_is_infinite_at_their_critical_compactness(self):
        # Their (3.6): r_s/R = 1 - 3^(-1/Q), where -g_tt vanishes at the centre.
        sp = self.sp
        p = self.reader.parameters
        r = self.reader.symbol["r"]
        for exponent in ("1/2", "1/4", "1/3"):
            Q = sp.Rational(exponent)
            at = {p["R"]: 1, p["r_s"]: 1 - sp.Integer(3) ** (-1 / Q), p["Q"]: Q, r: 0}
            self.assertEqual(sp.simplify(self.g[0, 0].subs(at, simultaneous=True)), 0)
        self.assertEqual(1 - sp.Integer(3) ** (-2), sp.Rational(8, 9))
        self.assertEqual(1 - sp.Integer(3) ** (-4), sp.Rational(80, 81))

    def test_their_numbers(self):
        # (3.10): at one density the mass goes as (r_s/R)^(3/2), so the star that fills its own
        # Schwarzschild radius holds (9/8)^(3/2) = 1.19 of Schwarzschild's heaviest.
        self.assertAlmostEqual((9 / 8) ** 1.5, 1.19, places=2)
        # (3.12): z_crit = 3^(1/(2Q)) - 1, which is 2 at Q = 1/2 and 2.358 at Q = 0.453.
        self.assertAlmostEqual(3 ** (1 / (2 * 0.5)) - 1, 2.0, places=12)
        Q = math.log(3) / (2 * math.log(3.358))
        self.assertAlmostEqual(Q, 0.453, places=3)
        # Bondi's u < 0.485 is a surface redshift of 4.77.
        self.assertAlmostEqual((1 - 2 * 0.485) ** -0.5 - 1, 4.77, places=2)
        # Andreasson's bound 1 - 1/(1 + 2 Omega)^2 at Omega = 1 and 3.
        self.assertAlmostEqual(1 - 1 / 3 ** 2, 8 / 9, places=12)
        self.assertAlmostEqual(1 - 1 / 7 ** 2, 48 / 49, places=12)

    def test_the_star_the_diagrams_draw(self):
        # Q = 1/4 at Schwarzschild's limit R = 9 r_s/8.
        sp = self.sp
        # The Einstein tensor in the metric functions divides by r, so the centre is taken at r = 1e-9 R.
        density, radial, across = self.stresses("1/4", "8/9", "1/1000000000")
        self.assertAlmostEqual(float(radial / density), 1 / math.sqrt(3), places=9)
        self.assertAlmostEqual(float(across / density), 1 / math.sqrt(3), places=9)
        density, radial, across = self.stresses("1/4", "8/9", 1)
        self.assertAlmostEqual(float(radial), 0, places=20)
        self.assertAlmostEqual(float(across / density), 1, places=12)
        # The stresses add up to more than the energy density, as Andreasson's theorem asks of a
        # star at 8/9 whose stresses are not all at his bound.
        self.assertGreater(float(radial + 2 * across), float(density))
        r, lapse, grr = self.star("1/4", "8/9")
        self.assertAlmostEqual(float(sp.N(lapse.subs(r, 0), 20)), 0.0179, places=4)
        self.assertAlmostEqual(float(sp.N(lapse.subs(r, 1), 20)), 1 / 9, places=12)
        # A redshift of 2 at the surface, and a cap of the three sphere 70.5 degrees wide.
        self.assertAlmostEqual(float(sp.N(lapse.subs(r, 1), 20)) ** -0.5 - 1, 2.0, places=12)
        self.assertAlmostEqual(math.degrees(math.asin(math.sqrt(8 / 9))), 70.5, places=1)
        # Light takes 8.28 r_s/c of t from the centre to the surface: R = 9/8 in units of r_s.
        speed = sp.lambdify(r, sp.sqrt(grr / lapse), "math")
        steps = 20000
        crossing = sum(speed((i + 0.5) / steps) for i in range(steps)) / steps * 9 / 8
        self.assertAlmostEqual(crossing, 8.28, places=2)
        self.assertAlmostEqual(float(sp.N(sp.sqrt(lapse / grr).subs(r, 0), 20)), 0.134, places=3)


@unittest.skipUnless(HAS_SYMPY, "sympy is not installed")
class TheHeldNames(unittest.TestCase):
    def test_the_declared_slopes_are_those_of_the_definitions(self):
        import sympy as sp
        import verify_metrics as vm
        metric = json.loads((DATA / "metrics" / "bowers_liang.json").read_text(encoding="utf-8"))
        entry = metric["coordinates"][0]
        key = ("bowers_liang", "areal")
        # The reader refuses a declared slope that is not its definition's.
        reader = vm.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]], (), held=vm.HELD[key],
                           rates=vm.RATES[key])
        r = reader.symbol["r"]
        N, P = reader.parameters["N"], reader.parameters["P"]
        self.assertEqual(set(reader.rates), {N, P})
        # P is Z^Q/N, so Z^Q = N P and N' = Q r_s r Z^(Q - 1)/R^3.
        Z, Q = reader.parameters["Z"], reader.parameters["Q"]
        written = reader.held[P] - Z ** Q / reader.held[N]
        self.assertEqual(vm.norm(written), 0)
        self.assertEqual(vm.norm(reader.rates[N][r] - Q * (1 - Z) * N * P / (r * Z)), 0)


if __name__ == "__main__":
    unittest.main()
