"""The Morgan-Morgan discs, static discs of counterrotating dust in Weyl's class:
python3 -m unittest discover -s _tools

What its texts and drawings state, held on the published files and on the disc every drawing
draws, the first of the family at m = a/5 in units of the disc's radius a. The first class needs
nothing but the files. The second needs numpy, scipy and sympy and is skipped where they are
absent: it holds the published psi and gamma to Weyl's field equations, the disc to the junction
conditions of a layer of dust in two streams, the speed of the streams to the limit m/a = 2/(3 pi),
the family to Morgan and Morgan's bound on the redshift, and the functions Weyl's chart is drawn
with to the published ones.
"""
import importlib.util
import json
import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "MFS" / "assets" / "data"
HAS_NUMBERS = all(importlib.util.find_spec(name) is not None for name in ("sympy", "numpy", "scipy", "contourpy"))
sys.path.insert(0, str(Path(__file__).resolve().parent / "derivations"))

M = 0.2                         # m/a of every drawing
K = 3 * math.pi * M / 4         # rho d_rho psi = K rho^2 on the disc


def load(folder):
    return json.loads((DATA / folder / "morgan_morgan.json").read_text(encoding="utf-8"))


class Drawings(unittest.TestCase):
    def test_the_embedded_rim_is_the_circle_the_potential_gives_it(self):
        """On the disc psi = -(3 pi m/8a)(2 - rho^2/a^2), so the rim, rho = a, is a circle of radius
        a e^(3 pi m/8a), where the disc and the plane outside it meet."""
        disc, plane = load("embedding")["views"][0]["surfaces"][0]["pieces"]
        want = math.exp(3 * math.pi * M / 8)
        for point in (disc["points"][-1], plane["points"][0]):
            self.assertAlmostEqual(point[0], 1.0)
            self.assertAlmostEqual(point[1], want, places=6)
        self.assertAlmostEqual(disc["points"][-1][2], plane["points"][0][2], places=6)

    def test_the_embedded_cap_is_level_at_its_centre_and_steepest_at_the_rim(self):
        disc, plane = load("embedding")["views"][0]["surfaces"][0]["pieces"]

        def slopes(piece):
            pts = piece["points"]
            return [(b[2] - a[2]) / (b[1] - a[1]) for a, b in zip(pts, pts[1:])]
        inner, outer = slopes(disc), slopes(plane)
        self.assertLess(inner[0], 0.05)
        self.assertEqual(max(inner), inner[-1])
        self.assertEqual(max(outer), outer[0])
        self.assertAlmostEqual(inner[-1], outer[0], delta=0.05)

    def test_the_cones_are_narrowest_at_the_centre_of_the_disc(self):
        """On the axis c dt/dz = e^(-2 psi), greatest at z = 0, where it is e^(3 pi m/2a) = 2.57."""
        view = next(v for v in load("diagrams")["systems"]["weyl"] if v["id"] == "axis")
        X0, X1, Y0, Y1 = view["box"]
        steepest = 0.0
        for cone in view["cones"]:
            dx, dy = cone["b"][0] * (X1 - X0), cone["b"][1] * (Y1 - Y0)
            steepest = max(steepest, abs(dy / dx))
            z = X0 + cone["at"][0] * (X1 - X0)
            if abs(z) > 2.5:
                self.assertLess(abs(dy / dx), 1.5)
        self.assertLessEqual(steepest, math.exp(2 * K) + 1e-3)
        self.assertGreater(steepest, 2.0)

    def test_the_rim_is_marked_singular_where_the_oblate_chart_draws_it(self):
        views = {v["id"]: v for v in load("diagrams")["systems"]["oblate_spheroidal"]}
        for name in ("plane", "disc"):
            self.assertIn({"kind": "singular", "edges": ["left"]}, views[name]["markers"])
        self.assertFalse(any(m["kind"] == "singular" for m in views["axis"]["markers"]))
        conformal = {v["id"]: v for v in load("conformal")["views"]}
        for name, singular in (("weyl_axis", False), ("oblate_axis", False), ("weyl_plane", True),
                               ("oblate_plane", True), ("oblate_disc", True)):
            found = any(layer.get("class") == "singular" for layer in conformal[name]["layers"])
            self.assertEqual(found, singular, name)


@unittest.skipUnless(HAS_NUMBERS, "needs numpy, scipy and sympy")
class Physics(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import numpy as np
        import sympy as sp

        import null_rays as nr
        import verify_metrics as vm
        cls.np, cls.sp, cls.nr, cls.vm = np, sp, nr, vm
        entry = next(c for c in load("metrics")["coordinates"] if c["id"] == "oblate_spheroidal")
        cls.entry = entry
        cls.reader = vm.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]], (),
                               held=vm.HELD[("morgan_morgan", "oblate_spheroidal")])
        r = cls.reader
        cls.xi, cls.eta = r.symbol["\\xi"], r.symbol["\\eta"]
        cls.m, cls.a = r.parameters["m"], r.parameters["a"]
        cls.psi, cls.gamma = r.held[r.parameters["psi"]], r.held[r.parameters["gamma"]]

    def test_the_published_functions_solve_weyls_equations(self):
        sp, x, y, P, G = self.sp, self.xi, self.eta, self.psi, self.gamma
        laplace = sp.diff((1 + x ** 2) * sp.diff(P, x), x) + sp.diff((1 - y ** 2) * sp.diff(P, y), y)
        self.assertEqual(self.vm.norm(laplace), 0)
        px, py = sp.diff(P, x), sp.diff(P, y)
        along = (1 - y ** 2) / (x ** 2 + y ** 2) * (x * (1 + x ** 2) * px ** 2 - x * (1 - y ** 2) * py ** 2
                                                      - 2 * y * (1 + x ** 2) * px * py)
        across = (1 + x ** 2) / (x ** 2 + y ** 2) * (y * (1 + x ** 2) * px ** 2 - y * (1 - y ** 2) * py ** 2
                                                       + 2 * x * (1 - y ** 2) * px * py)
        self.assertEqual(self.vm.norm(sp.diff(G, x) - along), 0)
        self.assertEqual(self.vm.norm(sp.diff(G, y) - across), 0)
        self.assertEqual(self.vm.norm(G.subs(y, 1)), 0)

    def test_the_disc_is_dust_in_two_streams(self):
        """Israel's junction conditions across z = 0 of Weyl's metric, for a field even in z: with
        n = e^(psi - gamma) d_z and K_ab = n^z d_z g_ab / 2 on the upper face, the layer's stress is
        S^a_b = -([K^a_b] - delta^a_b [K])/(8 pi G/c^4) with [K] twice the upper face's. Its component
        along rho vanishes by gamma's quadrature, and the ratio of the stress round the axis to the
        energy density is V^2/c^2 = rho d_rho psi/(1 - rho d_rho psi), the speed of two streams of dust."""
        sp = self.sp
        rho, z = sp.symbols("rho z", positive=True)
        psi, gam = sp.Function("psi")(rho, z), sp.Function("gamma")(rho, z)
        g = {"t": -sp.exp(2 * psi), "rho": sp.exp(2 * gam - 2 * psi), "phi": rho ** 2 * sp.exp(-2 * psi)}
        normal = sp.exp(psi - gam)
        mixed = {k: normal * sp.diff(v, z) / (2 * v) for k, v in g.items()}
        trace = sum(mixed.values())
        stress = {k: -(v - trace) for k, v in mixed.items()}      # times c^4/(4 pi G), the jump being twice this
        on_shell = {sp.Derivative(gam, z): 2 * rho * sp.Derivative(psi, rho) * sp.Derivative(psi, z)}
        stress = {k: sp.simplify(v.subs(on_shell)) for k, v in stress.items()}
        self.assertEqual(stress["rho"], 0)
        u = rho * sp.Derivative(psi, rho)
        self.assertEqual(sp.simplify(stress["phi"] / (-stress["t"]) - u / (1 - u)), 0)
        # The first disc: d_z psi = 3 m eta/a^2 on the upper face and rho d_rho psi = (3 pi m/4a^3) rho^2.
        x, y, m, a = self.xi, self.eta, self.m, self.a
        self.assertEqual(self.vm.norm(sp.diff(self.psi, x).subs(x, 0) / (a * y) - 3 * m * y / a ** 2), 0)
        on_disc = self.psi.subs(x, 0)
        turn = -(1 - y ** 2) / y * sp.diff(on_disc, y)
        self.assertEqual(self.vm.norm(turn - 3 * sp.pi * m * (1 - y ** 2) / (4 * a)), 0)

    def test_the_dust_reaches_the_speed_of_light_at_two_over_three_pi(self):
        """V = c where rho d_rho psi = 1/2, first at the rim, at m/a = 2/(3 pi); there light from the
        centre has the redshift e^(1/2) - 1. At m = a/5 the rim moves at 0.94 c and the redshift is 0.60."""
        limit = 2 / (3 * math.pi)
        self.assertAlmostEqual(3 * math.pi * limit / 4, 0.5)
        self.assertAlmostEqual(math.exp(3 * math.pi * limit / 4) - 1, math.exp(0.5) - 1)
        self.assertLess(M, limit)
        self.assertAlmostEqual(math.sqrt(K / (1 - K)), 0.944, places=3)
        self.assertAlmostEqual(math.exp(K) - 1, 0.602, places=3)
        psi_centre = float(self.psi.subs({self.xi: 0, self.eta: 1, self.m: M, self.a: 1}))
        self.assertAlmostEqual(psi_centre, -K)

    def test_the_family_reaches_morgan_and_morgans_redshift(self):
        """The redshift at which a stream first reaches c grows with the order of the disc. As the order
        grows without bound the density is a Gaussian of width s, psi = -sqrt(pi/2)(m/s) e^(-y) I_0(y)
        with y = rho^2/(4 s^2), and rho d_rho psi = 1/2 at its greatest caps the central redshift at
        z = 1.5803, the bound of Morgan and Morgan's abstract. The first disc alone stops at 0.6487."""
        from scipy.optimize import minimize_scalar
        from scipy.special import ive
        best = minimize_scalar(lambda y: -2 * y * (ive(0, y) - ive(1, y)), bounds=(0.1, 3), method="bounded")
        self.assertAlmostEqual(math.exp(1 / (2 * -best.fun)) - 1, 1.5803, places=3)
        self.assertAlmostEqual(math.exp(0.5) - 1, 0.6487, places=4)

    def test_the_kretschmann_scalar_diverges_at_the_rim_as_the_inverse_distance(self):
        """K (xi^2 + eta^2) settles on one value at the rim from the plane outside and from the disc,
        and rho - a = a xi^2/2 outside, so K grows as the inverse of the distance; on the axis it is finite."""
        sp, r = self.sp, self.reader
        K = r.surface(r(self.entry["kretschmann"].split("=", 1)[1])).subs({r.c: 1, self.m: sp.Rational(1, 5), self.a: 1})
        f = sp.lambdify((self.xi, self.eta), K, "mpmath")
        outside = [float(f(sp.Float(10) ** -k, 0) * sp.Float(10) ** (-2 * k)) for k in (3, 5)]
        on_disc = [float(f(0, sp.Float(10) ** -k) * sp.Float(10) ** (-2 * k)) for k in (3, 5)]
        self.assertAlmostEqual(outside[0] / outside[1], 1.0, places=2)
        self.assertAlmostEqual(on_disc[1] / outside[1], 1.0, places=2)
        self.assertGreater(outside[1], 1.0)
        self.assertLess(abs(float(f(sp.Float("1e-6"), 1))), 1.0)

    def test_weyls_chart_is_drawn_with_the_published_functions(self):
        """The strings the drawings of Weyl's chart declare, in rho and z, are the published psi and
        gamma at rho = a sqrt((1 + xi^2)(1 - eta^2)), z = a xi eta, and the plane's own strings and
        the closed forms in floats agree with them."""
        sp, np, nr = self.sp, self.np, self.nr
        rho, z = sp.symbols("rho z", real=True)
        at = {self.m: sp.Rational(1, 5), self.a: 1}
        general = {k: sp.lambdify((rho, z), sp.sympify(v, locals={"rho": rho, "z": z}), "numpy")
                   for k, v in nr._morgan_morgan_weyl().items()}
        plane = {k: sp.lambdify(rho, sp.sympify(v, locals={"rho": rho}), "numpy") for k, v in nr._morgan_morgan_plane().items()}
        published = {"psi": sp.lambdify((self.xi, self.eta), self.psi.subs(at), "numpy"),
                     "gamma": sp.lambdify((self.xi, self.eta), self.gamma.subs(at), "numpy")}
        rng = np.random.default_rng(3)
        for x, y in zip(rng.uniform(0.0, 3.0, 40), rng.uniform(-1.0, 1.0, 40)):
            r, h = math.sqrt((1 + x * x) * (1 - y * y)), x * y
            floats = dict(zip(("psi", "gamma"), nr._morgan_morgan_numbers(x, abs(y))))
            for name in ("psi", "gamma"):
                self.assertAlmostEqual(float(general[name](r, h)), float(published[name](x, y)), places=10)
                self.assertAlmostEqual(floats[name], float(published[name](x, y)), places=12)
        for r in (0.0, 0.3, 0.99, 1.0, 1.01, 2.5):
            x, y = (0.0, math.sqrt(1 - r * r)) if r <= 1 else (math.sqrt(r * r - 1), 0.0)
            for name in ("psi", "gamma"):
                self.assertAlmostEqual(float(plane[name](r)), float(published[name](x, y)), places=10)


if __name__ == "__main__":
    unittest.main()
