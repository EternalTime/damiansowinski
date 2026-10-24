"""The universe on a wall moving past a black hole in five dimensions:
.venv.noindex/bin/python -m unittest discover -s _tools

What its texts state about the wall, held on the published files. The Gaussian normal chart stands
the wall at w = 0, so Israel's junction conditions, with the mirror symmetry, read the matter on the
wall off its published metric there: 2K^i_i = -kappa^2 rho_b/3 and 2K^t_t = kappa^2 (2 rho_b + 3 p_b)/3,
with K^i_i = d_w A/A and K^t_t = d_w N/N at w = 0 and rho_b = sigma + rho, p_b = p - sigma. The first is
Q = 1 + rho/sigma at sigma = 6/(kappa^2 l), which with the published Q is the Friedmann equation the
History states, and the second, with the first, keeps the matter on the wall conserved. The static
chart's wall moves by the junction condition its parameter states, and the closed forms of the wall
every diagram draws, brane_wall_universe.py, are held to both charts. The tests need sympy and numpy
and are skipped where either is absent.
"""
import importlib.util
import json
import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "MFS" / "assets" / "data"
READY = all(importlib.util.find_spec(name) is not None for name in ("sympy", "numpy", "scipy"))
sys.path.insert(0, str(Path(__file__).resolve().parent / "derivations"))


def chart(system):
    metric = json.loads((DATA / "metrics" / "brane_wall_universe.json").read_text(encoding="utf-8"))
    return next(c for c in metric["coordinates"] if c["id"] == system)


@unittest.skipUnless(READY, "sympy, numpy and scipy are not installed")
class Wall(unittest.TestCase):
    def setUp(self):
        import sympy
        import verify_metrics
        self.sp, self.vm = sympy, verify_metrics
        entry = chart("gaussian_normal")
        self.reader = verify_metrics.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]], ())
        self.g = {tuple(e["indices"]): self.reader(e["value"]) for e in entry["metric_components"]}
        self.t, self.w = self.reader.symbol["t"], self.reader.symbol["w"]
        self.a = self.reader.functions["a"][0]
        self.ell, self.mu, self.k = (self.reader.parameters[n] for n in ("ell", "mu", "k"))

    def written(self, expression):
        return self.reader.surface(self.sp.sympify(expression))

    def on_wall(self, expression):
        return self.vm.norm(self.written(expression).subs(self.w, 0))

    def curvatures(self):
        """K^i_i and K^t_t on the side w > 0 of the wall, from the published g_chichi and g_tt."""
        sp = self.sp
        A = sp.sqrt(self.written(self.g[("\\chi", "\\chi")]) * (1 - self.k * self.reader.symbol["\\chi"] ** 2))
        N = sp.sqrt(-self.written(self.g[("t", "t")]))
        return self.vm.norm((sp.diff(A, self.w) / A).subs(self.w, 0)), self.vm.norm((sp.diff(N, self.w) / N).subs(self.w, 0))

    def test_the_wall_is_the_surface_w_zero_with_its_own_proper_time(self):
        self.assertEqual(self.on_wall(self.g[("t", "t")]), -1)
        A2 = self.on_wall(self.g[("\\theta", "\\theta")] / self.reader.symbol["\\chi"] ** 2)
        self.assertEqual(self.vm.norm(A2 - self.a ** 2), 0)

    def test_the_junction_condition_is_q_equal_to_one_plus_rho_over_sigma(self):
        sp = self.sp
        spatial, _ = self.curvatures()
        kappa2, rho = sp.symbols("kappa2 rho", positive=True)
        sigma = 6 / (kappa2 * self.ell)
        # 2 K^i_i = -kappa^2 (sigma + rho)/3, solved for rho.
        solved = sp.solve(sp.Eq(2 * spatial, -kappa2 * (sigma + rho) / 3), rho)[0]
        Q = self.written(self.reader.parameters["Q"])
        self.assertEqual(self.vm.norm(solved / sigma - (Q - 1)), 0)

    def test_the_friedmann_equation_of_the_history(self):
        """Q^2 = (1 + rho/sigma)^2 with sigma = 3c^4/(4 pi G_5 l) and G_4 = G_5/l is
        H^2 = (8 pi G_4/3c^2) rho (1 + rho/2 sigma) - k c^2/a^2 + mu c^2/a^4, H = (da/d tau)/a."""
        sp = self.sp
        c, G5, rho, a, adot = sp.symbols("c G_5 rho a adot", positive=True)
        ell, mu, k = self.ell, self.mu, self.k
        sigma = 3 * c ** 4 / (4 * sp.pi * G5 * ell)
        G4 = G5 / ell
        Q = self.written(self.reader.parameters["Q"])
        Q = Q.subs({sp.Derivative(self.a, self.t): adot}).subs(self.a, a)
        # adot = da/d(c tau), so H^2 = c^2 adot^2/a^2, solved from Q^2 = (1 + rho/sigma)^2.
        H2 = c ** 2 * (sp.solve(sp.Eq(Q ** 2, (1 + rho / sigma) ** 2), adot)[0] ** 2) / a ** 2
        stated = 8 * sp.pi * G4 / (3 * c ** 2) * rho * (1 + rho / (2 * sigma)) - k * c ** 2 / a ** 2 + mu * c ** 2 / a ** 4
        self.assertEqual(sp.simplify(H2 - stated), 0)

    def test_the_matter_on_the_wall_is_conserved(self):
        """rho_b' + 3 H (rho_b + p_b) = 0, with rho_b and p_b read off both junction conditions."""
        sp = self.sp
        spatial, temporal = self.curvatures()
        kappa2 = sp.Symbol("kappa2", positive=True)
        rho_b = -6 * spatial / kappa2
        p_b = (6 * temporal / kappa2 - 2 * rho_b) / 3
        a = self.a
        balance = sp.diff(rho_b, self.t) + 3 * sp.diff(a, self.t) / a * (rho_b + p_b)
        self.assertEqual(self.vm.norm(balance), 0)

    def test_the_lapse_is_the_rate_of_the_areal_radius(self):
        """N = d_t A/a', Binetruy, Deffayet, Ellwanger and Langlois's (38), as the parameter states."""
        A, N = (self.written(self.reader.parameters[n]) for n in ("A", "N"))
        self.assertEqual(self.vm.norm(N - self.sp.diff(A, self.t) / self.sp.diff(self.a, self.t)), 0)


@unittest.skipUnless(READY, "sympy, numpy and scipy are not installed")
class DrawnWall(unittest.TestCase):
    """The wall every diagram draws, at l = mu = 1, k = 0 and rho = 0."""

    def setUp(self):
        import numpy
        import brane_wall_universe
        self.np, self.bwu = numpy, brane_wall_universe

    def test_the_wall_moves_by_the_junction_condition(self):
        np, bwu = self.np, self.bwu
        a = np.linspace(1.2, 6, 2000)
        # (da/d(c tau))^2 = a^2/l^2 - f(a) = mu/a^2, so d(c tau)/da = a, and f(a) dt/d(c tau) = a/l.
        self.assertTrue(np.allclose(bwu.speed_squared(a, 1, 1, 0), 1 / a ** 2))
        self.assertTrue(np.allclose(np.gradient(bwu.proper_time(a), a, edge_order=2), a, rtol=1e-9))
        # dt/da = (dt/d(c tau))(d(c tau)/da), a times the static rate.
        self.assertTrue(np.allclose(np.gradient(bwu.static_time(a), a, edge_order=2), a * bwu.static_rate(a, 1, 1, 0),
                                    rtol=2e-3))
        self.assertTrue(np.allclose(bwu.static_time(bwu.radius_at_static_time(np.array([-1.0, 0.0, 2.0]))),
                                    [-1.0, 0.0, 2.0], atol=1e-10))

    def test_the_gaussian_normal_forms_of_the_drawn_wall(self):
        import sympy as sp
        import verify_metrics as vm
        np, bwu = self.np, self.bwu
        entry = chart("gaussian_normal")
        reader = vm.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]], ())
        t, w = reader.symbol["t"], reader.symbol["w"]
        drawn = {reader.parameters["ell"]: 1, reader.parameters["mu"]: 1, reader.parameters["k"]: 0}
        A, N = (reader.surface(reader.parameters[n]).subs(drawn) for n in ("A", "N"))
        A, N = (e.subs(reader.functions["a"][0], sp.sqrt(2 * t)).doit() for e in (A, N))
        for tt, ww in ((0.3, 0.1), (1.0, 0.4), (2.0, 0.9)):
            at = {t: tt, w: ww}
            self.assertTrue(math.isclose(float(A.subs(at)) ** 2, 2 * tt * bwu.phi(tt, ww), rel_tol=1e-12))
            self.assertTrue(math.isclose(float(N.subs(at)), bwu.psi(tt, ww) / math.sqrt(bwu.phi(tt, ww)), rel_tol=1e-12))
        self.assertTrue(np.allclose(bwu.psi(np.array([0.5, 2.0]), bwu.focal_distance(np.array([0.5, 2.0]))), 0, atol=1e-12))


if __name__ == "__main__":
    unittest.main()
