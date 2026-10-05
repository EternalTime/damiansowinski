"""Haggard and Rovelli's black hole fireworks, a shell of light that falls into a black hole and comes
back out of a white hole: .venv.noindex/bin/python -m unittest discover -s _tools

What its texts state, held on the construction black_to_white_hole.py draws and on the published
files. The shell of light is in no chart's mathematics, so the junction is taken here from its
definitions: flat space inside, ds^2 = -du dv + ((v - u)/2)^2 dOmega^2, Kruskal's chart outside,
(1 - r) e^r = UV at r_s = 1, the shell along v = 0 and V = V0, so that the outgoing ray u of the
interior crosses into U = (1 + u/2) e^(-u/2)/V0. The edge of the quantum region is the spacelike
geodesic from Delta, at r = 7/6 on U + V = 0, to the point E of the shell at r = 1/2, which keeps
the norm of its tangent and its energy along the Killing vector V d_V - U d_U. The drawings are held
to the same numbers: the shell crosses U + V = 0 at r = 1.12, crosses r_s at the Painleve-Gullstrand
ct = -1.41 and reaches E at -1.18, Delta is at ct = -1.10, and the drawing's q puts the surface of
time symmetry on T = 0 with the geodesic below it. The published charts are held to Schwarzschild's
metric pulled back. Tests that need numpy and scipy are skipped where they are absent, and those of
the published mathematics where sympy is.
"""
import importlib.util
import json
import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "MFS" / "assets" / "data"
HAS_NUMERIC = all(importlib.util.find_spec(m) is not None for m in ("numpy", "scipy"))
HAS_SYMPY = importlib.util.find_spec("sympy") is not None
sys.path.insert(0, str(Path(__file__).resolve().parent / "derivations"))


def metric():
    return json.loads((DATA / "metrics" / "black_to_white_hole.json").read_text(encoding="utf-8"))


def chart(system):
    return next(c for c in metric()["coordinates"] if c["id"] == system)


def captions(kind):
    data = json.loads((DATA / kind / "black_to_white_hole.json").read_text(encoding="utf-8"))
    if kind == "diagrams":
        return {f"{s}/{v['id']}": " ".join(v["caption"]) for s, views in data["systems"].items() for v in views}
    return {v["id"]: " ".join(v["caption"]) for v in data["views"]}


@unittest.skipUnless(HAS_NUMERIC, "numpy and scipy are not installed")
class Construction(unittest.TestCase):
    def setUp(self):
        import numpy as np
        import black_to_white_hole as b
        self.np, self.b = np, b

    def test_the_shell_has_one_radius_on_both_sides(self):
        np, b = self.np, self.b
        u = -np.geomspace(1e-3, 30, 200)
        self.assertLess(float(np.max(np.abs(b.radius(b.shell_U(u), b.V0) + u / 2))), 1e-10)
        self.assertLess(float(np.max(np.abs(b.interior_u(b.shell_U(u)) - u))), 1e-10)

    def test_E_and_Delta_are_where_haggard_and_rovelli_put_them(self):
        b = self.b
        self.assertAlmostEqual(float(b.radius(b.U_E, b.V0)), b.EPS, places=12)
        self.assertAlmostEqual(float(b.radius(-b.V_DELTA, b.V_DELTA)), 7 / 6, places=12)
        self.assertAlmostEqual(float(b.shell_U(-2 * b.EPS)), b.U_E, places=12)
        self.assertAlmostEqual(float(b.interior_u(0.0)), -2.0, places=12)

    def test_the_edge_of_the_quantum_region_is_a_spacelike_geodesic(self):
        np, b = self.np, self.b
        U, V, r = b.geodesic()
        self.assertAlmostEqual(U[0], -b.V_DELTA, places=12)
        self.assertAlmostEqual(V[0], b.V_DELTA, places=12)
        self.assertAlmostEqual(U[-1], b.U_E, places=9)
        self.assertAlmostEqual(V[-1], b.V0, places=9)
        dU, dV = np.gradient(U), np.gradient(V)
        self.assertTrue(np.all(dU > 0) and np.all(dV < 0))
        self.assertTrue(np.all(np.diff(r) < 0))
        F = 4 / r * np.exp(-r)
        norm = F * -dU * dV
        energy = F / 2 * (V * dU - U * dV)
        inner = slice(5, -5)
        for name, q in (("norm", norm), ("energy", energy)):
            spread = float(np.ptp(q[inner]) / np.mean(np.abs(q[inner])))
            self.assertLess(spread, 1e-5, name)

    def test_the_shell_crosses_the_surface_of_time_symmetry_inside_Delta(self):
        b = self.b
        a = b.shell_crossing()
        self.assertAlmostEqual(a, 1.1177, places=4)
        self.assertLess(a, b.R_DELTA)
        self.assertAlmostEqual(float(b.shell_radius_schwarzschild(0.0)), a, places=10)

    def test_the_painleve_gullstrand_chart_is_schwarzschilds_outside_r_s(self):
        np, b = self.np, self.b
        r = np.linspace(1.05, 6, 50)
        t = np.linspace(-3, 2, 50)
        T = t - 2 * np.sqrt(r) - np.log((np.sqrt(r) - 1) / (np.sqrt(r) + 1))
        for x, y in zip(b.painleve_UV(t, r), b.schwarzschild_UV(T, r)):
            self.assertLess(float(np.max(np.abs(x / y - 1))), 1e-12)
        self.assertLess(float(np.max(np.abs(b.painleve_time(*b.painleve_UV(t, r)) - t))), 1e-12)

    def test_the_drawings_place_the_shell_and_the_edge(self):
        np, b = self.np, self.b
        self.assertAlmostEqual(float(b.shell_painleve_time(1.0)), -1.408, places=3)
        self.assertAlmostEqual(float(b.shell_painleve_time(b.EPS)), -1.177, places=3)
        self.assertAlmostEqual(float(b.edge_painleve_time(b.R_DELTA)), -1.096, places=3)
        far = np.geomspace(b.V_DELTA, 60, 100)
        p, q = b.conformal_two(-far, far)
        self.assertLess(float(np.max(np.abs(p + q))), 1e-12)
        U, V, _ = b.geodesic()
        p, q = b.conformal_two(U, V)
        self.assertLessEqual(float(np.max(p + q)), 1e-12)
        self.assertAlmostEqual(float(p[-1]), -math.pi / 4, places=9)
        self.assertAlmostEqual(float(q[-1]), 0.0, places=12)
        u = -np.geomspace(1.0, 40, 100)
        p1, q1 = b.conformal_one(u, 0 * u)
        p2, q2 = b.conformal_two(b.shell_U(u), b.V0 + 0 * u)
        self.assertLess(float(np.max(np.abs(p1 - p2)) + np.max(np.abs(q1 - q2))), 1e-10)
        Vs = np.linspace(b.V0, 4, 400)
        self.assertTrue(np.all(np.diff(b.q_of_V(Vs)) > 0))

    def test_the_flaps_overlap_in_kruskals_spacetime_beside_Delta(self):
        np, b = self.np, self.b
        U, V = np.meshgrid(np.linspace(-1, 2, 301), np.linspace(0.5, 2, 301))
        both = (b.region_two(U, V) > 0) & (b.region_two(-V, -U) > 0)
        self.assertGreater(int(both.sum()), 0)
        r = b.radius(U[both], V[both])
        self.assertGreater(float(r.min()), b.shell_crossing() - 1e-3)
        self.assertLess(float(r.max()), b.R_DELTA)


class Words(unittest.TestCase):
    def test_the_captions_state_the_numbers_the_construction_gives(self):
        flat = captions("diagrams")
        self.assertIn("$r = 1.12\\,r_s$", flat["schwarzschild/radial"])
        self.assertIn("$ct = -1.41\\,r_s$", flat["painleve_gullstrand_ingoing/falling"])
        self.assertIn("$ct = -1.18\\,r_s$", flat["painleve_gullstrand_ingoing/falling"])
        self.assertIn("$ct = -1.10\\,r_s$", flat["painleve_gullstrand_ingoing/falling"])
        self.assertIn("$ct = 1.41\\,r_s$", flat["painleve_gullstrand_outgoing/rising"])
        self.assertIn("$r = 1.12\\,r_s$", flat["kruskal/rising"])
        for text in captions("conformal").values():
            self.assertIn("no event horizon", text)

    def test_every_chart_is_a_piece_of_minkowski_or_schwarzschild(self):
        ids = [c["id"] for c in metric()["coordinates"]]
        self.assertEqual(ids, ["interior", "kruskal", "schwarzschild", "painleve_gullstrand_ingoing",
                               "painleve_gullstrand_outgoing", "lemaitre"])
        for c in metric()["coordinates"]:
            ricci = c["ricci_tensor"]["variants"][c["ricci_tensor"]["default"]]["nonzero"]
            self.assertEqual(ricci, [], c["id"])
        self.assertEqual(chart("interior")["kretschmann"], "K = 0")
        for system in ("kruskal", "schwarzschild", "painleve_gullstrand_ingoing", "painleve_gullstrand_outgoing",
                       "lemaitre"):
            self.assertEqual(chart(system)["kretschmann"], "K = \\dfrac{12r_s^2}{r^6}", system)


@unittest.skipUnless(HAS_SYMPY, "sympy is not installed")
class Pullbacks(unittest.TestCase):
    def setUp(self):
        import sympy
        import verify_metrics
        self.sp, self.vm = sympy, verify_metrics

    def published(self, system, held=()):
        entry = chart(system)
        reader = self.vm.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]], (), held=held)
        g = self.sp.zeros(4, 4)
        for e in entry["metric_components"]:
            i, j = (entry["coords"].index(x) for x in e["indices"])
            g[i, j] = reader(e["value"])
        return entry, reader, g

    def test_the_painleve_gullstrand_charts_are_schwarzschild_pulled_back(self):
        sp = self.sp
        for system, sign in (("painleve_gullstrand_ingoing", 1), ("painleve_gullstrand_outgoing", -1)):
            _, reader, g = self.published(system)
            r, rs = reader.symbol["r"], reader.parameters["r_s"]
            f = 1 - rs / r
            # T = t - sign h(r), h' = sqrt(r_s/r)/f, in the chart x^0 = ct.
            J = sp.Matrix([[1, -sign * sp.sqrt(rs / r) / f], [0, 1]])
            pulled = J.T * sp.diag(-f, 1 / f) * J
            for i in range(2):
                for j in range(2):
                    self.assertEqual(self.vm.norm(pulled[i, j] - g[i, j]), 0, (system, i, j))

    def test_the_interior_is_flat_space_in_null_coordinates(self):
        sp = self.sp
        _, reader, g = self.published("interior")
        u, v = reader.symbol["u"], reader.symbol["v"]
        self.assertEqual(self.vm.norm(g[0, 1] + sp.Rational(1, 2)), 0)
        self.assertEqual(self.vm.norm(g[2, 2] - (v - u) ** 2 / 4), 0)

    def test_lemaitres_chart_holds_r_and_rho_minus_c_tau(self):
        entry = chart("lemaitre")
        values = [e["value"] for e in entry["christoffel"]["variants"]["ull"]["nonzero"]]
        self.assertTrue(any("\\rho - c\\,\\tau" in v for v in values))
        self.assertFalse(any("\\sqrt" in v for v in values))
        self.assertEqual([e["value"] for e in entry["metric_components"]], ["-1", "\\dfrac{r_s}{r}", "r^2",
                                                                               "r^2\\sin^2\\theta"])
