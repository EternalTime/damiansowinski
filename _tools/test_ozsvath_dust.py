"""Ozsvath's homogeneous universes of dust, in four charts:
.venv.noindex/bin/python -m unittest discover -s _tools

What its texts and drawings state, held on the published files. The first class needs nothing
but the files. The second needs sympy and numpy and is skipped where they are absent: it reads
each chart's published metric and Einstein tensor and holds every chart to dust with a
cosmological constant at the density and the constant its parameters name, class III to Godel's
universe at s = 0 and to Godel's ratio K/R^2 = 3 at every s, Farnsworth and Kerr's form to class
III, class I to Einstein's static universe at beta = 0, and class II to class III at gamma = pi/2.
"""
import importlib.util
import json
import math
import random
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "MFS" / "assets" / "data"
HAS_NUMBERS = all(importlib.util.find_spec(name) is not None for name in ("sympy", "numpy"))
sys.path.insert(0, str(Path(__file__).resolve().parent / "derivations"))
ME = "ozsvath_dust"
CHARTS = ("class_one", "class_two", "class_three", "farnsworth_kerr")


def load(folder, name=ME):
    return json.loads((DATA / folder / f"{name}.json").read_text(encoding="utf-8"))


def chart(system, name=ME):
    return next(c for c in load("metrics", name)["coordinates"] if c["id"] == system)


def density_and_constant(system, p):
    """8 pi G rho/c^2 and Lambda, as ozsvath_dust.md Step 2 gives them."""
    if system == "class_one":
        L, b = p["L"], p["beta"]
        return math.cos(2 * b) / (L * L * (1 + math.cos(b) ** 2)), 1 / (2 * L * L * (1 + math.cos(b) ** 2))
    if system == "class_two":
        a, g = p["a"], p["gamma"]
        return -math.cos(2 * g) / (a * a * (1 + math.cos(g) ** 2)), -1 / (2 * a * a * (1 + math.cos(g) ** 2))
    return 1 / p["a"] ** 2, -1 / (2 * p["a"] ** 2)


class Drawings(unittest.TestCase):
    def test_each_drawing_is_there(self):
        diagrams = load("diagrams")
        self.assertEqual({s: [v["id"] for v in views] for s, views in diagrams["systems"].items()},
                         {"class_one": ["tpsi"], "class_two": ["xy"], "class_three": ["xy"], "farnsworth_kerr": ["tx"]})
        self.assertEqual([v["id"] for v in load("embedding")["views"]], ["torus"])
        self.assertEqual([v["id"] for v in load("curvature")["views"]], ["class_one", "class_two"])
        self.assertFalse((DATA / "conformal" / f"{ME}.json").exists())

    def test_the_torus_swells_and_shrinks_twice_in_each_turn(self):
        """Its circles have the radius L sqrt(1 + epsilon cos 2psi), epsilon = sin(pi/6)/sqrt 2, widest at
        psi = 0 and pi and narrowest at pi/2."""
        view, = load("embedding")["views"]
        points = view["surfaces"][0]["pieces"][0]["points"]
        eps = math.sin(math.pi / 6) / math.sqrt(2)
        radii = [p[1] for p in points]
        self.assertAlmostEqual(radii[0], math.sqrt(1 + eps), places=6)
        self.assertAlmostEqual(radii[-1], math.sqrt(1 + eps), places=6)
        self.assertAlmostEqual(min(radii), math.sqrt(1 - eps), places=4)
        narrowest = min(points, key=lambda p: p[1])
        self.assertAlmostEqual(narrowest[0], math.pi / 2, places=1)
        n = 20000
        rise = sum(math.sqrt(1 + math.sin(math.pi / 6) ** 2 - (eps * math.sin(2 * q)) ** 2 / (1 + eps * math.cos(2 * q)))
                   for q in ((k + 0.5) * math.pi / n for k in range(n))) * math.pi / n
        self.assertAlmostEqual(points[-1][2], rise, places=5)

    def test_class_one_curvature_is_einsteins_at_beta_zero_and_vanishes_at_the_edges(self):
        view = next(v for v in load("curvature")["views"] if v["id"] == "class_one")
        K = [p for c in view["curves"] if c["panel"] == "K" for p in c["points"]]
        peak = max(K, key=lambda p: p[1])
        self.assertAlmostEqual(peak[0], 0.0, places=2)
        self.assertAlmostEqual(peak[1], 0.75, places=6)
        self.assertLess(K[0][1], 1e-4)
        self.assertLess(K[-1][1], 1e-4)

    def test_class_two_curvature_ends_on_class_three(self):
        view = next(v for v in load("curvature")["views"] if v["id"] == "class_two")
        last = {c["panel"]: c["points"][-1] for c in view["curves"]}
        self.assertAlmostEqual(last["K"][0], math.pi / 2, places=9)
        self.assertAlmostEqual(last["K"][1], 3.0, places=9)
        self.assertAlmostEqual(last["R"][1], -1.0, places=9)


@unittest.skipUnless(HAS_NUMBERS, "needs sympy and numpy")
class Physics(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import numpy as np
        import sympy as sp
        import verify_metrics as vm
        cls.np, cls.sp, cls.vm = np, sp, vm

    def published(self, system, field, variant=None, name=ME):
        entry = chart(system, name)
        reader = self.vm.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]], ())
        listed = entry[field]["variants"][variant]["nonzero"] if variant else entry[field]
        values = {tuple(e["indices"]): e["value"] for e in listed}
        n = len(entry["coords"])
        matrix = self.sp.Matrix(n, n, lambda i, j: reader(values.get((entry["coords"][i], entry["coords"][j]), "0")))
        return matrix.subs(reader.c, 1), reader

    def scalar(self, system, field, name=ME):
        entry = chart(system, name)
        reader = self.vm.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]], ())
        return reader(entry[field].partition("=")[2] or entry[field]), reader

    def test_every_chart_is_dust_with_its_cosmological_constant(self):
        """G^mu_nu + Lambda delta = kappa rho u^mu u_nu: one eigenvalue -kappa rho on a timelike u and
        three zeros, at random events and parameters inside each chart's range of positive density."""
        np, sp = self.np, self.sp
        rng = random.Random(1965)
        draws = {"class_one": lambda: {"L": rng.uniform(0.5, 2), "beta": rng.uniform(-0.7, 0.7)},
                 "class_two": lambda: {"a": rng.uniform(0.5, 2), "gamma": rng.uniform(0.82, 1.5)},
                 "class_three": lambda: {"a": rng.uniform(0.5, 2), "s": rng.uniform(-0.9, 0.9)},
                 "farnsworth_kerr": lambda: {"a": rng.uniform(0.5, 2), "v": rng.uniform(-2, 2)}}
        for system in CHARTS:
            G, reader = self.published(system, "einstein_tensor", "ul")
            g, _ = self.published(system, "metric_components")
            for _ in range(4):
                p = draws[system]()
                at = {reader.parameters[k]: v for k, v in p.items()}
                at.update({s: rng.uniform(0.1, 1.2) for s in reader.symbol.values()})
                krho, lam = density_and_constant(system, p)
                M = np.array(G.subs(at).evalf(), dtype=float) + lam * np.eye(4)
                values, vectors = np.linalg.eig(M)
                k = int(np.argmin(values.real))
                self.assertAlmostEqual(values[k].real, -krho, places=9, msg=system)
                self.assertTrue(np.allclose(np.delete(values, k), 0, atol=1e-9), system)
                u = vectors[:, k].real
                self.assertLess(u @ np.array(g.subs(at).evalf(), dtype=float) @ u, 0, f"{system}: u is timelike")
                self.assertGreater(krho, 0)

    def test_class_three_at_s_zero_is_godels_metric(self):
        sp = self.sp
        g, reader = self.published("class_three", "metric_components")
        x, y, z, u = (reader.symbol[c] for c in ("x", "y", "z", "u"))
        a, s = reader.parameters["a"], reader.parameters["s"]
        w1 = sp.Matrix([1, 0, -sp.sinh(y), 0])
        godel = a ** 2 * (-2 * w1 * w1.T + sp.diag(0, 1, sp.cosh(y) ** 2, 1))
        self.assertEqual(sp.simplify(g.subs(s, 0) - godel), sp.zeros(4, 4))

    def test_class_three_and_farnsworth_kerr_have_godels_ratio_k_over_r_squared(self):
        sp = self.sp
        godel_K, _ = self.scalar("cartesian", "kretschmann", "godel")
        godel_R, _ = self.scalar("cartesian", "ricci_scalar", "godel")
        self.assertEqual(sp.simplify(godel_K / godel_R ** 2), 3)
        for system in ("class_three", "farnsworth_kerr"):
            K, _ = self.scalar(system, "kretschmann")
            R, _ = self.scalar(system, "ricci_scalar")
            self.assertEqual(sp.simplify(K / R ** 2), 3, system)

    def test_class_one_at_beta_zero_is_einsteins_static_universe_of_radius_2L(self):
        sp = self.sp
        K, reader = self.scalar("class_one", "kretschmann")
        R, _ = self.scalar("class_one", "ricci_scalar")
        beta, L = reader.parameters["beta"], reader.parameters["L"]
        einstein_K, er = self.scalar("hyperspherical", "kretschmann", "einstein_static")
        einstein_R, _ = self.scalar("hyperspherical", "ricci_scalar", "einstein_static")
        radius = {er.parameters["R"]: 2 * L}
        self.assertEqual(sp.simplify(K.subs(beta, 0) - einstein_K.subs(radius)), 0)
        self.assertEqual(sp.simplify(R.subs(beta, 0) - einstein_R.subs(radius)), 0)
        entry = chart("class_one")
        weyl = [reader(e["value"]).subs(beta, 0) for e in entry["weyl_tensor"]["variants"]["ulll"]["nonzero"]]
        self.assertTrue(weyl)
        self.assertTrue(all(sp.simplify(w) == 0 for w in weyl))

    def test_class_two_at_gamma_half_pi_has_class_threes_scalars(self):
        sp = self.sp
        for field in ("kretschmann", "ricci_scalar"):
            two, r2 = self.scalar("class_two", field)
            three, r3 = self.scalar("class_three", field)
            self.assertEqual(sp.simplify(two.subs(r2.parameters["gamma"], sp.pi / 2)
                                         - three.subs(r3.parameters["a"], r2.parameters["a"])), 0, field)


if __name__ == "__main__":
    unittest.main()
