"""Global lattice, holonomy and modulus checks of the published torus cosmology.

Local curvature and tensor identities belong to the scoped metric checker.
These checks retain the independent global quotient data.
"""
import json
import math
import sys
import unittest
from pathlib import Path

import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "MFS" / "assets" / "data"
sys.path.insert(0, str(ROOT / "_tools" / "derivations"))
import verify_metrics as vm


def published(system):
    metric = json.loads((DATA / "metrics" / "torus_universe_2plus1.json").read_text())
    entry = next(c for c in metric["coordinates"] if c["id"] == system)
    reader = vm.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]], ())
    matrix = sp.zeros(3)
    for component in entry["metric_components"]:
        i, j = (entry["coords"].index(x) for x in component["indices"])
        matrix[i, j] = reader(component["value"]).subs(reader.c, 1)
    return entry, reader, matrix


def spatial_metric(t, a, b, lam, mu):
    _, R, matrix = published("proper_time")
    at = {R.symbol["t"]: t, **{R.parameters[n]: v for n, v in
                              zip(("a", "b", "lambda", "mu"), (a, b, lam, mu))}}
    return np.array(matrix[1:, 1:].subs(at), dtype=float)


class TorusModuli(unittest.TestCase):
    def test_area_and_modulus_from_the_published_gram_matrix(self):
        for a, b, lam, mu in ((1.2, .4, .8, 1.7), (2.1, -.7, -.3, .9), (.8, .6, 0, 1.4)):
            delta = a*mu-b*lam
            for t in (.3, .8, 1.7, 3.2):
                G = spatial_metric(t, a, b, lam, mu)
                area = math.sqrt(np.linalg.det(G))
                modulus = complex(G[0, 1], area)/G[0, 0]
                self.assertAlmostEqual(area, t*delta, places=11)
                self.assertAlmostEqual(abs(modulus-(b+1j*t*mu)/(a+1j*t*lam)), 0, places=12)
                self.assertGreater(modulus.imag, 0)

    def test_shape_follows_a_hyperbolic_geodesic(self):
        a, b, lam, mu = 1.2, .4, .8, 1.7
        centre = (b/a+mu/lam)/2
        radius = abs(b/a-mu/lam)/2
        for t in (.3, .8, 1.7, 3.2):
            G = spatial_metric(t, a, b, lam, mu)
            m = complex(G[0, 1], math.sqrt(np.linalg.det(G)))/G[0, 0]
            self.assertAlmostEqual(abs(m-centre), radius, places=12)
        self.assertGreater(abs((b+.3j*mu)/(a+.3j*lam)-(b+3.2j*mu)/(a+3.2j*lam)), .5)

    def test_integer_changes_of_basis_preserve_area_and_geometry(self):
        t, a, b, lam, mu = 1.4, 1.2, .4, .8, 1.7
        G = spatial_metric(t, a, b, lam, mu)
        modulus = (b+1j*t*mu)/(a+1j*t*lam)
        for M in (np.array([[1, 2], [0, 1]]), np.array([[0, -1], [1, 0]]), np.array([[2, 1], [1, 1]])):
            p, r, q, s = M.ravel()
            aa, bb = np.array([a, b]) @ M
            ll, mm = np.array([lam, mu]) @ M
            transformed = spatial_metric(t, aa, bb, ll, mm)
            np.testing.assert_allclose(transformed, M.T @ G @ M, rtol=1e-12, atol=1e-12)
            self.assertAlmostEqual(np.linalg.det(transformed), np.linalg.det(G), places=10)
            self.assertAlmostEqual(abs((bb+1j*t*mm)/(aa+1j*t*ll)-(r+s*modulus)/(p+q*modulus)), 0, places=12)

    def test_static_family_retains_the_free_lengths_and_angle(self):
        _, R, matrix = published("static")
        at = {R.parameters["r_1"]: 1.2, R.parameters["r_2"]: 1.6, R.parameters["phi"]: .9}
        G = np.array(matrix[1:, 1:].subs(at), dtype=float)
        area = math.sqrt(np.linalg.det(G))
        modulus = complex(G[0, 1], area)/G[0, 0]
        self.assertAlmostEqual(area, 1.2*1.6*math.sin(.9), places=12)
        self.assertAlmostEqual(abs(modulus-(1.6/1.2)*complex(math.cos(.9), math.sin(.9))), 0, places=12)
        self.assertFalse(matrix.free_symbols & set(R.symbol.values()))


class HolonomyAndBoundary(unittest.TestCase):
    def test_the_two_deck_generators_commute_and_match_the_development(self):
        def boost(P, rapidity, shift):
            T, X, Y = P
            return np.array([T*math.cosh(rapidity)+X*math.sinh(rapidity),
                             X*math.cosh(rapidity)+T*math.sinh(rapidity), Y+shift])

        a, b, lam, mu = 1.2, .4, .8, 1.7
        for t, x, y in ((.7, .2, .4), (1.3, .8, .1), (2.4, .3, .6)):
            for null_chart in (False, True):
                def developed(x, y):
                    q = lam*x+mu*y
                    if null_chart:
                        return np.array([(math.exp(q)+t*t*math.exp(-q))/2,
                                         (math.exp(q)-t*t*math.exp(-q))/2, a*x+b*y])
                    return np.array([t*math.cosh(q), t*math.sinh(q), a*x+b*y])
                P = developed(x, y)
                np.testing.assert_allclose(boost(P, lam, a), developed(x+1, y), rtol=1e-12, atol=1e-12)
                np.testing.assert_allclose(boost(P, mu, b), developed(x, y+1), rtol=1e-12, atol=1e-12)
                np.testing.assert_allclose(boost(boost(P, lam, a), mu, b),
                                           boost(boost(P, mu, b), lam, a), rtol=1e-12, atol=1e-12)

    def test_null_extension_has_timelike_closed_torus_circuits(self):
        entry, R, matrix = published("null_boundary")
        self.assertIn(r"u \in (0,\infty)", entry["domains"])
        at = {R.parameters["a"]: 1, R.parameters["b"]: sp.Rational(1, 3),
              R.parameters["lambda"]: 1, R.parameters["mu"]: sp.sqrt(2)}
        for u in (.2, 1.4, -.2, -1.4):
            g = np.array(matrix.subs(at).subs(R.symbol["u"], u), dtype=float)
            self.assertEqual(np.linalg.det(g[1:, 1:]) > 0, u > 0)
            if u < 0:
                winding = np.array([0, 1, -3])
                self.assertLess(winding @ g @ winding, 0)


class QuotientDrawings(unittest.TestCase):
    def test_every_chart_retains_both_identification_pairs(self):
        data = json.loads((DATA / "diagrams" / "torus_universe_2plus1.json").read_text())
        self.assertEqual(set(data["projections"]), {"proper_time", "mean_curvature", "minkowski", "null_boundary", "static"})
        for system, views in data["projections"].items():
            figure = views[0]
            for cls in ("above", "below"):
                lines = [np.array(L["points"]) for L in figure["turn"]["lines"] if L["class"] == cls]
                self.assertEqual(len(lines), 4, system)
                for P, Q in ((lines[0], lines[1]), (lines[2], lines[3])):
                    if system == "minkowski":
                        rapidity, shift = (1, 1) if cls == "above" else (math.sqrt(2), 1/3)
                        X, Y, T = P.T
                        image = np.column_stack([X*math.cosh(rapidity)+T*math.sinh(rapidity),
                                                 Y+shift, T*math.cosh(rapidity)+X*math.sinh(rapidity)])
                    else:
                        z = P[0, 2]
                        t = -1/z if system == "mean_curvature" else math.sqrt(2*z) if system == "null_boundary" else z
                        if system == "static":
                            shift = (1, 0, 0) if cls == "above" else (.55, 1.1*math.sin(math.pi/3), 0)
                        else:
                            shift = (1, t, 0) if cls == "above" else (1/3, t*math.sqrt(2), 0)
                        image = P+shift
                    np.testing.assert_allclose(Q, image, rtol=0, atol=5e-6, err_msg=system+" "+cls)


if __name__ == "__main__":
    unittest.main()
