"""The four metric motions and the proper subgroup inherited by the Maxwell field.

The chart printer checks the source and both Maxwell equations.
These tests independently differentiate the published metric and the reconstructed two-form
along the four motion generators, including the generator which induces a duality rotation.
"""
import json
from pathlib import Path
import sys
import unittest

import sympy as sp

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "_tools" / "derivations"))
import verify_metrics as vm


class SymmetryInheritance(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.metric = json.loads((ROOT / "MFS/assets/data/metrics/mclenaghan_tariq_tupper.json").read_text())
        entry = next(c for c in cls.metric["coordinates"] if c["id"] == "exponential")
        cls.reader = vm.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]], ())
        R = cls.reader
        cls.X = [R.symbol[c] for c in entry["coords"]]
        t, r, y, z = cls.X
        cls.a = a = R.parameters["a"]
        components = {tuple(c["indices"]): R(c["value"]) for c in entry["metric_components"]}
        cls.g = sp.Matrix(4, 4, lambda i, j: components.get((entry["coords"][i], entry["coords"][j]), 0))
        e = sp.Matrix([[1, 0, 0, -2*y/a], [0, 1, 0, 0],
                       [0, 0, sp.exp(-r/a), 0], [0, 0, 0, sp.exp(r/a)]])
        f = sp.zeros(4)
        f[0, 1], f[2, 3] = sp.sqrt(2)*sp.cos(2*r/a)/a, sp.sqrt(2)*sp.sin(2*r/a)/a
        f[1, 0], f[3, 2] = -f[0, 1], -f[2, 3]
        cls.F = e.T*f*e
        cls.generators = [sp.Matrix(v) for v in ((1, 0, 0, 0), (0, 0, 0, 1),
                                                (2*z/a, 0, 1, 0), (0, 1, y/a, -z/a))]
        cls.points = [dict(zip([*cls.X, a], p)) for p in
                      ((0.3, -0.8, 0.4, -0.6, 1.2), (-0.1, 0.5, -0.7, 0.3, 0.9), (0.9, 1.3, 1.1, -0.4, 1.7))]

    def lie(self, tensor, vector):
        X = self.X
        return sp.Matrix(4, 4, lambda i, j: sum(
            vector[k]*sp.diff(tensor[i, j], X[k]) + tensor[k, j]*sp.diff(vector[k], X[i])
            + tensor[i, k]*sp.diff(vector[k], X[j]) for k in range(4)))

    def assert_zero_at_points(self, tensor):
        for at in self.points:
            self.assertLess(max(abs(complex(sp.N(v.subs(at), 30))) for v in tensor), 1e-10)

    def test_four_independent_generators_preserve_the_published_metric(self):
        for vector in self.generators:
            self.assert_zero_at_points(self.lie(self.g, vector))
        matrix = sp.Matrix.hstack(*self.generators)
        for at in self.points:
            self.assertEqual(matrix.subs(at).rank(), 4)

    def test_exactly_three_generators_inherit_the_maxwell_field(self):
        for vector in self.generators[:3]:
            self.assert_zero_at_points(self.lie(self.F, vector))
        up = self.g.inv()*self.F*self.g.inv()
        dual = sp.Matrix(4, 4, lambda i, j: sum(
            sp.LeviCivita(i, j, k, l)*up[k, l]/2 for k in range(4) for l in range(4)))
        change = self.lie(self.F, self.generators[3])
        self.assert_zero_at_points(change + 2*dual/self.a)
        for at in self.points:
            self.assertGreater(max(abs(complex(sp.N(v.subs(at), 30))) for v in change), 0.1)

    def test_upper_half_plane_has_no_periodic_identification(self):
        chart = next(c for c in self.metric["coordinates"] if c["id"] == "upper_half_plane")
        self.assertIn(r"\phi \in (-\infty, \infty)", chart["domains"])
        self.assertIn("dimensionless real coordinate", chart["convention"])


if __name__ == "__main__":
    unittest.main()
