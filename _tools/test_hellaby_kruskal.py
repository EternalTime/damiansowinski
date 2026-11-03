"""The regular neck and the shared slice of the symmetric dust example."""
import json
from pathlib import Path
import unittest

import sympy as sp

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "MFS/assets/data"


class DustNeck(unittest.TestCase):
    def test_finite_radial_metric_and_density_at_maximum_expansion(self):
        r = sp.symbols("r", real=True)
        mass = (1 + r**2 / 10)**3
        radius = 2 * mass * sp.exp(r**2)
        energy = -sp.exp(-r**2)
        radial_metric = sp.diff(radius, r)**2 / (1 + energy)
        density_8pi = 2 * sp.diff(mass, r) / (radius**2 * sp.diff(radius, r))
        self.assertEqual(sp.limit(radial_metric, r, 0), sp.Rational(676, 25))
        self.assertEqual(sp.limit(density_8pi, r, 0), sp.Rational(3, 26))
        self.assertEqual(sp.diff(radius, r).subs(r, 0), 0)
        self.assertGreater(sp.diff(radius, r, 2).subs(r, 0), 0)

    def test_light_crosses_the_neck_and_slice_matches_embedding(self):
        diagram = json.loads((DATA / "diagrams/hellaby_kruskal.json").read_text())
        view = diagram["systems"]["comoving"][0]
        for family in view["rays"].values():
            self.assertTrue(any(min(p[0] for p in ray) < 0.5 < max(p[0] for p in ray)
                                for ray in family))
        self.assertEqual(view["slices"][0]["view"], "neck")
        embedding = json.loads((DATA / "embedding/hellaby_kruskal.json").read_text())
        pieces = embedding["views"][0]["surfaces"][0]["pieces"]
        self.assertEqual(len(pieces), 2)
        for piece in pieces:
            self.assertEqual(piece["points"][0][:2], [0.0, 2.0])


if __name__ == "__main__":
    unittest.main()
