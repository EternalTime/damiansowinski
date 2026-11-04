"""Expansion, recollapse, and curvature of the closed string cosmology."""

import json
from pathlib import Path
import sys
import unittest

import sympy as sp

ROOT = Path(__file__).resolve().parent.parent
METRIC_FILE = ROOT / "MFS/assets/data/metrics/nappi_witten.json"
sys.path.insert(0, str(ROOT / "_tools/derivations"))
import verify_metrics as vm


class NappiWitten(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        page = json.loads(METRIC_FILE.read_text())
        cls.chart = next(c for c in page["coordinates"] if c["id"] == "cosmological")
        cls.reader = vm.Reader(cls.chart["coords"], [p["symbol"] for p in cls.chart["parameters"]], ())
        cls.metric = sp.Matrix(vm.metric_from_line_element(cls.reader, cls.chart["line_element"], cls.chart["coords"]))
        cls.psi, cls.s, _, _ = [cls.reader.symbol[c] for c in cls.chart["coords"]]
        cls.L, cls.b = (cls.reader.parameters[k] for k in ("L", "b"))

    def scalar(self, field):
        return self.reader(self.chart[field].split("=", 1)[1])

    def test_opposite_circles_collapse_at_the_two_endpoints(self):
        for latitude in (sp.Rational(1, 3), sp.Rational(2, 3)):
            at = {self.s: latitude, self.L: 1, self.b: 1}
            self.assertEqual(self.metric[3, 3].subs(at).subs(self.psi, 0), 0)
            self.assertEqual(self.metric[2, 2].subs(at).subs(self.psi, sp.pi/2), 0)
            self.assertGreater(float(self.metric[2, 2].subs(at).subs(self.psi, 0)), 0)
            self.assertGreater(float(self.metric[3, 3].subs(at).subs(self.psi, sp.pi/2)), 0)

    def test_curvature_is_finite_when_recollapse_begins(self):
        c = sp.cos(2*self.s)
        for field, expected in (("ricci_scalar", 14*c*c/self.L**2),
                                ("kretschmann", (76*c**4+128*c*c-128)/self.L**4)):
            at_maximum = self.scalar(field).subs({self.psi: sp.pi/4, self.b: 1})
            self.assertEqual(vm.norm(at_maximum-expected), 0)
            for latitude in (0, sp.pi/8, sp.pi/4, 3*sp.pi/8, sp.pi/2):
                value = at_maximum.subs({self.s: latitude, self.L: 1})
                self.assertTrue(value.is_finite)

    def test_spacetime_curvature_varies_across_the_round_spatial_slice(self):
        ricci = self.scalar("ricci_scalar").subs({self.psi: sp.pi/4, self.b: 1, self.L: 1})
        self.assertEqual(vm.norm(ricci.subs(self.s, 0)), 14)
        self.assertEqual(vm.norm(ricci.subs(self.s, sp.pi/4)), 0)
        self.assertEqual(vm.norm(ricci.subs(self.s, sp.pi/2)), 14)


if __name__ == "__main__":
    unittest.main()
