"""The Einstein base, Ricci flat cone, and type IIB throat from stored tensors.

Run with .venv.noindex/bin/python _tools/test_klebanov_witten.py.
"""
import json
import sys
import unittest
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "_tools" / "derivations"))
import verify_metrics as vm


def matrix(reader, components, names):
    out = sp.zeros(len(names))
    for component in components:
        a, b = (names.index(i) for i in component["indices"])
        out[a, b] = reader(component["value"])
    return out


class Curvature(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.metric = json.loads((ROOT / "MFS/assets/data/metrics/klebanov_witten.json").read_text())
        cls.chart = cls.metric["coordinates"][0]
        cls.names = cls.chart["coords"]
        cls.reader = vm.Reader(cls.names, [p["symbol"] for p in cls.chart["parameters"]], ())
        cls.L = cls.reader.parameters["L"]
        cls.g = matrix(cls.reader, cls.chart["metric_components"], cls.names)
        cls.angles = [cls.reader.symbol[n] for n in cls.names[5:]]
        cls.h = (cls.g[5:, 5:]/cls.L**2).applyfunc(vm.norm)
        cls.base = vm.Geometry(cls.h, cls.angles, 60)
        cls.r = cls.reader.symbol["r"]
        cls.cone = vm.Geometry(sp.diag(1, cls.r**2*cls.h), [cls.r]+cls.angles, 60)

    def near_zero(self, value, at):
        self.assertLess(abs(complex(sp.N(value.subs(at), 40))), 1e-25)

    def angular_points(self):
        for a, b in ((sp.Rational(2, 5), sp.Rational(7, 8)),
                     (sp.Rational(7, 6), sp.Rational(11, 5)),
                     (sp.Rational(13, 7), sp.Rational(4, 9))):
            yield dict(zip(self.angles, (a, sp.Rational(1, 3), b, sp.Rational(2, 3), sp.Rational(3, 4))))

    def test_base_is_einstein_with_scalar_twenty(self):
        self.assertEqual(self.base.n, 5)
        ricci = sp.Matrix(self.base.ricci_ll())
        self.assertEqual((ricci-4*self.h).applyfunc(vm.norm), sp.zeros(5))
        self.assertEqual(self.base.ricci_scalar(), 20)

    def test_base_squared_riemann_is_136_at_three_polar_pairs(self):
        square = self.base.kretschmann()
        for at in self.angular_points():
            self.near_zero(square-136, at)

    def test_cone_is_ricci_flat_and_singular_at_its_tip(self):
        self.assertEqual(self.cone.n, 6)
        self.assertEqual(sp.Matrix(self.cone.ricci_ll()).applyfunc(vm.norm), sp.zeros(6))
        self.assertEqual(self.cone.ricci_scalar(), 0)
        square = self.cone.kretschmann()
        for n, at in enumerate(self.angular_points()):
            at[self.r] = sp.Rational(n+2, 3)
            self.near_zero(square-96/self.r**4, at)
        self.assertEqual(vm.norm(square*self.r**4), 96)
        self.assertEqual(sp.limit(square, self.r, 0, dir="+"), sp.oo)

    def test_throat_has_two_opposite_ricci_sectors_and_zero_scalar(self):
        self.assertEqual(len(self.names), 10)
        ricci = matrix(self.reader, self.chart["ricci_tensor"]["variants"]["ul"]["nonzero"], self.names)
        expected = sp.diag(*([-4/self.L**2]*5+[4/self.L**2]*5))
        self.assertEqual((ricci-expected).applyfunc(vm.norm), sp.zeros(10))
        scalar = self.reader(vm.scalar_parts(self.chart["ricci_scalar"])[0])
        self.assertEqual(scalar, 0)
        self.assertEqual(vm.norm(sp.trace(ricci)), scalar)
        self.assertNotEqual(ricci, sp.zeros(10))

    def test_throat_squared_riemann_is_40_plus_136(self):
        ads = vm.Geometry(self.g[:5, :5], [self.reader.symbol[n] for n in self.names[:5]], 60)
        self.assertEqual(vm.norm(ads.ricci_scalar()+20/self.L**2), 0)
        self.assertEqual(vm.norm(ads.kretschmann()-40/self.L**4), 0)
        square = self.reader(vm.scalar_parts(self.chart["kretschmann"])[0])
        self.assertEqual(vm.norm(square-176/self.L**4), 0)
        for at in self.angular_points():
            at[self.L] = sp.Rational(7, 5)
            self.near_zero(square-ads.kretschmann()-self.base.kretschmann()/self.L**4, at)
        for component in self.chart["riemann"]["variants"]["ulll"]["nonzero"]:
            sectors = {self.names.index(n)<5 for n in component["indices"]}
            self.assertEqual(len(sectors), 1, component["indices"])

    def test_the_embedded_fibre_closes_with_the_source_period(self):
        circumference = 4*sp.pi*sp.sqrt(self.h[4, 4])
        self.assertEqual(circumference, 2*sp.pi*sp.Rational(2, 3))
        self.assertEqual(vm.norm(4*self.h[4, 4]), sp.Rational(4, 9))

    def test_the_spatial_fibre_has_unit_chern_number_on_each_sphere(self):
        theta, phi, alpha, beta, psi = self.angles
        self.assertIn(r"\psi \in [0,4\pi)", self.chart["domains"])
        # psi/2 has period 2 pi, so the normalized circle connection is half the angular one.
        for polar, azimuth, index in ((theta, phi, 1), (alpha, beta, 3)):
            potential = vm.norm(self.h[4, index]/(2*self.h[4, 4]))
            self.assertEqual(vm.norm(potential-sp.cos(polar)/2), 0)
            chern = sp.integrate(sp.diff(potential, polar), (polar, 0, sp.pi))
            self.assertEqual(chern, -1)


if __name__ == "__main__":
    unittest.main()
