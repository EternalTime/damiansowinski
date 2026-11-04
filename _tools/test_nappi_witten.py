"""Expansion, recollapse, and curvature of the closed string cosmology."""

import json
from pathlib import Path
import sys
import unittest

import mpmath as mp
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

    def test_total_spatial_volume_expands_and_then_recollapses(self):
        # Integrate the stored spatial metric, including both angular periods.
        # The diagonal product is used directly, with no symbolic determinant.
        density = sp.sqrt(self.metric[1, 1]*self.metric[2, 2]*self.metric[3, 3])
        f = sp.lambdify([self.s, self.psi, self.b, self.L], density, "mpmath")
        with mp.workdps(40):
            for deformation, length in ((mp.mpf("0.5"), mp.mpf("0.8")),
                                         (mp.mpf(1), mp.mpf(1)),
                                         (mp.mpf(2), mp.mpf("1.25"))):
                volumes = []
                maximum = 2*mp.pi**2*length**3
                for y in (-4, -2, -1, 0, 1, 2, 4):
                    time = mp.atan(mp.exp(y)/deformation)
                    volume = 4*mp.pi**2*mp.quad(lambda s: f(s, time, deformation, length),
                                                [0, mp.pi/4, mp.pi/2])
                    expected = maximum*(mp.mpf(y)/mp.sinh(y) if y else 1)
                    self.assertLess(abs(volume-expected)/maximum, mp.mpf("1e-30"))
                    volumes.append(volume)
                self.assertTrue(all(a < b for a, b in zip(volumes[:3], volumes[1:4])))
                self.assertTrue(all(a > b for a, b in zip(volumes[3:6], volumes[4:])))
                for y in (-20, 20):
                    time = mp.atan(mp.exp(y)/deformation)
                    volume = 4*mp.pi**2*mp.quad(lambda s: f(s, time, deformation, length),
                                                [0, mp.pi/4, mp.pi/2])
                    self.assertLess(volume/maximum, mp.mpf("1e-6"))

    def test_local_expansion_is_not_uniform_at_the_volume_maximum(self):
        theta = sum(sp.diff(self.metric[i, i], self.psi)/(2*self.metric[i, i])
                    for i in (1, 2, 3))/self.L
        theta = theta.subs({self.psi: sp.pi/4, self.b: 1})
        self.assertEqual(vm.norm(theta+2*sp.cos(2*self.s)/self.L), 0)
        self.assertLess(float(theta.subs({self.s: sp.pi/8, self.L: 1})), 0)
        self.assertGreater(float(theta.subs({self.s: 3*sp.pi/8, self.L: 1})), 0)

    def test_the_later_wave_is_the_equal_negative_eigenvalue_member(self):
        page = json.loads((ROOT/"MFS/assets/data/metrics/cahen_wallach.json").read_text())
        chart = next(c for c in page["coordinates"] if c["id"] == "rotating")
        reader = vm.Reader(chart["coords"], [p["symbol"] for p in chart["parameters"]], ())
        metric = sp.Matrix(vm.metric_from_line_element(reader, chart["line_element"], chart["coords"]))
        u, v, a1, a2, b0 = sp.symbols("u v a1 a2 b0", real=True)
        transform = sp.Matrix([u, -v-b0*u/2, a1, a2])
        replacements = dict(zip((reader.symbol[c] for c in chart["coords"]), transform))
        replacements.update({reader.c: 1, reader.parameters["omega"]: -sp.Rational(1, 2)})
        jacobian = transform.jacobian([u, v, a1, a2])
        pulled = jacobian.T*metric.subs(replacements, simultaneous=True)*jacobian
        # Nappi and Witten 1993 (17), independently of the cosmological metric.
        source = sp.Matrix([[b0, 1, a2/2, -a1/2], [1, 0, 0, 0],
                            [a2/2, 0, 1, 0], [-a1/2, 0, 0, 1]])
        self.assertTrue(all(sp.expand(value) == 0 for value in pulled-source))


if __name__ == "__main__":
    unittest.main()
