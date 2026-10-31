"""Limits and frame interpretation of the published false vacuum background.

The chart printer owns the field equations and coordinate-map checks.
These checks pin the reference data in the de Sitter limit, the distinct
Einstein metric, the growing nucleation ratio and the null singular edge.
"""
import json
import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "_tools" / "derivations"))

import sympy as sp
import verify_metrics as vm


class Physics(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        metric = json.loads((ROOT / "MFS/assets/data/metrics/extended_inflation.json").read_text())
        cls.entry = next(c for c in metric["coordinates"] if c["id"] == "cosmic")
        cls.reader = vm.Reader(cls.entry["coords"], [p["symbol"] for p in cls.entry["parameters"]], ())
        cls.t = cls.reader.symbol["t"]
        cls.p, cls.t0 = (cls.reader.parameters[n] for n in ("p", "t_0"))
        rr = next(c["value"] for c in cls.entry["metric_components"] if c["indices"] == ["r", "r"])
        cls.a2 = cls.reader(rr)

    def scale(self, p, t0, time):
        return math.sqrt(float(self.a2.subs({self.p: p, self.t0: t0, self.t: time,
                                            self.reader.c: 1}).evalf(35)))

    def test_de_sitter_limit_keeps_reference_coupling_density_and_elapsed_time(self):
        # G_0=c=1 and rho_v=3/(8 pi) set the limiting Hubble rate to one.
        for elapsed in (-0.3, 0.5, 1.0):
            errors = []
            for omega in (1e3, 1e5, 1e7):
                reference = math.sqrt((6*omega+5)*(2*omega+3)/12)
                a = self.scale(omega+0.5, reference, reference+elapsed)
                Phi = (1+elapsed/reference)**2
                errors.append(abs(a/math.exp(elapsed)-1)+abs(Phi-1))
            self.assertLess(errors[1], errors[0]/50)
            self.assertLess(errors[2], errors[1]/50)
            self.assertLess(errors[2], 1e-6)

    def test_einstein_metric_has_its_own_cosmic_power(self):
        # g_E=Phi g_J, dt_E=sqrt(Phi)dt and t_E=t^2/(2t_0).
        for p in (1.1, 2.0, 17.0):
            first, second, t0 = 0.7, 1.4, 1.2
            aE1 = first/t0*self.scale(p, t0, first)
            aE2 = second/t0*self.scale(p, t0, second)
            tE_ratio = (second/first)**2
            self.assertAlmostEqual(math.log(aE2/aE1)/math.log(tE_ratio), (p+1)/2, places=11)
            # The Jordan and Einstein metrics differ at every Phi != 1.
            self.assertNotAlmostEqual(aE1, self.scale(p, t0, first))

    def test_constant_nucleation_ratio_grows_as_the_fourth_power_of_time(self):
        connection = self.entry["christoffel"]["variants"]["ull"]["nonzero"]
        H = self.reader.c*self.reader(next(c["value"] for c in connection
                                          if c["indices"] == ["r", "t", "r"]))
        at = {self.p: 2.3, self.t0: 1.7, self.reader.c: 1}
        rate1 = float(H.subs(at).subs(self.t, 0.8))
        rate2 = float(H.subs(at).subs(self.t, 2.4))
        self.assertAlmostEqual((rate1/rate2)**4, 3**4, places=10)

    def test_past_null_edge_is_curvature_singular(self):
        # Hold the advanced conformal coordinate eta+r=1 at p=2, t_0=c=1.
        K = self.reader(self.entry["kretschmann"].partition("=")[2])
        at = {self.p: 2, self.t0: 1, self.reader.c: 1}
        last_K, last_gap = 0, math.inf
        for time in (0.01, 0.001, 0.0001):
            eta, radius = -1/time, 1+1/time
            u, v = math.atan(eta-radius), math.atan(eta+radius)
            T, X = u+v, v-u
            gap = abs(T-X+math.pi)
            value = float(K.subs(at).subs(self.t, time))
            self.assertLess(gap, last_gap)
            self.assertGreater(value, last_K)
            last_K, last_gap = value, gap
        self.assertLess(last_gap, 0.0002)
        self.assertGreater(last_K, 1e18)


if __name__ == "__main__":
    unittest.main()
