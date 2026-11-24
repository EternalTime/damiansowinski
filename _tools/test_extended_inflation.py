"""Limits and frame interpretation of the published false vacuum background.

The chart printer owns the field equations and coordinate-map checks.
These checks pin the reference data in the de Sitter limit, the distinct
Einstein metric, the growing nucleation ratio and the null singular edge.
The embedding diagram's ring of free particles, the central observer's event
horizon and the ring's world tube are held to the metric at p = 2, from the
numbers in its file.
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

    def test_embedding_is_a_ring_of_free_particles_outgrowing_the_central_event_horizon(self):
        # At p = 2 and t_0 = c = 1 a ring at rest on r = c t_0 has the proper radius a r = t^2, and the
        # event horizon of the observer at the centre, of comoving radius c int_t^oo dt'/a = 1/t, the
        # proper radius t: the ring is inside it before t = t_0 and outside it after.
        embedding = json.loads((ROOT / "MFS/assets/data/embedding/extended_inflation.json").read_text())
        views = {view["id"]: view for view in embedding["views"]}
        self.assertEqual(list(views), ["tube", "ring"])
        ring = views["ring"]
        self.assertEqual([s["time"] for s in ring["surfaces"]], [0.5, 1, 1.5])
        frames = ring["movie"]["frames"]
        values = [frame["value"] for frame in frames]
        self.assertEqual((values[0], values[-1], len(values)), (0.5, 1.5, 41))
        for earlier, later in zip(values, values[1:]):
            self.assertAlmostEqual(later - earlier, 0.025, delta=1e-6)
        radii = []
        for frame in frames:
            t = frame["value"]
            curves = {curve["class"]: curve for curve in frame["curves"]}
            for name, want in (("particles", t * t), ("horizon", t)):
                for X, Y, Z in curves[name]["points"]:
                    self.assertAlmostEqual(math.hypot(X, Y), want, delta=1e-6, msg=(name, t))
                    self.assertEqual(Z, 0, (name, t))
            self.assertEqual(len(frame["dots"]), 12, t)
            radii.append(tuple(math.hypot(*curves[name]["points"][0][:2]) for name in ("particles", "horizon")))
            plane = frame["pieces"][0]
            for r, rho, z in plane["points"]:
                self.assertAlmostEqual(rho, t * t * r, delta=1e-6, msg=t)
                self.assertEqual(z, 0, t)
            self.assertGreater(plane["points"][-1][1], t * t, t)
        inside = [ring < horizon for ring, horizon in radii]
        self.assertEqual(inside, [t < 1 for t in values[:20]] + [False] * 21)
        self.assertAlmostEqual(radii[20][0], radii[20][1], delta=1e-12)
        grown = [later[0] - earlier[0] for earlier, later in zip(radii, radii[1:])]
        self.assertTrue(all(later > earlier > 0 for earlier, later in zip(grown, grown[1:])))
        tube = views["tube"]["surfaces"][0]
        grid = tube["pieces"][0]["grid"]
        for u, a, b, z in zip(grid["u"], grid["a"], grid["b"], grid["z"]):
            self.assertAlmostEqual(a, u * u, delta=1e-6, msg=u)
            self.assertAlmostEqual(b, u * u, delta=1e-6, msg=u)
            self.assertAlmostEqual(z, 3 * u, delta=1e-6, msg=u)
        self.assertEqual([c["time"] for c in tube["curves"]], [0.5, 1, 1.5])


if __name__ == "__main__":
    unittest.main()
