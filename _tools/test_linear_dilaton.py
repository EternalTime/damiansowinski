"""Frame transformation and curvature of the linear dilaton cosmology, and the embedding
diagram's ring of free particles and its world tube held to the published metric, from the
numbers in its file."""

import json
import math
from pathlib import Path
import sys
import unittest

import sympy as sp

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "_tools" / "derivations"))
import verify_metrics as vm


class LinearDilaton(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        page = json.loads((ROOT / "MFS/assets/data/metrics/linear_dilaton.json").read_text())
        cls.charts = {chart["id"]: chart for chart in page["coordinates"]}

    @staticmethod
    def read_chart(chart):
        reader = vm.Reader(chart["coords"], [p["symbol"] for p in chart["parameters"]], ())
        metric = vm.metric_from_line_element(reader, chart["line_element"], chart["coords"])
        return reader, sp.Matrix(metric)

    def test_ricci_scalar(self):
        chart = self.charts["cosmic"]
        reader, _ = self.read_chart(chart)
        published = reader(chart["ricci_scalar"].split("=", 1)[1])
        self.assertEqual(sp.simplify(published - 6 / reader.symbol[r"\tau"]**2), 0)

    def test_string_frame_factor_and_cosmic_time(self):
        cosmic_reader, cosmic = self.read_chart(self.charts["cosmic"])
        conformal_reader, conformal = self.read_chart(self.charts["conformal"])
        eta = conformal_reader.symbol[r"\eta"]
        q = conformal_reader.parameters["q"]
        tau = cosmic_reader.symbol[r"\tau"]
        cosmic = cosmic.subs({cosmic_reader.parameters["q"]: q, tau: sp.exp(q*eta)/q})
        jacobian = sp.diag(sp.exp(q*eta), 1, 1, 1)
        pulled = jacobian.T * cosmic * jacobian
        string_frame = sp.diag(-1, 1, 1, 1)
        expected = sp.exp(2*q*eta) * string_frame
        for entry in list(conformal - expected) + list(pulled - expected):
            self.assertEqual(sp.simplify(entry), 0)

    def test_embedding_is_a_ring_of_free_particles_of_proper_radius_tau(self):
        # Neither chart publishes a Gamma^i_00 of its own time, so a free particle at rest in the chart stays
        # at rest, and the ring at rest on x^2 + y^2 = 1/q^2 has the proper radius q tau / q = tau.
        for system, time in (("cosmic", r"\tau"), ("conformal", r"\eta")):
            published = self.charts[system]["christoffel"]["variants"]["ull"]["nonzero"]
            self.assertFalse([c for c in published if c["indices"][1:] == [time, time] and c["indices"][0] != time])
        reader, metric = self.read_chart(self.charts["cosmic"])
        tau, q = reader.symbol[r"\tau"], reader.parameters["q"]
        scale = sp.lambdify(tau, sp.sqrt(metric[1, 1].subs(q, 1)), "math")
        embedding = json.loads((ROOT / "MFS/assets/data/embedding/linear_dilaton.json").read_text())
        views = {view["id"]: view for view in embedding["views"]}
        self.assertEqual(list(views), ["tube", "ring"])
        for view in views.values():
            self.assertIn("$q = 1$", view["settings"])
            self.assertEqual(view["unit"], "$1/q$")
        ring = views["ring"]
        self.assertEqual([s["time"] for s in ring["surfaces"]], [0.5, 1, 1.5, 2])
        frames = ring["movie"]["frames"]
        values = [frame["value"] for frame in frames]
        self.assertEqual((len(frames), values[0], values[-1]), (61, 0.5, 2))
        radii = []
        for frame in frames:
            t = frame["value"]
            points = frame["curves"][0]["points"]
            for X, Y, Z in points:
                self.assertAlmostEqual(math.hypot(X, Y), scale(t), delta=1e-6, msg=t)
                self.assertAlmostEqual(math.hypot(X, Y), t, delta=1e-6, msg=t)
                self.assertEqual(Z, 0, t)
            self.assertEqual(len(frame["dots"]), 12, t)
            radii.append(points[0][0])
        # The frames stand a steady step of q tau apart, and the ring grows by the same step at each.
        for earlier, later, a, b in zip(values, values[1:], radii, radii[1:]):
            self.assertAlmostEqual(later - earlier, 0.025, delta=1e-6)
            self.assertAlmostEqual(b - a, later - earlier, delta=1e-6)
        # Each row of the world tube is the ring at its tau, at the height 2 tau: on the cone rho = z/2,
        # whose apex is the moment tau = 0 at the foot of the axis, where the curvature diverges.
        tube = views["tube"]["surfaces"][0]
        grid = tube["pieces"][0]["grid"]
        self.assertEqual((grid["u"][0], grid["u"][-1]), (0.5, 2))
        for u, a, b, z in zip(grid["u"], grid["a"], grid["b"], grid["z"]):
            self.assertAlmostEqual(a, u, delta=1e-6, msg=u)
            self.assertAlmostEqual(b, u, delta=1e-6, msg=u)
            self.assertAlmostEqual(z, 2 * u, delta=1e-6, msg=u)
        self.assertEqual([c["time"] for c in tube["curves"]], [0.5, 1, 1.5, 2])
        self.assertGreater(tube["axis"]["from"], 0)


if __name__ == "__main__":
    unittest.main()
