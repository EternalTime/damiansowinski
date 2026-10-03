"""The quantum Oppenheimer-Snyder black hole: python3 -m unittest discover -s _tools

What its texts state, held on the published files. Outside the dust, the static chart is
Lewandowski, Ma, Yang and Zhang's f = 1 - 2m/r + alpha m^2/r^4, whose mixed Einstein tensor is the
energy density 3 alpha m^2/r^6 that falls off as r^-6, whose Ricci and Kretschmann scalars are Kelly,
Santacruz and Wilson-Ewing's, and whose horizons merge at r = 3m/2 = 2 sqrt(alpha/3) at the least
mass 4 sqrt(alpha)/3 sqrt(3). Inside, the published Einstein tensor of the comoving chart with the
bouncing scale factor a^3 = 1 + 9 c^2 tau^2/alpha is the Friedmann equation of loop quantum cosmology,
and the surface chi_0 = r_b is the exterior's radial geodesic that falls from rest far away. The path
of the surface every diagram draws, quantum_os.py, is held to the same geodesic in each chart and to
the numbers the captions state. The tests need sympy and numpy and are skipped where either is absent.
"""
import importlib.util
import json
import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "MFS" / "assets" / "data"
READY = all(importlib.util.find_spec(name) is not None for name in ("sympy", "numpy", "scipy"))
sys.path.insert(0, str(Path(__file__).resolve().parent / "derivations"))


def metric():
    return json.loads((DATA / "metrics" / "quantum_oppenheimer_snyder.json").read_text(encoding="utf-8"))


def chart_of(system):
    return next(c for c in metric()["coordinates"] if c["id"] == system)


@unittest.skipUnless(READY, "sympy, numpy and scipy are not installed")
class Exterior(unittest.TestCase):
    def setUp(self):
        import sympy
        import verify_metrics
        self.sp, self.vm = sympy, verify_metrics
        self.r, self.m, self.alpha = sympy.symbols("r m alpha", positive=True)
        self.f = 1 - 2 * self.m / self.r + self.alpha * self.m ** 2 / self.r ** 4

    def read(self, system, field, variant=None):
        entry = chart_of(system)
        reader = self.vm.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]], ())
        mine = {reader.symbol["r"]: self.r, reader.parameters["m"]: self.m,
                reader.parameters["alpha"]: self.alpha, reader.c: 1}
        values = entry[field] if variant is None else entry[field]["variants"][variant]["nonzero"]
        if isinstance(values, str):
            return self.sp.sympify(reader(values.partition("=")[2])).subs(mine)
        return {tuple(e["indices"]): self.sp.sympify(reader(e["value"])).subs(mine) for e in values}

    def test_the_static_chart_is_the_deformed_schwarzschild_metric(self):
        g = self.read("static", "metric_components")
        self.assertEqual(self.sp.simplify(g[("t", "t")] + self.f), 0)
        self.assertEqual(self.sp.simplify(g[("r", "r")] - 1 / self.f), 0)

    def test_the_vacuum_carries_an_energy_density_falling_as_r_to_the_minus_six(self):
        for system in ("static", "painleve_gullstrand", "eddington_finkelstein_ingoing", "eddington_finkelstein_outgoing"):
            G = self.read(system, "einstein_tensor", "ul")
            density = -G[(chart_of(system)["coords"][0],) * 2]
            self.assertEqual(self.sp.simplify(density - 3 * self.alpha * self.m ** 2 / self.r ** 6), 0, system)
            self.assertEqual(self.sp.simplify(G[("\\theta", "\\theta")] - 6 * self.alpha * self.m ** 2 / self.r ** 6), 0)

    def test_the_scalars_are_kelly_santacruz_and_wilson_ewings(self):
        # Their R and K with R_S = 2m and gamma^2 Delta = alpha/4.
        R_S, gD = 2 * self.m, self.alpha / 4
        R = -6 * gD * R_S ** 2 / self.r ** 6
        K = 12 * R_S ** 2 / self.r ** 6 * (1 - 10 * gD * R_S / self.r ** 3 + 39 * gD ** 2 * R_S ** 2 / self.r ** 6)
        for system in ("static", "painleve_gullstrand"):
            self.assertEqual(self.sp.simplify(self.read(system, "ricci_scalar") - R), 0, system)
            self.assertEqual(self.sp.simplify(self.read(system, "kretschmann") - K), 0, system)

    def test_the_painleve_gullstrand_slices_fall_at_the_speed_of_a_fall_from_rest(self):
        g = self.read("painleve_gullstrand", "metric_components")
        self.assertEqual(self.sp.simplify(self.sp.expand(g[("\\tau", "r")] ** 2 - (1 - self.f))), 0)
        self.assertEqual(g[("r", "r")], 1)

    def test_the_least_mass_and_the_merged_horizon(self):
        least = {self.alpha: 27 * self.m ** 2 / 16}
        merged = 3 * self.m / 2
        self.assertEqual(self.sp.simplify(self.f.subs(least).subs(self.r, merged)), 0)
        self.assertEqual(self.sp.simplify(self.sp.diff(self.f, self.r).subs(least).subs(self.r, merged)), 0)
        a = self.sp.Symbol("a", positive=True)
        self.assertEqual(self.sp.simplify(4 * self.sp.sqrt(a) / (3 * self.sp.sqrt(3)) * 3 / 2 - 2 * self.sp.sqrt(a / 3)), 0)
        # Below the least mass f has no zero: at alpha = 2 m^2 its least value, at r^3 = 2 alpha m, is positive.
        F = self.f.subs({self.alpha: 2, self.m: 1})
        self.assertGreater(float(F.subs(self.r, 4 ** (1 / 3))), 0)

    def test_the_surface_turns_round_where_f_is_one(self):
        turn = (self.alpha * self.m / 2) ** self.sp.Rational(1, 3)
        self.assertEqual(self.sp.simplify(self.f.subs(self.r, turn) - 1), 0)


@unittest.skipUnless(READY, "sympy, numpy and scipy are not installed")
class Interior(unittest.TestCase):
    def test_the_bouncing_dust_obeys_the_friedmann_equation_of_loop_quantum_cosmology(self):
        import sympy as sp
        import verify_metrics as vm
        entry = chart_of("interior_comoving")
        reader = vm.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]], ())
        tau, a = reader.symbol["\\tau"], reader.parameters["a"]
        x, alpha = sp.symbols("x alpha", positive=True)
        scale = (1 + 9 * x ** 2 / alpha) ** sp.Rational(1, 3)
        bounce = {sp.Derivative(a, (tau, 2)): sp.diff(scale, x, 2), sp.Derivative(a, tau): sp.diff(scale, x), a: scale}
        G = {tuple(e["indices"]): reader(e["value"]) for e in entry["einstein_tensor"]["variants"]["ul"]["nonzero"]}
        density = 12 / (alpha * scale ** 3)              # 8 pi G rho/c^2, critical at a = 1
        self.assertEqual(sp.simplify(G[("\\tau", "\\tau")].subs(bounce) + density * (1 - density * alpha / 12)), 0)
        # Pressureless only where Einstein's equation would hold: the bounce needs G^chi_chi != 0.
        self.assertNotEqual(sp.simplify(G[("\\chi", "\\chi")].subs(bounce)), 0)

    def test_the_surface_is_the_geodesic_that_falls_from_rest_far_away(self):
        import quantum_os as qos
        import numpy as np
        tau = np.linspace(-4, 0.4, 23)
        R, h = qos.radius(tau), 1e-6
        speed = (qos.radius(tau + h) - qos.radius(tau - h)) / (2 * h)
        self.assertLess(float(np.max(np.abs(speed ** 2 - (1 - qos.f(R))))), 1e-8)
        self.assertAlmostEqual(float(qos.radius(0.0)), qos.RB, places=14)
        self.assertAlmostEqual(float(qos.radius(qos.TAU_M)), qos.RM, places=12)
        self.assertAlmostEqual(float(qos.scale_factor(1.0)) * qos.RB, float(qos.radius(1.0)), places=14)


@unittest.skipUnless(READY, "sympy, numpy and scipy are not installed")
class Drawn(unittest.TestCase):
    def setUp(self):
        import quantum_os
        import numpy
        self.q, self.np = quantum_os, numpy

    def test_the_numbers_the_captions_state(self):
        q = self.q
        self.assertAlmostEqual(q.RB, 0.855, places=3)
        self.assertAlmostEqual(q.RM, 1.127, places=3)
        self.assertAlmostEqual(q.RP, 1.777, places=3)
        self.assertAlmostEqual(q.kappa(q.RM) / q.kappa(q.RP), 3.34, places=2)
        self.assertLess(abs(float(q.f(q.RM))) + abs(float(q.f(q.RP))), 1e-12)
        self.assertGreater(1 / q.M_LEAST, 1)        # the hole drawn has more than the least mass

    def test_the_tortoise_coordinate(self):
        r, h = self.np.array([0.9, 1.0, 1.4, 2.5, 9.0]), 1e-6
        slope = (self.q.tortoise(r + h) - self.q.tortoise(r - h)) / (2 * h)
        self.assertLess(float(self.np.max(self.np.abs(slope * self.q.f(r) - 1))), 1e-6)
        self.assertEqual(float(self.q.tortoise(self.q.RB)), 0.0)

    def test_the_surface_in_each_chart(self):
        from scipy.integrate import quad
        q, np = self.q, self.np
        # The way in: dv/dtau = 1/(1 + N); the way out: dv/dtau = (1 + N)/f and dT/dtau = (1 + N^2)/f.
        v_in = -quad(lambda t: 1 / (1 + float(q.speed(q.radius(t)))), -2.0, 0)[0]
        self.assertAlmostEqual(float(q.advanced(-2.0)), v_in, places=7)
        v_out = quad(lambda t: (1 + float(q.speed(q.radius(t)))) / float(q.f(q.radius(t))), 0, 0.3)[0]
        self.assertAlmostEqual(float(q.advanced(0.3)), v_out, places=7)
        T_out = quad(lambda t: (1 + float(q.speed(q.radius(t))) ** 2) / float(q.f(q.radius(t))), 0, 0.3)[0]
        self.assertAlmostEqual(float(q.slice_time(0.3)), T_out, places=7)
        # The radius read back at those times, the way out the time reverse of the way in.
        self.assertAlmostEqual(float(q.radius_at_advanced(q.advanced(0.3))), float(q.radius(0.3)), places=6)
        self.assertAlmostEqual(float(q.radius_at_slice_time(-1.5)), float(q.radius(-1.5)), places=12)
        self.assertAlmostEqual(float(q.radius_at_retarded(0.7)), float(q.radius_at_advanced(-0.7)), places=12)
        # The way out reaches r_- only as v runs to infinity.
        self.assertGreater(float(q.advanced(q.TAU_M * (1 - 1e-6))), 9)

    def test_conformal_time(self):
        q, np = self.q, self.np
        tau, h = np.array([-5.0, -0.3, 0.0, 0.7, 4.0]), 1e-6
        rate = (q.conformal_time(tau + h) - q.conformal_time(tau - h)) / (2 * h)
        self.assertLess(float(np.max(np.abs(rate * q.scale_factor(tau) - 1))), 1e-7)

    def test_the_collapse_view_turns_round_on_r_b(self):
        views = json.loads((DATA / "conformal" / "quantum_oppenheimer_snyder.json").read_text(encoding="utf-8"))["views"]
        collapse = next(v for v in views if v["id"] == "collapse")
        surface = [p for layer in collapse["layers"] if layer["class"] == "surface" for p in layer["points"]]
        leftmost = min(surface, key=lambda p: p[0])
        self.assertAlmostEqual(leftmost[0], -math.pi / 2, places=3)
        self.assertAlmostEqual(leftmost[1], math.pi, places=2)
        self.assertEqual({v["id"] for v in views}, {"tower", "ingoing", "outgoing", "collapse"})

    def test_the_painleve_gullstrand_view_draws_the_surface_down_to_r_b(self):
        systems = json.loads((DATA / "diagrams" / "quantum_oppenheimer_snyder.json").read_text(encoding="utf-8"))["systems"]
        view = systems["painleve_gullstrand"][0]
        shell = [p for marker in view["markers"] if marker["kind"] == "shell" for line in marker["lines"] for p in line]
        x0, x1 = view["box"][0], view["box"][1]
        least = min(x0 + p[0] * (x1 - x0) for p in shell)
        self.assertAlmostEqual(least, self.q.RB, places=2)


if __name__ == "__main__":
    unittest.main()
