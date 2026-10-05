"""Hiscock's evaporating black hole, ingoing and outgoing Vaidya metrics joined on a timelike surface:
.venv.noindex/bin/python -m unittest discover -s _tools

What its texts state about the surface of pair creation, held on the published files. The surface
r = R is in neither chart's mathematics, so what it carries is taken here from the jump of the
extrinsic curvature between the published ingoing and outgoing metrics, with the published
Christoffel symbols, by Israel's junction conditions: [[K^theta_theta]] = -4 pi sigma and
[[K^tau_tau + K^theta_theta]] = 8 pi p. With the mass the same on both sides the surface has no
surface energy density and a surface pressure, a negative tension, as Hayward found for the surface
at a fixed radius, Phys. Rev. Lett. 96, 031103 (2006). The model the diagrams draw is held to the
same conditions and to Hiscock's, dM/du -> 0 as M -> 0. The tests need sympy, and the drawn model
numpy and scipy as well; each is skipped where its library is absent.
"""
import importlib.util
import json
import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "MFS" / "assets" / "data"
HAS_SYMPY = importlib.util.find_spec("sympy") is not None
HAS_DRAWING = all(importlib.util.find_spec(name) is not None for name in ("sympy", "numpy", "scipy", "contourpy"))
sys.path.insert(0, str(Path(__file__).resolve().parent / "derivations"))


@unittest.skipUnless(HAS_SYMPY, "sympy is not installed")
class Surface(unittest.TestCase):
    def setUp(self):
        import sympy
        import verify_metrics
        self.sp, self.vm = sympy, verify_metrics

    def curvatures(self, system):
        """(x^0, m, R, K_tau_tau, K^theta_theta) of the surface r = R(x^0) in one published chart,
        with the unit normal toward larger r: K_tau_tau = -n_mu a^mu for the surface's own unit
        tangent, and K^theta_theta = -Gamma^mu_theta_theta n_mu / R^2."""
        sp = self.sp
        metric = json.loads((DATA / "metrics" / "hiscock.json").read_text(encoding="utf-8"))
        entry = next(c for c in metric["coordinates"] if c["id"] == system)
        reader = self.vm.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]], ())
        g = {tuple(e["indices"]): reader(e["value"]) for e in entry["metric_components"]}
        gamma = {tuple(e["indices"]): reader(e["value"]) for e in entry["christoffel"]["variants"]["ull"]["nonzero"]}
        null = entry["coords"][0]
        x0, r = reader.symbol[null], reader.symbol["r"]
        m, R = reader.parameters["m"], reader.parameters["R"]
        names = [null, "r"]
        tangent = [sp.Integer(1), sp.diff(R, x0)]
        size = sp.sqrt(-(g[(null, null)] + 2 * g[(null, "r")] * tangent[1]))
        unit = [t / size for t in tangent]
        normal = [-tangent[1] / size, 1 / size]
        on = {r: R}
        accel = [(sp.diff(unit[k].subs(on), x0) / size.subs(on)
                  + sum(gamma.get((names[k], names[a], names[b]), 0) * unit[a] * unit[b]
                        for a in range(2) for b in range(2))).subs(on).doit() for k in range(2)]
        k_tt = -(normal[0] * accel[0] + normal[1] * accel[1]).subs(on)
        k_theta = -sum(gamma.get((names[k], "\\theta", "\\theta"), 0) * normal[k] for k in range(2)).subs(on) / R ** 2
        return x0, m, R, k_tt, k_theta

    def jumps(self, N, R, v):
        """[[K^tau_tau]] and [[K^theta_theta]], outside minus inside, on the surface r = R(v) with
        the mass N(v) inside and the same mass outside, as expressions in v. The outgoing chart's
        retarded time has du/dv = 1 - 2R'/(1 - 2N/R) on the surface, which makes the induced
        metric one."""
        sp = self.sp
        vi, mi, Ri, ktt_in, kth_in = self.curvatures("ingoing")
        uo, mo, Ro, ktt_out, kth_out = self.curvatures("outgoing")
        inside = lambda e: e.subs({mi: N, Ri: R}).doit().subs(vi, v)
        stretch = 1 - 2 * sp.diff(R, v) / (1 - 2 * N / R)
        d = lambda f: sp.diff(f, v) / stretch           # a derivative along u, on the surface
        outside = lambda e: (e.subs(sp.Derivative(Ro, (uo, 2)), d(d(R))).subs(sp.Derivative(Ro, uo), d(R))
                             .subs(sp.Derivative(mo, uo), d(N)).subs({Ro: R, mo: N}))
        # u.u = -1, so K^tau_tau = -K_tau_tau.
        return -(outside(ktt_out) - inside(ktt_in)), outside(kth_out) - inside(kth_in), stretch

    def test_the_induced_metric_is_one_from_both_sides(self):
        sp = self.sp
        v = sp.Symbol("v", real=True)
        N, R = sp.Function("N")(v), sp.Function("R")(v)
        F = 1 - 2 * N / R
        stretch = 1 - 2 * sp.diff(R, v) / F
        ingoing = -F + 2 * sp.diff(R, v)                                    # ds^2/dv^2 on the surface
        outgoing = (-F - 2 * sp.diff(R, v) / stretch) * stretch ** 2        # the same from the outgoing chart
        self.assertEqual(sp.simplify(ingoing - outgoing), 0)

    def test_a_surface_at_a_fixed_radius_has_no_density_and_haywards_tension(self):
        # p = -N'/(4 pi r_0 F^(3/2)) with F = 1 - 2N/r_0: a pressure while the mass falls.
        sp = self.sp
        v, r0 = sp.Symbol("v", real=True), sp.Symbol("r_0", positive=True)
        N = sp.Function("N")(v)
        k_tau, k_theta, stretch = self.jumps(N, r0 + 0 * v, v)
        self.assertEqual(stretch, 1)
        self.assertEqual(sp.simplify(k_theta), 0)
        F = 1 - 2 * N / r0
        stated = -sp.diff(N, v) / (4 * sp.pi * r0 * F ** sp.Rational(3, 2))
        self.assertEqual(sp.simplify((k_tau + k_theta) / (8 * sp.pi) - stated), 0)

    def test_the_drawn_surface_has_no_density_and_a_pressure(self):
        # R = 3N with N = cos^2(pi (v - 2)/12): g^rr = 1/3 on the surface and du/dv = 1 - 18 N'.
        sp = self.sp
        v = sp.Symbol("v", real=True)
        N = sp.cos(sp.pi * (v - 2) / 12) ** 2
        k_tau, k_theta, stretch = self.jumps(N, 3 * N, v)
        self.assertEqual(sp.simplify(stretch - (1 - 18 * sp.diff(N, v))), 0)
        for at in (sp.Rational(5, 2), 4, sp.Rational(11, 2), 7):
            self.assertLess(abs(complex(k_theta.subs(v, at).evalf(30))), 1e-20)
            self.assertGreater(float(((k_tau + k_theta) / (8 * sp.pi)).subs(v, at).evalf(30)), 0)

    def test_a_surface_with_unequal_masses_has_a_density(self):
        # The mass must be the same on both sides: at a fixed radius, sigma = (sqrt(F_in) - sqrt(F_out))/(4 pi r_0).
        sp = self.sp
        _, mi, Ri, _, kth_in = self.curvatures("ingoing")
        _, mo, Ro, _, kth_out = self.curvatures("outgoing")
        r0, a, b = sp.symbols("r_0 a b", positive=True)
        jump = kth_out.subs({Ro: r0, mo: b}).doit() - kth_in.subs({Ri: r0, mi: a}).doit()
        stated = (sp.sqrt(1 - 2 * a / r0) - sp.sqrt(1 - 2 * b / r0)) / (4 * sp.pi * r0)
        self.assertEqual(sp.simplify(-jump / (4 * sp.pi) - stated), 0)


@unittest.skipUnless(HAS_DRAWING, "numpy, scipy, sympy and contourpy are not all installed")
class DrawnModel(unittest.TestCase):
    def setUp(self):
        import numpy
        import null_rays
        self.np, self.nr = numpy, null_rays

    def test_the_mass_is_the_same_on_both_sides_of_the_surface(self):
        np, nr = self.np, self.nr
        v = np.linspace(2, 8, 601)
        N = nr.hiscock_mass(v)
        self.assertLess(float(np.max(np.abs(nr.hiscock_out(v - 18 * N) - N))), 1e-12)
        self.assertLess(float(np.max(np.abs(nr.hiscock_advanced(v - 18 * N) - v))), 1e-10)

    def test_the_retarded_time_runs_from_u1_to_u0(self):
        nr = self.nr
        self.assertEqual(float(nr.hiscock_out(-16.0)), 1.0)
        self.assertEqual(float(nr.hiscock_out(8.0)), 0.0)
        self.assertAlmostEqual(float(nr.hiscock_out(-16 + 1e-9)), 1.0, places=9)

    def test_hiscocks_condition_holds_at_the_end(self):
        # dM/du -> 0 as M -> 0, and dN/dv -> 0 as N -> 0.
        np, nr = self.np, self.nr
        u = 8 - np.array([1e-1, 1e-2, 1e-3, 1e-4])
        rate, mass = nr.hiscock_out(u, 1), nr.hiscock_out(u)
        self.assertTrue(np.all(np.diff(np.abs(rate)) < 0) and abs(rate[-1]) < 1e-4)
        self.assertTrue(np.all(mass > 0) and mass[-1] < 1e-9)
        self.assertTrue(np.all(rate < 0))

    def test_the_derivatives_of_the_outgoing_mass_are_its_own(self):
        np, nr = self.np, self.nr
        u, h = np.linspace(-15, 7, 23), 1e-5
        first = (nr.hiscock_out(u + h) - nr.hiscock_out(u - h)) / (2 * h)
        second = (nr.hiscock_out(u + h, 1) - nr.hiscock_out(u - h, 1)) / (2 * h)
        self.assertLess(float(np.max(np.abs(first - nr.hiscock_out(u, 1)))), 1e-8)
        self.assertLess(float(np.max(np.abs(second - nr.hiscock_out(u, 2)))), 1e-8)

    def test_the_outgoing_chart_begins_on_the_shell_and_then_on_the_surface(self):
        # Before u_1 the edge is the shell v = 0: u + 2r + 4 ln(r/2 - 1) + 12 + 4 ln 2 = 0 on it.
        np, nr = self.np, self.nr
        u = np.linspace(-40, -16.001, 50)
        r = nr.hiscock_edge(u)
        self.assertLess(float(np.max(np.abs(u + 2 * r + 4 * np.log(r / 2 - 1) + nr.HISCOCK_SHIFT))), 1e-10)
        late = np.linspace(-16, 8, 50)
        self.assertLess(float(np.max(np.abs(nr.hiscock_edge(late) - 3 * nr.hiscock_out(late)))), 1e-15)

    def test_the_slices_of_the_two_charts_meet_on_the_surface(self):
        nr = self.nr
        for T in (-1.0, 0.5, 3.0, 6.0, 7.9):
            v = nr.hiscock_meets(T)
            N = float(nr.hiscock_mass(v))
            self.assertAlmostEqual(v - 3 * N, T, places=10)                      # v - r = T on the surface
            self.assertAlmostEqual(nr.hiscock_outer_time(T), (v - 18 * N) + 3 * N, places=10)   # u + r there
        self.assertEqual(nr.hiscock_outer_time(9.0), 9.0)

    def test_the_event_horizon_is_inside_the_apparent_horizon(self):
        # The seed the spacetime diagram traces the event horizon from: inside r = 2N at v = 7.75.
        nr = self.nr
        v, r = (float(nr.number(x)) for x in nr.HISCOCK_HORIZON.values())
        self.assertLess(r, 2 * float(nr.hiscock_mass(v)))
        self.assertGreater(r, 0)
        self.assertAlmostEqual(math.cos(math.pi * (v - 2) / 12) ** 2, float(nr.hiscock_mass(v)), places=14)


if __name__ == "__main__":
    unittest.main()
