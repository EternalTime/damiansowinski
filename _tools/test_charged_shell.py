"""The charged shell of dust, flat space inside and Reissner-Nordstrom's field outside:
.venv.noindex/bin/python -m unittest discover -s _tools

What its texts state about the shell, held on the published files. The shell is in no chart, since
the time of one side is not the time of the other and g_rr jumps across it, so its motion is taken
here from the published metrics and Christoffel symbols of the charts by Israel's junction
conditions: the shell's world line has one proper time from both sides, the jump of K^theta_theta
is -4 pi G sigma/c^2 and the jump of K^tau_tau is +4 pi G sigma/c^2 for dust of surface density
sigma = m/(4 pi R^2), which together are the one first integral the History quotes,
M = m sqrt(1 + Rdot^2) + (Q^2 - m^2)/(2R) in G = c = 1. The shell at rest in the isotropic chart is
held to Arnowitt, Deser and Misner's relation between its masses and its radius. The motion of the
shell every diagram draws, charged_shell.py, is held to the same first integral and to the numbers
the captions state, and the embedding file to the surfaces its caption names. The tests need sympy
and numpy and are skipped where either is absent.
"""
import importlib.util
import json
import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "MFS" / "assets" / "data"
READY = all(importlib.util.find_spec(name) is not None for name in ("sympy", "numpy"))
sys.path.insert(0, str(Path(__file__).resolve().parent / "derivations"))


def chart_of(system):
    metric = json.loads((DATA / "metrics" / "charged_shell.json").read_text(encoding="utf-8"))
    return next(c for c in metric["coordinates"] if c["id"] == system)


@unittest.skipUnless(READY, "sympy and numpy are not installed")
class Shell(unittest.TestCase):
    def setUp(self):
        import sympy
        import verify_metrics
        import charged_shell
        self.sp, self.vm, self.shell = sympy, verify_metrics, charged_shell
        self.R, self.rs, self.rq, self.mu = sympy.symbols("R r_s r_q mu", positive=True)
        self.gamma = self.rs / (2 * self.mu) - (self.rq ** 2 - self.mu ** 2) / (2 * self.mu * self.R)
        self.beta = self.rs / (2 * self.mu) - (self.rq ** 2 + self.mu ** 2) / (2 * self.mu * self.R)
        self.f = 1 - self.rs / self.R + self.rq ** 2 / self.R ** 2

    def side(self, system, way):
        """One face of the shell in a published chart: the metric and the Christoffel symbols on the
        plane of the time and r at r = R, in this test's symbols and with c = 1, and the shell's
        velocity (dx^0/dtau, dR/dtau) there, on its way in (way = -1) or out (way = +1)."""
        entry = chart_of(system)
        reader = self.vm.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]], ())
        time, radius = entry["coords"][:2]
        mine = {reader.symbol["r"]: self.R, reader.parameters["r_s"]: self.rs, reader.parameters["r_q"]: self.rq,
                reader.c: 1}
        read = lambda text: self.sp.sympify(reader(text)).subs(mine)
        g = {tuple(e["indices"]): read(e["value"]) for e in entry["metric_components"]}
        gamma = {tuple(e["indices"]): read(e["value"]) for e in entry["christoffel"]["variants"]["ull"]["nonzero"]}
        names = (time, radius)
        h = self.sp.Matrix(2, 2, lambda i, j: g.get((names[i], names[j]), 0))
        G = [[[gamma.get((names[a], names[b], names[c]), 0) for c in range(2)] for b in range(2)] for a in range(2)]
        Rdot = way * self.sp.sqrt(self.gamma ** 2 - 1)
        rate = {"interior": self.gamma, "exterior": self.beta / self.f, "exterior_ingoing": 1 / (self.beta - Rdot),
                "exterior_outgoing": 1 / (self.beta + Rdot)}[system]
        return h, G, self.sp.Matrix([rate, Rdot]), self.sp.sqrt(g[("\\theta", "\\theta")]), Rdot

    def curvatures(self, system, way=-1):
        """K^theta_theta and K^tau_tau of the shell from one side, with the unit normal pointing to
        larger r where the shell is at rest: K^theta_theta = n^r/R and K^tau_tau = n_a a^a, a the
        acceleration of the dust along the shell's world line."""
        sp = self.sp
        h, G, u, areal, Rdot = self.side(system, way)
        root = sp.sqrt(-h.det())
        n = sp.Matrix([-root * u[1], root * u[0]])
        raised = h.inv() * n
        a = sp.Matrix([Rdot * sp.diff(u[i], self.R) + sum(G[i][b][c] * u[b] * u[c] for b in range(2) for c in range(2))
                       for i in range(2)])
        return raised[1] / self.R, (n.T * a)[0], (u.T * h * u)[0], (raised.T * n)[0], areal

    # Outside r_+, between the horizons and inside r_- of the shell drawn, a heavier hole, and a field
    # with no horizon at all.
    POINTS = ({"r_s": 1, "r_q": "12/25", "mu": "1/5", "R": 3}, {"r_s": 1, "r_q": "12/25", "mu": "1/5", "R": "1/2"},
              {"r_s": 1, "r_q": "12/25", "mu": "1/5", "R": "33/100"}, {"r_s": 2, "r_q": "9/10", "mu": "1/2", "R": 4},
              {"r_s": 1, "r_q": "3/5", "mu": "1/4", "R": 2})
    CHARTS = (("exterior", -1), ("exterior_ingoing", -1), ("exterior_outgoing", 1))

    def at(self, point):
        return {self.rs: self.sp.Rational(point["r_s"]), self.rq: self.sp.Rational(point["r_q"]),
                self.mu: self.sp.Rational(point["mu"]), self.R: self.sp.Rational(point["R"])}

    def zero(self, expression, values):
        self.assertLess(abs(complex(self.sp.sympify(expression).subs(values).evalf(30))), 1e-20)

    def test_the_shell_keeps_one_proper_time_in_every_chart(self):
        for system, way in (("interior", -1), ("interior", 1)) + self.CHARTS:
            _, _, speed, unit, areal = self.curvatures(system, way)
            self.assertEqual(areal, self.R, system)
            for point in self.POINTS:
                self.zero(speed + 1, self.at(point))
                self.zero(unit - 1, self.at(point))

    def test_the_jump_of_the_extrinsic_curvature_is_dust_of_constant_rest_mass(self):
        # [[K^theta_theta]] = -4 pi sigma and [[K^tau_tau]] = +4 pi sigma with 4 pi R^2 sigma = mu:
        # S_ab = sigma u_a u_b, no pressure in the shell, charged or not.
        for system, way in self.CHARTS:
            k_in, a_in = self.curvatures("interior", way)[:2]
            k_out, a_out = self.curvatures(system, way)[:2]
            for point in self.POINTS:
                self.zero(k_out - k_in + self.mu / self.R ** 2, self.at(point))
                self.zero(a_out - a_in - self.mu / self.R ** 2, self.at(point))

    def test_the_mass_outside_is_the_energy_of_the_shell(self):
        # M = m sqrt(1 + Rdot^2) + (Q^2 - m^2)/(2R), with r_s = 2M, r_q = Q and mu = m.
        energy = self.mu * self.gamma + (self.rq ** 2 - self.mu ** 2) / (2 * self.R)
        self.assertEqual(self.sp.simplify(energy - self.rs / 2), 0)
        # and beta^2 - f is the same Rdot^2 as gamma^2 - 1
        self.assertEqual(self.sp.simplify(self.beta ** 2 - self.f - self.gamma ** 2 + 1), 0)

    def test_without_the_charge_it_is_israels_shell(self):
        import israel_shell
        for R in (0.7, 1.5, 4.0):
            self.assertAlmostEqual(float(self.shell.gamma(R, 1.0, 0.0, 0.5)), float(israel_shell.gamma(R, 1.0, 0.5)))
            self.assertAlmostEqual(float(self.shell.beta(R, 1.0, 0.0, 0.5)), float(israel_shell.beta(R, 1.0, 0.5)))

    def test_a_shell_with_more_charge_than_mass_turns_inside_the_inner_horizon(self):
        # At R = (r_q^2 - mu^2)/(r_s - 2 mu) the shell is at rest, and for a black hole outside,
        # r_q < r_s/2, that radius is at most r_-, with equality at mu = r_-.
        sp = self.sp
        rest = (self.rq ** 2 - self.mu ** 2) / (self.rs - 2 * self.mu)
        self.assertEqual(sp.simplify((self.gamma ** 2 - 1).subs(self.R, rest)), 0)
        for rq in (sp.Rational(1, 10), sp.Rational(3, 10), sp.Rational(12, 25), sp.Rational(499, 1000)):
            inner = sp.Rational(1, 2) - sp.sqrt(sp.Rational(1, 4) - rq ** 2)
            for k in range(1, 20):
                mu = rq * k / 20
                self.assertLessEqual(float(rest.subs({self.rs: 1, self.rq: rq, self.mu: mu})), float(inner) + 1e-15)
            self.zero(rest.subs({self.rs: 1, self.rq: rq, self.mu: inner}) - inner, {})

    def test_the_shell_at_rest_in_the_isotropic_chart_is_arnowitt_deser_and_misners(self):
        # At rest gamma = 1 and beta = sqrt(f), so the junction condition is 1 - sqrt(-g_tt) = mu/R on
        # the shell rho = epsilon, with R the areal radius there. The published chart gives
        # mu = a + b + 2ab/epsilon, which is their (7.1), M = -eps + sqrt(eps^2 + 2 m_0 eps + Q^2), with
        # r_s = 2(a + b) and r_q = a - b.
        sp = self.sp
        entry = chart_of("exterior_isotropic")
        reader = self.vm.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]], ())
        g = {tuple(e["indices"]): sp.sympify(reader(e["value"])) for e in entry["metric_components"]}
        rho, a, b = reader.symbol["\\rho"], reader.parameters["a"], reader.parameters["b"]
        areal = (rho + a) * (rho + b) / rho
        lapse = (rho ** 2 - a * b) / ((rho + a) * (rho + b))
        self.assertEqual(sp.simplify(areal ** 2 - g[("\\theta", "\\theta")]), 0)
        self.assertEqual(sp.simplify(lapse ** 2 + g[("t", "t")]), 0)
        mass = sp.simplify((1 - lapse) * areal)
        self.assertEqual(sp.simplify(mass - (a + b + 2 * a * b / rho)), 0)
        rs, rq = 2 * (a + b), a - b
        self.assertEqual(sp.simplify((rs / 2 + rho) ** 2 - (rho ** 2 + 2 * mass * rho + rq ** 2)), 0)
        # The same shell by its areal radius: M = m + (Q^2 - m^2)/(2R).
        self.assertEqual(sp.simplify(mass + (rq ** 2 - mass ** 2) / (2 * areal) - rs / 2), 0)
        for eps in (2.0, 0.3, 0.01):
            r_s, R = self.shell.point_charge(eps, 0.5)
            self.assertAlmostEqual(float(self.shell.mass_at_rest(R, 1.0, 0.5)), float(r_s))
            # The balanced shell weighs its rest mass at every radius.
            self.assertAlmostEqual(float(self.shell.point_charge(eps, 1.0)[0]), 2.0)

    def test_the_shell_drawn_is_that_motion(self):
        import numpy as np
        shell = self.shell
        eta = np.array([-9.0, -4.0, -shell.ETA_P, -1.5, -shell.ETA_M, -0.4, 0.3, 0.9])
        R = shell.radius(eta)
        h = 1e-6
        per = lambda f: (f(eta + h) - f(eta - h)) / (shell.proper_time(eta + h) - shell.proper_time(eta - h))
        gamma, beta = shell.gamma(R, 1.0, shell.RQ, shell.MU), shell.beta(R, 1.0, shell.RQ, shell.MU)
        self.assertLess(float(np.max(np.abs(per(shell.radius) - shell.speed(eta)))), 1e-8)
        self.assertLess(float(np.max(np.abs(per(shell.inner_time) - gamma))), 1e-8)
        self.assertLess(float(np.max(np.abs(shell.speed(eta) ** 2 - gamma ** 2 + 1))), 1e-12)
        self.assertLess(float(np.max(np.abs(per(shell.advanced) * (beta - shell.speed(eta)) - 1))), 1e-7)
        # The way out is the way in turned over in time: u(eta) = -v(-eta), and du/dtau = 1/(beta + Rdot).
        self.assertLess(float(np.max(np.abs(shell.retarded(-eta) + shell.advanced(eta)))), 1e-12)
        eta = -eta
        self.assertLess(float(np.max(np.abs(per(shell.retarded) * (beta + shell.speed(eta)) - 1))), 1e-7)
        self.assertAlmostEqual(float(shell.outer_time(0.0)), 0.0)

    def test_the_numbers_the_captions_state(self):
        shell = self.shell
        self.assertAlmostEqual(shell.R_TURN, 119 / 375)
        self.assertAlmostEqual(shell.R_TURN, shell.turn(1.0, shell.RQ, shell.MU))
        self.assertAlmostEqual(float(shell.radius(-shell.ETA_P)), 16 / 25)
        self.assertAlmostEqual(float(shell.radius(shell.ETA_M)), 9 / 25)
        for value, stated in ((shell.R_TURN, 0.317), (shell.ROOT_A / shell.E, 0.92),
                              (-float(shell.proper_time(-shell.ETA_P)), 0.39),
                              (2 * float(shell.proper_time(shell.ETA_P)), 0.79),
                              (2 * float(shell.proper_time(shell.ETA_M)), 0.27),
                              (float(shell.inner_time(-shell.ETA_P)), -0.53),
                              (float(shell.inner_retarded(-shell.ETA_P)), -1.17),
                              (float(shell.inner_advanced(shell.ETA_M)), 0.50), (shell.CONFORMAL_LENGTH, 0.504),
                              (float(shell.advanced(-shell.ETA_P)), -0.02), (float(shell.advanced(-shell.ETA_M)), 0.11),
                              (float(shell.advanced(0.0)), 0.30), (float(shell.retarded(shell.ETA_P)), 0.02),
                              (float(shell.slice_time(0.0)), -0.015), (float(shell.slice_time(-shell.ETA_P)), -0.66),
                              (float(shell.slice_time(-shell.ETA_M)), -0.25)):
            digits = len(repr(stated).split(".")[1])
            self.assertEqual(round(value, digits), stated)

    def test_each_time_is_inverted_to_the_shell_it_came_from(self):
        import numpy as np
        shell = self.shell
        eta = np.array([-7.0, -3.0, -2.0, -1.0, -0.3, 0.2, 0.8])
        for forward, back in ((shell.inner_time, shell.eta_of_inner_time), (shell.advanced, shell.eta_of_advanced),
                              (shell.inner_advanced, shell.eta_of_inner_advanced),
                              (shell.inner_retarded, shell.eta_of_inner_retarded), (shell.slice_time, shell.eta_of_slice)):
            self.assertLess(float(np.max(np.abs(back(forward(eta)) - eta))), 1e-9)
        self.assertLess(float(np.max(np.abs(shell.eta_of_retarded(shell.retarded(-eta)) + eta))), 1e-9)
        for region, values in (("I", [-7.0, -3.0, -2.3]), ("III", [-0.8, 0.0, 0.5]), ("I'", [2.3, 3.0, 7.0])):
            values = np.array(values)
            self.assertLess(float(np.max(np.abs(shell.eta_of_outer_time(shell.outer_time(values), region) - values))), 1e-9)

    def test_the_conformal_map_is_one_on_the_shell_and_straight_on_its_horizons(self):
        import numpy as np
        shell = self.shell
        eta = np.array([-8.0, -3.0, -2.2, -1.5, -0.9, 0.0, 0.5, 0.9])
        inside = np.array(shell.conformal_inside(shell.inner_time(eta), shell.radius(eta)))
        self.assertLess(float(np.max(np.abs(inside - shell.conformal_ingoing(shell.advanced(eta), shell.radius(eta))))), 1e-8)
        mirror = np.array(shell.conformal_inside(shell.inner_time(-eta), shell.radius(eta)))
        self.assertLess(float(np.max(np.abs(mirror - shell.conformal_outgoing(shell.retarded(-eta), shell.radius(eta))))), 1e-8)
        v = np.array([0.5, 2.0, 9.0])
        p, _ = shell.conformal_ingoing(v, np.full(3, shell.RM))
        self.assertLess(float(np.max(np.abs(p + math.pi / 4))), 1e-12)              # the inner horizon, going in
        _, q = shell.conformal_ingoing(np.array([80.0, 300.0]), np.array([0.5, 3.0]))
        self.assertLess(float(np.max(np.abs(q - math.pi / 4))), 1e-9)               # the Cauchy horizon
        p, q = shell.conformal_inside(np.array([-5.0, 0.0, 5.0]), np.zeros(3))
        self.assertLess(float(np.max(np.abs(p - q))), 1e-15)                        # the centre, X = 0
        p, q = shell.conformal_beyond(np.array([-1.0, 0.0, 2.0]), np.full(3, 0.2))
        back = shell.conformal_beyond(-np.array([-1.0, 0.0, 2.0]), np.full(3, 0.2))
        self.assertLess(float(np.max(np.abs(p + back[1]))), 1e-12)                  # turned over in time


class Embedding(unittest.TestCase):
    """The embedding file, held to the surfaces its captions name, from its own numbers."""

    def setUp(self):
        self.views = {v["id"]: v for v in json.loads((DATA / "embedding" / "charged_shell.json").read_text(encoding="utf-8"))["views"]}

    def test_the_falling_shell_is_a_flat_disc_inside_and_the_slice_of_constant_v_minus_r_outside(self):
        rq = 12 / 25

        def height(r):
            root = math.sqrt(r - rq * rq)
            return 2 * root - 2 * rq * math.atan(root / rq)
        view = self.views["bounce"]
        self.assertEqual([s["time"] for s in view["surfaces"]], [-1.5, -0.5, 0.0, 1.0])
        radii = []
        for surface in view["surfaces"]:
            inside, outside = (next(p for p in surface["pieces"] if p["id"] == name)["points"] for name in ("inside", "outside"))
            self.assertTrue(all(z == inside[0][2] for _, _, z in inside))
            self.assertAlmostEqual(inside[-1][1], outside[0][1], places=6)
            self.assertAlmostEqual(inside[-1][2], outside[0][2], places=6)
            self.assertEqual(outside[-1][2], 0)
            shell = outside[0][0]
            radii.append(shell)
            for r, rho, z in outside:
                self.assertAlmostEqual(rho, r, places=6)
                self.assertAlmostEqual(z - outside[0][2], height(r) - height(shell), places=5)
        # Outside r_+, between the horizons, just past the turn at 119/375, and climbing back toward r_-.
        self.assertGreater(radii[0], 16 / 25)
        self.assertTrue(9 / 25 < radii[1] < 16 / 25)
        self.assertTrue(119 / 375 < radii[2] < radii[3] < 9 / 25)

    def test_the_balanced_shell_cuts_the_extreme_throat_off_with_a_flat_disc(self):
        def height(rho):
            root = math.sqrt(2 * (rho + 1) - 1)
            return 2 * root + math.log((root - 1) / (root + 1))
        (surface,) = self.views["point"]["surfaces"]
        pieces = {p["id"]: p for p in surface["pieces"]}
        inside, outside, throat = (pieces[name]["points"] for name in ("inside", "outside", "throat"))
        self.assertTrue(pieces["throat"]["reference"])
        self.assertTrue(all(z == inside[0][2] for _, _, z in inside))
        self.assertEqual((outside[0][0], inside[-1][0]), (1.0, 2.0))
        self.assertAlmostEqual(inside[-1][2], outside[0][2], places=6)
        for rho, radius, z in outside + throat:
            self.assertAlmostEqual(radius, rho + 1, places=6)                    # the areal radius is rho + a
            self.assertAlmostEqual(z, height(rho) - height(3.0), places=5)
        self.assertEqual(throat[0][0], 1 / 64)


if __name__ == "__main__":
    unittest.main()
