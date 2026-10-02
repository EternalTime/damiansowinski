"""Israel's collapsing shell of dust, flat space inside and Schwarzschild's vacuum outside:
python3 -m unittest discover -s _tools

What its texts state about the shell, held on the published files. The shell is in no chart, since
the time of one side is not the time of the other and g_rr jumps across it, so its motion is taken
here from the published metrics and Christoffel symbols of the three charts by Israel's junction
conditions: the shell's world line has one proper time from both sides, the jump of K^theta_theta
is -4 pi G sigma/c^2 and the jump of K^tau_tau is +4 pi G sigma/c^2 for dust of surface density
sigma = m/(4 pi R^2), which together are the one first integral the History quotes,
M = m sqrt(1 + Rdot^2) - m^2/(2R) in G = c = 1. The closed forms of the shell every diagram draws,
israel_shell.py, are held to the same motion and to the numbers the captions state. The tests need
sympy and numpy and are skipped where either is absent.
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


@unittest.skipUnless(READY, "sympy and numpy are not installed")
class Shell(unittest.TestCase):
    def setUp(self):
        import sympy
        import verify_metrics
        import israel_shell
        self.sp, self.vm, self.shell = sympy, verify_metrics, israel_shell
        self.R, self.rs, self.mu = sympy.symbols("R r_s mu", positive=True)
        self.gamma = self.rs / (2 * self.mu) + self.mu / (2 * self.R)
        self.beta = self.rs / (2 * self.mu) - self.mu / (2 * self.R)
        self.Rdot = -sympy.sqrt(self.gamma ** 2 - 1)

    def side(self, system):
        """One face of the shell in a published chart: the metric and the Christoffel symbols on the
        plane of the time and r at r = R, in this test's symbols and with c = 1, and the shell's
        velocity (dx^0/dtau, dR/dtau) there."""
        metric = json.loads((DATA / "metrics" / "israel_shell.json").read_text(encoding="utf-8"))
        entry = next(c for c in metric["coordinates"] if c["id"] == system)
        reader = self.vm.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]], ())
        time, radius = entry["coords"][:2]
        mine = {reader.symbol["r"]: self.R, reader.parameters["r_s"]: self.rs, reader.c: 1}
        read = lambda text: self.sp.sympify(reader(text)).subs(mine)
        g = {tuple(e["indices"]): read(e["value"]) for e in entry["metric_components"]}
        gamma = {tuple(e["indices"]): read(e["value"]) for e in entry["christoffel"]["variants"]["ull"]["nonzero"]}
        names = (time, radius)
        h = self.sp.Matrix(2, 2, lambda i, j: g.get((names[i], names[j]), 0))
        G = [[[gamma.get((names[a], names[b], names[c]), 0) for c in range(2)] for b in range(2)] for a in range(2)]
        f = 1 - self.rs / self.R
        rate = {"interior": self.gamma, "exterior": self.beta / f, "exterior_ingoing": 1 / (self.beta - self.Rdot)}[system]
        areal = self.sp.sqrt(g[("\\theta", "\\theta")])
        return h, G, self.sp.Matrix([rate, self.Rdot]), areal

    def curvatures(self, system):
        """K^theta_theta and K^tau_tau of the shell from one side, with the unit normal pointing to
        larger r where the shell is at rest: K^theta_theta = n^r/R and K^tau_tau = n_a a^a, a the
        acceleration of the dust along the shell's world line."""
        sp = self.sp
        h, G, u, areal = self.side(system)
        # n_a = sqrt(-det h) epsilon_ab u^b is the unit normal orthogonal to u.
        root = sp.sqrt(-h.det())
        n = sp.Matrix([-root * u[1], root * u[0]])
        raised = h.inv() * n
        a = sp.Matrix([self.Rdot * sp.diff(u[i], self.R) + sum(G[i][b][c] * u[b] * u[c] for b in range(2) for c in range(2))
                       for i in range(2)])
        return raised[1] / self.R, (n.T * a)[0], (u.T * h * u)[0], (raised.T * n)[0], areal

    POINTS = ({"r_s": 1, "mu": "1/2", "R": 3}, {"r_s": 1, "mu": "1/2", "R": "1/2"}, {"r_s": 2, "mu": "5/4", "R": 4},
              {"r_s": 1, "mu": "2/5", "R": "7/3"}, {"r_s": 1, "mu": "3/5", "R": "9/10"})

    def at(self, point):
        return {self.rs: self.sp.Rational(point["r_s"]), self.mu: self.sp.Rational(point["mu"]),
                self.R: self.sp.Rational(point["R"])}

    def zero(self, expression, values):
        self.assertLess(abs(complex(self.sp.sympify(expression).subs(values).evalf(30))), 1e-20)

    def test_the_shell_keeps_one_proper_time_in_every_chart(self):
        for system in ("interior", "exterior", "exterior_ingoing"):
            _, _, speed, unit, areal = self.curvatures(system)
            self.assertEqual(areal, self.R, system)
            for point in self.POINTS:
                self.zero(speed + 1, self.at(point))
                self.zero(unit - 1, self.at(point))

    def test_the_jump_of_the_extrinsic_curvature_is_dust_of_constant_rest_mass(self):
        # [[K^theta_theta]] = -4 pi sigma and [[K^tau_tau]] = +4 pi sigma with 4 pi R^2 sigma = mu:
        # S_ab = sigma u_a u_b, no pressure in the shell.
        k_in, a_in = self.curvatures("interior")[:2]
        for system in ("exterior", "exterior_ingoing"):
            k_out, a_out = self.curvatures(system)[:2]
            for point in self.POINTS:
                self.zero(k_out - k_in + self.mu / self.R ** 2, self.at(point))
                self.zero(a_out - a_in - self.mu / self.R ** 2, self.at(point))

    def test_the_two_schwarzschild_charts_see_one_shell(self):
        outer, ingoing = self.curvatures("exterior"), self.curvatures("exterior_ingoing")
        for point in self.POINTS:
            self.zero(outer[0] - ingoing[0], self.at(point))
            self.zero(outer[1] - ingoing[1], self.at(point))

    def test_the_mass_outside_is_the_energy_of_the_shell(self):
        # M = m sqrt(1 + Rdot^2) - m^2/(2R), with r_s = 2M and mu = m.
        energy = self.mu * self.sp.sqrt(1 + self.Rdot ** 2) - self.mu ** 2 / (2 * self.R)
        for point in self.POINTS:
            self.zero(energy - self.rs / 2, self.at(point))

    def test_a_shell_heavier_than_its_mass_outside_turns_back(self):
        # With M < m the shell is at rest at R = m^2/(2(m - M)).
        rest = self.mu ** 2 / (2 * (self.mu - self.rs / 2))
        values = {self.rs: 1, self.mu: self.sp.Rational(3, 5)}
        self.zero((self.gamma ** 2 - 1).subs(self.R, rest), values)
        self.assertGreater(float(rest.subs(values)), 1)

    def test_the_shell_drawn_is_that_motion_in_closed_form(self):
        sp = self.sp
        s = sp.Symbol("s", positive=True)
        drawn = {self.rs: 1, self.mu: sp.Rational(1, 2), self.R: (s ** 2 - 1) / 8}
        radius = (s ** 2 - 1) / 8
        proper = -(s - 1) ** 2 * (s + 2) / 24
        inner = -(s ** 3 + 3 * s - 4) / 24
        advanced = -(s ** 3 / 3 - s ** 2 + 5 * s - 15 - 16 * sp.log((s + 3) / 6)) / 8
        per = lambda e: sp.diff(e, s) / sp.diff(proper, s)
        for value in (sp.Rational(11, 10), 2, 3 + sp.Rational(1, 7), 9):
            at = {s: value}
            self.zero((per(radius) - self.Rdot.subs(drawn)).subs(at), {})
            self.zero((per(inner) - self.gamma.subs(drawn)).subs(at), {})
            self.zero((per(advanced) - (1 / (self.beta - self.Rdot)).subs(drawn)).subs(at), {})
        for name, form in (("radius", radius), ("proper_time", proper), ("inner_time", inner), ("advanced", advanced)):
            for value in (1.3, 2.0, 3.0, 7.5):
                self.assertAlmostEqual(float(getattr(self.shell, name)(value)), float(form.subs(s, value)), places=12)

    def test_the_numbers_the_captions_state(self):
        shell = self.shell
        self.assertAlmostEqual(float(shell.radius(3.0)), 1.0)
        self.assertAlmostEqual(float(shell.proper_time(3.0)), -5 / 6)
        self.assertAlmostEqual(float(shell.inner_time(3.0)), -4 / 3)
        self.assertAlmostEqual(float(shell.inner_retarded(3.0)), shell.HORIZON_U)
        self.assertAlmostEqual(shell.HORIZON_U, -7 / 3)
        self.assertAlmostEqual(float(shell.advanced(3.0)), 0.0)
        self.assertAlmostEqual(float(shell.speed(1.0, 1.0, 0.5)), -0.75)
        self.assertAlmostEqual(float(shell.advanced(1.0)), shell.V_END)
        self.assertEqual(round(shell.V_END, 2), 0.52)
        self.assertAlmostEqual(float(shell.radius(1.0)), 0.0)
        self.assertAlmostEqual(float(shell.inner_time(1.0)), 0.0)

    def test_each_time_is_inverted_to_the_shell_it_came_from(self):
        import numpy as np
        shell = self.shell
        s = np.array([1.2, 2.0, 3.5, 6.0, 20.0])
        for forward, back in ((shell.inner_time, shell.s_of_inner_time), (shell.advanced, shell.s_of_advanced),
                              (shell.inner_advanced, shell.s_of_inner_advanced),
                              (shell.inner_retarded, shell.s_of_inner_retarded),
                              (lambda x: shell.advanced(x) - shell.radius(x), shell.s_of_slice)):
            self.assertLess(float(np.max(np.abs(back(forward(s)) - s))), 1e-9)
        outside = s[s > 3]
        self.assertLess(float(np.max(np.abs(shell.s_of_outer_time(shell.outer_time(outside)) - outside))), 1e-9)

    def test_the_conformal_map_is_one_on_the_shell_and_straight_on_its_edges(self):
        import numpy as np
        shell = self.shell
        s = np.array([1.3, 2.0, 2.9, 3.0, 3.4, 8.0, 30.0])
        inside = shell.conformal_inside(shell.inner_time(s), shell.radius(s))
        outside = shell.conformal_outside(shell.advanced(s), shell.radius(s))
        self.assertLess(float(np.max(np.abs(np.array(inside) - np.array(outside)))), 1e-9)
        late = np.array([0.6, 1.0, 4.0, 20.0])
        p, q = shell.conformal_outside(late, np.full(4, 1e-10))
        self.assertLess(float(np.max(np.abs(p + q))), 1e-6)                     # the singularity, T = 0
        p, _ = shell.conformal_outside(late, np.ones(4))
        self.assertLess(float(np.max(np.abs(p + math.pi / 4))), 1e-12)          # r_s outside the shell
        p, q = shell.conformal_inside(np.array([-9.0, -7 / 3, -0.2]), np.zeros(3))
        self.assertLess(float(np.max(np.abs(p - q))), 1e-15)                    # the centre, X = 0
        self.assertAlmostEqual(float(p[1]), -math.pi / 4)                       # where the horizon starts


if __name__ == "__main__":
    unittest.main()
