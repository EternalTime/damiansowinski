"""Bonnor and Faulkes's oscillating sphere of uniform density with pressure:
.venv.noindex/bin/python -m unittest discover -s _tools

What its texts and drawings state, held on the published files and on Bonnor and Faulkes's own
equation (2.8), Fdot^2 = Q(F), in units of L with c = 1: the oscillation (4.7), alpha = 0.7,
r_0^2/L^2 = 4/9, m = 0.1084L, turns round at F = 0.84 and F = 1.20, as they found; its period, the
integral of 2 dF/sqrt(Q) over that swing, is 18.5L/c, where they printed 17.4L/c, and the history
states both; the small oscillations of the static star have the period 2 pi sqrt(6)L/c/sqrt((1 -
alpha)(1 + 2 alpha)), the limit of that integral as the swing shrinks; the drawings' solver follows
the same F; and the published Einstein tensor of the comoving chart is a fluid of uniform density
whose pressure vanishes on r = r_0. The tests need sympy and mpmath, and are skipped where absent.
"""
import importlib.util
import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "MFS" / "assets" / "data"
HAS_SYMPY = importlib.util.find_spec("sympy") is not None
sys.path.insert(0, str(Path(__file__).resolve().parent / "derivations"))

ALPHA, R0, M = "0.7", "2/3", "0.1084"


def metric():
    return json.loads((DATA / "metrics" / "oscillating_sphere.json").read_text(encoding="utf-8"))


def well(alpha, r0, m):
    """Q(F) of Bonnor and Faulkes's (2.8) with L = 1, as an mpmath function."""
    def Q(F):
        Y0 = F ** alpha - 1 + alpha + alpha * r0 ** 2 / 4
        return 2 * m * Y0 ** 3 / (r0 ** 3 * alpha ** 3 * F ** (3 * alpha + 1)) - (F ** -alpha - (1 - alpha) * F ** (-2 * alpha)) / alpha
    return Q


def period(mp, Q, lo, hi):
    """2 times the integral of dF/sqrt(Q) from lo to hi, with F = mid + half sin(theta), which takes
    the inverse square root of each simple zero of Q at the turning points into the measure; the
    absolute value keeps a root rounded to 30 digits from leaving a negative Q at the very ends."""
    mid, half = (lo + hi) / 2, (hi - lo) / 2
    return 2 * mp.quad(lambda th: half * mp.cos(th) / mp.sqrt(abs(Q(mid + half * mp.sin(th)))), [-mp.pi / 2, 0, mp.pi / 2])


@unittest.skipUnless(HAS_SYMPY, "sympy is not installed")
class Period(unittest.TestCase):
    def setUp(self):
        from mpmath import mp, mpf, fraction
        mp.dps = 30
        self.mp = mp
        self.alpha, self.r0, self.m = mpf(ALPHA), fraction(2, 3), mpf(M)
        self.Q = well(self.alpha, self.r0, self.m)

    def swing(self, Q, guesses):
        lo, hi = (self.mp.findroot(Q, g) for g in guesses)
        return lo, hi, period(self.mp, Q, lo, hi)

    def test_the_turning_points_are_theirs(self):
        lo, hi, _ = self.swing(self.Q, (0.84, 1.2))
        self.assertAlmostEqual(float(lo), 0.84, places=2)
        self.assertAlmostEqual(float(hi), 1.20, places=2)

    def test_the_period_is_18_5_and_the_history_says_so(self):
        _, _, T = self.swing(self.Q, (0.84, 1.2))
        self.assertAlmostEqual(float(T), 18.516, places=3)
        self.assertEqual(round(float(T), 1), 18.5)
        history = metric()["history"]
        self.assertIn("$18.5L/c$", history)
        self.assertIn("$17.4L/c$", history)
        # Their 7.6 years at a mean density of 1e-8 g/cm^3 scales with the period, and the mean
        # density came from the turning points alone, which agree.
        self.assertEqual(round(7.6 * float(T) / 17.4, 1), 8.1)
        self.assertIn("8.1 years", history)

    def test_small_swings_tend_to_the_printed_period(self):
        # Constants a little away from the static star's: with r_0 and alpha fixed, the static star
        # has 1/L^2 = (2m/r_0^3)(1 + mu_0^2)^3, their (3.2) with alpha = (1 - 2mu_0^2)/(1 + mu_0^2),
        # so mu_0^2 = (1 - alpha)/(2 + alpha) and L = 1 fixes r_0 and the static m.
        mp = self.mp
        alpha = self.alpha
        mu2 = (1 - alpha) / (2 + alpha)
        r0 = 2 * mp.sqrt(mu2)
        m_static = r0 ** 3 / (2 * (1 + mu2) ** 3)
        Q = well(alpha, r0, m_static * (1 + mp.mpf("1e-7")))
        grid = [mp.mpf(1) + k * mp.mpf("1e-5") for k in range(-200, 201)]
        signs = [(grid[i], grid[i + 1]) for i in range(len(grid) - 1) if Q(grid[i]) * Q(grid[i + 1]) < 0]
        lo = mp.findroot(Q, signs[0])
        hi = mp.findroot(Q, signs[-1])
        T = period(mp, Q, lo, hi)
        printed = 2 * mp.pi * mp.sqrt(6) / mp.sqrt((1 - alpha) * (1 + 2 * alpha))
        self.assertAlmostEqual(float(T), float(printed), places=3)

    def test_the_drawings_solver_follows_the_same_swing(self):
        import numpy as np
        import null_rays
        solver = null_rays.OscillationSolver("comoving", {"alpha": "7/10", "r_0": "2/3", "m": M})
        _, hi, T = self.swing(self.Q, (0.84, 1.2))
        self.assertAlmostEqual(solver.F_max, float(hi), places=9)
        F, _, _ = solver.values("F", np.array([0.0, float(T) / 2, float(T)]))
        self.assertAlmostEqual(F[0], float(hi), places=9)
        self.assertAlmostEqual(F[2], float(hi), places=6)
        self.assertLess(F[1], 0.85)


@unittest.skipUnless(HAS_SYMPY, "sympy is not installed")
class Fluid(unittest.TestCase):
    def test_the_density_is_uniform_and_the_pressure_vanishes_on_the_surface(self):
        import sympy as sp
        import verify_metrics as vm
        entry = next(c for c in metric()["coordinates"] if c["id"] == "comoving")
        reader = vm.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]], ())
        ul = entry["einstein_tensor"]["variants"]["ul"]["nonzero"]
        G = {tuple(c["indices"]): reader(c["value"]) for c in ul}
        t, r = reader.symbol["t"], reader.symbol["r"]
        al, L, m, r0, F = (reader.parameters[n] for n in ("alpha", "L", "m", "r_0", "F"))
        values = {al: sp.Rational(7, 10), L: 1, m: sp.Rational(1084, 10000), r0: sp.Rational(2, 3)}
        Fv = sp.Rational(21, 20)
        Q = well(*[sp.Float(x, 30) for x in (0.7, 2 / 3, 0.1084)])
        Qe = sp.sympify(Q(sp.Symbol("x", positive=True)))

        def at(value, radius):
            value = value.subs(reader.held).doit().subs(values) if reader.held else value.subs(values)
            x = sp.Symbol("x", positive=True)
            rate2 = Qe.subs(x, Fv)
            accel = sp.diff(Qe, x).subs(x, Fv) / 2
            value = value.subs(sp.Derivative(F, (t, 2)), accel).subs(sp.Derivative(F, t), sp.sqrt(rate2)).subs(F, Fv)
            return complex(sp.N(value.subs({r: radius, reader.c: 1}), 30))
        density = [at(-G[("t", "t")], x) for x in (sp.Rational(1, 10), sp.Rational(1, 3), sp.Rational(3, 5))]
        for d in density[1:]:
            self.assertAlmostEqual(abs(d - density[0]), 0.0, places=12)
        self.assertAlmostEqual(abs(at(G[("r", "r")], sp.Rational(2, 3))), 0.0, places=12)
        self.assertGreater(at(G[("r", "r")], sp.Rational(1, 3)).real, 0.0)


if __name__ == "__main__":
    unittest.main()
