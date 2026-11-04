"""Mukhanov and Brandenberger's universe with limiting curvature, held on the theory's action, on
the numbers limiting_curvature.py solves and on the published files:
.venv.noindex/bin/python -m unittest discover -s _tools -p test_limiting_curvature.py

What the page states. On a Robertson-Walker universe the invariants are I_2 = 12(dH/dt - k/a^2)^2
and, with the root taken with the sign of dH/dt - k/a^2, I_1 = 12(H^2 + k/a^2); with that root the
minisuperspace action gives the paper's constraints (2.25) and (2.26), and with the positive root
neither. The lapse gives the Hamiltonian equation with the coefficient 2 phi_1 + 1, which the paper
misprints, and the scale factor the second order equation the solver integrates. The potentials
(4.7) and (4.8) keep H^2 + k/a^2 below H_0^2, so a closed universe never falls below c/H_0, keep
|dH/dt - k/a^2| at most 3 sqrt(3) H_0^2/8, and give de Sitter space as |phi_2| grows. The closed
universe drawn bounces, H changing sign at a finite radius, has finite curvature at every time,
and approaches de Sitter space on both sides. Every comparison of the action's equations is made at
points, with numbers for the scale factor, the fields and their derivatives.
Needs sympy, numpy and scipy, and is skipped where they are absent.
"""
import importlib.util
import json
import math
import random
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "MFS" / "assets" / "data"
HAS_ALL = all(importlib.util.find_spec(name) is not None for name in ("sympy", "numpy", "scipy"))
sys.path.insert(0, str(Path(__file__).resolve().parent / "derivations"))

ACTION = {}


def action():
    """The Euler-Lagrange expressions of N a^3 [R + phi_1 I_1 - phi_2 S - V_1 - V_2] for N, a, phi_1
    and phi_2, with N = 1 after the variation, for each sign of the root S = sign sqrt(12)(dH/dt -
    k/a^2), in the site's signature: R = 6(dH/dt + 2H^2 + k/a^2) and I_1 = R - sqrt(3) S, every
    derivative along the time taken with the lapse, H = (da/dt)/(N a). The potentials are (4.7) and
    (4.8) at H_0 = 1. Returns the symbols and, for each sign, the four expressions."""
    if ACTION:
        return ACTION
    import sympy as sp
    t, k = sp.symbols("t k")
    N, a = sp.Function("N")(t), sp.Function("a")(t)
    p1, p2 = sp.Function("phi1")(t), sp.Function("phi2")(t)
    V1 = 12 * p1 ** 2 / (1 + p1) * (1 - sp.log(1 + p1) / (1 + p1))
    V2 = 2 * sp.sqrt(3) * p2 ** 2 / (1 + p2 ** 2)
    HN = sp.diff(a, t) / (a * N)
    HdN = sp.diff(HN, t) / N
    R = 6 * (HdN + 2 * HN ** 2 + k / a ** 2)
    out = {}
    for sign in (1, -1):
        S = sign * sp.sqrt(12) * (HdN - k / a ** 2)
        L = N * a ** 3 * (R + p1 * (R - sp.sqrt(3) * S) - p2 * S - V1 - V2)

        def euler(q):
            e = (sp.diff(L, q) - sp.diff(sp.diff(L, sp.diff(q, t)), t)
                 + sp.diff(sp.diff(L, sp.diff(q, t, 2)), t, 2))
            return e.subs(N, 1).doit()
        out[sign] = {name: euler(q) for name, q in (("N", N), ("a", a), ("phi1", p1), ("phi2", p2))}
    ACTION.update(symbols=(t, k, a, p1, p2), equations=out)
    return ACTION


def at(expr, values):
    """An expression of the action at a point: every derivative of a, phi_1 and phi_2 set to its
    number, the highest first, then the functions themselves."""
    import sympy as sp
    t, k, a, p1, p2 = action()["symbols"]
    out = expr
    for f, name in ((a, "a"), (p1, "p1"), (p2, "p2")):
        for n in range(4, 0, -1):
            if f"{name}{n}" in values:
                out = out.subs(sp.Derivative(f, (t, n)), values[f"{name}{n}"])
    out = out.subs({a: values["a0"], p1: values["p10"], p2: values["p20"], k: values["k"]})
    return complex(sp.N(out, 30)).real


def on_shell(rng, k, sign=1):
    """A point at which (2.25), (2.26) and (4.1) hold, at H_0 = 1: phi_1 and a chosen, da/dt from
    (2.25), d^2a/dt^2 from (2.26), and the rate of phi_1 from (4.1); `sign` picks the sign of H."""
    import limiting_curvature as lc
    while True:
        p1 = rng.uniform(0.2, 3.0)
        a0 = rng.uniform(1.2, 4.0)
        first = float(lc.V1p(p1)) / 12 - k / a0 ** 2
        if first > 0.01:
            break
    p2 = rng.uniform(-2.0, 2.0)
    H = sign * math.sqrt(first)
    Hd = k / a0 ** 2 - float(lc.V2p(p2)) / math.sqrt(12)
    p1d = -4 * math.sqrt(3) * H * float(lc.V2p(p2)) / float(lc.V1pp(p1))
    return {"a0": a0, "a1": a0 * H, "a2": a0 * (Hd + H * H), "p10": p1, "p11": p1d, "p20": p2,
            "p21": rng.uniform(-1.0, 1.0), "p22": rng.uniform(-1.0, 1.0), "k": k, "H": H, "Hd": Hd}


@unittest.skipUnless(HAS_ALL, "sympy, numpy and scipy are not all installed")
class TheAction(unittest.TestCase):
    def setUp(self):
        self.rng = random.Random(1993)

    def point(self, k, gap):
        """Any point, with dH/dt - k/a^2 = gap."""
        a0, a1 = self.rng.uniform(1.2, 3.0), self.rng.uniform(-1.5, 1.5)
        H = a1 / a0
        Hd = k / a0 ** 2 + gap
        return {"a0": a0, "a1": a1, "a2": a0 * (Hd + H * H), "a3": self.rng.uniform(-1, 1),
                "a4": self.rng.uniform(-1, 1), "p10": self.rng.uniform(0.1, 2.0), "p11": self.rng.uniform(-1, 1),
                "p12": self.rng.uniform(-1, 1), "p20": self.rng.uniform(-2, 2), "p21": self.rng.uniform(-1, 1),
                "p22": self.rng.uniform(-1, 1), "k": k, "H": H, "Hd": Hd}

    def test_the_signed_root_gives_the_constraints_and_the_positive_root_neither(self):
        import limiting_curvature as lc
        equations = action()["equations"]
        for k in (0, 1):
            for gap in (0.7, -0.4):
                v = self.point(k, gap)
                a0, H, Hd = v["a0"], v["H"], v["Hd"]
                # (2.25) and (2.26), each as the residual of its field's equation divided by a^3.
                first = 12 * (H * H + k / a0 ** 2) - float(lc.V1p(v["p10"]))
                second = -math.sqrt(12) * (Hd - k / a0 ** 2) - float(lc.V2p(v["p20"]))
                self.assertAlmostEqual(at(equations[1]["phi1"], v) / a0 ** 3, first, places=10)
                self.assertAlmostEqual(at(equations[1]["phi2"], v) / a0 ** 3, second, places=10)
                # Where dH/dt - k/a^2 < 0 the positive root is -sqrt(12)(dH/dt - k/a^2), and with it the
                # equations are no longer those two.
                if gap < 0:
                    self.assertGreater(abs(at(equations[-1]["phi1"], v) / a0 ** 3 - first), 1)
                    self.assertGreater(abs(at(equations[-1]["phi2"], v) / a0 ** 3 - second), 1)

    def test_the_lapse_gives_the_hamiltonian_equation_with_two_phi_1_plus_one(self):
        """On (2.25), (2.26) and (4.1) the lapse's equation is a fixed multiple of the residual of
        -(V_1 + V_2)/2 + 3H^2(1 - 2 phi_1) + 3(k/a^2)(2 phi_1 + 1) - sqrt(3)H(dphi_2/dt + 3H phi_2)
        + sqrt(3) k phi_2/a^2; the paper's 4 phi_1 + 1 and 4 phi_1 - 1 agree with it only at k = 0."""
        import limiting_curvature as lc
        equations = action()["equations"]

        def residual(v, c1, c0):
            p1, p2, a0, H = v["p10"], v["p20"], v["a0"], v["H"]
            return (-(float(lc.V1(p1)) + float(lc.V2(p2))) / 2 + 3 * H * H * (1 - 2 * p1)
                    + 3 * (v["k"] / a0 ** 2) * (c1 * p1 + c0) - math.sqrt(3) * H * (v["p21"] + 3 * H * p2)
                    + math.sqrt(3) * v["k"] * p2 / a0 ** 2)
        ratios = []
        for k in (0, 1):
            for sign in (1, -1):
                v = on_shell(self.rng, k, sign)
                lapse = at(equations[1]["N"], v) / v["a0"] ** 3
                ratios.append(lapse / residual(v, 2, 1))
                misprints = [residual(v, 4, 1), residual(v, 4, -1)]
                if k == 0:
                    self.assertAlmostEqual(misprints[0], residual(v, 2, 1), places=12)
                    self.assertAlmostEqual(misprints[1], residual(v, 2, 1), places=12)
                else:
                    self.assertTrue(all(abs(lapse / m - ratios[0]) > 1e-3 for m in misprints), (lapse, misprints))
        self.assertTrue(all(abs(r - ratios[0]) < 1e-10 for r in ratios), ratios)
        self.assertGreater(abs(ratios[0]), 1e-3)

    def test_the_scale_factor_gives_the_equation_the_solver_integrates(self):
        import numpy as np
        import sympy as sp
        import limiting_curvature as lc
        t, k, a, p1, p2 = action()["symbols"]
        e = action()["equations"][1]["a"]
        unknown = sp.Symbol("unknown")
        e = e.subs(sp.Derivative(p2, (t, 2)), unknown)
        # The coefficient of d^2phi_2/dt^2 is -2 sqrt(3) a^2, so the equation divides by nothing.
        coefficient = sp.diff(e, unknown)
        for sign in (1, -1):
            v = on_shell(self.rng, 1, sign)
            v.pop("p22")
            self.assertAlmostEqual(at(coefficient, v), -2 * math.sqrt(3) * v["a0"] ** 2, places=10)
            solved = -at(e.subs(unknown, 0), v) / at(coefficient, v)
            state = [v["a0"], v["H"], v["p10"], v["p20"], v["p21"], 0.0]
            rates = lc.rates(0.0, np.array(state))
            self.assertAlmostEqual(rates[4], solved, places=9)
            self.assertAlmostEqual(rates[1], v["Hd"], places=12)
            self.assertAlmostEqual(rates[2], v["p11"], places=12)


@unittest.skipUnless(HAS_ALL, "sympy, numpy and scipy are not all installed")
class ThePotentials(unittest.TestCase):
    def test_the_slope_of_v_1_rises_from_zero_toward_twelve_and_never_reaches_it(self):
        import numpy as np
        import sympy as sp
        import limiting_curvature as lc
        x = sp.Symbol("x", nonnegative=True)
        V1 = 12 * x ** 2 / (1 + x) * (1 - sp.log(1 + x) / (1 + x))
        slope = 12 * x * (x ** 2 + 2 * x + 2 - 2 * sp.log(1 + x)) / (1 + x) ** 3
        self.assertEqual(sp.cancel(sp.expand(sp.diff(V1, x) - slope)), 0)
        grid = np.unique(np.concatenate([np.linspace(0, 10, 20001), np.geomspace(10, 1e9, 2001)]))
        values = lc.V1p(grid) / 12
        self.assertEqual(values[0], 0)
        # Far out the slope sits within a rounding of 1, so it is held strictly rising up to phi_1 = 10^4.
        self.assertTrue(np.all(np.diff(values[grid <= 1e4]) > 0))
        self.assertTrue(np.all(np.diff(values) >= 0))
        self.assertTrue(np.all(values < 1))
        self.assertLess(1 - values[-1], 1e-8)
        self.assertTrue(np.all(lc.V1pp(grid[grid < 1e6]) > 0))

    def test_a_closed_universe_is_never_smaller_than_c_over_h_0(self):
        # (2.25) with k = 1: 1/a^2 = V_1'/12 - H^2 < 1.
        import numpy as np
        import limiting_curvature as lc
        self.assertTrue(np.all(lc.V1p(np.geomspace(1e-6, 1e9, 4001)) / 12 < 1))
        with self.assertRaises(ValueError):
            lc.bounce_data(1.0)

    def test_the_slope_of_v_2_bounds_dh_dt_minus_k_over_a_squared(self):
        import numpy as np
        import sympy as sp
        import limiting_curvature as lc
        x = sp.Symbol("x", real=True)
        V2 = 2 * sp.sqrt(3) * x ** 2 / (1 + x ** 2)
        slope = sp.diff(V2, x)
        self.assertEqual(sp.solve(sp.diff(slope, x), x), [-1 / sp.sqrt(3), 1 / sp.sqrt(3)])
        top = slope.subs(x, 1 / sp.sqrt(3))
        self.assertEqual(sp.nsimplify(top), sp.Rational(9, 4))
        self.assertEqual(sp.nsimplify(top / sp.sqrt(12) - 3 * sp.sqrt(3) / 8), 0)
        grid = np.linspace(-50, 50, 200001)
        self.assertLessEqual(np.max(np.abs(lc.V2p(grid))), 9 / 4 + 1e-12)
        # As |phi_2| grows the slope falls as 4 sqrt(3)/phi_2^3, and dH/dt - k/a^2 with it.
        far = np.array([1e2, 1e3, -1e4])
        self.assertTrue(np.allclose(lc.V2p(far) * far ** 3, 4 * math.sqrt(3), rtol=1e-3))

    def test_the_bounce_is_a_regular_point_of_the_hamiltonian_equation(self):
        """At H = 0 the Hamiltonian equation, divided by sqrt(3) H, has the coefficient lambda = 1 in
        front of phi_2/(t - t_b), so every solution crosses the bounce with phi_2 finite."""
        import sympy as sp
        p1, p2, a, k = sp.symbols("phi1 phi2 a k", positive=True)
        V1, V2 = sp.Function("V1"), sp.Function("V2")
        numerator = -(V1(p1) + V2(p2)) / 2 + 3 * (k / a ** 2) * (2 * p1 + 1) + sp.sqrt(3) * k * p2 / a ** 2
        rate = k / a ** 2 - sp.diff(V2(p2), p2) / sp.sqrt(12)
        self.assertEqual(sp.cancel(sp.diff(numerator, p2) / (sp.sqrt(3) * rate) - 1), 0)


@unittest.skipUnless(HAS_ALL, "sympy, numpy and scipy are not all installed")
class TheBounce(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import numpy as np
        import limiting_curvature as lc
        cls.np, cls.lc = np, lc
        cls.universe = lc.bounce()
        cls.t = np.linspace(-lc.T_FAR, lc.T_FAR, 24001)
        cls.jet = cls.universe.jet(cls.t)

    def test_the_constraints_hold_along_the_whole_solution(self):
        self.assertLess(self.universe.residual(self.t), 1e-9)

    def test_the_hubble_rate_changes_sign_at_a_finite_radius(self):
        H, a = self.jet["H"], self.jet["a"]
        middle = len(self.t) // 2
        self.assertEqual(self.t[middle], 0)
        self.assertEqual(H[middle], 0)
        self.assertTrue(self.np.all(H[:middle] < 0) and self.np.all(H[middle + 1:] > 0))
        self.assertTrue(self.np.all(self.np.diff(H) > 0))
        self.assertEqual(a[middle], 2.0)
        self.assertTrue(self.np.all(a >= 2.0))
        self.assertGreater(self.jet["d2a"][middle], 0)

    def test_the_bounds_hold_at_every_time(self):
        a, H, Hd = self.jet["a"], self.jet["H"], self.jet["dH"]
        self.assertTrue(self.np.all(H ** 2 + 1 / a ** 2 < 1))
        self.assertTrue(self.np.all(a > 1))
        self.assertTrue(self.np.all(self.np.abs(Hd - 1 / a ** 2) <= 3 * math.sqrt(3) / 8))

    def test_the_curvature_is_finite_and_the_universe_approaches_de_sitter_space(self):
        a, da, d2a, Hd = self.jet["a"], self.jet["da"], self.jet["d2a"], self.jet["dH"]
        K = 12 * ((d2a / a) ** 2 + ((da ** 2 + 1) / a ** 2) ** 2)
        self.assertTrue(self.np.all(self.np.isfinite(K)))
        self.assertLess(float(self.np.max(K)), 24)
        # dH/dt - k/a^2 = -V_2'/sqrt(12), and with it I_2 = 12(dH/dt - k/a^2)^2, falls toward zero.
        gap = self.np.abs(Hd - 1 / a ** 2)
        at = {t: float(gap[self.np.argmin(self.np.abs(self.t - t))]) for t in (2.0, 10.0, self.lc.T_FAR)}
        self.assertTrue(0.006 < at[2.0] < 0.009 and 0.002 < at[10.0] < 0.004 and 4e-4 < at[self.lc.T_FAR] < 6e-4, at)
        self.assertTrue(self.np.all(self.np.diff(gap[self.t > 3]) < 0))
        self.assertTrue(self.np.all(self.jet["phi2"][self.t > 3] < -7))

    def test_the_symmetry_in_time(self):
        np = self.np
        for key in ("a", "phi1", "phi2", "d2a"):
            self.assertTrue(np.allclose(self.jet[key], self.jet[key][::-1], rtol=1e-12, atol=1e-12), key)
        for key in ("H", "da", "eta"):
            self.assertTrue(np.allclose(self.jet[key], -self.jet[key][::-1], rtol=1e-12, atol=1e-12), key)

    def test_the_conformal_time_spans_less_than_pi(self):
        from scipy.integrate import quad
        end = self.universe.eta_end
        self.assertAlmostEqual(end, 0.8430, places=4)
        self.assertLess(2 * end, math.pi)
        # The conformal time again, by quadrature of 1/a, which the solver does not use.
        value, _ = quad(lambda s: 1 / float(self.universe.a(s)), 0, 2.0, epsabs=1e-12, epsrel=1e-12, limit=200)
        self.assertAlmostEqual(value, float(self.universe.eta(2.0)), places=9)


@unittest.skipUnless(HAS_ALL, "sympy, numpy and scipy are not all installed")
class ThePublishedCharts(unittest.TestCase):
    def setUp(self):
        import sympy
        import verify_metrics
        self.sp, self.vm = sympy, verify_metrics
        self.metric = json.loads((DATA / "metrics" / "limiting_curvature.json").read_text(encoding="utf-8"))

    def invariants(self, system):
        """R and I_2 = 4 R^mu_nu R^nu_mu - R^2 from the chart's published mixed Ricci tensor."""
        entry = next(c for c in self.metric["coordinates"] if c["id"] == system)
        reader = self.vm.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]], ())
        mixed = self.sp.zeros(4, 4)
        names = entry["coords"]
        for e in entry["ricci_tensor"]["variants"]["ul"]["nonzero"]:
            i, j = (names.index(n) for n in e["indices"])
            mixed[i, j] = reader(e["value"])
        R = reader(entry["ricci_scalar"].partition("=")[2])
        square = sum(mixed[i, j] * mixed[j, i] for i in range(4) for j in range(4))
        return R, 4 * square - R ** 2, reader

    def test_the_comoving_chart_has_i_1_only_with_the_signed_root(self):
        sp = self.sp
        R, I2, reader = self.invariants("comoving")
        a = reader.parameters["a"]
        x0 = reader.symbol["t"]
        rng = random.Random(1992)
        signs = set()
        for k in range(8):
            A, A1 = sp.Rational(rng.randint(50, 300), 100), sp.Rational(rng.randint(-200, 200), 100)
            A2 = (A1 ** 2 + 1) / A + sp.Rational(rng.randint(10, 200), 100) * (1 if k % 2 else -1)
            values = [(sp.Derivative(a, (x0, 2)), A2), (sp.Derivative(a, x0), A1), (a, A)]
            r, i2 = (float(sp.N(e.subs(values), 30)) for e in (R, I2))
            D = float(A2 / A - (A1 ** 2 + 1) / A ** 2)
            signs.add(D > 0)
            first = float(12 * (A1 ** 2 + 1) / A ** 2)
            self.assertAlmostEqual(i2, 12 * D * D, places=10)
            self.assertAlmostEqual(r - math.sqrt(3) * math.sqrt(12) * D, first, places=10)
            positive = r - math.sqrt(3) * math.sqrt(i2)
            if D > 0:
                self.assertAlmostEqual(positive, first, places=10)
            else:
                self.assertAlmostEqual(positive - first, 12 * D, places=10)
        self.assertEqual(signs, {True, False})

    def test_the_de_sitter_bounce_has_i_2_zero_and_the_limit_on_its_radius(self):
        sp = self.sp
        R, I2, reader = self.invariants("de_sitter_bounce")
        ell, x0 = reader.parameters["ell"], reader.symbol["t"]
        for at in ((1, 0), (sp.Rational(3, 2), sp.Rational(7, 10)), (2, -3)):
            values = {ell: at[0], x0: at[1]}
            self.assertLess(abs(float(sp.N(I2.subs(values), 30))), 1e-25)
            self.assertLess(abs(float(sp.N((R - 12 / ell ** 2).subs(values), 30))), 1e-25)
        # a = l cosh(ct/l) keeps H^2 + c^2/a^2 = c^2/l^2, below H_0^2 only for l > c/H_0.
        entry = next(c for c in self.metric["coordinates"] if c["id"] == "de_sitter_bounce")
        self.assertIn("\\ell > c/H_0", entry["parameters"][0]["description"])
        a = ell * sp.cosh(x0 / ell)
        first = sp.diff(a, x0) ** 2 / a ** 2 + 1 / a ** 2 - 1 / ell ** 2
        for at in ((1, 0), (sp.Rational(3, 2), sp.Rational(7, 10)), (2, -3)):
            self.assertLess(abs(float(sp.N(first.subs({ell: at[0], x0: at[1]}), 30))), 1e-25)

    def test_the_embedding_draws_the_solver_s_radius_and_none_below_c_over_h_0(self):
        import limiting_curvature as lc
        view = json.loads((DATA / "embedding" / "limiting_curvature.json").read_text(encoding="utf-8"))["views"][0]
        universe = lc.bounce()
        radii = []
        for frame, value in zip(view["movie"]["frames"], [f["value"] for f in view["movie"]["frames"]]):
            equator = [r for r in frame["rings"] if r["class"] == "chartedge"]
            self.assertEqual(len(equator), 1)
            radii.append(equator[0]["rho"])
            self.assertAlmostEqual(equator[0]["rho"], float(universe.a(value)), places=3)
        self.assertEqual(min(radii), 2.0)
        self.assertGreater(len(radii), 30)


if __name__ == "__main__":
    unittest.main()
