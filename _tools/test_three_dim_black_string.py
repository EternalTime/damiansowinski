"""The charged black string of three dimensions:
.venv.noindex/bin/python -m unittest _tools/test_three_dim_black_string.py

What its charts and texts state, held on the published files. The static chart is held to Horne and
Horowitz's line element as the History prints it, and its published Ricci tensor, Christoffel symbols and
Ricci scalar to the three field equations of the page's convention, with the dilaton
Phi = ln r + (1/2) ln(k/2) and H_rtx = Q/r^2: the metric equation at every k and in both of its forms, the
dilaton's equation with 8/k and with no other constant, and the field's. The two horizons are held to
r = M and r = Q^2/M, the singularity to being timelike, every chart to the static one pulled back, the
temperature and the period of x to the History's formulas, and each chart's own check to failing on a
metric that is not the solution. The tests need sympy and are skipped where it is absent.
"""
import importlib.util
import json
import math
import re
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "MFS" / "assets" / "data"
HAS_SYMPY = importlib.util.find_spec("sympy") is not None
sys.path.insert(0, str(Path(__file__).resolve().parent / "derivations"))
ID = "three_dim_black_string"
METRIC = json.loads((DATA / "metrics" / f"{ID}.json").read_text(encoding="utf-8"))
M, Q = 2.0, 1.0                      # Q^2/M = 1/2: the inner horizon at r = 1/2, the event horizon at r = 2
LEVEL = 52 / 23                      # k for the string alone
RADII = (0.2, 0.4, 0.7, 1.3, 1.9, 2.4, 5.0)


def chart(system):
    return next(c for c in METRIC["coordinates"] if c["id"] == system)


@unittest.skipUnless(HAS_SYMPY, "needs sympy")
class FieldEquations(unittest.TestCase):
    """The published static chart against the equations of the page's convention."""

    @classmethod
    def setUpClass(cls):
        import sympy as sp
        import verify_metrics as vm
        cls.sp = sp
        c = chart("static")
        cls.R = R = vm.Reader(c["coords"], [p["symbol"] for p in c["parameters"]], ())
        cls.t, cls.x, cls.r = (R.symbol[name] for name in c["coords"])
        cls.M, cls.Q, cls.k = (R.parameters[name] for name in ("M", "Q", "k"))
        cls.coords = c["coords"]
        cls.g = cls.matrix(c["metric_components"])
        cls.ginv = cls.matrix(c["inverse_metric_components"])
        cls.ricci = cls.matrix(c["ricci_tensor"]["variants"]["ll"]["nonzero"])
        cls.scalar = R(c["ricci_scalar"].split("=", 1)[1])
        cls.kretschmann = R(c["kretschmann"].split("=", 1)[1])
        cls.gamma = [[[0] * 3 for _ in range(3)] for _ in range(3)]
        for entry in c["christoffel"]["variants"]["ull"]["nonzero"]:
            a, b, d = (c["coords"].index(name) for name in entry["indices"])
            cls.gamma[a][b][d] = R(entry["value"])
        cls.Phi = sp.log(cls.r) + sp.log(cls.k / 2) / 2

    @classmethod
    def matrix(cls, entries):
        g = cls.sp.zeros(3)
        for entry in entries:
            i, j = (cls.coords.index(name) for name in entry["indices"])
            g[i, j] = cls.R(entry["value"])
        return g

    def point(self, r, k=LEVEL, mass=M, charge=Q):
        exact = lambda v: self.sp.nsimplify(v, rational=True)
        return {self.t: exact(0.3), self.x: exact(0.7), self.r: exact(r), self.M: exact(mass),
                self.Q: exact(charge), self.k: exact(k)}

    def value(self, expression, at):
        return complex(self.sp.N(self.sp.sympify(expression).subs(at), 30))

    def residuals(self, r, k=LEVEL, dilaton=None, constant=None, mass=M, charge=Q):
        """The metric equation R_ab - nabla_a nabla_b Phi - HH/4 at a point, and the dilaton's equation
        R - 2 nabla^2 Phi - (nabla Phi)^2 - H^2/12 + constant, from the published tensors."""
        sp = self.sp
        at = self.point(r, k, mass, charge)
        Phi = self.Phi if dilaton is None else dilaton
        constant = 8 / self.k if constant is None else constant
        s = (self.t, self.x, self.r)
        dPhi = [sp.diff(Phi, a) for a in s]
        hess = sp.Matrix(3, 3, lambda m, n: sp.diff(Phi, s[m], s[n])
                         - sum(self.gamma[c][m][n] * dPhi[c] for c in range(3)))
        g, gi = self.g.subs(at), self.ginv.subs(at)
        H = [[[0] * 3 for _ in range(3)] for _ in range(3)]
        # H_rtx = Q/r^2, antisymmetric in its three indices.
        for (a, b, c), sign in (((2, 0, 1), 1), ((0, 1, 2), 1), ((1, 2, 0), 1),
                                ((2, 1, 0), -1), ((1, 0, 2), -1), ((0, 2, 1), -1)):
            H[a][b][c] = sign * self.Q / self.r ** 2
        H = [[[sp.sympify(H[a][b][c]).subs(at) for c in range(3)] for b in range(3)] for a in range(3)]
        hh = sp.Matrix(3, 3, lambda m, n: sum(H[m][a][b] * H[n][c][d] * gi[a, c] * gi[b, d]
                                                for a in range(3) for b in range(3)
                                                for c in range(3) for d in range(3)))
        H2 = sum(H[a][b][c] * H[d][e][f] * gi[a, d] * gi[b, e] * gi[c, f]
                 for a in range(3) for b in range(3) for c in range(3)
                 for d in range(3) for e in range(3) for f in range(3))
        hess_at = hess.subs(at)
        grad = sp.Matrix([d.subs(at) for d in dPhi])
        metric = self.ricci.subs(at) - hess_at - hh / 4
        box = sum(gi[m, n] * hess_at[m, n] for m in range(3) for n in range(3))
        scalar = (self.scalar.subs(at) - 2 * box - (grad.T * gi * grad)[0, 0] - H2 / 12
                  + sp.sympify(constant).subs(at))
        return [self.value(x, {}) for x in metric], self.value(scalar, {})

    def test_the_metric_equation_vanishes_at_every_k(self):
        for k in (LEVEL, 0.5, 3.0, 40.0):
            for r in RADII:
                metric, _ = self.residuals(r, k)
                self.assertLess(max(abs(x) for x in metric), 1e-20, msg=f"k = {k}, r = {r}")

    def test_the_metric_equation_holds_for_phi_minus_half_Phi_and_not_plus(self):
        # R_ab + 2 nabla_a nabla_b phi - HH/4 = 0 is the same equation for phi = -Phi/2; phi = +Phi/2
        # turns the Hessian's sign, which no longer balances.
        for r in RADII:
            metric, _ = self.residuals(r)
            self.assertLess(max(abs(x) for x in metric), 1e-20)
            flipped, _ = self.residuals(r, dilaton=-self.Phi)
            self.assertGreater(max(abs(x) for x in flipped), 1e-3, msg=f"r = {r}")

    def test_the_dilatons_equation_needs_8_over_k(self):
        for k in (LEVEL, 0.5, 3.0):
            for r in RADII:
                _, scalar = self.residuals(r, k)
                self.assertLess(abs(scalar), 1e-20, msg=f"k = {k}, r = {r}")
                for wrong in (0, 4 / self.k, 9 / self.k):
                    _, off = self.residuals(r, k, constant=wrong)
                    self.assertGreater(abs(off), 1e-3)

    def test_e_to_the_Phi_times_the_field_is_constant_so_its_equation_holds(self):
        # sqrt(-g) e^Phi H^rtx is e^Phi Q sqrt(8/k)/r up to sign, so nabla_a(e^Phi H^abc) = 0.
        sp = self.sp
        root = sp.sqrt(-self.g.det())
        carried = sp.exp(self.Phi) * (self.Q / self.r ** 2) / root
        for r in RADII:
            at = self.point(r)
            self.assertAlmostEqual(abs(self.value(carried, at)), 2 * Q, places=12)
            self.assertLess(abs(self.value(sp.diff(carried, self.r), at)), 1e-20)

    def test_the_ricci_scalar_and_kretschmann_scalar_are_horne_and_horowitzs(self):
        for r in RADII:
            at = self.point(r)
            R = 4 * (2 * M * M * r - 7 * M * Q * Q + 2 * Q * Q * r) / (M * LEVEL * r * r)
            K = 16 * (4 * M ** 4 * r * r - 12 * M ** 3 * Q * Q * r + 19 * M * M * Q ** 4
                      - 12 * M * Q ** 4 * r + 4 * Q ** 4 * r * r) / (M * M * LEVEL ** 2 * r ** 4)
            self.assertAlmostEqual(self.value(self.scalar, at).real, R, places=10)
            self.assertAlmostEqual(self.value(self.kretschmann, at).real, K, places=8)

    def test_the_horizons_are_at_M_and_Q_squared_over_M(self):
        grr = self.ginv[2, 2]
        for horizon in (M, Q * Q / M):
            self.assertLess(abs(self.value(grr, self.point(horizon))), 1e-25)
        roots = self.sp.solve(self.sp.numer(self.sp.together(grr)), self.r)
        self.assertEqual(set(roots), {self.M, self.Q ** 2 / self.M})

    def test_the_singularity_is_timelike_and_x_is_time_inside_the_inner_horizon(self):
        for r in (0.05, 0.2, 0.45):
            at = self.point(r)
            self.assertGreater(self.value(self.ginv[2, 2], at).real, 0)     # r spacelike, so r = 0 is timelike
            self.assertLess(self.value(self.g[1, 1], at).real, 0)           # d/dx timelike
            self.assertGreater(self.value(self.g[0, 0], at).real, 0)        # d/dt spacelike
        for r in (0.7, 1.9):
            self.assertLess(self.value(self.ginv[2, 2], self.point(r)).real, 0)   # r is time between the horizons
        for r in (0.01, 0.001):
            self.assertGreater(abs(self.value(self.kretschmann, self.point(r))), 1e6)

    def test_the_temperature_is_the_surface_gravity_over_two_pi(self):
        # kappa^2 = -(1/2) nabla^a xi^b nabla_a xi_b on the horizon for xi = d/dt, which for this diagonal
        # metric is g^rr (d_r g_tt)^2 / (-4 g_tt), taken in the limit r -> M.
        sp = self.sp
        expr = sp.cancel(self.ginv[2, 2] * sp.diff(self.g[0, 0], self.r) ** 2 / (-4 * self.g[0, 0]))
        for k, mass, charge in ((LEVEL, M, Q), (3.0, 3.0, 2.0)):
            kappa2 = expr.subs({self.r: mass, self.M: mass, self.Q: charge, self.k: k})
            T = math.sqrt((mass ** 2 - charge ** 2) / (2 * k)) / (math.pi * mass)
            self.assertAlmostEqual(math.sqrt(float(kappa2)) / (2 * math.pi), T, places=12)
        self.assertIn("$T = \\sqrt{(M^2 - Q^2)/2k}\\,/\\pi M$", METRIC["history"])


@unittest.skipUnless(HAS_SYMPY, "needs sympy")
class Charts(unittest.TestCase):
    """Every published chart against the static chart pulled back, and each check against a wrong metric."""

    @classmethod
    def setUpClass(cls):
        import sympy as sp
        import verify_metrics as vm
        import three_dim_black_string as module
        cls.sp, cls.vm, cls.module = sp, vm, module

    def published(self, system, values, point):
        sp = self.sp
        c = chart(system)
        R = self.vm.Reader(c["coords"], [p["symbol"] for p in c["parameters"]], ())
        at = {R.symbol[name]: v for name, v in zip(c["coords"], point)}
        at.update({R.parameters[name]: v for name, v in values.items() if name in R.parameters})
        g = sp.zeros(3)
        for entry in c["metric_components"]:
            i, j = (c["coords"].index(name) for name in entry["indices"])
            g[i, j] = sp.N(R(entry["value"]).subs(at), 30)
        return g

    def test_every_chart_is_the_static_chart_pulled_back(self):
        sp, module = self.sp, self.module
        for system in module.SYSTEMS:
            for (mass, charge, k), point in module.POINTS[system]:
                mass, charge, k = (sp.sympify(v) for v in (mass, charge, k))
                point = [sp.sympify(v) for v in point]
                g = self.published(system, {"M": mass, "Q": charge, "k": k}, point)
                r = module.radius(system, point, mass, charge)
                J = module.jacobian(system, mass, charge, k, *point)
                pulled = J.T * module.static_metric(mass, charge, k, r) * J
                worst = max(abs(complex(sp.N(x, 30))) for x in g - pulled)
                self.assertLess(worst, 1e-20, msg=f"{system} at {point}")

    def test_the_superextremal_x_closes_without_a_cone_at_the_historys_period(self):
        # Near rho = 0 the circle of constant rho has circumference P sqrt(g_xx) and radius
        # integral sqrt(g_rho rho) d rho; no cone needs P d sqrt(g_xx)/d rho = 2 pi sqrt(g_rho rho) at rho = 0.
        sp = self.sp
        for mass, charge, k in ((1.0, 2.0, LEVEL), (2.0, 3.0, 5.0), (1.0, -1.5, 0.7)):
            eps = 1e-8
            g0 = self.published("superextremal", {"M": mass, "Q": charge, "k": k}, (0, 0, eps))
            slope = math.sqrt(float(g0[1, 1])) / eps
            period = 2 * math.pi * math.sqrt(float(g0[2, 2])) / slope
            self.assertAlmostEqual(period, math.pi * abs(charge) * math.sqrt(2 * k / (charge ** 2 - mass ** 2)),
                                   places=5)
        self.assertIn("period $\\pi|Q|\\sqrt{2k/(Q^2 - M^2)}$", METRIC["history"])

    def test_each_charts_check_fails_on_a_metric_that_is_not_the_solution(self):
        import chart_printer as cp
        module = self.module
        for system in module.SYSTEMS:
            line = module.LINES[system].replace("k", "3k", 1)
            wrong = cp.Chart(module.COORDS[system], module.PARAMETERS[system], line)
            with self.assertRaises(AssertionError, msg=system):
                module.check(wrong, system)


@unittest.skipUnless(HAS_SYMPY, "needs sympy")
class Texts(unittest.TestCase):

    def test_the_history_prints_the_static_chart(self):
        import sympy as sp
        import verify_metrics as vm
        shown = re.search(r"\$\$(ds\^2 = .*?),\$\$", METRIC["history"]).group(1)
        c = chart("static")
        R = vm.Reader(c["coords"], [p["symbol"] for p in c["parameters"]], ())
        difference = R(shown.split("=", 1)[1]) - R(c["line_element"].split("=", 1)[1])
        for r in RADII:
            at = {R.symbol["r"]: r, R.parameters["M"]: M, R.parameters["Q"]: Q, R.parameters["k"]: LEVEL}
            self.assertLess(abs(complex(sp.N(sp.expand(difference).subs(at), 30))), 1e-20)

    def test_the_history_states_the_horizons_and_the_dilaton(self):
        history = METRIC["history"]
        for text in ("$r = M$", "$r = Q^2/M$", "$\\Phi = \\ln r + \\tfrac{1}{2}\\ln(k/2)$", "$H_{rtx} = Q/r^2$"):
            self.assertIn(text, history)

    def test_the_convention_states_the_field_equation_and_the_contraction(self):
        convention = METRIC["convention"]
        for text in ("\\alpha' = 1/2", "R_{\\mu\\nu} - \\nabla_\\mu\\nabla_\\nu\\Phi - \\tfrac{1}{4}H_{\\mu\\alpha\\beta}H_\\nu{}^{\\alpha\\beta} = 0",
                     "\\phi = -\\Phi/2", "8/k", "R_{\\mu\\nu} = R^\\alpha{}_{\\mu\\alpha\\nu}"):
            self.assertIn(text, convention)


if __name__ == "__main__":
    unittest.main()
