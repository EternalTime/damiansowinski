"""The field equations, the curvature and the two frames of the pre-big bang dilaton vacuum.

The geometry is built here from the line element alone, with t a negative symbol so that every
power of -t is a power of a positive number, and held to what the History and the two charts
state: both string frame equations hold on the sphere sum_i beta_i^2 = 1 and fail off it, the
curvature grows without bound as t rises to zero for every exponents but the flat ones, the
shifted dilaton is the same for every exponents and under the inversion of the scale factors,
and the Einstein frame is Kasner's universe with a massless scalar field.
"""

import json
from pathlib import Path
import sys
import unittest

import sympy as sp

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "_tools" / "derivations"))
import verify_metrics as vm  # noqa: E402

t = sp.Symbol("t", negative=True)
SPACE = sp.symbols("x y z", real=True)
X = [t, *SPACE]
b1, b2, b3 = BETAS = sp.symbols("beta_1 beta_2 beta_3", real=True)
ROOT3 = sp.sqrt(3)

# Points of the sphere sum_i beta_i^2 = 1: two anisotropic, the expanding and the contracting
# isotropic members, one with exponents of either sign, and (-1, 0, 0), the inversion of the
# flat member (1, 0, 0), which is curved.
ON = [(sp.Rational(2, 3), sp.Rational(2, 3), sp.Rational(1, 3)),
      (sp.Rational(3, 5), sp.Rational(4, 5), 0),
      (-1 / ROOT3,) * 3,
      (1 / ROOT3,) * 3,
      (-sp.Rational(2, 3), sp.Rational(1, 3), -sp.Rational(2, 3)),
      (-sp.Integer(1), sp.Integer(0), sp.Integer(0))]
FLAT = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
OFF = (sp.Rational(1, 2), sp.Rational(1, 2), sp.Rational(1, 2))


def geometry(g):
    """Christoffel symbols, Ricci tensor, Ricci scalar and Kretschmann scalar of a diagonal metric."""
    n = len(X)
    gi = sp.diag(*[1 / g[a, a] for a in range(n)])

    def clean(e):
        return sp.simplify(sp.powsimp(sp.together(e), force=True))
    gamma = [[[clean(gi[a, a] * (sp.diff(g[a, b], X[c]) + sp.diff(g[a, c], X[b]) - sp.diff(g[b, c], X[a])) / 2)
               for c in range(n)] for b in range(n)] for a in range(n)]

    def riemann(a, b, c, d):
        r = sp.diff(gamma[a][b][d], X[c]) - sp.diff(gamma[a][b][c], X[d])
        r += sum(gamma[a][c][e] * gamma[e][b][d] - gamma[a][d][e] * gamma[e][b][c] for e in range(n))
        return clean(r)
    R = {(a, b, c, d): riemann(a, b, c, d) for a in range(n) for b in range(n) for c in range(n) for d in range(n)}
    ricci = sp.Matrix(n, n, lambda b, d: clean(sum(R[a, b, a, d] for a in range(n))))
    scalar = clean(sum(gi[a, a] * ricci[a, a] for a in range(n)))
    kretschmann = clean(sum(g[a, a] * gi[b, b] * gi[c, c] * gi[d, d] * R[a, b, c, d] ** 2 for (a, b, c, d) in R))
    return gi, gamma, ricci, scalar, kretschmann


def string_frame(betas):
    return sp.diag(-1, *[(-t) ** (2 * b) for b in betas])


def dilaton(betas):
    return (sum(betas) - 1) * sp.log(-t)


def equations(betas):
    """R_mu nu + nabla_mu nabla_nu phi and R + 2 box phi - (nabla phi)^2 of the line element."""
    g = string_frame(betas)
    gi, gamma, ricci, scalar, _ = geometry(g)
    phi = dilaton(betas)
    grad = [sp.diff(phi, x) for x in X]
    hessian = sp.Matrix(4, 4, lambda m, k: sp.diff(phi, X[m], X[k]) - sum(gamma[c][m][k] * grad[c] for c in range(4)))
    box = sum(gi[m, m] * hessian[m, m] for m in range(4))
    square = sum(gi[m, m] * grad[m] ** 2 for m in range(4))
    return (ricci + hessian).applyfunc(sp.simplify), sp.simplify(scalar + 2 * box - square)


class PreBigBang(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        page = json.loads((ROOT / "MFS/assets/data/metrics/pre_big_bang.json").read_text())
        cls.page = page
        cls.charts = {chart["id"]: chart for chart in page["coordinates"]}
        cls.free = geometry(string_frame(BETAS))

    def published(self, system, text):
        chart = self.charts[system]
        reader = vm.Reader(chart["coords"], [p["symbol"] for p in chart["parameters"]], ())
        return reader, reader(text.split("=", 1)[-1])

    def test_both_field_equations_hold_on_the_sphere_at_every_sampled_exponents(self):
        for betas in ON:
            self.assertEqual(sp.simplify(sum(b ** 2 for b in betas)), 1, betas)
            first, second = equations(betas)
            self.assertEqual(first, sp.zeros(4, 4), betas)
            self.assertEqual(second, 0, betas)

    def test_both_field_equations_fail_off_the_sphere(self):
        first, second = equations(OFF)
        self.assertEqual(sp.simplify(first[0, 0] * t ** 2), sp.Rational(1, 4))
        self.assertEqual(sp.simplify(second * t ** 2), -sp.Rational(1, 4))

    def test_the_sphere_is_exactly_what_the_equations_ask(self):
        # For free exponents the space equations hold and the other two are -+(sum beta_i^2 - 1)/t^2.
        first, second = equations(BETAS)
        sphere = sum(b ** 2 for b in BETAS) - 1
        self.assertEqual(sp.simplify(first[0, 0] + sphere / t ** 2), 0)
        for m in range(4):
            for k in range(4):
                if (m, k) != (0, 0):
                    self.assertEqual(sp.simplify(first[m, k]), 0, (m, k))
        self.assertEqual(sp.simplify(second - sphere / t ** 2), 0)

    def test_curvature_grows_without_bound_as_t_rises_to_zero(self):
        _, kretschmann = self.published("synchronous", self.charts["synchronous"]["kretschmann"])
        _, _, _, _, free = self.free
        ours = {"t": t, **{b.name: b for b in BETAS}}
        kretschmann = kretschmann.subs({s: ours[s.name] for s in kretschmann.free_symbols})
        self.assertEqual(sp.simplify(kretschmann - free), 0)
        for betas in ON:
            K = sp.simplify(free.subs(dict(zip(BETAS, betas))))
            strength = sp.simplify(K * t ** 4)
            self.assertTrue(strength.is_positive, betas)
            values = [K.subs(t, -sp.Rational(1, 10) ** n) for n in range(1, 7)]
            self.assertTrue(all(later > 100 * earlier for earlier, later in zip(values, values[1:])), betas)
            self.assertGreater(values[-1], 10 ** 23)
        for betas in FLAT:
            self.assertEqual(sp.simplify(free.subs(dict(zip(BETAS, betas)))), 0, betas)

    def test_ricci_scalar_and_tensor_are_those_of_the_dilaton(self):
        _, _, ricci, scalar, _ = self.free
        for betas in ON:
            at = dict(zip(BETAS, betas))
            B = sum(betas)
            self.assertEqual(sp.simplify(scalar.subs(at) - (B - 1) ** 2 / t ** 2), 0, betas)
            self.assertEqual(sp.simplify(ricci[0, 0].subs(at) - (B - 1) / t ** 2), 0, betas)

    def test_string_coupling_and_expansion_rates(self):
        # g_s^2 = e^phi vanishes, settles or grows as t rises to zero as sum beta_i is above, at or below one,
        # and each H_i = beta_i/t.
        for B, limit in ((2, 0), (1, 1), (sp.Rational(1, 2), sp.oo)):
            self.assertEqual(sp.limit(sp.exp((B - 1) * sp.log(-t)), t, 0, "-"), limit)
        for i, b in enumerate(BETAS, 1):
            a = (-t) ** b
            self.assertEqual(sp.simplify(sp.diff(a, t) / a - b / t), 0)

    def test_the_expanding_isotropic_member_superinflates(self):
        a = (-t) ** (-1 / ROOT3)
        H = sp.diff(a, t) / a
        for value in (sp.diff(a, t), sp.diff(a, t, 2), sp.diff(H, t)):
            for at in (-sp.Rational(1, 100), -1, -100):
                self.assertGreater(value.subs(t, at).evalf(), 0)
        self.assertEqual(sp.simplify(dilaton((-1 / ROOT3,) * 3) + (1 + ROOT3) * sp.log(-t)), 0)

    def test_shifted_dilaton_is_the_same_for_every_exponents_and_under_inversion(self):
        logs = sum(b * sp.log(-t) for b in BETAS)
        shifted = sp.expand(dilaton(BETAS) - logs)
        self.assertEqual(sp.simplify(shifted + sp.log(-t)), 0)
        # a_i -> 1/a_i is beta_i -> -beta_i, with phi -> phi - 2 sum_i ln a_i.
        inverted = dilaton(BETAS) - 2 * logs
        self.assertEqual(sp.simplify(inverted - dilaton([-b for b in BETAS])), 0)
        self.assertEqual(sp.simplify(sp.expand(inverted + logs) - shifted), 0)
        for betas in ON[:2]:
            first, second = equations([-b for b in betas])
            self.assertEqual((first, second), (sp.zeros(4, 4), 0), betas)

    def test_einstein_frame_is_kasner_with_a_massless_scalar_field(self):
        for betas in ON:
            B = sum(betas)
            phi = dilaton(betas)
            g = string_frame(betas) * (-t) ** (1 - B)
            _, _, ricci, _, _ = geometry(g)
            grad = [sp.diff(phi, x) for x in X]
            # In the Einstein frame the dilaton is a massless scalar, R_mu nu = (1/2) d_mu phi d_nu phi.
            stress = sp.Matrix(4, 4, lambda m, k: grad[m] * grad[k] / 2)
            self.assertEqual((ricci - stress).applyfunc(sp.simplify), sp.zeros(4, 4), betas)
            p = [(1 - B + 2 * b) / (3 - B) for b in betas]
            self.assertEqual(sp.simplify(sum(p)), 1, betas)
            self.assertEqual(sp.simplify(sum(q ** 2 for q in p) - 1 + 2 * (1 - B) ** 2 / (3 - B) ** 2), 0, betas)
        isotropic = [sp.simplify((1 - B + 2 * b) / (3 - B)) for B, b in [(-ROOT3, -1 / ROOT3)]]
        self.assertEqual(isotropic, [sp.Rational(1, 3)])

    def test_conformal_chart_is_the_isotropic_members_in_conformal_time(self):
        chart = self.charts["conformal"]
        reader, _ = self.published("conformal", "0")
        eta, k = reader.symbol["\\eta"], reader.parameters["k"]
        _, scalar = self.published("conformal", chart["ricci_scalar"])
        for value, beta in ((-1 / (1 + ROOT3), -1 / ROOT3), (1 / (ROOT3 - 1), 1 / ROOT3)):
            self.assertEqual(sp.simplify(value / (1 + value) - beta), 0)
            self.assertEqual(sp.simplify(2 * value ** 2 - 2 * value - 1), 0)
            # At eta = -(1 + k) s, with -t = s^(1 + k), the scalar is (3 beta - 1)^2/t^2.
            for s in (sp.Rational(1, 3), 2):
                at = {eta: -(1 + value) * s, k: value}
                minus_t = s ** (1 + value)
                self.assertLess(abs(sp.N(scalar.subs(at) - (3 * beta - 1) ** 2 / minus_t ** 2, 30)), 1e-20)


if __name__ == "__main__":
    unittest.main()
