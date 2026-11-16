"""The joining of the fireball to the free expansion, and what the drawings rest on.

The chart printer holds each chart to its field equations and to the flat background.
The collection checker holds every stored tensor and drawing, so those are not repeated here.
"""
import sys
import unittest
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "_tools" / "derivations"))
import brout_englert_gunzig as beg


def vanishes(expression, **values):
    """An expression that holds an inverse hyperbolic tangent inside a hyperbolic function, which
    sympy leaves standing: it is zero at each of three points, to thirty digits."""
    points = ((sp.Rational(5, 2), sp.Rational(4, 5), sp.Rational(1, 3)),
              (sp.Rational(7, 4), sp.Rational(3, 7), sp.Rational(9, 10)),
              (sp.Integer(1), sp.Rational(19, 20), sp.Rational(1, 2)))
    for end, x, part in points:
        at = {s: v for s, v in ((values.get("end"), end), (values.get("x"), x),
                                (values.get("tau"), part*end)) if s is not None}
        if abs(sp.N(expression.subs(at), 40)) > 1e-30:
            return False
    return True


class Joining(unittest.TestCase):
    """The free expansion that meets the fireball at tau_1 with no jump in the conformal factor
    or in its first derivative."""

    def setUp(self):
        self.tau, self.end, self.x = sp.symbols("tau tau_infty x", positive=True)
        self.tau_1 = self.x*self.end
        self.k = beg.joined(self.end, self.tau_1)
        self.fireball = 1/(1 - self.tau**2/self.end**2)
        self.expansion = self.k["Omega_f"]*(1 - self.k["tau_b"]**2/self.tau**2)

    def gaps(self, expansion):
        return [sp.simplify(sp.diff(self.fireball - expansion, self.tau, order).subs(self.tau, self.tau_1))
                for order in (0, 1)]

    def test_the_factor_and_its_slope_are_continuous(self):
        self.assertEqual(self.gaps(self.expansion), [0, 0])

    def test_another_constant_leaves_a_jump(self):
        wrong = self.k["Omega_1"]*(1 - self.k["tau_b"]**2/self.tau**2)
        self.assertNotEqual(self.gaps(wrong)[0], 0)
        shifted = self.k["Omega_f"]*(1 - (self.k["tau_b"]/2)**2/self.tau**2)
        self.assertNotEqual(self.gaps(shifted)[1], 0)

    def test_the_final_factor_is_the_square_of_the_factor_at_the_joining(self):
        # The authors' figure: e^lambda = 10^16 where creation ends and 10^32 at late times.
        k = beg.joined(1, sp.sqrt(1 - sp.Rational(1, 10**8)))
        self.assertEqual(sp.simplify(k["Omega_1"]), 10**8)
        self.assertEqual(sp.simplify(k["Omega_f"]), 10**16)

    def test_the_constants_in_proper_and_in_conformal_time(self):
        H, T_1 = self.k["H"], self.k["T_1"]
        names = {"end": self.end, "x": self.x}
        self.assertTrue(vanishes(sp.tanh(H*T_1/2) - self.x, **names))
        self.assertTrue(vanishes(self.k["b"] - 2*self.k["tau_b"]*self.k["Omega_f"], **names))
        self.assertTrue(vanishes(self.k["b"] - 1/(H*sp.sinh(self.k["eta_1"])**2), **names))
        self.assertTrue(vanishes(self.k["eta_b"] - sp.log(self.k["tau_b"]/self.end), **names))
        self.assertFalse(vanishes(self.k["b"] - self.k["tau_b"]*self.k["Omega_f"], **names))


class OneScaleFactor(unittest.TestCase):
    """The three coordinate systems state one scale factor, a = c tau Omega, at each stage."""

    def test_the_fireball(self):
        tau, end = sp.symbols("tau tau_infty", positive=True)
        H, a = 2/end, tau/(1 - tau**2/end**2)
        proper = end*sp.atanh(tau/end)
        names = {"end": end, "tau": tau}
        self.assertTrue(vanishes(sp.sinh(H*proper)/H - a, **names))
        self.assertTrue(vanishes(sp.diff(proper, tau) - a/tau, **names))
        eta = sp.log(tau/end)
        self.assertTrue(vanishes(-1/(H*sp.sinh(eta)) - a, **names))
        self.assertFalse(vanishes(sp.sinh(H*proper/2)*2/H - a, **names))

    def test_the_free_expansion(self):
        for end, tau_1 in ((sp.Rational(5, 2), sp.Integer(2)), (sp.Rational(7, 4), sp.Rational(3, 4))):
            k = beg.joined(end, tau_1)
            tau = sp.Symbol("tau", positive=True)
            a = k["Omega_f"]*(tau - k["tau_b"]**2/tau)
            # dT = Omega dtau from T_1 at tau_1
            proper = k["T_1"] + k["Omega_f"]*((tau - tau_1) + k["tau_b"]**2*(1/tau - 1/tau_1))
            eta = sp.log(tau/end)
            for value in (tau_1, 2*tau_1, 7*tau_1):
                at = {tau: value}
                self.assertLess(abs(sp.N((sp.sqrt((proper - k["T_b"])**2 - k["b"]**2) - a).subs(at), 40)), 1e-30)
                self.assertLess(abs(sp.N((k["b"]*sp.sinh(eta - k["eta_b"]) - a).subs(at), 40)), 1e-30)


class Causality(unittest.TestCase):
    def test_the_conformal_time_diverges_at_the_origin(self):
        # eta = ln tanh(HT/2) has no lower bound, so no comoving world line has a particle horizon.
        T, H = sp.symbols("T H", positive=True)
        self.assertEqual(sp.limit(sp.log(sp.tanh(H*T/2)), T, 0, "+"), -sp.oo)

    def test_the_fireball_alone_would_take_an_unbounded_proper_time_to_reach_its_pole(self):
        tau, end = sp.symbols("tau tau_infty", positive=True)
        self.assertEqual(sp.limit(end*sp.atanh(tau/end), tau, end, "-"), sp.oo)
        self.assertEqual(sp.limit(1/(1 - tau**2/end**2), tau, end, "-"), sp.oo)


class EmbeddedMoment(unittest.TestCase):
    """The plane of one background time t: rho = Omega r, and a surface of revolution in flat
    space carries it while g_rr >= (d rho/dr)^2."""

    def setUp(self):
        self.t, self.r, self.end = sp.symbols("t r tau_infty", positive=True)

    def room(self, Omega):
        return sp.simplify(Omega**2 - sp.diff(Omega*self.r, self.r)**2)

    def test_the_fireball_embeds_until_the_background_time_tau_infinity(self):
        t, r, end = self.t, self.r, self.end
        Omega = 1/(1 - (t*t - r*r)/end**2)
        expected = 4*r*r*Omega**3*(end*end - t*t)/end**4*Omega
        self.assertEqual(sp.simplify(self.room(Omega) - expected), 0)

    def test_the_free_expansion_embeds_until_the_same_time(self):
        # Its room vanishes where tau^4 = t^2 tau_b^2, and tau >= tau_1 with tau_b = tau_1^2/tau_infinity.
        t, r, b = self.t, self.r, sp.Symbol("tau_b", positive=True)
        tau2 = t*t - r*r
        Omega = 1 - b*b/tau2
        expected = 4*r*r*b*b*(tau2**2 - t*t*b*b)/tau2**4
        self.assertEqual(sp.simplify(self.room(Omega) - expected), 0)

    def test_the_crease_on_the_light_cone(self):
        # d(rho)/dr over sqrt(g_rr) at r = t is 1 - 2 t^2/tau_infinity^2 inside the cone and 1 outside.
        t, r, end = self.t, self.r, self.end
        Omega = 1/(1 - (t*t - r*r)/end**2)
        cosine = (sp.diff(Omega*r, r)/Omega).subs(r, t)
        self.assertEqual(sp.simplify(cosine - (1 - 2*t*t/end**2)), 0)
        self.assertLessEqual(beg.MOMENT, beg.END)
        self.assertGreater(beg.MOMENT, beg.RATIO*beg.END)


class Curvature(unittest.TestCase):
    def test_the_kretschmann_scalar_is_continuous_where_the_ricci_scalar_drops(self):
        end, x = sp.symbols("tau_infty x", positive=True)
        k = beg.joined(end, x*end)
        H = k["H"]
        after = 24*k["b"]**4/((k["T_1"] - k["T_b"])**2 - k["b"]**2)**4
        self.assertTrue(vanishes(after/H**4 - 24, end=end, x=x))


if __name__ == "__main__":
    unittest.main()
