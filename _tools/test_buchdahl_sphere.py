"""Buchdahl's fluid sphere, held on the published file with the checker's own reader:
.venv.noindex/bin/python -m unittest _tools/test_buchdahl_sphere.py

The published Einstein tensor is a perfect fluid at rest with 3p/(rho c^2) = f/(1 - f),
Buchdahl's (2.12); its pressure and density obey his equation of state (3.3),
p/p_c = (rho/rho_c)^(6/5)/(1 + 2q_c(1 - (rho/rho_c)^(1/5))), which goes over into the law of the
polytrope of index 5 as q_c goes to zero; his f is Schwarzschild's r_s/(4r) far out, so the mass is
the r_s of the chart; and light has circular orbits once GM/(c^2 R) reaches Igata and Takamori's
1.3629, the root of 3a^6 - 1512a^4 + 10240a^2 - 13824 = 0. Needs sympy, and is skipped where it is
absent.
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


@unittest.skipUnless(HAS_SYMPY, "sympy is not installed")
class TheSphere(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import sympy
        import verify_metrics
        cls.sp = sympy
        chart = json.loads((DATA / "metrics" / "buchdahl_sphere.json").read_text())["coordinates"][0]
        reader = verify_metrics.Reader(chart["coords"], [p["symbol"] for p in chart["parameters"]], ())
        cls.r, cls.rs, cls.R = reader.symbol["r"], reader.parameters["r_s"], reader.parameters["R"]
        mixed = {tuple(e["indices"]): reader(e["value"]) for e in chart["einstein_tensor"]["variants"]["ul"]["nonzero"]}
        # 8 pi G rho/c^2 and 8 pi G p/c^4, as functions of r at r_s = 1.
        cls.density = staticmethod(sympy.lambdify((cls.r, cls.R), -mixed[("t", "t")].subs(cls.rs, 1), "mpmath"))
        cls.pressure = staticmethod(sympy.lambdify((cls.r, cls.R), mixed[("r", "r")].subs(cls.rs, 1), "mpmath"))
        cls.angular = staticmethod(sympy.lambdify((cls.r, cls.R), mixed[("\\theta", "\\theta")].subs(cls.rs, 1), "mpmath"))
        g = {tuple(e["indices"]): reader(e["value"]) for e in chart["metric_components"]}
        cls.gtt, cls.gpp = g[("t", "t")], g[("\\theta", "\\theta")]

    def f(self, r, R):
        import mpmath
        return 1 / (4 * mpmath.sqrt(r * r + R * R))

    def test_the_fluid_is_at_rest_with_buchdahls_ratio_of_pressure_to_density(self):
        import mpmath
        mpmath.mp.dps = 30
        for r, R in (("0", "1"), ("0.7", "1"), ("3", "2"), ("11", "0.6")):
            r, R = mpmath.mpf(r), mpmath.mpf(R)
            f = self.f(r, R)
            p, rho = self.pressure(r, R), self.density(r, R)
            self.assertLess(abs(self.angular(r, R) / p - 1), 1e-20)
            self.assertLess(abs(3 * p / rho - f / (1 - f)), 1e-20)

    def test_pressure_and_density_obey_buchdahls_equation_of_state(self):
        import mpmath
        mpmath.mp.dps = 30
        R = mpmath.mpf(1)
        fc = self.f(0, R)
        qc = fc / (1 - fc)
        for r in (0.3, 1, 2.5, 7):
            x = self.density(r, R) / self.density(0, R)
            law = x ** (mpmath.mpf(6) / 5) / (1 + 2 * qc * (1 - x ** (mpmath.mpf(1) / 5)))
            self.assertLess(abs(self.pressure(r, R) / self.pressure(0, R) - law), 1e-20)

    def test_a_weak_centre_follows_the_law_of_index_five(self):
        import mpmath
        mpmath.mp.dps = 30
        R = mpmath.mpf(400)
        x = self.density(800, R) / self.density(0, R)
        ratio = self.pressure(800, R) / self.pressure(0, R)
        self.assertLess(abs(ratio / x ** (mpmath.mpf(6) / 5) - 1), 2e-3)

    def test_the_mass_is_the_schwarzschild_radius_of_the_chart(self):
        # sqrt(-g_tt) = 1 - r_s/(2r) + ... far out, Schwarzschild's lapse of the mass r_s c^2/(2G).
        sp = self.sp
        lapse = sp.sqrt(-self.gtt).subs({self.rs: 1, self.R: 1})
        self.assertEqual(sp.limit(self.r * (lapse - 1), self.r, sp.oo), -sp.Rational(1, 2))

    def test_light_circles_once_the_core_is_compact_enough(self):
        import mpmath
        mpmath.mp.dps = 30
        a0 = mpmath.findroot(lambda a: 3 * a ** 6 - 1512 * a ** 4 + 10240 * a ** 2 - 13824, 1.36)
        self.assertLess(abs(a0 - mpmath.mpf("1.3629")), 1e-4)

        def turns(a):
            # V = -g_tt/g_thth along the radius, R = 1 and r_s = 2a; light circles where V' = 0.
            V = self.sp.lambdify(self.r, (-self.gtt / self.gpp).subs({self.R: 1, self.rs: 2 * a}), "mpmath")
            slopes = [mpmath.diff(V, mpmath.mpf(k) / 100) for k in range(10, 400)]
            return sum(1 for u, v in zip(slopes, slopes[1:]) if u * v < 0)
        self.assertEqual(turns(a0 - mpmath.mpf("0.01")), 0)
        self.assertEqual(turns(a0 + mpmath.mpf("0.01")), 2)


if __name__ == "__main__":
    unittest.main()
