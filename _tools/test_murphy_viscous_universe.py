"""Murphy's universe with viscosity, the flat Friedmann universe of a fluid whose bulk viscosity is
proportional to its density:
.venv.noindex/bin/python -m unittest discover -s _tools

What its texts and drawings state, held on the published files at points: each chart is the viscous
fluid at rest, with Friedmann's equation for H = H_s/(1 + a^q) and the pressure w rho c^2 lowered by
the viscous pressure -3 zeta H; the cosmic chart's scale factor is a^q = W(e^(q H_s t)) and solves the
same rate; the strong energy condition breaks while a^q < 1/(q - 1) and the weak one holds; both
scalars are finite and tend to de Sitter's as a -> 0; light reaches a -> 0 at a finite affine
distance, and the conformal time a - 1/a of radiation runs over the whole line. The tests need sympy,
and are skipped where it is absent.
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
HAS_SYMPY = importlib.util.find_spec("sympy") is not None
sys.path.insert(0, str(Path(__file__).resolve().parent / "derivations"))


def chart(system):
    metric = json.loads((DATA / "metrics" / "murphy_viscous_universe.json").read_text(encoding="utf-8"))
    return next(c for c in metric["coordinates"] if c["id"] == system)


def points(n=8, seed=1973):
    rng = random.Random(seed)
    return [(rng.uniform(0.05, 6), rng.uniform(0.2, 4), rng.uniform(0.3, 3), rng.uniform(0.5, 3)) for _ in range(n)]


@unittest.skipUnless(HAS_SYMPY, "sympy is not installed")
class Physics(unittest.TestCase):
    def setUp(self):
        import sympy
        import verify_metrics
        self.sp, self.vm = sympy, verify_metrics

    def reader(self, entry):
        return self.vm.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]], ())

    def mixed(self, entry, reader, i):
        ul = entry["einstein_tensor"]["variants"]["ul"]["nonzero"]
        return reader(next(e["value"] for e in ul if e["indices"] == [i, i]))

    def test_each_scale_factor_chart_is_the_viscous_fluid_at_rest(self):
        for system, space in (("scale_factor", "x"), ("scale_factor_spherical", "r")):
            entry = chart(system)
            R = self.reader(entry)
            a, q, Hs, c = R.symbol["a"], R.parameters["q"], R.parameters["H_s"], R.c
            G00, Gii = self.mixed(entry, R, "a"), self.mixed(entry, R, space)
            for av, qv, hv, cv in points():
                at = {a: av, q: qv, Hs: hv, c: cv}
                H = hv / (1 + av ** qv)
                w = 2 * qv / 3 - 1
                self.assertAlmostEqual(float(G00.subs(at)), -3 * H ** 2 / cv ** 2, delta=1e-12 * H ** 2)
                self.assertAlmostEqual(float(Gii.subs(at)), (w - (1 + w) * H / hv) * 3 * H ** 2 / cv ** 2,
                                       delta=1e-12 * H ** 2, msg=system)

    def test_cosmic_scale_factor_is_lambert_w_and_solves_the_rate(self):
        sp = self.sp
        t = sp.Symbol("t", real=True)
        for qv in (sp.Rational(3, 2), 2, 3):
            a = sp.LambertW(sp.exp(qv * t)) ** (1 / sp.S(qv))
            for tv in (-3, -0.5, 0, 1.2, 4):
                av = float(a.subs(t, tv))
                self.assertAlmostEqual(math.log(av) + av ** float(qv) / float(qv), tv, places=10)
                rate = float(sp.diff(a, t).subs(t, tv)) / av
                self.assertAlmostEqual(rate, 1 / (1 + av ** float(qv)), places=10)

    def test_the_strong_energy_condition_breaks_early_and_the_weak_one_holds(self):
        for qv in (1.5, 2.0, 3.0):
            w = 2 * qv / 3 - 1
            for av in (0.05, 0.3, 0.9, 1.1, 2.0, 5.0):
                u = 1 / (1 + av ** qv)
                pressure = w - (1 + w) * u          # in units of rho c^2
                self.assertGreaterEqual(1 + pressure, 0)
                accelerates = av ** qv < 1 / (qv - 1)
                self.assertEqual(1 + 3 * pressure < 0, accelerates, (qv, av))

    def test_both_scalars_are_finite_and_tend_to_de_sitter(self):
        entry = chart("scale_factor")
        R = self.reader(entry)
        a, q, Hs, c = R.symbol["a"], R.parameters["q"], R.parameters["H_s"], R.c
        ricci = R(entry["ricci_scalar"].split("=", 1)[1])
        kretschmann = R(entry["kretschmann"].split("=", 1)[1])
        for qv in (1.5, 2.0, 3.0):
            at = {q: qv, Hs: 1, c: 1}
            self.assertAlmostEqual(float(ricci.subs({**at, a: 1e-9})), 12, places=6)
            self.assertAlmostEqual(float(kretschmann.subs({**at, a: 1e-9})), 24, places=6)
            for av in (0.01, 0.5, 1, 3, 30):
                self.assertTrue(math.isfinite(float(kretschmann.subs({**at, a: av}))))
        # Radiation: R = 12/(1 + a^2)^3.
        self.assertAlmostEqual(float(ricci.subs({q: 2, Hs: 1, c: 1, a: 1})), 1.5, places=12)

    def test_light_reaches_the_steady_state_at_a_finite_affine_distance(self):
        # d lambda is proportional to a dt = (1 + a^q) da/H_s, whose integral from 0 to 1 is 1 + 1/(q + 1).
        for qv in (1.5, 2.0, 3.0):
            n = 20000
            h = 1 / n
            total = sum((1 + ((k + 0.5) * h) ** qv) * h for k in range(n))
            self.assertAlmostEqual(total, 1 + 1 / (qv + 1), places=6)

    def test_conformal_time_of_radiation_runs_over_the_whole_line(self):
        # c d eta = (1 + a^2) da/a^2 for q = 2, so eta = a - 1/a, unbounded at both ends.
        def eta(a):
            return a - 1 / a
        for lo, hi in ((0.2, 0.7), (1, 3)):
            n = 20000
            h = (hi - lo) / n
            integral = sum((1 + (lo + (k + 0.5) * h) ** 2) / (lo + (k + 0.5) * h) ** 2 * h for k in range(n))
            self.assertAlmostEqual(integral, eta(hi) - eta(lo), places=6)
        self.assertLess(eta(1e-6), -1e5)
        self.assertGreater(eta(1e6), 1e5)


class Drawings(unittest.TestCase):
    def test_the_hubble_sphere_comes_nearest_the_centre_where_the_acceleration_ends(self):
        diagram = json.loads((DATA / "diagrams" / "murphy_viscous_universe.json").read_text(encoding="utf-8"))
        view = diagram["systems"]["scale_factor_spherical"]
        view = view["views"][0] if isinstance(view, dict) and "views" in view else view[0]
        marker = next(m for m in view["markers"] if m["kind"] == "apparent")
        X0, X1, Y0, Y1 = view["box"]
        X, Y = min((p for line in marker["lines"] for p in line), key=lambda p: p[0])
        # r = c(1 + a^2)/H_s a is least, 2c/H_s, at a = 1.
        self.assertAlmostEqual(X0 + X * (X1 - X0), 2.0, delta=0.02)
        self.assertAlmostEqual(Y0 + Y * (Y1 - Y0), 1.0, delta=0.05)


if __name__ == "__main__":
    unittest.main()
