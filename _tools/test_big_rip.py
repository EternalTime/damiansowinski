"""The Big Rip, the flat Friedmann universe of phantom energy alone:
.venv.noindex/bin/python -m unittest discover -s _tools

What its texts and drawings state, held on the published files: the fluid's pressure is w times its
density and violates the null energy condition, the density grows with the scale factor, the time
left before the rip is Caldwell, Kamionkowski and Weinberg's, every bound system comes apart about
0.3 of its period before the end at w = -3/2, the event horizon's proper radius is Chiba, Takahashi
and Sugiyama's 3(1 + w)c(-t)/(1 + 3w), the limit w -> -1 is de Sitter's exponential, and the drawn
horizons sit where that radius puts them. The tests need sympy, and are skipped where it is absent.
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
sys.path.insert(0, str(Path(__file__).resolve().parent / "derivations"))


def chart(system):
    metric = json.loads((DATA / "metrics" / "big_rip.json").read_text(encoding="utf-8"))
    return next(c for c in metric["coordinates"] if c["id"] == system)


@unittest.skipUnless(HAS_SYMPY, "sympy is not installed")
class Physics(unittest.TestCase):
    def setUp(self):
        import sympy
        import verify_metrics
        self.sp, self.vm = sympy, verify_metrics
        entry = chart("comoving")
        self.reader = verify_metrics.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]], ())
        self.entry = entry
        R = self.reader
        self.t, self.w, self.t0, self.c = R.symbol["t"], R.parameters["w"], R.parameters["t_0"], R.c

    def mixed_einstein(self, i):
        ul = self.entry["einstein_tensor"]["variants"]["ul"]["nonzero"]
        return self.reader(next(e["value"] for e in ul if e["indices"] == [i, i]))

    def test_pressure_is_w_times_density_and_breaks_the_null_condition(self):
        sp = self.sp
        # G^t_t = -8 pi G rho/c^2 and G^r_r = 8 pi G p/c^4, in the chart x^0 = ct.
        rho, p = -self.mixed_einstein("t"), self.mixed_einstein("r")
        self.assertEqual(sp.simplify(p / rho - self.w), 0)
        at = {self.w: sp.Rational(-3, 2), self.t: -1, self.c: 1, self.t0: 1}
        self.assertGreater(float(rho.subs(at)), 0)
        self.assertLess(float((rho + p).subs(at)), 0)

    def test_density_grows_as_a_to_minus_three_one_plus_w(self):
        sp = self.sp
        rho = -self.mixed_einstein("t")
        a = (-self.t / self.t0) ** (2 / (3 * (1 + self.w)))
        slope = sp.simplify(sp.diff(sp.log(rho), self.t) / sp.diff(sp.log(a), self.t))
        self.assertEqual(sp.simplify(slope + 3 * (1 + self.w)), 0)

    def test_hubble_rate_and_time_left(self):
        sp = self.sp
        ull = self.entry["christoffel"]["variants"]["ull"]["nonzero"]
        # Gamma^r_tr = H/c in the chart x^0 = ct.
        H = self.c * self.reader(next(e["value"] for e in ull if e["indices"] == ["r", "t", "r"]))
        H0 = H.subs(self.t, -self.t0)
        self.assertEqual(sp.simplify(H0 * self.t0 + 2 / (3 * (1 + self.w))), 0)
        # Caldwell, Kamionkowski and Weinberg's 22 billion years for w = -3/2, H_0 = 70 km/s/Mpc and
        # Omega_m = 0.3, their (2/3)|1 + w|^-1 H_0^-1 (1 - Omega_m)^-1/2.
        hubble_time = 3.0857e19 / 70 / 3.156e16
        self.assertAlmostEqual((2 / 3) / 0.5 * hubble_time / math.sqrt(0.7), 22.2, delta=0.3)

    def test_bound_systems_come_apart_about_a_third_of_a_period_before_the_end(self):
        w = -1.5
        self.assertAlmostEqual(math.sqrt(2 * abs(1 + 3 * w)) / (6 * math.pi * abs(1 + w)), 0.3, delta=0.02)

    def test_event_horizon_radius(self):
        sp = self.sp
        t, w, t0, c = self.t, self.w, self.t0, self.c
        u = sp.Symbol("u", positive=True)
        g_rr = self.reader(next(e["value"] for e in self.entry["metric_components"] if e["indices"] == ["r", "r"]))
        a = sp.sqrt(g_rr)
        # Light from the centre covers the comoving distance int_t^0 c dt'/a before the rip.
        at = {w: sp.Rational(-3, 2)}
        reach = sp.integrate((c / a).subs(at).subs(t, -u), (u, 0, -t))
        proper = sp.simplify(sp.powsimp(sp.powdenest((a.subs(at) * reach), force=True), force=True))
        want = (3 * (1 + w) * c * (-t) / (1 + 3 * w)).subs(at)
        self.assertEqual(sp.simplify(proper.subs(t, -sp.Rational(5, 3)) - want.subs(t, -sp.Rational(5, 3))), 0)

    def test_de_sitter_is_the_limit_w_to_minus_one(self):
        sp = self.sp
        H, tau, eps = sp.symbols("H tau epsilon", positive=True)
        # w = -1 - eps, t_0 = 2/(3 eps H), and tau = t + t_0 the time since a = 1.
        t0 = 2 / (3 * eps * H)
        a = (-(tau - t0) / t0) ** (2 / (3 * (-eps)))
        self.assertEqual(sp.simplify(sp.limit(sp.log(a), eps, 0, "+") - H * tau), 0)

    def test_conformal_scale_factor_is_the_comoving_one(self):
        sp = self.sp
        entry = chart("conformal")
        R = self.vm.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]], ())
        eta, eta0, w = R.symbol["\\eta"], R.parameters["eta_0"], R.parameters["w"]
        g_rr = R(next(e["value"] for e in entry["metric_components"] if e["indices"] == ["r", "r"]))
        at = {w: sp.Rational(-3, 2), eta0: sp.Rational(3, 7)}
        # At w = -3/2 and t_0 = c = 1, eta = -(3/7)(-t)^(7/3), and a^2 = (-t)^(-8/3).
        t = sp.Rational(-1, 2)
        value = g_rr.subs(at).subs(eta, -sp.Rational(3, 7) * (-t) ** sp.Rational(7, 3))
        self.assertAlmostEqual(float(value), float((-t) ** sp.Rational(-8, 3)), places=12)


class Drawings(unittest.TestCase):
    def test_spacetime_diagrams_mark_the_event_horizon(self):
        diagrams = json.loads((DATA / "diagrams" / "big_rip.json").read_text(encoding="utf-8"))
        views = {(s, v["id"]): v for s, vs in diagrams["systems"].items() for v in vs}
        # Unit square of each view: the comoving box (0, 2, -2, 0) in r/ct_0 and t/t_0, the conformal
        # one (0, 3, -3, 0) in r/eta_0 and eta/eta_0; the horizon is r = (3/7)(-t)^(7/3) and r = -eta.
        for key, size, horizon in ((("comoving", "radial"), 2, lambda t: 3 / 7 * (-t) ** (7 / 3)),
                                   (("conformal", "radial"), 3, lambda eta: -eta)):
            self.assertEqual(views[key]["box"], [0, size, -size, 0])
            event = next(m for m in views[key]["markers"] if m["kind"] == "event")
            for x, y in event["lines"][0]:
                r, t = size * x, size * (y - 1)
                if 0.05 < r < size - 0.05 and t < -0.05:
                    self.assertAlmostEqual(r, horizon(t), delta=0.01)

    def test_embedding_horizon_shrinks_as_three_sevenths_of_c_times_the_time_left(self):
        embedding = json.loads((DATA / "embedding" / "big_rip.json").read_text(encoding="utf-8"))
        ring = next(v for v in embedding["views"] if v["id"] == "ring")
        for frame in ring["movie"]["frames"]:
            t = frame["value"]
            curves = {c["class"]: c for c in frame["curves"]}
            radius = [math.hypot(x, y) for x, y, _ in curves["horizon"]["points"]]
            self.assertAlmostEqual(max(radius), 3 * (-t) / 7, delta=1e-4)
            self.assertAlmostEqual(min(radius), 3 * (-t) / 7, delta=1e-4)
            galaxies = [math.hypot(x, y) for x, y, _ in curves["particles"]["points"]]
            self.assertAlmostEqual(max(galaxies), 0.5 * (-t) ** (-4 / 3), delta=1e-4)

    def test_conformal_diagram_is_the_lower_half_of_minkowski_space(self):
        conformal = json.loads((DATA / "conformal" / "big_rip.json").read_text(encoding="utf-8"))
        view = conformal["views"][0]
        region = next(layer for layer in view["layers"] if layer["class"] == "region")
        self.assertEqual(sorted(map(tuple, region["points"])), sorted([(0, 0), (round(math.pi, 4), 0), (0, -round(math.pi, 4))]))
        rip = next(layer for layer in view["layers"] if layer["class"] == "singular")
        self.assertTrue(all(abs(T) < 1e-9 for _, T in rip["points"]))
        horizon = next(layer for layer in view["layers"] if layer["class"] == "horizon")
        self.assertTrue(all(abs(X + T) < 1e-9 for X, T in horizon["points"]))


if __name__ == "__main__":
    unittest.main()
