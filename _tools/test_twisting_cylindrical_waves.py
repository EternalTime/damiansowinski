"""The cylindrical waves with two polarisations of Kompaneets and of Jordan, Ehlers and Kundt:
.venv.noindex/bin/python -m unittest discover -s _tools

What the page and its drawings state of the physics, held on the published files and on the pulse
of Piran, Safier, and Katz at a = 1 and alpha = 2 that every drawing draws. The first class needs
nothing but the files. The second needs numpy and sympy and is skipped where they are absent: it
holds the pulse to the field equations each chart's parameters state, to a regular axis, to the
symmetry in time, and to the limits far out that the captions state, and the null chart's
equations to the cylindrical chart's.
"""
import importlib.util
import json
import random
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "MFS" / "assets" / "data"
HAS_NUMBERS = all(importlib.util.find_spec(name) is not None for name in ("sympy", "numpy"))
sys.path.insert(0, str(Path(__file__).resolve().parent / "derivations"))


def load(folder):
    return json.loads((DATA / folder / "twisting_cylindrical_waves.json").read_text(encoding="utf-8"))


class Files(unittest.TestCase):
    def test_both_charts_carry_the_three_functions(self):
        charts = {c["id"]: c for c in load("metrics")["coordinates"]}
        self.assertEqual(set(charts), {"cylindrical", "null"})
        for chart in charts.values():
            names = [p["symbol"].split(" = ")[0] for p in chart["parameters"]]
            self.assertEqual(names, ["\\psi", "\\omega", "\\gamma"])

    def test_omega_is_the_cross_term(self):
        for chart in load("metrics")["coordinates"]:
            cross = {tuple(c["indices"]): c["value"] for c in chart["metric_components"]}
            self.assertIn(("\\phi", "z"), cross)
            self.assertIn("\\omega", cross[("\\phi", "z")])

    def test_every_drawing_takes_the_same_pulse(self):
        import null_rays as nr
        inputs = [v["input"] for vs in load("diagrams")["systems"].values() for v in vs]
        inputs += [v["input"] for v in load("conformal")["views"]]
        inputs += [v["input"] for v in load("embedding")["views"]]
        inputs += [v["input"] for v in load("curvature")["views"]]
        self.assertEqual(set(inputs), {nr.PSK_INPUT})
        self.assertEqual(load("curvature")["views"][0]["map"]["functions"], nr._psk_pulse("t", "rho"))


@unittest.skipUnless(HAS_NUMBERS, "needs numpy and sympy")
class Pulse(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import sympy as sp
        import null_rays as nr
        cls.sp = sp
        cls.t, cls.r = sp.symbols("t rho", real=True)
        cls.f = {k: sp.sympify(v, locals={"t": cls.t, "rho": cls.r}) for k, v in nr._psk_pulse("t", "rho").items()}

    def at(self, expr, t, r, digits=40):
        return self.sp.N(expr.subs({self.t: self.sp.Rational(t), self.r: self.sp.Rational(r)}), digits)

    def test_the_pulse_solves_the_field_equations_of_the_cylindrical_chart(self):
        """The two wave equations and the two quadratures the cylindrical chart's parameters state."""
        sp, t, r = self.sp, self.t, self.r
        P, W, G = self.f["psi"], self.f["omega"], self.f["gamma"]
        d = sp.diff
        e = sp.exp(4 * P)
        equations = [
            d(P, t, 2) - d(P, r, 2) - d(P, r) / r - e * (d(W, t) ** 2 - d(W, r) ** 2) / (2 * r ** 2),
            d(W, t, 2) - d(W, r, 2) + d(W, r) / r - 4 * (d(P, r) * d(W, r) - d(P, t) * d(W, t)),
            d(G, r) - r * (d(P, r) ** 2 + d(P, t) ** 2) - e * (d(W, r) ** 2 + d(W, t) ** 2) / (4 * r),
            d(G, t) - 2 * r * d(P, r) * d(P, t) - e * d(W, r) * d(W, t) / (2 * r),
        ]
        rng = random.Random(1986)
        for _ in range(4):
            t0, r0 = f"{rng.randint(-60, 60)}/10", f"{rng.randint(1, 80)}/10"
            for eq in equations:
                self.assertLess(abs(self.at(eq, t0, r0)), 1e-30)

    def test_the_axis_is_regular(self):
        """gamma = omega = 0 on the axis, where lambda_u lambda_v = 1, so there is no conical singularity."""
        for t0 in ("-3", "0", "1/2", "5"):
            self.assertLess(abs(self.at(self.f["gamma"], t0, "1e-30")), 1e-25)
            self.assertLess(abs(self.at(self.f["omega"], t0, "1e-30")), 1e-25)

    def test_the_pulse_is_even_in_time(self):
        for name in ("psi", "omega", "gamma"):
            for t0, r0 in (("3/2", "1"), ("4", "7/2"), ("1/3", "6")):
                self.assertLess(abs(self.at(self.f[name], t0, r0) - self.at(self.f[name], "-" + t0, r0)), 1e-30)

    def test_far_out_space_is_a_cone_of_deficit_pi(self):
        """psi -> 0 and gamma -> ln alpha = ln 2, so the plane z = 0 tends to 4 drho^2 + rho^2 dphi^2,
        and omega tends to its constant sqrt(3)(2 + sqrt(3))."""
        sp = self.sp
        far = "1000000000000"
        self.assertLess(abs(self.at(self.f["psi"], "1", far)), 1e-10)
        self.assertLess(abs(self.at(self.f["gamma"], "1", far) - sp.log(2)), 1e-10)
        self.assertLess(abs(self.at(self.f["omega"], "1", far) - sp.sqrt(3) * (2 + sp.sqrt(3))), 1e-10)

    def test_the_pulse_is_deepest_on_the_axis_at_t_zero(self):
        """The diagrams' captions: psi = -1.32 on the axis at t = 0."""
        self.assertAlmostEqual(float(self.at(self.f["psi"], "0", "1e-30")), -1.317, places=3)

    def test_alpha_one_is_flat(self):
        """With alpha = 1 the pulse of d'Inverno, Dubal and Sarkies's (17) to (22) is Minkowski space."""
        sp = self.sp
        a = 1
        u, v = self.t - self.r, self.t + self.r
        lu, lv = sp.sqrt(1 + u ** 2) - u, sp.sqrt(1 + v ** 2) + v
        xi = 1 + lu * lv
        top = a ** 2 * (1 - lu * lv) ** 2 + (lu + lv) ** 2
        for expr in (top / (a ** 2 * xi ** 2 + (lu - lv) ** 2), top / ((1 + lu ** 2) * (1 + lv ** 2))):
            for t0, r0 in (("1", "2"), ("-3", "1/2"), ("7", "5")):
                self.assertLess(abs(self.at(expr, t0, r0) - 1), 1e-30)


if __name__ == "__main__":
    unittest.main()
