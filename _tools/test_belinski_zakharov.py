"""Belinski and Zakharov's gravitational solitons, the wave of two solitons on the Kasner universe
that expands alike in two directions:
.venv.noindex/bin/python -m unittest discover -s _tools

What its texts and drawings state, held on the published files at the values every drawing takes,
w = 1 and cosh(beta) = 5/4. The first class needs nothing but the files. The second needs sympy
and is skipped where it is absent: it reads the published charts through the checker's Reader and
holds the pole chart to Belinski and Zakharov's own formula, to the Kasner universe on both sides
of the pulses, and the canonical chart's published Einstein tensor to vanishing on the wave.
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

P1, P2 = 5 / 4, 3 / 4
MOMENTS = (0.25, 0.5, 1.0, 2.0)


def load(folder):
    return json.loads((DATA / folder / "belinski_zakharov.json").read_text(encoding="utf-8"))


def semi_axes(tau):
    """sqrt(g_xx) and sqrt(g_yy) on xi = 0, in units of w."""
    up, down = P1 * math.cosh(tau) + 1, P1 * math.cosh(tau) - 1
    return math.sqrt(math.sinh(tau) * up / down), math.sqrt(math.sinh(tau) * down / up)


class Drawings(unittest.TestCase):
    def test_the_rays_of_both_planes_are_at_45_degrees(self):
        """Both planes are conformally flat, so every ray keeps tau -+ xi, or ct -+ z."""
        systems = load("diagrams")["systems"]
        for system in ("pole", "canonical"):
            view = systems[system][0]
            X0, X1, Y0, Y1 = view["box"]
            for family in view["rays"].values():
                for ray in family:
                    xs = [X0 + u * (X1 - X0) for u, _ in ray]
                    ys = [Y0 + v * (Y1 - Y0) for _, v in ray]
                    plus = [y + x for x, y in zip(xs, ys)]
                    minus = [y - x for x, y in zip(xs, ys)]
                    self.assertLess(min(max(plus) - min(plus), max(minus) - min(minus)), 2e-3 * (X1 - X0))

    def test_the_marked_rays_are_the_light_cone_of_the_origin(self):
        """xi = +-tau in the pole chart and z = +-ct in the canonical one, the rays the pulses run along."""
        systems = load("diagrams")["systems"]
        for system in ("pole", "canonical"):
            view = systems[system][0]
            X0, X1, Y0, Y1 = view["box"]
            lines = [line for m in view["markers"] if m["kind"] == "shell" for line in m["lines"]]
            self.assertEqual(len(lines), 2)
            sides = set()
            for line in lines:
                self.assertGreater(len(line), 1)
                for u, v in line:
                    x, y = X0 + u * (X1 - X0), Y0 + v * (Y1 - Y0)
                    self.assertAlmostEqual(abs(x), y, delta=2e-3 * (X1 - X0))
                sides.add(line[-1][0] > 0.5 or line[0][0] > 0.5)
            self.assertEqual(sides, {True, False})

    def test_the_ring_is_the_ellipse_of_the_published_block(self):
        """On xi = 0 the block is diagonal: the ring of the tube at each named moment has the
        semi-axes sqrt(g_xx) and sqrt(g_yy), whose product is sinh(tau), the area Kasner's universe
        gives the ring, and whose ratio is (p cosh(tau) + 1)/(p cosh(tau) - 1), p = cosh(beta)."""
        tube = next(v for v in load("embedding")["views"] if v["id"] == "tube")
        curves = tube["surfaces"][0]["curves"]
        self.assertEqual([c["time"] for c in curves], list(MOMENTS))
        for curve in curves:
            tau = curve["time"]
            a = max(abs(p[0]) for p in curve["points"])
            b = max(abs(p[1]) for p in curve["points"])
            want_a, want_b = semi_axes(tau)
            self.assertAlmostEqual(a, want_a, delta=1e-5)
            self.assertAlmostEqual(b, want_b, delta=1e-5)
            self.assertAlmostEqual(a * b, math.sinh(tau), delta=1e-5)
            self.assertAlmostEqual(a / b, (P1 * math.cosh(tau) + 1) / (P1 * math.cosh(tau) - 1), delta=1e-4)

    def test_the_ring_rounds_off_as_the_pulses_leave(self):
        """The ratio of the ring's axes falls at every frame of the movie and its area grows."""
        ring = next(v for v in load("embedding")["views"] if v["id"] == "ring")
        times = [frame["value"] for frame in ring["movie"]["frames"]]
        self.assertEqual(times, sorted(times))
        ratios = [semi_axes(t)[0] / semi_axes(t)[1] for t in times]
        self.assertTrue(all(later < earlier for earlier, later in zip(ratios, ratios[1:])))
        self.assertGreater(ratios[0], 7)
        self.assertLess(ratios[-1], 1.6)

    def test_each_moment_stands_on_the_line_midway_between_the_pulses(self):
        """xi = 0 at tau is z = 0 at ct = w sinh(tau): X = 0 and T = 2 arctan(sinh(tau)) in the conformal
        diagram, and the same event in both of its views."""
        views = load("conformal")["views"]
        self.assertEqual([v["id"] for v in views], ["canonical", "pole"])
        for view in views:
            marked = sorted(tuple(s["points"][0]) for s in view["slices"] if s["view"] == "tube")
            self.assertEqual(len(marked), len(MOMENTS))
            for (X, T), tau in zip(marked, MOMENTS):
                self.assertAlmostEqual(X, 0.0, delta=1e-4)
                self.assertAlmostEqual(T, 2 * math.atan(math.sinh(tau)), delta=1e-4)


@unittest.skipUnless(HAS_SYMPY, "needs sympy")
class Physics(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import sympy as sp
        import verify_metrics as vm
        cls.sp, cls.vm = sp, vm
        cls.charts = {c["id"]: c for c in load("metrics")["coordinates"]}

    def reader(self, system):
        chart = self.charts[system]
        return self.vm.Reader(chart["coords"], [p["symbol"] for p in chart["parameters"]], ())

    def components(self, system, field="metric_components"):
        R = self.reader(system)
        return R, {tuple(e["indices"]): R(e["value"]) for e in self.charts[system][field]}

    def test_the_pole_chart_is_belinski_and_zakharovs_formula(self):
        """Their (5.10) with sigma = tanh^2(tau/2) and sin^2(phi) = 1/cosh^2(xi), (5.16), written
        out here apart from the generator: the block, and the conformal factor
        alpha^(3/2) Q/(4 sigma)."""
        sp = self.sp
        R, g = self.components("pole")
        tau, xi = R.symbol["\\tau"], R.symbol["\\xi"]
        beta, w = R.parameters["beta"], R.parameters["w"]
        for at in ({tau: 0.7, xi: 0.4, beta: 0.9, w: 1.3}, {tau: 2.1, xi: -1.6, beta: -0.5, w: 0.8}):
            t, x, b, scale = (at[s] for s in (tau, xi, beta, w))
            p1, p2 = math.cosh(b), math.sinh(b)
            alpha, sigma = math.sinh(t) * math.cosh(x), math.tanh(t / 2) ** 2
            sin2 = 1 / math.cosh(x) ** 2
            cos_2phi, sin_2phi = 1 - 2 * sin2, 2 * math.tanh(x) / math.cosh(x)
            H = 1 + sigma ** 2 - 2 * sigma * cos_2phi
            Q = p1 ** 2 * H - (1 - sigma) ** 2
            want = {("\\xi", "\\xi"): scale ** 2 * alpha ** 1.5 * Q / (4 * sigma),
                    ("x", "x"): alpha / Q * (p1 ** 2 * H - (1 - sigma) ** 2 * cos_2phi + 2 * p1 * (1 - sigma ** 2) * sin2),
                    ("y", "y"): alpha / Q * (p1 ** 2 * H - (1 - sigma) ** 2 * cos_2phi - 2 * p1 * (1 - sigma ** 2) * sin2),
                    ("x", "y"): -alpha / Q * p2 * (1 - sigma) ** 2 * sin_2phi}
            want[("\\tau", "\\tau")] = -want[("\\xi", "\\xi")]
            for index, value in want.items():
                self.assertAlmostEqual(float(sp.sympify(g[index]).xreplace(at)) / value, 1.0, places=10)

    def test_the_wave_is_a_vacuum_with_a_singularity_at_the_start(self):
        """No Ricci component is published, and the Kretschmann scalar grows without bound toward tau = 0."""
        chart = self.charts["pole"]
        for variant in chart["ricci_tensor"]["variants"].values():
            self.assertEqual(variant["nonzero"], [])
        R = self.reader("pole")
        K = R(chart["kretschmann"].partition("=")[2])
        tau, xi = R.symbol["\\tau"], R.symbol["\\xi"]
        fixed = {R.parameters["beta"]: math.log(2), R.parameters["w"]: 1}
        for where in (-1.5, 0.0, 0.8):
            near = abs(float(K.xreplace({**fixed, tau: 1e-2, xi: where})))
            nearer = abs(float(K.xreplace({**fixed, tau: 1e-3, xi: where})))
            self.assertGreater(nearer, 1e8)
            self.assertGreater(nearer, 50 * near)

    def test_the_universe_is_kasners_on_both_sides_of_the_pulses(self):
        """Far inside the light cone and far outside it the block tends to (ct/w) diag(1, 1) and
        the Kretschmann scalar to that of the Kasner universe with exponents (-1/3, 2/3, 2/3),
        K = (64/27)/T^4 in the proper time T = (4/3) k w (ct/w)^(3/4), where f tends to
        k^2 (ct/w)^(-1/2): k = cosh(beta) inside and sinh(beta) outside."""
        R, g = self.components("pole")
        K = R(self.charts["pole"]["kretschmann"].partition("=")[2])
        tau, xi = R.symbol["\\tau"], R.symbol["\\xi"]
        fixed = {R.parameters["beta"]: math.log(2), R.parameters["w"]: 1}
        for at, k in (({tau: 12.0, xi: 0.0}, P1), ({tau: 1e-4, xi: 10.0}, P2)):
            at = {**fixed, **at}
            alpha = math.sinh(at[tau]) * math.cosh(at[xi])
            for index, want in ((("x", "x"), alpha), (("y", "y"), alpha)):
                self.assertAlmostEqual(float(g[index].xreplace(at)) / want, 1.0, delta=2e-3)
            self.assertLess(abs(float(g[("x", "y")].xreplace(at))) / alpha, 2e-3)
            T = 4 / 3 * k * alpha ** 0.75
            self.assertAlmostEqual(float(K.xreplace(at)) * T ** 4, 64 / 27, delta=2e-2)

    def test_the_canonical_charts_einstein_tensor_vanishes_on_the_wave(self):
        """The published G_mu_nu of the free chart, with the wave's f, P and Q as the spacetime diagram
        declares them, every derivative taken in sympy: each component is zero to rounding."""
        sp = self.sp
        import null_rays as nr
        chart = self.charts["canonical"]
        R = self.reader("canonical")
        t, z = R.symbol["t"], R.symbol["z"]
        wave = {name: sp.sympify(text, locals={"t": t, "z": z}) for name, text in nr._bz_wave().items()}
        at = {t: sp.Rational(17, 10), z: sp.Rational(-9, 10), R.parameters["w"]: 1, R.c: 1}
        for entry in chart["einstein_tensor"]["variants"]["ll"]["nonzero"]:
            value = R(entry["value"])
            for name, expression in wave.items():
                value = value.replace(R.parameters[name].func, sp.Lambda((t, z), expression))
            self.assertLess(abs(value.doit().xreplace(at).evalf(30)), 1e-18, entry["indices"])


if __name__ == "__main__":
    unittest.main()
