"""The Giddings-Strominger wormhole, a Riemannian space of four dimensions with no time in it:
.venv.noindex/bin/python -m unittest discover -s _tools

What its charts, texts and drawings state, held on the published files. Each chart's published
curvature is read through the checker's own Reader at points and held to the solution written out
here: the curvature scalar -6a^4/r^6, the Einstein tensor -3a^4/r^6 along the radius and 3a^4/r^6
along each direction on the sphere, and the Ricci tensor (-6, 0, 0, 0)a^4/r^6, whose eigenvalues
are non-positive with one negative; the History to stating those values; the embedding diagram to
the funnel's closed form and its finite height sqrt(pi) Gamma(1/4) a/(4 Gamma(3/4)); the curvature
figure to its values at the throat; and the drawings that cannot exist to their reasons.
The tests that read the mathematics need sympy and mpmath and are skipped where they are absent.
"""
import importlib.util
import json
import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "MFS" / "assets" / "data"
HAS_SYMPY = importlib.util.find_spec("sympy") is not None and importlib.util.find_spec("mpmath") is not None
sys.path.insert(0, str(Path(__file__).resolve().parent / "derivations"))
METRIC = json.loads((DATA / "metrics" / "giddings_strominger.json").read_text(encoding="utf-8"))
RADII = ("1.37", "2.21", "4.05")
HEIGHT = math.sqrt(math.pi) * math.gamma(0.25) / (4 * math.gamma(0.75))


def chart(system):
    return next(c for c in METRIC["coordinates"] if c["id"] == system)


def funnel(r):
    """The funnel's height above the throat at a = 1, by Simpson's rule on dz/du = 2/sqrt((2 + u^2)((1 + u^2)^2 + 1))
    with r = 1 + u^2, which has no singular end."""
    top, steps = math.sqrt(r - 1), 2000
    f = lambda u: 2 / math.sqrt((2 + u * u) * ((1 + u * u) ** 2 + 1))  # noqa: E731
    h = top / steps
    return h / 3 * (f(0) + f(top) + sum((4 if i % 2 else 2) * f(i * h) for i in range(1, steps)))


@unittest.skipUnless(HAS_SYMPY, "needs sympy and mpmath")
class Curvature(unittest.TestCase):
    """Every chart's published curvature at points, against the solution."""

    @classmethod
    def setUpClass(cls):
        import mpmath
        import sympy as sp
        import verify_metrics as vm
        mpmath.mp.dps = 40
        cls.sp, cls.mpmath, cls.vm = sp, mpmath, vm

    def reader(self, system):
        c = chart(system)
        return c, self.vm.Reader(c["coords"], [p["symbol"] for p in c["parameters"]], ())

    def diagonal(self, entry, field, variant):
        return {m["indices"][0]: m["value"] for m in entry[field]["variants"][variant]["nonzero"]
                if m["indices"][0] == m["indices"][1]}

    def published(self, system, r):
        """The published R, G^mu_mu and R^mu_mu of the chart at the point of areal radius r, a = 1."""
        sp, mpmath = self.sp, self.mpmath
        c, R = self.reader(system)
        x, chi, theta, phi = (R.symbol[name] for name in c["coords"])
        r = mpmath.mpf(r)
        point = {chi: sp.Rational(7, 10), theta: sp.Rational(2, 5), phi: sp.Rational(3, 10), R.parameters["a"]: 1}
        carried = {}
        if system == "areal":
            point[x] = sp.Float(r, 40)
        elif system == "proper":
            # The free r(rho) at a point of the wormhole, its slopes carried by (dr/drho)^2 = 1 - a^4/r^4.
            fn = R.parameters["r"]
            carried = {sp.Derivative(fn, (x, 2)): sp.Float(2 / r ** 5, 40),
                       sp.Derivative(fn, x): sp.Float(mpmath.sqrt(1 - 1 / r ** 4), 40), fn: sp.Float(r, 40)}
            point[x] = sp.Float(mpmath.quad(lambda s: 1 / mpmath.sqrt(1 - 1 / s ** 4), [1, r]), 40)
        else:
            point[x] = sp.Float(mpmath.acosh(r ** 2) / 2, 40)

        def at(value):
            return complex(sp.N(R(value).xreplace(carried).xreplace(point), 40)).real
        scalar = at(c["ricci_scalar"].split("=", 1)[1])
        einstein = {k: at(v) for k, v in self.diagonal(c, "einstein_tensor", "ul").items()}
        ricci = {k: at(v) for k, v in self.diagonal(c, "ricci_tensor", "ul").items()}
        return c, scalar, einstein, ricci

    def test_the_curvature_scalar_is_minus_6a4_over_r6_in_every_chart(self):
        self.assertEqual(chart("areal")["ricci_scalar"], "R = -\\dfrac{6a^4}{r^6}")
        for system in ("areal", "proper", "conformal"):
            for r in RADII:
                _, scalar, _, _ = self.published(system, r)
                self.assertLess(abs(scalar + 6 / float(r) ** 6), 1e-12, f"{system} at r = {r}")

    def test_the_einstein_tensor_is_the_axions_stress_in_every_chart(self):
        for system in ("areal", "proper", "conformal"):
            for r in RADII:
                c, _, einstein, _ = self.published(system, r)
                q = 1 / float(r) ** 6
                radial, *sphere = c["coords"]
                self.assertEqual(set(einstein), set(c["coords"]))
                self.assertLess(abs(einstein[radial] + 3 * q), 1e-12, f"{system} G along the radius at r = {r}")
                for angle in sphere:
                    self.assertLess(abs(einstein[angle] - 3 * q), 1e-12, f"{system} G^{angle} at r = {r}")

    def test_the_ricci_tensor_has_one_negative_eigenvalue_and_no_positive_one(self):
        for system in ("areal", "proper", "conformal"):
            for r in RADII:
                c, _, _, ricci = self.published(system, r)
                q = 1 / float(r) ** 6
                values = [ricci.get(name, 0.0) for name in c["coords"]]
                self.assertLess(abs(values[0] + 6 * q), 1e-12, f"{system} at r = {r}")
                self.assertTrue(all(abs(v) < 1e-12 for v in values[1:]), f"{system} at r = {r}")
                self.assertEqual(sum(v < -1e-12 for v in values), 1)
                self.assertFalse(any(v > 1e-12 for v in values))

    def test_the_space_is_conformally_flat_and_k_is_60a8_over_r12(self):
        for system in ("areal", "proper", "conformal"):
            self.assertEqual(chart(system)["weyl_tensor"]["variants"]["ulll"]["nonzero"], [])
        self.assertEqual(chart("areal")["kretschmann"], "K = \\dfrac{60a^8}{r^{12}}")
        self.assertEqual(chart("conformal")["kretschmann"], "K = \\dfrac{60}{a^4\\cosh^6(2\\eta)}")

    def test_the_proper_chart_leaves_the_radius_free(self):
        c = chart("proper")
        self.assertEqual([p["symbol"] for p in c["parameters"]], ["a", "r = r(\\rho)"])
        self.assertIn("no component assumes the equation $(\\partial_\\rho r)^2 = 1 - a^4/r^4$", c["convention"])


class History(unittest.TestCase):
    def test_the_history_states_the_curvature_the_published_charts_hold(self):
        history = METRIC["history"]
        self.assertIn("$R = -6a^4/r^6$", history)
        self.assertIn("$G^r{}_r = -3a^4/r^6$ along the radius and $3a^4/r^6$ along each direction on the sphere",
                      history)
        self.assertIn("$$ds^2 = \\frac{dr^2}{1 - a^4/r^4} + r^2\\,d\\Omega_3^2,$$", history)

    def test_the_history_states_the_funnels_height(self):
        self.assertIn("$\\sqrt{\\pi}\\,\\Gamma(1/4)\\,a/(4\\Gamma(3/4)) \\approx 1.311a$", METRIC["history"])
        self.assertEqual(f"{HEIGHT:.3f}", "1.311")


class Drawings(unittest.TestCase):
    def test_no_null_ray_or_conformal_diagram_is_written_and_both_say_why(self):
        self.assertFalse((DATA / "diagrams" / "giddings_strominger.json").exists())
        self.assertFalse((DATA / "conformal" / "giddings_strominger.json").exists())
        import giddings_strominger_drawings as drawings
        self.assertEqual(set(drawings.NOT_DRAWN), {"areal", "proper", "conformal"})
        self.assertTrue(all("Positive definite" in reason for reason in drawings.NOT_DRAWN.values()))

    def test_the_funnel_height_is_finite_and_in_closed_form(self):
        # z_inf = int_a^inf a^2 dr/sqrt(r^4 - a^4) = sqrt(pi) Gamma(1/4) a/(4 Gamma(3/4)) = 1.3110287771 a.
        self.assertLess(abs(HEIGHT - 1.3110287771), 1e-10)
        # What is left to climb beyond r is (a^2/r) 2F1(1/2, 1/4; 5/4; a^4/r^4), summed term by term.
        for r in (2.0, 4.0):
            x, term, tail = r ** -4, 1.0, 0.0
            for n in range(60):
                tail += term / (1 + 4 * n)
                term *= (0.5 + n) / (n + 1) * x
            self.assertLess(abs(funnel(r) + tail / r - HEIGHT), 1e-10, f"r = {r}")
        if HAS_SYMPY:
            import mpmath
            with mpmath.workdps(30):
                quadrature = mpmath.quad(lambda u: 2 / mpmath.sqrt((2 + u ** 2) * ((1 + u ** 2) ** 2 + 1)),
                                         [0, 1, mpmath.inf])
                closed = mpmath.sqrt(mpmath.pi) * mpmath.gamma(mpmath.mpf(1) / 4) / (4 * mpmath.gamma(mpmath.mpf(3) / 4))
                self.assertLess(abs(quadrature - closed), mpmath.mpf(10) ** -25)
        # Ellis and Bronnikov's catenoid, z = l arccosh(r/l), passes that height at r = 2l and climbs without bound.
        self.assertGreater(math.acosh(2.0), HEIGHT)

    def test_the_funnel_is_drawn_to_its_closed_form_through_the_throat(self):
        views = json.loads((DATA / "embedding" / "giddings_strominger.json").read_text(encoding="utf-8"))["views"]
        self.assertEqual([v["id"] for v in views], ["funnel"])
        pieces = {p["id"]: p["points"] for s in views[0]["surfaces"] for p in s["pieces"]}
        self.assertEqual(set(pieces), {"near", "far"})
        for name, sense in (("near", 1), ("far", -1)):
            points = pieces[name]
            self.assertEqual(points[0][:3], [1, 1, 0])
            self.assertEqual(points[-1][0], 4)
            for r, rho, z in points[::5] + points[-1:]:
                self.assertLess(abs(rho - r), 1e-6, f"{name}: rho at r = {r}")
                self.assertLess(abs(z - sense * funnel(r)), 1e-5, f"{name}: z at r = {r}")
                self.assertLess(abs(z), HEIGHT)

    def test_the_curvature_figure_peaks_at_the_throat(self):
        view, = json.loads((DATA / "curvature" / "giddings_strominger.json").read_text(encoding="utf-8"))["views"]
        self.assertEqual((view["x"]["from"], view["x"]["to"]), (-1.5, 1.5))
        for cls, extreme in (("K", 60.0), ("R", -6.0)):
            points = [p for curve in view["curves"] if curve["class"] == cls for p in curve["points"]]
            values = [v for _, v in points]
            peak = max(values, key=abs)
            self.assertLess(abs(peak - extreme), 1e-9 * abs(extreme), cls)
            self.assertLess(abs(min((abs(x) for x, v in points if v == peak))), 1e-3, cls)
            for x, v in points:
                # Even in eta: the stored scalar's value at -eta is the same, 60/cosh^6 or -6/cosh^3 of 2 eta.
                want = extreme / math.cosh(2 * x) ** (6 if cls == "K" else 3)
                self.assertLess(abs(v - want), 1e-9 * abs(extreme) + 1e-12, f"{cls} at eta = {x}")


if __name__ == "__main__":
    unittest.main()
