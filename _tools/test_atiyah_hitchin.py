"""The Atiyah-Hitchin space, a Riemannian space of four dimensions with no time in it:
.venv.noindex/bin/python -m unittest discover -s _tools

What its charts, texts and drawings state, held on the published files. The modulus chart's
coefficients are held to the three products its convention states and to Weinberg and Yi's
products of Gibbons and Manton's solution, evaluated here from mpmath's elliptic integrals and
nothing of the chart; the core to a 2-sphere of radius pi n; the far region to Taub-NUT with a
negative mass; the published Riemann tensor to being self dual up to orientation; the radius
chart to being the modulus chart carried along r = 2nK; the embedding diagram's surfaces to the
lengths a/2, b/2 and |c|/2 at their radii; and the drawings that cannot exist to their reasons.
The tests that read the mathematics need sympy and mpmath and are skipped where they are absent.
"""
import importlib.util
import itertools
import json
import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "MFS" / "assets" / "data"
HAS_SYMPY = importlib.util.find_spec("sympy") is not None and importlib.util.find_spec("mpmath") is not None
sys.path.insert(0, str(Path(__file__).resolve().parent / "derivations"))
METRIC = json.loads((DATA / "metrics" / "atiyah_hitchin.json").read_text(encoding="utf-8"))


def chart(system):
    return next(c for c in METRIC["coordinates"] if c["id"] == system)


def products(k):
    """Weinberg and Yi's ab, bc and ca at n = 1, in Gibbons and Manton's angle beta = 2 arcsin k,
    with r = 2K(sin(beta/2)) and its slope taken numerically."""
    import mpmath
    beta = 2 * mpmath.asin(k)
    radius = lambda b: 2 * mpmath.ellipk(mpmath.sin(b / 2) ** 2)  # noqa: E731
    r, slope = radius(beta), mpmath.diff(radius, beta)
    s, c = mpmath.sin(beta), mpmath.cos(beta)
    return -s * r * slope + (1 - c) * r ** 2 / 2, -s * r * slope - (1 + c) * r ** 2 / 2, -s * r * slope


@unittest.skipUnless(HAS_SYMPY, "needs sympy and mpmath")
class Modulus(unittest.TestCase):
    """The modulus chart, read through the checker's own Reader at points."""

    @classmethod
    def setUpClass(cls):
        import mpmath
        import sympy as sp
        import verify_metrics as vm
        mpmath.mp.dps = 40
        cls.sp, cls.mpmath = sp, mpmath
        c = chart("modulus")
        cls.entry = c
        cls.reader = vm.Reader(c["coords"], [p["symbol"] for p in c["parameters"]], (),
                               held=vm.HELD[("atiyah_hitchin", "modulus")], rates=vm.RATES[("atiyah_hitchin", "modulus")])
        R = cls.reader
        cls.k, cls.theta, cls.phi, cls.psi = (R.symbol[x] for x in c["coords"])
        cls.K, cls.E, cls.n = R.parameters["K"], R.parameters["E"], R.parameters["n"]
        cls.g = {tuple(m["indices"]): R(m["value"]) for m in c["metric_components"]}

    def at(self, value, k, theta=0.7, phi=0.4, psi=0.0):
        sp, mpmath = self.sp, self.mpmath
        k = mpmath.mpf(k)
        names = {self.K: sp.Float(mpmath.ellipk(k ** 2), 40), self.E: sp.Float(mpmath.ellipe(k ** 2), 40)}
        point = {self.k: sp.Float(k, 40), self.theta: sp.Float(theta, 40), self.phi: sp.Float(phi, 40),
                 self.psi: sp.Float(psi, 40), self.n: 1}
        return complex(sp.N(sp.sympify(value).xreplace(names).xreplace(point), 40)).real

    def test_the_coefficients_are_gibbons_and_mantons_products(self):
        # At psi = 0 the three forms are -cos(theta)... : g_thth = b^2, g_psipsi = c^2, and
        # g_phiphi - c^2 cos^2(theta) = a^2 sin^2(theta).
        for k in ("0.2", "0.55", "0.9", "0.99"):
            ab, bc, ca = products(self.mpmath.mpf(k))
            a2, b2, c2 = ab * ca / bc, ab * bc / ca, bc * ca / ab
            theta = 0.7
            got_b2 = self.at(self.g[("\\theta", "\\theta")], k, theta)
            got_c2 = self.at(self.g[("\\psi", "\\psi")], k, theta)
            got_a2 = (self.at(self.g[("\\phi", "\\phi")], k, theta) - got_c2 * math.cos(theta) ** 2) / math.sin(theta) ** 2
            for name, got, want in (("a^2", got_a2, a2), ("b^2", got_b2, b2), ("c^2", got_c2, c2)):
                # The products take dr/dbeta numerically, which limits them to about a part in 10^12.
                self.assertLess(abs(got / float(want) - 1), 1e-9, f"{name} at k = {k}")

    def test_the_convention_states_the_products_the_line_element_holds(self):
        convention = self.entry["convention"]
        for text in ("bc = -4n^2KE", "ca = -4n^2K(E - (1 - k^2)K)", "ab = 4n^2K(K - E)", "r = 2nK"):
            self.assertIn(text, convention)

    def test_the_core_is_a_sphere_of_radius_pi_n(self):
        # b and -c tend to pi n as k -> 0, the radius of the 2-sphere the circle of sigma_1 closes on,
        # b as pi n + (r - pi n)/2 and -c as pi n - (r - pi n)/2, with r - pi n = pi n k^2/4 to first order.
        for k in (0.01, 0.001):
            b = math.sqrt(self.at(self.g[("\\theta", "\\theta")], k))
            c = math.sqrt(self.at(self.g[("\\psi", "\\psi")], k))
            self.assertLess(abs(b / math.pi - 1 - k * k / 8), 10 * k ** 4)
            self.assertLess(abs(c / math.pi - 1 + k * k / 8), 10 * k ** 4)
        self.assertIn("k = 0 \\;\\text{(the core, a real projective plane)}", self.entry["domains"])
        self.assertIn("\\psi \\in [0, \\pi)", self.entry["domains"])

    def test_far_away_it_is_taub_nut_with_a_negative_mass(self):
        # r = 2nK; a^2 and b^2 tend to r^2(1 - 2n/r) and c^2 to 4n^2/(1 - 2n/r), to terms that
        # fall off as e^(-r/n): at n = 1, Tong's b^2 = r^2(1 - 2/r) + 8r^3 e^(-r), 3.3e-7 of r^2 at r = 20.
        mpmath = self.mpmath
        k = mpmath.sqrt(1 - 16 * mpmath.exp(-20))   # where K(k^2) = 10 to leading order
        r = 2 * mpmath.ellipk(k ** 2)
        self.assertGreater(r, 19.9)
        c2 = self.at(self.g[("\\psi", "\\psi")], k)
        b2 = self.at(self.g[("\\theta", "\\theta")], k)
        self.assertLess(abs(c2 * (1 - 2 / float(r)) / 4 - 1), 1e-6)
        r = float(r)
        self.assertLess(abs((b2 - r ** 2 * (1 - 2 / r)) / (8 * r ** 3 * math.exp(-r)) - 1), 0.15)
        # The mass is negative: c^2 exceeds its far value 4n^2, as 4n^2/(1 + 2m/r) would with m = -n.
        self.assertGreater(c2, 4)

    def test_the_riemann_tensor_is_self_dual_up_to_orientation(self):
        sp = self.sp
        R = self.reader
        coords = [R.symbol[x] for x in self.entry["coords"]]
        names = self.entry["coords"]
        riemann = {}
        for m in self.entry["riemann"]["variants"]["llll"]["nonzero"]:
            riemann[tuple(names.index(i) for i in m["indices"])] = R(m["value"])
        g = sp.Matrix(4, 4, lambda i, j: self.g.get((names[i], names[j]), 0))
        k, theta, psi = "0.6", 0.9, 0.5
        G = sp.Matrix(4, 4, lambda i, j: self.at(g[i, j], k, theta, 0.4, psi))
        Ginv, root = G.inv(), math.sqrt(float(G.det()))
        low = {key: self.at(v, k, theta, 0.4, psi) for key, v in riemann.items()}

        def Rl(a, b, c, d):
            return low.get((a, b, c, d), 0.0)
        signs = []
        for a, b in itertools.combinations(range(4), 2):
            for c, d in itertools.combinations(range(4), 2):
                left = Rl(a, b, c, d)
                # (*R)_abcd = (1/2) eps_ab^ef R_efcd, with eps_abef = sqrt(g) [abef].
                dual = 0.0
                for e, f in itertools.permutations(range(4), 2):
                    for p, q in itertools.permutations(range(4), 2):
                        sign = sp.LeviCivita(a, b, p, q)
                        if sign:
                            dual += 0.5 * float(sign) * root * float(Ginv[p, e]) * float(Ginv[q, f]) * Rl(e, f, c, d)
                if abs(left) > 1e-8:
                    signs.append(dual / left)
                else:
                    self.assertLess(abs(dual), 1e-8)
        self.assertTrue(signs)
        orientation = round(signs[0])
        self.assertIn(orientation, (1, -1))
        for s in signs:
            self.assertLess(abs(s - orientation), 1e-8)

    def test_the_ricci_tensor_and_scalar_vanish(self):
        self.assertEqual(self.entry["ricci_tensor"]["variants"]["ll"]["nonzero"], [])
        self.assertEqual(self.entry["ricci_scalar"], "R = 0")


@unittest.skipUnless(HAS_SYMPY, "needs sympy and mpmath")
class Radius(unittest.TestCase):
    def test_the_radius_chart_states_the_equations_and_leaves_a_b_and_c_free(self):
        c = chart("radius")
        self.assertEqual([p["symbol"] for p in c["parameters"]], ["n", "a = a(r)", "b = b(r)", "c = c(r)"])
        for text in ("2rc\\,a' = a^2 - (b - c)^2", "2rac\\,b' = b\\left(b^2 - (c - a)^2\\right)",
                     "2rab\\,c' = b\\left(c^2 - (a - b)^2\\right)", "no component assumes"):
            self.assertIn(text, c["convention"])

    def test_the_equations_hold_for_gibbons_and_mantons_solution(self):
        # a, b and c of r from the products at k, with r = 2K(k^2), against the stated rates.
        import mpmath
        mpmath.mp.dps = 40

        def lengths(k):
            ab, bc, ca = products(k)
            return mpmath.sqrt(ab * ca / bc), mpmath.sqrt(ab * bc / ca), -mpmath.sqrt(bc * ca / ab)
        for k in ("0.3", "0.7", "0.95"):
            k = mpmath.mpf(k)
            a, b, c = lengths(k)
            r = 2 * mpmath.ellipk(k ** 2)
            drdk = mpmath.diff(lambda q: 2 * mpmath.ellipk(q ** 2), k)
            da, db, dc = (mpmath.diff(lambda q, i=i: lengths(q)[i], k) / drdk for i in range(3))
            self.assertLess(abs(2 * r * c * da - (a ** 2 - (b - c) ** 2)), 1e-25)
            self.assertLess(abs(2 * r * a * c * db - b * (b ** 2 - (c - a) ** 2)), 1e-25)
            self.assertLess(abs(2 * r * a * b * dc - b * (c ** 2 - (a - b) ** 2)), 1e-25)


class Drawings(unittest.TestCase):
    def test_no_null_ray_or_conformal_diagram_is_written_and_both_say_why(self):
        self.assertFalse((DATA / "diagrams" / "atiyah_hitchin.json").exists())
        self.assertFalse((DATA / "conformal" / "atiyah_hitchin.json").exists())
        import atiyah_hitchin_drawings as drawings
        self.assertEqual(set(drawings.NOT_DRAWN), {"modulus", "radius"})
        self.assertTrue(all("Positive definite" in reason for reason in drawings.NOT_DRAWN.values()))

    def test_the_embedding_surfaces_have_the_lengths_of_the_solution(self):
        views = {v["id"]: v for v in json.loads((DATA / "embedding" / "atiyah_hitchin.json")
                                                 .read_text(encoding="utf-8"))["views"]}
        self.assertEqual(set(views), {"cone", "trumpet"})
        if not HAS_SYMPY:
            return
        import atiyah_hitchin_drawings as drawings

        def pieces(view):
            return {p["id"]: [tuple(q) for q in p["points"]] for s in views[view]["surfaces"] for p in s["pieces"]}
        cone = pieces("cone")["cone"]
        self.assertEqual(cone[0][1], 0.0)
        for r, rho, _ in cone[1::7]:
            self.assertLess(abs(rho - float(drawings.solution(r)[0][0]) / 2), 1e-5, f"the cone's rho at r = {r}")
        trumpet = pieces("trumpet")
        for piece, index in (("tube", 2), ("flare", 1)):
            for r, rho, _ in trumpet[piece][1::7]:
                self.assertLess(abs(rho - abs(float(drawings.solution(r)[0][index])) / 2), 1e-5,
                                f"the trumpet's {piece} at r = {r}")
        # The two halves meet on the core's circle, rho = pi n/2.
        self.assertLess(abs(trumpet["tube"][0][1] - math.pi / 2), 1e-6)
        self.assertLess(abs(trumpet["flare"][0][1] - math.pi / 2), 1e-6)

    def test_the_curvature_figure_falls_from_the_core(self):
        figure = json.loads((DATA / "curvature" / "atiyah_hitchin.json").read_text(encoding="utf-8"))
        points = [p for curve in figure["views"][0]["curves"] for p in curve["points"]]
        values = [v for _, v in points]
        self.assertTrue(all(v > 0 for v in values))
        self.assertGreater(values[0], 0.55)
        self.assertTrue(all(b <= a * (1 + 1e-12) for a, b in zip(values, values[1:])))


if __name__ == "__main__":
    unittest.main()
