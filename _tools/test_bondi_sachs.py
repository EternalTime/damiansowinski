"""Bondi and Sachs's radiating metric, -(V/r) e^(2 beta) du^2 - 2 e^(2 beta) du dr
+ r^2 e^(2 gamma) (dtheta - U du)^2 + r^2 e^(-2 gamma) sin^2 theta dphi^2:
python3 -m unittest discover -s _tools

What its texts and drawings state, held on the published files. The first class needs nothing
but the files: the cones of the spacetime diagrams, the edge that is future null infinity, the
conformal diagram's cones ending on it, the sphere of the embedding diagram, and the relations each side
answers. The second reads
the published tensors through the checker's Reader and holds the burst the drawings declare to
them: Schwarzschild before it, the vacuum equations to the order Bondi's expansion keeps, the
mass loss formula, and the numbers the captions state; it needs sympy and is skipped where it
is absent.
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

CHARTS = ["bondi", "compactified"]


def load(folder, metric_id="bondi_sachs"):
    return json.loads((DATA / folder / f"{metric_id}.json").read_text(encoding="utf-8"))


class Drawings(unittest.TestCase):
    def test_each_chart_has_its_spacetime_diagrams(self):
        systems = load("diagrams")["systems"]
        self.assertEqual(sorted(systems), sorted(CHARTS))
        self.assertEqual([v["id"] for v in systems["bondi"]], ["equator", "axis"])
        self.assertEqual([v["id"] for v in systems["compactified"]], ["equator"])

    def test_the_drawn_windows_are_between_one_to_two_and_two_to_one(self):
        for system, views in load("diagrams")["systems"].items():
            for view in views:
                X0, X1, Y0, Y1 = view["box"]
                ratio = (X1 - X0) / (Y1 - Y0)
                self.assertTrue(0.5 <= ratio <= 2.0, (system, view["id"], ratio))

    def test_the_outgoing_rays_of_bondis_chart_run_at_45_degrees(self):
        """u constant is the line Y = X + cu on the axes r and cu + r."""
        for view in load("diagrams")["systems"]["bondi"]:
            X0, X1, Y0, Y1 = view["box"]
            for cone in view["cones"]:
                slopes = [b * (Y1 - Y0) / (a * (X1 - X0)) for a, b in (cone["a"], cone["b"])]
                self.assertEqual(sum(abs(s - 1) < 2e-3 for s in slopes), 1, (view["id"], slopes))

    def test_the_burst_is_marked_on_its_first_and_last_cone(self):
        """u = 0 and cu = 20 m_0: lines of slope one on Bondi's axes, level lines on the chart of 1/r."""
        systems = load("diagrams")["systems"]
        for view in systems["bondi"]:
            X0, X1, Y0, Y1 = view["box"]
            got = []
            for mark in view["markers"]:
                (a, b), = mark["lines"]
                at = [(X0 + x * (X1 - X0), Y0 + y * (Y1 - Y0)) for x, y in (a, b)]
                got.append([round(Y - X, 6) for X, Y in at])
            self.assertEqual(got, [[0.0, 0.0], [20.0, 20.0]], view["id"])
        (view,) = systems["compactified"]
        Y0, Y1 = view["box"][2:]
        levels = [[round(Y0 + y * (Y1 - Y0), 3) for _, y in mark["lines"][0]] for mark in view["markers"]]
        for got, want in zip(levels, ([0.0, 0.0], [20.0, 20.0])):
            for a, b in zip(got, want):
                self.assertAlmostEqual(a, b, delta=0.01)

    def test_the_future_cone_closes_onto_null_infinity(self):
        """In the chart of 1/r the outgoing ray is level, and the ingoing one, dl/d(cu) = l^3 V/2,
        stands more nearly upright the nearer the edge l = 0."""
        (view,) = load("diagrams")["systems"]["compactified"]
        lean = {}
        for cone in view["cones"]:
            flat = [v for v in (cone["a"], cone["b"]) if abs(v[1]) < 1e-9]
            steep = [v for v in (cone["a"], cone["b"]) if abs(v[1]) > 0.5]
            self.assertEqual((len(flat), len(steep)), (1, 1), cone)
            self.assertLess(flat[0][0], 0)          # toward smaller l, outward
            lean.setdefault(round(cone["at"][1], 3), []).append((cone["at"][0], steep[0][0] / steep[0][1]))
        for row in lean.values():
            row.sort()
            self.assertLess(row[0][1], row[-1][1])
            self.assertLess(abs(row[0][1]), 0.02)

    def test_the_cones_of_the_conformal_diagram_end_on_future_null_infinity(self):
        """Every cone drawn, the two of the burst among them, runs from the world tube to q = pi/2,
        the line T = pi - X."""
        (view,) = load("conformal")["views"]
        self.assertEqual(view["id"], "axis")
        cones = [layer["points"] for layer in view["layers"] if layer["class"] in ("null", "surface")]
        self.assertEqual(len(cones), 7)
        for points in cones:
            X, T = points[-1]
            self.assertAlmostEqual(X + T, math.pi, delta=2e-3)
            self.assertGreater(points[0][0], 0.05)      # it starts on the tube, off the left edge
        tube = next(layer["points"] for layer in view["layers"] if layer["class"] == "boundary")
        self.assertAlmostEqual(tube[0][1], -math.pi, delta=0.05)
        self.assertAlmostEqual(tube[-1][1], math.pi, delta=0.05)

    def test_the_sphere_of_the_world_tube_is_oblate_by_half_a_percent(self):
        """At cu = 10 m_0 and r = 10 m_0 the equator has the radius 1.0035 r and the surface keeps the
        area 4 pi r^2, as the caption states."""
        (view,) = load("embedding")["views"]
        self.assertEqual(view["id"], "sphere")
        (surface,) = view["surfaces"]
        (piece,) = surface["pieces"]
        points = piece["points"]
        self.assertAlmostEqual(max(p[1] for p in points) / 10, 1.0035, places=3)
        area = sum(math.pi * (a[1] + b[1]) * math.hypot(b[1] - a[1], b[2] - a[2]) for a, b in zip(points, points[1:]))
        self.assertAlmostEqual(area / (400 * math.pi), 1.0, places=3)
        height = max(p[2] for p in points) - min(p[2] for p in points)
        self.assertLess(height, 2 * 10 * 0.999)        # flattened along the axis

    def test_each_relation_is_answered_with_the_kind_that_goes_with_it(self):
        answers = {"special_case": "generalisation", "programme": "programme"}
        related = load("metrics")["related"]
        self.assertEqual([r["id"] for r in related],
                         ["schwarzschild", "vaidya", "robinson_trautman", "photon_rocket", "minkowski", "kerr", "pp_wave"])
        for entry in related:
            back = [r for r in load("metrics", entry["id"])["related"] if r["id"] == "bondi_sachs"]
            self.assertEqual([r["kind"] for r in back], [answers[entry["kind"]]], entry["id"])


@unittest.skipUnless(HAS_SYMPY, "needs sympy")
class Burst(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import mpmath
        import sympy as sp
        import null_rays as nr
        cls.sp, cls.nr, cls.mp = sp, nr, mpmath
        _, cls.entry, cls.reader = nr.load("bondi_sachs", "bondi")
        R = cls.reader
        cls.u, cls.r, cls.theta = (R.symbol[n] for n in ("u", "r", "\\theta"))
        names = {"u": cls.u, "r": cls.r, "theta": cls.theta}
        cls.member = {k: sp.sympify(v, locals=names) for k, v in nr.BONDI_BURST.items()}
        cls.facts = {k: sp.sympify(str(nr._BURST[k]), locals=names) for k in ("sigma", "M")}

    def on_burst(self, text):
        """A published value with the declared V, beta, U and gamma put in, as a function of
        (u, r, theta) in thirty digits."""
        sp, R = self.sp, self.reader
        e = R(text).subs(R.c, 1)
        for name, rep in self.member.items():
            fn = R.parameters[name]
            e = e.replace(fn.func, sp.Lambda(fn.args, rep))
        return sp.lambdify((self.u, self.r, self.theta), e.doit(), "mpmath")

    def test_the_charts_are_the_two_the_literature_uses(self):
        self.assertEqual([c["id"] for c in load("metrics")["coordinates"]], CHARTS)

    def test_before_the_burst_the_member_is_schwarzschilds(self):
        sp = self.sp
        before = {k: v.subs(self.u, -1) for k, v in self.member.items()}
        self.assertEqual(sp.simplify(before["V"] - (self.r - 2)), 0)
        for k in ("beta", "U", "gamma"):
            self.assertEqual(before[k], 0, k)

    def test_the_expansion_solves_the_vacuum_equations_to_the_order_it_keeps(self):
        """Every published Ricci component, in an orthonormal measure, against the square root of
        the Kretschmann scalar: under 1.5% at r = 10 m_0, as the drawings' input states, and
        falling as 1/r, by a factor between 1.6 and 2.4 for each doubling of r."""
        mp = self.mp
        mp.mp.dps = 30
        ricci = {tuple(c["indices"]): self.on_burst(c["value"])
                 for c in self.entry["ricci_tensor"]["variants"]["ll"]["nonzero"]}
        K = self.on_burst(self.entry["kretschmann"].split("=", 1)[1])
        worst = {}
        for r in (10, 20, 40):
            worst[r] = 0.0
            for u in (-4.0, 3.0, 7.0, 10.0, 13.0, 17.0, 24.0):
                for theta in (0.3, math.pi / 2 - 1e-9):
                    scale = {"u": 1.0, "r": 1.0, "\\theta": r, "\\phi": r * math.sin(theta)}
                    size = float(mp.sqrt(abs(K(u, r, theta))))
                    for index, f in ricci.items():
                        value = abs(float(f(mp.mpf(u), mp.mpf(r), mp.mpf(theta)))) / (scale[index[0]] * scale[index[1]])
                        if u < 0:
                            self.assertLess(value, 1e-25, (index, r, theta))
                        worst[r] = max(worst[r], value / size)
        self.assertLess(worst[10], 0.015)
        self.assertGreater(worst[10], 0.010)
        for r in (10, 20):
            self.assertTrue(1.6 < worst[r] / worst[2 * r] < 2.4, worst)

    def test_the_bondi_mass_falls_by_the_integral_of_the_news_squared(self):
        """dm/du = -(1/2) int (d_u sigma)^2 sin(theta) dtheta: the mean of the mass aspect over the
        sphere after the burst is m_0 less 8/15 of the integral of the news on the equator squared."""
        sp, nr = self.sp, self.nr
        x = sp.Symbol("x")
        M_after = self.facts["M"].subs(self.u, 25)
        mean = sp.integrate(M_after * sp.sin(self.theta), (self.theta, 0, sp.pi)) / 2
        news = sp.diff(self.facts["sigma"], self.u)
        during = news.subs(self.u, x).subs(self.theta, sp.pi / 2)
        radiated = sp.Rational(8, 15) * sp.Integral(during ** 2, (x, 0, nr.BURST["T"]))
        self.assertAlmostEqual(float(1 - mean), float(radiated.evalf(20)), places=12)
        self.assertAlmostEqual(float(1 - mean), float(nr.BURST_FACTS["radiated"]), places=12)
        self.assertEqual(f"{float(mean):.4f}", "0.9990")
        # The loss at each angle goes as sin^4, the square of the news there.
        flux = sp.integrate(news ** 2, (self.u, 0, nr.BURST["T"]))
        self.assertEqual(sp.simplify(flux - flux.subs(self.theta, sp.pi / 2) * sp.sin(self.theta) ** 4), 0)

    def test_the_numbers_the_captions_state(self):
        sp = self.sp
        M = sp.lambdify((self.u, self.theta), self.facts["M"], "math")
        sigma = sp.lambdify((self.u, self.theta), self.facts["sigma"], "math")
        news = sp.lambdify((self.u, self.theta), sp.diff(self.facts["sigma"], self.u), "math")
        times = [20 * (i + 0.5) / 4000 for i in range(4000)]
        equator, axis = math.pi / 2, 0.0
        self.assertEqual(f"{M(25.0, equator):.4f}", "0.9981")
        self.assertEqual(f"{M(25.0, axis):.4f}", "1.0000")
        on_axis = [M(t, axis) for t in times]
        self.assertEqual((f"{min(on_axis):.2f}", f"{max(on_axis):.2f}"), ("0.85", "1.08"))
        shear = [sigma(t, equator) for t in times]
        self.assertEqual((f"{min(shear):.3f}", f"{max(shear):.3f}"), ("-0.037", "0.019"))
        self.assertEqual(f"{max(abs(news(t, equator)) for t in times):.3f}", "0.019")
        self.assertEqual(f"{1000 * max(abs(s) for s in shear) / 10:.0f}", "4")      # parts in a thousand at r = 10
        self.assertTrue(all(sigma(t, axis) == 0 for t in times[::400]))

    def test_the_chart_of_the_inverse_distance_draws_the_same_burst(self):
        sp, nr = self.sp, self.nr
        ell = sp.Symbol("ell", real=True)
        names = {"u": self.u, "ell": ell, "theta": self.theta}
        for k, text in nr.BONDI_BURST_INVERSE.items():
            there = sp.sympify(text, locals=names).subs(ell, 1 / self.r)
            point = {self.u: sp.Rational(7), self.r: 13, self.theta: sp.Rational(9, 10)}
            self.assertAlmostEqual(float(there.subs(point)), float(self.member[k].subs(point)), places=14, msg=k)


if __name__ == "__main__":
    unittest.main()
