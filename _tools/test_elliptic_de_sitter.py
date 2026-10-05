"""Elliptic de Sitter space: .venv.noindex/bin/python -m unittest discover -s _tools

What its texts state, held on the hyperboloid -X_0^2 + X_1^2 + ... + X_4^2 = l^2 and on the published
files. An event and its antipode -X are spacelike separated and share no light ray, so the
identification makes no closed timelike curve; the antipodal map reverses the direction of time, so
the quotient has none; any two events of the quotient are joined by a geodesic, which fails on the
hyperboloid; and an observer sees exactly one of every antipodal pair, Schrodinger's theorem. Every
published chart is the hyperboloid pulled back, the identification each chart states is X -> -X, the
static patch and the planar chart hold no antipodal pair, and de Sitter's elliptical space, opposite
points of the sphere of his static chart made one, is the same identification. The drawings hold
what their captions say. The tests of the published mathematics need sympy and are skipped where it
is absent; the others read the files alone or need nothing.
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

HALF = math.pi / 2


def dot(X, Y):
    """The product of five dimensional Minkowski space, time first."""
    return -X[0] * Y[0] + sum(a * b for a, b in zip(X[1:], Y[1:]))


def event(rng, reach=2.0):
    """A random event of the hyperboloid at l = 1, in its global chart."""
    t = rng.uniform(-reach, reach)
    n = [rng.gauss(0, 1) for _ in range(4)]
    size = math.sqrt(sum(x * x for x in n))
    return [math.sinh(t)] + [math.cosh(t) * x / size for x in n]


def tangent(rng, X, kind):
    """A random tangent vector at X: timelike with V.V = -1, spacelike with V.V = 1, or null."""
    while True:
        V = [rng.gauss(0, 1) for _ in range(5)]
        k = dot(V, X)
        V = [v - k * x for v, x in zip(V, X)]
        norm = dot(V, V)
        if kind == "timelike" and norm < -1e-3:
            return [v / math.sqrt(-norm) for v in V]
        if kind == "spacelike" and norm > 1e-3:
            return [v / math.sqrt(norm) for v in V]
        if kind == "null" and norm > 1e-3:
            # A unit spacelike vector plus the unit timelike one orthogonal to it and to X.
            S = [v / math.sqrt(norm) for v in V]
            T = tangent(rng, X, "timelike")
            k = dot(T, S)
            T = [a - k * b for a, b in zip(T, S)]
            size = math.sqrt(-dot(T, T))
            return [a / size + b for a, b in zip(T, S)]


class Hyperboloid(unittest.TestCase):
    def setUp(self):
        self.rng = random.Random(1956)

    def test_an_event_and_its_antipode_are_spacelike_separated_and_share_no_light_ray(self):
        # (X - (-X))^2 = 4 l^2 > 0. A point Y = X + s N on the light cone of X, with N null and tangent,
        # has X.Y = l^2, and one on the light cone of -X has X.Y = -l^2, so no event lies on both.
        for _ in range(200):
            X = event(self.rng)
            anti = [-x for x in X]
            gap = [a - b for a, b in zip(X, anti)]
            self.assertAlmostEqual(dot(gap, gap), 4.0, places=9)
            N = tangent(self.rng, X, "null")
            self.assertAlmostEqual(dot(N, N), 0.0, places=9)
            s = self.rng.uniform(-3, 3)
            Y = [x + s * n for x, n in zip(X, N)]
            self.assertAlmostEqual(dot(Y, Y), 1.0, places=8)
            self.assertAlmostEqual(dot(X, Y), 1.0, places=8)
            to_antipode = [a - b for a, b in zip(Y, anti)]
            self.assertAlmostEqual(dot(to_antipode, to_antipode), 4.0, places=7)

    def test_the_antipodal_map_reverses_the_direction_of_time(self):
        # On the hyperboloid a timelike vector is future directed when its X_0 component is positive,
        # and X -> -X carries the vector V at X to -V at -X: tangent there, timelike, and past directed.
        for _ in range(200):
            X = event(self.rng)
            V = tangent(self.rng, X, "timelike")
            if V[0] < 0:
                V = [-v for v in V]
            image_at, image = [-x for x in X], [-v for v in V]
            self.assertAlmostEqual(dot(image, image_at), 0.0, places=9)
            self.assertAlmostEqual(dot(image, image), -1.0, places=9)
            self.assertGreater(V[0], 0)
            self.assertLess(image[0], 0)

    def joined(self, X, Y):
        """Whether a geodesic of the hyperboloid joins X to Y, built and followed to its end. The
        geodesics through X are X cos s + V sin s, X cosh s + V sinh s and X + s V for a unit
        spacelike, unit timelike and null tangent V, on which X.Y is cos s, cosh s and 1."""
        Z = dot(X, Y)
        if Z <= -1 + 1e-12:
            return max(abs(a + b) for a, b in zip(X, Y)) < 1e-9
        V = [y - Z * x for x, y in zip(X, Y)]
        norm = dot(V, V)
        if abs(Z) < 1:
            s = math.acos(Z)
            V = [v / math.sqrt(norm) for v in V]
            end = [x * math.cos(s) + v * math.sin(s) for x, v in zip(X, V)]
        elif Z > 1:
            s = math.acosh(Z)
            V = [v / math.sqrt(-norm) for v in V]
            end = [x * math.cosh(s) + v * math.sinh(s) for x, v in zip(X, V)]
        else:
            end = [x + v for x, v in zip(X, V)]
        scale = 1 + max(abs(y) for y in Y)
        return max(abs(a - b) for a, b in zip(end, Y)) < 1e-7 * scale

    def test_any_two_events_of_the_quotient_are_joined_by_a_geodesic_and_not_of_the_hyperboloid(self):
        # Calabi and Markus, as Hawking and Ellis state it: X.Y < -l^2 leaves no geodesic from X to Y, and
        # then X.(-Y) > l^2, a timelike geodesic to the antipode of Y, which is the same event.
        alone = 0
        for _ in range(400):
            X, Y = event(self.rng), event(self.rng)
            direct = self.joined(X, Y)
            through = self.joined(X, [-y for y in Y])
            self.assertTrue(direct or through)
            self.assertEqual(direct, dot(X, Y) > -1)
            alone += not direct
        self.assertGreater(alone, 20)

    def test_an_observer_sees_exactly_one_of_every_antipodal_pair(self):
        # The observer at the pole, X = (sinh tau, 0, 0, 0, cosh tau). An event Y is in the past of the
        # world line when some later point of it has X.Y >= l^2, which happens exactly when Y_4 > Y_0.
        def seen(Y):
            for k in range(1, 800):
                tau = -20 + 0.05 * k
                if math.sinh(tau) > Y[0] and -math.sinh(tau) * Y[0] + math.cosh(tau) * Y[4] >= 1:
                    return True
            return False
        for _ in range(200):
            Y = event(self.rng)
            if abs(Y[4] - Y[0]) < 1e-3:
                continue
            anti = [-y for y in Y]
            self.assertEqual(seen(Y), Y[4] > Y[0])
            self.assertNotEqual(seen(Y), seen(anti))

    def test_de_sitters_elliptical_space_is_this_identification(self):
        # His static chart with r = l sin chi: X_0 = l cos chi sinh(ct/l), X_4 = l cos chi cosh(ct/l) and
        # X_i = l sin chi n_i. Opposite points of the sphere of chi and n at one t, chi -> pi - chi and
        # n -> -n, are X and -X.
        for _ in range(200):
            t, chi = self.rng.uniform(-2, 2), self.rng.uniform(0, math.pi)
            n = [self.rng.gauss(0, 1) for _ in range(3)]
            size = math.sqrt(sum(x * x for x in n))
            n = [x / size for x in n]

            def at(chi, n):
                return ([math.cos(chi) * math.sinh(t)] + [math.sin(chi) * x for x in n]
                        + [math.cos(chi) * math.cosh(t)])
            X, opposite = at(chi, n), at(math.pi - chi, [-x for x in n])
            self.assertAlmostEqual(dot(X, X), 1.0, places=9)
            for a, b in zip(X, opposite):
                self.assertAlmostEqual(a, -b, places=9)


def published():
    return json.loads((DATA / "metrics" / "elliptic_de_sitter.json").read_text(encoding="utf-8"))


@unittest.skipUnless(HAS_SYMPY, "sympy is not installed")
class Published(unittest.TestCase):
    def setUp(self):
        import sympy
        import print_charts
        import verify_metrics
        self.sp, self.pc, self.vm = sympy, print_charts, verify_metrics
        self.rng = random.Random(2003)

    def chart(self, system):
        """The published metric of a chart as a matrix in the chart x^0 = ct at l = 1, its coordinate
        symbols, and its embedding in the hyperboloid."""
        entry = next(c for c in published()["coordinates"] if c["id"] == system)
        reader = self.vm.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]], ())
        x = [reader.symbol[n] for n in entry["coords"]]
        ell = reader.parameters["ell"]
        g = self.sp.zeros(4, 4)
        for e in entry["metric_components"]:
            i, j = (entry["coords"].index(k) for k in e["indices"])
            g[i, j] = reader(e["value"]).subs({reader.c: 1, ell: 1})
        X = [self.sp.sympify(f).subs(ell, 1) for f in self.pc.eds_embedding(system, x, ell)]
        return entry, x, g, X

    def point(self, system):
        r = self.rng.uniform
        return {"global": lambda: [r(-1.5, 1.5), r(0.1, 1.5), r(0.2, 2.9), r(0, 6)],
                "conformal": lambda: [r(-1.3, 1.3), r(0.1, 1.5), r(0.2, 2.9), r(0, 6)],
                "kruskal": lambda: [r(-0.9, 0.9), r(-0.9, 0.9), r(0.2, 2.9), r(0, 6)],
                "static": lambda: [r(-1.5, 1.5), r(0.05, 0.95), r(0.2, 2.9), r(0, 6)],
                "planar": lambda: [r(-1.5, 1.5), r(-2, 2), r(-2, 2), r(-2, 2)]}[system]()

    def test_every_chart_is_the_hyperboloid_pulled_back(self):
        sp = self.sp
        for system in self.pc.EDS_CHARTS:
            entry, x, g, X = self.chart(system)
            J = sp.Matrix([[sp.diff(f, v) for v in x] for f in X])
            pulled = J.T * sp.diag(-1, 1, 1, 1, 1) * J
            on = -X[0] ** 2 + sum(f ** 2 for f in X[1:]) - 1
            for _ in range(5):
                at = dict(zip(x, self.point(system)))
                self.assertAlmostEqual(float(on.subs(at)), 0.0, places=9, msg=system)
                for i in range(4):
                    for j in range(4):
                        self.assertAlmostEqual(float((pulled[i, j] - g[i, j]).subs(at)), 0.0, places=8,
                                               msg=f"{system} {i}{j}")

    def test_the_identification_each_chart_states_is_the_antipodal_map(self):
        # The last domain of the global, conformal and Kruskal charts glues the edge of the chart to
        # itself: the map it states, taken over the whole chart, is X -> -X.
        for system, image in self.pc.EDS_ANTIPODE.items():
            entry, x, g, X = self.chart(system)
            self.assertIn("\\sim", entry["domains"][-1])
            self.assertIn("antipodal identification", entry["domains"][-1])
            there = dict(zip(x, image(*x)))
            for _ in range(5):
                at = dict(zip(x, self.point(system)))
                for f in X:
                    self.assertAlmostEqual(float(f.subs(at)), -float(f.subs(there, simultaneous=True).subs(at)),
                                           places=9, msg=system)

    def test_the_edge_each_chart_glues_is_carried_onto_itself(self):
        # chi = pi/2 goes to chi = pi/2 with the time reversed, and U = V to itself with U -> -U.
        sp = self.sp
        for system in ("global", "conformal"):
            image = self.pc.EDS_ANTIPODE[system](sp.Symbol("a"), sp.pi / 2, sp.Symbol("b"), sp.Symbol("c"))
            self.assertEqual(image[0], -sp.Symbol("a"))
            self.assertEqual(image[1], sp.pi / 2)
        U = sp.Symbol("U")
        image = self.pc.EDS_ANTIPODE["kruskal"](U, U, sp.Symbol("b"), sp.Symbol("c"))
        self.assertEqual(image[:2], (-U, -U))

    def test_the_static_patch_and_the_planar_chart_hold_no_antipodal_pair(self):
        # X -> -X carries X_4 > |X_0| and X_0 + X_4 > 0 off themselves.
        for system, inside in (("static", lambda X: X[4] - abs(X[0])), ("planar", lambda X: X[0] + X[4])):
            entry, x, g, X = self.chart(system)
            for _ in range(20):
                at = dict(zip(x, self.point(system)))
                here = [float(f.subs(at)) for f in X]
                self.assertGreater(inside(here), 0, system)
                self.assertLess(inside([-v for v in here]), 0, system)

    def test_every_chart_has_the_curvature_of_de_sitter_space(self):
        for entry in published()["coordinates"]:
            self.assertEqual(entry["ricci_scalar"], "R = \\dfrac{12}{\\ell^2}")
            self.assertEqual(entry["kretschmann"], "K = \\dfrac{24}{\\ell^4}")
            self.assertEqual(entry["weyl_tensor"]["variants"]["llll"]["nonzero"], [])
            mixed = entry["einstein_tensor"]["variants"]["ul"]["nonzero"]
            self.assertEqual({e["value"] for e in mixed}, {"-\\dfrac{3}{\\ell^2}"})
            self.assertEqual(len(mixed), 4)


class Drawings(unittest.TestCase):
    def test_the_spacetime_diagrams_of_the_global_and_conformal_charts_are_a_band_with_one_glued_edge(self):
        diagrams = json.loads((DATA / "diagrams" / "elliptic_de_sitter.json").read_text(encoding="utf-8"))
        for system in ("global", "conformal"):
            view, = diagrams["systems"][system]
            self.assertTrue(view["mirror"])
            self.assertAlmostEqual(view["box"][1], HALF)
            edge, = [m for m in view["markers"] if m["kind"] == "surface"]
            self.assertEqual(edge["lines"], [[[1.0, 0.0], [1.0, 1.0]]])
            self.assertIn("glued to itself", edge["legend"])

    def test_the_kruskal_diagram_is_hatched_off_the_spacetime_and_marks_both_horizons(self):
        view, = json.loads((DATA / "diagrams" / "elliptic_de_sitter.json").read_text(encoding="utf-8"))["systems"]["kruskal"]
        X0, X1, T0, T1 = view["box"]
        self.assertEqual(X0, 0)
        horizons, = [m for m in view["markers"] if m["kind"] == "event"]
        # U = 0 and V = 0 leave the middle of the left edge at 45 degrees.
        for line in horizons["lines"]:
            ends = sorted(line)
            self.assertEqual(ends[0], [0.0, 0.5])
            dX, dT = (ends[1][0] - ends[0][0]) * (X1 - X0), (ends[1][1] - ends[0][1]) * (T1 - T0)
            self.assertAlmostEqual(abs(dT), dX, places=6)
        # A cone stands only where -1 < UV < 1.
        for cone in view["cones"]:
            X, T = X0 + cone["at"][0] * (X1 - X0), T0 + cone["at"][1] * (T1 - T0)
            U, V = T - X, T + X
            self.assertLess(abs(U * V), 1)

    def test_the_conformal_diagram_is_half_of_de_sitters_square_glued_along_its_edge(self):
        conformal = json.loads((DATA / "conformal" / "elliptic_de_sitter.json").read_text(encoding="utf-8"))
        self.assertEqual([v["id"] for v in conformal["views"]], ["global", "conformal", "kruskal", "static", "planar"])
        for view in conformal["views"]:
            region, = [layer for layer in view["layers"] if layer["class"] == "region"]
            xs, ts = [p[0] for p in region["points"]], [p[1] for p in region["points"]]
            self.assertEqual((min(xs), max(xs)), (0.0, round(HALF, 4)))
            self.assertEqual((min(ts), max(ts)), (-round(HALF, 4), round(HALF, 4)))
            edge, = [layer for layer in view["layers"] if layer["class"] == "surface"]
            self.assertTrue(all(p[0] == round(HALF, 4) for p in edge["points"]))
            # The two marks on the edge are one event: (chi, eta) -> (pi - chi, -eta).
            E, = [layer["at"] for layer in view["layers"] if layer["class"] == "mark" and layer["at"][0] < 1.5]
            P = [layer["at"] for layer in view["layers"] if layer["class"] == "mark" and layer["at"][0] > 1.5]
            self.assertEqual(len(P), 2)
            self.assertAlmostEqual(P[0][0], math.pi - P[1][0], places=3)
            self.assertAlmostEqual(P[0][1], -P[1][1], places=6)
            # E lies beyond the observer's future horizon eta = pi/2 - chi.
            self.assertGreater(E[1], HALF - E[0])
            # The ray runs at 45 degrees from E to one mark and from the other to the observer, which it
            # reaches before eta = pi/2.
            out, back = [layer["points"] for layer in view["layers"] if layer["class"] == "null"]
            for (a, b) in (out, back):
                self.assertAlmostEqual(abs(b[1] - a[1]), abs(b[0] - a[0]), places=3)
            self.assertEqual(out[0], E)
            self.assertIn(out[1], P)
            self.assertIn(back[0], P)
            self.assertNotEqual(out[1], back[0])
            self.assertEqual(back[1][0], 0.0)
            self.assertLess(back[1][1], HALF)
            labels = [label["text"] for label in view["labels"]]
            self.assertEqual(labels.count("$P$"), 2)
            self.assertEqual(labels.count("$E$"), 1)

    def test_the_planar_chart_covers_the_whole_half_square(self):
        # The half above the past horizon directly, and the triangle below it through the identification.
        conformal = json.loads((DATA / "conformal" / "elliptic_de_sitter.json").read_text(encoding="utf-8"))
        view = next(v for v in conformal["views"] if v["id"] == "planar")

        def area(points):
            return abs(sum(a[0] * b[1] - b[0] * a[1] for a, b in zip(points, points[1:] + points[:1]))) / 2
        direct, = [layer["points"] for layer in view["layers"] if layer["class"] == "cover"]
        carried, = [layer["points"] for layer in view["layers"] if layer["class"] == "cover2"]
        self.assertAlmostEqual(area(direct) + area(carried), HALF * math.pi, places=3)
        self.assertAlmostEqual(area(carried), HALF ** 2 / 2, places=3)
        edge, = [layer["points"] for layer in view["layers"] if layer["class"] == "chartedge"]
        self.assertEqual(edge, [[0.0, -round(HALF, 4)], [round(HALF, 4), 0.0]])

    def test_the_embedding_is_a_hemisphere_with_a_closed_geodesic_of_length_pi(self):
        embedding = json.loads((DATA / "embedding" / "elliptic_de_sitter.json").read_text(encoding="utf-8"))
        view, = embedding["views"]
        piece, = view["surfaces"][0]["pieces"]
        points = piece["points"]
        self.assertEqual((piece["system"], piece["coordinate"]), ("global", "\\chi"))
        self.assertEqual(points[0][1:], [0.0, -1.0])
        self.assertAlmostEqual(points[-1][0], HALF)
        self.assertEqual(points[-1][1:], [1.0, 0.0])
        for chi, rho, z in points:
            self.assertAlmostEqual(rho, math.sin(chi), places=5)
            self.assertAlmostEqual(z, -math.cos(chi), places=5)
        # The marked geodesic is the two opposite meridians, from rim to rim through the pole. Its two
        # ends on the rim are one point, so it closes after the length of half a great circle.
        marks = view["figure"]["turn"]["marks"]
        self.assertEqual(sorted(m["phi"] for m in marks), [0.0, math.pi])
        meridian = sum(math.hypot(b[1] - a[1], b[2] - a[2]) for a, b in zip(points, points[1:]))
        self.assertAlmostEqual(2 * meridian, math.pi, places=3)
        self.assertEqual(piece["end"]["kind"], "edge")


if __name__ == "__main__":
    unittest.main()
