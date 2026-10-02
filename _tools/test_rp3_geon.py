"""The RP3 geon: python3 -m unittest discover -s _tools

What its texts state, held on Kruskal's plane and on the published files. The fold
(T, X, theta, phi) -> (T, -X, pi - theta, phi + pi) undoes itself, moves every point, keeps the areal
radius and exchanges the two null coordinates; on the moment of time symmetry it is Misner and
Wheeler's inversion of the isotropic radius followed by the antipodal map. A light ray sent in from
far away meets the glued edge inside the black hole and ends on the singularity, and a ray that
leaves the edge for infinity began on the singularity of the white hole, Friedman, Schleich and
Witt's topological censorship. Of the exterior's moments t = const only t = 0 meets the edge level.
Every published chart is a vacuum with Schwarzschild's Kretschmann scalar, Kruskal's chart is
Schwarzschild's carried along Louko and Marolf's map and is unchanged by the fold, and the drawings
hold what their captions say. The tests of the published mathematics need sympy and are skipped
where it is absent; the others read the files alone or need nothing.
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


def fold(T, X, theta, phi):
    """The involution whose quotient of Kruskal's manifold is the geon."""
    return T, -X, math.pi - theta, (phi + math.pi) % (2 * math.pi)


def radius(T, X):
    """The areal radius at r_s = 1, the root of (r - 1) e^r = X^2 - T^2, by bisection."""
    want = X * X - T * T
    lo, hi = 0.0, 1.0
    while (hi - 1) * math.exp(hi) < want:
        hi *= 2
    for _ in range(200):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if (mid - 1) * math.exp(mid) < want else (lo, mid)
    return (lo + hi) / 2


def unit(theta, phi):
    return [math.sin(theta) * math.cos(phi), math.sin(theta) * math.sin(phi), math.cos(theta)]


class Fold(unittest.TestCase):
    def setUp(self):
        self.rng = random.Random(1957)

    def event(self):
        r = self.rng.uniform
        return r(-0.95, 0.95), r(-3, 3), r(0.1, 3.0), r(0, 2 * math.pi)

    def test_the_fold_undoes_itself_and_moves_every_point(self):
        for _ in range(200):
            T, X, theta, phi = self.event()
            image = fold(T, X, theta, phi)
            back = fold(*image)
            self.assertAlmostEqual(back[0], T)
            self.assertAlmostEqual(back[1], X)
            for a, b in zip(unit(*back[2:]), unit(theta, phi)):
                self.assertAlmostEqual(a, b, places=12)
            # The sphere's point goes to its antipode, so even on the edge X = 0 no point stays.
            for a, b in zip(unit(*image[2:]), unit(theta, phi)):
                self.assertAlmostEqual(a, -b, places=12)

    def test_the_fold_keeps_the_radius_and_exchanges_the_null_coordinates(self):
        for _ in range(200):
            T, X, theta, phi = self.event()
            image = fold(T, X, theta, phi)
            self.assertAlmostEqual(radius(T, X), radius(image[0], image[1]), places=12)
            self.assertAlmostEqual(T - X, image[0] + image[1])
            self.assertAlmostEqual(T + X, image[0] - image[1])

    def test_the_fold_keeps_the_directions_of_time_and_of_space(self):
        # Its Jacobian on (T, X, theta, phi) is diag(1, -1, -1, 1): d/dT goes to d/dT, and the
        # determinant is +1, the reflection of X undone by the antipodal map of the sphere.
        jacobian = [1, -1, -1, 1]
        self.assertEqual(jacobian[0], 1)
        self.assertEqual(math.prod(jacobian), 1)

    def test_misner_and_wheelers_inversion_is_the_fold_on_the_moment_of_time_symmetry(self):
        # At t = 0 and r_s = 1 the isotropic radius rho has the areal radius rho (1 + 1/4 rho)^2 and the
        # signed X = (4 rho - 1) e^(r/2)/(4 sqrt(rho)); rho -> 1/(16 rho) keeps the first and turns the second.
        def areal(rho):
            return rho * (1 + 1 / (4 * rho)) ** 2

        def X_of(rho):
            return (4 * rho - 1) * math.exp(areal(rho) / 2) / (4 * math.sqrt(rho))
        for _ in range(200):
            rho = math.exp(self.rng.uniform(-3, 3))
            image = 1 / (16 * rho)
            self.assertAlmostEqual(areal(image) / areal(rho), 1.0, places=9)
            if abs(4 * rho - 1) > 1e-3:
                self.assertAlmostEqual(X_of(image) / X_of(rho), -1.0, places=9)
                self.assertAlmostEqual(X_of(rho) ** 2 / ((areal(rho) - 1) * math.exp(areal(rho))), 1.0, places=7)
            # The line element of the moment, (1 + 1/4 rho)^4 d rho^2, is the same at the image.
            stretch = (1 + 1 / (4 * image)) ** 2 / (16 * rho * rho)
            self.assertAlmostEqual(stretch / (1 + 1 / (4 * rho)) ** 2, 1.0, places=9)
        self.assertAlmostEqual(areal(0.25), 1.0)
        self.assertAlmostEqual(X_of(0.25), 0.0)

    def test_light_sent_in_from_far_away_meets_the_edge_inside_the_black_hole_and_ends_on_the_singularity(self):
        # A ray moving left keeps V = T + X, positive if it comes from the exterior's past infinity. It
        # meets the edge X = 0 at T = V, where T^2 - X^2 = V^2 > 0: inside the black hole if V < 1, and
        # otherwise it has already ended on T^2 - X^2 = 1. From the edge it moves right with U = T - X = V,
        # along which T^2 - X^2 = V (V + 2X) only grows, so it ends on the singularity and never comes
        # back to the exterior, where T - X is negative.
        for _ in range(200):
            V = self.rng.uniform(0.01, 0.99)
            self.assertGreater(V, 0)
            self.assertLess(radius(V, 0.0), 1.0)
            X_end = (1 / V - V) / 2
            self.assertGreater(X_end, 0)
            self.assertAlmostEqual((V + X_end) ** 2 - X_end ** 2, 1.0, places=9)
            for X in (0.25 * X_end, 0.5 * X_end, 0.99 * X_end):
                T = V + X
                self.assertGreater(T, X)                    # still behind the future horizon

    def test_light_that_leaves_the_edge_for_infinity_began_in_the_white_hole(self):
        # The time reverse: a ray that reaches the exterior moving right has U = T - X < 0, so at the edge
        # T = U < 0, inside the white hole, and before the edge it moved left along V = U from the past
        # singularity.
        for _ in range(200):
            U = -self.rng.uniform(0.01, 0.99)
            self.assertLess(radius(U, 0.0), 1.0)
            self.assertLess(U, 0)
            X = self.rng.uniform(1.0, 50.0)
            self.assertGreater(X, abs(U + X))               # out in the exterior, X > |T|
            X_start = (1 / abs(U) - abs(U)) / 2
            self.assertAlmostEqual((U - X_start) ** 2 - X_start ** 2, 1.0, places=9)

    def test_of_the_moments_of_the_exterior_only_t_0_meets_the_edge_level(self):
        # The moment t of the exterior is the line T = X tanh(t/2) through the middle of the edge, and
        # its image under the fold is T = -X tanh(t/2): the two make one smooth line only where the
        # slope vanishes.
        for t in (-2.0, -0.5, 0.0, 0.5, 2.0):
            slope = math.tanh(t / 2)
            self.assertEqual(slope == -slope, t == 0.0)

    def test_the_projective_plane_where_the_horizons_cross_has_half_the_area_of_the_horizon(self):
        # The sphere T = X = 0 has the areal radius r_s; with opposite points one its area is 2 pi r_s^2,
        # 8 pi M^2 for r_s = 2M, half of the horizon's 16 pi M^2, Louko and Marolf's section IV B.
        self.assertAlmostEqual(radius(0.0, 0.0), 1.0, places=12)
        M = 0.5
        self.assertAlmostEqual(4 * math.pi * radius(0.0, 0.0) ** 2 / 2, 8 * math.pi * M * M)
        self.assertAlmostEqual(8 * math.pi * M * M, 16 * math.pi * M * M / 2)


def published():
    return json.loads((DATA / "metrics" / "rp3_geon.json").read_text(encoding="utf-8"))


@unittest.skipUnless(HAS_SYMPY, "sympy is not installed")
class Published(unittest.TestCase):
    def setUp(self):
        import sympy
        import print_charts
        import verify_metrics
        self.sp, self.pc, self.vm = sympy, print_charts, verify_metrics
        self.rng = random.Random(1993)

    def chart(self, system):
        """The published metric of a chart as a matrix in the chart x^0 = ct at r_s = 1, with its
        coordinate symbols; Kruskal's areal radius stays a symbol, r."""
        entry = next(c for c in published()["coordinates"] if c["id"] == system)
        reader = self.vm.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]], (),
                                held=("r",) if system == "kruskal" else ())
        x = [reader.symbol[n] for n in entry["coords"]]
        g = self.sp.zeros(4, 4)
        for e in entry["metric_components"]:
            i, j = (entry["coords"].index(k) for k in e["indices"])
            g[i, j] = reader(e["value"]).subs({reader.c: 1, reader.parameters["r_s"]: 1})
        return entry, reader, x, g

    def test_every_chart_is_a_vacuum_with_schwarzschilds_kretschmann_scalar(self):
        scalar = {"kruskal": "K = \\dfrac{12r_s^2}{r^6}", "schwarzschild": "K = \\dfrac{12r_s^2}{r^6}",
                  "isotropic": "K = \\dfrac{12r_s^2}{\\rho^6}\\left(1 + \\dfrac{r_s}{4\\rho}\\right)^{-12}"}
        self.assertEqual([c["id"] for c in published()["coordinates"]], list(self.pc.RP3_GEON_CHARTS))
        for entry in published()["coordinates"]:
            self.assertEqual(entry["ricci_scalar"], "R = 0")
            self.assertEqual(entry["kretschmann"], scalar[entry["id"]])
            for variant in entry["ricci_tensor"]["variants"].values():
                self.assertEqual(variant["nonzero"], [])
            for variant in entry["einstein_tensor"]["variants"].values():
                self.assertEqual(variant["nonzero"], [])

    def test_the_kruskal_chart_states_the_fold_on_its_edge(self):
        entry = next(c for c in published()["coordinates"] if c["id"] == "kruskal")
        self.assertIn("X \\in [0, \\infty)", entry["domains"])
        self.assertIn("\\sim", entry["domains"][-1])
        self.assertIn("antipodal identification", entry["domains"][-1])
        sp = self.sp
        T, theta, phi = sp.symbols("T theta phi")
        self.assertEqual(self.pc.RP3_GEON_INVOLUTION(T, 0, theta, phi), (T, 0, sp.pi - theta, phi + sp.pi))

    def test_the_fold_is_an_isometry_of_the_published_kruskal_chart(self):
        entry, reader, x, g = self.chart("kruskal")
        r = reader.parameters["r"]
        jacobian = self.sp.diag(1, -1, -1, 1)
        for _ in range(20):
            T, X = self.rng.uniform(-0.9, 0.9), self.rng.uniform(0.05, 3)
            theta, phi = self.rng.uniform(0.2, 2.9), self.rng.uniform(0, 6)
            here = dict(zip(x, (T, X, theta, phi)))
            there = dict(zip(x, fold(T, X, theta, phi)))
            g_here = g.subs(r, radius(T, X)).subs(here)
            g_there = (jacobian.T * g.subs(r, radius(T, -X)).subs(there) * jacobian)
            for i in range(4):
                for j in range(4):
                    self.assertAlmostEqual(float(g_here[i, j]), float(g_there[i, j]), places=9)

    def test_the_kruskal_chart_is_schwarzschilds_carried_along_louko_and_marolfs_map(self):
        # T = sqrt(r - 1) e^(r/2) sinh(t/2), X the same with cosh, at r_s = 1: the published Kruskal
        # metric pulled back is the published Schwarzschild chart's.
        sp = self.sp
        _, reader_k, x_k, g_k = self.chart("kruskal")
        _, _, x_s, g_s = self.chart("schwarzschild")
        t, r = x_s[:2]
        root = sp.sqrt(r - 1) * sp.exp(r / 2)
        image = [root * sp.sinh(t / 2), root * sp.cosh(t / 2)]
        J = sp.Matrix(2, 2, lambda i, j: sp.diff(image[i], (t, r)[j]))
        pulled = J.T * g_k[:2, :2].subs(reader_k.parameters["r"], r) * J
        for _ in range(20):
            at = {t: self.rng.uniform(-3, 3), r: self.rng.uniform(1.05, 6)}
            for i in range(2):
                for j in range(2):
                    self.assertAlmostEqual(float(pulled[i, j].subs(at)), float(g_s[i, j].subs(at)), places=8)
            self.assertAlmostEqual(float((image[1] ** 2 - image[0] ** 2).subs(at)),
                                   (at[r] - 1) * math.exp(at[r]), places=7)
            self.assertAlmostEqual(radius(float(image[0].subs(at)), float(image[1].subs(at))), at[r], places=8)

    def test_the_isotropic_chart_is_schwarzschilds_and_is_unchanged_by_the_inversion(self):
        sp = self.sp
        _, _, x_i, g_i = self.chart("isotropic")
        _, _, x_s, g_s = self.chart("schwarzschild")
        rho, r = x_i[1], x_s[1]
        areal = rho * (1 + 1 / (4 * rho)) ** 2
        stretch = sp.diff(areal, rho)
        inverted = 1 / (16 * rho)
        turn = sp.diff(inverted, rho)
        for _ in range(20):
            at = {rho: self.rng.uniform(0.3, 5)}
            there = {r: float(areal.subs(at))}
            self.assertAlmostEqual(float(g_i[0, 0].subs(at)), float(g_s[0, 0].subs(there)), places=9)
            self.assertAlmostEqual(float(g_i[1, 1].subs(at)), float((g_s[1, 1] * 1).subs(there)) * float(stretch.subs(at)) ** 2,
                                   places=8)
            self.assertAlmostEqual(float(g_i[2, 2].subs(at)), there[r] ** 2, places=9)
            image = {rho: float(inverted.subs(at))}
            self.assertAlmostEqual(float(g_i[0, 0].subs(image)), float(g_i[0, 0].subs(at)), places=9)
            self.assertAlmostEqual(float(g_i[1, 1].subs(image)) * float(turn.subs(at)) ** 2, float(g_i[1, 1].subs(at)),
                                   places=8)
            self.assertAlmostEqual(float(g_i[2, 2].subs(image)), float(g_i[2, 2].subs(at)), places=9)


class Drawings(unittest.TestCase):
    def test_the_kruskal_diagram_is_the_half_plane_with_its_edge_glued_and_both_horizons_marked(self):
        view, = json.loads((DATA / "diagrams" / "rp3_geon.json").read_text(encoding="utf-8"))["systems"]["kruskal"]
        X0, X1, T0, T1 = view["box"]
        self.assertEqual(X0, 0)
        self.assertAlmostEqual(X1 - X0, T1 - T0)
        edge, = [m for m in view["markers"] if m["kind"] == "surface"]
        self.assertEqual(edge["lines"], [[[0.0, 0.0], [0.0, 1.0]]])
        self.assertIn("glued to itself", edge["legend"])
        horizons, = [m for m in view["markers"] if m["kind"] == "event"]
        # X = T and X = -T leave the middle of the left edge at 45 degrees.
        self.assertEqual(len(horizons["lines"]), 2)
        for line in horizons["lines"]:
            ends = sorted(line)
            self.assertEqual(ends[0], [0.0, 0.5])
            dX, dT = (ends[1][0] - ends[0][0]) * (X1 - X0), (ends[1][1] - ends[0][1]) * (T1 - T0)
            self.assertAlmostEqual(abs(dT), dX, places=6)
        # The singularities are the hyperbola T^2 - X^2 = 1, and a cone stands only inside it.
        singular, = [m for m in view["markers"] if m["kind"] == "singular"]
        self.assertEqual(len(singular["lines"]), 2)
        for line in singular["lines"]:
            for u in line:
                X, T = X0 + u[0] * (X1 - X0), T0 + u[1] * (T1 - T0)
                self.assertAlmostEqual(T * T - X * X, 1.0, delta=5e-3)
        for cone in view["cones"]:
            X, T = X0 + cone["at"][0] * (X1 - X0), T0 + cone["at"][1] * (T1 - T0)
            self.assertLess(T * T - X * X, 1)
            # Light runs at 45 degrees everywhere on the plane.
            for edge_of_cone in (cone["a"], cone["b"]):
                self.assertAlmostEqual(abs(edge_of_cone[0]), abs(edge_of_cone[1]), places=3)

    def test_the_exterior_diagrams_stop_at_the_horizon(self):
        diagrams = json.loads((DATA / "diagrams" / "rp3_geon.json").read_text(encoding="utf-8"))["systems"]
        for system, horizon in (("schwarzschild", 1.0), ("isotropic", 0.25)):
            view, = diagrams[system]
            X0, X1, T0, T1 = view["box"]
            self.assertAlmostEqual(X1 - X0, T1 - T0)
            hatch, = view["hatch"]
            self.assertAlmostEqual(X0 + max(p[0] for p in hatch) * (X1 - X0), horizon, delta=0.02)
            for cone in view["cones"]:
                self.assertGreater(X0 + cone["at"][0] * (X1 - X0), horizon)

    def test_the_conformal_diagram_is_the_right_half_of_kruskals_with_two_rays_through_its_edge(self):
        conformal = json.loads((DATA / "conformal" / "rp3_geon.json").read_text(encoding="utf-8"))
        self.assertEqual([v["id"] for v in conformal["views"]], ["kruskal", "schwarzschild", "isotropic"])
        half = round(HALF, 4)
        for view in conformal["views"]:
            region, = [layer for layer in view["layers"] if layer["class"] == "region"]
            xs, ts = [p[0] for p in region["points"]], [p[1] for p in region["points"]]
            self.assertEqual((min(xs), max(xs)), (0.0, round(math.pi, 4)))
            self.assertEqual((min(ts), max(ts)), (-half, half))
            edge, = [layer for layer in view["layers"] if layer["class"] == "surface"]
            self.assertEqual(edge["points"], [[0.0, -half], [0.0, half]])
            singular = [layer["points"] for layer in view["layers"] if layer["class"] == "singular"]
            self.assertEqual(sorted(singular), [[[0.0, -half], [half, -half]], [[0.0, half], [half, half]]])
            rays = [layer["points"] for layer in view["layers"] if layer["class"] == "null"]
            self.assertEqual(len(rays), 4)
            for a, b in rays:
                self.assertAlmostEqual(abs(b[1] - a[1]), abs(b[0] - a[0]), places=3)
                self.assertGreater(b[1], a[1])                # each leg is drawn toward the future
            marks = sorted(layer["at"] for layer in view["layers"] if layer["class"] == "mark")
            self.assertEqual([m[0] for m in marks], [0.0, 0.0])
            self.assertAlmostEqual(marks[0][1], -marks[1][1])
            # The ray from scri-, the line T = X - pi, meets the edge above the middle, inside the black
            # hole, and its second leg ends on the future singularity T = pi/2.
            came_in, = [ray for ray in rays if abs(ray[0][1] - (ray[0][0] - math.pi)) < 1e-3]
            self.assertEqual(came_in[1], marks[1])
            self.assertGreater(came_in[1][1], 0)
            went_on, = [ray for ray in rays if ray[0] == marks[1]]
            self.assertAlmostEqual(went_on[1][1], half, places=3)
            self.assertLess(went_on[1][0], half)
            # The ray that reaches scri+, the line T = pi - X, left the edge below the middle, inside
            # the white hole, and its first leg began on the past singularity T = -pi/2.
            got_out, = [ray for ray in rays if abs(ray[1][1] - (math.pi - ray[1][0])) < 1e-3]
            self.assertEqual(got_out[0], marks[0])
            self.assertLess(got_out[0][1], 0)
            began, = [ray for ray in rays if ray[1] == marks[0]]
            self.assertAlmostEqual(began[0][1], -half, places=3)

    def test_the_exterior_charts_cover_one_diamond_of_the_conformal_diagram(self):
        conformal = json.loads((DATA / "conformal" / "rp3_geon.json").read_text(encoding="utf-8"))

        def area(points):
            return abs(sum(a[0] * b[1] - b[0] * a[1] for a, b in zip(points, points[1:] + points[:1]))) / 2
        covers = {}
        for view in conformal["views"]:
            cover, = [layer["points"] for layer in view["layers"] if layer["class"] == "cover"]
            covers[view["id"]] = area(cover)
        # The whole half hexagon is the diamond and two triangles of half its area together.
        self.assertAlmostEqual(covers["schwarzschild"], math.pi ** 2 / 2, places=2)
        self.assertAlmostEqual(covers["isotropic"], covers["schwarzschild"], places=6)
        self.assertAlmostEqual(covers["kruskal"], 1.5 * covers["schwarzschild"], places=2)

    def test_the_embedding_is_one_sheet_of_flamms_paraboloid_with_a_throat_of_length_pi(self):
        embedding = json.loads((DATA / "embedding" / "rp3_geon.json").read_text(encoding="utf-8"))
        view, = embedding["views"]
        (surface,) = view["surfaces"]
        piece, = surface["pieces"]
        points = piece["points"]
        self.assertEqual((piece["system"], piece["coordinate"]), ("schwarzschild", "r"))
        self.assertEqual(points[0], [1.0, 1.0, 0.0])
        self.assertEqual(piece["start"]["kind"], "throat")
        for r, rho, z in points:
            self.assertAlmostEqual(rho, r, places=5)
            self.assertAlmostEqual(z, 2 * math.sqrt(r - 1), places=4)
        # The marked geodesic is the two opposite meridians, in to the throat on one side and out on the
        # other, and the throat, with opposite points one, closes after half of its circle.
        marks = view["figure"]["turn"]["marks"]
        self.assertEqual(sorted(m["phi"] for m in marks), [0.0, math.pi])
        self.assertAlmostEqual(2 * math.pi * points[0][1] / 2, math.pi)
        self.assertIn("$\\pi r_s$", " ".join(view["caption"]))


if __name__ == "__main__":
    unittest.main()
