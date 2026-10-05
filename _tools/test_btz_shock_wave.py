"""Shenker and Stanford's shock wave in the BTZ black hole: .venv.noindex/bin/python -m unittest discover -s _tools

What its texts and drawings state, held on the published files and on the geometry of the black
hole's Kruskal plane. Every drawing is at l = R = 1 and alpha = 1. The first class holds the
numbers the History and the captions lean on to Shenker and Stanford's geodesic, their (17) and
(18): it crosses the shock halfway between the horizons and grows longer by 2 l ln(1 + alpha/2).
The second holds the moment the embedding diagram draws to its closed form, the third the
drawings to the plane's geometry, and the fourth the charts and the relations.
"""
import json
import math
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "MFS" / "assets" / "data"
ME = "btz_shock_wave"
CHARTS = ["kruskal", "discontinuous", "exterior"]
ALPHA = 1.0


def load(folder, metric_id=ME):
    return json.loads((DATA / folder / f"{metric_id}.json").read_text(encoding="utf-8"))


def geodesic(alpha, r=1e6, n=200001):
    """Shenker and Stanford's (17) and (18) at t_L = t_R = 0, l = R = 1: the length from the left
    boundary at the radius r to the right one through the point (0, v) of the shock, extremized
    over v, and the v of the extremum. d_1 + d_2 is concave in v, so the extremum is its greatest
    value."""
    s = math.sqrt(r * r - 1)
    v = np.linspace(-alpha - 0.999, 0.999, n)
    left, right = r + s * (v + alpha), r - s * v
    ok = (left >= 1) & (right >= 1)
    d = np.full_like(v, -np.inf)
    d[ok] = np.arccosh(left[ok]) + np.arccosh(right[ok])
    i = int(np.argmax(d))
    return float(d[i]), float(v[i])


def rho(w, alpha=ALPHA):
    """The radius of the circle on the moment u + v = -alpha/2 at w = |u|."""
    w = np.abs(np.asarray(w, dtype=float))
    return (1 - alpha * w / 2 + w * w) / (1 + alpha * w / 2 - w * w)


class Geodesic(unittest.TestCase):
    def test_the_bridge_grows_longer_by_twice_the_log_of_one_plus_half_the_shift(self):
        flat, _ = geodesic(0.0)
        for alpha in (0.5, 1.0, 2.0):
            d, _ = geodesic(alpha)
            self.assertAlmostEqual(d - flat, 2 * math.log(1 + alpha / 2), delta=1e-4)

    def test_the_geodesic_crosses_the_shock_halfway_between_the_horizons(self):
        for alpha in (0.5, 1.0, 2.0):
            _, v = geodesic(alpha)
            self.assertAlmostEqual(v, -alpha / 2, delta=1e-4)


class Moment(unittest.TestCase):
    def test_the_two_halves_meet_on_the_shock_at_an_angle(self):
        # d rho/ds at w -> 0 with ds = 2 dw/(1 + alpha w/2 - w^2): -alpha R/2l on either side.
        h = 1e-7
        self.assertAlmostEqual((rho(h) - 1) / (2 * h), -ALPHA / 2, places=5)

    def test_each_half_narrows_to_a_neck_inside_its_horizon(self):
        w = np.linspace(0, 0.5, 50001)
        i = int(np.argmin(rho(w)))
        self.assertAlmostEqual(w[i], ALPHA / 4, places=4)
        self.assertAlmostEqual(float(rho(w[i])), (16 - ALPHA ** 2) / (16 + ALPHA ** 2), places=9)
        self.assertAlmostEqual(float(rho(w[i])), 15 / 17, places=9)
        # It crosses the horizon at w = 1/2 with the horizon's radius R.
        self.assertAlmostEqual(float(rho(0.5)), 1.0, places=12)

    def test_the_moment_lies_level_where_the_circles_grow_as_fast_as_the_distance(self):
        level = (math.sqrt(33) - 3) / 4
        slope = (2 * level - ALPHA / 2) / (1 + ALPHA * level / 2 - level ** 2)
        self.assertAlmostEqual(slope, 1.0, places=12)

    def test_the_embedding_draws_that_moment(self):
        view, = load("embedding")["views"]
        surface, = view["surfaces"]
        pieces = {p["id"]: np.array(p["points"]) for p in surface["pieces"]}
        self.assertEqual(sorted(pieces), ["ahead", "ahead_minkowski", "behind", "behind_minkowski"])
        for name, P in pieces.items():
            np.testing.assert_allclose(P[:, 1], rho(P[:, 0]), atol=1e-6, err_msg=name)
        # The two halves are each other's mirror image through the circle of the shock.
        for one, other in (("behind", "ahead"), ("behind_minkowski", "ahead_minkowski")):
            A, B = pieces[one], pieces[other][::-1]
            np.testing.assert_allclose(A[:, 0], -B[:, 0], atol=1e-12)
            np.testing.assert_allclose(A[:, 1], B[:, 1], atol=1e-6)
            np.testing.assert_allclose(A[:, 2], -B[:, 2], atol=1e-6)
        self.assertAlmostEqual(float(pieces["behind"][:, 1].min()), 15 / 17, places=4)
        self.assertAlmostEqual(float(pieces["behind_minkowski"][-1, 1]), 3.0, places=6)


class Drawings(unittest.TestCase):
    def test_every_chart_has_a_spacetime_diagram_within_the_allowed_shapes(self):
        systems = load("diagrams")["systems"]
        self.assertEqual(sorted(systems), sorted(CHARTS))
        for name, views in systems.items():
            self.assertEqual(len(views), 1, name)
            X0, X1, Y0, Y1 = views[0]["box"]
            self.assertTrue(0.5 <= (X1 - X0) / (Y1 - Y0) <= 2.0, name)

    def null_points(self, view, unit):
        """(u, v) of points of the unit square of a view drawn with (v - u)/2 across and (u + v)/2 up."""
        X0, X1, Y0, Y1 = view["box"]
        X, Y = X0 + unit[:, 0] * (X1 - X0), Y0 + unit[:, 1] * (Y1 - Y0)
        return Y - X, X + Y

    def test_the_kruskal_view_marks_two_horizons_that_miss_by_alpha(self):
        view, = load("diagrams")["systems"]["kruskal"]
        horizons = {m["legend"]: np.array(m["lines"][0]) for m in view["markers"] if m["kind"] == "event"}
        self.assertEqual(len(horizons), 2)
        for legend, line in horizons.items():
            u, v = self.null_points(view, line)
            if "right" in legend:
                np.testing.assert_allclose(v, 0.0, atol=1e-3)
                self.assertTrue(np.all(u <= 1e-3))
            else:
                np.testing.assert_allclose(v, -ALPHA, atol=1e-3)
                self.assertTrue(np.all(u >= -1e-3))
        shock, = [m for m in view["markers"] if m["kind"] == "shell"]
        u, _ = self.null_points(view, np.array(shock["lines"][0]))
        np.testing.assert_allclose(u, 0.0, atol=1e-3)

    def test_each_ray_moving_left_comes_out_of_the_shock_later_by_alpha(self):
        view, = load("diagrams")["systems"]["discontinuous"]
        crossed = 0
        for ray in view["rays"]["P"]:
            u, v = self.null_points(view, np.array(ray))
            if u[0] < -0.3 and u[-1] > 0.3 and abs(v[0]) < 1:
                crossed += 1
                self.assertAlmostEqual(v[-1] - v[0], ALPHA, delta=0.05)
        self.assertGreater(crossed, 2)

    def test_the_moment_is_marked_where_the_embedding_cuts_it(self):
        kruskal, = load("diagrams")["systems"]["kruskal"]
        for line in kruskal["slices"][0]["lines"]:
            u, v = self.null_points(kruskal, np.array(line))
            np.testing.assert_allclose(u + v, -ALPHA / 2, atol=1e-3)
        disc, = load("diagrams")["systems"]["discontinuous"]
        for line in disc["slices"][0]["lines"]:
            U, V = self.null_points(disc, np.array(line))
            behind = float(np.mean(U)) > 0
            np.testing.assert_allclose(U + V - ALPHA * behind, -ALPHA / 2, atol=1e-3)

    def test_the_conformal_diagram_is_symmetric_and_puts_the_moment_on_t_zero(self):
        for view in load("conformal")["views"]:
            for line in view["slices"][0]["lines"]:
                np.testing.assert_allclose(np.array(line)[:, 1], 0.0, atol=1e-4)
            lines = {}
            for layer in view["layers"]:
                if layer["kind"] in ("line", "zig"):
                    lines.setdefault(layer["class"], []).append(np.array(layer["points"]))
            # The past horizon of the right outside and the future horizon of the left are the lines
            # X + T = 2q with q = +-arctan(alpha/2).
            sums = sorted(round(float(np.mean(h.sum(axis=1))), 3) for h in lines["horizon"])
            self.assertEqual(sums, [round(-2 * math.atan(ALPHA / 2), 3), round(2 * math.atan(ALPHA / 2), 3)])
            for h in lines["horizon"]:
                self.assertLess(float(np.ptp(h.sum(axis=1))), 1e-3)
            # Each edge of one half is an edge of the other turned through half a turn about the centre.
            for cls in ("boundary", "singular"):
                first, second = lines[cls]
                image = -first
                gaps = [float(np.min(np.hypot(*(second - p).T))) for p in image]
                self.assertLess(max(gaps), 5e-3, cls)


class Charts(unittest.TestCase):
    def setUp(self):
        self.metric = load("metrics")
        self.charts = {c["id"]: c for c in self.metric["coordinates"]}

    def test_every_chart_states_the_curvature_of_the_black_hole(self):
        self.assertEqual(list(self.charts), CHARTS)
        for name, chart in self.charts.items():
            self.assertEqual(chart["ricci_scalar"], "R = -\\dfrac{6}{\\ell^2}", name)
            self.assertEqual(chart["kretschmann"], "K = \\dfrac{12}{\\ell^4}", name)
            for variant in chart["weyl_tensor"]["variants"].values():
                self.assertEqual(variant["nonzero"], [], name)

    def test_the_shock_carries_shenker_and_stanfords_shell(self):
        """R_uu + (2/l^2) g_uu = 2 alpha delta(u), their T_uu = alpha delta(u)/(4 pi G)."""
        ricci = {tuple(e["indices"]): e["value"] for e in self.charts["kruskal"]["ricci_tensor"]["variants"]["ll"]["nonzero"]}
        self.assertEqual(ricci[("u", "u")], "2\\alpha\\,\\delta(u)")
        metric = {tuple(e["indices"]): e["value"] for e in self.charts["kruskal"]["metric_components"]}
        self.assertNotIn(("u", "u"), metric)
        ricci = {tuple(e["indices"]): e["value"]
                 for e in self.charts["discontinuous"]["ricci_tensor"]["variants"]["ll"]["nonzero"]}
        metric = {tuple(e["indices"]): e["value"] for e in self.charts["discontinuous"]["metric_components"]}
        # -6 alpha delta(U) + (2/l^2)(4 alpha l^2 delta(U)) = 2 alpha delta(U).
        self.assertEqual(ricci[("U", "U")], "-6\\alpha\\,\\delta(U)")
        self.assertEqual(metric[("U", "U")], "4\\alpha\\,\\ell^2\\,\\delta(U)")

    def test_the_spacetime_is_tagged_with_no_symmetry_in_time(self):
        """The exterior chart is static, and a shift of its Killing time rescales alpha, so the
        chart that covers one outside alone does not speak for the spacetime."""
        tags = set(self.metric["tags"])
        self.assertFalse(tags & {"static", "stationary", "Einstein space", "vacuum"})
        self.assertTrue({"conformally flat", "three-dimensional", "shock wave", "BTZ black hole"} <= tags)

    def test_each_relation_is_answered(self):
        for relation in self.metric["related"]:
            other = load("metrics", relation["id"])
            back = [r for r in other["related"] if r["id"] == ME]
            self.assertEqual(len(back), 1, relation["id"])
        kinds = {r["id"]: r["kind"] for r in self.metric["related"]}
        self.assertEqual(kinds["btz"], "piece")


if __name__ == "__main__":
    unittest.main()
