"""A shell of light falling in anti-de Sitter space: .venv.noindex/bin/python -m unittest discover -s _tools

What its History and its drawings state, held on the published files and on closed forms. Every
drawing is at l = 1 with the final horizon at r = l. The first class holds the numbers the History
quotes to the formulas of their sources, the next four hold the spacetime diagrams, the conformal
diagrams, the embedding diagram and the curvature maps to the geometry of the shell.
"""
import json
import math
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "MFS" / "assets" / "data"
ME = "ads_vaidya"
ROUND = ("ingoing_three", "ingoing_four", "ingoing_five")
POINCARE = ("poincare_three", "poincare_four")
CHARTS = ROUND + POINCARE + ("flat_front",)


def load(folder):
    return json.loads((DATA / folder / f"{ME}.json").read_text(encoding="utf-8"))


def chart_points(view, marker_kind):
    """The points of one kind of marker of a spacetime diagram, in the box's own units."""
    X0, X1, Y0, Y1 = view["box"]
    marker = next(m for m in view["markers"] if m["kind"] == marker_kind)
    return [(X0 + u[0] * (X1 - X0), Y0 + u[1] * (Y1 - Y0)) for line in marker["lines"] for u in line]


class History(unittest.TestCase):
    def test_the_entanglement_tsunami_runs_at_liu_and_suhs_speeds(self):
        # v_E = (eta - 1)^((eta - 1)/2)/eta^(eta/2) with eta = 2(d - 1)/d, their (4): 1, 0.687 and 0.620.
        for d, quoted in ((2, 1.0), (3, 0.687), (4, 0.620)):
            eta = 2 * (d - 1) / d
            speed = (eta - 1) ** ((eta - 1) / 2) / eta ** (eta / 2) if eta > 1 else 1.0
            self.assertAlmostEqual(speed, quoted, delta=5e-4)

    def test_half_of_hbar_over_kT_is_some_three_tenths_of_a_fermi_over_c_in_a_collider(self):
        # hbar c = 197.327 MeV fm, and Balasubramanian and coauthors take T from 300 to 400 MeV.
        times = [0.5 * 197.327 / T for T in (300, 400)]
        self.assertTrue(all(abs(t - 0.3) < 0.06 for t in times), times)

    def test_light_of_mass_lambda_v_leaves_a_ray_from_the_singularity_for_lambda_up_to_an_eighth(self):
        # Near the centre the r^2 of f is negligible, and a ray r = v/X of dr/dv = (1 + r^2 - lambda v/r)/2,
        # at l = 1 with 2Gm/c^3 = lambda v, needs lambda X^2 - X + 2 = 0, which has a root for lambda <= 1/8.
        for lam in (0.05, 0.1, 0.125):
            X = (1 - math.sqrt(1 - 8 * lam)) / (2 * lam)
            v = np.geomspace(1e-9, 2.0, 200001)
            r = np.empty_like(v)
            r[0] = v[0] / X
            for k in range(len(v) - 1):
                h = v[k + 1] - v[k]
                slope = lambda vv, rr: (1 + rr * rr - lam * vv / rr) / 2  # noqa: E731
                k1 = slope(v[k], r[k])
                k2 = slope(v[k] + h, r[k] + h * k1)
                r[k + 1] = r[k] + h * (k1 + k2) / 2
            # The ray stays off the centre and is still moving outward when it leaves the drawing.
            self.assertTrue(np.all(r > 0), lam)
            self.assertGreater(r[-1], 0.5, lam)
            self.assertGreater(r[-1] - r[-2], 0, lam)
        for lam in (0.13, 0.2):
            X = np.linspace(0, 100, 100001)
            self.assertTrue(np.all(lam * X * X - X + 2 > 0), lam)


class SpacetimeDiagrams(unittest.TestCase):
    def setUp(self):
        self.views = {system: views[0] for system, views in load("diagrams")["systems"].items()}

    def test_every_chart_has_its_diagram(self):
        self.assertEqual(sorted(self.views), sorted(CHARTS))

    def test_the_event_horizon_leaves_the_centre_at_minus_half_pi_and_stands_on_l_after_the_shell(self):
        for system in ROUND:
            points = chart_points(self.views[system], "event")
            first = min(points, key=lambda p: p[1])
            self.assertLess(first[0], 0.02, system)
            self.assertAlmostEqual(first[1], -math.pi / 2, delta=0.01, msg=system)
            late = [r for r, v in points if v > 0.3]
            self.assertTrue(late and max(abs(r - 1) for r in late) < 2e-3, system)

    def test_behind_a_flat_front_the_event_horizon_is_two_l_squared_over_two_l_minus_cv(self):
        for r, v in chart_points(self.views["flat_front"], "event"):
            if v < -0.1:
                self.assertAlmostEqual(r, 2 / (2 - v), delta=3e-3)

    def test_in_the_poincare_charts_the_event_horizon_is_l_minus_cv_over_two_before_the_shell(self):
        for system in POINCARE:
            for z, w in chart_points(self.views[system], "event"):
                v = w - z
                if v < -0.1:
                    self.assertAlmostEqual(z, 1 - v / 2, delta=3e-3, msg=system)
                elif v > 0.3:
                    self.assertAlmostEqual(z, 1.0, delta=3e-3, msg=system)

    def test_only_the_round_chart_of_three_dimensions_shows_the_moments(self):
        for system, view in self.views.items():
            self.assertEqual(len(view.get("slices", [])), 4 if system == "ingoing_three" else 0, system)


class ConformalDiagrams(unittest.TestCase):
    def setUp(self):
        self.views = {view["id"]: view for view in load("conformal")["views"]}

    def lines(self, view, cls):
        return [layer["points"] for layer in view["layers"] if layer["class"] == cls]

    def test_one_view_for_each_chart(self):
        self.assertEqual(sorted(self.views), sorted(CHARTS))
        for vid, view in self.views.items():
            self.assertEqual(view["system"], vid)

    def test_the_event_horizon_is_the_line_T_minus_X_equals_minus_half_pi_in_every_round_chart(self):
        for system in ROUND:
            for X, T in (p for line in self.lines(self.views[system], "event") for p in line):
                self.assertAlmostEqual(T - X, -math.pi / 2, delta=2e-4, msg=system)

    def test_the_boundary_is_the_vertical_line_X_equals_half_pi(self):
        for vid, view in self.views.items():
            for X, _ in (p for line in self.lines(view, "boundary") for p in line):
                self.assertAlmostEqual(X, math.pi / 2, delta=2e-4, msg=vid)

    def test_after_the_shell_r_equals_zero_is_the_straight_line_T_equals_zero_in_three_dimensions(self):
        points = [p for line in self.lines(self.views["ingoing_three"], "singular") for p in line]
        self.assertTrue(points)
        self.assertTrue(all(abs(T) < 2e-4 for _, T in points))

    def test_the_shell_is_a_light_ray(self):
        for vid, view in self.views.items():
            (a, b), = self.lines(view, "surface")
            self.assertAlmostEqual(abs(b[0] - a[0]), abs(b[1] - a[1]), delta=2e-4, msg=vid)

    def test_only_the_round_chart_of_three_dimensions_shows_the_moments(self):
        for vid, view in self.views.items():
            self.assertEqual(len(view.get("slices", [])), 4 if vid == "ingoing_three" else 0, vid)


class EmbeddingDiagram(unittest.TestCase):
    def setUp(self):
        self.view, = load("embedding")["views"]
        self.frames = self.view["movie"]["frames"]

    def test_the_movie_steps_through_thirty_six_moments_a_twentieth_apart(self):
        values = [f["value"] for f in self.frames]
        self.assertEqual(len(values), 36)
        self.assertTrue(np.allclose(np.diff(values), 0.05, atol=1e-9))
        self.assertAlmostEqual(values[0], -1.25)
        self.assertEqual(len({f["label"] for f in self.frames}), 36)

    def test_inside_the_shell_each_moment_is_the_hyperboloid(self):
        for frame in self.frames:
            for piece in frame["pieces"]:
                if piece["id"] == "inside":
                    r, rho, z = np.array(piece["points"]).T
                    self.assertTrue(np.allclose(rho, r, atol=1e-6))
                    self.assertTrue(np.allclose(z - z[0], np.sqrt(1 + r * r) - 1, atol=2e-5), frame["label"])

    def test_the_rim_stands_at_z_equals_zero_and_the_horizon_at_l_over_R_or_l(self):
        for frame in self.frames:
            T = frame["value"]
            self.assertAlmostEqual(frame["pieces"][-1]["points"][-1][2], 0.0, places=6)
            horizon = [ring["x"] for ring in frame["rings"] if ring["class"] == "horizon"]
            R = math.tan(-T) if T < 0 else 0.0
            self.assertEqual(len(horizon), 1, frame["label"])
            self.assertAlmostEqual(horizon[0], 1 / R if R > 1 else 1.0, places=9, msg=frame["label"])

    def test_after_the_shell_reaches_the_centre_the_moment_ends_in_the_tip_of_a_cone(self):
        for frame in self.frames:
            if frame["value"] >= 0:
                first = frame["pieces"][0]
                (r0, _, z0), (r1, _, z1) = first["points"][:2]
                self.assertEqual(first["id"], "between")
                self.assertAlmostEqual((z1 - z0) / (r1 - r0), math.sqrt(2), delta=2e-3)


class CurvatureMaps(unittest.TestCase):
    def setUp(self):
        self.views = {view["id"]: view for view in load("curvature")["views"]}

    def test_a_map_for_each_chart_of_four_and_five_dimensions(self):
        self.assertEqual(sorted(self.views), sorted(("ingoing_four", "flat_front", "poincare_four", "ingoing_five")))

    def test_K_is_anti_de_sitters_before_the_shell_and_the_holes_after_it(self):
        after = {"ingoing_four": (24, lambda x: 48 / x ** 6), "flat_front": (24, lambda x: 12 / x ** 6),
                 "poincare_four": (24, lambda x: 12 * x ** 6), "ingoing_five": (40, lambda x: 288 / x ** 8)}
        for vid, view in self.views.items():
            grid = view["map"]
            cells = np.array(grid["cells"], dtype=float)
            nx, ny = grid["nx"], grid["ny"]
            xs = view["x"]["from"] + (np.arange(nx) + 0.5) * (view["x"]["to"] - view["x"]["from"]) / nx
            ys = view["y"]["from"] + (np.arange(ny) + 0.5) * (view["y"]["to"] - view["y"]["from"]) / ny
            ads, hole = after[vid]
            self.assertTrue(np.allclose(cells[ys < -0.2], ads, rtol=1e-5), vid)
            X = np.broadcast_to(xs, (int(np.sum(ys > 0.2)), nx))
            self.assertTrue(np.allclose(cells[ys > 0.2], ads + hole(X), rtol=1e-5), vid)
