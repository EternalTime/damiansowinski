"""The generalised Vaidya spacetime: .venv.noindex/bin/python -m unittest discover -s _tools

What its History, its charts and its drawings state, held on the published files and on closed forms.
The first class holds the relations the History and the related texts quote, the second the published
charts to Wang and Wu's Einstein tensor and Husain's fluid, the next three the spacetime diagrams, the
embedding diagram and the curvature maps to the mass functions they declare. Every drawing is in units
of the length L.
"""
import json
import math
import sys
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "MFS" / "assets" / "data"
ME = "generalised_vaidya"
sys.path.insert(0, str(ROOT / "_tools" / "derivations"))
import generalised_vaidya_drawings as drawings  # noqa: E402


def load(folder, name=ME):
    return json.loads((DATA / folder / f"{name}.json").read_text(encoding="utf-8"))


def chart(system):
    return next(c for c in load("metrics")["coordinates"] if c["id"] == system)


def marker_points(view, kind):
    """The points of one kind of marker of a spacetime diagram, in the box's own units."""
    X0, X1, Y0, Y1 = view["box"]
    marker = next(m for m in view["markers"] if m["kind"] == kind)
    return [(X0 + u[0] * (X1 - X0), Y0 + u[1] * (Y1 - Y0)) for line in marker["lines"] for u in line]


def view_of(system, view_id):
    return next(v for v in load("diagrams")["systems"][system] if v["id"] == view_id)


class History(unittest.TestCase):
    def test_husains_static_term_is_kiselevs_at_w_equal_to_2k_minus_1_over_3(self):
        # Husain's lapse holds r^(-2k), Kiselev's (r_q/r)^(3w + 1): the powers agree at w = (2k - 1)/3,
        # Heydarzade and Darabi's map, and radiation's k = 1 is Kiselev's w = 1/3.
        for k in (1 / 3, 3 / 4, 1.0, 2.0):
            w = (2 * k - 1) / 3
            self.assertAlmostEqual(3 * w + 1, 2 * k, places=12)
        self.assertAlmostEqual((2 * 1 - 1) / 3, 1 / 3, places=12)

    def test_husains_k_equal_to_one_with_2g_equal_to_q_squared_is_reissner_nordstrom(self):
        # 1 - 2f/r + 2g/((2k - 1) r^(2k)) at k = 1, f = M and 2g = Q^2 is 1 - 2M/r + Q^2/r^2.
        M, Q = 0.7, 0.4
        for r in (0.5, 1.3, 4.0):
            g = Q * Q / 2
            self.assertAlmostEqual(1 - 2 * M / r + 2 * g / ((2 * 1 - 1) * r ** 2), 1 - 2 * M / r + Q * Q / r ** 2,
                                   places=14)

    def test_brown_and_husains_q_is_husains_g(self):
        # Brown and Husain's 1 - 2M/r + Q^(2k)/r^(2k) is Husain's lapse with Q^(2k) = 2g/(2k - 1).
        for k in (0.75, 1.5, 3.0):
            g, r = 0.3, 1.7
            Q2k = 2 * g / (2 * k - 1)
            self.assertAlmostEqual(Q2k / r ** (2 * k), 2 * g / ((2 * k - 1) * r ** (2 * k)), places=14)

    def test_the_drawn_collapse_settles_on_husains_cubic(self):
        # At late times 2m/r = 2A/r + 2B r^(-2/3) with A = L/2, B = L^(2/3)/8, and the apparent horizon
        # solves Husain's (r - 2A)^3 = (2B)^3 r.
        r, A, B = drawings.COLLAPSE_HORIZON, 0.5, 1 / 8
        self.assertAlmostEqual((r - 2 * A) ** 3, (2 * B) ** 3 * r, places=12)
        self.assertAlmostEqual(r, 2 * (A + B * r ** (1 / 3)), places=12)

    def test_the_drawn_hair_has_two_horizons(self):
        # Husain's eq. (20) at k = 3/4, A = L, B = 1/2, late: 1 - 2/r + (1 - B)/r^(3/2) vanishes twice.
        for r in drawings.HAIR_HORIZONS:
            self.assertAlmostEqual(1 - 2 / r + 0.5 / r ** 1.5, 0.0, places=12)
        early = [1 + 1.5 / r ** 1.5 for r in np.linspace(0.01, 4, 200)]
        self.assertTrue(min(early) > 1)

    def test_every_related_spacetime_answers(self):
        mine = {e["id"]: e["kind"] for e in load("metrics")["related"]}
        answer = {"special_case": "generalisation", "family": "family"}
        for other, kind in mine.items():
            theirs = [e for e in load("metrics", other)["related"] if e["id"] == ME]
            self.assertEqual(len(theirs), 1, other)
            self.assertEqual(theirs[0]["kind"], answer[kind], other)


class Charts(unittest.TestCase):
    def einstein(self, system):
        return {tuple(e["indices"]): e["value"] for e in chart(system)["einstein_tensor"]["variants"]["ul"]["nonzero"]}

    def test_the_free_charts_have_wang_and_wus_einstein_tensor(self):
        for system, w in (("ingoing", "v"), ("outgoing", "u")):
            G = self.einstein(system)
            self.assertEqual(G[(w, w)], "-\\dfrac{2\\partial_r m}{r^2}")
            self.assertEqual(G[("r", "r")], "-\\dfrac{2\\partial_r m}{r^2}")
            self.assertEqual(G[("r", w)], f"\\dfrac{{2\\partial_{w} m}}{{r^2}}")
            self.assertEqual(G[("\\theta", "\\theta")], "-\\dfrac{\\partial_r^2 m}{r}")
            self.assertEqual(len(G), 5)

    def test_husains_chart_has_his_fluid(self):
        G = self.einstein("husain")
        self.assertEqual(G[("v", "v")], "-\\dfrac{2h}{r^2}")
        self.assertEqual(G[("\\theta", "\\theta")], "\\dfrac{2k\\,h}{r^2}")
        self.assertEqual(chart("husain")["ricci_scalar"], "R = -\\dfrac{4h\\left(k - 1\\right)}{r^2}")

    def test_every_chart_marks_the_trapped_spheres(self):
        for system in ("ingoing", "outgoing", "husain"):
            self.assertIn("r = 2m \\;\\text{(the marginally trapped spheres)}", chart(system)["domains"])


class Diagrams(unittest.TestCase):
    def test_the_collapse_horizon_ends_on_husains_radius(self):
        view = view_of("ingoing", "collapse")
        event = marker_points(view, "event")
        top = max(event, key=lambda p: p[1])
        self.assertAlmostEqual(top[0], drawings.COLLAPSE_HORIZON, delta=2e-3)
        self.assertLess(min(p[0] for p in event), 0.02)

    def test_the_hairy_hole_horizon_ends_on_its_outer_radius(self):
        view = view_of("husain", "hair")
        event = marker_points(view, "event")
        top = max(event, key=lambda p: p[1])
        self.assertAlmostEqual(top[0], drawings.HAIR_HORIZONS[0], delta=2e-3)
        trapped = [p[0] for p in marker_points(view, "grr")]
        self.assertAlmostEqual(max(trapped), drawings.HAIR_HORIZONS[0], delta=5e-3)
        self.assertAlmostEqual(min(trapped), drawings.HAIR_HORIZONS[1], delta=5e-3)

    def test_the_strings_trap_light_inside_twice_their_mass(self):
        # m = 1/2 + (r^2 - u)/40: the apparent horizon r = 2m at each u of the box.
        view = view_of("outgoing", "strings")
        for X, Y in marker_points(view, "grr"):
            u = Y - X
            self.assertAlmostEqual(X, 2 * (0.5 + (X * X - u) / 40), delta=5e-3)


class Embedding(unittest.TestCase):
    def test_the_movie_is_evenly_spaced_with_distinct_labels(self):
        view = load("embedding")["views"][0]
        frames = view["movie"]["frames"]
        self.assertGreater(len(frames), 30)
        times = [f["value"] for f in frames]
        self.assertTrue(np.allclose(np.diff(times), 1 / 16))
        self.assertEqual(len({f["label"] for f in frames}), len(frames))

    def test_each_moment_marks_its_trapped_sphere(self):
        for T in (-1.0, 0.0, 1.5):
            r = drawings.trapped_on(T)
            self.assertAlmostEqual(r, 2 * drawings.collapse_mass(T + r, r), places=10)
        self.assertAlmostEqual(drawings.trapped_on(1.5), drawings.COLLAPSE_HORIZON, delta=1e-3)


class Curvature(unittest.TestCase):
    def test_no_conformal_diagram_is_drawn(self):
        self.assertFalse((DATA / "conformal" / f"{ME}.json").exists())

    def test_the_ricci_map_is_the_fluids_trace(self):
        view = next(v for v in load("curvature")["views"] if v["id"] == "ricci_scalar")
        m = view["map"]
        nx, ny = m["nx"], m["ny"]
        x0, x1, y0, y1 = view["x"]["from"], view["x"]["to"], view["y"]["from"], view["y"]["to"]
        cells = np.array(m["cells"])
        for i, j in ((0, 0), (nx // 2, ny // 3), (nx - 1, ny - 1)):
            r = x0 + (i + 0.5) * (x1 - x0) / nx
            v = y0 + (j + 0.5) * (y1 - y0) / ny
            want = (1 + math.tanh(2 * v)) * r ** (-8 / 3) / 18
            self.assertAlmostEqual(cells[j, i] / want, 1.0, delta=1e-3)


if __name__ == "__main__":
    unittest.main()
