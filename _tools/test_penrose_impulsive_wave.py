"""Penrose's spherical impulsive wave for a cosmic string that snaps:
python3 -m unittest discover -s _tools

What its texts and drawings state, held on the published files. The continuous charts are held to
a curvature that is a delta on the wave front and nothing else, and, where sympy and numpy are
present, to being the string's cone ahead of the front and Minkowski space behind it, through
Podolsky and Griffiths's transformations (3) to (5) of Class. Quantum Grav. 17, 1401 (2000) written
out here afresh for the warp h = Z^beta. The movie of the embedding diagram is held to the flat disc
and the cone it draws, to their meeting on the front, and to the ring of free particles, which rests
on R = 2 until the front reaches it and then falls at 3c/5, Podolsky and Steinbauer's
delta (1 - delta/2)/(1 - delta) = 3/4 per unit proper time at delta = 1/2.
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
READY = all(importlib.util.find_spec(name) is not None for name in ("sympy", "numpy"))
sys.path.insert(0, str(Path(__file__).resolve().parent / "derivations"))
BETA = 0.5


def chart(system):
    metric = json.loads((DATA / "metrics" / "penrose_impulsive_wave.json").read_text(encoding="utf-8"))
    return next(c for c in metric["coordinates"] if c["id"] == system)


class Curvature(unittest.TestCase):
    def test_the_curvature_of_each_continuous_chart_is_a_delta_on_the_front(self):
        for system, front in (("null", "\\delta(U)"), ("retarded", "\\delta(u)")):
            entry = chart(system)
            for tensor in ("riemann", "weyl_tensor"):
                for variant in entry[tensor]["variants"].values():
                    self.assertEqual(len(variant["nonzero"]), 8, (system, tensor))
                    for component in variant["nonzero"]:
                        self.assertIn(front, component["value"], (system, tensor, component["indices"]))
            for tensor in ("ricci_tensor", "einstein_tensor"):
                for variant in entry[tensor]["variants"].values():
                    self.assertEqual(variant["nonzero"], [], (system, tensor))
            self.assertEqual(entry["ricci_scalar"], "R = 0")
            self.assertEqual(entry["kretschmann"], "K = 0")

    def test_the_flat_charts_have_no_curvature(self):
        for system in ("behind", "ahead"):
            entry = chart(system)
            for tensor in ("riemann", "ricci_tensor", "einstein_tensor", "weyl_tensor"):
                for variant in entry[tensor]["variants"].values():
                    self.assertEqual(variant["nonzero"], [], (system, tensor))


@unittest.skipUnless(READY, "sympy and numpy are not installed")
class Gluing(unittest.TestCase):
    """Each side of each continuous chart against the flat space it is."""

    def flat(self, epsilon, U, V, rho, phi, beta, ahead):
        """Inertial (cT, X, Y, Z) of the event (U, V, rho e^{i phi}): Podolsky and Griffiths's (3) behind
        the wave, and (4) and (5) with h = Z^beta ahead of it, where the angle of eta is beta phi."""
        p = 1 + epsilon * rho ** 2
        if not ahead:
            u, v, eta, angle = rho ** 2 * V / p - U, V / p - epsilon * U, V * rho / p, phi
        else:
            k = (1 - beta ** 2) / 4
            A, B, C = rho ** (1 - beta) / (beta * p), rho ** (1 + beta) / (beta * p), rho / (beta * p)
            D = rho ** (1 - beta) * (p * (1 - beta) ** 2 / (4 * rho ** 2) + epsilon * beta) / beta
            E = rho ** (1 + beta) * (p * (1 + beta) ** 2 / (4 * rho ** 2) - epsilon * beta) / beta
            F = p * k / (beta * rho)
            u, v, eta, angle = B * V - E * U, A * V - D * U, C * V - F * U, beta * phi
        root = math.sqrt(2)
        return [(u + v) / root, root * eta * math.cos(angle), root * eta * math.sin(angle), (u - v) / root]

    def published(self, system, point, beta):
        import sympy as sp
        import verify_metrics as vm
        entry = chart(system)
        reader = vm.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]], ())
        at = dict(zip([reader.symbol[c] for c in entry["coords"]], point))
        at[reader.parameters["k"]] = (1 - beta ** 2) / 4
        g = [[0.0] * 4 for _ in range(4)]
        for e in entry["metric_components"]:
            i, j = (entry["coords"].index(x) for x in e["indices"])
            g[i][j] = float(sp.sympify(reader(e["value"])).subs(at))
        return g

    def pulled_back(self, image, point, h=1e-6):
        columns = []
        for j in range(4):
            up, down = list(point), list(point)
            up[j] += h
            down[j] -= h
            columns.append([(a - b) / (2 * h) for a, b in zip(image(up), image(down))])
        eta = (-1, 1, 1, 1)
        return [[sum(eta[m] * columns[i][m] * columns[j][m] for m in range(4)) for j in range(4)] for i in range(4)]

    def test_each_side_of_the_null_chart_is_flat_space_and_the_cone(self):
        rng = random.Random(3)
        for ahead in (False, True):
            for _ in range(5):
                beta = rng.uniform(0.3, 0.9)
                V, rho, phi = rng.uniform(2, 5), rng.uniform(1, 3), rng.uniform(0.1, 6)
                U = rng.uniform(0.2, 2) * (1 if ahead else -1)
                mine = self.published("null", [U, V, rho, phi], beta)
                flat = self.pulled_back(lambda x: self.flat(0, *x, beta, ahead), [U, V, rho, phi])
                for i in range(4):
                    for j in range(4):
                        self.assertAlmostEqual(mine[i][j], flat[i][j], places=6, msg=(ahead, i, j))

    def test_each_side_of_the_retarded_chart_is_flat_space_and_the_cone(self):
        rng = random.Random(4)

        def image(x, beta, ahead):
            u, r, theta, phi = x
            return self.flat(1, -u / math.sqrt(2), math.sqrt(2) * r, 1 / math.tan(theta / 2), phi, beta, ahead)

        for ahead in (False, True):
            for _ in range(5):
                beta = rng.uniform(0.3, 0.9)
                r, theta, phi = rng.uniform(3, 6), rng.uniform(0.9, 2.2), rng.uniform(0.1, 6)
                u = rng.uniform(0.2, 2) * (-1 if ahead else 1)
                mine = self.published("retarded", [u, r, theta, phi], beta)
                flat = self.pulled_back(lambda x: image(x, beta, ahead), [u, r, theta, phi])
                for i in range(4):
                    for j in range(4):
                        self.assertAlmostEqual(mine[i][j], flat[i][j], places=6, msg=(ahead, i, j))

    def test_the_retarded_chart_reads_the_inertial_time_and_radius_behind_the_wave(self):
        for u, r, theta in ((0.5, 2.0, 1.0), (1.5, 0.7, 2.0)):
            cT, X, Y, Z = self.flat(1, -u / math.sqrt(2), math.sqrt(2) * r, 1 / math.tan(theta / 2), 0.0, BETA, False)
            self.assertAlmostEqual(cT, u + r)
            self.assertAlmostEqual(math.hypot(X, Z), r)
            self.assertAlmostEqual(Z / r, math.cos(theta))

    def test_on_the_front_the_two_frames_meet_where_T_is_t_over_beta(self):
        for ct in (0.5, 1.0, 2.5):
            cT, X, Y, Z = self.flat(1, 0.0, math.sqrt(2) * ct, 1.0, 0.0, BETA, True)
            self.assertAlmostEqual(cT, ct / BETA)
            self.assertAlmostEqual(X, ct / BETA)
            self.assertAlmostEqual(Z, 0.0)


class Movie(unittest.TestCase):
    def setUp(self):
        self.view = json.loads((DATA / "embedding" / "penrose_impulsive_wave.json").read_text(encoding="utf-8"))["views"][0]
        self.frames = self.view["movie"]["frames"]

    def test_every_frame_is_a_level_disc_joined_to_the_strings_cone_on_the_front(self):
        slope = math.sqrt(1 - BETA ** 2)
        for frame in self.frames:
            ct = frame["value"]
            inside, outside = frame["pieces"]
            self.assertEqual((inside["system"], outside["system"]), ("behind", "ahead"))
            level = inside["points"][0][2]
            for r, rho, z in inside["points"]:
                self.assertAlmostEqual(rho, r, places=6)
                self.assertAlmostEqual(z, level, places=6)
            self.assertAlmostEqual(inside["points"][-1][0], ct, places=6)
            self.assertAlmostEqual(outside["points"][0][0], ct / BETA, places=6)
            for R, rho, z in outside["points"]:
                self.assertAlmostEqual(rho, BETA * R, places=6)
                self.assertAlmostEqual(z, slope * (R - 8), places=6)
            # One circle from both sides: the front has the circumference 2 pi ct.
            self.assertAlmostEqual(inside["points"][-1][1], outside["points"][0][1], places=6)
            self.assertAlmostEqual(inside["points"][-1][2], outside["points"][0][2], places=6)

    def test_the_ring_rests_until_the_front_reaches_it_and_then_falls_at_three_fifths_of_c(self):
        speed = (1 - BETA ** 2) / (1 + BETA ** 2)
        self.assertAlmostEqual(speed, 0.6)
        delta = 1 - BETA
        self.assertAlmostEqual(speed / math.sqrt(1 - speed ** 2), delta * (1 - delta / 2) / (1 - delta))
        for frame in self.frames:
            ct = frame["value"]
            rings = [r for r in frame["rings"] if r["class"] == "particles"]
            self.assertEqual(len(rings), 1, ct)
            ring = rings[0]
            if ct <= 1:
                self.assertEqual(ring["piece"], "outside")
                self.assertAlmostEqual(ring["x"], 2.0)
                self.assertAlmostEqual(ring["rho"], 1.0, places=6)
            else:
                self.assertEqual(ring["piece"], "inside")
                self.assertAlmostEqual(ring["x"], 1 - speed * (ct - 1), places=6)
                self.assertGreater(ring["x"], 0)

    def test_the_movie_runs_from_a_quarter_to_five_halves(self):
        self.assertAlmostEqual(self.frames[0]["value"], 0.25)
        self.assertAlmostEqual(self.frames[-1]["value"], 2.5)


if __name__ == "__main__":
    unittest.main()
