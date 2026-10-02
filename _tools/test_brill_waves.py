"""Brill's time-symmetric gravitational waves, one moment of a spacetime:
python3 -m unittest discover -s _tools

The physics its texts and drawings state, held on the solver the drawings use and on the published
files. The conformal factor of Holz, Miller, Wakano and Wheeler's wave solves Brill's equation; its
mass is Brill's integral, positive for either sign of the amplitude and of second order in a weak
wave; the masses are those Alcubierre and others tabulate; a minimal surface first surrounds the
centre at the amplitude 11.82, below the Penrose bound of the mass; and the potential of Brill's
equation integrates to zero. The embedding diagram's two views are then held to the metric from the
numbers written. The solver needs numpy and scipy, and its tests are skipped where they are absent.
"""
import json
import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "MFS" / "assets" / "data"
sys.path.insert(0, str(Path(__file__).resolve().parent / "derivations"))
try:
    import numpy as np
    import brill_wave
except ImportError:
    brill_wave = None

CRITICAL = 11.816       # the amplitude at which a minimal surface first surrounds the centre
NECK = 12.19            # and at which the circles of the plane z = 0 first have a narrowest
# Alcubierre and others, Class. Quantum Grav. 17, 2159 (2000), Table V: the amplitude, the mass
# of their axisymmetric solver and its stated error.
TABLE_V = ((1, 0.0338, 0.0004), (2, 0.1262, 0.0009), (5, 0.696, 0.003), (10, 2.912, 0.008), (12, 4.67, 0.01))


@unittest.skipIf(brill_wave is None, "needs numpy and scipy")
class TheConstraint(unittest.TestCase):
    def test_psi_solves_brills_equation(self):
        for a in (1.0, 5.0, 12.0, -3.0):
            wave = brill_wave.wave(a)
            for rho, z in ((0.7, 0.4), (1.5, 1.0), (0.3, 2.2), (2.5, 0.1)):
                self.assertLess(abs(wave.residual(rho, z)), 2e-5 * max(1.0, abs(a)), f"a = {a} at ({rho}, {z})")

    def test_the_potential_integrates_to_zero(self):
        # V = (d_rho^2 q + d_z^2 q)/4 is a quarter of the base metric's scalar curvature density,
        # and its integral over flat space vanishes: positive and negative regions balance.
        r = np.linspace(0, 9, 1801)
        mu = np.linspace(-1, 1, 801)
        R, MU = np.meshgrid(r, mu, indexing="ij")
        V = brill_wave.potential(1.0, R * np.sqrt(1 - MU ** 2), R * MU) * R ** 2
        total = np.trapezoid(np.trapezoid(V, mu, axis=1), r)
        scale = np.trapezoid(np.trapezoid(np.abs(V), mu, axis=1), r)
        self.assertLess(abs(total) / scale, 1e-5)

    def test_the_solution_has_converged(self):
        for a in (5.0, 12.0):
            coarse, fine = brill_wave.Wave(a, 64, 12), brill_wave.Wave(a, 128, 20)
            self.assertLess(abs(coarse.mass / fine.mass - 1), 1e-7)
            self.assertLess(abs(brill_wave.wave(a).mass / fine.mass - 1), 1e-8)
            for rho, z in ((0.5, 0.5), (1.2, 0.0), (0.0, 1.5)):
                self.assertLess(abs(coarse.psi(rho, z)[0] - fine.psi(rho, z)[0]), 1e-7)

    def test_with_no_wave_the_slice_is_flat(self):
        wave = brill_wave.wave(0.0)
        self.assertEqual(wave.mass, 0.0)
        self.assertEqual(wave.psi(0.8, 0.3), (1.0, 0.0, 0.0))


@unittest.skipIf(brill_wave is None, "needs numpy and scipy")
class TheMass(unittest.TestCase):
    def test_the_mass_is_brills_positive_integral(self):
        for a in (-5.0, -1.0, 1.0, 5.0, 12.0):
            wave = brill_wave.wave(a)
            self.assertGreater(wave.mass, 0, f"a = {a}")
            self.assertLess(abs(wave.brill_mass(160, 24) / wave.mass - 1), 1e-7, f"a = {a}")

    def test_a_weak_wave_has_a_mass_of_second_order(self):
        small = [brill_wave.wave(a).mass / a ** 2 for a in (0.01, -0.01, 0.02)]
        self.assertLess(abs(small[0] / small[1] - 1), 2e-2)
        self.assertLess(abs(small[0] / small[2] - 1), 2e-2)
        self.assertGreater(small[0], 0)

    def test_the_masses_are_those_alcubierre_and_others_tabulate(self):
        for a, mass, error in TABLE_V:
            self.assertLess(abs(brill_wave.wave(a).mass - mass), 2 * error, f"a = {a}")

    def test_far_away_psi_is_one_plus_half_the_mass_over_r(self):
        wave = brill_wave.wave(5.0)
        for rho, z in ((300.0, 0.0), (0.0, 300.0), (200.0, 200.0)):
            r = math.hypot(rho, z)
            self.assertLess(abs(2 * r * (wave.psi(rho, z)[0] - 1) / wave.mass - 1), 1e-4)


@unittest.skipIf(brill_wave is None, "needs numpy and scipy")
class TheHorizon(unittest.TestCase):
    def test_a_minimal_surface_first_appears_at_11_82(self):
        self.assertLess(abs(brill_wave.critical() - CRITICAL), 2e-3)
        # Alcubierre and others bracket it between 11.81 and 11.82.
        self.assertIsNone(brill_wave.wave(11.81).horizons())
        self.assertIsNotNone(brill_wave.wave(11.82).horizons())
        self.assertIsNone(brill_wave.wave(4.0).horizons())

    def test_the_horizon_is_minimal_and_obeys_the_penrose_inequality(self):
        wave = brill_wave.wave(12.0)
        inner, outer = wave.horizons()
        self.assertLess(inner, outer)
        for z0 in (inner, outer):
            self.assertLess(abs(wave.miss(z0)), 1e-8)
        area = wave.area(outer)
        # Alcubierre and others give the area as about 1.1e3 at a = 12.
        self.assertLess(abs(area / 1.1e3 - 1), 0.02)
        self.assertLess(area, 16 * math.pi * wave.mass ** 2)
        self.assertLess(16 * math.pi * wave.mass ** 2 / area, 1.01)
        # The outer of the two is a minimum of the area among the curves shot from the axis nearby,
        # once each is closed at the plane z = 0, and the inner one a maximum.
        self.assertGreater(wave.area(inner), area)

    def test_the_marked_circle_is_where_the_horizon_cuts_the_plane(self):
        wave = brill_wave.wave(12.0)
        outer = wave.horizons()[1]
        path, end = wave.curve(outer)
        rho, z, theta = path(end)
        self.assertLess(abs(z), 1e-9)
        self.assertLess(abs(math.cos(theta)), 1e-8)
        self.assertLess(abs(wave.waist(outer) - rho), 1e-12)


class Embedding(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.views = {v["id"]: v for v in json.loads(
            (DATA / "embedding" / "brill_waves.json").read_text(encoding="utf-8"))["views"]}

    def test_the_views_are_these(self):
        self.assertEqual(list(self.views), ["strong", "weak"])
        values = [f["value"] for f in self.views["strong"]["movie"]["frames"]]
        self.assertEqual((values[0], values[-1], len(values)), (11.0, 13.0, 41))
        self.assertEqual(self.views["strong"]["movie"]["variable"], "$a$")
        self.assertNotIn("movie", self.views["weak"])

    def test_the_horizon_is_marked_from_11_82_on(self):
        for frame in self.views["strong"]["movie"]["frames"]:
            a = frame["value"]
            horizon = [ring for ring in frame["rings"] if ring["class"] == "horizon"]
            self.assertEqual(len(horizon), 1 if a > CRITICAL else 0, f"a = {a}")
            if horizon and brill_wave is not None:
                wave = brill_wave.wave(a)
                self.assertLess(abs(wave.waist(wave.horizons()[1]) - horizon[0]["x"]), 1e-7)

    def test_the_plane_has_a_neck_from_12_19_on(self):
        for frame in self.views["strong"]["movie"]["frames"]:
            plane, = frame["pieces"]
            radii = [p[1] for p in plane["points"] if p[0] > 1.5]
            shrinks = any(b < a for a, b in zip(radii, radii[1:]))
            self.assertEqual(shrinks, frame["value"] > NECK, f"a = {frame['value']}")

    @unittest.skipIf(brill_wave is None, "needs numpy and scipy")
    def test_each_frame_is_drawn_to_its_metric(self):
        for frame in self.views["strong"]["movie"]["frames"][::10]:
            a = frame["value"]
            wave = brill_wave.wave(a)
            plane, = frame["pieces"]
            points = np.array(plane["points"])
            self.assertEqual(plane["points"][0], [0.0, 0.0, 0.0])
            x, rho, z = points.T
            psi = wave.equator(x)[0]
            self.assertLess(float(np.max(np.abs(rho - x * psi ** 2))), 2e-6, f"a = {a}")
            # Each chord is as long as the metric says, psi^2 e^q drho along the meridian, by Simpson's rule.
            def stretch(r):
                return wave.equator(r)[0] ** 2 * np.exp(a * r ** 2 * np.exp(-r ** 2))
            proper = (x[1:] - x[:-1]) * (stretch(x[:-1]) + 4 * stretch((x[:-1] + x[1:]) / 2) + stretch(x[1:])) / 6
            chords = np.hypot(np.diff(rho), np.diff(z))
            self.assertLess(float(np.max(np.abs(chords / proper - 1))), 3e-4, f"a = {a}")

    @unittest.skipIf(brill_wave is None, "needs numpy and scipy")
    def test_the_stalk_is_as_long_as_the_captions_say(self):
        frames = self.views["strong"]["movie"]["frames"]
        for frame, height, mass, stretch in ((frames[0], 38.0, 3.69, 57.0), (frames[-1], 89.1, 5.94, 119.0)):
            a = frame["value"]
            self.assertLess(abs(frame["pieces"][0]["points"][-1][2] - height), 0.1)
            self.assertLess(abs(brill_wave.wave(a).mass - mass), 5e-3)
            self.assertLess(abs(math.exp(a / math.e) - stretch), 0.5)

    @unittest.skipIf(brill_wave is None, "needs numpy and scipy")
    def test_the_weaker_wave_stops_where_its_circles_outgrow_the_distance(self):
        a = 4.0
        wave = brill_wave.wave(a)
        self.assertLess(abs(wave.mass - 0.460), 5e-4)
        disc, outside = self.views["weak"]["surfaces"][0]["pieces"]
        inner, outer = disc["points"][-1][0], outside["points"][0][0]
        self.assertLess(abs(inner - 1.59), 5e-3)
        self.assertLess(abs(outer - 2.43), 5e-3)

        def spare(r):
            psi, dpsi = (float(v) for v in wave.equator(r))
            return math.exp(2 * a * r * r * math.exp(-r * r)) - (1 + 2 * r * dpsi / psi) ** 2
        for r in np.linspace(inner + 1e-3, outer - 1e-3, 50):
            self.assertLess(spare(r), 0)
        for r in (0.5, 1.0, 1.5, 2.5, 4.0, 6.0):
            self.assertGreater(spare(r), 0)
        # The two circles where the surface stops are set at one height.
        self.assertLess(abs(disc["points"][-1][2] - outside["points"][0][2]), 1e-6)


if __name__ == "__main__":
    unittest.main()
