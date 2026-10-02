"""The T-models of Datt and Ruban, dust on a tube of spheres that all share one radius:
python3 -m unittest discover -s _tools

What its texts and its drawings state, held on the published files. The scale factor of
Ruban's chart is held to the numbers the History and the captions print: positive from the
bang to the crunch exactly where the dust is thick enough, mu >= 1/(2 pi) at epsilon = 1, and
growing toward both ends as the captions say. The embedding diagram's moments are held to the
closed forms of the geometry: the dust a tube of the radius of its cycloid whose marked shells
draw apart as a grows, the outside Flamm's paraboloid at the greatest expansion, the two
meeting in one point, and the rest mass of the stretch drawn 4/pi times the mass seen from
outside. The conformal diagram's corners are held to the conformal time of the dust. The
field equations, the density and the junction need sympy and are skipped where it is absent.
"""
import importlib.util
import json
import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "MFS" / "assets" / "data"
NAME = "datt_ruban_t_models"
HAS_SYMPY = importlib.util.find_spec("sympy") is not None
sys.path.insert(0, str(Path(__file__).resolve().parent / "derivations"))

MU = 1 / math.pi                    # the dust on every shell of the T-sphere drawn, with epsilon = 1


def load(folder):
    return json.loads((DATA / folder / f"{NAME}.json").read_text(encoding="utf-8"))


def scale(eta, mu, epsilon=1):
    """Ruban's a = epsilon cot(eta/2) + 2 mu (1 - (eta/2) cot(eta/2))."""
    return (epsilon - mu * eta) / math.tan(eta / 2) + 2 * mu


def cycloid(tau):
    """e in [0, pi] with (e + sin e)/2 = tau: the dust's parameter since its greatest expansion."""
    lo, hi = 0.0, math.pi
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if (mid + math.sin(mid)) / 2 < tau else (lo, mid)
    return 0.5 * (lo + hi)


class ScaleFactor(unittest.TestCase):
    ETAS = [2 * math.pi * k / 4000 for k in range(1, 4000)]

    def test_the_shells_keep_apart_exactly_where_the_dust_is_thick_enough(self):
        # The History: with epsilon = 1 the shells keep apart from the bang to the crunch only
        # where mu > 1/(2 pi). Below that a reaches zero before the crunch.
        for mu in (1 / (2 * math.pi) + 1e-6, 0.2, MU, 0.5, 3.0):
            self.assertTrue(all(scale(e, mu) > 0 for e in self.ETAS), f"mu = {mu}")
        for mu in (0.01, 0.1, 1 / (2 * math.pi) - 1e-3):
            self.assertTrue(any(scale(e, mu) < 0 for e in self.ETAS), f"mu = {mu}")

    def test_without_the_vacuum_part_the_tube_starts_from_a_point_and_never_crosses(self):
        for mu in (0.05, MU, 2.0):
            self.assertTrue(all(scale(e, mu, 0) > 0 for e in self.ETAS))
            self.assertAlmostEqual(scale(1e-3, mu, 0) / (mu * 1e-6 / 6), 1.0, delta=1e-3)

    def test_the_stretch_grows_toward_the_bang_and_the_crunch_as_the_caption_says(self):
        # a -> 2/eta at the bang and 2 (2 pi mu - 1)/(2 pi - eta) at the crunch.
        caption = " ".join(load("diagrams")["systems"]["ruban"][0]["caption"])
        self.assertIn("2/\\eta", caption)
        self.assertIn("2\\left(2\\pi\\mu - 1\\right)/\\left(2\\pi - \\eta\\right)", caption)
        for mu in (0.2, MU, 0.45):
            self.assertAlmostEqual(scale(1e-4, mu) * 1e-4 / 2, 1.0, delta=1e-3)
            x = 1e-4
            self.assertAlmostEqual(scale(2 * math.pi - x, mu) * x / (2 * (2 * math.pi * mu - 1)), 1.0, delta=1e-3)

    def test_the_drawn_tube_is_above_the_threshold_on_every_shell(self):
        # mu = (2 + tanh r)/(2 pi), the dust of every spacetime diagram of the tube.
        for k in range(-60, 61):
            mu = (2 + math.tanh(k / 10)) / (2 * math.pi)
            self.assertGreater(mu, 1 / (2 * math.pi))
            self.assertTrue(all(scale(e, mu) > 0 for e in self.ETAS[::40]))

    def test_the_t_sphere_drawn_is_symmetric_in_time(self):
        for e in (0.3, 1.0, 2.5):
            self.assertAlmostEqual(scale(math.pi + e, MU), scale(math.pi - e, MU), places=12)
        self.assertAlmostEqual(scale(math.pi, MU), 2 / math.pi, places=12)


class Embedding(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.view = load("embedding")["views"][0]
        cls.frames = cls.view["movie"]["frames"]

    @staticmethod
    def piece(frame, pid):
        return next(p for p in frame["pieces"] if p["id"] == pid)["points"]

    def test_the_dust_is_a_tube_of_its_cycloids_radius_whose_shells_draw_apart(self):
        lengths = []
        for frame in self.frames:
            tau = frame["value"]
            e = cycloid(tau)
            a, b = scale(math.pi + e, MU), math.cos(e / 2) ** 2
            dust = self.piece(frame, "dust")
            self.assertEqual(dust[0][0], -2)
            self.assertEqual(dust[-1][0], 0)
            for r, rho, z in dust:
                self.assertAlmostEqual(rho, b, delta=2e-5, msg=f"rho at {r}, c tau = {tau}")
                self.assertAlmostEqual(z - dust[-1][2], a * r, delta=2e-4, msg=f"z at {r}, c tau = {tau}")
            lengths.append(dust[-1][2] - dust[0][2])
        self.assertTrue(all(y > x for x, y in zip(lengths, lengths[1:])), "the tube does not lengthen")
        self.assertAlmostEqual(lengths[0], 4 / math.pi, delta=2e-4)

    def test_at_the_greatest_expansion_the_outside_is_flamms_paraboloid_on_the_tube(self):
        outside = self.piece(self.frames[0], "exterior")
        dust = self.piece(self.frames[0], "dust")
        self.assertEqual(outside[0][0], 0)
        for s, rho, z in outside:
            self.assertAlmostEqual(rho, s * s + 1, delta=2e-5, msg=f"rho at {s}")
            self.assertAlmostEqual(z - outside[0][2], 2 * s, delta=2e-4, msg=f"z at {s}")
        self.assertAlmostEqual(dust[-1][1], 1.0, delta=1e-6)       # as wide as the horizon

    def test_the_surface_is_the_smallest_circle_outside_and_the_two_meet_there(self):
        radii = []
        for frame in self.frames:
            tau = frame["value"]
            outside, dust = self.piece(frame, "exterior"), self.piece(frame, "dust")
            self.assertAlmostEqual(outside[0][1], math.cos(cycloid(tau) / 2) ** 2, delta=2e-5, msg=f"c tau = {tau}")
            self.assertAlmostEqual(outside[0][1], min(p[1] for p in outside), delta=1e-9)
            self.assertAlmostEqual(outside[0][1], dust[-1][1], delta=2e-5)
            self.assertAlmostEqual(outside[0][2], dust[-1][2], delta=2e-5)
            radii.append(outside[0][1])
        self.assertTrue(all(y < x for x, y in zip(radii, radii[1:])), "the tube does not narrow")
        self.assertLess(self.frames[-1]["value"], math.pi / 2)     # it would close at c tau = pi r_s/2

    def test_the_stretch_drawn_holds_four_over_pi_of_the_mass_seen_outside(self):
        # The rest mass between two shells is (c^2/G) mu dr, and M = r_s c^2/(2G): the ratio is 2 mu dr/r_s.
        caption = " ".join(self.view["caption"])
        self.assertIn("rest mass $4M/\\pi$", caption)
        self.assertIn("adds $2M/\\pi$", caption)
        self.assertAlmostEqual(2 * MU * 2.0, 4 / math.pi, places=12)
        self.assertAlmostEqual(2 * MU * 1.0, 2 / math.pi, places=12)


class Conformal(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.views = {v["id"]: v for v in load("conformal")["views"]}

    @staticmethod
    def sigma_m():
        """The conformal time of the dust from its greatest expansion to the crunch, in r_s, by Simpson's rule."""
        n = 200000
        h = math.pi / n
        f = lambda e: math.sin(e / 2) ** 2 / scale(e, MU) if e < 2 * math.pi else 0.0
        total = f(math.pi) + f(2 * math.pi)
        for k in range(1, n):
            total += (4 if k % 2 else 2) * f(math.pi + k * h)
        return total * h / 3

    def test_the_corners_stand_at_the_conformal_time_of_the_dust(self):
        sigma = self.sigma_m()
        self.assertAlmostEqual(sigma, 1.9147, delta=1e-4)
        top = 2 * math.atan(sigma)
        for view in self.views.values():
            points = [layer["at"] for layer in view["layers"] if layer["kind"] == "point"]
            self.assertIn([round(top, 4), round(top, 4)], points)             # i+
            self.assertIn([round(2 * top, 4), 0], points)                     # i0
            self.assertIn([round(-math.pi, 4), 0], points)                    # the far end of the tube

    def test_the_captions_numbers_are_the_dusts(self):
        sigma = self.sigma_m()
        text = " ".join(self.views["ruban"]["caption"])
        self.assertIn("$3.83\\,r_s$", text)
        self.assertIn("$r = -1.91\\,r_s$", text)
        self.assertAlmostEqual(2 * sigma, 3.83, delta=5e-3)
        self.assertAlmostEqual(sigma, 1.91, delta=5e-3)

    def test_the_horizons_cross_on_the_surface_at_the_greatest_expansion(self):
        for view in self.views.values():
            lines = {layer["class"]: layer["points"] for layer in view["layers"] if layer["kind"] == "line"
                     and layer["class"] in ("event", "horizon")}
            for cls, sign in (("event", 1), ("horizon", -1)):
                (x0, t0), (x1, t1) = lines[cls][0], lines[cls][-1]
                self.assertAlmostEqual(abs((t1 - t0) / (x1 - x0)), 1.0, places=6)       # a light ray
                self.assertAlmostEqual(t0 - sign * x0, 0.0, places=4)                   # through the origin


@unittest.skipUnless(HAS_SYMPY, "sympy is not installed")
class FieldEquations(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import sympy as sp
        import verify_metrics as vm
        cls.sp, cls.vm = sp, vm
        cls.charts = {c["id"]: c for c in json.loads(
            (DATA / "metrics" / f"{NAME}.json").read_text(encoding="utf-8"))["coordinates"]}

    def chart(self, system):
        entry = self.charts[system]
        reader = self.vm.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]], ())
        mixed = {tuple(e["indices"]): reader(e["value"]) for e in entry["einstein_tensor"]["variants"]["ul"]["nonzero"]}
        return entry, reader, mixed

    def test_each_explicit_chart_is_dust_at_rest_of_density_two_mu_over_a_b_squared(self):
        sp = self.sp
        for system, radius in (("ruban", lambda R: R.parameters["r_s"] * sp.sin(R.symbol["\\eta"] / 2) ** 2),
                               ("areal", lambda R: R.symbol["T"])):
            entry, R, mixed = self.chart(system)
            time = entry["coords"][0]
            self.assertEqual(set(mixed), {(time, time)}, system)
            a, mu = R.parameters["a"], R.parameters["mu"]
            self.assertEqual(sp.simplify(mixed[(time, time)] + 2 * mu / (a * radius(R) ** 2)), 0, system)

    def test_the_de_sitter_chart_is_dust_and_a_cosmological_constant(self):
        sp = self.sp
        entry, R, mixed = self.chart("de_sitter")
        ell, a, mu, t = R.parameters["ell"], R.parameters["a"], R.parameters["mu"], R.symbol["t"]
        b = ell * sp.cosh(R.c * t / ell)
        for name in entry["coords"]:
            want = -3 / ell ** 2 - (2 * mu / (a * b ** 2) if name == "t" else 0)
            self.assertEqual(sp.simplify((mixed[(name, name)] - want).rewrite(sp.exp)), 0, name)

    def test_rubans_scale_factor_carries_twice_mu(self):
        # Krasinski's factor of 2: with mu alone in front of 1 - (eta/2) cot(eta/2), as Ruban's (18)
        # prints it, the published G^r_r of the comoving chart does not vanish.
        sp = self.sp
        eta, mu = sp.symbols("eta mu", positive=True)
        b = sp.sin(eta / 2) ** 2
        d = lambda f: sp.diff(f, eta) / b                              # d/d(ct) at r_s = 1
        for factor, vanishes in ((2, True), (1, False)):
            a = sp.cot(eta / 2) + factor * mu * (1 - eta / 2 * sp.cot(eta / 2))
            density = sp.simplify((b * d(b) ** 2 * a + a * b + 2 * b ** 2 * d(a) * d(b)) / (a * b ** 3) - 2 * mu / (a * b ** 2))
            pressure = sp.simplify(d(d(a)) / a + d(d(b)) / b + d(a) * d(b) / (a * b))
            self.assertEqual(pressure, 0)
            self.assertEqual(density == 0, vanishes, f"factor {factor}")

    def test_a_surface_of_constant_r_is_the_vacuums_whatever_the_dust(self):
        # The junction of a T-sphere: in the areal chart the block of T, theta and phi holds no a,
        # mu or epsilon, and does not change along r, so the induced metric is that of the inside
        # of Schwarzschild's horizon and the extrinsic curvature (1/2a) d_r g vanishes.
        entry, R, _ = self.chart("areal")
        r = R.symbol["r"]
        for e in entry["metric_components"]:
            if "r" in e["indices"]:
                continue
            value = R(e["value"])
            self.assertFalse(value.has(R.parameters["a"], R.parameters["mu"], R.parameters["epsilon"]), e["indices"])
            self.assertEqual(self.sp.diff(value, r), 0)


if __name__ == "__main__":
    unittest.main()
