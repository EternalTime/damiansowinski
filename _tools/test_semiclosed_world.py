"""The semiclosed world, more than half of a closed universe of dust behind Schwarzschild's throat:
.venv.noindex/bin/python -m unittest discover -s _tools

What its texts and its drawings state, held on the published files. The embedding diagram's
moments are held to the closed forms of the geometry: the dust a sphere of the radius its cycloid
gives, the outside Flamm's paraboloid on both sheets at the moment of greatest expansion, the
throat on its own cycloid afterwards, and the two meeting in one point with one tangent. The
conformal diagram's horizons and corners are held to the angles chi_0 fixes, and the numbers the
captions and the History print are computed again. The junction and the dust's field equation
need sympy and are skipped where it is absent.
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

CHI0 = 3 * math.pi / 4
AM = 2 * math.sqrt(2.0)             # a_m in r_s: a_m sin^3 chi_0 = r_s


def load(folder):
    return json.loads((DATA / folder / "semiclosed_world.json").read_text(encoding="utf-8"))


def cycloid(target):
    """eta in [0, pi] with eta + sin eta = target."""
    lo, hi = 0.0, math.pi
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if mid + math.sin(mid) < target else (lo, mid)
    return 0.5 * (lo + hi)


class Mass(unittest.TestCase):
    def test_the_drawn_world_has_the_mass_and_the_surface_its_texts_state(self):
        # r_s = a_m sin^3 chi_0, and the surface reaches a_m sin chi_0 = 2 r_s.
        self.assertAlmostEqual(AM * math.sin(CHI0) ** 3, 1.0, places=12)
        self.assertAlmostEqual(AM * math.sin(CHI0), 2.0, places=12)

    def test_the_mass_seen_outside_falls_past_the_equator_and_vanishes_as_the_world_closes(self):
        mass = lambda chi: math.sin(chi) ** 3 / 2                      # G M / (c^2 a_m)
        rest = lambda chi: 0.75 * (chi - math.sin(chi) * math.cos(chi))  # G M_0 / (c^2 a_m)
        angles = [math.pi / 2 + k * (math.pi / 2) / 50 for k in range(51)]
        for a, b in zip(angles, angles[1:]):
            self.assertLess(mass(b), mass(a), "the mass seen from outside does not fall")
            self.assertGreater(rest(b), rest(a), "the dust within does not grow")
        self.assertAlmostEqual(mass(math.pi), 0.0, places=12)
        # The embedding's caption: a twelfth of the mass of the dust counted grain by grain.
        self.assertAlmostEqual(rest(CHI0) / mass(CHI0), 12.1, delta=0.05)
        caption = " ".join(load("embedding")["views"][0]["caption"])
        self.assertIn("a twelfth of the mass of the dust", caption)


class Embedding(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.view = load("embedding")["views"][0]
        cls.frames = cls.view["movie"]["frames"]

    @staticmethod
    def piece(frame, pid):
        return next(p for p in frame["pieces"] if p["id"] == pid)["points"]

    def test_the_dust_is_a_sphere_of_its_cycloids_radius_kept_past_its_equator(self):
        for frame in self.frames:
            tau = frame["value"]
            eta = cycloid(2 * tau / AM)
            a = AM * (1 + math.cos(eta)) / 2
            dust = self.piece(frame, "dust")
            self.assertEqual(dust[0][0], 0)
            self.assertAlmostEqual(dust[-1][0], CHI0, places=12)
            for chi, rho, z in dust:
                self.assertAlmostEqual(rho, a * math.sin(chi), delta=2e-5, msg=f"rho at {chi}, c tau = {tau}")
                self.assertAlmostEqual(z - dust[0][2], a * (1 - math.cos(chi)), delta=2e-4, msg=f"z at {chi}, c tau = {tau}")
            # More than half of the three sphere: the widest circle is the equator, inside the dust.
            widest = max(dust, key=lambda p: p[1])
            self.assertAlmostEqual(widest[1], a, delta=1e-3)
            self.assertLess(widest[0], CHI0)

    def test_at_the_greatest_expansion_the_outside_is_flamms_paraboloid_on_both_sheets(self):
        outside = self.piece(self.frames[0], "exterior")
        throat = next(p for p in outside if p[0] == 0)
        self.assertAlmostEqual(outside[0][0], 1 / math.tan(CHI0), places=12)   # the surface, s = cot chi_0 = -1
        self.assertLess(outside[0][0], 0)
        self.assertGreater(outside[-1][0], 0)
        for s, rho, z in outside:
            self.assertAlmostEqual(rho, s * s + 1, delta=2e-5, msg=f"rho at {s}")
            self.assertAlmostEqual(z - throat[2], 2 * s, delta=2e-4, msg=f"z at {s}")
        self.assertAlmostEqual(throat[1], 1.0, delta=1e-6)

    def test_the_throat_follows_its_own_cycloid_and_is_the_smallest_circle_outside(self):
        radii = []
        for frame in self.frames:
            tau = frame["value"]
            outside = self.piece(frame, "exterior")
            throat = next(p for p in outside if p[0] == 0)
            self.assertAlmostEqual(throat[1], math.cos(cycloid(2 * tau) / 2) ** 2, delta=2e-5, msg=f"c tau = {tau}")
            self.assertAlmostEqual(throat[1], min(p[1] for p in outside), delta=1e-9)
            radii.append(throat[1])
        self.assertTrue(all(b < a for a, b in zip(radii, radii[1:])), "the throat does not shrink")
        # It would close at c tau = pi r_s/2, after the last moment drawn.
        self.assertLess(self.frames[-1]["value"], math.pi / 2)
        self.assertAlmostEqual(cycloid(2 * (math.pi / 2)), math.pi, delta=1e-4)

    def test_the_dust_meets_the_outside_in_one_point_with_one_tangent(self):
        for frame in self.frames:
            dust, outside = self.piece(frame, "dust"), self.piece(frame, "exterior")
            self.assertAlmostEqual(dust[-1][1], outside[0][1], delta=2e-5)
            self.assertAlmostEqual(dust[-1][2], outside[0][2], delta=2e-5)
            # Both sides leave the join along (cos chi_0, sin chi_0): the circles shrink outward.
            for a, b in ((dust[-2], dust[-1]), (outside[0], outside[1])):
                length = math.hypot(b[1] - a[1], b[2] - a[2])
                self.assertAlmostEqual((b[1] - a[1]) / length, math.cos(CHI0), delta=0.03)
                self.assertAlmostEqual((b[2] - a[2]) / length, math.sin(CHI0), delta=0.03)

    def test_the_rim_stands_still_and_the_bag_hangs_below_it(self):
        for frame in self.frames:
            outside, dust = self.piece(frame, "exterior"), self.piece(frame, "dust")
            self.assertEqual(outside[-1][2], 0, "the rim of the drawing stands at z = 0")
            self.assertTrue(all(b[2] > a[2] for a, b in zip(outside, outside[1:])), "the outside does not rise to its rim")
            self.assertLess(dust[-1][2], 0)

    def test_the_marked_circles_are_the_marginally_trapped_spheres(self):
        for frame in self.frames[1:]:
            tau = frame["value"]
            eta = cycloid(2 * tau / AM)
            marked = [r for r in frame["rings"] if r["class"] == "horizon"]
            outside = sorted(r["x"] for r in marked if r["piece"] == "exterior")
            inside = sorted(r["x"] for r in marked if r["piece"] == "dust")
            self.assertEqual(len(outside), 2)
            self.assertAlmostEqual(outside[0], -outside[1], places=9)
            for r in marked:
                if r["piece"] == "exterior":
                    self.assertAlmostEqual(r["rho"], 1.0, delta=2e-5, msg="a marked circle outside is not at r_s")
            self.assertAlmostEqual(inside[0], math.pi / 2 - eta / 2, places=6)
            self.assertAlmostEqual(inside[1], math.pi / 2 + eta / 2, places=6)
            # There the areal radius a sin chi is 2GM/c^2 for the mass within, a_m sin^3 chi.
            a = AM * (1 + math.cos(eta)) / 2
            for chi in inside:
                self.assertAlmostEqual(a * math.sin(chi), AM * math.sin(chi) ** 3, places=6)

    def test_the_numbers_the_caption_prints(self):
        caption = " ".join(self.view["caption"])
        last = self.frames[-1]["value"]
        self.assertAlmostEqual(last, 1.5, places=9)
        self.assertIn("$c\\tau = 1.50\\,r_s$", caption)
        a = (1 + math.cos(cycloid(2 * last / AM))) / 2
        self.assertEqual(f"{a:.2f}", "0.93")
        self.assertIn("$0.93\\,a_m$", caption)
        self.assertIn("$c\\tau = \\pi r_s/2$", caption)


class Conformal(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.views = {v["id"]: v for v in load("conformal")["views"]}

    def lines(self, view, cls):
        return [layer["points"] for layer in self.views[view]["layers"] if layer["class"] == cls and "points" in layer
                and layer["kind"] != "fill"]

    def test_one_view_for_each_chart(self):
        metric = json.loads((DATA / "metrics" / "semiclosed_world.json").read_text(encoding="utf-8"))
        self.assertEqual(list(self.views), [c["id"] for c in metric["coordinates"]])
        for vid, view in self.views.items():
            self.assertEqual(view["system"], vid)

    def test_the_horizons_and_corners_stand_where_chi_0_puts_them(self):
        near = lambda got, want: self.assertLess(math.dist(got, want), 2e-4, f"{got} is not {want}")
        for vid in self.views:
            event, = self.lines(vid, "event")
            near(event[0], [3 * CHI0 - 2 * math.pi, -math.pi])      # the bang, at chi = 3 chi_0 - 2 pi
            near(event[-1], [3 * CHI0, math.pi])                    # i+
            # A null line, and it meets the surface chi = chi_0 at eta = pi - 2 chi_0.
            self.assertAlmostEqual(event[-1][1] - event[0][1], event[-1][0] - event[0][0], delta=2e-4)
            self.assertAlmostEqual(event[0][1] + (CHI0 - event[0][0]), math.pi - 2 * CHI0, delta=2e-4)
            white, = self.lines(vid, "horizon")
            near(white[0], [3 * CHI0 - 2 * math.pi, math.pi])
            near(white[-1], [3 * CHI0, -math.pi])
            scri = self.lines(vid, "scri")
            for line in scri:
                near(line[-1], [math.pi + 3 * CHI0, 0])             # i0
            surface, = self.lines(vid, "surface")
            self.assertTrue(all(abs(p[0] - CHI0) < 2e-4 for p in surface))
            up, down = self.lines(vid, "apparent")
            for line, sign in ((up, 1), (down, -1)):
                near(line[0], [0, sign * math.pi])
                near(line[1], [math.pi / 2, 0])                     # they cross on the equator at eta = 0
                near(line[2], [CHI0, sign * (2 * CHI0 - math.pi)])  # and meet the surface where it crosses r_s

    def test_the_centre_of_the_dust_lies_behind_the_event_horizon(self):
        # The event horizon is eta - chi = pi - 3 chi_0, and it reaches the bang at chi = pi/4 > 0,
        # so no event of the centre lies below it.
        self.assertGreater(3 * CHI0 - 2 * math.pi, 0)
        self.assertAlmostEqual(3 * CHI0 - 2 * math.pi, math.pi / 4, places=12)

    def test_each_moment_is_a_line_across_the_dust_and_a_curve_out_through_the_throat(self):
        moments = load("embedding")["views"][0]["surfaces"]
        for vid, view in self.views.items():
            self.assertEqual(len(view["slices"]), len(moments))
            for mark, moment in zip(view["slices"], moments):
                eta = cycloid(2 * moment["time"] / AM)
                inside, outside = mark["lines"]
                self.assertEqual(inside, [[0, round(eta, 4)], [round(CHI0, 4), round(eta, 4)]])
                self.assertLess(math.dist(outside[0], [CHI0, eta]), 2e-4)
                xs = [p[0] for p in outside]
                self.assertTrue(all(b > a for a, b in zip(xs, xs[1:])), "the moment does not run outward")
                # It ends on the far sheet, beyond the bifurcation sphere at X = 3 chi_0 - pi.
                self.assertGreater(outside[-1][0], 3 * CHI0 - math.pi)
                self.assertTrue(all(abs(p[1]) < math.pi for p in outside))


class Diagrams(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.systems = load("diagrams")["systems"]

    def view(self, system):
        view, = self.systems[system]
        return view

    def marker(self, system, kind):
        return next(m for m in self.view(system)["markers"] if m["kind"] == kind)

    def test_the_event_horizon_is_the_ray_through_the_surface_at_r_s(self):
        # In eta and chi the ray is the straight line eta - chi = -5 pi/4, from the bang at
        # chi = pi/4 to the surface at eta = -pi/2.
        view = self.view("conformal")
        X0, X1, Y0, Y1 = view["box"]
        line, = self.marker("conformal", "event")["lines"]
        for u in line:
            chi, eta = X0 + u[0] * (X1 - X0), Y0 + u[1] * (Y1 - Y0)
            if abs(eta) < math.pi - 0.02:
                self.assertAlmostEqual(eta - chi, -5 * math.pi / 4, delta=3e-3)
        white, = self.marker("conformal", "past")["lines"]
        for u in white:
            chi, eta = X0 + u[0] * (X1 - X0), Y0 + u[1] * (Y1 - Y0)
            if abs(eta) < math.pi - 0.02:
                self.assertAlmostEqual(eta + chi, 5 * math.pi / 4, delta=3e-3)
        # In the proper time the same ray reaches the surface at c tau = -(pi/2 + 1) a_m/2.
        self.assertEqual(f"{(math.pi / 2 + 1) / 2:.2f}", "1.29")
        caption = " ".join(self.view("comoving")["caption"])
        self.assertIn("$c\\tau = -1.29\\,a_m$", caption)

    def test_the_isotropic_planes_surface_starts_where_its_caption_says(self):
        # The areal radius 2 r_s behind the throat is the isotropic radius (3/2 - sqrt 2)/2.
        r = (1.5 - math.sqrt(2)) / 2
        self.assertAlmostEqual(r * (1 + 1 / (4 * r)) ** 2, 2.0, places=12)
        self.assertEqual(f"{r:.3f}", "0.043")
        view = self.view("isotropic")
        self.assertIn("$r_d = 0.043\\,r_s$", " ".join(view["caption"]))
        X0, X1, Y0, Y1 = view["box"]
        surface, = self.marker("isotropic", "surface")["lines"]
        start = min(surface, key=lambda u: u[1])
        self.assertAlmostEqual(X0 + start[0] * (X1 - X0), r, delta=2e-3)
        # It stays behind the throat, r < r_s/4, and nears it as t grows.
        xs = [X0 + u[0] * (X1 - X0) for u in sorted(surface, key=lambda u: u[1])]
        self.assertTrue(all(x < 0.25 for x in xs))
        self.assertTrue(all(b >= a for a, b in zip(xs, xs[1:])))

    def test_the_cones_behind_the_throat_point_down_the_page(self):
        view = self.view("isotropic")
        X0, X1, _, _ = view["box"]
        behind = [c for c in view["cones"] if X0 + c["at"][0] * (X1 - X0) < 0.25]
        beyond = [c for c in view["cones"] if X0 + c["at"][0] * (X1 - X0) > 0.25]
        self.assertTrue(behind and beyond)
        self.assertTrue(all(c["a"][1] < 0 and c["b"][1] < 0 for c in behind))
        self.assertTrue(all(c["a"][1] > 0 and c["b"][1] > 0 for c in beyond))

    def test_the_dust_is_drawn_from_its_bang_to_its_crunch(self):
        # Both dust planes are hatched along their bottom and top edges, the two singularities.
        for system in ("comoving", "conformal"):
            view = self.view(system)
            edges = sorted(min(p[1] for p in ring) for ring in view["hatch"])
            self.assertEqual(len(view["hatch"]), 2, system)
            self.assertEqual(edges[0], 0, system)
            self.assertGreater(edges[1], 0.99, system)
        # The life of the dust in its proper time is pi a_m/c, the height of the comoving plane.
        X0, X1, Y0, Y1 = self.view("comoving")["box"]
        self.assertAlmostEqual(Y1 - Y0, math.pi, places=12)


@unittest.skipUnless(HAS_SYMPY, "sympy is not installed")
class Junction(unittest.TestCase):
    def test_the_surface_is_a_geodesic_of_the_exterior_and_no_shell_stands_on_it(self):
        import print_charts
        print_charts.semiclosed_world_matching()

    def test_the_conformal_chart_is_dust_at_rest(self):
        import sympy as sp
        import print_charts
        import verify_metrics as vm
        spec = print_charts.semiclosed_world("conformal")
        chart = print_charts.cp.Chart(spec["system"]["coords"], spec["system"]["parameters"], spec["chart_line_element"],
                                      spec["printer"], spec.get("pretty"))
        print_charts.semiclosed_world_check(chart, "conformal")
        eta = chart.symbols[0]
        am = chart.reader.parameters["a_m"]
        a = am * (1 + sp.cos(eta)) / 2
        G = chart.geo.raise_indices(chart.geo.einstein_ll(), 2, (0,))
        value = (vm._at(G, (0, 0)) + 3 * am / a ** 3).subs({eta: sp.Rational(7, 10), am: 3})
        self.assertLess(abs(complex(value.evalf(30))), 1e-20)


if __name__ == "__main__":
    unittest.main()
