"""Lewis's stationary cylinders, the vacuum fields with cylindrical symmetry:
python3 -m unittest discover -s _tools

What the texts and drawings of `lewis` state, held on the published files. The parameters the
Lewis chart and the canonical chart are drawn at are held to being one spacetime by Costa,
Natario and Santos's map. The published components are held to the determinant -r^2 on the plane
of t and phi, to the radii the captions name, where Lewis's f and l and van Stockum's F and L
vanish, to the junction with the published dust of van Stockum at r = R, to the circular path of
light round the critical cylinder at every radius, to the Sagnac delay 4 pi C/c of the canonical
chart at every radius, and to the sign of the Lewis class's Kretschmann scalar on the two sides of
m = sqrt(3). The embedding diagram's surface is held, from the numbers written, to the circle's
radius sqrt(r^(3/2) - r^(1/2)), to leaving the axis along a light cone at the null circle and to
lying level where its two pieces meet. The tests of the published components need sympy and are
skipped where it is absent.
"""
import importlib.util
import json
import math
import sys
import unittest
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "MFS" / "assets" / "data"
HAS_SYMPY = importlib.util.find_spec("sympy") is not None
sys.path.insert(0, str(Path(__file__).resolve().parent / "derivations"))

# The values each chart is drawn at, as null_rays.py has them.
WEYL = {"n": Fraction(1, 2), "a": Fraction(4, 9), "b": Fraction(3, 2), "q": Fraction(1, 9), "ell": Fraction(1)}
CANONICAL = {"sigma": Fraction(1, 8), "j": Fraction(1, 8), "alpha": Fraction(1), "ell": Fraction(1)}
JOINED = math.exp(-0.5)             # ell = R/sqrt(e), at R = 1
DRAWN = {
    "lewis": {k: float(v) for k, v in WEYL.items()},
    "canonical": {k: float(v) for k, v in CANONICAL.items()},
    "lewis_class": {"m": 1.0, "a_1": 1.0, "b_1": 0.0, "a_2": 0.0, "ell": 1.0},
    "stockum_light": {"n": 0.5, "R": 1.0, "ell": JOINED},
    "stockum_critical": {"R": 1.0, "ell": JOINED},
    "stockum_heavy": {"m": 1.5, "R": 1.0, "ell": JOINED},
}


def published(metric_id, system):
    metric = json.loads((DATA / "metrics" / f"{metric_id}.json").read_text(encoding="utf-8"))
    return next(c for c in metric["coordinates"] if c["id"] == system)


class Parameters(unittest.TestCase):
    def test_the_two_charts_of_the_weyl_class_are_drawn_as_one_spacetime(self):
        """sigma = (1 - n)/4, alpha = (n ell - b q)^2/(n^2 ell^2 a) and j = b(n ell - b q)/(4 ell)."""
        n, a, b, q, ell = (WEYL[k] for k in ("n", "a", "b", "q", "ell"))
        self.assertEqual((1 - n) / 4, CANONICAL["sigma"])
        self.assertEqual((n * ell - b * q) ** 2 / (n ** 2 * ell ** 2 * a), CANONICAL["alpha"])
        self.assertEqual(b * (n * ell - b * q) / (4 * ell), CANONICAL["j"])
        # Lewis's coordinates turn at Omega = q/(n ell - b q) = 1/(3 ell), as the caption says.
        self.assertEqual(q / (n * ell - b * q), Fraction(1, 3))

    def test_the_generator_draws_at_these_values(self):
        import ast
        source = (ROOT / "_tools" / "derivations" / "null_rays.py").read_text(encoding="utf-8")
        for name, want in (("LEWIS_WEYL", WEYL), ("LEWIS_CANONICAL", CANONICAL)):
            line = next(row for row in source.splitlines() if row.startswith(name + " = "))
            drawn = ast.literal_eval(line.partition(" = ")[2])
            self.assertEqual({k: Fraction(str(v)) for k, v in drawn.items()}, want, name)


@unittest.skipUnless(HAS_SYMPY, "sympy is not installed")
class Components(unittest.TestCase):
    """The published metric, inverse and Christoffel symbols of each chart as numbers."""

    @classmethod
    def setUpClass(cls):
        import sympy as sp
        import verify_metrics as vm
        cls.sp = sp
        cls.charts = {}
        for system, values in DRAWN.items():
            entry = published("lewis", system)
            reader = vm.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]], ())
            at = {s: values[s.name] for s in reader.allowed if s.name in values}
            cls.charts[system] = (entry, reader, at)

    def value(self, system, text, r):
        entry, reader, at = self.charts[system]
        return float(reader(text).subs(at).subs(reader.symbol["r"], r))

    def block(self, system, r):
        """(g_tt, g_tphi, g_phiphi) at the radius r."""
        entry = self.charts[system][0]
        found = {tuple(e["indices"]): e["value"] for e in entry["metric_components"]}
        return tuple(self.value(system, found[pair], r) for pair in (("t", "t"), ("t", "\\phi"), ("\\phi", "\\phi")))

    def gamma(self, system, lower, r):
        entry = self.charts[system][0]
        found = {tuple(e["indices"]): e["value"] for e in entry["christoffel"]["variants"]["ull"]["nonzero"]}
        text = found.get(("r",) + lower)
        return self.value(system, text, r) if text else 0.0

    def zero(self, f, lo, hi):
        self.assertLess(f(lo) * f(hi), 0)
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if f(lo) * f(mid) > 0 else (lo, mid)
        return 0.5 * (lo + hi)

    def test_the_plane_of_t_and_phi_has_the_determinant_minus_r_squared(self):
        for system in DRAWN:
            for r in (1.0, 1.7, 3.0, 8.0):
                gtt, gtp, gpp = self.block(system, r)
                self.assertLess(abs(gtt * gpp - gtp ** 2 + r * r), 1e-9 * r * r, f"{system} at r = {r}")

    def test_the_radii_the_captions_name(self):
        gtt = lambda system: (lambda r: self.block(system, r)[0])       # noqa: E731
        gpp = lambda system: (lambda r: self.block(system, r)[2])       # noqa: E731
        # Lewis's chart: the null circle at ell, and d_t spacelike beyond 4 ell.
        self.assertLess(abs(self.zero(gpp("lewis"), 0.5, 2.0) - 1), 1e-9)
        self.assertLess(abs(self.zero(gtt("lewis"), 2.0, 6.0) - 4), 1e-9)
        # The canonical chart: the same null circle, and d_t timelike at every radius.
        self.assertLess(abs(self.zero(gpp("canonical"), 0.5, 2.0) - 1), 1e-9)
        self.assertTrue(all(gtt("canonical")(r) < 0 for r in (0.01, 0.5, 1, 4, 100, 1e4)))
        # The Lewis class at m = 1: f and l change sign together at r = e^(pi/2) ell.
        edge = math.exp(math.pi / 2)
        self.assertLess(abs(self.zero(gtt("lewis_class"), 2.0, 10.0) - edge), 1e-8)
        self.assertLess(abs(self.zero(gpp("lewis_class"), 2.0, 10.0) - edge), 1e-8)
        # van Stockum's light and critical cylinders: F = 0 at 3R and at eR, and L > 0 throughout.
        self.assertLess(abs(self.zero(gtt("stockum_light"), 1.0, 5.0) - 3), 1e-9)
        self.assertLess(abs(self.zero(gtt("stockum_critical"), 1.0, 5.0) - math.e), 1e-9)
        for system in ("stockum_light", "stockum_critical"):
            self.assertTrue(all(gpp(system)(r) > 0 for r in (1, 1.5, 3, 10, 100, 1e4)), system)
        # The heavy cylinder at m = 3/2: L < 0 from 1.14 R to 9.24 R, a band that repeats by e^(2 pi/m).
        first = self.zero(gpp("stockum_heavy"), 1.0, 1.5)
        last = self.zero(gpp("stockum_heavy"), 5.0, 12.0)
        self.assertEqual((round(first, 2), round(last, 2)), (1.14, 9.24))
        self.assertLess(abs(last / first - math.exp(math.pi / 1.5)), 1e-8)
        self.assertLess(gpp("stockum_heavy")(3.0), 0)
        self.assertGreater(gtt("stockum_heavy")(3.0), 0)

    def test_van_stockums_exteriors_join_the_published_dust_at_the_surface(self):
        """The metric and its first derivative along r at r = R = 1, with ell = R/sqrt(e), against
        the published chart of the dust, whose own R is 1/w."""
        import verify_metrics as vm
        sp = self.sp
        dust = published("stockum_dust", "cylindrical")
        reader = vm.Reader(dust["coords"], [p["symbol"] for p in dust["parameters"]], ())
        found = {tuple(e["indices"]): reader(e["value"]) for e in dust["metric_components"]}
        r = reader.symbol["r"]
        for system, w in (("stockum_light", math.sqrt(3) / 4), ("stockum_critical", 0.5),
                          ("stockum_heavy", math.sqrt(13) / 4)):
            entry, there, at = self.charts[system]
            ours = {tuple(e["indices"]): there(e["value"]).subs(at) for e in entry["metric_components"]}
            for pair in (("t", "t"), ("t", "\\phi"), ("r", "r"), ("\\phi", "\\phi"), ("z", "z")):
                inside = found[pair].subs(reader.parameters["R"], 1 / w)
                outside = ours[pair]
                for order in (0, 1):
                    left = float(sp.diff(inside, r, order).subs(r, 1))
                    right = float(sp.diff(outside, there.symbol["r"], order).subs(there.symbol["r"], 1))
                    self.assertLess(abs(left - right), 1e-9, f"{system}: g_{pair}, derivative {order}")

    def test_light_circles_the_critical_cylinder_at_every_radius(self):
        """c dt = (R/2) dphi is null, and Gamma^r_ab k^a k^b vanishes along it."""
        for r in (1.0, 2.0, 5.0, 40.0):
            gtt, gtp, gpp = self.block("stockum_critical", r)
            k = (0.5, 1.0)
            self.assertLess(abs(gtt * k[0] ** 2 + 2 * gtp * k[0] * k[1] + gpp * k[1] ** 2), 1e-9 * r)
            terms = [self.gamma("stockum_critical", ("t", "t"), r) * k[0] ** 2,
                     2 * self.gamma("stockum_critical", ("t", "\\phi"), r) * k[0] * k[1],
                     self.gamma("stockum_critical", ("\\phi", "\\phi"), r) * k[1] ** 2]
            self.assertLess(abs(sum(terms)), 1e-9 * sum(abs(x) for x in terms), f"r = {r}")

    def test_dust_at_rest_on_the_light_cylinders_surface_is_in_free_fall(self):
        self.assertLess(abs(self.gamma("stockum_light", ("t", "t"), 1.0)), 1e-12)
        self.assertGreater(abs(self.gamma("stockum_light", ("t", "t"), 2.0)), 1e-3)

    def test_the_sagnac_delay_is_the_same_at_every_radius(self):
        """In the canonical chart the laps of light round the axis take 2 pi (alpha r/u -+ C)/c, which
        differ by 4 pi C/c with C = 4j/(1 - 4 sigma) = ell, whatever the radius."""
        for r in (1.5, 3.0, 10.0, 200.0):
            gtt, gtp, gpp = self.block("canonical", r)
            root = math.sqrt(gtp ** 2 - gtt * gpp)
            slopes = sorted(((-gtp + root) / gtt, (-gtp - root) / gtt))
            # A lap toward +phi takes 2 pi slopes[1], and one toward -phi takes -2 pi slopes[0].
            self.assertLess(abs(-slopes[0] - slopes[1] - 2.0), 1e-9 * r, f"r = {r}")

    def test_the_lewis_classs_kretschmann_scalar_changes_sign_at_m_equal_root_three(self):
        entry, reader, at = self.charts["lewis_class"]
        m = next(s for s in reader.allowed if s.name == "m")
        K = reader(entry["kretschmann"].partition("=")[2])
        for value, sign in ((1.0, 1), (1.7, 1), (1.8, -1)):
            got = float(K.subs({**at, m: value}).subs(reader.symbol["r"], 2.0))
            self.assertGreater(sign * got, 0, f"m = {value}")


class Embedding(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        view, = json.loads((DATA / "embedding" / "lewis.json").read_text(encoding="utf-8"))["views"]
        surface, = view["surfaces"]
        cls.near, cls.far = surface["pieces"]

    def test_the_moment_leaves_the_axis_along_the_light_cone_and_lies_level_where_its_pieces_meet(self):
        """The moment t = 0 of the canonical chart at sigma = 1/8, j = ell/8 and alpha = 1: a piece in
        Minkowski space from the null circle r = ell, and a piece in flat space from the level circle
        r = 1.472 ell, the root of 16 s^17 - 9 s^16 - 16 s^9 + 6 s^8 - 1 in s = r^(1/8), out to 7 ell."""
        near, far = self.near, self.far
        self.assertEqual(near.get("space"), "minkowski")
        self.assertNotIn("space", far)
        level = near["points"][-1][0]
        s = level ** 0.125
        self.assertLess(abs(16 * s ** 17 - 9 * s ** 16 - 16 * s ** 9 + 6 * s ** 8 - 1), 1e-9)
        self.assertEqual(round(level, 3), 1.472)
        self.assertEqual(level, far["points"][0][0])
        self.assertLess(near["points"][0][0] - 1, 1e-4)
        self.assertEqual(far["points"][-1][0], 7.0)
        for piece in (near, far):
            for r, rho, z in piece["points"]:
                self.assertLess(abs(rho - math.sqrt(max(r ** 1.5 - math.sqrt(r), 0))), 2e-6, f"rho at r = {r}")
        P, Q = near["points"][:2]
        self.assertLess(abs((Q[2] - P[2]) / (Q[1] - P[1]) - 1), 1e-3, "along the light cone at the null circle")
        for (P, Q), where in ((near["points"][-2:], "the Minkowski piece"), (far["points"][:2], "the flat piece")):
            self.assertLess(abs((Q[2] - P[2]) / (Q[1] - P[1])), 0.1, f"{where} lies level at the join")

    def test_every_chord_is_the_metric_distance(self):
        """h dr^2 with h = r^(-3/8): the proper distance from r_0 is (16/13)(r^(13/16) - r_0^(13/16)),
        and a chord's length is taken with dX^2 - dZ^2 in Minkowski space and dX^2 + dZ^2 in flat."""
        proper = lambda r: 16 / 13 * r ** (13 / 16)       # noqa: E731
        for piece, sign in ((self.near, -1), (self.far, 1)):
            points = piece["points"]
            # Near the null circle the two terms of a chord cancel to the digits written, so the
            # sum along the piece is held.
            total = 0.0
            for (r0, x0, z0), (r1, x1, z1) in zip(points, points[1:]):
                total += math.sqrt(max((x1 - x0) ** 2 + sign * (z1 - z0) ** 2, 0.0))
            want = proper(points[-1][0]) - proper(points[0][0])
            self.assertLess(abs(total - want), 2e-3 * want, piece["id"])


if __name__ == "__main__":
    unittest.main()
