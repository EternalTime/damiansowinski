"""Bonnor's stars of charged dust, held on the published files:
python3 -m unittest discover -s _tools

Every chart is -U^{-2}c^2dt^2 + U^2 times flat space. The published potentials have the densities
their sources state, 4 pi G rho/c^2 = -(flat Laplacian of U)/U^3: Bonnor's (3.11) of 1965, Bonnor
and Wickramasuriya's (3.3) and (I.1) of 1975 and Lemos and Weinberg's (3.2), each positive. Each
star's U and its slope meet the exterior's at the surface, the sphere of 1975 has the central
redshift 3m/(2r_0) and keeps the radius 4m/3 and the central density 2/(9 pi m^2) as r_0 -> 0, the
spheroid has the redshift (4.8) on its disc, the cloud is the single hole at b = 0, and the two
exterior charts are one field, the published single hole of majumdar_papapetrou and the published
rn_metric at r_s = 2m, r_q = m. Those need sympy and are skipped where it is absent. The drawings
are held from their numbers alone: the embedded surfaces to r U and to their joins, the spheroid's
disc to being flat, and each ray of the spacetime diagrams to what ct keeps along it.
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
ID = "bonnor_charged_dust"


def read(folder):
    return json.loads((DATA / folder / f"{ID}.json").read_text(encoding="utf-8"))


@unittest.skipUnless(HAS_SYMPY, "sympy is not installed")
class ThePotentials(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import sympy
        import verify_metrics
        cls.sp, cls.vm = sympy, verify_metrics
        cls.charts = {}
        for entry in read("metrics")["coordinates"]:
            reader = cls.vm.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]], ())
            cls.charts[entry["id"]] = (reader, cls.vm.metric_from_line_element(reader, entry["line_element"],
                                                                               entry["coords"]))

    def potential(self, chart):
        """(reader, U, the radial coordinate) of a chart that names its potential, U written out."""
        reader, g = self.charts[chart]
        U = reader.parameters["U"]
        x = reader.symbol[reader.coords[1]]
        at = {s: self.sp.Rational(7, 10) for s in U.free_symbols}
        self.close(g[0, 0] * U ** 2 / reader.c ** 2 + 1, at)
        return reader, U, x

    def density(self, U, r):
        """4 pi G rho/c^2 of a spherical potential: -(U'' + 2U'/r)/U^3."""
        sp = self.sp
        return -(sp.diff(U, r, 2) + 2 * sp.diff(U, r) / r) / U ** 3

    def close(self, e, values, places=25):
        self.assertAlmostEqual(float(self.sp.N(self.sp.sympify(e).subs(values), 40)), 0, places=places)

    def test_every_chart_is_the_same_line_element_in_its_potential(self):
        sp = self.sp
        for chart in ("sphere_1965", "sphere_1975", "quasi_black_hole"):
            reader, g = self.charts[chart]
            U, r = reader.parameters["U"], reader.symbol["r"]
            self.assertEqual(sp.simplify(g[1, 1] - U ** 2), 0, chart)
            self.assertEqual(sp.simplify(g[2, 2] - U ** 2 * r ** 2), 0, chart)

    def test_bonnors_sphere_of_1965_has_his_density_and_meets_the_exterior(self):
        sp = self.sp
        reader, U, r = self.potential("sphere_1965")
        m, r0 = reader.parameters["m"], reader.parameters["r_0"]
        for at in ({m: sp.Rational(3, 2), r0: sp.Rational(7, 10), r: sp.Rational(2, 5)},
                   {m: sp.Rational(1, 3), r0: 2, r: sp.Rational(9, 5)}):
            self.close(self.density(U, r) - 3 * m / r0 ** 3 / (1 + m / r0) ** 3 / (1 + m * r ** 2 / r0 ** 3), at)
            self.assertGreater(sp.N(self.density(U, r).subs(at)), 0)
            self.assertLess(sp.N(sp.diff(self.density(U, r), r).subs(at)), 0)
            edge = {k: v for k, v in at.items() if k != r}
            self.close(U.subs(r, r0) - 1 - m / r0, edge)
            self.close(sp.diff(U, r).subs(r, r0) + m / r0 ** 2, edge)
            self.close(U.subs(r, 0) - (1 + m / r0) ** sp.Rational(3, 2), edge)

    def test_the_sphere_of_1975_has_its_density_its_redshift_and_its_limits(self):
        sp = self.sp
        reader, U, r = self.potential("sphere_1975")
        m, r0 = reader.parameters["m"], reader.parameters["r_0"]
        self.assertEqual(sp.simplify(self.density(U, r) - 3 * m / (r0 ** 3 * U ** 3)), 0)
        self.assertEqual(sp.simplify(U.subs(r, 0) - 1 - 3 * m / (2 * r0)), 0)
        self.assertEqual(sp.simplify(U.subs(r, r0) - 1 - m / r0), 0)
        self.assertEqual(sp.simplify(sp.diff(U, r).subs(r, r0) + m / r0 ** 2), 0)
        at = {m: 1, r0: sp.Rational(1, 2), r: sp.Rational(1, 4)}
        self.assertGreater(sp.N(sp.diff(self.density(U, r), r).subs(at)), 0)
        # As r_0 -> 0: the proper radius 4m/3, the area 4 pi m^2 and the central density 2/(9 pi m^2).
        self.assertEqual(sp.limit(sp.integrate(U, (r, 0, r0)), r0, 0, "+"), 4 * m / 3)
        self.assertEqual(sp.limit((r0 * U.subs(r, r0)) ** 2, r0, 0, "+"), m ** 2)
        self.assertEqual(sp.limit(self.density(U, r).subs(r, 0) / (4 * sp.pi), r0, 0, "+"), 2 / (9 * sp.pi * m ** 2))

    def test_lemos_and_weinbergs_cloud_has_their_density_and_is_the_hole_at_b_zero(self):
        sp = self.sp
        reader, U, r = self.potential("quasi_black_hole")
        m, b = reader.parameters["m"], reader.parameters["b"]
        at = {m: sp.Rational(6, 5), b: sp.Rational(3, 10), r: sp.Rational(7, 10)}
        root = sp.sqrt(r ** 2 + b ** 2)
        self.close(self.density(U, r) - 3 * m * b ** 2 / ((r ** 2 + b ** 2) * (m + root) ** 3), at)
        self.close(U.subs(b, 0) - 1 - m / r, at)
        self.close(U.subs(r, 0) - 1 - m / b, at)

    def test_the_spheroid_is_matched_at_its_surface_and_has_its_redshift_on_the_disc(self):
        sp = self.sp
        reader, inside, u = self.potential("spheroid_interior")
        other, outside, w = self.potential("spheroid_exterior")
        m, a, u0 = (reader.parameters[k] for k in ("m", "a", "u_0"))
        outside = outside.subs({other.parameters["m"]: m, other.parameters["a"]: a, w: u})
        at = {m: sp.Rational(5, 4), a: sp.Rational(3, 4), u0: sp.Rational(9, 10)}
        self.close((inside - outside).subs(u, u0), at)
        self.close((sp.diff(inside, u) - sp.diff(outside, u)).subs(u, u0), at)
        alpha = sp.atan(1 / sp.sinh(u0))
        self.close(inside.subs(u, 0) - 1 - m / a * (alpha + u0 / (4 * sp.cosh(u0))), at)
        # Outside, U is harmonic: (cosh u U')' = 0. Inside, the density is their (I.1) on the axis.
        self.close(sp.diff(sp.cosh(u) * sp.diff(outside, u), u), {**at, u: sp.Rational(7, 5)})
        laplace = sp.diff(sp.cosh(u) * sp.diff(inside, u), u) / (a ** 2 * sp.cosh(u) * (sp.sinh(u) ** 2 + 1))
        wanted = m * u ** 2 * (3 + u * sp.tanh(u)) / (sp.cosh(u0) * a ** 3 * u0 ** 3 * (sp.sinh(u) ** 2 + 1) * inside ** 3)
        self.close(-laplace / inside ** 3 - wanted, {**at, u: sp.Rational(2, 5)})
        # Flattened to a disc, u_0 -> 0, the redshift on it tends to pi m/(2a).
        flat = {**at, u0: sp.Rational(1, 10 ** 9)}
        self.close(inside.subs(u, 0) - 1 - sp.pi * m / (2 * a), flat, places=7)

    def test_the_two_exterior_charts_are_one_field_and_the_fields_of_their_neighbours(self):
        sp = self.sp
        reader, g = self.charts["exterior"]
        other, h = self.charts["exterior_areal"]
        r, m, R = reader.symbol["r"], reader.parameters["m"], other.symbol["R"]
        same = {other.parameters["m"]: m, other.symbol["\\theta"]: reader.symbol["\\theta"], other.c: reader.c}
        pulled = h.subs(same).subs(R, r + m)
        for i in range(4):
            self.assertEqual(sp.simplify(pulled[i, i] - g[i, i]), 0, i)
        for metric_id, chart, mine, values in (
                ("majumdar_papapetrou", "isotropic", g, lambda q: {q.parameters["m"]: m, q.symbol["r"]: r}),
                ("rn_metric", "spherical", h.subs(same).subs(reader.c, 1),
                 lambda q: {q.parameters["r_s"]: 2 * m, q.parameters["r_q"]: m, q.symbol["r"]: R})):
            there = next(c for c in json.loads((DATA / "metrics" / f"{metric_id}.json").read_text(encoding="utf-8"))
                         ["coordinates"] if c["id"] == chart)
            q = self.vm.Reader(there["coords"], [p["symbol"] for p in there["parameters"]], ())
            line = self.vm.metric_from_line_element(q, there["line_element"], there["coords"])
            at = {**values(q), q.symbol["\\theta"]: reader.symbol["\\theta"], q.c: reader.c}
            for i in range(4):
                self.assertEqual(sp.simplify(line[i, i].subs(at) - mine[i, i]), 0, f"{metric_id} {i}")


class TheDrawings(unittest.TestCase):
    """From the numbers in the files alone, at m = 1."""

    @staticmethod
    def sphere_1975(r):
        return 1.75 - r * r / 16

    @staticmethod
    def sphere_1965(r):
        return 27 ** 0.5 / math.sqrt(8 + r * r)

    @staticmethod
    def cloud(r):
        return 1 + 1 / math.sqrt(r * r + 0.25)

    def pieces(self, view_id):
        view = next(v for v in read("embedding")["views"] if v["id"] == view_id)
        return {p["id"]: p["points"] for p in view["surfaces"][0]["pieces"]}

    def test_each_embedded_circle_has_the_radius_r_times_its_potential(self):
        for view_id, piece, U in (("star", "star", self.sphere_1975), ("star_1965", "star", self.sphere_1965),
                                  ("cloud", "cloud", self.cloud), ("star", "exterior", lambda r: 1 + 1 / r)):
            for r, rho, _ in self.pieces(view_id)[piece]:
                self.assertAlmostEqual(rho, r * U(r), places=5, msg=f"{view_id} {piece} at {r}")

    def test_each_sphere_meets_its_exterior_in_one_circle_of_radius_three(self):
        for view_id in ("star", "star_1965"):
            pieces = self.pieces(view_id)
            self.assertEqual(pieces["star"][-1][0], 2.0)
            self.assertEqual(pieces["star"][-1][1:], pieces["exterior"][0][1:], view_id)
            self.assertAlmostEqual(pieces["exterior"][0][1], 3.0, places=6)

    def test_the_sphere_of_1975_has_the_proper_radius_of_its_potential(self):
        # The length of the profile from the centre to the surface is the integral of U, 7/2 - 1/6.
        points = self.pieces("star")["star"]
        length = sum(math.hypot(b[1] - a[1], b[2] - a[2]) for a, b in zip(points, points[1:]))
        self.assertAlmostEqual(length, 3.5 - 1 / 6, places=3)

    def test_the_spheroids_disc_is_flat_and_meets_the_body_at_the_focal_ring(self):
        pieces = self.pieces("spheroid")
        on_disc = 1 + math.atan(1 / math.sinh(1)) + 1 / (4 * math.cosh(1))
        heights = {round(z, 6) for _, _, z in pieces["disc"]}
        self.assertEqual(len(heights), 1)
        for theta, rho, _ in pieces["disc"]:
            self.assertAlmostEqual(rho, on_disc * math.cos(theta), places=5)
        self.assertAlmostEqual(pieces["body"][0][1], on_disc, places=5)
        self.assertAlmostEqual(pieces["body"][0][2], pieces["disc"][0][2], places=6)
        for mine, theirs in zip(pieces["body"][-1][1:], pieces["exterior"][0][1:]):
            self.assertAlmostEqual(mine, theirs, places=6)
        surface = (1 + math.atan(1 / math.sinh(1))) * math.cosh(1)
        self.assertAlmostEqual(pieces["exterior"][0][1], surface, places=5)

    def test_every_ray_of_the_spherical_planes_keeps_ct_plus_or_minus_the_integral_of_the_potential_squared(self):
        A, B = 1.75, 1 / 16
        keeps = {
            ("sphere_1975", "radial"): lambda r: A * A * r - 2 * A * B * r ** 3 / 3 + B * B * r ** 5 / 5,
            ("sphere_1965", "radial"): lambda r: 27 / math.sqrt(8) * math.atan(r / math.sqrt(8)),
            ("exterior", "radial"): lambda r: r + 2 * math.log(r) - 1 / r,
            ("exterior_areal", "radial"): lambda R: R + 2 * math.log(R - 1) - 1 / (R - 1),
            ("quasi_black_hole", "radial"): lambda r: r + 2 * math.asinh(2 * r) + 2 * math.atan(2 * r),
        }
        systems = read("diagrams")["systems"]
        checked = 0
        for (system, view_id), F in keeps.items():
            view = next(v for v in systems[system] if v["id"] == view_id)
            X0, X1, Y0, Y1 = view["box"]
            for family in view["rays"].values():
                for ray in family:
                    points = [(X0 + u * (X1 - X0), Y0 + w * (Y1 - Y0)) for u, w in ray]
                    inside = [(x, t) for x, t in points if X0 + 1e-6 < x <= X1]
                    if len(inside) < 3:
                        continue
                    # The file rounds a point to 1e-4 of the box, and F is steepest where U is greatest.
                    spreads = [max(k) - min(k) for k in ([t + sign * F(x) for x, t in inside] for sign in (1, -1))]
                    self.assertLess(min(spreads), 5e-3 * (Y1 - Y0), f"{system}/{view_id}")
                    checked += 1
        self.assertGreater(checked, 40)


if __name__ == "__main__":
    unittest.main()
