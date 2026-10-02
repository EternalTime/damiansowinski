"""Bonnor's rotating dust cloud, the dipole of van Stockum's class of rigidly rotating dust:
python3 -m unittest discover -s _tools

What its texts state, held on the published files: in both charts the Einstein tensor is that of
dust at rest, with the density the Ricci scalar; the spherical chart is the cylindrical one
carried along rho = r sin(theta), z = r cos(theta); the dust outside a sphere of radius r has the
mass a^4/3r^3 in units of c^2/G and the Komar mass inside the sphere is minus that, so the total
is zero and the centre holds a negative mass without bound; the circles about the axis are null
on r^3 = a^2 rho and timelike inside it, while the signature holds through that surface; the
singularity is directional, the density diverging along the axis and dying away in the equatorial
plane; and the drawings hold what their captions say. The tests of the published mathematics need
sympy and are skipped where it is absent; the tests of the drawings read the files alone.
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


@unittest.skipUnless(HAS_SYMPY, "sympy is not installed")
class Published(unittest.TestCase):
    def setUp(self):
        import sympy
        import verify_metrics
        self.sp, self.vm = sympy, verify_metrics

    def read(self, system):
        """The reader of a published chart, its coordinates, and a function from a field's name to
        the matrix of its published components, with c = 1."""
        metric = json.loads((DATA / "metrics" / "bonnor_rotating_dust.json").read_text(encoding="utf-8"))
        entry = next(c for c in metric["coordinates"] if c["id"] == system)
        reader = self.vm.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]], ())
        names = entry["coords"]
        x = [reader.symbol[n] for n in names]

        def matrix(field, variant=None):
            block = entry[field] if variant is None else entry[field]["variants"][variant]["nonzero"]
            M = self.sp.zeros(len(x), len(x))
            for e in block:
                i, j = (names.index(k) for k in e["indices"])
                M[i, j] = reader(e["value"]).subs(reader.c, 1)
            return M

        def scalar(field):
            return reader(entry[field].partition("=")[2]).subs(reader.c, 1)
        return reader, x, matrix, scalar

    def number(self, expression, at):
        return float(self.sp.sympify(expression).subs(at).evalf(30))

    POINTS = ((0.7, 0.4), (1.3, -0.9), (0.2, 1.1), (2.5, 0.3))    # (rho, z) at a = 1

    def cylindrical_points(self, reader, x):
        a = reader.parameters["a"]
        return [{a: 1, x[1]: self.sp.Rational(str(rho)), x[3]: self.sp.Rational(str(z))} for rho, z in self.POINTS]

    def test_the_source_is_dust_at_rest_with_the_ricci_scalar_for_its_density(self):
        sp = self.sp
        for system in ("cylindrical", "spherical"):
            reader, x, matrix, scalar = self.read(system)
            a = reader.parameters["a"]
            G, R = matrix("einstein_tensor", "uu"), scalar("ricci_scalar")
            if system == "cylindrical":
                points = self.cylindrical_points(reader, x)
                density = lambda at: self.number(a ** 4 * (x[1] ** 2 + 4 * x[3] ** 2) * sp.exp(
                    -a ** 4 * x[1] ** 2 * (x[1] ** 2 - 8 * x[3] ** 2) / (8 * (x[1] ** 2 + x[3] ** 2) ** 4))
                    / (x[1] ** 2 + x[3] ** 2) ** 4, at)
            else:
                points = [{a: 1, x[1]: sp.Rational(13, 10), x[2]: sp.Rational(7, 10)},
                          {a: 1, x[1]: sp.Rational(1, 2), x[2]: sp.Rational(1, 4)}]
                density = lambda at: self.number(a ** 4 * (1 + 3 * sp.cos(x[2]) ** 2) * sp.exp(
                    -a ** 4 * sp.sin(x[2]) ** 2 * (9 * sp.sin(x[2]) ** 2 - 8) / (8 * x[1] ** 4)) / x[1] ** 6, at)
            for at in points:
                D = density(at)
                self.assertGreater(D, 0)
                self.assertAlmostEqual(self.number(R, at) / D, 1.0, places=12)
                for i in range(4):
                    for j in range(4):
                        want = D if i == j == 0 else 0.0
                        self.assertAlmostEqual(self.number(G[i, j], at) / D, want / D, places=12, msg=f"{system} G^{i}{j}")

    def test_the_spherical_chart_is_the_cylindrical_one_carried_along(self):
        sp = self.sp
        cyl, xc, mc, _ = self.read("cylindrical")
        sph, xs, ms, _ = self.read("spherical")
        t, r, th, ph = xs
        images = [t, r * sp.sin(th), ph, r * sp.cos(th)]
        J = sp.Matrix(4, 4, lambda i, j: sp.diff(images[i], xs[j]))
        g = mc("metric_components").subs({cyl.parameters["a"]: sph.parameters["a"]})
        pulled = J.T * g.subs(dict(zip(xc, images)), simultaneous=True) * J
        published = ms("metric_components")
        for values in ((1.3, 0.7), (0.5, 0.25), (3.0, 2.0)):
            at = {sph.parameters["a"]: 1, r: sp.Rational(str(values[0])), th: sp.Rational(str(values[1]))}
            for i in range(4):
                for j in range(4):
                    self.assertAlmostEqual(self.number(pulled[i, j], at), self.number(published[i, j], at), places=12)

    def test_the_total_mass_is_zero_and_every_sphere_holds_a_negative_mass(self):
        # In units of c^2/G: the dust outside the sphere of radius r_0 has the mass
        # (1/8 pi) int R sqrt(h) d^3x over its rest space, sqrt(h) = e^mu r^2 sin(theta), and the Komar
        # mass inside the sphere is Bratek, Jalocha and Kutschera's (1/8 pi) int K d_r K/sin(theta)
        # dtheta dphi with K = g_tphi. The two cancel at every r_0, and both go to zero far out.
        import mpmath as mp
        sp = self.sp
        reader, x, matrix, scalar = self.read("spherical")
        t, r, th, ph = x
        at = {reader.parameters["a"]: 1}
        g = matrix("metric_components").subs(at)
        rest = sp.sqrt(g[1, 1] * g[2, 2] * (g[3, 3] + g[0, 3] ** 2))
        integrand = sp.lambdify((r, th), scalar("ricci_scalar").subs(at) * rest, "mpmath")
        K = g[0, 3]
        flux = sp.lambdify((r, th), K * sp.diff(K, r) / sp.sin(th), "mpmath")
        last = None
        for r0 in (0.5, 1.0, 3.0):
            dust = mp.quad(lambda s, u: integrand(s, u), [r0, mp.inf], [0, mp.pi]) / 4
            komar = mp.quad(lambda u: flux(r0, u), [0, mp.pi]) / 4
            self.assertAlmostEqual(float(dust) * 3 * r0 ** 3, 1.0, places=8)
            self.assertAlmostEqual(float(komar + dust) * r0 ** 3, 0.0, places=8)
            self.assertLess(float(komar), 0)
            if last is not None:
                self.assertLess(float(dust), last)          # less dust outside a larger sphere
            last = float(dust)

    def test_far_away_the_field_is_the_frame_dragging_of_an_angular_momentum_with_no_mass(self):
        # g_tt = -1 everywhere, so no mass term, and g_tphi r/sin^2(theta) = a^2 = 2GJ/c^3, Kerr's
        # frame dragging term for the angular momentum J.
        reader, x, matrix, _ = self.read("spherical")
        g = matrix("metric_components")
        self.assertEqual(g[0, 0], -1)
        self.assertEqual(self.sp.simplify(g[0, 3] * x[1] / self.sp.sin(x[2]) ** 2 - reader.parameters["a"] ** 2), 0)

    def test_the_circles_turn_null_on_r_cubed_equal_a_squared_rho_and_the_signature_holds(self):
        sp = self.sp
        reader, x, matrix, _ = self.read("cylindrical")
        g = matrix("metric_components")
        for at, (rho, z) in zip(self.cylindrical_points(reader, x), self.POINTS):
            inside = (rho * rho + z * z) ** 1.5 < rho
            self.assertEqual(self.number(g[2, 2], at) < 0, inside, (rho, z))
            self.assertLess(self.number(g.det(), at), 0)
        # On the surface itself, here at z = 0 and rho = a, g_phiphi vanishes and the determinant does not.
        on = {reader.parameters["a"]: 1, x[1]: 1, x[3]: 0}
        self.assertAlmostEqual(self.number(g[2, 2], on), 0.0, places=14)
        self.assertAlmostEqual(self.number(g.det(), on), -math.exp(0.25), places=12)
        # Inside it the circle run toward -phi is future directed: g(-d_phi, d_t) = -g_tphi < 0.
        self.assertGreater(self.number(g[0, 2], {reader.parameters["a"]: 1, x[1]: sp.Rational(1, 2), x[3]: 0}), 0)

    def test_the_singularity_is_directional(self):
        sp = self.sp
        reader, x, _, scalar = self.read("cylindrical")
        a, rho, z = reader.parameters["a"], x[1], x[3]
        R, K = scalar("ricci_scalar"), scalar("kretschmann")
        # Along the axis the density grows as 4a^4/z^6 and the Kretschmann scalar is 4a^4(8a^4 - 27z^4)/z^12.
        for value in (sp.Rational(1, 2), sp.Rational(1, 10)):
            at = {a: 1, rho: 0, z: value}
            self.assertAlmostEqual(self.number(R * z ** 6, at), 4.0, places=12)
            self.assertAlmostEqual(self.number(K * z ** 12 / (4 * (8 - 27 * z ** 4)), at), 1.0, places=12)
        # In the equatorial plane both carry e^(-a^4/8 rho^4) or its square and die away toward the centre.
        values = [self.number(R, {a: 1, rho: v, z: 0}) for v in (sp.Rational(1, 2), sp.Rational(1, 4), sp.Rational(1, 8))]
        self.assertTrue(values[0] > values[1] > values[2] >= 0)
        self.assertLess(values[2], 1e-200)
        self.assertLess(abs(self.number(K, {a: 1, rho: sp.Rational(1, 8), z: 0})), 1e-200)
        self.assertEqual(sp.simplify(R.subs(z, 0) - a ** 4 * sp.exp(-a ** 4 / (8 * rho ** 4)) / rho ** 6), 0)


class Drawings(unittest.TestCase):
    """The published diagrams against the metric, from the numbers in the files, in units of a."""

    @classmethod
    def setUpClass(cls):
        data = json.loads((DATA / "diagrams" / "bonnor_rotating_dust.json").read_text(encoding="utf-8"))
        cls.views = {(system, v["id"]): v for system, views in data["systems"].items() for v in views}
        cls.figure = data["projections"]["cylindrical"][0]
        cls.embedding = json.loads((DATA / "embedding" / "bonnor_rotating_dust.json").read_text(encoding="utf-8"))

    @staticmethod
    def chart_points(view, line):
        """A ray's points in the view's own axes, from the unit square the file holds them in."""
        X0, X1, Y0, Y1 = view["box"]
        return [(X0 + u * (X1 - X0), Y0 + w * (Y1 - Y0)) for u, w in line]

    def end_to_end(self, view):
        found = {}
        for family, lines in view["rays"].items():
            for line in lines:
                (xa, ya), (xb, yb) = (self.chart_points(view, line)[k] for k in (0, -1))
                if abs(xb - xa) > 0.1 * (view["box"][1] - view["box"][0]):
                    found.setdefault(family, []).append((yb - ya) / (xb - xa))
        return found

    def test_the_rays_on_the_axis_run_at_45_degrees(self):
        for system in ("cylindrical", "spherical"):
            found = self.end_to_end(self.views[system, "axis"])
            self.assertEqual(len(found), 2)
            for slopes in found.values():
                for k in slopes:
                    self.assertAlmostEqual(abs(k), 1.0, delta=2e-3)

    def test_the_rays_of_each_circle_are_the_null_lines_of_the_metric(self):
        # On the circle at rho the metric is -c^2dt^2 + (2/rho) c dt dphi + (rho^2 - 1/rho^2) dphi^2,
        # null on c dt = (1/rho +- rho) dphi; the view draws rho phi across.
        for view, rho in (("inside", 0.5), ("outside", 1.5)):
            found = self.end_to_end(self.views["cylindrical", view])
            self.assertEqual(len(found), 2, view)
            want = {(1 / rho + rho) / rho, (1 / rho - rho) / rho}
            for family, slopes in found.items():
                nearest = min(want, key=lambda k: abs(k - slopes[0]))
                want.discard(nearest)
                for k in slopes:
                    self.assertAlmostEqual(k, nearest, delta=3e-3, msg=f"{view} {family}")
                lam = nearest * rho
                self.assertAlmostEqual(-lam * lam + 2 * lam / rho + rho * rho - rho ** -2, 0.0, places=12)
            self.assertEqual(want, set(), view)

    def test_inside_the_null_circle_both_families_climb_toward_plus_phi(self):
        inside = self.end_to_end(self.views["cylindrical", "inside"])
        outside = self.end_to_end(self.views["cylindrical", "outside"])
        self.assertTrue(all(k[0] > 0 for k in inside.values()))
        self.assertEqual(sorted(k[0] > 0 for k in outside.values()), [False, True])

    def test_the_rays_of_no_angular_momentum_keep_their_slope(self):
        # c dt/drho = +- e^(1/(16 rho^4)) sqrt(1 - 1/rho^4) in the plane z = 0, from chord to chord.
        for system in ("cylindrical", "spherical"):
            view = self.views[system, "equator"]
            checked = 0
            for lines in view["rays"].values():
                for line in lines:
                    pts = self.chart_points(view, line)
                    for (ra, ta), (rb, tb) in zip(pts, pts[1:]):
                        if min(ra, rb) < 1.5 or abs(rb - ra) < 0.02:
                            continue
                        mid = 0.5 * (ra + rb)
                        want = math.exp(1 / (16 * mid ** 4)) * math.sqrt(1 - mid ** -4)
                        self.assertAlmostEqual(abs((tb - ta) / (rb - ra)), want, delta=0.03)
                        checked += 1
            self.assertGreater(checked, 50, system)

    def test_the_figure_draws_the_null_circle_at_a_and_every_cone_null(self):
        turn = self.figure["turn"]
        critical = [line["points"] for line in turn["lines"] if line["class"] == "critical"]
        self.assertEqual(len(critical), 1)
        for X, Y, T in critical[0]:
            self.assertAlmostEqual(math.hypot(X, Y), 1.0, delta=1e-5)
            self.assertAlmostEqual(T, 0.0, delta=1e-9)
        self.assertEqual(len(turn["cones"]), 12)
        for cone in turn["cones"]:
            X, Y, T = cone["apex"]
            rho = math.hypot(X, Y)
            for x, y, t in cone["rim"]:
                dX, dY, dT = x - X, y - Y, t - T
                drho, dphi = (X * dX + Y * dY) / rho, (X * dY - Y * dX) / rho ** 2
                size = dX * dX + dY * dY + dT * dT
                null = (-(dT - dphi / rho) ** 2 + rho ** 2 * dphi ** 2 + math.exp(1 / (8 * rho ** 4)) * drho ** 2)
                self.assertAlmostEqual(null / size, 0.0, delta=2e-3)

    def test_inside_the_null_circle_the_cones_hold_the_circle_run_clockwise(self):
        # The cone meets the cylinder of its own rho in two generators, its edges there, and the
        # clockwise tangent of the circle, -d_phi, lies inside the cone when it is a combination of
        # the two with positive weights. At rho = a/2 the edges are (c dt, dphi) along (5/2, 1) and
        # (-3/2, -1), and -d_phi is 3/2 of the first plus 5/2 of the second.
        for cone in self.figure["turn"]["cones"]:
            X, Y, T = cone["apex"]
            rho = math.hypot(X, Y)
            sides = []
            for x, y, t in cone["rim"]:
                dX, dY, dT = x - X, y - Y, t - T
                sides.append((abs(X * dX + Y * dY) / rho, dT, (X * dY - Y * dX) / rho ** 2))
            flat = sorted(sides)[:4]                         # the generators nearest the cylinder of this rho
            (_, t1, p1), (_, t2, p2) = max(flat, key=lambda g: g[1]), min(flat, key=lambda g: g[1])
            det = t1 * p2 - t2 * p1
            alpha, beta = t2 / det, -t1 / det                # alpha e_1 + beta e_2 = (0, -1)
            if rho < 0.99:
                self.assertTrue(alpha > 0 and beta > 0, rho)
                self.assertAlmostEqual(t1 / p1, 2.5, delta=0.3)
                self.assertAlmostEqual(t2 / p2, 1.5, delta=0.3)
            elif rho > 1.01:
                self.assertFalse(alpha > 0 and beta > 0, rho)

    def test_the_embedded_moment_is_the_plane_outside_the_null_circle_in_minkowski_space(self):
        view, = self.embedding["views"]
        surface, = view["surfaces"]
        sheet, = surface["pieces"]
        self.assertEqual(sheet["space"], "minkowski")
        points = sheet["points"]
        self.assertAlmostEqual(points[0][0], 1.0002, places=12)
        self.assertAlmostEqual(points[-1][0], 5.0, places=12)
        for rho, R, Z in points:
            self.assertAlmostEqual(R, math.sqrt(rho * rho - rho ** -2), delta=1e-9)
        slopes = []
        for (xa, Ra, Za), (xb, Rb, Zb) in zip(points, points[1:]):
            self.assertGreater(Zb, Za)
            self.assertGreater(Rb - Ra, Zb - Za)             # every chord is spacelike
            slopes.append((Zb - Za) / (Rb - Ra))
            # dZ/dR = sqrt(1 - g_rhorho/(dR/drho)^2) with (dR/drho)^2 = (1 + x)^2/(1 - x), x = 1/rho^4.
            mid = 0.5 * (xa + xb)
            x = mid ** -4
            self.assertAlmostEqual(slopes[-1], math.sqrt(1 - math.exp(x / 8) * (1 - x) / (1 + x) ** 2), delta=2e-3)
        self.assertGreater(slopes[0], 0.999)                 # it leaves the axis along the light cone
        self.assertLess(slopes[-1], 0.1)                     # and flattens far out


if __name__ == "__main__":
    unittest.main()
