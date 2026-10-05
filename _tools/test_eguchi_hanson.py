"""The Eguchi-Hanson space, a Riemannian space of four dimensions with no time in it:
.venv.noindex/bin/python -m unittest discover -s _tools

What its texts and drawings state, held on the published files. The embedding diagram's three
surfaces are held to their closed forms from the numbers written and nothing else: the fibre over
a point of the bolt, rho = (r/2) sqrt(1 - a^4/r^4) rising at dz/dr = sqrt(3 + a^4/r^4)/2, smooth on
the bolt and a cone of half angle 30 degrees far away, whose circles are half as long as a
plane's; the bolt, a sphere of radius a/2; and the surface through the bolt's equator, rho = r/2
rising at dz/dr = sqrt((3r^4 + a^4)/(4(r^4 - a^4))). The published curvature is held to being Ricci
flat and self-dual up to orientation, to the Kretschmann scalar 384 a^8/r^12, and to the integral
of the Euler density, 3/2, which with the boundary's 1/2 is the Euler number 2. Each chart after
Eguchi and Hanson's is held to being theirs carried along its map. The tests of the curvature need
sympy and are skipped where it is absent.
"""
import importlib.util
import itertools
import json
import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "MFS" / "assets" / "data"
HAS_SYMPY = importlib.util.find_spec("sympy") is not None
sys.path.insert(0, str(Path(__file__).resolve().parent / "derivations"))


def published(system):
    metric = json.loads((DATA / "metrics" / "eguchi_hanson.json").read_text(encoding="utf-8"))
    return next(c for c in metric["coordinates"] if c["id"] == system)


class Embedding(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.views = {v["id"]: v for v in json.loads(
            (DATA / "embedding" / "eguchi_hanson.json").read_text(encoding="utf-8"))["views"]}

    def piece(self, view, piece_id):
        pieces = self.views[view]["surfaces"][0]["pieces"]
        return [tuple(p) for p in next(p for p in pieces if p["id"] == piece_id)["points"]]

    def test_the_fibre_closes_smoothly_on_the_bolt_and_opens_into_a_cone_of_thirty_degrees(self):
        fibre = self.piece("fibre", "fibre")
        self.assertEqual(fibre[0], (1.0, 0.0, 0.0))
        for r, rho, z in fibre:
            self.assertLess(abs(rho - 0.5 * r * math.sqrt(max(1 - r ** -4, 0))), 2e-6, f"the fibre's rho at {r}")
        for (r0, rho0, z0), (r1, rho1, z1) in zip(fibre, fibre[1:]):
            mid = 0.5 * (r0 + r1)
            self.assertLess(abs((z1 - z0) / (r1 - r0) - 0.5 * math.sqrt(3 + mid ** -4)), 2e-3,
                            f"the fibre's dz/dr at {mid}")
        # On the bolt the circle grows as fast as the distance from it, as a plane's does about a point.
        (r0, rho0, z0), (r1, rho1, z1) = fibre[0], fibre[1]
        self.assertLess(abs((rho1 - rho0) / math.hypot(rho1 - rho0, z1 - z0) - 1), 1e-2)
        # At the rim, r = 4a, it grows (1 + a^4/r^4)/2 as fast: a cone of half angle 30 degrees to 0.4%.
        (r0, rho0, z0), (r1, rho1, z1) = fibre[-2], fibre[-1]
        grows = (rho1 - rho0) / math.hypot(rho1 - rho0, z1 - z0)
        self.assertLess(abs(grows - (1 + 4.0 ** -4) / 2), 1e-3)
        self.assertLess(abs(math.degrees(math.asin(grows)) - 30), 0.3)

    def test_the_bolt_is_a_sphere_of_radius_half_a(self):
        bolt = self.piece("bolt", "bolt")
        self.assertEqual(bolt[0][0], 0.0)
        self.assertLess(abs(bolt[-1][0] - math.pi), 1e-9)
        for theta, rho, z in bolt:
            self.assertLess(abs(rho - 0.5 * math.sin(theta)), 2e-6, f"the bolt's rho at {theta}")
            self.assertLess(abs(z + 0.5 * math.cos(theta)), 2e-6, f"the bolt's z at {theta}")

    def test_the_surface_through_the_bolt_has_the_equator_for_its_smallest_circle(self):
        for sign, pid in ((1, "near"), (-1, "far")):
            half = self.piece("equator", pid)
            self.assertEqual(half[0], (1.0, 0.5, 0.0))
            self.assertEqual(min(rho for _, rho, _ in half), 0.5)
            for r, rho, z in half:
                self.assertLess(abs(rho - 0.5 * r), 2e-6, f"the half {pid}'s rho at {r}")
                self.assertGreaterEqual(sign * z, 0.0)
            for (r0, _, z0), (r1, _, z1) in zip(half[1:], half[2:]):
                mid = 0.5 * (r0 + r1)
                if mid < 1.05:
                    continue
                want = sign * math.sqrt((3 * mid ** 4 + 1) / (4 * (mid ** 4 - 1)))
                self.assertLess(abs((z1 - z0) / (r1 - r0) - want), 2e-3 * abs(want), f"the half {pid}'s dz/dr at {mid}")
        near, far = self.piece("equator", "near"), self.piece("equator", "far")
        for (r, rho, z), (r2, rho2, z2) in zip(near, far):
            self.assertEqual((r, rho, z), (r2, rho2, -z2))


@unittest.skipUnless(HAS_SYMPY, "sympy is not installed")
class Published(unittest.TestCase):
    def geometry(self, system):
        import verify_metrics as vm
        entry = published(system)
        reader = vm.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]], ())
        g = vm.metric_from_line_element(reader, entry["line_element"], entry["coords"])
        symbols = [reader.symbol[c] for c in entry["coords"]]
        return entry, reader, vm.Geometry(g, symbols, 600), symbols

    def test_the_curvature_is_ricci_flat_and_self_dual_up_to_orientation(self):
        import sympy as sp
        import verify_metrics as vm
        _, reader, geo, symbols = self.geometry("eguchi_hanson")
        ricci = geo.ricci_ll()
        self.assertTrue(all(vm._at(ricci, (i, j)) == 0 for i in range(4) for j in range(4)))
        riemann = geo.riemann_llll()
        r, theta = reader.symbol["r"], reader.symbol["\\theta"]
        root = r ** 3 * sp.sin(theta) / 8
        self.assertEqual(sp.simplify(geo.g.det() - root ** 2), 0)
        # *R_{abcd} = (1/2) sqrt(g) eps_{cdef} g^{ep} g^{fq} R_{abpq}, and R = s *R with one sign s in every slot.
        signs = set()
        for a, b, c, d in itertools.product(range(4), repeat=4):
            if a >= b or c >= d:
                continue
            dual = sum(sp.LeviCivita(c, d, e, f) * geo.ginv[e, p] * geo.ginv[f, q] * vm._at(riemann, (a, b, p, q))
                       for e, f, p, q in itertools.product(range(4), repeat=4)) * root / 2
            value = vm._at(riemann, (a, b, c, d))
            if vm.norm(value) == 0:
                self.assertEqual(vm.norm(dual), 0)
                continue
            signs.add(sp.simplify(dual / value))
        self.assertEqual(len(signs), 1)
        self.assertIn(signs.pop(), (1, -1))

    def test_the_kretschmann_scalar_integrates_to_the_euler_number_less_the_boundary_term(self):
        import sympy as sp
        _, reader, geo, symbols = self.geometry("eguchi_hanson")
        r, theta = reader.symbol["r"], reader.symbol["\\theta"]
        a = reader.parameters["a"]
        K = sp.simplify(geo.kretschmann())
        self.assertEqual(sp.simplify(K - 384 * a ** 8 / r ** 12), 0)
        self.assertEqual(K.subs(r, a), 384 / a ** 4)
        # For a Ricci flat metric the Euler density is K/(32 pi^2), and psi has the period 2 pi.
        volume = r ** 3 * sp.sin(theta) / 8
        self.assertEqual(sp.simplify(geo.g.det() - volume ** 2), 0)
        bulk = sp.integrate(sp.integrate(K * volume, (r, a, sp.oo)), (theta, 0, sp.pi)) * (2 * sp.pi) ** 2
        self.assertEqual(sp.simplify(bulk / (32 * sp.pi ** 2)), sp.Rational(3, 2))
        # The bolt, r = a, is a sphere of area pi a^2.
        element = a ** 2 * sp.sin(theta) / 4
        self.assertEqual(sp.simplify((geo.g[1, 1] * geo.g[2, 2]).subs(r, a) - element ** 2), 0)
        area = sp.integrate(element, (theta, 0, sp.pi)) * 2 * sp.pi
        self.assertEqual(sp.simplify(area - sp.pi * a ** 2), 0)

    def test_the_kahler_chart_is_eguchi_and_hansons_at_rho_to_the_fourth_r_to_the_fourth_less_a_to_the_fourth(self):
        import sympy as sp
        _, reader, geo, symbols = self.geometry("eguchi_hanson")
        _, other, there, coords = self.geometry("kahler")
        a, b = reader.parameters["a"], other.parameters["a"]
        rho = coords[0]
        image = [(rho ** 4 + b ** 4) ** sp.Rational(1, 4)] + coords[1:]
        J = sp.Matrix(4, 4, lambda i, j: sp.diff(image[i], coords[j]))
        pulled = J.T * geo.g.subs({**dict(zip(symbols, image)), a: b}, simultaneous=True) * J
        own = there.g.subs(other.held) if getattr(other, "held", None) else there.g
        point = {rho: sp.Rational(7, 5), coords[1]: sp.Rational(9, 10), b: sp.Rational(6, 5)}
        for i in range(4):
            for j in range(4):
                self.assertLess(abs((pulled[i, j] - own[i, j]).subs(point).evalf(30)), 1e-25)

    def test_the_two_centre_chart_is_eguchi_and_hansons_with_the_bolt_between_the_centres(self):
        import sympy as sp
        _, reader, geo, symbols = self.geometry("eguchi_hanson")
        _, other, there, coords = self.geometry("two_centre")
        a, b = reader.parameters["a"], other.parameters["a"]
        r, theta, phi, psi = symbols
        image = [sp.sqrt(r ** 4 - a ** 4) * sp.sin(theta) / a, r ** 2 * sp.cos(theta) / a, psi, a * phi / 4]
        J = sp.Matrix(4, 4, lambda i, j: sp.diff(image[i], symbols[j]))
        own = there.g.subs(other.held) if getattr(other, "held", None) else there.g
        pulled = J.T * own.subs({**dict(zip(coords, image)), b: a}, simultaneous=True) * J
        for at in ({r: sp.Rational(3, 2), theta: sp.Rational(7, 10), a: 1},
                   {r: sp.Rational(5, 2), theta: sp.Rational(21, 10), a: sp.Rational(6, 5)}):
            for i in range(4):
                for j in range(4):
                    self.assertLess(abs((pulled[i, j] - geo.g[i, j]).subs(at).evalf(30)), 1e-25)
        # The bolt r = a lands on the axis between the centres, rho = 0 and |z| <= a, its poles on the centres.
        self.assertEqual(image[0].subs(r, a), 0)
        self.assertEqual([image[1].subs({r: a, theta: t}) for t in (0, sp.pi)], [a, -a])
        # The period of tau is that of phi, 2 pi, times a/4.
        self.assertIn("\\tau \\in [0, \\pi a/2)", published("two_centre")["domains"])


if __name__ == "__main__":
    unittest.main()
