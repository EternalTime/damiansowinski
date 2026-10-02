"""The gravitational instantons of 1977 and 1978, Riemannian spaces of four dimensions with no time:
python3 -m unittest discover -s _tools

What the texts and drawings state, held on the published files. The embedding diagram's surfaces
are held to their closed forms from the numbers written and nothing else: the cigar of the Euclidean
Schwarzschild solution, rho = 2 r_s sqrt(1 - r_s/r), smooth on the bolt and a cylinder of radius
2 r_s far away, the same surface in the regular chart, where rho = x, and Flamm's paraboloid through
the bolt; the nut's surface and Taub-bolt's fibre, cylinders of radius 4n far away; the axis of two
centres, a sphere between them and a cigar beyond each; the complex line of CP^2, a round sphere of
radius L/2; and the fibre of Page's space, closed and smooth on both bolts. The published curvature
is held to what each chart claims: Ricci flat or Einstein, self dual or not, the Kretschmann scalar,
the period that the surface gravity of each nut and bolt asks for, and the Euler number, the integral
of the Euler density, which counts one for a nut and two for a bolt. Page's chart is held to
regularity exactly on the quartic and to Gibbons and Hawking's form there. The tests of the
curvature need sympy and are skipped where it is absent.
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

NU = 0.281701557908774005923342651117


def published(system):
    metric = json.loads((DATA / "metrics" / "gravitational_instantons.json").read_text(encoding="utf-8"))
    return next(c for c in metric["coordinates"] if c["id"] == system)


class Embedding(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.views = {v["id"]: v for v in json.loads(
            (DATA / "embedding" / "gravitational_instantons.json").read_text(encoding="utf-8"))["views"]}

    def piece(self, view, piece_id):
        pieces = self.views[view]["surfaces"][0]["pieces"]
        return [tuple(p) for p in next(p for p in pieces if p["id"] == piece_id)["points"]]

    def rate(self, a, b):
        """How fast the circle's radius grows with distance along the profile between two points."""
        return (b[1] - a[1]) / math.hypot(b[1] - a[1], b[2] - a[2])

    def test_each_view_belongs_to_its_chart(self):
        self.assertEqual({k: v["system"] for k, v in self.views.items()},
                         {"cigar": "schwarzschild", "bridge": "schwarzschild", "cigar_x": "regular", "nut": "taub_nut",
                          "bolt_fibre": "taub_bolt", "centres": "multi_centre", "line": "cp2_distance",
                          "line_r": "cp2", "page_fibre": "page"})

    def test_the_cigar_closes_smoothly_on_the_bolt_and_widens_to_a_cylinder_of_radius_two_r_s(self):
        cigar = self.piece("cigar", "cigar")
        self.assertEqual(cigar[0], (1.0, 0.0, 0.0))
        for r, rho, z in cigar:
            self.assertLess(abs(rho - 2 * math.sqrt(max(1 - 1 / r, 0))), 2e-6, f"the cigar's rho at {r}")
            self.assertLess(rho, 2.0)
        for (r0, _, z0), (r1, _, z1) in zip(cigar[1:], cigar[2:]):
            mid = 0.5 * (r0 + r1)
            want = math.sqrt((mid + 1) * (mid * mid + 1) / mid ** 3)
            self.assertLess(abs((z1 - z0) / (r1 - r0) - want), 3e-3 * want, f"the cigar's dz/dr at {mid}")
        # On the bolt the circle grows as fast as the distance from it, and at r it grows r_s^2/r^2 as fast.
        self.assertLess(abs(self.rate(cigar[0], cigar[1]) - 1), 1e-2)
        self.assertLess(abs(self.rate(cigar[-2], cigar[-1]) - 1 / (0.5 * (cigar[-2][0] + cigar[-1][0])) ** 2), 1e-3)

    def test_the_regular_chart_draws_the_same_cigar_with_the_radius_x(self):
        cigar, plane = self.piece("cigar", "cigar"), self.piece("cigar_x", "cigar")
        self.assertEqual(plane[0], (0.0, 0.0, 0.0))
        for x, rho, z in plane:
            self.assertLess(abs(rho - x), 2e-6)
        # The same surface: at the areal radius r = 4 r_s^3/(4 r_s^2 - x^2) of a marked circle the heights agree.
        heights = {round(r, 9): z for r, _, z in cigar}
        for x, _, z in plane:
            r = round(4 / (4 - x * x), 9)
            if r in heights:
                self.assertLess(abs(z - heights[r]), 2e-4)
        self.assertLess(abs(4 / (4 - plane[-1][0] ** 2) - 256 / 31), 1e-12)

    def test_the_surface_through_the_bolt_is_flamms_paraboloid_on_both_sides(self):
        for sign, pid in ((1, "near"), (-1, "far")):
            half = self.piece("bridge", pid)
            self.assertEqual(half[0], (1.0, 1.0, 0.0))
            for r, rho, z in half:
                self.assertLess(abs(rho - r), 2e-6)
                self.assertLess(abs(z - sign * 2 * math.sqrt(max(r - 1, 0))), 5e-6, f"Flamm's height at {r}")

    def test_the_nut_and_the_bolt_of_taub_nut_close_smoothly_and_widen_to_cylinders_of_radius_four_n(self):
        nut = self.piece("nut", "nut")
        self.assertEqual(nut[0], (1.0, 0.0, 0.0))
        for r, rho, z in nut:
            self.assertLess(abs(rho - 4 * math.sqrt(max((r - 1) / (r + 1), 0))), 2e-6)
        self.assertLess(abs(self.rate(nut[0], nut[1]) - 1), 1e-2)
        fibre = self.piece("bolt_fibre", "fibre")
        self.assertEqual(fibre[0], (2.0, 0.0, 0.0))
        for r, rho, z in fibre:
            self.assertLess(abs(rho - 4 * math.sqrt(max((r - 2) * (2 * r - 1) / (2 * (r * r - 1)), 0))), 2e-6)
            self.assertLess(rho, 4.0)
        self.assertLess(abs(self.rate(fibre[0], fibre[1]) - 1), 1e-2)
        (r0, rho0, z0), (r1, rho1, z1) = fibre[-2], fibre[-1]
        mid = 0.5 * (r0 + r1)
        self.assertLess(abs(self.rate(fibre[-2], fibre[-1]) - (5 * mid * mid - 8 * mid + 5) / (mid * mid - 1) ** 2), 1e-3)

    def test_the_axis_of_two_centres_is_a_sphere_between_them_and_a_cigar_beyond_each(self):
        below, between, above = (self.piece("centres", k) for k in ("below", "between", "above"))
        potential = {"below": lambda z: 1 + 2 / (-2 - z) + 2 / (2 - z), "between": lambda z: 1 + 2 / (z + 2) + 2 / (2 - z),
                     "above": lambda z: 1 + 2 / (z + 2) + 2 / (z - 2)}
        for name, piece in (("below", below), ("between", between), ("above", above)):
            for z, rho, height in piece:
                want = 0.0 if abs(abs(z) - 2) < 1e-12 else 4 / math.sqrt(potential[name](z))
                self.assertLess(abs(rho - want), 2e-6, f"{name}: rho at z = {z}")
        # The three parts touch at the two centres, where each circle has closed to a point.
        self.assertEqual((below[-1][1], between[0][1], between[-1][1], above[0][1]), (0.0, 0.0, 0.0, 0.0))
        self.assertEqual(below[-1][2], between[0][2])
        self.assertEqual(between[-1][2], above[0][2])
        widest = max(between, key=lambda p: p[1])
        self.assertEqual(widest[0], 0.0)
        self.assertLess(abs(widest[1] - 4 / math.sqrt(3)), 2e-6)
        for a, b in ((between[0], between[1]), (above[0], above[1])):
            self.assertLess(abs(self.rate(a, b) - 1), 1e-2)

    def test_the_complex_line_is_a_round_sphere_of_radius_half_l(self):
        for view in ("line", "line_r"):
            line = self.piece(view, "line")
            self.assertEqual(line[0][:2], (0.0, 0.0))
            self.assertLess(abs(line[-1][0] - math.pi / 2), 1e-12)
            for chi, rho, z in line:
                self.assertLess(abs(rho - 0.5 * math.sin(2 * chi)), 2e-6)
                self.assertLess(abs(z + 0.5 * math.cos(2 * chi)), 2e-6)
        # The circles marked for the chart in r are at r = L tan(chi/L) = 1/2, 1 and 2 L.
        marked = sorted(math.tan(ring["x"]) for ring in self.views["line_r"]["surfaces"][0]["rings"])
        for got, want in zip(marked, (0.5, 1.0, 2.0)):
            self.assertLess(abs(got - want), 1e-12)

    def test_the_fibre_of_pages_space_is_closed_and_smooth_on_both_bolts(self):
        fibre = self.piece("page_fibre", "fibre")
        n = 3 + 6 * NU ** 2 - NU ** 4
        self.assertEqual((fibre[0][0], fibre[0][1]), (0.0, 0.0))
        self.assertLess(abs(fibre[-1][0] - math.pi), 1e-12)
        self.assertLess(fibre[-1][1], 1e-6)
        for chi, rho, z in fibre:
            p = 3 - NU ** 2 - NU ** 2 * (1 + NU ** 2) * math.cos(chi) ** 2
            q = 1 - NU ** 2 * math.cos(chi) ** 2
            self.assertLess(abs(rho - 4 * NU * math.sqrt(3 * (1 + NU ** 2) * p / q) * math.sin(chi) / n), 2e-6)
        self.assertLess(abs(self.rate(fibre[0], fibre[1]) - 1), 1e-2)
        self.assertLess(abs(self.rate(fibre[-1], fibre[-2]) - 1), 1e-2)
        # The numbers its caption states: 3.27 long from bolt to bolt, 1.00 at its widest.
        length = sum(math.hypot(b[1] - a[1], b[2] - a[2]) for a, b in zip(fibre, fibre[1:]))
        self.assertEqual(f"{length:.2f}", "3.27")
        self.assertEqual(f"{max(rho for _, rho, _ in fibre):.2f}", "1.00")
        self.assertEqual(f"{12 * math.pi * (1 - NU ** 4) / n:.2f}", "10.80")
        self.assertIn("$3.27/\\sqrt{\\Lambda}$ long", self.views["page_fibre"]["caption"][1])
        self.assertIn("area $10.80/\\Lambda$", self.views["page_fibre"]["caption"][1])


@unittest.skipUnless(HAS_SYMPY, "sympy is not installed")
class Published(unittest.TestCase):
    def geometry(self, system):
        import verify_metrics as vm
        entry = published(system)
        reader = vm.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]], ())
        g = vm.metric_from_line_element(reader, entry["line_element"], entry["coords"])
        symbols = [reader.symbol[c] for c in entry["coords"]]
        return entry, reader, vm.Geometry(g, symbols, 600), symbols

    def duality(self, geo, tensor, root):
        """The set of ratios *T/T over the independent slots of a tensor with Riemann's symmetries:
        {1} or {-1} where it is self dual for one orientation, more than one value where it is neither."""
        import sympy as sp
        import verify_metrics as vm
        ratios = set()
        for a, b, c, d in itertools.product(range(4), repeat=4):
            if a >= b or c >= d:
                continue
            dual = sum(sp.LeviCivita(c, d, e, f) * geo.ginv[e, p] * geo.ginv[f, q] * vm._at(tensor, (a, b, p, q))
                       for e, f, p, q in itertools.product(range(4), repeat=4)) * root / 2
            value = vm._at(tensor, (a, b, c, d))
            if vm.norm(value) == 0 and vm.norm(dual) == 0:
                continue
            ratios.add(sp.zoo if vm.norm(value) == 0 else sp.simplify(dual / value))
        return ratios

    def test_the_euclidean_schwarzschild_chart_is_schwarzschilds_with_the_sign_of_g_tt_reversed(self):
        import sympy as sp
        import verify_metrics as vm
        _, reader, geo, symbols = self.geometry("schwarzschild")
        there = next(c for c in json.loads((DATA / "metrics" / "schwarzschild.json").read_text(encoding="utf-8"))
                     ["coordinates"] if c["id"] == "spherical")
        other = vm.Reader(there["coords"], [p["symbol"] for p in there["parameters"]], ())
        theirs = {tuple(e["indices"]): other(e["value"]) for e in there["metric_components"]}
        at = {other.symbol[c]: s for c, s in zip(("r", "\\theta", "\\phi"), symbols[1:])}
        at[other.parameters["r_s"]] = reader.parameters["r_s"]
        self.assertEqual(sp.simplify(geo.g[0, 0] + theirs[("t", "t")].subs(at)), 0)
        for k, name in enumerate(("r", "\\theta", "\\phi"), start=1):
            self.assertEqual(sp.simplify(geo.g[k, k] - theirs[(name, name)].subs(at)), 0)

    def test_each_period_is_two_pi_over_the_surface_gravity_of_its_nut_or_bolt(self):
        import sympy as sp
        # For f d tau^2 + dr^2/f the surface gravity where f vanishes is f'/2, and tau closes smoothly at 2 pi over it.
        for system, at, period in (("schwarzschild", "r_s", "4\\pi r_s"), ("taub_nut", "n", "8\\pi n"),
                                   ("taub_bolt", "2*n", "8\\pi n")):
            entry, reader, geo, symbols = self.geometry(system)
            r, tau = reader.symbol["r"], reader.symbol["\\tau"]
            k = entry["coords"].index("\\tau")
            f = geo.g[k, k]
            scale = reader.parameters["r_s" if system == "schwarzschild" else "n"]
            root = sp.sympify(at, locals={"r_s": scale, "n": scale})
            self.assertEqual(sp.simplify(f.subs(r, root)), 0)
            self.assertEqual(sp.simplify(geo.g[entry["coords"].index("r"), entry["coords"].index("r")] * f - 1), 0)
            kappa = sp.simplify(sp.diff(f, r).subs(r, root) / 2)
            self.assertEqual(sp.simplify(2 * sp.pi / kappa - reader(period)), 0)
            self.assertIn("\\tau \\in [0, " + period + ")", entry["domains"])

    def test_the_ricci_flat_charts_count_one_for_a_nut_and_two_for_a_bolt(self):
        import sympy as sp
        import verify_metrics as vm
        # For a Ricci flat metric the Euler density is K/(32 pi^2); the boundary terms of these three vanish.
        for system, start, period, euler in (("schwarzschild", "r_s", "4*pi*r_s", 2), ("taub_nut", "n", "8*pi*n", 1),
                                             ("taub_bolt", "2*n", "8*pi*n", 2)):
            entry, reader, geo, symbols = self.geometry(system)
            ricci = geo.ricci_ll()
            self.assertTrue(all(vm._at(ricci, (i, j)) == 0 for i in range(4) for j in range(4)), system)
            r, theta = reader.symbol["r"], reader.symbol["\\theta"]
            scale = reader.parameters["r_s" if system == "schwarzschild" else "n"]
            names = {"r_s": scale, "n": scale, "pi": sp.pi}
            K = sp.simplify(geo.kretschmann())
            self.assertEqual(vm.norm(K - reader(entry["kretschmann"].removeprefix("K = "))), 0)
            volume = r ** 2 if system == "schwarzschild" else r ** 2 - scale ** 2
            self.assertEqual(sp.simplify(geo.g.det() - volume ** 2 * sp.sin(theta) ** 2), 0)
            radial = sp.integrate(sp.simplify(K * volume), (r, sp.sympify(start, locals=names), sp.oo))
            bulk = radial * 2 * 2 * sp.pi * sp.sympify(period, locals=names)
            self.assertEqual(sp.simplify(bulk / (32 * sp.pi ** 2)), euler, system)

    def test_hawkings_taub_nut_is_self_dual_and_pages_is_not(self):
        import sympy as sp
        _, reader, geo, symbols = self.geometry("taub_nut")
        r, theta, n = reader.symbol["r"], reader.symbol["\\theta"], reader.parameters["n"]
        root = (r ** 2 - n ** 2) * sp.sin(theta)
        self.assertEqual(sp.simplify(geo.g.det() - root ** 2), 0)
        ratios = self.duality(geo, geo.riemann_llll(), root)
        self.assertEqual(len(ratios), 1)
        self.assertIn(ratios.pop(), (1, -1))
        _, reader, geo, symbols = self.geometry("taub_bolt")
        r, theta, n = reader.symbol["r"], reader.symbol["\\theta"], reader.parameters["n"]
        self.assertGreater(len(self.duality(geo, geo.riemann_llll(), (r ** 2 - n ** 2) * sp.sin(theta))), 1)
        # Page's bolt, r = 2n, is a sphere of area 12 pi n^2.
        k = [symbols.index(reader.symbol[c]) for c in ("\\theta", "\\phi")]
        element = sp.sqrt(sp.simplify((geo.g[k[0], k[0]] * geo.g[k[1], k[1]]).subs(r, 2 * n) / sp.sin(theta) ** 2))
        self.assertEqual(sp.simplify(element * 4 * sp.pi - 12 * sp.pi * n ** 2), 0)

    def test_the_multi_centre_chart_is_a_vacuum_for_two_centres_and_their_twist(self):
        import sympy as sp
        entry, reader, geo, symbols = self.geometry("multi_centre")
        rho, z = symbols[0], symbols[1]
        V, omega = reader.parameters["V"], reader.parameters["omega"]
        R1, R2 = sp.sqrt(rho ** 2 + (z + 2) ** 2), sp.sqrt(rho ** 2 + (z - 2) ** 2)
        potential = 1 + 2 / R1 + 2 / R2
        twist = 2 * (z + 2) / R1 + 2 * (z - 2) / R2
        # The field equation the chart's convention states.
        self.assertEqual(sp.simplify(sp.diff(twist, rho) - rho * sp.diff(potential, z)), 0)
        self.assertEqual(sp.simplify(sp.diff(twist, z) + rho * sp.diff(potential, rho)), 0)
        point = {rho: sp.Rational(7, 5), z: sp.Rational(3, 10)}
        for e in entry["ricci_tensor"]["variants"]["ll"]["nonzero"]:
            value = reader(e["value"]).subs({V: potential, omega: twist}, simultaneous=True).doit()
            self.assertLess(abs(value.subs(point).evalf(30)), 1e-25, e["indices"])
        # Without its twist the same potential is no vacuum.
        some = [reader(e["value"]).subs({V: potential, omega: 0}, simultaneous=True).doit().subs(point).evalf(30)
                for e in entry["ricci_tensor"]["variants"]["ll"]["nonzero"]]
        self.assertGreater(max(abs(v) for v in some), 1e-3)

    def test_the_projective_plane_is_einstein_and_half_flat_with_euler_number_three(self):
        import sympy as sp
        import verify_metrics as vm
        entry, reader, geo, symbols = self.geometry("cp2")
        r, theta, lam = reader.symbol["r"], reader.symbol["\\theta"], reader.parameters["Lambda"]
        ricci = geo.ricci_ll()
        self.assertTrue(all(vm.norm(vm._at(ricci, (i, j)) - lam * geo.g[i, j]) == 0 for i in range(4) for j in range(4)))
        w = 1 + lam * r ** 2 / 6
        root = r ** 3 * sp.sin(theta) / (8 * w ** 3)
        self.assertEqual(sp.simplify(geo.g.det() - root ** 2), 0)
        ratios = self.duality(geo, geo.weyl_llll(), root)
        self.assertEqual(len(ratios), 1)
        self.assertIn(ratios.pop(), (1, -1))
        # The Riemann tensor itself is not half flat, since the Ricci tensor does not vanish.
        self.assertGreater(len(self.duality(geo, geo.riemann_llll(), root)), 1)
        K = sp.simplify(geo.kretschmann())
        self.assertEqual(sp.simplify(K - 16 * lam ** 2 / 3), 0)
        # psi has the period 4 pi: the volume is 18 pi^2/Lambda^2, and for an Einstein space the Euler density is K/(32 pi^2).
        positive = sp.Symbol("lam", positive=True)
        volume = sp.integrate((r ** 3 / (8 * w ** 3)).subs(lam, positive), (r, 0, sp.oo)) * 2 * 2 * sp.pi * 4 * sp.pi
        self.assertEqual(sp.simplify(volume - 18 * sp.pi ** 2 / positive ** 2), 0)
        self.assertEqual(sp.simplify(K.subs(lam, positive) * volume / (32 * sp.pi ** 2)), 3)
        self.assertIn("$18\\pi^2/\\Lambda^2$", entry["parameters"][0]["description"])

    def test_the_distance_chart_is_the_chart_in_r_along_r_equals_l_tan_chi_over_l(self):
        import sympy as sp
        _, reader, geo, symbols = self.geometry("cp2")
        _, other, there, coords = self.geometry("cp2_distance")
        L, chi = other.parameters["L"], coords[0]
        image = [L * sp.tan(chi / L)] + coords[1:]
        J = sp.Matrix(4, 4, lambda i, j: sp.diff(image[i], coords[j]))
        pulled = J.T * geo.g.subs({**dict(zip(symbols, image)), reader.parameters["Lambda"]: 6 / L ** 2},
                                  simultaneous=True) * J
        point = {chi: sp.Rational(7, 10), coords[1]: sp.Rational(9, 10), L: sp.Rational(6, 5)}
        for i in range(4):
            for j in range(4):
                self.assertLess(abs((pulled[i, j] - there.g[i, j]).subs(point).evalf(30)), 1e-25)
        # The bolt, chi = pi L/2, is a sphere of radius L/2.
        self.assertEqual(sp.simplify(there.g[1, 1].subs(chi, sp.pi * L / 2) - L ** 2 / 4), 0)

    def test_pages_chart_is_einstein_for_every_nu_and_regular_exactly_on_the_quartic(self):
        import sympy as sp
        import verify_metrics as vm
        entry, reader, geo, symbols = self.geometry("page")
        chi, theta = reader.symbol["\\chi"], reader.symbol["\\theta"]
        lam, nu = reader.parameters["Lambda"], reader.parameters["nu"]
        ricci = geo.ricci_ll()
        self.assertTrue(all(vm.norm(vm._at(ricci, (i, j)) - lam * geo.g[i, j]) == 0 for i in range(4) for j in range(4)))
        quartic = nu ** 4 + 4 * nu ** 3 - 6 * nu ** 2 + 12 * nu - 3
        N = 3 + 6 * nu ** 2 - nu ** 4
        # With psi of period 4 pi the circle over a point of a bolt has the radius 2 sqrt(g_psipsi), and it
        # grows with the distance from the bolt at the rate 4 nu (3 + nu^2)/N on both.
        k = entry["coords"].index("\\psi")
        # Near a bolt g_psipsi = c chi^2, so the radius is 2 sqrt(c) chi at the distance sqrt(g_chichi) chi.
        ratio = sp.simplify(4 * geo.g[k, k] / (geo.g[0, 0] * sp.sin(chi) ** 2))
        for end in (0, sp.pi):
            self.assertEqual(sp.simplify(ratio.subs(chi, end) - (4 * nu * (3 + nu ** 2) / N) ** 2), 0)
        self.assertEqual(sp.expand(4 * nu * (3 + nu ** 2) - N - quartic), 0)
        # At the root the chart is Gibbons and Hawking's (3.25), which writes 1/(3 + nu^2)^2 for 16 nu^2/N^2.
        root = sp.nsolve(quartic, nu, 0.28, prec=40)
        self.assertLess(abs(root - NU), 1e-15)
        theirs = (3 * (1 + nu ** 2) / lam * (3 - nu ** 2 - nu ** 2 * (1 + nu ** 2) * sp.cos(chi) ** 2) * sp.sin(chi) ** 2
                  / (4 * (3 + nu ** 2) ** 2 * (1 - nu ** 2 * sp.cos(chi) ** 2)))
        at = {nu: root, lam: 1, chi: sp.Rational(11, 10)}
        self.assertLess(abs((geo.g[k, k] - theirs).subs(at).evalf(40)), 1e-30)
        self.assertGreater(abs((geo.g[k, k] - theirs).subs({**at, nu: sp.Rational(1, 2)}).evalf(40)), 1e-3)
        # Each bolt is a sphere of area 12 pi (1 - nu^4)/(N Lambda), their (3.27).
        self.assertEqual(sp.simplify(4 * sp.pi * geo.g[1, 1].subs(chi, 0) - 12 * sp.pi * (1 - nu ** 4) / (N * lam)), 0)
        self.assertEqual(sp.simplify(geo.g[1, 1].subs(chi, 0) - geo.g[1, 1].subs(chi, sp.pi)), 0)

    def test_pages_space_has_euler_number_four(self):
        import mpmath
        import sympy as sp
        entry, reader, geo, symbols = self.geometry("page")
        chi, theta = reader.symbol["\\chi"], reader.symbol["\\theta"]
        lam, nu = reader.parameters["Lambda"], reader.parameters["nu"]
        K = reader(entry["kretschmann"].removeprefix("K = "))
        volume = sp.sqrt(sp.simplify(geo.g.det() / sp.sin(theta) ** 2))
        density = sp.lambdify(chi, (K * volume).subs({lam: 1, nu: sp.Float(NU, 30)}), "mpmath")
        bulk = mpmath.quad(density, [0, mpmath.pi]) * 2 * 2 * mpmath.pi * 4 * mpmath.pi
        self.assertLess(abs(bulk / (32 * mpmath.pi ** 2) - 4), 1e-9)


if __name__ == "__main__":
    unittest.main()
