"""Nordstrom's scalar theory of gravity as Einstein and Fokker wrote it, a conformal factor on
Minkowski's metric: python3 -m unittest discover -s _tools

What the texts and drawings of `nordstrom_scalar` state, held on the published files. The free
chart is held to R = -6 box Phi/Phi^3 and to no Weyl tensor, and each other chart to being that
chart at its own Phi, a solution of the field equation R = 24 pi G T/c^4: the point mass and the
uniform field vacua, R = 0, with a Ricci tensor that does not vanish, and the dust universe with
R Phi^3 constant. The point mass is held to Deruelle's post-Newtonian numbers gamma = -1 and
beta = 1/2, to the published metric of `ppn_metric` at those numbers, to the redshift 1/(1 - m/r),
and, by running a planet with the published Christoffel symbols, to a perihelion that turns back
by a sixth of Einstein's advance. A body dropped in the uniform field is run the same way to the
fall of Giulini's (26a), pi c/4a on the body's own clock, and a galaxy's clock in the dust universe to 4L/3c. The
spacetime diagrams are held to cones at 45 degrees in every chart, the conformal diagrams to their
singular edges, and the embedding diagram to the closed form of the point mass's surface in
Minkowski space and to the discs of radius Phi r. The tests that read the published components
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

CHARTS = ("conformal", "spherical", "uniform", "dust")


def published(system, metric_id="nordstrom_scalar"):
    metric = json.loads((DATA / "metrics" / f"{metric_id}.json").read_text(encoding="utf-8"))
    return next(c for c in metric["coordinates"] if c["id"] == system)


def drawn(kind):
    return json.loads((DATA / kind / "nordstrom_scalar.json").read_text(encoding="utf-8"))


def rk4(f, y, h, steps, stop=None):
    """The classical Runge-Kutta steps of y' = f(y), handed back as the list of states."""
    out = [list(y)]
    for _ in range(steps):
        k1 = f(y)
        k2 = f([a + h / 2 * b for a, b in zip(y, k1)])
        k3 = f([a + h / 2 * b for a, b in zip(y, k2)])
        k4 = f([a + h * b for a, b in zip(y, k3)])
        y = [a + h / 6 * (p + 2 * q + 2 * r + s) for a, p, q, r, s in zip(y, k1, k2, k3, k4)]
        out.append(list(y))
        if stop and stop(y):
            break
    return out


@unittest.skipUnless(HAS_SYMPY, "sympy is not installed")
class Theory(unittest.TestCase):
    """The published components, read as the checker reads them, in the chart x^0 = ct with c = 1."""

    @classmethod
    def setUpClass(cls):
        import sympy as sp
        import verify_metrics as vm
        cls.sp, cls.vm = sp, vm
        cls.entry, cls.reader, cls.g, cls.x, cls.P = {}, {}, {}, {}, {}
        for system in CHARTS:
            entry = published(system)
            reader = vm.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]], ())
            names = entry["coords"]
            there = {tuple(e["indices"]): e["value"] for e in entry["metric_components"]}
            g = sp.Matrix(4, 4, lambda i, j: reader(there.get((names[i], names[j]), "0")))
            cls.entry[system], cls.reader[system] = entry, reader
            cls.g[system] = g.subs(reader.c, 1)
            cls.x[system] = [reader.symbol[n] for n in names]
            cls.P[system] = reader.parameters

    def scalar(self, system):
        text = self.entry[system]["ricci_scalar"].partition("=")[2]
        return self.reader[system](text).subs(self.reader[system].c, 1)

    def tensor(self, system, field, variant):
        reader, names = self.reader[system], self.entry[system]["coords"]
        out = {}
        for c in self.entry[system][field]["variants"][variant]["nonzero"]:
            out[tuple(names.index(n) for n in c["indices"])] = reader(c["value"]).subs(reader.c, 1)
        return out

    def test_every_chart_is_a_conformal_factor_times_minkowskis_metric(self):
        sp = self.sp
        flat = {"conformal": sp.diag(-1, 1, 1, 1), "uniform": sp.diag(-1, 1, 1, 1)}
        for system in CHARTS:
            g, x = self.g[system], self.x[system]
            eta = flat.get(system, sp.diag(-1, 1, x[1] ** 2, x[1] ** 2 * sp.sin(x[2]) ** 2))
            factor = sp.simplify(-g[0, 0])
            self.assertEqual((g - factor * eta).applyfunc(sp.simplify), sp.zeros(4, 4), system)
            self.assertEqual(self.entry[system]["weyl_tensor"]["variants"]["llll"]["nonzero"], [], system)

    def test_the_free_charts_curvature_scalar_is_minus_six_box_phi_over_phi_cubed(self):
        sp = self.sp
        t, x, y, z = self.x["conformal"]
        Phi = self.P["conformal"]["Phi"]
        box = sum(sp.diff(Phi, s, 2) for s in (x, y, z)) - sp.diff(Phi, t, 2)
        self.assertEqual(sp.simplify(self.scalar("conformal") + 6 * box / Phi ** 3), 0)

    def test_the_point_mass_and_the_uniform_field_solve_the_wave_equation_and_are_not_flat(self):
        sp = self.sp
        t, x, y, z = self.x["conformal"]
        Phi = self.P["conformal"]["Phi"]
        m, a = sp.symbols("m a", positive=True)
        general = self.scalar("conformal")
        for factor in (1 - m / sp.sqrt(x ** 2 + y ** 2 + z ** 2), a * z):
            self.assertEqual(sp.simplify(general.subs(Phi, factor).doit()), 0)
        for system in ("spherical", "uniform"):
            self.assertEqual(self.scalar(system), 0, system)
            self.assertTrue(self.tensor(system, "ricci_tensor", "ll"), system)

    def test_the_dust_universe_keeps_r_phi_cubed_constant_and_its_galaxies_fall_freely(self):
        sp = self.sp
        t = self.x["dust"][0]
        L = self.P["dust"]["L"]
        factor = 1 - t ** 2 / L ** 2
        self.assertEqual(sp.simplify(self.scalar("dust") * factor ** 3 + 12 / L ** 2), 0)
        gamma = self.tensor("dust", "christoffel", "ull")
        self.assertFalse([k for k in gamma if k[0] != 0 and k[1:] == (0, 0)])
        # The free chart at Phi = 1 - c^2t^2/L^2 has the dust universe's curvature scalar.
        Phi = self.P["conformal"]["Phi"]
        T = self.x["conformal"][0]
        general = self.scalar("conformal").subs(Phi, 1 - T ** 2 / L ** 2).doit()
        self.assertEqual(sp.simplify(general.subs(T, t) - self.scalar("dust")), 0)

    def test_a_galaxys_clock_runs_four_thirds_of_l_over_c_from_the_bang_to_the_crunch(self):
        sp = self.sp
        t = self.x["dust"][0]
        L = self.P["dust"]["L"]
        rate = sp.sqrt(sp.factor(-self.g["dust"][0, 0]))
        self.assertEqual(sp.simplify(sp.integrate(sp.simplify(rate.subs(L, 1)), (t, -1, 1))), sp.Rational(4, 3))

    def test_the_point_mass_has_gamma_minus_one_and_beta_one_half(self):
        sp = self.sp
        r = self.x["spherical"][1]
        m = self.P["spherical"]["m"]
        U = sp.Symbol("U", positive=True)
        g = self.g["spherical"]
        lapse = sp.expand((-g[0, 0]).subs(r, m / U))
        space = sp.expand(g[1, 1].subs(r, m / U))
        # -g_tt = 1 - 2U + 2 beta U^2 and g_rr = 1 + 2 gamma U + ..., U = m/r.
        self.assertEqual(lapse, 1 - 2 * U + 2 * sp.Rational(1, 2) * U ** 2)
        self.assertEqual(space.coeff(U, 1), 2 * (-1))
        # Mercury's perihelion moves by (2 + 2 gamma - beta)/3 of Einstein's advance.
        self.assertEqual((2 + 2 * (-1) - sp.Rational(1, 2)) / 3, -sp.Rational(1, 6))

    def test_it_is_the_post_newtonian_metric_at_those_numbers(self):
        sp, vm = self.sp, self.vm
        there = published("isotropic", "ppn_metric")
        reader = vm.Reader(there["coords"], [p["symbol"] for p in there["parameters"]], ())
        ppn = {tuple(e["indices"]): reader(e["value"]) for e in there["metric_components"]}
        at = {reader.parameters["beta"]: sp.Rational(1, 2), reader.parameters["gamma"]: -1,
              reader.parameters["m"]: self.P["spherical"]["m"], reader.symbol["r"]: self.x["spherical"][1]}
        g = self.g["spherical"]
        m, r = self.P["spherical"]["m"], self.x["spherical"][1]
        self.assertEqual(sp.simplify(ppn[("t", "t")].subs(at) - g[0, 0]), 0)
        self.assertEqual(sp.simplify(g[1, 1] - ppn[("r", "r")].subs(at) - m ** 2 / r ** 2), 0)

    def test_a_clock_at_rest_outside_the_body_ticks_at_one_minus_m_over_r(self):
        sp = self.sp
        r = self.x["spherical"][1]
        m = self.P["spherical"]["m"]
        self.assertEqual(sp.simplify(-self.g["spherical"][0, 0] - (1 - m / r) ** 2), 0)

    def test_a_body_at_rest_in_the_uniform_field_weighs_c_to_the_fourth_over_a_z_squared(self):
        sp = self.sp
        z = self.x["uniform"][3]
        a = self.P["uniform"]["a"]
        gamma = self.tensor("uniform", "christoffel", "ull")
        factor = a * z
        self.assertEqual(sp.simplify(factor ** 2 + self.g["uniform"][0, 0]), 0)
        # Gamma^z_tt (u^t)^2 with u^t = 1/Phi, measured with sqrt(g_zz) = Phi.
        weight = sp.simplify(gamma[(3, 0, 0)] / factor ** 2 * factor)
        self.assertEqual(sp.simplify(weight - 1 / (a * z ** 2)), 0)
        self.assertEqual(sp.simplify(weight.subs(z, 1 / a) - a), 0)

    def runner(self, system, values):
        """The geodesic equations of the published Christoffel symbols, as a function of the state
        (x^mu, dx^mu/d tau), at the parameter values given."""
        sp = self.sp
        gamma = self.tensor(system, "christoffel", "ull")
        subs = {self.P[system][k]: v for k, v in values.items()}
        fn = {k: sp.lambdify(self.x[system], v.subs(subs), "math") for k, v in gamma.items()}

        def step(y):
            x, v = y[:4], y[4:]
            acc = [0.0] * 4
            for (k, i, j), f in fn.items():
                acc[k] -= f(*x) * v[i] * v[j]
            return list(v) + acc
        return step

    def test_a_body_dropped_from_where_phi_is_one_reaches_the_singular_plane_after_pi_c_over_4a(self):
        """In units of c^2/a a body released at z = 1 has dz/d tau = -sqrt(1 - z^2)/z^2 on its own
        clock, so it reads tau = pi/4 - (arcsin z - z sqrt(1 - z^2))/2 at the height z and pi/4 on
        the singular plane. Giulini's (26a) and (29a) are the same fall in the proper time of the
        flat background, s with d tau = Phi ds: z = sqrt(1 - s^2), which ends at s = c/a."""
        step = self.runner("uniform", {"a": 1})
        h = 1e-4
        path = rk4(step, [0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0], h, 7700)
        flat = 0.0
        for n, (a, b) in enumerate(zip(path, path[1:]), 1):
            flat += h * 2 / (a[3] + b[3])
            if n in (2000, 5000, 7700):
                z = b[3]
                self.assertAlmostEqual(n * h, math.pi / 4 - (math.asin(z) - z * math.sqrt(1 - z * z)) / 2, places=7)
                self.assertAlmostEqual(z, math.sqrt(1 - flat * flat), places=5)
        self.assertLess(path[7700][3], 0.36)
        self.assertLess(7700 * h, math.pi / 4)

    def test_a_planets_perihelion_turns_back_by_a_sixth_of_einsteins_advance(self):
        """A planet run with the published Christoffel symbols at m = 1: between two perihelia it
        goes round by 2 pi/sqrt(1 + m^2/l^2), with l its angular momentum per unit mass, short of a
        full turn by pi m^2/l^2 to leading order, where Einstein's theory adds 6 pi m^2/l^2."""
        step = self.runner("spherical", {"m": 1})
        r0 = 400.0                      # the perihelion, in units of m
        factor = 1 - 1 / r0
        speed = 1.2 * math.sqrt(1 / r0)   # faster than circular, so r0 is the perihelion
        gamma = 1 / math.sqrt(1 - speed ** 2)
        ut, uphi = gamma / factor, gamma * speed / (factor * r0)
        ell = factor ** 2 * r0 ** 2 * uphi
        h = 2.0
        y = [0.0, r0, math.pi / 2, 0.0, ut, 0.0, 0.0, uphi]
        path = rk4(step, y, h, 400000, stop=lambda s: s[5] > 0 and s[1] < r0 * 1.02 and s[3] > 5.0)
        # Past the aphelion, the next perihelion is where dr/d tau changes sign from - to +.
        turn = None
        for a, b in zip(path, path[1:]):
            if a[3] > 3.0 and a[5] < 0 <= b[5]:
                w = -a[5] / (b[5] - a[5])
                turn = a[3] + w * (b[3] - a[3])
                break
        self.assertIsNotNone(turn)
        self.assertAlmostEqual(turn, 2 * math.pi / math.sqrt(1 + 1 / ell ** 2), places=7)
        shift, einstein = turn - 2 * math.pi, 6 * math.pi / ell ** 2
        self.assertLess(shift, 0)
        self.assertAlmostEqual(shift / einstein, -1 / 6, places=3)


class SpacetimeDiagrams(unittest.TestCase):
    def test_every_cone_in_every_chart_is_minkowskis(self):
        for system, views in drawn("diagrams")["systems"].items():
            for view in views:
                x0, x1, y0, y1 = view["box"]
                self.assertTrue(view["cones"], system)
                for cone in view["cones"]:
                    for edge in (cone["a"], cone["b"]):
                        across, up = edge[0] * (x1 - x0), edge[1] * (y1 - y0)
                        self.assertGreater(up, 0, system)
                        self.assertAlmostEqual(abs(across) / up, 1.0, places=3, msg=system)

    def test_every_ray_is_a_straight_line_at_45_degrees(self):
        for system, views in drawn("diagrams")["systems"].items():
            for view in views:
                x0, x1, y0, y1 = view["box"]
                for family, sign in (("P", None), ("M", None)):
                    for ray in view["rays"][family]:
                        (ax, ay), (bx, by) = ray[0], ray[-1]
                        across, up = (bx - ax) * (x1 - x0), (by - ay) * (y1 - y0)
                        self.assertAlmostEqual(abs(across), abs(up), delta=2e-3 * (x1 - x0), msg=system)

    def test_the_singular_edges_are_where_the_conformal_factor_vanishes(self):
        views = {k: v[0] for k, v in drawn("diagrams")["systems"].items()}
        marked = {k: sorted(e for m in v["markers"] if m["kind"] == "singular" for e in m["edges"])
                  for k, v in views.items()}
        self.assertEqual(marked, {"conformal": [], "spherical": ["left"], "uniform": ["left"],
                                  "dust": ["bottom", "top"]})
        self.assertEqual(views["spherical"]["box"][0], 1)      # r = m
        self.assertEqual(views["uniform"]["box"][0], 0)        # z = 0
        self.assertEqual(views["dust"]["box"][2:], [-1, 1])    # ct = -L and L

    def test_the_free_chart_is_drawn_for_any_factor_with_a_plane_wave(self):
        view = drawn("diagrams")["systems"]["conformal"][0]
        self.assertTrue(view["input"].startswith("Any $\\Phi$"))
        source = (ROOT / "_tools" / "derivations" / "null_rays.py").read_text(encoding="utf-8")
        self.assertIn('functions={"Phi": "1 + cos(x - t)/2"}, any_factor="Phi"', source)
        # The wave solves the wave equation: it depends on x - ct alone.
        wave = lambda t, x: 1 + math.cos(x - t) / 2
        for t, x in ((0.3, 1.1), (-2.0, 0.4), (1.7, -3.0)):
            h = 1e-4
            d2t = (wave(t + h, x) - 2 * wave(t, x) + wave(t - h, x)) / h ** 2
            d2x = (wave(t, x + h) - 2 * wave(t, x) + wave(t, x - h)) / h ** 2
            self.assertAlmostEqual(d2x - d2t, 0.0, places=6)
            self.assertGreater(wave(t, x), 0)


class ConformalDiagrams(unittest.TestCase):
    def setUp(self):
        self.views = {v["id"]: v for v in drawn("conformal")["views"]}

    def classes(self, view):
        return [(layer["kind"], layer["class"]) for layer in self.views[view]["layers"]]

    def test_the_point_mass_and_the_uniform_field_have_a_timelike_singularity_and_null_infinity(self):
        for view in ("spherical", "uniform"):
            zigs = [layer for layer in self.views[view]["layers"] if layer["kind"] == "zig"]
            self.assertEqual(len(zigs), 1)
            self.assertEqual({round(p[0], 6) for p in zigs[0]["points"]}, {0.0})        # the edge X = 0
            self.assertIn(("line", "scri"), self.classes(view))

    def test_the_dust_universe_lies_between_two_spacelike_singularities_that_meet_at_spatial_infinity(self):
        layers = self.views["dust"]["layers"]
        zigs = [layer["points"] for layer in layers if layer["kind"] == "zig"]
        self.assertEqual(len(zigs), 2)
        self.assertNotIn(("line", "scri"), self.classes("dust"))
        for curve in zigs:
            # ct = +-L: tan p tan q = ..., p = arctan(+-1 - r), q = arctan(+-1 + r), in (X, T) = (q - p, p + q).
            sign = 1 if curve[len(curve) // 2][1] > 0 else -1
            for X, T in curve:
                p, q = (T - X) / 2, (T + X) / 2
                if abs(X - math.pi) > 1e-3:
                    self.assertAlmostEqual((math.tan(p) + math.tan(q)) / 2, sign, delta=0.03 * (1 + abs(math.tan(q))))
            ends = sorted(curve, key=lambda at: at[0])
            self.assertAlmostEqual(ends[0][0], 0.0, places=3)
            self.assertAlmostEqual(abs(ends[0][1]), math.pi / 2, places=3)
            self.assertAlmostEqual(ends[-1][0], math.pi, places=2)
            self.assertAlmostEqual(ends[-1][1], 0.0, places=2)
            # Spacelike: the curve never climbs faster than it runs out.
            for (xa, ta), (xb, tb) in zip(curve, curve[1:]):
                self.assertLessEqual(abs(tb - ta), abs(xb - xa) + 1e-3)

    def test_the_plane_wave_fills_minkowskis_diamond(self):
        self.assertNotIn("zig", {kind for kind, _ in self.classes("conformal")})
        self.assertEqual(self.classes("conformal").count(("line", "scri")), 4)


class EmbeddingDiagram(unittest.TestCase):
    def setUp(self):
        self.views = {v["id"]: v for v in drawn("embedding")["views"]}

    @staticmethod
    def height(r):
        root = math.sqrt(2 * r - 1)
        return 2 * root - 2 * math.atan(root)

    def test_the_point_masss_surface_stands_in_minkowski_space_at_its_closed_form(self):
        sheet = next(p for p in self.views["point_mass"]["surfaces"][0]["pieces"] if p["id"] == "sheet")
        self.assertEqual(sheet["space"], "minkowski")
        points = sheet["points"]
        r0, _, z0 = points[0]
        for r, rho, z in points:
            self.assertAlmostEqual(rho, r - 1, places=9)
            self.assertAlmostEqual(z - z0, self.height(r) - self.height(r0), places=7)

    def test_every_chord_of_it_has_its_proper_length_in_minkowski_space(self):
        sheet = next(p for p in self.views["point_mass"]["surfaces"][0]["pieces"] if p["id"] == "sheet")
        for (ra, rhoa, za), (rb, rhob, zb) in zip(sheet["points"], sheet["points"][1:]):
            chord = math.sqrt((rhob - rhoa) ** 2 - (zb - za) ** 2)
            proper = (rb - ra) - math.log(rb / ra)          # the integral of (1 - m/r) dr
            self.assertAlmostEqual(chord / proper, 1.0, delta=2e-3)

    def test_far_out_it_is_flamms_height_for_the_radius_two_m(self):
        """Z -> 2 sqrt(2mr) - pi m: Flamm's 2 sqrt(r_s (r - r_s)) at r_s = 2m to a constant."""
        for r in (1e6, 1e8):
            self.assertAlmostEqual(self.height(r) - (2 * math.sqrt(2 * r) - math.pi), 0.0, delta=3 / math.sqrt(r))

    def test_it_leaves_the_singular_point_along_the_light_cone(self):
        for rho in (1e-2, 1e-3):
            rise = self.height(1 + rho) - self.height(1)
            self.assertAlmostEqual(rise / rho, 1.0, delta=rho)
        cone = next(p for p in self.views["point_mass"]["surfaces"][0]["pieces"] if p["id"] == "cone")
        self.assertTrue(cone["reference"])
        (_, rho_a, z_a), (_, rho_b, z_b) = cone["points"][0], cone["points"][-1]
        self.assertAlmostEqual((z_b - z_a) / (rho_b - rho_a), 1.0, places=9)

    def test_the_dust_universes_moments_are_flat_discs_of_radius_phi_r(self):
        view = self.views["dust"]
        frames = view["movie"]["frames"]
        self.assertEqual([f["value"] for f in frames][0], -0.9)
        self.assertEqual([f["value"] for f in frames][-1], 0.9)
        widest = max(frames, key=lambda f: f["pieces"][0]["points"][-1][1])
        self.assertEqual(widest["value"], 0.0)
        for frame in frames:
            factor = 1 - frame["value"] ** 2
            for r, rho, z in frame["pieces"][0]["points"]:
                self.assertAlmostEqual(rho, factor * r, places=7)
                self.assertAlmostEqual(z, 0.0, places=7)
        self.assertEqual([s["time"] for s in view["surfaces"]], [-0.9, -0.6, 0.0, 0.6, 0.9])


if __name__ == "__main__":
    unittest.main()
