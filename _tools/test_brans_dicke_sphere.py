"""The static sphere of Brans and Dicke's theory, Brans's class I:
.venv.noindex/bin/python -m unittest discover -s _tools

What its texts and drawings state, held on the published files. The first two classes need
nothing but the files: the constants the three charts name for one solution and the coupling
constant each stands for, the numbers the History and the captions quote, the circles of the
embedding diagram and its level circle, the triangle of the conformal diagrams, and the relations
each side answers. The third reads the published tensors through the checker's Reader and holds
them to the vacuum equations of the theory, to Fisher's published metric in Dicke's units, to the
post-Newtonian limit, and to Schwarzschild's metric without the field; it needs sympy and is
skipped where it is absent.
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

CHARTS = ["isotropic", "spherical", "harmonic"]
# The values every diagram is drawn at, in each chart's own letters.
C, LAMBDA = -1 / 4, 1.0
M, N = 0.0, 1 / 4
K, B_HARMONIC, S = 1.0, 7 / 8, -1 / 4


def load(folder, metric_id="brans_dicke_sphere"):
    return json.loads((DATA / folder / f"{metric_id}.json").read_text(encoding="utf-8"))


def omega_isotropic(c, lam):
    """Brans and Dicke's (33), lambda^2 = (C + 1)^2 - C (1 - omega C/2), solved for omega."""
    return 2 * (lam ** 2 - c ** 2 - c - 1) / c ** 2


def omega_spherical(m, n):
    return -2 * (m ** 2 + n ** 2 + m * n + m - n) / (m + n) ** 2


def omega_harmonic(k, b, s):
    return 2 * (k ** 2 - b ** 2) / s ** 2 - 3 / 2


def tortoise(r):
    """r_* in units of r_0 at m = 0 and n = 1/4, zero at r = 2: the integral of A^(-7/8) dr, by
    Simpson's rule in w = A^(1/8), where it is the integral of 16/(1 - w^8)^2 dw from 0."""
    w = (1 - 2 / r) ** 0.125
    steps = 2000
    h = w / steps
    f = lambda x: 16 / (1 - x ** 8) ** 2
    return h / 3 * (f(0) + f(w) + sum((4 if i % 2 else 2) * f(i * h) for i in range(1, steps)))


class Constants(unittest.TestCase):
    def test_the_charts_are_the_three_the_literature_uses(self):
        self.assertEqual([c["id"] for c in load("metrics")["coordinates"]], CHARTS)

    def test_the_three_sets_of_letters_name_one_solution(self):
        """m = 1/lambda - 1 and n = 1 - (C + 1)/lambda; b = (C + 2) k/2 lambda and s = k C/lambda."""
        self.assertAlmostEqual(1 / LAMBDA - 1, M)
        self.assertAlmostEqual(1 - (C + 1) / LAMBDA, N)
        self.assertAlmostEqual((C + 2) * K / (2 * LAMBDA), B_HARMONIC)
        self.assertAlmostEqual(K * C / LAMBDA, S)

    def test_every_chart_is_drawn_at_omega_equal_six(self):
        self.assertAlmostEqual(omega_isotropic(C, LAMBDA), 6.0)
        self.assertAlmostEqual(omega_spherical(M, N), 6.0)
        self.assertAlmostEqual(omega_harmonic(K, B_HARMONIC, S), 6.0)

    def test_the_three_formulas_for_omega_agree_off_the_drawn_values(self):
        for c, lam in ((-0.125, math.sqrt(15) / 4), (-0.5, 0.9), (0.3, 1.4)):
            m, n = 1 / lam - 1, 1 - (c + 1) / lam
            b, s = (c + 2) / (2 * lam), c / lam
            want = omega_isotropic(c, lam)
            self.assertAlmostEqual(omega_spherical(m, n), want, places=9)
            self.assertAlmostEqual(omega_harmonic(1.0, b, s), want, places=9)

    def test_a_body_of_weak_gravity_has_the_constants_the_captions_quote(self):
        """C = -1/(2 + omega) and lambda = sqrt((2 omega + 3)/(2 omega + 4)): at omega = 6 the power
        of the isotropic cones is 0.94, n is 0.096, the level circle is r = 2.005 r_0 and its
        radius 1.5 r_0."""
        omega = 6
        c, lam = -1 / (2 + omega), math.sqrt((2 * omega + 3) / (2 * omega + 4))
        self.assertAlmostEqual(omega_isotropic(c, lam), omega, places=9)
        self.assertAlmostEqual(c + 1, (1 + omega) / (2 + omega))            # the post-Newtonian gamma
        self.assertAlmostEqual((c + 2) / lam - 1, 0.94, places=2)
        n = 1 - (c + 1) / lam
        self.assertAlmostEqual(n, 0.096, places=3)
        level = (2 - n) ** 2 / (2 * (1 - n))
        self.assertAlmostEqual(level, 2.005, places=3)
        self.assertAlmostEqual(level * (1 - 2 / level) ** (n / 2), 1.5, places=1)

    def test_the_numbers_of_the_history(self):
        """Light bends by (3 + 2 omega)/(4 + 2 omega) = (1 + gamma)/2 of Einstein's angle; at omega = 6
        the perihelion advances at 11/12 of Einstein's rate, 39.4 of his 43 seconds, which is 8
        percent short; and the Sun's J_2 moves it by 3e-4 J_2/1e-7 of that rate, under a tenth of
        a percent."""
        history = load("metrics")["history"]
        for omega in (6, 40000, -1):
            gamma = (1 + omega) / (2 + omega)
            self.assertAlmostEqual((3 + 2 * omega) / (4 + 2 * omega), (1 + gamma) / 2)
        share = (4 + 3 * 6) / (6 + 3 * 6)
        self.assertAlmostEqual(share * 43.0, 39.4, places=1)
        self.assertAlmostEqual(1 - share, 0.08, places=2)
        self.assertLess(3e-4 * 2.2, 1e-3)
        for quoted in ("$39.4$", "$42.6 \\pm 0.9$", "$\\omega \\geq 6$", "$(4 + 3\\omega)/(6 + 3\\omega)$",
                       "$(3 + 2\\omega)/(4 + 2\\omega)$", "8 percent"):
            self.assertIn(quoted, history)
        self.assertLessEqual(history.count("?"), 3)

    def test_the_relations_are_answered_from_the_other_side(self):
        kinds = {"conformal": "conformal", "special_case": "generalisation", "generalisation": "special_case",
                 "family": "family"}
        related = load("metrics")["related"]
        self.assertEqual([r["id"] for r in related],
                         ["fisher_jnw", "schwarzschild", "ppn_metric", "nordstrom_scalar", "dilaton_black_hole",
                          "morris_thorne"])
        for r in related:
            back = next(x for x in load("metrics", r["id"])["related"] if x["id"] == "brans_dicke_sphere")
            self.assertEqual(back["kind"], kinds[r["kind"]], r["id"])


class Drawings(unittest.TestCase):
    def test_the_tortoise_coordinate_gives_the_times_the_captions_quote(self):
        """18 r_0/c from r = 4 r_0; 29 B/c from rho = 3B, which is r = 8 r_0/3 with times doubled;
        13 k/c from ku = 1; and 6.2 r_0/c across the last thousandth of r_0."""
        self.assertAlmostEqual(tortoise(4.0), 17.51, places=2)
        self.assertEqual(round(tortoise(4.0)), 18)
        self.assertEqual(round(2 * tortoise(3 * (1 + 1 / 3) ** 2 / 2)), 29)
        self.assertEqual(round(tortoise(2 / (1 - math.exp(-2)))), 13)
        self.assertAlmostEqual(tortoise(2.001), 6.2, places=1)
        # The power of A is 7/8: at r = 3 r_0 the cones are 3^(1/8) = 1.15 of Schwarzschild's.
        self.assertAlmostEqual((1 / 3) ** (7 / 8 - 1), 1.15, places=2)

    def test_every_spacetime_diagram_keeps_its_window_between_one_to_two_and_two_to_one(self):
        systems = load("diagrams")["systems"]
        self.assertEqual(sorted(systems), sorted(CHARTS))
        for system, (view,) in systems.items():
            X0, X1, Y0, Y1 = view["box"]
            self.assertTrue(0.5 <= (X1 - X0) / (Y1 - Y0) <= 2.0, system)
            self.assertEqual(len(view["slices"]), 1, system)

    def test_the_faint_lines_of_the_harmonic_plane_are_the_spheres_the_caption_names(self):
        """Areal radius 2w/(1 - w^8) with w = e^(-u/4): 4k, 2k, 1.5 k and k at ku = 0.31, 0.83, 1.4, 2.8."""
        radius = lambda u: 2 * math.exp(-u / 4) / (1 - math.exp(-2 * u))
        for u, want in ((0.3106, 4.0), (0.8347, 2.0), (1.4011, 1.5), (2.7878, 1.0)):
            self.assertAlmostEqual(radius(u), want, places=3)

    def test_the_embedded_circles_have_the_radius_of_their_spheres(self):
        """Every point written is (u, rho, z) with rho = r A^(1/8) = 2w/(1 - w^8), w = e^(-u/4)."""
        (view,) = load("embedding")["views"]
        (surface,) = view["surfaces"]
        for piece in surface["pieces"]:
            self.assertEqual((piece["system"], piece["coordinate"]), ("harmonic", "u"))
            for u, rho, _ in piece["points"]:
                w = math.exp(-u / 4)
                self.assertAlmostEqual(rho, 2 * w / (1 - w ** 8), delta=2e-6 if piece["id"] == "far" else 1e-12)

    def test_the_surface_lies_level_on_the_circle_r_equal_49_r0_over_24(self):
        """A = 1/49 there, u = ln 7, and the radius is (49/24) 49^(-1/8) r_0; the far piece ends and
        the near piece, in Minkowski space, begins on it at one height."""
        (surface,) = load("embedding")["views"][0]["surfaces"]
        far, near = surface["pieces"]
        self.assertEqual(near.get("space"), "minkowski")
        self.assertIsNone(far.get("space"))
        want = 49 / 24 * 49 ** -0.125
        n = N
        self.assertAlmostEqual((2 - n) ** 2 / (2 * (1 - n)), 49 / 24)
        for point in (far["points"][-1], near["points"][0]):
            self.assertAlmostEqual(point[0], math.log(7), places=12)
            self.assertAlmostEqual(point[1], want, places=6)
            self.assertAlmostEqual(point[2], 0.0, places=6)

    def test_the_surface_nears_the_singular_point_along_a_light_cone(self):
        """In Minkowski space dZ/d rho -> 1: over the last chords written the two differ by 64 A/2."""
        near = load("embedding")["views"][0]["surfaces"][0]["pieces"][1]
        (u0, rho0, z0), (u1, rho1, z1) = near["points"][-2:]
        slope = (z1 - z0) / (rho1 - rho0)
        self.assertAlmostEqual(slope, 1.0, places=6)
        self.assertLess(slope, 1.0 + 1e-9)
        self.assertAlmostEqual(u1, 12.0)
        self.assertAlmostEqual(rho1, 2 * math.exp(-3), places=9)

    def test_each_conformal_view_is_minkowskis_triangle_with_the_singularity_on_its_axis(self):
        views = load("conformal")["views"]
        self.assertEqual([v["id"] for v in views], CHARTS)
        for v in views:
            region = next(l for l in v["layers"] if l["kind"] == "fill" and l["class"] == "region")
            self.assertEqual([[round(x, 4), round(y, 4)] for x, y in region["points"]],
                             [[0.0, -3.1416], [3.1416, 0.0], [0.0, 3.1416]])
            singular = [l for l in v["layers"] if l.get("class") == "singular"]
            self.assertEqual(len(singular), 1, v["id"])
            self.assertEqual(sum(1 for l in v["layers"] if l.get("class") == "r"), 7, v["id"])

    def test_the_first_line_of_the_spherical_view_stands_where_its_caption_says(self):
        """r = 2.001 r_0 crosses t = 0 at X = 2 arctan(r_*/l) with l = 16 r_0."""
        view = next(v for v in load("conformal")["views"] if v["id"] == "spherical")
        first = next(l for l in view["layers"] if l.get("class") == "r")
        X = min(first["points"], key=lambda p: abs(p[1]))[0]
        self.assertAlmostEqual(X, 2 * math.atan(tortoise(2.001) / 16), places=2)


@unittest.skipUnless(HAS_SYMPY, "the Reader needs sympy: .venv.noindex/bin/python -m unittest discover -s _tools")
class Geometry(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import sympy
        import verify_metrics
        cls.sp, cls.vm = sympy, verify_metrics
        cls.charts = {c["id"]: c for c in load("metrics")["coordinates"]}

    def components(self, system, field, variant=None, metric_id=None):
        chart = self.charts[system] if metric_id is None else next(
            c for c in load("metrics", metric_id)["coordinates"] if c["id"] == system)
        reader = self.vm.Reader(chart["coords"], [p["symbol"] for p in chart["parameters"]], ())
        block = chart[field]
        entries = block["variants"][variant]["nonzero"] if variant else block
        return reader, {tuple(e["indices"]): reader(e["value"]) for e in entries}

    def test_the_spherical_chart_solves_the_vacuum_equations_of_the_theory(self):
        """R_ab = omega d_a phi d_b phi/phi^2 + nabla_a nabla_b phi/phi with phi = A^(-(m + n)/2) and
        omega = -2 (m^2 + n^2 + mn + m - n)/(m + n)^2, at exact rational values of m, n and r."""
        sp = self.sp
        reader, ricci = self.components("spherical", "ricci_tensor", "ll")
        _, gamma = self.components("spherical", "christoffel", "ull")
        r, r0, m, n = reader.symbol["r"], reader.parameters["r_0"], reader.parameters["m"], reader.parameters["n"]
        theta = reader.symbol["\\theta"]
        field = (1 - 2 * r0 / r) ** (-(m + n) / 2)
        omega = -2 * (m ** 2 + n ** 2 + m * n + m - n) / (m + n) ** 2
        log, second = sp.diff(field, r) / field, sp.diff(field, r, 2) / field
        at = {r0: 1, m: sp.Rational(2, 7), n: sp.Rational(1, 3), r: sp.Rational(17, 5), theta: sp.Rational(6, 5)}
        for a in ("t", "r", "\\theta", "\\phi"):
            source = (omega * log ** 2 + second if a == "r" else 0) - gamma.get(("r", a, a), 0) * log
            got = ricci.get((a, a), sp.Integer(0))
            self.assertAlmostEqual(float((got - source).subs(at)), 0.0, places=12, msg=a)

    def test_in_dickes_units_the_metric_is_fishers(self):
        """phi/phi_0 times the published spherical metric is the published spherical metric of
        fisher_jnw with b = 2 r_0 and gamma = 1 + (m - n)/2, component by component."""
        sp = self.sp
        reader, mine = self.components("spherical", "metric_components")
        theirs, fisher = self.components("spherical", "metric_components", metric_id="fisher_jnw")
        r, r0, m, n = reader.symbol["r"], reader.parameters["r_0"], reader.parameters["m"], reader.parameters["n"]
        field = (1 - 2 * r0 / r) ** (-(m + n) / 2)
        names = {theirs.symbol["r"]: r, theirs.symbol["\\theta"]: reader.symbol["\\theta"],
                 theirs.parameters["b"]: 2 * r0, theirs.parameters["gamma"]: 1 + (m - n) / 2}
        at = {r0: 1, m: sp.Rational(2, 7), n: sp.Rational(1, 3), r: sp.Rational(17, 5),
              reader.symbol["\\theta"]: sp.Rational(6, 5)}
        self.assertEqual(sorted(mine), sorted(fisher))
        for slot, value in mine.items():
            ratio = (field * value / fisher[slot].subs(names, simultaneous=True)).subs(at)
            self.assertAlmostEqual(float(ratio), 1.0, places=12, msg=str(slot))

    def test_far_away_the_isotropic_chart_is_the_post_newtonian_metric_at_gamma_equal_c_plus_one(self):
        """-g_tt = 1 - 2M/rho and g_rho rho = 1 + 2 gamma M/rho with M = 2B/lambda and gamma = C + 1."""
        sp = self.sp
        reader, g = self.components("isotropic", "metric_components")
        rho, B = reader.symbol["\\rho"], reader.parameters["B"]
        c, lam = reader.parameters["C"], reader.parameters["lambda"]
        e = sp.Symbol("e", positive=True)
        values = {c: -sp.Rational(1, 8), lam: sp.Rational(9, 10)}
        mass = 2 / values[lam]
        lapse = sp.series((-g[("t", "t")]).subs(values).subs(B, e * rho), e, 0, 2).removeO()
        space = sp.series(g[("\\rho", "\\rho")].subs(values).subs(B, e * rho), e, 0, 2).removeO()
        self.assertEqual(sp.simplify(sp.re(lapse) - (1 - 2 * mass * e)), 0)
        self.assertEqual(sp.simplify(sp.re(space) - (1 + 2 * (values[c] + 1) * mass * e)), 0)

    def test_without_the_field_each_chart_is_a_vacuum(self):
        """C = 0 with lambda = 1, m = n = 0, and s = 0 with b = k: every Ricci component vanishes."""
        sp = self.sp
        for system, vacuum in (("isotropic", {"C": 0, "lambda": 1}), ("spherical", {"m": 0, "n": 0}),
                               ("harmonic", {"s": 0, "b": "k"})):
            reader, ricci = self.components(system, "ricci_tensor", "ll")
            P = reader.parameters
            at = {P[k]: (P[v] if isinstance(v, str) else v) for k, v in vacuum.items()}
            for slot, value in ricci.items():
                self.assertEqual(sp.simplify(value.subs(at)), 0, f"{system} {slot}")

    def test_the_kretschmann_scalar_diverges_on_the_singularity_as_the_derivation_says(self):
        """K (r - 2 r_0)^(2 + 2n) has a finite limit that is not zero, 81 sqrt(2)/4096 at m = 0, n = 1/4."""
        sp = self.sp
        chart = self.charts["spherical"]
        reader = self.vm.Reader(chart["coords"], [p["symbol"] for p in chart["parameters"]], ())
        K = reader(chart["kretschmann"].partition("=")[2])
        r, P = reader.symbol["r"], reader.parameters
        K = K.subs({P["r_0"]: 1, P["m"]: 0, P["n"]: sp.Rational(1, 4)})
        x = sp.Symbol("x", positive=True)
        scaled = sp.simplify(K.subs(r, 2 + x) * x ** sp.Rational(5, 2))
        self.assertEqual(sp.simplify(sp.limit(scaled, x, 0, "+") - 81 * sp.sqrt(2) / 4096), 0)


if __name__ == "__main__":
    unittest.main()
