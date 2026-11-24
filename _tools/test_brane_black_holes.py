"""The black holes on a wall of Emparan, Horowitz and Myers, hep-th/9911043 and hep-th/9912135:
.venv.noindex/bin/python -m unittest discover -s _tools

What its charts, texts and drawings state, held on the published files. Each chart's published metric
and Ricci tensor are read through the checker's own Reader and held at points to the field equations
R_ab = -(3/l_4^2) g_ab of the sources; the flat wall's chart to being the sequel's at lambda = 0 and
k = -1, to the metric on the wall, the first paper's \\bhbrane, and its walls to the extrinsic curvature
of \\extcur; the spinning hole on a flat wall to the equator of the published Kerr black hole; the static
charts' Weyl tensor to Petrov type D; the mass of the cone to the thermodynamic mass, the first paper's
\\tmass, with its first law and temperature; the History to stating its formulas; the references to first
citation and the bibliography; the related spacetimes to answering; and the drawings at mu = 1/8 to the
surfaces, curves and moments they draw. The tests that read the mathematics need sympy and are skipped
where it is absent.
"""
import importlib.util
import json
import math
import re
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "MFS" / "assets" / "data"
HAS_SYMPY = importlib.util.find_spec("sympy") is not None
sys.path.insert(0, str(Path(__file__).resolve().parent / "derivations"))
ME = "brane_black_holes"
METRIC = json.loads((DATA / "metrics" / f"{ME}.json").read_text(encoding="utf-8"))
SYSTEMS = ("flat_wall", "two_walls", "rotating")
POINTS = ((0.2, -0.75), (1 / 3, -5.0), (1 / 7, -1 / 9))


def chart(system, metric=METRIC):
    return next(c for c in metric["coordinates"] if c["id"] == system)


def drawing(kind):
    return json.loads((DATA / kind / f"{ME}.json").read_text(encoding="utf-8"))


def cited(text):
    keys = []
    for group in re.findall(r"\[([a-z][^\]]*)\]", text):
        for key in (k.strip() for k in group.split(",")):
            if key not in keys:
                keys.append(key)
    return keys


def positive_root(f, lo, hi):
    """The root of f between lo and hi, where f changes sign, by bisection to the last digit."""
    a, b = lo, hi
    for _ in range(200):
        m = 0.5 * (a + b)
        if (f(a) > 0) == (f(m) > 0):
            a = m
        else:
            b = m
    return 0.5 * (a + b)


@unittest.skipUnless(HAS_SYMPY, "needs sympy")
class Charts(unittest.TestCase):
    """Each chart at points against the sources' field equations and statements."""

    @classmethod
    def setUpClass(cls):
        import sympy as sp
        import verify_metrics as vm
        cls.sp, cls.vm = sp, vm

    def published(self, system, metric=METRIC):
        """The chart's reader and its published metric, Ricci tensor and Weyl tensor as sympy."""
        sp = self.sp
        c = chart(system, metric)
        reader = self.vm.Reader(c["coords"], [p["symbol"] for p in c["parameters"]], ())
        n = len(c["coords"])
        index = {name: i for i, name in enumerate(c["coords"])}
        g = sp.zeros(n, n)
        for e in c["metric_components"]:
            g[index[e["indices"][0]], index[e["indices"][1]]] = reader(e["value"])
        ricci = sp.zeros(n, n)
        for e in c["ricci_tensor"]["variants"]["ll"]["nonzero"]:
            ricci[index[e["indices"][0]], index[e["indices"][1]]] = reader(e["value"])
        weyl = {}
        for e in c["weyl_tensor"]["variants"].get("llll", {"nonzero": []})["nonzero"]:
            weyl[tuple(index[i] for i in e["indices"])] = reader(e["value"])
        return c, reader, g, ricci, weyl

    def values(self, reader, **numbers):
        sp = self.sp
        out = {}
        for name, value in numbers.items():
            symbol = reader.symbol.get(name) if name in ("x", "y", "t") else reader.parameters[name]
            out[symbol] = sp.nsimplify(value)
        return out

    def test_each_chart_is_an_einstein_space_of_the_sources_radius(self):
        sp = self.sp
        cases = {"flat_wall": [{"ell": sp.Rational(5, 2), "mu": sp.Rational(1, 8)}],
                 "two_walls": [{"A": sp.Rational(3, 2), "m": sp.Rational(1, 5), "lambda": sp.Rational(1, 3), "k": k}
                               for k in (-1, 0, 1)],
                 "rotating": [{"A": sp.Rational(3, 2), "m": sp.Rational(1, 5), "lambda": sp.Rational(1, 3), "k": k,
                               "a": sp.Rational(2, 7)} for k in (-1, 0, 1)]}
        for system in SYSTEMS:
            _, R, g, ricci, _ = self.published(system)
            P = R.parameters
            inverse_square = 1 / P["ell"] ** 2 if system == "flat_wall" else P["A"] ** 2 * (P["lambda"] + 1)
            for case in cases[system]:
                for x, y in POINTS:
                    at = self.values(R, x=x, y=y, **case)
                    for i in range(4):
                        for j in range(i, 4):
                            with self.subTest(system=system, case=case, point=(x, y), slot=(i, j)):
                                self.assertEqual(sp.simplify((ricci[i, j] + 3 * inverse_square * g[i, j]).subs(at)), 0)

    def test_the_first_paper_s_chart_is_the_sequel_s_at_lambda_zero(self):
        sp = self.sp
        _, F, gf, _, _ = self.published("flat_wall")
        _, W, gw, _, _ = self.published("two_walls")
        ell, mu = sp.Rational(5, 2), sp.Rational(1, 8)
        for x, y in POINTS:
            here = self.values(F, x=x, y=y, ell=ell, mu=mu)
            there = self.values(W, x=x, y=y, A=1 / ell, m=mu * ell, k=-1)
            there[W.parameters["lambda"]] = 0
            for i in range(4):
                with self.subTest(point=(x, y), slot=i):
                    self.assertEqual(sp.simplify(gf[i, i].subs(here) - gw[i, i].subs(there)), 0)

    def test_the_metric_on_the_flat_wall_is_schwarzschild_s_plane_and_a_cone(self):
        sp = self.sp
        _, R, g, _, _ = self.published("flat_wall")
        ell, mu, x, y = R.parameters["ell"], R.parameters["mu"], R.symbol["x"], R.symbol["y"]
        r = sp.Symbol("r", positive=True)
        on = {x: 0, y: -ell / r}
        f = 1 - 2 * mu * ell / r
        self.assertEqual(sp.simplify(g[0, 0].subs(on) / ell ** 2 + f), 0)
        self.assertEqual(sp.simplify(g[1, 1].subs(on) * (ell / r ** 2) ** 2 - 1 / f), 0)
        self.assertEqual(sp.simplify(g[3, 3].subs(on) - r ** 2), 0)
        # Far out the circle at r has the circumference r Delta varphi, short of 2 pi r by the deficit angle.
        self.assertIn("$\\Delta\\varphi = 2\\pi/(x_2(1 + 3\\mu x_2))$", chart("flat_wall")["convention"])

    def test_each_wall_bends_by_its_metric_times_a_positive_number(self):
        sp = self.sp
        _, F, gf, _, _ = self.published("flat_wall")
        ell, x, y = F.parameters["ell"], F.symbol["x"], F.symbol["y"]
        G = 1 - x ** 2 - 2 * F.parameters["mu"] * x ** 3
        for i in (0, 1, 3):
            # The first paper's \extcur, K = -((x - y) sqrt(G)/2 l) dg/dx at x = 0, is g/l there.
            K = -((x - y) * sp.sqrt(G) / (2 * ell)) * sp.diff(gf[i, i], x)
            self.assertEqual(sp.simplify((K - gf[i, i] / ell).subs(x, 0)), 0)
        _, W, gw, _, _ = self.published("two_walls")
        A, lam, x, y = W.parameters["A"], W.parameters["lambda"], W.symbol["x"], W.symbol["y"]
        k, m = W.parameters["k"], W.parameters["m"]
        G = 1 + k * x ** 2 - 2 * m * A * x ** 3
        H = -lam + k * y ** 2 - 2 * m * A * y ** 3
        for i in (0, 1, 3):
            # The outward normal of the kept quadrant points to smaller x on the wall x = 0.
            K = -A * (x - y) * sp.sqrt(G) / 2 * sp.diff(gw[i, i], x)
            self.assertEqual(sp.simplify((K - A * gw[i, i]).subs(x, 0)), 0)
        for i in (0, 2, 3):
            # and to larger y on the wall y = 0, where g^yy = A^2 x^2 lambda.
            K = A * (x - y) * sp.sqrt(-H) / 2 * sp.diff(gw[i, i], y)
            self.assertEqual(sp.simplify((K - A * sp.sqrt(lam) * gw[i, i]).subs(y, 0)), 0)

    def test_the_spinning_hole_on_a_flat_wall_is_kerr_s_equator(self):
        sp = self.sp
        _, R, g, _, _ = self.published("rotating")
        kerr = json.loads((DATA / "metrics" / "kerr.json").read_text(encoding="utf-8"))
        _, K, gk, _, _ = self.published("boyer_lindquist", kerr)
        A, m, a = sp.Rational(3, 2), sp.Rational(1, 9), sp.Rational(1, 10)
        here = {R.parameters["A"]: A, R.parameters["m"]: m, R.parameters["lambda"]: 0, R.parameters["k"]: -1,
                R.parameters["a"]: a, R.symbol["x"]: 0}
        that, rho = sp.symbols("that rho", positive=True)
        phi = R.symbol["\\varphi"]
        # The sequel's \branegeom: t-hat = (t + a varphi)/A and rho-hat = rho/A = -1/(A y), so the chart's
        # t, y and varphi are these functions of t-hat, rho-hat and varphi.
        old = [A * that - a * phi, -1 / (A * rho), phi]
        new = [that, rho, phi]
        J = sp.Matrix(3, 3, lambda i, j: sp.diff(old[i], new[j]))
        block = g.extract([0, 1, 3], [0, 1, 3]).subs(here).subs(R.symbol["y"], old[1])
        pulled = J.T * block * J
        there = {K.parameters["M"]: m, K.parameters["a"]: a / A, K.parameters["G"]: 1, K.c: 1,
                 K.symbol["\\theta"]: sp.pi / 2, K.symbol["r"]: rho}
        kblock = gk.extract([0, 1, 3], [0, 1, 3]).subs(there)
        for rv in (sp.Rational(1, 2), 2, 7):
            for i in range(3):
                for j in range(3):
                    with self.subTest(rho=rv, slot=(i, j)):
                        self.assertEqual(sp.simplify((pulled[i, j] - kblock[i, j]).subs(rho, rv)), 0)

    def test_the_static_charts_are_of_petrov_type_d(self):
        sp = self.sp
        for system, case in (("flat_wall", {"ell": sp.Rational(5, 2), "mu": sp.Rational(1, 8)}),
                             ("two_walls", {"A": sp.Rational(3, 2), "m": sp.Rational(1, 5), "lambda": sp.Rational(1, 3),
                                            "k": 1})):
            _, R, g, _, weyl = self.published(system)
            P = R.parameters
            for x, y in ((sp.Rational(1, 5), -sp.Rational(3, 4)), (sp.Rational(1, 7), -sp.Rational(1, 9))):
                at = self.values(R, x=x, y=y, **case)
                gn = g.subs(at)
                if not (gn[0, 0] < 0 < gn[1, 1]):
                    continue
                e = [sp.zeros(4, 1) for _ in range(4)]
                for i in range(4):
                    e[i][i] = 1 / sp.sqrt(abs(gn[i, i]))
                l, n = (e[0] + e[1]) / sp.sqrt(2), (e[0] - e[1]) / sp.sqrt(2)
                mm, mb = (e[2] + sp.I * e[3]) / sp.sqrt(2), (e[2] - sp.I * e[3]) / sp.sqrt(2)
                C = {key: value.subs(at) for key, value in weyl.items()}

                def psi(u, v, w, z):
                    return sp.nsimplify(sp.expand(sum(c * u[a] * v[b] * w[d] * z[f] for (a, b, d, f), c in C.items())))
                psis = [psi(l, mm, l, mm), psi(l, n, l, mm), psi(l, mm, mb, n), psi(l, n, mb, n), psi(n, mb, n, mb)]
                if system == "flat_wall":
                    want = -case["mu"] * (x - y) ** 3 / case["ell"] ** 2
                else:
                    want = -case["m"] * case["A"] ** 3 * (x - y) ** 3
                with self.subTest(system=system, point=(x, y)):
                    self.assertEqual([sp.simplify(p) for p in psis[:2] + psis[3:]], [0, 0, 0, 0])
                    self.assertEqual(sp.simplify(psis[2] - want), 0)
        self.assertIn("Petrov type D", METRIC["tags"])

    def test_the_axis_of_the_two_wall_chart_holds_its_null_curves(self):
        """On the axis x = x_2 no Christoffel symbol turns a curve of the plane of t and y out of it."""
        c = chart("two_walls")
        R = self.vm.Reader(c["coords"], [p["symbol"] for p in c["parameters"]], ())
        sp = self.sp
        at = {R.parameters["A"]: 1, R.parameters["m"]: 1, R.parameters["lambda"]: sp.Rational(1, 2),
              R.parameters["k"]: 1, R.symbol["x"]: 1}
        for e in c["christoffel"]["variants"]["ull"]["nonzero"]:
            up, a, b = e["indices"]
            if up in ("x", "\\varphi") and a in ("t", "y") and b in ("t", "y"):
                with self.subTest(indices=e["indices"]):
                    self.assertEqual(sp.simplify(R(e["value"]).subs(at)), 0)


class Thermodynamics(unittest.TestCase):
    """The first paper's masses: the mass of the cone, 8 pi G_3 M_3 = 2 pi - Delta varphi with
    G_3 = G_4/2l, and the thermodynamic mass \\tmass, in units l = G_4 = 1."""

    @staticmethod
    def hole(mu):
        x2 = positive_root(lambda x: 1 - x * x - 2 * mu * x ** 3, 0.0, 1.0)
        fold = 1 / (x2 * (1 + 3 * mu * x2))
        xh = 2 * mu * x2
        return x2, fold, xh

    def test_the_mass_of_the_cone_equals_the_thermodynamic_mass(self):
        for mu in (1e-3, 0.05, 0.125, 1 / (3 * math.sqrt(3)), 1.0, 10.0, 1e4):
            x2, fold, xh = self.hole(mu)
            m3 = (1 - fold) / 2                     # G_4 M_3 / l, from 8 pi (G_4/2l) M_3 = 2 pi (1 - fold)
            m4 = (1 - math.sqrt(1 + xh) / (1 + 1.5 * xh)) / 2
            with self.subTest(mu=mu):
                self.assertAlmostEqual(xh * xh * (1 + xh), 4 * mu * mu, delta=1e-12 * max(1, 4 * mu * mu))
                self.assertAlmostEqual(m3, m4, delta=1e-12)

    def test_the_first_law_holds_with_the_temperature_on_the_wall(self):
        for mu in (0.01, 0.125, 2.0, 50.0):
            def entropy(m):
                _, _, xh = self.hole(m)
                return math.pi * xh * xh / (1 + 1.5 * xh)

            def mass(m):
                _, fold, _ = self.hole(m)
                return (1 - fold) / 2
            h = 1e-5 * mu
            temperature = (mass(mu + h) - mass(mu - h)) / (entropy(mu + h) - entropy(mu - h))
            _, _, xh = self.hole(mu)
            with self.subTest(mu=mu):
                self.assertAlmostEqual(temperature * 8 * math.pi * mu, 1, delta=1e-6)
                self.assertAlmostEqual(1 / (4 * math.pi * xh * math.sqrt(1 + xh)), 1 / (8 * math.pi * mu), delta=1e-12)

    def test_the_mass_stays_below_a_quarter_while_the_area_grows_without_bound(self):
        last_m3, last_area = 0.0, 0.0
        for mu in (0.01, 0.1, 1.0, 10.0, 100.0, 1e3, 1e4, 1e5):
            x2, fold, _ = self.hole(mu)
            g3m3 = (1 - fold) / 4                 # G_3 M_3, in G_3 = G_4/2l = 1/2
            y0 = -1 / (2 * mu)
            area = 2 * (2 * math.pi * fold) * (1 / -y0 - 1 / (x2 - y0))   # the first paper's \bharea, l = 1
            with self.subTest(mu=mu):
                self.assertLess(g3m3, 0.25)
                self.assertGreater(g3m3, last_m3)
                self.assertGreater(area, last_area)
            last_m3, last_area = g3m3, area
        self.assertGreater(last_area, 1e3)
        self.assertAlmostEqual(last_area / ((8 * math.pi / 3) * (2e5) ** (2 / 3)), 1, delta=1e-3)


class Texts(unittest.TestCase):
    """The History states the sources' formulas, and every citation is held."""

    def test_history_states_the_solution(self):
        history = METRIC["history"]
        for formula in ("$1/(2\\pi G_4\\ell)$", "$G_3 = G_4/(2\\ell)$", "$M_3 = M_4$", "$1/(4G_3)$",
                        "$\\ell\\log(2\\mu)$", "$2\\mu\\ell$", "$x = 0$", "$y = 0$", "$M_3 \\approx 1/(32G_3)$"):
            with self.subTest(formula=formula):
                self.assertIn(formula, history)
        self.assertEqual(len(history.split("¶")), 6)

    def test_the_conventions_define_the_period_and_the_wall_s_radius(self):
        flat = chart("flat_wall")["convention"]
        self.assertIn("$r = -\\ell/y$", flat)
        self.assertIn("$\\hat t = \\ell t$", flat)
        self.assertIn("$\\Delta\\varphi = 2\\pi/(x_2(3mAx_2 - k))$", chart("two_walls")["convention"])
        self.assertIn("$t + ax_2^2\\varphi$", chart("rotating")["convention"])
        self.assertNotIn("carrying", METRIC["convention"] + "".join(c["convention"] for c in METRIC["coordinates"]))

    def test_references_follow_first_citation_and_are_in_the_bibliography(self):
        keys = cited(METRIC["history"])
        self.assertEqual(METRIC["references"][:len(keys)], keys)
        for entry in METRIC["related"]:
            for key in cited(entry["text"]):
                self.assertIn(key, METRIC["references"])
        bib = (ROOT / "assets" / "data" / "references.bib").read_text(encoding="utf-8")
        for key in METRIC["references"]:
            self.assertEqual(len(re.findall(r"^@\w+\{" + re.escape(key) + ",", bib, re.M)), 1, key)

    def test_every_related_spacetime_answers(self):
        answers = {"piece": "composite", "family": "family"}
        self.assertEqual(len(METRIC["related"]), 10)
        for entry in METRIC["related"]:
            other = json.loads((DATA / "metrics" / f"{entry['id']}.json").read_text(encoding="utf-8"))
            answer = [e for e in other["related"] if e["id"] == ME]
            with self.subTest(other=entry["id"]):
                self.assertEqual(entry["kind"], "piece" if entry["id"] == "anti_de_sitter" else "family")
                self.assertEqual(len(answer), 1)
                self.assertEqual(answer[0]["kind"], answers[entry["kind"]])
                for key in cited(answer[0]["text"]):
                    self.assertIn(key, other["references"])

    def test_the_tags_hold_what_the_charts_print(self):
        for tag in ("rotating black hole", "angular momentum", "Einstein space", "negative cosmological constant",
                    "Petrov type D", "braneworld", "four-dimensional"):
            with self.subTest(tag=tag):
                self.assertIn(tag, METRIC["tags"])
        self.assertEqual([c["id"] for c in METRIC["coordinates"]], list(SYSTEMS))


class Drawings(unittest.TestCase):
    """The drawings at mu = 1/8 in units of l, whose horizon meets the wall at r = l/4."""

    S = None

    @classmethod
    def setUpClass(cls):
        x2 = positive_root(lambda x: 1 - x * x - x ** 3 / 4, 0.0, 1.0)
        cls.S = 1 / (x2 * (1 + 3 * x2 / 8))

    def test_embedding_is_the_wall_from_the_horizon_opening_into_a_cone(self):
        (view,) = drawing("embedding")["views"]
        self.assertEqual(view["id"], "wall")
        (piece,) = view["surfaces"][0]["pieces"]
        points = piece["points"]
        for y, rho, _ in points:
            self.assertAlmostEqual(rho, self.S / abs(y), places=5)
        self.assertAlmostEqual(points[0][0], -4.0, places=12)
        self.assertAlmostEqual(points[-1][0], -1 / 3, places=12)
        # The slope dz/drho = sqrt(1/(1 - r_s/r) - s^2)/s steepens toward the horizon and tends far out to the
        # cone's sqrt(1 - s^2)/s, which it passes above: at r = 3 l it stands within 15% of the cone's.
        (y1, r1, z1), (y2, r2, z2) = points[-2], points[-1]
        r = -2 / (y1 + y2)
        self.assertAlmostEqual((z2 - z1) / (r2 - r1), math.sqrt(1 / (1 - 0.25 / r) - self.S ** 2) / self.S, delta=0.01)
        slopes = [(b[2] - a[2]) / (b[1] - a[1]) for a, b in zip(points, points[1:])]
        self.assertTrue(all(u > v for u, v in zip(slopes, slopes[1:])))
        cone = math.sqrt(1 - self.S ** 2) / self.S
        self.assertLess(slopes[-1] / cone - 1, 0.15)
        self.assertIn(f"{360 * (1 - self.S):.0f}°", view["caption"][0])

    def test_conformal_view_is_the_exterior_and_the_black_hole_with_the_moment(self):
        (view,) = drawing("conformal")["views"]
        self.assertEqual(view["id"], "flat_wall")
        classes = {layer["class"] for layer in view["layers"]}
        for cls in ("horizon", "singular", "scri", "infinity", "cover"):
            self.assertIn(cls, classes)
        self.assertEqual([s["view"] for s in view["slices"]], ["wall"])

    def test_curvature_figure_is_schwarzschild_s_plus_the_constant_of_anti_de_sitter_space(self):
        (view,) = drawing("curvature")["views"]
        self.assertEqual([p["id"] for p in view["panels"]], ["K"])
        (own,) = [c for c in view["curves"] if c["class"] == "K"]
        for (r, K), y in zip(own["points"], own["horizontal"]["values"]):
            self.assertAlmostEqual(y, -1 / r, places=12)
            self.assertAlmostEqual(K / (24 + 12 * 0.25 ** 2 / r ** 6), 1, places=9)
        self.assertTrue(all(c["metric"] == "schwarzschild" for c in view["curves"] if c["class"] == "reference"))
        self.assertEqual(view["marks"], [{"class": "horizon", "x": 0.25}])
        self.assertIn("The Ricci scalar is $-12/\\ell^2$.", view["caption"][0])

    def test_spacetime_diagrams_mark_the_moment_on_the_flat_wall_alone(self):
        systems = drawing("diagrams")["systems"]
        self.assertEqual(sorted(systems), sorted(SYSTEMS))
        (flat,) = systems["flat_wall"]
        (axis,) = systems["two_walls"]
        (spin,) = systems["rotating"]
        self.assertEqual(flat["box"], [-6, 0, -3, 3])
        self.assertEqual([s["view"] for s in flat["slices"]], ["wall"])
        self.assertNotIn("slices", axis)
        self.assertNotIn("slices", spin)
        hidden = json.loads((ROOT / "_tools" / "derivations" / "exceptions" / f"{ME}.json").read_text(encoding="utf-8"))
        self.assertEqual(sorted(hidden["hidden"][0]["drawings"]), [f"{ME}/rotating/wall", f"{ME}/two_walls/axis"])


if __name__ == "__main__":
    unittest.main()
