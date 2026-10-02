"""Einstein's static field of 1912, flat space with a speed of light c N(x, y, z):
python3 -m unittest discover -s _tools

What the texts and drawings of `einstein_1912_static` state, held on the published files. Every
chart is held to flat space and one lapse, the free chart to R = -2 Laplacian(N)/N with no tt
component of the Einstein tensor, and each other chart to being that chart at its own N: the
uniform field flat and the same chart as Minkowski's published Rindler chart, the body of
February a solution of Laplace's equation for N, and the body of March a solution of Einstein's
second equation, Laplace's for sqrt(N). Both bodies are held to gamma = 0, to beta = 1/2 and 3/4,
and to the published metric of `ppn_metric` at those numbers. A ray and a planet are run with
the published Christoffel symbols: the ray bends by half of general relativity's angle past
either body, and the planet's perihelion advances by 1/2 of Einstein's advance in February and
5/12 in March. A body dropped toward the sphere where N vanishes reaches it in a finite proper
time. The declared star is held to Giulini's interior solution and his bound. The spacetime
diagrams are held to cones of slope N, the conformal diagrams to their singular null edges, and
the embedding diagram to a flat plane under Flamm's paraboloid. The tests that read the
published components need sympy and are skipped where it is absent.
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
HAS_NUMPY = importlib.util.find_spec("numpy") is not None and importlib.util.find_spec("scipy") is not None
sys.path.insert(0, str(Path(__file__).resolve().parent / "derivations"))

CHARTS = ("static", "uniform", "february", "march")
LAPSE = {"february": lambda r: 1 - 1 / r, "march": lambda r: (1 - 1 / (2 * r)) ** 2}
EDGE = {"february": 1.0, "march": 0.5}


def published(system, metric_id="einstein_1912_static"):
    metric = json.loads((DATA / "metrics" / f"{metric_id}.json").read_text(encoding="utf-8"))
    return next(c for c in metric["coordinates"] if c["id"] == system)


def drawn(kind):
    return json.loads((DATA / kind / "einstein_1912_static.json").read_text(encoding="utf-8"))


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

    def lapse(self, system):
        sp = self.sp
        if system == "static":
            return self.P["static"]["N"]
        if system == "uniform":
            return self.P["uniform"]["a"] * self.x["uniform"][3]
        r, m = self.x[system][1], self.P[system]["m"]
        return 1 - m / r if system == "february" else (1 - m / (2 * r)) ** 2

    def test_every_chart_is_flat_space_with_one_lapse(self):
        sp = self.sp
        for system in CHARTS:
            g, x = self.g[system], self.x[system]
            space = sp.diag(1, 1, 1) if system in ("static", "uniform") else sp.diag(1, x[1] ** 2, x[1] ** 2 * sp.sin(x[2]) ** 2)
            self.assertEqual((g[1:, 1:] - space).applyfunc(sp.simplify), sp.zeros(3, 3), system)
            self.assertEqual(sp.simplify(g[0, 0] + self.lapse(system) ** 2), 0, system)
            self.assertEqual([g[0, k] for k in (1, 2, 3)], [0, 0, 0], system)

    def test_flat_space_leaves_the_einstein_tensor_without_a_tt_component(self):
        for system in CHARTS:
            self.assertNotIn((0, 0), self.tensor(system, "einstein_tensor", "ll"), system)

    def test_the_free_charts_curvature_scalar_is_minus_twice_the_laplacian_of_n_over_n(self):
        sp = self.sp
        t, x, y, z = self.x["static"]
        N = self.P["static"]["N"]
        laplacian = sum(sp.diff(N, s, 2) for s in (x, y, z))
        self.assertEqual(sp.simplify(self.scalar("static") + 2 * laplacian / N), 0)
        self.assertEqual(sp.simplify(self.tensor("static", "ricci_tensor", "ll")[(0, 0)] - N * laplacian), 0)

    def test_the_uniform_field_is_flat_and_is_minkowskis_rindler_chart(self):
        sp, vm = self.sp, self.vm
        for field in ("riemann", "ricci_tensor", "einstein_tensor", "weyl_tensor"):
            variants = self.entry["uniform"][field]["variants"]
            self.assertEqual([v["nonzero"] for v in variants.values()], [[]] * len(variants), field)
        self.assertEqual(self.entry["uniform"]["kretschmann"], "K = 0")
        there = published("rindler", "minkowski")
        reader = vm.Reader(there["coords"], [p["symbol"] for p in there["parameters"]], ())
        rindler = {tuple(e["indices"]): reader(e["value"]).subs(reader.c, 1) for e in there["metric_components"]}
        at = {reader.symbol["X"]: self.x["uniform"][3], reader.parameters["a"]: self.P["uniform"]["a"]}
        self.assertEqual(sp.simplify(rindler[("T", "T")].subs(at) - self.g["uniform"][0, 0]), 0)
        self.assertEqual([rindler[(s, s)] for s in "XYZ"], [1, 1, 1])

    def test_the_body_of_february_solves_laplaces_equation_and_the_body_of_march_einsteins_second(self):
        sp = self.sp
        t, x, y, z = self.x["static"]
        N = self.P["static"]["N"]
        m = sp.Symbol("m", positive=True)
        r = sp.sqrt(x ** 2 + y ** 2 + z ** 2)
        laplacian = lambda f: sum(sp.diff(f, s, 2) for s in (x, y, z))
        gradient = lambda f: sum(sp.diff(f, s) ** 2 for s in (x, y, z))
        february, march = 1 - m / r, (1 - m / (2 * r)) ** 2
        self.assertEqual(sp.simplify(laplacian(february)), 0)
        # N Laplacian(N) - (grad N)^2/2 = 0 in a vacuum, which is Laplace's equation for sqrt(N).
        self.assertEqual(sp.simplify(march * laplacian(march) - gradient(march) / 2), 0)
        self.assertEqual(sp.simplify(laplacian(1 - m / (2 * r))), 0)
        self.assertNotEqual(sp.simplify(laplacian(march)), 0)
        # The published scalars are the free chart's at those N.
        general = self.scalar("static")
        self.assertEqual(sp.simplify(general.subs(N, february).doit()), 0)
        self.assertEqual(self.scalar("february"), 0)
        rr, mm = self.x["march"][1], self.P["march"]["m"]
        far = sp.Symbol("far", positive=True)     # a point of the x axis, at the distance far
        on_axis = sp.simplify(general.subs(N, march).doit().subs({x: far, y: 0, z: 0}))
        self.assertEqual(sp.simplify(on_axis - self.scalar("march").subs({rr: far, mm: m})), 0)
        for system in ("february", "march"):
            self.assertTrue(self.tensor(system, "riemann", "ulll"), system)

    def test_the_uniform_field_solves_the_equation_of_february_and_fails_that_of_march(self):
        sp = self.sp
        z, a = self.x["uniform"][3], self.P["uniform"]["a"]
        N = a * z
        self.assertEqual(sp.diff(N, z, 2), 0)
        self.assertEqual(sp.simplify(N * sp.diff(N, z, 2) - sp.diff(N, z) ** 2 / 2), -a ** 2 / 2)

    def test_both_bodies_have_gamma_zero_and_beta_one_half_and_three_quarters(self):
        sp = self.sp
        U = sp.Symbol("U", positive=True)
        for system, beta in (("february", sp.Rational(1, 2)), ("march", sp.Rational(3, 4))):
            r, m, g = self.x[system][1], self.P[system]["m"], self.g[system]
            # -g_tt = 1 - 2U + 2 beta U^2 + ..., U = m/r, and g_rr = 1 + 2 gamma U = 1.
            series = sp.series((-g[0, 0]).subs(r, m / U), U, 0, 3).removeO()
            self.assertEqual(sp.expand(series - (1 - 2 * U + 2 * beta * U ** 2)), 0, system)
            self.assertEqual(g[1, 1], 1, system)
            # Light bends by (1 + gamma)/2 of general relativity's angle, a perihelion advances by
            # (2 + 2 gamma - beta)/3 of its advance.
            self.assertEqual(sp.Rational(1 + 0, 2), sp.Rational(1, 2))
            self.assertEqual((2 - beta) / 3, {"february": sp.Rational(1, 2), "march": sp.Rational(5, 12)}[system])

    def test_they_are_the_post_newtonian_metric_at_those_numbers(self):
        sp, vm = self.sp, self.vm
        there = published("isotropic", "ppn_metric")
        reader = vm.Reader(there["coords"], [p["symbol"] for p in there["parameters"]], ())
        ppn = {tuple(e["indices"]): reader(e["value"]) for e in there["metric_components"]}
        U = sp.Symbol("U", positive=True)
        for system, beta in (("february", sp.Rational(1, 2)), ("march", sp.Rational(3, 4))):
            m, r = self.P[system]["m"], self.x[system][1]
            at = {reader.parameters["beta"]: beta, reader.parameters["gamma"]: 0, reader.parameters["m"]: m,
                  reader.symbol["r"]: r}
            gap = sp.series((ppn[("t", "t")].subs(at) - self.g[system][0, 0]).subs(r, m / U), U, 0, 3).removeO()
            self.assertEqual(sp.simplify(gap), 0, system)
            self.assertEqual(sp.simplify(ppn[("r", "r")].subs(at) - self.g[system][1, 1]), 0, system)

    def test_the_kretschmann_scalar_diverges_where_the_lapse_vanishes(self):
        sp = self.sp
        for system, edge, power, limit in (("february", 1, 2, 24), ("march", sp.Rational(1, 2), 4, 16)):
            r, m = self.x[system][1], self.P[system]["m"]
            text = self.entry[system]["kretschmann"].partition("=")[2]
            K = self.reader[system](text).subs(self.reader[system].c, 1)
            self.assertEqual(sp.limit(K.subs(m, 1) * (r - edge) ** power, r, edge), limit, system)
            self.assertEqual(sp.simplify(self.lapse(system).subs({m: 1, r: edge})), 0, system)

    def runner(self, system, values):
        """The geodesic equations of the published Christoffel symbols, as a function of the state
        (x^mu, dx^mu/d lambda), at the parameter values given."""
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

    def test_a_ray_bends_by_half_of_general_relativitys_angle(self):
        """A ray run with the published Christoffel symbols at m = 1 from its nearest point r_0 out
        to 2000 r_0, where what is left of its turn is that of a straight line, arcsin(b/r) with b
        its impact parameter r_0/N(r_0). It turns through pi/2 + m/b on the way out, so the whole
        deflection is 2m/b, half of general relativity's 4m/b."""
        for system in ("february", "march"):
            step = self.runner(system, {"m": 1})
            r0 = 500.0
            lapse = LAPSE[system](r0)
            b = r0 / lapse
            # At the nearest point the ray moves along phi: N dt = r dphi.
            y = [0.0, r0, math.pi / 2, 0.0, 1 / lapse, 0.0, 0.0, 1 / r0]
            path = rk4(step, y, 25.0, 60000, stop=lambda s: s[1] > 2000 * r0)
            t, r, _, phi, _, vr, _, vphi = path[-1]
            self.assertGreater(r, 2000 * r0)
            # Far out the ray is a straight line whose direction makes the angle arctan(r vphi/vr) with the radius.
            swept = phi + math.atan2(r * vphi, vr)
            deflection = 2 * swept - math.pi
            self.assertAlmostEqual(deflection / (4 / b), 0.5, delta=5e-3, msg=system)

    def test_a_planets_perihelion_advances_by_a_half_and_by_five_twelfths_of_einsteins_advance(self):
        """A planet run with the published Christoffel symbols at m = 1 from its perihelion round
        to the next: it overshoots a full turn by (2 - beta)/3 of general relativity's
        6 pi m^2/l^2, with l its angular momentum per unit mass."""
        for system, share in (("february", 1 / 2), ("march", 5 / 12)):
            step = self.runner(system, {"m": 1})
            r0 = 400.0
            lapse = LAPSE[system](r0)
            speed = 1.2 * math.sqrt(1 / r0)   # faster than circular, so r0 is the perihelion
            boost = 1 / math.sqrt(1 - speed ** 2)
            ut, uphi = boost / lapse, boost * speed / r0
            ell = r0 ** 2 * uphi
            y = [0.0, r0, math.pi / 2, 0.0, ut, 0.0, 0.0, uphi]
            path = rk4(step, y, 2.0, 400000, stop=lambda s: s[5] > 0 and s[1] < r0 * 1.02 and s[3] > 5.0)
            turn = None
            for a, b in zip(path, path[1:]):
                if a[3] > 3.0 and a[5] < 0 <= b[5]:
                    w = -a[5] / (b[5] - a[5])
                    turn = a[3] + w * (b[3] - a[3])
                    break
            self.assertIsNotNone(turn, system)
            shift, einstein = turn - 2 * math.pi, 6 * math.pi / ell ** 2
            self.assertGreater(shift, 0, system)
            self.assertAlmostEqual(shift / einstein, share, delta=6e-3, msg=system)

    def test_a_dropped_body_reaches_the_singular_sphere_in_a_finite_proper_time(self):
        """Dropped from rest at r = 3m, a body has N^2 dt/d tau = N(3m) and falls at
        dr/d tau = -sqrt(N(3m)^2/N^2 - 1), which grows without bound where N vanishes, so its
        clock reads a finite time on arrival while the time t runs on without end."""
        for system in ("february", "march"):
            step = self.runner(system, {"m": 1})
            edge, lapse = EDGE[system], LAPSE[system]
            y = [0.0, 3.0, math.pi / 2, 0.0, 1 / lapse(3.0), 0.0, 0.0, 0.0]
            h = 5e-5
            # Down to where the lapse is a twentieth, r = 1.053 m in February and 0.644 m in March.
            near = next(edge + k * 1e-4 for k in range(1, 20000) if lapse(edge + k * 1e-4) > 0.05)
            path = rk4(step, y, h, 400000, stop=lambda s: s[1] < near)
            t, r, _, _, ut, ur, _, _ = path[-1]
            tau = h * (len(path) - 1)
            self.assertLess(r, near, system)
            self.assertLess(tau, 12.0, system)
            # The energy N^2 dt/d tau is kept, and the fall is as fast as the lapse says.
            self.assertAlmostEqual(lapse(r) ** 2 * ut, lapse(3.0), places=4, msg=system)
            self.assertAlmostEqual(ur, -math.sqrt(lapse(3.0) ** 2 / lapse(r) ** 2 - 1), delta=1e-3 * abs(ur), msg=system)
            # What is left of the fall takes less proper time than the distance left: dr/d tau < -1.
            self.assertLess(ur, -10.0, system)
            self.assertGreater(t, 2 * tau, system)


@unittest.skipUnless(HAS_NUMPY, "numpy and scipy are not installed")
class Star(unittest.TestCase):
    """The declared star of uniform density, Giulini's interior solution at omega R = 3/2."""

    @classmethod
    def setUpClass(cls):
        import null_rays as nr
        import numpy as np
        cls.nr, cls.np = nr, np

    def root(self, r, order=0):
        return float(self.nr._e12_root(self.np.array([r * r]), order)[0])

    def test_it_is_giulinis_interior_solution_joined_to_his_exterior(self):
        w, rg = self.nr.E12_OMEGA, self.nr.E12_RG
        self.assertAlmostEqual(rg, 1 - math.tanh(w) / w, places=15)
        for r in (0.1, 0.5, 0.9):
            self.assertAlmostEqual(self.root(r), math.sinh(w * r) / (w * r * math.cosh(w)), places=13)
        for r in (1.5, 4.0, 30.0):
            self.assertAlmostEqual(self.root(r), 1 - rg / r, places=13)
        self.assertAlmostEqual(self.root(0.0), 1 / math.cosh(w), places=13)
        # sqrt(N) and its slope are continuous across the surface.
        eps = 1e-7
        self.assertAlmostEqual(self.root(1 - eps), self.root(1 + eps), places=6)
        self.assertAlmostEqual(self.root(1 - eps, 1), self.root(1 + eps, 1), places=5)

    def test_inside_it_solves_einsteins_second_equation_for_a_uniform_density(self):
        """Laplacian(sqrt N) = omega^2 sqrt(N): in q = r^2 the Laplacian of f(q) is 4 q f'' + 6 f'."""
        w = self.nr.E12_OMEGA
        for r in (0.0, 0.05, 0.3, 0.7, 0.95):
            laplacian = 4 * r * r * self.root(r, 2) + 6 * self.root(r, 1)
            self.assertAlmostEqual(laplacian, w * w * self.root(r), places=10)
        for r in (1.2, 3.0):
            self.assertAlmostEqual(4 * r * r * self.root(r, 2) + 6 * self.root(r, 1), 0.0, places=12)

    def test_its_lapse_and_two_derivatives_are_those_of_its_square_root(self):
        np = self.np
        for q in (0.0, 0.2, 0.8, 2.0):
            at = np.array([q])
            f, f1, f2 = (float(self.nr._e12_root(at, k)[0]) for k in (0, 1, 2))
            self.assertAlmostEqual(float(self.nr._e12_star(at, 0)[0]), f * f, places=14)
            self.assertAlmostEqual(float(self.nr._e12_star(at, 1)[0]), 2 * f * f1, places=14)
            self.assertAlmostEqual(float(self.nr._e12_star(at, 2)[0]), 2 * f1 * f1 + 2 * f * f2, places=14)
        h = 1e-5
        for q in (0.3, 1.7):
            slope = (float(self.nr._e12_star(np.array([q + h]), 0)[0]) - float(self.nr._e12_star(np.array([q - h]), 0)[0])) / (2 * h)
            self.assertAlmostEqual(slope, float(self.nr._e12_star(np.array([q]), 1)[0]), places=8)

    def test_its_surface_lies_outside_the_sphere_where_the_exterior_lapse_vanishes(self):
        """Giulini's theorem, R_g < R, so m = 2 R_g < 2R; and N is positive throughout."""
        rg = self.nr.E12_RG
        self.assertLess(rg, 1.0)
        self.assertAlmostEqual(2 * rg, 0.793, places=3)
        self.assertAlmostEqual(self.root(0.0) ** 2, 0.181, places=3)
        self.assertAlmostEqual(self.root(1.0) ** 2, 0.364, places=3)
        self.assertTrue(all(self.root(r) > 0 for r in (0.0, 0.25, 0.5, 0.75, 1.0, 2.0, 10.0)))

    def test_light_crosses_it_in_8_89_of_its_radius_over_c(self):
        self.assertAlmostEqual(2 * float(self.nr.e12_tortoise(1.0)), 8.89, delta=5e-3)
        # Outside the star the tortoise coordinate grows at dx_*/dx = 1/N.
        h = 1e-5
        slope = (float(self.nr.e12_tortoise(2.0 + h)) - float(self.nr.e12_tortoise(2.0 - h))) / (2 * h)
        self.assertAlmostEqual(slope, 1 / self.root(2.0) ** 2, places=6)


class SpacetimeDiagrams(unittest.TestCase):
    def setUp(self):
        self.views = {k: v[0] for k, v in drawn("diagrams")["systems"].items()}

    def lapse_at(self, system, x):
        if system == "uniform":
            return x
        if system in LAPSE:
            return LAPSE[system](x)
        import null_rays as nr
        import numpy as np
        return float(nr._e12_star(np.array([x * x]), 0)[0])

    @unittest.skipUnless(HAS_NUMPY, "numpy and scipy are not installed")
    def test_every_cone_opens_at_the_speed_of_light_of_its_place(self):
        for system, view in self.views.items():
            x0, x1, y0, y1 = view["box"]
            self.assertTrue(view["cones"], system)
            for cone in view["cones"]:
                x = x0 + cone["at"][0] * (x1 - x0)
                for edge in (cone["a"], cone["b"]):
                    across, up = edge[0] * (x1 - x0), edge[1] * (y1 - y0)
                    self.assertGreater(up, 0, system)
                    self.assertAlmostEqual(abs(across) / up, self.lapse_at(system, x), delta=2e-3, msg=system)

    def test_the_singular_edges_are_where_the_lapse_vanishes(self):
        marked = {k: sorted(e for m in v["markers"] if m["kind"] == "singular" for e in m["edges"])
                  for k, v in self.views.items()}
        self.assertEqual(marked, {"static": [], "uniform": [], "february": ["left"], "march": ["left"]})
        self.assertEqual(self.views["february"]["box"][0], 1)      # r = m
        self.assertEqual(self.views["march"]["box"][0], 0.5)       # r = m/2
        self.assertEqual(self.views["uniform"]["box"][0], 0)       # z = 0, a horizon of flat spacetime

    def test_the_free_chart_is_drawn_through_the_declared_star(self):
        view = self.views["static"]
        self.assertTrue(view["input"].startswith("Any $N$"))
        shells = [m for m in view["markers"] if m["kind"] == "shell"]
        x0, x1 = view["box"][:2]
        self.assertEqual(sorted(round(x0 + m["lines"][0][0][0] * (x1 - x0), 3) for m in shells), [-1.0, 1.0])

    def test_the_moment_of_the_embedding_diagram_is_marked_on_the_body_of_march_alone(self):
        for system, view in self.views.items():
            self.assertEqual(bool(view.get("slices")), system == "march", system)


class ConformalDiagrams(unittest.TestCase):
    def setUp(self):
        self.views = {v["id"]: v for v in drawn("conformal")["views"]}

    def classes(self, view):
        return [(layer["kind"], layer["class"]) for layer in self.views[view]["layers"]]

    def test_each_body_is_a_full_diamond_whose_left_edges_are_singular_null_lines(self):
        for view in ("february", "march"):
            zigs = [layer["points"] for layer in self.views[view]["layers"] if layer["kind"] == "zig"]
            self.assertEqual(len(zigs), 2, view)
            for (xa, ta), (xb, tb) in zigs:
                self.assertAlmostEqual(abs(tb - ta), abs(xb - xa), places=3)      # null
                for X, T in ((xa, ta), (xb, tb)):
                    self.assertAlmostEqual(abs(T), math.pi + X, places=3)         # the left edges
            self.assertEqual(self.classes(view).count(("line", "scri")), 2)

    def test_the_star_fills_minkowskis_diamond_with_no_singularity(self):
        self.assertNotIn("zig", {kind for kind, _ in self.classes("static")})
        self.assertEqual(self.classes("static").count(("line", "scri")), 4)
        self.assertIn(("fill", "star"), self.classes("static"))
        self.assertEqual(self.classes("static").count(("line", "surface")), 2)

    def test_the_uniform_field_is_rindlers_wedge_with_a_horizon(self):
        self.assertNotIn("zig", {kind for kind, _ in self.classes("uniform")})
        self.assertEqual(self.classes("uniform").count(("line", "horizon")), 4)
        cover = next(layer for layer in self.views["uniform"]["layers"] if layer["class"] == "cover")
        corners = sorted(tuple(round(c, 3) for c in at) for at in cover["points"])
        half = round(math.pi / 2, 3)
        self.assertEqual(corners, sorted([(0.0, 0.0), (half, -half), (round(math.pi, 3), 0.0), (half, half)]))


class EmbeddingDiagram(unittest.TestCase):
    def setUp(self):
        self.pieces = {p["id"]: p for p in drawn("embedding")["views"][0]["surfaces"][0]["pieces"]}

    def test_the_equator_outside_the_body_is_a_flat_plane(self):
        sheet = self.pieces["sheet"]
        self.assertEqual((sheet["metric"], sheet["system"]), ("einstein_1912_static", "march"))
        self.assertEqual(sheet["points"][0][0], 0.5)
        for r, rho, z in sheet["points"]:
            self.assertAlmostEqual(rho, r, places=9)
            self.assertAlmostEqual(z, 0.0, places=9)

    def test_flamms_paraboloid_of_the_same_mass_rises_over_it_from_its_throat(self):
        flamm = self.pieces["flamm"]
        self.assertTrue(flamm["reference"])
        self.assertEqual((flamm["metric"], flamm["system"]), ("schwarzschild", "spherical"))
        self.assertEqual(flamm["points"][0][:2], [2.0, 2.0])
        for r, rho, z in flamm["points"]:
            self.assertAlmostEqual(rho, r, places=6)
            self.assertAlmostEqual(z, 2 * math.sqrt(2 * (r - 2)), places=5)


if __name__ == "__main__":
    unittest.main()
