"""Kerr's black hole in Melvin's magnetic universe, Ernst and Wild's solution:
python3 -m unittest discover -s _tools

What its texts and drawings state, held on the published files. The first class needs nothing
but the files. The second needs numpy, scipy and sympy and is skipped where they are absent: it
writes the published definitions out and holds them to Gibbons, Mujtaba and Pope's polynomials
and to Bicak and Hejda's complex factor, to the Einstein-Maxwell equations at a point no chart
was printed at, to Wald's charge 2amB by the flux of the dual field through a sphere, to the
regular axis, to the rate of the horizon, to the curvature at r = 0 on the axis and on the
ring, and to the ergoregion that runs to infinity beside the axis for every turning azimuth.
"""
import importlib.util
import json
import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "MFS" / "assets" / "data"
HAS_NUMBERS = all(importlib.util.find_spec(name) is not None for name in ("sympy", "numpy", "scipy", "mpmath"))
sys.path.insert(0, str(Path(__file__).resolve().parent / "derivations"))


def load(folder):
    return json.loads((DATA / folder / "kerr_melvin.json").read_text(encoding="utf-8"))


class Drawings(unittest.TestCase):
    def test_each_chart_draws_the_axis_and_the_equator(self):
        systems = load("diagrams")["systems"]
        for system in ("boyer_lindquist", "rotating"):
            self.assertEqual({v["id"] for v in systems[system]}, {"radial", "equator"})
        self.assertEqual({v["id"] for v in load("conformal")["views"]}, {"axis", "equator"})
        self.assertEqual([v["id"] for v in load("embedding")["views"]], ["equator"])
        figures = load("diagrams")["projections"]["boyer_lindquist"]
        self.assertEqual([f["id"] for f in figures], ["dragging", "tube"])

    def test_the_spacetime_diagrams_are_square(self):
        for views in load("diagrams")["systems"].values():
            for view in views:
                X0, X1, Y0, Y1 = view["box"]
                self.assertAlmostEqual((X1 - X0) / (Y1 - Y0), 1.0)

    def test_only_the_chart_at_rest_marks_an_ergosurface_on_the_equator(self):
        """In the chart that turns with the horizon g_tt is negative from r_+ out on the equator,
        so no line of g_tt = 0 is drawn there, and in the chart at rest it is drawn at 1.942 m."""
        systems = load("diagrams")["systems"]
        marked = {system: "the ergosurface" in json.dumps({k: v for k, v in view.items() if k != "caption"})
                  for system, views in systems.items() for view in views if view["id"] == "equator"}
        self.assertEqual(marked, {"boyer_lindquist": True, "rotating": False})

    def test_the_equator_widens_to_k_over_B_and_narrows(self):
        """The circle of r on the equator has radius k sqrt(P): 2mk/(1 + B^2 m^2) = 1.8871 m at the
        throat r_+ = 8m/5, greatest, k/B = 4.01 m, at r = 7.950 m, and smaller beyond."""
        view, = load("embedding")["views"]
        points = view["surfaces"][0]["pieces"][0]["points"]
        self.assertAlmostEqual(points[0][0], 1.6, places=6)
        self.assertAlmostEqual(points[0][1], 2 * 1.0025 / (1 + 1 / 16), places=5)
        radius, at = max((p[1], p[0]) for p in points)
        self.assertAlmostEqual(radius, 4.01, places=4)
        self.assertAlmostEqual(at, 7.95, delta=0.1)
        self.assertLess(points[-1][1], radius)
        self.assertAlmostEqual(points[-1][0], 10.0, places=6)

    def test_the_tube_is_drawn_apart_from_the_region_round_the_hole(self):
        """The figure of the ergoregion fills the region round the hole and one tube each way
        along the axis, on both sides of it, and covers the hole itself."""
        figure = next(f for f in load("diagrams")["projections"]["boyer_lindquist"] if f["id"] == "tube")
        fills = [layer for layer in figure["layers"] if layer["kind"] == "fill"]
        self.assertEqual([layer["class"] for layer in fills].count("double"), 6)
        self.assertEqual([layer["class"] for layer in fills].count("wedge"), 1)
        self.assertNotIn("turn", figure)
        X0, X1, Y0, Y1 = figure["box"]
        self.assertAlmostEqual((X1 - X0) / (Y1 - Y0), 1.0, places=2)
        tubes = [layer for layer in fills if layer["class"] == "double"
                 and max(abs(p[1]) for p in layer["points"]) > 11.9]
        self.assertEqual(len(tubes), 4)
        for layer in tubes:
            # A tube starts about four m from the hole and never touches the axis.
            self.assertGreater(min(abs(p[1]) for p in layer["points"]), 3.5)
            self.assertGreater(min(abs(p[0]) for p in layer["points"]), 0.1)


@unittest.skipUnless(HAS_NUMBERS, "needs numpy, scipy and sympy")
class Physics(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import sympy as sp
        import verify_metrics as vm
        cls.sp, cls.vm = sp, vm
        metric = load("metrics")
        cls.metric = metric
        cls.entry = {c["id"]: c for c in metric["coordinates"]}
        entry = cls.entry["boyer_lindquist"]
        cls.reader = R = vm.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]], (),
                                   held=vm.HELD[("kerr_melvin", "boyer_lindquist")])
        cls.x = [R.symbol[name] for name in entry["coords"]]
        cls.r, cls.th = cls.x[1], cls.x[2]
        cls.m, cls.a, cls.B, cls.k = (R.parameters[name] for name in ("m", "a", "B", "k"))
        g = sp.zeros(4, 4)
        for component in entry["metric_components"]:
            i, j = (entry["coords"].index(name) for name in component["indices"])
            g[i, j] = R(component["value"])
        cls.g = cls.written(g).subs(R.c, 1)
        cls.H, cls.omega, cls.A = (cls.written(sp.sympify(R.parameters[name])) for name in ("H", "omega", "A"))

    @classmethod
    def written(cls, value):
        """A value with the held names written out by the published definitions."""
        for _ in cls.reader.held:
            value = value.subs(cls.reader.held).doit()
        return value

    def point(self, r, sine, cosine, a, m, B):
        sp = self.sp
        return ({self.r: sp.Rational(r), self.a: sp.Rational(a), self.m: sp.Rational(m), self.B: sp.Rational(B)},
                {sp.sin(self.th): sp.Rational(sine), sp.cos(self.th): sp.Rational(cosine)})

    def geometry(self, at, trig, potential=None):
        """The Riemann tensor, the Ricci tensor and Maxwell's stress at a point of rational
        coordinates and parameters, exactly, from the written out metric and its derivatives."""
        sp, x, n = self.sp, self.x, 4
        number = lambda e: sp.cancel(e.subs(at).subs(trig))
        G = sp.Matrix(n, n, lambda i, j: number(self.g[i, j]))
        d1 = [[[number(sp.diff(self.g[i, j], x[k])) for k in range(n)] for j in range(n)] for i in range(n)]
        d2 = [[[[number(sp.diff(self.g[i, j], x[k], x[l])) if {k, l} <= {1, 2} else sp.Integer(0)
                 for l in range(n)] for k in range(n)] for j in range(n)] for i in range(n)]
        Gi = G.inv()
        first = [[[(d1[c][a][b] + d1[c][b][a] - d1[a][b][c]) / 2 for b in range(n)] for a in range(n)] for c in range(n)]
        gamma = [[[sum(Gi[i, c] * first[c][a][b] for c in range(n)) for b in range(n)] for a in range(n)] for i in range(n)]
        riemann = [[[[(d2[a][d][b][c] + d2[b][c][a][d] - d2[a][c][b][d] - d2[b][d][a][c]) / 2
                      + sum(gamma[e][b][c] * first[e][a][d] - gamma[e][b][d] * first[e][a][c] for e in range(n))
                      for d in range(n)] for c in range(n)] for b in range(n)] for a in range(n)]
        ricci = sp.Matrix(n, n, lambda b, d: sum(Gi[a, c] * riemann[a][b][c][d] for a in range(n) for c in range(n)))
        out = {"g": G, "inverse": Gi, "riemann": riemann, "ricci": ricci,
               "scalar": sp.cancel(sum(Gi[i, j] * ricci[i, j] for i in range(n) for j in range(n)))}
        if potential is not None:
            F = sp.Matrix(n, n, lambda i, j: number(sp.diff(potential[j], x[i]) - sp.diff(potential[i], x[j])))
            mixed = F * Gi
            square = sum((Gi * F * Gi)[i, j] * F[i, j] for i in range(n) for j in range(n))
            out["stress"] = sp.Matrix(n, n, lambda i, j: sum(F[i, k] * mixed[j, k] for k in range(n)) - G[i, j] * square / 4)
        return out

    def potential(self):
        """Phi_0 dx^0 + Phi_3 (k dphi - omega dx^0), Phi_3 = (dH/dB)/(2H), Phi_0 = (2/B)(omega - 2mra/A)."""
        sp = self.sp
        azimuthal = sp.diff(self.H, self.B) / (2 * self.H)
        timelike = 2 * (self.omega - 2 * self.m * self.r * self.a / self.A) / self.B
        return [timelike - azimuthal * self.omega, 0, 0, self.k * azimuthal]

    def test_the_line_elements_are_the_printer_s(self):
        import print_charts as pc
        for system, entry in self.entry.items():
            self.assertEqual(entry["line_element"], pc.kerr_melvin(system)["system"]["line_element"])
            self.assertEqual([p["symbol"] for p in entry["parameters"]], pc.kerr_melvin_parameters(system))

    def test_H_and_omega_are_gibbons_mujtaba_and_pope_s(self):
        """Their appendix B at q = p = 0: H = 1 + (H_2 B^2 + H_4 B^4)/R^2 and
        omega = (2mra + omega_4 B^4)/Sigma, with their R^2 the Sigma here and their Sigma the A."""
        sp, r, a, m, B = self.sp, self.r, self.a, self.m, self.B
        s, c = sp.sin(self.th), sp.cos(self.th)
        R2 = r ** 2 + a ** 2 * c ** 2
        delta = r ** 2 + a ** 2 - 2 * m * r
        big = (r ** 2 + a ** 2) ** 2 - a ** 2 * delta * s ** 2
        H2 = big * s ** 2 / 2
        H4 = ((r ** 2 + a ** 2) ** 2 * R2 * s ** 4 / 16 + m * a ** 2 * r * (r ** 2 + a ** 2) * s ** 6 / 4
              + m ** 2 * a ** 2 * (r ** 2 * (c ** 2 - 3) ** 2 * c ** 2 + a ** 2 * (1 + c ** 2) ** 2) / 4)
        w4 = (a ** 3 * m ** 3 * r * (3 + c ** 4) / 2
              + a * m ** 2 * (r ** 4 * (3 - 6 * c ** 2 + c ** 4) + 2 * a ** 2 * r ** 2 * (3 * s ** 2 - 2 * c ** 4)
                              - a ** 4 * (1 + c ** 4)) / 4
              + a * m * r * (r ** 2 + a ** 2) * (r ** 2 * (3 + 6 * c ** 2 - c ** 4)
                                                 - a ** 2 * (1 - 6 * c ** 2 - 3 * c ** 4)) / 8)
        theirs = {"H": 1 + (H2 * B ** 2 + H4 * B ** 4) / R2, "omega": (2 * m * r * a + w4 * B ** 4) / big}
        # Bicak and Hejda's complex factor for the uncharged seed, whose squared modulus H is.
        real = 1 + B ** 2 * big * s ** 2 / (4 * R2)
        imaginary = -B ** 2 * m * a * c * (3 - c ** 2 + a ** 2 * s ** 4 / R2) / 2
        for at, trig in (self.point("23/10", "3/5", "4/5", "3/5", 1, "2/7"),
                         self.point("7/2", "12/13", "-5/13", "9/10", "6/5", "1/2"),
                         self.point("1/3", "8/17", "15/17", "1/2", "3/4", "5/3")):
            number = lambda e: sp.cancel(e.subs(at).subs(trig))
            self.assertEqual(number(self.H - theirs["H"]), 0)
            self.assertEqual(number(self.omega - theirs["omega"]), 0)
            self.assertEqual(number(self.H - real ** 2 - imaginary ** 2), 0)

    def test_the_einstein_maxwell_equations_hold(self):
        """G = 2 (F F - g F^2/4) and R = 0, exactly, at a point no chart was printed at."""
        sp = self.sp
        at, trig = self.point("17/5", "5/13", "12/13", "7/10", "11/10", "3/8")
        found = self.geometry(at, trig, self.potential())
        einstein = found["ricci"] - found["g"] * found["scalar"] / 2
        self.assertEqual(found["scalar"], 0)
        self.assertEqual((einstein - 2 * found["stress"]).applyfunc(sp.cancel), sp.zeros(4, 4))
        self.assertNotEqual(einstein, sp.zeros(4, 4))

    def test_the_published_curvature_is_the_metric_s(self):
        """The Kretschmann scalar and the Ricci scalar as published, in the held names and their
        derivatives, against the curvature taken from the written out metric at one point."""
        sp, R = self.sp, self.reader
        entry = self.entry["boyer_lindquist"]
        at, trig = self.point("17/5", "5/13", "12/13", "7/10", "11/10", "3/8")
        found = self.geometry(at, trig)
        n = 4
        lift = found["inverse"]
        riemann = found["riemann"]
        raised = [[[[sum(lift[a, e] * lift[b, f] * riemann[e][f][c][d] for e in range(n) for f in range(n))
                     for d in range(n)] for c in range(n)] for b in range(n)] for a in range(n)]
        up = [[[[sum(lift[c, e] * lift[d, f] * raised[a][b][e][f] for e in range(n) for f in range(n))
                 for d in range(n)] for c in range(n)] for b in range(n)] for a in range(n)]
        kretschmann = sp.cancel(sum(riemann[a][b][c][d] * up[a][b][c][d] for a in range(n) for b in range(n)
                                    for c in range(n) for d in range(n)))
        published = R(entry["kretschmann"].split("=", 1)[1])
        values = {}
        for atom in published.atoms(sp.Derivative) | (published.atoms(sp.core.function.AppliedUndef) & set(R.held)):
            values[atom] = sp.cancel(self.written(atom).subs(at).subs(trig))
        self.assertEqual(sp.cancel(published.xreplace(values).subs(at).subs(R.c, 1) - kretschmann), 0)
        # The Ricci scalar is published in the held names, a sum of their second derivatives, and
        # for the published definitions it is zero.
        scalar = R(entry["ricci_scalar"].split("=", 1)[1])
        for atom in scalar.atoms(sp.Derivative) | (scalar.atoms(sp.core.function.AppliedUndef) & set(R.held)):
            values[atom] = sp.cancel(self.written(atom).subs(at).subs(trig))
        self.assertNotEqual(scalar, 0)
        self.assertEqual(sp.cancel(scalar.xreplace(values).subs(at).subs(R.c, 1)), 0)
        self.assertEqual(found["scalar"], 0)

    def test_the_hole_carries_wald_s_charge(self):
        """Q = (1/4 pi) times the flux of the dual field through a sphere of constant r is 2amB
        in size, the same through every sphere outside the horizon."""
        import mpmath
        sp, x = self.sp, self.x
        values = {self.a: sp.Rational(4, 5), self.m: 1, self.B: sp.Rational(1, 4)}
        g = self.g.subs(values)
        potential = [sp.sympify(p).subs(values) for p in self.potential()]
        F = sp.Matrix(4, 4, lambda i, j: sp.diff(potential[j], x[i]) - sp.diff(potential[i], x[j]))
        for radius in (sp.Rational(5, 2), 6):
            at = {self.r: radius}
            G = g.subs(at)
            Gi_tt, Gi_tp, Gi_rr = (sp.lambdify(self.th, e, "mpmath") for e in (
                (G[3, 3] / (G[0, 0] * G[3, 3] - G[0, 3] ** 2)), (-G[0, 3] / (G[0, 0] * G[3, 3] - G[0, 3] ** 2)), 1 / G[1, 1]))
            F_tr, F_pr = (sp.lambdify(self.th, F[i, 1].subs(at), "mpmath") for i in (0, 3))
            root = sp.lambdify(self.th, -G[1, 1] * G[2, 2] * (G[0, 0] * G[3, 3] - G[0, 3] ** 2), "mpmath")

            def density(theta):
                # sqrt(-g) F^{tr}, the theta phi component of the dual two-form.
                return mpmath.sqrt(root(theta)) * Gi_rr(theta) * (Gi_tt(theta) * F_tr(theta) + Gi_tp(theta) * F_pr(theta))
            charge = mpmath.quad(density, [0, mpmath.pi / 2, mpmath.pi]) / 2
            self.assertAlmostEqual(abs(float(charge)), 2 * 0.8 * 1 * 0.25, places=9)

    def test_the_axis_is_regular_and_the_horizon_turns_rigidly(self):
        sp, r, th, a, m, B = self.sp, self.r, self.th, self.a, self.m, self.B
        for pole in (0, sp.pi):
            self.assertEqual(sp.simplify(self.H.subs(th, pole) - self.k), 0)
        # g_phiphi / (theta^2 g_thetatheta) -> 1 on the axis, a circle's radius over its distance out.
        at = {r: sp.Rational(3), a: sp.Rational(4, 5), m: 1, B: sp.Rational(1, 4)}
        ratio = sp.limit((self.g[3, 3] / self.g[2, 2]).subs(at) / th ** 2, th, 0)
        self.assertEqual(sp.nsimplify(ratio), 1)
        rate = sp.Rational(4, 5) / (2 * sp.Rational(8, 5)) + sp.Rational(4, 5) * sp.Rational(1, 256) * sp.Rational(13, 5) / 2
        self.assertEqual(rate, sp.Rational(813, 3200))
        for sine, cosine in ((0, 1), ("3/5", "4/5"), (1, 0), ("5/13", "-12/13")):
            value = self.omega.subs({r: sp.Rational(8, 5), a: sp.Rational(4, 5), m: 1, B: sp.Rational(1, 4)})
            value = value.subs({sp.sin(th): sp.Rational(sine), sp.cos(th): sp.Rational(cosine)})
            self.assertEqual(sp.cancel(value), rate)
        # The second chart's diagrams turn at Omega = omega_+/k.
        self.assertEqual(rate / sp.Rational(401, 400), sp.Rational(813, 3208))

    def test_the_ring_is_singular_and_the_centre_of_its_disc_is_not(self):
        """Toward r = 0 on the equator the Kretschmann scalar grows without bound, and toward
        r = 0 on the axis it settles."""
        ring = []
        for r, t in (("1/10", "1/20"), ("1/100", "1/200"), ("1/1000", "1/2000")):
            q = self.sp.Rational(t)
            at, trig = self.point(r, (1 - q ** 2) / (1 + q ** 2), 2 * q / (1 + q ** 2), "4/5", 1, "1/4")
            ring.append(float(self.kretschmann(at, trig)))
        self.assertGreater(abs(ring[1]), 10 * abs(ring[0]))
        self.assertGreater(abs(ring[2]), 10 * abs(ring[1]))
        axis = []
        for r, t in (("1/10", "1/1000"), ("1/100", "1/1000"), ("1/1000", "1/1000")):
            q = self.sp.Rational(t)
            at, trig = self.point(r, 2 * q / (1 + q ** 2), (1 - q ** 2) / (1 + q ** 2), "4/5", 1, "1/4")
            axis.append(float(self.kretschmann(at, trig)))
        self.assertLess(abs(axis[2] - axis[1]), 0.05 * abs(axis[1]))
        self.assertLess(abs(axis[1]), 1e4)

    def kretschmann(self, at, trig):
        sp, n = self.sp, 4
        found = self.geometry(at, trig)
        lift, riemann = found["inverse"], found["riemann"]
        pairs = [(a, b) for a in range(n) for b in range(a + 1, n)]
        up = {(p, q): lift[p[0], q[0]] * lift[p[1], q[1]] - lift[p[0], q[1]] * lift[p[1], q[0]] for p in pairs for q in pairs}
        half = {(p, q): sum(up[p, s] * riemann[s[0]][s[1]][q[0]][q[1]] for s in pairs) for p in pairs for q in pairs}
        raised = {(p, q): sum(half[p, s] * up[s, q] for s in pairs) for p in pairs for q in pairs}
        return sp.cancel(4 * sum(riemann[p[0]][p[1]][q[0]][q[1]] * raised[p, q] for p in pairs for q in pairs))

    def test_the_ergoregion_reaches_infinity_beside_the_axis(self):
        """At a fixed distance from the axis g_tt grows as the square of the height, and for an
        azimuth turning at any one constant rate Omega the time translation K = d_t + Omega d_phi
        is spacelike far enough up, since the rate that holds it off, omega/k, grows there. On the
        equator g_tt is negative from the ergosurface out, and in the chart turning with the
        horizon it is negative from the horizon out."""
        import numpy as np
        sp = self.sp
        for field in ("1/4", "3/4"):
            values = {self.a: sp.Rational(4, 5), self.m: 1, self.B: sp.Rational(field)}
            g_tt, g_tp, g_pp = (sp.lambdify((self.r, self.th), self.g[i, j].subs(values), "numpy")
                                for i, j in ((0, 0), (0, 3), (3, 3)))
            rho = 2.0
            at = lambda z: (np.hypot(rho, z), np.arctan2(rho, z))
            far = 1e5
            self.assertGreater(g_tt(*at(far)), 0)
            self.assertAlmostEqual(g_tt(*at(2 * far)) / g_tt(*at(far)), 4.0, places=3)
            for turning in (-1.0, 0.3, 1.0, 10.0):
                up = 1e3 * far * max(1.0, abs(turning))
                self.assertGreater(g_tt(*at(up)) + 2 * turning * g_tp(*at(up)) + turning ** 2 * g_pp(*at(up)), 0)
            self.assertGreater(-g_tp(*at(2 * far)) / g_pp(*at(2 * far)), -1.5 * g_tp(*at(far)) / g_pp(*at(far)))
            r = np.linspace(2.2, 2000, 200001)
            self.assertTrue(np.all(g_tt(r, np.pi / 2) < 0))
            k = float(self.k.subs(values))
            rate = (0.8 / 3.2 + 0.8 * float(sp.Rational(field)) ** 4 * 2.6 / 2) / k
            r = np.linspace(1.6001, 2000, 200001)
            turning = g_tt(r, np.pi / 2) + 2 * rate * g_tp(r, np.pi / 2) + rate ** 2 * g_pp(r, np.pi / 2)
            self.assertTrue(np.all(turning < 0))
        # In the weak field the tube begins beyond 120 m, and in the strong one about 4 m up.
        values = {self.a: sp.Rational(4, 5), self.m: 1}
        for field, lowest in (("1/4", (120, 130)), ("3/4", (3.9, 4.2))):
            g_tt = sp.lambdify((self.r, self.th), self.g[0, 0].subs({**values, self.B: sp.Rational(field)}), "numpy")
            x, z = np.meshgrid(np.linspace(0.05, 12, 400), np.linspace(2.5, 200, 4000))
            positive = g_tt(np.hypot(x, z), np.arctan2(x, z)) > 0
            self.assertTrue(positive.any())
            self.assertTrue(lowest[0] < z[positive].min() < lowest[1], z[positive].min())

    def test_the_numbers_the_texts_state(self):
        """The ergosurface on the equator at 1.942 m, the dragging rate least at 7.950 m, where
        the circle of phi is widest, and growing as 3amB^4 r/8 beyond; Kerr's own light cylinder
        for the horizon's rate at 2.87 m; and the angular momentum am(1 - a^2 m^2 B^4) with
        Gibbons, Mujtaba and Pope's charge for the seed q = -amB half of Wald's in a weak field."""
        import numpy as np
        from scipy.optimize import brentq, minimize_scalar
        sp = self.sp
        values = {self.a: sp.Rational(4, 5), self.m: 1, self.B: sp.Rational(1, 4)}
        g_tt = sp.lambdify(self.r, self.g[0, 0].subs(values).subs(self.th, sp.pi / 2), "numpy")
        self.assertAlmostEqual(brentq(g_tt, 1.7, 2.5), 1.9424, places=4)
        omega = sp.lambdify(self.r, self.omega.subs(values).subs(self.th, sp.pi / 2), "numpy")
        least = minimize_scalar(omega, bounds=(2, 40), method="bounded").x
        self.assertAlmostEqual(least, 7.950, places=2)
        radius = sp.lambdify(self.r, sp.sqrt(self.g[3, 3]).subs(values).subs(self.th, sp.pi / 2), "numpy")
        self.assertAlmostEqual(minimize_scalar(lambda r: -radius(r), bounds=(2, 40), method="bounded").x, least, places=3)
        self.assertAlmostEqual(omega(2000.0) / (3 * 0.8 * 0.25 ** 4 * 2000 / 8), 1.0, places=2)
        kerr = {self.a: sp.Rational(4, 5), self.m: 1, self.B: 0}
        cylinder = sp.lambdify(self.r, (self.g[0, 0] + 2 * sp.Rational(1, 4) * self.g[0, 3]
                                        + sp.Rational(1, 16) * self.g[3, 3]).subs(kerr).subs(self.th, sp.pi / 2), "numpy")
        self.assertAlmostEqual(brentq(cylinder, 2.0, 10.0), 2.866, places=3)
        j, B = 0.8, 1e-3
        self.assertAlmostEqual(j * B * (1 + j ** 2 * B ** 4 / 4) / (2 * j * B), 0.5, places=9)
        self.assertIn("J = am(1 - a^2m^2B^4)", self.metric["history"])
        self.assertIn("Q = jB(1 + \\tfrac{1}{4}j^2B^4)", self.metric["history"])


if __name__ == "__main__":
    unittest.main()
