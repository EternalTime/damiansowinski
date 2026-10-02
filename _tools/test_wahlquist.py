"""Wahlquist's rotating perfect fluid: python3 -m unittest discover -s _tools

What its texts state, held on the published files: in every chart the published Einstein tensor
is a perfect fluid's at rest in the chart, with rho + 3p the stated constant; Mars's chart at
beta = 0, mu_0 = 0 and a_1 = 0 is the published metric of Kerr; the body the diagrams draw has
its axis at the root of the published h_2, the gamma that keeps that axis regular, its surface
on the equator where Whittaker's sphere has its own, and is longer through its poles than across
its equator; and the embedded equatorial plane holds the ring and the surface where the caption
puts them. The tests of the published mathematics need sympy and mpmath and are skipped where
sympy is absent; the tests of the drawings read the files alone.
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

# The body of the diagrams, in units of r_0.
K, B = "3/10", "4/5"
ETA_0, GAMMA, XI_S, X_S = 1.02460115214, 1.02881004732, 2.80977761761, 1.00271218011


def chart(system, metric_id="wahlquist"):
    metric = json.loads((DATA / "metrics" / f"{metric_id}.json").read_text(encoding="utf-8"))
    return next(c for c in metric["coordinates"] if c["id"] == system)


@unittest.skipUnless(HAS_SYMPY, "sympy is not installed")
class Published(unittest.TestCase):
    def setUp(self):
        import sympy
        import verify_metrics
        self.sp, self.vm = sympy, verify_metrics

    def read(self, system, metric_id="wahlquist"):
        """The reader of a published chart, with the functions the checker holds held, and a
        function from a field to the matrix of its published components, with c = 1 and every
        held function written out as the file defines it."""
        entry = chart(system, metric_id)
        held = self.vm.HELD.get((metric_id, system), ())
        reader = self.vm.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]], (), held=held)
        names = entry["coords"]

        def matrix(field, variant=None):
            block = entry[field] if variant is None else entry[field]["variants"][variant]["nonzero"]
            M = self.sp.zeros(4, 4)
            for e in block:
                i, j = (names.index(k) for k in e["indices"])
                M[i, j] = reader(e["value"]).subs(reader.c, 1).subs(reader.held).doit()
            return M
        return reader, [reader.symbol[n] for n in names], matrix

    def number(self, expression, at):
        return complex(self.sp.sympify(expression).subs(at).evalf(40))

    def fluid(self, system, at, constant):
        """G^mu_nu = (rho + p) u^mu u_nu + p delta^mu_nu with u along the time, and rho + 3p."""
        reader, x, matrix = self.read(system)
        at = {**{reader.parameters[k]: self.sp.Rational(v) for k, v in at["parameters"].items()},
              **{s: self.sp.Rational(v) for s, v in zip(x, at["point"])}}
        G, g = matrix("einstein_tensor", "ul"), matrix("metric_components")
        value = [[self.number(G[i, j], at) for j in range(4)] for i in range(4)]
        metric = [[self.number(g[i, j], at) for j in range(4)] for i in range(4)]
        pressure = value[2][2]
        density = 3 * pressure - sum(value[i][i] for i in range(4))
        for i in range(4):
            for j in range(4):
                expected = pressure * (i == j) - ((density + pressure) * metric[0][j] / metric[0][0] if i == 0 else 0)
                self.assertLess(abs(value[i][j] - expected), 1e-12, (system, i, j))
        self.assertLess(abs(density + 3 * pressure - complex(self.sp.Rational(constant))), 1e-12)
        return density.real, pressure.real

    def test_wahlquists_chart_is_a_perfect_fluid_with_rho_plus_3p_constant(self):
        # 2 k^2/(b^2 r_0^2) at k = 3/10, b = 4/5 and r_0 = 3/2.
        density, pressure = self.fluid(
            "wahlquist", {"parameters": {"r_0": "3/2", "k": K, "b": B, "eta_0": "41/40", "gamma": "103/100"},
                          "point": ("0", "7/5", "2/5", "0")}, "1/8")
        self.assertGreater(pressure, 0)
        self.assertGreater(density, 0)

    def test_marss_charts_are_a_perfect_fluid_with_rho_plus_3p_twice_mu_0(self):
        constants = {"Q_0": "1", "nu_0": "9/10", "mu_0": "1/7", "a_1": "1/5", "a_2": "-3/10", "beta": "2/5"}
        for system in ("mars", "mars_ingoing"):
            self.fluid(system, {"parameters": constants, "point": ("0", "6/5", "7/10", "0")}, "2/7")

    def test_whittakers_sphere_is_the_fluid_at_rest(self):
        density, pressure = self.fluid("whittaker", {"parameters": {"R_0": "1", "b": B},
                                                     "point": ("0", "1/2", "11/10", "0")}, "25/8")
        reader, x, matrix = self.read("whittaker")
        g = matrix("metric_components")
        self.assertEqual([(i, j) for i in range(4) for j in range(4) if i != j and g[i, j] != 0], [])
        self.assertGreater(pressure, 0)

    def test_the_pressure_vanishes_where_minus_g_tt_is_one_over_b_squared(self):
        reader, x, matrix = self.read("whittaker")
        at = {reader.parameters["R_0"]: 1, reader.parameters["b"]: self.sp.Rational(B),
              x[1]: self.sp.Float(str(X_S), 30), x[2]: self.sp.Rational(1, 2)}
        self.assertLess(abs(self.number(matrix("einstein_tensor", "ul")[2, 2], at)), 1e-9)
        self.assertLess(abs(self.number(matrix("metric_components")[0, 0], at) + 25 / 16), 1e-9)

    def test_marss_chart_without_the_fluid_is_kerr(self):
        """At beta -> 0, mu_0 = 0, a_1 = 0, nu_0 = 1, Q_0 = a^2 and a_2 = -2M, along y = r,
        z = a cos(theta), tau = t - a phi and sigma = -phi/a, Mars's (10)."""
        sp = self.sp
        reader, x, matrix = self.read("mars")
        kerr, kx, kmatrix = self.read("boyer_lindquist", "kerr")
        t, r, theta, phi = kx
        M, a = sp.Rational(1), sp.Rational(3, 5)
        p = reader.parameters
        constants = {p["Q_0"]: a ** 2, p["nu_0"]: 1, p["mu_0"]: 0, p["a_1"]: 0, p["a_2"]: -2 * M,
                     p["beta"]: sp.Rational(1, 10 ** 15)}
        image = dict(zip(x, [t - a * phi, r, a * sp.cos(theta), -phi / a]))
        jacobian = sp.Matrix(4, 4, lambda i, j: sp.diff(image[x[i]], kx[j]))
        pulled = jacobian.T * matrix("metric_components").subs(constants).subs(image, simultaneous=True) * jacobian
        theirs = kmatrix("metric_components").subs({kerr.parameters["M"]: M, kerr.parameters["a"]: a,
                                                     kerr.parameters["G"]: 1})
        at = {t: 0, r: sp.Rational(5, 2), theta: sp.Rational(7, 10), phi: 0}
        for i in range(4):
            for j in range(4):
                self.assertLess(abs(self.number(pulled[i, j] - theirs[i, j], at)), 1e-12, (i, j))

    def test_the_body_of_the_diagrams(self):
        """eta_0 is the root of the published h_2, gamma keeps the axis regular, the surface
        crosses the equator where k xi = sin(X_s), and the pole is farther from the centre than
        the equator is."""
        import mpmath
        sp = self.sp
        reader, x, matrix = self.read("wahlquist")
        xi, eta = x[1], x[2]
        p = reader.parameters
        constants = {p["k"]: sp.Rational(K), p["b"]: sp.Rational(B)}
        by_name = {str(f.func): definition.subs(constants) for f, definition in reader.held.items()}
        h1 = sp.lambdify(xi, by_name["h_1"], "mpmath")
        h2 = sp.lambdify(eta, by_name["h_2"], "mpmath")
        with mpmath.workdps(30):
            k, b = mpmath.mpf(3) / 10, mpmath.mpf(4) / 5
            eta0 = mpmath.findroot(h2, 1.02)
            self.assertLess(abs(eta0 - ETA_0), 1e-10)
            gamma = 2 / (mpmath.sqrt(1 + k ** 2 * eta0 ** 2) * abs(mpmath.diff(h2, eta0)))
            self.assertLess(abs(gamma - GAMMA), 1e-10)
            equator = mpmath.findroot(lambda s: (h1(s) - 1) / s ** 2 - 1 / b ** 2, 2.8)
            self.assertLess(abs(equator - XI_S), 1e-10)
            self.assertLess(abs(k * equator - mpmath.sin(mpmath.findroot(lambda s: s * mpmath.cot(s) - b ** 2, 1.0))), 1e-25)
            pole = mpmath.findroot(lambda s: h1(s) / (s ** 2 + eta0 ** 2) - 1 / b ** 2, 2.9)
            # The proper distances from the centre, in units of r_0: along the axis to the pole, and
            # across the disc and the plane to the equator.
            polar = mpmath.quad(lambda s: mpmath.sqrt((s ** 2 + eta0 ** 2) / ((1 - k ** 2 * s ** 2) * h1(s))), [0, pole])
            across = (mpmath.quad(lambda s: s / mpmath.sqrt((1 + k ** 2 * s ** 2) * h2(s)), [0, eta0]).real
                      + mpmath.quad(lambda s: s / mpmath.sqrt((1 - k ** 2 * s ** 2) * h1(s)), [0, equator]))
            self.assertAlmostEqual(float(polar), 3.332, places=3)
            self.assertAlmostEqual(float(across), 3.310, places=3)
            self.assertGreater(polar, across)


class Drawn(unittest.TestCase):
    def test_the_embedded_plane_of_the_rotating_body(self):
        views = {v["id"]: v for v in json.loads((DATA / "embedding" / "wahlquist.json").read_text(encoding="utf-8"))["views"]}
        pieces = {p["id"]: p for p in views["rotating"]["surfaces"][0]["pieces"]}
        disc, plane = pieces["disc"]["points"], pieces["plane"]["points"]
        # The disc and the plane meet on the ring, whose circle has the radius 1.028 r_0.
        self.assertEqual(disc[0][1:], plane[0][1:])
        self.assertAlmostEqual(plane[0][1], 1.0275, places=4)
        # The disc ends on the axis, at eta_0, and the plane on the surface, whose circle is 2.959 r_0.
        self.assertAlmostEqual(disc[-1][0], ETA_0, places=6)
        self.assertLess(disc[-1][1], 1e-4)
        self.assertAlmostEqual(plane[-1][0], XI_S, places=9)
        self.assertAlmostEqual(plane[-1][1], 2.959, places=3)
        # A bowl: the centre lowest, the surface highest, and the circles growing all the way out.
        self.assertLess(disc[-1][2], 0)
        self.assertGreater(plane[-1][2], 0)
        radii = [p[1] for p in reversed(disc)] + [p[1] for p in plane[1:]]
        self.assertEqual(radii, sorted(radii))

    def test_whittakers_sphere_ends_where_its_pressure_vanishes(self):
        views = {v["id"]: v for v in json.loads((DATA / "embedding" / "wahlquist.json").read_text(encoding="utf-8"))["views"]}
        (ball,) = views["static"]["surfaces"][0]["pieces"]
        self.assertAlmostEqual(ball["points"][-1][0], X_S, places=9)
        self.assertAlmostEqual(ball["points"][-1][1], math.sin(X_S), places=6)
        self.assertLess(abs(X_S / math.tan(X_S) - 0.64), 1e-9)

    def test_the_surface_is_marked_on_each_plane_that_reaches_it(self):
        systems = json.loads((DATA / "diagrams" / "wahlquist.json").read_text(encoding="utf-8"))["systems"]
        marked = {(system, v["id"]): [m["kind"] for m in v["markers"]] for system, views in systems.items() for v in views}
        for where in (("wahlquist", "equator"), ("mars", "equator"), ("mars_ingoing", "equator"),
                      ("whittaker", "radial"), ("whittaker", "through")):
            self.assertIn("surface", marked[where], where)
        self.assertNotIn("surface", marked[("wahlquist", "disc")])


if __name__ == "__main__":
    unittest.main()
