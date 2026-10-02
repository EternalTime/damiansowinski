"""Hotta and Tanaka's shock wave in de Sitter space, a delta on a curved background:
python3 -m unittest discover -s _tools

The reading of the delta its charts rest on, verify_metrics.on_the_shock, held to the rules of
distributions; what the published charts state, held to Einstein's equations with Lambda = 3/a^2,
to Sfetsos's angle, to Aichelburg and Sexl's logarithm as a goes to infinity and to the jump of a
ray across the shock; and the drawings' data, the ring of particles of the embedding diagram
against Podolsky and Ortaggio's geodesics and its moments against each chart's map into the
hyperboloid. The first two need sympy and the slices numpy; each is skipped where its library is
absent.
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
HAS_NUMPY = importlib.util.find_spec("numpy") is not None
sys.path.insert(0, str(Path(__file__).resolve().parent / "derivations"))

CHARTS = ("conformally_flat", "global", "kruskal", "null_cylindrical", "kundt")


def chart(system):
    metric = json.loads((DATA / "metrics" / "hotta_tanaka.json").read_text(encoding="utf-8"))
    return next(c for c in metric["coordinates"] if c["id"] == system)


@unittest.skipUnless(HAS_SYMPY, "sympy is not installed")
class TheDelta(unittest.TestCase):
    """on_the_shock: a pulse times a function that varies across it, in the limit of a delta."""

    def setUp(self):
        import sympy
        import verify_metrics
        self.sp, self.vm = sympy, verify_metrics
        self.u, self.v = sympy.symbols("u v", real=True)

    def pulse(self, order=0, argument=None):
        return self.vm.Pulse(self.u if argument is None else argument, self.sp.Integer(order))

    def test_a_delta_keeps_the_coefficient_on_the_shock(self):
        sp, u, v = self.sp, self.u, self.v
        self.assertEqual(self.vm.on_the_shock((3 + u * v + u ** 2) * self.pulse()), 3 * sp.DiracDelta(u))

    def test_a_derivative_of_the_delta_follows_leibniz(self):
        sp, u, v = self.sp, self.u, self.v
        # f delta' = f(0) delta' - f'(0) delta, and f delta'' = f(0) delta'' - 2 f'(0) delta' + f''(0) delta.
        self.assertEqual(sp.expand(self.vm.on_the_shock((1 + u * v) * self.pulse(1))),
                         sp.DiracDelta(u, 1) - v * sp.DiracDelta(u))
        self.assertEqual(sp.expand(self.vm.on_the_shock(sp.exp(v * u) * self.pulse(2))),
                         sp.DiracDelta(u, 2) - 2 * v * sp.DiracDelta(u, 1) + v ** 2 * sp.DiracDelta(u))

    def test_a_square_goes_with_a_coefficient_that_vanishes_on_the_shock(self):
        sp, u, v = self.sp, self.u, self.v
        self.assertEqual(self.vm.on_the_shock(u * v * self.pulse() ** 2 + sp.cos(u)), sp.cos(u))
        self.assertEqual(self.vm.on_the_shock(u ** 2 * self.pulse() * self.pulse(1)), 0)

    def test_a_product_no_reading_fixes_is_left_standing(self):
        sp, u = self.sp, self.u
        self.assertEqual(self.vm.on_the_shock(self.pulse() ** 2), sp.DiracDelta(u) ** 2)
        # u delta delta' diverges as the pulse narrows: its coefficient has to vanish twice.
        self.assertNotEqual(self.vm.on_the_shock(u * self.pulse() * self.pulse(1)), 0)

    def test_an_argument_of_two_coordinates_is_read_on_the_first(self):
        sp = self.sp
        eta, chi = sp.symbols("eta chi", real=True)
        found = self.vm.on_the_shock(sp.sin(eta) * self.pulse(1, eta - chi))
        self.assertEqual(sp.expand(found - (sp.sin(chi) * sp.DiracDelta(eta - chi, 1)
                                            - sp.cos(chi) * sp.DiracDelta(eta - chi))), 0)

    def test_the_limit_changes_nothing_already_in_it(self):
        sp, u, v = self.sp, self.u, self.v
        once = self.vm.on_the_shock((1 + u * v) * self.pulse(1) + u * self.pulse())
        self.assertEqual(sp.expand(self.vm.on_the_shock(once) - once), 0)


@unittest.skipUnless(HAS_SYMPY, "sympy is not installed")
class Published(unittest.TestCase):
    def setUp(self):
        import sympy
        import verify_metrics
        self.sp, self.vm = sympy, verify_metrics

    def read(self, system):
        entry = chart(system)
        reader = self.vm.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]], ())
        return entry, reader

    def test_every_chart_is_an_einstein_space_shock_and_all(self):
        """R_mu_nu = (3/a^2) g_mu_nu, the delta of the metric with it, so the shock adds no source
        off the two particles."""
        for system in CHARTS:
            entry, reader = self.read(system)
            a = reader.parameters["a"]
            g = {tuple(e["indices"]): reader(e["value"]) for e in entry["metric_components"]}
            ricci = {tuple(e["indices"]): reader(e["value"]) for e in entry["ricci_tensor"]["variants"]["ll"]["nonzero"]}
            self.assertEqual(set(g), set(ricci), system)
            for index, value in g.items():
                self.assertEqual(self.vm.norm(ricci[index] - 3 * value / a ** 2), 0, f"{system} {index}")
            self.assertEqual(entry["ricci_scalar"], "R = \\dfrac{12}{a^2}")
            self.assertEqual(entry["kretschmann"], "K = \\dfrac{24}{a^4}")

    def test_the_weyl_tensor_lives_on_the_shock(self):
        sp = self.sp
        for system in CHARTS:
            entry, reader = self.read(system)
            for e in entry["weyl_tensor"]["variants"]["llll"]["nonzero"]:
                value = reader(e["value"])
                self.assertTrue(value.has(sp.DiracDelta), f"{system} {e['indices']}")
                self.assertEqual(sp.simplify(value.replace(lambda x: isinstance(x, sp.DiracDelta), lambda x: 0)), 0)

    def test_the_shift_changes_sign_at_sfetsos_angle(self):
        """cos(theta) ln((1 + cos theta)/(1 - cos theta)) = 2 at 33.53 degrees from a particle, which
        Sfetsos prints as 33.52 and the History gives as 33.5."""
        sp = self.sp
        entry, reader = self.read("kruskal")
        theta = reader.symbol["\\theta"]
        g_uu = next(reader(e["value"]) for e in entry["metric_components"] if e["indices"] == ["u", "u"])
        profile = g_uu.coeff(sp.DiracDelta(reader.symbol["u"]))
        root = sp.nsolve(profile.subs({reader.parameters["G"]: 1, reader.parameters["E"]: 1, reader.c: 1}), theta, 0.6)
        self.assertAlmostEqual(math.degrees(float(root)), 33.534, places=2)

    def test_a_ray_crossing_the_shock_jumps_by_the_profile(self):
        """v'' = -Gamma^v_uu u'^2 integrated twice across u = 0: Delta v = (2GE/c^4) P(theta), half
        the coefficient of delta(u) du^2 over the 2 of -4 du dv."""
        sp = self.sp
        entry, reader = self.read("kruskal")
        u = reader.symbol["u"]
        gamma = next(reader(e["value"]) for e in entry["christoffel"]["variants"]["ull"]["nonzero"]
                     if e["indices"] == ["v", "u", "u"])
        g_uu = next(reader(e["value"]) for e in entry["metric_components"] if e["indices"] == ["u", "u"])
        jump = -sp.expand(gamma).coeff(sp.DiracDelta(u, 1))
        self.assertEqual(sp.simplify(jump - g_uu.coeff(sp.DiracDelta(u)) / 4), 0)
        at = {reader.parameters["G"]: 1, reader.parameters["E"]: sp.Rational(1, 8), reader.c: 1,
              reader.symbol["\\theta"]: sp.pi / 2}
        self.assertEqual(jump.subs(at), -sp.Rational(1, 2))

    def test_the_null_chart_reaches_aichelburg_and_sexl(self):
        """As a goes to infinity the conformal factor goes to 1 and the profile to
        -(4GE/c^4) ln(rho^2/rho_0^2) with rho_0 = 2a/e, Hotta and Tanaka's (21)."""
        sp = self.sp
        entry, reader = self.read("null_cylindrical")
        a, rho = reader.parameters["a"], reader.symbol["\\rho"]
        G, E, c = reader.parameters["G"], reader.parameters["E"], reader.c
        g = {tuple(e["indices"]): reader(e["value"]) for e in entry["metric_components"]}
        profile = g[("u", "u")].coeff(sp.DiracDelta(reader.symbol["u"]))
        shock = -4 * G * E / c ** 4 * sp.log(rho ** 2 * sp.E ** 2 / (4 * a ** 2))
        self.assertEqual(sp.limit(sp.expand_log(profile - shock, force=True), a, sp.oo), 0)
        self.assertEqual(sp.limit(g[("u", "v")], a, sp.oo), -sp.Rational(1, 2))
        self.assertEqual(sp.limit(g[("\\rho", "\\rho")], a, sp.oo), 1)


class Drawings(unittest.TestCase):
    def test_the_ring_follows_podolsky_and_ortaggio(self):
        """Each frame of the embedding's movie: the sphere's radius and the ring's latitude against
        Z = sin(theta_0) sqrt(1 + s^2) - g s/sin(theta_0) and
        Z_4 = cos(theta_0) sqrt(1 + s^2) + g ln(cot(theta_0/2)) s behind the shock, s = sinh(tau),
        at g = 4GE/(c^4 a) = 1/2 and theta_0 = 30 degrees, to twice the width of the drawn pulse."""
        view, = json.loads((DATA / "embedding" / "hotta_tanaka.json").read_text(encoding="utf-8"))["views"]
        theta0, g = math.pi / 6, 0.5
        frames = view["movie"]["frames"]
        self.assertGreater(len(frames), 20)
        latitudes = []
        for frame in frames:
            tau = frame["value"]
            s = math.sinh(tau)
            behind = max(s, 0.0)
            Z = math.sin(theta0) * math.sqrt(1 + s * s) - g / math.sin(theta0) * behind
            Z4 = math.cos(theta0) * math.sqrt(1 + s * s) + g * math.log(1 / math.tan(theta0 / 2)) * behind
            north = next(p for p in frame["pieces"] if p["id"] == "north")
            radius = north["points"][0][2]
            self.assertEqual(north["points"][0][1], 0)
            x, y, z = frame["curves"][0]["points"][0]
            self.assertAlmostEqual(radius, math.hypot(Z, Z4), delta=0.01)
            self.assertAlmostEqual(math.hypot(x, y), Z, delta=0.01)
            self.assertAlmostEqual(z, Z4, delta=0.01)
            latitudes.append(math.atan2(math.hypot(x, y), z))
        ahead = [lat for frame, lat in zip(frames, latitudes) if frame["value"] < -0.03]
        self.assertTrue(all(abs(lat - theta0) < 1e-6 for lat in ahead))
        after = [lat for frame, lat in zip(frames, latitudes) if frame["value"] > 0.03]
        self.assertTrue(all(b < a for a, b in zip(after, after[1:])), "behind the shock the ring slides toward the pole")
        self.assertLess(after[-1], 0.05)

    @unittest.skipUnless(HAS_NUMPY, "numpy is not installed")
    def test_each_moment_lies_on_its_sphere_in_every_chart(self):
        """slices.hotta_tanaka_event against each chart's map into the hyperboloid, written out here
        at a = 1: the event has the sphere's Z_0 + Z_1 and Z_0 - Z_1."""
        import slices

        def null_pair(system, x0, r, fixed):
            c, s = math.cos(fixed), math.sin(fixed)
            if system == "kruskal":
                return 2 * x0 / (1 - x0 * r), 2 * r / (1 - x0 * r)
            if system == "global":
                return (math.cos(r) - math.cos(x0)) / math.sin(x0), -(math.cos(x0) + math.cos(r)) / math.sin(x0)
            if system == "null_cylindrical":
                omega = 1 + (fixed ** 2 - x0 * r) / 4
                return x0 / omega, r / omega
            if system == "kundt":
                w, u = x0, r
                return u * (2 + u * w) * s / 2, -2 * c + 2 * w * s + u * s + u * u * w * s / 2
            eta, rho = x0, r
            return (eta * eta - rho * rho) / (2 * eta), (eta * eta - rho * rho - 4 + 4 * rho * c) / (2 * eta)
        met = 0
        for system, planes in (("kruskal", (math.pi / 2, math.pi / 12)), ("global", (math.pi / 2, math.pi / 12)),
                               ("conformally_flat", (math.pi / 2, math.pi / 12)), ("kundt", (math.pi / 2, math.pi / 12)),
                               ("null_cylindrical", (2.0, 0.25))):
            for fixed in planes:
                for tau in (-0.6, 0.0, 0.25, 0.5):
                    U, V = slices.hotta_tanaka_sphere(tau)
                    event = slices.hotta_tanaka_event(system, fixed, U, V)
                    if event is None:
                        continue
                    met += 1
                    found = null_pair(system, *event, fixed)
                    self.assertAlmostEqual(found[0], U, places=9, msg=f"{system} {fixed} {tau}")
                    self.assertAlmostEqual(found[1], V, places=9, msg=f"{system} {fixed} {tau}")
        self.assertEqual(met, 35)

    def test_the_conformal_diagram_is_taller_than_de_sitters_square(self):
        """Behind the shock future infinity, drawn with the map that keeps each crossing ray one
        line, stands above T = pi/2, and ahead of it the drawing is de Sitter's square."""
        view, = json.loads((DATA / "conformal" / "hotta_tanaka.json").read_text(encoding="utf-8"))["views"]
        scri = [layer["points"] for layer in view["layers"] if layer["class"] == "scri"]
        top = max(point[1] for line in scri for point in line)
        self.assertAlmostEqual(top, math.pi / 2 + math.atan(0.5), places=3)
        self.assertAlmostEqual(min(point[1] for line in scri for point in line), -math.pi / 2, places=4)
        self.assertGreater(view["box"][3], math.pi / 2 + math.atan(0.5))


if __name__ == "__main__":
    unittest.main()
