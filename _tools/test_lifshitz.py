"""Lifshitz spacetime, anti-de Sitter space with time scaling as a power of space:
.venv.noindex/bin/python -m unittest discover -s _tools

What its texts state, held on the published files: the anisotropic scaling is an isometry of
every chart's published metric, the published Einstein tensor is the stress of Taylor's massive
vector field with Kachru, Liu and Mulligan's cosmological constant, the null energy condition
holds exactly for z >= 1, at z = 1 the inverse radius chart is anti-de Sitter space's published
Poincare chart and the Weyl tensor vanishes, the tidal force on a falling observer grows as
r^(-2z), and the drawings hold what their captions say: each ray of the figure a catenary that
turns back at u = 2L sin(alpha), and the embedded moment a pseudosphere below r = L. The tests
of the published mathematics need sympy and are skipped where it is absent; the tests of the
drawings read the files alone.
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


def chart(system):
    metric = json.loads((DATA / "metrics" / "lifshitz_spacetime.json").read_text(encoding="utf-8"))
    return next(c for c in metric["coordinates"] if c["id"] == system)


@unittest.skipUnless(HAS_SYMPY, "sympy is not installed")
class Published(unittest.TestCase):
    def setUp(self):
        import sympy
        import verify_metrics
        self.sp, self.vm = sympy, verify_metrics

    def read(self, system, metric_id="lifshitz_spacetime"):
        """The reader of a published chart, its coordinates, and a function from a field's name to
        the matrix of its published components, with c = 1."""
        metric = json.loads((DATA / "metrics" / f"{metric_id}.json").read_text(encoding="utf-8"))
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
        return reader, x, matrix, entry

    def zero(self, expression, at):
        self.assertLess(abs(complex(self.sp.sympify(expression).subs(at).evalf(30))), 1e-20)

    def point(self, reader, x, z=None):
        """A point of a chart with positive coordinates, at a dynamical exponent that is no integer."""
        at = {s: self.sp.Rational(7 + 2 * i, 5 + i) for i, s in enumerate(x)}
        at[reader.parameters["z"]] = self.sp.Rational(7, 3) if z is None else z
        at[reader.parameters["L"]] = self.sp.Rational(3, 2)
        return at

    def test_the_anisotropic_scaling_is_an_isometry_of_every_chart(self):
        # t -> lambda^z t, x -> lambda x, y -> lambda y and r -> r/lambda, Kachru, Liu and Mulligan's
        # (2.2), as a vector field in each chart: the Lie derivative of the published metric vanishes.
        sp = self.sp
        for system, radial in (("kachru_liu_mulligan", lambda r, z, L: -r), ("poincare", lambda u, z, L: u),
                               ("proper_distance", lambda rho, z, L: -L), ("tortoise", lambda w, z, L: z * w),
                               ("eddington_finkelstein", lambda r, z, L: -r), ("affine", lambda s, z, L: -z * s)):
            reader, x, matrix, _ = self.read(system)
            z, L = reader.parameters["z"], reader.parameters["L"]
            g = matrix("metric_components")
            xi = [z * x[0], x[1], x[2], radial(x[3], z, L)]
            at = self.point(reader, x)
            for a in range(4):
                for b in range(a, 4):
                    lie = (sum(xi[c] * sp.diff(g[a, b], x[c]) for c in range(4))
                           + sum(g[c, b] * sp.diff(xi[c], x[a]) + g[a, c] * sp.diff(xi[c], x[b]) for c in range(4)))
                    with self.subTest(system=system, slot=(a, b)):
                        self.zero(lie, at)

    def test_the_einstein_tensor_is_the_massive_vector_fields_stress(self):
        # G_ab + Lambda g_ab = (F_ac F_b^c - g_ab F^2/4)/2 + m^2 (A_a A_b - g_ab A^2/2)/2 with
        # A = sqrt(2(z - 1)/z) (r/L)^z c dt, m^2 = 2z/L^2 and Lambda = -(z^2 + z + 4)/(2L^2).
        sp = self.sp
        reader, x, matrix, _ = self.read("kachru_liu_mulligan")
        z, L, r = reader.parameters["z"], reader.parameters["L"], x[3]
        g, gi, G = matrix("metric_components"), matrix("inverse_metric_components"), matrix("einstein_tensor", "ll")
        A = [sp.sqrt(2 * (z - 1) / z) * (r / L) ** z, 0, 0, 0]
        F = sp.Matrix(4, 4, lambda a, b: sp.diff(A[b], x[a]) - sp.diff(A[a], x[b]))
        F2 = sum(F[a, b] * (gi * F * gi)[a, b] for a in range(4) for b in range(4))
        A2 = gi[0, 0] * A[0] ** 2
        at = self.point(reader, x)
        for a in range(4):
            for b in range(4):
                stress = ((sum(F[a, c] * F[b, d] * gi[c, d] for c in range(4) for d in range(4)) - g[a, b] * F2 / 4) / 2
                          + (2 * z / L ** 2) * (A[a] * A[b] - g[a, b] * A2 / 2) / 2)
                with self.subTest(slot=(a, b)):
                    self.zero(G[a, b] - (z ** 2 + z + 4) / (2 * L ** 2) * g[a, b] - stress, at)

    def test_the_null_energy_condition_holds_exactly_for_z_of_one_or_more(self):
        # G_ab k^a k^b on the null vectors along r and along x: 2(z - 1)/L^2 and (z - 1)(z + 2)/L^2 times
        # the square of k^t in the static frame, both of the sign of z - 1.
        reader, x, matrix, _ = self.read("kachru_liu_mulligan")
        z, L = reader.parameters["z"], reader.parameters["L"]
        mixed = matrix("einstein_tensor", "ul") if "ul" in chart("kachru_liu_mulligan")["einstein_tensor"]["variants"] else \
            matrix("inverse_metric_components") * matrix("einstein_tensor", "ll")
        at = self.point(reader, x)
        self.zero(mixed[3, 3] - mixed[0, 0] - 2 * (z - 1) / L ** 2, at)
        self.zero(mixed[1, 1] - mixed[0, 0] - (z - 1) * (z + 2) / L ** 2, at)
        for value, sign in ((self.sp.Rational(1, 2), -1), (1, 0), (2, 1)):
            for combination in (mixed[3, 3] - mixed[0, 0], mixed[1, 1] - mixed[0, 0]):
                got = float(combination.subs(self.point(reader, x, value)))
                self.assertEqual((got > 1e-12) - (got < -1e-12), sign)

    def test_at_z_of_one_it_is_anti_de_sitter_space(self):
        reader, x, matrix, entry = self.read("poincare")
        other, y, there, _ = self.read("poincare", "anti_de_sitter")
        at = self.point(reader, x, 1)
        ads = there("metric_components").subs(dict(zip(y, x))).subs(other.parameters["L"], reader.parameters["L"])
        mine = matrix("metric_components")
        for a in range(4):
            for b in range(4):
                self.zero(mine[a, b] - ads[a, b], at)
        for e in entry["weyl_tensor"]["variants"]["llll"]["nonzero"]:
            self.zero(reader(e["value"]), at)
        # And for no other exponent: the Weyl tensor is z(z - 1) times a tensor that does not vanish.
        some = entry["weyl_tensor"]["variants"]["llll"]["nonzero"][0]
        self.assertGreater(abs(float(reader(some["value"]).subs(self.point(reader, x)))), 1e-3)

    def test_the_tidal_force_on_a_falling_observer_grows_as_a_power_of_r(self):
        # Along the timelike geodesic of energy E with no momentum along x or y, Copsey and Mann's (2.12):
        # the tidal tensor's component along x in the falling frame is (1 + (z - 1) E^2 (L/r)^(2z))/L^2.
        sp = self.sp
        reader, x, matrix, entry = self.read("kachru_liu_mulligan")
        z, L, r = reader.parameters["z"], reader.parameters["L"], x[3]
        names = entry["coords"]
        riemann = {tuple(names.index(k) for k in e["indices"]): reader(e["value"]).subs(reader.c, 1)
                   for e in entry["riemann"]["variants"]["llll"]["nonzero"]}
        E = sp.Rational(5, 4)
        u = [E * (L / r) ** (2 * z), 0, 0, -(r / L) * sp.sqrt(E ** 2 * (L / r) ** (2 * z) - 1)]
        g = matrix("metric_components")
        at = self.point(reader, x)
        at[r] = sp.Rational(1, 3)
        self.zero(sum(g[a, a] * u[a] ** 2 for a in range(4)) + 1, at)
        tidal = sum(value * u[c] * u[d] for (a, c, b, d), value in riemann.items() if (a, b) == (1, 1)) / g[1, 1]
        self.zero(tidal - (1 + (z - 1) * E ** 2 * (L / r) ** (2 * z)) / L ** 2, at)


class Drawings(unittest.TestCase):
    def test_every_ray_of_the_figure_turns_back_where_its_caption_says(self):
        # The page's own coordinates are x across and -u up. A ray at z = 2 that leaves u = 2L at the
        # angle alpha is the catenary u = u_0 cosh((x - x_0)/u_0) with u_0 = 2L sin(alpha), and the one
        # sent straight up reaches the boundary; at z = 1 every ray is a straight line that arrives.
        diagrams = json.loads((DATA / "diagrams" / "lifshitz_spacetime.json").read_text(encoding="utf-8"))
        figure, = diagrams["projections"]["poincare"]
        rays = [[(p[0], -p[1]) for p in layer["points"]] for layer in figure["layers"] if layer["class"] == "above"]
        straight = [[(p[0], -p[1]) for p in layer["points"]] for layer in figure["layers"] if layer["class"] == "below"]
        self.assertEqual((len(rays), len(straight)), (9, 8))
        turned = []
        for ray in rays:
            self.assertAlmostEqual(ray[0][0], 0.0, places=6)
            self.assertAlmostEqual(ray[0][1], 2.0, places=6)
            u0 = min(u for _, u in ray)
            if all(abs(x) < 1e-6 for x, _ in ray):
                self.assertLess(u0, 2e-3)
                continue
            side = 1 if ray[-1][0] > 0 else -1
            x0 = side * u0 * math.acosh(2.0 / u0)
            for x, u in ray:
                self.assertAlmostEqual(u, u0 * math.cosh((x - x0) / u0), delta=2e-3)
            self.assertGreater(ray[-1][1], 2.0)
            turned.append(round(math.degrees(math.asin(u0 / 2.0))))
        self.assertEqual(sorted(turned), [15, 15, 30, 30, 45, 45, 60, 60])
        for ray in straight:
            (xa, ua), (xb, ub) = ray[0], ray[-1]
            self.assertLess(ub, 2e-3)
            for x, u in ray:
                self.assertAlmostEqual((x - xa) * (ub - ua), (u - ua) * (xb - xa), delta=2e-3)

    def test_the_embedded_moment_is_a_pseudosphere_below_r_of_L(self):
        # Below rho = 0 the circle of radius e^rho stands at the height of the tractrix,
        # -(arcosh(e^-rho) - sqrt(1 - e^(2 rho))), and above it in Minkowski space at
        # sqrt(e^(2 rho) - 1) - arctan(sqrt(e^(2 rho) - 1)).
        embedding = json.loads((DATA / "embedding" / "lifshitz_spacetime.json").read_text(encoding="utf-8"))
        view, = embedding["views"]
        inner, outer = view["surfaces"][0]["pieces"]
        for rho, radius, height in inner["points"]:
            self.assertAlmostEqual(radius, math.exp(rho), delta=1e-6)
            self.assertAlmostEqual(height, -(math.acosh(math.exp(-rho)) - math.sqrt(max(0.0, 1 - math.exp(2 * rho)))), delta=1e-5)
        for rho, radius, height in outer["points"]:
            root = math.sqrt(max(0.0, math.expm1(2 * rho)))
            self.assertAlmostEqual(radius, math.exp(rho), delta=1e-6)
            self.assertAlmostEqual(height, root - math.atan(root), delta=1e-5)
        self.assertAlmostEqual(inner["points"][-1][1], 1.0, delta=1e-9)
        self.assertAlmostEqual(outer["points"][-1][1], 3.0, delta=1e-6)


if __name__ == "__main__":
    unittest.main()
