"""Witten's black hole in two dimensions, the first spacetime of two dimensions in the collection:
.venv.noindex/bin/python -m unittest discover -s _tools

What its dimension alone fixes, held on the checker's own Geometry: the Weyl and Einstein tensors
vanish and the Kretschmann scalar is the square of the Ricci scalar. And what its texts state, held
on the published files: the surface gravity, the longest proper time inside the horizon, and that
nothing is drawn beyond the singularity UV = 1 of the Kruskal plane. The first tests need sympy and
are skipped where it is absent; the rest read the files alone.
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


@unittest.skipUnless(HAS_SYMPY, "sympy is not installed")
class TwoDimensions(unittest.TestCase):
    def setUp(self):
        import sympy
        import verify_metrics
        self.sp, self.vm = sympy, verify_metrics

    def geometry(self, system):
        metric = json.loads((DATA / "metrics" / "witten_black_hole.json").read_text(encoding="utf-8"))
        entry = next(c for c in metric["coordinates"] if c["id"] == system)
        declared = self.vm.DIMENSIONS[("witten_black_hole", system)]
        reader = self.vm.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]],
                                self.vm.time_coordinates(declared, entry["coords"]))
        g = self.vm.metric_from_line_element(reader, entry["line_element"], entry["coords"])
        return reader, self.vm.Geometry(g, [reader.symbol[name] for name in entry["coords"]], 60)

    def test_the_weyl_and_einstein_tensors_of_two_dimensions_vanish(self):
        for system in ("witten", "kruskal", "eddington_finkelstein_ingoing"):
            _, geo = self.geometry(system)
            for tensor, rank in ((geo.weyl_llll(), 4), (geo.einstein_ll(), 2)):
                for index in self.vm._indices(2, rank):
                    self.assertEqual(self.vm._at(tensor, index), 0, f"{system} {index}")

    def test_the_kretschmann_scalar_is_the_square_of_the_ricci_scalar(self):
        for system in ("witten", "schwarzschild_gauge", "dilaton", "conformal", "kruskal"):
            _, geo = self.geometry(system)
            self.assertEqual(self.vm.norm(geo.kretschmann() - geo.ricci_scalar() ** 2), 0, system)

    def test_the_riemann_tensor_is_all_trace(self):
        # R_abcd = (R/2)(g_ac g_bd - g_ad g_bc), which is what leaves no Weyl tensor.
        _, geo = self.geometry("schwarzschild_gauge")
        riemann, R, g = geo.riemann_llll(), geo.ricci_scalar(), geo.g
        for a, b, c, d in self.vm._indices(2, 4):
            trace = R * (g[a, c] * g[b, d] - g[a, d] * g[b, c]) / 2
            self.assertEqual(self.vm.norm(riemann[a][b][c][d] - trace), 0, (a, b, c, d))

    def test_the_surface_gravity_is_lambda_for_every_mass(self):
        # f = 1 - m e^(-2 lambda x) vanishes at x_h = ln(m)/2 lambda, and f'(x_h)/2 = lambda.
        reader, geo = self.geometry("schwarzschild_gauge")
        x, lam, m = reader.symbol["x"], reader.parameters["lambda"], reader.parameters["m"]
        f = -geo.g[0, 0] / reader.c ** 2
        horizon = self.sp.log(m) / (2 * lam)
        self.assertEqual(self.sp.simplify(f.subs(x, horizon)), 0)
        self.assertEqual(self.sp.simplify(self.sp.diff(f, x).subs(x, horizon) / 2 - lam), 0)


class TheTimeInside(unittest.TestCase):
    def test_the_proper_time_from_the_horizon_to_the_singularity_is_pi_over_two_lambda(self):
        # Along constant t inside the horizon d tau = dx/sqrt(m e^(-2 lambda x) - 1), from the horizon
        # x_h = ln(m)/2 lambda down to x -> -infinity. With x = x_h - q^2 the integrand is
        # 2q/sqrt(e^(2 q^2) - 1), which tends to sqrt 2 at the horizon and falls off as 2q e^(-q^2);
        # Simpson's rule in q, at lambda = 1 and m = 3.
        m, n, far = 3.0, 4000, 8.0
        horizon = math.log(m) / 2

        def rate(q):
            if q == 0:
                return math.sqrt(2)
            x = horizon - q * q
            return 2 * q / math.sqrt(m * math.exp(-2 * x) - 1)
        h = far / n
        values = [rate(k * h) for k in range(n + 1)]
        total = h / 3 * (values[0] + values[-1] + 4 * sum(values[1:-1:2]) + 2 * sum(values[2:-1:2]))
        self.assertAlmostEqual(total, math.pi / 2, places=6)


class TheKruskalPlane(unittest.TestCase):
    def setUp(self):
        diagrams = json.loads((DATA / "diagrams" / "witten_black_hole.json").read_text(encoding="utf-8"))
        (self.view,) = diagrams["systems"]["kruskal"]

    def chart(self, u):
        """A point of the unit square as (U, V): the view draws (V - U)/2 across and (U + V)/2 up."""
        X0, X1, Y0, Y1 = self.view["box"]
        X, Y = X0 + u[0] * (X1 - X0), Y0 + u[1] * (Y1 - Y0)
        return Y - X, Y + X

    def test_no_ray_and_no_cone_lies_beyond_the_singularity(self):
        points = [u for family in self.view["rays"].values() for line in family for u in line]
        points += [cone["at"] for cone in self.view["cones"]]
        self.assertGreater(len(points), 100)
        for u in points:
            U, V = self.chart(u)
            self.assertLess(U * V, 1 + 2e-2, u)

    def test_the_singularity_is_drawn_on_the_hyperbola(self):
        (singular,) = [m for m in self.view["markers"] if m["kind"] == "singular"]
        self.assertEqual(len(singular["lines"]), 2)
        for line in singular["lines"]:
            for u in line:
                U, V = self.chart(u)
                self.assertAlmostEqual(U * V, 1, delta=2e-2, msg=str(u))

    def test_every_ray_keeps_its_null_coordinate(self):
        # Only g_UV is nonzero, so a ray moving left keeps V and one moving right keeps U.
        for name, keeps in (("P", 1), ("M", 0)):
            for line in self.view["rays"][name]:
                held = [self.chart(u)[keeps] for u in line]
                self.assertLess(max(held) - min(held), 2e-3, name)


class TheCigar(unittest.TestCase):
    def test_the_profile_is_the_closed_form(self):
        # rho = tanh r and z = arsinh(cosh r) - sqrt(1 + sech^2 r) + sqrt 2 - arsinh 1, at lambda = 1.
        embedding = json.loads((DATA / "embedding" / "witten_black_hole.json").read_text(encoding="utf-8"))
        (view,) = embedding["views"]
        (piece,) = view["surfaces"][0]["pieces"]
        self.assertEqual((piece["system"], piece["coordinate"]), ("witten", "r"))
        for r, rho, z in piece["points"]:
            self.assertAlmostEqual(rho, math.tanh(r), places=6)
            want = math.asinh(math.cosh(r)) - math.sqrt(1 + 1 / math.cosh(r) ** 2) + math.sqrt(2) - math.asinh(1)
            self.assertAlmostEqual(z, want, places=5)
        self.assertEqual(piece["points"][0][:2], [0, 0])


if __name__ == "__main__":
    unittest.main()
