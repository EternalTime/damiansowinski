"""Kopczynski and Trautman's universe with torsion: python3 -m unittest discover -s _tools

What the texts and drawings of `kopczynski_trautman` state, held on the published files. The
comoving charts are held to the modified Friedmann equation of Trautman's article, to the
Einstein tensor of dust plus the repulsion of its spins, and to the curvature the captions state
at the turn; the conformal chart to the same Einstein tensor where its scale factor obeys the
equation its parameter states. The numbers of the history, Trautman's centimetre and the density
there, are computed from the least radius. The spacetime diagrams are held to the Hubble sphere
their caption places and the conformal diagram to filling Minkowski's half diamond. The tests that
read the published components need sympy and are skipped where it is absent.
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

COMOVING = ("comoving_cartesian", "comoving_spherical")


def published(system):
    metric = json.loads((DATA / "metrics" / "kopczynski_trautman.json").read_text(encoding="utf-8"))
    return next(c for c in metric["coordinates"] if c["id"] == system)


def drawn(kind):
    return json.loads((DATA / kind / "kopczynski_trautman.json").read_text(encoding="utf-8"))


@unittest.skipUnless(HAS_SYMPY, "sympy is not installed")
class Theory(unittest.TestCase):
    """The published components, read as the checker reads them, in the chart x^0 = ct with c = 1."""

    def read(self, system):
        import sympy as sp
        import verify_metrics as vm
        chart = published(system)
        reader = vm.Reader(chart["coords"], [p["symbol"] for p in chart["parameters"]], ())

        def tensor(field, variant=None):
            block = chart[field] if variant is None else chart[field]["variants"][variant]["nonzero"]
            out = sp.zeros(4, 4)
            for entry in block:
                i, j = (chart["coords"].index(x) for x in entry["indices"])
                out[i, j] = reader.surface(reader(entry["value"])).subs(reader.c, 1)
            return out
        return chart, reader, tensor

    def test_the_scale_factor_solves_trautmans_equation(self):
        import sympy as sp
        for system in COMOVING:
            chart, reader, tensor = self.read(system)
            t, ell = reader.symbol["t"], reader.parameters["ell"]
            a = sp.sqrt(tensor("metric_components")[1, 1])
            # (da/d(ct))^2 = (4/9 l^2)(1/a - 1/a^4), Trautman's (33) counted from the least radius.
            self.assertEqual(sp.simplify(sp.diff(a, t) ** 2 - 4 * (1 / a - 1 / a ** 4) / (9 * ell ** 2)), 0, system)
            self.assertEqual(sp.simplify(a.subs(t, 0)), 1, system)
            # Far from the turn it is Friedmann's flat dust, a = (ct/l)^(2/3).
            self.assertEqual(sp.limit((a / t ** sp.Rational(2, 3)).subs(ell, 1), t, sp.oo), 1, system)

    def test_the_einstein_tensor_is_dust_and_the_repulsion_of_its_spins(self):
        import sympy as sp
        for system in COMOVING:
            chart, reader, tensor = self.read(system)
            ell = reader.parameters["ell"]
            a = sp.sqrt(tensor("metric_components")[1, 1])
            G = tensor("einstein_tensor", "ul")
            density, pressure = 4 * (a ** -3 - a ** -6) / (3 * ell ** 2), -4 * a ** -6 / (3 * ell ** 2)
            self.assertEqual(sp.simplify(G[0, 0] + density), 0, system)
            for i in range(1, 4):
                self.assertEqual(sp.simplify(G[i, i] - pressure), 0, system)
            # The strong energy condition fails near the turn: density + 3 pressure < 0 while a^3 < 4.
            self.assertEqual(sp.simplify((density + 3 * pressure) * 3 * ell ** 2 * a ** 6 / 4 - (a ** 3 - 4)), 0)

    def test_the_curvature_is_finite_and_greatest_at_the_turn(self):
        import sympy as sp
        for system in COMOVING:
            chart, reader, _ = self.read(system)
            t, ell = reader.symbol["t"], reader.parameters["ell"]
            K = reader.surface(reader(chart["kretschmann"].partition("=")[2])).subs(reader.c, 1)
            self.assertEqual(sp.simplify(K.subs(t, 0) - sp.Rational(16, 3) / ell ** 4), 0, system)
            # K = 16(5u^2 - 16u + 20)/27 l^4 u^4 with u = a^3 >= 1 falls as u grows.
            u = sp.Symbol("u", positive=True)
            slope = sp.diff((5 * u ** 2 - 16 * u + 20) / u ** 4, u)
            self.assertTrue(all(slope.subs(u, value) < 0 for value in (1, 2, 5, 50)))

    def test_the_conformal_chart_on_its_equation(self):
        import sympy as sp
        chart, reader, tensor = self.read("conformal")
        eta, a = reader.symbol["\\eta"], reader.parameters["a"]
        G = tensor("einstein_tensor", "ul")
        first, second = 4 * (a - a ** -2) / 9, 2 * (1 + 2 / a ** 3) / 9
        on = lambda e: sp.simplify(e.subs(sp.Derivative(a, (eta, 2)), second).subs(sp.Derivative(a, eta) ** 2, first))  # noqa: E731
        self.assertEqual(on(G[0, 0] + 4 * (a ** -3 - a ** -6) / 3), 0)
        self.assertEqual(on(G[1, 1] + 4 * a ** -6 / 3), 0)


class Numbers(unittest.TestCase):
    """The figures the history and the captions state."""

    def test_trautmans_centimetre(self):
        G, c, hbar, m, N = 6.674e-11, 2.998e8, 1.0546e-34, 1.6726e-27, 1e80
        cube = 3 * G * (N * hbar / 2) ** 2 / (2 * N * m * c ** 4)
        radius = cube ** (1 / 3)
        density = N * m / (4 / 3 * math.pi * cube)
        self.assertAlmostEqual(radius * 100, 1.27, places=2)
        self.assertEqual(round(math.log10(density * 1e-3)), 55)
        self.assertEqual(round(math.log10(c ** 5 / (hbar * G ** 2) / density)), 38)
        # The turn is where rho c^4 = 2 pi G sigma^2.
        sigma = (N * hbar / 2) / (4 / 3 * math.pi * cube)
        self.assertAlmostEqual(density * c ** 4 / (2 * math.pi * G * sigma ** 2), 1.0, places=12)

    def test_the_hubble_sphere_comes_nearest_at_root_three(self):
        def r(t):
            return 3 * (1 + t * t) ** (2 / 3) / (2 * abs(t))
        best = min((r(0.001 * k), 0.001 * k) for k in range(1, 6000))
        self.assertAlmostEqual(best[1], math.sqrt(3), places=2)
        self.assertAlmostEqual(best[0], 2.18, places=2)
        self.assertAlmostEqual(best[0], math.sqrt(3) / 2 * 4 ** (2 / 3), places=6)


class Drawings(unittest.TestCase):
    def test_every_chart_has_a_square_spacetime_diagram(self):
        systems = drawn("diagrams")["systems"]
        self.assertEqual(sorted(systems), sorted(COMOVING + ("conformal",)))

    def test_the_conformal_diagram_fills_the_half_diamond(self):
        views = drawn("conformal")["views"]
        self.assertEqual([v["id"] for v in views], ["comoving", "conformal"])
        for view in views:
            region = next(layer for layer in view["layers"] if layer["kind"] == "fill" and layer["class"] == "region")
            self.assertEqual([[round(x, 3), round(y, 3)] for x, y in region["points"]],
                             [[0, -3.142], [3.142, 0], [0, 3.142]])
            self.assertFalse([layer for layer in view["layers"] if layer["class"] == "singular"])

    def test_the_embedded_disc_is_smallest_at_the_turn(self):
        view, = drawn("embedding")["views"]
        radii = {}
        for surface in view["surfaces"]:
            piece, = surface["pieces"]
            radii[surface["time"]] = max(point[1] for point in piece["points"])
            # A flat plane: the surface stays at one height.
            self.assertLess(max(abs(point[2]) for point in piece["points"]), 1e-9)
        self.assertEqual(sorted(radii), [-3.0, -1.5, 0.0, 1.5, 3.0])
        for t, radius in radii.items():
            self.assertAlmostEqual(radius, (1 + t * t) ** (1 / 3), places=6)
        self.assertEqual(min(radii, key=radii.get), 0.0)


if __name__ == "__main__":
    unittest.main()
