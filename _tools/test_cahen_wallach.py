"""The symmetric plane waves of Cahen and Wallach:
.venv.noindex/bin/python -m unittest discover -s _tools

What the page and its drawings state of the physics, held on the published files: every chart has
four dimensions and vanishing Ricci and Kretschmann scalars; the Brinkmann chart's one Ricci
component is R_uu = -(lambda_1 + lambda_2), so the vacuum member is lambda_1 + lambda_2 = 0, and its
Weyl tensor vanishes at lambda_1 = lambda_2; the rotating chart is the conformally flat member, with
R_uu = 2 omega^2 and no Weyl tensor, and the Rosen chart the vacuum member, with no Ricci tensor; the
Riemann tensor of the Brinkmann chart is covariantly constant; light leaving the axis where the wave
pulls meets again at cu = pi L; and the drawings are the ones the captions name. The tests of the
published mathematics need sympy and are skipped where it is absent; the tests of the drawings read
the files alone.
"""
import importlib.util
import itertools
import json
import random
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "MFS" / "assets" / "data"
HAS_SYMPY = importlib.util.find_spec("sympy") is not None
sys.path.insert(0, str(Path(__file__).resolve().parent / "derivations"))
ID = "cahen_wallach"
CHARTS = ["brinkmann", "rotating", "rosen"]


def published():
    return json.loads((DATA / "metrics" / f"{ID}.json").read_text(encoding="utf-8"))


def drawing(kind):
    return json.loads((DATA / kind / f"{ID}.json").read_text(encoding="utf-8"))


def chart(system):
    return next(c for c in published()["coordinates"] if c["id"] == system)


class Files(unittest.TestCase):
    def test_three_charts_of_four_dimensions(self):
        self.assertEqual([c["id"] for c in published()["coordinates"]], CHARTS)
        for system in CHARTS:
            self.assertEqual(len(chart(system)["coords"]), 4, system)

    def test_every_scalar_vanishes(self):
        for system in CHARTS:
            self.assertEqual(chart(system)["ricci_scalar"], "R = 0", system)
            self.assertEqual(chart(system)["kretschmann"], "K = 0", system)

    def test_the_ricci_tensor_lies_along_the_wave(self):
        want = {"brinkmann": "-\\left(\\lambda_1 + \\lambda_2\\right)", "rotating": "2\\omega^2"}
        for system in CHARTS:
            nonzero = chart(system)["ricci_tensor"]["variants"]["ll"]["nonzero"]
            if system == "rosen":
                self.assertEqual(nonzero, [], "the Rosen chart is the vacuum member")
            else:
                self.assertEqual([(e["indices"], e["value"]) for e in nonzero], [(["u", "u"], want[system])])

    def test_the_rotating_chart_is_conformally_flat_and_the_rosen_chart_is_not(self):
        self.assertEqual(chart("rotating")["weyl_tensor"]["variants"]["llll"]["nonzero"], [])
        self.assertNotEqual(chart("rosen")["weyl_tensor"]["variants"]["llll"]["nonzero"], [])

    def test_the_drawings_are_the_ones_named(self):
        diagrams = drawing("diagrams")
        self.assertFalse(diagrams["systems"])
        self.assertEqual({s: [v["id"] for v in vs] for s, vs in diagrams["projections"].items()},
                         {"brinkmann": ["push", "pull"], "rosen": ["push", "pull"], "rotating": ["turning"]})
        self.assertEqual([v["id"] for v in drawing("conformal")["views"]], ["brinkmann", "rotating"])
        self.assertEqual([v["id"] for v in drawing("embedding")["views"]], ["tube", "ring"])
        self.assertFalse((DATA / "curvature" / f"{ID}.json").exists(), "both stored scalars vanish")

    def test_the_ring_is_drawn_at_its_six_moments(self):
        ring = next(v for v in drawing("embedding")["views"] if v["id"] == "ring")
        self.assertEqual([s["time"] for s in ring["surfaces"]], [0.0, 0.3, 0.6, 0.9, 1.2, 1.5])


@unittest.skipUnless(HAS_SYMPY, "needs sympy")
class Mathematics(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import sympy as sp
        import chart_printer as cp
        import print_charts as pc
        cls.sp, cls.cp, cls.pc = sp, cp, pc
        cls.brink = pc.cw_chart("brinkmann")

    def at(self, symbols, seed):
        rng = random.Random(seed)
        return {s: self.sp.Rational(rng.choice([-1, 1]) * rng.randint(3, 40), rng.randint(7, 19)) for s in symbols}

    def test_the_riemann_tensor_is_covariantly_constant(self):
        geo, x = self.brink.geo, self.brink.symbols
        params = list(self.brink.reader.parameters.values())
        gamma, riemann = geo.christoffel_ull(), geo.riemann_llll()
        point = self.at(x + params, 1970)
        for e, a, b, c, d in itertools.product(range(4), repeat=5):
            value = self.sp.diff(riemann[a][b][c][d], x[e])
            for f in range(4):
                value -= (gamma[f][e][a] * riemann[f][b][c][d] + gamma[f][e][b] * riemann[a][f][c][d]
                          + gamma[f][e][c] * riemann[a][b][f][d] + gamma[f][e][d] * riemann[a][b][c][f])
            self.assertEqual(self.sp.sympify(value).xreplace(point), 0, (e, a, b, c, d))

    def test_light_from_the_axis_meets_again_where_the_wave_pulls(self):
        # y'' = lambda_2 y u'^2 with u affine: at lambda_2 = -1/L^2 the ray y = s L sin(cu/L) is back on the
        # axis at cu = pi L, and at lambda_1 = 1/L^2 the ray x = s L sinh(cu/L) never is.
        sp = self.sp
        u, s = sp.symbols("u s", positive=True)
        lam1, lam2 = self.brink.reader.parameters.values()
        gamma = self.brink.geo.christoffel_ull()
        x, y = self.brink.symbols[2], self.brink.symbols[3]
        for k, lam, path in ((2, lam1, s * sp.sinh(u)), (3, lam2, s * sp.sin(u))):
            other = {lam1: 1, lam2: -1}
            on = {x: path, y: 0} if k == 2 else {x: 0, y: path}
            equation = sp.diff(path, u, 2) + gamma[k][0][0].subs(other).subs(on)
            self.assertEqual(sp.expand(equation), 0)
        self.assertEqual((s * sp.sin(u)).subs(u, sp.pi), 0)
        self.assertNotEqual((s * sp.sinh(u)).subs(u, sp.pi), 0)

    def test_the_rosen_and_rotating_charts_are_the_brinkmann_chart_pulled_back(self):
        for system in ("rotating", "rosen"):
            self.pc.cw_check(self.pc.cw_chart(system), system)
