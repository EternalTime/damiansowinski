"""The naked black holes of Horowitz and Ross:
.venv.noindex/bin/python -m unittest discover -s _tools

What the published files state, held to the articles. Each chart's published line element is read
through the checker's own Reader at points and held to its source, HR97's (gensol), (uuf) and (R),
(uf) and (ur), or (nsmet) and (dd); its published Riemann tensor, turned into the frame of a radial
fall, to HR97's tidal component (bhtide), or for the string chart (gdevb); the History to stating
the line element and the areal radius the dilaton chart prints; and every drawing of the dilaton
chart to the one member it is drawn at. The tests that read the mathematics need sympy and are
skipped where it is absent.
"""
import importlib.util
import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "MFS" / "assets" / "data"
HAS_SYMPY = importlib.util.find_spec("sympy") is not None
sys.path.insert(0, str(Path(__file__).resolve().parent / "derivations"))
METRIC = json.loads((DATA / "metrics" / "naked_black_holes.json").read_text(encoding="utf-8"))
SYSTEMS = ("dilaton", "two_charges", "string_ns")


def chart(system):
    return next(c for c in METRIC["coordinates"] if c["id"] == system)


def drawing(kind):
    return json.loads((DATA / kind / "naked_black_holes.json").read_text(encoding="utf-8"))


@unittest.skipUnless(HAS_SYMPY, "needs sympy")
class Charts(unittest.TestCase):
    """Each published chart at points against its source."""

    @classmethod
    def setUpClass(cls):
        import sympy as sp
        import naked_black_holes as nbh
        import verify_metrics as vm
        cls.sp, cls.nbh, cls.vm = sp, nbh, vm

    def read(self, system):
        reader = self.nbh.reader(system)
        c = chart(system)
        g = self.vm.metric_from_line_element(reader, c["line_element"].replace("c^2", ""), c["coords"])
        return c, reader, g

    def number(self, reader, expr, at):
        expr = self.sp.sympify(expr).xreplace(reader.held) if reader.held else self.sp.sympify(expr)
        return complex(self.sp.N(expr.doit().subs(at), 30))

    def sources(self, reader):
        """F, R and the point of each chart, and the time's factor in g_tt where the chart is not of
        HR97's form (gensol)."""
        sp, P, r, th = self.sp, reader.parameters, reader.symbol["r"], reader.symbol["\\theta"]
        if "a" in P:
            a, rp, rm = P["a"], P["r_+"], P["r_-"]
            R = r * (1 - rm / r) ** (a ** 2 / (1 + a ** 2))
            return (r - rp) * (r - rm) / R ** 2, R, 1, {a: sp.Rational(3, 4), rp: 2, rm: 1, r: 3, th: sp.Rational(1, 2)}
        if "Sigma" in P:
            S, rp, rm = P["Sigma"], P["r_+"], P["r_-"]
            R = sp.sqrt(r ** 2 - S ** 2)
            return (r - rp) * (r - rm) / R ** 2, R, 1, {S: sp.Rational(1, 2), rp: 2, rm: 1, r: 3, th: 1}
        r0, r1, rq = P["r_0"], P["r_1"], P["r_p"]
        delta = (1 + r1 / r) * (1 + rq / r)
        return 1 - r0 / r, r, delta, {r0: 1, r1: sp.Rational(1, 3), rq: 2, r: sp.Rational(5, 2), th: 1}

    def test_each_line_element_is_its_sources(self):
        for system in SYSTEMS:
            with self.subTest(system=system):
                c, reader, g = self.read(system)
                F, R, G, at = self.sources(reader)
                th = reader.symbol["\\theta"]
                for i, value in enumerate((-F / G, 1 / F, R ** 2, R ** 2 * self.sp.sin(th) ** 2)):
                    self.assertLess(abs(self.number(reader, g[i, i] - value, at)), 1e-25)

    def test_the_published_riemann_tensor_gives_the_tidal_force_of_a_radial_fall(self):
        """R_0'20'2 = u^a u^c R_a2c2 / g_22 for u^t = EG/F, u^r = -(E^2 G - F)^{1/2}, which is HR97's
        -[R''(E^2 G - F) + R'(E^2 G' - F')/2]/R, (gdevb), and (bhtide) where G = 1."""
        sp = self.sp
        E = sp.Rational(3, 2)
        for system in SYSTEMS:
            with self.subTest(system=system):
                c, reader, g = self.read(system)
                F, R, G, at = self.sources(reader)
                r = reader.symbol["r"]
                comps = {tuple(e["indices"]): e["value"] for e in c["riemann"]["variants"]["llll"]["nonzero"]}
                value = {k: self.number(reader, reader(comps.get(k, "0")), at).real
                         for k in (("t", "\\theta", "t", "\\theta"), ("t", "\\theta", "r", "\\theta"),
                                   ("r", "\\theta", "r", "\\theta"))}
                ut = float((E * G / F).subs(at))
                ur = -float(sp.sqrt(E ** 2 * G - F).subs(at))
                g22 = self.number(reader, g[2, 2], at).real
                falling = (ut ** 2 * value[("t", "\\theta", "t", "\\theta")]
                           + 2 * ut * ur * value[("t", "\\theta", "r", "\\theta")]
                           + ur ** 2 * value[("r", "\\theta", "r", "\\theta")]) / g22
                tide = -(sp.diff(R, r, 2) * (E ** 2 * G - F) + sp.diff(R, r) * (E ** 2 * sp.diff(G, r) - sp.diff(F, r)) / 2) / R
                self.assertAlmostEqual(falling, float(sp.N(tide.subs(at), 30)), places=12)


class History(unittest.TestCase):
    def test_it_states_the_line_element_and_the_radius_the_dilaton_chart_prints(self):
        history = METRIC["history"]
        self.assertIn("$ds^2 = -F\\,dt^2 + dr^2/F + R^2 d\\Omega^2$", history)
        self.assertIn("$R = r(1 - r_-/r)^{a^2/(1+a^2)}$", history)
        self.assertIn("$R_{0'20'2} = -[R''(E^2 - F) - F'R'/2]/R$", history)
        self.assertEqual(METRIC["references"][0], "horowitz1997naked")


class Drawings(unittest.TestCase):
    """Every drawing of the dilaton chart is of the member a = 1/2, r_- = r_+/2, the one the embedding
    diagram draws and whose moment the other drawings carry."""

    def test_every_drawing_of_the_dilaton_chart_is_of_one_member(self):
        settings = "$a = 1/2$, $r_- = r_+/2$, $r_+ = 1$."
        embedding = drawing("embedding")
        self.assertEqual([v["settings"] for v in embedding["views"]], [settings])
        conformal = next(v for v in drawing("conformal")["views"] if v["id"] == "dilaton")
        self.assertEqual(conformal["settings"], settings)
        figure = drawing("curvature")["views"]
        self.assertEqual([(v["system"], v["settings"]) for v in figure], [("dilaton", settings)])
        diagrams = drawing("diagrams")
        self.assertEqual([v["settings"].split(", $\\theta")[0] for v in diagrams["systems"]["dilaton"]],
                         ["$r_+ = 1$, $r_- = \\frac{1}{2}$, $a = \\frac{1}{2}$"])


if __name__ == "__main__":
    unittest.main()
