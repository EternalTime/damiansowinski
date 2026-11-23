"""The Koutras-McIntosh spacetimes: .venv.noindex/bin/python -m unittest discover -s _tools

What the History, the charts and the drawings state, held on the published files. The first class
holds the curvature the History quotes, the second the relations between the charts, the third the
drawings and the two absences koutras_mcintosh.md gives the reasons for.
"""
import ast
import json
import sys
import unittest
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "MFS" / "assets" / "data"
ME = "koutras_mcintosh"
sys.path.insert(0, str(ROOT / "_tools" / "derivations"))
import verify_metrics as vm  # noqa: E402


def load(folder, name=ME):
    return json.loads((DATA / folder / f"{name}.json").read_text(encoding="utf-8"))


def chart(system):
    return next(c for c in load("metrics")["coordinates"] if c["id"] == system)


def matrix(system):
    """A published chart's metric as a sympy matrix, with its reader."""
    c = chart(system)
    reader = vm.Reader(c["coords"], [p["symbol"] for p in c["parameters"]], ())
    values = {tuple(e["indices"]): e["value"] for e in c["metric_components"]}
    n = len(c["coords"])
    return reader, sp.Matrix(n, n, lambda i, j: reader(values.get((c["coords"][i], c["coords"][j]), "0")))


class Curvature(unittest.TestCase):
    def test_every_chart_is_conformally_flat_with_every_scalar_zero(self):
        # Koutras and McIntosh: no scalar invariant; Edgar and Ludwig: conformally flat pure radiation.
        for c in load("metrics")["coordinates"]:
            self.assertEqual(c["weyl_tensor"]["variants"]["llll"]["nonzero"], [], c["id"])
            self.assertEqual(c["ricci_scalar"], "R = 0", c["id"])
            self.assertEqual(c["kretschmann"], "K = 0", c["id"])

    def test_the_ricci_tensor_lies_along_the_rays(self):
        # Griffiths and Podolsky's Phi_22 = 2f/(ax + b) is R_uu/2(ax + b)^2: R_uu = 4(ax + b)f alone.
        km = chart("koutras_mcintosh")["ricci_tensor"]["variants"]["ll"]["nonzero"]
        self.assertEqual([(e["indices"], e["value"]) for e in km], [(["u", "u"], "4f\\left(a\\,x + b\\right)")])
        el = chart("edgar_ludwig")["ricci_tensor"]["variants"]["ll"]["nonzero"]
        self.assertEqual([(e["indices"], e["value"]) for e in el], [(["u", "u"], "4f\\,x")])


class Charts(unittest.TestCase):
    def test_wils_metric_is_koutras_and_mcintosh_at_a_one_b_zero(self):
        # Edgar and Ludwig's eq. (8) at g = h = 0 is eq. (3) at a = 1, b = 0 pulled back along w = v/x.
        km_reader, km = matrix("koutras_mcintosh")
        el_reader, el = matrix("edgar_ludwig")
        u, v, x, y = (el_reader.symbol[n] for n in ("u", "v", "x", "y"))
        F = sp.Rational(1, 2) + u / 3
        image = [u, v / x, x, y]
        J = sp.Matrix(4, 4, lambda i, j: sp.diff(image[i], (u, v, x, y)[j]))
        there = km.subs({km_reader.parameters["a"]: 1, km_reader.parameters["b"]: 0,
                         km_reader.parameters["f"]: F.subs(u, km_reader.symbol["u"])}).doit()
        there = there.subs(dict(zip((km_reader.symbol[n] for n in ("u", "w", "x", "y")), image)), simultaneous=True)
        here = el.subs({el_reader.parameters["f"]: F, el_reader.parameters["g"]: 0, el_reader.parameters["h"]: 0}).doit()
        pulled = J.T * there * J
        for at in ({u: sp.Rational(1, 3), v: sp.Rational(5, 4), x: sp.Rational(7, 5), y: sp.Rational(-2, 3)},
                   {u: -1, v: sp.Rational(-1, 2), x: sp.Rational(1, 4), y: 2}):
            for i in range(4):
                for j in range(4):
                    self.assertAlmostEqual(complex(sp.N((pulled[i, j] - here[i, j]).subs(at))), 0, places=20)

    def test_the_plane_wave_is_koutras_and_mcintosh_at_a_zero(self):
        # At a = 0 the metric is dx^2 + dy^2 - 2b du dw - 2fb(x^2 + y^2) du^2, a plane wave.
        reader, km = matrix("koutras_mcintosh")
        a, b, f = (reader.parameters[n] for n in ("a", "b", "f"))
        u, w, x, y = (reader.symbol[n] for n in ("u", "w", "x", "y"))
        plane = sp.Matrix([[-2 * f * b * (x ** 2 + y ** 2), -b, 0, 0], [-b, 0, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
        self.assertEqual(sp.expand(km.subs(a, 0) - plane), sp.zeros(4, 4))


class Drawings(unittest.TestCase):
    def test_every_chart_has_a_spacetime_diagram_of_the_allowed_shape(self):
        diagrams = load("diagrams")
        views = {(s, v["id"]): v for s, vs in diagrams["systems"].items() for v in vs}
        self.assertEqual(set(views), {("koutras_mcintosh", "front"), ("edgar_ludwig", "front")})
        for view in views.values():
            X0, X1, Y0, Y1 = view["box"]
            self.assertTrue(0.5 <= (X1 - X0) / (Y1 - Y0) <= 2.0)
        self.assertIn("kerr_schild", json.dumps(diagrams["projections"]))

    def test_no_conformal_and_no_embedding_diagram_by_design(self):
        # The free amplitude and the singular envelope leave no global completion, and no slice is
        # spacelike across the radiation; koutras_mcintosh.md gives both reasons.
        conformal = (ROOT / "_tools" / "derivations" / "conformal.py").read_text(encoding="utf-8")
        self.assertIn('NOT_DRAWN.add("koutras_mcintosh")', conformal)
        tree = ast.parse((ROOT / "_tools" / "derivations" / "embedding.py").read_text(encoding="utf-8"))
        not_drawn = next(ast.literal_eval(node.value) for node in tree.body if isinstance(node, ast.Assign)
                         and [getattr(t, "id", None) for t in node.targets] == ["NOT_DRAWN"])
        self.assertIn(ME, not_drawn)
        self.assertFalse((DATA / "embedding" / f"{ME}.json").exists())
        self.assertFalse((DATA / "curvature" / f"{ME}.json").exists())


if __name__ == "__main__":
    unittest.main()
