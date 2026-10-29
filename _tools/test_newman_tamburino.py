"""The Newman-Tamburino metrics, held on the published file with the checker's own reader:
.venv.noindex/bin/python -m unittest _tools/test_newman_tamburino.py

Both charts are vacua whose rays d_r are null geodesics, affinely parametrised by r; Steele's
vectors are symmetries, the Killing vector -u d_u + r d_r + 2x d_x + 2y d_y and the homothety
r d_r + x d_x + y d_y of the spherical class at A_0 = 0, and d_u, d_y and 2r d_r - x d_x + 2y d_y of
the cylindrical limit; and the cylindrical limit's Kretschmann scalar changes sign on
r x^2 = e^(3/2) ell. Needs sympy, and is skipped where it is absent.
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
    import verify_metrics
    entry = next(c for c in json.loads((DATA / "metrics" / "newman_tamburino.json").read_text())["coordinates"]
                 if c["id"] == system)
    reader = verify_metrics.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]], ())
    if reader.held:
        def read(text):
            return reader(text).subs(reader.held).doit()
    else:
        read = reader
    return entry, reader, read


@unittest.skipUnless(HAS_SYMPY, "sympy is not installed")
class TheCharts(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import sympy
        cls.sp = sympy
        cls.charts = {}
        for system in ("spherical", "cylindrical"):
            entry, reader, read = chart(system)
            X = [reader.symbol[c] for c in entry["coords"]]
            g = sympy.zeros(4, 4)
            for e in entry["metric_components"]:
                i, j = (entry["coords"].index(c) for c in e["indices"])
                g[i, j] = read(e["value"])
            cls.charts[system] = (entry, reader, read, X, g)

    def lie(self, system, xi, at):
        """The Lie derivative of the metric along xi, at a point, to thirty digits."""
        sp = self.sp
        _, _, _, X, g = self.charts[system]
        out = sp.zeros(4, 4)
        for a in range(4):
            for b in range(4):
                out[a, b] = (sum(xi[c] * sp.diff(g[a, b], X[c]) for c in range(4))
                             + sum(g[c, b] * sp.diff(xi[c], X[a]) + g[a, c] * sp.diff(xi[c], X[b]) for c in range(4)))
        return out.subs(at).evalf(30), g.subs(at).evalf(30)

    def point(self, system):
        sp = self.sp
        entry, reader, _, X, _ = self.charts[system]
        at = dict(zip(X, (sp.Rational(1, 3), sp.Rational(9, 2), sp.Rational(3, 5), sp.Rational(1, 4))))
        at[reader.c] = 1
        if system == "spherical":
            at.update({reader.parameters["A_0"]: 0, reader.parameters["b"]: sp.Rational(1, 2)})
        else:
            at[reader.parameters["ell"]] = sp.Rational(7, 5)
        return at

    def test_both_charts_publish_a_vanishing_ricci_tensor_and_scalar(self):
        for system in ("spherical", "cylindrical"):
            entry = self.charts[system][0]
            self.assertEqual(entry["ricci_scalar"], "R = 0")
            self.assertEqual(json.loads(entry["ricci_tensor"])["variants"]["ll"]["nonzero"]
                             if isinstance(entry["ricci_tensor"], str) else entry["ricci_tensor"]["variants"]["ll"]["nonzero"],
                             [])

    def test_the_rays_are_null_and_affinely_parametrised_by_r(self):
        # g_rr = 0, and Gamma^mu_rr = 0, so d_r is a null geodesic with r affine.
        for system in ("spherical", "cylindrical"):
            entry = self.charts[system][0]
            self.assertNotIn(["r", "r"], [e["indices"] for e in entry["metric_components"]])
            christoffel = entry["christoffel"]
            christoffel = json.loads(christoffel) if isinstance(christoffel, str) else christoffel
            ull = christoffel["variants"]["ull"]["nonzero"]
            self.assertFalse([e for e in ull if e["indices"][1:] == ["r", "r"]])

    def test_steeles_vectors_are_symmetries(self):
        sp = self.sp
        for system, vectors in (
                ("spherical", (("killing", lambda u, r, x, y: [-u, r, 2 * x, 2 * y], 0),
                               ("homothety", lambda u, r, x, y: [0, r, x, y], 0.5))),
                ("cylindrical", (("d_u", lambda u, r, x, y: [1, 0, 0, 0], 0),
                                 ("d_y", lambda u, r, x, y: [0, 0, 0, 1], 0),
                                 ("homothety", lambda u, r, x, y: [0, 2 * r, -x, 2 * y], 1)))):
            X = self.charts[system][3]
            for name, field, scale in vectors:
                with self.subTest(system=system, vector=name):
                    lie, g = self.lie(system, field(*X), self.point(system))
                    worst = max(abs(lie[a, b] - 2 * scale * g[a, b]) for a in range(4) for b in range(4))
                    self.assertLess(worst, sp.Float(10) ** -25)

    def test_the_cylindrical_kretschmann_scalar_changes_sign_on_r_x2_equal_e32_ell(self):
        sp = self.sp
        entry, reader, read, X, _ = self.charts["cylindrical"]
        K = read(entry["kretschmann"].partition("=")[2])
        u, r, x, y = X
        ell = reader.parameters["ell"]
        for k in (sp.Rational(1, 2), 1, 3):
            edge = sp.exp(sp.Rational(3, 2)) / k ** 2
            at = {x: k, ell: 1, reader.c: 1}
            self.assertLess(abs(sp.N(K.subs({**at, r: edge}), 30)), 1e-25)
            self.assertGreater(sp.N(K.subs({**at, r: edge * sp.Rational(9, 10)})), 0)
            self.assertLess(sp.N(K.subs({**at, r: edge * sp.Rational(11, 10)})), 0)
        self.assertAlmostEqual(math.exp(1.5), 4.4817, places=4)


if __name__ == "__main__":
    unittest.main()
