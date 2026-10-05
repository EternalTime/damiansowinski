"""Kiselev's black hole, Schwarzschild's with a power of the radius added for matter around it:
.venv.noindex/bin/python -m unittest discover -s _tools

What its texts and drawings state, held on the published files. The first class needs nothing
but the files: the horizons of the example w = -2/3 at r_q = 8 r_s where every drawing puts them,
the cycloid the matter alone embeds as, and the null infinity and the straight centres of the
conformal diagrams. The second reads the published Einstein tensor through the checker's Reader
and holds it to the pressures the History and kiselev.md state, and the published metric to the
temperatures; it needs sympy and is skipped where it is absent.
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

INNER, OUTER = 4 - 2 * math.sqrt(2), 4 + 2 * math.sqrt(2)


def load(folder):
    return json.loads((DATA / folder / "kiselev.json").read_text(encoding="utf-8"))


class Drawings(unittest.TestCase):
    def test_the_spacetime_diagrams_mark_both_horizons_of_the_example(self):
        """g^rr = 1 - r_s/r - r/r_q vanishes at r = (4 -+ 2 sqrt 2) r_s for r_q = 8 r_s, and each
        view of t, v or u against r marks those two lines and no other."""
        systems = load("diagrams")["systems"]
        for system in ("static", "linear", "eddington_finkelstein_ingoing", "eddington_finkelstein_outgoing"):
            for view in systems[system]:
                X0, X1 = view["box"][:2]
                marked = sorted(X0 + line[0][0] * (X1 - X0)
                                for m in view["markers"] if m["kind"] == "grr" for line in m["lines"])
                self.assertEqual(len(marked), 2, f"{system}/{view['id']}")
                for got, want in zip(marked, (INNER, OUTER)):
                    self.assertAlmostEqual(got, want, delta=2e-3, msg=f"{system}/{view['id']}")

    def test_the_embedding_runs_from_the_throat_to_the_widest_circle(self):
        view = next(v for v in load("embedding")["views"] if v["id"] == "black_hole")
        for piece in view["surfaces"][0]["pieces"]:
            radii = [p[1] for p in piece["points"]]
            self.assertAlmostEqual(min(radii), INNER, places=6)
            self.assertAlmostEqual(max(radii), OUTER, places=6)

    def test_the_matter_alone_embeds_as_a_cycloid(self):
        """dz/dr = sqrt(r/(r_q - r)) integrates to z = r_q(asin sqrt(r/r_q) - sqrt((r/r_q)(1 - r/r_q))),
        which with r = r_q(1 - cos psi)/2 is z = r_q(psi - sin psi)/2. Each sheet leaves its centre at
        z = -+pi r_q/2 and the two meet on the horizon r = r_q at z = 0."""
        view = next(v for v in load("embedding")["views"] if v["id"] == "free")
        pieces = view["surfaces"][0]["pieces"]
        self.assertEqual(len(pieces), 2)
        checked = 0
        for piece in pieces:
            points = piece["points"]
            side = math.copysign(1.0, points[0][2])
            self.assertAlmostEqual(abs(points[0][2]), math.pi / 2, places=6)
            for r, rho, z in points:
                self.assertAlmostEqual(rho, r, places=6)
                psi = math.acos(max(-1.0, min(1.0, 1 - 2 * r)))
                self.assertAlmostEqual(z, side * (math.pi / 2 - (psi - math.sin(psi)) / 2), delta=2e-6)
                checked += 1
        self.assertGreater(checked, 100)

    def test_infinity_is_null_and_the_centres_of_the_matter_alone_are_straight(self):
        """Beyond the outer horizon r* falls as -r_q ln r, so r = infinity is drawn on straight edges
        at 45 degrees that meet at (pi, +-pi), and without the black hole each centre r = 0 is a
        vertical line, since one surface gravity serves both null coordinates."""
        views = {v["id"]: v for v in load("conformal")["views"]}
        self.assertEqual(sorted(views), ["conformally_flat", "hyperbolic", "ingoing", "linear", "outgoing", "static"])
        for vid, view in views.items():
            scri = [layer["points"] for layer in view["layers"] if layer["class"] == "scri"]
            self.assertEqual(len(scri), 4, vid)
            tops = set()
            for line in scri:
                self.assertEqual(len(line), 2, vid)
                (xa, ta), (xb, tb) = line
                self.assertAlmostEqual(abs(xa - xb), abs(ta - tb), places=3, msg=vid)
                tops |= {(round(x, 3), round(t, 3)) for x, t in line}
            self.assertIn((round(math.pi, 3), round(math.pi, 3)), tops, vid)
            self.assertIn((round(math.pi, 3), -round(math.pi, 3)), tops, vid)
        for vid in ("hyperbolic", "conformally_flat"):
            centres = [layer["points"] for layer in views[vid]["layers"] if layer["class"] == "singular"]
            self.assertEqual(len(centres), 2, vid)
            for line, where in zip(sorted(centres), (math.pi / 2, 3 * math.pi / 2)):
                for x, _ in line:
                    self.assertAlmostEqual(x, where, places=3, msg=vid)


@unittest.skipUnless(HAS_SYMPY, "sympy is not installed")
class Matter(unittest.TestCase):
    def setUp(self):
        import sympy
        import verify_metrics
        self.sp, self.vm = sympy, verify_metrics
        metric = json.loads((DATA / "metrics" / "kiselev.json").read_text(encoding="utf-8"))
        self.charts = {c["id"]: c for c in metric["coordinates"]}

    def read(self, system, field, variant=None):
        entry = self.charts[system]
        reader = self.vm.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]], ())
        block = entry[field] if variant is None else entry[field]["variants"][variant]["nonzero"]
        return reader, {tuple(e["indices"]): reader(e["value"]) for e in block}

    def test_the_matter_has_two_pressures_whose_average_is_w_times_the_density(self):
        """With G = 8 pi T the published mixed Einstein tensor gives the density -G^t_t/8 pi, the radial
        pressure G^r_r/8 pi and the tangential pressure G^theta_theta/8 pi: p_r = -rho,
        p_t = (3w + 1) rho/2, and (p_r + 2 p_t)/3 = w rho, Visser's (3) to (5)."""
        sp = self.sp
        reader, G = self.read("static", "einstein_tensor", "ul")
        w = reader.parameters["w"]
        rho, p_r, p_t = -G[("t", "t")], G[("r", "r")], G[("\\theta", "\\theta")]
        self.assertEqual(sp.simplify(p_r + rho), 0)
        self.assertEqual(sp.simplify(p_t - (3 * w + 1) * rho / 2), 0)
        self.assertEqual(sp.simplify((p_r + 2 * p_t) / 3 - w * rho), 0)
        self.assertEqual(sp.simplify(G[("\\phi", "\\phi")] - p_t), 0)
        # The two pressures agree only for the cosmological constant.
        self.assertEqual(sp.solve(sp.simplify((p_t - p_r) / rho), w), [-1])
        # The density is positive for w < 0, at w = -2/3 and r = 3 r_s with r_q = 8 r_s.
        at = {w: sp.Rational(-2, 3), reader.symbol["r"]: 3, reader.parameters["r_s"]: 1, reader.parameters["r_q"]: 8}
        self.assertGreater(float(rho.subs(at)), 0)

    def test_the_matter_cools_the_black_hole(self):
        """The surface gravity of the black hole horizon of the example, f'(r_-)/2 with
        f = 1 - r_s/r - r/r_q, is below Schwarzschild's 1/2r_s for every r_q > 4 r_s, and it vanishes
        at r_q = 4 r_s, where the two horizons meet at 2 r_s."""
        sp = self.sp
        reader, g = self.read("linear", "inverse_metric_components")
        r, rs, rq = reader.symbol["r"], reader.parameters["r_s"], reader.parameters["r_q"]
        f = g[("r", "r")]
        inner = (rq - sp.sqrt(rq ** 2 - 4 * rq * rs)) / 2
        self.assertEqual(sp.simplify(f.subs(r, inner)), 0)
        kappa = sp.diff(f, r).subs(r, inner) / 2
        for ratio in (5, 8, 100):
            self.assertLess(float(kappa.subs({rs: 1, rq: ratio})), 0.5)
            self.assertGreater(float(kappa.subs({rs: 1, rq: ratio})), 0)
        self.assertEqual(sp.simplify(kappa.subs(rq, 4 * rs)), 0)
        self.assertEqual(sp.simplify(inner.subs(rq, 4 * rs) - 2 * rs), 0)

    def test_the_matter_alone_is_conformally_flat_with_a_singular_centre(self):
        for system in ("hyperbolic", "conformally_flat"):
            self.assertEqual(self.charts[system]["weyl_tensor"]["variants"]["ulll"]["nonzero"], [], system)
        reader, _ = self.read("conformally_flat", "metric_components")
        K = reader(self.charts["conformally_flat"]["kretschmann"].partition("=")[2])
        tau, rho = reader.symbol["\\tau"], reader.symbol["\\rho"]
        at = {reader.parameters["r_q"]: 1, tau: 1}
        self.assertGreater(float(K.subs(at).subs(rho, 1e-6)), 1e11)


if __name__ == "__main__":
    unittest.main()
