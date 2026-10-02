"""Cremmer and Scherk's spontaneous compactification, Minkowski space times a sphere in six dimensions:
python3 -m unittest discover -s _tools

What its texts state, held on the published files: the Einstein tensor is -1/a^2 on the four flat
dimensions and nothing on the sphere, which with the cosmological constant 1/2a^2 the parameter
states is the stress of a monopole's field of one strength over the sphere; Horvath, Palla, Cremmer
and Scherk's two relations, V_0 = e^2/(k (8 pi G)^2) and R_0^2 = 8 pi G k/e^2, are that
cosmological constant and that field; the flat plane times a sphere of Plebanski and Hacyan's page
has the same cosmological constant for its radius, as the two relations say; and the drawings hold
what their captions say. The tests of the published mathematics need sympy and are skipped where it
is absent; the tests of the drawings read the files alone.
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


def published(metric_id="cremmer_scherk"):
    return json.loads((DATA / "metrics" / f"{metric_id}.json").read_text(encoding="utf-8"))


def drawing(kind, metric_id="cremmer_scherk"):
    return json.loads((DATA / kind / f"{metric_id}.json").read_text(encoding="utf-8"))


@unittest.skipUnless(HAS_SYMPY, "sympy is not installed")
class Published(unittest.TestCase):
    def setUp(self):
        import sympy
        import verify_metrics
        self.sp, self.vm = sympy, verify_metrics

    def read(self, system, metric_id="cremmer_scherk"):
        """The reader of a published chart and a function from a field's name to the matrix of its
        published components, with c = 1."""
        entry = next(c for c in published(metric_id)["coordinates"] if c["id"] == system)
        reader = self.vm.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]], ())
        names = entry["coords"]

        def matrix(field, variant=None):
            block = entry[field] if variant is None else entry[field]["variants"][variant]["nonzero"]
            M = self.sp.zeros(len(names), len(names))
            for e in block:
                i, j = (names.index(k) for k in e["indices"])
                M[i, j] = reader(e["value"]).subs(reader.c, 1)
            return M
        return reader, matrix, entry

    def test_the_einstein_tensor_is_the_monopole_and_the_cosmological_constant(self):
        # G^mu_nu + Lambda delta = 8 pi G T^mu_nu, with Lambda = 1/2a^2 and the stress of a field of
        # one strength B over the sphere, (B^2/2) diag(-1, -1, -1, -1, 1, 1), at 8 pi G B^2 = 1/a^2.
        sp = self.sp
        reader, matrix, entry = self.read("cartesian")
        a = reader.parameters["a"]
        mixed = matrix("einstein_tensor", "ul")
        self.assertEqual(mixed, sp.diag(*([-1 / a ** 2] * 4 + [0, 0])))
        stress = sp.diag(*([sp.Rational(-1, 2)] * 4 + [sp.Rational(1, 2)] * 2)) / a ** 2
        self.assertEqual(sp.simplify(mixed + sp.eye(6) / (2 * a ** 2) - stress), sp.zeros(6, 6))
        self.assertIn("\\Lambda = 1/2a^2", entry["parameters"][0]["description"])

    def test_the_curvature_is_the_sphere_alone(self):
        sp = self.sp
        reader, matrix, entry = self.read("cartesian")
        a = reader.parameters["a"]
        theta = reader.symbol["\\theta"]
        ricci = matrix("ricci_tensor", "ll")
        self.assertEqual(sp.simplify(ricci - sp.diag(0, 0, 0, 0, 1, sp.sin(theta) ** 2)), sp.zeros(6, 6))
        self.assertEqual(sp.simplify(reader(entry["ricci_scalar"].partition("=")[2]) - 2 / a ** 2), 0)
        self.assertEqual(sp.simplify(reader(entry["kretschmann"].partition("=")[2]) - 4 / a ** 4), 0)
        # A product with a factor that is not conformally flat on its own account has a Weyl tensor.
        self.assertTrue(entry["weyl_tensor"]["variants"]["llll"]["nonzero"])

    def test_the_four_authors_relations_are_the_radius_and_the_cosmological_constant(self):
        # Their (4): V_0 = e^2/(k (8 pi G)^2) and R_0^2 = 8 pi G k/e^2. The field of their (3) has the
        # strength B = sqrt(k)/(e R_0^2), so 8 pi G B^2 = 1/R_0^2, and the term V_0/2 of their action is
        # Lambda = 8 pi G V_0/2 = 1/2R_0^2.
        sp = self.sp
        G, e, k = sp.symbols("G e k", positive=True)
        V0 = e ** 2 / (k * (8 * sp.pi * G) ** 2)
        R2 = 8 * sp.pi * G * k / e ** 2
        B2 = k / (e ** 2 * R2 ** 2)
        self.assertEqual(sp.simplify(8 * sp.pi * G * B2 - 1 / R2), 0)
        self.assertEqual(sp.simplify(8 * sp.pi * G * V0 / 2 - 1 / (2 * R2)), 0)
        description = published()["coordinates"][0]["parameters"][0]["description"]
        self.assertIn("a^2 = 8\\pi Gk/e^2", description)

    def test_the_flat_plane_times_a_sphere_has_the_same_cosmological_constant_for_its_radius(self):
        # The relation with Plebanski and Hacyan's page: in four dimensions G^mu_nu = -1/b^2 on the flat
        # plane and nothing on the sphere, the same balance, with Lambda = 1/2b^2.
        sp = self.sp
        reader, matrix, entry = self.read("sphere", "plebanski_hacyan")
        b = reader.parameters["b"]
        self.assertEqual(matrix("einstein_tensor", "ul"), sp.diag(-1 / b ** 2, -1 / b ** 2, 0, 0))
        self.assertIn("\\Lambda = 1/2b^2", entry["parameters"][0]["description"])
        mine = next(r for r in published()["related"] if r["id"] == "plebanski_hacyan")
        self.assertIn("\\Lambda = 1/2b^2", mine["text"])


class Drawings(unittest.TestCase):
    def test_the_chart_has_two_spacetime_diagrams_and_a_conformal_diagram(self):
        diagrams = drawing("diagrams")
        self.assertEqual(list(diagrams["systems"]), ["cartesian"])
        self.assertEqual([v["id"] for v in diagrams["systems"]["cartesian"]], ["tx", "circle"])
        self.assertEqual([v["id"] for v in drawing("conformal")["views"]], ["cartesian"])

    def test_a_ray_goes_round_the_equator_in_the_time_the_caption_states(self):
        # On the plane of the time and phi at a = 1 the rays are at 45 degrees and phi = 0 and 2 pi are
        # one line, so a ray is back after c t = 2 pi a: the box is one turn wide and one turn high.
        view = next(v for v in drawing("diagrams")["systems"]["cartesian"] if v["id"] == "circle")
        x0, x1, y0, y1 = view["box"]
        self.assertAlmostEqual(x1 - x0, 2 * math.pi, places=6)
        self.assertAlmostEqual(y1 - y0, 2 * math.pi, places=6)
        self.assertEqual(view["hatch"], [])
        self.assertIn("$2\\pi a/c$", " ".join(view["caption"]))
        self.assertTrue(view["cones"])
        for cone in view["cones"]:
            # Each cone's two edges are the null directions d(ct) = +-a dphi.
            for edge in (cone["a"], cone["b"]):
                self.assertAlmostEqual(abs(edge[0]), abs(edge[1]), places=4)

    def test_the_embedding_is_a_cylinder_and_a_sphere_of_the_one_radius(self):
        views = {v["id"]: v for v in drawing("embedding")["views"]}
        self.assertEqual(sorted(views), ["equator", "sphere"])
        tube = views["equator"]["surfaces"][0]["pieces"][0]
        self.assertEqual((tube["system"], tube["coordinate"]), ("cartesian", "x"))
        for x, rho, height in tube["points"]:
            self.assertAlmostEqual(rho, 1.0, places=6)
            self.assertAlmostEqual(height, x, places=6)
        ball = views["sphere"]["surfaces"][0]["pieces"][0]
        for theta, rho, height in ball["points"]:
            self.assertAlmostEqual(rho ** 2 + (height - 1) ** 2, 1.0, places=5)
            self.assertAlmostEqual(rho, math.sin(theta), places=5)

    def test_the_conformal_diagram_is_the_whole_diamond(self):
        view = drawing("conformal")["views"][0]
        cover = next(layer for layer in view["layers"] if layer["kind"] == "fill" and layer["class"] == "cover")
        pi = round(math.pi, 4)
        self.assertEqual([[round(x, 4), round(t, 4)] for x, t in cover["points"]],
                         [[pi, 0], [0, pi], [-pi, 0], [0, -pi]])
        # The cylinder's moment t = 0 from x = -a to a: p, q = arctan(-+1), so X = q - p = +-pi/2.
        line = view["slices"][0]["lines"][0]
        self.assertAlmostEqual(line[0][0], -math.pi / 2, places=3)
        self.assertAlmostEqual(line[-1][0], math.pi / 2, places=3)
        self.assertAlmostEqual(line[0][1], 0.0, places=6)

    def test_the_sphere_is_not_marked_on_the_plane_it_does_not_meet(self):
        views = {v["id"]: v for v in drawing("diagrams")["systems"]["cartesian"]}
        self.assertEqual([m["view"] for m in views["tx"]["slices"]], ["equator", "sphere"])
        self.assertEqual([m["view"] for m in views["circle"]["slices"]], ["equator"])


if __name__ == "__main__":
    unittest.main()
