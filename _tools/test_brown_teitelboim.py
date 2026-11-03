"""The charged membrane's junction and the drawings made from its two bulk regions."""
import json
import math
import unittest
from pathlib import Path

from test_build_mfs_data import checks, conformal_file, diagram_file, embedding_file


@checks("brown_teitelboim")
class BrownTeitelboimBubble(unittest.TestCase):
    """The depicted bubble has de Sitter radii 2 and 1 and minimum wall radius 1/2."""

    def test_the_signed_junction_on_the_expanding_and_contracting_wall(self):
        beta = math.sqrt(15) / 2 - math.sqrt(3)
        for tau in (-0.5, 0.0, 0.3, 1.0):
            radius = 0.5 * math.cosh(2 * tau)
            velocity = math.sinh(2 * tau)
            inner = math.sqrt(1 + velocity**2 - radius**2 / 4)
            outer = math.sqrt(1 + velocity**2 - radius**2)
            self.assertAlmostEqual(inner - outer, beta * radius, places=12)

    def test_the_radius_view_shows_the_same_wall_in_proper_time(self):
        view = diagram_file("brown_teitelboim")["systems"]["radius_inside"][0]
        (wall,) = [m for m in view["markers"] if m["kind"] == "shell"]
        x0, x1, y0, y1 = view["box"]
        for line in wall["lines"]:
            for x, y in line:
                radius, tau = x0 + x * (x1-x0), y0 + y * (y1-y0)
                self.assertAlmostEqual(radius, math.cosh(2*tau)/2, delta=0.00015)
        self.assertTrue(view["hatch"])

    def test_the_embedding_pieces_meet_at_the_membrane(self):
        view = embedding_file("brown_teitelboim")["views"][0]
        inner, outer = view["surfaces"][0]["pieces"]
        self.assertEqual(inner["points"][-1], outer["points"][0])
        self.assertAlmostEqual(inner["points"][-1][1], 0.5, places=7)
        for piece, ell in ((inner, 2), (outer, 1)):
            r0, _, z0 = piece["points"][0]
            for r, rho, z in piece["points"][::max(1, len(piece["points"])//10)]:
                self.assertAlmostEqual(rho, r, delta=1e-7)
                expected = z0 + math.sqrt(ell**2-r0**2) - math.sqrt(ell**2-r**2)
                self.assertAlmostEqual(z, expected, delta=1e-7)

    def test_the_conformal_wall_joins_the_two_sliced_regions(self):
        view = conformal_file("brown_teitelboim")["views"][0]
        (wall,) = [layer["points"] for layer in view["layers"] if layer["class"] == "surface"]
        self.assertEqual(len(wall), 2)
        self.assertEqual(wall[0][0], 0)
        self.assertEqual(wall[1][0], 0)
        self.assertAlmostEqual(wall[1][1], math.pi, delta=1e-4)
        self.assertAlmostEqual(wall[0][1], -math.pi, delta=1e-4)

    def test_bulk_curvature_takes_the_de_sitter_values(self):
        path = Path(__file__).resolve().parents[1] / "MFS/assets/data/curvature/brown_teitelboim.json"
        view = json.loads(path.read_text())["views"][0]
        for curve in view["curves"]:
            lam = 0.75 if curve["system"] == "static_inside" else 3
            expected = 8 * lam**2 / 3 if curve["scalar"] == "kretschmann" else 4 * lam
            for _, value in curve["points"]:
                self.assertAlmostEqual(value, expected, places=8)


if __name__ == "__main__":
    unittest.main()
