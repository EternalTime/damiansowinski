"""Roberts's collapse of a scalar field, from the numbers in its files alone.

The readers of the published files and MFS_SPACETIME are those of test_build_mfs_data.py."""
import json
import math
import unittest

import build_mfs_data as build
from test_build_mfs_data import ROBERTS_P, checks


@checks("roberts")
class RobertsCollapse(unittest.TestCase):
    """The published drawings of Roberts's collapse against its closed forms: where the apparent
    horizon and the singularity stand in each chart, and the shape of each embedded moment."""

    @classmethod
    def setUpClass(cls):
        root = build.METRICS_DIR.parent
        cls.diagrams = json.loads((root / "diagrams" / "roberts.json").read_text(encoding="utf-8"))
        cls.conformal = json.loads((root / "conformal" / "roberts.json").read_text(encoding="utf-8"))
        cls.embedding = json.loads((root / "embedding" / "roberts.json").read_text(encoding="utf-8"))

    @staticmethod
    def null(chart, p, a, b):
        """u and v, at l = 1, of the point (a, b) of a chart's plane."""
        if chart == "double_null":
            return a, b
        if chart == "advanced":
            return (1 + p) * a - 2 * b, a
        if chart == "areal":
            return a - 2 * math.sqrt(b * b + p * p * a * a / 4), a
        if chart == "diagonal":
            s = math.sqrt(1 + p)
            return s * (a - b), (a + b) / s
        return -2 * math.exp(-a), math.exp(-a) * math.expm1(2 * b)

    def marked(self, kind):
        """Every point of every marker of one kind, as (chart, case, u, v)."""
        for chart, views in self.diagrams["systems"].items():
            for view in views:
                X0, X1, Y0, Y1 = view["box"]
                (m00, m01), (m10, m11) = view["to_display"]
                det = m00 * m11 - m01 * m10
                for marker in view["markers"]:
                    if marker["kind"] != kind:
                        continue
                    for line in marker.get("lines", []):
                        for ux, uy in line:
                            X, Y = X0 + ux * (X1 - X0), Y0 + uy * (Y1 - Y0)
                            a, b = (m11 * X - m01 * Y) / det, (m00 * Y - m10 * X) / det
                            yield chart, view["id"], *self.null(chart, ROBERTS_P[view["id"]], a, b)

    def test_the_apparent_horizon_is_u_equal_to_one_minus_p_squared_times_v(self):
        found = set()
        for kind in ("apparent", "grr"):
            for chart, case, u, v in self.marked(kind):
                p = ROBERTS_P[case]
                self.assertGreater(p, 1, f"{chart}/{case}")
                self.assertAlmostEqual(u, (1 - p * p) * v, delta=5e-3 * (1 + abs(u) + abs(v)), msg=f"{chart}/{case}")
                found.add(chart)
        self.assertEqual(found, {"double_null", "advanced", "areal", "diagonal", "scaling"})

    def test_the_singularity_is_u_equal_to_one_minus_p_times_v_and_only_from_the_threshold_up(self):
        found = set()
        for chart, case, u, v in self.marked("singular"):
            p = ROBERTS_P[case]
            self.assertGreaterEqual(p, 1, f"{chart}/{case}")
            self.assertGreaterEqual(v, -1e-3, f"{chart}/{case}")
            self.assertAlmostEqual(u, (1 - p) * v, delta=5e-3 * (1 + abs(u) + abs(v)), msg=f"{chart}/{case}")
            found.add((chart, case))
        # At p = 1 the scaling chart's singularity is tau and x to infinity together, off every box.
        self.assertEqual(found, {(c, k) for c in self.diagrams["systems"] for k in ("threshold", "collapses")}
                         - {("scaling", "threshold")})

    def test_each_embedded_moment_has_the_areal_radius_and_the_height_of_its_closed_form(self):
        """R = sqrt(rho (rho - p t)) and Z = p t ln(sqrt(rho) + sqrt(rho - p t)) with the rim at Z = 0, in
        Minkowski space, from the first ray before t = 0 and from the last ray or the singularity after."""
        frames = 0
        for view in self.embedding["views"]:
            p = ROBERTS_P[view["id"]]
            for frame in view["movie"]["frames"]:
                t = float(frame["value"])
                (piece,) = frame["pieces"]
                self.assertEqual((piece["system"], piece["coordinate"]), ("diagonal", "\\rho"))
                points = [[float(x) for x in point] for point in piece["points"]]
                top = points[-1][0]
                first = -t if t <= 0 else t if p < 1 else p * t
                self.assertAlmostEqual(points[0][0], first, delta=3e-4 * max(abs(first), 1e-9) + 1e-12, msg=view["id"])
                for rho, R, Z in points:
                    self.assertAlmostEqual(R, math.sqrt(max(rho * (rho - p * t), 0.0)), delta=1e-9, msg=view["id"])
                    want = 0.0 if t == 0 else p * t * (math.log(math.sqrt(rho) + math.sqrt(max(rho - p * t, 0.0)))
                                                         - math.log(math.sqrt(top) + math.sqrt(top - p * t)))
                    self.assertAlmostEqual(Z, want, delta=1e-9, msg=f"{view['id']} at ct = {t}")
                frames += 1
        self.assertEqual(frames, 99)

    def test_the_conformal_singularity_is_null_on_the_threshold_and_level_at_p_equal_to_two(self):
        for view in self.conformal["views"]:
            case = view["id"].rsplit("_", 1)[1]
            lines = [layer["points"] for layer in view["layers"] if layer["class"] == "singular"]
            if case == "disperses":
                self.assertEqual(lines, [], view["id"])
                continue
            (line,) = lines
            for X, T in line:
                self.assertAlmostEqual(T, X if case == "threshold" else 0.0, delta=2e-4, msg=view["id"])


if __name__ == "__main__":
    unittest.main()
