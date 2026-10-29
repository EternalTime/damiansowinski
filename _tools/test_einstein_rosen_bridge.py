"""Einstein and Rosen's bridges: the horizon every view marks.

The readers of the published files and MFS_SPACETIME are those of test_build_mfs_data.py."""
import math
import unittest

from test_build_mfs_data import checks, diagram_file


@checks("einstein_rosen_bridge")
class EinsteinRosenHorizons(unittest.TestCase):
    """Every view of the Einstein-Rosen bridge whose plotted region holds a horizon marks it. The
    bridge is the horizon of Schwarzschild's metric, and of Reissner and Nordström's for the
    charged bridge: r = r_s in the spherical chart and r_s/4 in the isotropic one, r_+ in the
    charged chart, and u = 0 in both of Einstein and Rosen's own. In the two spherical charts it
    is the edge of the chart, where g^rr vanishes, and the generator lost its marker there on
    2 October 2026, when zero_set began to blank everything outside a row's `where`."""

    # chart -> the kind of marker and the horizon's place on the horizontal axis, in the view's units.
    HORIZONS = {"bridge": ("throat", 0.0), "spherical": ("grr", 1.0), "isotropic": ("throat", 0.25),
                "charged_spherical": ("grr", (1 + math.sqrt(1 + 4 * 3 / 4)) / 2), "charged_bridge": ("throat", 0.0)}

    def setUp(self):
        self.views = [(system, view) for system, views in diagram_file("einstein_rosen_bridge")["systems"].items()
                      for view in views]

    def test_every_chart_is_known_to_the_test(self):
        self.assertEqual({system for system, _ in self.views}, set(self.HORIZONS))

    def test_every_view_that_holds_a_horizon_marks_it_from_bottom_to_top(self):
        marked = 0
        for system, view in self.views:
            kind, at = self.HORIZONS[system]
            X0, X1 = view["box"][:2]
            if not X0 < at < X1:
                continue
            where = f"einstein_rosen_bridge/{system}/{view['id']}"
            lines = [line for marker in view["markers"] if marker["kind"] == kind for line in marker["lines"]]
            on_it = [line for line in lines
                     if all(abs(X0 + u * (X1 - X0) - at) < 2e-3 * (X1 - X0) for u, _ in line)]
            self.assertTrue(on_it, f"{where} holds the horizon at {at:.3f} and draws no {kind} marker on it")
            heights = [v for line in on_it for _, v in line]
            self.assertLess(min(heights), 1e-3, where)
            self.assertGreater(max(heights), 1 - 1e-3, where)
            marked += 1
        self.assertEqual(marked, len(self.views), "a view of the bridge leaves its horizon out of the box")

    def test_the_charged_horizon_is_the_one_the_view_declares(self):
        view, = [view for system, view in self.views if system == "charged_spherical"]
        self.assertIn("$r_s = 1$", view["settings"])
        self.assertIn("$r_q = \\frac{\\sqrt{3}}{2}$", view["settings"])
        self.assertAlmostEqual(self.HORIZONS["charged_spherical"][1], 1.5)


if __name__ == "__main__":
    unittest.main()
