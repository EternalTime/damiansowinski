"""The markers of every view whose row carries a `where`, held to what the generator draws today:
python3 -m unittest discover -s _tools

A row's `where` is positive on the spacetime, and nothing is marked outside it. Since the website's
`fa93957` of 2 October 2026 `Plot.zero_set` blanks the whole grid there, so a curve lying on the very
edge of `where`, as the horizon r = r_s of the Einstein-Rosen bridge's spherical chart does, has no
change of sign left to find, and that view lost its horizon marker the next time it was redrawn.
Redrawing a view traces its rays, which takes minutes, but its markers take a second, so every such
row's markers are computed here and held to the published file: a change to the generator that
drops or adds a marker on one of them fails here and names the view. The tests need sympy, numpy,
scipy and contourpy and are skipped where any is absent.
"""
import importlib.util
import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DIAGRAMS = ROOT / "MFS" / "assets" / "data" / "diagrams"
READY = all(importlib.util.find_spec(name) is not None for name in ("sympy", "numpy", "scipy", "contourpy"))
sys.path.insert(0, str(Path(__file__).resolve().parent / "derivations"))


def published(spec):
    data = json.loads((DIAGRAMS / f"{spec.metric}.json").read_text(encoding="utf-8"))
    return next(view for view in data["systems"][spec.system] if view["id"] == spec.view)


@unittest.skipUnless(READY, "sympy, numpy, scipy and contourpy are not installed")
class MarkersAtTheEdge(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import null_rays
        cls.nr = null_rays
        cls.rows = [spec for spec in null_rays.DIAGRAMS if spec.where]

    def test_there_are_rows_to_hold(self):
        self.assertGreaterEqual(len(self.rows), 25)

    def test_every_view_with_a_where_draws_the_markers_its_file_carries(self):
        for spec in self.rows:
            place = f"{spec.metric}/{spec.system}/{spec.view}"
            with self.subTest(place):
                drawn = self.nr.Plot(self.nr.Chart(spec)).markers()
                self.assertEqual(json.loads(json.dumps(drawn)), published(spec)["markers"],
                                 f"{place}: the generator's markers are no longer the published ones")

    def test_a_horizon_on_the_edge_of_where_is_marked_only_where_the_row_says_so(self):
        declared = {(spec.metric, spec.system, spec.view) for spec in self.rows if spec.edge_horizon}
        self.assertEqual(declared, {("einstein_rosen_bridge", "spherical", "radial"),
                                    ("einstein_rosen_bridge", "charged_spherical", "radial"),
                                    ("btz_multi_holes_wormholes", "exterior", "radial")})
        for spec in self.rows:
            kinds = [marker["kind"] for marker in published(spec)["markers"]]
            if spec.edge_horizon:
                self.assertIn("grr", kinds, f"{spec.metric}/{spec.system}/{spec.view}")
                edge = self.nr.Plot(self.nr.Chart(spec)).zero_set_on_the_edge("girr")
                self.assertEqual(len(edge), 1)

    def test_a_row_that_declares_a_horizon_its_edge_does_not_have_is_refused(self):
        row = next(spec for spec in self.rows if (spec.metric, spec.system) == ("hiscock", "flat"))
        from dataclasses import replace
        with self.assertRaises(SystemExit) as refused:
            self.nr.Plot(self.nr.Chart(replace(row, edge_horizon=True))).markers()
        self.assertIn("does not vanish on the edge", str(refused.exception))


if __name__ == "__main__":
    unittest.main()
