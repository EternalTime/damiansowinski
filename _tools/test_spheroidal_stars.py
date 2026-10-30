"""Family identities and coverage beyond the generator's matter and junction audit."""
import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "MFS/assets/data"
sys.path.insert(0, str(ROOT / "_tools/derivations"))

import sympy as sp
import spheroidal_stars as stars
import verify_metrics as vm


class StellarFamilies(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.metric = json.loads((DATA / "metrics/spheroidal_stars.json").read_text())
        cls.charts = {c["id"]: c for c in cls.metric["coordinates"]}

    def reader(self, system):
        chart = self.charts[system]
        return vm.Reader(chart["coords"], [p["symbol"] for p in chart["parameters"]], ())

    def test_no_family_is_skipped_for_an_incomplete_dimension_declaration(self):
        for system in stars.SYSTEMS:
            declared = vm.DIMENSIONS[(stars.METRIC, system)]
            vm.Dimensions(self.reader(system), declared)

    def test_complete_models_and_the_credited_duplicates_share_one_chart(self):
        self.assertEqual(tuple(self.charts), stars.SYSTEMS)
        vt = self.charts["vaidya_tikekar"]
        self.assertIn("Durgapal-Bannerji", vt["name"])
        self.assertIn("Buchdahl", self.metric["history"])
        self.assertIn("Heintzmann", self.charts["durgapal_iii"]["name"])
        self.assertNotIn("durgapal_i", self.charts)
        self.assertNotIn("durgapal_ii", self.charts)

    def test_the_vt_radius_and_both_k_conventions_are_the_same_geometry(self):
        rd = self.reader("vaidya_tikekar")
        r, C = rd.symbol["r"], rd.parameters["C"]
        Z = rd.parameters["Z"]
        L = sp.Rational(3, 2)
        for value in (sp.Rational(1, 7), sp.Rational(1, 2), sp.Rational(3, 4)):
            grr = (1/Z).subs({r: value, C: 2/L**2})
            plus = (1+2*value**2/L**2)/(1-value**2/L**2)
            minus = (1-(-2)*value**2/L**2)/(1-value**2/L**2)
            self.assertEqual(grr, plus)
            self.assertEqual(grr, minus)

    def test_both_lapse_constants_are_independent(self):
        # The two basis functions span the entire second order isotropy solution space.
        # Their nonzero Wronskians also guard against retaining only the radial ansatz.
        for system in ("vaidya_tikekar", "finch_skea"):
            rd = self.reader(system)
            y, r = rd.parameters["y"], rd.symbol["r"]
            A, B = rd.parameters["A"], rd.parameters["B"]
            f, g = sp.diff(y, A), sp.diff(y, B)
            w = f*sp.diff(g, r)-g*sp.diff(f, r)
            for radius in (sp.Rational(1, 5), sp.Rational(4, 5)):
                value = w.subs({r: radius, rd.parameters["C"]: sp.Rational(1, 3)})
                self.assertGreater(abs(float(sp.N(value, 30))), 1e-3)

    def test_the_finite_buchdahl_star_and_the_1964_infinite_sphere_are_distinct(self):
        other = json.loads((DATA / "metrics/buchdahl_sphere.json").read_text())
        self.assertIn(r"r \in [0, r_b]".replace(" ", ""),
                      [d.replace(" ", "") for d in self.charts["vaidya_tikekar"]["domains"]])
        self.assertTrue(any(r"\infty" in d and d.startswith("r ") for d in other["coordinates"][0]["domains"]))

    def test_all_families_have_each_kind_of_diagram(self):
        for folder in ("diagrams", "conformal", "embedding", "curvature"):
            data = json.loads((DATA / folder / "spheroidal_stars.json").read_text())
            if folder == "diagrams":
                ids = set(data["systems"])
            else:
                ids = {v.get("system", v["id"]) for v in data["views"]}
            self.assertEqual(ids, set(stars.SYSTEMS), folder)


if __name__ == "__main__":
    unittest.main()
