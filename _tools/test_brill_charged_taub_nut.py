"""Brill's charged Taub-NUT, held on the published files:
python3 -m unittest discover -s _tools

Every chart writes Sigma = r^2 + l^2 and Delta = r^2 - 2mr - l^2 + r_q^2. The published first
chart is the published Taub-NUT metric at r_q = 0 and Reissner and Nordstrom's at l = 0, its
Kretschmann scalar is the square of the Ricci tensor plus the square of the Weyl tensor that
Clement, Gal'tsov and Guenouche give and stays finite at r = 0, and the way print_charts.py
regroups a sum changes no value; those need sympy and are skipped where it is absent. The
drawings are held from their numbers alone: the tortoise coordinate the rays are checked against
to dr_*/dr = Sigma/Delta, the marked horizons to the roots of Delta, the wormhole's views to
marking none, and the embedded surfaces to the radius sqrt(r^2 + l^2), to the throat's radius l
and to the wormhole being the same on both sides.
"""
import importlib.util
import json
import math
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "MFS" / "assets" / "data"
DERIVATIONS = Path(__file__).resolve().parent / "derivations"
HAS_SYMPY = importlib.util.find_spec("sympy") is not None and importlib.util.find_spec("numpy") is not None
ID = "brill_charged_taub_nut"


def read(folder, metric_id=ID):
    return json.loads((DATA / folder / f"{metric_id}.json").read_text(encoding="utf-8"))


@unittest.skipUnless(HAS_SYMPY, "sympy and numpy are not installed")
class TheCharts(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        sys.path.insert(0, str(DERIVATIONS))
        import print_charts
        import slices
        import sympy
        import verify_metrics
        cls.pc, cls.slices, cls.sp, cls.vm = print_charts, slices, sympy, verify_metrics

    @classmethod
    def tearDownClass(cls):
        sys.path.remove(str(DERIVATIONS))

    def published(self, metric_id, chart_id):
        """(reader, the metric as a matrix) of a published chart, from its metric components."""
        entry = next(c for c in read("metrics", metric_id)["coordinates"] if c["id"] == chart_id)
        reader = self.vm.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]], ())
        values = {tuple(e["indices"]): e["value"] for e in entry["metric_components"]}
        g = self.sp.Matrix(4, 4, lambda i, j: reader(values.get((entry["coords"][i], entry["coords"][j]), "0")))
        return entry, reader, g

    def same(self, a, b, readers, at):
        """Two published metrics agree at rational points, with the second's symbols named by the first's."""
        sp = self.sp
        first, second = readers
        names = dict(zip([second.symbol[x] for x in second.coords], [first.symbol[x] for x in first.coords]))
        for i in range(4):
            for j in range(i, 4):
                gap = (a[i, j] - b[i, j].subs(names)).subs(at)
                for point in ({first.symbol["r"]: sp.Rational(5, 2), first.symbol["\\theta"]: sp.Rational(7, 10)},
                              {first.symbol["r"]: sp.Rational(-4, 3), first.symbol["\\theta"]: sp.Rational(11, 5)}):
                    value = gap.subs(point).subs(first.c, 1).subs(second.c, 1)
                    self.assertEqual(sp.simplify(value), 0, (i, j))

    def test_without_charge_it_is_the_published_taub_nut_metric(self):
        sp = self.sp
        _, mine, g = self.published(ID, "spherical")
        _, theirs, h = self.published("taub_nut", "spherical")
        m, l = sp.Rational(3, 2), sp.Rational(2, 3)
        self.same(g, h, (mine, theirs), {mine.parameters["m"]: m, mine.parameters["l"]: l, mine.parameters["r_q"]: 0,
                                         theirs.parameters["m"]: m, theirs.parameters["l"]: l})

    def test_without_twist_it_is_the_published_reissner_nordstrom_metric(self):
        sp = self.sp
        _, mine, g = self.published(ID, "spherical")
        _, theirs, h = self.published("rn_metric", "spherical")
        m, q = sp.Rational(3, 2), sp.Rational(4, 5)
        self.same(g, h, (mine, theirs), {mine.parameters["m"]: m, mine.parameters["l"]: 0, mine.parameters["r_q"]: q,
                                         theirs.parameters["r_s"]: 2 * m, theirs.parameters["r_q"]: q})

    def test_the_kretschmann_scalar_is_finite_at_the_centre_while_there_is_a_twist(self):
        # At r = 0 the published scalar is (56 r_q^4 - 96 l^2 r_q^2 - 48 l^2 (m^2 - l^2))/l^8, over a
        # power of l alone, and with l = 0 it is Reissner and Nordstrom's, which diverges there.
        sp = self.sp
        entry, reader, _ = self.published(ID, "spherical")
        K = reader(entry["kretschmann"].partition("=")[2])
        r, m, l, q = reader.symbol["r"], *(reader.parameters[k] for k in ("m", "l", "r_q"))
        at_centre = sp.simplify(K.subs(r, 0))
        self.assertEqual(sp.simplify(at_centre - (56 * q ** 4 - 96 * l ** 2 * q ** 2 - 48 * l ** 2 * (m ** 2 - l ** 2))
                                     / l ** 8), 0)
        flat = sp.simplify(K.subs(l, 0))
        self.assertEqual(sp.simplify(flat - (48 * m ** 2 * r ** 2 - 96 * m * q ** 2 * r + 56 * q ** 4) / r ** 8), 0)
        # Reissner and Nordstrom's own published scalar, at r_s = 2m.
        theirs_entry, theirs, _ = self.published("rn_metric", "spherical")
        rn = theirs(theirs_entry["kretschmann"].removeprefix("K = ")).subs(
            {theirs.parameters["r_s"]: 2 * m, theirs.parameters["r_q"]: q, theirs.symbol["r"]: r})
        self.assertEqual(sp.simplify(flat - rn), 0)

    def test_regrouping_a_sum_changes_no_value(self):
        sp, pc = self.sp, self.pc
        spec = pc.brill_charged_taub_nut("spherical")
        reader = self.vm.Reader(spec["system"]["coords"], spec["system"]["parameters"], ())
        r = reader.symbol["r"]
        m, l, q = (reader.parameters[k] for k in ("m", "l", "r_q"))
        Sigma, Delta = sp.Symbol("Sigma", positive=True), sp.Symbol("Delta")
        forms = pc.BrillForms(reader, "r", [(Sigma, r ** 2 + l ** 2), (Delta, sp.expand(r ** 2 - 2 * m * r - l ** 2 + q ** 2))])
        names = {Sigma: r ** 2 + l ** 2, Delta: r ** 2 - 2 * m * r - l ** 2 + q ** 2}
        C = forms.C
        for base, terms in ((r ** 4 + 6 * l ** 2 * r ** 2 - 8 * m * l ** 2 * r - 3 * l ** 4 + 4 * l ** 2 * q ** 2, 2),
                            (m * r ** 2 + 2 * l ** 2 * r - q ** 2 * r - m * l ** 2, 3),
                            (sp.expand((r ** 2 + l ** 2) ** 2 * (C ** 2 - 1)
                                       + 4 * l ** 2 * (r ** 2 - 2 * m * r - l ** 2 + q ** 2) * C ** 2), 2)):
            number, grouped = forms.regrouped(sp.expand(base))
            self.assertEqual(len(sp.Add.make_args(grouped)), terms)
            written = (number * grouped).subs(names).subs(forms.s, sp.sqrt(1 - C ** 2))
            self.assertEqual(sp.expand(written - base), 0)
        # A sum that is no shorter in the names is left as it stands.
        self.assertIsNone(forms.regrouped(sp.expand(r ** 3 + m)))

    def test_the_tortoise_coordinate_has_the_slope_sigma_over_delta(self):
        for case in ("black_hole", "wormhole"):
            m, l, q = self.slices.BRILL[case]
            star = self.slices.brill_rstar(case)
            self.assertAlmostEqual(float(star(0.0)), 0.0, places=14)
            for r in (-3.0, -1.0, 1.0, 3.0, 7.5):
                slope = float(star(r + 1e-6) - star(r - 1e-6)) / 2e-6
                self.assertAlmostEqual(slope, (r * r + l * l) / (r * r - 2 * m * r - l * l + q * q), places=6)

    def test_brills_universe_prints_again_as_published(self):
        pc = self.pc
        published = pc.METRICS / f"{ID}.json"
        with tempfile.TemporaryDirectory() as folder:
            copy = Path(folder) / published.name
            copy.write_bytes(published.read_bytes())
            kept, pc.METRICS = pc.METRICS, Path(folder)
            try:
                pc.write(pc.brill_charged_taub_nut("taub"))
            finally:
                pc.METRICS = kept
            self.assertEqual(copy.read_bytes(), published.read_bytes())


class TheDrawings(unittest.TestCase):
    """The diagram files, from their numbers alone."""

    HOLE = (1.0, 0.75, 1.0)          # m, l and r_q of the black hole, in units of m
    WORMHOLE = (0.0, 1.0, 1.5)       # and of the wormhole, in units of l

    def horizons(self, case):
        m, l, q = case
        inside = m * m + l * l - q * q
        return [] if inside < 0 else [m - math.sqrt(inside), m + math.sqrt(inside)]

    def test_the_two_cases_stand_on_either_side_of_the_extreme_charge(self):
        self.assertEqual(self.horizons(self.HOLE), [0.25, 1.75])
        self.assertEqual(self.horizons(self.WORMHOLE), [])

    def test_every_view_marks_the_horizons_of_its_case_and_no_others(self):
        views = read("diagrams")["systems"]
        seen = 0
        for system, drawn in views.items():
            for view in drawn:
                if system == "taub":
                    continue                                    # its box ends on both horizons
                x0, x1 = view["box"][:2]
                marked = sorted(x0 + line[0][0] * (x1 - x0) for marker in view["markers"] if marker["kind"] == "grr"
                                for line in marker["lines"])
                wanted = self.horizons(self.WORMHOLE if view["id"] == "wormhole" else self.HOLE)
                self.assertEqual(len(marked), len(wanted), (system, view["id"]))
                for a, b in zip(marked, wanted):
                    self.assertAlmostEqual(a, b, delta=0.01, msg=(system, view["id"]))
                seen += 1
        self.assertEqual(seen, 7)

    def test_brills_universe_is_drawn_between_the_horizons(self):
        view, = read("diagrams")["systems"]["taub"]
        self.assertEqual(view["box"][2:], self.horizons(self.HOLE))
        # The angle is drawn as the length 2 l psi, at l = 3m/4.
        self.assertEqual(view["to_display"][0], [0.0, 2 * self.HOLE[1]])

    def test_the_embedded_circles_have_the_radius_of_the_metric(self):
        for view in read("embedding")["views"]:
            l = (self.WORMHOLE if view["id"] == "wormhole" else self.HOLE)[1]
            for surface in view["surfaces"]:
                self.assertTrue(surface["rings"])
                for ring in surface["rings"]:
                    self.assertAlmostEqual(ring["rho"], math.sqrt(ring["x"] ** 2 + l * l), places=5)

    def test_the_wormholes_throat_has_the_radius_l_and_its_two_sides_are_mirror_images(self):
        view = next(v for v in read("embedding")["views"] if v["id"] == "wormhole")
        rings = {ring["x"]: ring for ring in view["surfaces"][0]["rings"]}
        self.assertAlmostEqual(rings[0.0]["rho"], 1.0, places=6)
        self.assertAlmostEqual(rings[0.0]["z"], 0.0, places=6)
        self.assertEqual(min(ring["rho"] for ring in rings.values()), rings[0.0]["rho"])
        for x, ring in rings.items():
            self.assertAlmostEqual(ring["z"], -rings[-x]["z"], places=6)
            self.assertAlmostEqual(ring["rho"], rings[-x]["rho"], places=6)

    def test_the_black_holes_surface_starts_on_its_outer_horizon(self):
        view = next(v for v in read("embedding")["views"] if v["id"] == "equator")
        inner = min(ring["x"] for ring in view["surfaces"][0]["rings"])
        self.assertAlmostEqual(inner, 1.75, places=9)


if __name__ == "__main__":
    unittest.main()
