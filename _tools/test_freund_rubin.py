"""Freund and Rubin's anti-de Sitter space of four dimensions times a 7-sphere of radius 2L:
.venv.noindex/bin/python -m unittest discover -s _tools

What its texts state, held on the published files: every chart has eleven dimensions, the mixed
Ricci tensor is -3/L^2 on the four large dimensions and 3/(2L^2) on the sphere, the ratio -2 of
Freund and Rubin's (7a) for each dimension, which is their field equation for a four-form f eps
on the anti-de Sitter factor at 8 pi G f^2 = 3/(4L^2); the Ricci scalar is -3/(2L^2) and the
Kretschmann scalar 117/(4L^4); the static, conformal and proper distance charts are the global or
Poincare chart carried along r = L sinh(rho), tan(chi) = sinh(rho) and r = L e^(sigma/L); and the
drawings hold what their captions say. The tests of the published mathematics need sympy and are
skipped where it is absent; the tests of the drawings read the files alone.
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
CHARTS = ["global", "conformal", "static", "poincare", "proper"]
SEVEN = ["\\alpha", "\\beta", "\\gamma", "\\kappa", "\\xi", "\\omega", "\\psi"]


def published(metric_id="freund_rubin"):
    return json.loads((DATA / "metrics" / f"{metric_id}.json").read_text(encoding="utf-8"))


def drawing(kind, metric_id="freund_rubin"):
    return json.loads((DATA / kind / f"{metric_id}.json").read_text(encoding="utf-8"))


@unittest.skipUnless(HAS_SYMPY, "sympy is not installed")
class Published(unittest.TestCase):
    def setUp(self):
        import sympy
        import verify_metrics
        self.sp, self.vm = sympy, verify_metrics

    def read(self, system):
        """The reader of a published chart and a function from a field's name to the diagonal of
        its published components with c = 1."""
        entry = next(c for c in published()["coordinates"] if c["id"] == system)
        reader = self.vm.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]], ())
        names = entry["coords"]

        def diagonal(field, variant=None):
            block = entry[field] if variant is None else entry[field]["variants"][variant]["nonzero"]
            out = [self.sp.Integer(0)] * len(names)
            for e in block:
                i, j = (names.index(k) for k in e["indices"])
                self.assertEqual(i, j, f"{system}: {field} has a component off the diagonal")
                out[i] = reader(e["value"]).subs(reader.c, 1)
            return out
        return reader, diagonal, entry

    def at(self, reader, L=1.3):
        """Numbers for every coordinate and for L."""
        values = iter((0.3, 0.7, 1.1, 0.9, 1.3, 0.4, 0.8, 1.2, 0.6, 1.4, 0.5))
        out = {s: next(values) for s in reader.symbol.values()}
        out[reader.parameters["L"]] = L
        return out

    def number(self, value, at):
        return float(self.sp.N(self.sp.sympify(value).subs(at)))

    def test_every_chart_has_eleven_dimensions_and_a_7_sphere(self):
        for system in CHARTS:
            reader, diagonal, entry = self.read(system)
            self.assertEqual(len(entry["coords"]), 11, system)
            self.assertEqual(entry["coords"][4:], SEVEN, system)
            self.assertEqual(entry["ricci_scalar"], "R = -\\dfrac{3}{2L^2}", system)
            self.assertEqual(entry["kretschmann"], "K = \\dfrac{117}{4L^4}", system)

    def test_the_mixed_ricci_tensor_is_freund_and_rubin_s(self):
        L = 1.3
        for system in CHARTS:
            reader, diagonal, entry = self.read(system)
            at = self.at(reader, L)
            mixed = [self.number(v, at) for v in diagonal("ricci_tensor", "ul")]
            for i, value in enumerate(mixed):
                self.assertAlmostEqual(value, -3 / L ** 2 if i < 4 else 1.5 / L ** 2, 12, f"{system} slot {i}")
            # Freund and Rubin's (7a) at d = 11, s = 4: each large dimension carries -2 times what each
            # small one does.
            self.assertAlmostEqual(mixed[0] / mixed[10], -2, 12, system)

    def test_the_four_form_solves_the_field_equations(self):
        """R_MN = 8 pi G (F_MPQR F_N^PQR - g_MN F^2/12) with F_MPQR F_N^PQR = -6 f^2 g_MN on the
        four large dimensions and F^2 = -24 f^2, at 8 pi G f^2 = 3/(4 L^2)."""
        L = 1.3
        coupling = 3 / (4 * L ** 2)
        for system in CHARTS:
            reader, diagonal, entry = self.read(system)
            at = self.at(reader, L)
            ricci = [self.number(v, at) for v in diagonal("ricci_tensor", "ll")]
            metric = [self.number(v, at) for v in diagonal("metric_components")]
            for i in range(11):
                source = coupling * ((-6 * metric[i] if i < 4 else 0) + 24 * metric[i] / 12)
                self.assertAlmostEqual(ricci[i], source, 10, f"{system} slot {i}")

    def test_the_sphere_has_twice_the_radius_of_anti_de_sitter_space(self):
        for system in CHARTS:
            reader, diagonal, entry = self.read(system)
            at = self.at(reader, 1.0)
            metric = diagonal("metric_components")
            self.assertAlmostEqual(self.number(metric[4], at), 4.0, 12, system)

    def test_the_charts_are_one_another_carried_along(self):
        sp = self.sp
        L = 1.3
        _, glob, _ = self.read("global")
        greader, _, _ = self.read("global")
        rho = greader.symbol["\\rho"]
        for system, image in (("static", lambda x: sp.asinh(x / L)), ("conformal", lambda x: sp.asinh(sp.tan(x)))):
            reader, diagonal, entry = self.read(system)
            x = reader.symbol[entry["coords"][1]]
            mine = diagonal("metric_components")
            theirs = glob("metric_components")
            for value in (0.3, 0.9, 1.4):
                at = self.at(reader, L)
                at[x] = value
                there = {s: at.get(s, 0.5) for s in greader.symbol.values()}
                there.update({greader.symbol[n]: at[reader.symbol[n]] for n in entry["coords"] if n in greader.symbol})
                there[rho] = float(image(value))
                there[greader.parameters["L"]] = L
                slope = float(sp.diff(image(x), x).subs(x, value))
                for i in range(11):
                    want = self.number(theirs[i], there) * (slope ** 2 if i == 1 else 1)
                    self.assertAlmostEqual(self.number(mine[i], at), want, 10, f"{system} slot {i} at {value}")
        _, poincare, pentry = self.read("poincare")
        preader, _, _ = self.read("poincare")
        reader, diagonal, entry = self.read("proper")
        sigma = reader.symbol["\\sigma"]
        for value in (-1.0, 0.0, 0.8):
            at = self.at(reader, L)
            at[sigma] = value
            there = {preader.symbol[n]: at[reader.symbol[n]] for n in entry["coords"] if n in preader.symbol}
            there[preader.symbol["r"]] = L * math.exp(value / L)
            there[preader.parameters["L"]] = L
            slope = math.exp(value / L)
            mine, theirs = diagonal("metric_components"), poincare("metric_components")
            for i in range(11):
                want = self.number(theirs[i], there) * (slope ** 2 if i == 3 else 1)
                self.assertAlmostEqual(self.number(mine[i], at), want, 10, f"proper slot {i} at {value}")


class Drawings(unittest.TestCase):
    def test_every_chart_has_a_spacetime_diagram_and_a_conformal_diagram(self):
        diagrams, conformal = drawing("diagrams"), drawing("conformal")
        self.assertEqual(sorted(diagrams["systems"]), sorted(CHARTS))
        self.assertEqual(sorted(v["id"] for v in conformal["views"]), sorted(CHARTS))

    def test_the_embedding_has_the_hyperbolic_plane_a_cylinder_and_a_sphere(self):
        views = {v["id"] for v in drawing("embedding")["views"]}
        self.assertEqual(views, {"anti_de_sitter", "circle", "sphere"})

    def test_a_ray_round_the_sphere_takes_four_times_a_crossing_of_anti_de_sitter_space(self):
        """On the plane of t and psi at sigma = 0 the metric is -c^2dt^2 + 4L^2dpsi^2: once round
        the great circle is ct = 4 pi L, the box's height, and boundary to centre and back out of
        anti-de Sitter space is ct = pi L, which the captions state."""
        view = next(v for v in drawing("diagrams")["systems"]["proper"] if v["id"] == "circle")
        x0, x1, y0, y1 = view["box"]
        self.assertAlmostEqual(x1 - x0, 4 * math.pi, 9)
        self.assertAlmostEqual(y1 - y0, 4 * math.pi, 9)
        caption = " ".join(view["caption"])
        self.assertIn("$4\\pi L/c$", caption)
        self.assertIn("four times", caption)

    def test_the_moment_is_marked_at_the_radius_each_chart_gives_it(self):
        """The hyperbolic plane reaches rho = 2, which is chi = arctan(sinh 2) and r = sinh 2 L; the
        cylinder reaches sigma = -L to L, which is r = L/e to e L."""
        systems = drawing("diagrams")["systems"]
        want = {"global": (0.0, 2.0), "conformal": (0.0, math.atan(math.sinh(2))), "static": (0.0, math.sinh(2)),
                "poincare": (math.exp(-1), math.exp(1)), "proper": (-1.0, 1.0)}
        for system, (lo, hi) in want.items():
            view = next(v for v in systems[system] if v["id"] == "radial")
            x0, x1, y0, y1 = view["box"]
            # The lines are stored in the unit square of the box, the drawn coordinate across.
            ends = sorted(x0 + p[0] * (x1 - x0) for line in view["slices"][0]["lines"] for p in line)
            self.assertTrue(all(abs(y0 + p[1] * (y1 - y0)) < 1e-3 for line in view["slices"][0]["lines"] for p in line))
            self.assertAlmostEqual(ends[0], lo, 3, system)
            self.assertAlmostEqual(ends[-1], hi, 3, system)


if __name__ == "__main__":
    unittest.main()
