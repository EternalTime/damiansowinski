"""Visser's Roman ring of wormholes, Phys. Rev. D 55, 5212 (1997), gr-qc/9702043:
.venv.noindex/bin/python -m unittest discover -s _tools

What its chart, texts and drawings state, held on the published files. The chart's published metric
is read through the checker's own Reader and held at points to Minkowski's, with every curvature
component and scalar empty; Visser's van Vleck determinant at a throat to his Chebyshev form and its
limit for crossings long beside the mouths; the mouths of the drawn ring to a polygon of sides L and
ell; the History to its references, its shape and its people named in full; the related spacetimes to
answering; and the drawings to the one figure the ring has. The tests that read the mathematics need
sympy and numpy, and are skipped where they are absent.
"""
import importlib.util
import json
import re
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "MFS" / "assets" / "data"
DERIVATIONS = Path(__file__).resolve().parent / "derivations"
HAS_SYMPY = all(importlib.util.find_spec(m) is not None for m in ("sympy", "numpy", "mpmath"))
sys.path.insert(0, str(DERIVATIONS))
sys.path.insert(0, str(ROOT / "_tools"))
ME = "roman_ring"
METRIC = json.loads((DATA / "metrics" / f"{ME}.json").read_text(encoding="utf-8"))


def chart(system="rest"):
    return next(c for c in METRIC["coordinates"] if c["id"] == system)


def cited(text):
    keys = []
    for group in re.findall(r"\[([a-z][^\]]*)\]", text):
        for key in (k.strip() for k in group.split(",")):
            if key not in keys:
                keys.append(key)
    return keys


@unittest.skipUnless(HAS_SYMPY, "needs sympy, numpy and mpmath")
class Chart(unittest.TestCase):
    """The rest frame of the mouths is Minkowski space, flat outside the mouths."""

    def test_published_metric_is_minkowskis(self):
        import sympy as sp
        import verify_metrics as vm
        c = chart()
        reader = vm.Reader(c["coords"], [p["symbol"] for p in c["parameters"]], ())
        g = sp.zeros(4, 4)
        for comp in c["metric_components"]:
            i, j = (c["coords"].index(x) for x in comp["indices"])
            g[i, j] = g[j, i] = reader(comp["value"])
        self.assertEqual(g, sp.diag(-1, 1, 1, 1))

    def test_every_curvature_component_vanishes(self):
        c = chart()
        for field in ("christoffel", "riemann", "ricci_tensor", "einstein_tensor", "weyl_tensor"):
            for name, variant in c[field]["variants"].items():
                self.assertEqual(variant["nonzero"], [], (field, name))
        self.assertEqual((c["ricci_scalar"], c["kretschmann"]), ("R = 0", "K = 0"))
        self.assertIn("\\text{flat outside the mouths: every curvature tensor vanishes}", c["domains"])

    def test_identification_joins_each_entrance_at_t_plus_T_to_its_exit_at_t(self):
        c = chart()
        self.assertIn("\\mathcal{A}_k(t + T) \\sim \\mathcal{B}_k(t) \\;\\text{(each pair of spheres identified)}",
                      c["domains"])
        T = next(p for p in c["parameters"] if p["symbol"] == "T")
        self.assertIn("the entrance mouth at time $t + T$ is the exit mouth at time $t$", T["description"])


@unittest.skipUnless(HAS_SYMPY, "needs sympy, numpy and mpmath")
class VanVleck(unittest.TestCase):
    """Visser 1994 eqs. (65)-(66) and (115)-(116) against Visser 1997 eqs. (6)-(7)."""

    def setUp(self):
        import sympy as sp
        self.sp = sp

    def P(self, n, z):
        return sum(self.sp.binomial(n + j + 1, 2 * j + 1) * z ** j for j in range(n + 1))

    def test_the_polynomials_are_chebyshev_polynomials_at_one_plus_half_z(self):
        sp, z = self.sp, self.sp.Symbol("z")
        for n in range(10):
            self.assertEqual(sp.expand(self.P(n, z) - sp.chebyshevu(n, 1 + z / 2)), 0, n)

    def test_the_argument_printed_in_eq_66_fails_at_n_1(self):
        sp, z = self.sp, self.sp.Symbol("z")
        self.assertNotEqual(sp.expand(self.P(1, z) - sp.chebyshevu(1, (z + 4) / 2)), 0)

    def test_the_1994_determinant_at_a_throat_is_the_1997_one(self):
        sp = self.sp
        for N in range(1, 8):
            for q in (sp.Rational(1, 3), 2, sp.Rational(77, 5)):
                old = N ** 2 / self.P(N - 1, 2 * q) ** 2
                new = (N / sp.chebyshevu(N - 1, 1 + q)) ** 2
                self.assertEqual(sp.simplify(old - new), 0, (N, q))

    def test_eq_7_is_the_limit_of_eq_6_and_more_wormholes_make_it_smaller(self):
        sp = self.sp
        for N in range(2, 8):
            q = 10 ** 6
            exact = (N / sp.chebyshevu(N - 1, 1 + q)) ** 2
            self.assertAlmostEqual(float(exact / (N ** 2 * sp.Rational(1, 2 * q) ** (2 * (N - 1)))), 1.0, places=4)
        for q in (2, 5, 40):
            values = [float((N / sp.chebyshevu(N - 1, 1 + q)) ** 2) for N in range(1, 9)]
            self.assertEqual(values[0], 1.0)
            self.assertTrue(all(b < a for a, b in zip(values, values[1:])), q)


@unittest.skipUnless(HAS_SYMPY, "needs sympy, numpy and mpmath")
class Mouths(unittest.TestCase):
    """The ring the figure draws, and the polygon the chart's check builds."""

    def test_the_drawn_polygon_has_sides_L_and_ell_and_matches_the_charts(self):
        import numpy as np
        import projections as pj
        import print_charts as pc
        for N, L, ell in ((2, 5.0, 1.0), (4, 4.0, 1.0), (7, 9.0, 1.5)):
            a, b, rho = pj.roman_ring_mouths(N, L, ell)
            self.assertLess(np.abs(np.linalg.norm(b - a, axis=1) - L).max(), 1e-12)
            self.assertLess(np.abs(np.linalg.norm(np.roll(a, -1, axis=0) - b, axis=1) - ell).max(), 1e-12)
            self.assertLess(np.abs(np.linalg.norm(np.vstack([a, b]), axis=1) - rho).max(), 1e-12)
            sa, sb, srho = pc.roman_ring_polygon(N, L, ell)
            self.assertLess(np.abs(np.array(sa, dtype=float) - a).max(), 1e-12)
            self.assertLess(np.abs(np.array(sb, dtype=float) - b).max(), 1e-12)

    def test_the_drawn_ring_is_short_of_a_time_machine(self):
        import projections as pj
        v = {k: float(eval(str(x))) for k, x in pj.ROMAN_RING.items()}
        self.assertLess(v["T"], v["ell"])
        self.assertLess(v["ell"], v["L"])
        self.assertLess(v["R"], v["ell"] / 2)
        self.assertEqual(v["N"], 4)


class Texts(unittest.TestCase):
    """The page's fields, references and relations."""

    def test_references_in_order_of_first_citation_and_in_the_bibliography(self):
        order = cited(METRIC["history"])
        for entry in METRIC["related"]:
            for key in cited(entry["text"]):
                if key not in order:
                    order.append(key)
        self.assertEqual(order, METRIC["references"])
        bib = (ROOT / "assets" / "data" / "references.bib").read_text(encoding="utf-8")
        for key in order:
            self.assertEqual(len(re.findall(r"^@\w+\{" + re.escape(key) + ",", bib, re.M)), 1, key)

    def test_history_shape(self):
        history = METRIC["history"]
        self.assertEqual(len(history.split("¶")), 6)
        self.assertTrue(history.startswith("A traveller who steps into one wormhole of a Roman ring"))
        self.assertLessEqual(len(history.replace("¶", "")), 6900)

    def test_people_named_in_full_first(self):
        history = METRIC["history"]
        for full, short in (("Matt Visser", "Visser wrote"), ("Stephen Hawking", "Cassidy and Hawking"),
                            ("Thomas Roman", "Roman's spacetime"), ("Maxim Lyutikov", "Lyutikov"),
                            ("Michael Cassidy", "Cassidy and Hawking")):
            self.assertLessEqual(history.index(full), history.index(short), full)

    def test_tags_the_mouths_decide(self):
        tags = set(METRIC["tags"])
        self.assertTrue({"locally flat", "stationary", "wormhole", "time machine"} <= tags)
        self.assertFalse({"flat", "vacuum", "conformally flat", "static"} & tags)

    def test_related_spacetimes_answer(self):
        kinds = {"locally_same": "locally_same", "family": "family"}
        for entry in METRIC["related"]:
            other = json.loads((DATA / "metrics" / f"{entry['id']}.json").read_text(encoding="utf-8"))
            answer = [e for e in other["related"] if e["id"] == ME]
            self.assertEqual(len(answer), 1, entry["id"])
            self.assertEqual(answer[0]["kind"], kinds[entry["kind"]], entry["id"])
            for key in cited(answer[0]["text"]):
                self.assertIn(key, other["references"], (entry["id"], key))


class Drawings(unittest.TestCase):
    """The ring has one figure, and no conformal, embedding or curvature drawing."""

    def test_the_chart_has_its_figure_in_three_dimensions(self):
        data = json.loads((DATA / "diagrams" / f"{ME}.json").read_text(encoding="utf-8"))
        self.assertEqual(data.get("systems", {}), {})
        self.assertEqual([v["id"] for v in data["projections"]["rest"]], ["ring"])

    def test_no_conformal_embedding_or_curvature_file(self):
        for kind in ("conformal", "embedding", "curvature"):
            self.assertFalse((DATA / kind / f"{ME}.json").exists(), kind)

    def test_the_absent_drawings_are_named_with_their_reasons(self):
        for script in ("conformal.py", "embedding.py"):
            text = (DERIVATIONS / script).read_text(encoding="utf-8")
            self.assertIn('"roman_ring"', text, script)
            self.assertIn("roman_ring.md", text, script)
        build = (ROOT / "_tools" / "build_mfs_data.py").read_text(encoding="utf-8")
        self.assertRegex(build, r'WITHOUT_EMBEDDING = frozenset\(\{[^}]*"roman_ring"')


if __name__ == "__main__":
    unittest.main()
