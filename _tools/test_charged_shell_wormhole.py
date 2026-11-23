"""Schein and Aichelburg's wormhole between charged shells, gr-qc/9606069 and Acta Phys. Pol. B 29, 1025:
.venv.noindex/bin/python -m unittest discover -s _tools

What its charts, texts and drawings state, held on the published files. Each chart's published metric
is read through the checker's own Reader and held at points to the founders' (metric), (bisph),
(sphcoor) and (RN); the potential's image series to taking one value on both spheres, to the
bispherical series and to the maximum principle; the shells to one areal radius from both sides and
to a positive energy density at the drawn parameters, by the founders' formula for it; the History
to its references and the bibliography; the related spacetimes to answering; and the drawings to
the views, shells and moments they draw. The tests that read the mathematics need sympy and mpmath,
and are skipped where they are absent.
"""
import importlib.util
import json
import re
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "MFS" / "assets" / "data"
HAS_SYMPY = importlib.util.find_spec("sympy") is not None and importlib.util.find_spec("mpmath") is not None
sys.path.insert(0, str(Path(__file__).resolve().parent / "derivations"))
ME = "charged_shell_wormhole"
METRIC = json.loads((DATA / "metrics" / f"{ME}.json").read_text(encoding="utf-8"))
SYSTEMS = ("cartesian", "bispherical", "polar", "interior")
# The drawn configuration, R = 1: the spheres at z = +-2R, m_1 = R/2, the black hole m = R/4, q = 6R/25.
R, D1, M1, M, Q = 1, 2, 0.5, 0.25, 0.24


def chart(system, metric=METRIC):
    return next(c for c in metric["coordinates"] if c["id"] == system)


def drawing(kind):
    return json.loads((DATA / kind / f"{ME}.json").read_text(encoding="utf-8"))


def cited(text):
    keys = []
    for group in re.findall(r"\[([a-z][^\]]*)\]", text):
        for key in (k.strip() for k in group.split(",")):
            if key not in keys:
                keys.append(key)
    return keys


@unittest.skipUnless(HAS_SYMPY, "needs sympy and mpmath")
class Charts(unittest.TestCase):
    """Each chart's published metric at points against the founders' line elements."""

    @classmethod
    def setUpClass(cls):
        import sympy as sp
        import verify_metrics as vm
        cls.sp, cls.vm = sp, vm

    def published(self, system, metric=METRIC):
        c = chart(system, metric)
        reader = self.vm.Reader(c["coords"], [p["symbol"] for p in c["parameters"]], ())
        values = {tuple(e["indices"]): e["value"] for e in c["metric_components"]}
        n = len(c["coords"])
        g = self.sp.Matrix(n, n, lambda i, j: reader(values.get((c["coords"][i], c["coords"][j]), "0")))
        return reader, g

    def test_cartesian_is_majumdar_and_papapetrou(self):
        sp = self.sp
        reader, g = self.published("cartesian")
        V = reader.parameters["V"]
        want = sp.diag(-1 / V ** 2, V ** 2, V ** 2, V ** 2)
        for i in range(4):
            for j in range(4):
                self.assertEqual(self.vm.norm(g[i, j] - want[i, j]), 0, (i, j))

    def test_bispherical_is_flat_space_pulled_back(self):
        sp = self.sp
        reader, g = self.published("bispherical")
        mu, eta, phi = (reader.symbol[x] for x in ("\\mu", "\\eta", "\\phi"))
        a = reader.parameters["a"]
        W = sp.cosh(mu) - sp.cos(eta)
        image = [a * sp.sin(eta) * sp.cos(phi) / W, a * sp.sin(eta) * sp.sin(phi) / W, a * sp.sinh(mu) / W]
        J = sp.Matrix(3, 3, lambda i, j: sp.diff(image[i], (mu, eta, phi)[j]))
        flat = J.T * J
        spatial = g[1:, 1:].subs(reader.parameters["V"], 1)
        if "W" in reader.parameters:
            spatial = spatial.subs(reader.parameters["W"], W)
        for at in ({mu: sp.Rational(1, 3), eta: sp.Rational(7, 5), phi: sp.Rational(1, 2), a: sp.sqrt(3)},
                   {mu: sp.Rational(-9, 10), eta: sp.Rational(2, 5), phi: 3, a: sp.Rational(3, 2)}):
            for i in range(3):
                for j in range(3):
                    self.assertLess(abs(complex(sp.N((spatial[i, j] - flat[i, j]).subs(at), 30))), 1e-25)

    def test_polar_is_flat_space_about_one_centre(self):
        sp = self.sp
        reader, g = self.published("polar")
        rho, th, ph = (reader.symbol[x] for x in ("\\rho", "\\theta", "\\phi"))
        d1 = reader.parameters["d_1"]
        image = [rho * sp.sin(th) * sp.cos(ph), rho * sp.sin(th) * sp.sin(ph), -d1 + rho * sp.cos(th)]
        J = sp.Matrix(3, 3, lambda i, j: sp.diff(image[i], (rho, th, ph)[j]))
        flat = J.T * J
        spatial = g[1:, 1:].subs(reader.parameters["V"], 1)
        for i in range(3):
            for j in range(3):
                self.assertEqual(self.vm.norm(spatial[i, j] - flat[i, j]), 0, (i, j))

    def test_interior_is_reissner_and_nordstrom(self):
        sp = self.sp
        reader, g = self.published("interior")
        rn_metric = json.loads((DATA / "metrics" / "rn_metric.json").read_text(encoding="utf-8"))
        other, h = self.published("spherical", rn_metric)
        swap = {other.symbol[x]: reader.symbol[y] for x, y in zip(("t", "r", "\\theta", "\\phi"),
                                                                    ("t", "r", "\\theta", "\\phi"))}
        swap.update({other.parameters["r_s"]: 2 * reader.parameters["m"], other.parameters["r_q"]: reader.parameters["q"]})
        for i in range(4):
            for j in range(4):
                self.assertEqual(self.vm.norm(h[i, j].subs(swap) - g[i, j]), 0, (i, j))

    def test_interior_kretschmann_is_reissner_and_nordstroms(self):
        sp = self.sp
        c = chart("interior")
        reader = self.vm.Reader(c["coords"], [p["symbol"] for p in c["parameters"]], ())
        K = reader(re.sub(r"^K\s*=\s*", "", c["kretschmann"]))
        r, m, q = reader.symbol["r"], reader.parameters["m"], reader.parameters["q"]
        want = 48 * m ** 2 / r ** 6 - 96 * m * q ** 2 / r ** 7 + 56 * q ** 4 / r ** 8
        self.assertEqual(self.vm.norm(K - want), 0)


@unittest.skipUnless(HAS_SYMPY, "needs sympy and mpmath")
class Potential(unittest.TestCase):
    """The founders' image series, in forty digits."""

    @classmethod
    def setUpClass(cls):
        import mpmath as mp
        mp.mp.dps = 40
        cls.mp = mp

    def images(self, R, d1, m1, terms=120):
        mp = self.mp
        mu0, a = mp.acosh(mp.mpf(d1) / R), mp.sqrt(mp.mpf(d1) ** 2 - R ** 2)
        return [((-1) ** (n + 1) * m1 * mp.sinh(mu0) / mp.sinh(n * mu0), a * mp.coth(n * mu0)) for n in range(1, terms + 1)]

    def V(self, rho, z, R=R, d1=D1, m1=M1):
        mp = self.mp
        return 1 + sum(mn * (1 / mp.sqrt(rho ** 2 + (z - dn) ** 2) + 1 / mp.sqrt(rho ** 2 + (z + dn) ** 2))
                       for mn, dn in self.images(R, d1, m1))

    def test_one_value_on_both_spheres(self):
        mp = self.mp
        for R_, d1, m1 in ((1, 2, 0.5), (1, 1.2, 0.3), (0.5, 3, 1)):
            for th in (0, 0.7, mp.pi / 2, 2.4, mp.pi):
                for centre in (-d1, d1):
                    v = self.V(R_ * mp.sin(th), centre + R_ * mp.cos(th), R_, d1, m1)
                    self.assertLess(abs(v - (1 + mp.mpf(m1) / R_)), 1e-12, (R_, d1, m1, th, centre))

    def test_one_far_away_and_largest_on_the_shells(self):
        self.assertLess(abs(self.V(0, 1e6) - 1), 1e-5)
        for rho, z in ((0, 0), (0.5, 0.2), (3, 0), (0, 3.5), (1.5, 2.2), (6, -4)):
            v = self.V(rho, z)
            self.assertTrue(1 < v < 1 + M1 / R, (rho, z, v))

    def test_bispherical_series_is_the_cartesian_series(self):
        mp = self.mp
        mu0, a = mp.acosh(2), mp.sqrt(3)
        for mu, eta in ((0.3, 1.1), (-1.0, 2.9), (1.25, 0.4)):
            W = mp.cosh(mu) - mp.cos(eta)
            rho, z = a * mp.sin(eta) / W, a * mp.sinh(mu) / W
            series = 1 + M1 / R * (1 + mp.sqrt(W) * sum((-1) ** (n + 1) / mp.sqrt(mp.cosh(mu + 2 * n * mu0) - mp.cos(eta))
                                                       for n in range(-120, 121)))
            self.assertLess(abs(series - self.V(rho, z)), 1e-15, (mu, eta))

    def sigma(self, R, d1, m1, m, q, th):
        """The founders' surface energy density at the polar angle th of the shell about z = -d_1,
        th = 0 toward the other sphere: -(1/4 pi V_0^2) dV/dr_+ - (1/4 pi R V_0)(1 - sqrt(f(R V_0)))."""
        mp = self.mp
        V0 = 1 + mp.mpf(m1) / R
        px, pz = R * mp.sin(th), -d1 + R * mp.cos(th)
        nx, nz = mp.sin(th), mp.cos(th)
        dV = 0
        for mn, dn in self.images(R, d1, m1, 400):
            for c in (dn, -dn):
                r = mp.sqrt(px ** 2 + (pz - c) ** 2)
                dV -= mn * (px * nx + (pz - c) * nz) / r ** 3
        x = R * V0
        return -dV / (4 * mp.pi * V0 ** 2) - (1 - mp.sqrt(1 - 2 * m / x + q ** 2 / x ** 2)) / (4 * mp.pi * x)

    def test_energy_density_positive_at_the_drawn_parameters(self):
        mp = self.mp
        for th in (0, mp.pi / 4, mp.pi / 2, 3 * mp.pi / 4, mp.pi):
            self.assertGreater(self.sigma(R, D1, M1, M, Q, th), 0, th)

    def test_energy_density_positive_for_spheres_that_nearly_touch(self):
        """Aichelburg and Schein's 1998 result, at m = q = 0 and R/d_1 = 0.95."""
        mp = self.mp
        for th in (0, mp.pi / 2, mp.pi):
            self.assertGreater(self.sigma(0.95, 1, 0.5, 0, 0, th), 0, th)

    def test_shells_sit_at_the_areal_radius_R_plus_m1(self):
        c = chart("interior")
        self.assertIn("r = R + m_1 \\;\\text{(the shells)}", c["domains"])
        self.assertGreater(R + M1, M + (M ** 2 - Q ** 2) ** 0.5)


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
        paragraphs = METRIC["history"].split("¶")
        self.assertEqual(len(paragraphs), 6)
        self.assertTrue(METRIC["history"].startswith("At the Institute for Theoretical Physics"))
        self.assertLessEqual(len(METRIC["history"].replace("¶", "")), 6900)

    def test_founders_named_in_full_first(self):
        history = METRIC["history"]
        self.assertLess(history.index("Friedrich Schein"), history.index("Schein and Aichelburg"))
        self.assertLess(history.index("Peter Christian Aichelburg"), history.index("Aichelburg and Schein"))

    def test_related_spacetimes_answer(self):
        kinds = {"piece": "composite", "family": "family"}
        for entry in METRIC["related"]:
            other = json.loads((DATA / "metrics" / f"{entry['id']}.json").read_text(encoding="utf-8"))
            answer = [e for e in other["related"] if e["id"] == ME]
            self.assertEqual(len(answer), 1, entry["id"])
            self.assertEqual(answer[0]["kind"], kinds[entry["kind"]], entry["id"])
            for key in cited(answer[0]["text"]):
                self.assertIn(key, other["references"], (entry["id"], key))

    def test_every_chart_has_a_spacetime_diagram(self):
        drawn = set(drawing("diagrams")["systems"])
        self.assertEqual(drawn, set(SYSTEMS))


class Drawings(unittest.TestCase):
    """What the drawings mark."""

    def test_the_conformal_tower_marks_both_shells_and_what_they_cut_away(self):
        view, = drawing("conformal")["views"]
        self.assertEqual((view["id"], view["system"]), ("tower", "interior"))
        text = json.dumps(view)
        self.assertIn("\"surface\"", text)
        self.assertIn("\"removed\"", text)

    def test_the_embedding_is_the_midplane(self):
        views = drawing("embedding")["views"]
        self.assertEqual([v["id"] for v in views], ["midplane"])

    def test_the_curvature_figure_draws_inside_and_on_the_midplane(self):
        views = drawing("curvature")["views"]
        self.assertEqual([(v["id"], v["system"]) for v in views], [("interior", "interior"), ("midplane", "cartesian")])


if __name__ == "__main__":
    unittest.main()
