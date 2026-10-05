"""The quantum BTZ black hole of Emparan, Frassino and Way (arXiv:2007.15999):
.venv.noindex/bin/python -m unittest discover -s _tools

What its texts and drawings state, held on the published files. The first class needs nothing
but the files: the horizons where every drawing puts them and the throat of the embedded moment.
The second reads the published metric and Einstein tensor through the checker's Reader and holds
them to the paper: the stress tensor of the quantum fields, the masses and stress strengths of the
family of solutions, the rescaling between the brane chart and the static one, and the first law
for the generalised entropy; it needs sympy and is skipped where it is absent.
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


def load(folder):
    return json.loads((DATA / folder / "quantum_btz.json").read_text(encoding="utf-8"))


class Drawings(unittest.TestCase):
    def test_the_spacetime_diagrams_mark_the_horizons_of_each_example(self):
        """The hole of the static and Eddington-Finkelstein charts, M = F = 1/4 and l = 15 l_3/16, has
        H = (r - 3/4)(r^2 + 3r/4 + 5/16)/r and one horizon at 3 l_3/4; the dressed cone, kappa = +1,
        mu = 6 and l = l_3/3, has H = (r - 1)(r^2 + r + 2)/r and its horizon at l_3; the rotating hole
        has r^2 H = (r - 1)(r - 1/2)(r^2 + 3r/2 + 3/4) and horizons at l_3 and l_3/2."""
        want = {"static": [0.75], "eddington_finkelstein_ingoing": [0.75], "eddington_finkelstein_outgoing": [0.75],
                "brane": [1.0], "rotating": [0.5, 1.0]}
        systems = load("diagrams")["systems"]
        self.assertEqual(sorted(systems), sorted(want))
        for system, horizons in want.items():
            for view in systems[system]:
                X0, X1 = view["box"][:2]
                marked = sorted(X0 + line[0][0] * (X1 - X0)
                                for m in view["markers"] if m["kind"] == "grr" for line in m["lines"])
                self.assertEqual(len(marked), len(horizons), f"{system}/{view['id']}")
                for got, at in zip(marked, horizons):
                    self.assertAlmostEqual(got, at, delta=2e-3, msg=f"{system}/{view['id']}")

    def test_the_embedded_moment_has_its_throat_on_the_horizon(self):
        """The moment t = 0 runs from the bifurcation circle r_+ = 3 l_3/4, where the radius of the
        surface is least, out to r = 2.5 l_3 on both sheets."""
        view = next(v for v in load("embedding")["views"] if v["id"] == "throat")
        radii = [p[1] for piece in view["surfaces"][0]["pieces"] for p in piece["points"]]
        self.assertAlmostEqual(min(radii), 0.75, places=6)
        self.assertAlmostEqual(max(radii), 2.5, places=6)

    def test_every_chart_has_a_conformal_diagram(self):
        views = {v["id"]: v["system"] for v in load("conformal")["views"]}
        self.assertEqual(views, {"static": "static", "ingoing": "eddington_finkelstein_ingoing",
                                 "outgoing": "eddington_finkelstein_outgoing", "dressed": "brane",
                                 "rotating": "rotating"})


@unittest.skipUnless(HAS_SYMPY, "sympy is not installed")
class Physics(unittest.TestCase):
    def setUp(self):
        import sympy
        import verify_metrics
        self.sp, self.vm = sympy, verify_metrics
        metric = json.loads((DATA / "metrics" / "quantum_btz.json").read_text(encoding="utf-8"))
        self.charts = {c["id"]: c for c in metric["coordinates"]}

    def read(self, system, field, variant=None):
        entry = self.charts[system]
        reader = self.vm.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]], ())
        block = entry[field] if variant is None else entry[field]["variants"][variant]["nonzero"]
        return reader, {tuple(e["indices"]): reader(e["value"]) for e in block}

    def test_the_quantum_fields_carry_the_holographic_stress_tensor(self):
        """8 pi G_3 <T^a_b> = G^a_b - delta^a_b/l_3^2 is (l F/2r^3) diag(1, 1, -2), traceless, falling
        as 1/r^3, the paper's (2.46), and in the brane chart mu takes F's place."""
        sp = self.sp
        for system, strength in (("static", "F"), ("brane", "mu")):
            reader, G = self.read(system, "einstein_tensor", "ul")
            r, ell3, ell = reader.symbol["r"], reader.parameters["ell_3"], reader.parameters["ell"]
            q = ell * reader.parameters[strength] / (2 * r ** 3)
            for index, sign in ((("t", "t"), 1), (("r", "r"), 1), (("\\phi", "\\phi"), -2)):
                self.assertEqual(sp.simplify(G[index] - 1 / ell3 ** 2 - sign * q), 0, f"{system} {index}")

    def test_the_family_runs_from_anti_de_sitter_space_to_a_largest_mass(self):
        """With k = -kappa x_1^2 the brane chart's M = -kappa Delta^2 and F = mu Delta^3, Delta = 2x_1/(3 -
        kappa x_1^2) and mu = (1 - kappa x_1^2)/x_1^3, are the parameter descriptions' M = 4k/(3 + k)^2
        and F = 8(1 + k)/(3 + k)^3. At k = -1, M = -1 and F = 0, anti-de Sitter space; the largest mass
        is M = 1/3 at k = 3, where the two branches of positive mass meet."""
        sp = self.sp
        x, kappa, k = sp.symbols("x kappa k", real=True)
        Delta = 2 * x / (3 - kappa * x ** 2)
        mu = (1 - kappa * x ** 2) / x ** 3
        on = {kappa: -k / x ** 2}
        self.assertEqual(sp.simplify((-kappa * Delta ** 2).subs(on) - 4 * k / (3 + k) ** 2), 0)
        self.assertEqual(sp.simplify((mu * Delta ** 3).subs(on) - 8 * (1 + k) / (3 + k) ** 3), 0)
        M = 4 * k / (3 + k) ** 2
        self.assertEqual(M.subs(k, -1), -1)
        self.assertEqual((8 * (1 + k) / (3 + k) ** 3).subs(k, -1), 0)
        self.assertEqual(sp.solve(sp.diff(M, k), k), [3])
        self.assertEqual(M.subs(k, 3), sp.Rational(1, 3))
        # The drawn hole, k = 1, has M = F = 1/4 on the hotter branch 0 < k < 3.
        self.assertEqual((M.subs(k, 1), (8 * (1 + k) / (3 + k) ** 3).subs(k, 1)), (sp.Rational(1, 4), sp.Rational(1, 4)))

    def test_the_brane_chart_rescaled_is_the_static_chart(self):
        """t = Delta t-bar, r = r-bar/Delta and phi = Delta phi-bar, the paper's (2.36), take the brane
        chart's H = r^2/l_3^2 + kappa - mu l/r to the static chart's with M = -kappa Delta^2 and
        F = mu Delta^3, so that phi-bar has period 2 pi."""
        sp = self.sp
        reader, g = self.read("brane", "metric_components")
        _, h = self.read("static", "metric_components")
        static, _ = self.read("static", "inverse_metric_components")
        r = reader.symbol["r"]
        Delta = sp.Symbol("Delta", positive=True)
        P, Q = reader.parameters, static.parameters
        at = {Q["M"]: -P["kappa"] * Delta ** 2, Q["F"]: P["mu"] * Delta ** 3, static.symbol["r"]: Delta * r}
        self.assertEqual(sp.simplify(h[("t", "t")].subs(at) / Delta ** 2 - g[("t", "t")]), 0)
        self.assertEqual(sp.simplify(h[("r", "r")].subs(at) * Delta ** 2 - g[("r", "r")]), 0)
        self.assertEqual(sp.simplify(h[("\\phi", "\\phi")].subs(at) / Delta ** 2 - g[("\\phi", "\\phi")]), 0)

    def test_the_temperature_on_the_brane_is_the_papers(self):
        """With nu = l/l_3 and z = l_3/(r_+ x_1), the paper's (2.60) to (2.62) give x_1, r_+ and mu x_1
        for kappa = -1, and Delta H'(r_+)/4 pi on the published brane chart is then its (2.68),
        z(2 + 3 nu z + nu z^3)/(2 pi l_3 (1 + 3z^2 + 2 nu z^3)), at several points of the family."""
        sp = self.sp
        reader, ginv = self.read("brane", "inverse_metric_components")
        r, P = reader.symbol["r"], reader.parameters
        H = ginv[("r", "r")]
        for nu, z in ((sp.Rational(1, 3), sp.Rational(1, 2)), (sp.Rational(1, 10), 1), (sp.Rational(1, 2), sp.Rational(1, 3))):
            x1 = sp.sqrt((1 - nu * z ** 3) / (z ** 2 * (1 + nu * z)))
            rp = sp.sqrt((1 + nu * z) / (1 - nu * z ** 3))
            mu = (1 + z ** 2) / (1 - nu * z ** 3) / x1
            Delta = 2 * x1 / (3 + x1 ** 2)
            at = {P["ell_3"]: 1, P["kappa"]: -1, P["mu"]: mu, P["ell"]: nu}
            self.assertEqual(sp.simplify(H.subs(at).subs(r, rp)), 0)
            T = Delta * sp.diff(H, r).subs(at).subs(r, rp) / (4 * sp.pi)
            want = z * (2 + 3 * nu * z + nu * z ** 3) / (2 * sp.pi * (1 + 3 * z ** 2 + 2 * nu * z ** 3))
            self.assertEqual(sp.simplify(T - want), 0, (nu, z))

    def test_the_generalised_entropy_obeys_the_first_law(self):
        """The paper's (2.63), (2.68) and (2.70): M, T and S_gen as functions of z at fixed nu, with
        curly G_3 = G_3/sqrt(1 + nu^2), satisfy dM = T dS_gen for every nu, and the area of the
        horizon on the brane, S_cl = (1 + nu z) S_gen/sqrt(1 + nu^2), its (2.71), does not."""
        sp = self.sp
        z, nu, G3, l3 = sp.symbols("z nu G_3 ell_3", positive=True)
        cal = G3 / sp.sqrt(1 + nu ** 2)
        D = 1 + 3 * z ** 2 + 2 * nu * z ** 3
        M = z ** 2 * (1 - nu * z ** 3) * (1 + nu * z) / (2 * cal * D ** 2)
        T = z * (2 + 3 * nu * z + nu * z ** 3) / (2 * sp.pi * l3 * D)
        S = sp.pi * l3 * z * sp.sqrt(1 + nu ** 2) / (G3 * D)
        self.assertEqual(sp.simplify(sp.diff(M, z) - T * sp.diff(S, z)), 0)
        S_cl = (1 + nu * z) * S / sp.sqrt(1 + nu ** 2)
        self.assertNotEqual(sp.simplify(sp.diff(M, z) - T * sp.diff(S_cl, z)), 0)
        # At nu = 0 both are the BTZ entropy, pi l_3 sqrt(2M/G_3) in units where c = hbar = 1.
        self.assertEqual(sp.simplify((S - sp.pi * l3 * sp.sqrt(2 * M / G3)).subs(nu, 0)), 0)


if __name__ == "__main__":
    unittest.main()
