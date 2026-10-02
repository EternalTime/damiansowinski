"""The parametrised post-Newtonian metric, the weak field of a body with Eddington's numbers beta
and gamma in front of its terms: python3 -m unittest discover -s _tools

What its texts state, held on the published files: the bending of light by (1 + gamma)/2 times
Einstein's, the advance of the perihelion by (2 + 2 gamma - beta)/3 times his, Shapiro's logarithm,
the dragging of frames, the vacuum of general relativity at beta = gamma = 1 and nowhere else, and
the numbers the captions of the drawings give. It also holds the rule by which the checker keeps a
post-Newtonian component, PostNewtonian in verify_metrics.py: terms of the metric one order beyond
the post-Newtonian metric move no component that is kept, and do move one cut an order later. The
tests need sympy and are skipped where it is absent.
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


def flat(tensor):
    if isinstance(tensor, list):
        for item in tensor:
            yield from flat(item)
    else:
        yield tensor


@unittest.skipUnless(HAS_SYMPY, "sympy is not installed")
class Published(unittest.TestCase):
    """The published charts, read through the checker's Reader, with c = 1."""

    @classmethod
    def setUpClass(cls):
        import sympy
        import verify_metrics
        cls.sp, cls.vm = sympy, verify_metrics
        cls.metric = json.loads((DATA / "metrics" / "ppn_metric.json").read_text(encoding="utf-8"))

    def chart(self, system):
        entry = next(c for c in self.metric["coordinates"] if c["id"] == system)
        reader = self.vm.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]], ())
        return entry, reader

    def components(self, system, field="metric_components"):
        entry, reader = self.chart(system)
        return reader, {tuple(e["indices"]): self.sp.sympify(reader(e["value"])) for e in entry[field]}

    def lapse_and_space(self, system):
        """A = -g_tt, B = g_rr and C = g_phiphi on the equator, with the symbols r, m, beta, gamma."""
        reader, g = self.components(system)
        on_equator = {reader.symbol["\\theta"]: self.sp.pi / 2}
        r = reader.symbol["r"]
        m, beta, gamma = (reader.parameters[k] for k in ("m", "beta", "gamma"))
        return (-g[("t", "t")], g[("r", "r")], g[("\\phi", "\\phi")].subs(on_equator)), (r, m, beta, gamma)

    def test_light_grazing_the_body_is_bent_by_half_of_one_plus_gamma_times_einsteins_angle(self):
        """The path of light in a static isotropic metric is that of a ray in a medium of index
        n = sqrt(g_rr/(-g_tt)), which turns through 2 * integral from b of dr/(r sqrt((r n/(b n_b))^2 - 1))
        less pi. At m/b = 10^-4 that is (1 + gamma)/2 times 4m/b to the next order in m/b."""
        import mpmath
        (A, B, _), (r, m, beta, gamma) = self.lapse_and_space("isotropic")
        index = self.sp.sqrt(B / A)
        for value in (1, self.sp.Rational(1, 2), 0):
            n = self.sp.lambdify(r, index.subs({m: 1, beta: 1, gamma: value}), "mpmath")
            b = mpmath.mpf(10) ** 4
            turned = 2 * mpmath.quad(lambda x: 1 / (x * mpmath.sqrt((x * n(x) / (b * n(b))) ** 2 - 1)),
                                     [b, 2 * b, 100 * b, mpmath.inf]) - mpmath.pi
            einstein = 4 / b
            self.assertAlmostEqual(float(turned / einstein), (1 + float(value)) / 2, places=3)

    def test_the_perihelion_advances_by_a_third_of_two_plus_two_gamma_minus_beta_times_einsteins(self):
        """A nearly circular orbit of -A dt^2 + B dr^2 + C dphi^2 closes after its radial oscillation,
        of frequency kappa, and the perihelion advances by 2 pi (1 - kappa/Omega) a turn. To first
        order in m/r that is 6 pi m/r times (2 + 2 gamma - beta)/3, in the isotropic chart and in
        the areal one."""
        sp = self.sp
        for system in ("isotropic", "areal"):
            (A, B, C), (r, m, beta, gamma) = self.lapse_and_space(system)
            E2, L2 = sp.symbols("E2 L2")
            radial = (E2 / A - L2 / C - 1) / B                  # (dr/dtau)^2
            circle = sp.solve([radial, sp.diff(radial, r)], [E2, L2], dict=True)[0]
            kappa2 = (-sp.diff(radial, r, 2) / 2).subs(circle)
            omega2 = (L2 / C ** 2).subs(circle)
            eps = sp.Symbol("eps", positive=True)
            # kappa/Omega is 1 - (1 - kappa^2/Omega^2)/2 to first order.
            ratio = sp.series(sp.simplify(kappa2 / omega2).subs(m, eps * m), eps, 0, 2).removeO()
            advance = sp.simplify((1 - ratio).subs(eps, 1) / 2 * r / (3 * m))
            self.assertEqual(sp.simplify(advance - (2 + 2 * gamma - beta) / 3), 0, system)

    def test_a_ray_is_delayed_by_one_plus_gamma_times_the_logarithm(self):
        """dr_*/dr = sqrt(g_rr/(-g_tt)) is 1 + (1 + gamma) m/r to first order, whose integral is
        Shapiro's logarithm, and the round trip the caption gives from 26 m to the surface at 10 m is 3.9 m."""
        import mpmath
        sp = self.sp
        (A, B, _), (r, m, beta, gamma) = self.lapse_and_space("isotropic")
        eps = sp.Symbol("eps", positive=True)
        speed = sp.series(sp.sqrt(B / A).subs(m, eps * m), eps, 0, 2).removeO().subs(eps, 1)
        self.assertEqual(sp.simplify(speed - 1 - (1 + gamma) * m / r), 0)
        drawn = sp.lambdify(r, sp.sqrt(B / A).subs({m: 1, beta: 1, gamma: 1}), "mpmath")
        self.assertEqual(round(float(2 * (mpmath.quad(drawn, [10, 26]) - 16)), 1), 3.9)

    def test_the_numbers_of_the_captions(self):
        sp = self.sp
        (A, B, _), (r, m, beta, gamma) = self.lapse_and_space("isotropic")
        edge = sp.lambdify(r, sp.sqrt(A / B).subs({m: 1, beta: 1, gamma: 1}))
        self.assertEqual(round(edge(10.0), 3), 0.827)
        self.assertEqual(round(edge(26.0), 3), 0.927)
        (A, B, _), (r, m, beta, gamma) = self.lapse_and_space("areal")
        edge = sp.lambdify(r, sp.sqrt(A / B).subs({m: 1, beta: 1, gamma: 1}))
        self.assertEqual(round(edge(11.0), 3), 0.832)
        self.assertEqual(round(1 - 2 / 11, 3), 0.818)
        # At beta = gamma the term in m^2 leaves the areal chart's g_tt.
        self.assertEqual(sp.simplify(A.subs(beta, gamma) - (1 - 2 * m / r)), 0)

    def test_the_rotating_body_drags_frames_at_delta_times_lense_and_thirrings_rate(self):
        """d phi/d(ct) = -g_tphi/g_phiphi is 2 Delta a m/(r^3 (1 + 2 gamma m/r)), with
        Delta = (1 + gamma + alpha_1/4)/2, which is Lense and Thirring's published 2am/r^3 to first
        order at general relativity's values; at the surface drawn the speed around is 0.037 c."""
        sp = self.sp
        reader, g = self.components("rotating")
        r, th = reader.symbol["r"], reader.symbol["\\theta"]
        m, a, gamma, alpha = (reader.parameters[k] for k in ("m", "a", "gamma", "alpha_1"))
        rate = sp.simplify(-g[("t", "\\phi")] / g[("\\phi", "\\phi")])
        delta = (1 + gamma + alpha / 4) / 2
        self.assertEqual(sp.simplify(rate - 2 * delta * a * m / (r ** 3 * (1 + 2 * gamma * m / r))), 0)
        theirs = json.loads((DATA / "metrics" / "hartle_thorne.json").read_text(encoding="utf-8"))
        entry = next(c for c in theirs["coordinates"] if c["id"] == "lense_thirring")
        weak = self.vm.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]], ())
        lt = {tuple(e["indices"]): sp.sympify(weak(e["value"])) for e in entry["metric_components"]}
        same = {weak.symbol["r"]: r, weak.symbol["\\theta"]: th, weak.parameters["m"]: m, weak.parameters["a"]: a}
        self.assertEqual(sp.simplify(lt[("t", "\\phi")].subs(same) - g[("t", "\\phi")].subs({gamma: 1, alpha: 0})), 0)
        drawn = {m: 1, a: 2, gamma: 1, alpha: 0, r: 10, th: sp.pi / 2}
        around = sp.sqrt(g[("\\phi", "\\phi")]) * rate
        self.assertEqual(round(float(around.subs(drawn)), 3), 0.037)

    def test_the_ricci_tensor_vanishes_at_beta_and_gamma_one_and_nowhere_else(self):
        """R_tt = -(1 - 2 beta + gamma) m^2/r^4 and R_rr = 2 (1 - gamma) m/r^3 in the isotropic chart:
        both vanish only at beta = gamma = 1, general relativity's point of the plane."""
        sp = self.sp
        reader, ricci = None, None
        entry, reader = self.chart("isotropic")
        ricci = {tuple(e["indices"]): sp.sympify(reader(e["value"]))
                 for e in entry["ricci_tensor"]["variants"]["ll"]["nonzero"]}
        r = reader.symbol["r"]
        m, beta, gamma = (reader.parameters[k] for k in ("m", "beta", "gamma"))
        self.assertEqual(sp.simplify(ricci[("t", "t")] + (1 - 2 * beta + gamma) * m ** 2 / r ** 4), 0)
        self.assertEqual(sp.simplify(ricci[("r", "r")] - 2 * (1 - gamma) * m / r ** 3), 0)
        self.assertEqual(sp.solve([ricci[("t", "t")], ricci[("r", "r")]], [beta, gamma], dict=True),
                         [{beta: 1, gamma: 1}])
        self.assertEqual(sp.simplify(reader(entry["kretschmann"].split("=", 1)[1])
                                     - 24 * (1 + gamma ** 2) * m ** 2 / r ** 6), 0)

    def test_the_embedded_equator_rises_as_gamma_alone_sets(self):
        """A moment of t on the equator is B (dr^2 + r^2 dphi^2), whose circle has the radius
        sqrt(r (r + 2 gamma m)) and which rises as dz/dr = sqrt(gamma m (2r + 3 gamma m)/(r (r + 2 gamma m)));
        by 3R it stands 0.30 m below Flamm's paraboloid over the same circles."""
        import mpmath
        sp = self.sp
        (_, B, C), (r, m, beta, gamma) = self.lapse_and_space("isotropic")
        self.assertEqual(sp.simplify(C - r * (r + 2 * gamma * m)), 0)
        rise2 = sp.simplify(B - sp.diff(sp.sqrt(C), r) ** 2)
        self.assertEqual(sp.simplify(rise2 - gamma * m * (2 * r + 3 * gamma * m) / (r * (r + 2 * gamma * m))), 0)
        self.assertFalse(rise2.has(beta))
        height = mpmath.quad(lambda x: mpmath.sqrt((2 * x + 3) / (x * (x + 2))), [10, 30])
        flamm = lambda rho: 2 * mpmath.sqrt(2 * (rho - 2))
        short = flamm(mpmath.sqrt(30 * 32)) - flamm(mpmath.sqrt(10 * 12)) - height
        self.assertEqual(round(float(short), 2), 0.30)


@unittest.skipUnless(HAS_SYMPY, "sympy is not installed")
class Orders(unittest.TestCase):
    """The post-Newtonian orders the checker keeps."""

    @classmethod
    def setUpClass(cls):
        import sympy
        import verify_metrics
        cls.sp, cls.vm = sympy, verify_metrics

    def test_each_component_is_kept_as_far_as_its_time_indices_say(self):
        reader = self.vm.Reader(["t", "r", "\\theta", "\\phi"], ["m", "a"], (), kept=({"m": 2, "a": 1}, 2, "t"))
        for field, index, order in (("christoffel", (1, 0, 0), 4), ("christoffel", (0, 0, 1), 4),
                                    ("christoffel", (1, 0, 3), 3), ("christoffel", (1, 2, 2), 2),
                                    ("riemann", (0, 1, 0, 1), 4), ("riemann", (1, 2, 1, 2), 2),
                                    ("ricci_tensor", (0, 0), 4), ("ricci_tensor", (1, 1), 2),
                                    ("einstein_tensor", (0, 0), 2), ("einstein_tensor", (0, 3), 3),
                                    ("weyl_tensor", (0, 1, 0, 1), 2), ("weyl_tensor", (0, 1, 1, 3), 3),
                                    ("ricci_scalar", (), 2), ("kretschmann", (), 4)):
            self.assertEqual(reader.kept_order(field, index), order, (field, index))
        m, r = reader.parameters["m"], reader.symbol["r"]
        self.assertEqual(self.vm.norm(reader.through(1 / (1 - m / r), 4) - (1 + m / r + m ** 2 / r ** 2)), 0)
        self.assertEqual(self.vm.norm(reader.through(1 / (1 - m / r), 2) - (1 + m / r)), 0)
        with self.assertRaises(self.vm.LatexError):
            self.vm.Reader(["t", "r"], ["m"], (), kept=({"m": 2}, 2, "u"))

    def test_terms_of_the_next_order_in_the_metric_move_nothing_that_is_kept(self):
        """The post-Newtonian metric leaves g_tt open at the sixth order of v/c, g_ti at the fifth and
        g_ij at the fourth. With an arbitrary function of r and theta put in each of the ten components
        at that order, every tensor the checker hands over is what it was, so nothing kept depends
        on a term the metric does not have; the Einstein tensor's tt component cut two orders
        later, as the Ricci tensor's is, does move."""
        import print_charts
        sp, vm = self.sp, self.vm
        spec = print_charts.ppn_metric("isotropic")
        coords = spec["system"]["coords"]

        def geometry(open_terms):
            reader = vm.Reader(coords, spec["system"]["parameters"] + ["p", "q", "s"], (),
                               kept=({"m": 2, "p": 6, "q": 5, "s": 4}, 2, "t"))
            g = vm.metric_from_line_element(reader, spec["chart_line_element"], coords)
            symbols = [reader.symbol[name] for name in coords]
            if open_terms:
                strength = [reader.parameters[k] for k in "sqp"]        # by the number of time indices
                for i in range(4):
                    for j in range(i, 4):
                        g[i, j] += strength[(i == 0) + (j == 0)] * sp.Function(f"F{i}{j}")(symbols[1], symbols[2])
                        g[j, i] = g[i, j]
            return vm.geometry_of(g, symbols, 10 ** 6, reader)

        fixed, opened = geometry(False), geometry(True)
        for name in ("christoffel_ull", "christoffel_lll", "riemann_llll", "ricci_ll", "einstein_ll", "weyl_llll",
                     "ricci_scalar"):
            ours, theirs = getattr(fixed, name)(), getattr(opened, name)()
            self.assertTrue(all(vm.norm(a - b) == 0 for a, b in zip(flat(ours), flat(theirs))), name)
            if name in ("riemann_llll", "ricci_ll", "einstein_ll", "weyl_llll"):
                rank = 4 if name.endswith("llll") else 2
                raised = [g.raise_indices(t, rank, (0,)) for g, t in ((fixed, ours), (opened, theirs))]
                self.assertTrue(all(vm.norm(a - b) == 0 for a, b in zip(*map(flat, raised))), name + " raised")
        late = [g._geometry.einstein_ll()[0][0] for g in (fixed, opened)]
        self.assertNotEqual(vm.norm(late[0] - late[1]), 0)


if __name__ == "__main__":
    unittest.main()
