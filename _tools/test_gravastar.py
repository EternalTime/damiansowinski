"""The gravastar, a ball of de Sitter space under a thin shell with Schwarzschild's vacuum outside:
python3 -m unittest discover -s _tools

What its texts state about the shell, held on the published files. The shell is in no chart, since
g_rr jumps across it, so its surface density and tension are taken here from the jump of the
extrinsic curvature between the published interior and exterior metrics, by the junction conditions
of Israel as Visser and Wiltshire wrote them for a static shell, their (12) and (13), and compared
with the closed forms the History, the embedding diagram's caption and gravastar.md state. The
tests need sympy and are skipped where it is absent.
"""
import importlib.util
import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "MFS" / "assets" / "data"
HAS_SYMPY = importlib.util.find_spec("sympy") is not None
sys.path.insert(0, str(Path(__file__).resolve().parent / "derivations"))


@unittest.skipUnless(HAS_SYMPY, "sympy is not installed")
class Shell(unittest.TestCase):
    def setUp(self):
        import sympy
        import verify_metrics
        self.sp, self.vm = sympy, verify_metrics
        self.rs, self.R, self.L, self.r = sympy.symbols("r_s R L r", positive=True)

    def published(self, metric_id, system):
        """The published g_tt and g_rr of a chart of t and r, in the symbols of this test, G = c = 1."""
        metric = json.loads((DATA / "metrics" / f"{metric_id}.json").read_text(encoding="utf-8"))
        entry = next(c for c in metric["coordinates"] if c["id"] == system)
        reader = self.vm.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]], ())
        values = {tuple(e["indices"]): reader(e["value"]) for e in entry["metric_components"]}
        mine = {reader.symbol["r"]: self.r, **{reader.parameters[n]: getattr(self, a)
                                               for n, a in (("r_s", "rs"), ("R", "R"), ("L", "L"))
                                               if n in reader.parameters}}
        return tuple(self.sp.sympify(values[k]).subs(mine) for k in (("t", "t"), ("r", "r")))

    def curvatures(self, metric_id, system):
        """K^t_t and K^theta_theta of the sphere r = R, with the normal pointing outward:
        sqrt(h) f'/(2 f) and sqrt(h)/r, where f = -g_tt and h = 1/g_rr."""
        gtt, grr = self.published(metric_id, system)
        f, h = -gtt, 1 / grr
        at = lambda e: e.subs(self.r, self.R)
        return at(self.sp.sqrt(h) * self.sp.diff(f, self.r) / (2 * f)), at(self.sp.sqrt(h) / self.r)

    def shell(self):
        """The surface density and the surface tension, from [[K^theta_theta]] = -4 pi sigma and
        [[K^t_t + K^theta_theta]] = -8 pi tension."""
        (kt_in, kth_in), (kt_out, kth_out) = self.curvatures("gravastar", "interior"), self.curvatures("gravastar", "exterior")
        sigma = -(kth_out - kth_in) / (4 * self.sp.pi)
        tension = -((kt_out + kth_out) - (kt_in + kth_in)) / (8 * self.sp.pi)
        return sigma, tension

    def zero(self, expression, values):
        self.assertLess(abs(complex(expression.subs(values).evalf(30))), 1e-20)

    POINTS = ({"r_s": 1, "R": "5/4", "L": 2}, {"r_s": 1, "R": "21/20", "L": "3/2"}, {"r_s": 2, "R": 3, "L": 7})

    def at(self, point):
        return {self.rs: self.sp.Rational(point["r_s"]), self.R: self.sp.Rational(point["R"]),
                self.L: self.sp.Rational(point["L"])}

    def test_the_time_is_one_on_both_faces_of_the_shell(self):
        inside, outside = self.published("gravastar", "interior")[0], self.published("gravastar", "exterior")[0]
        self.assertEqual(self.sp.simplify((inside - outside).subs(self.r, self.R)), 0)

    def test_the_surface_density_is_the_jump_of_the_square_root_of_h(self):
        # sigma = (sqrt(1 - R^2/L^2) - sqrt(1 - r_s/R))/(4 pi R), Visser and Wiltshire's (12).
        sigma, _ = self.shell()
        R, L, rs, sp = self.R, self.L, self.rs, self.sp
        stated = (sp.sqrt(1 - R ** 2 / L ** 2) - sp.sqrt(1 - rs / R)) / (4 * sp.pi * R)
        for point in self.POINTS:
            self.zero(sigma - stated, self.at(point))

    def test_the_surface_tension_is_visser_and_wiltshires(self):
        # [[(1 - m/a - m')/(a sqrt(1 - 2m/a))]] = -8 pi tension, their (13), with m = r_s/2 outside
        # and m = r^3/(2 L^2) inside.
        _, tension = self.shell()
        R, L, rs, sp = self.R, self.L, self.rs, self.sp
        outside = (1 - rs / (2 * R)) / (R * sp.sqrt(1 - rs / R))
        inside = (1 - 2 * R ** 2 / L ** 2) / (R * sp.sqrt(1 - R ** 2 / L ** 2))
        for point in self.POINTS:
            self.zero(tension + (outside - inside) / (8 * sp.pi), self.at(point))

    def test_the_shell_has_no_mass_where_the_vacuum_accounts_for_all_of_it(self):
        # At L^2 = R^3/r_s the surface density vanishes, C = 1 and g_rr is continuous, the thin shell
        # gravastar of Pani and his collaborators; the shell is then under pressure, a negative tension.
        sigma, tension = self.shell()
        metric = json.loads((DATA / "metrics" / "gravastar.json").read_text(encoding="utf-8"))
        entry = next(c for c in metric["coordinates"] if c["id"] == "interior")
        reader = self.vm.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]], ())
        C = reader.parameters["C"].subs({reader.parameters["r_s"]: self.rs, reader.parameters["R"]: self.R,
                                         reader.parameters["L"]: self.L})
        for rs, R in ((1, self.sp.Rational(5, 4)), (1, self.sp.Rational(21, 20)), (2, 3)):
            values = {self.rs: rs, self.R: R, self.L: self.sp.sqrt(R ** 3 / self.sp.Integer(rs))}
            self.zero(sigma, values)
            self.assertEqual(self.sp.simplify(C.subs(values)), 1)
            jump = self.published("gravastar", "interior")[1] - self.published("gravastar", "exterior")[1]
            self.zero(jump.subs(self.r, self.R), values)
            self.assertLess(float(tension.subs(values)), 0)

    def test_the_drawn_gravastar_has_a_shell_of_positive_density(self):
        sigma, _ = self.shell()
        self.assertGreater(float(sigma.subs(self.at(self.POINTS[0]))), 0)

    def test_schwarzschilds_star_at_its_own_radius_is_the_interior_with_c_a_quarter(self):
        # Mazur and Mottola's limit: at R = r_s the star of uniform density has
        # -g_tt = (1 - r^2/r_s^2)/4 and 1/g_rr = 1 - r^2/r_s^2, the gravastar's interior with L = r_s
        # and C = 1/4.
        gtt, grr = self.published("interior_schwarzschild", "spherical")
        star = {self.R: self.rs}
        mine_tt, mine_rr = self.published("gravastar", "interior")
        h = 1 - self.r ** 2 / self.rs ** 2
        self.assertEqual(self.sp.simplify(gtt.subs(star) + h / 4), 0)
        self.assertEqual(self.sp.simplify(1 / grr.subs(star) - h), 0)
        self.assertEqual(self.sp.simplify(1 / mine_rr.subs(self.L, self.rs) - h), 0)


if __name__ == "__main__":
    unittest.main()
