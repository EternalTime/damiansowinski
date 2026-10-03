"""Vuorio's universe and warped anti-de Sitter space, the homogeneous vacua of topologically
massive gravity with a timelike Killing vector of constant twist:
python3 -m unittest discover -s _tools

What its texts state, held on the published files: every chart is -(c dt + A)^2 + h with dA twice
the twist Omega times the area of h; the Ricci scalar is 2(Omega^2 - m^2) and the Kretschmann
scalar 44 Omega^4 - 24 Omega^2 m^2 + 4 m^4, anti-de Sitter space's at m = 2 Omega; the field
equation G + Lambda g + C/mu = 0 holds with mu = 3 Omega and Lambda = (Omega^2 - m^2)/3, so that
Vuorio's member m = Omega has no cosmological constant; that member is Vuorio's (2.21) as Chow,
Pope and Sezgin quote it; the circle of constant t and r is null where tanh(mr/2) = m/(2 Omega)
and timelike beyond; and the drawings hold what their captions say. The tests of the published
mathematics need sympy and are skipped where it is absent; the tests of the drawings read the files
alone.
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
CHARTS = ["cylindrical", "disc", "fibred", "horospherical"]
# Vuorio's member, m = Omega = 1, at which every drawing is made, and its null circle.
RC = math.log(3)


def published():
    return json.loads((DATA / "metrics" / "vuorio_warped_ads.json").read_text(encoding="utf-8"))


def drawing(kind):
    return json.loads((DATA / kind / "vuorio_warped_ads.json").read_text(encoding="utf-8"))


def twist_and_spread(r):
    """H and D of the cylindrical chart at m = Omega = 1: -(c dt + H dphi)^2 + dr^2 + D^2 dphi^2."""
    return 4 * math.sinh(r / 2) ** 2, math.sinh(r)


@unittest.skipUnless(HAS_SYMPY, "sympy is not installed")
class Published(unittest.TestCase):
    def setUp(self):
        import sympy
        import verify_metrics
        self.sp, self.vm = sympy, verify_metrics

    def chart(self, system):
        """The published metric of a chart as a sympy matrix with c = 1, its coordinates, and its
        parameters Omega and m."""
        entry = next(c for c in published()["coordinates"] if c["id"] == system)
        reader = self.vm.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]], ())
        names = entry["coords"]
        g = self.sp.zeros(3, 3)
        for e in entry["metric_components"]:
            i, j = (names.index(k) for k in e["indices"])
            g[i, j] = reader(e["value"]).subs(reader.c, 1)
        X = [reader.symbol[n] for n in names]
        return g, X, reader.parameters["Omega"], reader.parameters["m"], entry, reader

    def scalar(self, entry, reader, field):
        return reader(entry[field].partition("=")[2]).subs(reader.c, 1)

    def test_every_chart_is_a_twisted_line_over_a_plane_of_curvature_minus_m_squared(self):
        sp = self.sp
        for system in CHARTS:
            g, X, Om, m, _, _ = self.chart(system)
            self.assertEqual(sp.simplify(g[0, 0] + 1), 0, system)
            A = [-g[0, 1], -g[0, 2]]
            h = sp.Matrix(2, 2, lambda i, j: g[i + 1, j + 1] + A[i] * A[j])
            self.assertTrue(all(sp.diff(v, X[0]) == 0 for v in list(g)), system)
            area = sp.sqrt(sp.factor(h.det()))
            curl = sp.diff(A[1], X[1]) - sp.diff(A[0], X[2])
            # Twice the twist times the area, and the plane's Gaussian curvature -m^2, at points of
            # the chart.
            for point in ({X[1]: sp.Rational(3, 10), X[2]: sp.Rational(7, 10), Om: sp.Rational(6, 5), m: sp.Rational(4, 5)},
                          {X[1]: sp.Rational(1, 2), X[2]: sp.Rational(-2, 5), Om: 1, m: 1}):
                self.assertAlmostEqual(float(sp.N((curl - 2 * Om * area).subs(point), 30)), 0.0, places=20, msg=system)
            E, F, G = h[0, 0], h[0, 1], h[1, 1]
            self.assertEqual(F, 0, system)
            # Brioschi's formula for an orthogonal metric E du^2 + G dv^2.
            u, v = X[1], X[2]
            K = -(sp.diff(sp.diff(G, u) / sp.sqrt(E * G), u) + sp.diff(sp.diff(E, v) / sp.sqrt(E * G), v)) / (2 * sp.sqrt(E * G))
            for point in ({u: sp.Rational(3, 10), v: sp.Rational(7, 10), Om: sp.Rational(6, 5), m: sp.Rational(4, 5)},
                          {u: sp.Rational(1, 2), v: sp.Rational(-2, 5), Om: 1, m: 1}):
                self.assertAlmostEqual(float(sp.N((K + m ** 2).subs(point), 30)), 0.0, places=20, msg=system)

    def test_the_scalars_are_those_of_the_family_and_of_anti_de_sitter_space_at_m_2_omega(self):
        sp = self.sp
        for system in CHARTS:
            _, _, Om, m, entry, reader = self.chart(system)
            R = self.scalar(entry, reader, "ricci_scalar")
            K = self.scalar(entry, reader, "kretschmann")
            self.assertEqual(sp.expand(R - 2 * (Om ** 2 - m ** 2)), 0, system)
            self.assertEqual(sp.expand(K - (44 * Om ** 4 - 24 * Om ** 2 * m ** 2 + 4 * m ** 4)), 0, system)
            # Anti-de Sitter space of radius l = 1/Omega: R = -6/l^2 and K = 12/l^4.
            self.assertEqual(sp.expand(R.subs(m, 2 * Om) + 6 * Om ** 2), 0, system)
            self.assertEqual(sp.expand(K.subs(m, 2 * Om) - 12 * Om ** 4), 0, system)
            # Vuorio's member has no curvature scalar, as a vacuum with no cosmological constant.
            self.assertEqual(R.subs(m, Om), 0, system)

    def test_the_horospherical_chart_solves_topologically_massive_gravity(self):
        """G_mu_nu + Lambda g_mu_nu + C_mu_nu/mu = 0 with mu = 3 Omega and Lambda = (Omega^2 - m^2)/3,
        the Cotton tensor C_mu_nu = eps_mu^ab nabla_a (R_b_nu - R g_b_nu/4) and eps^{txy} = 1/sqrt(-g)."""
        sp = self.sp
        g, X, Om, m, _, _ = self.chart("horospherical")
        gi = g.inv()
        n = 3
        gam = [[[sum(gi[a, d] * (sp.diff(g[d, b], X[c]) + sp.diff(g[d, c], X[b]) - sp.diff(g[b, c], X[d]))
                     for d in range(n)) / 2 for c in range(n)] for b in range(n)] for a in range(n)]

        def riemann(a, b, c, d):
            return (sp.diff(gam[a][b][d], X[c]) - sp.diff(gam[a][b][c], X[d])
                    + sum(gam[a][c][e] * gam[e][b][d] - gam[a][d][e] * gam[e][b][c] for e in range(n)))
        ric = sp.Matrix(n, n, lambda b, d: sp.simplify(sum(riemann(a, b, a, d) for a in range(n))))
        R = sp.simplify(sum(gi[a, b] * ric[a, b] for a in range(n) for b in range(n)))
        S = ric - R * g / 4
        root = sp.exp(m * X[1])
        self.assertEqual(sp.simplify(root ** 2 + g.det()), 0)

        def nabla(a, b, c):
            return sp.diff(S[b, c], X[a]) - sum(gam[e][a][b] * S[e, c] + gam[e][a][c] * S[b, e] for e in range(n))
        C = g * sp.Matrix(n, n, lambda i, c: sum(sp.LeviCivita(i, a, b) * nabla(a, b, c) for a in range(n)
                                                 for b in range(n)) / root)
        field = (ric - R * g / 2) + (Om ** 2 - m ** 2) / 3 * g + C / (3 * Om)
        self.assertEqual(field.applyfunc(sp.simplify), sp.zeros(3, 3))

    def test_the_member_m_equal_omega_is_vuorios_solution(self):
        """Chow, Pope and Sezgin's quotation of Vuorio's (2.21), (9/mu^2)[-(dt_V + 2 dtheta -
        2 cosh(sigma) dtheta)^2 + dsigma^2 + sinh^2(sigma) dtheta^2] at mu = 3 Omega, pulled back
        through t_V = Omega t, sigma = Omega r and theta = -phi."""
        sp = self.sp
        g, (t, r, phi), Om, m, _, _ = self.chart("cylindrical")
        tV, sigma, theta = Om * t, Om * r, -phi
        d = [sp.Matrix([[sp.diff(f, v) for v in (t, r, phi)]]) for f in (tV, sigma, theta)]
        fibre = d[0] + 2 * d[2] - 2 * sp.cosh(sigma) * d[2]
        vuorio = 9 / (3 * Om) ** 2 * (-(fibre.T * fibre) + d[1].T * d[1] + sp.sinh(sigma) ** 2 * (d[2].T * d[2]))
        diff = (vuorio - g.subs(m, Om)).applyfunc(lambda v: sp.simplify(v.rewrite(sp.exp)))
        self.assertEqual(diff, sp.zeros(3, 3))

    def test_the_circle_of_constant_t_and_r_is_null_where_tanh_of_mr_over_2_is_m_over_2_omega(self):
        sp = self.sp
        g, (t, r, phi), Om, m, _, _ = self.chart("cylindrical")
        for omega, mass in ((1, 1), (sp.Rational(3, 2), 1), (1, sp.Rational(3, 2))):
            at = {Om: omega, m: mass}
            rc = float(2 * sp.atanh(sp.Rational(mass, 2) / omega) / mass)
            gpp = sp.lambdify(r, g[2, 2].subs(at))
            self.assertAlmostEqual(gpp(rc), 0.0, places=12)
            self.assertGreater(gpp(rc / 2), 0)
            self.assertLess(gpp(1.5 * rc), 0)
        # At Vuorio's member r_c = ln 3/Omega, where sinh^2(Omega r/2) = 1/3.
        self.assertAlmostEqual(2 * math.atanh(0.5), RC, places=14)
        self.assertAlmostEqual(math.sinh(RC / 2) ** 2, 1 / 3, places=14)
        # At m = 2 Omega the circles never close in time: g_phiphi > 0 for every r.
        gpp = sp.lambdify(r, g[2, 2].subs({Om: 1, m: 2}))
        self.assertTrue(all(gpp(x) > 0 for x in (0.1, 1.0, 3.0, 8.0)))


class Drawings(unittest.TestCase):
    """The published diagrams against the metric at m = Omega = 1, from the numbers in the files."""

    @classmethod
    def setUpClass(cls):
        data = drawing("diagrams")
        cls.views = {(system, v["id"]): v for system, views in data["systems"].items() for v in views}
        cls.figure = data["projections"]["cylindrical"][0]

    def slopes(self, view):
        """d(ct)/dX of every ray from end to end, by family; a ray of a view whose metric is the same
        at every point is straight, and its points are rounded to a ten thousandth of the drawing,
        so a ray shorter than a tenth is passed over."""
        X0, X1, Y0, Y1 = view["box"]
        found = {}
        for family, lines in view["rays"].items():
            for line in lines:
                (ua, wa), (ub, wb) = line[0], line[-1]
                if abs(ub - ua) > 0.1:
                    found.setdefault(family, []).append((wb - wa) * (Y1 - Y0) / ((ub - ua) * (X1 - X0)))
        return found

    def test_the_rays_of_each_cylinder_are_the_null_lines_of_the_metric(self):
        for view, r in (("inside", RC / 2), ("beyond", 3 * RC / 2)):
            found = self.slopes(self.views["cylindrical", view])
            self.assertEqual(len(found), 2, view)
            H, D = twist_and_spread(r)
            # c dt = (-H +- D) dphi, drawn against r phi.
            want = {(-H + D) / r, (-H - D) / r}
            for family, slopes in found.items():
                nearest = min(want, key=lambda k: abs(k - slopes[0]))
                want.discard(nearest)
                for k in slopes:
                    self.assertAlmostEqual(k, nearest, delta=2e-3 * max(1.0, abs(nearest)), msg=f"{view} {family}")
                self.assertAlmostEqual(-(nearest * r + H) ** 2 + D * D, 0.0, places=12)
            self.assertEqual(want, set(), view)

    def test_the_slopes_the_captions_state(self):
        H, D = twist_and_spread(RC / 2)
        self.assertAlmostEqual(-H + D, 2 - math.sqrt(3), places=14)
        self.assertAlmostEqual(-H - D, -(5 / math.sqrt(3) - 2), places=14)
        H, D = twist_and_spread(3 * RC / 2)
        self.assertAlmostEqual(-H + D, -(5 / math.sqrt(3) - 2), places=14)
        self.assertAlmostEqual(-H - D, -(41 / (3 * math.sqrt(3)) - 2), places=14)

    def test_beyond_the_null_circle_both_families_run_down_in_t_toward_plus_phi(self):
        inside = self.slopes(self.views["cylindrical", "inside"])
        beyond = self.slopes(self.views["cylindrical", "beyond"])
        self.assertEqual(sorted(k[0] > 0 for k in inside.values()), [False, True])
        self.assertTrue(all(k[0] < 0 for k in beyond.values()))

    def test_the_planes_of_the_fibred_and_horospherical_charts_run_at_45_degrees(self):
        for key in (("fibred", "plane"), ("horospherical", "tx")):
            for slopes in self.slopes(self.views[key]).values():
                for k in slopes:
                    self.assertAlmostEqual(abs(k), 1.0, delta=2e-3, msg=str(key))

    def test_the_rays_of_poincares_disc_keep_t_plus_or_minus_2_artanh_R(self):
        view = self.views["disc", "plane"]
        # A mirrored view keeps its rays on the half R >= 0 of the box; the page reflects them.
        X0, X1, Y0, Y1 = view["box"]
        worst = 0.0
        for family, lines in view["rays"].items():
            for line in lines:
                kept = []
                for u, w in line:
                    R, t = X0 + u * (X1 - X0), Y0 + w * (Y1 - Y0)
                    if abs(R) < 0.9:
                        kept.append((t + 2 * math.atanh(R), t - 2 * math.atanh(R)))
                if len(kept) > 3:
                    spread = min(max(k[i] for k in kept) - min(k[i] for k in kept) for i in (0, 1))
                    worst = max(worst, spread)
        self.assertLess(worst, 2e-2)

    def test_the_figure_draws_the_null_circle_at_r_c_and_every_cone_null(self):
        turn = self.figure["turn"]
        critical = [line["points"] for line in turn["lines"] if line["class"] == "critical"]
        self.assertEqual(len(critical), 1)
        for X, Y, T in critical[0]:
            self.assertAlmostEqual(math.hypot(X, Y), RC, delta=2e-6)
            self.assertAlmostEqual(T, 0.0, delta=1e-9)
        # A generator (dX, dY, dT) at (X, Y) is null: with r dr = X dX + Y dY and r^2 dphi = X dY - Y dX,
        # -(dT + H dphi)^2 + dr^2 + D^2 dphi^2 = 0.
        self.assertEqual(len(turn["cones"]), 13)
        for cone in turn["cones"]:
            X, Y, T = cone["apex"]
            r = math.hypot(X, Y)
            H, D = twist_and_spread(r)
            for x, y, t in cone["rim"]:
                dX, dY, dT = x - X, y - Y, t - T
                size = dX * dX + dY * dY + dT * dT
                if r < 1e-4:
                    null = -dT * dT + dX * dX + dY * dY
                else:
                    dr, dphi = (X * dX + Y * dY) / r, (X * dY - Y * dX) / r ** 2
                    null = -(dT + H * dphi) ** 2 + dr * dr + D * D * dphi * dphi
                self.assertAlmostEqual(null / size, 0.0, delta=3e-4)

    def test_the_embedding_follows_the_circles_and_stops_where_the_caption_says(self):
        piece = drawing("embedding")["views"][0]["surfaces"][0]["pieces"][0]
        for r, rho, _ in piece["points"]:
            s2 = math.sinh(r / 2) ** 2
            self.assertAlmostEqual(rho, 2 * math.sqrt(s2 * max(1 - 3 * s2, 0.0)), delta=1e-6, msg=f"rho at {r}")
        self.assertAlmostEqual(math.sinh(piece["points"][-1][0] / 2) ** 2, (math.sqrt(3) - 1) / 3, places=10)
        self.assertAlmostEqual(max(p[1] for p in piece["points"]), 1 / math.sqrt(3), delta=1e-5)
        # The drawing stops inside the null circle, where a surface of constant t is still space.
        self.assertLess(piece["points"][-1][0], RC)


if __name__ == "__main__":
    unittest.main()
