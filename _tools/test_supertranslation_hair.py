"""Compère and Long's flat space with supertranslation hair:
python3 -m unittest discover -s _tools

What its texts and drawings state, held on the published files. The drawings take the field
C = ell (3 cos^2 theta - 1)/2 at ell = 1. The shells of constant rho are held to Compère and Long's
Cartesian map, X = (rho - C) n - C' e_theta, and their quoted sizes to it; the supertranslation
horizon on the equator is held to rho = 5/2 and to its ring of radius 3; the light rays of the Bondi
planes to the closed forms the captions give; and the dots of the conformal views to the events
the shells are cut at. Where sympy is installed, every published chart is held to a vanishing
Riemann tensor, the static chart to the flat metric pulled back through the map, and each Bondi
chart to the static chart pulled back through t = u + rho or t = v - rho; those tests are skipped
where it is absent.
"""
import importlib.util
import json
import math
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "MFS" / "assets" / "data"
HAS_SYMPY = importlib.util.find_spec("sympy") is not None


def C(th):
    return (3 * math.cos(th) ** 2 - 1) / 2


def dC(th):
    return -3 * math.sin(th) * math.cos(th)


def d2C(th):
    return -3 * math.cos(2 * th)


def cartesian(rho, th, ph):
    """Compère and Long's map, their 1602.05197 section 3.1 with C axisymmetric."""
    n = (math.sin(th) * math.cos(ph), math.sin(th) * math.sin(ph), math.cos(th))
    e = (math.cos(th) * math.cos(ph), math.cos(th) * math.sin(ph), -math.sin(th))
    return tuple((rho - C(th)) * n[i] - dC(th) * e[i] for i in range(3))


def load(kind):
    return json.loads((DATA / kind / "supertranslation_hair.json").read_text(encoding="utf-8"))


class Shells(unittest.TestCase):
    def test_the_map_carries_the_flat_metric_onto_the_static_chart(self):
        # Finite differences of the map give g_rhorho = 1, g_thth = (rho - C - C'')^2,
        # g_phph = B^2 and no cross terms.
        h = 1e-6
        for rho, th in ((3.0, 0.4), (4.2, 1.3), (6.0, 2.7)):
            def d(k):
                up, down = [rho, th, 0.3], [rho, th, 0.3]
                up[k] += h
                down[k] -= h
                return [(a - b) / (2 * h) for a, b in zip(cartesian(*up), cartesian(*down))]
            tangents = [d(k) for k in range(3)]
            g = [[sum(a * b for a, b in zip(tangents[i], tangents[j])) for j in range(3)] for i in range(3)]
            B = (rho - C(th)) * math.sin(th) - dC(th) * math.cos(th)
            expected = [[1, 0, 0], [0, (rho - C(th) - d2C(th)) ** 2, 0], [0, 0, B * B]]
            for i in range(3):
                for j in range(3):
                    self.assertAlmostEqual(g[i][j], expected[i][j], places=6)

    def test_each_shell_is_oblate_by_the_sizes_the_caption_quotes(self):
        for rho in (2.6, 3.0, 4.5, 6.0):
            equator = math.hypot(*cartesian(rho, math.pi / 2, 0.0)[:2])
            pole = cartesian(rho, 0.0, 0.0)[2]
            self.assertAlmostEqual(equator, rho + 0.5, places=12)
            self.assertAlmostEqual(pole, rho - 1, places=12)
            self.assertAlmostEqual(cartesian(rho, math.pi, 0.0)[2], -(rho - 1), places=12)

    def test_the_embedding_draws_each_shell_where_the_map_puts_it(self):
        view = load("embedding")["views"][0]
        self.assertEqual(view["id"], "shells")
        self.assertEqual(view["movie"]["variable"], "$\\rho$")
        for surface in view["surfaces"]:
            rho = surface["time"]
            for th, radius, z in surface["pieces"][0]["points"]:
                X = cartesian(rho, th, 0.0)
                self.assertAlmostEqual(radius, math.hypot(X[0], X[1]), delta=2e-3 * rho)
                self.assertAlmostEqual(z, -X[2], delta=2e-3 * rho)

    def test_the_meridian_bends_hardest_at_the_equator(self):
        # The principal radius rho - C - C'' = rho - (5 - 9 cos^2 theta)/2 is least at the equator,
        # rho - 5/2, and the other, B/sin(theta) = rho + (1 + 3 cos^2 theta)/2, stays positive.
        for rho in (2.6, 4.0):
            radii = [rho - C(th) - d2C(th) for th in [k * math.pi / 200 for k in range(1, 200)]]
            self.assertAlmostEqual(min(radii), rho - 2.5, places=3)
            self.assertAlmostEqual(rho - C(math.pi / 2) - d2C(math.pi / 2), rho - 2.5, places=12)


class Horizon(unittest.TestCase):
    def test_the_supertranslation_horizon_on_the_equator(self):
        # rho = C + max(C'', C' cot theta): on the equator 5/2, a ring of radius 3 in flat space,
        # where W = sqrt(r^2 + sigma^2) is 3/2 with sigma = -3/2 and r = 0.
        th = math.pi / 2
        edge = C(th) + max(d2C(th), 0.0)
        self.assertAlmostEqual(edge, 2.5, places=12)
        self.assertAlmostEqual(math.hypot(*cartesian(edge, th, 0.0)[:2]), 3.0, places=12)
        sigma = (dC(th) * math.cos(th) / math.sin(th) - d2C(th)) / 2
        self.assertAlmostEqual(sigma, -1.5, places=12)
        P = C(th) + (d2C(th) + dC(th) * math.cos(th) / math.sin(th)) / 2
        self.assertAlmostEqual(abs(sigma) + P, edge, places=12)

    def test_the_static_view_hatches_inside_the_horizon(self):
        view = load("diagrams")["systems"]["static"][0]
        self.assertLess(view["box"][0], 2.5)


class Rays(unittest.TestCase):
    def test_the_ingoing_rays_of_the_retarded_plane_keep_u_plus_2W(self):
        # On the equator g_uu = -1 and g_ur = -r/W, so a null curve off u = const has
        # du/dr = -2r/W, and u + 2W is constant along it.
        for r in (0.1, 1.0, 3.7):
            W = math.sqrt(r * r + 2.25)
            h = 1e-6
            dW = (math.sqrt((r + h) ** 2 + 2.25) - math.sqrt((r - h) ** 2 + 2.25)) / (2 * h)
            self.assertAlmostEqual(-1 * (-2 * dW) ** 2 - 2 * (r / W) * (-2 * dW), 0.0, places=8)

    def test_the_dots_are_the_shells_on_every_plane(self):
        moments = [s["time"] for s in load("embedding")["views"][0]["surfaces"]]
        diagrams = load("diagrams")["systems"]
        for system, place in (("static", lambda rho: (rho, 0.0)),
                              ("bondi_retarded", lambda rho: (math.sqrt((rho - 1) ** 2 - 2.25), -rho)),
                              ("bondi_advanced", lambda rho: (math.sqrt((rho - 1) ** 2 - 2.25), rho))):
            view = diagrams[system][0]
            X0, X1, Y0, Y1 = view["box"]
            for rho, mark in zip(moments, view["slices"]):
                (u,) = mark["points"]
                X, Y = place(rho)
                self.assertAlmostEqual(X0 + u[0] * (X1 - X0), X, delta=2e-3)
                self.assertAlmostEqual(Y0 + u[1] * (Y1 - Y0), Y, delta=2e-3)

    def test_the_conformal_dots_are_the_shells(self):
        # p, q = arctan(t -+ (rho - 5/2)) at t = 0, drawn at X = q - p, T = p + q.
        moments = [s["time"] for s in load("embedding")["views"][0]["surfaces"]]
        for view in load("conformal")["views"]:
            for rho, mark in zip(moments, view["slices"]):
                (at,) = mark["points"]
                self.assertAlmostEqual(at[0], 2 * math.atan(rho - 2.5), delta=1e-3)
                self.assertAlmostEqual(at[1], 0.0, delta=1e-3)


@unittest.skipUnless(HAS_SYMPY, "sympy is not installed")
class Charts(unittest.TestCase):
    def setUp(self):
        import sys
        sys.path.insert(0, str(ROOT / "_tools" / "derivations"))
        import print_charts as pc
        import chart_printer as cp
        import verify_metrics as vm
        self.pc, self.cp, self.vm = pc, cp, vm

    def test_every_published_chart_is_flat_and_checked_against_its_source(self):
        published = {c["id"]: c for c in json.loads(
            (DATA / "metrics" / "supertranslation_hair.json").read_text(encoding="utf-8"))["coordinates"]}
        for system in self.pc.STH_CHARTS:
            spec = self.pc.supertranslation_hair(system)
            self.assertEqual(published[system]["line_element"], spec["system"]["line_element"])
            self.assertEqual(published[system]["riemann"]["variants"]["ulll"]["nonzero"], [])
            chart = self.cp.Chart(spec["system"]["coords"], spec["system"]["parameters"], spec["chart_line_element"],
                                  spec["printer"], spec.get("pretty"), None, None, spec.get("reduce"))
            # Flatness, the Cartesian map, each Bondi chart pulled back from the static one, and the
            # retarded chart as a member of Bondi's published metric.
            spec["check"](chart)


if __name__ == "__main__":
    unittest.main()
