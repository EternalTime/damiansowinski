"""Lin, Lunin and Maldacena's bubbling anti-de Sitter space:
python3 -m unittest discover -s _tools

What its texts state, held on the published files: the chart of concentric droplets is their
metric (2.4) on the polar coordinates of the plane, and at the disc of radius L^2 it is the global
chart of anti-de Sitter space times a 5-sphere carried along y = L^2 sinh(rho) sin(theta),
r = L^2 cosh(rho) cos(theta), phi = psi - t; the global chart has R^M_N = -4/L^2 on anti-de Sitter
space and +4/L^2 on the sphere and the Kretschmann scalar 80/L^4; the plane wave's one Ricci component
is R_tt = 8; the disc, the black ring drawn and the half plane solve the field equation of type IIB
supergravity with Lin, Lunin and Maldacena's 5-form, and the ring's G and V obey the first order
equations the chart's parameters state; and the drawings hold what their captions say. The tests of
the published mathematics need sympy and are skipped where it is absent; the tests of the drawings
read the files alone.
"""
import importlib.util
import json
import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "MFS" / "assets" / "data"
HAS_SYMPY = importlib.util.find_spec("sympy") is not None and importlib.util.find_spec("scipy") is not None
sys.path.insert(0, str(Path(__file__).resolve().parent / "derivations"))


def published():
    return json.loads((DATA / "metrics" / "bubbling_ads.json").read_text(encoding="utf-8"))


def drawing(kind):
    return json.loads((DATA / kind / "bubbling_ads.json").read_text(encoding="utf-8"))


def chart(system):
    return next(c for c in published()["coordinates"] if c["id"] == system)


@unittest.skipUnless(HAS_SYMPY, "sympy is not installed")
class Published(unittest.TestCase):
    def setUp(self):
        import sympy
        import bubbling_ads
        import verify_metrics
        self.sp, self.ba, self.vm = sympy, bubbling_ads, verify_metrics

    def read(self, system):
        entry = chart(system)
        reader = self.vm.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]], ())
        n = len(entry["coords"])
        g = self.sp.zeros(n)
        for e in entry["metric_components"]:
            i, j = (entry["coords"].index(k) for k in e["indices"])
            g[i, j] = g[j, i] = reader(e["value"])
        return reader, entry, g

    def test_the_chart_of_concentric_droplets_is_lin_lunin_and_maldacena_s_metric(self):
        reader, entry, g = self.read("rings")
        self.assertEqual(entry["coords"][:4], ["t", "r", "\\phi", "y"])
        symbols = [reader.symbol[c] for c in entry["coords"]]
        G, V = reader.parameters["G"], reader.parameters["V"]
        theirs = self.ba.metric("polar", symbols[:4], G, V).subs(dict(zip(self.ba.ANGLES, symbols[4:])))
        plain = {G: self.sp.Symbol("_G"), V: self.sp.Symbol("_V")}
        for i in range(10):
            for j in range(10):
                self.assertEqual(self.sp.simplify((g[i, j] - theirs[i, j]).subs(plain).rewrite(self.sp.exp)), 0,
                                 f"g_{entry['coords'][i]}{entry['coords'][j]}")

    def test_the_disc_is_the_global_chart_pulled_back(self):
        sp = self.sp
        reader, entry, g = self.read("rings")
        greader, gentry, gg = self.read("global")
        s = [reader.symbol[c] for c in entry["coords"]]
        gs = [greader.symbol[c] for c in gentry["coords"]]
        t, r, phi, y = s[:4]
        T, rho, theta, psi = gs[0], gs[1], gs[5], gs[6]
        L = sp.Rational(13, 10)
        G, V = self.ba.concentric_polar(r, y, [L ** 2])
        disc = g.subs({reader.parameters["G"]: G, reader.parameters["V"]: V}, simultaneous=True)
        image = {t: T, r: L ** 2 * sp.cosh(rho) * sp.cos(theta), phi: psi - T, y: L ** 2 * sp.sinh(rho) * sp.sin(theta)}
        image.update(dict(zip(s[4:], gs[2:5] + gs[7:])))
        order = [0, 1, 5, 6, 2, 3, 4, 7, 8, 9]
        plane = [T, rho, theta, psi]
        J = sp.zeros(10)
        for i in range(4):
            for a in range(4):
                J[i, a] = sp.diff(image[s[i]], plane[a])
        for k in range(4, 10):
            J[k, k] = 1
        at = {T: sp.Rational(3, 10), rho: sp.Rational(7, 10), theta: sp.Rational(6, 10), psi: sp.Rational(11, 10),
              **{a: sp.Rational(k + 5, 10) for k, a in enumerate(gs[2:5] + gs[7:])}}
        pulled = (J.T * disc.subs(image, simultaneous=True) * J).subs(at).evalf(30)
        want = gg.subs(greader.parameters["L"], L).subs(at).evalf(30)
        for i in range(10):
            for j in range(10):
                self.assertLess(abs(pulled[i, j] - want[order[i], order[j]]), 1e-25)

    def test_the_global_chart_is_anti_de_sitter_space_times_a_sphere_of_one_radius(self):
        entry = chart("global")
        self.assertEqual(entry["kretschmann"], "K = \\dfrac{80}{L^4}")
        self.assertEqual(entry["ricci_scalar"], "R = 0")
        reader, _, g = self.read("global")
        L = reader.parameters["L"]
        coords = entry["coords"]
        for e in entry["ricci_tensor"]["variants"]["ul"]["nonzero"]:
            i, j = (coords.index(k) for k in e["indices"])
            self.assertEqual(i, j)
            self.assertEqual(self.sp.simplify(reader(e["value"]) - (-4 if i < 5 else 4) / L ** 2), 0)

    def test_the_plane_wave_has_one_ricci_component(self):
        entry = chart("plane_wave")
        self.assertEqual(entry["kretschmann"], "K = 0")
        ricci = entry["ricci_tensor"]["variants"]["ll"]["nonzero"]
        self.assertEqual([(e["indices"], e["value"]) for e in ricci], [(["t", "t"], "8")])

    def test_the_members_solve_the_field_equation(self):
        sp, ba = self.sp, self.ba
        r, y = ba.generic("polar")[0][1], ba.generic("polar")[0][3]
        for radii in ([sp.Rational(3, 2)], list(ba.RING)):
            residual = ba.field_equation_residual("polar", ba.concentric_polar(r, y, radii), [(0.4, 0.3, 0.7)])
            self.assertLess(residual, 1e-9, radii)
        x, w, y = ba.generic("cartesian")[0][1:]
        self.assertLess(ba.field_equation_residual("cartesian", ba.half_plane(x, w, y), [(0.4, 0.9, 0.7)]), 1e-9)
        # The opposite twist is no solution, so the test can fail.
        G, V = ba.concentric_polar(r, y, [sp.Rational(3, 2)])
        self.assertGreater(ba.field_equation_residual("polar", (G, -V), [(0.4, 0.3, 0.7)]), 1e-2)

    def test_the_ring_obeys_the_equations_its_parameters_state(self):
        sp, ba = self.sp, self.ba
        r, y = sp.symbols("r y", positive=True)
        G, V = ba.concentric_polar(r, y, ba.RING)
        z = sp.tanh(G) / 2
        at = {r: sp.Rational(7, 10), y: sp.Rational(4, 10)}
        for law in (sp.diff(z, r, 2) + sp.diff(z, r) / r + y * sp.diff(sp.diff(z, y) / y, y),
                    y * sp.diff(V, y) - r * sp.diff(z, r), y * sp.diff(V, r) + r * sp.diff(z, y)):
            self.assertLess(abs(sp.N(law.subs(at), 30)), 1e-20)
        # On the plane y = 0 the ring is white inside r = L^2, black out to sqrt(2) L^2, and white beyond.
        for radius, colour in ((0.5, 0.5), (1.2, -0.5), (2.0, 0.5)):
            self.assertAlmostEqual(float(z.subs({r: radius, y: sp.Rational(1, 10 ** 6)})), colour, places=9)


class Drawings(unittest.TestCase):
    def test_a_ray_on_the_ring_s_axis_reaches_the_boundary_when_the_captions_say(self):
        if not HAS_SYMPY:
            self.skipTest("scipy is not installed")
        import bubbling_ads
        edge = bubbling_ads.axis_tortoise(math.inf)
        self.assertAlmostEqual(edge, 1.2526, places=4)
        self.assertAlmostEqual(bubbling_ads.axis_tortoise(math.inf, (1,)), math.pi / 2, places=10)
        views = drawing("diagrams")["systems"]["rings"]
        self.assertIn("$t = 1.253$", " ".join(views[0]["caption"]))
        rings = next(v for v in drawing("conformal")["views"] if v["id"] == "rings")
        self.assertIn("$0 \\le y_* < 1.253$", " ".join(rings["caption"]))

    def test_the_plane_of_the_disc_is_a_hemisphere_on_a_cylinder_of_radius_l(self):
        view = drawing("embedding")["views"][0]
        disc, white = view["surfaces"][0]["pieces"]
        for x, rho, z in disc["points"]:
            self.assertAlmostEqual(rho, math.cos(x), places=4)
            self.assertAlmostEqual(z, math.sin(x), places=4)
        for x, rho, z in white["points"]:
            self.assertAlmostEqual(rho, 1.0, places=6)
            self.assertAlmostEqual(z, -x, places=4)
        self.assertEqual((disc["points"][0][1], disc["points"][0][2]), (white["points"][0][1], white["points"][0][2]))


if __name__ == "__main__":
    unittest.main()
