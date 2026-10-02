"""The white hole, Oppenheimer and Snyder's ball of dust with the time reversed:
python3 -m unittest discover -s _tools

What its texts state, held on the published files. The surface of the core is in no chart's
mathematics, so the junction is taken here from the published metrics and Christoffel symbols of
the dust and of Schwarzschild's exterior: along the cycloid a = a_m sin^2(eta/2),
c tau = (a_m/2)(eta - sin eta), the surface chi = chi_0 has the areal radius R = a sin chi_0, is a
radial geodesic of the exterior of energy cos chi_0 where r_s = a_m sin^3 chi_0, and carries no
surface layer, K^theta_theta being cos chi_0 / R from both sides. The same files give the moments
the captions name: the surface crosses r_s at eta = 2 chi_0 and 2 pi - 2 chi_0, where
|grad R|^2 = 0 inside the dust, and for chi_0 = pi/4 it is at rest at u = (pi + ln 2) r_s of the
retarded time that counts from the first crossing, and leaves r = 0 at u = (2 + ln 2 - pi) r_s.
Novikov's chart is held to Tolman's density, to Schwarzschild's vacuum where the mass is constant,
and to R = F where q = 2F/3, the Schwarzschild sphere. The drawings are held to the same numbers.
The tests of the published mathematics need sympy and are skipped where it is absent.
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

CHI0 = math.pi / 4                                   # the core every drawing declares
U_REST = math.pi + math.log(2)                       # u of the moment of rest, in r_s
U_BANG = 2 + math.log(2) - math.pi                   # u at which the surface leaves r = 0, in r_s


def chart(system):
    metric = json.loads((DATA / "metrics" / "white_hole.json").read_text(encoding="utf-8"))
    return next(c for c in metric["coordinates"] if c["id"] == system)


@unittest.skipUnless(HAS_SYMPY, "sympy is not installed")
class Junction(unittest.TestCase):
    def setUp(self):
        import sympy
        import verify_metrics
        self.sp, self.vm = sympy, verify_metrics

    def read(self, system, held=()):
        entry = chart(system)
        reader = self.vm.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]], (), held=held)
        return entry, reader

    def components(self, entry, reader, field, variant=None):
        block = entry[field] if variant is None else entry[field]["variants"][variant]["nonzero"]
        return {tuple(e["indices"]): reader(e["value"]) for e in block}

    def cycloid(self, reader, eta):
        """a, da/d(c tau) and d^2a/d(c tau)^2 along a = a_m sin^2(eta/2), c dtau = a d eta."""
        sp = self.sp
        tau, a, am = reader.symbol["\\tau"], reader.parameters["a"], reader.parameters["a_m"]
        return {sp.Derivative(a, (tau, 2)): -1 / (2 * am * sp.sin(eta / 2) ** 4),
                sp.Derivative(a, tau): sp.cos(eta / 2) / sp.sin(eta / 2), a: am * sp.sin(eta / 2) ** 2}

    def test_the_core_is_dust_of_the_mass_the_exterior_has(self):
        sp = self.sp
        entry, reader = self.read("interior_comoving")
        eta = sp.Symbol("eta", positive=True)
        on = self.cycloid(reader, eta)
        einstein = self.components(entry, reader, "einstein_tensor", "ul")
        for name in ("\\chi", "\\theta", "\\phi"):
            self.assertEqual(sp.simplify(einstein[(name, name)].subs(on)), 0, f"a pressure along {name}")
        # 8 pi G rho / c^2 = -G^tau_tau = 3 a_m / a^3, and (8 pi G rho / 3 c^2) R^3 = r_s = a_m sin^3 chi_0.
        am, chi0 = reader.parameters["a_m"], reader.parameters["chi_0"]
        density = -einstein[("\\tau", "\\tau")].subs(on)
        radius = (am * sp.sin(eta / 2) ** 2) * sp.sin(chi0)
        self.assertEqual(sp.simplify(density * radius ** 3 / 3 - am * sp.sin(chi0) ** 3), 0)

    def test_the_conformal_chart_is_the_same_dust(self):
        sp = self.sp
        entry, reader = self.read("interior_conformal")
        eta, am = reader.symbol["\\eta"], reader.parameters["a_m"]
        g = self.components(entry, reader, "metric_components")
        self.assertEqual(sp.simplify(g[("\\eta", "\\eta")] + (am * sp.sin(eta / 2) ** 2) ** 2), 0)
        einstein = self.components(entry, reader, "einstein_tensor", "ul")
        self.assertEqual(list(einstein), [("\\eta", "\\eta")])
        self.assertEqual(sp.simplify(einstein[("\\eta", "\\eta")] + 3 / (am ** 2 * sp.sin(eta / 2) ** 6)), 0)
        self.assertEqual(entry["weyl_tensor"]["variants"]["llll"]["nonzero"], [])

    def surface(self):
        """The surface in the exterior's chart x^0 = ct, as functions of eta: its areal radius, its
        velocity (dt/dtau, dr/dtau) from the energy cos chi_0, and the exterior's f = 1 - r_s/R."""
        sp = self.sp
        eta = sp.Symbol("eta", positive=True)
        am, chi0 = sp.Symbol("a_m", positive=True), sp.Symbol("chi_0", positive=True)
        a = am * sp.sin(eta / 2) ** 2
        R, rs, E = a * sp.sin(chi0), am * sp.sin(chi0) ** 3, sp.cos(chi0)
        f = 1 - rs / R
        rdot = sp.diff(R, eta) / a                       # c dtau = a d eta
        return eta, a, R, rs, E, f, E / f, rdot

    def test_the_surface_is_a_radial_geodesic_of_the_exterior(self):
        sp = self.sp
        entry, reader = self.read("exterior_schwarzschild")
        eta, a, R, rs, E, f, tdot, rdot = self.surface()
        at = {reader.symbol["r"]: R, reader.parameters["r_s"]: rs}
        g = {k: v.subs(at) for k, v in self.components(entry, reader, "metric_components").items()}
        gamma = {k: v.subs(at) for k, v in self.components(entry, reader, "christoffel", "ull").items()}
        speed = g[("t", "t")] * tdot ** 2 + g[("r", "r")] * rdot ** 2
        self.assertEqual(sp.simplify(speed + 1), 0, "the surface does not keep unit speed with the energy cos chi_0")
        radial = (sp.diff(rdot, eta) / a + gamma[("r", "t", "t")] * tdot ** 2 + gamma[("r", "r", "r")] * rdot ** 2)
        self.assertEqual(sp.simplify(radial), 0, "the surface is not a geodesic of the published exterior")

    def test_the_surface_carries_no_layer(self):
        """K^theta_theta = -Gamma^mu_theta_theta n_mu / R^2 with the unit normal toward the outside:
        cos chi_0 / R from the dust, where n = d chi / a, and from the vacuum, where n_mu is
        (-dr/dtau, dt/dtau)."""
        sp = self.sp
        inside, reader_in = self.read("interior_comoving")
        outside, reader_out = self.read("exterior_schwarzschild")
        eta, a, R, rs, E, f, tdot, rdot = self.surface()
        chi0 = sp.Symbol("chi_0", positive=True)
        gamma_in = self.components(inside, reader_in, "christoffel", "ull")
        on = self.cycloid(reader_in, eta)
        on[reader_in.symbol["\\chi"]] = reader_in.parameters["chi_0"]
        k_in = (-gamma_in[("\\chi", "\\theta", "\\theta")] * reader_in.parameters["a"]).subs(on) / R ** 2
        k_in = k_in.subs({reader_in.parameters["chi_0"]: chi0, reader_in.parameters["a_m"]: sp.Symbol("a_m", positive=True)})
        gamma_out = self.components(outside, reader_out, "christoffel", "ull")
        at = {reader_out.symbol["r"]: R, reader_out.parameters["r_s"]: rs}
        k_out = -(gamma_out[("r", "\\theta", "\\theta")].subs(at) * tdot) / R ** 2
        self.assertEqual(sp.simplify(k_in - sp.cos(chi0) / R), 0)
        self.assertEqual(sp.simplify(k_out - sp.cos(chi0) / R), 0)

    def test_the_surface_crosses_the_schwarzschild_sphere_where_the_caption_says(self):
        sp = self.sp
        eta, a, R, rs, *_ = self.surface()
        chi0 = sp.Symbol("chi_0", positive=True)
        for crossing in (2 * chi0, 2 * sp.pi - 2 * chi0):
            self.assertEqual(sp.simplify((R - rs).subs(eta, crossing)), 0)
        # Inside the dust the areal radius a sin chi has a null gradient on eta = 2 chi and 2 pi - 2 chi.
        entry, reader = self.read("interior_conformal")
        e, chi = reader.symbol["\\eta"], reader.symbol["\\chi"]
        inverse = self.components(entry, reader, "inverse_metric_components")
        square = self.components(entry, reader, "metric_components")[("\\theta", "\\theta")]
        # |grad R|^2 = g^ab d_a(R^2) d_b(R^2) / 4R^2, which takes no root.
        gradient = (inverse[("\\eta", "\\eta")] * sp.diff(square, e) ** 2
                    + inverse[("\\chi", "\\chi")] * sp.diff(square, chi) ** 2) / (4 * square)
        for line in (2 * chi, 2 * sp.pi - 2 * chi):
            self.assertEqual(sp.simplify(sp.expand_trig(gradient.subs(e, line))), 0)

    def test_the_retarded_time_counts_from_the_crossing(self):
        """Along the surface on its way out du/dtau = 1/(E + sqrt(E^2 - f)) by the published outgoing
        chart, regular at r_s, so u runs from its value at r = 0 through 0 at the crossing to the
        moment of rest: for chi_0 = pi/4, (2 + ln 2 - pi) r_s and (pi + ln 2) r_s."""
        sp = self.sp
        entry, reader = self.read("exterior_eddington_finkelstein")
        g = self.components(entry, reader, "metric_components")
        r, rs = reader.symbol["r"], reader.parameters["r_s"]
        E, udot, rdot = sp.symbols("E udot rdot", positive=True)
        f = 1 - rs / r
        rate = 1 / (E + sp.sqrt(E ** 2 - f))
        # Unit speed and the energy E = -(g_uu du/dtau + g_ur dr/dtau), with dr/dtau = sqrt(E^2 - f).
        out = {udot: rate, rdot: sp.sqrt(E ** 2 - f)}
        speed = g[("u", "u")] * udot ** 2 + 2 * g[("u", "r")] * udot * rdot
        point = {r: sp.Rational(3, 2), rs: 1, E: 1 / sp.sqrt(2)}
        self.assertAlmostEqual(float((speed.subs(out) + 1).subs(point)), 0, places=12)
        self.assertAlmostEqual(float((-(g[("u", "u")] * udot + g[("u", "r")] * rdot).subs(out) - E).subs(point)), 0, places=12)
        am, energy = 1 / math.sin(CHI0) ** 3, math.cos(CHI0)        # r_s = 1

        def du(eta):
            a = am * math.sin(eta / 2) ** 2
            radius = a * math.sin(CHI0)
            return a / (energy + math.sqrt(max(energy ** 2 - (1 - 1 / radius), 0.0)))

        def simpson(lo, hi, n=20000):
            h = (hi - lo) / n
            return h / 3 * (du(lo) + du(hi) + sum((4 if k % 2 else 2) * du(lo + k * h) for k in range(1, n)))
        self.assertAlmostEqual(simpson(2 * CHI0, math.pi), U_REST, places=7)
        self.assertAlmostEqual(-simpson(1e-9, 2 * CHI0), U_BANG, places=6)

    def test_kruskals_radius_is_the_schwarzschild_radius_on_the_horizons_and_zero_on_the_singularity(self):
        sp = self.sp
        entry, reader = self.read("exterior_kruskal", held=("r",))
        U, V, rs = reader.symbol["U"], reader.symbol["V"], reader.parameters["r_s"]
        radius = reader.parameters["r"].subs(reader.held).doit()
        self.assertEqual(sp.simplify(radius.subs(V, 0) - rs), 0)
        self.assertEqual(sp.simplify(radius.subs(U, 0) - rs), 0)
        self.assertEqual(sp.simplify(radius.subs({U: 1, V: 1})), 0)
        self.assertEqual(entry["ricci_tensor"]["variants"]["ll"]["nonzero"], [])
        # The declared derivatives of r, which its description states, are its definition's.
        E = sp.exp(-radius / rs)
        for variable, other in ((U, V), (V, U)):
            self.assertEqual(self.vm.norm(sp.diff(radius, variable) + rs ** 2 * other * E / radius), 0)

    def test_novikovs_chart_is_tolmans_dust_and_schwarzschilds_vacuum(self):
        sp = self.sp
        entry, reader = self.read("novikov_comoving")
        t, r = reader.symbol["t"], reader.symbol["r"]
        F, b, q, R = (reader.parameters[n] for n in ("F", "b", "q", "R"))
        dF, db = sp.Derivative(F, r), sp.Derivative(b, r)
        x = sp.Symbol("x", positive=True)
        explicit = {R: (sp.Rational(9, 4) * F) ** sp.Rational(1, 3) * x ** sp.Rational(2, 3), q: x}
        # Tolman's density, F'/(R^2 d_r R), with d_r R as the chart's description states it.
        einstein = self.components(entry, reader, "einstein_tensor", "ll")
        self.assertEqual(list(einstein), [("t", "t")])
        slope = R * (q * dF - 2 * F * db) / (3 * F * q)
        self.assertEqual(sp.simplify((einstein[("t", "t")] - dF / (R ** 2 * slope)).subs(explicit)), 0)
        # A vacuum where F is constant, with Schwarzschild's Kretschmann scalar at the areal radius R.
        K = reader(entry["kretschmann"].split("=", 1)[1])
        self.assertEqual(sp.simplify((K.subs(dF, 0) - 12 * F ** 2 / R ** 6).subs(explicit)), 0)
        # The areal radius has a null gradient where R = F, which is q = 2F/3: -(d_0 R)^2 + 1 = 0.
        gradient = 1 - (2 * R / (3 * q)) ** 2
        self.assertEqual(sp.simplify(gradient.subs(explicit).subs(x, 2 * F / 3)), 0)
        self.assertEqual(sp.simplify((R - F).subs(explicit).subs(x, 2 * F / 3)), 0)
        # With no delay and uniform dust it is the Einstein-de Sitter universe, K = 80/(27 q^4).
        self.assertEqual(sp.simplify(K.subs(db, 0) - sp.Rational(80, 27) / q ** 4), 0)


class Drawings(unittest.TestCase):
    """The numbers the captions state, read back from the drawn files."""

    def view(self, kind, system=None, view_id=None):
        data = json.loads((DATA / kind / "white_hole.json").read_text(encoding="utf-8"))
        if kind == "diagrams":
            return next(v for v in data["systems"][system] if v["id"] == view_id)
        return data["views"][0]

    def marker(self, view, kind, legend=None):
        return next(m for m in view["markers"] if m["kind"] == kind and (legend is None or legend in m.get("legend", "")))

    def chart_point(self, view, unit):
        """A point of a flat view's unit square in its display coordinates."""
        X0, X1, Y0, Y1 = view["box"]
        return X0 + unit[0] * (X1 - X0), Y0 + unit[1] * (Y1 - Y0)

    def near(self, line, x=None, y=None):
        """The point of a polyline, in display coordinates, at the given x or y, by interpolation."""
        best = None
        for (xa, ya), (xb, yb) in zip(line, line[1:]):
            lo, hi, value = (xa, xb, x) if x is not None else (ya, yb, y)
            if (lo - value) * (hi - value) <= 0 and lo != hi:
                s = (value - lo) / (hi - lo)
                best = (xa + s * (xb - xa), ya + s * (yb - ya))
        return best

    def test_the_outgoing_chart_draws_the_surface_through_the_crossing_and_the_rest(self):
        view = self.view("diagrams", "exterior_eddington_finkelstein", "finkelstein")
        line = [self.chart_point(view, p) for p in self.marker(view, "surface")["lines"][0]]
        # X = r, Y = u + r: the crossing r = r_s at u = 0, on the way out, and the rest r = 2 r_s at u = pi + ln 2.
        rising = line[:max(range(len(line)), key=lambda k: line[k][0]) + 1]
        crossing = self.near(rising, x=1.0)
        self.assertAlmostEqual(crossing[1], 1.0, delta=0.03)
        rest = max(line, key=lambda p: p[0])
        self.assertAlmostEqual(rest[0], 2.0, delta=0.01)
        self.assertAlmostEqual(rest[1], 2.0 + U_REST, delta=0.15)
        self.assertAlmostEqual(line[0][1] - line[0][0], U_BANG, delta=0.05)

    def test_kruskals_chart_draws_the_surface_out_of_the_singularity_and_through_the_past_horizon(self):
        view = self.view("diagrams", "exterior_kruskal", "kruskal")
        surface = [self.chart_point(view, p) for p in self.marker(view, "surface")["lines"][0]]
        horizon = [self.chart_point(view, p) for p in self.marker(view, "event")["lines"][0]]
        # Drawn with X = (V - U)/2 and Y = (U + V)/2: the past horizon V = 0 is Y = -X, and it ends
        # on the surface at U = -1, the point (1/2, -1/2).
        for X, Y in horizon:
            self.assertAlmostEqual(X + Y, 0.0, delta=0.01)
        end = min(horizon, key=lambda p: p[0])
        self.assertAlmostEqual(end[0], 0.5, delta=0.02)
        # The surface leaves the singularity UV = 1 at U = -e^(pi/2 - 1)/sqrt(2).
        X, Y = surface[0]
        U, V = Y - X, X + Y
        self.assertAlmostEqual(U, -math.exp(math.pi / 2 - 1) / math.sqrt(2), delta=0.02)
        self.assertAlmostEqual(U * V, 1.0, delta=0.05)

    def test_novikovs_vacuole_draws_the_schwarzschild_sphere_and_the_delayed_singularity(self):
        view = self.view("diagrams", "novikov_comoving", "vacuole")
        sphere = [self.chart_point(view, p) for p in self.marker(view, "apparent")["lines"][0]]
        # R = F: ct = 3 + (2/3) r^3 in the core, r_2 - r + 2/3 in the vacuum, (2/3)(r/4)^3 beyond.
        for r, ct in ((1.0, 3 + 2 / 3), (2.5, 1.5 + 2 / 3), (4.0, 2 / 3), (5.0, (2 / 3) * (5 / 4) ** 3)):
            self.assertAlmostEqual(self.near(sphere, x=r)[1], ct, delta=0.06)
        bang = [self.chart_point(view, p) for p in self.marker(view, "singular")["lines"][0]]
        for r, ct in ((0.5, 3.0), (2.5, 1.5), (3.5, 0.5)):
            self.assertAlmostEqual(self.near(bang, x=r)[1], ct, delta=0.05)

    def test_the_conformal_diagram_draws_both_horizons_from_the_centre(self):
        view = self.view("conformal")
        lines = {layer["class"]: layer["points"] for layer in view["layers"] if layer.get("kind") == "line"}
        past, event = lines["horizon"], lines["event"]
        self.assertAlmostEqual(past[0][0], 0.0, places=3)
        self.assertAlmostEqual(past[0][1], -(math.pi - 3 * CHI0), places=3)
        self.assertAlmostEqual(past[-1][0], 3 * CHI0, places=3)
        self.assertAlmostEqual(past[-1][1], -math.pi, places=3)
        # The event horizon is the past horizon mirrored in the moment of rest, T = 0.
        for (xa, ta), (xb, tb) in zip(past, event):
            self.assertAlmostEqual(xa, xb, places=3)
            self.assertAlmostEqual(ta, -tb, places=3)

    def test_the_embedding_runs_from_the_first_moment_to_the_rest(self):
        view = self.view("embedding")
        am = 2 / math.sin(CHI0)
        times = [surface["time"] for surface in view["surfaces"]]
        self.assertEqual(times, sorted(times))
        self.assertAlmostEqual(times[-1], math.pi * am / 2, places=6)
        for surface, fraction in zip(view["surfaces"], (0.2, 0.4, 0.7, 1.0)):
            eta = fraction * math.pi
            self.assertAlmostEqual(surface["time"], (am / 2) * (eta - math.sin(eta)), places=6)
            # The dust's cap reaches chi_0 and the areal radius of its rim is a sin chi_0.
            dust = next(p for p in surface["pieces"] if p["id"] == "dust")
            self.assertAlmostEqual(dust["points"][-1][0], CHI0, places=6)


if __name__ == "__main__":
    unittest.main()
