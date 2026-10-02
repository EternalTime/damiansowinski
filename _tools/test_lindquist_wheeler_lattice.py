"""The lattice universe of Lindquist and Wheeler:
python3 -m unittest discover -s _tools

What its texts state, held on the published files: the angular radius of a cell of equal volume
for each of the six regular tilings, and the lattice's largest radius against Friedmann's, 1.43
times at five masses and 1.2 percent apart at 600; the comoving chart a vacuum with
Schwarzschild's Kretschmann scalar on the shells' cycloids; the comparison hypersphere dust at
rest of the density the lattice's equation states; the cosmological time chart's comoving shells
unit normals of its slices; and the drawings: every frame of the embedding's movie tangent where
a cell meets the hypersphere, Flamm's paraboloid at the moment of rest, the boundary's turn and
its proper time on the spacetime and conformal diagrams. The tests of the published mathematics
need sympy and numpy and are skipped where they are absent; the tests of the drawings read the
files alone.
"""
import importlib.util
import json
import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "MFS" / "assets" / "data"
NAME = "lindquist_wheeler_lattice"
HAS_SYMPY = all(importlib.util.find_spec(m) is not None for m in ("sympy", "numpy", "scipy", "contourpy"))
sys.path.insert(0, str(Path(__file__).resolve().parent / "derivations"))


def psi_of(cells):
    """The angular radius of one of N cells of equal volume on the three sphere."""
    lo, hi = 0.0, math.pi / 2
    for _ in range(200):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if cells * (2 * mid - math.sin(2 * mid)) < 2 * math.pi else (lo, mid)
    return (lo + hi) / 2


PSI = psi_of(8)
RM = 1 / math.sin(PSI) ** 2
AM = 1 / math.sin(PSI) ** 3


def published(kind="metrics"):
    return json.loads((DATA / kind / f"{NAME}.json").read_text(encoding="utf-8"))


class Lattice(unittest.TestCase):
    def test_the_ratio_to_friedmann_is_what_the_history_states(self):
        ratio = {n: 3 * math.pi / (2 * n * math.sin(psi_of(n)) ** 3) for n in (5, 8, 16, 24, 120, 600)}
        self.assertAlmostEqual(ratio[5], 1.43, places=2)
        self.assertAlmostEqual(100 * (ratio[600] - 1), 1.2, places=1)
        values = list(ratio.values())
        self.assertTrue(all(a > b > 1 for a, b in zip(values, values[1:])))
        # As N grows the cell is a small ball, N (4/3) psi^3 = 2 pi, and the ratio tends to 1.
        self.assertAlmostEqual(3 * math.pi / (2 * 10 ** 9 * math.sin(psi_of(10 ** 9)) ** 3), 1.0, places=5)

    def test_the_boundary_is_tangent_to_the_hypersphere(self):
        self.assertAlmostEqual(math.sqrt(1 - 1 / RM), math.cos(PSI), places=14)
        self.assertAlmostEqual(AM * math.sin(PSI), RM, places=14)
        # The boundary's proper time from rest to r = 0 is the hypersphere's.
        self.assertAlmostEqual(math.pi / 2 * RM ** 1.5, math.pi / 2 * AM, places=13)

    def test_the_history_keeps_its_shape_and_its_sources(self):
        metric = published()
        history = metric["history"]
        self.assertLessEqual(history.count("?"), 3)
        self.assertEqual(len(history.split("¶")), 5)
        for key in ("lindquist1957", "lindquist2011", "clifton2009", "clifton2012", "bentivegna2012", "wheeler1983"):
            self.assertIn(key, metric["references"])


class Drawings(unittest.TestCase):
    def test_every_frame_of_the_movie_is_tangent_at_the_boundary(self):
        view = published("embedding")["views"][0]
        frames = view["movie"]["frames"]
        self.assertEqual(len(view["surfaces"]), 4)
        self.assertGreater(len(frames), 40)
        for frame in frames:
            cell, sphere, far = frame["pieces"]
            self.assertEqual((cell["system"], sphere["system"]), ("lindquist_wheeler", "comparison_hypersphere"))
            # One point: the boundary shell r_m of the cell is the circle chi = psi of the hypersphere.
            self.assertAlmostEqual(cell["points"][-1][0], RM, places=12)
            self.assertAlmostEqual(sphere["points"][0][0], PSI, places=12)
            for a, b in ((cell["points"][-1], sphere["points"][0]), (far["points"][-1], sphere["points"][-1])):
                self.assertAlmostEqual(a[1], b[1], places=6)
                self.assertAlmostEqual(a[2], b[2], places=6)
            # The hypersphere is a sphere about the origin, of the radius its rim fixes.
            a = sphere["points"][0][1] / math.sin(PSI)
            for chi, rho, z in sphere["points"]:
                self.assertAlmostEqual(rho, a * math.sin(chi), places=5)
                self.assertAlmostEqual(z, -a * math.cos(chi), places=5)
            # One tangent: the last chord of the cell has dr/dl = cos(psi), to the length of a chord.
            (_, r1, z1), (_, r2, z2) = cell["points"][-2:]
            self.assertAlmostEqual((r2 - r1) / math.hypot(r2 - r1, z2 - z1), math.cos(PSI), delta=0.05)

    def test_the_moment_of_rest_is_flamms_paraboloid(self):
        first = published("embedding")["views"][0]["surfaces"][0]
        cell = first["pieces"][0]
        self.assertEqual(cell["start"]["kind"], "throat")
        z0 = cell["points"][0][2]
        for rho, r, z in cell["points"]:
            self.assertAlmostEqual(r, rho, places=6)
            self.assertAlmostEqual(z - z0, 2 * math.sqrt(max(rho - 1, 0)), places=4)

    def test_the_funnel_ends_in_a_point_once_the_throat_has_gone(self):
        for surface in published("embedding")["views"][0]["surfaces"][1:]:
            cell = surface["pieces"][0]
            tau = surface["time"]
            self.assertGreater(tau, math.pi / 2)
            self.assertEqual(cell["start"]["kind"], "apex")
            # The movie keeps a moment's time to six places.
            self.assertAlmostEqual(cell["points"][0][0], (2 * tau / math.pi) ** (2 / 3), places=6)
            self.assertAlmostEqual(cell["points"][0][1], 0.0, places=6)

    def test_the_boundary_turns_at_r_m_on_the_spacetime_diagrams(self):
        systems = published("diagrams")["systems"]
        cell = systems["schwarzschild_cell"][0]
        (X0, X1, Y0, Y1) = cell["box"]
        line = next(m for m in cell["markers"] if m["kind"] == "surface")["lines"][0]
        turn = max(line, key=lambda p: p[0])
        self.assertAlmostEqual(X0 + turn[0] * (X1 - X0), RM, places=3)
        self.assertAlmostEqual(Y0 + turn[1] * (Y1 - Y0), 0.0, places=2)
        # In the cosmological time the boundary leaves r = 0 at tau = 0 and is at rest after (pi/2) r_m^(3/2).
        cosmo = systems["cosmological_time"][0]
        (X0, X1, Y0, Y1) = cosmo["box"]
        line = next(m for m in cosmo["markers"] if m["kind"] == "surface")["lines"][0]
        self.assertEqual(line[0], [0.0, 0.0])
        self.assertAlmostEqual(X0 + line[-1][0] * (X1 - X0), RM, places=3)
        self.assertAlmostEqual(Y0 + line[-1][1] * (Y1 - Y0), math.pi / 2 * RM ** 1.5, places=2)
        self.assertIn("left", next(m for m in cosmo["markers"] if m["kind"] == "singular")["edges"])

    def test_the_comoving_view_closes_on_a_bang_and_a_crunch(self):
        view = published("diagrams")["systems"]["lindquist_wheeler"][0]
        (X0, X1, Y0, Y1) = view["box"]
        curves = next(m for m in view["markers"] if m["kind"] == "singular")["lines"]
        self.assertEqual(len(curves), 2)
        for curve in curves:
            for u, v in curve:
                rho, tau = X0 + u * (X1 - X0), Y0 + v * (Y1 - Y0)
                if abs(tau) < Y1 - 0.05:
                    self.assertAlmostEqual(abs(tau), math.pi / 2 * rho ** 1.5, places=2)
        self.assertEqual(sorted(1 if c[0][1] > 0.5 else -1 for c in curves), [-1, 1])

    def test_the_conformal_cell_lies_between_the_throat_and_the_boundary(self):
        views = {v["id"]: v for v in published("conformal")["views"]}
        self.assertEqual(sorted(views), ["cell", "expanding", "hypersphere", "shells"])
        rest = 2 * math.atan(math.sqrt((RM - 1) * math.exp(RM)))
        for vid in ("cell", "expanding", "shells"):
            layers = views[vid]["layers"]
            edge = next(l for l in layers if l["class"] == "surface")["points"]
            self.assertAlmostEqual(max(p[0] for p in edge), rest, places=3)
            self.assertAlmostEqual(edge[0][1], -math.pi / 2, places=3)
            self.assertAlmostEqual(edge[-1][1], math.pi / 2, places=3)
            throat = next(l for l in layers if l["class"] == "throat")["points"]
            self.assertTrue(all(p[0] == 0 for p in throat))
            # The boundary crosses the horizon U = 0, the line T = X, at eta = pi - 2 psi.
            horizon = [l["points"] for l in layers if l["class"] == "horizon"]
            eta = math.pi - 2 * PSI
            k = math.sqrt(RM - 1)
            V = (k * math.cos(eta / 2) + math.sin(eta / 2)) * math.exp((1 + k * (eta + RM / 2 * (eta + math.sin(eta)))) / 2)
            self.assertTrue(any(abs(h[1][0] - math.atan(V)) < 1e-3 and abs(h[1][1] - math.atan(V)) < 1e-3 for h in horizon))
        # The moments embedded lie on the cell and on the hypersphere, and none on the expanding chart.
        self.assertEqual([len(views[v].get("slices", [])) for v in ("cell", "expanding", "shells", "hypersphere")], [4, 0, 4, 4])


@unittest.skipUnless(HAS_SYMPY, "sympy, numpy, scipy and contourpy are not all installed")
class Published(unittest.TestCase):
    def setUp(self):
        import numpy
        import sympy
        import null_rays
        import verify_metrics
        self.np, self.sp, self.nr, self.vm = numpy, sympy, null_rays, verify_metrics

    def chart(self, system):
        entry = next(c for c in published()["coordinates"] if c["id"] == system)
        return entry, self.vm.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]], ())

    def test_the_cycloid_solver_inverts_eta_plus_sin_eta(self):
        np = self.np
        want = np.linspace(-math.pi, math.pi, 2001)
        eta = np.array([self.nr._lw_eta_one(w) for w in want])
        self.assertLess(float(np.max(np.abs(eta + np.sin(eta) - want))), 1e-14)
        self.assertTrue(math.isnan(self.nr._lw_eta_one(3.5)))
        self.assertEqual(self.nr._lw_eta_one(math.pi), math.pi)

    def test_the_comoving_chart_is_schwarzschilds_vacuum_on_the_cycloids(self):
        np, sp = self.np, self.sp
        entry, reader = self.chart("lindquist_wheeler")
        tau, rho = reader.symbol["\\tau"], reader.symbol["\\rho"]
        r = reader.parameters["r"]
        rng = np.random.default_rng(0)
        rhos = rng.uniform(1.05, RM, 40)
        taus = rng.uniform(-0.9, 0.9, 40) * math.pi / 2 * rhos ** 1.5
        eta = self.nr._lw_eta(taus, rhos)
        names = {}

        def numbers(text):
            e = reader(text).subs({reader.c: 1, reader.parameters["r_s"]: 1})
            for d in sorted(e.atoms(sp.Derivative), key=lambda d: -len(d.variables)):
                i, j = d.variables.count(tau), d.variables.count(rho)
                names.setdefault((i, j), sp.Symbol(f"r{i}{j}"))
                e = e.subs(d, names[(i, j)])
            names.setdefault((0, 0), sp.Symbol("r00"))
            e = e.subs(r, names[(0, 0)])
            keys = sorted(names)
            f = sp.lambdify([rho] + [names[k] for k in keys], e, "numpy")
            return f(rhos, *[self.nr._lw_jet(i, j)(eta, rhos) for i, j in keys])
        radius = self.nr._lw_jet(0, 0)(eta, rhos)
        for component in entry["einstein_tensor"]["variants"]["ul"]["nonzero"]:
            self.assertLess(float(np.max(np.abs(numbers(component["value"])) * radius ** 3)), 1e-8, component["indices"])
        K = numbers(entry["kretschmann"].partition("=")[2])
        self.assertLess(float(np.max(np.abs(K * radius ** 6 / 12 - 1))), 1e-8)

    def test_the_hypersphere_is_dust_of_the_lattices_density(self):
        np, sp = self.np, self.sp
        entry, reader = self.chart("comparison_hypersphere")
        tau, a = reader.symbol["\\tau"], reader.parameters["a"]
        e = np.linspace(-2.8, 2.8, 41)
        values = {"a": AM * (1 + np.cos(e)) / 2, "a1": -np.tan(e / 2), "a2": -1 / (2 * AM * np.cos(e / 2) ** 4)}

        def numbers(text):
            x = reader(text).subs({reader.c: 1, reader.parameters["r_s"]: 1, reader.parameters["psi"]: PSI})
            x = x.subs(sp.Derivative(a, (tau, 2)), sp.Symbol("a2")).subs(sp.Derivative(a, tau), sp.Symbol("a1")).subs(a, sp.Symbol("a"))
            x = x.subs({reader.symbol["\\chi"]: 0.7, reader.symbol["\\theta"]: 1.1})
            return sp.lambdify(sp.symbols("a a1 a2"), x, "numpy")(values["a"], values["a1"], values["a2"]) * np.ones_like(e)
        mixed = {tuple(c["indices"]): c["value"] for c in entry["einstein_tensor"]["variants"]["ul"]["nonzero"]}
        # (da/d(c tau))^2 = r_s/(a sin^3 psi) - 1 on the cycloid, and G^tau_tau = -3 r_s/(a^3 sin^3 psi).
        self.assertLess(float(np.max(np.abs(values["a1"] ** 2 - (AM / values["a"] - 1)))), 1e-12)
        density = numbers(mixed[("\\tau", "\\tau")])
        self.assertLess(float(np.max(np.abs(density * values["a"] ** 3 / (-3 * AM) - 1))), 1e-10)
        for index in ("\\chi", "\\theta", "\\phi"):
            self.assertLess(float(np.max(np.abs(numbers(mixed[(index, index)])) * values["a"] ** 2)), 1e-9)

    def test_the_comoving_shells_are_normal_to_the_slices_of_cosmological_time(self):
        sp = self.sp
        entry, reader = self.chart("cosmological_time")
        r, rs, E = reader.symbol["r"], reader.parameters["r_s"], reader.parameters["E"]
        g = {tuple(c["indices"]): reader(c["value"]).subs(reader.c, 1) for c in entry["metric_components"]}
        W = sp.sqrt(E - 1 + rs / r)
        at = {r: sp.Rational(7, 5), rs: 1, E: sp.Rational(2, 5)}
        low_tau = (g[("\\tau", "\\tau")] + g[("\\tau", "r")] * W).subs(at)
        low_r = (g[("\\tau", "r")] + g[("r", "r")] * W).subs(at)
        self.assertAlmostEqual(float(low_tau), -1.0, places=12)
        self.assertAlmostEqual(float(low_r), 0.0, places=12)
        # The root vanishes at r_s/(1 - E), the boundary's turn: r_m at E = cos^2(psi).
        self.assertAlmostEqual(1 / (1 - math.cos(PSI) ** 2), RM, places=12)

    def test_the_boundary_is_the_cycloid_in_both_cell_charts(self):
        np, nr = self.np, self.nr
        specs = {(s.metric, s.system, s.view): s for s in nr.DIAGRAMS}
        bang = nr.Chart(specs[(NAME, "cosmological_time", "radial")]).surface
        eta = np.linspace(0.02, math.pi - 0.05, 60)
        r = RM / 2 * (1 - np.cos(eta))
        tau = RM ** 1.5 / 2 * (eta - np.sin(eta))
        # To the straight lines between the traced points, a twentieth of a pixel of the drawing.
        self.assertLess(float(np.max(np.abs(bang(tau) - r))), 1e-4)
        self.assertTrue(np.isnan(bang(np.array([-0.1, 3.6]))).all())
        rest = nr.Chart(specs[(NAME, "schwarzschild_cell", "radial")]).surface
        self.assertAlmostEqual(float(rest(np.array([0.0]))[0]), RM, places=12)
        t = np.array([0.5, 1.5, 3.0])
        self.assertLess(float(np.max(np.abs(rest(t) - rest(-t)))), 1e-15)
        # Novikov's t of the boundary shell at the cycloid parameter 0.6, against the traced geodesic.
        k, e = math.sqrt(RM - 1), 0.6
        novikov = math.log(abs((k + math.tan(e / 2)) / (k - math.tan(e / 2)))) + k * (e + RM / 2 * (e + math.sin(e)))
        self.assertAlmostEqual(float(rest(np.array([novikov]))[0]), RM * math.cos(e / 2) ** 2, places=5)


if __name__ == "__main__":
    unittest.main()
