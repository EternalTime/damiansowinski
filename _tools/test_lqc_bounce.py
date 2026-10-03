"""The bounce of loop quantum cosmology, the flat universe of a massless scalar field in its
effective dynamics: python3 -m unittest discover -s _tools

What its texts state, held on the published files: the scale factor a = (1 + t^2/t_b^2)^(1/6)
solves the modified Friedmann equation H^2 = (8 pi G/3) rho (1 - rho/rho_c) with rho = rho_c/a^6 and
24 pi G rho_c t_b^2 = 1; the Ricci and Kretschmann scalars are finite, the Kretschmann scalar
greatest at the bounce, 4/3 c^4 t_b^4; the rate of expansion is greatest at t = t_b, where the
density is half the critical density; the harmonic chart is the cosmic one pulled back through
t = t_b sinh(tau/t_b); and the drawings hold what their captions say: the rays of the spacetime
diagrams keep x -+ eta, the observers of the embedding stand at radius a r, and the Hubble sphere
comes nearest the centre at t = sqrt(3/2) t_b. The tests of the published mathematics need sympy
and are skipped where it is absent; the tests of the drawings read the files alone.
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


def published():
    return json.loads((DATA / "metrics" / "lqc_bounce.json").read_text(encoding="utf-8"))


def drawing(kind):
    return json.loads((DATA / kind / "lqc_bounce.json").read_text(encoding="utf-8"))


def scale(t):
    """The scale factor at the proper time t, in units of t_b."""
    return (1 + t * t) ** (1 / 6)


def eta(t, n=4000):
    """The conformal time int_0^t dt/a in units of t_b, by Simpson's rule, which no drawing uses."""
    h = t / n
    total = 1 / scale(0) + 1 / scale(t)
    for k in range(1, n):
        total += (4 if k % 2 else 2) / scale(k * h)
    return total * h / 3


@unittest.skipUnless(HAS_SYMPY, "sympy is not installed")
class Published(unittest.TestCase):
    def setUp(self):
        import sympy
        import verify_metrics
        self.sp, self.vm = sympy, verify_metrics

    def chart(self, system):
        entry = next(c for c in published()["coordinates"] if c["id"] == system)
        reader = self.vm.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]], ())
        names = entry["coords"]
        g = self.sp.zeros(4, 4)
        for e in entry["metric_components"]:
            i, j = (names.index(k) for k in e["indices"])
            g[i, j] = reader(e["value"])
        return g, [reader.symbol[n] for n in names], reader.parameters["t_b"], entry, reader

    def scalar(self, entry, reader, field):
        return reader(entry[field].partition("=")[2])

    def test_the_scale_factor_solves_the_modified_friedmann_equation(self):
        sp = self.sp
        for system in ("cosmic", "comoving_spherical"):
            g, X, tb, _, reader = self.chart(system)
            t = X[0]
            a = sp.sqrt(g[1, 1])
            H = sp.diff(a, t) / a
            # (8 pi G/3) rho_c = 1/9 t_b^2 and rho/rho_c = a^-6.
            self.assertEqual(self.vm.norm(H ** 2 - a ** -6 * (1 - a ** -6) / (9 * tb ** 2)), 0, system)
            self.assertEqual(self.vm.norm(g[0, 0] + 1), 0, system)
            self.assertEqual(self.vm.norm(a.subs(t, 0) - 1), 0, system)

    def test_curvature_is_finite_and_greatest_at_the_bounce(self):
        sp = self.sp
        g, X, tb, entry, reader = self.chart("cosmic")
        t, c = X[0], reader.c
        a = sp.sqrt(g[1, 1])
        R = self.scalar(entry, reader, "ricci_scalar")
        # A flat Friedmann universe has R = 6(a'' / a + (a'/a)^2)/c^2.
        expected = 6 * (sp.diff(a, t, 2) / a + (sp.diff(a, t) / a) ** 2) / c ** 2
        self.assertEqual(self.vm.norm(R - expected), 0)
        K = sp.lambdify(t, self.scalar(entry, reader, "kretschmann").subs({c: 1, tb: 1}))
        self.assertAlmostEqual(K(0.0), 4 / 3, places=12)
        self.assertTrue(all(0 < K(0.01 * k) <= 4 / 3 for k in range(-1000, 1001)))

    def test_the_expansion_rate_is_greatest_at_half_the_critical_density(self):
        H = [t / (3 * (1 + t * t)) for t in (0.999, 1.0, 1.001)]
        self.assertGreater(H[1], H[0])
        self.assertGreater(H[1], H[2])
        self.assertAlmostEqual(scale(1.0) ** -6, 0.5)

    def test_the_harmonic_chart_is_the_cosmic_one_pulled_back(self):
        sp = self.sp
        gc, Xc, tbc, _, _ = self.chart("cosmic")
        gh, Xh, tbh, entry, reader = self.chart("harmonic")
        tau = Xh[0]
        t_of = tbh * sp.sinh(tau / tbh)
        image = {Xc[0]: t_of, tbc: tbh}
        self.assertEqual(self.vm.norm(gh[1, 1] - gc[1, 1].subs(image)), 0)
        self.assertEqual(self.vm.norm(gh[0, 0] + sp.diff(t_of, tau) ** 2), 0)
        # The harmonic time solves the wave equation: d_tau(sqrt(-g) g^tautau) = 0.
        self.assertEqual(self.vm.norm(sp.diff(sp.sqrt(-gh.det()) / gh[0, 0], tau)), 0)


class Drawings(unittest.TestCase):
    def test_every_ray_keeps_x_plus_or_minus_its_conformal_time(self):
        views = drawing("diagrams")["systems"]
        for system, view, time in (("cosmic", "tx", lambda y: y), ("harmonic", "taux", math.sinh)):
            v = next(w for w in views[system] if w["id"] == view)
            x0, x1, y0, y1 = v["box"]
            for family, sign in (("P", 1), ("M", -1)):
                for ray in v["rays"][family]:
                    values = [(x0 + (x1 - x0) * X) + sign * eta(time(y0 + (y1 - y0) * Y))
                              for X, Y in ray if 0 <= X <= 1 and 0 <= Y <= 1]
                    if len(values) > 1:
                        self.assertLess(max(values) - min(values), 0.02, f"{system} {family} {ray[:2]}")

    def test_each_comoving_observer_stands_at_radius_a_r(self):
        view = drawing("embedding")["views"][0]
        for surface in view["surfaces"]:
            a = scale(surface["time"])
            for r, rho, z in surface["pieces"][0]["points"]:
                self.assertAlmostEqual(rho, a * r, delta=1e-4)
                self.assertAlmostEqual(z, 0.0, delta=1e-9)

    def test_the_hubble_sphere_comes_nearest_the_centre_at_t_root_three_halves(self):
        # r = 3 c (t_b^2 + t^2)/(a |t|) in units of c t_b, the caption's 5.26 at sqrt(3/2).
        r = [3 * (1 + t * t) / (scale(t) * t) for t in (1.2, math.sqrt(1.5), 1.25)]
        self.assertLess(r[1], r[0])
        self.assertLess(r[1], r[2])
        self.assertAlmostEqual(r[1], 5.26, delta=0.005)

    def test_the_conformal_time_runs_over_the_whole_line(self):
        # eta grows as (3/2) t^(2/3), so the universe fills Minkowski's whole half diamond.
        self.assertAlmostEqual(eta(1e4, 200000) / (1.5 * 1e4 ** (2 / 3)), 1, delta=0.01)


if __name__ == "__main__":
    unittest.main()
