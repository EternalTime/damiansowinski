"""The point charge of Born and Infeld's electrodynamics with its own gravity, Hoffmann's
solution of 1935: .venv.noindex/bin/python -m unittest discover -s _tools

What its texts and drawings state, held on the published files. The mass function, which the
charts write with an incomplete elliptic integral of the first kind, is computed again here as
the whole mass less the energy of the field outside each radius, by quadrature of its slope, and
the numbers the History and the captions quote are computed from it: Born and Infeld's 1.2360,
the mass of Hoffmann's particle, the cone at its centre, the black hole's horizon and surface
gravity, and the cone of the electron. The embedding diagram's surfaces are held to the
quadrature of their slopes, from the numbers written and nothing else, and the other drawings to
the horizon and the singularities they mark. Where sympy is installed, the published metric of
every chart is held to Einstein's equations with the field of a point charge of Born and Infeld's
Lagrangian for its source, and the reader's elliptic integral to its integrand; those tests are
skipped where it is absent.
"""
import importlib.util
import json
import math
import sys
import types
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "MFS" / "assets" / "data"
HAS_SYMPY = importlib.util.find_spec("sympy") is not None
sys.path.insert(0, str(Path(__file__).resolve().parent / "derivations"))

# Born and Infeld's number, the energy of the field of a point charge in units of Q^2/r_0.
NUMBER = math.gamma(0.25) ** 2 / (6 * math.sqrt(math.pi))
# The charge as its diagrams draw it, in units of r_0: r_q = r_0/2, Hoffmann's particle, whose
# whole mass is the energy of its field, and the black hole.
RQ = 0.5
PARTICLE = 2 * NUMBER * RQ ** 2
HOLE = 2.0
HORIZON = 1.8666065519401185


def simpson(of, a, b, n=2000):
    h = (b - a) / n
    return h / 3 * (of(a) + of(b) + sum((4 if k % 2 else 2) * of(a + k * h) for k in range(1, n)))


def elliptic_f(phi, m):
    """The incomplete elliptic integral of the first kind, F(phi | m), by Simpson's rule."""
    return simpson(lambda t: 1 / math.sqrt(1 - m * math.sin(t) ** 2), 0.0, phi)


def mass(r, rs, rq=RQ, r0=1.0):
    """The mass function as the charts define it."""
    W = math.sqrt(r ** 4 + r0 ** 4)
    return rs / 2 + rq ** 2 * r / (3 * (r * r + W)) - rq ** 2 / (3 * r0) * elliptic_f(2 * math.atan2(r0, r), 0.5)


def slope(r, rq=RQ, r0=1.0):
    """dm/dr, the energy of the field in a shell."""
    return rq ** 2 / (r * r + math.sqrt(r ** 4 + r0 ** 4))


def energy_outside(r, rq=RQ):
    """The energy of the field outside r, at r_0 = 1: the quadrature of the slope, in 1/x beyond x = 1."""
    def in_w(w):
        return rq ** 2 / (1 + math.sqrt(1 + w ** 4))      # slope(1/w)/w^2
    if r >= 1:
        return simpson(in_w, 0.0, 1 / r)
    return simpson(lambda x: slope(x, rq), r, 1.0) + simpson(in_w, 0.0, 1.0)


def f(r, rs):
    return 1 - 2 * mass(r, rs) / r


def bisect(of, lo, hi):
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if of(lo) * of(mid) > 0 else (lo, mid)
    return 0.5 * (lo + hi)


class Numbers(unittest.TestCase):
    def test_born_and_infelds_number_is_the_energy_of_the_field(self):
        # (2/3) int_0^infinity dy/sqrt(y^4 + 1) = Gamma(1/4)^2/(6 sqrt(pi)) = 1.2360, which their (8.6) prints as 1.2360.
        integral = simpson(lambda y: 1 / math.sqrt(y ** 4 + 1), 0, 1) + simpson(lambda w: 1 / math.sqrt(1 + w ** 4), 0, 1)
        self.assertAlmostEqual(2 * integral / 3, NUMBER, places=10)
        self.assertEqual(round(NUMBER, 4), 1.2360)
        self.assertEqual(round(NUMBER, 3), 1.236)
        self.assertAlmostEqual(energy_outside(0.0, 1.0), NUMBER, places=10)

    def test_the_mass_function_is_the_whole_mass_less_the_energy_of_the_field_outside(self):
        for rs in (PARTICLE, HOLE, 0.7):
            for r in (0.01, 0.3, 1.0, 1.9, 5.0, 40.0):
                self.assertAlmostEqual(mass(r, rs), rs / 2 - energy_outside(r), places=9, msg=f"at {r}")

    def test_the_mass_function_far_away_and_at_the_centre(self):
        # r_s/2 - r_q^2/(2r) far away, Reissner and Nordstrom's, and r_s/2 - 1.2360 r_q^2/r_0 at the centre.
        self.assertAlmostEqual((mass(1e3, HOLE) - HOLE / 2) * 1e3, -RQ ** 2 / 2, places=7)
        self.assertAlmostEqual(mass(1e-9, HOLE), HOLE / 2 - NUMBER * RQ ** 2, places=8)
        self.assertAlmostEqual(mass(1e-9, 1.0, 0.3, 0.7), 0.5 - NUMBER * 0.09 / 0.7, places=8)

    def test_maxwells_limit_is_reissner_and_nordstroms_metric(self):
        # As r_0 -> 0 the metric function falls short of theirs by r_q^2 r_0^4/(20 r^6).
        r0, r, rs = 0.05, 1.3, HOLE
        ours = 1 - 2 * mass(r, rs, RQ, r0) / r
        theirs = 1 - rs / r + RQ ** 2 / r ** 2
        self.assertAlmostEqual((ours - theirs) / (-RQ ** 2 * r0 ** 4 / (20 * r ** 6)), 1.0, places=3)

    def test_hoffmanns_particle(self):
        # Its whole mass is the field's, r_s = 0.618 r_0 at r_q = r_0/2, the mass function vanishes at
        # the centre as r_q^2 r/r_0^2, g^rr tends to 1 - 2 r_q^2/r_0^2 = 1/2 there and is positive at
        # every radius, and a small circle about the centre is sqrt(1/2) of 2 pi times its radius round.
        self.assertEqual(round(PARTICLE, 3), 0.618)
        self.assertAlmostEqual(PARTICLE, float("0.61802489243379063947795011573"), places=14)
        self.assertAlmostEqual(mass(1e-9, PARTICLE), 0.0, places=8)
        self.assertAlmostEqual(f(1e-3, PARTICLE), 0.5, places=6)
        self.assertTrue(all(f(r, PARTICLE) > 0 for r in [0.001 * 1.02 ** k for k in range(700)]))
        proper = simpson(lambda x: 1 / math.sqrt(f(x, PARTICLE)), 1e-9, 1e-3, 20)
        self.assertAlmostEqual(1e-3 / proper, math.sqrt(1 - 2 * RQ ** 2), places=5)
        self.assertAlmostEqual(proper / 1e-3, math.sqrt(2), places=5)
        # The outgoing rays of the ingoing chart, dv/dr = 2/f: 4 at the centre and 2 far away.
        self.assertAlmostEqual(2 / f(1e-3, PARTICLE), 4.0, places=4)
        self.assertAlmostEqual(2 / f(1e6, PARTICLE), 2.0, places=5)

    def test_the_particle_has_a_horizon_only_past_two_r_q_squared_equal_to_r_0_squared(self):
        # Rasheed's "one non-degenerate horizon for Q > b/2 and none otherwise".
        def least(rq):
            rs = 2 * NUMBER * rq ** 2
            return min(1 - 2 * mass(r, rs, rq) / r for r in [0.002 * 1.03 ** k for k in range(400)])
        self.assertGreater(least(0.70), 0)
        self.assertGreater(least(1 / math.sqrt(2)), 0)
        self.assertLess(least(0.72), 0)

    def test_the_black_hole_drawn(self):
        # One horizon, at 1.8666 r_0, with the surface gravity 0.2490/r_0; Reissner and Nordstrom's
        # of the same mass and charge has horizons at 1.8660 and 0.1340 r_0; and the field holds
        # 0.309 r_0 in all, less than r_s/2 = r_0, so the mass at the centre is positive.
        rh = bisect(lambda r: f(r, HOLE), 1.5, 2.2)
        self.assertAlmostEqual(rh, HORIZON, places=9)
        self.assertEqual(round(rh, 3), 1.867)
        radii = [0.001 * 1.01 ** k for k in range(1200)]
        self.assertEqual(sum(1 for a, b in zip(radii, radii[1:]) if f(a, HOLE) * f(b, HOLE) < 0), 1)
        kappa = (1 - 2 * slope(rh)) / (2 * rh)
        h = 1e-5
        self.assertAlmostEqual((f(rh + h, HOLE) - f(rh - h, HOLE)) / (4 * h), kappa, places=6)
        self.assertEqual(round(kappa, 3), 0.249)
        self.assertEqual(round(1 + math.sqrt(1 - RQ ** 2), 4), 1.8660)
        self.assertEqual(round(1 - math.sqrt(1 - RQ ** 2), 3), 0.134)
        self.assertEqual(round(rh, 4), 1.8666)
        self.assertEqual(round(NUMBER * RQ ** 2, 3), 0.309)
        self.assertGreater(mass(1e-9, HOLE), 0)

    def test_the_kretschmann_scalar_at_the_centre(self):
        # K = (16/r^6)((m - r_q^2 r/W)^2 + (m - r_q^2 r/(r^2 + W))^2 + m^2): 16 r_q^4/(r_0^4 r^4) at
        # the particle's centre and Schwarzschild's 48 m(0)^2/r^6 at the black hole's.
        def K(r, rs):
            W, m = math.sqrt(r ** 4 + 1), mass(r, rs)
            return 16 / r ** 6 * ((m - RQ ** 2 * r / W) ** 2 + (m - RQ ** 2 * r / (r * r + W)) ** 2 + m ** 2)
        self.assertAlmostEqual(K(1e-3, PARTICLE) * 1e-12 / (16 * RQ ** 4), 1.0, places=5)
        self.assertAlmostEqual(K(1e-4, HOLE) * 1e-24 / (48 * mass(1e-9, HOLE) ** 2), 1.0, places=3)

    def test_the_cone_of_the_electron_is_very_nearly_flat(self):
        # 2 r_q^2/r_0^2 with r_q^2 = G e^2/c^4 in Gaussian units and Born and Infeld's
        # r_0 = 1.236 e^2/(m_0 c^2), their (8.7), 1.236 times the classical radius of the electron.
        G, e, c, m0 = 6.6743e-8, 4.80320471e-10, 2.99792458e10, 9.1093837e-28
        r0 = NUMBER * e ** 2 / (m0 * c ** 2)
        self.assertEqual(round(r0 * 1e13, 2), 3.48)
        deficit = 2 * G * e ** 2 / c ** 4 / r0 ** 2
        self.assertEqual(f"{deficit:.0e}", "3e-43")


class Surfaces(unittest.TestCase):
    def setUp(self):
        self.views = {v["id"]: v for v in json.loads(
            (DATA / "embedding" / "born_infeld_charge.json").read_text(encoding="utf-8"))["views"]}

    def pieces(self, view):
        return {p["id"]: p for p in self.views[view]["surfaces"][0]["pieces"]}

    def test_the_particles_slice_leaves_the_centre_as_a_cone_of_slope_one(self):
        piece = self.pieces("particle")["particle"]
        self.assertEqual(piece["start"]["kind"], "apex")
        points = piece["points"]
        self.assertEqual(points[0], [0, 0, 0])
        r, _, z = points[1]
        self.assertLess(abs(z / r - 1), 2e-2)

        def rise(x):
            m = mass(x, PARTICLE)
            return math.sqrt(2 * m / (x - 2 * m))
        for r, rho, z in points[4::6] + [points[-1]]:
            self.assertLess(abs(rho - r), 2e-6)
            self.assertLess(abs(z - simpson(rise, 1e-9, r, 200)), 2e-5, f"at {r}")

    def test_the_black_holes_slice_stands_at_the_quadrature_of_its_slope(self):
        first, second = slope(HORIZON), -2 * RQ ** 2 * HORIZON / (
            math.sqrt(HORIZON ** 4 + 1) * (HORIZON ** 2 + math.sqrt(HORIZON ** 4 + 1)))

        def in_q(q):
            # dz/dr = sqrt(2m/(r - 2m)) with r = r_h + q^2, where it is finite at the throat; next to
            # the throat r - 2m is its series there, (1 - 2m') q^2 - m'' q^4.
            r = HORIZON + q * q
            m = mass(r, HOLE)
            if q < 1e-2:
                return 2 * math.sqrt(2 * m / (1 - 2 * first - second * q * q))
            return 2 * q * math.sqrt(2 * m / (r - 2 * m))
        for sign, pid in ((1, "exterior"), (-1, "other_exterior")):
            points = self.pieces("hole")[pid]["points"]
            self.assertAlmostEqual(points[0][0], HORIZON, places=12)
            for r, rho, z in points[::9] + [points[-1]]:
                self.assertLess(abs(rho - r), 2e-6)
                self.assertLess(abs(z - sign * simpson(in_q, 0.0, math.sqrt(r - HORIZON), 100)), 4e-5, f"{pid} at {r}")


class Drawings(unittest.TestCase):
    def test_the_spacetime_diagrams_mark_the_horizon_of_the_black_hole_alone(self):
        diagrams = json.loads((DATA / "diagrams" / "born_infeld_charge.json").read_text(encoding="utf-8"))["systems"]
        seen = 0
        for views in diagrams.values():
            for view in views:
                kinds = {m["kind"]: m for m in view["markers"]}
                self.assertEqual(kinds["singular"]["edges"], ["left"])
                if view["id"] == "particle":
                    self.assertNotIn("grr", kinds)
                    self.assertIn("\\approx 0.618", view["settings"])
                else:
                    (line,) = kinds["grr"]["lines"]
                    self.assertTrue(all(abs(x - HORIZON / 4) < 1e-4 for x, _ in line))
                seen += 1
        self.assertEqual(seen, 8)

    def test_the_conformal_diagrams_put_the_singularities_where_the_captions_say(self):
        views = {v["id"]: v for v in json.loads(
            (DATA / "conformal" / "born_infeld_charge.json").read_text(encoding="utf-8"))["views"]}
        self.assertEqual(sorted(views), ["hole", "hole_ingoing", "hole_outgoing",
                                         "particle", "particle_ingoing", "particle_outgoing"])
        for vid, view in views.items():
            singular = [layer["points"] for layer in view["layers"] if layer["class"] == "singular"]
            if vid.startswith("particle"):
                # Timelike, on X = 0, with no horizon drawn.
                self.assertEqual(len(singular), 1)
                self.assertTrue(all(x == 0 for x, _ in singular[0]))
                self.assertFalse([layer for layer in view["layers"] if layer["class"] == "horizon"])
            else:
                # Spacelike, on T = +-pi/2, behind the two null lines of the horizon.
                self.assertEqual(sorted(round(points[0][1], 4) for points in singular),
                                 [-round(math.pi / 2, 4), round(math.pi / 2, 4)])
                for points in singular:
                    self.assertEqual(len({t for _, t in points}), 1)
                self.assertEqual(len([layer for layer in view["layers"] if layer["class"] == "horizon"]), 2)


@unittest.skipUnless(HAS_SYMPY, "sympy is not installed")
class FieldEquations(unittest.TestCase):
    def test_the_readers_elliptic_integral_is_the_integral_of_its_integrand(self):
        import sympy as sp
        import verify_metrics as vm
        reader = vm.Reader(["t", "x"], ["a"], ())
        x = reader.symbol["x"]
        read = reader("\\mathrm{F}\\left(x \\mid \\dfrac{1}{2}\\right)")
        self.assertEqual(read, vm.EllipticF(x, sp.Rational(1, 2)))
        self.assertEqual(sp.simplify(sp.diff(read, x) - 1 / sp.sqrt(1 - sp.sin(x) ** 2 / 2)), 0)
        for phi in (0.3, 1.2, math.pi):
            number = vm.EllipticF(sp.Float(phi), sp.Rational(1, 2))
            self.assertAlmostEqual(float(number.evalf(30)), elliptic_f(phi, 0.5), places=10)
            self.assertAlmostEqual(float(sp.lambdify(x, read, "numpy")(phi)), elliptic_f(phi, 0.5), places=10)

    def test_every_published_chart_solves_einsteins_equations_with_born_and_infelds_field(self):
        import print_charts
        import verify_metrics as vm
        metric = json.loads((DATA / "metrics" / "born_infeld_charge.json").read_text(encoding="utf-8"))
        self.assertEqual([c["id"] for c in metric["coordinates"]], print_charts.BORN_INFELD_CHARTS)
        for entry in metric["coordinates"]:
            key = ("born_infeld_charge", entry["id"])
            # The checker's own reader: it checks the declared slope of the mass function against
            # the elliptic integral that defines it, and raises where the two disagree.
            reader = vm.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]], (),
                               held=vm.HELD[key], rates=vm.RATES[key])
            symbols = [reader.symbol[c] for c in entry["coords"]]
            n = len(symbols)
            published = {tuple(e["indices"]): reader(e["value"]) for e in entry["metric_components"]}
            g = vm.sp.Matrix(n, n, lambda i, j: published.get((entry["coords"][i], entry["coords"][j]), 0))
            geometry = print_charts.cp.Reduced(vm.Geometry(g, symbols, 600), reader.surface)
            chart = types.SimpleNamespace(geo=geometry, reader=reader, symbols=symbols, coords_tex=entry["coords"])
            # Raises where the published metric misses the stress of the field of a point charge
            # in Born and Infeld's electrodynamics, Schwarzschild's metric at r_q = 0, Reissner and
            # Nordstrom's as r_0 -> 0, the mass at the centre, or the cone of Hoffmann's particle.
            print_charts.born_infeld_check(chart, entry["id"])


if __name__ == "__main__":
    unittest.main()
