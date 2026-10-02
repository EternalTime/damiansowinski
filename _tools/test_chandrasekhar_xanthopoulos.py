"""Chandrasekhar and Xanthopoulos's colliding waves, whose region of interaction is Kerr's metric
between its horizons: python3 -m unittest discover -s _tools

What its texts state, held on the published files. The published chart of eta and mu is carried
onto the published Boyer-Lindquist chart by r = m(1 - p eta), cos(theta) = mu, t = x - 2qy/p and
phi = -y/(mp) with a = mq, as the Boyer-Lindquist chart's convention says; with the polarizations
aligned, q = 0, it is the published metric of Schwarzschild inside its horizon; its Kretschmann
scalar is finite on the horizon eta = 1; and an observer at rest on mu = 0 reaches the horizon
after the proper time m(pi/2 - p). The ring's world tube is held to the published metric of the
wave front, to the area pi cos(psi), to being carried by parallel transport, which is that the
rate of its frame has no rotation in it, and to the turn of its long axis. The tests of the
published charts need sympy and are skipped where it is absent; the tube's need nothing.
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
P, Q = 0.6, 0.8      # the p and q the diagrams draw, at alpha = atan(4/3)


@unittest.skipUnless(HAS_SYMPY, "sympy is not installed")
class Published(unittest.TestCase):
    def setUp(self):
        import sympy
        import verify_metrics
        self.sp, self.vm = sympy, verify_metrics

    def chart(self, metric_id, system):
        """(reader, coordinate symbols, g) of one published chart."""
        metric = json.loads((DATA / "metrics" / f"{metric_id}.json").read_text(encoding="utf-8"))
        entry = next(c for c in metric["coordinates"] if c["id"] == system)
        declared = self.vm.DIMENSIONS[(metric_id, system)]
        reader = self.vm.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]],
                                self.vm.time_coordinates(declared, entry["coords"]))
        n = len(entry["coords"])
        g = self.sp.zeros(n, n)
        for component in entry["metric_components"]:
            i, j = (entry["coords"].index(x) for x in component["indices"])
            g[i, j] = reader(component["value"])
        return reader, [reader.symbol[c] for c in entry["coords"]], g, entry

    def points(self):
        sp = self.sp
        for eta, share, alpha in ((sp.Rational(1, 5), sp.Rational(1, 2), sp.Rational(9, 10)),
                                  (sp.Rational(7, 10), sp.Rational(-3, 4), sp.Rational(3, 10)),
                                  (sp.Rational(19, 20), sp.Rational(1, 10), sp.Rational(6, 5))):
            yield eta, eta * share, alpha

    def test_the_chart_of_eta_and_mu_is_kerr_between_its_horizons(self):
        sp = self.sp
        waves, (eta, mu, x, y), g, _ = self.chart("chandrasekhar_xanthopoulos", "prolate")
        kerr, (t, r, theta, phi), k, _ = self.chart("chandrasekhar_xanthopoulos", "boyer_lindquist")
        m, alpha = waves.parameters["m"], waves.parameters["alpha"]
        p, q = sp.cos(alpha), sp.sin(alpha)
        image = [x - 2 * q * y / p, m * (1 - p * eta), sp.acos(mu), -y / (m * p)]
        J = sp.Matrix(4, 4, lambda i, j: sp.diff(image[i], (eta, mu, x, y)[j]))
        at = {r: image[1], theta: image[2], kerr.parameters["m"]: m, kerr.parameters["a"]: m * q}
        difference = J.T * k.subs(at, simultaneous=True) * J - g
        for e, w, a in self.points():
            point = {eta: e, mu: w, alpha: a, m: sp.Rational(7, 5)}
            worst = max(abs(sp.N(entry.subs(point), 30)) for entry in difference)
            self.assertLess(worst, 1e-20, point)
            # Both Killing vectors are spacelike there, and r lies between Kerr's horizons.
            self.assertGreater(sp.N(g[2, 2].subs(point)), 0)
            self.assertGreater(sp.N((g[2, 2] * g[3, 3] - g[2, 3] ** 2).subs(point)), 0)
            self.assertLess(sp.N((image[1] ** 2 - 2 * m * image[1] + m ** 2 * q ** 2).subs(point)), 0)

    def test_with_the_polarizations_aligned_it_is_schwarzschild_inside_its_horizon(self):
        sp = self.sp
        waves, (eta, mu, x, y), g, _ = self.chart("chandrasekhar_xanthopoulos", "prolate")
        hole, (t, r, theta, phi), s, _ = self.chart("schwarzschild", "spherical")
        m, alpha = waves.parameters["m"], waves.parameters["alpha"]
        # q = 0 and p = 1: r = m(1 - eta) runs from m down to the singularity, r_s = 2m, ct = x, phi = -y/m.
        image = [x / hole.c, m * (1 - eta), sp.acos(mu), -y / m]
        J = sp.Matrix(4, 4, lambda i, j: sp.diff(image[i], (eta, mu, x, y)[j]))
        # The published components are those of the chart x^0 = ct.
        J[0, :] = J[0, :] * hole.c
        at = {r: image[1], theta: image[2], hole.parameters["r_s"]: 2 * m}
        difference = (J.T * s.subs(at, simultaneous=True) * J - g).subs(alpha, 0)
        for e, w, _ in self.points():
            point = {eta: e, mu: w, m: sp.Rational(7, 5)}
            self.assertLess(max(abs(sp.N(entry.subs(point), 30)) for entry in difference), 1e-20, point)

    def test_the_curvature_is_finite_on_the_horizon_and_the_proper_time_to_it_is_m_times_pi_over_two_less_p(self):
        sp = self.sp
        waves, (eta, mu, x, y), g, entry = self.chart("chandrasekhar_xanthopoulos", "prolate")
        m, alpha = waves.parameters["m"], waves.parameters["alpha"]
        at = {m: 1, alpha: sp.atan(sp.Rational(4, 3))}
        K = waves(entry["kretschmann"].partition("=")[2]).subs(at)
        # On eta = 1 and mu = 0 Kerr's radius is m(1 - p) and the scalar is 48/(1 - p)^6.
        self.assertAlmostEqual(float(K.subs({eta: 1, mu: 0})), 48 / (1 - P) ** 6, delta=1e-6)
        for w in (-0.9, -0.3, 0.5, 0.99):
            self.assertTrue(math.isfinite(float(K.subs({eta: 1, mu: w}))))
        lapse = sp.sqrt(-g[0, 0].subs(at).subs(mu, 0))
        self.assertAlmostEqual(float(sp.Integral(lapse, (eta, 0, 1)).evalf(20)), math.pi / 2 - P, delta=1e-12)


class Tube(unittest.TestCase):
    """The ring's world tube as it is written: each row the matrix ((a, sx), (sy, b)) that carries
    the unit circle of the chart to the ring at its psi."""

    def setUp(self):
        views = json.loads((DATA / "embedding" / "chandrasekhar_xanthopoulos.json").read_text(encoding="utf-8"))["views"]
        self.grid = next(v for v in views if v["id"] == "tube")["surfaces"][0]["pieces"][0]["grid"]
        self.rows = list(zip(self.grid["u"], self.grid["a"], self.grid["sx"], self.grid["sy"], self.grid["b"]))

    @staticmethod
    def front(psi):
        """The metric of the wave front on lambda = 0, in x and y, from the line element."""
        eta = math.sin(psi)
        rho = 1 - P * eta
        return ((1 + P * eta) / rho, -2 * Q * eta / rho,
                (4 * Q * Q * eta * eta + (1 - eta * eta) * rho * rho) / (rho * (1 + P * eta)))

    def test_every_row_is_an_orthonormal_frame_of_the_wave_front(self):
        for psi, a, sx, sy, b in self.rows:
            hxx, hxy, hyy = self.front(psi)
            self.assertAlmostEqual(a * a + sy * sy, hxx, delta=2e-6, msg=psi)
            self.assertAlmostEqual(a * sx + sy * b, hxy, delta=2e-6, msg=psi)
            self.assertAlmostEqual(sx * sx + b * b, hyy, delta=2e-6, msg=psi)
            # The area the ring encloses, in units of pi l^2.
            self.assertAlmostEqual(a * b - sx * sy, math.cos(psi), delta=2e-6, msg=psi)

    def test_the_tube_starts_as_the_circle(self):
        psi, a, sx, sy, b = self.rows[0]
        self.assertEqual((psi, a, sx, sy, b), (0, 1, 0, 0, 1))

    def test_the_frame_is_carried_by_parallel_transport(self):
        """The frame's covectors obey dE/dpsi = E h^-1 h'/2 along the world line of the ring's centre,
        with E = 1 at the collision: run here by the classical Runge-Kutta rule, a thousandth of
        psi a step, and held against every row."""
        def rate(psi, E):
            (hxx, hxy, hyy), d = self.front(psi), 1e-6
            up, down = self.front(psi + d), self.front(psi - d)
            dxx, dxy, dyy = ((u - w) / (2 * d) for u, w in zip(up, down))
            det = hxx * hyy - hxy * hxy
            # G = h^-1 h'/2, and the rate is E G.
            g = ((hyy * dxx - hxy * dxy) / (2 * det), (hyy * dxy - hxy * dyy) / (2 * det),
                 (hxx * dxy - hxy * dxx) / (2 * det), (hxx * dyy - hxy * dxy) / (2 * det))
            a, sx, sy, b = E
            return (a * g[0] + sx * g[2], a * g[1] + sx * g[3], sy * g[0] + b * g[2], sy * g[1] + b * g[3])

        E, psi, step = (1.0, 0.0, 0.0, 1.0), 0.0, 1e-3
        for target, a, sx, sy, b in self.rows:
            while psi < target - 1e-12:
                h = min(step, target - psi)
                k1 = rate(psi, E)
                k2 = rate(psi + h / 2, tuple(e + h / 2 * k for e, k in zip(E, k1)))
                k3 = rate(psi + h / 2, tuple(e + h / 2 * k for e, k in zip(E, k2)))
                k4 = rate(psi + h, tuple(e + h * k for e, k in zip(E, k3)))
                E = tuple(e + h / 6 * (p + 2 * q + 2 * r + s) for e, p, q, r, s in zip(E, k1, k2, k3, k4))
                psi += h
            for got, want in zip(E, (a, sx, sy, b)):
                self.assertAlmostEqual(got, want, delta=5e-6, msg=target)

    def test_the_long_axis_turns_through_sixty_three_degrees(self):
        def long_axis(a, sx, sy, b):
            # The direction of the largest stretch: the leading eigenvector of E E^T.
            m00, m01, m11 = a * a + sx * sx, a * sy + sx * b, sy * sy + b * b
            return math.degrees(0.5 * math.atan2(2 * m01, m00 - m11)) % 180
        first = long_axis(*self.rows[1][1:])
        last = long_axis(*self.rows[-1][1:])
        # It starts along -alpha/2 from the x axis, 153.4 degrees, and reaches 90 on the horizon, so
        # by the tube's last row, at psi = 1.2, it has turned through most of those 63.4 degrees.
        self.assertAlmostEqual(first, 180 - math.degrees(math.atan2(Q, P)) / 2, delta=2.0)
        self.assertGreater(first - last, 35)
        self.assertLess(first - last, 180 - math.degrees(math.atan2(Q, P)) / 2 - 90)


if __name__ == "__main__":
    unittest.main()
