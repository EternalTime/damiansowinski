#!/usr/bin/env python3
"""Neugebauer and Meinel's rigidly rotating disc of dust, as numbers, in units of the coordinate
radius rho_0 of the disc, G = c = 1.

The disc is zeta = 0, rho <= 1 in Weyl's coordinates, and the metric is

    ds^2 = e^(-2U) (e^(2k) (drho^2 + dzeta^2) + rho^2 dphi^2) - e^(2U) (dt + a dphi)^2

(G. Neugebauer and R. Meinel, Phys. Rev. Lett. 75, 3046, 1995, their (1)). Everything follows
from the Ernst potential f = e^(2U) + i b, which they give in closed form: with X_1^2 = (i - mu)/mu,
X_2 = -conj(X_1), and the curve

    W^2 = mu^2 (X - X_1)(X - conj X_1)(X - X_2)(X - conj X_2)(X - zeta - i rho)(X - zeta + i rho),

of genus two, f is a quotient of two of its theta functions times an exponential (G. Neugebauer,
A. Kleinwaechter and R. Meinel, Helv. Phys. Acta 69, 472, 1996, their (2.44); G. Neugebauer and
R. Meinel, J. Math. Phys. 44, 3407, 2003, their (113)),

    f = theta(a - C)/theta(a + C) exp(-(gamma_0 u + gamma_1 v + mu w)),
    theta(x) = sum over m, n of (-1)^(m + n) exp(B_11 m^2 + 2 B_12 m n + B_22 n^2 + 2 m x_1 + 2 n x_2),

where u, v and w are the integrals of h, h X and h X^2 over W_1 along the imaginary axis from -i
to i, h = arsinh(mu (1 + X^2))/(pi i sqrt(1 + mu^2 (1 + X^2)^2)) and W_1^2 = (X - zeta)^2 + rho^2
with Re W_1 < 0; d omega_1 and d omega_2 are the differentials (c_0 + c_1 X) dX/W whose periods
round the cuts of X_1 and of X_2 are pi i and 0, and 0 and pi i; B is their matrix of periods
along the two curves that run from those cuts to the cut of zeta -+ i rho; a_j is the j-th
differential's combination of u and v; C_j is minus its integral from zeta - i rho to infinity on
the sheet W -> mu X^3; and gamma_0 + gamma_1 X + mu X^2 over W has no period round either cut.

Point evaluates all of this at one (rho, zeta) with zeta >= 0, by Gauss-Legendre quadrature along
polylines between the branch points, the square root followed continuously along each and the
inverse square roots at the ends taken out by a change of variable. Three arrangements of the six
branch points occur, and Point.ways holds the curves of each: the cut of zeta -+ i rho to the right
of both others, their picture; over and under the cut of X_2, rho > Im X_2, where nothing is
crossed; and between its ends, rho < Im X_2 and zeta < Re X_2, where that cut has let the pair
through and now runs round it. The potential is analytic across all three, which main() checks.

The appendix of the 2003 review gives the metric in the same theta functions, with theta* the
theta function at x + (i pi/2, i pi/2). As printed its second formula has theta*(0) where
theta*(c) belongs; with that put right,

    e^(2U)                           = theta(c) theta*(c) theta(a) theta*(a)/(theta(0) theta*(0) theta(a + c) theta*(a + c)) E,
    1 + (1 + Omega a) e^(2U)/(Omega rho) = theta(0) theta*(0) theta(a + 2c) theta*(a)/(theta(c) theta*(c) theta(a + c) theta*(a + c)),
    1 - (1 + Omega a) e^(2U)/(Omega rho) = theta(0) theta*(0) theta*(a + 2c) theta(a)/(theta(c) theta*(c) theta(a + c) theta*(a + c)),

with E the exponential above and Omega^2 = mu e^(2 V_0)/2, V_0 = U at the centre of the disc. The
third is the second with theta and theta* exchanged in two places, and the two sum to 2, an
identity of the theta functions that main() checks to ten digits, and with it the correction.
Their product over e^(2U) is regular where e^(2U) vanishes, so the metric is written here in
components that stay finite on the ergosurface, where U and a do not:

    F = e^(2U) = -g_tt,   S = F + Omega A = Omega rho (P_+ - P_-)/2,   A = a e^(2U) = -g_t phi,
    e^(2U') = -Omega^2 rho^2 P_+ P_-/F,   g_phi phi = (2 S - F - e^(2U'))/Omega^2,
    g_rho rho = e^(2k)/F = kappa/(kappa_far F),

with U' the potential in the frame that turns with the disc and kappa the review's, which holds
the double integral 2 k_0 over the square of the imaginary axis; kappa_far is its value far away,
where e^(2k) = 1. In Bardeen and Wagoner's form, ds^2 = e^(2 alpha)(drho^2 + dzeta^2) +
rho^2 e^(-2 nu)(dphi - omega dt)^2 - e^(2 nu) dt^2, these are e^(2 nu) = rho^2/g_phi phi,
omega = A/g_phi phi and e^(2 alpha) = g_rho rho.

On the disc itself the conditions Neugebauer and Meinel solved hold: e^(2U') = e^(2 V_0), the
same at every radius, which gives S^2 = e^(2 V_0)(F + mu rho^2/2), and the turning frame's a' has
no derivative across the disc. The theta functions are never told so; disc_conditions works the
metric of the disc out from the two conditions, and the theta functions' own values agree. main() prints every check; the tests in _tools/test_neugebauer_meinel.py
hold the same numbers.

    .venv.noindex/bin/python _tools/derivations/nm_disc.py
"""
import math
from functools import lru_cache

import numpy as np
from numpy.polynomial import chebyshev as cheb
from numpy.polynomial.legendre import leggauss

MU_0 = 4.62966184347          # the first zero of the denominator of V_0(mu), their (20) of 1993
MU_ERGO = 1.68849             # where the ergoregion appears (Neugebauer, Kleinwaechter and Meinel 1996)
NODES = 400                   # Gauss-Legendre nodes on each leg of a curve


@lru_cache(maxsize=None)
def gauss(n):
    """Gauss-Legendre nodes and weights on (0, 1)."""
    x, w = leggauss(n)
    return 0.5 * (x + 1), 0.5 * w


def follow(values):
    """The square root of each value in turn, continuous along the list."""
    root = np.sqrt(values.astype(complex))
    flips = np.real(np.conj(root[:-1]) * root[1:]) < 0
    root[1:] *= np.where(np.cumsum(flips) % 2 == 1, -1.0, 1.0)
    return root


def theta(x, B, star=False, terms=18):
    """The theta function of the docstring; with star, at x + (i pi/2, i pi/2), which drops the signs."""
    m = np.arange(-terms, terms + 1)
    M, N = np.meshgrid(m, m, indexing="ij")
    sign = 1.0 if star else (-1.0) ** (M + N)
    return np.sum(sign * np.exp(B[0, 0] * M * M + 2 * B[0, 1] * M * N + B[1, 1] * N * N + 2 * M * x[0] + 2 * N * x[1]))


def theta_hessian(B, star, terms=18):
    """The second derivatives of ln theta at x = 0, where theta is even."""
    m = np.arange(-terms, terms + 1)
    M, N = np.meshgrid(m, m, indexing="ij")
    sign = np.ones_like(M, dtype=float) if star else (-1.0) ** (M + N)
    e = sign * np.exp(B[0, 0] * M * M + 2 * B[0, 1] * M * N + B[1, 1] * N * N)
    return 4 * np.array([[(e * M * M).sum(), (e * M * N).sum()], [(e * M * N).sum(), (e * N * N).sum()]]) / e.sum()


class Point:
    """Neugebauer and Meinel's solution at one point (rho, zeta), zeta >= 0, for the disc of parameter mu."""

    def __init__(self, mu, rho, zeta, nodes=NODES):
        if not 0 < mu < MU_0:
            raise ValueError(f"mu = {mu} is outside (0, mu_0)")
        if zeta < 0 or rho <= 0:
            raise ValueError("Point wants rho > 0 and zeta >= 0; below the disc f is the complex conjugate")
        self.mu, self.rho, self.zeta, self.n = float(mu), float(rho), float(zeta), nodes
        X1 = np.sqrt(complex(-1, 1 / mu))
        X1 = -X1 if X1.real > 0 else X1
        self.X1, self.X1b, self.X2, self.X2b = X1, np.conj(X1), -np.conj(X1), -X1
        self.eps, self.delta = -X1.real, abs(X1.imag)
        self.Ez, self.Ezb = complex(zeta, -rho), complex(zeta, rho)
        self.points = (self.X1, self.X1b, self.X2, self.X2b, self.Ez, self.Ezb)
        self._solve()

    # ------------------------------------------------------------- the curve
    def _radicand(self, X):
        out = self.mu ** 2 * np.ones_like(X, dtype=complex)
        for e in self.points:
            out = out * (X - e)
        return out

    def _path(self, way, far):
        """Nodes and weights along a polyline that starts on a branch point and ends on one, or, with
        far, runs out to infinity along the direction given last."""
        s, w = gauss(self.n)
        Xs, dXs = [], []
        legs = len(way) - 1
        for i in range(legs):
            A, B = way[i], way[i + 1]
            if far and i == legs - 1:
                t = s / (1 - s)
                if legs == 1:           # from the branch point A itself
                    X, dX = A + B * t * t, B * 2 * t / (1 - s) ** 2
                else:
                    X, dX = A + B * t, B / (1 - s) ** 2
            elif i == 0 and i == legs - 1:
                X, dX = A + (B - A) * np.sin(0.5 * np.pi * s) ** 2, (B - A) * 0.5 * np.pi * np.sin(np.pi * s)
            elif i == 0:
                X, dX = A + (B - A) * s * s, (B - A) * 2 * s
            elif i == legs - 1:
                X, dX = B + (A - B) * (1 - s) ** 2, (B - A) * 2 * (1 - s)
            else:
                X, dX = A + (B - A) * s, (B - A) * np.ones_like(s)
            Xs.append(X), dXs.append(dX * w)
        return np.concatenate(Xs), np.concatenate(dXs)

    def _moments(self, way, far=False, powers=(0, 1, 2)):
        """The integrals of X^k dX/W along the polyline."""
        X, dX = self._path(way, far)
        radicand = self._radicand(X)
        if far:         # the sheet W -> mu X^3, followed inward from the far end
            root = follow(radicand[::-1])[::-1]
            if np.real(np.conj(root[-1]) * self.mu * X[-1] ** 3) < 0:
                root = -root
        else:
            root = follow(radicand)
        return np.array([np.sum(X ** k * dX / root) for k in powers])

    def ways(self):
        """The curves: round the cut of X_1, round the cut of X_2, from each of those cuts to the cut of
        zeta -+ i rho, and from zeta - i rho to infinity."""
        e, d, z, r = self.eps, self.delta, self.zeta, self.rho
        a1 = [self.X1, self.X1b]
        b2 = [self.X2b, self.Ezb]
        if r > d or z > e + 0.3:
            a2 = [self.X2, self.X2b]
            b1 = [self.X1, complex(e, -d - 0.5), self.Ez] if z > e else [self.X1, self.Ez]
            tail = [self.Ez, -1j] if r > d else [self.Ez, 1.0 + 0j]
        else:
            # The pair lies between the ends of the cut of X_2, or just to the right of it. Once it is
            # through, that cut runs round the left of it, and so does the curve round the cut; the
            # same curve serves while the pair is still just to the right, where a straight one would
            # graze it. The curve from the cut of X_1 comes round below X_2 and in from the right.
            left = min(0.5 * (z - e), 0.0)
            a2 = [self.X2, complex(left, -d - 0.4), complex(left, d + 0.4), self.X2b]
            b1 = [self.X1, complex(e, -d - 0.4), complex(e + 0.9, -d - 0.4), complex(e + 0.9, 0.0), self.Ez]
            tail = [self.Ez, complex(e + 1.3, 0.0), 1.0 + 0j]
        return a1, a2, b1, b2, tail

    def _contour(self):
        """Nodes X and weights dX from -i to i for u, v, w and for kappa. Far from the disc the path is
        the imaginary axis. Over the disc, rho < 1 and zeta < 0.6, the integrands' branch points
        zeta -+ i rho lie on the path or beside it, so the path is moved to the left of the axis, which
        crosses nothing, since h is analytic to the right of the cuts of X_1; on the disc itself,
        zeta = 0, that is the limit from above."""
        s, w = gauss(self.n)
        if self.rho < 1 and self.zeta < 0.6:
            c = min(0.5 * self.eps, 0.3)
            way = [-1j, complex(-c, -0.6), complex(-c, 0.6), 1j]
            Xs, dXs = [], []
            # At -i and i the branch points zeta -+ i rho come close when rho is near 1, and the
            # integrands vary as a square root there, which a square in the variable takes out.
            Xs.append(way[0] + (way[1] - way[0]) * s * s), dXs.append((way[1] - way[0]) * 2 * s * w)
            Xs.append(way[1] + (way[2] - way[1]) * s), dXs.append((way[2] - way[1]) * w)
            Xs.append(way[3] + (way[2] - way[3]) * (1 - s) ** 2), dXs.append((way[3] - way[2]) * 2 * (1 - s) * w)
            return np.concatenate(Xs), np.concatenate(dXs)
        y = -1 + 2 * np.sin(0.5 * np.pi * s) ** 2
        return 1j * y, 1j * np.pi * np.sin(np.pi * s) * w

    def _w1(self, X):
        """W_1, the branch with Re W_1 < 0 on the imaginary axis above the disc."""
        return -np.sqrt((X - self.zeta) ** 2 + self.rho ** 2 + 0j)

    def _h(self, X):
        g = self.mu * (1 + X * X)
        return np.arcsinh(g) / (np.pi * 1j * np.sqrt(1 + g * g))

    def _solve(self):
        a1, a2, b1, b2, tail = self.ways()
        rounds = 2 * np.array([self._moments(a1), self._moments(a2)])
        first = np.linalg.solve(rounds[:, :2], np.pi * 1j * np.eye(2))         # column n: d omega_n
        third = np.linalg.solve(rounds[:, :2], -self.mu * rounds[:, 2])        # gamma_0, gamma_1
        across = 2 * np.array([self._moments(b1, powers=(0, 1)), self._moments(b2, powers=(0, 1))])
        B = across @ first
        for i in (0, 1):
            if B[i, i].real > 0:
                across[i] = -across[i]
        B = across @ first
        self.asymmetry = abs(B[0, 1] - B[1, 0])
        B[0, 1] = B[1, 0] = 0.5 * (B[0, 1] + B[1, 0])
        self.B = B
        self.C = -(self._moments(tail, far=True, powers=(0, 1)) @ first)
        X, dX = self._contour()
        weight = self._h(X) * dX / self._w1(X)
        u, v, w = (np.sum(weight * X ** k) for k in range(3))
        self.a = np.array([first[0, 0] * u + first[1, 0] * v, first[0, 1] * u + first[1, 1] * v])
        self.E = np.exp(-(third[0] * u + third[1] * v + self.mu * w))
        t = lambda x, star=False: theta(x, B, star)
        a, c, o = self.a, self.C, np.zeros(2)
        self._t = {name: (t(x), t(x, True)) for name, x in (("0", o), ("a", a), ("c", c), ("a+c", a + c), ("a+2c", a + 2 * c))}
        self._f = t(a - c) / t(a + c) * self.E

    # ------------------------------------------------------------ the fields
    def ernst(self):
        """f = e^(2U) + i b."""
        return self._f

    def F(self):
        """e^(2U), from the review's own formula, which agrees with the real part of f."""
        t = self._t
        return (t["c"][0] * t["c"][1] * t["a"][0] * t["a"][1] / (t["0"][0] * t["0"][1] * t["a+c"][0] * t["a+c"][1]) * self.E).real

    def P(self):
        """P_+ and P_-, 1 +- (1 + Omega a) e^(2U)/(Omega rho)."""
        t = self._t
        den = t["c"][0] * t["c"][1] * t["a+c"][0] * t["a+c"][1]
        return ((t["0"][0] * t["0"][1] * t["a+2c"][0] * t["a"][1] / den).real,
                (t["0"][0] * t["0"][1] * t["a+2c"][1] * t["a"][0] / den).real)

    def turning(self, Omega):
        """e^(2U') in the frame that turns with the disc, -Omega^2 rho^2 P_+ P_-/F with the zero of F
        cancelled against the theta functions of a."""
        t = self._t
        num = (t["0"][0] * t["0"][1]) ** 3 * t["a+2c"][0] * t["a+2c"][1]
        den = (t["c"][0] * t["c"][1]) ** 3 * t["a+c"][0] * t["a+c"][1] * self.E
        return -(Omega * self.rho) ** 2 * (num / den).real

    def metric(self, Omega):
        """(F, A, g_phi phi) with F = -g_tt and A = -g_t phi, regular on the ergosurface."""
        F = self.F()
        plus, minus = self.P()
        S = Omega * self.rho * (plus - minus) / 2
        return F, (S - F) / Omega, (2 * S - F - self.turning(Omega)) / Omega ** 2

    def kappa_over_F(self):
        """kappa/F of the review's appendix, regular where F vanishes: g_rho rho times kappa far away."""
        X, dX = self._contour()
        lam = self._w1(X) / (X - self.Ez)
        h = self._h(X)
        left = h * (X - self.X1) * (X - self.X2) * dX
        right = h * (X + self.X1) * (X + self.X2) * dX
        L, Lp = lam[:, None], lam[None, :]
        with np.errstate(all="ignore"):
            kernel = (L - Lp) ** 2 / (L * Lp) / (X[:, None] - X[None, :]) ** 2
        slope = 0.5 * (1 / (X - self.Ezb) - 1 / (X - self.Ez))          # lam'/lam, for the diagonal
        kernel[np.diag_indices_from(kernel)] = slope ** 2
        k0 = self.mu ** 2 / 4 * np.sum(kernel * left[:, None] * right[None, :])
        H = theta_hessian(self.B, False) + theta_hessian(self.B, True)
        t = self._t
        ratio = t["a+c"][0] * t["a+c"][1] / (t["c"][0] * t["c"][1] * self.E)
        return (ratio * np.exp(k0 - 0.5 * self.a @ H @ self.a)).real


# ------------------------------------------------------------------ one disc
@lru_cache(maxsize=None)
def centre(mu):
    """f at the centre of the disc, from above: f is even in rho, so two radii give it to fourth order."""
    near, far = Point(mu, 0.01, 0.0).ernst(), Point(mu, 0.02, 0.0).ernst()
    return (4 * near - far) / 3


def e2V0(mu):
    """e^(2 V_0), the real part of f at the centre of the disc."""
    return centre(mu).real


def omega_disc(mu):
    """Omega rho_0, from mu = 2 Omega^2 rho_0^2 e^(-2 V_0)."""
    return math.sqrt(mu * e2V0(mu) / 2)


def weierstrass(x, g2, g3):
    """Weierstrass's function, the p with x = integral from p to infinity of dt/sqrt(4t^3 - g2 t - g3),
    for x below the first half period, by bisection on the integral."""
    from scipy.integrate import quad
    from scipy.optimize import brentq
    e1 = max(np.roots([4, 0, -g2, -g3]).real)

    def integral(p):
        # t = p + s^2/(1 - s)^2 takes out nothing but the tail; the integrand is regular for p > e1
        return quad(lambda s: 2 * s / (1 - s) ** 3 / math.sqrt(max(4 * (p + (s / (1 - s)) ** 2) ** 3
                                                                   - g2 * (p + (s / (1 - s)) ** 2) - g3, 1e-300)), 0, 1, limit=400)[0]
    return brentq(lambda p: integral(p) - x, e1 + 1e-12, 1e8, xtol=1e-14, rtol=1e-14)


def V0_closed(mu):
    """V_0(mu) in closed form (Neugebauer and Meinel 1994, as their 1996 paper gives it, (3.1) and (3.2)):
    V_0 = -arsinh(mu + (1 + mu^2)/(p - 2 mu/3))/2 with p Weierstrass's function at I(mu)."""
    from scipy.integrate import quad
    I = quad(lambda s: 2 * math.asinh(mu - s * s) / math.sqrt(1 + (mu - s * s) ** 2), 0, math.sqrt(mu))[0] / math.pi
    p = weierstrass(I, 4 * mu ** 2 / 3 - 4, 8 * mu * (1 + mu ** 2 / 9) / 3)
    return -0.5 * math.asinh(mu + (1 + mu ** 2) / (p - 2 * mu / 3))


class Fit:
    """A function known at Chebyshev points of [lo, hi], with its derivatives."""

    def __init__(self, f, lo, hi, degree):
        x = 0.5 * (lo + hi) + 0.5 * (hi - lo) * np.cos(np.pi * (np.arange(degree + 1) + 0.5) / (degree + 1))
        self.x, self.values = x, np.array([f(v) for v in x])
        self.series = cheb.Chebyshev.fit(x, self.values, degree, domain=[lo, hi])

    def __call__(self, x, order=0):
        return (self.series.deriv(order) if order else self.series)(x)


class Table:
    """A table of several functions of one variable at Chebyshev points, by rows."""

    def __init__(self, rows, lo, hi, degree, names):
        x = 0.5 * (lo + hi) + 0.5 * (hi - lo) * np.cos(np.pi * (np.arange(degree + 1) + 0.5) / (degree + 1))
        values = np.array([rows(v) for v in x])
        self.x = x
        self.series = {n: cheb.Chebyshev.fit(x, values[:, i], degree, domain=[lo, hi]) for i, n in enumerate(names)}

    def __call__(self, name, x, order=0):
        s = self.series[name]
        return (s.deriv(order) if order else s)(x)


DISC_DEGREE, OUT_DEGREE = 56, 72


@lru_cache(maxsize=None)
def kappa_far(mu):
    """kappa far away, where e^(2k) = 1."""
    p = Point(mu, 3000.0, 0.0)
    return p.kappa_over_F() * p.F()


@lru_cache(maxsize=None)
def disc_table(mu):
    """F, A/rho^2, g_phi phi/rho^2 and g_rho rho on the disc, zeta = 0 from above, as functions of
    eta = sqrt(1 - rho^2), from the theta functions. They are smooth in eta and not in eta^2: odd
    powers of eta from the third on enter at the rim, as odd powers of xi do outside it."""
    Om, far = omega_disc(mu), kappa_far(mu)

    def rows(eta):
        rho = math.sqrt(1 - eta * eta)
        p = Point(mu, rho, 0.0)
        F, A, gpp = p.metric(Om)
        return F, A / rho ** 2, gpp / rho ** 2, p.kappa_over_F() / far
    return Table(rows, 0.0, 1.0, DISC_DEGREE, ("F", "A", "gpp", "grr"))


@lru_cache(maxsize=None)
def plane_table(mu):
    """F, A, g_phi phi/rho^2 and g_rho rho on the plane zeta = 0 outside the disc, as functions of
    s = xi/(1 + xi), rho = sqrt(1 + xi^2)."""
    Om, far = omega_disc(mu), kappa_far(mu)

    def rows(s):
        xi = s / (1 - s)
        rho = math.sqrt(1 + xi * xi)
        p = Point(mu, rho, 0.0)
        F, A, gpp = p.metric(Om)
        return F, A, gpp / rho ** 2, p.kappa_over_F() / far
    return Table(rows, 0.0, 1.0, OUT_DEGREE, ("F", "A", "gpp", "grr"))


RIM = 0.01          # the reach, in s or in eta, of the correction that levels each fit at the rim


def _levelled(table, name, s, order=0):
    """A function of the plane from its Chebyshev fit in s outside the disc, or in eta on it, or its
    slope in that variable. Each has no slope at the rim, s = 0 or eta = 0, where d rho/ds and
    d rho/d eta vanish; the fit's own slope there, a few parts in a million, is taken off over the
    first hundredth of the variable, which moves no value by a millionth."""
    c = table(name, 0.0, 1)
    if order == 0:
        return table(name, s) - c * s * np.exp(-s / RIM)
    return table(name, s, 1) - c * (1 - s / RIM) * np.exp(-s / RIM)


def _over(table, name, s):
    """The slope of a levelled fit over its variable, finite at the rim, where it is the second
    derivative."""
    close = s < 1e-6
    return np.where(close, table(name, 0.0, 2) + 2 * table(name, 0.0, 1) / RIM, _levelled(table, name, s, 1) / np.where(close, 1.0, s))


def plane(mu, rho, order=0):
    """(F, A, G, g_rho rho) on the plane zeta = 0 at rho >= 0, on the disc and outside it, with
    F = -g_tt, A = -g_t phi and G = g_phi phi/rho^2, in units of rho_0; with order 1, their
    derivatives along rho, by the chain rule through eta = sqrt(1 - rho^2) on the disc, where
    d rho = -eta d eta/rho, and through s = xi/(1 + xi), xi = sqrt(rho^2 - 1), outside it, where
    d rho = xi (1 + xi)^2 ds/rho. Each is finite everywhere, the ergosurface and the centre included."""
    rho = np.asarray(rho, dtype=float)
    inside = rho <= 1
    d, t = disc_table(mu), plane_table(mu)
    r = np.where(inside, rho, 0.5)
    eta = np.sqrt(np.maximum(1 - r * r, 0.0))
    R = np.where(inside, 1.5, rho)
    xi = np.sqrt(np.maximum(R * R - 1, 0.0))
    s = xi / (1 + xi)
    names = ("F", "A", "gpp", "grr")
    if order == 0:
        disc = tuple(_levelled(d, name, eta) * (r ** 2 if name == "A" else 1.0) for name in names)
        out = tuple(_levelled(t, name, s) for name in names)
    else:
        along = {name: -r * _over(d, name, eta) for name in names}
        along["A"] = along["A"] * r ** 2 + 2 * r * _levelled(d, "A", eta)
        disc = tuple(along[name] for name in names)
        out = tuple(_over(t, name, s) * R / (1 + xi) ** 3 for name in names)
    return tuple(np.where(inside, x, y) for x, y in zip(disc, out))


def disc_conditions(mu, rho):
    """(A, g_phi phi, g_rho rho) on the disc from the conditions Neugebauer and Meinel solved, with
    only F and b on the disc taken from the theta functions, out to rho = 0.97 rho_0: e^(2U') =
    e^(2 V_0) gives S = sqrt(e^(2 V_0)(F + mu rho^2/2)), A = (S - F)/Omega and g_phi phi =
    (2S - F - e^(2 V_0))/Omega^2. The turning frame's a' = (1 - S e^(-2 V_0))/Omega has no derivative
    across the disc, which with a_zeta = -rho b_rho/F^2 gives the derivatives across it,

        (F_z + Omega A_z)(2 S^2 - e^(2 V_0) F) = e^(2 V_0) S F_z,      F A_z - A F_z = -rho b_rho,

    and k', which vanishes at the centre, by k'_rho = -rho (U'_zeta^2 + e^(4 V_0)(a'_rho/rho)^2/4) with
    U'_zeta = (F_z + Omega A_z)/(2S); g_rho rho = e^(2k' - 2 V_0). The tests hold the theta functions'
    own values to these."""
    V, Om = e2V0(mu), omega_disc(mu)
    table = _disc_potential(mu)

    def slope(r):
        eta = math.sqrt(1 - r * r)
        F = table("F", eta)
        Fr, br = (-r * table(name, eta, 1) / eta for name in ("F", "b"))
        S = math.sqrt(V * (F + mu * r ** 2 / 2))
        Sr = V * (Fr + mu * r) / (2 * S)
        A = (S - F) / Om
        c = 2 * S ** 2 - V * F
        Fz, Az = np.linalg.solve(np.array([[c - V * S, Om * c], [-A, F]]), np.array([0.0, -r * br]))
        return -r * ((Fz + Om * Az) / (2 * S)) ** 2 - V ** 2 * (Sr / (V * Om)) ** 2 / (4 * r)
    kprime = Fit(slope, 0.0, 0.97, DISC_DEGREE).series.integ(lbnd=0.0)     # b's slope along rho has no bound at the rim
    rho = np.asarray(rho, dtype=float)
    F = table("F", np.sqrt(1 - rho * rho))
    S = np.sqrt(V * (F + mu * rho ** 2 / 2))
    return (S - F) / Om, (2 * S - F - V) / Om ** 2, np.exp(2 * kprime(rho)) / V


@lru_cache(maxsize=None)
def _disc_potential(mu):
    return Table(lambda eta: (lambda f: (f.real, f.imag))(Point(mu, math.sqrt(1 - eta * eta), 0.0).ernst()), 0.0, 1.0,
                 DISC_DEGREE, ("F", "b"))


def profile(mu, name, rho, order=0):
    """One function of the plane zeta = 0 at rho >= 0, or with order 1 or 2 its derivatives along rho:
    F = -g_tt, A = -g_t phi, G = g_phi phi/rho^2 and grr of plane(); omega = A/(rho^2 G), the dragging
    of Bardeen and Wagoner's form; S = F + Omega A; and E = F + 2 Omega A - Omega^2 rho^2 G, the
    e^(2U') of the frame that turns with the disc. The second derivative is the difference of the
    first over 2e-4 of rho_0."""
    rho = np.asarray(rho, dtype=float)
    if order == 2:
        h = 1e-4
        return (profile(mu, name, rho + h, 1) - profile(mu, name, np.abs(rho - h), 1) * np.where(rho < h, -1.0, 1.0)) / (2 * h)
    Om = omega_disc(mu)
    F, A, G, grr = plane(mu, rho)
    if order == 1:
        dF, dA, dG, dgrr = plane(mu, rho, 1)
    if name in ("F", "A", "G", "grr"):
        index = ("F", "A", "G", "grr").index(name)
        return plane(mu, rho, order)[index]
    if name == "S":
        return F + Om * A if order == 0 else dF + Om * dA
    if name == "E":
        if order == 0:
            return F + 2 * Om * A - Om ** 2 * rho ** 2 * G
        return dF + 2 * Om * dA - Om ** 2 * (2 * rho * G + rho ** 2 * dG)
    if name == "omega":
        inside = rho <= 1
        r = np.where(inside, rho, 0.5)          # on the disc A/rho^2 over G, finite at the centre
        R = np.where(inside, 1.5, rho)
        d = disc_table(mu)
        eta = np.sqrt(np.maximum(1 - r * r, 0.0))
        a, g = _levelled(d, "A", eta), _levelled(d, "gpp", eta)
        if order == 0:
            return np.where(inside, a / g, A / (R ** 2 * G))
        return np.where(inside, -r * (_over(d, "A", eta) * g - a * _over(d, "gpp", eta)) / g ** 2,
                        (dA - A * (2 / R + dG / G)) / (R ** 2 * G))
    raise KeyError(name)


def light_cylinder(mu):
    """The radius in the plane at which a point at rest in the turning frame moves at the speed of
    light, the zero of e^(2U')."""
    from scipy.optimize import brentq
    return brentq(lambda r: float(profile(mu, "E", r)), 1.0, 200.0, xtol=1e-13)


def ergosurface(mu):
    """The radii at which the ergosurface crosses the plane, the zeros of F, the inner one on the disc;
    none below mu = 1.68849."""
    from scipy.optimize import brentq
    r = np.linspace(1e-3, 60, 6000)
    F = profile(mu, "F", r)
    return [brentq(lambda x: float(profile(mu, "F", x)), a, b, xtol=1e-13) for a, b, fa, fb in zip(r[:-1], r[1:], F[:-1], F[1:]) if fa * fb < 0]


AXIS_DEGREE = 64


@lru_cache(maxsize=None)
def axis_table(mu):
    """The real and imaginary parts of f on the axis rho = 0 above the disc, as functions of
    s = zeta/(1 + zeta): f is even in rho, so two radii give it to fourth order."""
    def rows(s):
        zeta = s / (1 - s)
        if zeta == 0:
            f = centre(mu)
        else:
            near = 0.01 * (1 + zeta)             # the radii grow with the height, as the field's own scale does
            f = (4 * Point(mu, near, zeta).ernst() - Point(mu, 2 * near, zeta).ernst()) / 3
        return f.real, f.imag
    return Table(rows, 0.0, 1.0, AXIS_DEGREE, ("F", "b"))


def axis(mu, zeta, order=0):
    """e^(2U) on the axis at the height zeta, of either sign, or with order 1 or 2 its derivatives along
    the axis, the second a difference of the first; a = 0 and k = 0 there."""
    zeta = np.asarray(zeta, dtype=float)
    z = np.abs(zeta)
    s = z / (1 + z)
    if order == 0:
        return axis_table(mu)("F", s)
    if order == 2:
        h = 1e-4
        return (axis(mu, z + h, 1) - np.where(z < h, -1.0, 1.0) * axis(mu, np.abs(z - h), 1)) / (2 * h)
    return np.sign(zeta) * axis_table(mu)("F", s, 1) / (1 + z) ** 2


def axis_dragging(mu, zeta, order=0):
    """omega on the axis in Bardeen and Wagoner's form, half the derivative of b along the axis:
    near the axis a = rho^2 b_zeta e^(-4U)/2 and omega = a e^(4U)/rho^2."""
    z = np.abs(np.asarray(zeta, dtype=float))
    if order == 1:
        h = 1e-4
        return np.sign(zeta) * (axis_dragging(mu, z + h) - axis_dragging(mu, np.abs(z - h))) / (2 * h)
    return 0.5 * axis_table(mu)("b", z / (1 + z), 1) / (1 + z) ** 2


def far_field(mu):
    """The mass M and the angular momentum J in units of rho_0, from F = 1 - 2M/r + 2M^2/r^2 and
    A = 2J/rho on the plane far away, good to a few parts in ten thousand up to mu = 3. Nearer
    mu_0 the mass still is, and J, which enters the binding energy against a small e^(V_0), is not."""
    Om = omega_disc(mu)
    out = []
    for rho in (50.0 / Om, 100.0 / Om):           # far out in units of the mass, which is of the order of 1/Omega
        F, A, _ = Point(mu, rho, 0.0).metric(Om)
        M = rho * (1 - math.sqrt(max(2 * F - 1, 0.0))) / 2     # 1 - 2x + 2x^2 = F for x = M/rho
        out.append((M, A * rho / 2))
    (M1, J1), (M2, J2) = out
    return 2 * M2 - M1, 2 * J2 - J1                          # the next order falls as 1/rho


def main():
    print("Neugebauer and Meinel's disc: checks")
    for mu in (0.5, 1.0, 3.0, 4.5):
        V, Om = e2V0(mu), omega_disc(mu)
        f0 = centre(mu)
        print(f"  mu = {mu}: e^(2V_0) = {V:.10f}, closed form {math.exp(2 * V0_closed(mu)):.10f}; Omega rho_0 = {Om:.8f}; "
              f"|f_0|^2 + 4 Omega^2 - 1 = {abs(f0) ** 2 + 4 * Om ** 2 - 1:.1e}")
        worst = {"asymmetry": 0.0, "P sum": 0.0, "F": 0.0, "turning": 0.0}
        for rho, zeta in ((2.0, 3.0), (0.5, 1.0), (1.5, 0.8), (2.0, 0.0), (1.05, 0.0), (0.5, 0.0), (0.3, 0.5 * Point(mu, 1, 1).eps), (0.9, 0.0)):
            p = Point(mu, rho, zeta)
            plus, minus = p.P()
            worst["asymmetry"] = max(worst["asymmetry"], p.asymmetry)
            worst["P sum"] = max(worst["P sum"], abs(plus + minus - 2))
            worst["F"] = max(worst["F"], abs(p.F() - p.ernst().real))
            if zeta == 0 and rho < 1:
                worst["turning"] = max(worst["turning"], abs(p.turning(Om) / V - 1))
        print("     B_12 - B_21 {asymmetry:.1e}, P_+ + P_- - 2 {P sum:.1e}, e^(2U) against Re f {F:.1e}, "
              "e^(2U')/e^(2V_0) - 1 on the disc {turning:.1e}".format(**worst))
        r = np.array([0.3, 0.7, 0.95])
        theta_values, held = plane(mu, r), disc_conditions(mu, r)
        print("     on the disc, theta functions against the disc's conditions: A {:.1e}, g_phiphi {:.1e}, g_rhorho {:.1e}".format(
            *(float(np.max(np.abs(x / y - 1))) for x, y in zip((theta_values[1], theta_values[2] * r ** 2, theta_values[3]), held))))
        inner, outer = plane(mu, 1 - 1e-9), plane(mu, 1 + 1e-9)
        print(f"     at the rim F = {float(inner[0]):.8f} and {float(outer[0]):.8f} (1 - mu/2 = {1 - mu / 2}), "
              f"g_rhorho = {float(inner[3]):.6f} and {float(outer[3]):.6f}")
        M, J = far_field(mu)
        M0 = (M - 2 * Om * J) / math.sqrt(V)
        print(f"     M = {M:.6f}, J = {J:.6f}, 2 Omega M = {2 * Om * M:.6f}, M^2/J = {M * M / J:.6f}, binding energy {(M0 - M) / M0:.6f}")


if __name__ == "__main__":
    main()
