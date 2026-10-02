"""Brill's time-symmetric gravitational wave, solved: the conformal factor, the mass and the
apparent horizon of Holz, Miller, Wakano and Wheeler's wave.

The slice is psi^4 (e^(2q) (drho^2 + dz^2) + rho^2 dphi^2) in cylindrical coordinates (rho, z, phi),
with q = a rho^2 exp(-r^2), r^2 = rho^2 + z^2, every length in units of the wave's width lambda. At a
moment of time symmetry the one constraint is R = 0, which is Brill's equation (Ann. Phys. 7, 466,
1959),

    Laplacian(psi) + V psi = 0,    V = (d_rho^2 q + d_z^2 q)/4 = (a/4) e^(-r^2) (2 - 12 rho^2 + 4 rho^2 r^2),

with the flat Laplacian and psi -> 1 far away: a scattering problem at zero energy off a potential
that integrates to zero. In the Legendre polynomials of mu = cos(theta), V = V_0(r) + V_2(r) P_2(mu),
so psi = 1 + sum over even l of u_l(r) P_l(mu) with

    u_l'' + 2 u_l'/r - l(l + 1) u_l/r^2 + V_0 u_l + V_2 sum_l' G_ll' u_l' = -V_0 delta_l0 - V_2 delta_l2,

G_ll' = (2l + 1)/2 times the integral of P_l P_2 P_l'. Each u_l is collocated on Chebyshev's points
of [0, EDGE], regular at the centre, and beyond the edge, where V is below 1e-27, it is the
multipole u_l(EDGE) (EDGE/r)^(l + 1), which is the condition u_l' + (l + 1) u_l/r = 0 at the edge.

The mass, as a length GM/c^2, is twice the coefficient of 1/r in psi, M = 2 EDGE u_0(EDGE), and
Brill's theorem is that it equals (1/2 pi) times the integral of (grad psi/psi)^2 over flat space,
which is positive. Alcubierre and others (Class. Quantum Grav. 17, 2159, 2000, Table V) give M =
0.0338, 0.1262, 0.696, 2.912 and 4.67 at a = 1, 2, 5, 10 and 12.

A surface of revolution about the axis is a curve (rho(s), z(s)) of the half plane, s its flat arc
length and theta the angle of its tangent. Its area is the integral of 2 pi psi^4 e^q rho ds, so it
is minimal exactly when the curve is a geodesic of (psi^4 e^q rho)^2 (drho^2 + dz^2),

    dtheta/ds = n . grad ln(psi^4 e^q rho),    n = (-sin theta, cos theta),

and on the axis, where q vanishes and the surface closes smoothly, dtheta/ds = 2 d_z psi/psi. At a
moment of time symmetry a minimal surface is marginally trapped, and the outermost is the apparent
horizon. A curve shot from the axis at right angles closes into a surface when it crosses the
plane z = 0 at right angles. Alcubierre and others find the first such surface at an amplitude
between 11.81 and 11.82; critical() finds 11.82.

It needs numpy and scipy and nothing else, so the tests run it as it stands.
"""
import functools
import math

import numpy as np
from numpy.polynomial import legendre
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, minimize_scalar

EDGE = 8.0          # the radius beyond which the potential is dropped, e^(-64) of its size
NODES = 96          # Chebyshev's points along the radius
MODES = 16          # even Legendre polynomials, l = 0, 2, ..., 30
RTOL, ATOL = 1e-10, 1e-12
FAR = 40.0          # a curve that strays this far has missed

# What the drawings say they were drawn with.
INPUT = ("The conformal factor $\\psi$ of Holz, Miller, Wakano, and Wheeler's wave, the solution of Brill's equation "
         "$\\nabla^2\\psi + \\tfrac{1}{4}\\psi\\left(\\partial_\\rho^2 q + \\partial_z^2 q\\right) = 0$ that goes to $1$ "
         "far away, found numerically.")


def q(a, rho, z):
    """Holz, Miller, Wakano and Wheeler's q and its two first derivatives at (rho, z)."""
    e = a * math.exp(-rho * rho - z * z)
    return e * rho * rho, 2 * e * rho * (1 - rho * rho), -2 * e * rho * rho * z


def potential(a, rho, z):
    """V = (d_rho^2 q + d_z^2 q)/4, the potential of Brill's equation."""
    r2 = rho * rho + z * z
    return a / 4 * np.exp(-r2) * (2 - 12 * rho * rho + 4 * rho * rho * r2)


def _chebyshev(n):
    """Chebyshev's points on [0, EDGE], from the edge down to the centre, and the matrix that
    differentiates a function given on them (Trefethen, Spectral Methods in MATLAB, chapter 6)."""
    k = np.arange(n + 1)
    x = np.cos(math.pi * k / n)
    c = np.where((k == 0) | (k == n), 2.0, 1.0) * (-1.0) ** k
    X = np.tile(x, (n + 1, 1)).T
    D = np.outer(c, 1 / c) / (X - X.T + np.eye(n + 1))
    D -= np.diag(D.sum(axis=1))
    return EDGE * (x + 1) / 2, D * 2 / EDGE


@functools.lru_cache(maxsize=None)
def _grid(nodes, modes):
    r, D = _chebyshev(nodes)
    # (2l + 1)/2 times the integral of P_l P_2 P_l' over mu, by Gauss's rule, exact for these degrees.
    mu, w = legendre.leggauss(2 * modes + 4)
    P = np.array([legendre.legval(mu, [0] * (2 * k) + [1]) for k in range(modes)])
    P2 = (3 * mu * mu - 1) / 2
    G = np.array([[(4 * i + 1) / 2 * np.sum(w * P[i] * P2 * P[j]) for j in range(modes)] for i in range(modes)])
    weights = np.where((np.arange(nodes + 1) == 0) | (np.arange(nodes + 1) == nodes), 0.5, 1.0) * (-1.0) ** np.arange(nodes + 1)
    return r, D, G, weights


class Wave:
    """The wave of amplitude a: psi = 1 + sum_l u_l(r) P_l(mu) on [0, EDGE] and its multipoles beyond."""

    def __init__(self, a, nodes=NODES, modes=MODES):
        self.a, self.nodes, self.modes = float(a), nodes, modes
        r, D, G, self.weights = _grid(nodes, modes)
        self.r = r
        n = nodes + 1
        V0 = a / 4 * np.exp(-r * r) * (2 + (4 * r ** 4 - 12 * r * r) * 2 / 3)
        V2 = -a / 4 * np.exp(-r * r) * (4 * r ** 4 - 12 * r * r) * 2 / 3
        inside = slice(1, nodes)
        ri = r[inside]
        radial = (D @ D)[inside] + (2 / ri)[:, None] * D[inside]
        A = np.zeros((modes * n, modes * n))
        b = np.zeros(modes * n)
        for i in range(modes):
            l = 2 * i
            rows = slice(i * n + 1, i * n + nodes)
            A[rows, i * n:(i + 1) * n] = radial
            A[rows, i * n + 1:i * n + nodes] += np.diag(-l * (l + 1) / ri ** 2 + V0[inside])
            for j in range(modes):
                if G[i, j] != 0:
                    A[rows, j * n + 1:j * n + nodes] += np.diag(V2[inside] * G[i, j])
            if i == 0:
                b[rows] = -V0[inside]
            elif i == 1:
                b[rows] = -V2[inside]
            # The edge, r[0]: the multipole falls off as r^-(l + 1).
            A[i * n, i * n:(i + 1) * n] = D[0]
            A[i * n, i * n] += (l + 1) / EDGE
            # The centre, r[-1]: u_0' = 0, and u_l = 0 for l > 0.
            if i == 0:
                A[i * n + nodes, i * n:(i + 1) * n] = D[nodes]
            else:
                A[i * n + nodes, i * n + nodes] = 1.0
        self.u = np.linalg.solve(A, b).reshape(modes, n)
        self.du = self.u @ D.T
        self.mass = 2 * EDGE * self.u[0, 0]
        self.l = 2 * np.arange(modes)
        # P_l(0) for the even l.
        self.at_equator = np.array([legendre.legval(0.0, [0] * (2 * k) + [1]) for k in range(modes)])

    def _radial(self, r):
        """u_l(r) and u_l'(r) for every l at one radius."""
        if r >= EDGE:
            f = self.u[:, 0] * (EDGE / r) ** (self.l + 1)
            return f, -(self.l + 1) * f / r
        d = r - self.r
        k = int(np.argmin(np.abs(d)))
        if abs(d[k]) < 1e-14:
            return self.u[:, k], self.du[:, k]
        w = self.weights / d
        return (self.u @ w) / w.sum(), (self.du @ w) / w.sum()

    def psi(self, rho, z):
        """psi and its two first derivatives at (rho, z)."""
        r = math.hypot(rho, z)
        if r < 1e-12:
            return 1 + self.u[0, -1], 0.0, 0.0
        mu = z / r
        f, df = self._radial(r)
        c = np.zeros(2 * self.modes)
        c[::2] = f
        value = 1 + legendre.legval(mu, c)
        along_mu = legendre.legval(mu, legendre.legder(c))
        c[::2] = df
        along_r = legendre.legval(mu, c)
        # d(mu)/d(rho) = -rho z/r^3 and d(mu)/dz = rho^2/r^3.
        return value, along_r * rho / r - along_mu * rho * z / r ** 3, along_r * z / r + along_mu * rho * rho / r ** 3

    def equator(self, rho):
        """psi and d(psi)/d(rho) on the plane z = 0, for an array of rho: there mu = 0, so each is
        the sum of u_l or u_l' times P_l(0)."""
        shape = np.shape(rho)
        r = np.atleast_1d(np.asarray(rho, dtype=float)).ravel()
        f, df = np.empty((len(r), self.modes)), np.empty((len(r), self.modes))
        out = r >= EDGE
        if out.any():
            f[out] = self.u[:, 0] * (EDGE / r[out, None]) ** (self.l + 1)
            df[out] = -(self.l + 1) * f[out] / r[out, None]
        if (~out).any():
            d = r[~out, None] - self.r[None, :]
            on = np.abs(d) < 1e-14
            w = self.weights / np.where(on, 1.0, d)
            w[on.any(axis=1)] = on[on.any(axis=1)]
            total = w.sum(axis=1, keepdims=True)
            f[~out], df[~out] = w @ self.u.T / total, w @ self.du.T / total
        return (1 + f @ self.at_equator).reshape(shape), (df @ self.at_equator).reshape(shape)

    def brill_mass(self, n_r=400, n_mu=48):
        """Brill's integral for the mass, (1/2 pi) times the integral of (grad psi/psi)^2 over flat
        space: by Gauss's rule inside the edge, and beyond it term by term in the multipoles, where
        psi - 1 is small and the integral of (grad psi)^2/psi^2 is taken with psi's own values."""
        mu, w = legendre.leggauss(n_mu)
        total = 0.0
        for lo, hi, n in ((0.0, EDGE, n_r), (EDGE, 40 * EDGE, n_r)):
            # Beyond the edge the integrand falls as 1/r^2, so the points are spread evenly in 1/r.
            x, wx = legendre.leggauss(n)
            if lo == 0.0:
                r, wr = (hi - lo) / 2 * x + (hi + lo) / 2, (hi - lo) / 2 * wx
            else:
                s = (1 / lo - 1 / hi) / 2 * x + (1 / lo + 1 / hi) / 2
                r, wr = 1 / s, (1 / lo - 1 / hi) / 2 * wx / s ** 2
            for ri, wi in zip(r, wr):
                for m, wm in zip(mu, w):
                    if m < 0:
                        continue
                    rho, z = ri * math.sqrt(1 - m * m), ri * m
                    p, pr, pz = self.psi(rho, z)
                    total += 2 * wi * wm * ri * ri * (pr * pr + pz * pz) / (p * p)
        # The rest, beyond 40 EDGE, is the monopole's: (M^2/2)/(r + M/2).
        far = 40 * EDGE
        return total + self.mass ** 2 / 2 / (far + self.mass / 2)

    def residual(self, rho, z, h=1e-3):
        """Laplacian(psi) + V psi at (rho, z) by central differences, which vanishes for a solution."""
        def f(dr=0.0, dz=0.0):
            return self.psi(rho + dr, z + dz)[0]
        lap = ((f(h) + f(-h) - 2 * f()) / h ** 2 + (f(h) - f(-h)) / (2 * h * rho) + (f(0, h) + f(0, -h) - 2 * f()) / h ** 2)
        return lap + potential(self.a, rho, z) * f()

    # -- minimal surfaces ---------------------------------------------------------------

    def turning(self, rho, z, theta):
        """dtheta/ds of the minimal surface through (rho, z) with tangent angle theta."""
        p, pr, pz = self.psi(rho, z)
        if rho < 1e-9:
            return 2 * pz / p * (1 if math.cos(theta) > 0 else -1)
        _, qr, qz = q(self.a, rho, z)
        s, c = math.sin(theta), math.cos(theta)
        return 4 * (-s * pr + c * pz) / p + (-s * qr + c * qz) - s / rho

    def shoot(self, z0, dense=False):
        """The minimal surface that leaves the axis at z0 at right angles, until it crosses the
        plane z = 0 or strays."""
        def rhs(s, y):
            return [math.cos(y[2]), math.sin(y[2]), self.turning(y[0], y[1], y[2])]

        def plane(s, y):
            return y[1] if y[0] > 1e-6 else 1.0

        def lost(s, y):
            return FAR - math.hypot(y[0], y[1])
        plane.terminal = lost.terminal = True
        return solve_ivp(rhs, [0, 4 * FAR], [0.0, z0, 0.0], events=[plane, lost], rtol=RTOL, atol=ATOL,
                         dense_output=dense, first_step=1e-6)

    def miss(self, z0):
        """The cosine of the tangent's angle where the curve from z0 crosses the plane z = 0, which
        is zero for a surface that the reflection z -> -z closes."""
        sol = self.shoot(z0)
        if not len(sol.t_events[0]):
            return float("nan")
        return math.cos(float(sol.y_events[0][0][2]))

    def deepest(self, lo=0.4, hi=6.0, n=57):
        """Where on the axis miss is least, and its value there. Far up the axis the curve is
        nearly a sphere of flat space and the miss is positive; it is negative only between an inner
        and an outer minimal surface, so they exist exactly when the least value is negative."""
        zs = np.linspace(lo, hi, n)
        v = np.array([self.miss(z) for z in zs])
        v[~np.isfinite(v)] = np.inf
        k = int(np.argmin(v))
        found = minimize_scalar(self.miss, bounds=(zs[max(k - 1, 0)], zs[min(k + 1, n - 1)]), method="bounded",
                                options={"xatol": 1e-11})
        return float(found.x), float(found.fun)

    def horizons(self):
        """The points z_0 of the axis that the inner and the outer minimal surface leave from, the
        outer being the apparent horizon, or None where the wave has none."""
        z, least = self.deepest()
        if least >= 0:
            return None
        inner = brentq(self.miss, 0.4, z, xtol=1e-12, rtol=1e-13) if self.miss(0.4) > 0 else None
        return inner, brentq(self.miss, z, 6.0, xtol=1e-12, rtol=1e-13)

    def curve(self, z0):
        """The curve from z0 as a dense solution, with the arc length at which it meets z = 0."""
        sol = self.shoot(z0, dense=True)
        return sol.sol, float(sol.t_events[0][0])

    def waist(self, z0):
        """rho of the circle in which the minimal surface from z0 cuts the plane z = 0."""
        path, end = self.curve(z0)
        return float(path(end)[0])

    def area(self, z0, n=2001):
        """The area of the minimal surface from z0, both halves: 4 pi times the integral of
        psi^4 e^q rho ds from the axis to the plane z = 0."""
        path, end = self.curve(z0)
        s = np.linspace(0, end, n)
        y = path(s)
        f = np.array([self.psi(r, z)[0] ** 4 * math.exp(q(self.a, r, z)[0]) * r for r, z in zip(y[0], y[1])])
        h = s[1] - s[0]
        return 4 * math.pi * h / 3 * (f[0] + f[-1] + 4 * f[1:-1:2].sum() + 2 * f[2:-1:2].sum())


@functools.lru_cache(maxsize=None)
def wave(a, nodes=NODES, modes=MODES):
    return Wave(a, nodes, modes)


def critical(lo=11.5, hi=12.2):
    """The amplitude at which a minimal surface first appears in the wave, 11.82."""
    return brentq(lambda a: Wave(a).deepest()[1], lo, hi, xtol=1e-6)
