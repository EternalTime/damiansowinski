#!/usr/bin/env python3
"""Bartnik and McKinnon's solitons, as numbers, in units of the length ell, G = c = 1.

The static, spherically symmetric SU(2) Einstein-Yang-Mills equations with the purely magnetic
potential of amplitude w(r) (R. Bartnik and J. McKinnon, Phys. Rev. Lett. 61, 141, 1988, as
M. S. Volkov and D. V. Gal'tsov restate them, Phys. Rep. 319, 1, 1999), on

    ds^2 = -sigma^2 N dt^2 + dr^2/N + r^2 dOmega^2,      N = 1 - 2m/r,

are, with a prime for d/dr and every length in units of ell,

    m' = N w'^2 + (1 - w^2)^2/(2 r^2),
    r^2 N w'' + (2m - (1 - w^2)^2/r) w' + w (1 - w^2) = 0,
    sigma'/sigma = 2 w'^2/r.

No solution is known in closed form. A regular centre leaves one free number, b:

    w = 1 - b r^2 + (3b^2/10 + 4b^3/5) r^4 + ...,      m = 2 b^2 r^3 - (8/5) b^3 r^5 + ...,

and an asymptotically flat end two, the mass M and a:

    w = +-(1 - a/r + 3a(a - 2M)/(4r^2) + ...),          m = M - a^2/r^3 + ...

A solution regular at both ends exists only for a discrete ladder of b, one for each number
n = 1, 2, 3, ... of zeros of w. For any other b the amplitude leaves the strip |w| < 1 and the
solution ends, and the number of zeros it has by then steps up by one as b passes each b_n, which
is how bracket() finds them.

A soliton is then solved from both ends at once, since the end at infinity is unstable when it is
integrated outward (the perturbation of w that grows is r^2): from the centre outward with b, from
far out inward with M and a, the two meeting in w, w' and m just beyond the last zero of w, which
Newton's method brings about in (b, M, a). Beyond the far radius the series stands for the
solution, and inside R_NEAR the series at the centre.

Three more functions ride along: delta = -ln sigma, zero at infinity, from delta' = -2 w'^2/r;
the isotropic radius rho, with d(ln rho)/dr = 1/(r sqrt N) and rho/r -> 1 at infinity, so that
ds^2 = -sigma^2 N dt^2 + (r/rho)^2 (drho^2 + rho^2 dOmega^2); and the tortoise coordinate x, with
dx/dr = 1/(sigma N) and x = 0 at the centre, so that the plane of t and x is conformally flat.

    /tmp/mfs-venv/bin/python _tools/derivations/bartnik_mckinnon.py
"""
import math

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve
from scipy.interpolate import CubicSpline

R_NEAR = 1e-3       # inside it the series at the centre stands for the solution: its next term is 1e-24
MEET, FAR = 1.5, 200.0      # where the two integrations meet, and where the series at infinity takes
                            # over, as multiples of the radius of the last zero of w: the zeros of the
                            # higher solitons stand a factor of about 40 apart, and the solution is
                            # near its flat end only beyond the last. The series' next term there is 1e-14.
GRID = (1e-6, 1e10)         # the radii the tables cover
TOL = dict(rtol=1e-13, atol=1e-15, method="DOP853", dense_output=True)


def centre(b, r):
    """(w, w', m) from the series at a regular centre, through r^6 in w and r^7 in m."""
    c4 = 3 * b ** 2 / 10 - 4 * b ** 3 / 5
    c6 = -b ** 3 * (128 * b ** 2 - 80 * b + 7) / 70
    m5 = -8 * b ** 3 / 5
    m7 = 24 * b ** 4 * (12 * b ** 2 - 14 * b + 7) / 175
    w = 1 - b * r ** 2 + c4 * r ** 4 + c6 * r ** 6
    dw = -2 * b * r + 4 * c4 * r ** 3 + 6 * c6 * r ** 5
    m = 2 * b ** 2 * r ** 3 + m5 * r ** 5 + m7 * r ** 7
    return w, dw, m


def far(M, a, r):
    """(w, w', m) from the series at an asymptotically flat end, through 1/r^5 in w and 1/r^7 in
    m, for the branch w -> +1; the other is its mirror image, w -> -w."""
    x = 1 / r
    cw = [1, -a, 3 * a * (a - 2 * M) / 4, -a * (48 * M ** 2 - 42 * M * a + 11 * a ** 2) / 20,
          a * (-1920 * M ** 3 + 2244 * M ** 2 * a - 1076 * M * a ** 2 + 193 * a ** 3) / 480,
          -a * (7680 * M ** 4 - 10800 * M ** 3 * a + 7172 * M ** 2 * a ** 2 - 2408 * M * a ** 3
                + 329 * a ** 4 - 400 * a ** 2) / 1120]
    cm = [M, 0, 0, -a ** 2, a ** 2 * (4 * a - 5 * M) / 2,
          -3 * a ** 2 * (68 * M ** 2 - 100 * M * a + 37 * a ** 2) / 40,
          a ** 2 * (-588 * M ** 3 + 1200 * M ** 2 * a - 833 * M * a ** 2 + 196 * a ** 3) / 60,
          -a ** 2 * (77376 * M ** 4 - 195888 * M ** 3 * a + 192140 * M ** 2 * a ** 2 - 85876 * M * a ** 3
                     + 14679 * a ** 4 - 1800 * a ** 2) / 4200]
    w = sum(c * x ** k for k, c in enumerate(cw))
    dw = -sum(k * c * x ** (k + 1) for k, c in enumerate(cw))
    m = sum(c * x ** k for k, c in enumerate(cm))
    return w, dw, m


def second(r, w, dw, m):
    """w'' from the Yang-Mills equation."""
    return -((2 * m - (1 - w * w) ** 2 / r) * dw + w * (1 - w * w)) / (r * (r - 2 * m))


def mass_rate(r, w, dw, m):
    """m' from the Hamiltonian constraint, the energy density of the field times r^2."""
    return (1 - 2 * m / r) * dw * dw + (1 - w * w) ** 2 / (2 * r * r)


def field(r, y):
    """The equations for (w, w', m), and for delta, ln(rho/r) and x where y carries them."""
    w, dw, m = y[:3]
    N = 1 - 2 * m / r
    out = [dw, second(r, w, dw, m), mass_rate(r, w, dw, m)]
    if len(y) > 3:
        out += [-2 * dw * dw / r, (1 / math.sqrt(N) - 1) / r, math.exp(y[3]) / N]
    return out


def zeros(b, r_max=1e4):
    """How many zeros w has before it leaves the strip |w| < 1 or a horizon forms, for the
    solution regular at the centre with w = 1 - b r^2: the number that steps at each b_n."""
    def leaves(r, y):
        return abs(y[0]) - (1 + 1e-7)

    def horizon(r, y):
        return 1 - 2 * y[2] / r - 1e-5

    def node(r, y):
        return y[0]
    leaves.terminal = horizon.terminal = True
    run = solve_ivp(field, (R_NEAR, r_max), centre(b, R_NEAR), events=[leaves, horizon, node],
                    rtol=1e-11, atol=1e-13, method="DOP853")
    return len(run.t_events[2])


def bracket(n):
    """b_n to a part in 10^12, by bisection on the number of zeros: at most n below it, more above."""
    lo, hi = 0.05, 0.706
    for _ in range(44):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if zeros(mid) <= n else (lo, mid)
    return 0.5 * (lo + hi)


def _centre_tail(b, r):
    """(delta - delta(0), ln(rho/r) - its value at the centre, x e^(-delta(0))) from the series at
    the centre: delta' = -2 w'^2/r, d ln(rho/r)/dr = (1/sqrt(N) - 1)/r with N = 1 - 4 b^2 r^2 +
    (16/5) b^3 r^4, and dx/dr = e^delta/N, whose first correction is of order r^4."""
    c4 = 3 * b ** 2 / 10 - 4 * b ** 3 / 5
    return (-4 * b ** 2 * r ** 2 + 8 * b * c4 * r ** 4,
            b ** 2 * r ** 2 + (6 * b ** 4 - 8 * b ** 3 / 5) * r ** 4 / 4,
            r + 0 * r ** 5)


class Soliton:
    """The soliton whose w has n zeros: b, M and a, and every function of it at any radius."""

    def __init__(self, n):
        self.n, self.sign = n, (-1) ** n
        self._kept = {}
        b = bracket(n)
        # The solution shot outward follows the soliton well past its last zero, which sets the
        # two radii, and there the series at infinity, solved for M and a, gives Newton its start.

        def node(r, y):
            return y[0]
        out = solve_ivp(field, (R_NEAR, 1e4), centre(b, R_NEAR), events=node, **TOL)
        last = float(out.t_events[0][n - 1])
        self.r_meet, self.r_far = MEET * last, FAR * last
        w, dw, m = out.sol(8 * last)
        M, a = fsolve(lambda p: [far(p[0], p[1], 8 * last)[0] - abs(w), far(p[0], p[1], 8 * last)[2] - m],
                      [m, (1 - abs(w)) * 8 * last], xtol=1e-13)
        self.b, self.M, self.a = self._newton(np.array([b, M, a]))
        self._solve()

    def _gap(self, p):
        b, M, a = p
        inner = solve_ivp(field, (R_NEAR, self.r_meet), centre(b, R_NEAR), **TOL).y[:, -1]
        w, dw, m = far(M, a, self.r_far)
        outer = solve_ivp(field, (self.r_far, self.r_meet), [self.sign * w, self.sign * dw, m], **TOL).y[:, -1]
        return inner - outer

    def _newton(self, p):
        for _ in range(40):
            gap = self._gap(p)
            J = np.empty((3, 3))
            for k in range(3):
                h = 1e-7 * max(abs(p[k]), 1e-3)
                step = np.zeros(3)
                step[k] = h
                J[:, k] = (self._gap(p + step) - self._gap(p - step)) / (2 * h)
            move = np.linalg.solve(J, gap)
            p = p - move * min(1.0, 0.1 / float(np.max(np.abs(move / p))))
            if np.max(np.abs(move)) < 1e-13:
                break
        self.gap = float(np.max(np.abs(self._gap(p))))
        return p

    def _solve(self):
        """The two halves with delta, ln(rho/r) and x riding along. Beyond the far radius delta
        is the integral of 2 w'^2/r, a^2/(2 r^4) to leading order, and the isotropic radius and
        the tortoise coordinate are Schwarzschild's; delta and ln(rho/r) take their constants
        from there. The inner half starts from the series at the centre with those constants
        left out, and x from the centre, where it is r/sigma(0)."""
        M, R = self.M, self.r_far
        w, dw, m = far(M, self.a, R)
        self.outer = solve_ivp(field, (R, self.r_meet),
                               [self.sign * w, self.sign * dw, m, self.a ** 2 / (2 * R ** 4),
                                math.log((R - M + math.sqrt(R * (R - 2 * M))) / (2 * R)), 0.0], **TOL)
        self.inner = solve_ivp(field, (R_NEAR, self.r_meet),
                               [*centre(self.b, R_NEAR), *_centre_tail(self.b, R_NEAR)], **TOL)
        at_in, at_out = self.inner.sol(self.r_meet), self.outer.sol(self.r_meet)
        self.delta0 = float(at_out[3] - at_in[3])           # delta at the centre, -ln sigma(0)
        self.iso0 = float(at_out[4] - at_in[4])             # ln(rho/r) at the centre
        self.x_meet = math.exp(self.delta0) * float(at_in[5])
        self.x_out = float(at_out[5])
        self.x_far = self.x_meet + float(self.outer.sol(R)[5]) - self.x_out
        self._tables()

    def exact_state(self, r):
        """(w, w', m, delta, ln(rho/r), x) at radii r, from the two integrations and the two series."""
        r = np.atleast_1d(np.asarray(r, dtype=float))
        out = np.empty((6, r.size))
        near, mid = r < R_NEAR, (r >= R_NEAR) & (r < self.r_meet)
        outer, beyond = (r >= self.r_meet) & (r <= self.r_far), r > self.r_far
        if near.any():
            rr = r[near]
            d, iso, x = _centre_tail(self.b, rr)
            out[:, near] = [*centre(self.b, rr), self.delta0 + d, self.iso0 + iso, math.exp(self.delta0) * x]
        if mid.any():
            y = self.inner.sol(r[mid])
            # The inner half carries delta without its constant, so its x lacks the factor 1/sigma(0).
            out[:, mid] = [y[0], y[1], y[2], y[3] + self.delta0, y[4] + self.iso0, math.exp(self.delta0) * y[5]]
        if outer.any():
            y = self.outer.sol(r[outer])
            y[5] += self.x_meet - self.x_out
            out[:, outer] = y
        if beyond.any():
            rr, M, R = r[beyond], self.M, self.r_far
            w, dw, m = far(M, self.a, rr)
            star = rr - R + 2 * M * np.log((rr - 2 * M) / (R - 2 * M))
            out[:, beyond] = [self.sign * w, self.sign * dw, m, self.a ** 2 / (2 * rr ** 4),
                              np.log((rr - M + np.sqrt(rr * (rr - 2 * M))) / (2 * rr)), self.x_far + star]
        return out

    def exact_jet(self, r):
        """m, delta and N = 1 - 2m/r with their first two derivatives along r, and w and w', at
        radii r: a dict of arrays. The derivatives are the field equations' own, m' = N w'^2 +
        (1 - w^2)^2/(2 r^2) and delta' = -2 w'^2/r, differentiated once more with w'' from the
        Yang-Mills equation; inside R_NEAR they are the derivatives of the series, term by term,
        so that nothing is divided by a small radius."""
        r = np.atleast_1d(np.asarray(r, dtype=float))
        w, dw, m, delta = self.exact_state(r)[:4]
        b = self.b
        near = r < R_NEAR
        rr = np.where(near, 1.0, r)
        N = 1 - 2 * m / rr
        V = (1 - w * w) ** 2 / (2 * rr * rr)
        d2w = second(rr, w, dw, m)
        dm = N * dw * dw + V
        dN = 2 * (m - rr * dm) / (rr * rr)
        dV = -2 * w * dw * (1 - w * w) / (rr * rr) - 2 * V / rr
        d2m = dN * dw * dw + 2 * N * dw * d2w + dV
        d2N = -2 * d2m / rr - 2 * dN / rr
        dd = -2 * dw * dw / rr
        d2d = -4 * dw * d2w / rr + 2 * dw * dw / (rr * rr)
        if near.any():
            x = r[near]
            c4 = 3 * b ** 2 / 10 - 4 * b ** 3 / 5
            m5 = -8 * b ** 3 / 5
            m7 = 24 * b ** 4 * (12 * b ** 2 - 14 * b + 7) / 175
            dm[near] = 6 * b ** 2 * x ** 2 + 5 * m5 * x ** 4 + 7 * m7 * x ** 6
            d2m[near] = 12 * b ** 2 * x + 20 * m5 * x ** 3 + 42 * m7 * x ** 5
            # N = 1 - 2m/r = 1 - 4 b^2 r^2 - 2 m5 r^4 - 2 m7 r^6.
            N[near] = 1 - 4 * b ** 2 * x ** 2 - 2 * m5 * x ** 4 - 2 * m7 * x ** 6
            dN[near] = -8 * b ** 2 * x - 8 * m5 * x ** 3 - 12 * m7 * x ** 5
            d2N[near] = -8 * b ** 2 - 24 * m5 * x ** 2 - 60 * m7 * x ** 4
            dd[near] = -8 * b ** 2 * x + 32 * b * c4 * x ** 3
            d2d[near] = -8 * b ** 2 + 96 * b * c4 * x ** 2
        return {"w": w, "dw": dw, "m": m, "dm": dm, "d2m": d2m, "delta": delta, "ddelta": dd, "d2delta": d2d,
                "N": N, "dN": dN, "d2N": d2N}

    # ---- the functions a drawing names, each at radii r

    def w(self, r):
        return self.state(r)[0]

    def dw(self, r):
        return self.state(r)[1]

    def m(self, r):
        return self.state(r)[2]

    def delta(self, r):
        return self.state(r)[3]

    def lapse(self, r):
        """sigma^2 N = -g_tt."""
        j = self.jet(r)
        return np.exp(-2 * j["delta"]) * j["N"]

    def isotropic(self, r):
        """The isotropic radius rho of the sphere of areal radius r."""
        return np.asarray(r, dtype=float) * np.exp(self.state(r)[4])

    def tortoise(self, r):
        """The tortoise coordinate x of the sphere of areal radius r."""
        return self.state(r)[5]

    def _tables(self):
        """Every function on a grid even in ln r, from 1e-6 to 1e10 ell, as cubic splines: a drawing
        asks for them point by point, hundreds of thousands of times, and each of the integrations'
        answers costs several times a spline's. The spacing is a part in 900, which leaves the
        splines within a part in 1e7 of the integrations, as the tests hold; below the grid the series at the centre answers, and the two inverse maps
        start from splines of r against rho and against x."""
        u = np.linspace(math.log(GRID[0]), math.log(GRID[1]), 33001)
        # The neck about r = ell narrows as zeros are added, to a width of a few thousandths of
        # ell for three, so the grid is twenty times finer from ell/2 to 4 ell.
        u = np.unique(np.concatenate([u, np.linspace(math.log(0.5), math.log(4.0), 40001)]))
        r = np.exp(u)
        y, j = self.exact_state(r), self.exact_jet(r)
        self._names = ("w", "dw", "m", "dm", "d2m", "delta", "ddelta", "d2delta", "N", "dN", "d2N")
        self._spline = CubicSpline(u, np.array([j[k] for k in self._names] + [y[4], y[5]]), axis=1)
        self._from_rho = CubicSpline(u + y[4], u)
        self._from_x = CubicSpline(np.log(y[5]), u)

    def _read(self, r):
        """The eleven functions of the jet, ln(rho/r) and x at radii r, from the tables."""
        r = np.atleast_1d(np.asarray(r, dtype=float))
        # One point of a drawing asks for several of these functions at one radius, so the
        # last answer is kept.
        kept = self._kept.get("read")
        if kept is not None and kept[0].shape == r.shape and np.array_equal(kept[0], r):
            return kept[1]
        out = np.empty((13, r.size))
        self._kept["read"] = (r.copy(), out)
        low = r < GRID[0]
        if low.any():
            j, y = self.exact_jet(r[low]), self.exact_state(r[low])
            out[:, low] = [j[k] for k in self._names] + [y[4], y[5]]
        if not low.all():
            out[:, ~low] = self._spline(np.log(r[~low]))
        return out

    def state(self, r):
        """(w, w', m, delta, ln(rho/r), x) at radii r."""
        v = self._read(r)
        return v[[0, 1, 2, 5, 11, 12]]

    def jet(self, r):
        """What exact_jet returns, read from the tables."""
        return dict(zip(self._names, self._read(r)[:11]))

    def from_isotropic(self, rho):
        """The areal radius of the sphere of isotropic radius rho, polished by Newton's method."""
        rho = np.atleast_1d(np.asarray(rho, dtype=float))
        kept = self._kept.get("from_isotropic")
        if kept is not None and kept[0].shape == rho.shape and np.array_equal(kept[0], rho):
            return kept[1]
        low = rho < GRID[0] * math.exp(self.iso0)
        r = np.where(low, rho * math.exp(-self.iso0), np.exp(self._from_rho(np.log(np.where(low, 1.0, rho)))))
        for _ in range(2):
            v = self._read(r)
            r = np.maximum(r - (r * np.exp(v[11]) - rho) * np.exp(-v[11]) * np.sqrt(v[8]), 0.0)
        self._kept["from_isotropic"] = (rho.copy(), r)
        return r

    def from_tortoise(self, x):
        """The areal radius of the sphere of tortoise coordinate x, polished by Newton's method."""
        x = np.atleast_1d(np.asarray(x, dtype=float))
        kept = self._kept.get("from_tortoise")
        if kept is not None and kept[0].shape == x.shape and np.array_equal(kept[0], x):
            return kept[1]
        low = x < GRID[0] * math.exp(self.delta0)
        r = np.where(low, x * math.exp(-self.delta0), np.exp(self._from_x(np.log(np.where(low, 1.0, x)))))
        for _ in range(2):
            v = self._read(r)
            r = np.maximum(r - (v[12] - x) * np.exp(-v[5]) * v[8], 0.0)
        self._kept["from_tortoise"] = (x.copy(), r)
        return r

    def stretch(self, rho, order=0):
        """k = r/rho along the isotropic radius, which is finite at the centre, or with order 1
        or 2 its first or second derivative along rho: k' = k (sqrt(N) - 1)/rho, written
        k' = -2 (m/r^2) k^2/(1 + sqrt(N)) so that nothing is divided by a small radius."""
        rho = np.atleast_1d(np.asarray(rho, dtype=float))
        r = self.from_isotropic(rho)
        k = np.exp(-self.state(r)[4])
        if order == 0:
            return k
        j = self.jet(r)
        root = np.sqrt(j["N"])
        g = (1 - j["N"]) / 2            # m/r
        safe = np.where(r < R_NEAR, 1.0, r)
        u = np.where(r < R_NEAR, self._m_over(r, 2), g / safe)                      # m/r^2
        du = np.where(r < R_NEAR, self._m_over(r, 2, 1), (j["dm"] - 2 * g) / safe ** 2)  # (m/r^2)'
        dk = -2 * u * k * k / (1 + root)
        if order == 1:
            return dk
        # d/drho = (dr/drho) d/dr with dr/drho = k sqrt(N).
        return -2 * (k * root * du * k * k / (1 + root) + 2 * u * k * dk / (1 + root)
                     - u * k * k * k * root * j["dN"] / (2 * root * (1 + root) ** 2))

    def _m_over(self, r, power, order=0):
        """m/r^power from the series at the centre, power 2 or 3, or its derivative along r."""
        b = self.b
        c = [2 * b ** 2, -8 * b ** 3 / 5, 24 * b ** 4 * (12 * b ** 2 - 14 * b + 7) / 175]
        e = [3 - power, 5 - power, 7 - power]
        if order == 0:
            return sum(ck * r ** ek for ck, ek in zip(c, e))
        return sum(ck * ek * r ** (ek - 1) for ck, ek in zip(c, e) if ek > 0)

    # ---- what the field equations leave to check

    def residual(self, r):
        """The Yang-Mills equation, the two constraints and w' along the solution, by central
        differences of the solved functions over a step of r/10^4: each vanishes, to the square
        of that step, if the functions solve what they were built from."""
        r = np.asarray(r, dtype=float)
        h = 1e-4 * r
        up, at, down = self.exact_state(r + h), self.exact_state(r), self.exact_state(r - h)
        d = (up - down) / (2 * h)
        w, dw, m, delta = at[:4]
        return (d[1] - second(r, w, dw, m), d[2] - mass_rate(r, w, dw, m), d[3] + 2 * dw * dw / r, d[0] - dw)

    def facts(self):
        """The numbers the page quotes of this soliton."""
        r = np.geomspace(1e-3, 1e3 * self.r_meet, 400001)
        y = self.exact_state(r)
        N = 1 - 2 * y[2] / r
        k = int(np.argmin(N))
        nodes = r[:-1][np.sign(y[0][:-1]) != np.sign(y[0][1:])]
        return {"n": self.n, "b": float(self.b), "M": float(self.M), "a": float(self.a),
                "sigma0": math.exp(-self.delta0), "N_min": float(N[k]), "r_N_min": float(r[k]),
                "nodes": [float(v) for v in nodes], "rho_over_r_centre": math.exp(self.iso0)}


_SOLITONS = {}


def soliton(n=1):
    """The soliton with n zeros, solved once."""
    if n not in _SOLITONS:
        _SOLITONS[n] = Soliton(int(n))
    return _SOLITONS[n]


if __name__ == "__main__":
    for n in (1, 2, 3):
        s = soliton(n)
        f = s.facts()
        print(f"n = {n}: b = {f['b']:.10f}, M = {f['M']:.10f}, a = {f['a']:.8f}, sigma(0) = {f['sigma0']:.6f}, "
              f"least N = {f['N_min']:.6f} at r = {f['r_N_min']:.4f}, zeros of w at "
              + ", ".join(f"{v:.4f}" for v in f["nodes"]) + f", mismatch where the two halves meet {s.gap:.1e}")
        r = np.geomspace(2e-3, 2e3, 400)
        print("   largest residuals:", ", ".join(f"{float(np.max(np.abs(v))):.1e}" for v in s.residual(r)))
