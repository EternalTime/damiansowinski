"""The droplets of Lin, Lunin and Maldacena's bubbling anti-de Sitter space in closed form, and the
field equation of type IIB supergravity they solve.

Lin, Lunin and Maldacena (hep-th/0409174), section 2: the metric
    ds^2 = -h^-2 (dt + V_i dx^i)^2 + h^2 (dy^2 + dx^i dx^i) + y e^G dOmega_3^2 + y e^-G dtildeOmega_3^2,
    h^-2 = 2 y cosh G,   z = tanh(G)/2,
is fixed by z, which solves a linear equation with z = -1/2 on the black droplets and +1/2 on the
white of the plane y = 0. Their x_1 and x_2 are written x and w here, and on the plane's polar
coordinates r and phi their V_phi is written V. The members drawn and checked:

- the disc of radius r0, their (2.13) and (2.14): z = (r^2 - r0^2 + y^2)/(2S) with
  S = sqrt((r^2 + r0^2 + y^2)^2 - 4 r^2 r0^2), and V = -((r^2 + y^2 + r0^2)/S - 1)/2;
- concentric circles, their (2.17): z = 1/2 + sum_i (-1)^(i+1) (z_disc(r0_i) - 1/2), the outermost
  circle first, and V the same alternating sum, so radii (r2, r1) with r2 > r1 are a black ring
  around a white hole;
- the half plane, their (2.9) and (2.10): z = w/(2 sqrt(w^2 + y^2)), V_1 = 1/(2 sqrt(w^2 + y^2)), V_2 = 0.

V_1 dx + V_2 dw = V dphi, so V_1 = -V w/r^2 and V_2 = V x/r^2. With these signs y d_y V_1 = -d_w z
and y d_y V_2 = d_x z, the orientation their epsilon_ij takes for their own solutions, and the
disc is anti-de Sitter space times a 5-sphere in the global chart with their tilde phi = phi + t;
their text prints tilde phi = phi - t, which goes with the opposite sign of V.
field_equation_residual holds a member to R_MN = F_MPQRS F_N^PQRS/96 with the 5-form built from
their (2.5) and (2.6).
"""
import functools
import itertools
import math

import sympy as sp

ANGLES = sp.symbols("_a1 _b1 _c1 _a2 _b2 _c2", real=True)
SPHERE_POINT = (0.7, 1.1, 0.3, 0.9, 0.5, 2.0)


# -- the members --------------------------------------------------------------------------

def _S(r2, y, r0):
    return sp.sqrt((r2 + r0 ** 2 + y ** 2) ** 2 - 4 * r2 * r0 ** 2)


def disc_z(r2, y, r0):
    """z of the black disc of radius r0, with r2 the square of the radius on the plane."""
    return (r2 - r0 ** 2 + y ** 2) / (2 * _S(r2, y, r0))


def disc_twist(r2, y, r0):
    """V/r^2 of the disc, -((r^2 + y^2 + r0^2)/S - 1)/(2r^2) written so that it is finite on the
    axis r = 0: the difference of the two squares is 4 r^2 r0^2."""
    S = _S(r2, y, r0)
    return -2 * r0 ** 2 / (S * (r2 + y ** 2 + r0 ** 2 + S))


def _concentric(r2, y, radii):
    z = sp.Rational(1, 2) + sum((-1) ** i * (disc_z(r2, y, r) - sp.Rational(1, 2)) for i, r in enumerate(radii))
    twist = sum((-1) ** i * disc_twist(r2, y, r) for i, r in enumerate(radii))
    return sp.atanh(2 * z), twist


def concentric_polar(r, y, radii):
    """(G, V) of black and white rings with the circles of the given radii, the outermost first and
    the region inside it black, on the plane's polar coordinates."""
    G, twist = _concentric(r ** 2, y, list(radii))
    return G, r ** 2 * twist


def concentric(x, w, y, radii):
    """(G, V_1, V_2) of the same rings on the plane's Cartesian coordinates."""
    G, twist = _concentric(x ** 2 + w ** 2, y, list(radii))
    return G, -twist * w, twist * x


def half_plane(x, w, y):
    """(G, V_1, V_2) of the black half plane w < 0."""
    s = sp.sqrt(w ** 2 + y ** 2)
    return sp.atanh(w / s), 1 / (2 * s), sp.Integer(0)


RING = (sp.sqrt(2), sp.Integer(1))      # the black ring drawn: radii sqrt(2) L^2 and L^2 at L = 1


def functions_text(radii=RING):
    """G and V of concentric circles as text in the plain names r and y, the form a drawing's row
    declares its functions in."""
    r, y = sp.symbols("r y", positive=True)
    G, V = concentric_polar(r, y, radii)
    return {"G": str(G), "V": str(V)}


def axis_inverse_speed(y, radii=RING):
    """dt/dy of a light ray along the axis r = 0, 1/(2y cosh G), from the alternating sum
    c = sum_i (-1)^i/(r_i^2 + y^2), the outermost circle first, which keeps it exact at y = 0:
    1 - 2z = 2 - 2y^2 c for a white centre (an even number of circles) and 2y^2 c for a black one."""
    radii = [float(r) for r in radii]
    c = sum((-1) ** i / (r * r + y * y) for i, r in enumerate(radii))
    if len(radii) % 2 == 0:
        return math.sqrt(-c * (1 + y * y * c))
    return math.sqrt(c * (1 - y * y * c))


def axis_tortoise(y, radii=RING):
    """y_* = the integral of axis_inverse_speed from 0 to y, so that a ray on the axis keeps
    t -+ y_* constant."""
    from scipy.integrate import quad
    return quad(axis_inverse_speed, 0, y, args=(tuple(radii),), epsabs=1e-13, epsrel=1e-13, limit=200)[0]


# -- the metric and the field equation ----------------------------------------------------

def metric(kind, base, G, V):
    """The ten-dimensional metric on `base` (t and three coordinates) and the two 3-spheres of
    ANGLES: Cartesian, base (t, x, w, y) and V the pair (V_1, V_2); polar, base (t, r, phi, y) and
    V the one function V."""
    t, a, b, y = base
    H = 2 * y * sp.cosh(G)
    shift = [1, V[0], V[1], 0] if kind == "cartesian" else [1, 0, V, 0]
    flat = [0, 1, 1 if kind == "cartesian" else a ** 2, 1]
    g = sp.zeros(10)
    for i in range(4):
        for j in range(4):
            g[i, j] = -H * shift[i] * shift[j]
        g[i, i] += flat[i] / H
    for k, (u, radius) in enumerate(((ANGLES[:3], y * sp.exp(G)), (ANGLES[3:], y * sp.exp(-G)))):
        unit = [1, sp.sin(u[0]) ** 2, sp.sin(u[0]) ** 2 * sp.sin(u[1]) ** 2]
        for m in range(3):
            g[4 + 3 * k + m, 4 + 3 * k + m] = radius * unit[m]
    return g


@functools.lru_cache(maxsize=None)
def generic(kind):
    """The metric, its inverse and its Ricci tensor R_mn = R^a_man with G and V left free, unsimplified."""
    base = sp.symbols("t a b y", real=True)
    args = base[1:]
    G = sp.Function("G")(*args)
    V = (sp.Function("V1")(*args), sp.Function("V2")(*args)) if kind == "cartesian" else sp.Function("V")(*args)
    g = metric(kind, base, G, V)
    X = list(base) + list(ANGLES)
    n = 10
    gi = sp.zeros(n)
    gi[:4, :4] = g[:4, :4].inv()
    for k in range(4, n):
        gi[k, k] = 1 / g[k, k]
    dg = [[[sp.diff(g[i, j], X[k]) for k in range(n)] for j in range(n)] for i in range(n)]
    low = [[[(dg[a][b][c] + dg[a][c][b] - dg[b][c][a]) / 2 for c in range(n)] for b in range(n)] for a in range(n)]
    up = [[[sum(gi[a, d] * low[d][b][c] for d in range(n)) for c in range(n)] for b in range(n)] for a in range(n)]
    ricci = sp.zeros(n)
    for b in range(n):
        for d in range(b, n):
            ricci[b, d] = ricci[d, b] = sum(
                sp.diff(up[a][b][d], X[a]) - sp.diff(up[a][b][a], X[d])
                + sum(up[a][a][k] * up[k][b][d] - up[a][d][k] * up[k][b][a] for k in range(n)) for a in range(n))
    names = (G,) + (V if kind == "cartesian" else (V,))
    return base, names, g, ricci


def two_forms(kind, base, G, V):
    """Lin, Lunin and Maldacena's two-forms (2.5) and (2.6) on the base, F = dB_t ^ (dt + V)
    + B_t dV + d hat B with B_t = -y^2 e^(2G)/4 and d hat B = -(y^3/4) *_3 d((z + 1/2)/y^2), and
    F tilde the same with B_t = -y^2 e^(-2G)/4 and z - 1/2, as antisymmetric matrices. The flat
    star *_3 of dy^2 + dx^2 + dw^2 is taken in the orientation dw ^ dx ^ dy, the one their
    epsilon_ij takes for the signs of V in their (2.10) and (2.14); with dx ^ dw = r dr ^ dphi that
    is -r dr ^ dphi ^ dy on the polar coordinates."""
    t, a, b, y = base
    shift = [1, V[0], V[1], 0] if kind == "cartesian" else [1, 0, V, 0]
    z = sp.tanh(G) / 2
    # *_3 of da, db and dy as (coefficient, first index, second index), in the orientation above.
    if kind == "cartesian":
        star = {1: (1, 3, 2), 2: (1, 1, 3), 3: (1, 2, 1)}
    else:
        star = {1: (-a, 2, 3), 2: (-1 / a, 3, 1), 3: (-a, 1, 2)}
    out = []
    for B, shifted in ((-y ** 2 * sp.exp(2 * G) / 4, z + sp.Rational(1, 2)),
                       (-y ** 2 * sp.exp(-2 * G) / 4, z - sp.Rational(1, 2))):
        F = sp.zeros(4)
        dB = [sp.diff(B, c) for c in base]
        for m in range(4):
            for n in range(4):
                F[m, n] += (dB[m] * shift[n] - dB[n] * shift[m]
                            + B * (sp.diff(shift[n], base[m]) - sp.diff(shift[m], base[n])))
        phi = shifted / y ** 2
        for k, (factor, i, j) in star.items():
            value = -y ** 3 / 4 * factor * sp.diff(phi, base[k])
            F[i, j] += value
            F[j, i] -= value
        out.append(F)
    return out


def field_equation_residual(kind, funcs, points):
    """The largest |R_MN - F_MPQRS F_N^PQRS/96| over every component at each point (a, b, y) of the
    base and SPHERE_POINT on the spheres, for the member whose G and V are `funcs` in the plain
    symbols of generic(kind)'s base, with the 5-form 4 (F ^ Omega_3 + F tilde ^ Omega tilde_3) of
    the two-forms above and the unit spheres' volume forms."""
    import numpy as np
    base, names, g, ricci = generic(kind)
    exprs = dict(zip(names, funcs))
    atoms = set().union(*(e.atoms(sp.Derivative) for e in ricci))
    V = tuple(funcs[1:]) if kind == "cartesian" else funcs[1]
    forms = two_forms(kind, base, funcs[0], V)
    worst = 0.0
    for at in points:
        point = {s: sp.Float(v, 40) for s, v in zip(base[1:], at)}
        values = {d: sp.diff(exprs[d.expr], *d.variables).subs(point).evalf(40) for d in atoms}
        values.update({f: sp.sympify(e).subs(point).evalf(40) for f, e in exprs.items()})
        values.update(point)
        values.update({s: sp.Float(v, 40) for s, v in zip(ANGLES, SPHERE_POINT)})
        R = np.array(ricci.xreplace(values).evalf(30).tolist(), dtype=float)
        gn = np.array(g.xreplace(values).evalf(30).tolist(), dtype=float)
        gi = np.linalg.inv(gn)
        F5 = np.zeros((10,) * 5)
        for F2, block, (s1, s2) in zip(forms, ((4, 5, 6), (7, 8, 9)), (SPHERE_POINT[:2], SPHERE_POINT[3:5])):
            F2 = np.array(F2.subs(point).evalf(30).tolist(), dtype=float)
            volume = 4 * math.sin(s1) ** 2 * math.sin(s2)
            for m in range(4):
                for n in range(m + 1, 4):
                    if F2[m, n] == 0:
                        continue
                    slots = (m, n) + block
                    for perm in itertools.permutations(range(5)):
                        sign = round(np.linalg.det(np.eye(5)[list(perm)]))
                        F5[tuple(slots[p] for p in perm)] += sign * F2[m, n] * volume
        up = np.einsum("abcde,bB,cC,dD,eE->aBCDE", F5, gi, gi, gi, gi)
        square = np.einsum("mbcde,nbcde->mn", F5, up)
        worst = max(worst, float(np.max(np.abs(R - square / 96))))
    return worst
