"""Elvang and Figueras's black Saturn as numbers: the four functions of its polar chart, with
their first and second partial derivatives along r and theta.

The chart publishes its Kretschmann scalar in V, Omega, W and nu and their derivatives, and
leaves the four free, so a figure of that scalar needs them as numbers. Each is built from the
names Weyl's chart defines, e^(2V) = (G_y H_x^2 - (omega_psi + q H_y)^2)/(H_x H_y),
Omega = (omega_psi + q H_y) e^(-2V)/H_x, e^(2W) = G_y and e^(2nu) = k^2 H_x P, whose own
definitions reach down to the solitons mu_i of the rod ends. The polar chart's r and theta enter
through Weyl's rho = r^2 sin(theta) cos(theta)/L^2 and z = r^2 (cos^2 theta - sin^2 theta)/(2 L^2).
print_charts.black_saturn_functions reads the same four definitions for the chart printer's own
checks, and the curvature figure holds these functions to it.

Each derivative is carried through that arithmetic as a second order jet, a value together with
its two first and three second partial derivatives, so every one is exact to the precision the
arithmetic is done in. Jet.through chains an elementary function, and the arithmetic runs over
numpy arrays for a whole grid at once or over mpmath numbers at any precision.

black_saturn.md beside this file derives the solution, and null_rays.py declares these functions
to the readers of the published scalars as bs_V, bs_O, bs_W and bs_nu.
"""
from functools import lru_cache

import mpmath
import numpy as np
import sympy as sp

# Which component of a jet each order of derivative is, by the digits of the arguments it is taken
# along: 1 for the chart's first coordinate and 2 for its second.
ORDERS = ("", "1", "2", "11", "12", "22")


class Jet:
    """A value and its first and second partial derivatives along two variables, carried through
    arithmetic so that a composite's derivatives come out exact.

    v is the value, a and b the derivatives along the first and the second variable, and aa, ab
    and bb the second ones. Each part is a number or an array of them, so one jet carries a whole
    grid."""

    __slots__ = ("v", "a", "b", "aa", "ab", "bb")

    def __init__(self, v, a=0, b=0, aa=0, ab=0, bb=0):
        self.v, self.a, self.b, self.aa, self.ab, self.bb = v, a, b, aa, ab, bb

    def parts(self):
        return (self.v, self.a, self.b, self.aa, self.ab, self.bb)

    def __add__(self, other):
        if not isinstance(other, Jet):
            return Jet(self.v + other, self.a, self.b, self.aa, self.ab, self.bb)
        return Jet(*[x + y for x, y in zip(self.parts(), other.parts())])

    __radd__ = __add__

    def __neg__(self):
        return Jet(*[-x for x in self.parts()])

    def __sub__(self, other):
        return self + (-other)

    def __rsub__(self, other):
        return (-self) + other

    def __mul__(self, other):
        if not isinstance(other, Jet):
            return Jet(*[x * other for x in self.parts()])
        return Jet(self.v * other.v,
                   self.a * other.v + self.v * other.a,
                   self.b * other.v + self.v * other.b,
                   self.aa * other.v + 2 * self.a * other.a + self.v * other.aa,
                   self.ab * other.v + self.a * other.b + self.b * other.a + self.v * other.ab,
                   self.bb * other.v + 2 * self.b * other.b + self.v * other.bb)

    __rmul__ = __mul__

    def through(self, g, dg, ddg):
        """g of this jet, from g and its first two derivatives at the value."""
        v, first, second = g(self.v), dg(self.v), ddg(self.v)
        return Jet(v, first * self.a, first * self.b,
                   second * self.a * self.a + first * self.aa,
                   second * self.a * self.b + first * self.ab,
                   second * self.b * self.b + first * self.bb)

    def power(self, p):
        return self.through(lambda u: u ** p, lambda u: p * u ** (p - 1),
                            lambda u: p * (p - 1) * u ** (p - 2))

    def __pow__(self, p):
        return self.power(p)

    def __truediv__(self, other):
        return self * other.power(-1) if isinstance(other, Jet) else self * (1 / other)

    def __rtruediv__(self, other):
        return self.power(-1) * other


def jet_module(sqrt, log, exp, sin, cos, sign):
    """The elementary functions of a backend, each taking a jet as well as a number, as the
    modules lambdify evaluates an expression with."""
    def one(g, dg, ddg):
        return lambda u: u.through(g, dg, ddg) if isinstance(u, Jet) else g(u)
    return {
        "sqrt": one(sqrt, lambda u: 1 / (2 * sqrt(u)), lambda u: -1 / (4 * u * sqrt(u))),
        "log": one(log, lambda u: 1 / u, lambda u: -1 / u ** 2),
        "exp": one(exp, exp, exp),
        "sin": one(sin, cos, lambda u: -sin(u)),
        "cos": one(cos, lambda u: -sin(u), lambda u: -cos(u)),
        # The absolute value of a jet away from zero, where it is smooth: sympy writes the
        # sqrt(G_y/rho^2) of omega_psi with the absolute value of Weyl's rho.
        "Abs": lambda u: u * sign(u.v) if isinstance(u, Jet) else abs(u),
    }


NUMPY = jet_module(np.sqrt, np.log, np.exp, np.sin, np.cos, np.sign)
MPMATH = jet_module(mpmath.sqrt, mpmath.log, mpmath.exp, mpmath.sin, mpmath.cos, mpmath.sign)


@lru_cache(maxsize=4)
def polar(kappa_1, kappa_2, kappa_3, beta):
    """The polar chart's V, Omega, W and nu at the given parameters, each a pair of jet functions
    of Weyl's rho and z, over numpy arrays and over mpmath numbers."""
    import null_rays as nr
    _, _, reader = nr.load("black_saturn", "weyl")
    at = {reader.parameters[name]: sp.Rational(value) for name, value in
          (("kappa_1", kappa_1), ("kappa_2", kappa_2), ("kappa_3", kappa_3), ("beta", beta))}
    Hx, Hy, Gy, P, om, q, k = (sp.sympify(reader.parameters[name]).subs(at)
                               for name in ("H_x", "H_y", "G_y", "P", "omega_psi", "q", "k"))
    shift = om + q * Hy
    e2V = (Gy * Hx ** 2 - shift ** 2) / (Hx * Hy)
    exprs = {"V": sp.log(e2V) / 2, "Omega": shift / (Hx * e2V), "W": sp.log(Gy) / 2,
             "nu": sp.log(k ** 2 * Hx * P) / 2}
    held = {"H_x": Hx, "H_y": Hy, "G_y": Gy, "P": P, "omega_psi": om, "e2V": e2V, **exprs}
    rho, z = sp.Symbol("rho", real=True), sp.Symbol("z", real=True)
    return ({name: sp.lambdify((rho, z), expr, [NUMPY, "numpy"]) for name, expr in exprs.items()},
            {name: sp.lambdify((rho, z), expr, [MPMATH, "mpmath"]) for name, expr in exprs.items()},
            held)


def weyl_jets(r, theta, module, one, zero):
    """Weyl's rho and z as jets of the polar chart's r and theta, in units of L^2 at L = 1."""
    radius = Jet(r, one, zero, zero, zero, zero)
    angle = Jet(theta, zero, one, zero, zero, zero)
    sin, cos = module["sin"](angle), module["cos"](angle)
    square = radius * radius
    return square * sin * cos, square * (cos * cos - sin * sin) * (one / 2)


# The last two points or grids the four functions were worked out at. One reading of the published
# scalar asks for twenty-four of their values, and each costs the whole arithmetic, so the whole of
# it is kept by the point it was done at.
_LAST = {}


def values(parameters, r, theta, digits=None):
    """The polar chart's four functions and their first and second partial derivatives along r and
    theta, by name and then by the order ORDERS lists, at one point in `digits` digits or over a
    whole array of points in floating point."""
    kappas = tuple(str(parameters[name]) for name in ("kappa_1", "kappa_2", "kappa_3", "beta"))
    if digits is None:
        r, theta = np.asarray(r, dtype=float), np.asarray(theta, dtype=float)
        key = (kappas, digits, r.shape, r.tobytes(), theta.shape, theta.tobytes())
    else:
        key = (kappas, digits, repr(float(r)), repr(float(theta)))
    if key in _LAST:
        return _LAST[key]
    arrays, exact, _ = polar(*kappas)
    if digits is None:
        rho, z = weyl_jets(r, theta, NUMPY, 1.0, 0.0)
        got = {name: dict(zip(ORDERS, f(rho, z).parts())) for name, f in arrays.items()}
    else:
        with mpmath.workdps(digits):
            rho, z = weyl_jets(mpmath.mpf(r), mpmath.mpf(theta), MPMATH, mpmath.mpf(1), mpmath.mpf(0))
            got = {name: dict(zip(ORDERS, f(rho, z).parts())) for name, f in exact.items()}
    if len(_LAST) >= 2:
        _LAST.clear()
    _LAST[key] = got
    return got
