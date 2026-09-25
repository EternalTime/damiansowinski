#!/usr/bin/env python3
"""Check every published coordinate system in MFS/assets/data/metrics against sympy.

For each system the script reads the line element, builds g_{mu nu} from it, and then
computes the inverse metric, the Christoffel symbols in both variants, the Riemann
tensor, the Ricci tensor, the Ricci scalar, the Kretschmann scalar, the Einstein
tensor and the Weyl tensor. Every value the file publishes is compared against what
sympy got, and every component the file leaves out is required to vanish, so an
omission is caught as well as a wrong number.

Every published expression is also checked for dimensional consistency, which needs no
sympy algebra and so runs in a moment over the whole collection. See below.

Nothing passes silently. A system whose values cannot be parsed, or whose declaration
is missing, or a tensor sympy cannot finish inside the time budget, is reported as
UNCHECKED with the reason. The script exits non-zero if anything disagreed.

    python3 -m venv /tmp/mfs-venv && /tmp/mfs-venv/bin/pip install sympy
    /tmp/mfs-venv/bin/python _tools/derivations/verify_metrics.py

Pass --system <metric_id>/<system_id> to check one system, repeatable.
Pass --budget <seconds> to change the per tensor time budget, which defaults to 120.
Pass --dimensions-only to run the dimensional pass alone and skip the sympy algebra.


The chart convention
--------------------

The collection writes x^0 = cT. A coordinate that carries dimensions of time is
therefore not itself the chart coordinate: the chart coordinate is c times it, and
every published component is a component in that chart even though the index is
printed with the bare name. Schwarzschild is the clearest example, publishing
g_{tt} = -(1 - r_s/r) against a line element whose time term is -(1 - r_s/r)c^2dt^2.

Because the rescaling x^0 -> c x^0 is linear with constant coefficients, a component
in the chart is the component computed with the bare coordinate multiplied by

    c^(number of upper time indices - number of lower time indices),

and the Christoffel symbols follow the same rule as the tensors, since the
inhomogeneous term in their transformation law carries a second derivative of the
coordinate change and so vanishes for a linear one.

DIMENSIONS below declares, per system, the dimension of every coordinate and every
parameter. It is the one thing the script cannot read off the file, because telling a
time from a length needs the dimensions of the parameters as well, and those are prose.
A system that is not declared there is reported UNCHECKED rather than guessed at. A
coordinate declared as a time is a coordinate the chart multiplies by c; every other
coordinate is already its own chart coordinate.


The dimensional pass
--------------------

Every published expression has a dimension its left hand side fixes, and every term of
it has to carry that same dimension. With x^0 = cT the chart coordinates are the ones
DIMENSIONS declares, except that a time is multiplied by c and so becomes a length, and
then

    [g_{mu nu}] = L^2 / ([x^mu][x^nu]),

with an upper index contributing [x^mu]/L and a lower one L/[x^mu], on top of the
dimension the field carries when every coordinate is a length: 1 for the metric, 1/L
for either Christoffel variant, 1/L^2 for Riemann, Ricci, Einstein, Weyl and the Ricci
scalar, and 1/L^4 for Kretschmann. A geodesic equation is measured against its own
second derivative, so its terms carry [x^mu] over the affine parameter squared, and the
dots in it are chart velocities like everywhere else.

An argument of exp, log or a trigonometric function has to be dimensionless, which is
what pins Godel's coordinates down. A power whose exponent is not a number contributes
only the dimension of its numeric part, because such a power is read with its base in a
fixed unit; the Kasner entry says so of its own t^{2p_i} in as many words.


The curvature conventions
-------------------------

The signature is (-,+,+,+) and the Riemann tensor is

    R^mu_{nu rho sigma} = d_rho Gamma^mu_{nu sigma} - d_sigma Gamma^mu_{nu rho}
                          + Gamma^mu_{rho lam} Gamma^lam_{nu sigma}
                          - Gamma^mu_{sigma lam} Gamma^lam_{nu rho},

which is what the published Riemann components are in.

The Ricci tensor is the standard contraction R_{mu nu} = R^a_{mu a nu}, on the first
lower index, which is the one the signature (-,+,+,+) asks for: it makes the Einstein
tensor of ordinary matter positive where the energy density is, so FRW publishes
G_{tt} = 3(adot^2 + k)/a^2 and the scalar +6(addot/a + adot^2/a^2 + k/a^2). The
collection once contracted on the last index instead, which is minus this, and the
change of convention was settled on 2026-09-18.

The Weyl tensor is built from the same contraction, as it always was, because it is
defined by removing the traces of Riemann and those traces are Riemann's own.


Constrained parameters
----------------------

Some entries carry parameters that are not free. Kasner prints three exponents bound
by sum p_i = sum p_i^2 = 1, and claims its values only on the surface those two
equations cut out: the Ricci tensor it publishes as zero is not zero for arbitrary
exponents. Checking such an entry against free symbols would test a stronger claim
than it makes and report a disagreement that is not one.

PARAMETER_RELATIONS below carries, per system, a rational parametrisation of that
surface. Every constrained parameter is replaced by its parametrised value on both
sides of each comparison, so what is checked is an identity along the surface. That is
exact rather than a sample: a rational parametrisation of an irreducible variety
covers a dense subset of it, so an identity in the parameter is an identity on the
whole surface. The parametrising symbol must not be a name the system already uses.
"""

import argparse
import json
import re
import signal
import sys
from pathlib import Path

import sympy as sp
from sympy.polys.polyerrors import BasePolynomialError
from sympy.polys.rings import PolyRing
from sympy.parsing.sympy_parser import (
    implicit_multiplication,
    parse_expr,
    split_symbols_custom,
    standard_transformations,
)

ROOT = Path(__file__).resolve().parents[2]
METRICS_DIR = ROOT / "MFS" / "assets" / "data" / "metrics"

DEFAULT_BUDGET_SECONDS = 120

LENGTH = sp.Symbol("L", positive=True)
TIME = sp.Symbol("T", positive=True)
MASS = sp.Symbol("M", positive=True)
AFFINE = sp.Symbol("lambda", positive=True)
BASE_DIMENSIONS = {"L": LENGTH, "T": TIME, "M": MASS, "1": sp.Integer(1)}

# The dimension of every coordinate and every parameter, per system. A coordinate
# declared T is the one the chart multiplies by c; the rest are their own chart
# coordinates. A system absent from this table is UNCHECKED.
DIMENSIONS = {
    # Both parameters are dimensionless: v_s is the bubble velocity in units of c, so that
    # c v_s f is the shift, and f is the shape function. The solution names no length of its
    # own, because the scale f varies on is left to f.
    ("alcubierre", "cartesian"): {
        "t": "T", "x": "L", "y": "L", "z": "L", "v_s": "1", "f": "1",
    },
    # The one length of the solution is the anti-de Sitter radius, which the entry calls
    # L and the dimensional pass calls L as well; the two never meet in one expression.
    ("anti_de_sitter", "static_global"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "L": "L",
    },
    ("anti_de_sitter", "poincare"): {
        "t": "T", "x": "L", "y": "L", "z": "L", "L": "L",
    },
    # b is the one length of the solution, the common radius of the two factors, and
    # neither r nor x is an areal radius: both run along the AdS_2 factor.
    ("bertotti_robinson", "static"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "b": "L",
    },
    ("bertotti_robinson", "poincare"): {
        "t": "T", "x": "L", "\\theta": "1", "\\phi": "1", "b": "L",
    },
    ("bianchi", "type_i_cartesian"): {
        "t": "T", "x": "L", "y": "L", "z": "L", "a_1": "1", "a_2": "1", "a_3": "1",
    },
    # The exterior keeps G and a mass per unit length explicit, so that the deficit is
    # the dimensionless 4G mu/c^2 the entry prints; delta is the deficit angle itself,
    # which the entry declares and states rather than uses. The interior is written with
    # a dimensionless polar angle chi, so its one length is the radius of the cap.
    ("cosmic_string", "conical"): {
        "t": "T", "r": "L", "\\phi": "1", "z": "L",
        "\\mu": "M/L", "G": "L**3/(M*T**2)", "\\delta": "1",
    },
    ("cosmic_string", "interior_cap"): {
        "t": "T", "\\chi": "1", "\\phi": "1", "z": "L",
        "\\ell": "L", "\\chi_0": "1", "\\rho": "M/L**3", "\\mu": "M/L",
        "G": "L**3/(M*T**2)",
    },
    # The static patch carries the cosmological constant itself, a curvature; the flat
    # slicing carries the Hubble rate instead, a frequency, with 3H^2/c^2 = Lambda.
    ("de_sitter", "static_spherical"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "\\Lambda": "1/L**2",
    },
    ("de_sitter", "flat_slicing"): {
        "t": "T", "x": "L", "y": "L", "z": "L", "H": "1/T",
    },
    ("ellis_bronnikov", "spherical"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "\\ell": "L",
    },
    # r is the comoving radial coordinate the entry calls it, a length, and a(t) is
    # dimensionless, which is the normalisation its published curvature obeys; k is then a
    # curvature, carrying 1/L^2, and the values -1, 0 and +1 the entry lists for it are in
    # units of the curvature radius. That is why the closed case's domain ends at
    # r = 1/sqrt(k) rather than at 1.
    ("frw", "comoving_spherical"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "a": "1", "k": "1/L**2",
    },
    ("frw", "conformal_spherical"): {
        "\\eta": "L", "r": "L", "\\theta": "1", "\\phi": "1", "a": "1", "k": "1/L**2",
    },
    # e^x forces x dimensionless, and with it the other three coordinates, so the whole
    # length of the Godel solution sits in 1/omega.
    ("godel", "cartesian"): {
        "t": "1", "x": "1", "y": "1", "z": "1", "\\omega": "1/L",
    },
    ("interior_schwarzschild", "spherical"): {
        "t": "L", "r": "L", "\\theta": "1", "\\phi": "1", "r_s": "L", "R": "L",
    },
    ("kasner", "cartesian"): {
        "t": "T", "x": "L", "y": "L", "z": "L", "p_1": "1", "p_2": "1", "p_3": "1",
    },
    # The other entry that keeps G and a mass explicit rather than folding them into a
    # length. The spin per unit mass a = J/(Mc) is a length, which is what makes
    # r^2 + a^2cos^2(theta) and r^2 - 2GMr/c^2 + a^2 areas.
    ("kerr", "boyer_lindquist"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1",
        "M": "M", "a": "L", "G": "L**3/(M*T**2)",
    },
    # Kerr with the charge folded into a length beside the spin, as Reissner-Nordstrom
    # folds it: r_Q^2 = GQ^2/(4 pi epsilon_0 c^4) is an area, which is what lets it sit
    # in Delta beside r^2 and a^2. The charge itself never appears, so no dimension for
    # it is needed, and the mass stays a mass with G beside it as it does for Kerr.
    ("kerr_newman", "boyer_lindquist"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1",
        "M": "M", "a": "L", "G": "L**3/(M*T**2)", "r_Q": "L",
    },
    # The one parameter is the dimensionless shape function of the tube, and the chart is
    # cylindrical about the axis the tube is laid along, so r is the distance from that axis
    # rather than an areal radius. The solution names no length of its own: the radius of the
    # tube and the thickness of its wall are both left to k.
    ("krasnikov", "cylindrical"): {
        "t": "T", "x": "L", "r": "L", "\\phi": "1", "k": "1",
    },
    ("minkowski", "cartesian"): {"t": "T", "x": "L", "y": "L", "z": "L"},
    ("minkowski", "spherical"): {"t": "T", "r": "L", "\\theta": "1", "\\phi": "1"},
    ("minkowski", "double_null"): {"u": "T", "v": "T", "y": "L", "z": "L"},
    ("minkowski", "spherical_null"): {"u": "T", "v": "T", "\\theta": "1", "\\phi": "1"},
    ("minkowski", "rindler"): {"T": "T", "X": "L", "Y": "L", "Z": "L", "a": "L/T**2"},
    # The redshift function sits inside an exponential and so is dimensionless, and the
    # shape function is a length beside r, which is what leaves 1 - b/r dimensionless.
    # In the proper distance chart l is the radial coordinate and the areal radius r is
    # a declared function of it rather than a coordinate.
    ("morris_thorne", "spherical"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1",
        "\\Phi": "1", "b": "L", "b_0": "L",
    },
    ("morris_thorne", "proper_radial"): {
        "t": "T", "l": "L", "\\theta": "1", "\\phi": "1",
        "\\Phi": "1", "r": "L", "b_0": "L",
    },
    # The slices are flat space and the whole of the geometry is the flow field carried on them,
    # so its components are velocities and the chart components of the metric are powers of V/c.
    ("natario", "cartesian_flow"): {
        "t": "T", "x": "L", "y": "L", "z": "L",
        "u": "L/T", "v": "L/T", "w": "L/T",
    },
    ("natario", "plane_flow"): {
        "t": "T", "x": "L", "y": "L", "z": "L", "u": "L/T",
    },
    # The collapse is two charts. Inside, the comoving polar angle chi is dimensionless
    # and the scale factor carries the length, so an areal radius is a sin(chi) and a dot
    # on a is dimensionless; chi_0 marks the surface and a_m is the scale factor at
    # release, both declared because the entry states the matching in them. Outside it is
    # Schwarzschild, whose one length is r_s = 2GM/c^2 = a_m sin^3(chi_0).
    ("oppenheimer_snyder", "interior_comoving"): {
        "\\tau": "T", "\\chi": "1", "\\theta": "1", "\\phi": "1",
        "a": "L", "\\chi_0": "1", "a_m": "L",
    },
    ("oppenheimer_snyder", "exterior_schwarzschild"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "r_s": "L",
    },
    # u is the retarded time and v the affine parameter along the rays, which is a
    # length, so the wave profile H is dimensionless and the amplitudes of the exact
    # plane wave, multiplying x^2, are curvatures.
    ("pp_wave", "brinkmann"): {
        "u": "T", "v": "L", "x": "L", "y": "L", "H": "1",
    },
    ("pp_wave", "exact_plane_wave"): {
        "u": "T", "v": "L", "x": "L", "y": "L", "A": "1/L**2", "B": "1/L**2",
    },
    ("rn_metric", "spherical"): {
        "t": "L", "r": "L", "\\theta": "1", "\\phi": "1", "r_s": "L", "r_q": "L",
    },
    ("schwarzschild", "spherical"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "r_s": "L",
    },
    ("schwarzschild", "eddington_finkelstein_outgoing"): {
        "u": "L", "r": "L", "\\theta": "1", "\\phi": "1", "r_s": "L",
    },
    ("schwarzschild", "eddington_finkelstein_ingoing"): {
        "v": "L", "r": "L", "\\theta": "1", "\\phi": "1", "r_s": "L",
    },
    ("stockum_dust", "cylindrical"): {
        "t": "L", "r": "L", "\\phi": "1", "z": "L", "R": "L",
    },
    # The mass is folded into the length m = GM/c^2 and the NUT parameter is a length
    # beside it, so the cross term 2l cos(theta) carries the one length g_{t phi} wants.
    ("taub_nut", "spherical"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "m": "L", "l": "L",
    },
    # The areal radius is a length and the comoving shell label r is a length beside it,
    # so \partial_r R is dimensionless and the energy function has to be dimensionless
    # as well, which is what leaves 1 + 2E a pure number. Every published derivative of R
    # is taken along a chart coordinate, so each order of it carries one inverse length.
    ("tolman_bondi", "comoving_synchronous"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "R": "L", "E": "1",
    },
    # The one entry that keeps G and a mass explicit rather than folding them into a
    # length like r_s, and so the only one whose declarations need a mass at all.
    ("vaidya", "eddington_finkelstein_outgoing"): {
        "u": "T", "r": "L", "\\theta": "1", "\\phi": "1",
        "G": "L**3/(M*T**2)", "m": "M",
    },
    ("vaidya", "eddington_finkelstein_ingoing"): {
        "v": "T", "r": "L", "\\theta": "1", "\\phi": "1",
        "G": "L**3/(M*T**2)", "m": "M",
    },
}

# What a field carries when every coordinate is a length; the indices supply the rest.
FIELD_DIMENSIONS = {
    "metric_components": sp.Integer(1),
    "inverse_metric_components": sp.Integer(1),
    "christoffel": 1 / LENGTH,
    "riemann": 1 / LENGTH ** 2,
    "ricci_tensor": 1 / LENGTH ** 2,
    "ricci_scalar": 1 / LENGTH ** 2,
    "einstein_tensor": 1 / LENGTH ** 2,
    "kretschmann": 1 / LENGTH ** 4,
    "weyl_tensor": 1 / LENGTH ** 2,
}


def time_coordinates(declared, coords):
    """The coordinates the chart multiplies by c, which are exactly the times."""
    return {name for name in coords
            if sp.sympify(declared[name], locals=BASE_DIMENSIONS) == TIME}

# A rational parametrisation of the surface an entry's constrained parameters live on,
# written as {parameter name: expression in a fresh symbol}. See the header.
PARAMETER_RELATIONS = {
    # The Kasner circle, where the plane sum p_i = 1 cuts the sphere sum p_i^2 = 1.
    # Every point of it is reached, and the exponent ordering p_1 <= p_2 <= p_3 holds
    # on u >= 1, which is the range Belinskii, Khalatnikov and Lifshitz bounce within.
    ("kasner", "cartesian"): {
        "p_1": "-u/(1 + u + u**2)",
        "p_2": "(1 + u)/(1 + u + u**2)",
        "p_3": "u*(1 + u)/(1 + u + u**2)",
    },
}

GREEK = [
    "theta", "phi", "eta", "omega", "Omega", "ell", "pi", "lambda", "mu", "nu",
    "rho", "sigma", "tau", "chi", "psi", "alpha", "beta", "gamma", "delta",
    "epsilon", "kappa", "xi", "zeta", "Lambda", "Phi", "Theta", "Psi", "Sigma",
]

TRIG = ["sinh", "cosh", "tanh", "coth", "sin", "cos", "tan", "cot", "sec", "csc"]

FUNCTIONS = {
    "sin": sp.sin, "cos": sp.cos, "tan": sp.tan, "cot": sp.cot,
    "sec": sp.sec, "csc": sp.csc, "sinh": sp.sinh, "cosh": sp.cosh,
    "tanh": sp.tanh, "coth": sp.coth, "exp": sp.exp, "log": sp.log,
    "ln": sp.log, "sqrt": sp.sqrt,
}


def norm(expression):
    """A canonical form of the expression, in which one that vanishes is exactly zero.

    Every tensor is built through this, so it has to be fast as well as exact. What it
    does, and why, measured on Kerr in Boyer-Lindquist coordinates on 2026-09-22:

    - sympy's general simplify, which this used to call twice, ran Kerr's Riemann
      tensor for over ninety minutes of processor time without finishing.
    - sin and cos in terms of symbols, and then sympy's cancel, is exact for the
      zero test but takes 40s on a single Riemann component, because it multiplies
      every denominator of a sum together and then takes a gcd of the product.
    - The same in a field of rational functions, which cancels a gcd after every
      addition, takes 2 to 9s per component, and 192s for the Kretschmann scalar.
    - What is here keeps every denominator as a product of irreducible factors, so a
      common denominator is a maximum over multiplicities and cancelling is trial
      division by those few factors; no gcd is ever taken. Kerr's Riemann tensor takes
      0.5s, its Kretschmann scalar 1.4s, and the whole system 17s, including the
      comparison of every published value. sympy's inv() on the metric was then the
      slowest step left, which is why Geometry takes the adjugate instead.

    sin and cos are reduced by cos^2 = 1 - sin^2 rather than by a Weierstrass
    substitution or by rewriting through exp, because both of those raise the degree of
    every polynomial and neither prints back as anything a reader would recognise. A
    radical is split into the roots of the irreducible factors of its radicand, which
    assumes each factor is positive; _factored_root says why that is the right reading
    for this collection, and it is what lets the interior Schwarzschild radicals, which
    sympy would never combine unaided, cancel. Anything _canonical has no generators
    for, which is nothing in the collection as it stands, falls back to simplify.
    """
    if isinstance(expression, sp.MatrixBase):
        return expression.applyfunc(norm)
    expression = sp.sympify(expression)
    if expression.is_Number:
        return expression
    try:
        return _canonical(expression)
    except (BasePolynomialError, NotImplementedError, ZeroDivisionError):
        return _simplified(expression)


def _simplified(expression):
    """The general simplifier, for the rare expression _canonical has no generators for.

    sympy's simplify leaves forms such as sin(2x)tan(x) + cos(2x) - 1 standing, which is
    why it is followed by a second pass through expand_trig.
    """
    simplified = sp.simplify(expression)
    if simplified == 0:
        return sp.Integer(0)
    harder = sp.simplify(sp.expand_trig(sp.expand(simplified)))
    if harder == 0:
        return sp.Integer(0)
    return harder if sp.count_ops(harder) < sp.count_ops(simplified) else simplified


# Every function the reader knows that is not sin, cos or exp, in terms of those three.
_IN_SIN_COS_EXP = {
    sp.tan: lambda u: sp.sin(u) / sp.cos(u),
    sp.cot: lambda u: sp.cos(u) / sp.sin(u),
    sp.sec: lambda u: 1 / sp.cos(u),
    sp.csc: lambda u: 1 / sp.sin(u),
    sp.sinh: lambda u: (sp.exp(u) - sp.exp(-u)) / 2,
    sp.cosh: lambda u: (sp.exp(u) + sp.exp(-u)) / 2,
    sp.tanh: lambda u: (sp.exp(u) - sp.exp(-u)) / (sp.exp(u) + sp.exp(-u)),
    sp.coth: lambda u: (sp.exp(u) + sp.exp(-u)) / (sp.exp(u) - sp.exp(-u)),
}


def _canonical(expression):
    """The expression written as N / (f_1^k_1 ... f_m^k_m), in the one way it can be.

    The expression is read as a rational function of generators: every symbol, every
    sin(u) and cos(u), every exp(v) and power with a symbolic exponent after it has been
    split into integer powers of one exp or power per term of the exponent, the square
    root of every irreducible factor of a radicand, and every other function application.
    Two kinds of generator are algebraic rather than free: cos(u) wherever sin(u) is also
    present, since cos^2 = 1 - sin^2, and a square root w of b, since w^2 = b. The
    numerator N is a polynomial in the generators, and the denominator is kept as a
    product of powers of irreducible, normalised factors. _Fraction explains how that is
    kept canonical. A root above the square is left to the general simplifier.
    """
    for function, rewrite in _IN_SIN_COS_EXP.items():
        if expression.has(function):
            expression = expression.replace(function, rewrite)
    if any(not atom.args[0].is_Symbol for atom in expression.atoms(sp.sin, sp.cos)):
        expression = sp.expand_trig(expression)
    # The reader writes a mixed partial in the order the coordinates are listed and diff
    # writes it in its own; doit puts both in diff's, so they become one generator.
    if expression.has(sp.Derivative):
        expression = expression.replace(lambda x: isinstance(x, sp.Derivative), lambda x: x.doit())
    # Innermost first, so a radicand is itself in this form before it is factored.
    expression = expression.replace(
        lambda x: x.is_Pow and x.exp.is_Rational and not x.exp.is_Integer and not x.base.is_Number,
        _factored_root)

    generators = set()
    radicals = {}
    _collect_generators(expression, generators, radicals)
    algebraic = []
    for cosine in [g for g in generators if g.func is sp.cos]:
        if sp.sin(cosine.args[0]) in generators:
            algebraic.append((cosine, 1 - sp.sin(cosine.args[0]) ** 2))
    for base, q in radicals.items():
        if q != 2:
            raise NotImplementedError(f"a root of {base} above the square")
        root = sp.I if base == -1 else sp.sqrt(base)
        if not (root.is_Pow or root is sp.I):
            raise NotImplementedError(f"sqrt({base}) does not stay a power")
        generators.add(root)
        algebraic.append((root, base))
    # A radicand that holds another algebraic generator is reduced before that one is,
    # so the conjugates it leaves behind can still be reduced on the inner one.
    roots = [root for root, _ in algebraic]
    algebraic.sort(key=lambda pair: -sum(1 for root in roots if sp.sympify(pair[1]).has(root)))

    ring = PolyRing(sorted(generators, key=sp.default_sort_key), sp.QQ, "lex")
    index = {symbol: i for i, symbol in enumerate(ring.symbols)}
    reader = _FractionReader(ring, index)
    value = reader(expression)
    relations = [(index[root], reader(base)) for root, base in algebraic]
    for _ in range(2 * len(relations) + 2):
        if not any(value.holds(i) for i, _ in relations):
            break
        for i, base in relations:
            if value.holds(i):
                value = value.reduce(i, base)
    else:
        raise NotImplementedError("the algebraic generators do not settle")
    value = value.cancelled()
    if not value.numerator:
        return sp.Integer(0)
    return value.as_expr()


def _factored_root(power):
    """b^e as the product of f^(k e) over the irreducible factors f^k of b.

    That is an identity only where every factor is positive, and it is what lets
    sqrt((R - r_s)(R^3 - r^2 r_s)) and R^2 sqrt(1 - r_s/R) sqrt(1 - r^2 r_s/R^3) meet as
    one product of two generators. Every radicand the collection prints is positive factor
    by factor on the region its entry describes, the interior of a star or the outside of
    a horizon, so reading them that way is reading them where the entry is claimed.
    A factor normalised to the opposite sign of the one the entry means comes out with a
    factor of i, and so disagrees rather than agreeing by accident.
    """
    out = sp.Integer(1)
    for part, sign in zip(sp.fraction(norm(power.base)), (1, -1)):
        content, factors = sp.factor_list(part)
        out *= sp.Pow(content, sign * power.exp)
        for factor, multiplicity in factors:
            out *= sp.Pow(factor, sign * multiplicity * power.exp)
    return out


def _exponent_terms(power):
    """exp(2x + y) as [(exp(x), 2), (exp(y), 1)], and t**(2p - 1) as [(t**p, 2), (t, -1)]."""
    base, exponent = (sp.E, power.args[0]) if power.func is sp.exp else power.as_base_exp()
    terms = []
    for term in sp.Add.make_args(sp.expand(exponent)):
        coefficient, rest = term.as_coeff_Mul()
        if not coefficient.is_Rational:
            coefficient, rest = sp.Integer(1), term
        if rest == 1:
            terms.append((base, coefficient))
        else:
            terms.append((sp.exp(rest) if base is sp.E else sp.Pow(base, rest), coefficient))
    return terms


def _is_transcendental_power(expression):
    return expression.func is sp.exp or (expression.is_Pow and not expression.exp.is_Number)


def _collect_generators(expression, generators, radicals):
    """The generators a rational function of the expression is taken in."""
    if expression.is_Rational:
        return
    if expression.is_Add or expression.is_Mul:
        for argument in expression.args:
            _collect_generators(argument, generators, radicals)
    elif expression.is_Pow and expression.exp.is_Integer:
        _collect_generators(expression.base, generators, radicals)
    elif expression.is_Pow and expression.exp.is_Rational:
        base = expression.base
        radicals[base] = sp.ilcm(radicals.get(base, 1), expression.exp.q)
        _collect_generators(base, generators, radicals)
    elif expression is sp.I:
        radicals[sp.Integer(-1)] = 2
    elif _is_transcendental_power(expression):
        for generator, coefficient in _exponent_terms(expression):
            if not coefficient.is_Integer:
                raise NotImplementedError(f"{expression} is not an integer power of a generator")
            generators.add(generator)
    else:
        generators.add(expression)


class _FractionReader:
    """An expression, read into a _Fraction over the given ring."""

    def __init__(self, ring, index):
        self.ring = ring
        self.index = index
        self.factored = {}

    def generator(self, symbol):
        polynomial = self.ring.gens[self.index[symbol]]
        return _Fraction(self, polynomial, {})

    def constant(self, value):
        return _Fraction(self, self.ring(value), {})

    def __call__(self, expression):
        expression = sp.sympify(expression)
        if expression in self.index:
            return self.generator(expression)
        if expression.is_Rational:
            return self.constant(sp.QQ(expression.p, expression.q))
        if expression.is_Add:
            return _Fraction.sum([self(argument) for argument in expression.args], self)
        if expression.is_Mul:
            out = self.constant(1)
            for argument in expression.args:
                out = out * self(argument)
            return out
        if expression.is_Pow and expression.exp.is_Integer:
            return self(expression.base) ** int(expression.exp)
        if expression.is_Pow and expression.exp.is_Rational:
            root = sp.I if expression.base == -1 else sp.sqrt(expression.base)
            return self.generator(root) ** int(expression.exp * 2)
        if expression is sp.I:
            return self.generator(sp.I)
        if _is_transcendental_power(expression):
            out = self.constant(1)
            for generator, coefficient in _exponent_terms(expression):
                out = out * self(generator) ** int(coefficient)
            return out
        raise NotImplementedError(f"{expression} is not a generator")

    def factor(self, polynomial):
        """(content, {normalised irreducible factor: multiplicity}), memoised."""
        if polynomial not in self.factored:
            content, factors = polynomial.factor_list()
            normalised = {}
            for factor, multiplicity in factors:
                scale, factor = _normalised(factor)
                content *= scale ** multiplicity
                normalised[factor] = normalised.get(factor, 0) + multiplicity
            self.factored[polynomial] = (content, normalised)
        return self.factored[polynomial]


def _normalised(polynomial):
    """(scale, p) with polynomial = scale * p, p primitive with a positive leading term."""
    content, primitive = polynomial.primitive()
    if primitive.LC < 0:
        content, primitive = -content, -primitive
    return content, primitive


class _Fraction:
    """numerator / prod(factor**multiplicity), the factors irreducible and normalised.

    Keeping the denominator factored is what makes this fast: a common denominator is
    the maximum multiplicity of each factor, and cancelling is trial division by the
    few factors there are, so no polynomial gcd is ever taken. Once every factor is free
    of the algebraic generators and the numerator has been reduced by their relations,
    removing each factor the numerator is divisible by leaves a form that is unique, and
    so an expression that vanishes has the numerator zero.
    """

    def __init__(self, reader, numerator, denominator):
        self.reader = reader
        self.numerator = numerator
        self.denominator = {f: k for f, k in denominator.items() if k}

    @staticmethod
    def sum(terms, reader):
        common = {}
        for term in terms:
            for factor, multiplicity in term.denominator.items():
                common[factor] = max(common.get(factor, 0), multiplicity)
        numerator = reader.ring.zero
        for term in terms:
            if not term.numerator:
                continue
            cofactor = term.numerator
            for factor, multiplicity in common.items():
                missing = multiplicity - term.denominator.get(factor, 0)
                if missing:
                    cofactor = cofactor * factor ** missing
            numerator += cofactor
        return _Fraction(reader, numerator, common)

    def __add__(self, other):
        return _Fraction.sum([self, other], self.reader)

    def __sub__(self, other):
        return self + other * self.reader.constant(-1)

    def __mul__(self, other):
        denominator = dict(self.denominator)
        for factor, multiplicity in other.denominator.items():
            denominator[factor] = denominator.get(factor, 0) + multiplicity
        return _Fraction(self.reader, self.numerator * other.numerator, denominator)

    def __pow__(self, exponent):
        if exponent < 0:
            return self.inverse() ** -exponent
        return _Fraction(self.reader, self.numerator ** exponent,
                         {f: k * exponent for f, k in self.denominator.items()})

    def inverse(self):
        if not self.numerator:
            raise ZeroDivisionError("the inverse of zero")
        content, factors = self.reader.factor(self.numerator)
        numerator = self.reader.ring(1 / content)
        for factor, multiplicity in self.denominator.items():
            numerator = numerator * factor ** multiplicity
        return _Fraction(self.reader, numerator, factors)

    def holds(self, i):
        """Whether generator i is still anywhere it should not be."""
        return (self.numerator.degree(i) > 1
                or any(factor.degree(i) > 0 for factor in self.denominator))

    def reduce(self, i, base):
        """Fold generator i, w with w^2 = base, out of every denominator factor and down
        to at most its first power in the numerator."""
        free = {f: k for f, k in self.denominator.items() if f.degree(i) <= 0}
        out = _Fraction(self.reader, self.reader.ring.one, free)
        for factor, multiplicity in self.denominator.items():
            if factor.degree(i) <= 0:
                continue
            folded = _fold(self.reader, factor, i, base)
            even, odd = _split(folded.numerator, i)
            unfolded = _Fraction(self.reader, self.reader.ring.one, folded.denominator)
            if odd:
                # 1/(A + Bw) = (A - Bw)/(A^2 - B^2 base), and the right side is free of w.
                square = (_Fraction(self.reader, even ** 2, {})
                          - _Fraction(self.reader, odd ** 2, {}) * base)
                conjugate = _Fraction(self.reader, even - odd * self.reader.ring.gens[i], {})
                reciprocal = unfolded.inverse() * conjugate * square.inverse()
            else:
                reciprocal = unfolded.inverse() * _Fraction(self.reader, even, {}).inverse()
            out = out * reciprocal ** multiplicity
        numerator = _fold(self.reader, self.numerator, i, base)
        return out * numerator

    def cancelled(self):
        numerator = self.numerator
        denominator = {}
        for factor, multiplicity in self.denominator.items():
            while multiplicity and numerator:
                (quotient,), remainder = numerator.div([factor])
                if remainder:
                    break
                numerator, multiplicity = quotient, multiplicity - 1
            denominator[factor] = multiplicity
        return _Fraction(self.reader, numerator, denominator)

    def as_expr(self):
        denominator = sp.Mul(*(factor.as_expr() ** multiplicity for factor, multiplicity
                               in sorted(self.denominator.items(), key=lambda item: str(item[0]))))
        return self.numerator.as_expr() / denominator


def _split(polynomial, i):
    """polynomial = A + B w, for w the generator i, which appears at most linearly."""
    parts = [{}, {}]
    for monomial, coefficient in polynomial.terms():
        k = monomial[i]
        parts[k][monomial[:i] + (0,) + monomial[i + 1:]] = coefficient
    ring = polynomial.ring
    return ring.from_dict(parts[0]), ring.from_dict(parts[1])


def _fold(reader, polynomial, i, base):
    """polynomial with every w^k, w the generator i, replaced by w^(k mod 2) base^(k div 2)."""
    by_degree = {}
    for monomial, coefficient in polynomial.terms():
        free = monomial[:i] + (0,) + monomial[i + 1:]
        by_degree.setdefault(monomial[i], {})[free] = coefficient
    terms = []
    w = reader.ring.gens[i]
    for k, part in by_degree.items():
        term = _Fraction(reader, polynomial.ring.from_dict(part) * w ** (k % 2), {})
        terms.append(term * base ** (k // 2) if k > 1 else term)
    return _Fraction.sum(terms, reader)


class LatexError(Exception):
    pass


class DimensionError(Exception):
    pass


class Timeout(Exception):
    pass


class budget:
    """Abandon a computation that runs past `seconds` rather than let it hang."""

    def __init__(self, seconds):
        self.seconds = seconds

    def __enter__(self):
        signal.signal(signal.SIGALRM, self._fire)
        signal.setitimer(signal.ITIMER_REAL, self.seconds)

    def __exit__(self, *_):
        signal.setitimer(signal.ITIMER_REAL, 0)
        return False

    def _fire(self, *_):
        raise Timeout(f"sympy did not finish inside {self.seconds}s")


def matching_brace(text, opening):
    depth = 0
    for position in range(opening, len(text)):
        if text[position] == "{":
            depth += 1
        elif text[position] == "}":
            depth -= 1
            if depth == 0:
                return position
    raise LatexError(f"unbalanced braces in {text!r}")


def expand_fractions(text):
    while True:
        found = re.search(r"\\[dt]?frac\s*\{", text)
        if not found:
            return text
        first = text.index("{", found.start())
        first_end = matching_brace(text, first)
        second = text.index("{", first_end)
        second_end = matching_brace(text, second)
        numerator = text[first + 1:first_end]
        denominator = text[second + 1:second_end]
        text = f"{text[:found.start()]}(({numerator})/({denominator})){text[second_end + 1:]}"


def expand_braced_call(text, command, replacement):
    """Turn \\command{body} into replacement(body), innermost first."""
    while True:
        found = re.search(r"\\" + command + r"\s*\{", text)
        if not found:
            return text
        opening = text.index("{", found.start())
        closing = matching_brace(text, opening)
        body = text[opening + 1:closing]
        text = f"{text[:found.start()]}{replacement}({body}){text[closing + 1:]}"


PARTIAL = re.compile(
    r"\\partial_\s*\{?\s*(?:\\([A-Za-z]+)|([A-Za-z]))\s*\}?\s*(?:\^\s*\{?\s*(\d+)\s*\}?)?\s*")
PARTIAL_TARGET = re.compile(r"\\?([A-Za-z]+)")
# The single token expand_partials writes, split into the function and the coordinates.
PARTIAL_NAME = re.compile(r"\b([A-Za-z]+)_partial_([A-Za-z]+(?:_[A-Za-z]+)*)\b")


def expand_partials(text):
    """A run of \\partial factors into the one name the reader declares for it.

    \\partial_x^2 H becomes H_partial_x_x and \\partial_x\\partial_y H becomes
    H_partial_x_y, so a published partial derivative is a single token by the time the
    parser sees it. The run binds to the one function name that follows it.
    """
    out = []
    position = 0
    while True:
        found = PARTIAL.search(text, position)
        if not found:
            out.append(text[position:])
            return "".join(out)
        out.append(text[position:found.start()])
        variables = []
        at = found.start()
        while True:
            factor = PARTIAL.match(text, at)
            if factor is None:
                break
            variables += [factor.group(1) or factor.group(2)] * int(factor.group(3) or 1)
            at = factor.end()
        target = PARTIAL_TARGET.match(text, at)
        if target is None:
            raise LatexError(f"the \\partial in {text!r} names no function")
        out.append(f" {target.group(1)}_partial_{'_'.join(variables)} ")
        position = target.end()


def expand_superscript_braces(text):
    """^{body} into **(body), so that ^{-1} and e^{-r^2} survive the parser."""
    out = []
    position = 0
    while position < len(text):
        if text[position] == "^" and position + 1 < len(text) and text[position + 1] == "{":
            closing = matching_brace(text, position + 1)
            out.append("**(" + expand_superscript_braces(text[position + 2:closing]) + ")")
            position = closing + 1
        else:
            out.append(text[position])
            position += 1
    return "".join(out)


class Reader:
    """The LaTeX dialect the metric files are written in, translated into sympy.

    Every name the reader is willing to produce has to be declared up front, so a
    typo in a published value becomes an error here rather than a silent new symbol.
    """

    def __init__(self, coords, parameters, time_coords=frozenset(), relations=None):
        self.coords = list(coords)
        self.time_coords = set(time_coords)
        self.symbol = {}
        self.differential = {}
        self.dot = {}
        self.ddot = {}
        for name in self.coords:
            plain = self._plain(name)
            self.symbol[name] = sp.Symbol(plain, real=True)
            self.differential[name] = sp.Symbol("d" + plain)
            self.dot[name] = sp.Symbol(plain + "_dot")
            self.ddot[name] = sp.Symbol(plain + "_ddot")
        self.c = sp.Symbol("c", positive=True)
        self.parameters = {}
        self.parameter_names = []
        self.functions = {}
        self.primed = set()
        for declaration in parameters:
            self._declare_parameter(declaration)
        self.allowed = (
            set(self.symbol.values())
            | set(self.differential.values())
            | set(self.dot.values())
            | set(self.ddot.values())
            | {self.c}
        )
        for value in self.parameters.values():
            self.allowed |= value.free_symbols
        self.local = {str(s): s for s in self.allowed}
        self.local.update(self.parameters)
        self.local.update(FUNCTIONS)
        self.known = set(self.local)
        self.transforms = standard_transformations + (
            split_symbols_custom(lambda name, _=None: name not in self.known),
            implicit_multiplication,
        )
        self.relations = {}
        for name, value in (relations or {}).items():
            if name not in self.parameters:
                raise LatexError(f"a relation is declared for {name!r}, which is not a parameter")
            self.relations[self.parameters[name]] = sp.sympify(value)

    def surface(self, expression):
        """The expression on the surface the entry's constrained parameters live on."""
        return expression.subs(self.relations)

    @staticmethod
    def _plain(name):
        return name.replace("\\", "").strip()

    def _declare_parameter(self, declaration):
        """`a` is a constant; `a = a(t)` is a function of the coordinates it names.

        For the second kind the printed rate is a derivative with respect to the chart
        coordinate, which is c times the named one when that one is a time, so the rate
        carries the matching power of c. A function of one coordinate answers to a dot
        and to a prime, up to the second derivative, and a function of any number of
        coordinates answers to \\partial at any order, which _declare_partials resolves
        when a published value names one.
        """
        plain = self._plain(declaration.split("=")[0])
        self.parameter_names.append(plain)
        argument = re.search(r"\(([^()]*)\)", declaration)
        if argument is None:
            self.parameters[plain] = sp.Symbol(plain, real=True)
            return
        names = [name.strip() for name in argument.group(1).split(",")]
        unknown = [name for name in names if name not in self.symbol]
        if unknown:
            raise LatexError(f"parameter {declaration!r} varies with {unknown}, "
                             f"which are not coordinates of this system")
        function = sp.Function(plain, real=True)(*(self.symbol[name] for name in names))
        self.parameters[plain] = function
        self.functions[plain] = (function, names)
        if len(names) > 1:
            return
        variable = self.symbol[names[0]]
        scale = self._scale(names[0])
        first = sp.Derivative(function, variable) / scale
        second = sp.Derivative(function, variable, 2) / scale ** 2
        # Both spellings appear: a dot in the comoving chart, a prime in the conformal one.
        for suffix, value in (("_dot", first), ("_prime", first),
                              ("_ddot", second), ("_pprime", second)):
            self.parameters[plain + suffix] = value
        self.primed.add(plain)

    def _scale(self, name):
        """What a derivative along this coordinate is divided by to be one along the chart."""
        return self.c if name in self.time_coords else sp.Integer(1)

    def _declare_partials(self, text):
        """Declare every partial derivative the preprocessed text names, at any order.

        No order is fixed in advance, because the order a curvature reaches depends on
        the entry: a curvature is two derivatives of the metric, so it reaches the
        third derivative of a function whose first derivative the line element already
        carries, as Tolman-Bondi's g_rr carries \\partial_r R. A partial derivative is
        still only declared for a function the entry declares and along coordinates it
        declares that function of, so a typo is an error here rather than a new symbol.
        """
        for found in PARTIAL_NAME.finditer(text):
            name = found.group(0)
            if name in self.local:
                continue
            plain, variables = found.group(1), found.group(2).split("_")
            if plain not in self.functions:
                raise LatexError(f"{text!r} differentiates {plain}, "
                                 f"which is not declared as a function of the coordinates")
            function, arguments = self.functions[plain]
            by_plain = {self._plain(argument): argument for argument in arguments}
            stray = [variable for variable in variables if variable not in by_plain]
            if stray:
                raise LatexError(f"{text!r} differentiates {plain} along {stray}, "
                                 f"and it is declared a function of {arguments} only")
            # One object for every spelling, since mixed partials commute.
            run = sorted((by_plain[variable] for variable in variables), key=self.coords.index)
            value = sp.Derivative(function, *(self.symbol[argument] for argument in run))
            for argument in run:
                value /= self._scale(argument)
            self.parameters[name] = value
            self.local[name] = value
            self.known.add(name)

    def _preprocess(self, latex):
        text = latex
        text = text.replace("\\left", " ").replace("\\right", " ")
        text = re.sub(r"\\[,;:!>]", " ", text)
        text = re.sub(r"\\ ", " ", text)
        text = text.replace("\\cdot", "*")
        for name in sorted(self.primed, key=len, reverse=True):
            text = re.sub(re.escape(name) + r"''", f" {name}_pprime ", text)
            text = re.sub(re.escape(name) + r"'", f" {name}_prime ", text)
        text = expand_partials(text)
        text = expand_fractions(text)
        text = expand_braced_call(text, "sqrt", "sqrt")
        text = expand_braced_call(text, "ddot", "DDOT")
        text = expand_braced_call(text, "dot", "DOT")
        text = re.sub(r"DDOT\s*\(\s*\\?([A-Za-z]+)\s*\)", r" \1_ddot ", text)
        text = re.sub(r"DOT\s*\(\s*\\?([A-Za-z]+)\s*\)", r" \1_dot ", text)
        # d\Omega^2 is the unit two sphere, written out so the reader sees differentials.
        text = re.sub(r"d\\Omega\s*\^\s*2", "(dtheta**2 + sin(theta)**2*dphi**2)", text)
        # Differentials become single atoms before anything is allowed to pad with spaces.
        for name in sorted(self.coords, key=len, reverse=True):
            text = text.replace("d" + name, "d" + self._plain(name))
        for command in sorted(GREEK + TRIG, key=len, reverse=True):
            text = text.replace("\\" + command, " " + command + " ")
        text = text.replace("\\exp", " exp ").replace("\\ln", " log ").replace("\\log", " log ")
        text = expand_superscript_braces(text)
        text = text.replace("^", "**")
        # A trig call written bare, as \sin^2\theta or \cot\theta rather than sin(theta).
        names = "|".join(sorted(TRIG, key=len, reverse=True))
        text = re.sub(r"\b(" + names + r")\s*\*\*\s*\(?\s*(-?\d+)\s*\)?\s*([A-Za-z]\w*)",
                      r"(\1(\3))**\2", text)
        text = re.sub(r"\b(" + names + r")\s+([A-Za-z]\w*)", r"\1(\2)", text)
        text = re.sub(r"(?<![A-Za-z_])e\s*\*\*", " E**", text)
        # 3R\sqrt{..} leaves 3Rsqrt(..), whose one name token would be split letter by letter.
        for name in FUNCTIONS:
            text = re.sub(r"(?<=[A-Za-z0-9_])(" + name + r")\s*\(", r" \1(", text)
        if "\\" in text:
            raise LatexError(f"unhandled LaTeX in {latex!r}: {text!r}")
        return text

    def __call__(self, latex):
        text = self._preprocess(latex)
        self._declare_partials(text)
        try:
            expression = parse_expr(
                text, local_dict=dict(self.local), transformations=self.transforms, evaluate=True
            )
        except Exception as error:
            raise LatexError(f"cannot parse {latex!r} (as {text!r}): {error}") from error
        unknown = expression.free_symbols - self.allowed
        if unknown:
            raise LatexError(f"{latex!r} produced undeclared symbols {sorted(map(str, unknown))}")
        return expression


class Dimensions:
    """What every symbol of a system carries, and the dimension of an expression in them.

    The declared dimension of a coordinate is the one its own letter carries; the chart
    coordinate is c times it when it is a time, so a chart coordinate is never a time.
    Differentials in a line element are written on the bare letter and carry the bare
    dimension, while the dots in a geodesic equation are chart velocities and carry the
    chart dimension over the affine parameter.
    """

    def __init__(self, reader, declared):
        declared = {Reader._plain(name): value for name, value in declared.items()}
        wanted = [Reader._plain(name) for name in reader.coords] + reader.parameter_names
        missing = [name for name in wanted if name not in declared]
        if missing:
            raise DimensionError(f"no declared dimension for {missing}")
        extra = [name for name in declared if name not in wanted]
        if extra:
            raise DimensionError(
                f"a dimension is declared for {extra}, which the system does not use")
        self.of_symbol = {reader.c: LENGTH / TIME}
        self.of_function = {}
        self.chart = {}
        for name in reader.coords:
            bare = sp.sympify(declared[Reader._plain(name)], locals=BASE_DIMENSIONS)
            self.chart[name] = LENGTH if bare == TIME else bare
            self.of_symbol[reader.symbol[name]] = bare
            self.of_symbol[reader.differential[name]] = bare
            self.of_symbol[reader.dot[name]] = self.chart[name] / AFFINE
            self.of_symbol[reader.ddot[name]] = self.chart[name] / AFFINE ** 2
        for name in reader.parameter_names:
            value = reader.parameters[name]
            dimension = sp.sympify(declared[name], locals=BASE_DIMENSIONS)
            if value.is_Symbol:
                self.of_symbol[value] = dimension
            else:
                self.of_function[name] = dimension

    def __call__(self, expression):
        if expression.is_Number or isinstance(expression, sp.NumberSymbol):
            return sp.Integer(1)
        if expression.is_Symbol:
            if expression in self.of_symbol:
                return self.of_symbol[expression]
            raise DimensionError(f"{expression} has no declared dimension")
        if isinstance(expression, sp.core.function.AppliedUndef):
            name = expression.func.__name__
            if name in self.of_function:
                return self.of_function[name]
            raise DimensionError(f"{name} has no declared dimension")
        if expression.is_Add:
            first = self(expression.args[0])
            for term in expression.args[1:]:
                other = self(term)
                if other != first:
                    raise DimensionError(
                        f"{expression} adds {term}, which carries {other}, "
                        f"to {expression.args[0]}, which carries {first}")
            return first
        if expression.is_Mul:
            out = sp.Integer(1)
            for factor in expression.args:
                out *= self(factor)
            return out
        if expression.is_Pow:
            base, exponent = expression.args
            if self(exponent) != 1:
                raise DimensionError(f"the exponent of {expression} carries {self(exponent)}")
            # A symbolic exponent is read with its base in a fixed unit, as Kasner says
            # of its t^{2p_i}, so only the numeric part of the exponent counts.
            numeric, _ = exponent.as_coeff_Add()
            return self(base) ** numeric
        if isinstance(expression, sp.Derivative):
            out = self(expression.expr)
            for variable, order in expression.variable_count:
                out /= self(variable) ** order
            return out
        if isinstance(expression, sp.Function):
            for argument in expression.args:
                if self(argument) != 1:
                    raise DimensionError(
                        f"{expression} takes {argument}, which carries {self(argument)}")
            return sp.Integer(1)
        raise DimensionError(f"cannot take the dimension of {expression}")

    def index_weight(self, variance, coords, index):
        """An upper index carries [x^mu]/L and a lower one L/[x^mu]."""
        out = sp.Integer(1)
        for slot, position in zip(variance, index):
            scale = self.chart[coords[position]] / LENGTH
            out *= scale if slot == "u" else 1 / scale
        return out


def check_dimensions(report, reader, dimensions, where, entry, coords):
    """Every published expression against the dimension its left hand side fixes."""

    def term_by_term(label, expression, expected):
        for term in sp.Add.make_args(expression):
            if term == 0:
                continue
            try:
                carried = dimensions(term)
            except DimensionError as error:
                report.dimension(label, str(error))
                continue
            if carried != expected:
                report.dimension(label, f"the term {term} carries {carried}, not {expected}")

    def read(label, latex):
        try:
            return reader(latex)
        except LatexError as error:
            report.skip(label, str(error))
            return None

    _, _, right = entry["line_element"].partition("=")
    form = read(f"{where}.line_element", right)
    if form is not None:
        term_by_term(f"{where}.line_element", sp.expand(form), LENGTH ** 2)

    for field, variance, published in published_blocks(entry):
        base = FIELD_DIMENSIONS[field]
        for component in published:
            names = component["indices"]
            if any(name not in coords for name in names) or len(names) != len(variance):
                continue
            label = f"{where}.{field}.{variance} {names}"
            value = read(label, component["value"])
            if value is None:
                continue
            index = [coords.index(name) for name in names]
            term_by_term(label, sp.expand(value),
                         base * dimensions.index_weight(variance, coords, index))

    for field in ("ricci_scalar", "kretschmann"):
        if field not in entry:
            continue
        label = f"{where}.{field}"
        text = entry[field].split("=", 1)[1] if "=" in entry[field] else entry[field]
        value = read(label, text)
        if value is not None:
            term_by_term(label, sp.expand(value), FIELD_DIMENSIONS[field])

    for equation in entry.get("geodesics", []):
        label = f"{where}.geodesics"
        left, _, right = equation.partition("=")
        residual = read(f"{label} {equation!r}", left + "-(" + (right or "0") + ")")
        if residual is None:
            continue
        carried = [name for name in coords if residual.has(reader.ddot[name])]
        if len(carried) != 1:
            continue
        expected = dimensions.chart[carried[0]] / AFFINE ** 2
        term_by_term(f"{label} {equation!r}", sp.expand(residual), expected)


def published_blocks(entry):
    """Every published tensor block of an entry, as (field, variance, components)."""
    for field in ("metric_components", "inverse_metric_components"):
        if field in entry:
            yield field, "uu" if field.startswith("inverse") else "ll", entry[field]
    for field in ("christoffel", "riemann", "ricci_tensor", "einstein_tensor", "weyl_tensor"):
        for variance, block in entry.get(field, {}).get("variants", {}).items():
            yield field, variance, block["nonzero"]


def metric_from_line_element(reader, line_element, coords):
    """Read g_{mu nu} off the line element, in the chart its own coords name."""
    _, _, right = line_element.partition("=")
    form = sp.expand(reader(right))
    differentials = [reader.differential[name] for name in coords]
    size = len(coords)
    g = sp.zeros(size, size)
    for i in range(size):
        for j in range(size):
            if i == j:
                g[i, j] = form.coeff(differentials[i], 2)
            else:
                g[i, j] = form.coeff(differentials[i], 1).coeff(differentials[j], 1) / 2
    remainder = sp.expand(
        form - sum(g[i, j] * differentials[i] * differentials[j] for i in range(size) for j in range(size))
    )
    if norm(remainder) != 0:
        raise LatexError(f"{line_element!r} is not a quadratic form in its own differentials")
    return norm(g)


class Geometry:
    """Every tensor the files publish, computed from g in the chart the coords name."""

    def __init__(self, g, coords, seconds):
        self.n = len(coords)
        self.coords = coords
        self.seconds = seconds
        self.g = norm(g)
        # The adjugate over the determinant, which sympy's own inv() takes far longer to reach.
        self.ginv = norm(self.g.adjugate() / self.g.det())
        self._cache = {}
        self.unavailable = {}

    def _timed(self, name, build):
        if name in self._cache:
            return self._cache[name]
        if name in self.unavailable:
            raise Timeout(self.unavailable[name])
        try:
            with budget(self.seconds):
                self._cache[name] = build()
        except Timeout as error:
            self.unavailable[name] = str(error)
            raise
        return self._cache[name]

    def _zeros(self, rank):
        if rank == 1:
            return [sp.Integer(0)] * self.n
        return [self._zeros(rank - 1) for _ in range(self.n)]

    def christoffel_ull(self):
        def build():
            out = self._zeros(3)
            for mu in range(self.n):
                for nu in range(self.n):
                    for rho in range(nu, self.n):
                        value = norm(sum(
                            self.ginv[mu, alpha] * (
                                sp.diff(self.g[alpha, rho], self.coords[nu])
                                + sp.diff(self.g[alpha, nu], self.coords[rho])
                                - sp.diff(self.g[nu, rho], self.coords[alpha])
                            )
                            for alpha in range(self.n)
                        ) / 2)
                        out[mu][nu][rho] = value
                        out[mu][rho][nu] = value
            return out
        return self._timed("christoffel_ull", build)

    def christoffel_lll(self):
        def build():
            gamma = self.christoffel_ull()
            out = self._zeros(3)
            for mu in range(self.n):
                for nu in range(self.n):
                    for rho in range(self.n):
                        out[mu][nu][rho] = norm(
                            sum(self.g[mu, alpha] * gamma[alpha][nu][rho] for alpha in range(self.n))
                        )
            return out
        return self._timed("christoffel_lll", build)

    def riemann_ulll(self):
        def build():
            gamma = self.christoffel_ull()
            out = self._zeros(4)
            for mu in range(self.n):
                for nu in range(self.n):
                    for rho in range(self.n):
                        for sigma in range(rho + 1, self.n):
                            value = norm(
                                sp.diff(gamma[mu][nu][sigma], self.coords[rho])
                                - sp.diff(gamma[mu][nu][rho], self.coords[sigma])
                                + sum(
                                    gamma[mu][rho][lam] * gamma[lam][nu][sigma]
                                    - gamma[mu][sigma][lam] * gamma[lam][nu][rho]
                                    for lam in range(self.n)
                                )
                            )
                            out[mu][nu][rho][sigma] = value
                            out[mu][nu][sigma][rho] = -value
            return out
        return self._timed("riemann_ulll", build)

    def riemann_llll(self):
        def build():
            upper = self.riemann_ulll()
            out = self._zeros(4)
            for mu in range(self.n):
                for nu in range(self.n):
                    for rho in range(self.n):
                        for sigma in range(self.n):
                            out[mu][nu][rho][sigma] = norm(sum(
                                self.g[mu, alpha] * upper[alpha][nu][rho][sigma]
                                for alpha in range(self.n)
                            ))
            return out
        return self._timed("riemann_llll", build)

    def ricci_ll(self):
        """R_{mu nu} = R^a_{mu a nu}, the standard contraction on the first lower index."""
        def build():
            riemann = self.riemann_ulll()
            out = self._zeros(2)
            for nu in range(self.n):
                for sigma in range(nu, self.n):
                    value = norm(sum(riemann[mu][nu][mu][sigma] for mu in range(self.n)))
                    out[nu][sigma] = value
                    out[sigma][nu] = value
            return out
        return self._timed("ricci_ll", build)

    def ricci_scalar(self):
        def build():
            ricci = self.ricci_ll()
            return norm(sum(
                self.ginv[a, b] * ricci[a][b] for a in range(self.n) for b in range(self.n)
            ))
        return self._timed("ricci_scalar", build)

    def einstein_ll(self):
        def build():
            ricci = self.ricci_ll()
            scalar = self.ricci_scalar()
            out = self._zeros(2)
            for a in range(self.n):
                for b in range(self.n):
                    out[a][b] = norm(ricci[a][b] - scalar * self.g[a, b] / 2)
            return out
        return self._timed("einstein_ll", build)

    def kretschmann(self):
        def build():
            lower = self.riemann_llll()
            # One index at a time, so each component is a sum of n terms rather than n^4.
            upper = self.raise_indices(lower, 4, (0, 1, 2, 3))
            return norm(sum(
                lower[a][b][cc][d] * upper[a][b][cc][d]
                for a in range(self.n) for b in range(self.n)
                for cc in range(self.n) for d in range(self.n)
            ))
        return self._timed("kretschmann", build)

    def weyl_llll(self):
        def build():
            n = self.n
            riemann = self.riemann_llll()
            ricci = self.ricci_ll()
            scalar = self.ricci_scalar()
            out = self._zeros(4)
            for a in range(n):
                for b in range(n):
                    for cc in range(n):
                        for d in range(n):
                            correction = (
                                self.g[a, cc] * ricci[d][b] - self.g[a, d] * ricci[cc][b]
                                - self.g[b, cc] * ricci[d][a] + self.g[b, d] * ricci[cc][a]
                            ) / (n - 2)
                            trace = scalar * (
                                self.g[a, cc] * self.g[d, b] - self.g[a, d] * self.g[cc, b]
                            ) / ((n - 1) * (n - 2))
                            out[a][b][cc][d] = norm(riemann[a][b][cc][d] - correction + trace)
            return out
        return self._timed("weyl_llll", build)

    def raise_indices(self, tensor, rank, positions):
        """Raise the given index positions of a fully lowered tensor."""
        out = tensor
        for position in positions:
            source = out
            out = self._zeros(rank)
            for index in _indices(self.n, rank):
                total = 0
                for alpha in range(self.n):
                    swapped = list(index)
                    swapped[position] = alpha
                    total += self.ginv[index[position], alpha] * _at(source, swapped)
                _put(out, index, norm(total))
        return out


def _indices(n, rank):
    if rank == 0:
        yield []
        return
    for head in range(n):
        for tail in _indices(n, rank - 1):
            yield [head] + tail


def _at(tensor, index):
    for i in index:
        tensor = tensor[i]
    return tensor


def _put(tensor, index, value):
    for i in index[:-1]:
        tensor = tensor[i]
    tensor[index[-1]] = value


class Report:
    def __init__(self):
        self.disagreements = []
        self.dimensional = []
        self.unchecked = []
        self.checked_systems = 0
        self.systems = 0

    def disagree(self, where, message):
        self.disagreements.append(f"{where}: {message}")

    def dimension(self, where, message):
        self.dimensional.append(f"{where}: {message}")

    def skip(self, where, reason):
        self.unchecked.append(f"{where}: {reason}")

    def guarded(self, where, seconds, work):
        """Run a comparison under the clock, so a slow one is named rather than hung on."""
        try:
            with budget(seconds):
                work()
        except Timeout as error:
            self.skip(where, f"comparing the published values: {error}")


def variance_weight(variance, coords, index, time_coords):
    """c^(upper time indices - lower time indices) for the x^0 = cT chart."""
    exponent = 0
    for slot, position in zip(variance, index):
        if coords[position] in time_coords:
            exponent += 1 if slot == "u" else -1
    return exponent


def compare_block(report, reader, where, published, computed, variance, coords, time_coords, c):
    """Every published component against sympy, and every omitted one against zero."""
    rank = len(variance)
    seen = set()
    for entry in published:
        names = entry["indices"]
        if len(names) != rank:
            report.disagree(where, f"{names} has {len(names)} indices, expected {rank}")
            continue
        try:
            index = [coords.index(name) for name in names]
        except ValueError:
            report.disagree(where, f"{names} names an index that is not one of {coords}")
            continue
        seen.add(tuple(index))
        try:
            value = reader.surface(reader(entry["value"]))
        except LatexError as error:
            report.skip(f"{where} {names}", str(error))
            continue
        expected = reader.surface(
            _at(computed, index) * c ** variance_weight(variance, coords, index, time_coords))
        if norm(value - expected) != 0:
            report.disagree(where, f"{names} published as {entry['value']} "
                                   f"({norm(value)}), sympy says {norm(expected)}")
    for index in _indices(len(coords), rank):
        if tuple(index) in seen:
            continue
        names = [coords[i] for i in index]
        weighted = reader.surface(
            _at(computed, index) * c ** variance_weight(variance, coords, index, time_coords))
        if norm(weighted) != 0:
            report.disagree(where, f"{names} is missing, sympy says {norm(weighted)}")


def compare_scalar(report, reader, where, published, computed):
    text = published.split("=", 1)[1] if "=" in published else published
    try:
        value = reader.surface(reader(text))
    except LatexError as error:
        report.skip(where, str(error))
        return
    computed = reader.surface(computed)
    if norm(value - computed) != 0:
        report.disagree(where, f"published as {published.strip()}, sympy says {norm(computed)}")


def compare_geodesics(report, reader, where, published, gamma, coords, time_coords, c):
    """Each equation against xddot^mu + Gamma^mu_{nu rho} xdot^nu xdot^rho = 0.

    The printed dots are the chart velocities, matching the printed Christoffels, so
    the weighted symbols are the ones the residual is built from.
    """
    if len(published) != len(coords):
        report.disagree(where, f"{len(published)} equations for {len(coords)} coordinates")
    for equation in published:
        left, _, right = equation.partition("=")
        try:
            residual = reader.surface(sp.expand(reader(left) - reader(right)))
        except LatexError as error:
            report.skip(f"{where} {equation!r}", str(error))
            continue
        carried = [name for name in coords if residual.has(reader.ddot[name])]
        if len(carried) != 1:
            report.disagree(where, f"{equation!r} carries second derivatives of {carried}, expected one")
            continue
        name = carried[0]
        mu = coords.index(name)
        expected = reader.surface(reader.ddot[name] + sum(
            gamma[mu][nu][rho]
            * c ** variance_weight("ull", coords, [mu, nu, rho], time_coords)
            * reader.dot[coords[nu]] * reader.dot[coords[rho]]
            for nu in range(len(coords)) for rho in range(len(coords))
        ))
        if norm(residual - expected) != 0 and norm(residual + expected) != 0:
            report.disagree(where, f"{name} equation {equation!r} is not the geodesic equation, "
                                   f"sympy makes the residual {norm(expected)}")


VECTOR_VARIANTS = {"ll": (0, "ll"), "ul": ((0,), "ul"), "uu": ((0, 1), "uu")}
RANK4_VARIANTS = {"llll": (None, "llll"), "ulll": ((0,), "ulll")}


def check_system(report, metric_id, entry, seconds, dimensions_only=False):
    where = f"{metric_id}/{entry['id']}"
    coords = entry["coords"]
    report.systems += 1

    declared = DIMENSIONS.get((metric_id, entry["id"]))
    if declared is None:
        report.skip(where, "no dimension declaration in DIMENSIONS")
        return
    unknown = [name for name in coords if name not in declared]
    if unknown:
        report.skip(where, f"no declared dimension for the coordinates {unknown}")
        return
    declaration = time_coordinates(declared, coords)

    parameters = [p["symbol"] for p in entry.get("parameters", [])]
    relations = PARAMETER_RELATIONS.get((metric_id, entry["id"]), {})
    try:
        reader = Reader(coords, parameters, declaration, relations)
    except LatexError as error:
        report.skip(where, f"parameters unreadable: {error}")
        return
    c = reader.c
    try:
        dimensions = Dimensions(reader, declared)
    except DimensionError as error:
        report.skip(where, f"dimensions undeclared: {error}")
        return
    check_dimensions(report, reader, dimensions, where, entry, coords)
    if dimensions_only:
        report.checked_systems += 1
        print(f"  {where}")
        return
    try:
        g = metric_from_line_element(reader, entry["line_element"], coords)
    except LatexError as error:
        report.skip(where, f"line element unreadable: {error}")
        return
    if norm(g.det()) == 0:
        report.skip(where, "the line element gives a degenerate metric")
        return

    # Everything is computed with the bare coordinates and then weighted into the
    # x^0 = cT chart by compare_block, which is exact because the rescaling is linear.
    symbols = [reader.symbol[name] for name in coords]
    geometry = Geometry(g, symbols, seconds)
    report.checked_systems += 1
    print(f"  {where}")

    published_metric = [(f"{where}.metric_components", entry.get("metric_components"), "ll", g)]
    if "inverse_metric_components" in entry:
        published_metric.append((
            f"{where}.inverse_metric_components", entry["inverse_metric_components"],
            "uu", geometry.ginv,
        ))
    else:
        # Every system in the collection publishes its inverse metric, since the page
        # prints one beside the metric; Ellis-Bronnikov was the last without one and now
        # has one too, so this is not a slot the collection leaves empty on purpose. A
        # system that omits it is missing something the page would show, and saying so
        # here keeps that visible rather than passing it in silence.
        report.skip(f"{where}.inverse_metric_components", "the entry does not publish one")
    for label, published, variance, matrix in published_metric:
        as_lists = [[matrix[i, j] for j in range(len(coords))] for i in range(len(coords))]
        report.guarded(label, seconds, lambda label=label, published=published,
                       as_lists=as_lists, variance=variance: compare_block(
            report, reader, label, published, as_lists, variance, coords, declaration, c))
    if norm(g * geometry.ginv) != sp.eye(len(coords)):
        report.disagree(where, "the published metric and inverse are not inverse to each other")

    blocks = [
        ("christoffel", "ull", lambda: geometry.christoffel_ull()),
        ("christoffel", "lll", lambda: geometry.christoffel_lll()),
        ("riemann", "ulll", lambda: geometry.raise_indices(geometry.riemann_llll(), 4, (0,))),
        ("riemann", "llll", lambda: geometry.riemann_llll()),
        ("ricci_tensor", "ll", lambda: geometry.ricci_ll()),
        ("ricci_tensor", "ul", lambda: geometry.raise_indices(geometry.ricci_ll(), 2, (0,))),
        ("ricci_tensor", "uu", lambda: geometry.raise_indices(geometry.ricci_ll(), 2, (0, 1))),
        ("einstein_tensor", "ll", lambda: geometry.einstein_ll()),
        ("einstein_tensor", "ul", lambda: geometry.raise_indices(geometry.einstein_ll(), 2, (0,))),
        ("einstein_tensor", "uu", lambda: geometry.raise_indices(geometry.einstein_ll(), 2, (0, 1))),
        ("weyl_tensor", "llll", lambda: geometry.weyl_llll()),
        ("weyl_tensor", "ulll", lambda: geometry.raise_indices(geometry.weyl_llll(), 4, (0,))),
    ]
    for field, variance, build in blocks:
        variants = entry.get(field, {}).get("variants", {})
        if variance not in variants:
            continue
        label = f"{where}.{field}.{variance}"
        try:
            computed = build()
        except Timeout as error:
            report.skip(label, str(error))
            continue
        report.guarded(label, seconds, lambda label=label, computed=computed,
                       variance=variance, published=variants[variance]["nonzero"]:
                       compare_block(report, reader, label, published, computed,
                                     variance, coords, declaration, c))

    for field, build in (("ricci_scalar", geometry.ricci_scalar), ("kretschmann", geometry.kretschmann)):
        if field not in entry:
            continue
        label = f"{where}.{field}"
        try:
            computed = build()
        except Timeout as error:
            report.skip(label, str(error))
            continue
        report.guarded(label, seconds, lambda label=label, computed=computed, field=field:
                       compare_scalar(report, reader, label, entry[field], computed))

    if "geodesics" in entry:
        label = f"{where}.geodesics"
        try:
            gamma = geometry.christoffel_ull()
        except Timeout as error:
            report.skip(label, str(error))
        else:
            report.guarded(label, seconds, lambda: compare_geodesics(
                report, reader, label, entry["geodesics"], gamma, coords, declaration, c))
    else:
        report.skip(f"{where}.geodesics", "the entry does not publish any")


CONFLICT_COPY = re.compile(r" \d+$")


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--system", action="append", default=[],
                        help="check only <metric_id>/<system_id>, repeatable")
    parser.add_argument("--budget", type=float, default=DEFAULT_BUDGET_SECONDS,
                        help="seconds sympy may spend on one tensor")
    parser.add_argument("--dimensions-only", action="store_true",
                        help="run the dimensional pass alone and skip the sympy algebra")
    arguments = parser.parse_args()

    report = Report()
    files = 0
    against = "for dimensional balance" if arguments.dimensions_only else "against sympy"
    print(f"Checking every published coordinate system {against}, in the x^0 = cT chart.\n")
    for path in sorted(METRICS_DIR.glob("*.json")):
        if CONFLICT_COPY.search(path.stem):
            continue
        files += 1
        metric = json.loads(path.read_text(encoding="utf-8"))
        for entry in metric.get("coordinates", []):
            name = f"{metric['id']}/{entry['id']}"
            if arguments.system and name not in arguments.system:
                continue
            check_system(report, metric["id"], entry, arguments.budget,
                         arguments.dimensions_only)

    print(f"\n{files} metric files, {report.systems} coordinate systems, "
          f"{report.checked_systems} checked.")
    if report.unchecked:
        print(f"\n{len(report.unchecked)} UNCHECKED:")
        for line in report.unchecked:
            print(f"  UNCHECKED {line}")
    if report.dimensional:
        print(f"\n{len(report.dimensional)} terms whose dimensions do not balance:",
              file=sys.stderr)
        for line in report.dimensional:
            print(f"  {line}", file=sys.stderr)
    if report.disagreements:
        print(f"\n{len(report.disagreements)} disagreements:", file=sys.stderr)
        for line in report.disagreements:
            print(f"  {line}", file=sys.stderr)
    if report.dimensional or report.disagreements:
        return 1
    print(f"\nEvery published value balances dimensionally"
          f"{'' if arguments.dimensions_only else ' and agrees with sympy'}.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
