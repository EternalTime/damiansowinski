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

    _tools/setup-env.sh
    .venv.noindex/bin/python _tools/derivations/verify_metrics.py

Pass --system <metric_id>/<system_id> to check one system, repeatable.
Pass --metric <metric_id> to check every system that metric publishes, repeatable.
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


Systems kept to an order
------------------------

A slowly rotating star's exterior is known as a series in its spin, and Hartle and Thorne's
line element solves the vacuum equations through the second order of it and no further. Such
a system claims its connection and curvature to that order only, and its Ricci tensor, zero
to that order, is not zero beyond it.

ORDERS below names, per system, the order each small parameter counts as and the highest
order kept. Every tensor of such a system is built as a Taylor polynomial to that order,
cut after each product, which gives the same polynomial as cutting the exact tensor and
keeps every step small, and each comparison is made between the two sides' polynomials.
A published metric component or inverse component may stand as the line element writes it,
uncut, as Hartle and Thorne's (1 - 2j_2 P_2)/F does, and the drawings, which read both,
then read a metric and its exact inverse. A published connection or curvature that carries a
term beyond the order is a disagreement, so nothing of higher order is printed as if it were
known.

A post-Newtonian metric is kept to an order that depends on the component: with the velocity
v/c counted as first order, g_tt is known through the fourth, g_ti through the third and g_ij
through the second. Such a system names its time coordinate as a third item of its ORDERS
line, and each tensor is then cut component by component: a connection, Riemann or Ricci
component with n time indices, wherever they stand, is kept n orders further than the order
the line names, an Einstein or Weyl component only where n is odd, since both take a trace
with the metric, the Ricci scalar at the order named and the Kretschmann scalar at twice the
lowest order of the curvature. Those are the terms the metric fixes: one order further, each
would need a term of the metric that the post-Newtonian metric does not have, which
_tools/test_ppn_metric.py checks by adding such terms and finding the kept ones unmoved.

The small quantities of such a system may be declared functions of the coordinates, as the
potentials of Einstein, Infeld and Hoffmann's field of many bodies are: its ORDERS line gives
each function its order, and a derivative of one along the named time counts one order more
each time it is taken, since the bodies that make the potentials move slowly. A line whose third
item is None declares its potentials the same way with no time slow, as the perturbed Friedmann
universe's does, and then every derivative of a potential counts as the potential does; a line of
two items counts constants alone and refuses a function.
"""

import argparse
import functools
import inspect
import json
import re
import signal
import sys
import time
from pathlib import Path

import sympy as sp
from sympy.polys import compatibility, factortools
from sympy.polys.polyerrors import BasePolynomialError
from sympy.polys.rings import PolyRing
from sympy.parsing.sympy_parser import (
    implicit_multiplication,
    parse_expr,
    split_symbols_custom,
    standard_transformations,
)

# The seed every factorisation of a polynomial in several variables draws its evaluation points
# from. sympy factors such a polynomial by Wang's algorithm, which lifts the factors of the
# polynomial at a point chosen at random, and draws that point from a generator seeded from the
# operating system when sympy is imported, so the time one factorisation takes changed from run
# to run, whatever PYTHONHASHSEED was. On 2 October 2026 the same tensor of Brill's charged
# Taub-NUT took 2 seconds in one run and over 800 in another, in the Christoffel symbols, the
# Riemann tensor or the geodesics. Drawn from this seed, each polynomial is lifted at the same
# points in every run, whatever was factored before it, so a chart prints in the same time each
# time; the factors are unique, so no printed value changes.
WANG_SEED = 0


def _seeded(wang):
    """Wang's algorithm with its points drawn from WANG_SEED wherever the caller names no seed.
    sympy calls it again with a larger modulus when a choice of points fails, and that call is
    seeded the same way, through the module's name."""
    if "seed" not in inspect.signature(wang).parameters:
        raise ImportError("sympy's dmp_zz_wang takes no seed, so its factorisations cannot be seeded")

    @functools.wraps(wang)
    def seeded(f, u, K, mod=None, seed=None):
        return wang(f, u, K, mod=mod, seed=WANG_SEED if seed is None else seed)
    seeded.seeded = True
    return seeded


if not getattr(factortools.dmp_zz_wang, "seeded", False):
    factortools.dmp_zz_wang = _seeded(factortools.dmp_zz_wang)
    # The ring's own method calls the name it imported.
    compatibility.dmp_zz_wang = factortools.dmp_zz_wang

ROOT = Path(__file__).resolve().parents[2]
METRICS_DIR = ROOT / "MFS" / "assets" / "data" / "metrics"

DEFAULT_BUDGET_SECONDS = 120
# The systems whose slowest tensor needs more than the default, each with what it takes. The
# toroidal chart of Bach and Weyl's ring writes psi with the radical sqrt(cosh(zeta) - cos(sigma))
# and e^(zeta/2) beside both elliptic integrals, and its Kretschmann scalar took 100 to 130 seconds
# on 2 October 2026, as the machine was loaded.
SYSTEM_BUDGET_SECONDS = {("bach_weyl_ring", "toroidal"): 480}

LENGTH = sp.Symbol("L", positive=True)
TIME = sp.Symbol("T", positive=True)
MASS = sp.Symbol("M", positive=True)
AFFINE = sp.Symbol("lambda", positive=True)
# The name the parser reads a parameter \lambda under, since `lambda` is a Python keyword.
KEYWORD_LAMBDA = "lambda_"
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
    # Three dimensions, where Newton's constant carries L^2/(M T^2): the entry's M = 8Gm/c^2
    # is dimensionless and J = 8Gj/c^3 a length, which is what lets N^2 = r^2/l^2 - M + J^2/(4r^2)
    # be a sum of numbers. The advanced and retarded times of the Eddington-Finkelstein charts
    # are lengths, as Schwarzschild's are.
    ("btz", "stationary"): {
        "t": "T", "r": "L", "\\phi": "1", "\\ell": "L", "M": "1", "J": "L",
    },
    ("btz", "eddington_finkelstein_ingoing"): {
        "v": "L", "r": "L", "\\tilde\\phi": "1", "\\ell": "L", "M": "1", "J": "L",
    },
    ("btz", "eddington_finkelstein_outgoing"): {
        "u": "L", "r": "L", "\\tilde\\phi": "1", "\\ell": "L", "M": "1", "J": "L",
    },
    # The quantum BTZ black hole keeps the BTZ entry's reading of the mass, M = 8 G_3 m/c^2 a number,
    # with the brane's anti-de Sitter radius l_3 and the backreaction length l both lengths and
    # Emparan, Frassino and Way's F, kappa and mu numbers, so that l F/r and mu l/r are numbers
    # beside r^2/l_3^2. The rotation a is a length, which makes a^2/r^2 a number and a/r^2 times
    # c dt an angle. The advanced and retarded times are lengths.
    # The Big Rip: w is a number below -1, t_0 the time left at a = 1 and the scale factor a a
    # number, so the comoving coordinates x, y, z and r are lengths, the proper lengths at
    # t = -t_0. The conformal time eta and its scale eta_0 are lengths, d(eta) = c dt/a.
    ("big_rip", "cartesian"): {
        "t": "T", "x": "L", "y": "L", "z": "L", "w": "1", "t_0": "T", "a": "1",
    },
    ("big_rip", "comoving"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "w": "1", "t_0": "T", "a": "1",
    },
    ("big_rip", "conformal"): {
        "\\eta": "L", "r": "L", "\\theta": "1", "\\phi": "1", "w": "1", "\\eta_0": "L", "a": "1",
    },
    ("quantum_btz", "static"): {
        "t": "T", "r": "L", "\\phi": "1", "\\ell_3": "L", "M": "1", "\\ell": "L", "F": "1",
    },
    ("quantum_btz", "eddington_finkelstein_ingoing"): {
        "v": "L", "r": "L", "\\phi": "1", "\\ell_3": "L", "M": "1", "\\ell": "L", "F": "1",
    },
    ("quantum_btz", "eddington_finkelstein_outgoing"): {
        "u": "L", "r": "L", "\\phi": "1", "\\ell_3": "L", "M": "1", "\\ell": "L", "F": "1",
    },
    ("quantum_btz", "brane"): {
        "t": "T", "r": "L", "\\phi": "1", "\\ell_3": "L", "\\kappa": "1", "\\mu": "1", "\\ell": "L",
    },
    ("quantum_btz", "rotating"): {
        "t": "T", "r": "L", "\\phi": "1", "\\ell_3": "L", "\\kappa": "1", "\\mu": "1", "\\ell": "L", "a": "L",
    },
    # Witten's black hole in two dimensions: lambda is an inverse length and the mass parameter m,
    # the value of e^(-2 Phi) on the horizon, a number. The dilaton chart's w = e^(-2 Phi) and the
    # Kruskal coordinates U and V are numbers, and the advanced and retarded times are lengths.
    ("witten_black_hole", "witten"): {"t": "T", "r": "L", "\\lambda": "1/L"},
    ("witten_black_hole", "schwarzschild_gauge"): {"t": "T", "x": "L", "\\lambda": "1/L", "m": "1"},
    ("witten_black_hole", "dilaton"): {"t": "T", "w": "1", "\\lambda": "1/L", "m": "1"},
    ("witten_black_hole", "conformal"): {"t": "T", "\\sigma": "L", "\\lambda": "1/L", "m": "1"},
    ("witten_black_hole", "kruskal"): {"U": "1", "V": "1", "\\lambda": "1/L", "m": "1"},
    ("witten_black_hole", "eddington_finkelstein_ingoing"): {"v": "L", "x": "L", "\\lambda": "1/L", "m": "1"},
    ("witten_black_hole", "eddington_finkelstein_outgoing"): {"u": "L", "x": "L", "\\lambda": "1/L", "m": "1"},
    # The black hole of Jackiw and Teitelboim's gravity: the curvature radius L and the horizon r_h are
    # lengths in every chart; the Kruskal and global coordinates are numbers, and the Poincare z a length.
    ("jackiw_teitelboim_black_hole", "static"): {"t": "T", "r": "L", "L": "L", "r_h": "L"},
    ("jackiw_teitelboim_black_hole", "proper_distance"): {"t": "T", "\\rho": "L", "L": "L", "r_h": "L"},
    ("jackiw_teitelboim_black_hole", "kruskal"): {"U": "1", "V": "1", "L": "L", "r_h": "L"},
    ("jackiw_teitelboim_black_hole", "global"): {"\\tau": "1", "\\sigma": "1", "L": "L", "r_h": "L"},
    ("jackiw_teitelboim_black_hole", "poincare"): {"T": "T", "z": "L", "L": "L", "r_h": "L"},
    # Maldacena and Qi's eternal traversable wormhole: every chart is L^2 times theirs, so every
    # coordinate is a number and L the one length; the dilaton's least value phi_h is a number.
    ("eternal_traversable_wormhole", "global"): {"t": "1", "\\sigma": "1", "L": "L", "\\phi_h": "1"},
    ("eternal_traversable_wormhole", "poincare"): {"t_P": "1", "z": "1", "L": "L", "\\phi_h": "1"},
    ("eternal_traversable_wormhole", "thermal"): {"t_R": "1", "\\rho": "1", "L": "L", "\\phi_h": "1"},
    # The wave amplitude psi and gamma sit in exponentials and are dimensionless; the null chart's
    # u = ct - rho and v = ct + rho are lengths, as the Eddington-Finkelstein times of BTZ are.
    ("einstein_rosen_waves", "cylindrical"): {
        "t": "T", "\\rho": "L", "\\phi": "1", "z": "L", "\\psi": "1", "\\gamma": "1",
    },
    ("einstein_rosen_waves", "null"): {
        "u": "L", "v": "L", "\\phi": "1", "z": "L", "\\psi": "1", "\\gamma": "1",
    },
    # The second polarisation omega multiplies dphi beside dz, so it is a length; psi and gamma are as
    # Einstein and Rosen's.
    ("twisting_cylindrical_waves", "cylindrical"): {
        "t": "T", "\\rho": "L", "\\phi": "1", "z": "L", "\\psi": "1", "\\omega": "L", "\\gamma": "1",
    },
    ("twisting_cylindrical_waves", "null"): {
        "u": "L", "v": "L", "\\phi": "1", "z": "L", "\\psi": "1", "\\omega": "L", "\\gamma": "1",
    },
    # Gowdy's coordinates are all pure numbers, the areal time among them, and the one length L
    # multiplies the whole line element, so no coordinate is a time the chart multiplies by c.
    ("gowdy", "areal"): {
        "t": "1", "\\theta": "1", "\\sigma": "1", "\\delta": "1", "L": "L", "P": "1", "Q": "1", "\\lambda": "1",
    },
    ("gowdy", "logarithmic"): {
        "\\tau": "1", "\\theta": "1", "\\sigma": "1", "\\delta": "1", "L": "L", "P": "1", "Q": "1", "\\lambda": "1",
    },
    ("gowdy", "sphere"): {
        "t": "1", "\\theta": "1", "\\sigma": "1", "\\delta": "1", "L": "L", "P": "1", "Q": "1", "a": "1",
    },
    # Belinski and Zakharov's solitons: the pole chart's tau and xi are pure numbers and the length
    # w carries the scale, x and y lengths along the Killing vectors; the canonical chart's t is a
    # time and its three functions pure numbers.
    ("belinski_zakharov", "pole"): {
        "\\tau": "1", "\\xi": "1", "x": "L", "y": "L", "w": "L", "\\beta": "1", "N": "1",
    },
    ("belinski_zakharov", "canonical"): {
        "t": "T", "z": "L", "x": "L", "y": "L", "w": "L", "f": "1", "P": "1", "Q": "1",
    },
    # B = sqrt(G) B_0/c^2 folds the field in Gaussian units into an inverse length, so B rho and
    # B r sin(theta) are pure numbers; Ernst's hole keeps Schwarzschild's r_s = 2GM/c^2.
    ("melvin", "cylindrical"): {
        "t": "T", "\\rho": "L", "\\phi": "1", "z": "L", "B": "1/L",
    },
    ("melvin", "ernst"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "r_s": "L", "B": "1/L",
    },
    # Astorino, Martelli and Vigano's Schwarzschild hole in a swirling universe: m = GM/c^2 is a
    # length and the rotation parameter j an inverse area, so j r^2 and j rho^2 are pure numbers and
    # the dragging 4j(r - 2m)cos(theta), 4jz in the universe with no hole, an inverse length.
    ("swirling_universe_black_hole", "spherical"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "m": "L", "j": "1/L**2",
    },
    ("swirling_universe_black_hole", "cylindrical"): {
        "t": "T", "\\rho": "L", "\\phi": "1", "z": "L", "j": "1/L**2",
    },
    # Ernst and Wild's Kerr hole in Melvin's universe: m = GM/c^2 and a are lengths and B an
    # inverse length; k, the value of H on the axis, H and the squared lapse N are pure numbers,
    # Sigma, Delta, F and P areas, A the square of an area, and the dragging rate omega an inverse
    # length, as the rate Omega at which the second chart's azimuth turns is.
    ("kerr_melvin", "boyer_lindquist"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "m": "L", "a": "L", "B": "1/L", "k": "1",
        "\\Sigma": "L**2", "\\Delta": "L**2", "A": "L**4", "H": "1", "\\omega": "1/L",
        "N": "1", "F": "L**2", "P": "L**2",
    },
    ("kerr_melvin", "rotating"): {
        "t": "T", "r": "L", "\\theta": "1", "\\tilde\\phi": "1", "m": "L", "a": "L", "B": "1/L", "k": "1",
        "\\Sigma": "L**2", "\\Delta": "L**2", "A": "L**4", "H": "1", "\\omega": "1/L",
        "N": "1", "F": "L**2", "P": "L**2", "\\Omega": "1/L",
    },
    # Podolsky and Ovcharenko's Kerr hole in Bertotti and Robinson's field: m = GM/c^2 and a are
    # lengths and B an inverse length; I_1, I_2, the conicity C, P, the conformal factor Omega and
    # the squared lapse N are pure numbers, Delta, Sigma, Q, F and W areas, A the square of an area,
    # and the dragging rate omega an inverse length. The static chart's f is a pure number.
    ("kerr_bertotti_robinson", "boyer_lindquist"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "m": "L", "a": "L", "B": "1/L", "I_1": "1", "I_2": "1",
        "C": "1", "\\Delta": "L**2", "\\Sigma": "L**2", "P": "1", "Q": "L**2", "\\Omega": "1", "A": "L**4",
        "\\omega": "1/L", "N": "1", "F": "L**2", "W": "L**2",
    },
    ("kerr_bertotti_robinson", "static"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "m": "L", "B": "1/L", "C": "1", "P": "1", "f": "1",
        "\\Omega": "1",
    },
    # sigma = G lambda/c^2, the mass per unit length lambda as a pure number, and the conicity C
    # are dimensionless. A power of the radius whose exponent holds sigma, or one of Kasner's
    # exponents, is read with the radius in a fixed unit, so rho^{2 - 4 sigma} is an area and
    # ell, the radius of the circle at unit r, is a length.
    ("levi_civita", "weyl"): {
        "t": "T", "\\rho": "L", "\\phi": "1", "z": "L", "\\sigma": "1", "C": "1",
    },
    ("levi_civita", "kasner"): {
        "t": "T", "r": "L", "\\phi": "1", "z": "L", "p_0": "1", "p_2": "1", "p_3": "1", "\\ell": "L",
    },
    # f and gamma are pure numbers and omega is a length: with x^0 = ct the line element's
    # c dt - omega dphi is a length term by term.
    ("double_kerr", "weyl"): {
        "t": "T", "\\rho": "L", "\\phi": "1", "z": "L", "f": "1", "\\omega": "L", "\\gamma": "1",
    },
    # Neugebauer and Meinel's disc: U and k are pure numbers and a is a length, so that c dt + a dphi
    # is a length term by term, in the frame at rest far away and in the one that turns with the
    # disc. In Bardeen and Wagoner's form nu and alpha are pure numbers and omega is an inverse
    # length, the dragging per unit of ct. The spheroidal xi and eta are pure numbers and rho_0
    # carries the length; the limit's m = GM/c^2 is a length.
    ("neugebauer_meinel", "weyl"): {
        "t": "T", "\\rho": "L", "z": "L", "\\phi": "1", "U": "1", "a": "L", "k": "1",
    },
    ("neugebauer_meinel", "corotating"): {
        "t": "T", "\\rho": "L", "z": "L", "\\varphi": "1", "U": "1", "a": "L", "k": "1",
    },
    ("neugebauer_meinel", "bardeen_wagoner"): {
        "t": "T", "\\rho": "L", "z": "L", "\\phi": "1", "\\nu": "1", "\\omega": "1/L", "\\alpha": "1",
    },
    ("neugebauer_meinel", "spheroidal"): {
        "t": "T", "\\xi": "1", "\\eta": "1", "\\phi": "1", "\\rho_0": "L", "\\nu": "1", "\\omega": "1/L", "\\alpha": "1",
    },
    ("neugebauer_meinel", "black_hole_limit"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "m": "L",
    },
    # Weyl's two functions are pure numbers. In the oblate spheroidal chart xi and eta are pure
    # numbers too, the disc's radius a carries the length, and alpha = arccot(xi) is an angle.
    ("morgan_morgan", "weyl"): {
        "t": "T", "\\rho": "L", "\\phi": "1", "z": "L", "\\psi": "1", "\\gamma": "1",
    },
    ("morgan_morgan", "oblate_spheroidal"): {
        "t": "T", "\\xi": "1", "\\eta": "1", "\\phi": "1", "m": "L", "a": "L", "\\alpha": "1", "\\psi": "1",
        "\\gamma": "1",
    },
    # m = GM/c^2 is a length, and so is R = sqrt(rho^2 + z^2), the name Weyl's chart defines.
    ("curzon_chazy", "weyl"): {
        "t": "T", "\\rho": "L", "\\phi": "1", "z": "L", "m": "L", "R": "L",
    },
    ("curzon_chazy", "spherical"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "m": "L",
    },
    # m is half the length of Weyl's rod; f = 1 - 2m/r and h are the names both charts define.
    ("zipoy_voorhees", "spherical"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "m": "L", "q": "1", "f": "1", "h": "1",
    },
    ("zipoy_voorhees", "prolate_spheroidal"): {
        "t": "T", "x": "1", "y": "1", "\\phi": "1", "m": "L", "\\delta": "1", "f": "1", "h": "1",
    },
    # m = GM/c^2 is a length and the quadrupole parameter q a pure number; L is the logarithm both
    # charts name, and Weyl's two functions are pure numbers.
    ("erez_rosen", "prolate_spheroidal"): {
        "t": "T", "x": "1", "y": "1", "\\phi": "1", "m": "L", "q": "1", "L": "1", "\\psi": "1", "\\gamma": "1",
    },
    ("erez_rosen", "spherical"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "m": "L", "q": "1", "L": "1", "\\psi": "1",
        "\\gamma": "1",
    },
    # Bach and Weyl's ring: m = GM/c^2 and the ring's radius a in Weyl's coordinates are lengths,
    # l_2 is the greatest distance to the ring and kappa the parameter of the elliptic integrals.
    ("bach_weyl_ring", "weyl"): {
        "t": "T", "\\rho": "L", "\\phi": "1", "z": "L", "m": "L", "a": "L", "l_2": "L", "\\kappa": "1",
        "\\psi": "1", "\\gamma": "1",
    },
    ("bach_weyl_ring", "toroidal"): {
        "t": "T", "\\zeta": "1", "\\sigma": "1", "\\phi": "1", "m": "L", "a": "L", "\\kappa": "1",
        "\\psi": "1", "\\gamma": "1",
    },
    ("bach_weyl_ring", "oblate_spheroidal"): {
        "t": "T", "\\xi": "1", "\\eta": "1", "\\phi": "1", "m": "L", "a": "L", "\\kappa": "1", "\\psi": "1",
        "\\gamma": "1",
    },
    # m = GM/c^2 is a length and the quadrupole parameter q of the external field a pure number; U and
    # V, the distortion of Weyl's two functions, are pure numbers, and so are the prolate spheroidal
    # x and y that Weyl's chart names.
    ("distorted_schwarzschild", "prolate_spheroidal"): {
        "t": "T", "x": "1", "y": "1", "\\phi": "1", "m": "L", "q": "1", "U": "1", "V": "1",
    },
    ("distorted_schwarzschild", "spherical"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "m": "L", "q": "1", "U": "1", "V": "1",
    },
    ("distorted_schwarzschild", "weyl"): {
        "t": "T", "\\rho": "L", "\\phi": "1", "z": "L", "m": "L", "q": "1", "x": "1", "y": "1", "\\psi": "1",
        "\\gamma": "1",
    },
    # Sen and Banerjee's black hole in a scalar universe, found again by Cardoso and Natario: m = GM/c^2
    # is a length and the gradient b of the field an inverse length, so b r and b rho are pure numbers,
    # and so is F, the factor e^{-b^2 (r^2 - 2mr) sin^2(theta)} the spherical chart names.
    ("scalar_field_black_hole", "spherical"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "m": "L", "b": "1/L", "F": "1",
    },
    ("scalar_field_black_hole", "scalar_universe"): {
        "t": "T", "\\rho": "L", "\\phi": "1", "z": "L", "b": "1/L",
    },
    ("scalar_field_black_hole", "weyl"): {
        "t": "T", "\\rho": "L", "\\phi": "1", "z": "L", "m": "L", "b": "1/L", "x": "1", "y": "1", "\\psi": "1",
        "\\gamma": "1",
    },
    # m is half the mass as a length and b the dipole moment over the mass, a length; P, Q, Y
    # and Z are Bonnor's four polynomials, each an area.
    ("bonnor_magnetic_dipole", "spheroidal"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "m": "L", "b": "L", "P": "L^2", "Q": "L^2", "Y": "L^2", "Z": "L^2",
    },
    # Fisher, Janis, Newman and Winicour's b is a length and gamma a pure number; f and h are the
    # ratios the charts name. Bronnikov's harmonic coordinate u is an inverse length, e^{-2ku} = 1 - b/r.
    **{("static_five_soliton", system): {
        "t": "T", radius: dimension, "\\theta": "1", "\\phi": "1", "y": "L",
        "r_s": "L", "a": "1", "b": "1", **names,
    } for system, radius, dimension, names in (
        ("spherical", "r", "L", {"f": "1"}),
        ("isotropic", "\\rho", "L", {"h": "1"}),
        ("lake", "h", "1", {}))},
    ("fisher_jnw", "spherical"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "b": "L", "\\gamma": "1", "f": "1",
    },
    ("fisher_jnw", "jnw"): {
        "t": "T", "R": "L", "\\theta": "1", "\\phi": "1", "b": "L", "\\gamma": "1", "f": "1",
    },
    ("fisher_jnw", "isotropic"): {
        "t": "T", "\\rho": "L", "\\theta": "1", "\\phi": "1", "b": "L", "\\gamma": "1", "h": "1",
    },
    ("fisher_jnw", "harmonic"): {
        "t": "T", "u": "L**(-1)", "\\theta": "1", "\\phi": "1", "m": "L", "k": "L",
    },
    # Brans and Dicke's static sphere: B is a length and C and lambda pure numbers in Brans's
    # isotropic chart; r_0 is a length and m and n pure numbers in Campanelli and Lousto's; and
    # Bronnikov's harmonic coordinate u is an inverse length, with k, b and s lengths. h and A are
    # the ratios the charts name.
    ("brans_dicke_sphere", "isotropic"): {
        "t": "T", "\\rho": "L", "\\theta": "1", "\\phi": "1", "B": "L", "C": "1", "\\lambda": "1", "h": "1",
    },
    ("brans_dicke_sphere", "spherical"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "r_0": "L", "m": "1", "n": "1", "A": "1",
    },
    ("brans_dicke_sphere", "harmonic"): {
        "t": "T", "u": "L**(-1)", "\\theta": "1", "\\phi": "1", "k": "L", "b": "L", "s": "L",
    },
    # The exponential metric of Papapetrou and Yilmaz: m = GM/c^2 is a length, the Cartesian chart
    # names its isotropic radius, the areal chart holds the isotropic radius r(R) as a length, and
    # Bronnikov's harmonic coordinate u = 1/r is an inverse length.
    ("exponential_metric", "isotropic"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "m": "L",
    },
    ("exponential_metric", "cartesian"): {
        "t": "T", "x": "L", "y": "L", "z": "L", "m": "L", "r": "L",
    },
    ("exponential_metric", "areal"): {
        "t": "T", "R": "L", "\\theta": "1", "\\phi": "1", "m": "L", "r": "L",
    },
    ("exponential_metric", "harmonic"): {
        "t": "T", "u": "L**(-1)", "\\theta": "1", "\\phi": "1", "m": "L",
    },
    # Bondi and Sachs's radiating metric: r is the luminosity distance and l = 1/r; V stands beside r
    # and is a length, U beside d(theta)/d(cu) and is an inverse length, beta and gamma are exponents.
    ("bondi_sachs", "bondi"): {
        "u": "T", "r": "L", "\\theta": "1", "\\phi": "1", "V": "L", "\\beta": "1", "U": "1/L", "\\gamma": "1",
    },
    ("bondi_sachs", "compactified"): {
        "u": "T", "\\ell": "1/L", "\\theta": "1", "\\phi": "1", "V": "L", "\\beta": "1", "U": "1/L", "\\gamma": "1",
    },
    # The moving mirror of Fulling and Davies: the null coordinates are lengths, the mirror's position
    # z, its ray-tracing function p and the inverse f are lengths, and kappa is an inverse length.
    ("moving_mirror", "inertial"): {"t": "T", "x": "L", "z": "L"},
    ("moving_mirror", "null"): {"u": "L", "v": "L", "p": "L"},
    ("moving_mirror", "mirror_rest"): {"U": "L", "v": "L", "f": "L"},
    ("moving_mirror", "thermal"): {"T": "T", "X": "L", "\\kappa": "1/L"},
    ("moving_mirror", "collapse"): {"T": "T", "X": "L", "\\kappa": "1/L"},
    ("moving_mirror", "rindler"): {"\\eta": "T", "\\xi": "L", "\\kappa": "1/L"},
    # Roberts's collapsing scalar field: the null coordinates u and v are lengths and p a pure number.
    # Roberts's lambda, which the areal chart names, is a length, and Frolov's scaling coordinates are
    # pure numbers, counted in a length ell.
    ("roberts", "double_null"): {"u": "L", "v": "L", "\\theta": "1", "\\phi": "1", "p": "1"},
    ("roberts", "advanced"): {"v": "L", "r": "L", "\\theta": "1", "\\phi": "1", "p": "1"},
    ("roberts", "areal"): {"v": "L", "R": "L", "\\theta": "1", "\\phi": "1", "p": "1", "\\lambda": "L"},
    ("roberts", "diagonal"): {"t": "T", "\\rho": "L", "\\theta": "1", "\\phi": "1", "p": "1"},
    ("roberts", "scaling"): {"\\tau": "1", "x": "1", "\\theta": "1", "\\phi": "1", "p": "1", "\\ell": "L"},
    ("bianchi", "type_i_cartesian"): {
        "t": "T", "x": "L", "y": "L", "z": "L", "a_1": "1", "a_2": "1", "a_3": "1",
    },
    # The exterior keeps G and a mass per unit length explicit, so that the deficit is
    # the dimensionless 4G mu/c^2 the entry prints; delta is the deficit angle itself,
    # which the entry declares and states rather than uses. The interior is written with
    # a dimensionless polar angle chi, so its one length is the radius of the cap.
    # The mass is folded into the length m = GM/c^2, as Taub-NUT folds it, and the
    # acceleration alpha is an inverse length, so alpha r and 2 alpha m are pure numbers; C
    # scales the angle about the axis. The Hong-Teo chart is dimensionless throughout, with
    # tau = alpha ct, y = 1/(alpha r) and x = cos(theta), so its one length is 1/alpha.
    ("c_metric", "spherical"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "m": "L", "\\alpha": "1/L", "C": "1",
    },
    ("c_metric", "hong_teo"): {
        "\\tau": "1", "y": "1", "x": "1", "\\phi": "1", "m": "L", "\\alpha": "1/L", "C": "1",
    },
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
    # The one length is the radius R of the three sphere; the hyperspherical angles are pure
    # numbers, and the areal r and Einstein's projected x, y and z are lengths below R.
    ("einstein_static", "hyperspherical"): {
        "t": "T", "\\chi": "1", "\\theta": "1", "\\phi": "1", "R": "L",
    },
    ("einstein_static", "static_areal"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "R": "L",
    },
    ("einstein_static", "einstein_cartesian"): {
        "t": "T", "x": "L", "y": "L", "z": "L", "R": "L",
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
    # The global monopole keeps Schwarzschild's length r_s = 2GM/c^2 beside the solid angle
    # deficit Delta, a pure number, so 1 - Delta - r_s/r is one too. The Eddington-Finkelstein
    # times u = ct - r_* and v = ct + r_* are lengths.
    ("global_monopole", "static"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "\\Delta": "1", "r_s": "L",
    },
    ("global_monopole", "conical"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "\\Delta": "1",
    },
    ("global_monopole", "eddington_finkelstein_outgoing"): {
        "u": "L", "r": "L", "\\theta": "1", "\\phi": "1", "\\Delta": "1", "r_s": "L",
    },
    ("global_monopole", "eddington_finkelstein_ingoing"): {
        "v": "L", "r": "L", "\\theta": "1", "\\phi": "1", "\\Delta": "1", "r_s": "L",
    },
    # The black hole on a cosmic string keeps Schwarzschild's length r_s beside b = 1 - 4G mu/c^2,
    # a pure number, and the wedge chart's angle b phi is as dimensionless as phi.
    ("string_black_hole", "static"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "r_s": "L", "b": "1",
    },
    ("string_black_hole", "wedge"): {
        "t": "T", "r": "L", "\\theta": "1", "\\tilde\\phi": "1", "r_s": "L", "b": "1",
    },
    ("string_black_hole", "eddington_finkelstein_outgoing"): {
        "u": "L", "r": "L", "\\theta": "1", "\\phi": "1", "r_s": "L", "b": "1",
    },
    ("string_black_hole", "eddington_finkelstein_ingoing"): {
        "v": "L", "r": "L", "\\theta": "1", "\\phi": "1", "r_s": "L", "b": "1",
    },
    # The dilaton black hole keeps Schwarzschild's length r_s beside r_d = Q^2/M in units
    # G = c = 1, a length too, so r(r - r_d) is an area. The string charts are the same chart
    # with the metric multiplied by a power of 1 - r_d/r, a pure number.
    ("dilaton_black_hole", "static"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "r_s": "L", "r_d": "L",
    },
    ("dilaton_black_hole", "eddington_finkelstein_outgoing"): {
        "u": "L", "r": "L", "\\theta": "1", "\\phi": "1", "r_s": "L", "r_d": "L",
    },
    ("dilaton_black_hole", "eddington_finkelstein_ingoing"): {
        "v": "L", "r": "L", "\\theta": "1", "\\phi": "1", "r_s": "L", "r_d": "L",
    },
    ("dilaton_black_hole", "string_magnetic"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "r_s": "L", "r_d": "L",
    },
    ("dilaton_black_hole", "string_electric"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "r_s": "L", "r_d": "L",
    },
    # e^x forces x dimensionless, and with it the other three coordinates, so the whole
    # length of the Godel solution sits in 1/omega.
    ("godel", "cartesian"): {
        "t": "1", "x": "1", "y": "1", "z": "1", "\\omega": "1/L",
    },
    ("godel", "cylindrical"): {
        "t": "1", "r": "1", "\\phi": "1", "z": "1", "\\omega": "1/L",
    },
    ("interior_schwarzschild", "spherical"): {
        "t": "L", "r": "L", "\\theta": "1", "\\phi": "1", "r_s": "L", "R": "L",
    },
    ("kasner", "cartesian"): {
        "t": "T", "x": "L", "y": "L", "z": "L", "p_1": "1", "p_2": "1", "p_3": "1",
    },
    # Kasner's universe with a scalar field: the exponents and the strength q of the field are pure
    # numbers, and a power of the time whose exponent holds one is read with the time in a fixed
    # unit, as Kasner's is. The logarithmic time tau = -ln(t/t_0) is a pure number and ell = c t_0
    # a length, and the chart of five dimensions has the extra coordinate w for a length.
    ("kasner_scalar", "synchronous"): {
        "t": "T", "x": "L", "y": "L", "z": "L", "p_1": "1", "p_2": "1", "p_3": "1", "q": "1",
    },
    ("kasner_scalar", "logarithmic"): {
        "\\tau": "1", "x": "L", "y": "L", "z": "L", "p_1": "1", "p_2": "1", "p_3": "1", "q": "1", "\\ell": "L",
    },
    ("kasner_scalar", "kaluza_klein"): {
        "T": "T", "x": "L", "y": "L", "z": "L", "w": "L", "s_1": "1", "s_2": "1", "s_3": "1", "s_5": "1",
    },
    # Kasner's universe with a magnetic field: the exponents and the strength b of the field are
    # pure numbers, and the powers of the time are read with the time in a fixed unit, as Kasner's
    # are, so that 1 + b^2 t^(2 p_3) balances. Rosen's time eta is a pure number and ell a length.
    ("kasner_magnetic", "kasner_time"): {
        "t": "T", "x": "L", "y": "L", "z": "L", "p_1": "1", "p_2": "1", "p_3": "1", "b": "1",
    },
    ("kasner_magnetic", "rosen"): {
        "\\eta": "1", "x": "L", "y": "L", "z": "L", "\\ell": "L",
    },
    # Many black holes and wormholes in three dimensions: the one length is the anti-de Sitter
    # radius. The disc's rho and the stereographic x and y are lengths, as Brill's lectures write
    # them, and M is the pure number of the Banados-Teitelboim-Zanelli chart.
    ("btz_multi_holes_wormholes", "sausage"): {"t": "T", "\\rho": "L", "\\phi": "1", "\\ell": "L"},
    ("btz_multi_holes_wormholes", "stereographic"): {"\\tau": "T", "x": "L", "y": "L", "\\ell": "L"},
    ("btz_multi_holes_wormholes", "free_fall"): {"T": "T", "\\rho": "L", "\\phi": "1", "\\ell": "L"},
    ("btz_multi_holes_wormholes", "exterior"): {"t": "T", "r": "L", "\\phi": "1", "\\ell": "L", "M": "1"},
    # The shock wave in the BTZ black hole: Shenker and Stanford's Kruskal u and v and their
    # discontinuous U and V are pure numbers, the anti-de Sitter radius l and the radius R of the
    # horizon are lengths, and the shift alpha is a pure number, as is the step Theta.
    ("btz_shock_wave", "kruskal"): {"u": "1", "v": "1", "\\phi": "1", "\\ell": "L", "R": "L", "\\alpha": "1",
                                    "\\Theta": "1"},
    ("btz_shock_wave", "discontinuous"): {"U": "1", "V": "1", "\\phi": "1", "\\ell": "L", "R": "L", "\\alpha": "1"},
    ("btz_shock_wave", "exterior"): {"t": "T", "r": "L", "\\phi": "1", "\\ell": "L", "R": "L"},
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
    # Khan and Penrose's null coordinates, and tau and sigma, are pure numbers, measured in units
    # of the waves' focal length L, which multiplies the part of the metric along them.
    ("khan_penrose", "double_null"): {"u": "1", "v": "1", "x": "L", "y": "L", "L": "L"},
    ("khan_penrose", "cosmological"): {"\\tau": "1", "\\sigma": "1", "x": "L", "y": "L", "L": "L"},
    # Bell and Szekeres's null coordinates are lengths and the two wave strengths inverse lengths;
    # xi, eta and the coordinates of the regular, global and Bertotti-Robinson charts are pure
    # numbers, with 1/(2ab) carrying the length squared.
    ("bell_szekeres", "double_null"): {"u": "L", "v": "L", "x": "L", "y": "L", "a": "1/L", "b": "1/L"},
    ("bell_szekeres", "time_space"): {"\\xi": "1", "\\eta": "1", "x": "L", "y": "L", "a": "1/L", "b": "1/L"},
    ("bell_szekeres", "regular"): {"T": "1", "Z": "1", "X": "1", "Y": "1", "a": "1/L", "b": "1/L"},
    ("bell_szekeres", "global"): {"\\chi": "1", "\\rho": "1", "\\theta": "1", "\\phi": "1", "a": "1/L", "b": "1/L"},
    ("bell_szekeres", "kruskal_szekeres"): {"U": "L", "V": "L", "\\eta": "1", "x": "L", "a": "1/L", "b": "1/L"},
    ("bell_szekeres", "bertotti_robinson"): {"t": "1", "r": "1", "\\theta": "1", "\\phi": "1", "a": "1/L", "b": "1/L"},
    # Chandrasekhar and Xanthopoulos's eta and mu, and psi and lambda, are pure
    # numbers; x and y are lengths along the wave fronts and m is the length that multiplies the
    # part of the metric across them. In Kerr's charts m and a are lengths, and t and v are
    # lengths that stand along a spacelike Killing vector between the horizons.
    ("chandrasekhar_xanthopoulos", "prolate"): {
        "\\eta": "1", "\\mu": "1", "x": "L", "y": "L", "m": "L", "\\alpha": "1", "p": "1", "q": "1",
        "\\rho": "1", "X": "1", "Y": "1",
    },
    ("chandrasekhar_xanthopoulos", "angular"): {
        "\\psi": "1", "\\lambda": "1", "x": "L", "y": "L", "m": "L", "\\alpha": "1", "p": "1", "q": "1",
        "\\rho": "1", "X": "1", "Y": "1",
    },
    ("chandrasekhar_xanthopoulos", "boyer_lindquist"): {
        "t": "L", "r": "L", "\\theta": "1", "\\phi": "1", "m": "L", "a": "L", "\\Sigma": "L**2", "\\Delta": "L**2",
    },
    ("chandrasekhar_xanthopoulos", "kerr_ingoing"): {
        "v": "L", "r": "L", "\\theta": "1", "\\tilde\\phi": "1", "m": "L", "a": "L", "\\Sigma": "L**2",
        "\\Delta": "L**2",
    },
    ("krasnikov", "cylindrical"): {
        "t": "T", "x": "L", "r": "L", "\\phi": "1", "k": "1",
    },
    # The three Euler angles of the three sphere are pure numbers, so the lengths of the
    # closed universe sit in the three scale factors, as a radius does.
    ("mixmaster", "euler_angles"): {
        "t": "T", "\\psi": "1", "\\theta": "1", "\\phi": "1", "a_1": "L", "a_2": "L", "a_3": "L",
    },
    # The Kantowski-Sachs cylinder: r is a length along its axis, so a is a pure number, and
    # the spheres' radius b carries the length. Inside Schwarzschild's horizon the time T is
    # the areal radius, a length, and the dust's parametric time eta is a pure number.
    ("kantowski_sachs", "comoving"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "a": "1", "b": "L",
    },
    ("kantowski_sachs", "schwarzschild_interior"): {
        "T": "L", "r": "L", "\\theta": "1", "\\phi": "1", "r_s": "L",
    },
    ("kantowski_sachs", "dust"): {
        "\\eta": "1", "r": "L", "\\theta": "1", "\\phi": "1", "b_0": "L", "\\kappa": "1",
    },
    # Minkowski's inertial chart with its metric multiplied by the square of a conformal
    # factor, which multiplies proper time and so has to be a pure number.
    # Aichelburg and Sexl keep G and the energy E of the source explicit, so 8GE/c^4 is the one
    # length of the geometry and rho_0 an arbitrary second one; u = ct - z and v = ct + z are
    # lengths, and each Dirac delta carries the inverse of its argument.
    ("aichelburg_sexl", "cartesian"): {
        "t": "T", "x": "L", "y": "L", "z": "L", "G": "L**3/(M*T**2)", "E": "M*L**2/T**2", "\\rho_0": "L",
    },
    ("aichelburg_sexl", "null_cartesian"): {
        "u": "L", "v": "L", "x": "L", "y": "L", "G": "L**3/(M*T**2)", "E": "M*L**2/T**2", "\\rho_0": "L",
    },
    ("aichelburg_sexl", "null_cylindrical"): {
        "u": "L", "v": "L", "\\rho": "L", "\\phi": "1", "G": "L**3/(M*T**2)", "E": "M*L**2/T**2",
        "\\rho_0": "L",
    },
    # Hotta and Tanaka's shock in de Sitter space of radius a, a length: 8GE/c^4 is the other
    # length. The conformal time and the radius of the conformally flat chart are lengths, the
    # global chart's are angles, and the null coordinates of the others are lengths.
    ("hotta_tanaka", "conformally_flat"): {
        "\\eta": "L", "\\rho": "L", "\\theta": "1", "\\phi": "1", "G": "L**3/(M*T**2)", "E": "M*L**2/T**2", "a": "L",
    },
    ("hotta_tanaka", "global"): {
        "\\eta": "1", "\\chi": "1", "\\theta": "1", "\\phi": "1", "G": "L**3/(M*T**2)", "E": "M*L**2/T**2", "a": "L",
    },
    ("hotta_tanaka", "kruskal"): {
        "u": "L", "v": "L", "\\theta": "1", "\\phi": "1", "G": "L**3/(M*T**2)", "E": "M*L**2/T**2", "a": "L",
    },
    ("hotta_tanaka", "null_cylindrical"): {
        "u": "L", "v": "L", "\\rho": "L", "\\phi": "1", "G": "L**3/(M*T**2)", "E": "M*L**2/T**2", "a": "L",
    },
    ("hotta_tanaka", "kundt"): {
        "u": "L", "w": "L", "\\theta": "1", "\\phi": "1", "G": "L**3/(M*T**2)", "E": "M*L**2/T**2", "a": "L",
    },
    # Bonnor's profile A multiplies (c dt - dz)^2 and is a pure number; his u and v, with
    # sqrt(2) u = ct - z, are lengths, and 8 pi G epsilon/c^4, with epsilon an energy density,
    # is an inverse area, so the uniform beam's A is a pure number as well.
    ("light_beam", "cartesian"): {"t": "T", "x": "L", "y": "L", "z": "L", "A": "1"},
    ("light_beam", "null_cartesian"): {"u": "L", "v": "L", "x": "L", "y": "L", "A": "1"},
    ("light_beam", "null_cylindrical_interior"): {
        "u": "L", "v": "L", "\\rho": "L", "\\phi": "1", "G": "L**3/(M*T**2)", "\\epsilon": "M/(L*T**2)", "R": "L",
    },
    ("light_beam", "null_cylindrical_exterior"): {
        "u": "L", "v": "L", "\\rho": "L", "\\phi": "1", "G": "L**3/(M*T**2)", "\\epsilon": "M/(L*T**2)", "R": "L",
    },
    # k = 2 pi G sigma / c^4 is an inverse length, 1/k the radius of the wall when it stops.
    ("domain_wall", "planar"): {"t": "T", "x": "L", "y": "L", "z": "L", "k": "1/L"},
    ("domain_wall", "global"): {"t": "T", "z": "L", "\\theta": "1", "\\phi": "1", "k": "1/L"},
    ("domain_wall", "conformal"): {"t": "T", "w": "L", "\\theta": "1", "\\phi": "1", "k": "1/L"},
    ("domain_wall", "inertial"): {"T": "T", "R": "L", "\\theta": "1", "\\phi": "1", "k": "1/L"},
    # Coleman and De Luccia's bubble: the rapidities psi and chi and the conformal time eta are pure
    # numbers, the proper distance xi and the proper time tau carry the dimensions, and rho and
    # the scale factor a are lengths. Each static chart has the cosmological constant of its
    # side and the wall's areal radius, a length.
    ("coleman_de_luccia", "wall"): {"\\psi": "1", "\\xi": "L", "\\theta": "1", "\\phi": "1", "\\rho": "L"},
    ("coleman_de_luccia", "open"): {"\\tau": "T", "\\chi": "1", "\\theta": "1", "\\phi": "1", "a": "L"},
    ("coleman_de_luccia", "open_conformal"): {"\\eta": "1", "\\chi": "1", "\\theta": "1", "\\phi": "1", "a": "L"},
    ("coleman_de_luccia", "static_inside"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "\\Lambda_T": "1/L**2", "r_w": "L",
    },
    ("coleman_de_luccia", "static_outside"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "\\Lambda_F": "1/L**2", "r_w": "L",
    },
    ("malament_hogarth", "cartesian"): {
        "t": "T", "x": "L", "y": "L", "z": "L", "\\Omega": "1",
    },
    # The shift is the gradient of a potential and a flow in units of c, so the potential is
    # a length and each first derivative of it a pure number, as Lentz's N_i is with c = 1.
    ("lentz", "cartesian"): {
        "t": "T", "x": "L", "y": "L", "z": "L", "\\phi": "L",
    },
    # McVittie's mass is folded into r_s = 2GM/c^2, and his comoving r is a length with the
    # scale factor a pure number, as FRW's are, so r_s/4ar is one too. The areal chart carries
    # the Hubble rate H = (da/dt)/a, a frequency as de Sitter's flat slicing carries it, so HR/c
    # is a pure number.
    ("mcvittie", "isotropic"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "r_s": "L", "a": "1",
    },
    ("mcvittie", "areal"): {
        "t": "T", "R": "L", "\\theta": "1", "\\phi": "1", "r_s": "L", "H": "1/T",
    },
    # The Milne universe has no length of its own: the comoving chi and r are pure numbers and
    # ct carries the length, and the logarithmic time tau needs the time t_0 at which it is zero.
    ("milne", "comoving_hyperbolic"): {
        "t": "T", "\\chi": "1", "\\theta": "1", "\\phi": "1",
    },
    ("milne", "comoving_spherical"): {
        "t": "T", "r": "1", "\\theta": "1", "\\phi": "1",
    },
    ("milne", "logarithmic_time"): {
        "\\tau": "T", "\\chi": "1", "\\theta": "1", "\\phi": "1", "t_0": "T",
    },
    ("milne", "inertial"): {
        "T": "T", "R": "L", "\\theta": "1", "\\phi": "1",
    },
    ("minkowski", "cartesian"): {"t": "T", "x": "L", "y": "L", "z": "L"},
    ("minkowski", "spherical"): {"t": "T", "r": "L", "\\theta": "1", "\\phi": "1"},
    ("minkowski", "double_null"): {"u": "T", "v": "T", "y": "L", "z": "L"},
    ("minkowski", "spherical_null"): {"u": "T", "v": "T", "\\theta": "1", "\\phi": "1"},
    ("minkowski", "rindler"): {"T": "T", "X": "L", "Y": "L", "Z": "L", "a": "L/T**2"},
    # Misner space names no length of its own, since Minkowski space identified under a boost
    # is identified under every dilation too. Misner's T is then an area and psi a pure number;
    # the rapidity chi of the Milne chart and eta of the Rindler chart are pure numbers as well.
    ("misner", "misner"): {"T": "L**2", "\\psi": "1", "y": "L", "z": "L", "\\psi_0": "1"},
    ("misner", "milne"): {"t": "T", "\\chi": "1", "y": "L", "z": "L", "\\psi_0": "1"},
    ("misner", "rindler"): {"\\eta": "1", "\\xi": "L", "y": "L", "z": "L", "\\psi_0": "1"},
    # Gott's two strings: the centre of momentum chart and one string's conical chart carry the
    # strings' mass per length, speed and distance from the plane between them; the Rindler and
    # Milne charts of Grant's covering space carry the rapidity a of the boost round both strings,
    # a pure number, and the shift b along its axis, a length.
    ("gott_time_machine", "centre_of_momentum"): {
        "t": "T", "x": "L", "y": "L", "z": "L",
        "\\mu": "M/L", "G": "L**3/(M*T**2)", "v": "L/T", "d": "L", "\\alpha": "1", "\\gamma": "1",
    },
    ("gott_time_machine", "string_rest"): {
        "t": "T", "r": "L", "\\phi": "1", "z": "L", "\\mu": "M/L", "G": "L**3/(M*T**2)", "d": "L",
    },
    ("gott_time_machine", "grant_rindler"): {"\\eta": "1", "\\xi": "L", "Y": "L", "z": "L", "a": "1", "b": "L"},
    ("gott_time_machine", "grant_milne"): {"\\tau": "T", "\\chi": "1", "Y": "L", "z": "L", "a": "1", "b": "L"},
    # The spinning string's two parameters are the length a = 4GJ/c^3 and the pure number
    # b = 1 - 4G mu/c^2; the helical chart's angle b phi is as dimensionless as phi, and the
    # extended source's M(r) and rho(r) are lengths, as is the proper radius r_0 of its surface.
    # Ori's z is a pure number, periodic with period L, so dz dT makes T an area, as Misner's T is,
    # and f - T makes f one too; the numbers a and e of his example are then pure, which
    # e > (2e + a)^2 needs. In the Brinkmann chart u = -2 exp(-z/2) is pure and v = T exp(z/2) an area.
    ("ori_time_machine", "vacuum_core"): {"T": "L**2", "x": "L", "y": "L", "z": "1", "f": "L**2", "L": "1"},
    ("ori_time_machine", "foliation"): {"t": "L**2", "x": "L", "y": "L", "z": "1", "a": "1", "e": "1", "L": "1"},
    ("ori_time_machine", "brinkmann"): {"u": "1", "v": "L**2", "x": "L", "y": "L", "a": "1", "L": "1"},
    ("spinning_string", "proper_radius"): {"t": "T", "r": "L", "\\phi": "1", "z": "L", "a": "L", "b": "1"},
    ("spinning_string", "rescaled_radius"): {"t": "T", "\\rho": "L", "\\phi": "1", "z": "L", "a": "L", "b": "1"},
    ("spinning_string", "circumference_radius"): {"t": "T", "R": "L", "\\phi": "1", "z": "L", "a": "L", "b": "1"},
    ("spinning_string", "helical"): {"\\tau": "T", "r": "L", "\\tilde\\phi": "1", "z": "L", "a": "L", "b": "1"},
    ("spinning_string", "extended_source"): {
        "t": "T", "r": "L", "\\phi": "1", "z": "L", "M": "L", "\\rho": "L", "r_0": "L",
    },
    # The throat of extreme Kerr with one length, r_0^2 = 2GJ/c^3. Bardeen and Horowitz's
    # Poincare-type chart keeps a time and a length, the inverse radius x = r_0^2/r is a length
    # too, their global chart's tau and y are pure numbers with r_0^2 out in front, and the
    # near-NHEK chart adds the length k = pi r_0^2 T of its horizon's temperature.
    ("near_horizon_extreme_kerr", "poincare"): {"t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "r_0": "L"},
    ("near_horizon_extreme_kerr", "inverse_radius"): {
        "t": "T", "x": "L", "\\theta": "1", "\\phi": "1", "r_0": "L",
    },
    ("near_horizon_extreme_kerr", "global"): {"\\tau": "1", "y": "1", "\\theta": "1", "\\phi": "1", "r_0": "L"},
    ("near_horizon_extreme_kerr", "near_nhek"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "r_0": "L", "k": "L",
    },
    # Senovilla's a is an inverse length, so act and 3a rho are pure numbers; g_phiphi carries the
    # area 1/(9a^2) and z is a length.
    ("senovilla", "cylindrical"): {
        "t": "T", "\\rho": "L", "\\phi": "1", "z": "L", "a": "1/L",
    },
    # The bounce of loop quantum cosmology keeps c. Its bounce time t_b is a time, its scale factor
    # a, named for (1 + t^2/t_b^2)^(1/6) or cosh^(1/3)(tau/t_b), a pure number, and x, y and z
    # comoving lengths.
    ("lqc_bounce", "cosmic"): {"t": "T", "x": "L", "y": "L", "z": "L", "t_b": "T", "a": "1"},
    ("lqc_bounce", "comoving_spherical"): {"t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "t_b": "T", "a": "1"},
    ("lqc_bounce", "harmonic"): {"\\tau": "T", "x": "L", "y": "L", "z": "L", "t_b": "T", "a": "1"},
    # The scale factor is the time of two of Murphy's charts, a pure number, so no coordinate of
    # theirs is a time the chart multiplies by c.
    ("murphy_viscous_universe", "scale_factor"): {"a": "1", "x": "L", "y": "L", "z": "L", "q": "1", "H_s": "1/T"},
    ("murphy_viscous_universe", "scale_factor_spherical"): {
        "a": "1", "r": "L", "\\theta": "1", "\\phi": "1", "q": "1", "H_s": "1/T",
    },
    ("murphy_viscous_universe", "cosmic"): {"t": "T", "x": "L", "y": "L", "z": "L", "q": "1", "H_s": "1/T", "a": "1"},
    # The travelling wave on a string: u = ct - z and v = ct + z are lengths, so no coordinate is a
    # time and the profile F is a pure number, as is b = 1 - 4G mu/c^2. The string's displacements
    # A(u) and B(u) are lengths, as are the isotropic x and y, their distance rho from the string
    # and the arbitrary length ell that rho is measured in.
    ("string_wave", "null_conical"): {"u": "L", "v": "L", "r": "L", "\\phi": "1", "b": "1", "F": "1"},
    ("string_wave", "isotropic"): {
        "u": "L", "v": "L", "x": "L", "y": "L", "b": "1", "\\ell": "L", "A": "L", "B": "L", "\\rho": "L",
    },
    ("string_wave", "moving_string"): {
        "u": "L", "V": "L", "X": "L", "Y": "L", "b": "1", "\\ell": "L", "A": "L", "B": "L", "\\rho": "L",
    },
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
    # Teo's example with the throat radius restored: the spin a = GJ/(c^3 b_0^2) is a pure
    # number, and in the proper distance chart r is a declared function of l.
    ("teo_wormhole", "spherical"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "b_0": "L", "a": "1",
    },
    ("teo_wormhole", "proper_radial"): {
        "t": "T", "l": "L", "\\theta": "1", "\\phi": "1", "b_0": "L", "a": "1", "r": "L",
    },
    # Morris, Thorne and Yurtsever's wormhole with one mouth accelerating: the flat space outside the
    # mouths in its Lorentz chart, their own chart through the wormhole, whose acceleration g of
    # the right mouth is a function of t and whose lapse N = 1 + g l F cos(theta)/c^2 is a name
    # the chart defines, and the mouth of Friedman and his coauthors' throat of zero length.
    ("wormhole_time_machine", "lorentz"): {"T": "T", "X": "L", "Y": "L", "Z": "L", "b": "L"},
    ("wormhole_time_machine", "wormhole"): {
        "t": "T", "l": "L", "\\theta": "1", "\\phi": "1",
        "g": "L/T**2", "F": "1", "\\Phi": "1", "r": "L", "N": "1",
    },
    ("wormhole_time_machine", "short_throat"): {"t": "T", "l": "L", "\\theta": "1", "\\phi": "1", "b": "L"},
    # Tippett and Tsang's bubble: the top hat function h is a pure number, and lambda, the polar
    # angle of the plane of ct and x, is the time of the two charts that use it.
    ("tippett_tsang", "cartesian"): {"t": "T", "x": "L", "y": "L", "z": "L", "h": "1"},
    ("tippett_tsang", "polar"): {"\\lambda": "1", "\\xi": "L", "y": "L", "z": "L", "h": "1"},
    ("tippett_tsang", "interior"): {"t": "T", "x": "L", "y": "L", "z": "L"},
    ("tippett_tsang", "rindler"): {"\\lambda": "1", "\\xi": "L", "y": "L", "z": "L"},
    # Nordstrom's conformal factor is a pure number, 1 far from every mass. The point mass is the
    # length m = GM/c^2, the uniform field an acceleration, and the dust universe the length L,
    # the light travel distance of half its span of the inertial time.
    ("nordstrom_scalar", "conformal"): {"t": "T", "x": "L", "y": "L", "z": "L", "\\Phi": "1"},
    ("nordstrom_scalar", "spherical"): {"t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "m": "L"},
    ("nordstrom_scalar", "uniform"): {"t": "T", "x": "L", "y": "L", "z": "L", "a": "L/T**2"},
    ("nordstrom_scalar", "dust"): {"t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "L": "L"},
    # Einstein's speed of light of 1912 is c N, with N a pure number that is 1 far from every mass.
    # The uniform field is an acceleration and a body the length m = GM/c^2.
    ("einstein_1912_static", "static"): {"t": "T", "x": "L", "y": "L", "z": "L", "N": "1"},
    ("einstein_1912_static", "uniform"): {"t": "T", "x": "L", "y": "L", "z": "L", "a": "L/T**2"},
    ("einstein_1912_static", "february"): {"t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "m": "L"},
    ("einstein_1912_static", "march"): {"t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "m": "L"},
    # The slices are flat space and the whole of the geometry is the flow field carried on them,
    # so its components are velocities and the chart components of the metric are powers of V/c.
    ("natario", "cartesian_flow"): {
        "t": "T", "x": "L", "y": "L", "z": "L",
        "u": "L/T", "v": "L/T", "w": "L/T",
    },
    ("natario", "plane_flow"): {
        "t": "T", "x": "L", "y": "L", "z": "L", "u": "L/T",
    },
    # Plebanski and Hacyan's products and anti-Nariai: a and b are the radii of the factors, tau a
    # boost angle or a dimensionless time, and f and g, which multiply x and y in g_uu, inverse lengths.
    # Lindquist and Wheeler's lattice keeps Schwarzschild's length r_s = 2GM/c^2 for the mass of
    # one cell. psi, the angular radius of a cell on the comparison hypersphere, and E = cos^2(psi)
    # are pure numbers; the comoving label rho is the largest radius a shell reaches, a length, so
    # that d_rho r is a pure number; and the hypersphere's radius a carries the length, chi an angle.
    ("lindquist_wheeler_lattice", "schwarzschild_cell"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "r_s": "L", "\\psi": "1"},
    ("lindquist_wheeler_lattice", "cosmological_time"): {
        "\\tau": "T", "r": "L", "\\theta": "1", "\\phi": "1", "r_s": "L", "E": "1"},
    ("lindquist_wheeler_lattice", "lindquist_wheeler"): {
        "\\tau": "T", "\\rho": "L", "\\theta": "1", "\\phi": "1", "r_s": "L", "r": "L"},
    ("lindquist_wheeler_lattice", "comparison_hypersphere"): {
        "\\tau": "T", "\\chi": "1", "\\theta": "1", "\\phi": "1", "a": "L", "r_s": "L", "\\psi": "1"},
    ("plebanski_hacyan", "sphere"): {"t": "T", "z": "L", "\\theta": "1", "\\phi": "1", "b": "L"},
    ("plebanski_hacyan", "sphere_rindler"): {"\\tau": "1", "\\chi": "L", "\\theta": "1", "\\phi": "1", "b": "L"},
    ("plebanski_hacyan", "plane"): {"u": "L", "w": "L", "x": "L", "y": "L", "a": "L"},
    ("plebanski_hacyan", "plane_null"): {"u": "L", "v": "L", "x": "L", "y": "L", "a": "L"},
    ("plebanski_hacyan", "plane_static"): {"\\tau": "1", "\\chi": "1", "x": "L", "y": "L", "a": "L"},
    ("plebanski_hacyan", "anti_nariai"): {"\\tau": "1", "\\chi": "1", "\\theta": "1", "\\phi": "1", "a": "L"},
    ("plebanski_hacyan", "anti_nariai_static"): {"t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "a": "L"},
    ("plebanski_hacyan", "exceptional"): {"u": "L", "w": "L", "x": "L", "y": "L", "a": "L", "f": "1/L", "g": "1/L"},
    # Cremmer and Scherk's Minkowski space times a sphere: a is the sphere's radius.
    ("cremmer_scherk", "cartesian"): {"t": "T", "x": "L", "y": "L", "z": "L", "\\theta": "1", "\\phi": "1", "a": "L"},
    # The three-brane: L is the radius of its throat, and H, Gibbons, Horowitz and Townsend's w and
    # the five angles of the sphere are pure numbers.
    **{("three_brane_throat", chart): {
        "t": "T", "x": "L", "y": "L", "z": "L", radial: unit, "\\alpha": "1", "\\beta": "1", "\\psi": "1",
        "\\theta": "1", "\\phi": "1", "L": "L", **extra}
       for chart, radial, unit, extra in (("isotropic", "\\rho", "L", {"H": "1"}), ("areal", "r", "L", {}),
                                          ("horizon", "w", "1", {"r": "L"}), ("throat", "r", "L", {}),
                                          ("throat_proper", "\\sigma", "L", {}))},
    # Freund and Rubin's anti-de Sitter space times a seven-sphere: L is the radius of anti-de Sitter
    # space, the conformal chart's chi and every angle are pure numbers, and the global chart's rho too.
    # Lin, Lunin and Maldacena's bubbling anti-de Sitter space: the time t is a pure number in every
    # chart and nothing is multiplied by c; the plane of the droplets and y carry an area, G is a
    # pure number and so is the twist V; the plane wave's radii r_1 and r_2 are lengths.
    ("bubbling_ads", "rings"): {
        "t": "1", "r": "L**2", "\\phi": "1", "y": "L**2", "\\alpha": "1", "\\beta": "1", "\\gamma": "1",
        "\\kappa": "1", "\\xi": "1", "\\omega": "1", "G": "1", "V": "1",
    },
    ("bubbling_ads", "global"): {
        "t": "1", "\\rho": "1", "\\alpha": "1", "\\beta": "1", "\\gamma": "1", "\\theta": "1", "\\psi": "1",
        "\\kappa": "1", "\\xi": "1", "\\omega": "1", "L": "L",
    },
    ("bubbling_ads", "plane_wave"): {
        "t": "1", "x": "L**2", "r_1": "L", "r_2": "L", "\\alpha": "1", "\\beta": "1", "\\gamma": "1",
        "\\kappa": "1", "\\xi": "1", "\\omega": "1",
    },
    **{("freund_rubin", chart): {
        "t": "T", **{c: "L" for c in large}, radial: unit, **({} if large else {"\\theta": "1", "\\phi": "1"}), "\\alpha": "1",
        "\\beta": "1", "\\gamma": "1", "\\kappa": "1", "\\xi": "1", "\\omega": "1", "\\psi": "1", "L": "L"}
       for chart, radial, unit, large in (("global", "\\rho", "1", ()), ("conformal", "\\chi", "1", ()),
                                          ("static", "r", "L", ()), ("poincare", "r", "L", ("x", "y")),
                                          ("proper", "\\sigma", "L", ("x", "y")))},
    # The maximally supersymmetric plane wave: mu is an inverse length; the light cone time u and
    # the conformally flat chart's U are times, every other coordinate a length.
    ("maximally_supersymmetric_plane_wave", "brinkmann"): {
        "u": "T", "v": "L", **{f"x_{i}": "L" for i in range(1, 9)}, "\\mu": "1/L"},
    ("maximally_supersymmetric_plane_wave", "rosen"): {
        "u": "T", "v": "L", **{f"y_{i}": "L" for i in range(1, 9)}, "\\mu": "1/L"},
    ("maximally_supersymmetric_plane_wave", "conformally_flat"): {
        "U": "T", "V": "L", **{f"X_{i}": "L" for i in range(1, 9)}, "\\mu": "1/L"},
    # Cahen and Wallach's waves: the eigenvalues are inverse lengths squared, omega and k inverse
    # lengths; u is a time, every other coordinate a length.
    ("cahen_wallach", "brinkmann"): {
        "u": "T", "v": "L", "x": "L", "y": "L", "\\lambda_1": "1/L**2", "\\lambda_2": "1/L**2"},
    ("cahen_wallach", "rotating"): {"u": "T", "v": "L", "x": "L", "y": "L", "\\omega": "1/L"},
    ("cahen_wallach", "rosen"): {"u": "T", "v": "L", "p": "L", "q": "L", "k": "1/L"},
    # Hauser's coordinates and his function are pure numbers, and the one length L multiplies the
    # whole line element, so no coordinate is a time the chart multiplies by c.
    ("hauser_wave", "kerr_debney_schild"): {"t": "1", "r": "1", "x": "1", "y": "1", "L": "L", "f": "1", "H": "1"},
    # The conformal chart's eta and chi are angles; 1/Lambda carries the length squared.
    ("nariai", "static"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "\\Lambda": "1/L**2",
    },
    ("nariai", "global"): {
        "t": "T", "\\chi": "1", "\\theta": "1", "\\phi": "1", "\\Lambda": "1/L**2",
    },
    ("nariai", "conformal"): {
        "\\eta": "1", "\\chi": "1", "\\theta": "1", "\\phi": "1", "\\Lambda": "1/L**2",
    },
    # Elliptic de Sitter space: the one length is the de Sitter radius; the global chart's chi, the
    # conformal chart's eta and chi and the Kruskal chart's U and V are pure numbers.
    ("elliptic_de_sitter", "global"): {"t": "T", "\\chi": "1", "\\theta": "1", "\\phi": "1", "\\ell": "L"},
    ("elliptic_de_sitter", "conformal"): {"\\eta": "1", "\\chi": "1", "\\theta": "1", "\\phi": "1", "\\ell": "L"},
    ("elliptic_de_sitter", "kruskal"): {"U": "1", "V": "1", "\\theta": "1", "\\phi": "1", "\\ell": "L"},
    ("elliptic_de_sitter", "static"): {"t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "\\ell": "L"},
    ("elliptic_de_sitter", "planar"): {"t": "T", "x": "L", "y": "L", "z": "L", "\\ell": "L"},
    # Gott and Li's de Sitter space identified under a boost: beta is a length, the period of ct
    # in the static chart and of l in the Kantowski-Sachs chart, and r_0 the de Sitter radius.
    ("self_creating_universe", "static"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "r_0": "L", "\\beta": "L"},
    ("self_creating_universe", "kantowski_sachs"): {
        "\\tau": "T", "l": "L", "\\theta": "1", "\\phi": "1", "r_0": "L", "\\beta": "L"},
    ("self_creating_universe", "steady_state"): {
        "\\tau": "T", "x": "L", "y": "L", "z": "L", "r_0": "L", "\\beta": "L"},
    ("self_creating_universe", "conformal"): {
        "\\eta": "L", "\\rho": "L", "\\theta": "1", "\\phi": "1", "r_0": "L", "\\beta": "L"},
    # The universe from nothing: l is the radius of the 4-sphere and of de Sitter's waist. The
    # 4-sphere's tau and the scale factor a are lengths, with no c, and the conformal eta is a pure number.
    ("universe_from_nothing", "closed"): {"t": "T", "\\chi": "1", "\\theta": "1", "\\phi": "1", "\\ell": "L"},
    ("universe_from_nothing", "four_sphere"): {"\\tau": "L", "\\chi": "1", "\\theta": "1", "\\phi": "1", "\\ell": "L"},
    ("universe_from_nothing", "scale_factor"): {"a": "L", "\\chi": "1", "\\theta": "1", "\\phi": "1", "\\ell": "L"},
    ("universe_from_nothing", "conformal"): {"\\eta": "1", "\\chi": "1", "\\theta": "1", "\\phi": "1", "\\ell": "L"},
    ("universe_from_nothing", "lapse"): {"t": "T", "\\chi": "1", "\\theta": "1", "\\phi": "1", "\\ell": "L"},
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
    # Bonnor and Faulkes's sphere: alpha and F(t) are pure numbers, and so are Y and Z, since
    # r^2/4L^2 is; the mass m = GM/c^2 and the surface r_0 are lengths. Thompson and Whitrow's
    # B is a pure number and C an inverse area, so that B + Cr^2 is one.
    ("oscillating_sphere", "comoving"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1",
        "\\alpha": "1", "L": "L", "m": "L", "r_0": "L", "F": "1", "Y": "1", "Z": "1",
    },
    ("oscillating_sphere", "central_time"): {
        "\\tau": "T", "r": "L", "\\theta": "1", "\\phi": "1", "B": "1", "C": "1/L^2",
    },
    ("oscillating_sphere", "exterior_schwarzschild"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "m": "L",
    },
    # The semiclosed world is the same construction with the surface past the equator of the three
    # sphere. Its conformal time eta is an angle, with the cycloid written out in a_m, and the
    # isotropic radius is a length.
    ("semiclosed_world", "comoving"): {
        "\\tau": "T", "\\chi": "1", "\\theta": "1", "\\phi": "1",
        "a": "L", "\\chi_0": "1", "a_m": "L",
    },
    ("semiclosed_world", "conformal"): {
        "\\eta": "1", "\\chi": "1", "\\theta": "1", "\\phi": "1", "a_m": "L", "\\chi_0": "1",
    },
    ("semiclosed_world", "schwarzschild"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "r_s": "L",
    },
    ("semiclosed_world", "isotropic"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "r_s": "L",
    },
    # The white hole is that collapse with the time reversed, and its charts carry the same
    # dimensions: inside, chi is an angle and the scale factor the length, and the conformal time
    # eta of the cycloid a = a_m sin^2(eta/2) is an angle too; outside, Schwarzschild's chart and
    # the outgoing Eddington-Finkelstein chart, whose retarded time u is a length.
    ("white_hole", "interior_comoving"): {
        "\\tau": "T", "\\chi": "1", "\\theta": "1", "\\phi": "1",
        "a": "L", "\\chi_0": "1", "a_m": "L",
    },
    ("white_hole", "interior_conformal"): {
        "\\eta": "1", "\\chi": "1", "\\theta": "1", "\\phi": "1", "\\chi_0": "1", "a_m": "L",
    },
    ("white_hole", "exterior_schwarzschild"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "r_s": "L",
    },
    ("white_hole", "exterior_eddington_finkelstein"): {
        "u": "L", "r": "L", "\\theta": "1", "\\phi": "1", "r_s": "L",
    },
    # Ehlers and Kundt's A- and B-metrics, with their constant b a length. A time that stands
    # beside an angle is itself a number: tau in every B-metric is a boost parameter or the time
    # of a de Sitter or anti-de Sitter space of two dimensions with its radius factored out. Beyond the
    # horizon of AII the time is the length sigma. Kruskal's U and V are numbers, and the inertial
    # charts take their time T as a length and name sigma and r, the intervals from the tachyon's line.
    ("ab_metrics", "a2_static"): {"t": "T", "r": "L", "\\chi": "1", "\\phi": "1", "b": "L"},
    ("ab_metrics", "a2_cone"): {"\\sigma": "L", "z": "L", "\\chi": "1", "\\phi": "1", "b": "L"},
    ("ab_metrics", "a2_kruskal"): {"U": "1", "V": "1", "\\chi": "1", "\\phi": "1", "b": "L", "r": "L"},
    ("ab_metrics", "a2_cartesian"): {"T": "L", "X": "L", "Y": "L", "Z": "L", "b": "L", "\\sigma": "L"},
    ("ab_metrics", "b1_static"): {"\\tau": "1", "r": "L", "\\theta": "1", "z": "L", "b": "L"},
    ("ab_metrics", "b1_cone"): {"\\tau": "1", "r": "L", "\\phi": "1", "z": "L", "b": "L"},
    ("ab_metrics", "b1_neck"): {"\\tau": "1", "\\rho": "1", "\\phi": "1", "z": "L", "b": "L"},
    ("ab_metrics", "b1_cartesian"): {"T": "L", "X": "L", "Y": "L", "Z": "L", "b": "L", "r": "L"},
    ("ab_metrics", "a3"): {"t": "T", "r": "L", "\\chi": "1", "\\phi": "1", "b": "L"},
    ("ab_metrics", "b2_static"): {"\\tau": "1", "r": "L", "\\chi": "1", "z": "L", "b": "L"},
    ("ab_metrics", "b2_neck"): {"\\tau": "1", "\\rho": "1", "\\chi": "1", "z": "L", "b": "L"},
    ("ab_metrics", "b3"): {"\\tau": "1", "r": "L", "x": "1", "z": "L", "b": "L"},
    # Kruskal's null coordinates U and V are numbers, and the factor 4 r_s^3 e^(-r/r_s)/r carries
    # the area.
    ("white_hole", "exterior_kruskal"): {
        "U": "1", "V": "1", "\\theta": "1", "\\phi": "1", "r_s": "L", "r": "L",
    },
    # The RP3 geon: Kruskal's T and X are numbers, and the factor 4 r_s^3 e^(-r/r_s)/r carries the
    # area; Schwarzschild's chart of the one exterior; and the isotropic radius rho, a length.
    ("rp3_geon", "kruskal"): {
        "T": "1", "X": "1", "\\theta": "1", "\\phi": "1", "r_s": "L", "r": "L",
    },
    ("rp3_geon", "schwarzschild"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "r_s": "L",
    },
    ("rp3_geon", "isotropic"): {
        "t": "T", "\\rho": "L", "\\theta": "1", "\\phi": "1", "r_s": "L",
    },
    # The T-models of Datt and Ruban are the Kantowski-Sachs cylinder with its scale factor along
    # the axis free to differ from shell to shell: r is a length along the axis, a a pure number and
    # b the radius of the spheres. In Ruban's chart the cycloid's parameter eta is a pure number and
    # mu, the rest mass per unit of r as a length, is one too; the areal chart's time T is the
    # radius of the spheres, a length; and the de Sitter chart's ell is de Sitter's radius.
    ("datt_ruban_t_models", "comoving"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "a": "1", "b": "L",
    },
    ("datt_ruban_t_models", "ruban"): {
        "\\eta": "1", "r": "L", "\\theta": "1", "\\phi": "1", "r_s": "L", "\\epsilon": "1", "\\mu": "1", "a": "1",
    },
    ("datt_ruban_t_models", "areal"): {
        "T": "L", "r": "L", "\\theta": "1", "\\phi": "1", "r_s": "L", "\\epsilon": "1", "\\mu": "1", "a": "1",
    },
    ("datt_ruban_t_models", "de_sitter"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "\\ell": "L", "\\epsilon": "1", "\\mu": "1", "a": "1",
    },
    ("datt_ruban_t_models", "exterior_kruskal"): {
        "U": "1", "V": "1", "\\theta": "1", "\\phi": "1", "r_s": "L", "r": "L",
    },
    # Novikov's chart is Tolman's with no energy function: the shell label r is a length beside the
    # areal radius, F = 2GM(r)/c^2 is the Schwarzschild radius of the mass inside the shell, and b is
    # c times the moment the shell leaves the singularity, so q = ct - b is a length as well.
    ("white_hole", "novikov_comoving"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "F": "L", "b": "L", "q": "L", "R": "L",
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
    # Majumdar and Papapetrou's potential U multiplies proper length and divides proper time,
    # so it is a pure number, and so is 1 + m/r, whose mass parameter m = GM/c^2 is a length.
    ("majumdar_papapetrou", "cartesian"): {
        "t": "T", "x": "L", "y": "L", "z": "L", "U": "1",
    },
    ("majumdar_papapetrou", "cylindrical"): {
        "t": "T", "\\rho": "L", "\\phi": "1", "z": "L", "U": "1",
    },
    ("majumdar_papapetrou", "isotropic"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "m": "L",
    },
    # The draining bathtub. Here c is the speed of sound, and the chart is x^0 = ct with that c. The
    # strengths A and B of the radial flow and of the swirl are a speed times a length.
    ("draining_bathtub", "laboratory"): {"t": "T", "r": "L", "\\theta": "1", "A": "L**2/T", "B": "L**2/T"},
    ("draining_bathtub", "kerr_like"): {"T": "T", "r": "L", "\\phi": "1", "A": "L**2/T", "B": "L**2/T"},
    ("draining_bathtub", "vortex_filament"): {
        "t": "T", "r": "L", "\\theta": "1", "z": "L", "A": "L**2/T", "B": "L**2/T",
    },
    # Unruh's acoustic black hole, the spherical flow at constant density, with c the speed of sound and the
    # chart x^0 = ct with that c. The horizon radius r_0 is a length.
    ("unruh_acoustic_hole", "laboratory"): {"t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "r_0": "L"},
    ("unruh_acoustic_hole", "unruh"): {"\\tau": "T", "r": "L", "\\theta": "1", "\\phi": "1", "r_0": "L"},
    # Bonnor's stars of charged dust. The potential U is a pure number, as Majumdar and Papapetrou's
    # is; the mass parameter m = GM/c^2, the radius r_0, the focal radius a and Lemos and Weinberg's
    # core length b are lengths, and the spheroidal u and its surface u_0 are pure numbers.
    ("bonnor_charged_dust", "harmonic"): {
        "t": "T", "x": "L", "y": "L", "z": "L", "U": "1",
    },
    ("bonnor_charged_dust", "sphere_1965"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "m": "L", "r_0": "L", "U": "1",
    },
    ("bonnor_charged_dust", "sphere_1975"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "m": "L", "r_0": "L", "U": "1",
    },
    ("bonnor_charged_dust", "exterior"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "m": "L", "r_0": "L",
    },
    ("bonnor_charged_dust", "exterior_areal"): {
        "t": "T", "R": "L", "\\theta": "1", "\\phi": "1", "m": "L", "r_0": "L",
    },
    ("bonnor_charged_dust", "spheroid_interior"): {
        "t": "T", "u": "1", "\\theta": "1", "\\phi": "1", "m": "L", "a": "L", "u_0": "1", "U": "1",
    },
    ("bonnor_charged_dust", "spheroid_exterior"): {
        "t": "T", "u": "1", "\\theta": "1", "\\phi": "1", "m": "L", "a": "L", "u_0": "1", "U": "1",
    },
    ("bonnor_charged_dust", "quasi_black_hole"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "m": "L", "b": "L", "U": "1",
    },
    # The black holes of string theory with three and four charges. The radius r_0 of the horizon
    # and the charge radii r_i are lengths, r_0^2 sinh^2 alpha_i = r_i^2 in five dimensions and
    # r_0 sinh^2 alpha_i = r_i in four, and f and the harmonic functions H_i are pure numbers.
    ("string_bh_three_four_charges", "five_charges"): {
        "t": "T", "r": "L", "\\psi": "1", "\\theta": "1", "\\phi": "1", "r_0": "L", "r_1": "L", "r_2": "L",
        "r_3": "L", "f": "1", "H_1": "1", "H_2": "1", "H_3": "1",
    },
    ("string_bh_three_four_charges", "five_extreme"): {
        "t": "T", "r": "L", "\\psi": "1", "\\theta": "1", "\\phi": "1", "r_1": "L", "r_2": "L", "r_3": "L",
        "H_1": "1", "H_2": "1", "H_3": "1",
    },
    ("string_bh_three_four_charges", "five_areal"): {
        "t": "T", "\\rho": "L", "\\psi": "1", "\\theta": "1", "\\phi": "1", "r_0": "L", "r_q": "L",
    },
    ("string_bh_three_four_charges", "four_charges"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "r_0": "L", "r_1": "L", "r_2": "L", "r_3": "L",
        "r_4": "L", "f": "1", "H_1": "1", "H_2": "1", "H_3": "1", "H_4": "1",
    },
    ("string_bh_three_four_charges", "four_extreme"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "r_1": "L", "r_2": "L", "r_3": "L", "r_4": "L",
        "H_1": "1", "H_2": "1", "H_3": "1", "H_4": "1",
    },
    # Israel, Wilson and Perjes's W = |U| is a pure number, as Majumdar and Papapetrou's U is, and
    # omega stands beside c dt, so it is a length; the spin a and the NUT parameter l are lengths too.
    ("israel_wilson_perjes", "cylindrical"): {
        "t": "T", "\\rho": "L", "\\phi": "1", "z": "L", "W": "1", "\\omega": "L",
    },
    ("israel_wilson_perjes", "spheroidal"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "m": "L", "a": "L",
        "\\Sigma_0": "L**2", "\\Sigma": "L**2",
    },
    ("israel_wilson_perjes", "spherical"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "m": "L", "l": "L",
    },
    # The rocket's mass function is the length GM/c^2 and its acceleration alpha an inverse
    # length; in the Robinson-Trautman chart the four-velocity over c and p are pure numbers.
    ("photon_rocket", "rectilinear"): {
        "u": "T", "r": "L", "\\theta": "1", "\\phi": "1", "m": "L", "\\alpha": "1/L",
    },
    ("photon_rocket", "robinson_trautman"): {
        "u": "T", "r": "L", "\\theta": "1", "\\phi": "1", "m": "L", "U_1": "1", "U_2": "1", "U_3": "1", "p": "1",
    },
    # Robinson and Trautman's r is an affine length along the rays and the fronts' coordinates
    # are angles, so P, f and H are pure numbers.
    ("robinson_trautman", "stereographic"): {
        "u": "T", "r": "L", "x": "1", "y": "1", "P": "1", "H": "1",
    },
    ("robinson_trautman", "axisymmetric"): {
        "u": "T", "r": "L", "\\theta": "1", "\\phi": "1", "f": "1", "H": "1",
    },
    # Newman and Tamburino's r is an affine length along the rays and the fronts' x and y are pure
    # numbers in the spherical class, so A, and a = A s, are lengths and b is pure; in the cylindrical
    # limit y carries the length, and ell is the length that makes r^2 x^4 pure inside the logarithm.
    ("newman_tamburino", "spherical"): {
        "u": "T", "r": "L", "x": "1", "y": "1", "b": "1", "A_0": "L", "A": "L", "s": "1", "a": "L", "L": "1",
    },
    ("newman_tamburino", "cylindrical"): {"u": "T", "r": "L", "x": "1", "y": "L", "\\ell": "L"},
    # Penney's a and b are lengths, and gamma and the ratio power E are pure numbers.
    ("penney_scalar", "spherical"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "a": "L", "b": "L", "\\gamma": "1", "E": "1",
    },
    # Bronnikov's harmonic u is an inverse length, so k, h, the charge q and K = k coth(ku) and
    # H = h coth(h(u + u_1)) are lengths and u_1 an inverse length.
    ("penney_scalar", "harmonic"): {
        "t": "T", "u": "1/L", "\\theta": "1", "\\phi": "1", "k": "L", "h": "L", "u_1": "1/L", "q": "L", "K": "L",
        "H": "L",
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
    # Kottler's metric keeps the two lengths of its parents: Schwarzschild's r_s = 2GM/c^2 and
    # de Sitter's Lambda, a curvature, which is what leaves 1 - r_s/r - Lambda r^2/3 a pure
    # number. The Eddington-Finkelstein times u = ct - r_* and v = ct + r_* are lengths.
    ("schwarzschild_de_sitter", "static"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "r_s": "L", "\\Lambda": "1/L**2",
    },
    ("schwarzschild_de_sitter", "eddington_finkelstein_outgoing"): {
        "u": "L", "r": "L", "\\theta": "1", "\\phi": "1", "r_s": "L", "\\Lambda": "1/L**2",
    },
    ("schwarzschild_de_sitter", "eddington_finkelstein_ingoing"): {
        "v": "L", "r": "L", "\\theta": "1", "\\phi": "1", "r_s": "L", "\\Lambda": "1/L**2",
    },
    # Kiselev's black hole keeps Schwarzschild's length r_s beside the length r_q of the matter
    # around it and the state parameter w, a pure number, so h = (r_q/r)^{3w + 1}, the name the
    # charts that keep w define, is one too. The Eddington-Finkelstein times are lengths, and
    # Kiselev's hyperbolic and conformally flat coordinates are all pure numbers.
    ("kiselev", "static"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "r_s": "L", "r_q": "L", "w": "1", "h": "1",
    },
    ("kiselev", "eddington_finkelstein_ingoing"): {
        "v": "L", "r": "L", "\\theta": "1", "\\phi": "1", "r_s": "L", "r_q": "L", "w": "1", "h": "1",
    },
    ("kiselev", "eddington_finkelstein_outgoing"): {
        "u": "L", "r": "L", "\\theta": "1", "\\phi": "1", "r_s": "L", "r_q": "L", "w": "1", "h": "1",
    },
    ("kiselev", "linear"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "r_s": "L", "r_q": "L",
    },
    ("kiselev", "hyperbolic"): {
        "\\eta": "1", "\\chi": "1", "\\theta": "1", "\\phi": "1", "r_q": "L",
    },
    ("kiselev", "conformally_flat"): {
        "\\tau": "1", "\\rho": "1", "\\theta": "1", "\\phi": "1", "r_q": "L",
    },
    # Israel's disc source of Kerr's metric: Kerr's metric with the mass as the length m = GM/c^2,
    # the spin per unit mass a a length, Sigma and Delta areas, and the Kerr-Schild time tau a time.
    ("kerr_disc_source", "boyer_lindquist"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "m": "L", "a": "L", "\\Sigma": "L**2", "\\Delta": "L**2",
    },
    ("kerr_disc_source", "kerr_schild"): {
        "\\tau": "T", "r": "L", "\\theta": "1", "\\psi": "1", "m": "L", "a": "L", "\\Sigma": "L**2", "\\Delta": "L**2",
    },
    # Carter's rotating black hole with a cosmological constant of either sign. The mass enters as
    # the length r_s = 2GM/c^2, the rotation parameter a is a length as Kerr's is, and Lambda is a
    # curvature. Every chart defines four names: Xi and Delta_theta are pure numbers, rho is a
    # length and Delta_r an area. The Kerr times v and u are lengths.
    ("kerr_de_sitter", "boyer_lindquist"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "r_s": "L", "a": "L", "\\Lambda": "1/L**2",
        "\\Xi": "1", "\\rho": "L", "\\Delta_r": "L**2", "\\Delta_\\theta": "1",
    },
    ("kerr_de_sitter", "nonrotating"): {
        "t": "T", "r": "L", "\\theta": "1", "\\Phi": "1", "r_s": "L", "a": "L", "\\Lambda": "1/L**2",
        "\\Xi": "1", "\\rho": "L", "\\Delta_r": "L**2", "\\Delta_\\theta": "1",
    },
    ("kerr_de_sitter", "kerr_ingoing"): {
        "v": "L", "r": "L", "\\theta": "1", "\\tilde\\phi": "1", "r_s": "L", "a": "L", "\\Lambda": "1/L**2",
        "\\Xi": "1", "\\rho": "L", "\\Delta_r": "L**2", "\\Delta_\\theta": "1",
    },
    ("kerr_de_sitter", "kerr_outgoing"): {
        "u": "L", "r": "L", "\\theta": "1", "\\tilde\\phi": "1", "r_s": "L", "a": "L", "\\Lambda": "1/L**2",
        "\\Xi": "1", "\\rho": "L", "\\Delta_r": "L**2", "\\Delta_\\theta": "1",
    },
    ("kerr_de_sitter", "kerr_schild"): {
        "\\tau": "T", "r": "L", "\\theta": "1", "\\psi": "1", "r_s": "L", "a": "L", "\\Lambda": "1/L**2",
        "\\Xi": "1", "\\rho": "L", "\\Delta_r": "L**2", "\\Delta_\\theta": "1",
    },
    # Kerr-anti-de Sitter spun up to a = ell. The anti-de Sitter radius ell and the mass parameter m
    # are lengths, Sigma, Delta, Q and P areas, and the azimuth psi a pure number. In the
    # Carter-Plebanski chart q and p are lengths, tau a time, and sigma = -psi/ell an inverse length.
    ("super_entropic", "hennigar_kubiznak_mann"): {
        "t": "T", "r": "L", "\\theta": "1", "\\psi": "1", "\\ell": "L", "m": "L",
        "\\Sigma": "L**2", "\\Delta": "L**2",
    },
    ("super_entropic", "carter_plebanski"): {
        "\\tau": "T", "q": "L", "p": "L", "\\sigma": "1/L", "\\ell": "L", "m": "L",
        "Q": "L**2", "P": "L**2",
    },
    # Kerr's black hole with a NUT parameter. The mass enters as the length m = GM/c^2, as on the
    # Taub-NUT page, and a and l are lengths. Sigma and Delta are areas and chi a length. The Kerr
    # times v and u are lengths. In Plebanski and Demianski's chart q and p are lengths, tau is a
    # length that stands as the time, and sigma is an inverse length, so that p^2 sigma is a length.
    ("kerr_taub_nut", "boyer_lindquist"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "m": "L", "a": "L", "l": "L",
        "\\Sigma": "L**2", "\\Delta": "L**2", "\\chi": "L",
    },
    ("kerr_taub_nut", "one_string"): {
        "t_N": "T", "r": "L", "\\theta": "1", "\\phi": "1", "m": "L", "a": "L", "l": "L",
        "\\Sigma": "L**2", "\\Delta": "L**2", "\\chi": "L",
    },
    ("kerr_taub_nut", "kerr_ingoing"): {
        "v": "L", "r": "L", "\\theta": "1", "\\tilde\\phi": "1", "m": "L", "a": "L", "l": "L",
        "\\Sigma": "L**2", "\\Delta": "L**2", "\\chi": "L",
    },
    ("kerr_taub_nut", "kerr_outgoing"): {
        "u": "L", "r": "L", "\\theta": "1", "\\tilde\\phi": "1", "m": "L", "a": "L", "l": "L",
        "\\Sigma": "L**2", "\\Delta": "L**2", "\\chi": "L",
    },
    ("kerr_taub_nut", "plebanski"): {
        "\\tau": "L", "q": "L", "p": "L", "\\sigma": "1/L", "m": "L", "a": "L", "l": "L",
        "Q": "L**2", "P": "L**2",
    },
    # Brill's charged Taub-NUT. The mass enters as the length m = GM/c^2, as on the Taub-NUT page,
    # and the NUT parameter l and the charge radius r_q are lengths. Sigma and Delta are areas. The
    # times v and u are lengths. In the chart of Brill's universe tau is a length that stands as
    # the time, psi is an Euler angle, and U is an area.
    ("brill_charged_taub_nut", "spherical"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "m": "L", "l": "L", "r_q": "L",
        "\\Sigma": "L**2", "\\Delta": "L**2",
    },
    ("brill_charged_taub_nut", "one_string"): {
        "t_N": "T", "r": "L", "\\theta": "1", "\\phi": "1", "m": "L", "l": "L", "r_q": "L",
        "\\Sigma": "L**2", "\\Delta": "L**2",
    },
    ("brill_charged_taub_nut", "eddington_finkelstein_ingoing"): {
        "v": "L", "r": "L", "\\theta": "1", "\\phi": "1", "m": "L", "l": "L", "r_q": "L",
        "\\Sigma": "L**2", "\\Delta": "L**2",
    },
    ("brill_charged_taub_nut", "eddington_finkelstein_outgoing"): {
        "u": "L", "r": "L", "\\theta": "1", "\\phi": "1", "m": "L", "l": "L", "r_q": "L",
        "\\Sigma": "L**2", "\\Delta": "L**2",
    },
    ("brill_charged_taub_nut", "taub"): {
        "\\tau": "L", "\\psi": "1", "\\theta": "1", "\\phi": "1", "m": "L", "l": "L", "r_q": "L",
        "\\Sigma": "L**2", "U": "L**2",
    },
    # The Kaluza-Klein monopole has one length, m, with G nowhere in the line element. Gross
    # and Perry's fifth coordinate x_5 is a length of period 16 pi m, and the Hopf angle psi is a pure
    # number.
    ("kaluza_klein_monopole", "gross_perry"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "x_5": "L", "m": "L",
    },
    ("kaluza_klein_monopole", "hopf"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "\\psi": "1", "m": "L",
    },
    ("kaluza_klein_monopole", "taub_nut"): {
        "t": "T", "\\rho": "L", "\\theta": "1", "\\phi": "1", "\\psi": "1", "m": "L",
    },
    # The Eguchi-Hanson space is Riemannian: every coordinate is spatial and none is a time. Its one
    # parameter a is a length, and the Euler angles are pure numbers.
    ("eguchi_hanson", "eguchi_hanson"): {"r": "L", "\\theta": "1", "\\phi": "1", "\\psi": "1", "a": "L"},
    ("eguchi_hanson", "kahler"): {"\\rho": "L", "\\theta": "1", "\\phi": "1", "\\psi": "1", "a": "L", "W": "L^2"},
    # The two-centre chart's fibre coordinate tau is a length of period pi a/2; V is a pure number and
    # omega, the potential that stands beside d tau, a length.
    ("eguchi_hanson", "two_centre"): {
        "\\rho": "L", "z": "L", "\\psi": "1", "\\tau": "L", "a": "L", "R_1": "L", "R_2": "L", "V": "1", "\\omega": "L",
    },
    # The gravitational instantons of 1977 and 1978 are Riemannian, as the Eguchi-Hanson space is:
    # every coordinate is spatial and none is a time. The imaginary time tau is a length, r_s and
    # the nut parameter n are lengths, Lambda an inverse area, and Page's nu a pure number. The
    # regular chart's x is a length and its r a name for the areal radius. The multi-centre chart's
    # V is a pure number and omega, which multiplies d phi beside d tau, a length.
    ("gravitational_instantons", "schwarzschild"): {"\\tau": "L", "r": "L", "\\theta": "1", "\\phi": "1", "r_s": "L"},
    ("gravitational_instantons", "regular"): {
        "x": "L", "\\tau": "L", "\\theta": "1", "\\phi": "1", "r_s": "L", "r": "L",
    },
    ("gravitational_instantons", "taub_nut"): {"r": "L", "\\theta": "1", "\\phi": "1", "\\tau": "L", "n": "L"},
    ("gravitational_instantons", "taub_bolt"): {"r": "L", "\\theta": "1", "\\phi": "1", "\\tau": "L", "n": "L"},
    ("gravitational_instantons", "multi_centre"): {
        "\\rho": "L", "z": "L", "\\phi": "1", "\\tau": "L", "V": "1", "\\omega": "L",
    },
    ("gravitational_instantons", "cp2"): {"r": "L", "\\theta": "1", "\\phi": "1", "\\psi": "1", "\\Lambda": "1/L**2"},
    ("gravitational_instantons", "cp2_distance"): {"\\chi": "L", "\\theta": "1", "\\phi": "1", "\\psi": "1", "L": "L"},
    ("gravitational_instantons", "page"): {
        "\\chi": "1", "\\theta": "1", "\\phi": "1", "\\psi": "1", "\\Lambda": "1/L**2", "\\nu": "1",
        "P": "1", "Q": "1", "N": "1",
    },
    # Misner's and Brill and Lindquist's initial data are one moment, a Riemannian space of three
    # dimensions: every coordinate is spatial. The conformal factors psi, chi and Psi are pure numbers,
    # the bispherical coordinates are pure numbers, and their scale a and the one hole's r_s are lengths.
    ("misner_brill_lindquist", "cartesian"): {"x": "L", "y": "L", "z": "L", "\\psi": "1"},
    ("misner_brill_lindquist", "cylindrical"): {"\\rho": "L", "\\phi": "1", "z": "L", "\\psi": "1"},
    ("misner_brill_lindquist", "bispherical"): {"\\mu": "1", "\\eta": "1", "\\phi": "1", "a": "L", "\\Psi": "1"},
    ("misner_brill_lindquist", "isotropic"): {"r": "L", "\\theta": "1", "\\phi": "1", "r_s": "L"},
    ("misner_brill_lindquist", "charged"): {"x": "L", "y": "L", "z": "L", "\\chi": "1", "\\psi": "1"},
    ("bowen_york", "cartesian"): {"x": "L", "y": "L", "z": "L", "\\psi": "1"},
    ("bowen_york", "spherical"): {"r": "L", "\\theta": "1", "\\phi": "1", "\\psi": "1"},
    ("bowen_york", "cylindrical"): {"\\rho": "L", "z": "L", "\\phi": "1", "\\psi": "1"},
    # Brill's waves, one moment of a spacetime too: the conformal factor psi and the free function q
    # of the base metric, which stands in an exponential, are pure numbers.
    ("brill_waves", "cylindrical"): {"\\rho": "L", "z": "L", "\\phi": "1", "\\psi": "1", "q": "1"},
    ("brill_waves", "spherical"): {"r": "L", "\\theta": "1", "\\phi": "1", "\\psi": "1", "q": "1"},
    ("brill_waves", "three_dimensional"): {"\\rho": "L", "z": "L", "\\phi": "1", "\\psi": "1", "q": "1"},
    # Lewis's stationary cylinders. His coordinates are pure numbers, so each chart measures r in a
    # length, ell or van Stockum's radius R, and every power whose exponent holds a parameter is a
    # power of that ratio: the names u, h and H are pure numbers, and so are Lewis's constants n, a
    # and q, the last his c, while his b is a length beside c dt. f and F are numbers, k and M
    # lengths and l and L areas, as -f c^2dt^2 + 2k c dt dphi + l dphi^2 wants. The Komar mass per
    # unit length sigma is a number and the angular momentum per unit length j a length, each with
    # its power of G and c, and van Stockum's w, the dust's angular velocity over c, is an inverse length.
    ("lewis", "lewis"): {
        "t": "T", "r": "L", "\\phi": "1", "z": "L", "n": "1", "a": "1", "b": "L", "q": "1", "\\ell": "L",
        "u": "1", "h": "1", "f": "1", "k": "L", "l": "L**2",
    },
    ("lewis", "canonical"): {
        "t": "T", "r": "L", "\\phi": "1", "z": "L", "\\sigma": "1", "j": "L", "\\alpha": "1", "\\ell": "L",
        "u": "1", "h": "1",
    },
    ("lewis", "lewis_class"): {
        "t": "T", "r": "L", "\\phi": "1", "z": "L", "m": "1", "a_1": "1", "b_1": "1", "a_2": "1", "\\ell": "L",
        "b_2": "1", "\\psi": "1", "h": "1", "f": "1", "k": "L", "l": "L**2",
    },
    ("lewis", "stockum_light"): {
        "t": "T", "r": "L", "\\phi": "1", "z": "L", "n": "1", "R": "L", "\\ell": "L", "w": "1/L",
        "u": "1", "H": "1", "F": "1", "M": "L", "L": "L**2",
    },
    ("lewis", "stockum_critical"): {
        "t": "T", "r": "L", "\\phi": "1", "z": "L", "R": "L", "\\ell": "L", "H": "1", "F": "1", "M": "L",
        "L": "L**2",
    },
    ("lewis", "stockum_heavy"): {
        "t": "T", "r": "L", "\\phi": "1", "z": "L", "m": "1", "R": "L", "\\ell": "L", "w": "1/L",
        "\\psi": "1", "H": "1", "F": "1", "M": "L", "L": "L**2",
    },
    # Petrov's homogeneous vacuum. In his own chart every coordinate but the time is a length, phi
    # among them, which is no angle there, and ell is the one length of the solution. Outside the
    # dust cylinder phi is an angle, and F, M and L are a number, a length and an area, as in the
    # heavy cylinder of the Lewis family.
    ("petrov_homogeneous", "petrov"): {
        "t": "T", "r": "L", "\\phi": "L", "z": "L", "\\ell": "L", "\\psi": "1",
    },
    ("petrov_homogeneous", "cylinder"): {
        "t": "T", "r": "L", "\\phi": "1", "z": "L", "R": "L", "\\ell": "L", "\\psi": "1", "F": "1", "M": "L",
        "L": "L**2",
    },
    # The Kaluza-Klein black holes keep three lengths: the horizon radius r_s and the charge
    # lengths q and p. The fifth coordinate y is a length, and so are the advanced time v = ct + ...
    # and its companion w along the circle.
    ("kaluza_klein_black_hole", "electric"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "y": "L", "r_s": "L", "q": "L",
    },
    ("kaluza_klein_black_hole", "eddington_finkelstein_ingoing"): {
        "v": "L", "r": "L", "\\theta": "1", "\\phi": "1", "w": "L", "r_s": "L", "q": "L",
    },
    ("kaluza_klein_black_hole", "magnetic"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "y": "L", "r_s": "L", "p": "L",
    },
    ("kaluza_klein_black_hole", "dyonic"): {
        "t": "T", "\\rho": "L", "\\theta": "1", "\\phi": "1", "y": "L", "r_s": "L", "p": "L",
    },
    ("kaluza_klein_black_hole", "einstein"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "r_s": "L", "q": "L",
    },
    ("kaluza_klein_black_hole", "einstein_eddington_finkelstein"): {
        "v": "L", "r": "L", "\\theta": "1", "\\phi": "1", "r_s": "L", "q": "L",
    },
    # Lifshitz spacetime keeps one length, L, and one pure number, the dynamical exponent z. Every
    # radial coordinate is a length, and so is the advanced time v = ct - w of the two null charts.
    # Wahlquist's rotating fluid. In his own chart xi and eta are pure numbers and r_0 carries the
    # length; h_1, h_2 and f are numbers, and the Weyl scalars W_1 and W_2 are curvatures. In Mars's
    # charts y and z are lengths, sigma and varphi inverse lengths, Q_0, v_1, v_2, U and V areas,
    # mu_0 a curvature, beta an inverse length, and a_1, a_2, m and n lengths.
    ("wahlquist", "wahlquist"): {
        "t": "T", "\\xi": "1", "\\eta": "1", "\\phi": "1", "r_0": "L", "k": "1", "b": "1", "\\eta_0": "1",
        "\\gamma": "1", "h_1": "1", "h_2": "1", "W_1": "1/L**2", "W_2": "1/L**2",
    },
    ("wahlquist", "mars"): {
        "\\tau": "T", "y": "L", "z": "L", "\\sigma": "1/L", "Q_0": "L**2", "\\nu_0": "1", "\\mu_0": "1/L**2",
        "a_1": "L", "a_2": "L", "\\beta": "1/L", "v_1": "L**2", "v_2": "L**2", "n": "L", "m": "L", "U": "L**2",
        "V": "L**2", "W_1": "1/L**2", "W_2": "1/L**2",
    },
    ("wahlquist", "mars_ingoing"): {
        "v": "L", "y": "L", "z": "L", "\\phi": "1/L", "Q_0": "L**2", "\\nu_0": "1", "\\mu_0": "1/L**2",
        "a_1": "L", "a_2": "L", "\\beta": "1/L", "v_1": "L**2", "v_2": "L**2", "n": "L", "m": "L", "U": "L**2",
        "V": "L**2", "W_1": "1/L**2", "W_2": "1/L**2",
    },
    ("wahlquist", "whittaker"): {"t": "T", "X": "1", "\\theta": "1", "\\phi": "1", "R_0": "L", "b": "1", "F": "1"},
    ("lifshitz_spacetime", "kachru_liu_mulligan"): {"t": "T", "x": "L", "y": "L", "r": "L", "z": "1", "L": "L"},
    ("lifshitz_spacetime", "poincare"): {"t": "T", "x": "L", "y": "L", "u": "L", "z": "1", "L": "L"},
    ("lifshitz_spacetime", "proper_distance"): {"t": "T", "x": "L", "y": "L", "\\rho": "L", "z": "1", "L": "L"},
    ("lifshitz_spacetime", "tortoise"): {"t": "T", "x": "L", "y": "L", "w": "L", "z": "1", "L": "L"},
    ("lifshitz_spacetime", "eddington_finkelstein"): {"v": "L", "x": "L", "y": "L", "r": "L", "z": "1", "L": "L"},
    ("lifshitz_spacetime", "affine"): {"v": "L", "x": "L", "y": "L", "s": "L", "z": "1", "L": "L"},
    # Van Den Broeck's drive adds one more pure number to Alcubierre's two: B, the factor every
    # length of a slice is multiplied by. The comoving chart is the inside of the bubble, f = 1.
    ("van_den_broeck", "cartesian"): {
        "t": "T", "x": "L", "y": "L", "z": "L", "v_s": "1", "f": "1", "B": "1",
    },
    ("van_den_broeck", "pocket"): {"t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "B": "1", "R": "L"},
    # Krasnikov's chart of the pocket: l is the proper distance along a radius and r(l) the areal radius.
    ("van_den_broeck", "proper_radial"): {"t": "T", "l": "L", "\\theta": "1", "\\phi": "1", "r": "L", "l_0": "L"},
    # Tangherlini's black hole quotes its mass as the horizon radius r_h, a length in every
    # dimension, which leaves 1 - (r_h/r)^(D-3) a pure number; the charts are D = 5 and D = 6,
    # and the Eddington-Finkelstein times u = ct - r_* and v = ct + r_* are lengths.
    ("tangherlini", "spherical"): {
        "t": "T", "r": "L", "\\psi": "1", "\\theta": "1", "\\phi": "1", "r_h": "L",
    },
    ("tangherlini", "eddington_finkelstein_outgoing"): {
        "u": "L", "r": "L", "\\psi": "1", "\\theta": "1", "\\phi": "1", "r_h": "L",
    },
    ("tangherlini", "eddington_finkelstein_ingoing"): {
        "v": "L", "r": "L", "\\psi": "1", "\\theta": "1", "\\phi": "1", "r_h": "L",
    },
    ("tangherlini", "spherical_six"): {
        "t": "T", "r": "L", "\\chi": "1", "\\psi": "1", "\\theta": "1", "\\phi": "1", "r_h": "L",
    },
    # Boulware and Deser's black hole quotes its mass as Tangherlini's radius r_0 and the coupling
    # of the Gauss-Bonnet term as the length l, l^2 = 2 alpha; the radical W = sqrt(r^4 + 4 l^2 r_0^2)
    # is an area, and the Eddington-Finkelstein times u = ct - r_* and v = ct + r_* are lengths.
    ("boulware_deser", "spherical"): {
        "t": "T", "r": "L", "\\psi": "1", "\\theta": "1", "\\phi": "1", "r_0": "L", "\\ell": "L", "W": "L**2",
    },
    ("boulware_deser", "eddington_finkelstein_outgoing"): {
        "u": "L", "r": "L", "\\psi": "1", "\\theta": "1", "\\phi": "1", "r_0": "L", "\\ell": "L", "W": "L**2",
    },
    ("boulware_deser", "eddington_finkelstein_ingoing"): {
        "v": "L", "r": "L", "\\psi": "1", "\\theta": "1", "\\phi": "1", "r_0": "L", "\\ell": "L", "W": "L**2",
    },
    ("boulware_deser", "spherical_plus"): {
        "t": "T", "r": "L", "\\psi": "1", "\\theta": "1", "\\phi": "1", "r_0": "L", "\\ell": "L", "W": "L**2",
    },
    # Randall and Sundrum's wall has one inverse length, k, the curvature of the anti-de Sitter space
    # on each side; x_1, x_2 and x_3 run along the wall, and the fifth coordinate is a length in every
    # chart but the one between two walls, where it is the angle phi and r_c carries the length.
    ("randall_sundrum", "proper_distance"): {"t": "T", "x_1": "L", "x_2": "L", "x_3": "L", "y": "L", "k": "1/L"},
    ("randall_sundrum", "conformal"): {"t": "T", "x_1": "L", "x_2": "L", "x_3": "L", "w": "L", "k": "1/L"},
    ("randall_sundrum", "poincare"): {"t": "T", "x_1": "L", "x_2": "L", "x_3": "L", "z": "L", "k": "1/L"},
    ("randall_sundrum", "two_walls"): {
        "t": "T", "x_1": "L", "x_2": "L", "x_3": "L", "\\phi": "1", "k": "1/L", "r_c": "L",
    },
    # Bardeen's regular black hole keeps Schwarzschild's r_s; g, the monopole's charge, is a length.
    # The Eddington-Finkelstein times u = ct - r_* and v = ct + r_* are lengths.
    ("bardeen", "static"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "r_s": "L", "g": "L",
    },
    ("bardeen", "eddington_finkelstein_outgoing"): {
        "u": "L", "r": "L", "\\theta": "1", "\\phi": "1", "r_s": "L", "g": "L",
    },
    ("bardeen", "eddington_finkelstein_ingoing"): {
        "v": "L", "r": "L", "\\theta": "1", "\\phi": "1", "r_s": "L", "g": "L",
    },
    # Born and Infeld's point charge keeps Schwarzschild's r_s and Reissner and Nordstrom's r_q; r_0,
    # the radius at which the charge's Coulomb field would be Born and Infeld's greatest field, is a
    # length, the radical W an area and the mass function m = G M(r)/c^2 a length. The
    # Eddington-Finkelstein times u = ct - r_* and v = ct + r_* are lengths.
    ("born_infeld_charge", "static"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "r_s": "L", "r_q": "L", "r_0": "L", "W": "L^2", "m": "L",
    },
    ("born_infeld_charge", "eddington_finkelstein_outgoing"): {
        "u": "L", "r": "L", "\\theta": "1", "\\phi": "1", "r_s": "L", "r_q": "L", "r_0": "L", "W": "L^2", "m": "L",
    },
    ("born_infeld_charge", "eddington_finkelstein_ingoing"): {
        "v": "L", "r": "L", "\\theta": "1", "\\phi": "1", "r_s": "L", "r_q": "L", "r_0": "L", "W": "L^2", "m": "L",
    },
    # The black string is Schwarzschild's black hole with the length z along the string added, and
    # in six dimensions Tangherlini's of five with it; the advanced time v = ct + r_* is a length and
    # the Kerr-Schild time T = (v - r)/c a time.
    ("black_string", "static"): {"t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "z": "L", "r_s": "L"},
    ("black_string", "eddington_finkelstein_ingoing"): {
        "v": "L", "r": "L", "\\theta": "1", "\\phi": "1", "z": "L", "r_s": "L",
    },
    ("black_string", "kerr_schild"): {"T": "T", "r": "L", "\\theta": "1", "\\phi": "1", "z": "L", "r_s": "L"},
    ("black_string", "static_six"): {
        "t": "T", "r": "L", "\\psi": "1", "\\theta": "1", "\\phi": "1", "z": "L", "r_h": "L",
    },
    # Myers and Perry's rotating black hole keeps the mass parameter mu, which is r_h^(D-3) at no
    # spin: an area in D = 5 and a volume in D = 6. Each spin parameter is a length, as Kerr's is.
    # The ingoing time v = ct + r_* is a length, and the equal-spin chart's rho^2 = r^2 + a^2.
    ("myers_perry", "boyer_lindquist"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "\\psi": "1", "\\mu": "L**2", "a": "L",
    },
    ("myers_perry", "ingoing_kerr"): {
        "v": "L", "r": "L", "\\theta": "1", "\\tilde\\phi": "1", "\\psi": "1", "\\mu": "L**2", "a": "L",
    },
    ("myers_perry", "two_spins"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "\\psi": "1", "\\mu": "L**2", "a": "L", "b": "L",
    },
    ("myers_perry", "equal_spins"): {
        "t": "T", "\\rho": "L", "\\theta": "1", "\\phi": "1", "\\psi": "1", "\\mu": "L**2", "a": "L",
    },
    ("myers_perry", "boyer_lindquist_six"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "\\chi": "1", "\\psi": "1", "\\mu": "L**3", "a": "L",
    },
    # Black Saturn's Weyl coordinates are in units of L^2, Elvang and Figueras's rho/L^2 and z-bar, so
    # rho and z are pure numbers and the one length L carries the dimensions: every soliton mu_i, each
    # of their polynomials and functions, their constants in units of L and the four free functions of
    # the chart are pure numbers. The polar chart's r is a length, and rho and z are names it defines.
    # The ring chart keeps its radius R, and nu, lambda and C are pure numbers.
    ("black_saturn", "weyl"): {
        "t": "T", "\\rho": "1", "z": "1", "\\psi": "1", "\\phi": "1", "L": "L", "\\kappa_1": "1", "\\kappa_2": "1",
        "\\kappa_3": "1", "\\beta": "1", "c_1": "1", "c_2": "1", "q": "1", "k": "1", "\\mu_1": "1", "\\mu_2": "1",
        "\\mu_3": "1", "\\mu_4": "1", "\\mu_5": "1", "R_1": "1", "R_2": "1", "M_0": "1", "M_1": "1", "M_2": "1",
        "M_3": "1", "M_4": "1", "F": "1", "H_x": "1", "H_y": "1", "G_y": "1", "P": "1", "\\omega_\\psi": "1",
        "V": "1", "\\Omega": "1", "W": "1", "\\nu": "1",
    },
    ("black_saturn", "polar"): {
        "t": "T", "r": "L", "\\theta": "1", "\\psi": "1", "\\phi": "1", "L": "L", "\\rho": "1", "z": "1",
        "\\kappa_1": "1", "\\kappa_2": "1", "\\kappa_3": "1", "\\beta": "1", "c_1": "1", "c_2": "1", "q": "1",
        "k": "1", "\\mu_1": "1", "\\mu_2": "1", "\\mu_3": "1", "\\mu_4": "1", "\\mu_5": "1", "R_1": "1", "R_2": "1",
        "M_0": "1", "M_1": "1", "M_2": "1", "M_3": "1", "M_4": "1", "F": "1", "H_x": "1", "H_y": "1", "G_y": "1",
        "P": "1", "\\omega_\\psi": "1", "V": "1", "\\Omega": "1", "W": "1", "\\nu": "1",
    },
    ("black_saturn", "ring"): {
        "t": "T", "x": "1", "y": "1", "\\psi": "1", "\\phi": "1", "R": "L", "\\nu": "1", "\\lambda": "1", "C": "1",
    },
    # Kastor and Traschen's holes: the potential V of the holes is a pure number, as Majumdar and
    # Papapetrou's U is, and so are U = H tau + V, the scale factor a = e^{Ht} and Omega = 1 + V/a,
    # each a name its chart defines; H is a frequency with 3H^2/c^2 = Lambda, as de Sitter's is.
    ("kastor_traschen", "cartesian"): {
        "\\tau": "T", "x": "L", "y": "L", "z": "L", "H": "1/T", "V": "1", "U": "1",
    },
    ("kastor_traschen", "cylindrical"): {
        "\\tau": "T", "\\rho": "L", "\\phi": "1", "z": "L", "H": "1/T", "V": "1", "U": "1",
    },
    ("kastor_traschen", "isotropic"): {
        "\\tau": "T", "r": "L", "\\theta": "1", "\\phi": "1", "H": "1/T", "m": "L", "U": "1",
    },
    ("kastor_traschen", "comoving"): {
        "t": "T", "x": "L", "y": "L", "z": "L", "H": "1/T", "a": "1", "V": "1", "\\Omega": "1",
    },
    # Som and Raychaudhuri's universe keeps c, and its Omega is the angular velocity of the dust,
    # so that Omega r^2/c is a length beside c dt.
    ("som_raychaudhuri", "cylindrical"): {"t": "T", "r": "L", "\\phi": "1", "z": "L", "\\Omega": "1/T"},
    ("som_raychaudhuri", "cartesian"): {"t": "T", "x": "L", "y": "L", "z": "L", "\\Omega": "1/T"},
    # Vuorio's universe and warped anti-de Sitter space keep c, with the twist Omega of d_t and the
    # curvature scale m of the plane it turns about both inverse lengths.
    ("vuorio_warped_ads", "cylindrical"): {"t": "T", "r": "L", "\\phi": "1", "\\Omega": "1/L", "m": "1/L"},
    ("vuorio_warped_ads", "disc"): {"t": "T", "R": "1", "\\phi": "1", "\\Omega": "1/L", "m": "1/L"},
    ("vuorio_warped_ads", "fibred"): {"t": "T", "\\sigma": "1", "u": "1", "\\Omega": "1/L", "m": "1/L"},
    ("vuorio_warped_ads", "horospherical"): {"t": "T", "x": "L", "y": "L", "\\Omega": "1/L", "m": "1/L"},
    # Maitra's dust keeps c. Its a is the one length of the solution, 8 pi G rho/c^2 = 1/a^2 on the
    # axis; s and gamma are numbers and k is a length, so that k dphi stands beside c dt.
    ("maitra_dust", "cylindrical"): {
        "t": "T", "r": "L", "\\phi": "1", "z": "L", "a": "L", "s": "1", "k": "L", "\\gamma": "1",
    },
    # Bonnor's rotating dust cloud keeps c, and its a is a length, a^2 = 2h = 2GJ/c^3, so that
    # a^2 rho^2/r^3 is a length beside c dt; r = sqrt(rho^2 + z^2) is the name the cylindrical chart defines.
    ("bonnor_rotating_dust", "cylindrical"): {"t": "T", "\\rho": "L", "\\phi": "1", "z": "L", "a": "L", "r": "L"},
    ("bonnor_rotating_dust", "spherical"): {"t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "a": "L"},
    # Schwarzschild-anti-de Sitter keeps Schwarzschild's r_s and anti-de Sitter's radius, which the
    # entry calls L and the dimensional pass calls L as well, as for anti-de Sitter space itself.
    # The Eddington-Finkelstein times u = ct - r_* and v = ct + r_* are lengths.
    ("schwarzschild_ads", "static"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "r_s": "L", "L": "L",
    },
    ("schwarzschild_ads", "eddington_finkelstein_outgoing"): {
        "u": "L", "r": "L", "\\theta": "1", "\\phi": "1", "r_s": "L", "L": "L",
    },
    ("schwarzschild_ads", "eddington_finkelstein_ingoing"): {
        "v": "L", "r": "L", "\\theta": "1", "\\phi": "1", "r_s": "L", "L": "L",
    },
    # The topological black holes keep anti-de Sitter's radius L and a mass parameter mu, a length,
    # which is Schwarzschild-anti-de Sitter's r_s at k = 1. The curvature k of the horizon's own
    # metric and its radial coordinate rho carry no dimension. Lemos's string has an angle phi and
    # a length z; the brane's x, y and z are lengths, as in anti-de Sitter's Poincare patch, and
    # its horizon z_h = L^2/r_h is one too.
    ("topological_black_hole", "static"): {
        "t": "T", "r": "L", "\\rho": "1", "\\phi": "1", "\\mu": "L", "L": "L", "k": "1",
    },
    ("topological_black_hole", "eddington_finkelstein_ingoing"): {
        "v": "L", "r": "L", "\\rho": "1", "\\phi": "1", "\\mu": "L", "L": "L", "k": "1",
    },
    ("topological_black_hole", "eddington_finkelstein_outgoing"): {
        "u": "L", "r": "L", "\\rho": "1", "\\phi": "1", "\\mu": "L", "L": "L", "k": "1",
    },
    ("topological_black_hole", "black_string"): {
        "t": "T", "r": "L", "\\phi": "1", "z": "L", "\\mu": "L", "L": "L",
    },
    ("topological_black_hole", "brane"): {
        "t": "T", "x": "L", "y": "L", "z": "L", "z_h": "L", "L": "L",
    },
    ("topological_black_hole", "hyperbolic"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "\\mu": "L", "L": "L",
    },
    # The anti-de Sitter soliton keeps anti-de Sitter's radius L and the radius r_0 of its tip, both
    # lengths. The circle tau and the flat directions are lengths, as in Horowitz and Myers's
    # (3.14); the Poincare chart's z = L^2/r and its tip z_0 = L^2/r_0 are lengths too, and the
    # polar chart's rho is the proper distance from the tip and its phi an angle.
    ("ads_soliton", "horowitz_myers"): {"t": "T", "r": "L", "\\tau": "L", "x": "L", "r_0": "L", "L": "L"},
    ("ads_soliton", "poincare"): {"t": "T", "z": "L", "\\tau": "L", "x": "L", "z_0": "L", "L": "L"},
    ("ads_soliton", "polar"): {"t": "T", "\\rho": "L", "\\phi": "1", "x": "L", "r_0": "L", "L": "L"},
    ("ads_soliton", "five_dimensional"): {
        "t": "T", "r": "L", "\\tau": "L", "x": "L", "y": "L", "r_0": "L", "L": "L",
    },
    ("ads_soliton", "three_dimensional"): {"t": "T", "r": "L", "\\tau": "L", "r_0": "L", "L": "L"},
    # The topological star keeps two lengths, r_S and r_B, Bah and Heidmann's, and its fifth
    # coordinate y is a length. The chart about the bubble has their rho, a pure number, the angle
    # psi = y/R_y, and r as a name for r_B + (r_B - r_S) rho^2/4. The black string's advanced time
    # v = ct + r_* is a length, and the extremal string keeps one length m = r_S = r_B, with the
    # isotropic radius rho = r - m a length.
    ("topological_star", "bah_heidmann"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "y": "L", "r_S": "L", "r_B": "L",
    },
    ("topological_star", "bubble"): {
        "t": "T", "\\rho": "1", "\\theta": "1", "\\phi": "1", "\\psi": "1", "r_S": "L", "r_B": "L", "r": "L",
    },
    ("topological_star", "eddington_finkelstein_ingoing"): {
        "v": "L", "r": "L", "\\theta": "1", "\\phi": "1", "y": "L", "r_S": "L", "r_B": "L",
    },
    ("topological_star", "extremal"): {"t": "T", "\\rho": "L", "\\theta": "1", "\\phi": "1", "y": "L", "m": "L"},
    ("topological_star", "einstein"): {"t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "r_S": "L", "r_B": "L"},
    # The charged black hole in de Sitter space keeps the three lengths of its parents: r_s, the
    # charge radius r_q and Lambda, a curvature. Its cosmological chart, which exists at
    # r_q = r_s/2, carries the Hubble rate H, a frequency with 3H^2/c^2 = Lambda, as de Sitter's
    # flat slicing does, and a comoving length rho.
    ("reissner_nordstrom_de_sitter", "static"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "r_s": "L", "r_q": "L", "\\Lambda": "1/L**2",
    },
    ("reissner_nordstrom_de_sitter", "eddington_finkelstein_outgoing"): {
        "u": "L", "r": "L", "\\theta": "1", "\\phi": "1", "r_s": "L", "r_q": "L", "\\Lambda": "1/L**2",
    },
    ("reissner_nordstrom_de_sitter", "eddington_finkelstein_ingoing"): {
        "v": "L", "r": "L", "\\theta": "1", "\\phi": "1", "r_s": "L", "r_q": "L", "\\Lambda": "1/L**2",
    },
    ("reissner_nordstrom_de_sitter", "cosmological"): {
        "\\tau": "T", "\\rho": "L", "\\theta": "1", "\\phi": "1", "r_s": "L", "H": "1/T",
    },
    # The charged black hole in anti-de Sitter space keeps the three lengths of its parents: r_s,
    # the charge radius r_q and anti-de Sitter's radius L. The Eddington-Finkelstein times
    # u = ct - r_* and v = ct + r_* are lengths.
    ("reissner_nordstrom_ads", "static"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "r_s": "L", "r_q": "L", "L": "L",
    },
    ("reissner_nordstrom_ads", "eddington_finkelstein_outgoing"): {
        "u": "L", "r": "L", "\\theta": "1", "\\phi": "1", "r_s": "L", "r_q": "L", "L": "L",
    },
    ("reissner_nordstrom_ads", "eddington_finkelstein_ingoing"): {
        "v": "L", "r": "L", "\\theta": "1", "\\phi": "1", "r_s": "L", "r_q": "L", "L": "L",
    },
    # Visser's thin shell wormhole joins two Schwarzschild exteriors at the throat radius a, a
    # length beside r_s. Through the throat, ell = +-(r - a) is a length, |ell| carries what ell
    # carries and sgn(ell) none, and the delta at the throat an inverse length.
    ("thin_shell_wormhole", "spherical"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "r_s": "L", "a": "L",
    },
    ("thin_shell_wormhole", "throat"): {
        "t": "T", "\\ell": "L", "\\theta": "1", "\\phi": "1", "r_s": "L", "a": "L",
    },
    # Hiscock's evaporating black hole in his own units, G = c = 1: the null times v and u are
    # lengths and so is the mass m = GM/c^2, a function of one of them, so that a dot on m carries
    # no dimension. R, the areal radius of the surface of pair creation, enters the domains alone.
    ("hiscock", "ingoing"): {
        "v": "L", "r": "L", "\\theta": "1", "\\phi": "1", "m": "L", "R": "L",
    },
    ("hiscock", "outgoing"): {
        "u": "L", "r": "L", "\\theta": "1", "\\phi": "1", "m": "L", "R": "L",
    },
    ("hiscock", "flat"): {
        "u": "L", "v": "L", "\\theta": "1", "\\phi": "1",
    },
    # Penrose's spherical impulsive wave for a snapping string: the null coordinates U and V of
    # Podolsky and Griffiths are lengths and rho is a stereographic coordinate on the wave front, a
    # pure number; u and r of the retarded chart are lengths; k and beta are pure numbers, and so
    # is the step Theta each continuous chart defines.
    ("penrose_impulsive_wave", "null"): {
        "U": "L", "V": "L", "\\rho": "1", "\\phi": "1", "k": "1", "\\Theta": "1",
    },
    ("penrose_impulsive_wave", "retarded"): {
        "u": "L", "r": "L", "\\theta": "1", "\\phi": "1", "k": "1", "\\Theta": "1",
    },
    ("penrose_impulsive_wave", "behind"): {"t": "T", "r": "L", "\\theta": "1", "\\phi": "1"},
    ("penrose_impulsive_wave", "ahead"): {"T": "T", "R": "L", "\\phi": "1", "z": "L", "\\beta": "1"},
    # Israel's shell of dust: flat space inside in the time T of its centre, Schwarzschild's vacuum
    # outside in Schwarzschild's time t and in the advanced time v, a length. R, the areal radius of
    # the shell, and mu = G m/c^2, the rest mass of the dust as a length, enter the domains and the
    # shell's motion alone.
    ("shell_void", "exterior_proper"): {"t":"L", "\\chi":"1", "\\theta":"1", "\\phi":"1", "a":"L", "f":"1", "k":"1", "\\chi_s":"1"},
    ("shell_void", "interior"): {"T": "L", "r": "L", "\\theta": "1", "\\phi": "1", "R": "L"},
    ("shell_void", "exterior"): {"t": "L", "r": "L", "\\theta": "1", "\\phi": "1", "a": "1", "k": "L^-2", "R": "L"},
    ("israel_shell", "interior"): {
        "T": "T", "r": "L", "\\theta": "1", "\\phi": "1", "r_s": "L", "\\mu": "L", "R": "L",
    },
    ("israel_shell", "exterior"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "r_s": "L", "\\mu": "L", "R": "L",
    },
    ("israel_shell", "exterior_ingoing"): {
        "v": "L", "r": "L", "\\theta": "1", "\\phi": "1", "r_s": "L", "\\mu": "L", "R": "L",
    },
    # The charged shell of dust: flat space inside in the time T of its centre, Reissner-Nordstrom's
    # field outside in its static time t, in the advanced and retarded times v and u, both lengths,
    # and in the isotropic radius rho, where a = r_s/4 + r_q/2 and b = r_s/4 - r_q/2. R, the areal
    # radius of the shell, epsilon, its isotropic radius, and mu = G m/c^2, the rest mass of the dust as a length, enter the domains and the
    # shell's motion alone.
    ("charged_shell", "interior"): {
        "T": "T", "r": "L", "\\theta": "1", "\\phi": "1", "r_s": "L", "r_q": "L", "\\mu": "L", "R": "L",
    },
    ("charged_shell", "exterior"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "r_s": "L", "r_q": "L", "\\mu": "L", "R": "L",
    },
    ("charged_shell", "exterior_ingoing"): {
        "v": "L", "r": "L", "\\theta": "1", "\\phi": "1", "r_s": "L", "r_q": "L", "\\mu": "L", "R": "L",
    },
    ("charged_shell", "exterior_outgoing"): {
        "u": "L", "r": "L", "\\theta": "1", "\\phi": "1", "r_s": "L", "r_q": "L", "\\mu": "L", "R": "L",
    },
    ("charged_shell", "exterior_isotropic"): {
        "t": "T", "\\rho": "L", "\\theta": "1", "\\phi": "1", "a": "L", "b": "L", "\\mu": "L",
        "\\epsilon": "L",
    },
    # The gravastar: a ball of de Sitter space of radius L inside a thin shell at the areal radius
    # R, Schwarzschild's vacuum of radius r_s outside. C is a name for the number that makes g_tt
    # continuous across the shell, and the tortoise coordinate x of the interior is a length.
    ("gravastar", "interior"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "r_s": "L", "R": "L", "L": "L", "C": "1",
    },
    ("gravastar", "interior_tortoise"): {
        "t": "T", "x": "L", "\\theta": "1", "\\phi": "1", "r_s": "L", "R": "L", "L": "L", "C": "1",
    },
    ("gravastar", "exterior"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "r_s": "L", "R": "L",
    },
    # Damour and Solodukhin's lambda is a pure number beside Schwarzschild's r_s. Bueno and his
    # collaborators' rho is a pure number too, and Einstein and Rosen's u has u^2 = r - r_s, so it
    # carries the square root of a length.
    ("damour_solodukhin", "spherical"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "r_s": "L", "\\lambda": "1",
    },
    ("damour_solodukhin", "rescaled"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "r_s": "L", "\\lambda": "1",
    },
    ("damour_solodukhin", "throat"): {
        "t": "T", "\\rho": "1", "\\theta": "1", "\\phi": "1", "r_s": "L", "\\lambda": "1",
    },
    ("damour_solodukhin", "isotropic"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "r_s": "L", "\\lambda": "1",
    },
    ("damour_solodukhin", "einstein_rosen"): {
        "t": "T", "u": "L**(1/2)", "\\theta": "1", "\\phi": "1", "r_s": "L", "\\lambda": "1",
    },
    # Blazquez-Salcedo, Knoll and Radu's charge and mass are lengths beside the throat's radius r_0.
    # Bronnikov and Kim's u has u^2 = r - r_0, so it carries the square root of a length, and the
    # compact x, r = r_0/(1 - x^2), is a pure number.
    ("einstein_dirac_maxwell_wormhole", "areal"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "r_0": "L", "Q_e": "L", "M": "L",
    },
    ("einstein_dirac_maxwell_wormhole", "bronnikov_kim"): {
        "t": "T", "u": "L**(1/2)", "\\theta": "1", "\\phi": "1", "r_0": "L", "Q_e": "L", "M": "L",
    },
    ("einstein_dirac_maxwell_wormhole", "compact"): {
        "t": "T", "x": "1", "\\theta": "1", "\\phi": "1", "r_0": "L", "Q_e": "L", "M": "L",
    },
    # James, von Tunzelmann, Franklin and Thorne's Dneg wormhole: rho, a and M are lengths, and so
    # is the proper distance l; the flare's x = 2(l - a)/(pi M) is a pure number.
    ("interstellar_wormhole", "proper_distance"): {
        "t": "T", "\\ell": "L", "\\theta": "1", "\\phi": "1", "r": "L", "\\rho": "L", "a": "L", "M": "L",
    },
    ("interstellar_wormhole", "cylinder"): {
        "t": "T", "\\ell": "L", "\\theta": "1", "\\phi": "1", "\\rho": "L", "a": "L",
    },
    ("interstellar_wormhole", "flare"): {
        "t": "T", "\\ell": "L", "\\theta": "1", "\\phi": "1", "\\rho": "L", "a": "L", "M": "L", "x": "1", "r": "L",
    },
    # Einstein and Rosen's u has u^2 = r - r_s on the neutral bridge, the square root of a length, and
    # u^2 = r^2 - r_q^2 on the charged one, a length; r_q is the charge radius, their epsilon/sqrt(2).
    ("einstein_rosen_bridge", "bridge"): {
        "t": "T", "u": "L**(1/2)", "\\theta": "1", "\\phi": "1", "r_s": "L",
    },
    ("einstein_rosen_bridge", "spherical"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "r_s": "L",
    },
    ("einstein_rosen_bridge", "isotropic"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "r_s": "L",
    },
    ("einstein_rosen_bridge", "charged_spherical"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "r_s": "L", "r_q": "L",
    },
    ("einstein_rosen_bridge", "charged_bridge"): {
        "t": "T", "u": "L", "\\theta": "1", "\\phi": "1", "r_q": "L",
    },
    # Simpson and Visser's a is a length beside Schwarzschild's r_s, and rho = sqrt(r^2 + a^2) is the
    # areal radius, a name in the charts that keep their r and the coordinate of Tsukamoto's.
    ("simpson_visser", "spherical"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "r_s": "L", "a": "L", "\\rho": "L",
    },
    ("simpson_visser", "areal"): {
        "t": "T", "\\rho": "L", "\\theta": "1", "\\phi": "1", "r_s": "L", "a": "L",
    },
    ("simpson_visser", "eddington_finkelstein_outgoing"): {
        "u": "L", "r": "L", "\\theta": "1", "\\phi": "1", "r_s": "L", "a": "L", "\\rho": "L",
    },
    ("simpson_visser", "eddington_finkelstein_ingoing"): {
        "v": "L", "r": "L", "\\theta": "1", "\\phi": "1", "r_s": "L", "a": "L", "\\rho": "L",
    },
    ("stockum_dust", "cylindrical"): {
        "t": "L", "r": "L", "\\phi": "1", "z": "L", "R": "L",
    },
    # The mass is folded into the length m = GM/c^2 and the NUT parameter is a length
    # beside it, so the cross term 2l cos(theta) carries the one length g_{t phi} wants.
    ("taub_nut", "spherical"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "m": "L", "l": "L",
    },
    # The star of infinite central density. The length a sets the rate of t, r_b is the radius of
    # a star, Z = g^rr and Tolman's exponent n are pure numbers.
    ("misner_zapolsky", "areal"): {"t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "a": "L"},
    ("misner_zapolsky", "tolman_v"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "r_b": "L", "Z": "1",
    },
    ("misner_zapolsky", "tolman_vi"): {"t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "r_b": "L"},
    ("misner_zapolsky", "power_law"): {"t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "n": "1", "a": "L"},
    # Bowers and Liang's anisotropic sphere of uniform density. Their exponent Q is a pure number,
    # and so are Z = g^rr, N, whose power 1/Q is -g_tt, and P = 1 + 3p_r/(rho c^2).
    ("bowers_liang", "areal"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "R": "L", "r_s": "L", "Q": "1", "Z": "1", "N": "1",
        "P": "1",
    },
    # Tolman's solution VII. Lattimer and Prakash's compactness beta = GM/(Rc^2) is a pure number,
    # and so are Z = g^rr and the phase psi; Tolman's R and A are lengths and his B and C numbers.
    ("tolman_vii", "spherical"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "R": "L", "r_s": "L", "\\beta": "1", "Z": "1",
        "\\psi": "1",
    },
    ("tolman_vii", "tolman"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "R": "L", "A": "L", "B": "1", "C": "1", "Z": "1",
        "\\psi": "1",
    },
    # Buchdahl's fluid sphere: r_s = 2GM/c^2 and the core radius R = k^(-1/2) are lengths, and his
    # f = r_s/(4 sqrt(r^2 + R^2)) is a pure number.
    ("buchdahl_sphere", "isotropic"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "r_s": "L", "R": "L", "f": "1",
    },
    # The mass function is folded into a length, m = GM(r)/c^2, as Oppenheimer and Volkoff
    # folded it into their u, so 1 - 2m/r is dimensionless and the redshift function, sitting
    # in an exponential, is dimensionless too.
    ("tov", "spherical"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "\\Phi": "1", "m": "L",
    },
    # A boson star: the lapse and the radial function are pure numbers, as the alpha and a of
    # Liebling and Palenzuela's line element are.
    ("boson_star", "areal"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "\\alpha": "1", "a": "1",
    },
    ("boson_star", "isotropic"): {
        "t": "T", "R": "L", "\\theta": "1", "\\phi": "1", "\\alpha": "1", "\\psi": "1",
    },
    # Kodama's kink: the redshift functions eta and h and the radial function alpha sit in
    # exponentials and are pure numbers; the areal radius r and the throat radius r_0 are lengths.
    ("kodama_kink", "areal"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "\\eta": "1", "\\alpha": "1", "r_0": "L",
    },
    ("kodama_kink", "proper_radial"): {
        "t": "T", "R": "L", "\\theta": "1", "\\phi": "1", "h": "1", "r": "L", "r_0": "L",
    },
    # Herdeiro and Radu's black holes with scalar hair: F_0, F_1 and F_2 are pure numbers, and W is
    # the rate of dragging per unit of ct, an inverse length, since the line element writes W c dt.
    # In the Kerr member both constants r_H and b = -c_t are lengths.
    ("kerr_scalar_hair", "herdeiro_radu"): {
        "t": "T", "r": "L", "\\theta": "1", "\\varphi": "1", "r_H": "L", "F_0": "1", "F_1": "1", "F_2": "1",
        "W": "1/L",
    },
    ("kerr_scalar_hair", "kerr_member"): {
        "t": "T", "r": "L", "\\theta": "1", "\\varphi": "1", "r_H": "L", "b": "L",
    },
    # Bartnik and McKinnon's soliton: the mass function is a length, as TOV's is, the amplitude w
    # of the Yang-Mills field and the lapse's sigma are numbers, and ell is the one length of the
    # theory, which stands in the field equations and in no component. The isotropic chart's f
    # and h are numbers, and the tortoise chart's r is the areal radius as a function of xi.
    ("bartnik_mckinnon", "areal"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "\\sigma": "1", "m": "L", "w": "1", "\\ell": "L",
    },
    ("bartnik_mckinnon", "isotropic"): {
        "t": "T", "\\rho": "L", "\\theta": "1", "\\phi": "1", "f": "1", "h": "1",
    },
    ("bartnik_mckinnon", "flow"): {
        "t": "T", "\\tau": "1", "\\theta": "1", "\\phi": "1", "A": "1", "N": "1", "r": "L",
    },
    ("bartnik_mckinnon", "tortoise"): {
        "t": "T", "\\xi": "L", "\\theta": "1", "\\phi": "1", "F": "1", "r": "L",
    },
    # Einstein's cluster: the mass function is a length, as TOV's is, and here it is a name for
    # r^2 Phi'/(1 + 2r Phi'), the mass that leaves no pressure along the radius. V and sigma are
    # ratios, a speed to c and a mass to twice an isotropic radius, and chi is an angle on the
    # three sphere of radius a.
    ("einstein_cluster", "areal"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "\\Phi": "1", "m": "L",
    },
    ("einstein_cluster", "constant_speed"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "V": "1", "R": "L",
    },
    ("einstein_cluster", "isotropic"): {
        "t": "T", "\\rho": "L", "\\theta": "1", "\\phi": "1", "\\sigma": "1", "\\rho_0": "L",
        "a": "1", "b": "1",
    },
    ("einstein_cluster", "uniform"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "r_s": "L", "R": "L",
    },
    ("einstein_cluster", "hyperspherical"): {
        "t": "T", "\\chi": "1", "\\theta": "1", "\\phi": "1", "a": "L", "\\chi_0": "1",
    },
    # Szekeres's shells: the areal radius R and the shell label r are lengths, as Tolman-Bondi's
    # are, the stereographic coordinates p and q of a shell are pure numbers, and so are S, P
    # and Q, which place and scale them, the sign epsilon, the energy function f, and E.
    ("szekeres", "stereographic"): {
        "t": "T", "r": "L", "p": "1", "q": "1", "R": "L", "f": "1", "S": "1", "P": "1", "Q": "1",
        "\\epsilon": "1", "E": "1",
    },
    ("szekeres", "axisymmetric"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "R": "L", "f": "1", "S": "1",
    },
    # The areal radius is a length and the comoving shell label r is a length beside it,
    # so \partial_r R is dimensionless and the energy function has to be dimensionless
    # as well, which is what leaves 1 + 2E a pure number. Every published derivative of R
    # is taken along a chart coordinate, so each order of it carries one inverse length.
    ("tolman_bondi", "comoving_synchronous"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "R": "L", "E": "1",
    },
    # Hayward's regular black hole keeps his own two lengths, the mass m = GM/c^2 and the length
    # ell of the de Sitter core. The Eddington-Finkelstein times u = ct - r_* and v = ct + r_* are
    # lengths, and so is the advanced time of the chart whose mass is a function m(v), where a
    # dot on m is a derivative along v and carries no dimension.
    # The quantum Oppenheimer-Snyder black hole: the mass is the length m = GM/c^2 and alpha, the
    # quantum correction, is an area, so alpha m^2/r^4 is a pure number. Inside, the scale factor is a
    # pure number, 1 at the bounce, and chi is the areal radius of a shell then, a length.
    ("quantum_oppenheimer_snyder", "static"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "m": "L", "\\alpha": "L**2",
    },
    ("quantum_oppenheimer_snyder", "painleve_gullstrand"): {
        "\\tau": "T", "r": "L", "\\theta": "1", "\\phi": "1", "m": "L", "\\alpha": "L**2",
    },
    ("quantum_oppenheimer_snyder", "eddington_finkelstein_ingoing"): {
        "v": "L", "r": "L", "\\theta": "1", "\\phi": "1", "m": "L", "\\alpha": "L**2",
    },
    ("quantum_oppenheimer_snyder", "eddington_finkelstein_outgoing"): {
        "u": "L", "r": "L", "\\theta": "1", "\\phi": "1", "m": "L", "\\alpha": "L**2",
    },
    ("quantum_oppenheimer_snyder", "interior_comoving"): {
        "\\tau": "T", "\\chi": "L", "\\theta": "1", "\\phi": "1", "a": "1", "\\chi_0": "L",
    },
    # Flat space with supertranslation hair: Compere and Long's supertranslation field C and the
    # shear sigma it makes are lengths, and so are rho, the distance along each straight line of
    # constant angle, the luminosity distance r and the shorthand W = sqrt(r^2 + sigma^2).
    ("supertranslation_hair", "static"): {
        "t": "T", "\\rho": "L", "\\theta": "1", "\\phi": "1", "C": "L", "B": "L",
    },
    ("supertranslation_hair", "bondi_retarded"): {
        "u": "T", "r": "L", "\\theta": "1", "\\phi": "1", "\\sigma": "L", "W": "L",
    },
    ("supertranslation_hair", "bondi_advanced"): {
        "v": "T", "r": "L", "\\theta": "1", "\\phi": "1", "\\sigma": "L", "W": "L",
    },
    # The eye of the storm: Kerr's metric with the mass m = GM/c^2 suppressed by e^{-ell/r}, so
    # ell is a length and the exponent is a number.
    ("eye_of_the_storm", "boyer_lindquist"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "m": "L", "a": "L", "\\ell": "L",
        "\\Sigma": "L**2", "\\Delta": "L**2",
    },
    ("eye_of_the_storm", "kerr_ingoing"): {
        "v": "L", "r": "L", "\\theta": "1", "\\tilde\\phi": "1", "m": "L", "a": "L", "\\ell": "L",
        "\\Sigma": "L**2", "\\Delta": "L**2",
    },
    # Ozsvath's homogeneous dust: class I has a time t and the length L, the other charts are
    # dimensionless, Wylleman's t among them, with the length a carrying the scale.
    ("ozsvath_dust", "class_one"): {
        "t": "T", "\\theta": "1", "\\phi": "1", "\\psi": "1", "L": "L", "\\beta": "1", "\\epsilon": "1",
    },
    ("ozsvath_dust", "class_two"): {
        "x": "1", "y": "1", "z": "1", "u": "1", "a": "L", "\\gamma": "1", "s": "1", "\\lambda": "1",
    },
    ("ozsvath_dust", "class_three"): {
        "x": "1", "y": "1", "z": "1", "u": "1", "a": "L", "s": "1",
    },
    ("ozsvath_dust", "farnsworth_kerr"): {
        "t": "1", "x": "1", "y": "1", "u": "1", "a": "L", "v": "1",
    },
    ("hayward", "static"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "m": "L", "\\ell": "L",
    },
    ("hayward", "eddington_finkelstein_ingoing"): {
        "v": "L", "r": "L", "\\theta": "1", "\\phi": "1", "m": "L", "\\ell": "L",
    },
    ("hayward", "eddington_finkelstein_outgoing"): {
        "u": "L", "r": "L", "\\theta": "1", "\\phi": "1", "m": "L", "\\ell": "L",
    },
    ("hayward", "evaporating"): {
        "v": "L", "r": "L", "\\theta": "1", "\\phi": "1", "m": "L", "\\ell": "L",
    },
    # Bueno, Cano and Hennigar's regular black hole of five dimensions quotes its mass parameter
    # m = 8GM/(3 pi c^2), an area, as Tangherlini's radius r_0 with r_0^2 = m, and the coupling
    # alpha of the quasi-topological tower, an area, as the core's length l with l^2 = alpha; the
    # Eddington-Finkelstein times u = ct - r_* and v = ct + r_* are lengths.
    ("regular_black_hole_pure_gravity", "spherical"): {
        "t": "T", "r": "L", "\\psi": "1", "\\theta": "1", "\\phi": "1", "r_0": "L", "\\ell": "L",
    },
    ("regular_black_hole_pure_gravity", "eddington_finkelstein_ingoing"): {
        "v": "L", "r": "L", "\\psi": "1", "\\theta": "1", "\\phi": "1", "r_0": "L", "\\ell": "L",
    },
    ("regular_black_hole_pure_gravity", "eddington_finkelstein_outgoing"): {
        "u": "L", "r": "L", "\\psi": "1", "\\theta": "1", "\\phi": "1", "r_0": "L", "\\ell": "L",
    },
    # Mass inflation: the charged Vaidya metric in advanced time, a length, with the mass function
    # and the charge lengths too; the double null chart, where the areal radius r is a function of
    # the two null coordinates and sigma a pure number; and Brady and Smith's chart of advanced time
    # and radius, whose two functions are pure numbers.
    ("mass_inflation", "ingoing"): {
        "v": "L", "r": "L", "\\theta": "1", "\\phi": "1", "m": "L", "r_q": "L",
    },
    ("mass_inflation", "double_null"): {
        "u": "L", "v": "L", "\\theta": "1", "\\phi": "1", "r": "L", "\\sigma": "1",
    },
    ("mass_inflation", "advanced"): {
        "v": "L", "r": "L", "\\theta": "1", "\\phi": "1", "g": "1", "h": "1",
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
    # Bonnor and Vaidya's charged radiating star: the null time, the mass m and the charge q are
    # lengths, m = GM/c^2 and q^2 = GQ^2/(4 pi epsilon_0 c^4), so that r_s = 2m and r_q = q.
    ("bonnor_vaidya", "eddington_finkelstein_outgoing"): {
        "u": "L", "r": "L", "\\theta": "1", "\\phi": "1", "m": "L", "q": "L",
    },
    ("bonnor_vaidya", "eddington_finkelstein_ingoing"): {
        "v": "L", "r": "L", "\\theta": "1", "\\phi": "1", "m": "L", "q": "L",
    },
    # The homothetic chart: V and R are lengths, as the scales M and Q are, and the rate mu a number.
    ("bonnor_vaidya", "homothetic"): {
        "V": "L", "R": "L", "\\theta": "1", "\\phi": "1", "M": "L", "Q": "L", "\\mu": "1",
    },
    # Sultana and Dyer's black hole: the conformal time eta is a length, as FRW's is, with
    # a d(eta) = c dt, and so are eta_0, where the scale factor a = eta^2/eta_0^2 is one, and the
    # advanced time v. Schwarzschild's t is a time, and eta a name for ct + r_s ln(r/r_s - 1) there.
    ("sultana_dyer", "kerr_schild"): {
        "\\eta": "L", "r": "L", "\\theta": "1", "\\phi": "1", "r_s": "L", "\\eta_0": "L",
    },
    ("sultana_dyer", "schwarzschild_time"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "r_s": "L", "\\eta_0": "L", "\\eta": "L",
    },
    ("sultana_dyer", "eddington_finkelstein_ingoing"): {
        "v": "L", "r": "L", "\\theta": "1", "\\phi": "1", "r_s": "L", "\\eta_0": "L",
    },
    # Kopczynski and Trautman's universe: l = c/sqrt(6 pi G rho_0) is a length, the scale factor a pure
    # number that is one at the bounce, and the conformal time a length, as FRW's is.
    ("kopczynski_trautman", "comoving_cartesian"): {
        "t": "T", "x": "L", "y": "L", "z": "L", "\\ell": "L", "a": "1",
    },
    ("kopczynski_trautman", "comoving_spherical"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "\\ell": "L", "a": "1",
    },
    ("kopczynski_trautman", "conformal"): {
        "\\eta": "L", "r": "L", "\\theta": "1", "\\phi": "1", "a": "1",
    },
    # Ellis's small universes: on the torus x, y and z and their periods are comoving lengths and the
    # scale factor a pure number, as FRW's is, with the conformal time a length; in the hyperbolic
    # chart chi is a pure number and the scale factor is the radius of curvature of space, a length.
    ("small_universes", "torus"): {
        "t": "T", "x": "L", "y": "L", "z": "L", "a": "1", "L_1": "L", "L_2": "L", "L_3": "L",
    },
    ("small_universes", "torus_conformal"): {
        "\\eta": "L", "x": "L", "y": "L", "z": "L", "a": "1", "L_1": "L", "L_2": "L", "L_3": "L",
    },
    ("small_universes", "hyperbolic"): {"t": "T", "\\chi": "1", "\\theta": "1", "\\phi": "1", "a": "L"},
    # The horn's x, y and z are pure numbers, y and z flat coordinates on the horospheres of hyperbolic
    # space with pure numbers for periods, and the scale factor is again the radius of curvature.
    ("small_universes", "horn"): {"t": "T", "x": "1", "y": "1", "z": "1", "a": "L", "b_2": "1", "b_3": "1"},
    # The Poincare dodecahedral universe: the coordinates of the unit 3-sphere are pure numbers and the
    # scale factor is the radius of curvature of space, a length, as Aurich, Lustig and Steiner's R(t)
    # is; their conformal time, a d eta = c dt, is then a pure number too.
    ("poincare_dodecahedral", "comoving"): {"t": "T", "\\chi": "1", "\\theta": "1", "\\phi": "1", "a": "L"},
    ("poincare_dodecahedral", "conformal"): {"\\eta": "1", "\\chi": "1", "\\theta": "1", "\\phi": "1", "a": "L"},
    ("poincare_dodecahedral", "toroidal"): {"t": "T", "v": "1", "\\alpha": "1", "\\gamma": "1", "a": "L"},
    # The tilted universes of Bianchi's type V: r, y and z are pure numbers, as the horn's are, the
    # time u = ct + Cr of the surfaces of homogeneity is a length, and so are the tilt C, the two
    # scale factors X and Y, and W, which sets the density of Farnsworth's dust; his parameter eta
    # is a pure number. The flat model's inertial coordinates are a time and three lengths.
    ("tilted_universes", "homogeneous"): {"u": "L", "r": "1", "y": "1", "z": "1", "C": "L", "X": "L", "Y": "L"},
    ("tilted_universes", "farnsworth"): {"\\eta": "1", "r": "1", "y": "1", "z": "1", "C": "L", "W": "L", "X": "L"},
    ("tilted_universes", "flat_model"): {"u": "L", "r": "1", "y": "1", "z": "1", "C": "L"},
    ("tilted_universes", "inertial"): {"T": "T", "x": "L", "\\xi": "L", "\\zeta": "L", "C": "L"},
    # Hartle and Thorne's exterior: the mass and the spin per unit mass are lengths and the
    # quadrupole moment per unit mass an area, so that Kerr's value is q = a^2. R, the star's
    # radius, enters the domain alone, and the rest are names for the functions of r and theta the
    # line element is written in.
    ("hartle_thorne", "hartle_thorne"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "m": "L", "a": "L", "q": "L**2", "R": "L",
        "L": "1", "A": "1", "B": "1", "P_2": "1", "F": "1", "h_2": "1", "j_2": "1", "k_2": "1",
    },
    ("hartle_thorne", "lense_thirring"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "m": "L", "a": "L", "R": "L",
    },
    ("hartle_thorne", "painleve_gullstrand"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "m": "L", "a": "L", "R": "L",
    },
    # The parametrised post-Newtonian metric: the mass is the length m = GM/c^2, the spin the
    # length a = J/(Mc), and Eddington's numbers beta and gamma, the preferred frame parameter
    # alpha_1 and the dragging strength Delta are pure numbers.
    ("ppn_metric", "isotropic"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "m": "L", "\\beta": "1", "\\gamma": "1", "R": "L",
    },
    ("ppn_metric", "cartesian"): {
        "t": "T", "x": "L", "y": "L", "z": "L", "m": "L", "\\beta": "1", "\\gamma": "1", "R": "L", "r": "L",
    },
    ("ppn_metric", "areal"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "m": "L", "\\beta": "1", "\\gamma": "1", "R": "L",
    },
    ("ppn_metric", "rotating"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "m": "L", "a": "L", "\\beta": "1", "\\gamma": "1",
        "\\alpha_1": "1", "R": "L", "\\Delta": "1",
    },
    # Einstein, Infeld and Hoffmann's field of many bodies: the potentials U, psi and the three
    # components of V are pure numbers, each mass standing in them as the length Gm/c^2 over a
    # distance, and the superpotential chi, the sum of each such length times a distance, is an area.
    ("eih_many_bodies", "harmonic"): {
        "t": "T", "x": "L", "y": "L", "z": "L", "U": "1", "\\psi": "1", "\\chi": "L**2",
        "V_1": "1", "V_2": "1", "V_3": "1",
    },
    ("eih_many_bodies", "standard"): {
        "t": "T", "x": "L", "y": "L", "z": "L", "U": "1", "\\psi": "1", "\\chi": "L**2",
        "V_1": "1", "V_2": "1", "V_3": "1",
    },
    # The perturbed Friedmann universe: the conformal time eta, with a deta = c dt, and the comoving
    # x, y and z are lengths, the scale factor and the potentials pure numbers, the curvature k an
    # inverse area as in FRW's charts, and E, whose second derivatives along space perturb the metric
    # of space, an area.
    ("bardeen_perturbations", "longitudinal"): {
        "\\eta": "L", "x": "L", "y": "L", "z": "L", "a": "1", "\\Phi": "1", "\\Psi": "1",
    },
    ("bardeen_perturbations", "longitudinal_curved"): {
        "\\eta": "L", "x": "L", "y": "L", "z": "L", "a": "1", "k": "1/L**2", "\\Omega": "1", "\\Phi": "1", "\\Psi": "1",
    },
    ("bardeen_perturbations", "synchronous"): {
        "\\eta": "L", "x": "L", "y": "L", "z": "L", "a": "1", "\\psi": "1", "E": "L**2",
    },
    # Three dimensions, where 4Gm/c^2 is a pure number, so alpha = 1 - 4Gm/c^2 is one too. The
    # isotropic radius and the isotropic x and y are lengths, measured in the arbitrary length ell,
    # as are the coordinate distances rho_1 and rho_2 from two particles at rest and half their
    # separation d; the moving particle's v is a speed and its gamma a number.
    ("point_particle_2plus1", "conical"): {"t": "T", "r": "L", "\\phi": "1", "\\alpha": "1"},
    ("point_particle_2plus1", "wedge"): {"t": "T", "r": "L", "\\theta": "1", "\\alpha": "1"},
    ("point_particle_2plus1", "circumference"): {"t": "T", "R": "L", "\\phi": "1", "\\alpha": "1"},
    ("point_particle_2plus1", "planet"): {"t": "T", "\\chi": "1", "\\phi": "1", "a": "L", "\\chi_0": "1"},
    ("point_particle_2plus1", "isotropic"): {"t": "T", "\\rho": "L", "\\phi": "1", "\\alpha": "1", "\\ell": "L"},
    ("point_particle_2plus1", "two_bodies"): {
        "t": "T", "x": "L", "y": "L", "\\alpha_1": "1", "\\alpha_2": "1", "d": "L", "\\ell": "L",
        "\\rho_1": "L", "\\rho_2": "L", "\\Omega": "1",
    },
    ("point_particle_2plus1", "moving"): {
        "t": "T", "x": "L", "y": "L", "\\alpha": "1", "v": "L/T", "\\gamma": "1",
    },
    # Siklos's waves: every coordinate of his chart is a length, the null ones included, and the
    # profile H is a pure number, as Kaigorodov's x^3/L^3 is. The disc's xi and eta are lengths and
    # its profile h = LH/x is pure, as are the two functions p and q. The homogeneous form's Z is a
    # pure number and its k one too, and the Kundt form's x and y are pure, L times them Siklos's.
    ("siklos", "siklos"): {"u": "L", "v": "L", "x": "L", "y": "L", "L": "L", "H": "1"},
    ("siklos", "ozsvath_robinson_rozga"): {
        "u": "L", "v": "L", "\\xi": "L", "\\eta": "L", "L": "L", "h": "1", "p": "1", "q": "1",
    },
    ("siklos", "kaigorodov"): {"u": "L", "v": "L", "x": "L", "y": "L", "L": "L"},
    ("siklos", "kaigorodov_poincare"): {"t": "T", "x": "L", "y": "L", "z": "L", "L": "L"},
    ("siklos", "kaigorodov_horospheric"): {"u": "L", "v": "L", "y": "L", "\\rho": "L", "L": "L"},
    ("siklos", "kaigorodov_stationary"): {"u": "L", "v": "L", "y": "L", "\\rho": "L", "L": "L"},
    ("siklos", "kaigorodov_homogeneous"): {"U": "L", "X": "L", "y": "L", "Z": "1", "L": "L", "k": "1"},
    ("siklos", "kaigorodov_kundt"): {"U": "L", "V": "L", "x": "1", "y": "1", "L": "L"},
    # Kundt's waves: u and v are pure numbers, u the tangent of half the angle a wave front makes,
    # so the profile G is a length and the affine parameter w = 2x^2 v an area; the simplest wave's
    # l is a length. On flat space T is a time and the three names are x, a length, and the pure u
    # and v. With a cosmological constant the chart is the disc's of Siklos's waves, u and v
    # lengths, so alpha is pure, beta an inverse length and kappa, like Lambda, an inverse area.
    ("kundt_waves", "kundt"): {"u": "1", "w": "L^2", "x": "L", "y": "L", "G": "L"},
    ("kundt_waves", "podolsky_belan"): {"u": "1", "v": "1", "x": "L", "y": "L", "G": "L"},
    ("kundt_waves", "simplest_wave"): {"u": "1", "v": "1", "x": "L", "y": "L", "\\ell": "L"},
    ("kundt_waves", "kerr_schild"): {
        "T": "T", "X": "L", "Y": "L", "Z": "L", "\\ell": "L", "x": "L", "u": "1", "v": "1",
    },
    ("kundt_waves", "ozsvath_robinson_rozga"): {
        "u": "L", "v": "L", "\\xi": "L", "\\eta": "L", "\\Lambda": "1/L^2", "\\alpha": "1", "\\beta": "1/L",
        "h": "1", "p": "1", "q": "1", "\\kappa": "1/L^2",
    },
    # Schrodinger spacetime: t is the time of the nonrelativistic theory and every other
    # coordinate a length, the null xi included. beta is a length, so that beta^2/r^2 is pure, the
    # dynamical exponent z a pure number and the term h = (beta/r)^{2z - 2} it names one too, and
    # the trap's omega a frequency.
    ("schrodinger_spacetime", "poincare"): {
        "t": "T", "\\xi": "L", "x": "L", "r": "L", "L": "L", "\\beta": "L",
    },
    ("schrodinger_spacetime", "inverse_radius"): {
        "t": "T", "\\xi": "L", "x": "L", "\\rho": "L", "L": "L", "\\beta": "L",
    },
    ("schrodinger_spacetime", "global"): {
        "T": "T", "V": "L", "X": "L", "R": "L", "L": "L", "\\beta": "L", "\\omega": "1/T",
    },
    ("schrodinger_spacetime", "dynamical_exponent"): {
        "t": "T", "\\xi": "L", "x": "L", "r": "L", "L": "L", "\\beta": "L", "z": "1", "h": "1",
    },
    ("schrodinger_spacetime", "poincare_5d"): {
        "t": "T", "\\xi": "L", "x_1": "L", "x_2": "L", "r": "L", "L": "L", "\\beta": "L",
    },
    ("schrodinger_spacetime", "poincare_6d"): {
        "t": "T", "\\xi": "L", "x_1": "L", "x_2": "L", "x_3": "L", "r": "L", "L": "L", "\\beta": "L",
    },
    # Haggard and Rovelli's bounce of a shell of light: the null coordinates of the flat interior,
    # u = ct - r and v = ct + r, are lengths; Kruskal's U and V are numbers, with the factor
    # 4 r_s^3 e^(-r/r_s)/r carrying the area; Schwarzschild's and Painleve and Gullstrand's t is a
    # time; Lemaitre's tau is a time and rho a length.
    ("black_to_white_hole", "interior"): {"u": "L", "v": "L", "\\theta": "1", "\\phi": "1"},
    ("black_to_white_hole", "kruskal"): {
        "U": "1", "V": "1", "\\theta": "1", "\\phi": "1", "r_s": "L", "r": "L",
    },
    ("black_to_white_hole", "schwarzschild"): {"t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "r_s": "L"},
    ("black_to_white_hole", "painleve_gullstrand_ingoing"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "r_s": "L",
    },
    ("black_to_white_hole", "painleve_gullstrand_outgoing"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "r_s": "L",
    },
    ("black_to_white_hole", "lemaitre"): {
        "\\tau": "T", "\\rho": "L", "\\theta": "1", "\\phi": "1", "r_s": "L", "r": "L",
    },
    # The universe on a wall moving past a black hole in five dimensions: the anti-de Sitter radius
    # ell and the areal radii r, a and A are lengths, mu an area, as r^2 f is, and the curvature k,
    # the comoving radius chi, the junction factor Q and the lapse N pure numbers; w is the proper
    # distance from the wall, a length.
    ("brane_wall_universe", "static"): {
        "t": "T", "r": "L", "\\chi": "1", "\\theta": "1", "\\phi": "1", "\\ell": "L", "\\mu": "L**2", "k": "1",
        "a": "L",
    },
    ("brane_wall_universe", "gaussian_normal"): {
        "t": "T", "w": "L", "\\chi": "1", "\\theta": "1", "\\phi": "1", "\\ell": "L", "\\mu": "L**2", "k": "1",
        "a": "L", "Q": "1", "A": "L", "N": "1",
    },
    # Inside a holographic superconductor: the anti-de Sitter radius L, the radial z and the planar
    # x and y are lengths and t a time, and the blackening function f and chi are pure numbers. The
    # final Kasner chart reads its powers of the proper time tau with tau in a fixed unit, as
    # Kasner's chart does, and its t is a time; the inversion chart's rho and beta are pure numbers
    # and ell a length.
    ("holographic_superconductor", "planar"): {
        "t": "T", "z": "L", "x": "L", "y": "L", "L": "L", "f": "1", "\\chi": "1",
    },
    ("holographic_superconductor", "kasner"): {
        "\\tau": "T", "t": "T", "x": "L", "y": "L", "p_t": "1", "p_x": "1",
    },
    ("holographic_superconductor", "inversion"): {
        "\\rho": "1", "t": "T", "x": "L", "y": "L", "\\beta": "1", "\\ell": "L",
    },
    **{("stiff_fluid_universes", chart): {**{c: "1" for c in coords}, "a": "1", "b": "1", "\\ell": "L", "F": "1"}
       for chart, coords in (("areal", ("t", "x", "y", "z")), ("comoving", ("t", "\\chi", "y", "z")),
                             ("characteristic", ("u", "v", "y", "z")))},

    ("mclenaghan_tariq_tupper", "original"): {"u": "L", "r": "L", "y": "L", "z": "L", "a": "L"},
    ("mclenaghan_tariq_tupper", "exponential"): {"t": "T", "r": "L", "y": "L", "z": "L", "a": "L"},
    ("mclenaghan_tariq_tupper", "upper_half_plane"): {"t": "T", "x": "L", "y": "L", "\\phi": "1", "a": "L"},
}

# What a field carries when every coordinate is a length; the indices supply the rest.
DIMENSIONS.update({("charged_godel_universe", system): {
    "t": "T", "r": "L", "\\phi": "1", "z": "L", "\\omega": "1/L",
    **({parameter: "1/L"} if parameter else {})}
    for system, parameter in (("hyperbolic", "m"), ("linear", None), ("trigonometric", "k"))})

# Taub's vacuum under the ghost neutrinos: L a length, x and y lengths, mu and nu pure numbers.
DIMENSIONS.update({
    ("ghost_neutrinos", "taub"): {"t": "T", "x": "L", "y": "L", "z": "L", "L": "L"},
    ("ghost_neutrinos", "areal"): {"t": "T", "x": "L", "y": "L", "r": "L", "L": "L"},
    ("ghost_neutrinos", "null"): {"\\mu": "1", "\\nu": "1", "x": "L", "y": "L", "L": "L"},
})

# Bekenstein's universes: the conformal time and the comoving coordinates are pure numbers, and the
# scale factor a, named for its value in eta, and the constants l, A and a_b are lengths.
DIMENSIONS.update({("bekenstein_universes", system): {
    "\\eta": "1", radial: "1", "\\theta": "1", "\\phi": "1", "\\ell": "L", second: "L", "a": "L"}
    for system, radial, second in (("closed", "\\chi", "A"), ("flat", "r", "a_b"), ("open", "\\chi", "A"))})
# Kerr radiating and Kerr in a universe, in units c = G = 1: every coordinate but the angles is a
# length, as m, a, k and R are; A^2 = k^2 + a u makes Vaidya and Patel's a a length too, and
# m = C/A^3 makes C a length to the fourth. Vaidya's F is the logarithm of a scale factor.
DIMENSIONS.update({("vaidya_patel_kerr", chart): {
    **{c: ("1" if c in ("\\alpha", "\\beta", "\\phi") else "L") for c in coords},
    **({"a": "L", "k": "L", "C": "L**4", "A": "L", "m": "L"} if chart in ("vaidya_patel", "retarded")
       else {"m": "L", "a": "L", "R": "L"}),
    **({"q": "L"} if chart == "vaidya_patel" else {}),
    **({"\\Sigma": "L**2", "\\mu": "1/L", "\\Delta": "L**2"} if chart == "einstein_boyer_lindquist" else {}),
    **({"r": "L"} if chart.endswith("_vaidya") else {}),
    **({"F": "1"} if chart.startswith("expanding") else {})}
    for chart, coords in (("vaidya_patel", ("u", "\\alpha", "\\beta", "T")), ("retarded", ("u", "r", "\\alpha", "\\beta")),
                          ("einstein_vaidya", ("u", "\\alpha", "\\beta", "t")),
                          ("einstein_boyer_lindquist", ("t", "r", "\\alpha", "\\phi")),
                          ("expanding_vaidya", ("u", "\\alpha", "\\beta", "t")))})

DIMENSIONS.update({("charach_malin", system): {
    ("t" if system.endswith("areal") else r"\tau"): "1", r"\theta": "1", "x": "1", "y": "1",
    "L": "L", "P": "1", r"\lambda": "1", (r"\Phi" if system.startswith("scalar") else r"\chi"): "1",
} for system in ("scalar_areal", "scalar_logarithmic", "electromagnetic_areal", "electromagnetic_logarithmic")})

# Four periodic spatial lengths, cosmic time, and two constant reference times.
import chodos_detweiler as _chodos_detweiler
DIMENSIONS.update(_chodos_detweiler.DIMENSIONS)


import electron_bubble as _electron_bubble
DIMENSIONS.update(_electron_bubble.DIMENSIONS)
