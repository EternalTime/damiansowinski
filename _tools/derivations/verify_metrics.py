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
    # The wave amplitude psi and gamma sit in exponentials and are dimensionless; the null chart's
    # u = ct - rho and v = ct + rho are lengths, as the Eddington-Finkelstein times of BTZ are.
    ("einstein_rosen_waves", "cylindrical"): {
        "t": "T", "\\rho": "L", "\\phi": "1", "z": "L", "\\psi": "1", "\\gamma": "1",
    },
    ("einstein_rosen_waves", "null"): {
        "u": "L", "v": "L", "\\phi": "1", "z": "L", "\\psi": "1", "\\gamma": "1",
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
    # B = sqrt(G) B_0/c^2 folds the field in Gaussian units into an inverse length, so B rho and
    # B r sin(theta) are pure numbers; Ernst's hole keeps Schwarzschild's r_s = 2GM/c^2.
    ("melvin", "cylindrical"): {
        "t": "T", "\\rho": "L", "\\phi": "1", "z": "L", "B": "1/L",
    },
    ("melvin", "ernst"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "r_s": "L", "B": "1/L",
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
    # Fisher, Janis, Newman and Winicour's b is a length and gamma a pure number; f and h are the
    # ratios the charts name. Bronnikov's harmonic coordinate u is an inverse length, e^{-2ku} = 1 - b/r.
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
    # The slices are flat space and the whole of the geometry is the flow field carried on them,
    # so its components are velocities and the chart components of the metric are powers of V/c.
    ("natario", "cartesian_flow"): {
        "t": "T", "x": "L", "y": "L", "z": "L",
        "u": "L/T", "v": "L/T", "w": "L/T",
    },
    ("natario", "plane_flow"): {
        "t": "T", "x": "L", "y": "L", "z": "L", "u": "L/T",
    },
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
    # Visser's thin shell wormhole joins two Schwarzschild exteriors at the throat radius a, a
    # length beside r_s. Through the throat, ell = +-(r - a) is a length, |ell| carries what ell
    # carries and sgn(ell) none, and the delta at the throat an inverse length.
    ("thin_shell_wormhole", "spherical"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "r_s": "L", "a": "L",
    },
    ("thin_shell_wormhole", "throat"): {
        "t": "T", "\\ell": "L", "\\theta": "1", "\\phi": "1", "r_s": "L", "a": "L",
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
    # The mass function is folded into a length, m = GM(r)/c^2, as Oppenheimer and Volkoff
    # folded it into their u, so 1 - 2m/r is dimensionless and the redshift function, sitting
    # in an exponential, is dimensionless too.
    ("tov", "spherical"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "\\Phi": "1", "m": "L",
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
    # The same circle through Levi-Civita's mass parameter, written s here: the exponents of
    # the Kasner form along t, phi and z. Every point of the circle but (0, 0, 1), the limit of
    # large s, is reached.
    ("levi_civita", "kasner"): {
        "p_0": "2*s/(4*s**2 - 2*s + 1)",
        "p_2": "(1 - 2*s)/(4*s**2 - 2*s + 1)",
        "p_3": "2*s*(2*s - 1)/(4*s**2 - 2*s + 1)",
    },
}

# The order each small parameter of a system counts as, and the highest order the system keeps.
# See the header.
ORDERS = {
    # Hartle and Thorne's exterior: the spin is first order and the quadrupole moment second.
    ("hartle_thorne", "hartle_thorne"): ({"a": 1, "q": 2}, 2),
    # Lense and Thirring's field is linear in the source: the mass is first order, and the
    # angular momentum, which enters as the product a m, comes with it.
    ("hartle_thorne", "lense_thirring"): ({"m": 1}, 1),
}

GREEK = [
    "theta", "phi", "eta", "omega", "Omega", "ell", "pi", "lambda", "mu", "nu",
    "rho", "sigma", "tau", "chi", "psi", "alpha", "beta", "gamma", "delta",
    "epsilon", "kappa", "xi", "zeta", "Lambda", "Phi", "Theta", "Psi", "Sigma", "Delta", "Xi",
]

# The Dirac delta and its derivatives the reader reads, as \delta, \delta' and \delta''.
DIRAC_ORDERS = 3

TRIG = ["sinh", "cosh", "tanh", "coth", "sin", "cos", "tan", "cot", "sec", "csc"]

FUNCTIONS = {
    "sin": sp.sin, "cos": sp.cos, "tan": sp.tan, "cot": sp.cot,
    "sec": sp.sec, "csc": sp.csc, "sinh": sp.sinh, "cosh": sp.cosh,
    "tanh": sp.tanh, "coth": sp.coth, "exp": sp.exp, "log": sp.log,
    "ln": sp.log, "sqrt": sp.sqrt, "Abs": sp.Abs, "sign": sp.sign,
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
    kept canonical. A root above the square, w with w^q = b, as the cube root in Senovilla's
    cosh^(-2/3), is a generator of the same kind, reduced to powers below q in the numerator;
    below the line it may stand only as a factor of its own, 1/w = w^(q-1)/b, and a sum that
    holds one there is left to the general simplifier. Such roots are taken as independent of
    one another, which the roots of distinct irreducible factors are; were two of them not,
    a value that vanishes could be left standing, and none that does not could be lost.
    """
    for function, rewrite in _IN_SIN_COS_EXP.items():
        if expression.has(function):
            expression = expression.replace(function, rewrite)
    if expression.has(sp.Abs, sp.sign, sp.DiracDelta):
        expression = _on_a_kink(expression)
    if any(not atom.args[0].is_Symbol for atom in expression.atoms(sp.sin, sp.cos)):
        expression = sp.expand_trig(expression)
    # A logarithm's argument is put in one form, since expand spreads (x^2 + y^2)/rho_0^2 into two
    # fractions and the logarithm of each spelling would otherwise be a generator of its own.
    if expression.has(sp.log):
        expression = expression.replace(lambda x: isinstance(x, sp.log), lambda x: sp.log(sp.factor(x.args[0])))
    # The reader writes a mixed partial in the order the coordinates are listed and diff
    # writes it in its own; doit puts both in diff's, so they become one generator.
    if expression.has(sp.Derivative):
        expression = expression.replace(lambda x: isinstance(x, sp.Derivative), lambda x: x.doit())
    # Innermost first, so a radicand is itself in this form before it is factored.
    expression = expression.replace(
        lambda x: x.is_Pow and x.exp.is_Rational and not x.exp.is_Integer and not x.base.is_Number,
        _factored_root)
    # A power of a sum or a product with a symbolic exponent, as (1 - 2m/r)^q, is split the same way,
    # and so is one of a reciprocal, as the (1/ell)^(2b) sympy makes of (rho/ell)^(2b - 2).
    expression = expression.replace(
        lambda x: x.func is sp.Pow and not x.exp.is_Number
        and (x.base.is_Add or x.base.is_Mul or (x.base.is_Pow and x.base.exp.is_Integer)),
        _factored_root)

    generators = set()
    radicals = {}
    _collect_generators(expression, generators, radicals)
    algebraic = []
    for cosine in [g for g in generators if g.func is sp.cos]:
        if sp.sin(cosine.args[0]) in generators:
            algebraic.append((cosine, 1 - sp.sin(cosine.args[0]) ** 2))
    # sgn(x)^2 = 1 wherever x is not zero, which is where a function of x is read.
    for sign in [g for g in generators if g.func is sp.sign]:
        algebraic.append((sign, sp.Integer(1)))
    for base, q in radicals.items():
        if base == -1 and q != 2:
            raise NotImplementedError("a root of -1 above the square")
        root = sp.I if base == -1 else sp.root(base, q)
        if not (root is sp.I or (root.is_Pow and (q == 2 or root.base == base))):
            raise NotImplementedError(f"the root {q} of {base} does not stay a power")
        generators.add(root)
        algebraic.append((root, base, q))
    # A radicand that holds another algebraic generator is reduced before that one is,
    # so the conjugates it leaves behind can still be reduced on the inner one.
    roots = [pair[0] for pair in algebraic]
    algebraic.sort(key=lambda pair: -sum(1 for root in roots if sp.sympify(pair[1]).has(root)))

    ring = PolyRing(sorted(generators, key=sp.default_sort_key), sp.QQ, "lex")
    index = {symbol: i for i, symbol in enumerate(ring.symbols)}
    reader = _FractionReader(ring, index, radicals)
    value = reader(expression)
    relations = [(index[pair[0]], reader(pair[1]), pair[2] if len(pair) == 3 else 2) for pair in algebraic]
    for _ in range(2 * len(relations) + 2):
        if not any(value.holds(i, q) for i, _, q in relations):
            break
        for i, base, q in relations:
            if value.holds(i, q):
                value = value.reduce(i, base, q)
    else:
        raise NotImplementedError("the algebraic generators do not settle")
    value = value.cancelled()
    if not value.numerator:
        return sp.Integer(0)
    return value.as_expr()


def _on_a_kink(expression):
    """A metric with a kink, as the domain wall's (1 - k|z|)^2, read as a distribution.

    |x| of a coordinate is written x sgn(x), so that diff gives its derivative sgn(x) and
    the derivative of that, 2 delta(x), and sgn(x)^2 = 1 is one of _canonical's relations.
    Such a metric is continuous with a jump in its first derivative, so its curvature is a
    bounded function plus a delta of x times a function continuous at x = 0, which is the
    value of that function at x = 0 times the delta; that is the reduction made here, term
    by term in the power of delta(x), with sgn(x)^2 taken as 1 and x sgn(x) as |x|. A delta
    times a function that jumps at x = 0, which no reading as a distribution fixes, stops
    the reader. A square of a delta, which only the Kretschmann scalar holds, is left as it
    is; `off_support` compares such a scalar where x is not zero, as the entry states it.
    A delta whose argument is not a bare coordinate, or which carries a prime, is passed
    over, as the Aichelburg-Sexl shock's are, since nothing multiplying it varies across it.
    """
    expression = expression.replace(lambda e: isinstance(e, sp.Abs) and e.args[0].is_Symbol,
                                    lambda e: e.args[0] * sp.sign(e.args[0]))
    for delta in sorted(expression.atoms(sp.DiracDelta), key=sp.default_sort_key):
        if len(delta.args) != 1 or not delta.args[0].is_Symbol:
            continue
        x = delta.args[0]
        D, s = sp.Dummy("D"), sp.Dummy("s")
        e = expression.xreplace({delta: D, sp.sign(x): s})
        if not e.has(x) and not e.has(s):
            continue
        numerator, denominator = sp.fraction(sp.together(e))
        if denominator.has(D):
            raise NotImplementedError(f"{delta} stands below a fraction line")
        terms = []
        for (power,), coefficient in sp.Poly(sp.expand(numerator), D).as_dict().items():
            coefficient = coefficient if hasattr(coefficient, "free_symbols") else sp.sympify(coefficient)
            if power != 1:
                terms.append(D ** power * coefficient / denominator)
                continue
            on = [_even_in(sp.expand(side.subs(x, 0)), s) for side in (coefficient, denominator)]
            if any(side.has(s) for side in on) or on[1] == 0:
                raise NotImplementedError(f"{delta} multiplies a function that jumps at {x} = 0")
            terms.append(D * on[0] / on[1])
        expression = sp.Add(*terms).xreplace({D: delta, s: sp.sign(x)})
    return expression


def _even_in(expression, s):
    """The expression with every power of s, a sign, reduced by s^2 = 1."""
    return sp.expand(expression.replace(lambda e: e.is_Pow and e.base == s and e.exp.is_Integer,
                                        lambda e: s ** (int(e.exp) % 2)))


def _factored_root(power):
    """b^e as the product of f^(k e) over the irreducible factors f^k of b, for e a fraction
    or an exponent that holds a parameter, as the q of (1 - 2m/r)^q.

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
            if _is_transcendental_power(generator):
                generators.add(generator)
            else:
                # The whole part of the exponent leaves the base itself, which is read as it stands.
                _collect_generators(generator, generators, radicals)
    else:
        generators.add(expression)


class _FractionReader:
    """An expression, read into a _Fraction over the given ring."""

    def __init__(self, ring, index, radicals=None):
        self.ring = ring
        self.index = index
        self.radicals = radicals or {}
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
            q = self.radicals.get(expression.base, 2)
            root = sp.I if expression.base == -1 else sp.root(expression.base, q)
            return self.generator(root) ** int(expression.exp * q)
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

    def holds(self, i, q=2):
        """Whether generator i, a root of order q, is still anywhere it should not be."""
        return (self.numerator.degree(i) > q - 1
                or any(factor.degree(i) > 0 for factor in self.denominator))

    def reduce(self, i, base, q=2):
        """Fold generator i, w with w^q = base, out of every denominator factor and down
        to powers below q in the numerator."""
        if q != 2:
            return self._reduce_root(i, base, q)
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

    def _reduce_root(self, i, base, q):
        """The same for a root above the square, which may stand below the line only as a
        factor of its own: 1/w = w^(q-1)/base."""
        w = self.reader.ring.gens[i]
        free = {f: k for f, k in self.denominator.items() if f.degree(i) <= 0}
        out = _Fraction(self.reader, self.reader.ring.one, free)
        for factor, multiplicity in self.denominator.items():
            if factor.degree(i) <= 0:
                continue
            if factor != w:
                raise NotImplementedError(f"a root above the square stands in the sum {factor.as_expr()} "
                                          "below a fraction line")
            out = out * (_Fraction(self.reader, w ** (q - 1), {}) * base.inverse()) ** multiplicity
        return out * _fold(self.reader, self.numerator, i, base, q)

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


def _fold(reader, polynomial, i, base, q=2):
    """polynomial with every w^k, w the generator i, replaced by w^(k mod q) base^(k div q)."""
    by_degree = {}
    for monomial, coefficient in polynomial.terms():
        free = monomial[:i] + (0,) + monomial[i + 1:]
        by_degree.setdefault(monomial[i], {})[free] = coefficient
    terms = []
    w = reader.ring.gens[i]
    for k, part in by_degree.items():
        term = _Fraction(reader, polynomial.ring.from_dict(part) * w ** (k % q), {})
        terms.append(term * base ** (k // q) if k >= q else term)
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


def powered_trig_calls(text, names):
    """Turn trig ** n (argument) into (trig(argument))**n, the argument a balanced bracket,
    which is how \\cosh^2\\left(\\sqrt{\\Lambda}\\,t\\right) reads once the brackets are bare.
    The power may be negative or a fraction, as Senovilla's \\cosh^{-2/3}(3a\\rho) is."""
    pattern = re.compile(r"\b(" + names + r")\s*\*\*\s*\(?\s*(-?\d+(?:\s*/\s*\d+)?)\s*\)?\s*\(")
    while True:
        found = pattern.search(text)
        if not found:
            return text
        opening = found.end() - 1
        depth = 0
        for closing in range(opening, len(text)):
            depth += {"(": 1, ")": -1}.get(text[closing], 0)
            if depth == 0:
                break
        else:
            raise LatexError(f"unbalanced bracket after {found.group(1)} in {text!r}")
        body = text[opening + 1:closing]
        text = f"{text[:found.start()]}({found.group(1)}({body}))**({found.group(2)}){text[closing + 1:]}"


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

    def __init__(self, coords, parameters, time_coords=frozenset(), relations=None, kept=None):
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
        # A name the entry defines, as Weyl's R = \sqrt{\rho^2 + z^2}: read as the expression it names.
        self.defined = {}
        self.primed = set()
        # A parameter spelled with a command and a subscript, as \chi_0 is, is read whole: the
        # Greek letters are turned into words below, and \chi_0 would otherwise be read as chi
        # times a stray _0, which is zero. The subscript may be a command too, as in \Delta_\theta.
        self.spelled = {}
        for declaration in parameters:
            self._declare_parameter(declaration)
            spelling = declaration.split("=")[0].strip()
            if re.fullmatch(r"\\[A-Za-z]+_\\?\w+", spelling):
                self.spelled[spelling] = self._plain(spelling)
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
        # \delta(u) is the Dirac delta, and each prime on it one derivative with respect to its
        # argument, unless the system names a parameter \delta of its own, as the cosmic
        # string names its deficit angle.
        self.dirac = "delta" not in self.parameters
        if self.dirac:
            for order in range(DIRAC_ORDERS):
                self.local[f"DIRAC{order}"] = lambda argument, k=order: sp.DiracDelta(argument, k)
        self.known = set(self.local) | ({KEYWORD_LAMBDA} if "lambda" in self.local else set())
        self.transforms = standard_transformations + (
            split_symbols_custom(lambda name, _=None: name not in self.known),
            implicit_multiplication,
        )
        # A defined name whose definition holds a function of the coordinates, as Szekeres's E
        # holds S(r), P(r) and Q(r), is held as a function of the coordinates it varies with
        # while the tensors are built, and written out by surface() where two values are compared.
        self.held = {}
        for plain, definition in self.defined.items():
            value = self.defined[plain] = self(definition)
            if value.atoms(sp.core.function.AppliedUndef):
                names = [name for name in self.coords if value.has(self.symbol[name])]
                function = sp.Function(plain, real=True)(*(self.symbol[name] for name in names))
                self.functions[plain] = (function, names)
                self.held[function] = value
                value = function
            self.parameters[plain] = self.local[plain] = value
            self.known.add(plain)
        self.relations = {}
        for name, value in (relations or {}).items():
            if name not in self.parameters:
                raise LatexError(f"a relation is declared for {name!r}, which is not a parameter")
            self.relations[self.parameters[name]] = sp.sympify(value)
        # A system kept to an order counts each of its small parameters as a power of one symbol.
        self.order = None
        if kept:
            weights, highest = kept
            stray = [name for name in weights if name not in self.parameters or not self.parameters[name].is_Symbol]
            if stray:
                raise LatexError(f"an order is declared for {stray}, which are not constants of the system")
            small = sp.Dummy("epsilon", positive=True)
            self.order = ({self.parameters[name]: small ** weight * self.parameters[name]
                           for name, weight in weights.items()}, small, highest)

    def surface(self, expression):
        """The expression on the surface the entry's constrained parameters live on, and to
        the order the system keeps."""
        expression = expression.subs(self.relations)
        expression = expression.subs(self.held).doit() if self.held else expression
        return self.truncated(expression) if self.order else expression

    def truncated(self, expression):
        """The Taylor polynomial of the expression through the order the system keeps, in
        canonical form: each small parameter is scaled by its power of one symbol, and the
        polynomial in that symbol is cut. A matrix is cut entry by entry. The expression is
        differentiated as it stands and put in canonical form only once the symbol is set to
        zero, where its denominators are those of the zeroth order: the canonical form of the
        whole would factor every denominator, and one that holds the small parameters, as
        Hartle and Thorne's 1 + 2h_2 P_2 does, does not factor in ten minutes."""
        if self.order is None:
            return norm(expression)
        if isinstance(expression, sp.MatrixBase):
            return expression.applyfunc(self.truncated)
        scaled, small, highest = self.order
        term = sp.sympify(expression).xreplace(scaled)
        out = sp.Integer(0)
        for k in range(highest + 1):
            out += norm(term.subs(small, 0)) / sp.factorial(k)
            term = sp.diff(term, small)
        return norm(out)

    def zeroth(self, g):
        """A metric at zeroth order, every small parameter set to zero."""
        scaled, small, _ = self.order
        return norm(g.applyfunc(lambda e: e.xreplace(scaled).subs(small, 0)))

    def inverse(self, g):
        """The inverse of a metric to the order the system keeps, as the series
        g_0^-1 - g_0^-1 d g_0^-1 + ..., with g_0 the metric at zeroth order and d the rest of
        its polynomial, each term cut as it is formed."""
        highest = self.order[2]
        g = self.truncated(g)
        zeroth = self.zeroth(g)
        zeroth_inverse = norm(zeroth.adjugate() / zeroth.det())
        rest = norm(g - zeroth)
        out = term = zeroth_inverse
        for _ in range(highest):
            term = self.truncated(-term * rest * zeroth_inverse)
            out += term
        return norm(out)

    @staticmethod
    def _plain(name):
        return name.replace("\\", "").strip()

    def _declare_parameter(self, declaration):
        """`a` is a constant; `a = a(t)` is a function of the coordinates it names; and
        `R = \\sqrt{\\rho^2 + z^2}` is a name for the expression on its right, in the
        coordinates and the parameters declared before it, which every published value
        that writes the name is read as.

        For the second kind the printed rate is a derivative with respect to the chart
        coordinate, which is c times the named one when that one is a time, so the rate
        carries the matching power of c. A function of one coordinate answers to a dot
        and to a prime, a dot up to the second derivative and a prime up to the third,
        and a function of any number of coordinates answers to \\partial at any order,
        which _declare_partials resolves when a published value names one.
        """
        name, _, right = declaration.partition("=")
        plain = self._plain(name)
        self.parameter_names.append(plain)
        if right.strip() and not re.fullmatch(re.escape(name.strip()) + r"\s*\([^()]*\)", right.strip()):
            self.defined[plain] = right
            return
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
        # A third prime is as far as a chart goes: the travelling wave on a string writes its
        # g_uu with A'', and a Christoffel symbol differentiates that once more.
        third = sp.Derivative(function, variable, 3) / scale ** 3
        # Both spellings appear: a dot in the comoving chart, a prime in the conformal one.
        for suffix, value in (("_dot", first), ("_prime", first),
                              ("_ddot", second), ("_pprime", second), ("_ppprime", third)):
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
        for spelling in sorted(self.spelled, key=len, reverse=True):
            text = text.replace(spelling, f" {self.spelled[spelling]} ")
        text = re.sub(r"\\[,;:!>]", " ", text)
        text = re.sub(r"\\ ", " ", text)
        text = text.replace("\\cdot", "*")
        # |z| is the absolute value of a coordinate and \mathrm{sgn}(z) its sign, as a metric
        # with a kink writes them; _on_a_kink says how they are read.
        text = re.sub(r"\|\s*(\\?[A-Za-z]+)\s*\|", r" Abs(\1) ", text)
        text = text.replace("\\mathrm{sgn}", " sign ")
        if self.dirac:
            text = re.sub(r"\\delta\s*('*)\s*\(", lambda m: f" DIRAC{len(m.group(1))}(", text)
        for name in sorted(self.primed, key=len, reverse=True):
            text = re.sub(re.escape(name) + r"'''", f" {name}_ppprime ", text)
            text = re.sub(re.escape(name) + r"''", f" {name}_pprime ", text)
            text = re.sub(re.escape(name) + r"'", f" {name}_prime ", text)
        text = expand_partials(text)
        text = expand_fractions(text)
        text = expand_braced_call(text, "sqrt", "sqrt")
        text = expand_braced_call(text, "ddot", "DDOT")
        text = expand_braced_call(text, "dot", "DOT")
        # A dotted name may be accented, as \dot{\tilde\phi} is, and reads as the name tildephi, or
        # carry a subscript, as the Kaluza-Klein monopole's \dot{x_5}.
        text = re.sub(r"DDOT\s*\(\s*((?:\\?[A-Za-z]+)+(?:_[A-Za-z0-9]+)?)\s*\)",
                      lambda m: f" {self._plain(m.group(1))}_ddot ", text)
        text = re.sub(r"DOT\s*\(\s*((?:\\?[A-Za-z]+)+(?:_[A-Za-z0-9]+)?)\s*\)",
                      lambda m: f" {self._plain(m.group(1))}_dot ", text)
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
        # The same with a bracketed argument, as \cosh^2\left(\sqrt{\Lambda}\,t\right).
        text = powered_trig_calls(text, names)
        text = re.sub(r"\b(" + names + r")\s+([A-Za-z]\w*)", r"\1(\2)", text)
        # A power of e is the exponential, unless the system declares e as a name of its own, as
        # Ori's time machine does its parameter e, whose square is then that parameter's.
        if "e" not in self.local:
            text = re.sub(r"(?<![A-Za-z_])e\s*\*\*", " E**", text)
        # 3R\sqrt{..} leaves 3Rsqrt(..), whose one name token would be split letter by letter.
        for name in FUNCTIONS:
            text = re.sub(r"(?<=[A-Za-z0-9_])(" + name + r")\s*\(", r" \1(", text)
        # Python reads 2J as the imaginary number 2j, so a number is parted from a J it multiplies.
        text = re.sub(r"(?<![A-Za-z_])(\d+)([jJ])", r"\1 \2", text)
        if "\\" in text:
            raise LatexError(f"unhandled LaTeX in {latex!r}: {text!r}")
        return re.sub(r"[A-Za-z][A-Za-z0-9_]*", self._part_subscripted, text)

    def _part_subscripted(self, match):
        """A name written against the letters before it, as the b_0 of 2ab_0^2, parted from them.

        The parser splits a token it does not know letter by letter, and a subscript split so
        leaves a stray _0, which multiplies the whole term by zero without a word. So a token
        that is not a name the system knows and ends in one that carries a subscript is read as
        what stands before it times that name, and any other unknown token with an underscore
        is an error rather than a zero."""
        token = match.group(0)
        if "_" not in token or token in self.local or PARTIAL_NAME.fullmatch(token):
            return token
        for start in range(1, len(token)):
            if "_" in token[start:] and token[start:] in self.local and "_" not in token[:start]:
                return f"{token[:start]} {token[start:]}"
        raise LatexError(f"{match.string!r} holds {token!r}, a subscripted name the system does not declare")

    def __call__(self, latex):
        text = self._preprocess(latex)
        self._declare_partials(text)
        local = dict(self.local)
        if "lambda" in local:
            # A system may name a function \lambda, as Gowdy's cosmologies do, and Python keeps
            # that word for itself, so the parser is handed it under a spelling of its own.
            text = re.sub(r"(?<![A-Za-z0-9_])lambda(?![A-Za-z0-9_])", KEYWORD_LAMBDA, text)
            local[KEYWORD_LAMBDA] = local.pop("lambda")
        try:
            expression = parse_expr(
                text, local_dict=local, transformations=self.transforms, evaluate=True
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
            if name in reader.defined:
                if name in reader.functions:
                    self.of_function[name] = dimension
                continue
            if value.is_Symbol:
                self.of_symbol[value] = dimension
            else:
                self.of_function[name] = dimension
        # A defined name carries the dimension of its definition, which has to be the declared one.
        for name, value in reader.defined.items():
            dimension = sp.sympify(declared[name], locals=BASE_DIMENSIONS)
            if sp.simplify(self(value) / dimension) != 1:
                raise DimensionError(f"{name} is declared {dimension} and its definition carries {self(value)}")

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
        if isinstance(expression, sp.Abs):
            return self(expression.args[0])
        if isinstance(expression, sp.sign):
            self(expression.args[0])
            return sp.Integer(1)
        if isinstance(expression, sp.DiracDelta):
            # The delta carries the inverse of its argument, and each derivative one more.
            order = expression.args[1] if len(expression.args) > 1 else 0
            return self(expression.args[0]) ** -(order + 1)
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
        text, _ = scalar_parts(entry[field])
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


def _adjugate(matrix):
    """The adjugate and the determinant of a square matrix, by cofactors, every minor put
    through norm as it is formed and kept, so that no step holds an unreduced product.

    sympy's own adjugate() and det() reduce each product through cancel, which on a metric with
    three cross terms, as Kerr-de Sitter's Kerr-Schild chart has, did not finish in ten minutes;
    on a diagonal metric the two agree at once."""
    n = matrix.shape[0]
    minors = {}

    def minor(rows, columns):
        """The determinant of the submatrix of the given rows and columns."""
        if not rows:
            return sp.Integer(1)
        key = (rows, columns)
        if key not in minors:
            total = sp.Integer(0)
            for k, column in enumerate(columns):
                entry = matrix[rows[0], column]
                if entry != 0:
                    total += (-1) ** k * entry * minor(rows[1:], columns[:k] + columns[k + 1:])
            minors[key] = norm(total)
        return minors[key]

    every = tuple(range(n))
    cofactors = sp.Matrix(n, n, lambda i, j: (-1) ** (i + j) * minor(every[:j] + every[j + 1:], every[:i] + every[i + 1:]))
    return cofactors, minor(every, every)


class Geometry:
    """Every tensor the files publish, computed from g in the chart the coords name."""

    def __init__(self, g, coords, seconds, order=None):
        self.n = len(coords)
        self.coords = coords
        self.seconds = seconds
        # A system kept to an order passes `order`, its Reader: the metric and the inverse are
        # then their polynomials to that order, and every tensor is cut after each product.
        self.settle = order.truncated if order else norm
        self.g = self.settle(norm(g))
        # The adjugate over the determinant, which sympy's own inv() takes far longer to reach.
        if order:
            self.ginv = order.inverse(self.g)
        else:
            adjugate, determinant = _adjugate(self.g)
            self.ginv = norm(adjugate / determinant)
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
                        value = self.settle(sum(
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
                        out[mu][nu][rho] = self.settle(
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
                            value = self.settle(
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
                            out[mu][nu][rho][sigma] = self.settle(sum(
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
                    value = self.settle(sum(riemann[mu][nu][mu][sigma] for mu in range(self.n)))
                    out[nu][sigma] = value
                    out[sigma][nu] = value
            return out
        return self._timed("ricci_ll", build)

    def ricci_scalar(self):
        def build():
            ricci = self.ricci_ll()
            return self.settle(sum(
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
                    out[a][b] = self.settle(ricci[a][b] - scalar * self.g[a, b] / 2)
            return out
        return self._timed("einstein_ll", build)

    def kretschmann(self):
        def build():
            lower = self.riemann_llll()
            # One index at a time, so each component is a sum of n terms rather than n^4.
            upper = self.raise_indices(lower, 4, (0, 1, 2, 3))
            return self.settle(sum(
                lower[a][b][cc][d] * upper[a][b][cc][d]
                for a in range(self.n) for b in range(self.n)
                for cc in range(self.n) for d in range(self.n)
            ))
        return self._timed("kretschmann", build)

    def weyl_llll(self):
        def build():
            n = self.n
            if n < 3:
                # In two dimensions the Riemann tensor is (R/2)(g g - g g), all trace, and the Weyl
                # tensor, what is left of it with the traces removed, vanishes identically.
                return self._zeros(4)
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
                            out[a][b][cc][d] = self.settle(riemann[a][b][cc][d] - correction + trace)
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
                _put(out, index, self.settle(total))
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


def beyond_order(reader, value):
    """Whether a published value of a system kept to an order carries a term beyond it."""
    return reader.order is not None and norm(value - reader.truncated(value)) != 0


def compare_block(report, reader, where, published, computed, variance, coords, time_coords, c, cut=True):
    """Every published component against sympy, and every omitted one against zero. In a
    system kept to an order a value has to be cut at it as well, unless `cut` is off, as it
    is for the metric and its inverse, which may stand as the line element has them."""
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
            value = reader(entry["value"])
        except LatexError as error:
            report.skip(f"{where} {names}", str(error))
            continue
        if cut and beyond_order(reader, value):
            report.disagree(where, f"{names} published as {entry['value']} carries a term beyond the order kept")
        value = reader.surface(value)
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


WHERE = re.compile(r"\s*\\;\\text\{for\}\\;\s*(\\?[A-Za-z]+)\s*\\neq\s*0\s*$")


def scalar_parts(published):
    """A published scalar as (its value, the coordinate it holds away from the zero of, or None).

    The domain wall's Kretschmann scalar is written "K = 0 \\;\\text{for}\\; z \\neq 0": its
    Riemann tensor carries delta(z), and the square of a delta is no distribution, so the
    scalar is stated where z is not zero, and compared there by off_support.
    """
    text = published.split("=", 1)[1] if "=" in published.split("\\;")[0] else published
    found = WHERE.search(text)
    if found is None:
        return text, None
    return text[:found.start()], found.group(1)


def off_support(expression, x):
    """The expression where the coordinate x is not zero: every delta of x set to zero."""
    return expression.replace(lambda e: isinstance(e, sp.DiracDelta) and e.args[0].has(x), lambda e: sp.Integer(0))


def compare_scalar(report, reader, where, published, computed):
    text, away = scalar_parts(published)
    try:
        value = reader(text)
    except LatexError as error:
        report.skip(where, str(error))
        return
    if beyond_order(reader, value):
        report.disagree(where, f"published as {published.strip()} carries a term beyond the order kept")
    value = reader.surface(value)
    computed = reader.surface(computed)
    if away is not None:
        if away not in reader.symbol:
            report.disagree(where, f"{published!r} holds away from {away}, which is not a coordinate")
            return
        computed = off_support(computed, reader.symbol[away])
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
            residual = sp.expand(reader(left) - reader(right))
        except LatexError as error:
            report.skip(f"{where} {equation!r}", str(error))
            continue
        if beyond_order(reader, residual):
            report.disagree(where, f"{equation!r} carries a term beyond the order kept")
        residual = reader.surface(residual)
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
        reader = Reader(coords, parameters, declaration, relations, ORDERS.get((metric_id, entry["id"])))
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
    # A system kept to an order is a series about its zeroth order, which is what has to be a metric.
    if _adjugate(reader.zeroth(g) if reader.order else g)[1] == 0:
        report.skip(where, "the line element gives a degenerate metric")
        return

    # Everything is computed with the bare coordinates and then weighted into the
    # x^0 = cT chart by compare_block, which is exact because the rescaling is linear.
    symbols = [reader.symbol[name] for name in coords]
    geometry = Geometry(g, symbols, seconds, reader if reader.order else None)
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
            report, reader, label, published, as_lists, variance, coords, declaration, c, cut=False))
    if reader.truncated(g * geometry.ginv) != sp.eye(len(coords)):
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
