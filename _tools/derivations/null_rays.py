#!/usr/bin/env python3
"""Spacetime diagrams of null rays and future light cones, drawn from the published metrics.

Every diagram is computed from a coordinate system's published metric_components,
inverse_metric_components, kretschmann and domains, read through the Reader of
verify_metrics.py beside this file, in the x^0 = cT chart the collection writes, with
c = 1. Nothing is drawn by eye: the rays are integrated, the cones are the null
directions at a point, and every marker is the zero set of a published quantity.

    python3 -m venv /tmp/mfs-venv && /tmp/mfs-venv/bin/pip install sympy numpy scipy contourpy
    /tmp/mfs-venv/bin/python _tools/derivations/null_rays.py
    python3 _tools/build_mfs_data.py

writes MFS/assets/data/diagrams/<metric_id>.json, one file per spacetime that has a
diagram, and the second command stamps each file's version into the index. Pass
--metric <metric_id> to redraw one spacetime, repeatable, and --verify to check the
rays against closed forms the drawing never uses instead of writing anything.


The null condition
------------------

On a plane of the chart's time x^0 and one spatial coordinate r, every other coordinate
held fixed, the metric on that plane is the 2 by 2 matrix h = [[g_00, g_0r], [g_0r, g_rr]],
and a null direction (dx^0, dr) solves

    g_00 (dx^0)^2 + 2 g_0r dx^0 dr + g_rr dr^2 = 0,

which has two real roots exactly where D = g_0r^2 - g_00 g_rr > 0. They are written as

    P ~ (-g_0r + sqrt(D), g_00) ~ (g_rr, -g_0r - sqrt(D))
    M ~ (-g_0r - sqrt(D), g_00) ~ (g_rr, -g_0r + sqrt(D))

taking at each point whichever form does not vanish, which is what lets one formula
serve a diagonal h, the Eddington-Finkelstein charts' g_rr = 0 and a double null
chart's g_uu = g_vv = 0. h is divided by its largest entry first, which changes
no direction and keeps the numbers finite next to a horizon. P and M are continuous line
fields: P is the family that conserves an advanced coordinate, drawn as the ingoing or
left moving family, and M the retarded one. --verify confirms that labelling chart by
chart. Their integral curves are traced both ways by fourth order Runge-Kutta in the
drawing's own unit square, from seeds spaced evenly along its four edges.

Every curve drawn is a null curve, its tangent null. It is also a null geodesic, the path
a light ray takes, only if nothing accelerates it out of the plane, Gamma^A_ab k^a k^b = 0
for every fixed coordinate A. The planes in DIAGRAMS pass that test, except Godel's, whose
caption says its null curves are not null geodesics.


Rays that leave the plane
-------------------------

Off its axis, Kerr's light rays that run straight in and straight out stay in no plane of
two coordinates: they turn in phi as they go. They are its principal null congruence, and
a row with principal=True draws it from the published Weyl tensor rather than from a
formula for it.

At each point the published C_abcd, taken in a frame the published metric makes
orthonormal, is an operator C^AB_CD on the six bivectors, its entries the size of its
eigenvalues. Where the Weyl tensor is of Petrov type D, the two eigenvalues of largest
size are -2 times the other four, which fall into two equal pairs, and the eigenbivectors
of those two span a pair of simple planes, one timelike and one spacelike. The timelike
one is the principal plane, and its two null directions are the repeated principal null
directions, the k with C_abc[d k_e] k^b k^c = 0.

The principal plane is written in the basis U_0, U_r whose shadow on the drawn pair of
coordinates is their coordinate basis, and its metric h_ij = g(U_i, U_j) takes the place
of the coordinate plane's metric in the null condition above: the same quadratic, the same
P and M, the same tracing and the same cones, which are now the future cone of the
principal plane, whose edges are the two principal null directions. A ray traced so is the
shadow on the drawn pair of a ray that also moves in the coordinates the row names in
`leaves`. That is exact when the published metric does not depend on them, since the
principal plane then does not either, and the script requires it: every coordinate is
drawn, held fixed or left, and neither the published metric nor the Weyl tensor may depend
on one that is left. It refuses to draw, naming the point, where the Weyl tensor is not of
type D, where the principal plane has a component along a coordinate held fixed, so that
its rays would leave the surface drawn, or where the plane does not project one to one
onto the drawn pair. Beside a curvature singularity the published components lose to
rounding the digits the plane is read from, so it is not taken where the Kretschmann
scalar passes K_END, and a ray stops there.

The same rays seen from above are the pair (phi, r) drawn as X = r cos phi and
Y = r sin phi, which to_display=POLAR asks for. Their seeds are spaced evenly around the
circle inscribed in the box, so the rays of a family are copies of one another turned
about the axis, as the independence of phi makes them, and with inside=True a ray stops a
thousandth of the drawing short of the edge of the published domain, where phi winds
without end.

Before a principal view is written, its rays are checked to be null geodesics against the
published Christoffel symbols, k^b nabla_b k^a parallel to k^a, and its directions against
C_abc[d k_e] k^b k^c = 0 with the published Weyl tensor, and both fields are stamped.
--verify adds the closed forms of Kerr and Kerr-Newman, which the drawing never uses, and
shows that in Schwarzschild, Reissner-Nordstrom and Taub-NUT the principal plane is the
plane of t and r their radial views already draw, and that van Stockum's Weyl tensor is of
type I, with no principal congruence to draw.


Which way is the future
-----------------------

Every cone is a future cone, oriented by the rule the table names:

  'tau'       a time function: future where k.dtau > 0, with dtau timelike by the
              published inverse metric.
  'ingoing'   the ingoing family future directed toward smaller r everywhere. It agrees
              with t on the side a static chart is written for, and it carries the
              orientation through horizons the chart's t cannot, reading the region inside
              r_s as the black hole, as an ingoing Eddington-Finkelstein chart does.
  'outgoing'  the outgoing family future directed toward larger r or x.
  'vector'    the chart's time direction: future where g(k, d_0) < 0. Godel needs it,
              since its published g^tt is positive and no time function exists.

The other edge of a cone is the null direction whose product with the first is negative.


Markers
-------

  grr         g^rr = 0 (or g^xx), where r = const is null or the chart ends.
  g00         g^00 = 0, where x^0 = const stops being spacelike, on request.
  gtt         g_tt = 0 for the chart's own time t, on request: where d_t turns null, which
              for Kerr and Kerr-Newman is the ergosurface. It names itself in the legend.
  throat      d_r R = 0, with R^2 = g_theta_theta the areal radius.
  apparent    |grad R|^2 = 0 where grad R is not zero: marginally trapped spheres.
  singular    an edge along which the published Kretschmann scalar exceeds 1e8 and grows
              at least fiftyfold between 1e-4 and 1e-5 of the drawing, where h is
              Lorentzian.
  hatch       outside the entry's published domains, parsed by the same Reader; a domain
              ending at the undeclared r_+ is read at the outermost zero of g^rr.


Declared inputs
---------------

An entry that leaves a function free cannot be drawn without a choice. The choice is
written in the table and printed beside the diagram: either an expression for the
function, or dust, whose scale factors are solved from the entry's own published
G^i_i = 0. Parameter values, plot ranges and fixed coordinates are the other choices.


Output
------

Each view records the published fields it was drawn from and their version, as
build_mfs_data.diagram_source_version computes it; build_mfs_data.py --check fails when
a metric changes without its diagram being redrawn. It carries its ticks, and every
label as TeX in $...$, so the page and the application draw and set the same ones.
to_display is the linear map's rows, or "polar". A view whose cones are not the future
light cone says what they are in `cone`, and a marker that names itself carries `legend`.
"""

import argparse
import json
import math
import sys
from dataclasses import dataclass, field, replace
from decimal import Decimal
from fractions import Fraction
from pathlib import Path

import contourpy
import numpy as np
import sympy as sp
from scipy.integrate import quad, solve_ivp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import build_mfs_data as build  # noqa: E402
import verify_metrics as vm  # noqa: E402

DIAGRAMS_DIR = build.DIAGRAMS_DIR
BASE_FIELDS = ["coords", "parameters", "metric_components", "inverse_metric_components",
               "kretschmann", "domains"]

FINKELSTEIN_IN = ((0, 1), (1, -1))      # X = r, Y = v - r
FINKELSTEIN_OUT = ((0, 1), (1, 1))      # X = r, Y = u + r
NULL_TO_TR = ((-0.5, 0.5), (0.5, 0.5))  # (u, v) of Minkowski to X = (v - u)/2, Y = (u + v)/2
UV_TO_TZ = ((-0.5, 1), (0.5, 1))        # (u, v) of the plane wave to z = v - u/2, t = v + u/2
RADIAL = ("ingoing", "outgoing")
SIDEWAYS = ("moving left", "moving right")
EQUATOR = {"theta": "pi/2", "phi": "0"}
POLAR = "polar"                         # (phi, r) drawn from above: X = r cos phi, Y = r sin phi
PRINCIPAL_CONE = "future cone of the principal plane"


@dataclass
class Diagram:
    """One view of one coordinate system, and every choice its drawing makes."""

    metric: str
    system: str
    view: str                       # unique within the system
    label: str                      # the view's name on the page; every label is text with its
                                    # mathematics in $...$, as a caption is, so it is set as TeX
    plane: tuple                    # (time coordinate, spatial coordinate), as the entry spells them
    box: tuple                      # (Xmin, Xmax, Ymin, Ymax) in the drawn axes
    xlabel: str                     # the axes, as the line element spells their symbols
    ylabel: str
    params: dict = field(default_factory=dict)      # parameter -> value
    fixed: dict = field(default_factory=dict)       # every other coordinate -> value
    to_display: tuple = ((0, 1), (1, 0))            # rows: X = a.(x^0, r), Y = b.(x^0, r), or POLAR
    mirror: bool = False            # a line through the centre: x = r and its reflection x = -r
    orient: str = "tau"             # tau, ingoing, outgoing or vector; see the header
    tau: str = "t"                  # the time function, in the Reader's plain names
    families: tuple = RADIAL        # the names of P and M on the page
    cones: tuple = (7, 7)
    areal: bool = False             # spherical: R^2 = g_theta_theta
    areal_contours: tuple = ()
    functions: dict = field(default_factory=dict)   # a declared function -> its expression
    dust: dict = None               # scale factors solved as dust, see DustSolver
    reference: str = None           # label of the line where a dust solution starts, as TeX
    kretschmann: bool = True
    mark_g00: bool = False
    mark_gtt: str = None            # mark g_tt = 0, named in the legend by this prose
    input: str = None               # the declared input, in prose, printed beside the diagram
    principal: bool = False         # the principal null congruence; see "Rays that leave the plane"
    leaves: tuple = ()              # the coordinates its rays move in off the drawn pair
    ring: int = 0                   # a polar view's rays: this many of each family, evenly around
    inside: bool = False            # rays stop short of the edge of the published domain
    cone: str = None                # what a cone is, where it is not the future light cone
    periodic: tuple = ()            # drawn coordinates whose two ends are one line, never hatched
    marked: tuple = ()              # (kind, point, family, legend[, "past"]): a ray of one family,
                                    # or of both, through a point, drawn and named; see Plot.marked
    points: tuple = ()              # (kind, (x^0, r), legend): an event drawn and named
    lines: tuple = ()               # (kind, "x0" or "r", value, legend): a line of constant
                                    # coordinate drawn and named, such as where a declared function jumps
    surface: str = None             # r at which a star's surface is released from rest; see Surface
    singular_runs: bool = False     # mark a singular stretch of an edge, not only a whole edge
    star: dict = None               # a declared polytrope, {"K": ..., "rho_c": ...}; see StarSolver
    any_factor: str = None          # a declared conformal factor the drawing holds for every value of
    crunch: bool = False            # mark where the metric stops being finite as a singular curve,
                                    # checked on the Kretschmann scalar, and hatch what lies beyond it
    solves: tuple = ()              # published Einstein components the declared functions must zero


def _alcubierre_profile():
    return ("(tanh(4*(sqrt((x - 2*t)**2 + y**2 + z**2) + 1))"
            " - tanh(4*(sqrt((x - 2*t)**2 + y**2 + z**2) - 1)))/(2*tanh(4))")


def _natario_field():
    """Natario's zero expansion field for the same profile, n = f/2, as three strings.

    X = v_s [(2n + rho n') e_x - n' x_r (x_r, y, z)/rho], with x_r = x - v_s t. Its
    divergence is checked to vanish before it is used.
    """
    t, x, y, z = sp.symbols("t x y z", real=True)
    rs = sp.Symbol("rs", positive=True)
    xr = x - 2 * t
    rho = sp.sqrt(xr ** 2 + y ** 2 + z ** 2)
    n_of = (sp.tanh(4 * (rs + 1)) - sp.tanh(4 * (rs - 1))) / (4 * sp.tanh(4))
    n, dn = n_of.subs(rs, rho), sp.diff(n_of, rs).subs(rs, rho)
    u = 2 * (2 * n + rho * dn - dn * xr ** 2 / rho)
    v = -2 * dn * xr * y / rho
    w = -2 * dn * xr * z / rho
    div = sp.lambdify((t, x, y, z), sp.diff(u, x) + sp.diff(v, y) + sp.diff(w, z), "numpy")
    points = np.random.default_rng(1).uniform(-2, 2, (4000, 4))
    if np.nanmax(np.abs(div(*points.T))) > 1e-12:
        raise AssertionError("the declared Natario field is not divergence free")
    return {"u": str(u), "v": str(v), "w": str(w)}


_KRASNIKOV_TUBE = ("1 - (2 - 1/5)*(1 + tanh((1 - r**2)/(2*3/20)))/2*(1 + tanh((t - x)/(3/20)))/2"
                   "*(1 + tanh(x/(3/20)))/2*(1 + tanh((4 - x)/(3/20)))/2")

# A conformal factor of the Malament-Hogarth kind: 1 outside the unit ball about the removed
# event and growing as 1/rho toward it, so as 1/|ct| along the axis.
_MH_FACTOR = ("Piecewise((1 + exp(1 - 1/(1 - (t**2 + x**2 + y**2 + z**2)))/sqrt(t**2 + x**2 + y**2 + z**2),"
              " t**2 + x**2 + y**2 + z**2 < 1), (1, True))")

# Marginally bound dust, E = 0, whose density at t = 0 falls as 1 - r^2 to zero at r_b = 1, with
# 2GM/c^2 = r_b/2 and R(r, 0) = r: M(r) = (5r^3 - 3r^5)/8 inside and 1/4 beyond, and each shell
# falls as R^(3/2) = r^(3/2) - (3/2) sqrt(2M) t.
_TB_R = ("Piecewise((r*(1 - 3*sqrt(5 - 3*r**2)*t/4)**Rational(2, 3), r < 1),"
         " ((r**Rational(3, 2) - 3*t/(2*sqrt(2)))**Rational(2, 3), True))")

FRW_DUST = {"funcs": ["a"], "eqs": [["r", "r"]], "rates": [1.0], "params": {"k": 0}}
# The Oppenheimer-Snyder dust released from rest at a = a_m, which is the unit, at tau = 0.
OS_DUST = {"funcs": ["a"], "eqs": [["\\chi", "\\chi"]], "rates": [0.0], "start": [1.0], "origin": "reference"}

# The polytrope the conformal diagram declares, and what it makes of the star.
POLYTROPE = {"K": 100, "rho_c": "1.28e-3"}
POLYTROPE_INPUT = ("A polytrope, $p = K\\rho_0^2$ with rest mass density $\\rho_0$ and energy density "
                   "$\\rho c^2 = \\rho_0c^2 + p$, at $K = 100$ and a central $\\rho_0 = 1.28\\times10^{-3}$ in "
                   "units where $G = c = M_\\odot = 1$, solved from this spacetime's own $G^t{}_t$ and "
                   "$G^r{}_r$: a star of $M = 1.40\\,M_\\odot$ and $R = 14.2$ km.")
POLYTROPE_STATED = {"M": (1.40, 2), "R_km": (14.2, 1)}
KM = 1.4766250614                       # GM_sun/c^2 in km

BIANCHI_DUST = {"funcs": ["a_1", "a_2", "a_3"], "eqs": [["x", "x"], ["y", "y"], ["z", "z"]],
                "rates": [-0.5, 1.5, 2.0]}

# Godel's radius r_c = ln(1 + sqrt 2), where sinh r = 1 and the circles of constant t, r and z
# turn from spacelike to timelike; the views read it off the published g_phiphi as well.
GODEL_RC = math.asinh(1.0)

# Every view the page draws, in the order it shows them. Plot ranges are chosen with
# equal scales on both axes, so light in flat space runs at 45 degrees, and cone
# lattices so that no cone sits exactly on a line where the chart is singular.
DIAGRAMS = [
    Diagram("schwarzschild", "spherical", "radial", "$t$ and $r$", ("t", "r"), (0, 6, -3, 3),
            "$r/r_s$", "$ct/r_s$", {"r_s": 1}, EQUATOR, orient="ingoing", areal=True),
    Diagram("schwarzschild", "eddington_finkelstein_ingoing", "finkelstein", "against $v - r$",
            ("v", "r"), (0, 6, -3, 3), "$r/r_s$", "$(v - r)/r_s$", {"r_s": 1}, EQUATOR,
            to_display=FINKELSTEIN_IN, tau="v - r", areal=True),
    Diagram("schwarzschild", "eddington_finkelstein_ingoing", "chart", "against $v$",
            ("v", "r"), (0, 6, 0, 6), "$r/r_s$", "$v/r_s$", {"r_s": 1}, EQUATOR,
            tau="v - r", areal=True),
    Diagram("schwarzschild", "eddington_finkelstein_outgoing", "finkelstein", "against $u + r$",
            ("u", "r"), (0, 6, -3, 3), "$r/r_s$", "$(u + r)/r_s$", {"r_s": 1}, EQUATOR,
            to_display=FINKELSTEIN_OUT, tau="u + r", areal=True),
    Diagram("schwarzschild", "eddington_finkelstein_outgoing", "chart", "against $u$",
            ("u", "r"), (0, 6, -6, 0), "$r/r_s$", "$u/r_s$", {"r_s": 1}, EQUATOR,
            tau="u + r", areal=True),
    Diagram("frw", "comoving_spherical", "radial", "$t$ and $r$", ("t", "r"), (0, 3, 0, 2),
            "$r\\;[c/H_0]$", "$ct\\;[c/H_0]$", {"k": 0}, EQUATOR, areal=True, dust=FRW_DUST,
            reference="$a = 1$",
            input="Dust: $a(t)$ solved from this spacetime's own $G^r{}_r = 0$ with $k = 0$, "
                  "starting from $a = 1$ and $\\dot a = H_0$ at the dashed line."),
    Diagram("frw", "comoving_spherical", "through", "through the observer", ("t", "r"),
            (0, 1.5, 0, 3), "$x\\;[c/H_0]$", "$ct\\;[c/H_0]$", {"k": 0}, EQUATOR, mirror=True,
            families=SIDEWAYS, cones=(4, 8), areal=True, dust=FRW_DUST, reference="$a = 1$",
            input="Dust: $a(t)$ solved from this spacetime's own $G^r{}_r = 0$ with $k = 0$, "
                  "starting from $a = 1$ and $\\dot a = H_0$ at the dashed line."),
    Diagram("frw", "conformal_spherical", "radial", "$\\eta$ and $r$", ("\\eta", "r"), (0, 3, 0, 3),
            "$r\\;[c/H_0]$", "$\\eta\\;[c/H_0]$", {"k": 0}, EQUATOR, tau="eta", areal=True,
            dust=FRW_DUST, reference="$a = 1$",
            input="Dust, for the Hubble sphere and the Kretschmann scalar: $a(\\eta)$ solved from the "
                  "conformal chart's own $G^r{}_r = 0$ with $k = 0$, starting from $a = 1$ and "
                  "$a' = 1$ at the dashed line."),
    Diagram("ellis_bronnikov", "spherical", "radial", "$t$ and $r$", ("t", "r"), (-3, 3, -3, 3),
            "$r/\\ell$", "$ct/\\ell$", {"ell": 1}, EQUATOR, families=SIDEWAYS, areal=True,
            areal_contours=(1.5, 2.0, 3.0)),
    Diagram("morris_thorne", "spherical", "radial", "$t$ and $r$", ("t", "r"), (0, 4, -2, 2),
            "$r/b_0$", "$ct/b_0$", {"b_0": 1}, EQUATOR, areal=True,
            functions={"Phi": "0", "b": "b_0**2/r"},
            input="$\\Phi = 0$ and $b = b_0^2/r$, the member of the family that is the "
                  "Ellis-Bronnikov wormhole."),
    Diagram("minkowski", "spherical", "radial", "$t$ and $r$", ("t", "r"), (0, 4, -2, 2),
            "$r$", "$ct$", {}, EQUATOR, areal=True),
    Diagram("minkowski", "spherical_null", "radial", "$t$ and $r$", ("u", "v"), (0, 4, -2, 2),
            "$(v - u)/2$", "$(u + v)/2$", {}, EQUATOR, to_display=NULL_TO_TR, tau="u + v",
            areal=True),
    Diagram("minkowski", "cartesian", "tx", "$t$ and $x$", ("t", "x"), (-2, 2, -2, 2),
            "$x$", "$ct$", {}, {"y": "0", "z": "0"}, families=SIDEWAYS),
    Diagram("minkowski", "rindler", "tx", "$T$ and $X$", ("T", "X"), (0, 3, -1.5, 1.5),
            "$X\\;[c^2/a]$", "$cT\\;[c^2/a]$", {"a": 1}, {"Y": "0", "Z": "0"}, tau="T",
            families=SIDEWAYS),
    Diagram("de_sitter", "static_spherical", "radial", "$t$ and $r$", ("t", "r"), (0, 2, -1, 1),
            "$r\\sqrt{\\Lambda/3}$", "$ct\\sqrt{\\Lambda/3}$", {"Lambda": 3}, EQUATOR,
            orient="outgoing", cones=(8, 7), areal=True),
    Diagram("de_sitter", "static_spherical", "through", "through the observer", ("t", "r"),
            (0, 2, -2, 2), "$x\\sqrt{\\Lambda/3}$", "$ct\\sqrt{\\Lambda/3}$", {"Lambda": 3},
            EQUATOR, mirror=True, orient="outgoing", families=SIDEWAYS, cones=(4, 8), areal=True),
    Diagram("de_sitter", "flat_slicing", "tx", "$t$ and $x$", ("t", "x"), (-2, 2, -1, 3),
            "$x\\;[c/H]$", "$ct\\;[c/H]$", {"H": 1}, {"y": "0", "z": "0"}, families=SIDEWAYS),
    Diagram("anti_de_sitter", "static_global", "radial", "$t$ and $r$", ("t", "r"), (0, 4, -2, 2),
            "$r/L$", "$ct/L$", {"L": 1}, EQUATOR, areal=True),
    Diagram("anti_de_sitter", "static_global", "through", "through the centre", ("t", "r"),
            (0, 4, -4, 4), "$x/L$", "$ct/L$", {"L": 1}, EQUATOR, mirror=True,
            families=SIDEWAYS, cones=(4, 8), areal=True),
    Diagram("anti_de_sitter", "poincare", "tx", "$t$ and $x$", ("t", "x"), (-2, 2, -2, 2),
            "$x/L$", "$ct/L$", {"L": 1}, {"y": "0", "z": "1"}, families=SIDEWAYS),
    Diagram("rn_metric", "spherical", "radial", "$t$ and $r$", ("t", "r"), (0, 2, -1, 1),
            "$r/r_s$", "$t/r_s$", {"r_s": 1, "r_q": "2/5"}, EQUATOR, orient="ingoing",
            areal=True),
    Diagram("taub_nut", "spherical", "radial", "$t$ and $r$", ("t", "r"), (0, 6, -3, 3),
            "$r/m$", "$ct/m$", {"m": 1, "l": "1/2"}, EQUATOR, orient="ingoing"),
    Diagram("bertotti_robinson", "static", "radial", "$t$ and $r$", ("t", "r"), (0, 3, -1.5, 1.5),
            "$r/b$", "$ct/b$", {"b": 1}, EQUATOR),
    Diagram("bertotti_robinson", "poincare", "tx", "$t$ and $x$", ("t", "x"), (0, 4, -2, 2),
            "$x/b$", "$ct/b$", {"b": 1}, EQUATOR, families=SIDEWAYS),
    Diagram("interior_schwarzschild", "spherical", "radial", "$t$ and $r$", ("t", "r"),
            (0, 1.5, -0.75, 0.75), "$r/r_s$", "$t/r_s$", {"r_s": 1, "R": "3/2"}, EQUATOR,
            areal=True),
    Diagram("interior_schwarzschild", "spherical", "through", "through the centre", ("t", "r"),
            (0, 1.5, -1.5, 1.5), "$x/r_s$", "$t/r_s$", {"r_s": 1, "R": "3/2"}, EQUATOR,
            mirror=True, families=SIDEWAYS, cones=(4, 8), areal=True),
    Diagram("kerr", "boyer_lindquist", "radial", "$t$ and $r$ on the axis", ("t", "r"), (0, 4, -2, 2),
            "$r/(GM/c^2)$", "$ct/(GM/c^2)$", {"G": 1, "M": 1, "a": "9/10"},
            {"theta": "0", "phi": "0"}, orient="ingoing"),
    # The principal null congruence on the equator, where the ergoregion is widest and the
    # ring singularity lies, drawn twice: its shadow on t and r, which carries the cones and
    # the horizons as the axis does, and its shadow on the equatorial plane seen from above,
    # which carries the turning in phi. r changes monotonically along every ray, so the two
    # shadows, sharing r, fix the ray. Twelve rays of each family from above, one every 30
    # degrees, keep the two spirals apart at the ergosurface, and a box of 3 GM/c^2 each
    # way keeps the winding onto the horizon as large as the straight runs far out.
    Diagram("kerr", "boyer_lindquist", "principal", "principal null rays, $t$ and $r$", ("t", "r"),
            (0, 4, -2, 2), "$r/(GM/c^2)$", "$ct/(GM/c^2)$", {"G": 1, "M": 1, "a": "9/10"},
            {"theta": "pi/2"}, orient="ingoing", principal=True, leaves=("phi",),
            mark_gtt="the ergosurface", cone=PRINCIPAL_CONE),
    Diagram("kerr", "boyer_lindquist", "above", "principal null rays from above", ("\\phi", "r"),
            (-3, 3, -3, 3), "$r\\cos\\phi/(GM/c^2)$", "$r\\sin\\phi/(GM/c^2)$",
            {"G": 1, "M": 1, "a": "9/10"}, {"theta": "pi/2"}, to_display=POLAR, cones=(0, 0),
            principal=True, leaves=("t",), ring=12, inside=True, mark_gtt="the ergosurface"),
    Diagram("kerr_newman", "boyer_lindquist", "radial", "$t$ and $r$ on the axis", ("t", "r"),
            (0, 4, -2, 2), "$r/(GM/c^2)$", "$ct/(GM/c^2)$",
            {"G": 1, "M": 1, "a": "3/5", "r_Q": "1/2"}, {"theta": "0", "phi": "0"},
            orient="ingoing"),
    Diagram("kerr_newman", "boyer_lindquist", "principal", "principal null rays, $t$ and $r$",
            ("t", "r"), (0, 4, -2, 2), "$r/(GM/c^2)$", "$ct/(GM/c^2)$",
            {"G": 1, "M": 1, "a": "3/5", "r_Q": "1/2"}, {"theta": "pi/2"}, orient="ingoing",
            principal=True, leaves=("phi",), mark_gtt="the ergosurfaces", cone=PRINCIPAL_CONE),
    Diagram("kerr_newman", "boyer_lindquist", "above", "principal null rays from above",
            ("\\phi", "r"), (-3, 3, -3, 3), "$r\\cos\\phi/(GM/c^2)$", "$r\\sin\\phi/(GM/c^2)$",
            {"G": 1, "M": 1, "a": "3/5", "r_Q": "1/2"}, {"theta": "pi/2"}, to_display=POLAR,
            cones=(0, 0), principal=True, leaves=("t",), ring=12, inside=True,
            mark_gtt="the ergosurfaces"),
    Diagram("kasner", "cartesian", "tx", "$t$ and $x$", ("t", "x"), (-1, 1, 0, 2), "$x$", "$ct$",
            {"p_1": "-2/7", "p_2": "3/7", "p_3": "6/7"}, {"y": "0", "z": "0"},
            families=SIDEWAYS,
            input="Exponents $(p_1, p_2, p_3) = (-2/7, 3/7, 6/7)$, a point on the Kasner circle: "
                  "they sum to 1, and so do their squares."),
    Diagram("kasner", "cartesian", "tz", "$t$ and $z$", ("t", "z"), (-1, 1, 0, 2), "$z$", "$ct$",
            {"p_1": "-2/7", "p_2": "3/7", "p_3": "6/7"}, {"x": "0", "y": "0"},
            families=SIDEWAYS,
            input="Exponents $(p_1, p_2, p_3) = (-2/7, 3/7, 6/7)$, a point on the Kasner circle: "
                  "they sum to 1, and so do their squares."),
    Diagram("bianchi", "type_i_cartesian", "tx", "$t$ and $x$", ("t", "x"), (-1, 1, 0, 2),
            "$x\\;[c/\\bar H]$", "$ct\\;[c/\\bar H]$", {}, {"y": "0", "z": "0"}, families=SIDEWAYS,
            dust=BIANCHI_DUST, reference="$a_i = 1$",
            input="Dust: the three scale factors solved from this spacetime's own "
                  "$G^x{}_x = G^y{}_y = G^z{}_z = 0$, starting from $a_i = 1$ with rates "
                  "$(-0.5, 1.5, 2.0)\\,\\bar H$ at the dashed line, $\\bar H$ their mean."),
    Diagram("godel", "cartesian", "tx", "$t$ and $x$", ("t", "x"), (-2, 2, -2, 2), "$x$", "$t$",
            {"omega": 1}, {"y": "0", "z": "0"}, orient="vector", families=SIDEWAYS),
    Diagram("alcubierre", "cartesian", "tx", "$t$ and $x$ on the axis", ("t", "x"), (-3, 3, -2, 2),
            "$x/R$", "$ct/R$", {}, {"y": "0", "z": "0"}, families=SIDEWAYS, cones=(8, 7),
            functions={"v_s": "2", "f": _alcubierre_profile()},
            input="$v_s = 2$, and Alcubierre's own profile, "
                  "$f = [\\tanh\\sigma(r_s + R) - \\tanh\\sigma(r_s - R)]/(2\\tanh\\sigma R)$ "
                  "with $R = 1$ and $\\sigma = 4$."),
    Diagram("natario", "cartesian_flow", "tx", "$t$ and $x$ on the axis", ("t", "x"), (-3, 3, -2, 2),
            "$x/R$", "$ct/R$", {}, {"y": "0", "z": "0"}, families=SIDEWAYS, cones=(8, 7),
            functions=_natario_field(),
            input="$v_s = 2$, $n = f/2$ with Alcubierre's profile, and the zero expansion field "
                  "$X = v_s[(2n + \\rho n')\\,e_x - n'\\,x_r\\,(x_r, y, z)/\\rho]$, "
                  "$x_r = x - v_s t$, whose divergence vanishes."),
    Diagram("krasnikov", "cylindrical", "tx", "$t$ and $x$ on the axis", ("t", "x"), (-1, 5, -1, 5),
            "$x$", "$ct$", {}, {"r": "0", "phi": "0"}, orient="outgoing", families=SIDEWAYS,
            functions={"k": _KRASNIKOV_TUBE}, kretschmann=False, mark_g00=True,
            input="A tube along $x$ from $0$ to $D = 4$, built by a ship that left $x = 0$ at "
                  "$t = 0$ at the speed of light: $k = 1 - (2 - \\delta)\\,S(\\tfrac{\\rho_0^2 - r^2}"
                  "{2\\rho_0})\\,S(t - x)\\,S(x)\\,S(D - x)$ with $\\delta = 0.2$, $\\rho_0 = 1$ and "
                  "$S$ a step of width $0.15$ built from $\\tanh$."),
    Diagram("pp_wave", "exact_plane_wave", "tz", "$t$ and $z$ on the axis", ("u", "v"), (-2, 2, -2, 2),
            "$z$", "$ct$", {}, {"x": "0", "y": "0"}, to_display=UV_TO_TZ, tau="u + 2*v",
            families=SIDEWAYS, functions={"A": "exp(-u**2)", "B": "0"},
            input="The amplitudes $A(u)$ and $B(u)$ do not enter on the axis, so any profile gives "
                  "this diagram."),
    # The cylinders of t and phi at one r, unrolled, phi scaled by r/R so that the cones next to
    # the axis would stand at 45 degrees. The metric on each is the same at every point of it.
    Diagram("stockum_dust", "cylindrical", "inside", "$t$ and $\\phi$ at $r = R/2$", ("t", "\\phi"),
            (-math.pi / 2, math.pi / 2, -math.pi / 2, math.pi / 2), "$r\\phi/R$", "$t/R$", {"R": 1},
            {"r": "1/2", "z": "0"}, to_display=((0, 0.5), (1, 0)), orient="vector", families=SIDEWAYS,
            cones=(5, 5), periodic=("\\phi",)),
    Diagram("stockum_dust", "cylindrical", "beyond", "$t$ and $\\phi$ at $r = 3R/2$", ("t", "\\phi"),
            (-3 * math.pi / 2, 3 * math.pi / 2, -3 * math.pi / 2, 3 * math.pi / 2), "$r\\phi/R$", "$t/R$",
            {"R": 1}, {"r": "3/2", "z": "0"}, to_display=((0, 1.5), (1, 0)), orient="vector",
            families=SIDEWAYS, cones=(5, 5), periodic=("\\phi",)),
    # Godel's cylinders of t and phi about one world line of the dust, drawn as van Stockum's are,
    # phi scaled by r, at half and one and a half times r_c.
    Diagram("godel", "cylindrical", "inside", "$t$ and $\\phi$ at $r = r_c/2$", ("t", "\\phi"),
            tuple(s * math.pi * GODEL_RC / 2 for s in (-1, 1, -1, 1)), "$r\\phi$", "$t$", {"omega": 1},
            {"r": "asinh(1)/2", "z": "0"}, to_display=((0, GODEL_RC / 2), (1, 0)), orient="vector",
            families=SIDEWAYS, cones=(5, 5), periodic=("\\phi",)),
    Diagram("godel", "cylindrical", "beyond", "$t$ and $\\phi$ at $r = 3r_c/2$", ("t", "\\phi"),
            tuple(s * math.pi * 3 * GODEL_RC / 2 for s in (-1, 1, -1, 1)), "$r\\phi$", "$t$", {"omega": 1},
            {"r": "3*asinh(1)/2", "z": "0"}, to_display=((0, 3 * GODEL_RC / 2), (1, 0)), orient="vector",
            families=SIDEWAYS, cones=(5, 5), periodic=("\\phi",)),
    Diagram("tov", "spherical", "radial", "$t$ and $r$", ("t", "r"), (0, 16, -8, 8),
            "$r\\;[GM_\\odot/c^2]$", "$ct\\;[GM_\\odot/c^2]$", {}, EQUATOR, areal=True, star=POLYTROPE,
            input=POLYTROPE_INPUT),
    Diagram("tov", "spherical", "through", "through the centre", ("t", "r"), (0, 16, -16, 16),
            "$x\\;[GM_\\odot/c^2]$", "$ct\\;[GM_\\odot/c^2]$", {}, EQUATOR, mirror=True, families=SIDEWAYS,
            cones=(4, 8), areal=True, star=POLYTROPE, input=POLYTROPE_INPUT),
    Diagram("malament_hogarth", "cartesian", "tx", "$t$ and $x$", ("t", "x"), (-2, 2, -2, 2), "$x$", "$ct$", {},
            {"y": "0", "z": "0"}, families=SIDEWAYS, functions={"Omega": _MH_FACTOR}, any_factor="Omega",
            lines=(("world", "r", "0", "the computer's world line, up the axis into the removed event",
                    ("-1000000", "0")),),
            points=(("removed", ("0", "0"), "the removed event"), ("mark", ("1", "0"), "the event $p$")),
            marked=(("past", {"x0": "1", "r": "0"}, "both", "the past light cone of $p$", "past"),),
            input="Any $\\Omega$, since the rays and cones are the same for every one; drawn with "
                  "$\\Omega = 1 + e^{1 - 1/(1 - \\rho^2)}/\\rho$ for $\\rho^2 = c^2t^2 + x^2 + y^2 + z^2 < 1$ "
                  "and $\\Omega = 1$ beyond, which grows as $1/|ct|$ along the axis, and checked to give "
                  "the null directions of $\\Omega = 1$."),
    Diagram("tolman_bondi", "comoving_synchronous", "collapse", "a collapsing cloud", ("t", "r"),
            (0, 1.5, -0.5, 1.0), "$r/r_b$", "$ct/r_b$", {}, EQUATOR, areal=True,
            functions={"E": "0", "R": _TB_R}, solves=(("r", "r"),), crunch=True,
            lines=(("surface", "r", "1", "the surface of the cloud, $r = r_b$"),),
            marked=(("event", {"r": "11/10", "areal": "1/2"}, 1, "the event horizon"),),
            input="Marginally bound dust, $E = 0$, its density at $t = 0$ falling as $1 - r^2/r_b^2$ to "
                  "zero at $r_b$, with $R(r, 0) = r$ and $2GM/c^2 = r_b/2$: each shell falls as "
                  "$R^{3/2} = r^{3/2} - \\tfrac{3}{2}\\sqrt{2GM(r)/c^2}\\,ct$, checked to solve this "
                  "spacetime's own $G^r{}_r = 0$."),
    # The collapse the conformal diagram draws, released from rest at R_0 = 2 r_s, chi_0 = pi/4:
    # the dust in its own chart, and the vacuum outside it in Schwarzschild's.
    Diagram("oppenheimer_snyder", "interior_comoving", "through", "through the centre", ("\\tau", "\\chi"),
            (0, math.pi / 4, 0, math.pi / 2), "$\\chi$", "$c\\tau/a_m$", {"chi_0": "pi/4", "a_m": 1}, EQUATOR,
            mirror=True, families=SIDEWAYS, cones=(4, 8), tau="tau", areal=True, dust=OS_DUST,
            lines=(("surface", "r", "pi/4", "the surface of the star, $\\chi = \\chi_0$"),),
            marked=(("event", {"r": "pi/4", "areal": "sqrt(2)/4"}, 1, "the event horizon"),),
            input="Dust released from rest at $\\tau = 0$ with $a = a_m$, $a(\\tau)$ solved from this "
                  "spacetime's own $G^\\chi{}_\\chi = 0$, and $\\chi_0 = \\pi/4$, so that the star starts "
                  "at twice its Schwarzschild radius."),
    Diagram("oppenheimer_snyder", "exterior_schwarzschild", "radial", "$t$ and $r$", ("t", "r"), (0, 3, 0, 6),
            "$r/r_s$", "$ct/r_s$", {"r_s": 1}, EQUATOR, orient="ingoing", areal=True, surface="2",
            input="The surface released from rest at $r = 2r_s$ at $t = 0$, the star of the interior "
                  "coordinates with $\\chi_0 = \\pi/4$."),
    # A shell of null dust falls in along v = 0: flat inside, Schwarzschild outside, as the
    # conformal diagram declares it. The outgoing chart draws its time reverse.
    Diagram("vaidya", "eddington_finkelstein_ingoing", "shell", "an imploding shell", ("v", "r"),
            (0, 4, -2.5, 1.5), "$r/r_s$", "$(cv - r)/r_s$", {"G": 1}, EQUATOR, to_display=FINKELSTEIN_IN,
            tau="v - r", areal=True, functions={"m": "Heaviside(v)/2"}, lines=(("shell", "x0", "0", "the shell, $v = 0$"),),
            marked=(("event", "1/1000", 1, "the event horizon"),), singular_runs=True,
            input="$m(v) = 0$ for $v < 0$ and $M$ for $v > 0$, with $r_s = 2GM/c^2$: a shell of null "
                  "dust of mass $M$ falling in along $v = 0$."),
    Diagram("vaidya", "eddington_finkelstein_outgoing", "shell", "an exploding shell", ("u", "r"),
            (0, 4, -1.5, 2.5), "$r/r_s$", "$(cu + r)/r_s$", {"G": 1}, EQUATOR, to_display=FINKELSTEIN_OUT,
            tau="u + r", areal=True, functions={"m": "Heaviside(-u)/2"}, lines=(("shell", "x0", "0", "the shell, $u = 0$"),),
            marked=(("event", "-1/1000", 0, "the white hole's horizon"),), singular_runs=True,
            input="$m(u) = M$ for $u < 0$ and $0$ for $u > 0$, with $r_s = 2GM/c^2$: a shell of null "
                  "dust carrying off the whole mass $M$ along $u = 0$."),
]


# ---------------------------------------------------------------- the captions

CAPTIONS = {
    ("schwarzschild", "spherical", "radial"): [
        "This is the plane of $t$ and $r$ at $\\theta = \\pi/2$ and $\\phi = 0$, and by spherical "
        "symmetry every fixed angle gives the same picture. Outside $r_s$ the cones narrow toward "
        "the vertical as $r \\to r_s$, because $dt/dr = \\pm(1 - r_s/r)^{-1}$ diverges there. The "
        "ingoing family piles up against $r_s$ toward $t \\to +\\infty$, and the outgoing family "
        "peels away from it out of $t \\to -\\infty$; the two Eddington-Finkelstein charts carry "
        "them across.",
        "The domain of this chart stops at $r_s$. Evaluated inside it, the same components give "
        "cones lying on their side, because there it is $r$ that is the time. They cannot say "
        "whether that region is the black hole or the white hole; the ingoing "
        "Eddington-Finkelstein chart, which runs smoothly across $r_s$, makes it the black hole, "
        "and there every cone points to $r = 0$. The Kretschmann scalar $12r_s^2/r^6$ is finite "
        "at $r_s$ and diverges only at $r = 0$.",
    ],
    ("schwarzschild", "eddington_finkelstein_ingoing", "finkelstein"): [
        "This is the plane of $v$ and $r$ at $\\theta = \\pi/2$ and $\\phi = 0$, drawn with $v - "
        "r$ as the vertical axis so that the ingoing rays, $v = $ const, run at 45°. The chart "
        "crosses the horizon. The outgoing family has $dv/dr = 2(1 - r_s/r)^{-1}$, so it stands "
        "exactly vertical at $r_s$: the horizon is itself an outgoing ray that stays where it is.",
        "The cones cross $r_s$ smoothly and keep tipping. Inside, both edges of every future cone "
        "point to smaller $r$, so every future directed ray ends at $r = 0$.",
    ],
    ("schwarzschild", "eddington_finkelstein_ingoing", "chart"): [
        "This is the same plane of $v$ and $r$ at $\\theta = \\pi/2$ and $\\phi = 0$, drawn "
        "against the chart's own coordinates. The ingoing family is $v = $ const and runs "
        "horizontally here, since $v$ is itself a null coordinate. The outgoing family turns "
        "vertical at $r_s$ and leans back toward smaller $r$ inside it.",
    ],
    ("schwarzschild", "eddington_finkelstein_outgoing", "finkelstein"): [
        "This is the plane of $u$ and $r$ at $\\theta = \\pi/2$ and $\\phi = 0$, drawn with $u + "
        "r$ as the vertical axis so that the outgoing rays, $u = $ const, run at 45°. The retarded "
        "chart crosses the other horizon. Inside $r_s$ both edges of every future cone point to "
        "larger $r$: this is the white hole, which nothing from outside can enter.",
        "The chart is the time reverse of the ingoing one: its $g_{ur}$ is $-1$ where the ingoing "
        "chart's $g_{vr}$ is $+1$.",
    ],
    ("schwarzschild", "eddington_finkelstein_outgoing", "chart"): [
        "This is the same plane of $u$ and $r$ at $\\theta = \\pi/2$ and $\\phi = 0$, drawn "
        "against the chart's own coordinates. The outgoing family is $u = $ const and runs "
        "horizontally here, since $u$ is itself a null coordinate. The ingoing family turns "
        "vertical at $r_s$ and leans toward larger $r$ inside it.",
    ],
    ("frw", "comoving_spherical", "radial"): [
        "This is the plane of $t$ and the comoving $r$ at $\\theta = \\pi/2$ and $\\phi = 0$, and "
        "spherical symmetry makes it the same at every other angle. The line element leaves "
        "$a(t)$ free, and here it is the scale factor of pressureless dust, solved from "
        "$G^r{}_r = -2\\ddot a/a - \\dot a^2/a^2 - k/a^2 = 0$. Run back from the dashed line, it "
        "reaches $a = 0$ a time $2/(3H_0)$ before it, the Einstein-de Sitter age, and $t$ is "
        "counted from that instant, the big bang.",
        "The edges of the cones are $dt/dr = \\pm a(t)$, and they open out toward the bang, where "
        "the curvature diverges, so every ray leaves $t = 0$ almost flat and reaches only a finite "
        "comoving distance, the particle horizon. The dotted curve, where $|\\nabla R|^2 = 0$ for "
        "the areal radius $R = ar$, is the Hubble sphere $R = c/H$, which is the apparent horizon "
        "of an observer at $r = 0$.",
    ],
    ("frw", "comoving_spherical", "through"): [
        "This is the same universe along a line through the observer at $r = 0$, in the plane "
        "$\\theta = \\pi/2$: $x = r$ on the right is $\\phi = 0$ and $x = -r$ on the left is "
        "$\\phi = \\pi$, and spherical symmetry makes the two halves mirror images. The past light "
        "cone of an event on the observer's world line $x = 0$ flares out as it runs back toward "
        "the bang, and the Hubble sphere lies at the same distance on either side.",
    ],
    ("frw", "conformal_spherical", "radial"): [
        "This is the plane of the conformal time $\\eta$ and $r$ at $\\theta = \\pi/2$ and $\\phi "
        "= 0$. For $k = 0$ the metric on it is $a^2(-d\\eta^2 + dr^2)$; the scale factor "
        "multiplies both terms and drops out of the null condition, so the rays are straight 45° "
        "lines whatever $a(\\eta)$ is.",
        "The Hubble sphere and the Kretschmann scalar do depend on $a(\\eta)$, and for dust the "
        "Hubble sphere sits at $r = \\eta/2$, half the comoving radius of the particle horizon.",
    ],
    ("ellis_bronnikov", "spherical", "radial"): [
        "This is the plane of $t$ and $r$ at $\\theta = \\pi/2$ and $\\phi = 0$, with $r$ running "
        "from one mouth of the wormhole, $r < 0$, through the throat at $r = 0$ to the other "
        "mouth, $r > 0$. The metric on this plane is $-dt^2 + dr^2$, exactly flat, so the rays are "
        "straight 45° lines that pass through the throat without bending.",
        "The throat is where the spheres are smallest. The angular part of the metric is "
        "$g_{\\theta\\theta} = r^2 + \\ell^2$, so the areal radius $R = \\sqrt{r^2 + \\ell^2}$ "
        "takes its least value, $\\ell$, at $r = 0$, where $\\partial_r R = 0$. The faint vertical "
        "lines are the spheres $R = 1.5\\ell$, $2\\ell$ and $3\\ell$, one of each size on either "
        "side of the throat. The Kretschmann scalar $12\\ell^4/(r^2 + \\ell^2)^4$ is finite "
        "everywhere.",
    ],
    ("morris_thorne", "spherical", "radial"): [
        "This is the plane of $t$ and the areal radius $r$ at $\\theta = \\pi/2$ and $\\phi = 0$. "
        "The metric leaves $\\Phi(r)$ and $b(r)$ free. With $\\Phi = 0$ and $b = b_0^2/r$ it is the "
        "Ellis-Bronnikov wormhole, with $r^2 = r_{\\rm EB}^2 + \\ell^2$ and $b_0 = \\ell$, and "
        "these rays conserve $t \\mp \\sqrt{r^2 - b_0^2} = t \\mp r_{\\rm EB}$: they are the same "
        "rays as in the Ellis-Bronnikov chart.",
        "In this areal chart the cones close toward the throat at $r = b_0$, as they would at a "
        "horizon, because $g_{rr} = (1 - b_0^2/r^2)^{-1}$ diverges there. But $g_{tt} = -1$ stays "
        "finite, so $\\partial_t$ is timelike right up to the throat, and $r = b_0$ is only the "
        "edge of this chart. The rays reach the throat in finite $t$ and pass into the other "
        "mouth, which the Ellis-Bronnikov chart, running through the throat, covers in full. "
        "Below $b_0$ the formula gives a metric on this plane with no null directions at all.",
    ],
    ("minkowski", "spherical", "radial"): [
        "This is the plane of $t$ and $r$ at $\\theta = \\pi/2$ and $\\phi = 0$ in flat spacetime. "
        "The metric on it is $-c^2dt^2 + dr^2$, and every ray is at 45°. Each ingoing ray meets "
        "an outgoing one on the axis $r = 0$.",
    ],
    ("minkowski", "spherical_null", "radial"): [
        "This is the same plane of flat spacetime in the double null chart, $u$ and $v$ at "
        "$\\theta = \\pi/2$ and $\\phi = 0$, drawn against $(v - u)/2$ and $(u + v)/2$. Only "
        "$g_{uv}$ is nonzero on it, so a null direction has $du\\,dv = 0$, and the light rays are "
        "the coordinate lines $u = $ const and $v = $ const themselves.",
    ],
    ("minkowski", "cartesian", "tx"): [
        "This is the plane of $t$ and $x$ at $y = z = 0$ in flat spacetime. The metric on it is "
        "$-c^2dt^2 + dx^2$, and every ray is at 45°.",
    ],
    ("minkowski", "rindler", "tx"): [
        "This is the plane of $T$ and $X$ at $Y = Z = 0$ in the Rindler chart, which covers the "
        "wedge $X > 0$ seen by observers of constant proper acceleration $a$, with $g_{TT} = "
        "-a^2X^2/c^4$. The cones close toward $X = 0$, where $g_{TT}$ vanishes, and a ray takes "
        "infinite $T$ to get there, $cT = \\pm(c^2/a)\\ln X + $ const.",
        "That line is the Rindler horizon, which only the accelerated observers have. The "
        "spacetime is flat, its Kretschmann scalar is zero, and the rays carry on across $X = 0$ "
        "into the rest of Minkowski spacetime, which the Cartesian chart covers whole.",
    ],
    ("de_sitter", "static_spherical", "radial"): [
        "This is the plane of $t$ and $r$ at $\\theta = \\pi/2$ and $\\phi = 0$ in the static "
        "chart, and spherical symmetry makes it the same at every other angle. The cones close at "
        "the cosmological horizon $r = \\sqrt{3/\\Lambda}$, where $g^{rr} = 1 - \\Lambda r^2/3$ "
        "vanishes, on the far side from the observer at $r = 0$.",
        "The static chart covers only the observer's side of the horizon. Beyond it $t$ is "
        "spacelike, and the cones point along the outgoing rays to larger $r$, into the region "
        "the observer's own light goes on to reach.",
    ],
    ("de_sitter", "static_spherical", "through"): [
        "This is the static chart along a line through the observer in the plane $\\theta = "
        "\\pi/2$: $x = r$ on the right is $\\phi = 0$ and $x = -r$ on the left is $\\phi = \\pi$, "
        "and spherical symmetry makes the two halves mirror images. The horizon crosses the line "
        "on both sides, at $x = \\pm\\sqrt{3/\\Lambda}$, and beyond it the cones point away from "
        "the observer.",
    ],
    ("de_sitter", "flat_slicing", "tx"): [
        "This is the plane of $t$ and $x$ at $y = z = 0$ in the flat slicing, where $ds^2 = "
        "-c^2dt^2 + e^{2Ht}dx^2$, so the cones narrow as $e^{-Ht}$ toward the future and open out "
        "toward the past.",
        "A ray covers only a finite comoving distance however long it runs, $x = \\pm(c/H)e^{-Ht} "
        "+ $ const, so an observer at $x = 0$ has an event horizon.",
    ],
    ("anti_de_sitter", "static_global", "radial"): [
        "This is the plane of $t$ and $r$ at $\\theta = \\pi/2$ and $\\phi = 0$ in the global "
        "chart. The cones stay open everywhere, since $g^{rr} = 1 + r^2/L^2$ never vanishes. But "
        "$dt/dr = \\pm(1 + r^2/L^2)^{-1}$ falls off fast enough that a ray reaches $r \\to "
        "\\infty$, the conformal boundary, in the finite time $\\pi L/2$, and the rays flatten as "
        "they near it.",
    ],
    ("anti_de_sitter", "static_global", "through"): [
        "This is the global chart along a line through the centre in the plane $\\theta = \\pi/2$: "
        "$x = r$ on the right is $\\phi = 0$ and $x = -r$ on the left is $\\phi = \\pi$. A ray "
        "from the centre reaches the boundary on either side in the finite time $\\pi L/2$, and "
        "the rays flatten as they near it.",
    ],
    ("anti_de_sitter", "poincare", "tx"): [
        "This is the plane of $t$ and $x$ at $y = 0$ and $z = L$ in the Poincaré chart. The metric "
        "on it is $(L^2/z^2)(-c^2dt^2 + dx^2)$, conformal to flat, so the rays are exact 45° "
        "lines, as they are at every $z$. The conformal boundary, which a ray reaches in finite "
        "time, lies off this plane, at $z \\to 0$.",
    ],
    ("rn_metric", "spherical", "radial"): [
        "This is the plane of $t$ and $r$ at $\\theta = \\pi/2$ and $\\phi = 0$, drawn for "
        "$r_q = 0.4\\,r_s$, and spherical symmetry makes it the same at every other angle. There "
        "$g^{rr}$ vanishes twice, at $r_\\pm = (r_s \\pm \\sqrt{r_s^2 - 4r_q^2})/2$, which is "
        "$0.8\\,r_s$ and $0.2\\,r_s$, and the cones close at both. Between them $r$ is the time "
        "and the cones point to smaller $r$. Inside $r_-$, $t$ is a time again.",
        "The chart alone cannot say which way is future in the two inner regions. An ingoing "
        "chart runs smoothly through both horizons, and taking the future from it makes the "
        "region between the horizons the black hole. The Kretschmann scalar diverges at "
        "$r = 0$.",
    ],
    ("taub_nut", "spherical", "radial"): [
        "This is the plane of $t$ and $r$ at $\\theta = \\pi/2$ and $\\phi = 0$, drawn for $l = "
        "m/2$. There $g^{rr}$ vanishes at $r_+ = m + \\sqrt{m^2 + l^2}$, about $2.118\\,m$, and "
        "inside it is the Taub region, where $r$ is the time and the cones, following the ingoing "
        "family, point to smaller $r$.",
        "The twist the NUT parameter brings sits in the cross term $g_{t\\phi}$, which drops out "
        "of the metric on this plane. No Christoffel symbol turns these null curves out of the "
        "plane either, so they are null geodesics, the paths light takes, even though the "
        "solution is not spherically symmetric.",
    ],
    ("bertotti_robinson", "static", "radial"): [
        "This is the plane of $t$ and $r$ at $\\theta = \\pi/2$ and $\\phi = 0$, which is the AdS₂ "
        "factor of the product, $-(r^2/b^2)\\,c^2dt^2 + (b^2/r^2)\\,dr^2$. With $dt/dr = \\pm "
        "b^2/r^2$ the rays take infinite $t$ to reach $r = 0$, where $g^{rr} = r^2/b^2$ vanishes. "
        "That is a Poincaré horizon, the kind of edge that bounds the Poincaré chart of anti-de "
        "Sitter space. The other factor is the two sphere, of radius $b$ everywhere, and each "
        "point of this plane stands for one such sphere.",
    ],
    ("bertotti_robinson", "poincare", "tx"): [
        "This is the plane of $t$ and $x$ at $\\theta = \\pi/2$ and $\\phi = 0$ in the Poincaré "
        "chart, where the spacetime is the product of $(b^2/x^2)(-c^2dt^2 + dx^2)$ with a two "
        "sphere of radius $b$. The first factor is conformal to flat, so the rays are at 45°. "
        "Here $x$ is a coordinate on the AdS₂ factor, with the boundary at $x \\to 0$ and the "
        "Poincaré horizon at $x \\to \\infty$.",
    ],
    ("interior_schwarzschild", "spherical", "radial"): [
        "This is the plane of $t$ and $r$ at $\\theta = \\pi/2$ and $\\phi = 0$ over the whole "
        "domain of the chart, $r \\in [0, R]$, for a star with $R = 1.5\\,r_s$. The cones are "
        "narrowest at the centre, where $|g_{tt}|$ is least and the redshift greatest, and they "
        "stay open. They would close at the centre exactly when $3\\sqrt{1 - r_s/R} = 1$, which is "
        "Buchdahl's $R = 9r_s/8$. Beyond $R$ the spacetime is Schwarzschild's exterior, and the rays "
        "go on into it as they do in Schwarzschild's own chart.",
    ],
    ("interior_schwarzschild", "spherical", "through"): [
        "This is the line through the centre of the star in the plane $\\theta = \\pi/2$: $x = r$ "
        "on the right is $\\phi = 0$ and $x = -r$ on the left is $\\phi = \\pi$. Rays cross the "
        "centre smoothly, and the cones are narrowest there.",
    ],
    ("kerr", "boyer_lindquist", "radial"): [
        "This is the plane of $t$ and $r$ on the rotation axis, $\\theta = 0$, drawn for $a = "
        "0.9\\,GM/c^2$. The curves drawn are null, and on the axis they are also null geodesics, "
        "the paths light takes. Off the axis a light ray launched along a curve of fixed $\\theta$ "
        "and $\\phi$ is turned out of the plane by $\\Gamma^\\theta{}_{tt}$, "
        "$\\Gamma^\\theta{}_{rr}$ and $\\Gamma^\\phi{}_{tr}$. On the axis $g^{rr} = \\Delta/\\Sigma$ "
        "vanishes at both roots of $\\Delta$, $r_\\pm = GM/c^2 \\pm \\sqrt{(GM/c^2)^2 - a^2}$, "
        "and the cones close at both; between them they point to smaller $r$, following the "
        "ingoing family.",
        "The domain of the chart begins at $r_+$, the outer zero of $g^{rr}$. On the axis the "
        "ergosurface touches the horizon, so the metric on the plane stays Lorentzian. The "
        "Kretschmann scalar stays finite at $r = 0$ on the axis, because the ring singularity "
        "lies in the equatorial plane.",
    ],
    ("kerr", "boyer_lindquist", "principal"): [
        "This is the equatorial plane $\\theta = \\pi/2$ drawn in $t$ and $r$, with $\\phi$ left "
        "out, for $a = 0.9\\,GM/c^2$. Its rays are Kerr's principal null congruence, the light "
        "rays that run straight in and straight out. At each point they are the two null "
        "directions of the plane of $\\partial_r$ and $(r^2 + a^2)\\,\\partial_t + "
        "a\\,\\partial_\\phi$, which are the repeated principal null directions of the Weyl "
        "tensor. That plane tilts into $\\phi$, so every ray turns as it goes, at $d\\phi/dr = "
        "\\pm a/\\Delta$ with $\\Delta = r^2 - 2GMr/c^2 + a^2$, and the curves drawn are the rays' "
        "projections, $d(ct)/dr = \\pm(r^2 + a^2)/\\Delta$. The projections are the same at every "
        "$\\theta$, and on the axis they are the radial light rays themselves.",
        "The rays are null geodesics, the paths light takes, and each cone is the future cone of "
        "the principal plane, its edges the two principal directions. The dotted line is the "
        "ergosurface, $g_{tt} = 0$, at $r = 2GM/c^2$ on the equator. Between it and $r_+$ nothing, "
        "light included, can keep $\\phi$ fixed, and the plane of $t$ and $r$ at fixed $\\phi$ has "
        "no null direction, while the principal rays, already turning, cross it smoothly. The "
        "cones close at both horizons, where $\\Delta = 0$, and point to smaller $r$ between them. "
        "The ingoing rays end on the ring singularity at $r = 0$, which lies in this plane.",
    ],
    ("kerr", "boyer_lindquist", "above"): [
        "This is the equatorial plane $\\theta = \\pi/2$ seen from above, along the axis from "
        "$\\theta = 0$, with $r$ and $\\phi$ drawn as polar coordinates and $t$ left out, for "
        "$a = 0.9\\,GM/c^2$. Its rays are Kerr's principal null congruence, each turning at "
        "$d\\phi/dr = \\pm a/\\Delta$ with $\\Delta = r^2 - 2GMr/c^2 + a^2$, so both families wind "
        "counterclockwise, the way the hole turns: the ingoing rays as they fall and the outgoing "
        "rays as they climb. Far out the turning dies away as $a/r^2$ and the rays run nearly "
        "straight. The rate $d\\phi/dr$ is the same at every $\\theta$, and on the equator the rays "
        "stay in this plane.",
        "At the horizon $r_+$ the angle $\\phi$ runs to infinity along every ray, as $t$ does in "
        "the plane of $t$ and $r$, so the ingoing rays wind without end onto the horizon and the "
        "outgoing rays unwind off it. The winding is in the coordinate $\\phi$ alone: along an "
        "ingoing ray $\\tilde\\phi = \\phi + \\int a\\,dr/\\Delta$ stays fixed, and in it the ray "
        "crosses the horizon at a finite angle. The dotted circle is the ergosurface, "
        "$r = 2GM/c^2$ on the equator, and between it and $r_+$ nothing, light included, can keep "
        "$\\phi$ fixed.",
    ],
    ("kerr_newman", "boyer_lindquist", "radial"): [
        "This is the plane of $t$ and $r$ on the rotation axis, $\\theta = 0$, drawn for $a = "
        "0.6\\,GM/c^2$ and $r_Q = 0.5\\,GM/c^2$. The curves drawn are null, and on the axis they "
        "are also null geodesics, the paths light takes; off the axis a light ray launched along "
        "one is turned out of the plane, as it is in Kerr. With the charge, $g^{rr}$ vanishes at "
        "$r_\\pm = GM/c^2 \\pm \\sqrt{(GM/c^2)^2 - a^2 - r_Q^2}$, the cones close at both, and "
        "between them they point to smaller $r$.",
        "The domain of the chart begins at $r_+$, the outer zero of $g^{rr}$. On the axis the "
        "Kretschmann scalar stays finite at $r = 0$, because the ring singularity lies in the "
        "equatorial plane.",
    ],
    ("kerr_newman", "boyer_lindquist", "principal"): [
        "This is the equatorial plane $\\theta = \\pi/2$ drawn in $t$ and $r$, with $\\phi$ left "
        "out, for $a = 0.6\\,GM/c^2$ and $r_Q = 0.5\\,GM/c^2$. Its rays are the principal null "
        "congruence, the light rays that run straight in and straight out. As in Kerr, they are "
        "the two null directions of the plane of $\\partial_r$ and $(r^2 + a^2)\\,\\partial_t + "
        "a\\,\\partial_\\phi$, which are the repeated principal null directions of the Weyl "
        "tensor. Every ray turns as it goes, at $d\\phi/dr = \\pm a/\\Delta$ with $\\Delta = r^2 - "
        "2GMr/c^2 + a^2 + r_Q^2$, and the curves drawn are the rays' projections, $d(ct)/dr = "
        "\\pm(r^2 + a^2)/\\Delta$, which are the same at every $\\theta$ and on the axis are the "
        "radial light rays themselves.",
        "The rays are null geodesics, the paths light takes, and each cone is the future cone of "
        "the principal plane, its edges the two principal directions. The cones close at both "
        "horizons, $r_\\pm = GM/c^2 \\pm \\sqrt{(GM/c^2)^2 - a^2 - r_Q^2}$, and point to smaller "
        "$r$ between them. The dotted lines are the outer and inner ergosurfaces, where $g_{tt} = "
        "0$, at $r = GM/c^2 \\pm \\sqrt{(GM/c^2)^2 - r_Q^2}$ on the equator, $1.866$ and "
        "$0.134\\,GM/c^2$. Between them $\\partial_t$ is spacelike, and inside the inner one it is "
        "timelike again. The ingoing rays end on the ring singularity at $r = 0$, which lies in "
        "this plane.",
    ],
    ("kerr_newman", "boyer_lindquist", "above"): [
        "This is the equatorial plane $\\theta = \\pi/2$ seen from above, along the axis from "
        "$\\theta = 0$, with $r$ and $\\phi$ drawn as polar coordinates and $t$ left out, for "
        "$a = 0.6\\,GM/c^2$ and $r_Q = 0.5\\,GM/c^2$. Its rays are the principal null congruence, "
        "each turning at $d\\phi/dr = \\pm a/\\Delta$ with $\\Delta = r^2 - 2GMr/c^2 + a^2 + r_Q^2$, "
        "so both families wind counterclockwise, the way the hole turns, and far out the turning "
        "dies away as $a/r^2$. The charge enters only through $\\Delta$: it pulls the horizon in "
        "to $r_+ = 1.625\\,GM/c^2$, and the rays wind onto it without end, as $\\phi$ runs to "
        "infinity there.",
        "The winding is in the coordinate $\\phi$ alone: along an ingoing ray "
        "$\\tilde\\phi = \\phi + \\int a\\,dr/\\Delta$ stays fixed, and in it the ray crosses the "
        "horizon at a finite angle. The dotted circles are the ergosurfaces, "
        "$r = GM/c^2 \\pm \\sqrt{(GM/c^2)^2 - r_Q^2}$ on the equator, and between the outer one "
        "and $r_+$ nothing, light included, can keep $\\phi$ fixed.",
    ],
    ("kasner", "cartesian", "tx"): [
        "This is the plane of $t$ and $x$ at $y = z = 0$. Along $x$ the scale factor $t^{-2/7}$ "
        "grows toward the singularity, so the cones close up as $t \\to 0$: $dx/dt = \\pm "
        "t^{2/7}$. The Kretschmann scalar $-16p_1p_2p_3/t^4$ diverges at $t = 0$, the "
        "singularity.",
    ],
    ("kasner", "cartesian", "tz"): [
        "This is the plane of $t$ and $z$ at $x = y = 0$. Along $z$ the scale factor $t^{6/7}$ "
        "goes to zero at the singularity, and the cones open out flat: $dz/dt = \\pm t^{-6/7}$. "
        "In the plane of $t$ and $x$ the same singularity closes the cones instead, since the "
        "scale factor there, $t^{-2/7}$, grows as $t \\to 0$.",
    ],
    ("bianchi", "type_i_cartesian", "tx"): [
        "This is the plane of $t$ and $x$ at $y = z = 0$. The singularity comes about $0.378/\\bar "
        "H$ before the dashed line, and $t$ is counted from it. Near its singularity the dust "
        "universe is a Kasner spacetime, its scale factors running as powers of $t$ whose "
        "exponents lie on the Kasner circle.",
        "So along $x$ the cones close toward $t = 0$, as in Kasner's contracting direction, and "
        "open again later as $a_1$ turns round.",
    ],
    ("godel", "cartesian", "tx"): [
        "This is the plane of $t$ and $x$ at $y = z = 0$. The Gödel universe has no time function: "
        "its $g^{tt} = 2\\omega^2$ is positive, so the surfaces $t = $ const are not "
        "spacelike, and the cones are oriented instead by $\\partial_t$, which is timelike "
        "everywhere.",
        "The metric on this plane is $(-dt^2 + dx^2)/2\\omega^2$, so the curves drawn are null and "
        "run at 45°. They are not null geodesics, though. $\\Gamma^y{}_{tx} = -e^{-x}$ is not "
        "zero, so a light ray launched along one of these curves is turned out of the plane into "
        "$y$. The closed timelike curves of the Gödel universe circle each world line of the dust "
        "beyond a critical radius, through $y$ as well as $x$, so they cross this plane rather "
        "than lie in it.",
    ],
    ("alcubierre", "cartesian", "tx"): [
        "This is the plane of $t$ and $x$ on the bubble's axis of motion, $y = z = 0$, for a "
        "bubble moving at twice the speed of light along $x = 2ct$. There $\\partial_y f = "
        "\\partial_z f = 0$, so the null curves drawn are null geodesics, the paths light takes. "
        "Inside the bubble and far outside, the cones are Minkowski's; in the walls they tilt with "
        "the shift $v_s f$.",
        "The dashed lines are $g^{xx} = 1 - v_s^2 f^2 = 0$, at $v_s f = 1$. For $v_s = 2$ that is "
        "also where a forward ray keeps pace with the bubble, $f = 1 - 1/v_s$, so here they are "
        "the two horizons: forward rays from inside stall at the front wall, and forward rays from "
        "behind stall at the back wall. At other speeds the two places differ.",
    ],
    ("natario", "cartesian_flow", "tx"): [
        "This is the plane of $t$ and $x$ on the axis of motion, $y = z = 0$. There the field "
        "reduces to $u = 2nv_s = v_s f$, which is exactly Alcubierre's shift, so on this plane "
        "the two metrics are the same and so are their light rays.",
        "The drives differ only off the axis, where Natário's flow slides space sideways instead "
        "of compressing it.",
    ],
    ("krasnikov", "cylindrical", "tx"): [
        "This is the plane of $t$ and $x$ along the axis of the tube, $r = 0$. Outside the tube $k "
        "= 1$ and the cones are Minkowski's. Inside, $k$ comes close to $\\delta - 1$ and the edge "
        "moving left tips below the horizontal, so a ray going back toward $x = 0$ loses about "
        "$0.8$ in $ct$ for every unit of $x$ it covers. The edge moving right, $c\\,dt = dx$, is "
        "the same everywhere, so it points to the future inside the tube as it does outside, where "
        "$t$ is a time.",
        "The dash dot line is $g^{tt} = 0$, where $k = 0$ and the surfaces $t = $ const stop being "
        "spacelike. On the axis the expression for the Kretschmann scalar is $0/0$; the tube's "
        "curvature is concentrated in its thin walls.",
    ],
    ("pp_wave", "exact_plane_wave", "tz"): [
        "This is the plane the wave travels in, on its axis $x = y = 0$, drawn with $u = t - z$ "
        "and $v = (t + z)/2$ so that the axes read as $t$ and $z$; the chart's own $u$ and $v$ are "
        "both null. On the axis the profile $A(x^2 - y^2) + 2Bxy$ vanishes whatever $A$ and $B$ "
        "are, so the metric on this plane is flat and the rays are at 45°.",
        "Off the axis a light ray is pushed out of the plane, since $\\ddot x = (Ax + By)\\dot "
        "u^2$ and $\\ddot y = -(Ay - Bx)\\dot u^2$: the wave squeezes a beam toward the axis in "
        "one transverse direction and stretches it in the other, across this plane in $x$ and "
        "$y$. On the axis both accelerations vanish, so the null curves drawn there are null "
        "geodesics.",
    ],
    ("stockum_dust", "cylindrical", "inside"): [
        "This is the cylinder of $t$ and $\\phi$ at $r = R/2$ and $z = 0$, opened along the line "
        "$\\phi = \\pm\\pi$ and drawn with $r\\phi/R$ across, so that its left and right edges are that "
        "one line. The metric on it is $-dt^2 - (2r^2/R)\\,dt\\,d\\phi + r^2(1 - r^2/R^2)\\,d\\phi^2$, the "
        "same at every point, so its null curves are straight: $dt = r(1 - r/R)\\,d\\phi$ moving to "
        "$+\\phi$ and $dt = -r(1 + r/R)\\,d\\phi$ moving to $-\\phi$. The cross term tilts every cone "
        "toward $+\\phi$, and a curve moving that way covers three times the $\\phi$ in a given $t$ that "
        "one moving the other way does.",
        "At this radius the curve moving to $+\\phi$ is a null geodesic: $\\Gamma^r{}_{t\\phi}$ and "
        "$\\Gamma^r{}_{\\phi\\phi}$ cancel along it, and light sent that way circles the axis at "
        "$r = R/2$. The curve moving to $-\\phi$ is not a geodesic, and light launched along it is "
        "turned away from the axis. The horizontal lines, circles of constant $t$, lie outside every "
        "cone and are spacelike.",
    ],
    ("stockum_dust", "cylindrical", "beyond"): [
        "This is the cylinder of $t$ and $\\phi$ at $r = 3R/2$ and $z = 0$, opened along "
        "$\\phi = \\pm\\pi$ in the same way. Beyond $r = R$ the coefficient $g_{\\phi\\phi} = "
        "r^2(1 - r^2/R^2)$ is negative, and the cones have tipped over past the horizontal: the null "
        "curve moving to $+\\phi$, $dt = r(1 - r/R)\\,d\\phi$, goes down in $t$, while the one moving to "
        "$-\\phi$, $dt = -r(1 + r/R)\\,d\\phi$, climbs steeply. Every horizontal line, run toward "
        "$+\\phi$, points into the future cones, so the circle of constant $t$, $r$ and $z$ is a closed "
        "timelike curve.",
        "The curve moving to $+\\phi$ comes round to its own $\\phi$ at a $t$ earlier by $3\\pi R/2$ after "
        "each turn. None of the curves drawn here is a null geodesic: light launched along one moving "
        "to $+\\phi$ is turned toward the axis by $\\Gamma^r{}_{t\\phi}$ and $\\Gamma^r{}_{\\phi\\phi}$, and "
        "light launched along one moving to $-\\phi$ is turned away from it.",
    ],
    ("godel", "cylindrical", "inside"): [
        "This is the cylinder of $t$ and $\\phi$ at $r = r_c/2$ and $z = 0$, with $\\sinh r_c = 1$, about "
        "the axis $r = 0$, which is the world line of one particle of the dust. It is opened along the "
        "line $\\phi = \\pm\\pi$ and drawn with $r\\phi$ across, so that its left and right edges are that "
        "one line. The metric on it is the same at every point, so its null curves are straight: "
        "$dt = \\sinh r\\,(\\cosh r - \\sqrt{2}\\sinh r)\\,d\\phi$ moving to $+\\phi$, which is "
        "$dt = \\tfrac{1}{2}(\\sqrt{2} - 1)\\,d\\phi$ here, and $dt = -\\sinh r\\,(\\cosh r + "
        "\\sqrt{2}\\sinh r)\\,d\\phi$ moving to $-\\phi$, which is $dt = -\\tfrac{1}{2}(3 - \\sqrt{2})\\,d\\phi$. "
        "The cross term tilts every cone toward $+\\phi$, and a curve moving that way covers "
        "$1 + 2\\sqrt{2}$ times the $\\phi$ in a given $t$ that one moving the other way does.",
        "At this radius the curve moving to $+\\phi$ is a null geodesic: $\\Gamma^r{}_{t\\phi}$ and "
        "$\\Gamma^r{}_{\\phi\\phi}$ cancel along it, and light sent that way circles the axis at "
        "$r = r_c/2$. The curve moving to $-\\phi$ is not a geodesic, and light launched along it is "
        "turned away from the axis. The horizontal lines, circles of constant $t$, lie outside every "
        "cone and are spacelike.",
    ],
    ("godel", "cylindrical", "beyond"): [
        "This is the cylinder of $t$ and $\\phi$ at $r = 3r_c/2$ and $z = 0$, opened along "
        "$\\phi = \\pm\\pi$ in the same way. Beyond $r_c$ the coefficient $g_{\\phi\\phi} = "
        "2\\sinh^2 r\\,(1 - \\sinh^2 r)/\\omega^2$ is negative, and the cones have tipped over past the "
        "horizontal: the null curve moving to $+\\phi$, $dt = -\\tfrac{1}{2}(3 - \\sqrt{2})\\,d\\phi$, goes "
        "down in $t$, while the one moving to $-\\phi$, $dt = -\\tfrac{1}{2}(17 - \\sqrt{2})\\,d\\phi$, climbs "
        "steeply. Every horizontal line, run toward $+\\phi$, points into the future cones, so the circle "
        "of constant $t$, $r$ and $z$ is a closed timelike curve.",
        "The curve moving to $+\\phi$ comes round to its own $\\phi$ at a $t$ earlier by "
        "$(3 - \\sqrt{2})\\pi$ after each turn. None of the curves drawn here is a null geodesic: light "
        "launched along one moving to $+\\phi$ is turned toward the axis by $\\Gamma^r{}_{t\\phi}$ and "
        "$\\Gamma^r{}_{\\phi\\phi}$, and light launched along one moving to $-\\phi$ is turned away from it.",
    ],
    ("tov", "spherical", "radial"): [
        "This is the plane of $t$ and $r$ at $\\theta = \\pi/2$ and $\\phi = 0$ through a star of fluid "
        "with a polytrope for its equation of state, and spherical symmetry makes it the same at "
        "every other angle. Its mass and redshift functions come from $G^t{}_t$ and $G^r{}_r$ and its pressure "
        "from $\\partial_r p = -(\\rho c^2 + p)\\,\\partial_r\\Phi$, and beyond its surface, where the "
        "pressure falls to zero, the same coordinates carry on as Schwarzschild's exterior. The "
        "rays obey $c\\,dt = \\pm e^{-\\Phi}(1 - 2m/r)^{-1/2}\\,dr$.",
        "The cones are narrowest at the centre, where a clock runs at $e^\\Phi = 0.67$ of the rate "
        "$t$ counts, and they open toward the surface, where it runs at $0.84$, and on toward 45° far "
        "away. They stay open everywhere, since at $2GM/c^2R = 0.29$ the star is well inside "
        "Buchdahl's bound of $8/9$ and has no horizon.",
    ],
    ("tov", "spherical", "through"): [
        "This is the line through the centre of the star in the plane $\\theta = \\pi/2$: $x = r$ on the "
        "right is $\\phi = 0$ and $x = -r$ on the left is $\\phi = \\pi$, and spherical symmetry makes the "
        "two halves mirror images. Rays cross the centre smoothly, where the cones are narrowest and "
        "the Kretschmann scalar is finite, and the surface crosses the line on both sides.",
    ],
    ("malament_hogarth", "cartesian", "tx"): [
        "This is the plane of $t$ and $x$ at $y = z = 0$, through the removed event at the origin. "
        "The metric on it is $\\Omega^2(-c^2dt^2 + dx^2)$, and $\\Omega^2$ drops out of the null "
        "condition, so for every $\\Omega$ the light rays are Minkowski's, straight at 45°. The "
        "computer's world line runs up the axis $x = 0$ into the removed event and has no end in "
        "the spacetime.",
        "The event $p$ at $ct = 1$ on the axis has the whole of that world line in its past light "
        "cone, whose two edges run back from it, so a signal the computer sends at any moment can "
        "reach $p$ by passing round the removed event. Where $\\Omega$ grows at least as fast as "
        "$1/|t|$ toward the origin, the computer's proper time $\\int\\Omega\\,dt$ up to the removed "
        "event is infinite, while an observer who keeps away from the origin reaches $p$ in a "
        "finite time of their own.",
    ],
    ("oppenheimer_snyder", "interior_comoving", "through"): [
        "This is the line through the centre of the collapsing star in the plane $\\theta = \\pi/2$, in "
        "its own comoving coordinates: $\\chi$ on the right is $\\phi = 0$ and on the left $\\phi = \\pi$, "
        "and the surface is $\\chi_0 = \\pi/4$ on either side. The star is a closed Friedmann universe "
        "of dust, $-c^2d\\tau^2 + a^2(d\\chi^2 + \\sin^2\\chi\\,d\\Omega^2)$, released from rest at "
        "$\\tau = 0$, and its light rays obey $c\\,d\\tau = \\pm a\\,d\\chi$. As $a$ shrinks the cones open "
        "out flat, and every ray ends on the crunch at $c\\tau = \\pi a_m/2$, where the Kretschmann "
        "scalar diverges.",
        "The event horizon is the outgoing ray that leaves the centre at $c\\tau = 0.75\\,a_m$ and "
        "reaches the surface as the surface crosses $r_s = a_m\\sin^3\\chi_0$; light that leaves the "
        "centre after it never gets out. The dotted curve, $|\\nabla R|^2 = 0$ for the areal radius "
        "$R = a\\sin\\chi$, bounds the trapped spheres: it starts at the surface at the same moment "
        "and runs inward, reaching the centre only at the crunch.",
    ],
    ("oppenheimer_snyder", "exterior_schwarzschild", "radial"): [
        "This is the plane of $t$ and $r$ at $\\theta = \\pi/2$ and $\\phi = 0$ outside the collapsing "
        "star, where the metric is Schwarzschild's. The surface falls freely from rest at $r = 2r_s$ "
        "at $t = 0$, along the radial geodesic the Christoffel symbols give, and what lies inside it, "
        "$r < R(t)$, is the star, which these coordinates do not cover.",
        "The surface reaches $r_s$ only as $t \\to \\infty$, though its own clock reads a finite time "
        "there, and each outgoing ray it sends takes longer than the last to climb away. Outside the "
        "star the horizon and the black hole behind it lie beyond these coordinates; inside it the "
        "comoving coordinates carry on across both to the crunch.",
    ],
    ("tolman_bondi", "comoving_synchronous", "collapse"): [
        "This is the plane of $t$ and the comoving $r$ at $\\theta = \\pi/2$ and $\\phi = 0$, through a "
        "cloud of dust whose density falls from its centre to zero at its surface $r_b$, with vacuum "
        "outside. Every shell falls on its own clock, $R^{3/2} = r^{3/2} - \\tfrac{3}{2}"
        "\\sqrt{2GM(r)/c^2}\\,ct$, and reaches $R = 0$ at its own time: the centre first, at "
        "$ct = 0.60\\,r_b$, and the surface at $0.94\\,r_b$. The singularity, where the Kretschmann "
        "scalar diverges, is that curve, and beyond it there is no spacetime. The rays obey $c\\,dt = \\pm\\partial_r R\\,dr$, and outside the cloud the same "
        "coordinates are Lemaître's for Schwarzschild's exterior, carried by observers who fall "
        "freely from rest at infinity.",
        "The dotted curve, where $|\\nabla R|^2 = 1 - 2GM(r)/c^2R$ vanishes, bounds the trapped "
        "spheres. It meets the surface as the surface crosses $r_s = r_b/2$ and runs inward, dipping "
        "below the moment the centre is crushed, so trapped spheres form in the body of the cloud "
        "before its centre becomes singular. The event horizon leaves the centre at "
        "$ct = -0.42\\,r_b$, well before any of this, and outside the cloud it runs along the dotted "
        "curve, the sphere $R = 2GM/c^2$.",
    ],
    ("vaidya", "eddington_finkelstein_ingoing", "shell"): [
        "This is the plane of $v$ and $r$ at $\\theta = \\pi/2$ and $\\phi = 0$, drawn with $cv - r$ as "
        "the vertical axis so that the ingoing rays, $v$ constant, run at 45°. A shell of null dust of mass $M$ falls in "
        "along $v = 0$: before it $m = 0$ and the metric is flat, and after it $m = M$ and the metric "
        "is Schwarzschild's in ingoing coordinates, with $g^{rr} = 1 - r_s/r$. The outgoing rays run "
        "at 45° until the shell reaches them and bend away from $r_s$ after.",
        "The event horizon is the outgoing ray that reaches $r_s$ just as the shell does and stays "
        "there. It leaves the centre at $cv = -2r_s$, before the shell arrives, and crosses flat "
        "space, so an outgoing ray that leaves the centre after that moment ends at $r = 0$ even "
        "though its first stretch is flat. Behind the shell every future cone inside $r_s$ points to "
        "smaller $r$, and $r = 0$ is where the Kretschmann scalar $48G^2m^2/c^4r^6$ diverges.",
    ],
    ("vaidya", "eddington_finkelstein_outgoing", "shell"): [
        "This is the plane of $u$ and $r$ at $\\theta = \\pi/2$ and $\\phi = 0$, drawn with $cu + r$ as "
        "the vertical axis so that the outgoing rays, $u$ constant, run at 45°. It is the imploding shell run backward "
        "in time: a shell of null dust carries the whole mass $M$ out along $u = 0$, with "
        "Schwarzschild's metric in outgoing coordinates before it passes and flat space after.",
        "Before the shell, the region inside $r_s$ is a white hole: every future cone there points "
        "to larger $r$, and $r = 0$, where the Kretschmann scalar diverges, lies in its past. Its "
        "horizon is the ingoing ray that stays at $r_s$ until the shell leaves and then crosses flat "
        "space to the centre, arriving at $cu = 2r_s$. Every ingoing ray that reaches the centre "
        "before that moment came out of the white hole, and every one after came in from far away.",
    ],
}


# Why a coordinate system of a spacetime with no view at all has none, as prose a reader is
# given in place of the diagram. A system of a spacetime drawn elsewhere that has no view is
# accounted for in _tools/README.md instead.
NOT_DRAWN = {
    "lentz": {
        "cartesian": [
            "Lentz's soliton has no spacetime diagram. Its potential is known only as a numerical "
            "integral over the rhomboid sources Lentz laid out, so no soliton of the class can be "
            "written down to draw, and a potential written in its place would draw some other "
            "soliton of the class."],
    },
    "mixmaster": {
        "euler_angles": [
            "The Mixmaster universe has no spacetime diagram. Its scale factors are the solutions of "
            "its own vacuum field equations, which reach the singularity through an endless sequence "
            "of Kasner epochs packed ever closer toward it, so a drawing of any one solution would "
            "end on a handful of epochs and leave out the infinitely many that follow. The one plane "
            "of time and an Euler angle whose null curves are light rays, $t$ against $\\psi$ at "
            "$\\theta = \\pi/2$, has the cones $c\\,dt = \\pm a_3\\,d\\psi$ and would follow $a_3$ alone."],
    },
    "cosmic_string": {
        "interior_cap": [
            "Gott's core has no spacetime diagram of its own. On its plane of $t$ and $\\chi$ the metric "
            "is $-c^2dt^2 + \\ell^2d\\chi^2$, which is flat, and its light rays are Minkowski's. The core "
            "acts on light across that plane, in the cap of $\\chi$ and $\\phi$, "
            "$\\ell^2(d\\chi^2 + \\sin^2\\chi\\,d\\phi^2)$, a piece of a sphere whose total curvature, "
            "$2\\pi(1 - \\cos\\chi_0) = \\delta$, is the wedge the cone outside it is missing."],
    },
}


# ---------------------------------------------------------------- reading a chart

def strip_lhs(latex):
    """'K = ...' -> '...'; a bare value is returned unchanged."""
    head, eq, tail = latex.partition("=")
    return tail if eq and len(head.strip()) <= 3 else latex


def published_matrix(reader, entry, field_name):
    coords = entry["coords"]
    n = len(coords)
    g = sp.zeros(n, n)
    for comp in entry[field_name]:
        i, j = (coords.index(x) for x in comp["indices"])
        g[i, j] = reader(comp["value"])
    for i in range(n):
        for j in range(n):
            if g[i, j] == 0 and g[j, i] != 0:
                g[i, j] = g[j, i]
    return g


def load(metric_id, system_id):
    metric = json.loads((build.METRICS_DIR / f"{metric_id}.json").read_text(encoding="utf-8"))
    entry = next(e for e in metric["coordinates"] if e["id"] == system_id)
    declared = vm.DIMENSIONS[(metric_id, system_id)]
    reader = vm.Reader(entry["coords"], [p["symbol"] for p in entry.get("parameters", [])],
                       vm.time_coordinates(declared, entry["coords"]))
    return metric, entry, reader


def number(value):
    """A table value, "2/5", "pi/2" or 3, as an exact sympy number."""
    return sp.sympify(value)


# ---------------------------------------------------------------- dust

class DustSolver:
    """Scale factors from an entry's own published G^i_i set to zero, which is dust.

    Each published spatial diagonal G^i_i is read through the Reader, the declared
    functions and their first two derivatives become plain symbols, and the equations are
    solved for the second derivatives, which they are linear in. The ODE starts from
    a_i = 1 with the declared rates at a reference instant and runs back until a scale
    factor leaves (1e-7, 1e7), which is the singularity; the chart time is shifted so
    that the singularity is t = 0.
    """

    def __init__(self, metric_id, system_id, time_name, funcs, eqs, rates, params=None, start=None,
                 origin="bang"):
        _, entry, reader = load(metric_id, system_id)
        ul = entry["einstein_tensor"]["variants"]["ul"]["nonzero"]
        published = [next(c["value"] for c in ul if c["indices"] == list(ix)) for ix in eqs]
        t = reader.symbol[time_name]
        symbols = []
        exprs = [reader(text).subs(reader.c, 1) for text in published]
        for name in funcs:
            fn = reader.parameters[name]
            A0, A1, A2 = sp.symbols(f"{name}_0 {name}_1 {name}_2")
            symbols.append((A0, A1, A2))
            exprs = [e.subs(sp.Derivative(fn, (t, 2)), A2).subs(sp.Derivative(fn, t), A1).subs(fn, A0)
                     for e in exprs]
        exprs = [e.subs({reader.parameters[k]: v for k, v in (params or {}).items()}) for e in exprs]
        second = [A2 for _, _, A2 in symbols]
        solution = sp.solve(exprs, second, dict=True)[0]
        state = [A for A0, A1, _ in symbols for A in (A0, A1)]
        self.f = sp.lambdify(state, [solution[A2] for A2 in second], "numpy")
        self.funcs = list(funcs)
        n = len(funcs)

        def rhs(_, y):
            acc = self.f(*y)
            return [v for i in range(n) for v in (y[2 * i + 1], acc[i])]

        def singular(_, y):
            return min(np.min(y[0::2]) - 1e-7, 1e7 - np.max(y[0::2]))
        singular.terminal = True
        y0 = [v for a0, rate in zip(start or [1.0] * n, rates) for v in (float(a0), float(rate))]
        self.back = solve_ivp(rhs, (0, -20), y0, events=singular, rtol=1e-11, atol=1e-13, dense_output=True)
        self.fwd = solve_ivp(rhs, (0, 20), y0, events=singular, rtol=1e-11, atol=1e-13, dense_output=True)
        # The chart's time is shifted so that the singularity before the reference instant is
        # t = 0, or, with origin "reference", left with the reference instant at t = 0, as for
        # dust released from rest there and collapsing to its singularity after it.
        self.t_sing = self.back.t_events[0][0] if origin == "bang" else 0.0
        self.t_ref = -self.t_sing       # the reference instant, in time since the singularity

    def state(self, t):
        t = np.asarray(t, dtype=float)
        s = t + self.t_sing
        y = np.full((2 * len(self.funcs),) + t.shape, np.nan)
        back = (s <= 0) & (s >= self.t_sing)
        fwd = (s > 0) & (s <= self.fwd.t[-1])
        if back.any():
            y[:, back] = self.back.sol(s[back])
        if fwd.any():
            y[:, fwd] = self.fwd.sol(s[fwd])
        y[0::2] = np.where(y[0::2] > 0, y[0::2], np.nan)
        return y

    def at(self, name, x0, r):
        return self.values(name, x0)

    def values(self, name, t):
        """(a, a dot, a double dot) of one scale factor at chart times t."""
        i = self.funcs.index(name)
        y = self.state(t)
        acc = self.f(*y)
        return y[2 * i], y[2 * i + 1], np.asarray(acc[i], dtype=float) * np.ones_like(y[0])


class StarSolver:
    """A static star of fluid, its redshift and mass functions solved for a declared polytrope
    from an entry's own published Einstein tensor, G = c = M_sun = 1.

    G^t_t = -8 pi rho and G^r_r = 8 pi p, read from the published components with Phi and m
    as plain symbols, are linear in m' and Phi', and give them; the pressure follows from the
    conservation law p' = -(rho + p) Phi', which the Bianchi identity makes the same statement
    as G^theta_theta = G^r_r, so the published G^theta_theta is left for a check. The equation
    of state is p = K rho_0^2 with energy density rho = rho_0 + p. Outside the surface, where
    p = 0, m = M and e^(2 Phi) = 1 - 2M/r, the Schwarzschild exterior in the same coordinates,
    and Phi inside is shifted to meet it. conformal.py draws the same star.
    """

    funcs = ("Phi", "m")

    def __init__(self, metric_id, system_id, K, rho_c, fixed=None):
        _, entry, reader = load(metric_id, system_id)
        self.entry, self.kappa, self.rho_c = entry, K, rho_c
        ul = entry["einstein_tensor"]["variants"]["ul"]["nonzero"]
        r = reader.symbol["r"]
        by_plain = {reader._plain(n): sym for n, sym in reader.symbol.items()}
        held = {by_plain[k]: number(v) for k, v in (fixed or EQUATOR).items()}
        symbols = {}
        for name in self.funcs:
            symbols[name] = (reader.parameters[name], sp.symbols(f"_{name}0 _{name}1 _{name}2"))

        def published(index):
            e = reader(next(c["value"] for c in ul if c["indices"] == [index, index]))
            for fn, (F0, F1, F2) in symbols.values():
                e = e.subs(sp.Derivative(fn, (r, 2)), F2).subs(sp.Derivative(fn, r), F1).subs(fn, F0)
            return sp.simplify(e.subs(reader.c, 1).subs(held))
        Gtt, Grr = published("t"), published("r")
        self.Gthth = published("\\theta")
        (_, F1, F2), (M0, M1, _) = symbols["Phi"][1], symbols["m"][1]
        self.r, self.symbols = r, (F1, F2, M0, M1)
        rho, pres = sp.Symbol("rho"), sp.Symbol("p")
        solved = sp.solve([Gtt + 8 * sp.pi * rho, Grr - 8 * sp.pi * pres], [M1, F1], dict=True)[0]
        self.dm = sp.lambdify((r, M0, rho), solved[M1], "numpy")
        self.dPhi = sp.lambdify((r, M0, pres), solved[F1], "numpy")
        self.p_c = K * rho_c ** 2
        self.r0 = 1e-6

        def rhs(x, y):
            m, _, p = y
            e, _ = self.energy(p)
            f = self.dPhi(x, m, p)
            return [self.dm(x, m, e), f, -(e + max(p, 0)) * f]

        def surface(x, y):
            return y[2]
        surface.terminal = True
        r0 = self.r0
        self.ivp = solve_ivp(rhs, (r0, 100), [4 * np.pi / 3 * (rho_c + self.p_c) * r0 ** 3, 0.0, self.p_c],
                             events=surface, rtol=1e-12, atol=1e-16, dense_output=True, method="DOP853")
        self.R = float(self.ivp.t_events[0][0])
        self.M = float(self.ivp.sol(self.R)[0])
        self.shift = 0.5 * np.log(1 - 2 * self.M / self.R) - float(self.ivp.sol(self.R)[1])

    def energy(self, p):
        rho0 = np.sqrt(np.maximum(p, 0) / self.kappa)
        return rho0 + p, rho0

    def values(self, x):
        """(value, first, second derivative) of Phi and of m at radii x."""
        K, R, M, r0 = self.kappa, self.R, self.M, self.r0
        x = np.atleast_1d(np.asarray(x, dtype=float))
        inside = x < R
        out = {"Phi": [np.empty_like(x) for _ in range(3)], "m": [np.empty_like(x) for _ in range(3)]}
        if inside.any():
            xi = np.maximum(x[inside], r0)
            m, phi, p = self.ivp.sol(xi)
            p = np.maximum(p, 0)
            e, rho0 = self.energy(p)
            f = self.dPhi(xi, m, p)
            me = self.dm(xi, m, e)
            dp = -(e + p) * f
            de = -(1 + 2 * K * rho0) * f / (2 * K) + dp
            num, den = m + 4 * np.pi * xi ** 3 * p, xi * (xi - 2 * m)
            dnum = me + 12 * np.pi * xi ** 2 * p + 4 * np.pi * xi ** 3 * dp
            dden = 2 * xi - 2 * m - 2 * me * xi
            for store, v in zip(out["Phi"], (phi + self.shift, f, (dnum * den - num * dden) / den ** 2)):
                store[inside] = v
            for store, v in zip(out["m"], (m, me, 8 * np.pi * xi * e + 4 * np.pi * xi ** 2 * de)):
                store[inside] = v
        outer = ~inside
        if outer.any():
            xo = x[outer]
            for store, v in zip(out["Phi"], (0.5 * np.log(1 - 2 * M / xo), M / (xo * (xo - 2 * M)),
                                             -2 * M * (xo - M) / (xo ** 2 * (xo - 2 * M) ** 2))):
                store[outer] = v
            for store, v in zip(out["m"], (np.full_like(xo, M), 0 * xo, 0 * xo)):
                store[outer] = v
        return out

    def at(self, name, x0, r):
        """A declared function and its first two derivatives at chart points (x^0, r)."""
        shape = np.broadcast(np.asarray(x0), np.asarray(r)).shape
        return [v.reshape(shape) for v in self.values(np.broadcast_to(np.asarray(r, dtype=float), shape).ravel())[name]]

    def theta_theta(self, x):
        """The published G^theta_theta - 8 pi p along the solution, against 8 pi p_c: zero if
        the star solves the one field equation its construction did not use."""
        F1, F2, M0, M1 = self.symbols
        f = sp.lambdify((self.r, F1, F2, M0, M1), self.Gthth, "numpy")
        v = self.values(x)
        p = np.maximum(self.ivp.sol(x)[2], 0)
        return (f(x, v["Phi"][1], v["Phi"][2], v["m"][0], v["m"][1]) - 8 * np.pi * p) / (8 * np.pi * self.p_c)


_SOLVERS = {}


def star_solver(metric_id, system_id, star):
    key = (metric_id, system_id, json.dumps(star, sort_keys=True))
    if key not in _SOLVERS:
        _SOLVERS[key] = StarSolver(metric_id, system_id, float(number(star["K"])), float(number(star["rho_c"])))
    return _SOLVERS[key]


def dust_solver(metric_id, system_id, time_name, dust):
    key = (metric_id, system_id, json.dumps(dust, sort_keys=True))
    if key not in _SOLVERS:
        _SOLVERS[key] = DustSolver(metric_id, system_id, time_name, dust["funcs"], dust["eqs"],
                                   dust["rates"], dust.get("params"), dust.get("start"),
                                   dust.get("origin", "bang"))
    return _SOLVERS[key]


# ---------------------------------------------------------------- one chart, numerically

class Chart:
    """The plane of one Diagram, as numpy functions of (x^0, r)."""

    def __init__(self, spec):
        self.spec = spec
        metric, entry, reader = load(spec.metric, spec.system)
        self.metric, self.entry, self.reader = metric, entry, reader
        coords = entry["coords"]
        self.x0 = reader.symbol[spec.plane[0]]
        self.xr = reader.symbol[spec.plane[1]]
        a, b = coords.index(spec.plane[0]), coords.index(spec.plane[1])

        subs = {reader.c: 1}
        subs.update({reader.parameters[k]: number(v) for k, v in spec.params.items()})
        # The coordinates held fixed stay symbols and get their values when evaluated,
        # which is exact and never asks sympy for a limit of a long expression.
        by_plain = {reader._plain(name): symbol for name, symbol in reader.symbol.items()}
        self.fixed_syms = [by_plain[name] for name in spec.fixed]
        self.fixed_vals = [float(number(spec.fixed[name])) for name in spec.fixed]
        self.solver = (dust_solver(spec.metric, spec.system, spec.plane[0], spec.dust) if spec.dust
                       else star_solver(spec.metric, spec.system, spec.star) if spec.star else None)

        def prep(expr):
            expr = sp.sympify(expr)
            for name, rep in (spec.functions or {}).items():
                expr = expr.replace(reader.parameters[name].func, _as_lambda(reader, name, rep)).doit()
            return expr.subs(subs)

        g = published_matrix(reader, entry, "metric_components")
        gi = published_matrix(reader, entry, "inverse_metric_components")
        self.fn = {}
        # A principal view's null directions come from the principal plane instead of the
        # coordinate plane; its markers still come from the published inverse metric.
        drawn = (() if spec.principal else (("g00", g[a, a]), ("g0r", g[a, b]), ("grr", g[b, b])))
        for name, expr in drawn + (("gi00", gi[a, a]), ("gi0r", gi[a, b]), ("girr", gi[b, b])):
            self.fn[name] = self.lambdify(prep(expr))
        if spec.mark_gtt:
            self.fn["gtt"] = self.lambdify(prep(g[0, 0]))
        self.principal = PrincipalPlane(self, g, prep, a, b) if spec.principal else None
        K = prep(reader(strip_lhs(entry["kretschmann"]))) if spec.kretschmann else sp.Integer(0)
        self.fn["K"] = self.lambdify(K)
        tau = sp.sympify(spec.tau, locals={str(self.x0): self.x0, str(self.xr): self.xr})
        self.fn["dtau0"] = self.lambdify(sp.diff(tau, self.x0))
        self.fn["dtaur"] = self.lambdify(sp.diff(tau, self.xr))

        self.prep = prep
        self.surface = Surface(self, number(spec.surface)) if spec.surface else None
        self.same_as_grr = True
        if spec.areal:
            # R^2 = g_theta theta. Derivatives are taken of R^2 and divided by 2R, so that a
            # regular centre, where R^2 = r^2, is not read as a stationary point of R.
            ith = coords.index("\\theta")
            R2 = prep(g[ith, ith])
            dR2 = [sp.diff(R2, self.x0), sp.diff(R2, self.xr)]
            gi2 = [prep(gi[a, a]), prep(gi[a, b]), prep(gi[b, b])]
            grad2 = (gi2[0] * dR2[0] ** 2 + 2 * gi2[1] * dR2[0] * dR2[1] + gi2[2] * dR2[1] ** 2) / (4 * R2)
            self.fn["R"] = self.lambdify(sp.sqrt(R2))
            self.fn["dRr"] = self.lambdify(dR2[1] / (2 * sp.sqrt(R2)))
            self.fn["grad2"] = self.lambdify(grad2)
            self.same_as_grr = sp.simplify(grad2 - gi2[2]) == 0

    def lambdify(self, expr):
        """A numpy function of (x^0, r), with the fixed coordinates and any dust fed in."""
        expr = sp.sympify(expr)
        extra = []
        if self.solver:
            for name in self.solver.funcs:
                fn = self.reader.parameters[name]
                var = fn.args[0]
                A0, A1, A2 = sp.symbols(f"{name}_0 {name}_1 {name}_2")
                expr = expr.subs(sp.Derivative(fn, (var, 2)), A2).subs(sp.Derivative(fn, var), A1).subs(fn, A0)
                extra += [A0, A1, A2]
        f = sp.lambdify((self.x0, self.xr, *self.fixed_syms, *extra), expr, "numpy")

        def call(x0, r):
            x0 = np.asarray(x0, dtype=float)
            r = np.asarray(r, dtype=float)
            shape = np.broadcast(x0, r).shape
            args = [np.full(shape, v) for v in self.fixed_vals]
            if self.solver:
                for name in self.solver.funcs:
                    args += list(self.solver.at(name, x0, r))
            with np.errstate(all="ignore"):
                out = np.asarray(f(x0, r, *args))
            if np.iscomplexobj(out):
                real = np.abs(out.imag) <= 1e-9 * np.abs(out.real) + 1e-12
                out = np.where(real, out.real, np.nan)
            return np.broadcast_to(out.astype(float), shape)

        return call

    def block(self, x0, r):
        """The metric on the drawn plane, g_00, g_0r and g_rr: the coordinate plane's, or a
        principal view's principal plane in the basis whose shadow is (d_0, d_r)."""
        if self.principal:
            return self.principal.block(x0, r)
        return self.fn["g00"](x0, r), self.fn["g0r"](x0, r), self.fn["grr"](x0, r)

    def outside(self, x0, r):
        """Where the chart's published domain leaves off inside a star's surface."""
        if self.surface is None:
            return np.zeros(np.broadcast(np.asarray(x0), np.asarray(r)).shape, dtype=bool)
        with np.errstate(invalid="ignore"):
            return np.asarray(r) < self.surface(np.asarray(x0, dtype=float))

    def null_dirs(self, x0, r):
        """The directions P and M in (dx^0, dr), and D. NaN where there are none, and inside a
        star's surface, which the chart does not cover."""
        if self.surface is not None:
            P, M, D = self._null_dirs(x0, r)
            gone = self.outside(x0, r)
            return (np.where(gone[..., None], np.nan, P), np.where(gone[..., None], np.nan, M),
                    np.where(gone, np.nan, D))
        return self._null_dirs(x0, r)

    def _null_dirs(self, x0, r):
        g00, g0r, grr = self.block(x0, r)
        scale = np.maximum.reduce([np.abs(g00), np.abs(g0r), np.abs(grr)])
        with np.errstate(all="ignore"):
            g00, g0r, grr = g00 / scale, g0r / scale, grr / scale
            D = g0r ** 2 - g00 * grr
            sD = np.sqrt(np.where(D > 0, D, np.nan))
            P1, P2 = np.stack([-g0r + sD, g00], -1), np.stack([grr, -g0r - sD], -1)
            M1, M2 = np.stack([-g0r - sD, g00], -1), np.stack([grr, -g0r + sD], -1)
        P = np.where((np.linalg.norm(P1, axis=-1) >= np.linalg.norm(P2, axis=-1))[..., None], P1, P2)
        M = np.where((np.linalg.norm(M1, axis=-1) >= np.linalg.norm(M2, axis=-1))[..., None], M1, M2)
        return P, M, D

    def orient(self, x0, r, P, M):
        """Signs making P and M future directed at one point, or (0, 0) where none do."""
        at = (np.array([x0]), np.array([r]))
        g00, g0r, grr = (value[0] for value in self.block(*at))

        def dot(p, q):
            return g00 * p[0] * q[0] + g0r * (p[0] * q[1] + p[1] * q[0]) + grr * p[1] * q[1]

        mode = self.spec.orient
        if mode == "ingoing":
            sP = -np.sign(P[1])
            return sP, -np.sign(dot(sP * P, M))
        if mode == "outgoing":
            sM = np.sign(M[1])
            return -np.sign(dot(P, sM * M)), sM
        if mode == "vector":
            if g00 >= 0:
                return 0, 0
            return -np.sign(g00 * P[0] + g0r * P[1]), -np.sign(g00 * M[0] + g0r * M[1])
        d0, dr = self.fn["dtau0"](*at)[0], self.fn["dtaur"](*at)[0]
        gi00, gi0r, girr = (self.fn[k](*at)[0] for k in ("gi00", "gi0r", "girr"))
        if not gi00 * d0 ** 2 + 2 * gi0r * d0 * dr + girr * dr ** 2 < 0:
            return 0, 0
        return np.sign(P[0] * d0 + P[1] * dr), np.sign(M[0] * d0 + M[1] * dr)


class Surface:
    """The surface of a star of dust released from rest, as an exterior chart draws it: the
    radial timelike geodesic from rest at r0, integrated with the published Christoffel
    symbols in proper time and checked to keep unit speed and its energy -g_00 dx^0/dtau
    against the published metric. It starts at x^0 = 0, and R(x^0) is its radius; the chart's
    domain is what lies outside it, r >= R, and nothing is drawn inside."""

    legend = "the surface of the star, falling freely from rest"

    def __init__(self, chart, r0):
        entry, reader = chart.entry, chart.reader
        name_t, name_r = chart.spec.plane
        gamma = {tuple(c["indices"]): c["value"] for c in entry["christoffel"]["variants"]["ull"]["nonzero"]}

        def published(a, b, c):
            text = gamma.get((a, b, c))
            return chart.lambdify(chart.prep(reader(text))) if text else (lambda x0, r: 0.0)
        Gt_tr, Gr_tt, Gr_rr = published(name_t, name_t, name_r), published(name_r, name_t, name_t), \
            published(name_r, name_r, name_r)
        g00, grr = chart.fn["g00"], chart.fn["grr"]

        def one(f, t, r):
            return float(f(np.array([t]), np.array([r]))[0])

        def rhs(_, y):
            t, r, td, rd = y
            return [td, rd, -2 * one(Gt_tr, t, r) * td * rd,
                    -one(Gr_tt, t, r) * td * td - one(Gr_rr, t, r) * rd * rd]

        t_end = 2 * max(abs(v) for v in chart.spec.box)

        def far(_, y):
            return y[0] - t_end
        far.terminal = True

        def close(_, y):
            return -one(g00, y[0], y[1]) - 1e-10
        close.terminal = True
        r0 = float(r0)
        y0 = [0.0, r0, 1 / np.sqrt(-one(g00, 0.0, r0)), 0.0]
        sol = solve_ivp(rhs, (0, 1e4), y0, events=(far, close), rtol=1e-12, atol=1e-14, method="DOP853",
                        max_step=0.01)
        t, r, td, rd = sol.y
        speed = np.array([one(g00, a, b) for a, b in zip(t, r)]) * td ** 2 + \
            np.array([one(grr, a, b) for a, b in zip(t, r)]) * rd ** 2
        energy = -np.array([one(g00, a, b) for a, b in zip(t, r)]) * td
        # Against the size of the terms, which grow without bound as the surface nears r_s.
        speed = np.abs(speed + 1) / (1 + energy * td)
        if not (speed.max() < 1e-9 and np.ptp(energy) < 1e-9 and np.all(np.diff(t) > 0)):
            raise SystemExit(f"{key(chart.spec)}: the surface misses unit speed by {speed.max():.1e} "
                             f"or its energy drifts by {np.ptp(energy):.1e}")
        self.t, self.r, self.energy = t, r, float(energy[0])

    def __call__(self, x0):
        x0 = np.asarray(x0, dtype=float)
        return np.where((x0 >= 0) & (x0 <= self.t[-1]), np.interp(x0, self.t, self.r), np.nan)


def _as_lambda(reader, name, rep):
    fn = reader.parameters[name]
    body = sp.sympify(rep, locals={**reader.local, **{str(arg): arg for arg in fn.args}})
    return sp.Lambda(fn.args, body)


# ---------------------------------------------------------------- rays that leave the plane

PAIRS = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))    # the bivectors e_a ^ e_b, a < b
_PI, _PJ = np.array(PAIRS).T
TYPE_D = 1e-6       # how far, in the size of the eigenvalues, type D may miss
TANGENT = 1e-8      # how large, against the plane's own basis, a fixed coordinate's part may be
# Beside a curvature singularity the published components lose to rounding the digits the
# principal plane is read from: beside Kerr's ring, where K = 48/r^6, the directions miss
# the closed form by rounding over r^4, 1e-8 at r = 0.01 and 3e-4 at r = 0.001. So the plane
# is not taken where K passes K_END, r = 0.006 there, under a pixel from the ring, and a
# ray stops at it as it does at any point where the published metric is not finite.
K_END = 1e12


def simple(X, Y):
    """eps_abcd X^ab Y^cd / 4 for bivectors given by their six components: zero on X = Y
    exactly when X is simple, the wedge of two vectors."""
    return (X[..., 0] * Y[..., 5] + X[..., 5] * Y[..., 0] - X[..., 1] * Y[..., 4]
            - X[..., 4] * Y[..., 1] + X[..., 2] * Y[..., 3] + X[..., 3] * Y[..., 2])


class PrincipalPlane:
    """The principal plane of the published Weyl tensor at points of a view: the metric on
    it, in the basis U_0, U_r whose shadow on the drawn pair is their coordinate basis. The
    header's "Rays that leave the plane" gives the construction."""

    def __init__(self, chart, g, prep, a, b):
        spec, entry, reader = chart.spec, chart.entry, chart.reader
        coords = entry["coords"]
        self.name, self.a, self.b = key(spec), a, b
        plain = [reader._plain(c) for c in coords]
        drawn = [reader._plain(spec.plane[0]), reader._plain(spec.plane[1])]
        every = sorted(drawn + list(spec.fixed) + list(spec.leaves))
        if len(coords) != 4 or spec.dust or every != sorted(plain):
            raise SystemExit(f"{self.name}: a principal view needs four coordinates, no dust, and each "
                             "coordinate drawn, held fixed or left exactly once")
        self.fixed_index = [plain.index(name) for name in spec.fixed]
        by_plain = {reader._plain(name): symbol for name, symbol in reader.symbol.items()}
        left = {by_plain[name] for name in spec.leaves}

        weyl = entry["weyl_tensor"]["variants"]["llll"]["nonzero"]
        gamma = entry["christoffel"]["variants"]["ull"]["nonzero"]
        self.weyl_index = [tuple(coords.index(x) for x in comp["indices"]) for comp in weyl]
        self.gamma_index = [tuple(coords.index(x) for x in comp["indices"]) for comp in gamma]
        exprs = ([prep(g[i, j]) for i in range(4) for j in range(4)]
                 + [prep(reader(strip_lhs(entry["kretschmann"])))]
                 + [prep(reader(comp["value"])) for comp in weyl])
        gammas = [prep(reader(comp["value"])) for comp in gamma]
        for expr in exprs + gammas:
            if expr.free_symbols & left:
                raise SystemExit(f"{self.name}: the published metric depends on "
                                 f"{sorted(map(str, expr.free_symbols & left))}, which the rays leave "
                                 "the drawn plane in, so the plane's rays are not the shadow of one ray")
        args = (chart.x0, chart.xr, *chart.fixed_syms)
        self.f = sp.lambdify(args, exprs, "numpy")
        self.f_gamma = sp.lambdify(args, gammas, "numpy")
        self.fixed_vals = chart.fixed_vals
        self.worst = {"type D": 0.0, "tangent": 0.0}

    def _values(self, f, x0, r):
        x0, r = np.broadcast_arrays(np.asarray(x0, dtype=float), np.asarray(r, dtype=float))
        n = x0.size
        args = [np.full(n, v) for v in self.fixed_vals]
        with np.errstate(all="ignore"):
            out = f(x0.ravel(), r.ravel(), *args)
        return x0.shape, np.array([np.broadcast_to(np.asarray(v, dtype=float), (n,)) for v in out])

    def fields(self, x0, r):
        """The published g_ab, Kretschmann scalar and C_abcd at the points, stacked along the
        first axis."""
        shape, v = self._values(self.f, x0, r)
        n = v.shape[1]
        C = np.zeros((n, 4, 4, 4, 4))
        for (i, j, k, l), value in zip(self.weyl_index, v[17:]):
            C[:, i, j, k, l] = value
        return shape, v[:16].T.reshape(n, 4, 4), v[16], C

    def christoffel(self, x0, r):
        """The published Gamma^a_bc at the points, stacked along the first axis."""
        _, v = self._values(self.f_gamma, x0, r)
        G = np.zeros((v.shape[1], 4, 4, 4))
        for (i, j, k), value in zip(self.gamma_index, v):
            G[:, i, j, k] = value
        return G

    def plane(self, x0, r):
        """U, 4 by 2 at each point: its columns span the principal plane, and their shadow on
        the drawn pair is (d_0, d_r). NaN where a published quantity is not finite, or where
        the Kretschmann scalar passes K_END."""
        shape, g, K, C = self.fields(x0, r)
        U = np.full(g.shape[:1] + (4, 2), np.nan)
        good = np.isfinite(g).all((1, 2)) & np.isfinite(C).all((1, 2, 3, 4)) & (np.abs(K) < K_END)
        if good.any():
            U[good] = self._plane(g[good], C[good], np.asarray(x0), np.asarray(r))
        return shape, U, g

    def _plane(self, g, C, x0, r):
        n = g.shape[0]
        # The operator is taken in a frame orthonormal by g, e = O |Lambda|^(-1/2) from g's
        # eigenvectors, where its entries are the size of its eigenvalues; in the coordinate
        # basis, where g_theta theta = r^2 and g^theta theta = 1/r^2, they differ by powers of
        # r that cost the eigenvalues six digits beside the ring.
        lam_g, O = np.linalg.eigh(g)
        frame = O / np.sqrt(np.abs(lam_g))[:, None, :]
        eta = np.sign(lam_g)
        Cf = np.einsum("nabcd,naA,nbB,ncC,ndD->nABCD", C, frame, frame, frame, frame, optimize=True)
        Cuu = Cf * eta[:, :, None, None, None] * eta[:, None, :, None, None]
        W = Cuu[:, _PI[:, None], _PJ[:, None], _PI[None, :], _PJ[None, :]]
        lam, V = np.linalg.eig(W)
        order = np.argsort(-np.abs(lam), axis=1)
        lam = np.take_along_axis(lam, order, 1)
        V = np.take_along_axis(V, order[:, None, :], 2)
        # Type D: each of the four smaller eigenvalues is -1/2 of one of the two larger.
        miss = (np.abs(lam[:, 2:, None] + lam[:, None, :2] / 2).min(2).max(1)
                / np.abs(lam[:, :2]).max(1))
        self._require("type D", miss, TYPE_D, x0, r, "the Weyl tensor is not of Petrov type D")
        # The real span of the larger two's eigenbivectors holds l ^ n and its dual, the two
        # simple bivectors of that span, found as the null directions of eps on it.
        span = np.linalg.svd(np.concatenate([V[:, :, :2].real, V[:, :, :2].imag], 2))[0][:, :, :2]
        B1, B2 = span[..., 0], span[..., 1]
        s12 = simple(B1, B2)
        mu, Q = np.linalg.eigh(np.stack([np.stack([simple(B1, B1), s12], -1),
                                         np.stack([s12, simple(B2, B2)], -1)], -2))
        w = np.stack([np.sqrt(np.maximum(mu[:, 1], 0)), np.sqrt(np.maximum(-mu[:, 0], 0))], -1)
        planes = []
        for sign in (1.0, -1.0):
            c = np.einsum("nij,nj->ni", Q, w * np.array([1.0, sign]))
            B = c[:, :1] * B1 + c[:, 1:] * B2
            M = np.zeros((n, 4, 4))
            M[:, _PI, _PJ], M[:, _PJ, _PI] = B, -B
            E = np.einsum("nab,nbi->nai", frame, np.linalg.svd(M)[0][:, :, :2])
            planes.append((E, np.linalg.det(np.einsum("nai,nab,nbj->nij", E, g, E))))
        (E1, d1), (E2, d2) = planes
        if not np.all((d1 < 0) != (d2 < 0)):
            raise SystemExit(f"{self.name}: the principal bivectors are not one timelike plane and "
                             "one spacelike plane")
        E = np.where((d1 < 0)[:, None, None], E1, E2)
        shadow = E[:, [self.a, self.b], :]
        if np.any(np.abs(np.linalg.det(shadow)) < 1e-9):
            raise SystemExit(f"{self.name}: the principal plane does not project one to one onto the "
                             "drawn pair")
        U = E @ np.linalg.inv(shadow)
        part = (np.abs(U[:, self.fixed_index, :]).max((1, 2)) / np.abs(U).max((1, 2))
                if self.fixed_index else np.zeros(n))
        self._require("tangent", part, TANGENT, x0, r, "the principal plane leaves the surface of the "
                                                         "coordinates held fixed")
        return U

    def _require(self, what, values, tolerance, x0, r, failure):
        worst = float(np.max(values)) if values.size else 0.0
        self.worst[what] = max(self.worst[what], worst)
        if not worst <= tolerance:
            raise SystemExit(f"{self.name}: {failure}, by {worst:.1e}, among the points "
                             f"{np.ravel(x0)[:3]}, {np.ravel(r)[:3]}")

    def block(self, x0, r):
        shape, U, g = self.plane(x0, r)
        h = np.einsum("nai,nab,nbj->nij", U, g, U)
        return h[:, 0, 0].reshape(shape), h[:, 0, 1].reshape(shape), h[:, 1, 1].reshape(shape)


# ---------------------------------------------------------------- one view of a chart

class Plot:
    """A Diagram drawn: rays, cones and markers in the unit square of its axes."""

    def __init__(self, chart):
        spec = chart.spec
        self.c = chart
        self.polar = spec.to_display == POLAR
        if not self.polar:
            self.A = np.array(spec.to_display, float)
            self.Ai = np.linalg.inv(self.A)
        X0, X1, Y0, Y1 = spec.box
        self.lo = np.array([X0, Y0], float)
        self.span = np.array([X1 - X0, Y1 - Y0], float)
        self.edge = self.domain_edge() if spec.inside else None

    def to_chart(self, q):
        """Drawn axes (X, Y) to the chart's (x^0, r): linear with no offset, or for a polar
        view the angle phi in [0, 2 pi) and the radius r."""
        if self.polar:
            q = np.asarray(q, dtype=float)
            phi = np.mod(np.arctan2(q[..., 1], q[..., 0]), 2 * np.pi)
            return np.where(phi < 2 * np.pi, phi, 0.0), np.hypot(q[..., 0], q[..., 1])
        p = np.asarray(q) @ self.Ai.T
        return p[..., 0], p[..., 1]

    def from_unit(self, u):
        return np.asarray(u) * self.span + self.lo

    def push(self, k, x0, r):
        """Directions (dx^0, dr) of the chart at chart points, as directions of the drawn axes."""
        if self.polar:
            c, s = np.cos(x0)[..., None], np.sin(x0)[..., None]
            dphi, dr, r = k[..., :1], k[..., 1:], np.asarray(r)[..., None]
            return np.concatenate([c * dr - r * s * dphi, s * dr + r * c * dphi], -1)
        return k @ self.A.T

    def dirs_unit(self, u):
        """P and M at unit square points u, as unit vectors of the square, with D, P and M."""
        x0, r = self.to_chart(self.from_unit(u))
        P, M, D = self.c.null_dirs(x0, r)
        out = []
        for k in (P, M):
            d = self.push(k, x0, r) / self.span
            with np.errstate(all="ignore"):
                out.append(d / np.linalg.norm(d, axis=-1, keepdims=True))
        return out[0], out[1], D, P, M

    def domain_edge(self):
        """For inside=True, the published domain of the drawn r pulled in by a thousandth of
        the drawing: a ray is traced up to it and no further."""
        lo, hi, _, _ = parse_domains(self.c).get(self.c.spec.plane[1], (None, None, False, False))
        margin = 1e-3 * float(np.max(self.span))
        return (None if lo is None else lo + margin, None if hi is None else hi - margin)

    def beyond(self, u):
        if self.edge is None:
            return False
        _, r = self.to_chart(self.from_unit(u))
        lo, hi = self.edge
        return bool((lo is not None and r < lo) or (hi is not None and r > hi))

    # ---- rays

    def trace(self, u0, family, sign, h=0.0025, steps=5000):
        u = np.array(u0, float)
        points = [u.copy()]

        def direction(p, previous):
            d = self.dirs_unit(p[None, :])[family][0]
            if not np.all(np.isfinite(d)):
                return None
            return -d if previous is not None and np.dot(d, previous) < 0 else d

        d = direction(u, None)
        if d is None:
            return np.array(points)
        previous = sign * d
        for _ in range(steps):
            k1 = direction(u, previous)
            if k1 is None:
                break
            k2 = direction(u + 0.5 * h * k1, k1)
            if k2 is None:
                break
            k3 = direction(u + 0.5 * h * k2, k2)
            if k3 is None:
                break
            k4 = direction(u + h * k3, k3)
            if k4 is None:
                break
            step = (k1 + 2 * k2 + 2 * k3 + k4) / 6
            n = np.linalg.norm(step)
            if not np.isfinite(n) or n < 1e-6:
                break
            u = u + h * step / n
            if self.beyond(u):
                break
            previous = step / n
            points.append(u.copy())
            if np.any(u < -0.02) or np.any(u > 1.02):
                break
        return np.array(points)

    def ray_through(self, seed, family):
        return np.vstack([self.trace(seed, family, -1)[::-1], self.trace(seed, family, +1)[1:]])

    def rays(self, per_edge=15):
        """Both families, seeded evenly along the four edges, a seed passed over when a ray
        of its own family already runs within half a spacing of it."""
        if self.c.spec.ring:
            return self.ring_rays(self.c.spec.ring)
        edge = (np.arange(per_edge) + 0.5) / per_edge
        eps = 1e-3
        seeds = ([(x, eps) for x in edge] + [(eps, y) for y in edge]
                 + [(1 - eps, y) for y in edge] + [(x, 1 - eps) for x in edge])
        cell = 0.5 / per_edge
        families = {}
        for family in (0, 1):
            kept, grid = [], {}

            def near(p):
                i, j = int(p[0] / cell), int(p[1] / cell)
                return any((q[0] - p[0]) ** 2 + (q[1] - p[1]) ** 2 < cell ** 2
                           for di in (-1, 0, 1) for dj in (-1, 0, 1)
                           for q in grid.get((i + di, j + dj), ()))

            for seed in seeds:
                seed = np.array(seed)
                if near(seed):
                    continue
                line = self.ray_through(seed, family)
                inside = line[(line >= -0.02).all(1) & (line <= 1.02).all(1)]
                if len(line) < 3 or len(inside) < 3:
                    continue
                for p in line:
                    grid.setdefault((int(p[0] / cell), int(p[1] / cell)), []).append(p)
                kept.append(line)
            families[family] = kept
        if self.c.spec.mirror:
            self.join_at_axis(families)
        return families

    def ring_rays(self, n):
        """A polar view's rays, n of each family, through seeds spaced evenly around the
        circle inscribed in the box: copies of one another turned about the axis."""
        X0, X1, Y0, Y1 = self.c.spec.box
        rho = 0.98 * min(X1 - X0, Y1 - Y0) / 2
        centre = np.array([(X0 + X1) / 2, (Y0 + Y1) / 2])
        families = {}
        for family in (0, 1):
            kept = []
            for j in range(n):
                angle = 2 * np.pi * j / n
                seed = (centre + rho * np.array([np.cos(angle), np.sin(angle)]) - self.lo) / self.span
                line = self.ray_through(seed, family)
                if len(line) >= 3:
                    kept.append(line)
            families[family] = kept
        return families

    def join_at_axis(self, families, tol=0.004):
        """Every ray that meets r = 0 carries on through it as the other family, which the
        reflection draws on the far side, so a ray of that family is traced from the point."""
        for family in (0, 1):
            other = 1 - family
            for line in list(families[family]):
                for end in (line[0], line[-1]):
                    if end[0] > tol or not 0 <= end[1] <= 1:
                        continue
                    if any(np.min(np.hypot(l[:, 0] - end[0], l[:, 1] - end[1])) < tol
                           for l in families[other]):
                        continue
                    line2 = self.ray_through(np.array([1e-3, end[1]]), other)
                    if len(line2) > 3:
                        families[other].append(line2)

    # ---- cones

    def cones(self):
        """Future cones on a lattice. A point where the chart is singular gets none; a point
        where h has no null direction at all is recorded as such."""
        nx, ny = self.c.spec.cones
        out = []
        for i in range(nx):
            for j in range(ny):
                u = np.array([(i + 0.5) / nx, (j + 0.5) / ny])
                dP, dM, D, P, M = self.dirs_unit(u[None, :])
                if not np.isfinite(D[0]):
                    continue
                if not (np.all(np.isfinite(dP)) and np.all(np.isfinite(dM))):
                    if D[0] < 0:
                        out.append({"at": rounded(u), "kind": "none"})
                    continue
                x0, r = self.to_chart(self.from_unit(u))
                sP, sM = self.c.orient(float(x0), float(r), P[0], M[0])
                if sP == 0 or sM == 0:
                    raise AssertionError(f"{key(self.c.spec)}: no future at {u}; the orientation "
                                         "rule does not reach this cone")
                out.append({"at": rounded(u), "a": rounded(sP * dP[0]), "b": rounded(sM * dM[0])})
        return out

    # ---- markers

    def grid(self, n=321):
        U = np.linspace(0, 1, n)
        UU, VV = np.meshgrid(U, U)
        x0, r = self.to_chart(self.from_unit(np.stack([UU, VV], -1)))
        return UU, VV, x0, r

    def zero_set(self, name, keep=None, drop_edge=False):
        UU, VV, x0, r = self.grid()
        Z = self.c.fn[name](x0, r).astype(float)
        if keep is not None:
            Z = np.where(keep(x0, r), Z, np.nan)
        if self.c.surface is not None:
            Z = np.where(self.c.outside(x0, r), np.nan, Z)
        Z = np.where(np.isfinite(Z), Z, np.nan)
        lines = contourpy.contour_generator(UU, VV, Z, line_type="Separate").lines(0.0)
        return [rounded(thin(l, 0.0008)) for l in lines if len(l) > 3 and not (drop_edge and on_edge(l))]

    def level_sets(self, name, levels):
        UU, VV, x0, r = self.grid()
        generator = contourpy.contour_generator(UU, VV, self.c.fn[name](x0, r).astype(float),
                                                line_type="Separate")
        return [{"level": level, "lines": [rounded(thin(l, 0.0008)) for l in generator.lines(level) if len(l) > 3]}
                for level in levels]

    def singular_edges(self):
        """Edges along which the published Kretschmann scalar diverges."""
        t = np.linspace(0.01, 0.99, 99)
        edges = {"left": lambda e: np.stack([np.full_like(t, e), t], -1),
                 "right": lambda e: np.stack([np.full_like(t, 1 - e), t], -1),
                 "bottom": lambda e: np.stack([t, np.full_like(t, e)], -1),
                 "top": lambda e: np.stack([t, np.full_like(t, 1 - e)], -1)}
        out = []
        for name, at in edges.items():
            near, far = (np.abs(self.c.fn["K"](*self.to_chart(self.from_unit(at(e))))) for e in (1e-5, 1e-4))
            # A principal plane is timelike wherever it is taken, and is not taken this
            # close to a singularity.
            lorentzian = (True if self.c.principal
                          else self.c.null_dirs(*self.to_chart(self.from_unit(at(1e-4))))[2] > 0)
            with np.errstate(all="ignore"):
                if ((near > 1e8) & (near / far > 50) & lorentzian).mean() > 0.5:
                    out.append(name)
        return out

    def finite(self, x0, r):
        """Where the metric on the plane is finite."""
        with np.errstate(all="ignore"):
            return np.all([np.isfinite(v) for v in self.c.block(x0, r)], axis=0)

    def crunch_curves(self, n=241):
        """The curve past which the metric stops being finite, found up each vertical line of the
        drawing by bisection, and checked to be a curvature singularity: the Kretschmann scalar
        passes 1e8 and grows fiftyfold between 1e-4 and 1e-5 of the drawing short of it."""
        points, runs = [], []
        for u in np.linspace(0.002, 0.998, n):
            ends = [self.finite(*self.to_chart(self.from_unit(np.array([u, v])))) for v in (0.0, 1.0)]
            if not (ends[0] and not ends[1]):
                if points:
                    runs.append(points)
                points = []
                continue
            lo, hi = 0.0, 1.0
            for _ in range(60):
                mid = 0.5 * (lo + hi)
                lo, hi = (mid, hi) if self.finite(*self.to_chart(self.from_unit(np.array([u, mid])))) else (lo, mid)
            near, far = (abs(float(self.c.fn["K"](*self.to_chart(self.from_unit(np.array([u, lo - d])))))) for d in (1e-5, 1e-4))
            if not (near > 1e8 and near / far > 50):
                raise SystemExit(f"{key(self.c.spec)}: the metric stops being finite at {u, lo} where the "
                                 f"Kretschmann scalar does not diverge: {far:.1e}, {near:.1e}")
            points.append([u, lo])
        if points:
            runs.append(points)
        return [rounded(thin(np.array(p), 0.0006)) for p in runs if len(p) > 2]

    def singular_runs(self):
        """The singular_edges test point by point along each edge: an edge singular all the way
        is named, and a stretch of one is returned as a line in the unit square."""
        t = np.linspace(0.005, 0.995, 199)
        edges = {"left": lambda e: np.stack([np.full_like(t, e), t], -1),
                 "right": lambda e: np.stack([np.full_like(t, 1 - e), t], -1),
                 "bottom": lambda e: np.stack([t, np.full_like(t, e)], -1),
                 "top": lambda e: np.stack([t, np.full_like(t, 1 - e)], -1)}
        whole, runs = [], []
        for name, at in edges.items():
            near, far = (np.abs(self.c.fn["K"](*self.to_chart(self.from_unit(at(e))))) for e in (1e-5, 1e-4))
            with np.errstate(all="ignore"):
                hit = (near > 1e8) & (near / far > 50)
            if hit.all():
                whole.append(name)
                continue
            edge = at(0.0)
            for run in np.split(np.arange(t.size), np.flatnonzero(np.diff(hit.astype(int))) + 1):
                if hit[run[0]] and run.size > 2:
                    runs.append(rounded(edge[[run[0], run[-1]]]))
        return whole, runs

    def marked(self, at, family, toward=None):
        """The rays a view marks through one point: of `family`, 0 or 1, or of both, traced both
        ways or, with toward="past", only into the past. The point is an x^0, meaning the
        outermost zero of g^rr there; {"r": r, "areal": R}, the x^0 on that line where the areal
        radius is R; or {"x0": x^0, "r": r} itself."""
        fn = self.c.fn
        if isinstance(at, dict) and "areal" in at:
            r0, area = float(number(at["r"])), float(number(at["areal"]))
            x0 = root_between(lambda t: float(fn["R"](np.array([t]), np.array([r0]))[0]) - area, *self.x0_range())
        elif isinstance(at, dict):
            x0, r0 = float(number(at["x0"])), float(number(at["r"]))
        else:
            x0 = float(number(at))
            r0 = outer_root(self.c, x0)
        seed = self.to_unit(self.to_display(x0, r0))
        families = (0, 1) if family == "both" else (family,)
        if toward is None:
            return [self.ray_through(seed, f) for f in families]
        _, _, _, P, M = self.dirs_unit(seed[None, :])
        signs = self.c.orient(x0, r0, P[0], M[0])
        return [self.trace(seed, f, -signs[f]) for f in families]

    def x0_range(self):
        """The least and greatest x^0 over the drawing's box."""
        corners = self.to_chart(self.from_unit(np.array([[0, 0], [0, 1], [1, 0], [1, 1]], dtype=float)))
        return float(np.min(corners[0])), float(np.max(corners[0]))

    def surface_line(self):
        """The star's surface through the drawing, in the unit square."""
        t = np.linspace(*self.x0_range(), 2001)
        R = self.c.surface(t)
        keep = np.isfinite(R)
        return inside_unit(np.array([self.to_unit(self.to_display(a, b)) for a, b in zip(t[keep], R[keep])]))

    def to_display(self, x0, r):
        return np.asarray(self.A @ np.array([x0, r], dtype=float))

    def to_unit(self, q):
        return (np.asarray(q, dtype=float) - self.lo) / self.span

    def hatch(self):
        """Where a chart point lies outside the entry's published domains, as polygons."""
        domains = parse_domains(self.c)
        if not domains and self.c.surface is None and not self.c.spec.crunch:
            return []
        UU, VV, x0, r = self.grid(161)
        outside = np.zeros_like(UU, dtype=bool)
        for name, values in ((self.c.spec.plane[0], x0), (self.c.spec.plane[1], r)):
            if name in domains and name not in self.c.spec.periodic:
                lo, hi, lo_open, hi_open = domains[name]
                if lo is not None:
                    outside |= (values < lo) | ((values == lo) & lo_open)
                if hi is not None:
                    outside |= (values > hi) | ((values == hi) & hi_open)
        outside |= self.c.outside(x0, r)
        if self.c.spec.crunch:
            outside |= ~self.finite(x0, r)
        polygons, offsets = contourpy.contour_generator(
            UU, VV, outside.astype(float), fill_type="OuterOffset").filled(0.5, 1.5)
        out = []
        for polygon, offset in zip(polygons, offsets):
            for a, b in zip(offset[:-1], offset[1:]):
                ring = polygon[a:b]
                area = 0.5 * abs(np.dot(ring[:, 0], np.roll(ring[:, 1], 1))
                                 - np.dot(ring[:, 1], np.roll(ring[:, 0], 1)))
                if area > 1e-3:
                    out.append(rounded(thin(ring, 0.001)))
        return out

    def markers(self):
        spec, fn = self.c.spec, self.c.fn
        out = []
        if spec.mark_g00:
            lines = self.zero_set("gi00")
            if lines:
                out.append({"kind": "g00", "lines": lines})
        # Under a conformal factor g^rr vanishes only where the factor is infinite, at an event
        # the view marks itself.
        lines = [] if spec.any_factor else self.zero_set("girr")
        if lines:
            out.append({"kind": "grr", "lines": lines})
        if spec.mark_gtt:
            # Not on the ring itself, where cos(pi/2) is not zero in rounding and the published
            # g_tt reads -1: only where the Kretschmann scalar is below K_END.
            lines = self.zero_set("gtt", keep=lambda x0, r: np.abs(fn["K"](x0, r)) < K_END)
            if lines:
                t = self.c.entry["coords"][0]
                out.append({"kind": "gtt", "lines": lines, "legend": f"$g_{{{t}{t}}} = 0$, {spec.mark_gtt}"})
        if spec.areal:
            throat = self.zero_set("dRr", keep=lambda x0, r: fn["R"](x0, r) > 1e-6, drop_edge=True)
            if throat:
                out.append({"kind": "throat", "lines": throat})
            if not self.c.same_as_grr:
                apparent = self.zero_set("grad2", keep=lambda x0, r: np.abs(fn["dRr"](x0, r)) > 1e-6)
                apparent = [l for l in apparent if not same_line(l, throat)]
                if apparent:
                    out.append({"kind": "apparent", "lines": apparent})
        for kind, which, at, legend, *span in spec.lines:
            at = float(number(at))
            ends = [float(number(v)) for v in span[0]] if span else (-1e6, 1e6)
            pairs = [(at, v) for v in ends] if which == "x0" else [(v, at) for v in ends]
            line = clip_unit(np.array([self.to_unit(self.to_display(*pair)) for pair in pairs]))
            if line is not None:
                out.append({"kind": kind, "lines": [rounded(line)], "legend": legend})
        for kind, (x0, r), legend in spec.points:
            u = self.to_unit(self.to_display(float(number(x0)), float(number(r))))
            out.append({"kind": kind, "points": [rounded(u)], "legend": legend})
        if spec.surface:
            out.append({"kind": "surface", "lines": [rounded(thin(self.surface_line(), 0.0006))],
                        "legend": self.c.surface.legend})
        if spec.star:
            R = self.c.solver.R
            line = clip_unit(np.array([self.to_unit(self.to_display(x0, R)) for x0 in (-1e6, 1e6)]))
            out.append({"kind": "surface", "lines": [rounded(line)],
                        "legend": "the surface of the star, where the pressure falls to zero"})
        for kind, at, family, legend, *toward in spec.marked:
            lines = [inside_unit(line) for line in self.marked(at, family, toward[0] if toward else None)]
            out.append({"kind": kind, "lines": [rounded(thin(line, 0.0006)) for line in lines], "legend": legend})
        if spec.crunch:
            curves = self.crunch_curves()
            if curves:
                out.append({"kind": "singular", "edges": [], "lines": curves})
        elif spec.singular_runs:
            edges, runs = self.singular_runs()
            if edges or runs:
                out.append({"kind": "singular", "edges": edges, **({"lines": runs} if runs else {})})
        else:
            edges = self.singular_edges()
            if edges:
                out.append({"kind": "singular", "edges": edges})
        if spec.reference:
            y = (self.c.solver.t_ref - self.lo[1]) / self.span[1]
            out.append({"kind": "reference", "label": spec.reference, "y": round(float(y), 4)})
        return out


def rounded(points):
    return np.round(np.asarray(points, float), 4).tolist()


def on_edge(line, tol=2e-3):
    L = np.asarray(line)
    return bool(np.all(L[:, 0] < tol) or np.all(L[:, 0] > 1 - tol)
                or np.all(L[:, 1] < tol) or np.all(L[:, 1] > 1 - tol))


def same_line(a, others):
    A = np.asarray(a)
    for b in others:
        B = np.asarray(b)
        if np.median(np.min(np.linalg.norm(A[:, None, :] - B[None, :, :], axis=-1), axis=1)) < 0.01:
            return True
    return False


def thin(points, tol):
    """Ramer-Douglas-Peucker: the fewest points within tol of the curve."""
    P = np.asarray(points, float)
    if len(P) < 3:
        return P
    keep = np.zeros(len(P), bool)
    keep[0] = keep[-1] = True
    stack = [(0, len(P) - 1)]
    while stack:
        a, b = stack.pop()
        if b <= a + 1:
            continue
        chord = P[b] - P[a]
        length = np.linalg.norm(chord)
        rel = P[a + 1:b] - P[a]
        d = (np.linalg.norm(rel, axis=1) if length < 1e-12
             else np.abs(rel[:, 0] * chord[1] - rel[:, 1] * chord[0]) / length)
        i = int(np.argmax(d))
        if d[i] > tol:
            stack += [(a, a + 1 + i), (a + 1 + i, b)]
            keep[a + 1 + i] = True
    return P[keep]


def inside_unit(line):
    """A traced ray cut where it leaves the unit square, the cut made on the square's edge."""
    L = np.asarray(line, dtype=float)
    ok = np.all((L >= 0) & (L <= 1), axis=1)
    idx = np.flatnonzero(ok)
    a, b = idx[0], idx[-1]
    out = [L[a:b + 1]]
    before = clip_unit(L[[a - 1, a]]) if a > 0 else None
    after = clip_unit(L[[b, b + 1]]) if b < len(L) - 1 else None
    if before is not None:
        out.insert(0, before[:1])
    if after is not None:
        out.append(after[1:])
    return np.vstack(out)


def root_between(f, lo, hi):
    """A sign change of f on [lo, hi], by bisection on a grid and then on its bracket."""
    x = np.linspace(lo, hi, 4001)
    v = np.array([f(a) for a in x])
    i = np.flatnonzero(np.sign(v[:-1]) * np.sign(v[1:]) <= 0)
    if not i.size:
        raise SystemExit(f"no root between {lo} and {hi}")
    a, b = x[i[0]], x[i[0] + 1]
    for _ in range(100):
        m = 0.5 * (a + b)
        if np.sign(f(m)) == np.sign(f(a)):
            a = m
        else:
            b = m
    return 0.5 * (a + b)


def clip_unit(line):
    """The part of a straight segment inside the unit square, or None."""
    a, b = np.asarray(line, dtype=float)
    lo, hi = 0.0, 1.0
    d = b - a
    for k in range(2):
        if abs(d[k]) < 1e-15:
            if not 0 <= a[k] <= 1:
                return None
            continue
        s0, s1 = sorted(((0 - a[k]) / d[k], (1 - a[k]) / d[k]))
        lo, hi = max(lo, s0), min(hi, s1)
    return None if lo >= hi else np.array([a + lo * d, a + hi * d])


def outer_root(chart, x0=0.0):
    """The outermost zero of g^rr along the drawn radial range at x^0, where r_+ sits."""
    lo, hi = chart.spec.box[0], chart.spec.box[1]
    if chart.spec.to_display == POLAR:
        lo, hi = 0.0, float(np.hypot(np.max(np.abs(chart.spec.box[:2])), np.max(np.abs(chart.spec.box[2:]))))
    X = np.linspace(lo, hi, 20001)
    f = chart.fn["girr"](np.full_like(X, x0), X)
    crossings = np.where(np.sign(f[:-1]) * np.sign(f[1:]) < 0)[0]
    zeros = np.flatnonzero(f == 0)
    if not crossings.size and not zeros.size:
        return None
    if zeros.size and (not crossings.size or X[zeros[-1]] > X[crossings[-1]]):
        return float(X[zeros[-1]])
    i = crossings[-1]
    return float(X[i] - f[i] * (X[i + 1] - X[i]) / (f[i + 1] - f[i]))


DOMAIN_CONDITION = "\\;\\text{for}\\;"
CONDITION_OPERATORS = {"\\le": np.less_equal, "\\ge": np.greater_equal, "\\ne": np.not_equal,
                       "<": np.less, ">": np.greater, "=": np.equal}


def domain_holds(chart, condition, text):
    """Whether a domain's condition on the parameters, such as k \\le 0, holds at the view's
    values. A condition this cannot read, or one on a parameter the view leaves free, stops
    the run: a domain it silently passed over would draw the wrong region unhatched."""
    subs = {chart.reader.c: 1}
    subs.update({chart.reader.parameters[k]: number(v) for k, v in chart.spec.params.items()})
    for op in sorted(CONDITION_OPERATORS, key=len, reverse=True):
        lhs, found, rhs = condition.partition(op)
        if found:
            try:
                a, b = (float(sp.N(chart.reader(s.strip()).subs(subs))) for s in (lhs, rhs))
            except (vm.LatexError, TypeError, ValueError) as exc:
                raise SystemExit(f"{key(chart.spec)}: cannot evaluate the condition of the "
                                 f"domain {text!r} at {chart.spec.params}: {exc}")
            return bool(CONDITION_OPERATORS[op](a, b))
    raise SystemExit(f"{key(chart.spec)}: cannot read the condition of the domain {text!r}")


def parse_domains(chart):
    """{coordinate: (low, high, low open, high open)} for the two coordinates of the plane.

    A domain may end in a condition on the parameters, "r \\in [0, \\infty) \\;\\text{for}\\;
    k \\le 0", and is then read only for a view whose parameter values satisfy it; two domains
    for one coordinate that both hold are refused."""
    out = {}
    for text in chart.entry.get("domains", []):
        t, conditional, condition = text.partition(DOMAIN_CONDITION)
        if conditional and not domain_holds(chart, condition.strip(), text):
            continue
        t = t.replace("\\left", "").replace("\\right", "").replace("\\,", "")
        name, found, interval = t.partition("\\in")
        name, interval = name.strip(), interval.strip()
        if not found or name not in chart.spec.plane or len(interval) < 2:
            continue
        if name in out:
            raise SystemExit(f"{key(chart.spec)}: two domains for {name} hold at "
                             f"{chart.spec.params}")
        if interval[0] not in "([" or interval[-1] not in ")]":
            continue
        body, depth, cut = interval[1:-1], 0, None
        for i, ch in enumerate(body):
            depth += (ch in "({[") - (ch in ")}]")
            if ch == "," and depth == 0:
                cut = i
        if cut is None:
            continue
        ends = []
        for s in (body[:cut].strip(), body[cut + 1:].strip()):
            if s in ("\\infty", "+\\infty", "-\\infty"):
                ends.append(None)
            elif s == "r_+":
                ends.append(outer_root(chart))
            else:
                subs = {chart.reader.c: 1}
                subs.update({chart.reader.parameters[k]: number(v) for k, v in chart.spec.params.items()})
                try:
                    ends.append(float(sp.N(chart.reader(s).subs(subs))))
                except (vm.LatexError, TypeError, ValueError):
                    ends.append(None)
        out[name] = (ends[0], ends[1], interval[0] == "(", interval[-1] == ")")
    return out


# ---------------------------------------------------------------- the published file

def key(spec):
    return f"{spec.metric}/{spec.system}/{spec.view}"


def ticks(lo, hi, target):
    """Evenly spaced ticks on [lo, hi], about target of them at a step of 1, 2, 2.5 or 5
    times a power of ten, each with its label as TeX. Exact, so a tick is never 0.30000004."""
    lo, hi = Fraction(str(lo)), Fraction(str(hi))
    raw = (hi - lo) / target
    power = Fraction(10) ** math.floor(math.log10(raw))
    while power > raw:
        power /= 10
    while power * 10 <= raw:
        power *= 10
    step = next(s * power for s in (1, 2, Fraction(5, 2), 5, 10) if s * power >= raw)
    out, at = [], math.ceil(lo / step) * step
    while at <= hi:
        shown = Decimal(at.numerator) / Decimal(at.denominator)
        out.append({"at": float(at), "label": f"${shown.normalize():f}$"})
        at += step
    return out


def axes(spec):
    """The ticks of both axes. A mirrored view's horizontal axis runs from -Xmax to Xmax, and
    both axes are drawn at one scale, so the vertical axis takes as many ticks per unit."""
    x0, x1, y0, y1 = spec.box
    if spec.mirror:
        x0 = -x1
    return {"x": ticks(x0, x1, 6), "y": ticks(y0, y1, Fraction(6) * (Fraction(str(y1)) - Fraction(str(y0)))
                                             / (Fraction(str(x1)) - Fraction(str(x0))))}


def settings(spec, entry):
    """The parameter values and fixed coordinates, as LaTeX the page prints."""
    names = {vm.Reader._plain(p["symbol"].split("=")[0]): p["symbol"].split("=")[0].strip()
             for p in entry.get("parameters", [])}
    names.update({vm.Reader._plain(c): c for c in entry["coords"]})
    parts = []
    for plain, value in list(spec.params.items()) + list(spec.fixed.items()):
        shown = sp.latex(number(value), ln_notation=True)
        parts.append(f"${names[plain]} = {shown}$")
    return ", ".join(parts)


PRINCIPAL_FIELDS = ["weyl_tensor", "christoffel"]
# How far k^b nabla_b k^a may miss k^a, against the size of its two terms. Beside the horizons
# both terms grow as 1/Delta^2 and cancel, and rounding in their central differences leaves
# about 2e-6 at any step; a direction that is not a geodesic misses by order one.
GEODESIC = 1e-4
REPEATED = 1e-9     # how far C_abc[d k_e] k^b k^c may miss zero, against |C| |k|^2 |k_e|


def principal_checks(chart, n=241):
    """A principal view's rays against what the drawing does not use.

    At n points along the drawn r, away from where g^rr vanishes, each family's tangent k,
    scaled to k^r = 1, is checked to satisfy k^b nabla_b k^a = lambda k^a with the published
    Christoffel symbols, derivatives along both drawn coordinates taken by central
    differences, so that the rays are null geodesics; and C_abc[d k_e] k^b k^c = 0 with the
    published Weyl tensor, so that they are its repeated principal null directions. Returns
    the worst of each and refuses the view past GEODESIC or REPEATED.
    """
    spec, pp = chart.spec, chart.principal
    lo, hi = spec.box[0], spec.box[1]
    if spec.to_display == POLAR:
        lo, hi = 0.0, float(np.max(np.abs(spec.box)))
    r = np.linspace(lo, hi, n)[1:-1]
    x0 = np.full_like(r, 0.5)
    with np.errstate(all="ignore"):
        keep = np.abs(chart.fn["girr"](x0, r)) > 1e-2
    x0, r = x0[keep], r[keep]
    step = 1e-5 * np.maximum(1.0, np.abs(r))

    def tangents(t, rr):
        _, U, _ = pp.plane(t, rr)
        P, M, _ = chart.null_dirs(t, rr)
        return [k / k[:, pp.b:pp.b + 1] for k in (np.einsum("nai,ni->na", U, d) for d in (P, M))]

    here = tangents(x0, r)
    # d k / d x^0 and d k / d r by central differences, the other coordinates being ones the
    # metric does not depend on or ones k has no part along.
    along = [[(p - m) / (2 * step[:, None]) for p, m in zip(tangents(*plus), tangents(*minus))]
             for plus, minus in (((x0 + step, r), (x0 - step, r)), ((x0, r + step), (x0, r - step)))]
    _, g, _, C = pp.fields(x0, r)
    G = pp.christoffel(x0, r)
    geodesic, repeated, finite = [], [], np.ones(r.size, dtype=bool)
    for f, k in enumerate(here):
        dk = k[:, pp.a:pp.a + 1] * along[0][f] + k[:, pp.b:pp.b + 1] * along[1][f]
        Gkk = np.einsum("nabc,nb,nc->na", G, k, k)
        acc = dk + Gkk
        across = acc - (np.sum(acc * k, 1) / np.sum(k * k, 1))[:, None] * k
        size = np.linalg.norm(dk, axis=1) + np.linalg.norm(Gkk, axis=1)
        geodesic.append(np.linalg.norm(across, axis=1) / size)
        lower = np.einsum("nab,nb->na", g, k)
        Q = np.einsum("nabcd,nb,nc->nad", C, k, k)
        T = Q[:, :, :, None] * lower[:, None, None, :] - Q[:, :, None, :] * lower[:, None, :, None]
        scale = (np.abs(C).max((1, 2, 3, 4)) * np.abs(k).max(1) ** 2 * np.abs(lower).max(1))
        repeated.append(np.abs(T).max((1, 2, 3)) / scale)
        finite &= np.isfinite(geodesic[-1]) & np.isfinite(repeated[-1])
    # A point where the plane is not taken, beside a singularity, is left out; a NaN must never
    # pass as agreement, so the worst is taken over the finite points and too few of them fail.
    if finite.sum() < 0.8 * r.size:
        raise SystemExit(f"{key(spec)}: the principal checks are finite at {finite.sum()} of {r.size} points")
    geodesic = float(max(np.max(g[finite]) for g in geodesic))
    repeated = float(max(np.max(q[finite]) for q in repeated))
    if not (geodesic <= GEODESIC and repeated <= REPEATED):
        raise SystemExit(f"{key(spec)}: the principal rays miss the geodesic equation by {geodesic:.1e} "
                         f"or the principal condition by {repeated:.1e}")
    return {"points": int(finite.sum()), "geodesic": geodesic, "repeated": repeated}


def factor_check(spec, chart):
    """A view drawn for every value of a conformal factor: its null directions with the declared
    factor are those with the factor 1, to rounding, at points all over the drawing."""
    flat = Chart(replace(spec, functions={**spec.functions, spec.any_factor: "1"}, any_factor=None))
    u = np.random.default_rng(3).uniform(0.01, 0.99, (2000, 2))
    plot = Plot(chart)
    x0, r = plot.to_chart(plot.from_unit(u))
    gap = _same_directions(chart, flat, x0, r)
    if not gap < 1e-12:
        raise SystemExit(f"{key(spec)}: the null directions change with the conformal factor, by {gap:.1e}")


def solves_check(spec, chart):
    """Declared functions that are a solution: the published Einstein components the row names,
    mixed, vanish on them at points all over the drawing where the metric is finite."""
    ul = chart.entry["einstein_tensor"]["variants"]["ul"]["nonzero"]
    plot = Plot(chart)
    u = np.random.default_rng(5).uniform(0.01, 0.99, (1500, 2))
    x0, r = plot.to_chart(plot.from_unit(u))
    for index in spec.solves:
        value = chart.lambdify(chart.prep(chart.reader(next(c["value"] for c in ul if c["indices"] == list(index)))))
        size = chart.lambdify(chart.prep(chart.reader(next(c["value"] for c in ul
                                                           if c["indices"] == [chart.entry["coords"][0]] * 2))))
        with np.errstate(all="ignore"):
            v, scale = value(x0, r), np.abs(size(x0, r))
        ok = np.isfinite(v) & np.isfinite(scale)
        miss = float(np.max(np.abs(v[ok]) / (scale[ok] + 1)))
        if ok.sum() < 500 or not miss < 1e-9:
            raise SystemExit(f"{key(spec)}: the declared functions miss G^{index[0]}_{index[1]} = 0 by {miss:.1e} "
                             f"at {ok.sum()} finite points")


def star_checks(spec, star):
    """A declared star solves the one field equation its construction did not use, the
    published G^theta_theta = 8 pi p, and is the star its declared input says it is."""
    x = np.linspace(0.05 * star.R, 0.95 * star.R, 200)
    miss = float(np.max(np.abs(star.theta_theta(x))))
    stated = {"M": star.M, "R_km": star.R * KM}
    wrong = {k: v for k, v in stated.items() if round(v, POLYTROPE_STATED[k][1]) != POLYTROPE_STATED[k][0]}
    if not miss < 1e-7 or wrong or not 2 * star.M / star.R < 8 / 9:
        raise SystemExit(f"{key(spec)}: the star misses G^theta_theta by {miss:.1e}, or is not the star its "
                         f"input states: {wrong}")


def draw(spec):
    chart = Chart(spec)
    plot = Plot(chart)
    families = plot.rays()
    fields = (BASE_FIELDS + (["einstein_tensor"] if spec.dust or spec.star or spec.solves else [])
              + (PRINCIPAL_FIELDS if spec.principal else []))
    if spec.principal:
        principal_checks(chart)
    if spec.star:
        star_checks(spec, chart.solver)
    if spec.any_factor:
        factor_check(spec, chart)
    if spec.solves:
        solves_check(spec, chart)
    view = {
        "id": spec.view, "label": spec.label, "plane": list(spec.plane), "families": list(spec.families),
        "xlabel": spec.xlabel, "ylabel": spec.ylabel, "ticks": axes(spec), "box": list(spec.box),
        "to_display": POLAR if plot.polar else [list(map(float, row)) for row in spec.to_display],
        "mirror": spec.mirror,
        "rays": {name: [rounded(thin(l, 0.0006)) for l in families[i]] for i, name in ((0, "P"), (1, "M"))},
        "cones": plot.cones(), **({"cone": spec.cone} if spec.cone else {}),
        "markers": plot.markers(), "hatch": plot.hatch(),
        "settings": settings(spec, chart.entry), "input": spec.input,
        "caption": CAPTIONS[(spec.metric, spec.system, spec.view)],
        "source": {"fields": fields, "version": build.diagram_source_version(chart.entry, fields)},
    }
    if spec.areal_contours:
        view["areal"] = plot.level_sets("R", spec.areal_contours)
    return view


NONE_FIELDS = ["coords", "parameters", "metric_components"]


def not_drawn(metric_id):
    """The reasons a spacetime's coordinate systems carry no view, as its file carries them,
    each stamped with the fields it speaks of so that a changed metric asks for it again."""
    metric = json.loads((build.METRICS_DIR / f"{metric_id}.json").read_text(encoding="utf-8"))
    systems = {s["id"]: s for s in metric["coordinates"]}
    return {system_id: {"text": text, "source": {"fields": NONE_FIELDS, "version": build.diagram_source_version(
                systems[system_id], NONE_FIELDS)}}
            for system_id, text in NOT_DRAWN[metric_id].items()}


def by_metric():
    """Every spacetime with a diagram file: its flat views and its figures, in table order."""
    import projections as pj
    out = {}
    for spec in DIAGRAMS:
        out.setdefault(spec.metric, ([], []))[0].append(spec)
    for spec in pj.FIGURES:
        out.setdefault(spec.metric, ([], []))[1].append(spec)
    for metric_id in NOT_DRAWN:
        out.setdefault(metric_id, ([], []))
    return out


def write(metric_ids=None):
    import projections as pj
    if not DIAGRAMS_DIR.exists():
        DIAGRAMS_DIR.mkdir(parents=True)
    for metric_id, (specs, figures) in by_metric().items():
        if metric_ids and metric_id not in metric_ids:
            continue
        systems, projections = {}, {}
        for spec in specs:
            systems.setdefault(spec.system, []).append(draw(spec))
            print(f"  {key(spec)}", flush=True)
        for spec in figures:
            projections.setdefault(spec.system, []).append(pj.draw(spec))
            print(f"  {pj.key(spec)}", flush=True)
        data = {"metric": metric_id, "systems": systems}
        if projections:
            data["projections"] = projections
        if metric_id in NOT_DRAWN:
            data["none"] = not_drawn(metric_id)
        path = DIAGRAMS_DIR / f"{metric_id}.json"
        path.write_text(json.dumps(data, ensure_ascii=False, separators=(",", ":"), allow_nan=False) + "\n",
                        encoding="utf-8")
        print(f"wrote {path.relative_to(build.ROOT)}")


# ---------------------------------------------------------------- closed forms

def _rstar(r, horizons):
    """The tortoise coordinate of f = prod(1 - r_i/r) with simple roots r_i, up to a constant."""
    out = np.asarray(r, float).copy()
    for i, ri in enumerate(horizons):
        others = [rj for j, rj in enumerate(horizons) if j != i]
        coefficient = ri ** len(horizons) / np.prod([ri - rj for rj in others]) if others else ri
        out = out + coefficient * np.log(np.abs(r - ri))
    return out


# (metric, system, view): (what P conserves, what M conserves, where to compare). None
# where a family has no closed form. P moves toward smaller r or x, M toward larger.
CLOSED_FORMS = {
    ("schwarzschild", "spherical", "radial"):
        (lambda t, r: t + _rstar(r, [1]), lambda t, r: t - _rstar(r, [1]), lambda t, r: np.abs(r - 1) > 0.05),
    ("schwarzschild", "eddington_finkelstein_ingoing", "finkelstein"):
        (lambda v, r: v, lambda v, r: v - 2 * _rstar(r, [1]), lambda v, r: np.abs(r - 1) > 0.05),
    ("schwarzschild", "eddington_finkelstein_outgoing", "finkelstein"):
        (lambda u, r: u + 2 * _rstar(r, [1]), lambda u, r: u, lambda u, r: np.abs(r - 1) > 0.05),
    ("rn_metric", "spherical", "radial"):
        (lambda t, r: t + _rstar(r, [0.8, 0.2]), lambda t, r: t - _rstar(r, [0.8, 0.2]),
         lambda t, r: (np.abs(r - 0.8) > 0.05) & (np.abs(r - 0.2) > 0.05)),
    ("de_sitter", "static_spherical", "radial"):
        (lambda t, r: t + 0.5 * np.log(np.abs((1 + r) / (1 - r))),
         lambda t, r: t - 0.5 * np.log(np.abs((1 + r) / (1 - r))), lambda t, r: np.abs(r - 1) > 0.05),
    ("anti_de_sitter", "static_global", "radial"):
        (lambda t, r: t + np.arctan(r), lambda t, r: t - np.arctan(r), None),
    ("ellis_bronnikov", "spherical", "radial"): (lambda t, r: t + r, lambda t, r: t - r, None),
    ("morris_thorne", "spherical", "radial"):
        (lambda t, r: t + np.sqrt(r ** 2 - 1), lambda t, r: t - np.sqrt(r ** 2 - 1), lambda t, r: r > 1.0005),
    ("minkowski", "spherical", "radial"): (lambda t, r: t + r, lambda t, r: t - r, None),
    ("minkowski", "cartesian", "tx"): (lambda t, x: t + x, lambda t, x: t - x, None),
    ("minkowski", "rindler", "tx"):
        (lambda T, X: T + np.log(X), lambda T, X: T - np.log(X), lambda T, X: X > 1e-3),
    ("frw", "conformal_spherical", "radial"): (lambda e, r: e + r, lambda e, r: e - r, None),
    ("de_sitter", "flat_slicing", "tx"): (lambda t, x: x - np.exp(-t), lambda t, x: x + np.exp(-t), None),
    ("anti_de_sitter", "poincare", "tx"): (lambda t, x: t + x, lambda t, x: t - x, None),
    ("bertotti_robinson", "poincare", "tx"): (lambda t, x: t + x, lambda t, x: t - x, None),
    ("godel", "cartesian", "tx"): (lambda t, x: t + x, lambda t, x: t - x, None),
    ("kasner", "cartesian", "tx"):
        (lambda t, x: x + t ** (9 / 7) * 7 / 9, lambda t, x: x - t ** (9 / 7) * 7 / 9, lambda t, x: t > 1e-3),
    ("kasner", "cartesian", "tz"):
        (lambda t, z: z + 7 * t ** (1 / 7), lambda t, z: z - 7 * t ** (1 / 7), lambda t, z: t > 1e-3),
    ("pp_wave", "exact_plane_wave", "tz"): (lambda u, v: v, lambda u, v: u, None),
    ("krasnikov", "cylindrical", "tx"): (None, lambda t, x: t - x, None),
}
# The cylinders of t and phi, where each family runs straight, dt = k dphi, and conserves
# t - k phi: van Stockum's k = r(1 - r/R) moving right and -r(1 + r/R) moving left, and
# Godel's sinh r (cosh r - sqrt 2 sinh r) and -sinh r (cosh r + sqrt 2 sinh r), which are the
# numbers the captions state.
CYLINDERS = {
    ("stockum_dust", "cylindrical", "inside"): (-0.75, 0.25),
    ("stockum_dust", "cylindrical", "beyond"): (-3.75, -0.75),
    ("godel", "cylindrical", "inside"): (-(3 - math.sqrt(2)) / 2, (math.sqrt(2) - 1) / 2),
    ("godel", "cylindrical", "beyond"): (-(17 - math.sqrt(2)) / 2, -(3 - math.sqrt(2)) / 2),
}
CLOSED_FORMS.update({where: (lambda t, phi, k=left: t - k * phi, lambda t, phi, k=right: t - k * phi, None)
                     for where, (left, right) in CYLINDERS.items()})
# What light launched along each family of a cylinder does, as its caption says: stays on the
# cylinder as a null geodesic, or is turned toward or away from the axis.
TURNING = {
    ("stockum_dust", "cylindrical", "inside"): ("away", "geodesic"),
    ("stockum_dust", "cylindrical", "beyond"): ("away", "toward"),
    ("godel", "cylindrical", "inside"): ("away", "geodesic"),
    ("godel", "cylindrical", "beyond"): ("away", "toward"),
}


def _kerr_forms(a, rQ=0.0):
    """Kerr's and Kerr-Newman's principal null rays, with M = 1: t -+ r_* and phi -+ r_# are
    conserved on the ingoing and outgoing rays, dr_*/dr = (r^2 + a^2)/Delta and
    dr_#/dr = a/Delta, Delta = (r - r_+)(r - r_-)."""
    rp, rm = 1 + math.sqrt(1 - a * a - rQ * rQ), 1 - math.sqrt(1 - a * a - rQ * rQ)

    def rstar(r):
        return (r + (rp * rp + a * a) / (rp - rm) * np.log(np.abs(r - rp))
                - (rm * rm + a * a) / (rp - rm) * np.log(np.abs(r - rm)))

    def rsharp(r):
        return a / (rp - rm) * np.log(np.abs((r - rp) / (r - rm)))

    def away(x, r):
        return (np.abs(r - rp) > 0.05) & (np.abs(r - rm) > 0.05) & (r > 0.05)
    return ((lambda t, r: t + rstar(r), lambda t, r: t - rstar(r), away),
            (lambda phi, r: phi + rsharp(r), lambda phi, r: phi - rsharp(r), away))


for _metric, _forms in (("kerr", _kerr_forms(0.9)), ("kerr_newman", _kerr_forms(0.6, 0.5))):
    CLOSED_FORMS[(_metric, "boyer_lindquist", "principal")] = _forms[0]
    CLOSED_FORMS[(_metric, "boyer_lindquist", "above")] = _forms[1]


def verify():
    """Trace rays as the page does and measure how far each family's closed form drifts.
    Returns the number of failures."""
    failures = 0
    specs = {(s.metric, s.system, s.view): s for s in DIAGRAMS}
    frw = specs[("frw", "comoving_spherical", "radial")]
    solver = Chart(frw).solver

    def eta(t):
        return np.array([quad(lambda s: 1 / solver.values("a", np.array([s]))[0][0], 1e-12, x, limit=400)[0]
                         for x in np.atleast_1d(t)])
    forms = dict(CLOSED_FORMS)
    forms[("frw", "comoving_spherical", "radial")] = (lambda t, r: eta(t) + r, lambda t, r: eta(t) - r,
                                                      lambda t, r: t > 0.02)
    print(f"{'view':56s} {'P drift':>9s} {'M drift':>9s}  other family spread")
    traced = {}
    for where, (own_P, own_M, keep) in forms.items():
        plot = Plot(Chart(specs[where]))
        traced[where] = plot.c
        drift, spread = {0: 0.0, 1: 0.0}, {0: 0.0, 1: 0.0}
        for s in np.linspace(0.05, 0.95, 7):
            for seed in ((s, 0.001), (0.999, s), (s, 0.999), (0.001, s)):
                for family, own, other in ((0, own_P, own_M), (1, own_M, own_P)):
                    if own is None:
                        continue
                    line = plot.ray_through(np.array(seed), family)
                    line = line[(line >= 0).all(1) & (line <= 1).all(1)]
                    x0, r = plot.to_chart(plot.from_unit(line))
                    if plot.polar:
                        x0 = np.unwrap(x0)
                    mask = np.ones_like(r, bool) if keep is None else keep(x0, r)
                    if mask.sum() < 10:
                        continue
                    drift[family] = max(drift[family], float(np.ptp(own(x0[mask], r[mask]))))
                    if other is not None:
                        spread[family] = max(spread[family], float(np.ptp(other(x0[mask], r[mask]))))
        # A family conserves its own quantity and not the other's, which is what
        # confirms that P and M are labelled the way the page colours them.
        worst = max(drift.values())
        if own_P is None or own_M is None:
            separated = math.nan  # one family alone has a closed form: nothing to tell apart
        else:
            separated = min(spread.values())
        ok = worst < 1e-5 and (math.isnan(separated) or separated > 1.0)
        failures += not ok
        print(f"{'/'.join(where):56s} {drift[0]:9.1e} {drift[1]:9.1e}  {separated:6.2f}  {'ok' if ok else 'FAILED'}")
    t_sing = solver.t_sing
    t = np.linspace(0.01, 2.0, 400)
    eds = float(np.max(np.abs(solver.values("a", t)[0] / (t / -t_sing) ** (2 / 3) - 1)))
    ok = abs(t_sing + 2 / 3) < 1e-8 and eds < 1e-7
    failures += not ok
    print(f"FRW dust: the bang {t_sing:.9f} from a = 1, Einstein-de Sitter -2/3; "
          f"a(t) against (t/t0)^(2/3) to {eds:.1e}  {'ok' if ok else 'FAILED'}")
    bianchi = Chart(specs[("bianchi", "type_i_cartesian", "tx")]).solver
    y = bianchi.state(np.array([1e-4, 2e-4]))
    p = [float(np.log(y[2 * i][1] / y[2 * i][0]) / np.log(2)) for i in range(3)]
    ok = abs(sum(p) - 1) < 1e-3 and abs(sum(q * q for q in p) - 1) < 1e-3
    failures += not ok
    print(f"Bianchi I dust: exponents at the singularity {', '.join(f'{q:.4f}' for q in p)}, "
          f"on the Kasner circle  {'ok' if ok else 'FAILED'}")
    alcubierre = Chart(specs[("alcubierre", "cartesian", "tx")])
    natario = Chart(specs[("natario", "cartesian_flow", "tx")])
    T, X = np.meshgrid(np.linspace(-2, 2, 201), np.linspace(-3, 3, 301))
    gap = max(float(np.nanmax(np.abs(alcubierre.fn[k](T, X) - natario.fn[k](T, X)))) for k in ("g00", "g0r", "grr"))
    ok = gap < 1e-12
    failures += not ok
    print(f"Natario against Alcubierre on the axis: the metrics on the plane differ by {gap:.1e}  {'ok' if ok else 'FAILED'}")
    return failures + verify_turning(specs) + verify_principal(specs, traced)


def verify_turning(specs):
    """On each cylinder of t and phi, the acceleration off the cylinder of light launched along
    each family, -Gamma^r_ab k^a k^b from the published Christoffel symbols, against what the
    caption says it is. Returns the number of failures."""
    failures = 0
    print()
    for where, said in TURNING.items():
        spec = specs[where]
        _, entry, reader = load(spec.metric, spec.system)
        coords = entry["coords"]
        at = {reader.c: 1}
        at.update({reader.parameters[k]: number(v) for k, v in spec.params.items()})
        by_plain = {reader._plain(n): sym for n, sym in reader.symbol.items()}
        at.update({by_plain[k]: number(v) for k, v in spec.fixed.items()})
        at.update({reader.symbol[c]: 0 for c in spec.plane})
        g = np.array(published_matrix(reader, entry, "metric_components").subs(at), dtype=float)
        gamma = {tuple(c["indices"]): float(reader(c["value"]).subs(at))
                 for c in entry["christoffel"]["variants"]["ull"]["nonzero"] if c["indices"][0] == "r"}
        i, j = (coords.index(c) for c in spec.plane)
        found = []
        for sign in (-1, 1):
            # (dt, dphi) = (lam, sign), null on the cylinder and future along d_t.
            a, b, c = g[i, i], 2 * g[i, j] * sign, g[j, j]
            roots = [(-b + e * math.sqrt(b * b - 4 * a * c)) / (2 * a) for e in (-1, 1)]
            lam = next(x for x in roots if g[i, i] * x + g[i, j] * sign < 0)
            k = {spec.plane[0]: lam, spec.plane[1]: sign}
            terms = [gamma.get((("r",) + (m, n)), 0.0) * k[m] * k[n] for m in k for n in k]
            accel = -sum(terms)
            scale = sum(abs(x) for x in terms)
            found.append("geodesic" if abs(accel) < 1e-12 * scale else "toward" if accel < 0 else "away")
        ok = tuple(found) == said
        failures += not ok
        print(f"{'/'.join(where):56s} light moving left turned {found[0]}, moving right {found[1]}, "
              f"as the caption says  {'ok' if ok else 'FAILED'}")
    return failures


def _sine(A, B):
    """The sine of the angle between two stacks of plane directions, blind to their sign."""
    with np.errstate(all="ignore"):
        return np.abs(A[..., 0] * B[..., 1] - A[..., 1] * B[..., 0]) / (
            np.linalg.norm(A, axis=-1) * np.linalg.norm(B, axis=-1))


def _same_directions(first, second, x0, r):
    """The largest angle, as a sine, between the P of two charts and between their M."""
    one, two = first.null_dirs(x0, r), second.null_dirs(x0, r)
    return float(np.nanmax([_sine(one[0], two[0]), _sine(one[1], two[1])]))


def verify_principal(specs, traced):
    """The principal views' own checks where their rays ran, and the method run where it has
    nothing new to draw. Returns the number of failures."""
    failures = 0
    print()
    for spec in DIAGRAMS:
        if not spec.principal:
            continue
        chart = traced[(spec.metric, spec.system, spec.view)]
        c = principal_checks(chart)
        w = chart.principal.worst
        print(f"{key(spec):56s} type D to {w['type D']:.0e}, tangent to {w['tangent']:.0e} where "
              f"traced; at {c['points']} points geodesic to {c['geodesic']:.0e}, principal to "
              f"{c['repeated']:.0e}  ok")
    r = np.linspace(0.07, 3.97, 40)
    t = np.full_like(r, 0.3)
    for where in (("schwarzschild", "spherical", "radial"), ("rn_metric", "spherical", "radial"),
                  ("taub_nut", "spherical", "radial")):
        spec = specs[where]
        box = spec.box
        T, R = np.meshgrid(np.linspace(box[2], box[3], 9)[1:-1], np.linspace(box[0], box[1], 41)[1:-1])
        gap = _same_directions(Chart(spec), Chart(replace(spec, principal=True)), T, R)
        ok = gap < 1e-9
        failures += not ok
        print(f"{'/'.join(where):56s} the principal plane is the plane of t and r drawn: "
              f"directions agree to {gap:.0e}  {'ok' if ok else 'FAILED'}")
    for metric in ("kerr", "kerr_newman"):
        spec = specs[(metric, "boyer_lindquist", "principal")]
        at_equator = Chart(spec)
        gap_theta = _same_directions(at_equator, Chart(replace(spec, fixed={"theta": "pi/5"})), t, r)
        gap_axis = _same_directions(at_equator, Chart(specs[(metric, "boyer_lindquist", "radial")]), t, r)
        ok = gap_theta < 1e-9 and gap_axis < 1e-9
        failures += not ok
        print(f"{metric + '/boyer_lindquist/principal':56s} in t and r the same at theta = pi/5 to "
              f"{gap_theta:.0e}, and on the axis to {gap_axis:.0e}  {'ok' if ok else 'FAILED'}")
    stockum = Diagram("stockum_dust", "cylindrical", "principal", "", ("t", "r"), (0, 2, -1, 1), "", "",
                      {"R": 1}, {"z": "0"}, principal=True, leaves=("phi",))
    try:
        Chart(stockum).null_dirs(np.zeros(5), np.linspace(0.2, 1.8, 5))
        refused = ""
    except SystemExit as exc:
        refused = str(exc)
    ok = "Petrov type D" in refused
    failures += not ok
    print(f"{'stockum_dust/cylindrical':56s} refused, its Weyl tensor of type I: {refused or 'NOT refused'}"
          f"  {'ok' if ok else 'FAILED'}")
    return failures


def check_table():
    """Every view and figure has a caption, no two share a place, and a coordinate system with
    a stated reason for drawing nothing draws nothing."""
    import projections as pj
    flat = [(spec.metric, spec.system, spec.view) for spec in DIAGRAMS]
    figures = [(spec.metric, spec.system, spec.view) for spec in pj.FIGURES]
    if len(set(flat + figures)) != len(flat + figures):
        raise SystemExit("two views share a metric, system and view id")
    for table, places, captions in (("DIAGRAMS", flat, CAPTIONS), ("FIGURES", figures, pj.CAPTIONS)):
        missing = ["/".join(place) for place in places if place not in captions]
        stray = set(captions) - set(places)
        if missing or stray:
            raise SystemExit(f"{table}: views without a caption: {missing}; captions without a view: "
                             f"{sorted(stray)}")
    drawn = {(metric, system) for metric, system, _ in flat + figures}
    both = [f"{metric}/{system}" for metric, systems in NOT_DRAWN.items() for system in systems
            if (metric, system) in drawn]
    if both:
        raise SystemExit(f"NOT_DRAWN gives a reason for systems that are drawn: {both}")


def main(argv=None):
    check_table()
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--metric", action="append", default=[], help="redraw only this spacetime, repeatable")
    parser.add_argument("--verify", action="store_true",
                        help="check the rays against closed forms instead of writing")
    args = parser.parse_args(argv)
    unknown = set(args.metric) - set(by_metric())
    if unknown:
        parser.error(f"no diagram is drawn for {sorted(unknown)}")
    if args.verify:
        return 1 if verify() else 0
    write(set(args.metric))
    return 0


if __name__ == "__main__":
    sys.exit(main())
