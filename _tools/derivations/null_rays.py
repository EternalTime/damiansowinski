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
diagram, removing the file of any spacetime that has none, and the second command stamps
each file's version into the index. Pass
--metric <metric_id> to redraw one spacetime, repeatable, and --verify to check the
rays against closed forms the drawing never uses instead of writing anything, with --metric
checking only those spacetimes' closed forms. --slices
rewrites only the slices below, from the embedding files as they stand, tracing no ray,
which is what to run after an embedding diagram moves its moments.


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
left moving family, and M the retarded one. Where the drawn coordinate falls outward, as
the C-metric's Hong-Teo y = 1/(alpha r) does, the left moving family is the outgoing one and
P conserves the retarded coordinate, so the page's colour of P means moving left in every
view. A view drawn against a time that runs down the chart, as Gowdy's -tau is, mirrors the
plane, and there the chart's M is the family drawn moving left and is the one named P.
--verify confirms that labelling chart by chart. Their integral curves are traced
both ways by fourth order Runge-Kutta in the drawing's own unit square, from seeds spaced
evenly along its four edges.

Every curve drawn is a null curve, its tangent null. It is also a null geodesic, the path
a light ray takes, only if nothing accelerates it out of the plane, Gamma^A_ab k^a k^b = 0
for every fixed coordinate A. The planes in DIAGRAMS pass that test, except Godel's and Som and
Raychaudhuri's planes of t and x, whose captions say their null curves are not null geodesics.


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


Rays of no angular momentum
---------------------------

The rotating BTZ hole turns every light ray that runs straight in or out in phi, as Kerr's
does, and between its ergosurface and its outer horizon the plane of t and r at fixed phi
has no null direction at all. Its Weyl tensor vanishes, as every Weyl tensor in three
dimensions does, so it has no principal directions to draw them by. A row with
quotient='phi' divides the circles of phi out instead: where the published metric does not
depend on the coordinate k named, the plane's metric is

    h_ab = g_ab - g_ak g_bk / g_kk,

the part of g orthogonal to the circles, and it takes the place of the coordinate plane's
metric in the null condition, as the principal plane's does. A null direction of h, lifted
by dx^k = -(g_ka dx^a + g_kb dx^b)/g_kk, is null in the spacetime and carries no momentum
along k, and the rays of h are the shadows on the drawn pair of the null geodesics with no
angular momentum: those project to the geodesics of h, and in two dimensions every null
curve of h is one. The upper left block of the published inverse metric is the inverse of
h, so the markers read from it mark the zeros of h's own g^rr. Every coordinate is drawn,
held fixed or divided out, and nothing published may depend on the one divided out.

Before such a view is written, each family's lifted direction is checked null against the
published metric and to be a geodesic by the published Christoffel symbols, the orbits of
k are checked spacelike, and the Christoffel symbols are stamped with the view.


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
  'split'     'ingoing' below the radius the row names in `split`, and 'outgoing' above it,
              for a static region between a black hole horizon and a cosmological one: both
              agree with t inside that region, and together they carry the orientation
              through each horizon into the black hole and into the expanding region beyond.
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
              A row may declare a curve singular in `singular_zero` where the scalar
              diverges too slowly for that; the same test is then taken at 1e-20 and 1e-30
              of the chart's unit from the curve, in 60 digits.
  hatch       outside the entry's published domains, parsed by the same Reader; a domain
              ending at the undeclared r_+ is read at the outermost zero of g^rr, and one from
              r_+ to the undeclared r_c at the two outermost zeros.
              A domain that is an inequality between the plane's two coordinates, as UV < 1,
              is stated in the row's `where`, positive on the spacetime: what lies beyond
              is hatched, and no ray or cone is drawn there.



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


Slices
------

A view on which the spacetime's embedding diagram is cut from a moment that meets its plane
carries that moment under `slices`, one entry for each surface of each embedding view it
shows: the embedding view's id and the surface's place, the moment's label as TeX, a stamp
over the embedding surface, and the moment as lines, points and regions in the view's unit
square. slices.py declares what each moment is in the view's chart, and view_slices()
carries it into the square with the view's own map and box.
"""

import argparse
import json
import math
import sys
from dataclasses import dataclass, field, replace
from decimal import Decimal
from fractions import Fraction
from functools import lru_cache
from pathlib import Path

import contourpy
import numpy as np
import sympy as sp
from scipy.integrate import quad, solve_ivp
from scipy.special import erf as scipy_erf
from scipy.special import expi as scipy_expi

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import build_mfs_data as build  # noqa: E402
import ori_shell  # noqa: E402
import slices  # noqa: E402
import verify_metrics as vm  # noqa: E402

DIAGRAMS_DIR = build.DIAGRAMS_DIR
BASE_FIELDS = ["coords", "parameters", "metric_components", "inverse_metric_components",
               "kretschmann", "domains"]

FINKELSTEIN_IN = ((0, 1), (1, -1))      # X = r, Y = v - r
FINKELSTEIN_OUT = ((0, 1), (1, 1))      # X = r, Y = u + r
NULL_TO_TR = ((-0.5, 0.5), (0.5, 0.5))  # (u, v) of Minkowski to X = (v - u)/2, Y = (u + v)/2
NULL_TO_SUM = ((-1, 1), (1, 1))         # (u, v) of Khan and Penrose to X = v - u, Y = u + v
UV_TO_TZ = ((-0.5, 1), (0.5, 1))        # (u, v) of the plane wave to z = v - u/2, t = v + u/2
RADIAL = ("ingoing", "outgoing")
SIDEWAYS = ("moving left", "moving right")
EQUATOR = {"theta": "pi/2", "phi": "0"}
# Tangherlini's planes of the time and r hold every angle fixed, three in five dimensions and four in six.
TANGHERLINI_FIVE = {"psi": "pi/2", "theta": "pi/2", "phi": "0"}
TANGHERLINI_SIX = {"chi": "pi/2", **TANGHERLINI_FIVE}
# Boulware and Deser's black hole in units of its horizon radius, r_h^2 = r_0^2 - l^2 = 1, and the
# other branch in units of l at r_0 = l.
BD = {"r_0": "13/12", "ell": "5/12"}
BD_PLUS = {"r_0": 1, "ell": 1}
# The black string's planes of the time and r hold the angles and the length z along the string fixed.
STRING_FIVE = {"theta": "pi/2", "phi": "0", "z": "0"}
STRING_SIX = {"psi": "pi/2", **STRING_FIVE}
# The Kaluza-Klein monopole's planes of t and its radius lie on the half axis theta = 0, where Gross and
# Perry's potential vanishes and the embedding diagram's cigar stands.
KK_AXIS = {"theta": "0", "phi": "0", "x_5": "0"}
KK_HOPF = {"theta": "0", "phi": "0", "psi": "0"}
# The Kaluza-Klein black holes, in units of r_s: the electric hole at q = 2 r_s, a boost of
# cosh^2(alpha) = 2, the magnetic hole at p = 2 r_s, and the hole of equal charges at p = 3 r_s, whose
# horizons are 2 r_s and r_s. The electric and the equal charges have the circle divided out, and the
# magnetic hole's plane lies on the half axis theta = 0, where its potential vanishes.
KKBH_ELECTRIC = {"r_s": 1, "q": 2}
KKBH_MAGNETIC = {"r_s": 1, "p": 2}
KKBH_EQUAL = {"r_s": 1, "p": 3}
KKBH_CONE = "future cone of no momentum along the circle"
# X = r, Y = v - sqrt(q/r_s) r at q = 2 r_s, a time function of the ingoing chart at every r.
KKBH_FINKELSTEIN = ((0, 1), (1, -math.sqrt(2)))
# Myers and Perry's black hole: the plane transverse to the rotation, theta = 0, and the plane of
# rotation, theta = pi/2, in five dimensions and in six, and the hole of five dimensions with equal
# spins. Lengths are in sqrt(mu) in five dimensions and in the cube root of mu in six.
MYERS_ONE = {"mu": 1, "a": "3/5"}
MYERS_EQUAL = {"mu": 1, "a": "2/5"}
MYERS_SIX = {"mu": 1, "a": "3/2"}
MYERS_CONE = "future cone of no angular momentum"
# Randall and Sundrum's planes of t and the fifth coordinate hold the three coordinates along the wall fixed.
RS_WALL = {"x_1": "0", "x_2": "0", "x_3": "0"}
POLAR = "polar"                         # (phi, r) drawn from above: X = r cos phi, Y = r sin phi
PRINCIPAL_CONE = "future cone of the principal plane"
BTZ_CONE = "future cone of no angular momentum"
NHEK_CONE = "future cone of no angular momentum"


# Levi-Civita's cylinder as every one of its diagrams draws it: sigma = 1/4 and C = 1 in Weyl's
# coordinates, which is the point (2/3, 2/3, -1/3) of Kasner's circle, with ell = (3/4)^(2/3).
LC_WEYL = {"sigma": "1/4", "C": 1}
LC_KASNER = {"p_0": "2/3", "p_2": "2/3", "p_3": "-1/3", "ell": "(3/4)**(2/3)"}

# Bonnor's uniform beam of light as every one of its diagrams draws it: radius R, the unit, and
# pi G epsilon R^2/c^4 = 1/32, his m, so that A = rho^2/8 inside the beam and (1 + 2 ln rho)/8
# outside it, A = 1/8 at its edge. The free profile of the Cartesian chart is that beam, and the
# null Cartesian chart's is two of them side by side, their axes at x = +-2R, his section 6.
LB = {"G": 1, "epsilon": "1/(32*pi)", "R": 1}


def _lb_beam(x="x"):
    s = f"(({x})**2 + y**2)"
    return f"Piecewise(({s}/8, {s} <= 1), ((1 + log({s}))/8, True))"


LB_ONE = _lb_beam()
LB_TWO = _lb_beam("x - 2") + " + " + _lb_beam("x + 2")
LB_ONE_INPUT = ("A uniform beam of radius $R$ with $\\pi G\\epsilon R^2/c^4 = 1/32$, $\\epsilon$ its energy density: "
                "$A = \\rho^2/8R^2$ inside it and $A = (1 + 2\\ln(\\rho/R))/8$ outside, with $\\rho^2 = x^2 + y^2$.")
LB_TWO_INPUT = ("Two uniform beams of radius $R$ shining the same way, their axes at $x = \\pm 2R$, $y = 0$, each "
                "with $\\pi G\\epsilon R^2/c^4 = 1/32$: $A$ is the sum of the two profiles, each $\\rho^2/8R^2$ inside "
                "its beam and $(1 + 2\\ln(\\rho/R))/8$ outside, $\\rho$ the distance from that beam's axis.")
LB_FAMILIES = ("against the beam", "with the beam")
# The box holds the six wave fronts of the embedding diagram, u = 0 to 10 R, with sqrt(2) u = ct - z.
LB_BOX = (-8, 8, -2, 14)
LB_NULL_TO_TZ = ((-1 / math.sqrt(2), 1 / math.sqrt(2)), (1 / math.sqrt(2), 1 / math.sqrt(2)))  # z and ct of Bonnor's u and v

# Tolman's solution VII for a star of radius R = 2 r_s, in units of r_s: Lattimer and Prakash's
# parameters, and Tolman's constants for the same star.
TOLMAN_VII_STAR = {"r_s": 1, "R": 2}
TOLMAN_VII_CONSTANTS = {"R": "4/sqrt(5)", "A": "4/3**(1/4)", "B": "sqrt(7/12)",
                        "C": "sqrt(3)/2*(1/6 + sqrt(6)/3)*exp(2*atan(1/sqrt(6)) - pi)"}


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
    orient: str = "tau"             # tau, ingoing, outgoing, split or vector; see the header
    split: float = None             # the radius where orient="split" turns from ingoing to outgoing
    tau: str = "t"                  # the time function, in the Reader's plain names
    families: tuple = RADIAL        # the names of P and M on the page
    cones: tuple = (7, 7)
    areal: bool = False             # spherical: R^2 = g_theta_theta
    areal_contours: tuple = ()
    functions: dict = field(default_factory=dict)   # a declared function -> its expression
    delta: str = None               # the Dirac delta drawn as this pulse, an expression in s; see smoothed
    step: float = 0.0025            # the tracing step, in the unit square; a view whose metric changes
                                    # across a narrow pulse takes a finer one
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
    marked: tuple = ()              # (kind, point, family, legend[, "past" or "future"]): a ray of one family,
                                    # or of both, through a point, drawn and named; see Plot.marked
    points: tuple = ()              # (kind, (x^0, r), legend): an event drawn and named
    lines: tuple = ()               # (kind, "x0" or "r", value, legend): a line of constant
                                    # coordinate drawn and named, such as where a declared function jumps;
                                    # legend None for a marker the page names itself, as grr
    surface: str = None             # r at which a star's surface is released from rest; see Surface
    singular_runs: bool = False     # mark a singular stretch of an edge, not only a whole edge
    singular_where_claimed: bool = False  # judge a singular edge only inside the published domains
    singular_near: float = 1e-5     # how near an edge singular_runs takes the Kretschmann scalar, in the unit
                                    # square: nearer for a scalar whose coefficient falls to zero along the
                                    # edge, as Kiselev's 2(tau + rho)^2/rho^2 does toward tau = 0 on rho = 0
    star: dict = None               # a declared polytrope, {"K": ..., "rho_c": ...}; see StarSolver
    fronts: dict = None             # declared Robinson-Trautman initial data, {"epsilon": ...}; see FrontSolver
    any_factor: str = None          # a declared conformal factor the drawing holds for every value of
    crunch: bool = False            # mark where the metric stops being finite as a singular curve,
                                    # checked on the Kretschmann scalar, and hatch what lies beyond it
    solves: tuple = ()              # published Einstein components the declared functions must zero
    singular_zero: str = None       # an expression in the chart's plain names whose zero set the row
                                    # declares a curvature singularity, where the Kretschmann scalar
                                    # diverges too slowly for the test at 1e-5 of the drawing; drawn as a
                                    # singular curve and checked in 60 digits, see Plot.weak_singularity
    ends_on_singular: bool = False  # a ray ends on the row's singular_zero curve, where the expression
                                    # turns negative: Kastor and Traschen's metric is built on U^2, so the
                                    # null condition runs on through U = 0, beyond which there is no spacetime
    no_throat: bool = False         # the areal radius is stationary along the drawn radius on a curve that
                                    # is no throat, as on tau = 0 of the lukewarm hole's cosmological chart,
                                    # where every rho has the areal radius r_s/2; the curve is left unmarked
    where: str = None               # an expression in the chart's plain names that is positive on the
                                    # spacetime, for a published domain that is an inequality between the
                                    # plane's two coordinates, as UV < 1 in the Kruskal chart of Witten's
                                    # black hole: the view is hatched where it is not positive, and no ray
                                    # or cone is drawn there
    where_is_infinity: bool = False  # the edge of `where` is infinity, where a conformal factor diverges and
                                    # g^rr vanishes with it, as on tau + rho = 0 of Kiselev's conformally
                                    # flat chart: no zero of g^rr is marked on that edge
    null_radius: bool = False       # the drawn radius is a null coordinate, so the areal radius is stationary
                                    # along it exactly on the marginally trapped spheres: that curve is
                                    # marked as the apparent horizon, and no throat is
    quotient: str = None            # a coordinate the metric does not depend on, divided out: the
                                    # plane's metric is g_ab - g_ak g_bk/g_kk; see "Rays of no angular
                                    # momentum"


# Coleman and De Luccia's two bubbles, thin walled, in units of the curvature radius l: rho(xi) of
# the wall chart on both sides of the wall, continuous, with rho' dropping across it.
CDL_INTO_FLAT = "Piecewise((xi, xi < Rational(4, 5)), (sin(xi - Rational(4, 5) + asin(Rational(4, 5))), True))"
CDL_OUT_OF_FLAT = "Piecewise((sinh(xi), xi < log(3)), (xi - log(3) + Rational(4, 3), True))"
# The wall in the static chart of the anti-de Sitter space inside it, which it leaves at ct = pi l/2.
CDL_ADS_WALL = "Piecewise((sqrt(Rational(16, 9) + Rational(25, 9)*tan(t)**2), t < pi/2), (10**9, True))"
CDL_INPUT = {
    "into_flat": "A thin wall between a false vacuum of positive energy, de Sitter space of radius $\\ell$, "
                 "and a true vacuum of zero energy, born with the radius $\\bar\\rho = 4\\ell/5$, where "
                 "$\\partial_\\xi\\rho$ drops from $1$ to $3/5$.",
    "out_of_flat": "A thin wall between a false vacuum of zero energy and a true vacuum of negative energy, "
                   "anti-de Sitter space of radius $\\ell$, born with the radius $\\bar\\rho = 4\\ell/3$, where "
                   "$\\partial_\\xi\\rho$ drops from $5/3$ to $1$.",
}


def _alcubierre_profile():
    return ("(tanh(4*(sqrt((x - 2*t)**2 + y**2 + z**2) + 1))"
            " - tanh(4*(sqrt((x - 2*t)**2 + y**2 + z**2) - 1)))/(2*tanh(4))")


# Van Den Broeck's pocket as Krasnikov shapes it, a surface every circle of which grows no faster
# than the distance out to it, so that it stands in flat space, in units of R, the radius where
# the wall of the bubble begins. In the proper distance l from the middle of the neck the areal
# radius is r = l + 4 on the flat floor, 7/4 - (l + 2)^2 round the rim, -l on the flat lid,
# l^2 + 1/4 through the neck, Krasnikov's quadratic, and l outside. van_den_broeck.md derives the
# factor B of the comoving radius rho from it, through d ln(rho) = dl/r and B = r/rho.
_VDB_R = ("Piecewise((l + 4, l < -Rational(5, 2)), (Rational(7, 4) - (l + 2)**2, l < -Rational(3, 2)),"
          " (-l, l < -Rational(1, 2)), (l**2 + Rational(1, 4), l < Rational(1, 2)), (l, True))")
_VDB_RHO = {"neck": "exp(-pi)/2", "lid": "exp(-pi)/6",
            "floor": "exp(-pi - 4*atanh(1/sqrt(7))/sqrt(7))/6"}


def _vdb_factor(rho):
    """B of the declared pocket as a function of the comoving radius, written in `rho`."""
    return ("Piecewise((Rational(3, 2)/({floor}), {rho} < {floor}),"
            " (Rational(7, 4)/(({rho})*cosh(sqrt(7)/2*log(({rho})/({lid})) + atanh(1/sqrt(7)))**2), {rho} < {lid}),"
            " (exp(-pi)/(4*({rho})**2), {rho} < {neck}),"
            " (1/(2*({rho})*(1 - sin(log(2*({rho}))))), {rho} < Rational(1, 2)), (1, True))"
            ).format(rho=rho, **_VDB_RHO)


def _vdb_shape(rho):
    """The wall of the declared bubble, f falling from 1 at R to 0 at 3R/2 with two continuous
    derivatives, as a function of the comoving radius, written in `rho`."""
    s = f"(2*(({rho}) - 1))"
    return (f"Piecewise((1, {rho} < 1), (1 - 10*{s}**3 + 15*{s}**4 - 6*{s}**5, {rho} < Rational(3, 2)), (0, True))")


_VDB_RS = "sqrt((x - 2*t)**2 + y**2 + z**2)"
VDB_POCKET = ("Krasnikov's shape for the pocket, the areal radius $r(l)$ of the sphere at proper distance $l$ from "
              "the middle of the neck: $l + 4R$ on a flat floor, $7R/4 - (l + 2R)^2/R$ round the rim, $-l$ on a "
              "flat lid, $l^2/R + R/4$ through the neck, and $l$ outside it, which makes $1 + \\alpha = 380$")
VDB_WALL = "$f$ falls from $1$ at $r_s = R$ to $0$ at $3R/2$ as $1 - 10s^3 + 15s^4 - 6s^5$, $s = 2(r_s - R)/R$"


def _vdb_distance(rho):
    """The proper distance l from the middle of the neck at the comoving radius rho, the closed
    form the rays of the comoving chart are checked against."""
    rho = np.asarray(rho, dtype=float)
    neck, lid = math.exp(-math.pi) / 2, math.exp(-math.pi) / 6
    a = math.atanh(1 / math.sqrt(7))
    floor = lid * math.exp(-4 * a / math.sqrt(7))
    with np.errstate(all="ignore"):
        return np.select(
            [rho < floor, rho < lid, rho < neck, rho < 0.5],
            [-4 + 1.5 * rho / floor,
             -2 + math.sqrt(7) / 2 * np.tanh(math.sqrt(7) / 2 * np.log(rho / lid) + a),
             -math.exp(-math.pi) / (4 * rho),
             np.tan(np.log(2 * rho) / 2 + math.pi / 4) / 2], rho)


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

# The dipole of the Szekeres cloud: S'/S = 2r(1 - r^2) inside r_b = 1, and S constant beyond.
_SZ_S = "Piecewise((exp(r**2 - r**4/2), r < 1), (exp(Rational(1, 2)), True))"
SZ_INPUT = ("Marginally bound dust, $f = 0$, with the areal radius and the mass of the Tolman-Bondi cloud, "
            "$R(0, r) = r$ and $2GM(r)/c^2 = r^3\\left(5 - 3r^2/r_b^2\\right)/4r_b^2$ inside $r_b$, each shell "
            "falling as $R^{3/2} = r^{3/2} - \\tfrac{3}{2}\\sqrt{2GM(r)/c^2}\\,ct$, and with "
            "$S = \\exp\\left(r^2/r_b^2 - r^4/2r_b^4\\right)$ inside $r_b$ and constant beyond, checked to solve "
            "this spacetime's own $G^r{}_r = 0$.")


def _er_pulse(t, rho):
    """The pulse of Weber, Wheeler, and Bonnor at C = a = 1, psi and gamma as strings in the
    plain names of a chart whose ct and rho are the expressions t and rho: D_+ and D_-^2 are
    (1 + rho^2 - t^2)^2 +- 4t^2, and D_-^2 is written as the polynomial it is, never as a root."""
    q = f"(1 + ({rho})**2 - ({t})**2)"
    dp = f"sqrt({q}**2 + 4*({t})**2)"
    dm2 = f"({q}**2 - 4*({t})**2)"
    return {"psi": f"sqrt(2)*sqrt({dp} + {q})/{dp}",
            "gamma": f"(1 - 2*({rho})**2*{dm2}/{dp}**4 + (({rho})**2 - 1 - ({t})**2)/{dp})/2"}


ER_INPUT = ("The pulse of Weber, Wheeler, and Bonnor, of width $a$: $\\psi = \\sqrt{2}\\,C\\sqrt{D_+ + a^2 + \\rho^2 - c^2t^2}/D_+$ "
            "and $\\gamma = \\tfrac{C^2}{2a^2}\\left(1 - 2a^2\\rho^2D_-^2/D_+^4 + (\\rho^2 - a^2 - c^2t^2)/D_+\\right)$, "
            "with $D_\\pm^2 = (a^2 + \\rho^2 - c^2t^2)^2 \\pm 4a^2c^2t^2$, at $C = a$, checked to solve this "
            "spacetime's own field equations.")
ER_SOLVES = (("t", "t"), ("t", "\\rho"), ("\\phi", "\\phi"), ("z", "z"))

def _gowdy_wave(t):
    """A polarised Gowdy wave on the torus, as strings in the plain names of a chart whose areal
    time is the expression t: Q = 0, P = A Y_0(t) cos(theta) and
    lambda = A^2 (t^2 (Y_0^2 + Y_1^2)/2 - t Y_0 Y_1 cos^2(theta)), with Y_0 and Y_1 the Bessel
    functions of the second kind and A = -pi/4, so that P -> (cos(theta)/2) tau as t -> 0."""
    y0, y1 = f"bessely(0, {t})", f"bessely(1, {t})"
    return {"P": f"-pi/4*{y0}*cos(theta)", "Q": "0",
            "lambda": f"(pi/4)**2*(({t})**2*({y0}**2 + {y1}**2)/2 - ({t})*{y0}*{y1}*cos(theta)**2)"}


GOWDY_INPUT = ("A polarised wave once round the torus, $Q = 0$, $P = -\\tfrac{\\pi}{4}Y_0(t)\\cos\\theta$, and "
               "$\\lambda = \\tfrac{\\pi^2}{16}\\left(\\tfrac{1}{2}t^2\\left(Y_0^2 + Y_1^2\\right) - tY_0Y_1\\cos^2\\theta\\right)$, "
               "with $Y_0$ and $Y_1$ the Bessel functions of the second kind at $t$, checked to solve this "
               "spacetime's own field equations; as $t \\to 0$, $P \\to \\tfrac{1}{2}\\cos\\theta\\,\\tau$ with "
               "$\\tau = -\\ln t$.")


# The inside of Schwarzschild's horizon in the sphere chart. print_charts.py checks in sympy that it
# is Schwarzschild's metric and a vacuum, exactly, so its row names no `solves`: beside the
# singularity the terms of G^mu_nu are 10^10 and their sum is lost in the rounding.
GOWDY_HOLE = {"P": "log((1 - cos(t))**2*sin(theta)/sin(t))", "Q": "0", "a": "log(1 - cos(t))"}
GOWDY_HOLE_INPUT = ("The inside of Schwarzschild's horizon, $Q = 0$, $e^{2a} = (1 - \\cos t)^2$, and "
                    "$e^{P} = (1 - \\cos t)^2\\sin\\theta/\\sin t$, checked to solve this spacetime's own field "
                    "equations; $r = L(1 - \\cos t)$ is Schwarzschild's radius with $r_s = 2L$, and $L\\delta$ "
                    "his time.")


def _gowdy_solves(time):
    return ((time, time), (time, "\\theta"), ("\\theta", "\\theta"), ("\\sigma", "\\sigma"), ("\\delta", "\\delta"))


FRW_DUST = {"funcs": ["a"], "eqs": [["r", "r"]], "rates": [1.0], "params": {"k": 0}}
# The Oppenheimer-Snyder dust released from rest at a = a_m, which is the unit, at tau = 0.
OS_DUST = {"funcs": ["a"], "eqs": [["\\chi", "\\chi"]], "rates": [0.0], "start": [1.0], "origin": "reference"}

# The polytrope the conformal diagram declares, and what it makes of the star.
POLYTROPE = {"K": 100, "rho_c": "1.28e-3"}
POLYTROPE_INPUT = ("A polytrope, $p = K\\rho_0^2$ with rest mass density $\\rho_0$ and energy density "
                   "$\\rho c^2$, where $\\rho = \\rho_0 + p/c^2$, at $K = 100$ and a central $\\rho_0 = 1.28\\times10^{-3}$ in "
                   "units where $G = c = M_\\odot = 1$, solved from this spacetime's own $G^t{}_t$ and "
                   "$G^r{}_r$: a star of $M = 1.40\\,M_\\odot$ and $R = 14.2$ km.")
# Robinson and Trautman's fronts at u = 0: Macedo and Saa's prolate data, drawn at epsilon = 4/5.
# Kinnersley's photon rocket: a burn from u = 0 to cu = 10 m_0, in units of the mass m_0 it
# starts with, whose acceleration rises and falls as 2 sin^2(pi cu/10 m_0)/(25 m_0), losing mass
# at the least rate that keeps the density of its radiation positive in every direction,
# d_u m = -3 alpha m. Its rapidity w is the integral of alpha, 2/5 at the end, and the
# Robinson-Trautman chart's four-velocity is (cosh w, sinh w, 0, 0). alpha m stays below 1/16,
# 0.050 at most, so behind the rocket the two zeros of g^rr never meet.
_ROCKET_U = "Min(Max(u, 0), 10)"
_ROCKET_W = f"({_ROCKET_U} - 5*sin(pi*{_ROCKET_U}/5)/pi)/25"
ROCKET_BURN = {"alpha": f"2*sin(pi*{_ROCKET_U}/10)**2/25", "m": f"exp(-3*{_ROCKET_W})"}
ROCKET_FLIGHT = {"U_1": f"sinh({_ROCKET_W})", "U_2": "0*u", "U_3": "0*u", "m": f"exp(-3*{_ROCKET_W})"}
ROCKET_INPUT = ("A burn from $u = 0$ to $cu = 10\\,m_0$ along the first axis, with the acceleration "
                "$\\alpha = 2\\sin^2(\\pi cu/10m_0)/(25\\,m_0)$ and the rapidity $w = \\int\\alpha\\,c\\,du$, which ends at "
                "$2/5$, a speed of $0.38\\,c$. The mass is $m = m_0e^{-3w}$ from the mass $m_0$ at the start, the least loss that keeps the density "
                "of the radiation positive in every direction, and ends at $0.30\\,m_0$.")
RT_FRONTS = {"epsilon": "4/5"}
RT_INPUT = ("The first front $f(0, \\theta)^2 = f_0^2\\left(1 - \\epsilon^2\\cos^2\\theta\\right)$ at "
            "$\\epsilon = 4/5$, with $f_0^2 = \\ln\\left((1 + \\epsilon)/(1 - \\epsilon)\\right)/(2\\epsilon)$ so that "
            "every front has the area $4\\pi r^2$, and after it the solution of the Robinson-Trautman equation, "
            "with the vacuum $2H = K - 2r\\,\\partial_u\\ln f - 2m/r$.")

POLYTROPE_STATED = {"M": (1.40, 2), "R_km": (14.2, 1)}
KM = 1.4766250614                       # GM_sun/c^2 in km

BIANCHI_DUST = {"funcs": ["a_1", "a_2", "a_3"], "eqs": [["x", "x"], ["y", "y"], ["z", "z"]],
                "rates": [-0.5, 1.5, 2.0]}
# The Kantowski-Sachs dust at its moment of greatest b, where a = 1 and b = b_0, the unit, both at rest.
KS_DUST = {"funcs": ["a", "b"], "eqs": [["r", "r"], ["\\theta", "\\theta"]], "rates": [0.0, 0.0], "start": [1.0, 1.0]}
KS_INPUT = ("Dust: $a(t)$ and $b(t)$ solved from this spacetime's own $G^r{}_r = G^\\theta{}_\\theta = 0$, starting "
            "from $a = 1$ and $b = b_0$ at rest at the dashed line.")

# Godel's radius r_c = ln(1 + sqrt 2), where sinh r = 1 and the circles of constant t, r and z
# turn from spacelike to timelike; the views read it off the published g_phiphi as well.
GODEL_RC = math.asinh(1.0)

# Kottler's black hole at Lambda r_s^2 = 1/5, so that r_h = 1.085 r_s and r_c = 3.215 r_s, and the
# radius (3r_s/2Lambda)^(1/3) of its static observer in free fall, where f is greatest.
SDS = {"r_s": 1, "Lambda": "1/5"}
# Visser's thin shell wormhole with its throat at a = 5r_s/4, inside the photon sphere 3r_s/2.
TSW = {"r_s": 1, "a": "5/4"}
# The gravastar every diagram draws: the shell at R = 5 r_s/4 round a ball of de Sitter space of
# radius L = 2 r_s, where C = 64/195 and the shell has a positive surface density.
GRAVASTAR = {"r_s": 1, "R": "5/4", "L": 2}
GRAVASTAR_C = 64 / 195
GRAVASTAR_X = 2 / math.sqrt(GRAVASTAR_C) * math.atanh(5 / 8)
SDS_STATIC = 7.5 ** (1 / 3)
# Kiselev's black hole in his example w = -2/3 at r_q = 8 r_s, where f = 1 - r_s/r - r/r_q vanishes at
# r = (4 -+ 2 sqrt 2) r_s, the black hole horizon and the outer one, and is greatest at sqrt(r_s r_q),
# the radius of the static observer in free fall. KISELEV_FREE is the same matter with no black hole.
KISELEV = {"r_s": 1, "r_q": 8, "w": "-2/3"}
KISELEV_LINEAR = {"r_s": 1, "r_q": 8}
KISELEV_STATIC = 8 ** 0.5
KISELEV_ROOTS = (4 - 2 * 2 ** 0.5, 4 + 2 * 2 ** 0.5)
KISELEV_FREE = {"r_q": 1}
# Mass inflation as every diagram draws it, Ori's model at the numbers of ori_shell.py, in units of
# the final mass m_0: the charge 0.96 m_0 of Reissner-Nordstrom's own diagrams, Price's tail of a
# quadrupole switched on at v_0 = 10 m_0 with m_0/50 still to fall in, and a shell of mass m_0/50
# that crosses v_0 at r = m_0. Behind the shell the advanced time is zero on the Cauchy horizon.
MI = {"r_q": "24/25"}
MI_MASS = "1 - (10/Max(v, 10))**11/50"
MI_INPUT = ("$m(v) = m_0 - (m_0/50)(v_0/v)^{11}$ from $v_0 = 10\\,m_0$ on and $0.98\\,m_0$ before, an influx "
            "that falls off as $v^{-12}$, with $m_0$ the final mass and the unit of every length.")
MI_BEHIND_INPUT = ("Ori's shell of mass $m_0/50$, crossing the ingoing ray $v_0$ at $r = m_0$, in front of the tail "
                   "$m_0 - (m_0/50)(v_0/v)^{11}$: the mass function behind it, $m(v)$, is integrated along the "
                   "shell, with $v = 0$ the Cauchy horizon and $m_0$ the final mass outside.")
MI_V2_AT_V0 = float(ori_shell.shell().v2_at(ori_shell.V0))
# The lukewarm charged black hole in de Sitter space, r_q = r_s/2 and Lambda r_s^2 = 27/64, which is
# H r_s/c = 3/8: f = (1 - 1/(2r))^2 - 9r^2/64 vanishes at r_c = 2, r_+ = 2/3, r_- = (2 sqrt 7 - 4)/3 and
# -(2 sqrt 7 + 4)/3, and is greatest between r_+ and r_c at the root 1.298 of 9r^4 - 32r + 16.
RNDS = {"r_s": 1, "r_q": "1/2", "Lambda": "27/64"}
RNDS_ROOTS = (2.0, 2 / 3, (2 * math.sqrt(7) - 4) / 3, -(2 * math.sqrt(7) + 4) / 3)
RNDS_STATIC = 1.2979848366419
# A time function of its Eddington-Finkelstein charts down to r = 0: v - h(r) with h' = r^2/(r^2 + r_q^2),
# which lies between 0 and 2/f wherever f > 0, since f < 1 + r_q^2/r^2; v - r alone is spacelike where f > 2.
RNDS_TIME_IN = "v - r + atan(2*r)/2"
RNDS_TIME_OUT = "u + r - atan(2*r)/2"
# The charged black hole in anti-de Sitter space at L = 1, r_s = 27/8 and r_q^2 = 11/8, where
# r^2 f = (r - 1)(r - 1/2)(r^2 + 3r/2 + 11/4): the horizons are r_+ = L and r_- = L/2, the first the
# horizon Schwarzschild-anti-de Sitter is drawn with.
RNADS = {"r_s": "27/8", "r_q": "sqrt(22)/4", "L": 1}
# A time function of its Eddington-Finkelstein charts at every r > 0: v - h(r) with
# h' = r^2/((r^2 + 1)(r^2 + 2)), which lies between 0 and 2/f wherever f > 0, since
# f < 1 + r_q^2/r^2 + r^2 < (r^2 + 1)(r^2 + 2)/r^2; v - r alone is spacelike where f > 2, near the
# singularity and far outside.
RNADS_TIME_IN = "v + atan(r) - sqrt(2)*atan(r/sqrt(2))"
RNADS_TIME_OUT = "u - atan(r) + sqrt(2)*atan(r/sqrt(2))"
# Kerr-de Sitter at Kerr's spin and Kottler's cosmological constant, a = 0.45 r_s and
# Lambda = 0.2/r_s^2, so that r_- = 0.279, r_+ = 0.785 and r_c = 3.232 r_s, and Kerr-anti-de Sitter
# at r_s = 2 l and a = l/2 with l = sqrt(-3/Lambda) the unit, so that r_- = 0.137 and r_+ = 0.859 l.
# KDS_SPLIT is a radius between r_+ and r_c, where both orientations agree with t.
KDS = {"r_s": 1, "a": "9/20", "Lambda": "1/5"}
KADS = {"r_s": 2, "a": "1/2", "Lambda": -3}
KDS_SPLIT = 2.0
KDS_AXIS = {"theta": "0", "phi": "0"}
KDS_KERR_AXIS = {"theta": "0", "tildephi": "0"}

# The Aichelburg-Sexl shock on the plane of u and v at three distances from the source, in units
# of 8GE/c^4 with rho_0 = 8GE/c^4: each ray moving left jumps along the shock by -ln(rho/rho_0),
# ln 2, ln 8 and ln 32, the same step of ln 4 between neighbours. The delta is drawn as a pulse of
# width 1/20, which moves every such ray by the same amount, traced at a fifth of the usual step so
# that a ray turning into the pulse and out of it keeps its closed form to 1e-5. The box stands half a unit higher in ct
# than it is wide in z, so that no cone of the lattice falls within seven widths of the shock.
AS_PULSE = "20*exp(-400*s**2)/sqrt(pi)"
AS_INPUT = ("$\\delta(u)$ drawn as the pulse $e^{-u^2/w^2}/(w\\sqrt{\\pi})$ with $w = 0.05$, the field of a pulse "
            "of light of that length in $u$ along the axis, carrying the same energy $E$.")
AS_RHO = {"half": "1/2", "eighth": "1/8", "thirtysecond": "1/32"}

# The global monopole at Delta = 0.19, so that sqrt(1 - Delta) = 0.9 and the equatorial cone of its
# embedding diagram lacks the 36 degrees the cosmic string's does, and Letelier's black hole at
# r_s = 1 inside it, with its horizon at r_s/(1 - Delta) = 100/81 r_s.
GM = {"Delta": "19/100", "r_s": 1}
# McVittie's mass in a universe of dust and a cosmological constant, the expansion Lake and
# Abdelqader chose: H = H_0 coth(3 H_0 t/2), a = sinh^(2/3)(3 H_0 t/2), with H_0 = c/(sqrt(15) r_s),
# which is Lambda r_s^2 = 3 H_0^2 r_s^2/c^2 = 1/5, the value Schwarzschild-de Sitter is drawn at.
MCV_H = "coth(3*t/(2*sqrt(15)))/sqrt(15)"
MCV_A = "sinh(3*t/(2*sqrt(15)))**Rational(2, 3)"
MCV_INPUT = ("A universe of dust and a cosmological constant, $a = \\sinh^{2/3}(3H_0t/2)$ and "
             "$H = H_0\\coth(3H_0t/2)$, the expansion Kayll Lake and Majd Abdelqader chose, with "
             "$H_0 = c/(\\sqrt{15}\\,r_s)$, which is $\\Lambda r_s^2 = 1/5$ for the cosmological constant $\\Lambda$.")

# Gott's two strings at half deficit angle alpha = pi/3, 4 G mu/c^2 = 1/3, moving at v = 4c/5
# with d = l/2, where gamma sin(alpha) = 5/(2 sqrt 3) > 1: the boost round both strings has
# cosh(a/4) = gamma sin(alpha), and the shift is b = 4 gamma v d sin(alpha)/(c sinh(a/4)).
GOTT = {"a": "4*acosh(5/(2*sqrt(3)))", "b": "8/sqrt(13)"}
GOTT_A = 4 * math.acosh(5 / (2 * math.sqrt(3)))
# Ori's vacuum core as every one of its diagrams draws it: his example f = a(x^2 - y^2)/2 at
# a = 1/16, with e = 1/8 for the surfaces of constant t, so that e > a and e > (2e + a)^2 = 25/256,
# and z running once round 2 pi.
ORI = {"a": "1/16", "e": "1/8", "L": "2*pi"}
ORI_F = "(x**2 - y**2)/32"
ORI_INPUT = "Ori's example, $f = a(x^2 - y^2)/2$ at $a = 1/16$, a vacuum, with $L = 2\\pi$."
# The spinning string: b = 0.9, the cosmic string's deficit, and a = 0.9 in units of r_c = a/b.
SPINNING = {"a": "9/10", "b": "9/10"}
# The travelling wave on a string: b = 1/2, a string whose cone lacks half a turn, heavy enough for
# the wave's field to show, with ell the unit of length and the string displaced along one
# transverse direction by the pulse A = (ell/2) exp(-4u^2/ell^2), B = 0. In the null conical chart
# that is the profile F = -2 ell (b r/ell)^(1/b) A'' cos(phi), harmonic on the cone.
STRING_WAVE = {"b": "1/2", "ell": 1}
# A point particle in three dimensions of mass c^2/16G, whose deficit angle is a right angle, and
# two of them at rest at x = +-d, with ell = d the unit of length.
PARTICLE = {"alpha": "3/4"}
TWO_PARTICLES = {"alpha_1": "3/4", "alpha_2": "3/4", "d": 1, "ell": 1}
STRING_WAVE_PULSE = {"A": "exp(-4*u**2)/2", "B": "0"}
STRING_WAVE_PROFILE = {"F": "-2*(r/2)**2*(32*u**2 - 4)*exp(-4*u**2)*cos(phi)"}
STRING_WAVE_INPUT = ("The pulse $A = (\\ell/2)\\,e^{-4u^2/\\ell^2}$, $B = 0$ on a string with $b = 1/2$, whose cone "
                     "lacks half a turn.")
STRING_WAVE_PROFILE_INPUT = ("The string's own travelling wave, $F = -2\\ell(br/\\ell)^{1/b}(A''\\cos\\phi + "
                             "B''\\sin\\phi)$, for the pulse $A = (\\ell/2)\\,e^{-4u^2/\\ell^2}$, $B = 0$ on a "
                             "string with $b = 1/2$, whose cone lacks half a turn.")
STRING_WAVE_CREST = (("shell", "x0", "0", "the crest of the pulse, $u = 0$"),)

# Siklos's waves, in units of the anti-de Sitter radius L: the free charts are drawn for
# Kaigorodov's profile under a pulse, H = exp(-4u^2) x^3, which solves Siklos's equation since u
# never enters it, and which is h = LH/x = exp(-4u^2) p^2/q^2 on the disc; Kaigorodov's own
# charts are drawn with the Killing direction d/du divided out, d/dx in the Poincare chart.
SIKLOS = {"L": 1}
SIKLOS_PULSE = {"H": "exp(-4*u**2)*x**3"}
SIKLOS_DISC_PULSE = {"h": "exp(-4*u**2)*(1 - (xi**2 + eta**2)/4)**2/((1 + xi/2)**2 + eta**2/4)**2"}
SIKLOS_INPUT = ("Kaigorodov's profile switched on and off, $H = e^{-4u^2/L^2}x^3/L^3$, which solves "
                "$\\partial_x^2H + \\partial_y^2H - (2/x)\\,\\partial_xH = 0$.")
SIKLOS_DISC_INPUT = ("Kaigorodov's profile switched on and off, $h = e^{-4u^2/L^2}p^2/q^2$, which solves "
                     "$p^2(\\partial_\\xi^2h + \\partial_\\eta^2h) = 2h/L^2$.")
SIKLOS_CREST = (("shell", "x0", "0", "the crest of the pulse, $u = 0$"),)
SIKLOS_CONE = "future cone of no momentum along ${}$"
GOTT_STRINGS = {"mu": "1/12", "G": 1, "v": "4/5", "d": "1/2", "alpha": "pi/3", "gamma": "5/3"}

# Morris, Thorne and Yurtsever's round trip, in the throat radius r_0 and with c = 1: the right
# mouth's rapidity is eta = (3/2) sin^3(2 pi tau/P) over 0 <= tau <= P = 54 r_0, out for the first
# half and home for the second, so its acceleration g = d eta/d tau is greatest, 0.2015/r_0, at
# tau = (P/2 pi) arctan(sqrt 2) = 8.21 r_0 and least, the same with a minus sign, at P/2 - 8.21 r_0.
# The mouths start a distance D = 10 r_0 apart. wormhole_time_machine.md derives the rest.
WTM_ETA0, WTM_P, WTM_D = 1.5, 54.0, 10.0
WTM_G = "pi/6*sin(pi*t/27)**2*cos(pi*t/27)"
WTM_FUNCTIONS = {"g": WTM_G, "F": "Piecewise((0, l <= 0), (exp(-1/(4*l**2)), True))", "Phi": "0", "r": "sqrt(1 + l**2)"}
WTM_PEAK = WTM_P / (2 * math.pi) * math.atan(math.sqrt(2))


class WormholeTrip:
    """The world line of the right mouth on Morris, Thorne and Yurtsever's round trip, in the
    Lorentz coordinates (T, Z) of the flat space outside, with c = 1 and lengths in r_0.

    The left mouth rests at Z = 0 and its proper time is T. The right mouth starts at Z = D
    with the same proper time tau = T = 0, and moves with rapidity eta(tau) = eta_0
    sin^3(2 pi tau/P) for 0 <= tau <= P and none after, so that dT/dtau = cosh(eta) and
    dZ/dtau = sinh(eta), integrated here by Simpson's rule. eta is odd about P/2, so the mouth
    comes home, Z(P) = D, and it has aged less than the left one by shift = T(P) - P. The two
    mouths are identified at equal tau. A light ray leaving the left mouth at tau along +Z
    reaches the right one at T = tau + Z, and tau_c is the first tau at which that is the
    right mouth's own T(tau): the closed null geodesic, the only one, since after it the
    light arrives early at every tau. Every claim is checked as it is made."""

    def __init__(self, eta0=WTM_ETA0, period=WTM_P, distance=WTM_D, n=200001):
        from scipy.integrate import cumulative_simpson
        from scipy.optimize import brentq
        self.eta0, self.period, self.distance = eta0, period, distance
        tau = np.linspace(0.0, period, n)
        eta = eta0 * np.sin(2 * np.pi * tau / period) ** 3
        self._tau = tau
        self._T = cumulative_simpson(np.cosh(eta), x=tau, initial=0.0)
        self._Z = distance + cumulative_simpson(np.sinh(eta), x=tau, initial=0.0)
        self.shift = float(self._T[-1] - period)
        if abs(self._Z[-1] - distance) > 1e-9:
            raise SystemExit("the wormhole's right mouth does not come home")
        if not self.shift > distance:
            raise SystemExit("the wormhole's time shift is less than the distance between its mouths")
        gap = lambda t: self.right(t)[0] - t - self.right(t)[1]
        grid = np.linspace(0.0, period, 5401)
        values = np.array([gap(t) for t in grid])
        first = int(np.flatnonzero(values >= 0)[0])
        self.tau_c = float(brentq(gap, grid[first - 1], grid[first], xtol=1e-13))
        if not (values[:first] < 0).all():
            raise SystemExit("a closed causal curve threads the wormhole before the closed null geodesic")
        if not (values[first:] >= 0).all():
            raise SystemExit("the closed null geodesic is not the only one: light from the left mouth "
                             "arrives late again after it")

    def rapidity(self, tau):
        tau = np.asarray(tau, dtype=float)
        inside = (tau >= 0) & (tau <= self.period)
        return np.where(inside, self.eta0 * np.sin(2 * np.pi * tau / self.period) ** 3, 0.0)

    def right(self, tau):
        """(T, Z) of the right mouth at its proper time tau: at rest before the trip and after."""
        tau = np.asarray(tau, dtype=float)
        inside = np.clip(tau, 0.0, self.period)
        T = np.interp(inside, self._tau, self._T) + (tau - inside)
        return T, np.interp(inside, self._tau, self._Z)

    def left(self, tau):
        tau = np.asarray(tau, dtype=float)
        return tau, np.zeros_like(tau)
WTM_INPUT = ("$\\Phi = 0$ and $r = \\sqrt{r_0^2 + l^2}$, the Ellis-Bronnikov wormhole, with the form factor "
             "$F = e^{-r_0^2/4l^2}$ for $l > 0$ and the acceleration $g = d\\eta/dt$ of a round trip with rapidity "
             "$\\eta = \\tfrac{3}{2}\\sin^3(2\\pi t/P)$ and $P = 54\\,r_0/c$, greatest at $0.2\\,c^2/r_0$.")
GM_CONE = {"Delta": "19/100"}
GM_RH = 100 / 81
# The black hole on a cosmic string at the deficit the cosmic string is drawn at, 4G mu/c^2 = 0.1.
SBH = {"r_s": 1, "b": "9/10"}
DILATON = {"r_s": 1, "r_d": "1/2"}
# Schwarzschild-anti-de Sitter at r_s = 2L, where r^3 + L^2 r - L^2 r_s = (r - L)(r^2 + L r + 2L^2)
# and the horizon is r_h = L, the black hole of Hawking and Page's temperature T_1.
SADS = {"r_s": 2, "L": 1}
# The topological black holes as every one of their diagrams draws them, in units of L: the flat
# hole at mu = L, where r^3 - L^2 mu = (r - L)(r^2 + L r + L^2) and the horizon is r_h = L; the
# hyperbolic hole without mass, whose horizon is r_h = L too; and the hyperbolic hole of negative
# mass mu = -120L/343, where r^3 - L^2 r - L^2 mu = (r - 5L/7)(r - 3L/7)(r + 8L/7).
TBH_FLAT = {"mu": 1, "L": 1, "k": 0}
TBH_MASSLESS = {"mu": 0, "L": 1, "k": -1}
TBH_NEGATIVE = {"mu": "-120/343", "L": 1, "k": -1}
TBH_POINT = {"rho": "1/2", "phi": "0"}
TBH_THETA = {"theta": "1", "phi": "0"}
TBH_CASES = (("flat", "$k = 0$", TBH_FLAT, (0, 3, -1.5, 1.5)),
             ("massless", "$k = -1$, $\\mu = 0$", TBH_MASSLESS, (0, 3, -1.5, 1.5)),
             ("negative", "$k = -1$, $\\mu < 0$", TBH_NEGATIVE, (0, 1.5, -0.75, 0.75)))
# Bardeen's regular black hole at g = r_s/3, below the extremal 2 r_s/(3 sqrt 3) = 0.385 r_s: two
# horizons, r_- = 0.301 r_s and r_+ = 0.775 r_s, about a regular centre.
BARDEEN = {"r_s": 1, "g": "1/3"}

# Damour and Solodukhin's wormhole at lambda = 1/5, where the throat's clocks run five times slow.
DS = {"r_s": 1, "lambda": "1/5"}
# Einstein and Rosen's charged bridge with a mass, at r_q = sqrt(3) r_s/2, where r^2 - r_s r - r_q^2 =
# (r - 3r_s/2)(r + r_s/2) and the bridge is the sphere r_+ = 3r_s/2.
ERB_CHARGED = {"r_s": 1, "r_q": "sqrt(3)/2"}
# Fisher, Janis, Newman and Winicour's scalar field at gamma = 1/2, the value of Abdolrahimi and Shoom's
# figures, in units of b; the harmonic chart in units of k = b/2, where m = gamma k.
FJNW = {"b": 1, "gamma": "1/2"}
FJNW_HARMONIC = {"m": "1/2", "k": 1}
# Roberts's collapsing scalar field for its three outcomes: p = 9/10, where the field disperses, the
# threshold p = 1, and p = 2, where it makes a black hole. Nothing in it sets a scale, so lengths
# are in any unit ell.
ROBERTS_CASES = {"disperses": ("$p = 9/10$", "9/10"), "threshold": ("$p = 1$", "1"), "collapses": ("$p = 2$", "2")}

# Witten's black hole in two dimensions in units of 1/lambda, at m = 1, where the horizon is x = 0
# in the charts of x and w = 1 in the dilaton chart.
WITTEN = {"lambda": 1, "m": 1}

# Simpson and Visser's three geometries at r_s = 1: the black bounce, a = r_s/2, whose horizons are
# r = +-sqrt(3)/2, the one way wormhole, a = r_s, and the traversable wormhole, a = 2 r_s.
SV_CASES = {"bounce": ("$a = r_s/2$", "1/2"), "null": ("$a = r_s$", "1"), "wormhole": ("$a = 2\\,r_s$", "2")}

# Hayward's regular black hole at ell = 12m/(7 sqrt 7) = 0.648 m, where r^3 - 2m r^2 + 2m ell^2 =
# (r - 6m/7)(r - 12m/7)(r + 4m/7): the horizons are r_- = 6m/7 and r_+ = 12m/7, with surface
# gravities -5/(12m) and 1/(6m). The extremal hole has ell = 4m/(3 sqrt 3) = 0.770 m.
HAYWARD = {"m": 1, "ell": "12/(7*sqrt(7))"}
# Bonnor and Vaidya's charged shell: its charge in units of its mass M, the value of Reissner-Nordstrom's
# drawings, r_q = 0.48 r_s, and the radius q^2/2M at which the shell's energy density changes sign.
BV_CHARGE = "Rational(24, 25)"
BV_TURN = "288/625"
BV_INPUT = {
    "ingoing": "$m(v) = 0$ and $q(v) = 0$ for $v < 0$, and $m = M$ and $q = 0.96\\,M$ for $v > 0$: a shell of "
               "charged null dust of mass $M$ falling in along $v = 0$.",
    "outgoing": "$m(u) = M$ and $q(u) = 0.96\\,M$ for $u < 0$, and $m = 0$ and $q = 0$ for $u > 0$: a shell of "
                "charged null dust carrying off the whole mass $M$ and the whole charge along $u = 0$.",
}
# The same length with the mass a function of advanced time, in units of its greatest value m_0:
# it grows as sin^2 from v = 0 to 2, stays m_0 until v = 4 and falls as cos^2 to zero at v = 8.
# Trapped spheres exist while m > m_* = 3 sqrt(3) ell/4 = 0.842 m_0, from v = 1.479 to 5.042.
HAYWARD_MASS = "sin(pi*Min(Max(v, 0), 2)/4)**2*cos(pi*(Min(Max(v, 4), 8) - 4)/8)**2"
HAYWARD_INPUT = ("$m(v) = m_0\\sin^2(\\pi v/4m_0)$ from $v = 0$ to $2\\,m_0$, $m_0$ until $v = 4\\,m_0$, and "
                 "$m_0\\cos^2(\\pi(v - 4m_0)/8m_0)$ until $v = 8\\,m_0$, with no mass before or after: "
                 "ingoing radiation of positive energy forms the black hole, and ingoing radiation of "
                 "negative energy evaporates it.")

# Two of Majumdar and Papapetrou's holes, each of mass parameter m, the unit, at z = +-2m.
MP_TWO = "1 + 1/sqrt(x**2 + y**2 + (z - 2)**2) + 1/sqrt(x**2 + y**2 + (z + 2)**2)"
MP_TWO_CYLINDRICAL = "1 + 1/sqrt(rho**2 + (z - 2)**2) + 1/sqrt(rho**2 + (z + 2)**2)"
MP_TWO_INPUT = ("Two holes, each of mass parameter $m$, on the axis at $z = \\pm 2m$: "
                "$U = 1 + m/\\sqrt{x^2 + y^2 + (z - 2m)^2} + m/\\sqrt{x^2 + y^2 + (z + 2m)^2}$.")

# Israel, Wilson and Perjes's two sources of Hartle and Hawking's (4.21), U = 1 + (m - il)/r_1 + (m + il)/r_2
# with r_1 and r_2 the distances from z = +2m and z = -2m on the axis, at m = 1 and l = 1/2: W = |U|, and
# omega solves its equation with omega = 0 on the axis beyond the sources, where it is 4l(1 + m/d) = 3
# between them; israel_wilson_perjes.md derives it.
_IWP_R1, _IWP_R2 = "sqrt(rho**2 + (z - 2)**2)", "sqrt(rho**2 + (z + 2)**2)"
IWP_TWO = {"W": f"sqrt((1 + 1/{_IWP_R1} + 1/{_IWP_R2})**2 + (1/{_IWP_R2} - 1/{_IWP_R1})**2/4)",
           "omega": f"(z + 2)/{_IWP_R2} - (z - 2)/{_IWP_R1} - ((rho**2 + z**2 - 4)/({_IWP_R1}*{_IWP_R2}) - 1)/2"}
IWP_TWO_INPUT = ("Two sources on the axis, $U = 1 + (m - il)/r_1 + (m + il)/r_2$ with $l = m/2$ and $r_1$ and $r_2$ "
                 "the distances from $z = 2m$ and $z = -2m$, so that $W = |U|$ and $\\omega = 2l(z + 2m)/r_2 - "
                 "2l(z - 2m)/r_1 - (l/2)\\left((\\rho^2 + z^2 - 4m^2)/(r_1r_2) - 1\\right)$, which vanishes on "
                 "the axis beyond the sources.")

# Two of Kastor and Traschen's holes, each of mass parameter m, the unit, at z = +-2m, falling
# together at H = -3c/(32m), and one hole alone at H = -3c/(16m), the lukewarm hole; slices.py
# holds the rates and the last ray of each plane to reach infinity, the event horizon there.
KT_TWO = "1/sqrt(x**2 + y**2 + (z - 2)**2) + 1/sqrt(x**2 + y**2 + (z + 2)**2)"
KT_TWO_CYLINDRICAL = "1/sqrt(rho**2 + (z - 2)**2) + 1/sqrt(rho**2 + (z + 2)**2)"
KT_PARAMS = {"H": "-3/32"}
KT_TWO_INPUT = ("Two holes, each of mass parameter $m$, on the axis at $z = \\pm 2m$, falling together at "
                "$H = -3c/(32m)$: $V = m/\\sqrt{x^2 + y^2 + (z - 2m)^2} + m/\\sqrt{x^2 + y^2 + (z + 2m)^2}$.")
KT_HORIZON = "the event horizon, the last ray in the plane to reach infinity"


def _kt_point(plane, x, comoving=False):
    """An event of the last ray of a plane of the two holes to reach infinity, at the distance x."""
    tau = float(slices.kt_last_ray(plane)(abs(x)))
    return {"x0": repr(float(slices.kt_comoving_t(tau)) if comoving else tau), "r": repr(float(x))}


def _kt_rows():
    rows = []
    for system, time, ylabel, box, comoving in (("cartesian", "\\tau", "$c\\tau/m$", (-12, 12, -12, 12), False),
                                                ("comoving", "t", "$ct/m$", (-12, 12, 2, 26), True)):
        common = dict(families=SIDEWAYS, functions={"V": KT_TWO}, input=KT_TWO_INPUT,
                      tau="t" if comoving else "tau")
        if not comoving:
            common.update(singular_zero="H*tau + V", ends_on_singular=True)
        rows.append(Diagram("kastor_traschen", system, "tz", "the axis", (time, "z"), box, "$z/m$", ylabel,
                            KT_PARAMS, {"x": "0", "y": "0"}, **common,
                            lines=(("grr", "r", "2", None), ("grr", "r", "-2", None)),
                            marked=(("event", _kt_point("axis", 6, comoving), 1, KT_HORIZON),
                                    ("event", _kt_point("axis", -6, comoving), 0, KT_HORIZON))))
        rows.append(Diagram("kastor_traschen", system, "tx", "the midplane", (time, "x"), box, "$x/m$", ylabel,
                            KT_PARAMS, {"y": "0", "z": "0"}, **common,
                            marked=(("event", _kt_point("midplane", 0, comoving), "both", KT_HORIZON, "future"),)))
        if not comoving:
            rows.append(Diagram("kastor_traschen", "cylindrical", "radial", "$\\tau$ and $\\rho$", ("\\tau", "\\rho"),
                                (0, 12, -9, 3), "$\\rho/m$", "$c\\tau/m$", KT_PARAMS, {"phi": "0", "z": "0"},
                                functions={"V": KT_TWO_CYLINDRICAL}, input=KT_TWO_INPUT, tau="tau",
                                singular_zero="H*tau + V", ends_on_singular=True,
                                marked=(("event", _kt_point("midplane", 4), 1, KT_HORIZON),)))
            rows.append(Diagram("kastor_traschen", "isotropic", "radial", "$\\tau$ and $r$", ("\\tau", "r"),
                                (0, 8, -4, 4), "$r/m$", "$c\\tau/m$", {"m": 1, "H": "-3/16"}, EQUATOR, tau="tau",
                                areal=True, no_throat=True, singular_zero="H*tau*r + m", ends_on_singular=True))
    return rows


# Every view the page draws, in the order it shows them. Plot ranges are chosen with
# equal scales on both axes, so light in flat space runs at 45 degrees, and cone
# lattices so that no cone sits exactly on a line where the chart is singular.

def teo_inverse(x):
    """(rho, sigma) at x = l/b_0 for the wormhole of Teo's example, whose proper radial distance
    is l = +-b_0 (sqrt(rho(rho - 1)) + ln(sqrt(rho) + sqrt(rho - 1))) with rho = r/b_0, his eq.
    (28). With rho = cosh^2 w that is Kepler's kind of equation, x = w + sinh(2w)/2, solved for w
    by Newton's method from w = asinh(2x)/2, which it meets for large |x|; then rho = cosh^2 w
    and sigma = d rho/dx = tanh w, odd and smooth through the throat."""
    x = np.asarray(x, dtype=float)
    w = np.arcsinh(2 * x) / 2
    for _ in range(60):
        step = (w + np.sinh(2 * w) / 2 - x) / (2 * np.cosh(w) ** 2)
        w = w - step
        if np.max(np.abs(step), initial=0.0) < 1e-15:
            break
    return np.cosh(w) ** 2, np.tanh(w)


class teo_rho(sp.Function):
    """r/b_0 as a function of l/b_0 on Teo's wormhole, a declared function a row may name: sympy
    differentiates it by d rho/dx = sigma and d sigma/dx = 1/(2 rho^2), which are dr/dl =
    +-sqrt(1 - b_0/r) and its derivative, and lambdify evaluates it with teo_inverse."""
    nargs = 1
    is_real = True
    _imp_ = staticmethod(lambda x: teo_inverse(x)[0])

    def fdiff(self, argindex=1):
        return teo_sigma(self.args[0])


class teo_sigma(sp.Function):
    nargs = 1
    is_real = True
    _imp_ = staticmethod(lambda x: teo_inverse(x)[1])

    def fdiff(self, argindex=1):
        return 1 / (2 * teo_rho(self.args[0]) ** 2)


class ori_mass_behind(sp.Function):
    """The mass function behind Ori's shell as a function of the advanced time there, which is
    negative and reaches zero on the Cauchy horizon: ori_shell.py integrates it, and from zero on,
    where there is no spacetime, it is infinite."""
    nargs = 1
    is_real = True

    @staticmethod
    def _imp_(v):
        v = np.asarray(v, dtype=float)
        with np.errstate(all="ignore"):
            return np.where(v < 0, ori_shell.shell().mass_behind(np.minimum(v, -1e-300)), np.inf)

    def fdiff(self, argindex=1):
        return ori_influx_behind(self.args[0])


class ori_influx_behind(sp.Function):
    nargs = 1
    is_real = True

    @staticmethod
    def _imp_(v):
        v = np.asarray(v, dtype=float)
        with np.errstate(all="ignore"):
            return np.where(v < 0, ori_shell.shell().influx_behind(np.minimum(v, -1e-300)), np.inf)


class ori_shell_behind(sp.Function):
    """The radius of Ori's shell on the ingoing ray of advanced time v behind it."""
    nargs = 1
    is_real = True

    @staticmethod
    def _imp_(v):
        v = np.asarray(v, dtype=float)
        return ori_shell.shell().radius_behind(np.minimum(v, -1e-300))


# Hiscock's evaporating black hole as every one of its diagrams draws it, in units of the mass m_0
# of the shell of light that makes the hole, which falls in along v = 0. From v_1 = 2 to v_0 = 8
# negative energy falls in and the mass inside the surface of pair creation is N(v) =
# cos^2(pi (v - 2)/12), which starts and ends with no slope, so that dN/dv -> 0 as N -> 0, the
# condition Hiscock found necessary. The surface is r = R = 3 N(v), outside the apparent horizon
# r = 2 N by half its radius; on it g^rr = 1/3, and the retarded time of the outgoing chart beyond
# it, continuous in the induced metric, has du/dv = 1 - 2 R'/g^rr = 1 - 18 N', so u = v - 18 N(v),
# normalised to u = v - 2r in the flat space left at the end: u_1 = -16 and u_0 = v_0 = 8. Beyond
# the surface the mass is M(u) = N(v(u)), the same on both sides of it, so the surface has no
# surface energy density; hiscock.md derives the tension it has. Before u_1 the outgoing chart is
# Schwarzschild's, u = v - 2 r_* - 12 - 4 ln 2 with r_* = r + 2 ln(r/2 - 1), and begins on the shell.
HISCOCK_V1, HISCOCK_V0, HISCOCK_DELAY = 2.0, 8.0, 18.0
HISCOCK_SHIFT = 12 + 4 * math.log(2)        # v - 2 r_* - u in Schwarzschild's region
HISCOCK_MASS = "Heaviside(v)*cos(pi*(Min(Max(v, 2), 8) - 2)/12)**2"
HISCOCK_EDGE = "3*cos(pi*(Min(Max(v, 2), 8) - 2)/12)**2 + 1000*Heaviside(2 - v)"
HISCOCK_HORIZON = {"x0": "31/4", "r": "0.00766033576251"}       # the event horizon at v = 7.75
HISCOCK_UNITS = "with $m_0$ the mass of the shell and the unit of every length"
HISCOCK_INPUT = ("$m(v) = 0$ before $v = 0$, $m_0$ until $v_1 = 2\\,m_0$, and $m_0\\cos^2(\\pi(v - v_1)/12m_0)$ until "
                 "$v_0 = 8\\,m_0$: a shell of light of mass $m_0$ makes the hole, and negative energy falling in "
                 "evaporates it. The surface of pair creation is $R = 3m$.")
HISCOCK_OUT_INPUT = ("$m(u)$ is the mass of the ingoing chart on the surface of pair creation $R = 3m$, where "
                     "$u = v - 18m$: it is $m_0$ until $u_1 = -16\\,m_0$ and zero from $u_0 = 8\\,m_0$, and the "
                     "radiation between carries positive energy outward.")
HISCOCK_SHELLS_INPUT = ("$m(v) = m_0$ from $v = 0$ to $v_1 = 6\\,m_0$ and zero before and after: a shell of light of "
                        "mass $m_0$ makes the hole and a shell of mass $-m_0$ removes it.")
HISCOCK_FLAT_INPUT = ("The flat space left after the evaporation, beyond the last ray $u_0 = 8\\,m_0$ of the outgoing "
                      "chart, " + HISCOCK_UNITS + ".")


def hiscock_mass(v):
    """N(v), the mass inside the surface of pair creation, HISCOCK_MASS as numbers."""
    v = np.asarray(v, dtype=float)
    return np.where(v > 0, np.cos(np.pi * (np.clip(v, HISCOCK_V1, HISCOCK_V0) - HISCOCK_V1) / 12) ** 2, 0.0)


def _hiscock_table():
    v = np.linspace(HISCOCK_V1, HISCOCK_V0, 4001)
    return v - HISCOCK_DELAY * hiscock_mass(v), v


def hiscock_advanced(u):
    """The advanced time v at which the outgoing ray u leaves the surface of pair creation, the
    root of v - 18 N(v) = u between v_1 and v_0, by Newton's rule from a table of the map, which
    rises everywhere; u_1 = -16 and u_0 = 8 are its ends."""
    u = np.asarray(u, dtype=float)
    table = hiscock_advanced.__dict__.setdefault("table", _hiscock_table())
    v = np.interp(u, *table)
    for _ in range(6):
        angle = np.pi * (v - HISCOCK_V1) / 6
        v = np.clip(v - (v - HISCOCK_DELAY * np.cos(angle / 2) ** 2 - u) / (1 + HISCOCK_DELAY * np.pi / 12 * np.sin(angle)),
                    HISCOCK_V1, HISCOCK_V0)
    return v


def hiscock_out(u, order=0):
    """M(u), the mass of the outgoing chart, and with order 1 or 2 its first or second derivative:
    M = N(v(u)), M' = N'/(1 - 18 N') and M'' = N''/(1 - 18 N')^3."""
    u = np.asarray(u, dtype=float)
    v = hiscock_advanced(u)
    during = (u > HISCOCK_V1 - HISCOCK_DELAY) & (u < HISCOCK_V0)
    angle = np.pi * (v - HISCOCK_V1) / 6
    d1, d2 = -np.pi / 12 * np.sin(angle), -np.pi ** 2 / 72 * np.cos(angle)
    if order == 0:
        return np.where(u <= HISCOCK_V1 - HISCOCK_DELAY, 1.0, np.where(during, hiscock_mass(v), 0.0))
    stretch = 1 - HISCOCK_DELAY * d1
    return np.where(during, d1 / stretch if order == 1 else d2 / stretch ** 3, 0.0)


def hiscock_edge(u):
    """The inner edge of the outgoing chart at the retarded time u: the collapsing shell v = 0
    before u_1, r = 2(1 + W(exp(-(u + 12 + 4 ln 2)/4 - 1))) with W Lambert's function, and the
    surface of pair creation R = 3 M(u) from u_1 on."""
    from scipy.special import lambertw
    u = np.asarray(u, dtype=float)
    early = np.minimum(u, HISCOCK_V1 - HISCOCK_DELAY)
    shell = 2 * (1 + lambertw(np.exp(-(early + HISCOCK_SHIFT) / 4 - 1)).real)
    return np.where(u < HISCOCK_V1 - HISCOCK_DELAY, shell, 3 * hiscock_out(u))


def hiscock_meets(T):
    """The advanced time at which the slice v - r = T meets the surface of pair creation, the root
    of v - 3N(v) = T between v_1 and v_0, for T from -1 to 8."""
    lo, hi = HISCOCK_V1, HISCOCK_V0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if mid - 3 * float(hiscock_mass(mid)) < T else (lo, mid)
    return 0.5 * (lo + hi)


def hiscock_outer_time(T):
    """The u + r of the slice outside the surface of pair creation that meets the slice v - r = T
    on it: there u = v - 18N and r = 3N, so u + r = v - 15N(v) at v = hiscock_meets(T). From T = 8
    on, when the hole is gone, it is T itself."""
    if T >= HISCOCK_V0:
        return T
    v = hiscock_meets(T)
    return v - 15 * float(hiscock_mass(v))


class hiscock_m_out(sp.Function):
    """M(u) of Hiscock's outgoing chart, a declared function a row may name: sympy differentiates
    it by hiscock_out's own derivatives, and lambdify evaluates it with hiscock_out."""
    nargs = 1
    is_real = True
    _imp_ = staticmethod(hiscock_out)

    def fdiff(self, argindex=1):
        return hiscock_dm_out(self.args[0])


class hiscock_dm_out(sp.Function):
    nargs = 1
    is_real = True
    _imp_ = staticmethod(lambda u: hiscock_out(u, 1))

    def fdiff(self, argindex=1):
        return hiscock_d2m_out(self.args[0])


class hiscock_d2m_out(sp.Function):
    nargs = 1
    is_real = True
    _imp_ = staticmethod(lambda u: hiscock_out(u, 2))


class hiscock_r_out(sp.Function):
    """The inner edge of Hiscock's outgoing chart, hiscock_edge, as a declared function."""
    nargs = 1
    is_real = True
    _imp_ = staticmethod(hiscock_edge)


# Functions a row's `functions` may name beside the elementary ones, each a sympy function
# that carries its own derivative and its own numbers.
DECLARED_FUNCTIONS = {"teo_rho": teo_rho, "teo_sigma": teo_sigma, "ori_mass_behind": ori_mass_behind,
                      "ori_influx_behind": ori_influx_behind, "ori_shell_behind": ori_shell_behind,
                      "hiscock_m_out": hiscock_m_out, "hiscock_dm_out": hiscock_dm_out,
                      "hiscock_d2m_out": hiscock_d2m_out, "hiscock_r_out": hiscock_r_out}
TEO_RADIUS = "b_0*teo_rho(l/b_0)"
TEO_RADIUS_INPUT = ("$r(l)$ from Teo's $l = \\pm\\left(\\sqrt{r(r - b_0)} + b_0\\ln\\left(\\sqrt{r/b_0} + "
                    "\\sqrt{r/b_0 - 1}\\right)\\right)$, inverted by Newton's method.")
TEO_CONE = "future cone of no angular momentum"

# Hartle and Thorne's star as every one of its diagrams draws it, in units of m = GM/c^2: the
# surface at R = 6 m, 12.4 km for 1.4 solar masses, the spin a = m/4, and a quadrupole moment
# four times Kerr's, q = 4a^2, inside the range Laarakkers and Poisson found for neutron stars.
HT_STAR = {"m": 1, "a": "1/4", "q": "1/4", "R": 6}
HT_SIMPLE = {"m": 1, "a": "1/4", "R": 6}
HT_AXIS = {"theta": "0", "phi": "0"}
HT_CONE = "future cone of no angular momentum"


def ht_functions(r, a=0.25, q=0.25):
    """Hartle and Thorne's F, h_2, j_2 and k_2 at m = 1, written out here from their equation (A1)
    and never read from the file, so that the drawings can be checked against them."""
    r = np.asarray(r, dtype=float)
    L = np.log(r / (r - 2))
    A = 15 * r * (r - 2) / 16 * L - 5 * (r - 1) * (3 * r ** 2 - 6 * r - 2) / (8 * r * (r - 2))
    B = 5 * (3 * r ** 2 - 6 * r + 1) / (4 * r * (r - 2)) - 15 * (r - 1) / 8 * L
    F = 1 - 2 / r + 2 * a ** 2 / r ** 4
    h2 = a ** 2 / r ** 3 * (1 + 1 / r) + (q - a ** 2) * A
    j2 = a ** 2 / r ** 3 * (1 - 5 / r) + (q - a ** 2) * A
    k2 = -a ** 2 / r ** 3 * (1 + 2 / r) + (q - a ** 2) * (B - A)
    return F, h2, j2, k2


def _ht_speed(r, P2):
    F, h2, j2, _ = ht_functions(r)
    return np.sqrt((1 - 2 * j2 * P2) / (1 + 2 * h2 * P2)) / F


@lru_cache(maxsize=None)
def _ht_table(P2):
    r = 6 + np.concatenate([[0.0], np.geomspace(1e-9, 1e4, 400001)])
    speed = _ht_speed(r, P2)
    return r, speed, np.concatenate([[0.0], np.cumsum(0.5 * (speed[1:] + speed[:-1]) * np.diff(r))])


def ht_star(r, P2):
    """The tortoise coordinate of Hartle and Thorne's chart on the axis, P_2 = 1, or on the equator
    with phi divided out, P_2 = -1/2, at the declared star: the quadrature of sqrt(-g_rr/g_tt) =
    sqrt((1 - 2j_2 P_2)/(1 + 2h_2 P_2))/F from the surface r = 6m, where it vanishes, by the
    trapezoid rule on a table and on the last part of a cell. Past r = 10^4 m, where the terms of
    second order are below 10^-13, it runs on as Schwarzschild's."""
    r = np.asarray(r, dtype=float)
    grid, speed, table = _ht_table(P2)
    top = grid[-1]
    inside = np.clip(r, 6.0, top)
    i = np.clip(np.searchsorted(grid, inside, side="right") - 1, 0, grid.size - 1)
    near = table[i] + 0.5 * (speed[i] + _ht_speed(inside, P2)) * (inside - grid[i])
    far = table[-1] + (r - top) + 2 * np.log(np.maximum(r, top + 1) - 2) - 2 * np.log(top - 2)
    return np.where(r <= top, near, far)


def pg_out(r):
    """The integral of dr/(1 - sqrt(2m/r)) at m = 1, which t minus it keeps along an outgoing ray of
    the Painleve-Gullstrand chart, zero at the surface r = 6m."""
    f = lambda x: x + 2 * np.sqrt(2 * x) + 4 * np.log(np.sqrt(x / 2) - 1)
    return f(np.asarray(r, dtype=float)) - f(6.0)


def pg_in(r):
    """The integral of dr/(1 + sqrt(2m/r)) at m = 1, which t plus it keeps along an ingoing ray."""
    f = lambda x: x - 2 * np.sqrt(2 * x) + 4 * np.log(np.sqrt(x / 2) + 1)
    return f(np.asarray(r, dtype=float)) - f(6.0)


# Einstein's cluster as its diagrams draw it: the declared redshift function of the cluster with no
# surface, in units of its Schwarzschild radius with the core radius b = 2 r_s, and the isotropic
# radius and the polar angle of the surface of the clusters with R = 3 r_s.
EC_CORE = "log(1 - 1/sqrt(r**2 + 4))/2"
EC_CORE_INPUT = ("A cluster with no surface, $e^{2\\Phi} = 1 - r_s/\\sqrt{r^2 + b^2}$, with $r_s$ the Schwarzschild "
                 "radius of its whole mass and the core radius $b = 2\\,r_s$, so that $m \\to r_s/2$ far away and the "
                 "fastest particles, near $r = 2.3\\,r_s$, move at $0.37\\,c$.")
EC_RHO0 = 3 / (6 - 2 * math.sqrt(6)) ** 2
EC_CHI0 = math.asin(1 / math.sqrt(3))

DIAGRAMS = [
    *[Diagram("aichelburg_sexl", "null_cartesian", view, f"$\\rho = \\rho_0/{rho[2:]}$", ("u", "v"), (-3, 3, -2.5, 3.5),
              "$z\\;[8GE/c^4]$", "$ct\\;[8GE/c^4]$", {"G": 1, "E": "1/8", "rho_0": 1}, {"x": rho, "y": "0"},
              to_display=NULL_TO_TR, tau="u + v", families=SIDEWAYS, delta=AS_PULSE, step=0.0005,
              lines=(("shell", "x0", "0", "the shock, $u = 0$"),), input=AS_INPUT)
      for view, rho in AS_RHO.items()],
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
    Diagram("schwarzschild_de_sitter", "static", "radial", "$t$ and $r$", ("t", "r"), (0, 4, -2, 2),
            "$r/r_s$", "$ct/r_s$", SDS, EQUATOR, orient="split", split=SDS_STATIC, areal=True),
    Diagram("schwarzschild_de_sitter", "eddington_finkelstein_ingoing", "finkelstein", "against $v - r$",
            ("v", "r"), (0, 4, -2, 2), "$r/r_s$", "$(v - r)/r_s$", SDS, EQUATOR,
            to_display=FINKELSTEIN_IN, tau="v - r", areal=True),
    Diagram("schwarzschild_de_sitter", "eddington_finkelstein_ingoing", "chart", "against $v$",
            ("v", "r"), (0, 4, 0, 4), "$r/r_s$", "$v/r_s$", SDS, EQUATOR, tau="v - r", areal=True),
    Diagram("schwarzschild_de_sitter", "eddington_finkelstein_outgoing", "finkelstein", "against $u + r$",
            ("u", "r"), (0, 4, -2, 2), "$r/r_s$", "$(u + r)/r_s$", SDS, EQUATOR,
            to_display=FINKELSTEIN_OUT, tau="u + r", areal=True),
    Diagram("schwarzschild_de_sitter", "eddington_finkelstein_outgoing", "chart", "against $u$",
            ("u", "r"), (0, 4, -4, 0), "$r/r_s$", "$u/r_s$", SDS, EQUATOR, tau="u + r", areal=True),
    Diagram("kiselev", "static", "radial", "$t$ and $r$", ("t", "r"), (0, 9, -4.5, 4.5),
            "$r/r_s$", "$ct/r_s$", KISELEV, EQUATOR, orient="split", split=KISELEV_STATIC, areal=True),
    Diagram("kiselev", "eddington_finkelstein_ingoing", "finkelstein", "against $v - r$",
            ("v", "r"), (0, 9, -4.5, 4.5), "$r/r_s$", "$(v - r)/r_s$", KISELEV, EQUATOR,
            to_display=FINKELSTEIN_IN, tau="v - r", areal=True),
    Diagram("kiselev", "eddington_finkelstein_ingoing", "chart", "against $v$",
            ("v", "r"), (0, 9, 0, 9), "$r/r_s$", "$v/r_s$", KISELEV, EQUATOR, tau="v - r", areal=True),
    Diagram("kiselev", "eddington_finkelstein_outgoing", "finkelstein", "against $u + r$",
            ("u", "r"), (0, 9, -4.5, 4.5), "$r/r_s$", "$(u + r)/r_s$", KISELEV, EQUATOR,
            to_display=FINKELSTEIN_OUT, tau="u + r", areal=True),
    Diagram("kiselev", "eddington_finkelstein_outgoing", "chart", "against $u$",
            ("u", "r"), (0, 9, -9, 0), "$r/r_s$", "$u/r_s$", KISELEV, EQUATOR, tau="u + r", areal=True),
    Diagram("kiselev", "linear", "radial", "$t$ and $r$", ("t", "r"), (0, 9, -4.5, 4.5),
            "$r/r_s$", "$ct/r_s$", KISELEV_LINEAR, EQUATOR, orient="split", split=KISELEV_STATIC, areal=True),
    Diagram("kiselev", "hyperbolic", "radial", "$\\eta$ and $\\chi$", ("\\eta", "\\chi"), (0, 4, -2, 2),
            "$\\chi$", "$\\eta$", KISELEV_FREE, EQUATOR, tau="eta", areal=True),
    Diagram("kiselev", "conformally_flat", "radial", "$\\tau$ and $\\rho$", ("\\tau", "\\rho"), (0, 4, -2, 2),
            "$\\rho$", "$\\tau$", KISELEV_FREE, EQUATOR, tau="tau", areal=True, where="tau + rho",
            where_is_infinity=True, no_throat=True, singular_runs=True, singular_where_claimed=True,
            singular_near=1e-9),
    Diagram("reissner_nordstrom_de_sitter", "static", "radial", "$t$ and $r$", ("t", "r"), (0, 2.6, -1.3, 1.3),
            "$r/r_s$", "$ct/r_s$", RNDS, EQUATOR, orient="split", split=RNDS_STATIC, areal=True),
    Diagram("reissner_nordstrom_de_sitter", "eddington_finkelstein_ingoing", "finkelstein", "against $v - r$",
            ("v", "r"), (0, 2.6, -1.3, 1.3), "$r/r_s$", "$(v - r)/r_s$", RNDS, EQUATOR,
            to_display=FINKELSTEIN_IN, tau=RNDS_TIME_IN, areal=True),
    Diagram("reissner_nordstrom_de_sitter", "eddington_finkelstein_ingoing", "chart", "against $v$",
            ("v", "r"), (0, 2.6, 0, 2.6), "$r/r_s$", "$v/r_s$", RNDS, EQUATOR, tau=RNDS_TIME_IN, areal=True),
    Diagram("reissner_nordstrom_de_sitter", "eddington_finkelstein_outgoing", "finkelstein", "against $u + r$",
            ("u", "r"), (0, 2.6, -1.3, 1.3), "$r/r_s$", "$(u + r)/r_s$", RNDS, EQUATOR,
            to_display=FINKELSTEIN_OUT, tau=RNDS_TIME_OUT, areal=True),
    Diagram("reissner_nordstrom_de_sitter", "eddington_finkelstein_outgoing", "chart", "against $u$",
            ("u", "r"), (0, 2.6, -2.6, 0), "$r/r_s$", "$u/r_s$", RNDS, EQUATOR, tau=RNDS_TIME_OUT, areal=True),
    Diagram("reissner_nordstrom_de_sitter", "cosmological", "plane", "$\\tau$ and $\\rho$", ("\\tau", "\\rho"),
            (0, 3, -2, 4), "$\\rho/r_s$", "$c\\tau/r_s$", {"r_s": 1, "H": "3/8"}, EQUATOR, tau="tau", areal=True,
            no_throat=True, singular_zero="2*H*tau*rho + r_s"),
    Diagram("reissner_nordstrom_ads", "static", "radial", "$t$ and $r$", ("t", "r"), (0, 3, -1.5, 1.5),
            "$r/L$", "$ct/L$", RNADS, EQUATOR, orient="ingoing", areal=True),
    Diagram("reissner_nordstrom_ads", "eddington_finkelstein_ingoing", "finkelstein", "against $v - r$",
            ("v", "r"), (0, 3, -3, 0), "$r/L$", "$(v - r)/L$", RNADS, EQUATOR,
            to_display=FINKELSTEIN_IN, tau=RNADS_TIME_IN, areal=True),
    Diagram("reissner_nordstrom_ads", "eddington_finkelstein_ingoing", "chart", "against $v$",
            ("v", "r"), (0, 3, -2.5, 0.5), "$r/L$", "$v/L$", RNADS, EQUATOR, tau=RNADS_TIME_IN, areal=True),
    Diagram("reissner_nordstrom_ads", "eddington_finkelstein_outgoing", "finkelstein", "against $u + r$",
            ("u", "r"), (0, 3, 0, 3), "$r/L$", "$(u + r)/L$", RNADS, EQUATOR,
            to_display=FINKELSTEIN_OUT, tau=RNADS_TIME_OUT, areal=True),
    Diagram("reissner_nordstrom_ads", "eddington_finkelstein_outgoing", "chart", "against $u$",
            ("u", "r"), (0, 3, -0.5, 2.5), "$r/L$", "$u/L$", RNADS, EQUATOR, tau=RNADS_TIME_OUT, areal=True),
    Diagram("tangherlini", "spherical", "radial", "$t$ and $r$", ("t", "r"), (0, 6, -3, 3),
            "$r/r_h$", "$ct/r_h$", {"r_h": 1}, TANGHERLINI_FIVE, orient="ingoing", areal=True),
    Diagram("tangherlini", "eddington_finkelstein_ingoing", "finkelstein", "against $v - r$",
            ("v", "r"), (0, 6, -3, 3), "$r/r_h$", "$(v - r)/r_h$", {"r_h": 1}, TANGHERLINI_FIVE,
            to_display=FINKELSTEIN_IN, tau="v - r", areal=True),
    Diagram("tangherlini", "eddington_finkelstein_ingoing", "chart", "against $v$",
            ("v", "r"), (0, 6, 0, 6), "$r/r_h$", "$v/r_h$", {"r_h": 1}, TANGHERLINI_FIVE,
            tau="v - r", areal=True),
    Diagram("tangherlini", "eddington_finkelstein_outgoing", "finkelstein", "against $u + r$",
            ("u", "r"), (0, 6, -3, 3), "$r/r_h$", "$(u + r)/r_h$", {"r_h": 1}, TANGHERLINI_FIVE,
            to_display=FINKELSTEIN_OUT, tau="u + r", areal=True),
    Diagram("tangherlini", "eddington_finkelstein_outgoing", "chart", "against $u$",
            ("u", "r"), (0, 6, -6, 0), "$r/r_h$", "$u/r_h$", {"r_h": 1}, TANGHERLINI_FIVE,
            tau="u + r", areal=True),
    Diagram("tangherlini", "spherical_six", "radial", "$t$ and $r$", ("t", "r"), (0, 6, -3, 3),
            "$r/r_h$", "$ct/r_h$", {"r_h": 1}, TANGHERLINI_SIX, orient="ingoing", areal=True),
    Diagram("boulware_deser", "spherical", "radial", "$t$ and $r$", ("t", "r"), (0, 6, -3, 3),
            "$r/r_h$", "$ct/r_h$", BD, TANGHERLINI_FIVE, orient="ingoing", areal=True),
    Diagram("boulware_deser", "eddington_finkelstein_ingoing", "finkelstein", "against $v - r$",
            ("v", "r"), (0, 6, -3, 3), "$r/r_h$", "$(v - r)/r_h$", BD, TANGHERLINI_FIVE,
            to_display=FINKELSTEIN_IN, tau="v - r", areal=True),
    Diagram("boulware_deser", "eddington_finkelstein_ingoing", "chart", "against $v$",
            ("v", "r"), (0, 6, 0, 6), "$r/r_h$", "$v/r_h$", BD, TANGHERLINI_FIVE,
            tau="v - r", areal=True),
    Diagram("boulware_deser", "eddington_finkelstein_outgoing", "finkelstein", "against $u + r$",
            ("u", "r"), (0, 6, -3, 3), "$r/r_h$", "$(u + r)/r_h$", BD, TANGHERLINI_FIVE,
            to_display=FINKELSTEIN_OUT, tau="u + r", areal=True),
    Diagram("boulware_deser", "eddington_finkelstein_outgoing", "chart", "against $u$",
            ("u", "r"), (0, 6, -6, 0), "$r/r_h$", "$u/r_h$", BD, TANGHERLINI_FIVE,
            tau="u + r", areal=True),
    Diagram("boulware_deser", "spherical_plus", "radial", "$t$ and $r$", ("t", "r"), (0, 3, -1.5, 1.5),
            "$r/\\ell$", "$ct/\\ell$", BD_PLUS, TANGHERLINI_FIVE, areal=True),
    Diagram("black_string", "static", "radial", "$t$ and $r$", ("t", "r"), (0, 6, -3, 3),
            "$r/r_s$", "$ct/r_s$", {"r_s": 1}, STRING_FIVE, orient="ingoing", areal=True),
    Diagram("black_string", "eddington_finkelstein_ingoing", "finkelstein", "against $v - r$",
            ("v", "r"), (0, 6, -3, 3), "$r/r_s$", "$(v - r)/r_s$", {"r_s": 1}, STRING_FIVE,
            to_display=FINKELSTEIN_IN, tau="v - r", areal=True),
    Diagram("black_string", "eddington_finkelstein_ingoing", "chart", "against $v$",
            ("v", "r"), (0, 6, 0, 6), "$r/r_s$", "$v/r_s$", {"r_s": 1}, STRING_FIVE,
            tau="v - r", areal=True),
    Diagram("black_string", "kerr_schild", "radial", "$T$ and $r$", ("T", "r"), (0, 6, -3, 3),
            "$r/r_s$", "$cT/r_s$", {"r_s": 1}, STRING_FIVE, tau="T", areal=True),
    Diagram("black_string", "static_six", "radial", "$t$ and $r$", ("t", "r"), (0, 6, -3, 3),
            "$r/r_h$", "$ct/r_h$", {"r_h": 1}, STRING_SIX, orient="ingoing", areal=True),
    # Myers and Perry's black hole. One spin in five dimensions, a = 0.6 sqrt(mu), where the horizon is
    # 0.8 sqrt(mu): the plane transverse to the rotation, whose rays are null geodesics, and the plane
    # of rotation with phi divided out, in the Boyer-Lindquist chart and the ingoing one. Equal spins,
    # a = 0.4 sqrt(mu), with the Hopf fibre divided out, and one spin in six dimensions at a = 1.5 mu^(1/3),
    # past every spin a hole of five dimensions can have. The chart of two unequal spins has no view:
    # where one plane of rotation lies, the other plane's circle has shrunk to a point, and off both
    # planes a ray moves in theta.
    Diagram("myers_perry", "boyer_lindquist", "transverse", "the transverse plane", ("t", "r"), (0, 3, -1.5, 1.5),
            "$r/\\sqrt{\\mu}$", "$ct/\\sqrt{\\mu}$", MYERS_ONE, {"theta": "0", "phi": "0", "psi": "0"}, orient="ingoing"),
    Diagram("myers_perry", "boyer_lindquist", "rotation", "the plane of rotation", ("t", "r"), (0, 3, -1.5, 1.5),
            "$r/\\sqrt{\\mu}$", "$ct/\\sqrt{\\mu}$", MYERS_ONE, {"theta": "pi/2", "psi": "0"}, orient="ingoing",
            quotient="phi", mark_gtt="the ergosurface", cone=MYERS_CONE),
    Diagram("myers_perry", "ingoing_kerr", "transverse", "the transverse plane", ("v", "r"), (0, 3, -1.5, 1.5),
            "$r/\\sqrt{\\mu}$", "$(v - r)/\\sqrt{\\mu}$", MYERS_ONE, {"theta": "0", "tildephi": "0", "psi": "0"},
            to_display=FINKELSTEIN_IN, tau="v - r"),
    Diagram("myers_perry", "ingoing_kerr", "rotation", "the plane of rotation", ("v", "r"), (0, 3, -1.5, 1.5),
            "$r/\\sqrt{\\mu}$", "$(v - r)/\\sqrt{\\mu}$", MYERS_ONE, {"theta": "pi/2", "psi": "0"},
            to_display=FINKELSTEIN_IN, orient="ingoing", quotient="tildephi", mark_gtt="the ergosurface", cone=MYERS_CONE),
    Diagram("myers_perry", "equal_spins", "radial", "$t$ and $\\rho$", ("t", "\\rho"), (0, 3, -1.5, 1.5),
            "$\\rho/\\sqrt{\\mu}$", "$ct/\\sqrt{\\mu}$", MYERS_EQUAL, {"theta": "pi/2", "phi": "0"}, orient="ingoing",
            quotient="psi", mark_gtt="the ergosurface", cone=MYERS_CONE),
    Diagram("myers_perry", "boyer_lindquist_six", "transverse", "the transverse space", ("t", "r"), (0, 4, -2, 2),
            "$r/\\mu^{1/3}$", "$ct/\\mu^{1/3}$", MYERS_SIX, {"theta": "0", "phi": "0", "chi": "pi/2", "psi": "0"},
            orient="ingoing"),
    Diagram("myers_perry", "boyer_lindquist_six", "rotation", "the plane of rotation", ("t", "r"), (0, 4, -2, 2),
            "$r/\\mu^{1/3}$", "$ct/\\mu^{1/3}$", MYERS_SIX, {"theta": "pi/2", "chi": "pi/2", "psi": "0"}, orient="ingoing",
            quotient="phi", mark_gtt="the ergosurface", cone=MYERS_CONE),
    Diagram("schwarzschild_ads", "static", "radial", "$t$ and $r$", ("t", "r"), (0, 3, -1.5, 1.5),
            "$r/L$", "$ct/L$", SADS, EQUATOR, orient="ingoing", areal=True),
    Diagram("schwarzschild_ads", "eddington_finkelstein_ingoing", "finkelstein", "against $v - r$",
            ("v", "r"), (0, 3, -1.5, 1.5), "$r/L$", "$(v - r)/L$", SADS, EQUATOR,
            to_display=FINKELSTEIN_IN, orient="ingoing", areal=True),
    Diagram("schwarzschild_ads", "eddington_finkelstein_ingoing", "chart", "against $v$",
            ("v", "r"), (0, 3, -1, 2), "$r/L$", "$v/L$", SADS, EQUATOR, orient="ingoing", areal=True),
    Diagram("schwarzschild_ads", "eddington_finkelstein_outgoing", "finkelstein", "against $u + r$",
            ("u", "r"), (0, 3, -1.5, 1.5), "$r/L$", "$(u + r)/L$", SADS, EQUATOR,
            to_display=FINKELSTEIN_OUT, orient="outgoing", areal=True),
    Diagram("schwarzschild_ads", "eddington_finkelstein_outgoing", "chart", "against $u$",
            ("u", "r"), (0, 3, -2, 1), "$r/L$", "$u/L$", SADS, EQUATOR, orient="outgoing", areal=True),
    # Bardeen's regular black hole on its plane of the time and r, in each chart: the cones close at
    # two horizons, as Reissner-Nordstrom's do, and open again down to a centre that is regular.
    Diagram("bardeen", "static", "radial", "$t$ and $r$", ("t", "r"), (0, 2, -1, 1),
            "$r/r_s$", "$ct/r_s$", BARDEEN, EQUATOR, orient="ingoing", areal=True),
    Diagram("bardeen", "eddington_finkelstein_ingoing", "finkelstein", "against $v - r$",
            ("v", "r"), (0, 2, -1, 1), "$r/r_s$", "$(v - r)/r_s$", BARDEEN, EQUATOR,
            to_display=FINKELSTEIN_IN, orient="ingoing", areal=True),
    Diagram("bardeen", "eddington_finkelstein_ingoing", "chart", "against $v$",
            ("v", "r"), (0, 2, -0.5, 1.5), "$r/r_s$", "$v/r_s$", BARDEEN, EQUATOR, orient="ingoing", areal=True),
    Diagram("bardeen", "eddington_finkelstein_outgoing", "finkelstein", "against $u + r$",
            ("u", "r"), (0, 2, -1, 1), "$r/r_s$", "$(u + r)/r_s$", BARDEEN, EQUATOR,
            to_display=FINKELSTEIN_OUT, orient="outgoing", areal=True),
    Diagram("bardeen", "eddington_finkelstein_outgoing", "chart", "against $u$",
            ("u", "r"), (0, 2, -1.5, 0.5), "$r/r_s$", "$u/r_s$", BARDEEN, EQUATOR, orient="outgoing", areal=True),
    # The Kaluza-Klein monopole at m = 1 on its plane of t and the radius, in each chart, and through the
    # nut in Gross and Perry's: no horizon, and rays that slow toward the nut as sqrt(r/(r + 4m)).
    Diagram("kaluza_klein_monopole", "gross_perry", "radial", "$t$ and $r$", ("t", "r"), (0, 16, -8, 8),
            "$r/m$", "$ct/m$", {"m": 1}, KK_AXIS),
    Diagram("kaluza_klein_monopole", "gross_perry", "through", "through the nut", ("t", "r"), (0, 16, -16, 16),
            "$x/m$", "$ct/m$", {"m": 1}, KK_AXIS, mirror=True, families=SIDEWAYS, cones=(4, 8)),
    Diagram("kaluza_klein_monopole", "hopf", "radial", "$t$ and $r$", ("t", "r"), (0, 16, -8, 8),
            "$r/m$", "$ct/m$", {"m": 1}, KK_HOPF),
    Diagram("kaluza_klein_monopole", "taub_nut", "radial", "$t$ and $\\rho$", ("t", "\\rho"), (0, 16, -8, 8),
            "$\\rho/m$", "$ct/m$", {"m": 1}, KK_HOPF),
    # The Kaluza-Klein black holes: the rays with no momentum along the circle, which are the rays of
    # the Einstein metric of four dimensions, on the plane of the time and the radius of each chart.
    Diagram("kaluza_klein_black_hole", "electric", "radial", "$t$ and $r$", ("t", "r"), (0, 6, -3, 3),
            "$r/r_s$", "$ct/r_s$", KKBH_ELECTRIC, EQUATOR, orient="ingoing", quotient="y",
            mark_gtt="the ergosurface", cone=KKBH_CONE),
    Diagram("kaluza_klein_black_hole", "eddington_finkelstein_ingoing", "finkelstein",
            "against $v - \\sqrt{q/r_s}\\,r$", ("v", "r"), (0, 6, -3, 3), "$r/r_s$",
            "$(v - \\sqrt{q/r_s}\\,r)/r_s$", KKBH_ELECTRIC, EQUATOR, to_display=KKBH_FINKELSTEIN,
            tau="v - sqrt(2)*r", quotient="w", mark_gtt="the ergosurface", cone=KKBH_CONE),
    Diagram("kaluza_klein_black_hole", "magnetic", "radial", "$t$ and $r$", ("t", "r"), (0, 6, -3, 3),
            "$r/r_s$", "$ct/r_s$", KKBH_MAGNETIC, {"theta": "0", "phi": "0", "y": "0"}, orient="ingoing"),
    Diagram("kaluza_klein_black_hole", "dyonic", "radial", "$t$ and $\\rho$", ("t", "\\rho"), (0, 6, -3, 3),
            "$\\rho/r_s$", "$ct/r_s$", KKBH_EQUAL, EQUATOR, orient="ingoing", quotient="y", cone=KKBH_CONE),
    Diagram("kaluza_klein_black_hole", "einstein", "radial", "$t$ and $r$", ("t", "r"), (0, 6, -3, 3),
            "$r/r_s$", "$ct/r_s$", KKBH_ELECTRIC, EQUATOR, orient="ingoing", areal=True),
    Diagram("kaluza_klein_black_hole", "einstein_eddington_finkelstein", "finkelstein", "against $v - r$",
            ("v", "r"), (0, 6, -3, 3), "$r/r_s$", "$(v - r)/r_s$", KKBH_ELECTRIC, EQUATOR,
            to_display=FINKELSTEIN_IN, tau="v - r", areal=True),
    Diagram("hayward", "static", "radial", "$t$ and $r$", ("t", "r"), (0, 4, -2, 2),
            "$r/m$", "$ct/m$", HAYWARD, EQUATOR, orient="ingoing", areal=True),
    Diagram("hayward", "eddington_finkelstein_ingoing", "finkelstein", "against $v - r$",
            ("v", "r"), (0, 4, -2, 2), "$r/m$", "$(v - r)/m$", HAYWARD, EQUATOR,
            to_display=FINKELSTEIN_IN, tau="v - r", areal=True),
    Diagram("hayward", "eddington_finkelstein_ingoing", "chart", "against $v$",
            ("v", "r"), (0, 4, -1, 3), "$r/m$", "$v/m$", HAYWARD, EQUATOR, tau="v - r", areal=True),
    Diagram("hayward", "eddington_finkelstein_outgoing", "finkelstein", "against $u + r$",
            ("u", "r"), (0, 4, -2, 2), "$r/m$", "$(u + r)/m$", HAYWARD, EQUATOR,
            to_display=FINKELSTEIN_OUT, tau="u + r", areal=True),
    Diagram("hayward", "eddington_finkelstein_outgoing", "chart", "against $u$",
            ("u", "r"), (0, 4, -3, 1), "$r/m$", "$u/m$", HAYWARD, EQUATOR, tau="u + r", areal=True),
    # The black hole forming and evaporating: the mass rises, holds and falls back to zero, and the
    # curve g^rr = 0, which the script marks, is closed. Outside 0 < v < 8 the plane is Minkowski's.
    Diagram("hayward", "evaporating", "history", "forming and evaporating", ("v", "r"), (0, 6, -6.5, 8.5),
            "$r/m_0$", "$(v - r)/m_0$", {"ell": "12/(7*sqrt(7))"}, EQUATOR, to_display=FINKELSTEIN_IN,
            tau="v - r", areal=True, functions={"m": HAYWARD_MASS}, input=HAYWARD_INPUT, cones=(11, 14),
            lines=(("shell", "x0", "0", "the first radiation arrives, $v = 0$"),
                   ("shell", "x0", "8", "the last of the mass is gone, $v = 8\\,m_0$"))),
    # Mass inflation: the tail falling into the charged hole, with the ray Ori's shell runs along,
    # and the charged Vaidya metric behind that shell, whose mass function ori_shell.py integrates.
    Diagram("mass_inflation", "ingoing", "tail", "the tail falling in", ("v", "r"), (0, 3, 7.5, 16.5),
            "$r/m_0$", "$(v - r)/m_0$", MI, EQUATOR, to_display=FINKELSTEIN_IN, orient="ingoing", areal=True,
            functions={"m": MI_MASS}, input=MI_INPUT, cones=(9, 14),
            lines=(("surface", "x0", "10", "the tail begins, $v = v_0$"),),
            marked=(("event", {"x0": "16", "r": repr(float(ori_shell.event_horizon(16.0)))}, 1, "the event horizon"),
                    ("shell", {"x0": "10", "r": "1"}, 1, "the outgoing ray Ori's shell runs along"))),
    Diagram("mass_inflation", "ingoing", "behind", "behind the shell", ("v", "r"), (0, 1.5, -6, 0),
            "$r/m_0$", "$(v - r)/m_0$", MI, EQUATOR, to_display=FINKELSTEIN_IN, orient="ingoing", areal=True,
            functions={"m": "ori_mass_behind(v)"}, input=MI_BEHIND_INPUT, cones=(9, 14),
            where="Min(ori_shell_behind(v) - r, -v)",
            lines=(("event", "x0", "0", "the Cauchy horizon, $v = 0$, where the mass function diverges"),
                   ("surface", "x0", repr(MI_V2_AT_V0), "the tail begins")),
            marked=(("shell", {"x0": repr(MI_V2_AT_V0), "r": "1"}, 1, "the shell"),)),
    # Hiscock's evaporating black hole, HISCOCK_MASS: the ingoing chart inside the surface of pair
    # creation, hatched beyond it, with the event horizon traced back from the last point of the
    # singularity; the outgoing chart beyond the surface; the simplest model, two shells; and the
    # flat space left at the end. Outgoing rays converge on the event horizon toward the past, so it
    # is traced back from the end, by scipy from r = 0.0007 at v = 7.9, inside the apparent horizon
    # r = 2N = 0.0014, to HISCOCK_HORIZON at v = 7.75, where r = 0.0077, and drawn from there to the
    # past only: rays part from it toward the future, and its last stretch lies on the edge r = 0.
    Diagram("hiscock", "ingoing", "history", "forming and evaporating", ("v", "r"), (0, 5, -5, 8.5),
            "$r/m_0$", "$(v - r)/m_0$", {}, EQUATOR, to_display=FINKELSTEIN_IN, tau="v - r", areal=True,
            functions={"m": HISCOCK_MASS, "R": HISCOCK_EDGE}, input=HISCOCK_INPUT, cones=(9, 14),
            singular_runs=True,
            lines=(("shell", "x0", "0", "the shell of light, $v = 0$"),),
            where="3*cos(pi*(Min(Max(v, 2), 8) - 2)/12)**2 + 1000*Heaviside(2 - v) - r",
            marked=(("event", HISCOCK_HORIZON, 1, "the event horizon", "past"),)),
    Diagram("hiscock", "ingoing", "shells", "two shells", ("v", "r"), (0, 5, -5, 8.5),
            "$r/m_0$", "$(v - r)/m_0$", {}, EQUATOR, to_display=FINKELSTEIN_IN, tau="v - r", areal=True,
            functions={"m": "Heaviside(v) - Heaviside(v - 6)", "R": "1000"}, input=HISCOCK_SHELLS_INPUT,
            cones=(9, 14), singular_runs=True,
            lines=(("shell", "x0", "0", "the shells, $v = 0$ and $v = v_1$"),
                   ("shell", "x0", "6", "the shells, $v = 0$ and $v = v_1$")),
            marked=(("event", {"x0": "6 - 1/1000", "r": "sqrt(2/1000)"}, 1, "the event horizon", "past"),
                    ("mark", {"x0": "6 + 1/1000", "r": "1/2000"}, 1, "the Cauchy horizon", "future"))),
    Diagram("hiscock", "outgoing", "history", "the radiation leaving", ("u", "r"), (0, 9, -18, 12),
            "$r/m_0$", "$(u + r)/m_0$", {}, EQUATOR, to_display=FINKELSTEIN_OUT, tau="u + r", areal=True,
            functions={"m": "hiscock_m_out(u)", "R": "hiscock_r_out(u)"}, input=HISCOCK_OUT_INPUT, cones=(9, 15),
            where="r - hiscock_r_out(u)", singular_where_claimed=True,
            lines=(("shell", "x0", "-16", "the first and the last of the radiation, $u_1$ and $u_0$"),
                   ("shell", "x0", "8", "the first and the last of the radiation, $u_1$ and $u_0$"))),
    Diagram("hiscock", "flat", "after", "after the hole", ("u", "v"), (0, 4, 8, 12),
            "$(v - u)/2m_0$", "$(u + v)/2m_0$", {}, EQUATOR, to_display=NULL_TO_TR, tau="u + v", areal=True,
            where="u - 8", input=HISCOCK_FLAT_INPUT,
            lines=(("mark", "x0", "8", "the Cauchy horizon, $u = u_0$"),)),
    *[Diagram("topological_black_hole", "static", view, label, ("t", "r"), box, "$r/L$", "$ct/L$", params,
              TBH_POINT, orient="ingoing") for view, label, params, box in TBH_CASES],
    # Against v - r the static moment t = 0 of the flat hole, v = r_*, lies below -1.5 L, and against
    # u + r above 1.5 L, so those planes are drawn lower and higher by half their height, as the BTZ hole's are.
    *[Diagram("topological_black_hole", "eddington_finkelstein_ingoing", view, label, ("v", "r"),
              box if view == "negative" else (0, 3, -3, 0), "$r/L$",
              "$(v - r)/L$", params, TBH_POINT, to_display=FINKELSTEIN_IN, orient="ingoing")
      for view, label, params, box in TBH_CASES],
    *[Diagram("topological_black_hole", "eddington_finkelstein_outgoing", view, label, ("u", "r"),
              box if view == "negative" else (0, 3, 0, 3), "$r/L$",
              "$(u + r)/L$", params, TBH_POINT, to_display=FINKELSTEIN_OUT, orient="outgoing")
      for view, label, params, box in TBH_CASES],
    Diagram("topological_black_hole", "black_string", "radial", "$t$ and $r$", ("t", "r"), (0, 3, -1.5, 1.5),
            "$r/L$", "$ct/L$", {"mu": 1, "L": 1}, {"phi": "0", "z": "0"}, orient="ingoing"),
    # In z = L^2/r the black hole is the region beyond z_h, so the future there is toward larger z.
    Diagram("topological_black_hole", "brane", "tz", "$t$ and $z$", ("t", "z"), (0, 3, -1.5, 1.5),
            "$z/L$", "$ct/L$", {"z_h": 1, "L": 1}, {"x": "0", "y": "0"}, orient="outgoing", families=SIDEWAYS),
    *[Diagram("topological_black_hole", "hyperbolic", view, label.partition(", ")[2], ("t", "r"), box, "$r/L$",
              "$ct/L$", {k: v for k, v in params.items() if k != "k"}, TBH_THETA, orient="ingoing")
      for view, label, params, box in TBH_CASES[1:]],
    Diagram("global_monopole", "static", "radial", "$t$ and $r$", ("t", "r"), (0, 6, -3, 3),
            "$r/r_s$", "$ct/r_s$", GM, EQUATOR, orient="ingoing", areal=True),
    Diagram("global_monopole", "conical", "radial", "$t$ and $r$", ("t", "r"), (0, 4, -2, 2),
            "$r$", "$ct$", GM_CONE, EQUATOR),
    Diagram("global_monopole", "eddington_finkelstein_ingoing", "finkelstein", "against $v - r$",
            ("v", "r"), (0, 6, -3, 3), "$r/r_s$", "$(v - r)/r_s$", GM, EQUATOR,
            to_display=FINKELSTEIN_IN, tau="v - r", areal=True),
    Diagram("global_monopole", "eddington_finkelstein_ingoing", "chart", "against $v$",
            ("v", "r"), (0, 6, 0, 6), "$r/r_s$", "$v/r_s$", GM, EQUATOR, tau="v - r", areal=True),
    Diagram("global_monopole", "eddington_finkelstein_outgoing", "finkelstein", "against $u + r$",
            ("u", "r"), (0, 6, -3, 3), "$r/r_s$", "$(u + r)/r_s$", GM, EQUATOR,
            to_display=FINKELSTEIN_OUT, tau="u + r", areal=True),
    Diagram("global_monopole", "eddington_finkelstein_outgoing", "chart", "against $u$",
            ("u", "r"), (0, 6, -6, 0), "$r/r_s$", "$u/r_s$", GM, EQUATOR, tau="u + r", areal=True),
    Diagram("string_black_hole", "static", "radial", "$t$ and $r$", ("t", "r"), (0, 6, -3, 3),
            "$r/r_s$", "$ct/r_s$", SBH, EQUATOR, orient="ingoing", areal=True),
    Diagram("string_black_hole", "wedge", "radial", "$t$ and $r$", ("t", "r"), (0, 6, -3, 3),
            "$r/r_s$", "$ct/r_s$", SBH, {"theta": "pi/2", "tildephi": "0"}, orient="ingoing", areal=True),
    Diagram("string_black_hole", "eddington_finkelstein_ingoing", "finkelstein", "against $v - r$",
            ("v", "r"), (0, 6, -3, 3), "$r/r_s$", "$(v - r)/r_s$", SBH, EQUATOR,
            to_display=FINKELSTEIN_IN, tau="v - r", areal=True),
    Diagram("string_black_hole", "eddington_finkelstein_ingoing", "chart", "against $v$",
            ("v", "r"), (0, 6, 0, 6), "$r/r_s$", "$v/r_s$", SBH, EQUATOR, tau="v - r", areal=True),
    Diagram("string_black_hole", "eddington_finkelstein_outgoing", "finkelstein", "against $u + r$",
            ("u", "r"), (0, 6, -3, 3), "$r/r_s$", "$(u + r)/r_s$", SBH, EQUATOR,
            to_display=FINKELSTEIN_OUT, tau="u + r", areal=True),
    Diagram("string_black_hole", "eddington_finkelstein_outgoing", "chart", "against $u$",
            ("u", "r"), (0, 6, -6, 0), "$r/r_s$", "$u/r_s$", SBH, EQUATOR, tau="u + r", areal=True),
    Diagram("mcvittie", "isotropic", "radial", "$t$ and $r$", ("t", "r"), (0, 2, 0, 8),
            "$r/r_s$", "$ct/r_s$", {"r_s": 1}, EQUATOR, areal=True, functions={"a": MCV_A},
            singular_zero="4*a*r - r_s", input=MCV_INPUT),
    Diagram("mcvittie", "areal", "radial", "$t$ and $R$", ("t", "R"), (0, 5, 0, 12),
            "$R/r_s$", "$ct/r_s$", {"r_s": 1}, EQUATOR, areal=True, functions={"H": MCV_H},
            singular_zero="R - r_s", input=MCV_INPUT),
    Diagram("dilaton_black_hole", "static", "radial", "$t$ and $r$", ("t", "r"), (0.5, 6.5, -3, 3),
            "$r/r_s$", "$ct/r_s$", DILATON, EQUATOR, orient="ingoing", areal=True),
    Diagram("dilaton_black_hole", "eddington_finkelstein_ingoing", "finkelstein", "against $v - r$",
            ("v", "r"), (0.5, 6.5, -3, 3), "$r/r_s$", "$(v - r)/r_s$", DILATON, EQUATOR,
            to_display=FINKELSTEIN_IN, tau="v - r", areal=True),
    Diagram("dilaton_black_hole", "eddington_finkelstein_ingoing", "chart", "against $v$",
            ("v", "r"), (0.5, 6.5, 0, 6), "$r/r_s$", "$v/r_s$", DILATON, EQUATOR, tau="v - r", areal=True),
    Diagram("dilaton_black_hole", "eddington_finkelstein_outgoing", "finkelstein", "against $u + r$",
            ("u", "r"), (0.5, 6.5, -3, 3), "$r/r_s$", "$(u + r)/r_s$", DILATON, EQUATOR,
            to_display=FINKELSTEIN_OUT, tau="u + r", areal=True),
    Diagram("dilaton_black_hole", "eddington_finkelstein_outgoing", "chart", "against $u$",
            ("u", "r"), (0.5, 6.5, -6, 0), "$r/r_s$", "$u/r_s$", DILATON, EQUATOR, tau="u + r", areal=True),
    Diagram("dilaton_black_hole", "string_magnetic", "radial", "$t$ and $r$", ("t", "r"), (0.5, 6.5, -3, 3),
            "$r/r_s$", "$ct/r_s$", DILATON, EQUATOR, orient="ingoing", areal=True),
    Diagram("dilaton_black_hole", "string_electric", "radial", "$t$ and $r$", ("t", "r"), (0.5, 6.5, -3, 3),
            "$r/r_s$", "$ct/r_s$", DILATON, EQUATOR, orient="ingoing", areal=True),
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
    Diagram("thin_shell_wormhole", "throat", "radial", "$t$ and $\\ell$", ("t", "\\ell"), (-3, 3, -3, 3),
            "$\\ell/r_s$", "$ct/r_s$", TSW, EQUATOR, families=SIDEWAYS, areal=True, areal_contours=(1.5, 2.0, 3.0)),
    Diagram("thin_shell_wormhole", "spherical", "radial", "$t$ and $r$", ("t", "r"), (1.25, 4.25, -1.5, 1.5),
            "$r/r_s$", "$ct/r_s$", TSW, EQUATOR, areal=True,
            lines=(("shell", "r", "5/4", "the shell at the throat, $r = a$"),)),
    # The gravastar on its plane of the time and the radius: inside the shell in the areal radius and
    # in the tortoise coordinate, where every ray runs at 45 degrees, and outside it in Schwarzschild's
    # chart from the shell out.
    Diagram("gravastar", "interior", "radial", "$t$ and $r$", ("t", "r"), (0, 1.25, -1.25, 1.25),
            "$r/r_s$", "$ct/r_s$", GRAVASTAR, EQUATOR, areal=True,
            lines=(("shell", "r", "5/4", "the shell, $r = R$"),)),
    Diagram("gravastar", "interior", "through", "through the centre", ("t", "r"), (0, 1.25, -1.25, 1.25),
            "$x/r_s$", "$ct/r_s$", GRAVASTAR, EQUATOR, mirror=True, families=SIDEWAYS, cones=(4, 8), areal=True,
            lines=(("shell", "r", "5/4", "the shell, $r = R$"),)),
    Diagram("gravastar", "interior_tortoise", "radial", "$t$ and $x$", ("t", "x"),
            (0, 2.8, -1.4, 1.4), "$x/r_s$", "$ct/r_s$", GRAVASTAR, EQUATOR,
            areal=True, areal_contours=(0.5, 1.0),
            lines=(("shell", "r", repr(GRAVASTAR_X), "the shell, $r = R$"),)),
    Diagram("gravastar", "exterior", "radial", "$t$ and $r$", ("t", "r"), (1.25, 4.25, -1.5, 1.5),
            "$r/r_s$", "$ct/r_s$", {"r_s": 1, "R": "5/4"}, EQUATOR, areal=True,
            lines=(("shell", "r", "5/4", "the shell, $r = R$"),)),
    # The wormhole time machine on the axis of the right mouth's acceleration, theta = 0, where
    # N = 1 + g l F: a stretch of the trip about the greatest acceleration, and one about the
    # greatest deceleration, half a trip's turn later. And the mouth of the short throat.
    *[Diagram("wormhole_time_machine", "wormhole", view, label, ("t", "l"),
              (-1.5, 1.5, centre - 1.5, centre + 1.5), "$l/r_0$", "$ct/r_0$", {}, {"theta": "0", "phi": "0"},
              families=SIDEWAYS, functions=WTM_FUNCTIONS, input=WTM_INPUT,
              lines=(("surface", "r", "0", "the throat, $l = 0$"),))
      for view, label, centre in (("speeding", "$t$ and $l$ on the axis, speeding up", WTM_PEAK),
                                  ("slowing", "$t$ and $l$ on the axis, slowing down", WTM_P / 2 - WTM_PEAK))],
    Diagram("wormhole_time_machine", "short_throat", "radial", "$t$ and $l$", ("t", "l"), (-3, 3, -3, 3),
            "$l/b$", "$ct/b$", {"b": 1}, EQUATOR, families=SIDEWAYS, areal=True, areal_contours=(1.5, 2.0, 3.0)),
    # Damour and Solodukhin's wormhole from the throat r = r_s out, where its own chart and the chart
    # of the rescaled time end, and through the throat in Bueno and his collaborators' rho, in the
    # isotropic radius, whose throat is r_s/4, and in Einstein and Rosen's u, u^2 = r - r_s, in units
    # of sqrt(r_s).
    Diagram("damour_solodukhin", "spherical", "radial", "$t$ and $r$", ("t", "r"), (1, 6, -2.5, 2.5),
            "$r/r_s$", "$ct/r_s$", DS, EQUATOR, areal=True),
    Diagram("damour_solodukhin", "rescaled", "radial", "$t$ and $r$", ("t", "r"), (1, 6, -2.5, 2.5),
            "$r/r_s$", "$ct/r_s$", DS, EQUATOR, areal=True),
    Diagram("damour_solodukhin", "throat", "radial", "$t$ and $\\rho$", ("t", "\\rho"), (-6, 6, -6, 6),
            "$\\rho$", "$ct/r_s$", DS, EQUATOR, families=SIDEWAYS, areal=True, areal_contours=(1.5, 2.0, 3.0)),
    Diagram("damour_solodukhin", "isotropic", "radial", "$t$ and $r$", ("t", "r"), (0, 2, -4, 4),
            "$r/r_s$", "$ct/r_s$", DS, EQUATOR, areal=True, areal_contours=(1.5, 2.0)),
    Diagram("damour_solodukhin", "einstein_rosen", "radial", "$t$ and $u$", ("t", "u"), (-2, 2, -5, 5),
            "$u/\\sqrt{r_s}$", "$ct/r_s$", DS, EQUATOR, families=SIDEWAYS, areal=True,
            areal_contours=(1.5, 2.0, 3.0)),
    # Simpson and Visser's black bounce, one view for each of its three geometries in each chart: their
    # own t and r through r = 0, Tsukamoto's areal radius on the side r > 0, from the sphere of least
    # area out, and the two Eddington-Finkelstein charts of Simpson, Martin-Moruno and Visser. Where
    # a < r_s the region between the horizons takes its future from the ingoing chart.
    *[Diagram("simpson_visser", "spherical", case, label, ("t", "r"), (-4, 4, -4, 4), "$r/r_s$", "$ct/r_s$",
              {"r_s": 1, "a": a}, EQUATOR, families=SIDEWAYS, areal=True,
              orient="ingoing" if case == "bounce" else "tau")
      for case, (label, a) in SV_CASES.items()],
    *[Diagram("simpson_visser", "areal", case, label, ("t", "\\rho"), (float(Fraction(a)), float(Fraction(a)) + 6, -3, 3),
              "$\\rho/r_s$", "$ct/r_s$", {"r_s": 1, "a": a}, EQUATOR, areal=True,
              orient="ingoing" if case == "bounce" else "tau")
      for case, (label, a) in SV_CASES.items()],
    *[Diagram("simpson_visser", "eddington_finkelstein_ingoing", case, label, ("v", "r"), (-4, 4, -4, 4),
              "$r/r_s$", "$(v - r)/r_s$", {"r_s": 1, "a": a}, EQUATOR, to_display=FINKELSTEIN_IN, tau="v - r",
              families=SIDEWAYS, areal=True)
      for case, (label, a) in SV_CASES.items()],
    *[Diagram("simpson_visser", "eddington_finkelstein_outgoing", case, label, ("u", "r"), (-4, 4, -4, 4),
              "$r/r_s$", "$(u + r)/r_s$", {"r_s": 1, "a": a}, EQUATOR, to_display=FINKELSTEIN_OUT, tau="u + r",
              families=SIDEWAYS, areal=True)
      for case, (label, a) in SV_CASES.items()],
    # Fisher, Janis, Newman and Winicour's scalar field on its plane of the time and the radial
    # coordinate in each of its four charts, at gamma = 1/2: from the singularity out in Wyman's r, in
    # Janis, Newman and Winicour's R = r - 3b/4 and in the isotropic radius, whose singularity is b/4,
    # and from spatial infinity u = 0 in toward the singularity u = infinity in Bronnikov's harmonic u.
    Diagram("fisher_jnw", "spherical", "radial", "$t$ and $r$", ("t", "r"), (1, 5, -2, 2),
            "$r/b$", "$ct/b$", FJNW, EQUATOR, areal=True),
    Diagram("fisher_jnw", "jnw", "radial", "$t$ and $R$", ("t", "R"), (0.25, 4.25, -2, 2),
            "$R/b$", "$ct/b$", FJNW, EQUATOR, areal=True),
    Diagram("fisher_jnw", "isotropic", "radial", "$t$ and $\\rho$", ("t", "\\rho"), (0.25, 4.25, -2, 2),
            "$\\rho/b$", "$ct/b$", FJNW, EQUATOR, areal=True, inside=True),
    Diagram("fisher_jnw", "harmonic", "radial", "$t$ and $u$", ("t", "u"), (0, 4, -4, 4),
            "$ku$", "$ct/k$", FJNW_HARMONIC, EQUATOR, families=("outgoing", "ingoing"), areal=True,
            areal_contours=(0.5, 1.0, 2.0, 4.0)),
    # The slowly rotating star from its surface out. On the axis the dragging term vanishes and the
    # plane of t and r holds its rays; on the equator the rays of no angular momentum are those of
    # the plane with phi divided out. The weak field chart holds where m/r is small, and there
    # its planes are flat to the eye, so none is drawn.
    Diagram("hartle_thorne", "hartle_thorne", "axis", "$t$ and $r$ on the axis", ("t", "r"), (6, 14, -4, 4),
            "$r/m$", "$ct/m$", HT_STAR, HT_AXIS, lines=(("surface", "r", "6", "the surface of the star, $r = R$"),)),
    Diagram("hartle_thorne", "hartle_thorne", "equator", "$t$ and $r$ on the equator", ("t", "r"), (6, 14, -4, 4),
            "$r/m$", "$ct/m$", HT_STAR, {"theta": "pi/2"}, quotient="phi", cone=HT_CONE,
            lines=(("surface", "r", "6", "the surface of the star, $r = R$"),)),
    Diagram("hartle_thorne", "painleve_gullstrand", "axis", "$t$ and $r$ on the axis", ("t", "r"), (6, 14, -4, 4),
            "$r/m$", "$ct/m$", HT_SIMPLE, HT_AXIS, lines=(("surface", "r", "6", "the surface of the star, $r = R$"),)),
    Diagram("hartle_thorne", "painleve_gullstrand", "equator", "$t$ and $r$ on the equator", ("t", "r"),
            (6, 14, -4, 4), "$r/m$", "$ct/m$", HT_SIMPLE, {"theta": "pi/2"}, quotient="phi", cone=HT_CONE,
            lines=(("surface", "r", "6", "the surface of the star, $r = R$"),)),
    # Witten's black hole in two dimensions, the whole plane in each of its seven charts: outside the
    # horizon in his own proper distance and in the conformal chart, across it to the singularity in
    # the Schwarzschild gauge, in the dilaton chart and in both Eddington-Finkelstein charts, and
    # the Kruskal plane with both horizons and both branches of the singularity UV = 1.
    Diagram("witten_black_hole", "witten", "radial", "$t$ and $r$", ("t", "r"), (0, 3, -1.5, 1.5),
            "$\\lambda r$", "$\\lambda ct$", {"lambda": 1}),
    Diagram("witten_black_hole", "schwarzschild_gauge", "radial", "$t$ and $x$", ("t", "x"), (-1.5, 2.5, -2, 2),
            "$\\lambda x$", "$\\lambda ct$", WITTEN, orient="ingoing"),
    Diagram("witten_black_hole", "dilaton", "radial", "$t$ and $w$", ("t", "w"), (0, 4, -2, 2),
            "$w/m$", "$\\lambda ct$", WITTEN, orient="ingoing"),
    Diagram("witten_black_hole", "conformal", "radial", "$t$ and $\\sigma$", ("t", "\\sigma"), (-2.5, 2.5, -2.5, 2.5),
            "$\\lambda\\sigma$", "$\\lambda ct$", WITTEN),
    Diagram("witten_black_hole", "kruskal", "plane", "$U$ and $V$", ("U", "V"), (-2, 2, -2, 2),
            "$(V - U)/2$", "$(U + V)/2$", WITTEN, to_display=NULL_TO_TR, tau="U + V", families=SIDEWAYS,
            where="1 - U*V", singular_zero="1 - U*V",
            marked=(("shell", {"x0": "0", "r": "0"}, "both", "the horizons, $UV = 0$", "past"),)),
    Diagram("witten_black_hole", "eddington_finkelstein_ingoing", "finkelstein", "against $v - x$",
            ("v", "x"), (-1.5, 2.5, -2, 2), "$\\lambda x$", "$\\lambda(v - x)$", WITTEN,
            to_display=FINKELSTEIN_IN, tau="v - x"),
    Diagram("witten_black_hole", "eddington_finkelstein_outgoing", "finkelstein", "against $u + x$",
            ("u", "x"), (-1.5, 2.5, -2, 2), "$\\lambda x$", "$\\lambda(u + x)$", WITTEN,
            to_display=FINKELSTEIN_OUT, tau="u + x"),
    # Roberts's collapsing scalar field, one view for each of its three outcomes in each chart. The
    # double null chart is drawn against (v - u)/2 and (u + v)/2, Roberts's v and r against r and
    # (1 + p) v - r, which are sqrt(1 + p) times the diagonal chart's rho and ct, the areal chart
    # against R and (1 + p/2) v - R, and the scaling chart against x and tau + x, since neither tau nor
    # x is a time. The field arrives on v = 0, the lower edge of each wedge, with flat space below it.
    # Each box holds a stretch of every moment of Roberts's t that the embedding diagram marks.
    *[Diagram("roberts", "double_null", case, label, ("u", "v"), (0, 6, -5, 3), "$(v - u)/2\\ell$", "$(u + v)/2\\ell$",
              {"p": p}, EQUATOR, to_display=NULL_TO_TR, tau="u + v", areal=True, null_radius=True,
              singular_where_claimed=True, singular_runs=case == "disperses",
              singular_zero=None if case == "disperses" else f"(1 - {p})*v - u")
      for case, (label, p) in ROBERTS_CASES.items()],
    *[Diagram("roberts", "advanced", case, label, ("v", "r"), (0, 8, -4, 4), "$r/\\ell$", "$((1 + p)v - r)/\\ell$",
              {"p": p}, EQUATOR, to_display=((0, 1), (1 + float(Fraction(p)), -1)), tau=f"(1 + {p})*v - r", areal=True,
              singular_where_claimed=True, singular_zero=None if case == "disperses" else f"r - {p}*v")
      for case, (label, p) in ROBERTS_CASES.items()],
    *[Diagram("roberts", "areal", case, label, ("v", "R"), (0, 6, -4, 7), "$R/\\ell$", "$((1 + p/2)v - R)/\\ell$",
              {"p": p}, EQUATOR, to_display=((0, 1), (1 + float(Fraction(p)) / 2, -1)),
              tau=f"(1 + ({p})/2)*v - sqrt(R**2 + ({p})**2*v**2/4)", areal=True, singular_where_claimed=True,
              singular_runs=True, step=0.001, inside=case == "collapses")
      for case, (label, p) in ROBERTS_CASES.items()],
    *[Diagram("roberts", "diagonal", case, label, ("t", "\\rho"), (0, 6, -3, 3), "$\\rho/\\ell$", "$ct/\\ell$",
              {"p": p}, EQUATOR, areal=True, singular_where_claimed=True,
              singular_zero=None if case == "disperses" else f"rho - {p}*t")
      for case, (label, p) in ROBERTS_CASES.items()],
    *[Diagram("roberts", "scaling", case, label, ("\\tau", "x"), (0, 0.6 if case == "collapses" else 2, -2, 4),
              "$x$", "$\\tau + x$", {"p": p, "ell": 1}, EQUATOR, to_display=((0, 1), (1, 1)), tau="tau + x",
              families=("outgoing", "ingoing"), areal=True, null_radius=True, singular_where_claimed=True,
              singular_zero="x - atanh(1/2)" if case == "collapses" else None)
      for case, (label, p) in ROBERTS_CASES.items()],
    # Einstein and Rosen's bridges: the neutral one through the bridge in their u, in units of
    # sqrt(r_s), on one sheet in Schwarzschild's r, and through the bridge in the isotropic radius,
    # where it is r_s/4; the charged one with a mass on one sheet from r_+ = 3r_s/2 out, and the
    # charged one with no mass through the bridge in their u, in units of r_q.
    Diagram("einstein_rosen_bridge", "bridge", "radial", "$t$ and $u$", ("t", "u"), (-2, 2, -5, 5),
            "$u/\\sqrt{r_s}$", "$ct/r_s$", {"r_s": 1}, EQUATOR, families=SIDEWAYS, cones=(8, 7), areal=True,
            areal_contours=(1.5, 2.0, 3.0)),
    Diagram("einstein_rosen_bridge", "spherical", "radial", "$t$ and $r$", ("t", "r"), (0.5, 6, -2.75, 2.75),
            "$r/r_s$", "$ct/r_s$", {"r_s": 1}, EQUATOR, areal=True, where="r - 1"),
    Diagram("einstein_rosen_bridge", "isotropic", "radial", "$t$ and $r$", ("t", "r"), (0, 2, -4, 4),
            "$r/r_s$", "$ct/r_s$", {"r_s": 1}, EQUATOR, areal=True, areal_contours=(1.5, 2.0)),
    Diagram("einstein_rosen_bridge", "charged_spherical", "radial", "$t$ and $r$", ("t", "r"), (1, 6.5, -2.75, 2.75),
            "$r/r_s$", "$ct/r_s$", ERB_CHARGED, EQUATOR, areal=True, where="2*r - 3"),
    Diagram("einstein_rosen_bridge", "charged_bridge", "radial", "$t$ and $u$", ("t", "u"), (-3, 3, -3, 3),
            "$u/r_q$", "$ct/r_q$", {"r_q": 1}, EQUATOR, families=SIDEWAYS, cones=(8, 7), areal=True,
            areal_contours=(1.5, 2.0, 3.0)),
    Diagram("morris_thorne", "spherical", "radial", "$t$ and $r$", ("t", "r"), (0, 4, -2, 2),
            "$r/b_0$", "$ct/b_0$", {"b_0": 1}, EQUATOR, areal=True,
            functions={"Phi": "0", "b": "b_0**2/r"},
            input="$\\Phi = 0$ and $b = b_0^2/r$, the member of the family that is the "
                  "Ellis-Bronnikov wormhole."),
    # Teo's example: on the axis the dragging term vanishes and the plane of t and r holds its
    # rays, drawn at his a = 1/4; on the equator N = 1, the rays of no angular momentum are
    # those of the plane with phi divided out, and a = 1 puts the ergosurface at sqrt(2) b_0.
    Diagram("teo_wormhole", "spherical", "axis", "$t$ and $r$ on the axis", ("t", "r"), (1, 5, -2, 2),
            "$r/b_0$", "$ct/b_0$", {"b_0": 1, "a": "1/4"}, {"theta": "0", "phi": "0"}),
    Diagram("teo_wormhole", "spherical", "equator", "$t$ and $r$ on the equator", ("t", "r"), (1, 5, -2, 2),
            "$r/b_0$", "$ct/b_0$", {"b_0": 1, "a": 1}, {"theta": "pi/2"}, quotient="phi",
            mark_gtt="the ergosurface", cone=TEO_CONE),
    Diagram("teo_wormhole", "proper_radial", "axis", "$t$ and $l$ on the axis", ("t", "l"), (-3, 3, -3, 3),
            "$l/b_0$", "$ct/b_0$", {"b_0": 1, "a": "1/4"}, {"theta": "0", "phi": "0"}, families=SIDEWAYS,
            functions={"r": TEO_RADIUS}, input=TEO_RADIUS_INPUT, lines=(("throat", "r", "0", None),)),
    Diagram("teo_wormhole", "proper_radial", "equator", "$t$ and $l$ on the equator", ("t", "l"), (-3, 3, -3, 3),
            "$l/b_0$", "$ct/b_0$", {"b_0": 1, "a": 1}, {"theta": "pi/2"}, families=SIDEWAYS, quotient="phi",
            functions={"r": TEO_RADIUS}, input=TEO_RADIUS_INPUT, lines=(("throat", "r", "0", None),),
            mark_gtt="the ergosurface", cone=TEO_CONE),
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
    # Misner space at Li and Gott's psi_0 = 4 pi, a boost of rapidity 2 pi, where the Milne chart's
    # chi and the Rindler chart's eta each run once round 2 pi. Misner's own plane is drawn against
    # psi/2, that rapidity, so that it stands as wide as the Milne chart's.
    Diagram("misner", "misner", "plane", "$T$ and $\\psi$", ("T", "\\psi"), (0, 2 * math.pi, -2, 2),
            "$\\psi/2$", "$T$", {"psi_0": "4*pi"}, {"y": "0", "z": "0"}, to_display=((0, 0.5), (1, 0)),
            orient="outgoing", families=SIDEWAYS, mark_g00=True, periodic=("\\psi",)),
    Diagram("misner", "milne", "plane", "$t$ and $\\chi$", ("t", "\\chi"), (0, 2 * math.pi, -2.5, 0.5),
            "$\\chi$", "$ct$", {"psi_0": "4*pi"}, {"y": "0", "z": "0"}, families=SIDEWAYS,
            periodic=("\\chi",)),
    Diagram("misner", "rindler", "plane", "$\\eta$ and $\\xi$", ("\\eta", "\\xi"), (0, 2 * math.pi, 0, 2 * math.pi),
            "$\\xi$", "$\\eta$", {"psi_0": "4*pi"}, {"y": "0", "z": "0"}, tau="eta", families=SIDEWAYS,
            periodic=("\\eta",)),
    # Gott's two strings away from the strings, Grant's generalised Misner space, for strings of
    # half deficit angle pi/3 at v = 4c/5 and d = l/2: cosh(a/4) = gamma sin(alpha) = 5/(2 sqrt 3)
    # and b = 4 gamma v d sin(alpha)/(c sinh(a/4)) = 8 l/sqrt 13. The Rindler plane marks the first
    # two polarised hypersurfaces, xi_n = n b/(2 sinh(n a/2)).
    Diagram("gott_time_machine", "grant_rindler", "plane", "$\\eta$ and $\\xi$", ("\\eta", "\\xi"),
            (0, 2, 0, GOTT_A), "$\\xi/\\ell$", "$\\eta$", GOTT, {"Y": "0", "z": "0"}, tau="eta",
            families=SIDEWAYS, periodic=("\\eta",),
            lines=(("surface", "r", f"({GOTT['b']})/(2*sinh(({GOTT['a']})/2))",
                    "the first polarised hypersurface, $\\xi = b/(2\\sinh(a/2))$"),
                   ("surface", "r", f"({GOTT['b']})/sinh({GOTT['a']})",
                    "the second, $\\xi = b/\\sinh a$"))),
    Diagram("gott_time_machine", "grant_milne", "plane", "$\\tau$ and $\\chi$", ("\\tau", "\\chi"),
            (0, GOTT_A, -2.5, 0), "$\\chi$", "$c\\tau/\\ell$", GOTT, {"Y": "0", "z": "0"}, tau="tau",
            families=SIDEWAYS, periodic=("\\chi",)),
    # Ori's vacuum core with his example f = a(x^2 - y^2)/2 at a = 1/16, e = 1/8 and L = 2 pi: the
    # cylinder of the time and z through the central circle x = y = 0, which is totally geodesic and
    # Misner's, and the cylinder at x = 4 l, y = 0, where f = l^2/2 and e x^2 = 2 l^2, in each of
    # his two charts; and the Brinkmann plane of u and v through the central circle, flat, whose
    # line u = -2 a circuit of z carries onto u = -2 exp(-pi).
    Diagram("ori_time_machine", "vacuum_core", "centre", "the central circle", ("T", "z"),
            (0, 2 * math.pi, -2.5, 2), "$z$", "$T/\\ell^2$", {"L": "2*pi"}, {"x": "0", "y": "0"},
            orient="outgoing", families=SIDEWAYS, periodic=("z",), functions={"f": ORI_F}, input=ORI_INPUT,
            lines=(("surface", "x0", "0", "the closed null geodesic $N$, $T = 0$"),)),
    Diagram("ori_time_machine", "vacuum_core", "off_centre", "$x = 4\\,\\ell$", ("T", "z"),
            (0, 2 * math.pi, -4, 2), "$z$", "$T/\\ell^2$", {"L": "2*pi"}, {"x": "4", "y": "0"},
            orient="outgoing", families=SIDEWAYS, periodic=("z",), functions={"f": ORI_F}, input=ORI_INPUT,
            lines=(("surface", "x0", "1/2", "the null circle, $T = f = \\ell^2/2$"),)),
    Diagram("ori_time_machine", "foliation", "centre", "the central circle", ("t", "z"),
            (0, 2 * math.pi, -2.5, 2), "$z$", "$t/\\ell^2$", ORI, {"x": "0", "y": "0"},
            orient="outgoing", families=SIDEWAYS, periodic=("z",),
            lines=(("surface", "x0", "0", "the closed null geodesic $N$, $t = 0$"),)),
    Diagram("ori_time_machine", "foliation", "off_centre", "$x = 4\\,\\ell$", ("t", "z"),
            (0, 2 * math.pi, -2.5, 3.5), "$z$", "$t/\\ell^2$", ORI, {"x": "4", "y": "0"},
            orient="outgoing", families=SIDEWAYS, periodic=("z",),
            lines=(("surface", "x0", "2", "the null circle, $t = ex^2 = 2\\,\\ell^2$"),)),
    Diagram("ori_time_machine", "brinkmann", "plane", "$u$ and $v$", ("u", "v"), (-1, 2.5, -2.5, 1),
            "$(v/\\ell^2 - u)/2$", "$(u + v/\\ell^2)/2$", {"a": "1/16", "L": "2*pi"}, {"x": "0", "y": "0"},
            to_display=NULL_TO_TR, tau="u + v", families=SIDEWAYS,
            lines=(("surface", "r", "0", "the closed null geodesic $N$, $v = 0$"),
                   ("shell", "x0", "-2", "$u = -2$ and $u = -2e^{-L/2}$, one line of the core"),
                   ("shell", "x0", "-2*exp(-pi)", "$u = -2$ and $u = -2e^{-L/2}$, one line of the core"))),
    # The spinning string's cylinders of one time and one angle at one radius, unrolled, at b = 0.9,
    # the deficit the cosmic string is drawn at, and a = 0.9, so that the null circle r_c = a/b is
    # the unit of length: inside it and outside it in the proper radius, rescaled radius and helical
    # charts, and at R = a in the circumference radius chart, which ends on the null circle. The
    # angle is scaled by the circle's own radius, b r, rho, R or r, so that the cones of a string
    # without spin would stand at 45 degrees.
    *[Diagram("spinning_string", system, view, f"$t$ and $\\phi$ at ${name}$", ("t", "\\phi"),
              tuple(s * math.pi * width for s in (-1, 1, -1, 1)), across, "$ct/r_c$", SPINNING, {held: at, "z": "0"},
              to_display=((0, width), (1, 0)), orient="vector", families=SIDEWAYS, cones=(5, 5), periodic=("\\phi",))
      for system, view, name, across, held, at, width in (
          ("proper_radius", "inside", "r = r_c/2", "$br\\phi/r_c$", "r", "1/2", 0.45),
          ("proper_radius", "outside", "r = 3r_c/2", "$br\\phi/r_c$", "r", "3/2", 1.35),
          ("rescaled_radius", "inside", "\\rho = a/2", "$\\rho\\phi/r_c$", "rho", "9/20", 0.45),
          ("rescaled_radius", "outside", "\\rho = 3a/2", "$\\rho\\phi/r_c$", "rho", "27/20", 1.35),
          ("circumference_radius", "outside", "R = a", "$R\\phi/r_c$", "R", "9/10", 0.9))],
    *[Diagram("spinning_string", "helical", view, f"$\\tau$ and $\\tilde\\phi$ at ${name}$", ("\\tau", "\\tilde\\phi"),
              (0, 2 * math.pi * 0.9 * r, 0, 2 * math.pi * 0.9 * r), "$r\\tilde\\phi/r_c$", "$c\\tau/r_c$", SPINNING,
              {"r": at, "z": "0"}, to_display=((0, r), (1, 0)), tau="tau", families=SIDEWAYS, cones=(5, 5),
              periodic=("\\tilde\\phi",))
      for view, name, at, r in (("inside", "r = r_c/2", "1/2", 0.5), ("outside", "r = 3r_c/2", "3/2", 1.5))],
    # The travelling wave on a string, on the plane of u and v drawn against z and ct, at fixed
    # places across the string: on the two sides of it in the null conical chart, at r = ell, where
    # the isotropic x is ell/4; at that x in the isotropic chart; and at X = -ell/4 and 3 ell/4 in the
    # moving string chart, which the string draws away from and comes to within ell/4 of. The time
    # function is 3u + v, whose gradient is timelike wherever g_uu > -3, since g_uu falls to -2 here
    # and u + v stops being a time below -1.
    *[Diagram("string_wave", "null_conical", view, label, ("u", "v"), (-2, 2, -2, 2), "$z/\\ell$", "$ct/\\ell$",
              {"b": "1/2"}, {"r": "1", "phi": phi}, to_display=NULL_TO_TR, tau="3*u + v", families=SIDEWAYS,
              functions=STRING_WAVE_PROFILE, solves=(("v", "u"),), lines=STRING_WAVE_CREST,
              input=STRING_WAVE_PROFILE_INPUT)
      for view, label, phi in (("toward", "$\\phi = 0$", "0"), ("away", "$\\phi = \\pi$", "pi"))],
    Diagram("string_wave", "isotropic", "beside", "$x = \\ell/4$", ("u", "v"), (-2, 2, -2, 2), "$z/\\ell$", "$ct/\\ell$",
            STRING_WAVE, {"x": "1/4", "y": "0"}, to_display=NULL_TO_TR, tau="3*u + v", families=SIDEWAYS,
            functions=STRING_WAVE_PULSE, lines=STRING_WAVE_CREST, input=STRING_WAVE_INPUT),
    *[Diagram("string_wave", "moving_string", view, label, ("u", "V"), (-2, 2, -2, 2), "$z/\\ell$", "$ct/\\ell$",
              STRING_WAVE, {"X": X, "Y": "0"}, to_display=NULL_TO_TR, tau="3*u + V", families=SIDEWAYS,
              functions={"rho": "sqrt((X - A)**2 + (Y - B)**2)", **STRING_WAVE_PULSE}, lines=STRING_WAVE_CREST,
              input=STRING_WAVE_INPUT)
      for view, label, X in (("behind", "$X = -\\ell/4$", "-1/4"), ("ahead", "$X = 3\\ell/4$", "3/4"))],
    # Point particles in three dimensions, on the plane of t and the radius at alpha = 3/4: the
    # conical chart and the chart with the wedge removed, where every ray runs at 45 degrees, and
    # the isotropic chart, where ct -+ (ell/alpha)(rho/ell)^alpha are constant. Two particles at rest
    # are drawn on the plane of t and x through both, y = 0, a plane of symmetry, between them and
    # beyond one of them, each view ending on a particle, where a ray ends too.
    Diagram("point_particle_2plus1", "conical", "radial", "$t$ and $r$", ("t", "r"), (0, 4, -2, 2),
            "$r/\\ell$", "$ct/\\ell$", PARTICLE, {"phi": "0"}),
    Diagram("point_particle_2plus1", "wedge", "radial", "$t$ and $r$", ("t", "r"), (0, 4, -2, 2),
            "$r/\\ell$", "$ct/\\ell$", PARTICLE, {"theta": "0"}),
    Diagram("point_particle_2plus1", "circumference", "radial", "$t$ and $R$", ("t", "R"), (0, 4, -2, 2),
            "$R/\\ell$", "$ct/\\ell$", PARTICLE, {"phi": "0"}),
    # Gott and Alpert's planet of the same mass, cos(chi_0) = 3/4, from its centre to its edge.
    Diagram("point_particle_2plus1", "planet", "radial", "$t$ and $\\chi$", ("t", "\\chi"),
            (0, math.acos(0.75), -math.acos(0.75) / 2, math.acos(0.75) / 2),
            "$\\chi$", "$ct/a$", {"a": 1, "chi_0": "acos(3/4)"}, {"phi": "0"}),
    Diagram("point_particle_2plus1", "isotropic", "radial", "$t$ and $\\rho$", ("t", "\\rho"), (0, 4, -2, 2),
            "$\\rho/\\ell$", "$ct/\\ell$", {**PARTICLE, "ell": 1}, {"phi": "0"}),
    Diagram("point_particle_2plus1", "two_bodies", "between", "between the particles", ("t", "x"), (-1, 1, -1, 1),
            "$x/d$", "$ct/d$", TWO_PARTICLES, {"y": "0"}, families=SIDEWAYS),
    Diagram("point_particle_2plus1", "two_bodies", "beyond", "beyond them", ("t", "x"), (1, 5, -2, 2),
            "$x/d$", "$ct/d$", TWO_PARTICLES, {"y": "0"}, families=SIDEWAYS),
    Diagram("de_sitter", "static_spherical", "radial", "$t$ and $r$", ("t", "r"), (0, 2, -1, 1),
            "$r\\sqrt{\\Lambda/3}$", "$ct\\sqrt{\\Lambda/3}$", {"Lambda": 3}, EQUATOR,
            orient="outgoing", cones=(8, 7), areal=True),
    Diagram("de_sitter", "static_spherical", "through", "through the observer", ("t", "r"),
            (0, 2, -2, 2), "$x\\sqrt{\\Lambda/3}$", "$ct\\sqrt{\\Lambda/3}$", {"Lambda": 3},
            EQUATOR, mirror=True, orient="outgoing", families=SIDEWAYS, cones=(4, 8), areal=True),
    Diagram("de_sitter", "flat_slicing", "tx", "$t$ and $x$", ("t", "x"), (-2, 2, -1, 3),
            "$x\\;[c/H]$", "$ct\\;[c/H]$", {"H": 1}, {"y": "0", "z": "0"}, families=SIDEWAYS),
    Diagram("einstein_static", "hyperspherical", "radial", "$t$ and $\\chi$", ("t", "\\chi"),
            (0, math.pi, 0, math.pi), "$\\chi$", "$ct/R$", {"R": 1}, EQUATOR),
    Diagram("einstein_static", "hyperspherical", "through", "through the pole", ("t", "\\chi"),
            (0, math.pi, 0, 2 * math.pi), "$\\chi$", "$ct/R$", {"R": 1}, EQUATOR, mirror=True,
            families=SIDEWAYS, cones=(4, 8)),
    Diagram("einstein_static", "static_areal", "radial", "$t$ and $r$", ("t", "r"), (0, 1.2, -0.6, 0.6),
            "$r/R$", "$ct/R$", {"R": 1}, EQUATOR, areal=True),
    Diagram("milne", "comoving_hyperbolic", "through", "through the observer", ("t", "\\chi"),
            (0, 3, 0, 3), "$\\chi$", "$ct$", {}, EQUATOR, mirror=True, families=SIDEWAYS, cones=(4, 8)),
    Diagram("milne", "comoving_spherical", "radial", "$t$ and $r$", ("t", "r"), (0, 4, 0, 3),
            "$r$", "$ct$", {}, EQUATOR, areal=True),
    Diagram("milne", "logarithmic_time", "radial", "$\\tau$ and $\\chi$", ("\\tau", "\\chi"),
            (0, 3, -2, 1.5), "$\\chi$", "$c\\tau\\;[ct_0]$", {"t_0": 1}, EQUATOR, tau="tau"),
    Diagram("milne", "inertial", "through", "through the observer", ("T", "R"), (0, 4, 0, 5),
            "$x$", "$cT$", {}, EQUATOR, mirror=True, tau="T", families=SIDEWAYS, cones=(4, 8), areal=True),
    # The domain wall at k = 1: Rindler's plane on either side of the wall at z = 0, out to the
    # horizons z = +-1/k, and the inertial chart of one side, the inside of the hyperbola the wall is.
    Diagram("domain_wall", "planar", "tz", "$t$ and $z$", ("t", "z"), (-1, 1, -1.5, 1.5),
            "$kz$", "$kct$", {"k": 1}, {"x": "0", "y": "0"}, families=SIDEWAYS, cones=(6, 7),
            lines=(("shell", "r", "0", "the wall, $z = 0$"),)),
    Diagram("domain_wall", "inertial", "through", "through the centre", ("T", "R"), (0, 2, -2, 2),
            "$kx$", "$kcT$", {"k": 1}, EQUATOR, mirror=True, tau="T", families=SIDEWAYS, cones=(4, 8),
            areal=True, marked=(("shell", {"x0": "0", "r": "0"}, "both", "the rays $R = c|T|$ through the centre at "
                                "$T = 0$, the horizons of the planar and global charts"),)),
    # Randall and Sundrum's wall at k = 1, on the plane of t and the fifth coordinate in each chart, every
    # x_i held fixed: the proper distance chart out to 2.5/k on either side, where the cones close as
    # e^{-k|y|}; the conformally flat chart, whose rays run at 45 degrees; the Poincare chart of one side
    # from the wall at z = 1/k; and the two walls at k r_c = 1/2, whose edges phi = +-pi are one wall.
    Diagram("randall_sundrum", "proper_distance", "ty", "$t$ and $y$", ("t", "y"), (-2.5, 2.5, -6, 6),
            "$ky$", "$kct$", {"k": 1}, RS_WALL, families=SIDEWAYS, cones=(6, 7),
            lines=(("shell", "r", "0", "the wall, $y = 0$"),)),
    Diagram("randall_sundrum", "conformal", "tw", "$t$ and $w$", ("t", "w"), (-4, 4, -4, 4),
            "$kw$", "$kct$", {"k": 1}, RS_WALL, families=SIDEWAYS, cones=(6, 7),
            lines=(("shell", "r", "0", "the wall, $w = 0$"),)),
    Diagram("randall_sundrum", "poincare", "tz", "$t$ and $z$", ("t", "z"), (1, 5, -2, 2),
            "$kz$", "$kct$", {"k": 1}, RS_WALL, families=SIDEWAYS, cones=(5, 7)),
    Diagram("randall_sundrum", "two_walls", "tphi", "$t$ and $\\phi$", ("t", "\\phi"), (-math.pi, math.pi, -4, 4),
            "$\\phi$", "$kct$", {"k": 1, "r_c": "1/2"}, RS_WALL, families=SIDEWAYS, cones=(6, 7),
            periodic=("\\phi",), lines=(("shell", "r", "0", "the wall of positive tension, $\\phi = 0$"),)),
    # Coleman and De Luccia's bubble with a thin wall, in units of the curvature radius l of the
    # vacuum that has one. Decay into flat space: de Sitter space outside, Minkowski space inside,
    # CDL's rho_0 = l, so the wall is at rho_bar = 4l/5 and rho' drops from 1 to 3/5 across it.
    # Decay of flat space: Minkowski space outside, anti-de Sitter space inside, rho_0 = l, the
    # wall at rho_bar = 4l/3, where rho' drops from 5/3 to 1. coleman_de_luccia.md derives both.
    Diagram("coleman_de_luccia", "wall", "into_flat", "decay into flat space", ("\\psi", "\\xi"), (0, 3, -2, 2),
            "$\\xi/\\ell$", "$\\psi$", {}, EQUATOR, tau="psi", functions={"rho": CDL_INTO_FLAT},
            lines=(("shell", "r", "4/5", "the wall, $\\rho = \\bar\\rho$"),), input=CDL_INPUT["into_flat"]),
    Diagram("coleman_de_luccia", "wall", "out_of_flat", "decay of flat space", ("\\psi", "\\xi"), (0, 3, -2, 2),
            "$\\xi/\\ell$", "$\\psi$", {}, EQUATOR, tau="psi", functions={"rho": CDL_OUT_OF_FLAT},
            lines=(("shell", "r", "log(3)", "the wall, $\\rho = \\bar\\rho$"),), input=CDL_INPUT["out_of_flat"]),
    Diagram("coleman_de_luccia", "open", "zero", "zero vacuum energy", ("\\tau", "\\chi"), (0, 3, 0, 3),
            "$\\chi$", "$c\\tau/\\ell$", {}, EQUATOR, tau="tau", functions={"a": "tau"},
            input="A new vacuum of zero energy, $a = c\\tau$, with $\\ell$ any length."),
    Diagram("coleman_de_luccia", "open", "positive", "positive vacuum energy", ("\\tau", "\\chi"), (0, 3, 0, 3),
            "$\\chi$", "$c\\tau/\\ell$", {}, EQUATOR, tau="tau", functions={"a": "sinh(tau)"},
            input="A new vacuum of positive energy, $a = \\ell\\sinh(c\\tau/\\ell)$ with $\\ell = \\sqrt{3/\\Lambda}$."),
    Diagram("coleman_de_luccia", "open", "negative", "negative vacuum energy", ("\\tau", "\\chi"),
            (0, 3, 0, math.pi), "$\\chi$", "$c\\tau/\\ell$", {}, EQUATOR, tau="tau", functions={"a": "sin(tau)"},
            input="A new vacuum of negative energy, $a = \\ell\\sin(c\\tau/\\ell)$ with $\\ell = \\sqrt{-3/\\Lambda}$."),
    Diagram("coleman_de_luccia", "open_conformal", "negative", "negative vacuum energy", ("\\eta", "\\chi"),
            (0, 3, -3, 3), "$\\chi$", "$\\eta$", {}, EQUATOR, tau="eta", functions={"a": "1/cosh(eta)"},
            input="A new vacuum of negative energy, $a = \\ell/\\cosh\\eta$ with $\\ell = \\sqrt{-3/\\Lambda}$, "
                  "for the Kretschmann scalar."),
    Diagram("coleman_de_luccia", "static_inside", "into_flat", "decay into flat space", ("t", "r"), (0, 4, 0, 4),
            "$r/\\ell$", "$ct/\\ell$", {"Lambda_T": 0}, EQUATOR, areal=True,
            functions={"r_w": "sqrt(Rational(16, 25) + t**2)"}, input=CDL_INPUT["into_flat"]),
    Diagram("coleman_de_luccia", "static_inside", "out_of_flat", "decay of flat space", ("t", "r"),
            (0, 5, 0, math.pi), "$r/\\ell$", "$ct/\\ell$", {"Lambda_T": -3}, EQUATOR, areal=True,
            functions={"r_w": CDL_ADS_WALL}, input=CDL_INPUT["out_of_flat"]),
    Diagram("coleman_de_luccia", "static_outside", "into_flat", "decay into flat space", ("t", "r"), (0, 2, 0, 3),
            "$r/\\ell$", "$ct/\\ell$", {"Lambda_F": 3}, EQUATOR, orient="outgoing", areal=True,
            functions={"r_w": "sqrt(Rational(16, 25) + Rational(9, 25)*tanh(t)**2)"}, input=CDL_INPUT["into_flat"]),
    Diagram("coleman_de_luccia", "static_outside", "out_of_flat", "decay of flat space", ("t", "r"), (0, 5, 0, 4),
            "$r/\\ell$", "$ct/\\ell$", {"Lambda_F": 0}, EQUATOR, areal=True,
            functions={"r_w": "sqrt(Rational(16, 9) + t**2)"}, input=CDL_INPUT["out_of_flat"]),
    Diagram("anti_de_sitter", "static_global", "radial", "$t$ and $r$", ("t", "r"), (0, 4, -2, 2),
            "$r/L$", "$ct/L$", {"L": 1}, EQUATOR, areal=True),
    Diagram("anti_de_sitter", "static_global", "through", "through the centre", ("t", "r"),
            (0, 4, -4, 4), "$x/L$", "$ct/L$", {"L": 1}, EQUATOR, mirror=True,
            families=SIDEWAYS, cones=(4, 8), areal=True),
    Diagram("anti_de_sitter", "poincare", "tx", "$t$ and $x$", ("t", "x"), (-2, 2, -2, 2),
            "$x/L$", "$ct/L$", {"L": 1}, {"y": "0", "z": "1"}, families=SIDEWAYS),
    Diagram("rn_metric", "spherical", "radial", "$t$ and $r$", ("t", "r"), (0, 2, -1, 1),
            "$r/r_s$", "$t/r_s$", {"r_s": 1, "r_q": "12/25"}, EQUATOR, orient="ingoing",
            areal=True),
    # Majumdar and Papapetrou's two holes, each of mass parameter m = 1, at z = +-2 on the axis,
    # the pair the embedding diagram's midplane is cut between; and one hole alone, U = 1 + m/r.
    Diagram("majumdar_papapetrou", "cartesian", "tz", "the axis", ("t", "z"), (-5, 5, -5, 5),
            "$z/m$", "$ct/m$", {}, {"x": "0", "y": "0"}, families=SIDEWAYS,
            functions={"U": MP_TWO}, input=MP_TWO_INPUT,
            lines=(("grr", "r", "2", None), ("grr", "r", "-2", None))),
    Diagram("majumdar_papapetrou", "cartesian", "tx", "the midplane", ("t", "x"), (-5, 5, -5, 5),
            "$x/m$", "$ct/m$", {}, {"y": "0", "z": "0"}, families=SIDEWAYS,
            functions={"U": MP_TWO}, input=MP_TWO_INPUT),
    Diagram("majumdar_papapetrou", "cylindrical", "radial", "$t$ and $\\rho$", ("t", "\\rho"), (0, 5, -2.5, 2.5),
            "$\\rho/m$", "$ct/m$", {}, {"phi": "0", "z": "0"},
            functions={"U": MP_TWO_CYLINDRICAL}, input=MP_TWO_INPUT),
    Diagram("majumdar_papapetrou", "isotropic", "radial", "$t$ and $r$", ("t", "r"), (0, 4, -2, 2),
            "$r/m$", "$ct/m$", {"m": 1}, EQUATOR, areal=True),
    Diagram("israel_wilson_perjes", "cylindrical", "midplane", "the midplane", ("t", "\\rho"), (1, 7, -3, 3),
            "$\\rho/m$", "$ct/m$", {}, {"z": "0"}, functions=IWP_TWO, input=IWP_TWO_INPUT, quotient="phi",
            cone=BTZ_CONE),
    Diagram("israel_wilson_perjes", "spheroidal", "axis", "the axis, through the ring", ("t", "r"), (-4, 4, -4, 4),
            "$r/m$", "$ct/m$", {"m": 1, "a": 1}, {"theta": "0", "phi": "0"}, families=SIDEWAYS),
    Diagram("israel_wilson_perjes", "spheroidal", "principal", "principal null rays, $t$ and $r$", ("t", "r"),
            (-0.5, 4.5, -2.5, 2.5), "$r/m$", "$ct/m$", {"m": 1, "a": 1}, {"theta": "pi/2"},
            principal=True, leaves=("phi",), cone=PRINCIPAL_CONE),
    Diagram("israel_wilson_perjes", "spherical", "radial", "$t$ and $r$", ("t", "r"), (1, 6, -2.5, 2.5),
            "$r/m$", "$ct/m$", {"m": 1, "l": "1/2"}, EQUATOR),
    *_kt_rows(),
    Diagram("btz", "stationary", "static", "$J = 0$", ("t", "r"), (0, 3, -1.5, 1.5),
            "$r/\\ell$", "$ct/\\ell$", {"ell": 1, "M": 1, "J": 0}, {"phi": "0"}, orient="ingoing"),
    Diagram("btz", "stationary", "rotating", "$J = 4\\ell/5$", ("t", "r"), (0, 2, -1, 1),
            "$r/\\ell$", "$ct/\\ell$", {"ell": 1, "M": 1, "J": "4/5"}, orient="ingoing", quotient="phi",
            mark_gtt="the ergosurface", cone=BTZ_CONE),
    Diagram("btz", "eddington_finkelstein_ingoing", "static", "$J = 0$", ("v", "r"), (0, 3, -3, 0),
            "$r/\\ell$", "$(v - r)/\\ell$", {"ell": 1, "M": 1, "J": 0}, {"tildephi": "0"},
            to_display=FINKELSTEIN_IN, orient="ingoing"),
    Diagram("btz", "eddington_finkelstein_ingoing", "rotating", "$J = 4\\ell/5$", ("v", "r"), (0, 2, -1, 1),
            "$r/\\ell$", "$(v - r)/\\ell$", {"ell": 1, "M": 1, "J": "4/5"}, to_display=FINKELSTEIN_IN,
            orient="ingoing", quotient="tildephi", mark_gtt="the ergosurface", cone=BTZ_CONE),
    Diagram("btz", "eddington_finkelstein_outgoing", "static", "$J = 0$", ("u", "r"), (0, 3, 0, 3),
            "$r/\\ell$", "$(u + r)/\\ell$", {"ell": 1, "M": 1, "J": 0}, {"tildephi": "0"},
            to_display=FINKELSTEIN_OUT, orient="outgoing"),
    Diagram("btz", "eddington_finkelstein_outgoing", "rotating", "$J = 4\\ell/5$", ("u", "r"), (0, 2, -1, 1),
            "$r/\\ell$", "$(u + r)/\\ell$", {"ell": 1, "M": 1, "J": "4/5"}, to_display=FINKELSTEIN_OUT,
            orient="outgoing", quotient="tildephi", mark_gtt="the ergosurface", cone=BTZ_CONE),
    Diagram("taub_nut", "spherical", "radial", "$t$ and $r$", ("t", "r"), (0, 6, -3, 3),
            "$r/m$", "$ct/m$", {"m": 1, "l": "1/2"}, EQUATOR, orient="ingoing"),
    Diagram("bertotti_robinson", "static", "radial", "$t$ and $r$", ("t", "r"), (0, 3, -1.5, 1.5),
            "$r/b$", "$ct/b$", {"b": 1}, EQUATOR),
    Diagram("bertotti_robinson", "poincare", "tx", "$t$ and $x$", ("t", "x"), (0, 4, -2, 2),
            "$x/b$", "$ct/b$", {"b": 1}, EQUATOR, families=SIDEWAYS),
    # The throat of extreme Kerr at r_0 = 1, on the equator with phi divided out, where the plane's
    # metric is half that of the two dimensional anti-de Sitter space of each chart: the rays of no
    # angular momentum are the throat's principal null rays, the same on every cone of theta, which
    # print_charts.py checks. The near-NHEK chart is drawn at k = r_0/2, its horizon at r = r_0.
    Diagram("near_horizon_extreme_kerr", "poincare", "equator", "$t$ and $r$ on the equator", ("t", "r"),
            (0, 3, -1.5, 1.5), "$r/r_0$", "$ct/r_0$", {"r_0": 1}, {"theta": "pi/2"}, quotient="phi", cone=NHEK_CONE),
    Diagram("near_horizon_extreme_kerr", "inverse_radius", "equator", "$t$ and $x$ on the equator", ("t", "x"),
            (0, 4, -2, 2), "$x/r_0$", "$ct/r_0$", {"r_0": 1}, {"theta": "pi/2"}, families=SIDEWAYS, quotient="phi",
            cone=NHEK_CONE),
    Diagram("near_horizon_extreme_kerr", "global", "equator", "$\\tau$ and $y$ on the equator", ("\\tau", "y"),
            (-3, 3, -2, 2), "$y$", "$\\tau$", {"r_0": 1}, {"theta": "pi/2"}, tau="tau", families=SIDEWAYS,
            quotient="phi", mark_gtt="where $\\partial_\\tau$ turns spacelike", cone=NHEK_CONE),
    Diagram("near_horizon_extreme_kerr", "near_nhek", "equator", "$t$ and $r$ on the equator", ("t", "r"),
            (1, 4, -1.5, 1.5), "$r/r_0$", "$ct/r_0$", {"r_0": 1, "k": "1/2"}, {"theta": "pi/2"}, quotient="phi",
            cone=NHEK_CONE),
    # The Nariai universe's static patch reaches from one horizon, r = -1/sqrt(Lambda), to the other,
    # and its global chart runs round the whole circle of chi, the static patch the half 0 < chi < pi.
    Diagram("nariai", "static", "patch", "$t$ and $r$", ("t", "r"), (-1, 1, -2, 2),
            "$r\\sqrt{\\Lambda}$", "$ct\\sqrt{\\Lambda}$", {"Lambda": 1}, EQUATOR, families=SIDEWAYS),
    Diagram("nariai", "global", "circle", "$t$ and $\\chi$", ("t", "\\chi"), (0, 2 * math.pi, -2, 2),
            "$\\chi$", "$ct\\sqrt{\\Lambda}$", {"Lambda": 1}, EQUATOR, families=SIDEWAYS, periodic=("\\chi",)),
    Diagram("interior_schwarzschild", "spherical", "radial", "$t$ and $r$", ("t", "r"),
            (0, 1.5, -0.75, 0.75), "$r/r_s$", "$t/r_s$", {"r_s": 1, "R": "3/2"}, EQUATOR,
            areal=True),
    Diagram("interior_schwarzschild", "spherical", "through", "through the centre", ("t", "r"),
            (0, 1.5, -1.5, 1.5), "$x/r_s$", "$t/r_s$", {"r_s": 1, "R": "3/2"}, EQUATOR,
            mirror=True, families=SIDEWAYS, cones=(4, 8), areal=True),
    # Tolman's solution VII, the star whose density reaches zero at its surface, at R = 2 r_s, which
    # is beta = 1/4: in Lattimer and Prakash's form, and in Tolman's constants for the same star,
    # R^2 -> R^2/(5 beta), A^4 = 4R^4/(3 beta), B^2 = 1 - 5 beta/3 and the C that makes his phase
    # pi/2 less theirs, which print_charts.tolman_vii_check holds the two charts to.
    Diagram("tolman_vii", "spherical", "radial", "$t$ and $r$", ("t", "r"),
            (0, 2, -1, 1), "$r/r_s$", "$ct/r_s$", TOLMAN_VII_STAR, EQUATOR, areal=True),
    Diagram("tolman_vii", "spherical", "through", "through the centre", ("t", "r"),
            (0, 2, -2, 2), "$x/r_s$", "$ct/r_s$", TOLMAN_VII_STAR, EQUATOR,
            mirror=True, families=SIDEWAYS, cones=(4, 8), areal=True),
    Diagram("tolman_vii", "tolman", "radial", "$t$ and $r$", ("t", "r"),
            (0, 2, -1, 1), "$r/r_s$", "$ct/r_s$", TOLMAN_VII_CONSTANTS, EQUATOR, areal=True),
    Diagram("tolman_vii", "tolman", "through", "through the centre", ("t", "r"),
            (0, 2, -2, 2), "$x/r_s$", "$ct/r_s$", TOLMAN_VII_CONSTANTS, EQUATOR,
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
    # Kerr-de Sitter and Kerr-anti-de Sitter. Carter's chart: the axis, where the plane of t and r
    # holds its rays, and the principal null congruence on the equator, in t and r and from above.
    # The chart turned to the frame that does not rotate at infinity has the same t and r, so only
    # its view from above is new. The Kerr charts and the Kerr-Schild chart run through the
    # horizons on the axis.
    Diagram("kerr_de_sitter", "boyer_lindquist", "axis", "$\\Lambda > 0$, the axis", ("t", "r"), (0, 4, -2, 2),
            "$r/r_s$", "$ct/r_s$", KDS, KDS_AXIS, orient="split", split=KDS_SPLIT),
    Diagram("kerr_de_sitter", "boyer_lindquist", "principal", "$\\Lambda > 0$, principal null rays", ("t", "r"),
            (0, 4, -2, 2), "$r/r_s$", "$ct/r_s$", KDS, {"theta": "pi/2"}, orient="split", split=KDS_SPLIT,
            principal=True, leaves=("phi",), mark_gtt="the ergosurfaces", cone=PRINCIPAL_CONE),
    Diagram("kerr_de_sitter", "boyer_lindquist", "above", "$\\Lambda > 0$, from above", ("\\phi", "r"),
            (-3.4, 3.4, -3.4, 3.4), "$r\\cos\\phi/r_s$", "$r\\sin\\phi/r_s$", KDS, {"theta": "pi/2"},
            to_display=POLAR, cones=(0, 0), principal=True, leaves=("t",), ring=12, inside=True,
            mark_gtt="the ergosurfaces"),
    Diagram("kerr_de_sitter", "boyer_lindquist", "axis_ads", "$\\Lambda < 0$, the axis", ("t", "r"),
            (0, 3, -1.5, 1.5), "$r/\\ell$", "$ct/\\ell$", KADS, KDS_AXIS, orient="ingoing"),
    Diagram("kerr_de_sitter", "boyer_lindquist", "principal_ads", "$\\Lambda < 0$, principal null rays", ("t", "r"),
            (0.2, 3.2, -1.5, 1.5), "$r/\\ell$", "$ct/\\ell$", KADS, {"theta": "pi/2"}, orient="ingoing",
            principal=True, leaves=("phi",), mark_gtt="the ergosurface", cone=PRINCIPAL_CONE),
    Diagram("kerr_de_sitter", "nonrotating", "above", "$\\Lambda > 0$, from above", ("\\Phi", "r"),
            (-3.4, 3.4, -3.4, 3.4), "$r\\cos\\Phi/r_s$", "$r\\sin\\Phi/r_s$", KDS, {"theta": "pi/2"},
            to_display=POLAR, cones=(0, 0), principal=True, leaves=("t",), ring=12, inside=True),
    Diagram("kerr_de_sitter", "nonrotating", "above_ads", "$\\Lambda < 0$, from above", ("\\Phi", "r"),
            (-3, 3, -3, 3), "$r\\cos\\Phi/\\ell$", "$r\\sin\\Phi/\\ell$", KADS, {"theta": "pi/2"},
            to_display=POLAR, cones=(0, 0), principal=True, leaves=("t",), ring=12, inside=True),
    Diagram("kerr_de_sitter", "kerr_ingoing", "axis", "$\\Lambda > 0$, the axis", ("v", "r"), (0, 4, -2, 2),
            "$r/r_s$", "$(v - r)/r_s$", KDS, KDS_KERR_AXIS, to_display=FINKELSTEIN_IN, tau="v - r"),
    Diagram("kerr_de_sitter", "kerr_ingoing", "axis_ads", "$\\Lambda < 0$, the axis", ("v", "r"), (0, 3, -3, 0),
            "$r/\\ell$", "$(v - r)/\\ell$", KADS, KDS_KERR_AXIS, to_display=FINKELSTEIN_IN, orient="ingoing"),
    Diagram("kerr_de_sitter", "kerr_outgoing", "axis", "$\\Lambda > 0$, the axis", ("u", "r"), (0, 4, -2, 2),
            "$r/r_s$", "$(u + r)/r_s$", KDS, KDS_KERR_AXIS, to_display=FINKELSTEIN_OUT, tau="u + r"),
    Diagram("kerr_de_sitter", "kerr_outgoing", "axis_ads", "$\\Lambda < 0$, the axis", ("u", "r"), (0, 3, 0, 3),
            "$r/\\ell$", "$(u + r)/\\ell$", KADS, KDS_KERR_AXIS, to_display=FINKELSTEIN_OUT, orient="outgoing"),
    Diagram("kerr_de_sitter", "kerr_schild", "axis", "$\\Lambda > 0$, the axis", ("\\tau", "r"), (0, 4, -2, 2),
            "$r/r_s$", "$c\\tau/r_s$", KDS, {"theta": "0", "psi": "0"}, tau="tau"),
    Diagram("kerr_de_sitter", "kerr_schild", "axis_ads", "$\\Lambda < 0$, the axis", ("\\tau", "r"), (0, 3, -1.5, 1.5),
            "$r/\\ell$", "$c\\tau/\\ell$", KADS, {"theta": "0", "psi": "0"}, tau="tau"),
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
    # The Kantowski-Sachs universes: the dust universe symmetric in time, in its comoving chart with
    # a and b solved from rest and in the dust chart's own time eta, and the vacuum member, the
    # inside of Schwarzschild's horizon, whose future lies toward smaller T.
    Diagram("kantowski_sachs", "comoving", "tr", "$t$ and $r$", ("t", "r"), (-2, 2, 0, math.pi),
            "$r\\;[b_0]$", "$ct\\;[b_0]$", {}, EQUATOR, families=SIDEWAYS, dust=KS_DUST,
            reference="$b = b_0$", input=KS_INPUT),
    Diagram("kantowski_sachs", "dust", "etar", "$\\eta$ and $r$", ("\\eta", "r"), (-2, 2, -math.pi / 2, math.pi / 2),
            "$r/b_0$", "$\\eta$", {"b_0": 1, "kappa": 0}, EQUATOR, tau="eta", families=SIDEWAYS),
    Diagram("kantowski_sachs", "schwarzschild_interior", "Tr", "$T$ and $r$", ("T", "r"), (-2, 2, 0, 1),
            "$r/r_s$", "$T/r_s$", {"r_s": 1}, EQUATOR, tau="-T", families=SIDEWAYS),
    Diagram("godel", "cartesian", "tx", "$t$ and $x$", ("t", "x"), (-2, 2, -2, 2), "$x$", "$t$",
            {"omega": 1}, {"y": "0", "z": "0"}, orient="vector", families=SIDEWAYS),
    Diagram("alcubierre", "cartesian", "tx", "$t$ and $x$ on the axis", ("t", "x"), (-3, 3, -2, 2),
            "$x/R$", "$ct/R$", {}, {"y": "0", "z": "0"}, families=SIDEWAYS, cones=(8, 7),
            functions={"v_s": "2", "f": _alcubierre_profile()},
            input="$v_s = 2$, and Alcubierre's own profile, "
                  "$f = [\\tanh\\sigma(r_s + R) - \\tanh\\sigma(r_s - R)]/(2\\tanh\\sigma R)$ "
                  "with $R = 1$ and $\\sigma = 4$."),
    Diagram("van_den_broeck", "cartesian", "tx", "$t$ and $x$ on the axis", ("t", "x"), (-3, 3, -2, 2),
            "$x/R$", "$ct/R$", {}, {"y": "0", "z": "0"}, families=SIDEWAYS, cones=(8, 7), kretschmann=False,
            functions={"v_s": "2", "f": _vdb_shape(_VDB_RS), "B": _vdb_factor(_VDB_RS)},
            input="$v_s = 2$, a wall in which " + VDB_WALL + ", and " + VDB_POCKET + "."),
    Diagram("van_den_broeck", "pocket", "radial", "$t$ and $r$", ("t", "r"), (0, 1.25, -1.25, 1.25),
            "$r/R$", "$ct/R$", {"R": 1}, EQUATOR, functions={"B": _vdb_factor("r")},
            lines=(("shell", "r", "1/2", "the outer end of the neck, where $B$ returns to $1$"),),
            input=VDB_POCKET + "."),
    Diagram("van_den_broeck", "proper_radial", "radial", "$t$ and $l$", ("t", "l"), (-4, 1, -2.5, 2.5),
            "$l/R$", "$ct/R$", {"l_0": 4}, EQUATOR, families=SIDEWAYS, functions={"r": _VDB_R},
            lines=(("surface", "r", "-2", "the widest sphere of the pocket, $r = 7R/4$"),
                   ("shell", "r", "0", "the narrowest sphere of the neck, $r = R/4$")),
            input=VDB_POCKET + "."),
    Diagram("natario", "cartesian_flow", "tx", "$t$ and $x$ on the axis", ("t", "x"), (-3, 3, -2, 2),
            "$x/R$", "$ct/R$", {}, {"y": "0", "z": "0"}, families=SIDEWAYS, cones=(8, 7),
            functions=_natario_field(),
            input="$v_s = 2$, $n = f/2$ with Alcubierre's shape function $f$ for a bubble of radius $R$, and the zero expansion field "
                  "$X = v_s[(2n + \\rho n')\\,e_x - n'\\,x_r\\,(x_r, y, z)/\\rho]$, "
                  "$x_r = x - v_s t$, $\\rho^2 = x_r^2 + y^2 + z^2$, whose divergence vanishes."),
    Diagram("krasnikov", "cylindrical", "tx", "$t$ and $x$ on the axis", ("t", "x"), (-1, 5, -1, 5),
            "$x$", "$ct$", {}, {"r": "0", "phi": "0"}, orient="outgoing", families=SIDEWAYS,
            functions={"k": _KRASNIKOV_TUBE}, kretschmann=False, mark_g00=True,
            input="A tube along $x$ from $0$ to $D = 4$, built by a ship that left $x = 0$ at $t = 0$ "
                  "at the speed of light: $k = 1 - (2 - \\delta)\\,S(\\tfrac{\\rho_0^2 - "
                  "r^2}{2\\rho_0})\\,S(t - x)\\,S(x)\\,S(D - x)$ with $\\delta = 0.2$, $\\rho_0 = 1$, "
                  "and $S$ a step of width $0.15$ built from $\\tanh$."),
    # Where both waves have passed; the chart's domain ends on the two wave fronts, and above it
    # the metric stops being finite on the curvature singularity u^2 + v^2 = 1.
    Diagram("khan_penrose", "double_null", "plane", "$u$ and $v$", ("u", "v"), (-1, 1, 0, 1.5),
            "$v - u$", "$u + v$", {"L": 1}, {"x": "0", "y": "0"}, to_display=NULL_TO_SUM, tau="u + v",
            families=SIDEWAYS, crunch=True),
    Diagram("khan_penrose", "cosmological", "plane", "$\\tau$ and $\\sigma$", ("\\tau", "\\sigma"),
            (-math.pi / 2, math.pi / 2, 0, math.pi / 2), "$\\sigma$", "$\\tau$", {"L": 1}, {"x": "0", "y": "0"},
            tau="tau", families=SIDEWAYS, singular_where_claimed=True),
    # Bell and Szekeres's colliding electromagnetic waves at a = b = 1, where both have passed: each
    # chart's domain ends on the two wave fronts and on the Killing-Cauchy horizon au + bv = pi/2, and
    # the regular, global and Kruskal-Szekeres charts are drawn on through the horizon.
    Diagram("bell_szekeres", "double_null", "plane", "$u$ and $v$", ("u", "v"), (-1.6, 1.6, 0, 1.6),
            "$a(v - u)$", "$a(u + v)$", {"a": 1, "b": 1}, {"x": "0", "y": "0"}, to_display=NULL_TO_SUM,
            tau="u + v", families=SIDEWAYS, kretschmann=False),
    Diagram("bell_szekeres", "time_space", "plane", "$\\xi$ and $\\eta$", ("\\xi", "\\eta"),
            (-math.pi / 2, math.pi / 2, 0, math.pi / 2), "$\\eta$", "$\\xi$", {"a": 1, "b": 1},
            {"x": "0", "y": "0"}, tau="xi", families=SIDEWAYS, kretschmann=False),
    Diagram("bell_szekeres", "regular", "plane", "$T$ and $Z$", ("T", "Z"), (-1.5, 1.5, -1.5, 0.5),
            "$Z$", "$T$", {"a": 1, "b": 1}, {"X": "99/100", "Y": "0"}, tau="T/sqrt(1 + Z**2)", families=SIDEWAYS,
            kretschmann=False, marked=(("shell", {"x0": "0", "r": "0"}, "both",
                                        "the Killing-Cauchy horizon, $T = -|Z|$", "past"),)),
    Diagram("bell_szekeres", "global", "plane", "$\\chi$ and $\\rho$", ("\\chi", "\\rho"),
            (-2.5, 2.5, -math.pi / 2, 0.5), "$\\rho$", "$\\chi$", {"a": 1, "b": 1},
            {"theta": "pi/2", "phi": "0"}, tau="chi", families=SIDEWAYS, kretschmann=False,
            marked=(("shell", {"x0": "0", "r": "0"}, "both",
                     "the Killing-Cauchy horizon, $\\cos\\chi\\cosh\\rho = 1$", "past"),)),
    Diagram("bell_szekeres", "kruskal_szekeres", "plane", "$U$ and $V$", ("U", "V"), (-2, 2, -2.5, 0.5),
            "$a(V - U)$", "$a(U + V)$", {"a": 1, "b": 1}, {"eta": "0", "x": "0"}, to_display=NULL_TO_SUM,
            tau="U + V", families=SIDEWAYS, kretschmann=False,
            marked=(("shell", {"x0": "0", "r": "0"}, "both", "the Killing-Cauchy horizon, $UV = 0$", "past"),)),
    Diagram("bell_szekeres", "bertotti_robinson", "plane", "$t$ and $r$", ("t", "r"), (0, 3, 0, 3),
            "$r$", "$t$", {"a": 1, "b": 1}, {"theta": "pi/2", "phi": "0"}, families=SIDEWAYS, kretschmann=False),
    Diagram("pp_wave", "exact_plane_wave", "tz", "$t$ and $z$ on the axis", ("u", "v"), (-2, 2, -2, 2),
            "$z$", "$ct$", {}, {"x": "0", "y": "0"}, to_display=UV_TO_TZ, tau="u + 2*v",
            families=SIDEWAYS, functions={"A": "exp(-u**2)", "B": "0"},
            input="The amplitudes $A(u)$ and $B(u)$ do not enter on the axis, so the diagram is the "
                  "same for every profile."),
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
    # Som and Raychaudhuri's cylinders of t and phi about one world line of the dust, in units of
    # r_c = c/Omega, where the block of t and phi is van Stockum's with R = r_c.
    Diagram("som_raychaudhuri", "cylindrical", "inside", "$t$ and $\\phi$ at $r = r_c/2$", ("t", "\\phi"),
            (-math.pi / 2, math.pi / 2, -math.pi / 2, math.pi / 2), "$r\\phi/r_c$", "$ct/r_c$", {"Omega": 1},
            {"r": "1/2", "z": "0"}, to_display=((0, 0.5), (1, 0)), orient="vector", families=SIDEWAYS,
            cones=(5, 5), periodic=("\\phi",)),
    Diagram("som_raychaudhuri", "cylindrical", "beyond", "$t$ and $\\phi$ at $r = 3r_c/2$", ("t", "\\phi"),
            (-3 * math.pi / 2, 3 * math.pi / 2, -3 * math.pi / 2, 3 * math.pi / 2), "$r\\phi/r_c$", "$ct/r_c$",
            {"Omega": 1}, {"r": "3/2", "z": "0"}, to_display=((0, 1.5), (1, 0)), orient="vector",
            families=SIDEWAYS, cones=(5, 5), periodic=("\\phi",)),
    Diagram("som_raychaudhuri", "cartesian", "tx", "$t$ and $x$", ("t", "x"), (-2, 2, -2, 2), "$x/r_c$", "$ct/r_c$",
            {"Omega": 1}, {"y": "0", "z": "0"}, orient="vector", families=SIDEWAYS),
    Diagram("tov", "spherical", "radial", "$t$ and $r$", ("t", "r"), (0, 16, -8, 8),
            "$r\\;[GM_\\odot/c^2]$", "$ct\\;[GM_\\odot/c^2]$", {}, EQUATOR, areal=True, star=POLYTROPE,
            input=POLYTROPE_INPUT),
    Diagram("tov", "spherical", "through", "through the centre", ("t", "r"), (0, 16, -16, 16),
            "$x\\;[GM_\\odot/c^2]$", "$ct\\;[GM_\\odot/c^2]$", {}, EQUATOR, mirror=True, families=SIDEWAYS,
            cones=(4, 8), areal=True, star=POLYTROPE, input=POLYTROPE_INPUT),
    # Einstein's cluster, in units of its Schwarzschild radius: the areal chart for the declared cluster
    # with no surface that the embedding diagram draws, Einstein's own at V = 1/2, where m = r/6
    # and R = 3 r_s, in the areal radius and in his isotropic one, where sigma = 5 - 2 sqrt 6 and
    # rho_0 = R/(1 + sigma)^2, and Florides's of the same mass and radius in the areal radius and
    # on its three sphere, a = sqrt(R^3/r_s) and sin(chi_0) = R/a.
    Diagram("einstein_cluster", "areal", "radial", "$t$ and $r$", ("t", "r"), (0, 8, -4, 4),
            "$r/r_s$", "$ct/r_s$", {}, EQUATOR, areal=True, functions={"Phi": EC_CORE}, input=EC_CORE_INPUT),
    Diagram("einstein_cluster", "areal", "through", "through the centre", ("t", "r"), (0, 8, -8, 8),
            "$x/r_s$", "$ct/r_s$", {}, EQUATOR, mirror=True, families=SIDEWAYS, cones=(4, 8), areal=True,
            functions={"Phi": EC_CORE}, input=EC_CORE_INPUT),
    Diagram("einstein_cluster", "constant_speed", "radial", "$t$ and $r$", ("t", "r"), (0, 3, -1.5, 1.5),
            "$r/r_s$", "$ct/r_s$", {"V": "1/2", "R": 3}, EQUATOR, areal=True),
    Diagram("einstein_cluster", "isotropic", "radial", "$t$ and $\\rho$", ("t", "\\rho"),
            (0, EC_RHO0, -EC_RHO0 / 2, EC_RHO0 / 2), "$\\rho/r_s$", "$ct/r_s$",
            {"sigma": "5 - 2*sqrt(6)", "rho_0": "3/(6 - 2*sqrt(6))**2"}, EQUATOR, areal=True),
    Diagram("einstein_cluster", "uniform", "radial", "$t$ and $r$", ("t", "r"), (0, 3, -1.5, 1.5),
            "$r/r_s$", "$ct/r_s$", {"r_s": 1, "R": 3}, EQUATOR, areal=True),
    Diagram("einstein_cluster", "uniform", "through", "through the centre", ("t", "r"), (0, 3, -3, 3),
            "$x/r_s$", "$ct/r_s$", {"r_s": 1, "R": 3}, EQUATOR, mirror=True, families=SIDEWAYS, cones=(4, 8),
            areal=True),
    Diagram("einstein_cluster", "hyperspherical", "radial", "$t$ and $\\chi$", ("t", "\\chi"),
            (0, EC_CHI0, -1.6, 1.6), "$\\chi$", "$ct/r_s$", {"a": "sqrt(27)", "chi_0": "asin(1/sqrt(3))"},
            EQUATOR, areal=True),
    # The pulse of Weber, Wheeler, and Bonnor, whose metric on the plane of t and rho is conformally
    # flat, so its rays are at 45 degrees for every pulse; it comes in along rho = -ct and goes out
    # along rho = ct.
    Diagram("einstein_rosen_waves", "cylindrical", "radial", "$t$ and $\\rho$", ("t", "\\rho"), (0, 10, -1, 9),
            "$\\rho/a$", "$ct/a$", {}, {"phi": "0", "z": "0"}, functions=_er_pulse("t", "rho"), solves=ER_SOLVES,
            marked=(("shell", {"x0": "0", "r": "0"}, "both", "the rays $\\rho = |ct|$ through that event, just inside "
                     "the crest of the pulse"),),
            points=(("mark", ("0", "0"), "the event on the axis where the pulse is greatest, $\\psi = 2C/a$"),),
            input=ER_INPUT),
    Diagram("einstein_rosen_waves", "null", "radial", "$t$ and $\\rho$", ("u", "v"), (0, 10, -1, 9),
            "$(v - u)/2a$", "$(u + v)/2a$", {}, {"phi": "0", "z": "0"}, to_display=NULL_TO_TR, tau="u + v",
            functions=_er_pulse("(u + v)/2", "(v - u)/2"),
            marked=(("shell", {"x0": "0", "r": "0"}, "both", "the rays $v = 0$ and $u = 0$ through that event, just "
                     "inside the crest of the pulse"),),
            points=(("mark", ("0", "0"), "the event on the axis where the pulse is greatest, $\\psi = 2C/a$"),),
            input=ER_INPUT),
    # Gowdy's torus: the plane of the time and theta is conformally flat in the areal time, so its
    # rays are at 45 degrees for every wave, and in tau = -ln t they are theta -+ e^{-tau} = const.
    Diagram("gowdy", "areal", "plane", "$t$ and $\\theta$", ("t", "\\theta"), (0, 2 * math.pi, 0, 2 * math.pi),
            "$\\theta$", "$t$", {"L": 1}, {"sigma": "0", "delta": "0"}, families=SIDEWAYS, periodic=("\\theta",),
            functions=_gowdy_wave("t"), solves=_gowdy_solves("t"), input=GOWDY_INPUT),
    Diagram("gowdy", "logarithmic", "plane", "$\\tau$ and $\\theta$", ("\\tau", "\\theta"), (0, 2 * math.pi, -4, 2),
            "$\\theta$", "$-\\tau$", {"L": 1}, {"sigma": "0", "delta": "0"}, to_display=((0, 1), (-1, 0)), tau="-tau",
            families=SIDEWAYS, periodic=("\\theta",), functions=_gowdy_wave("exp(-tau)"),
            solves=_gowdy_solves("\\tau"), input=GOWDY_INPUT),
    # The sphere chart, for the inside of Schwarzschild's horizon, r = L(1 - cos t) with r_s = 2L:
    # the plane is conformally flat for every wave, so the rays are t -+ theta = const.
    Diagram("gowdy", "sphere", "plane", "$t$ and $\\theta$", ("t", "\\theta"), (0, math.pi, 0, math.pi),
            "$\\theta$", "$t$", {"L": 1}, {"sigma": "0", "delta": "0"}, families=SIDEWAYS,
            functions=GOWDY_HOLE, input=GOWDY_HOLE_INPUT),
    # Melvin's plane of t and rho, conformally flat, and Ernst's equator, conformal to Schwarzschild's
    # plane of t and r, at B r_s = 1/2, as the embedding diagram draws it.
    Diagram("melvin", "cylindrical", "radial", "$t$ and $\\rho$", ("t", "\\rho"), (0, 4, -2, 2),
            "$B\\rho$", "$Bct$", {"B": 1}, {"phi": "0", "z": "0"},
            lines=(("surface", "r", "2", "the Melvin radius $\\rho = 2/B$, where the circles about the axis are widest"),)),
    # Senovilla's plane of t and rho, conformally flat, so its rays are at 45 degrees through the
    # bounce t = 0, with the radius cosh(3a rho) = 2 of the one circular light ray and the event
    # where the density is greatest.
    Diagram("senovilla", "cylindrical", "radial", "$t$ and $\\rho$", ("t", "\\rho"), (0, 2, -1.5, 1.5),
            "$a\\rho$", "$act$", {"a": 1}, {"phi": "0", "z": "0"},
            lines=(("surface", "r", "acosh(2)/3", "the radius $\\cosh(3a\\rho) = 2$ of the one circular path of light"),),
            points=(("mark", ("0", "0"), "the event on the axis where the density is greatest, at the bounce $t = 0$"),)),
    Diagram("melvin", "ernst", "radial", "$t$ and $r$", ("t", "r"), (0, 6, -3, 3),
            "$r/r_s$", "$ct/r_s$", {"r_s": 1, "B": "1/2"}, EQUATOR, orient="ingoing",
            lines=(("surface", "r", "4", "$r = 2/B$, the widest circle of the equator"),)),
    # Levi-Civita's plane of t and its radius at sigma = 1/4, where the Kasner exponents are
    # (2/3, 2/3, -1/3), in Weyl's coordinates and in the Kasner form, whose r is the proper distance.
    Diagram("levi_civita", "weyl", "radial", "$t$ and $\\rho$", ("t", "\\rho"), (0, 4, -2, 2),
            "$\\rho$", "$ct$", LC_WEYL, {"phi": "0", "z": "0"}),
    Diagram("levi_civita", "kasner", "radial", "$t$ and $r$", ("t", "r"), (0, 4, -2, 2),
            "$r$", "$ct$", LC_KASNER, {"phi": "0", "z": "0"}),
    # The Curzon-Chazy particle on its two totally geodesic planes, the axis and the plane z = 0, in
    # each chart: on the axis the cones close toward R = 0, and in the plane they open.
    Diagram("curzon_chazy", "weyl", "axis", "$t$ and $z$ on the axis", ("t", "z"), (0, 4, -2, 2),
            "$z/m$", "$ct/m$", {"m": 1}, {"rho": "0", "phi": "0"}),
    Diagram("curzon_chazy", "weyl", "equator", "$t$ and $\\rho$ in the plane $z = 0$", ("t", "\\rho"), (0, 4, -2, 2),
            "$\\rho/m$", "$ct/m$", {"m": 1}, {"phi": "0", "z": "0"},
            lines=(("surface", "r", "1", "$\\rho = m$, the narrowest circle about the axis"),)),
    Diagram("curzon_chazy", "spherical", "axis", "$t$ and $r$ on the axis", ("t", "r"), (0, 4, -2, 2),
            "$r/m$", "$ct/m$", {"m": 1}, {"theta": "0", "phi": "0"}),
    Diagram("curzon_chazy", "spherical", "equator", "$t$ and $r$ in the plane $\\theta = \\pi/2$", ("t", "r"), (0, 4, -2, 2),
            "$r/m$", "$ct/m$", {"m": 1}, {**EQUATOR, "phi": "0"},
            lines=(("surface", "r", "1", "$r = m$, the narrowest circle about the axis"),)),
    # Zipoy and Voorhees's metric on its two totally geodesic planes, the axis and the equatorial
    # plane, in each chart, for the oblate q = 1 and the prolate q = -1/2. The prolate equator's
    # curvature diverges only as the 3/2 power, so its rows declare the edge singular.
    *[Diagram("zipoy_voorhees", "spherical", f"{plane}_{shape}", f"{name}, $q = {q}$", ("t", "r"), (2, 6, -2, 2),
              "$r/m$", "$ct/m$", {"m": 1, "q": q}, fixed,
              lines=((("surface", "r", "3", "$r = 3\\,m$, the narrowest circle about the axis"),)
                     if (plane, shape) == ("equator", "oblate") else ()),
              singular_zero="r - 2" if (plane, shape) == ("equator", "prolate") else None)
      for shape, q in (("oblate", "1"), ("prolate", "-1/2"))
      for plane, name, fixed in (("axis", "The axis", {"theta": "0", "phi": "0"}),
                                 ("equator", "The equatorial plane", {**EQUATOR}))],
    *[Diagram("zipoy_voorhees", "prolate_spheroidal", f"{plane}_{shape}", f"{name}, $\\delta = {d}$", ("t", "x"),
              (1, 5, -2, 2), "$x$", "$ct/m$", {"m": 1, "delta": d}, fixed,
              lines=((("surface", "r", "2", "$x = 2$, the narrowest circle about the axis"),)
                     if (plane, shape) == ("equator", "oblate") else ()),
              singular_zero="x - 1" if (plane, shape) == ("equator", "prolate") else None)
      for shape, d in (("oblate", "2"), ("prolate", "1/2"))
      for plane, name, fixed in (("axis", "The axis", {"y": "1", "phi": "0"}),
                                 ("equator", "The equatorial plane", {"y": "0", "phi": "0"}))],
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
                  "zero at $r_b$, with $R(r, 0) = r$ and $2GM/c^2 = r_b/2$ for the mass $M$ of the cloud: each shell falls as "
                  "$R^{3/2} = r^{3/2} - \\tfrac{3}{2}\\sqrt{2GM(r)/c^2}\\,ct$, checked to solve this "
                  "spacetime's own $G^r{}_r = 0$."),
    # Szekeres's quasispherical dust with an axis of symmetry, drawn on the two halves of that axis,
    # which its light rays never leave: the Tolman-Bondi cloud above with the centres of its shells
    # moved along the axis, S'/S = 2r(1 - r^2) inside the cloud and S constant outside it, where
    # the spacetime is Schwarzschild's. The density stays positive and no shells cross, since
    # S'/S < M'/(3M) = 5(1 - r^2)/(r(5 - 3r^2)) for r < 1.
    *[Diagram("szekeres", "axisymmetric", view, label, ("t", "r"), (0, 1.5, -0.5, 1.0), "$r/r_b$", "$ct/r_b$", {},
              {"theta": theta, "phi": "0"}, functions={"f": "0", "R": _TB_R, "S": _SZ_S},
              solves=(("r", "r"),), crunch=True,
              lines=(("surface", "r", "1", "the surface of the cloud, $r = r_b$"),),
              marked=(("event", {"x0": "2*sqrt(2)/3*((11/10)**Rational(3, 2) - (1/2)**Rational(3, 2))", "r": "11/10"},
                       1, "the last ray along this half of the axis to reach infinity"),),
              input=SZ_INPUT)
      for view, label, theta in (("north", "the axis, $\\theta = 0$", "0"), ("south", "the axis, $\\theta = \\pi$", "pi"))],
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
            (0, 4, -3, 1.5), "$r/r_s$", "$(cv - r)/r_s$", {"G": 1}, EQUATOR, to_display=FINKELSTEIN_IN,
            tau="v - r", areal=True, functions={"m": "Heaviside(v)/2"}, lines=(("shell", "x0", "0", "the shell, $v = 0$"),),
            marked=(("event", "1/1000", 1, "the event horizon"),), singular_runs=True,
            input="$m(v) = 0$ for $v < 0$ and $M$ for $v > 0$, with $r_s = 2GM/c^2$: a shell of null "
                  "dust of mass $M$ falling in along $v = 0$."),
    Diagram("vaidya", "eddington_finkelstein_outgoing", "shell", "an exploding shell", ("u", "r"),
            (0, 4, -1.5, 3), "$r/r_s$", "$(cu + r)/r_s$", {"G": 1}, EQUATOR, to_display=FINKELSTEIN_OUT,
            tau="u + r", areal=True, functions={"m": "Heaviside(-u)/2"}, lines=(("shell", "x0", "0", "the shell, $u = 0$"),),
            marked=(("event", "-1/1000", 0, "the white hole's horizon"),), singular_runs=True,
            input="$m(u) = M$ for $u < 0$ and $0$ for $u > 0$, with $r_s = 2GM/c^2$: a shell of null "
                  "dust carrying off the whole mass $M$ along $u = 0$."),
    # Bonnor and Vaidya's charged shell, as the conformal diagram declares it: a shell of charged
    # null dust of mass M and charge q = 24M/25 falls in along v = 0, flat inside and
    # Reissner-Nordstrom's metric outside, with r_+ = 1.28 M and r_- = 0.72 M. The shell's energy
    # density changes sign at r = q^2/2M, which the row marks. The outgoing chart draws its time
    # reverse. Inside r = sqrt(M^2 + q^2) - M the slices of constant v - r are timelike, so the
    # cones are oriented by the family that keeps the null time. On the shell itself both steps take
    # the flat side's value: with sympy's own 1/2 there, q^2 would be a quarter where m is a half.
    Diagram("bonnor_vaidya", "eddington_finkelstein_ingoing", "shell", "a charged shell falling in", ("v", "r"),
            (0, 4, -4, 2), "$r/M$", "$(v - r)/M$", {}, EQUATOR, to_display=FINKELSTEIN_IN,
            tau="v - r", orient="ingoing", areal=True, functions={"m": "Heaviside(v, 0)", "q": BV_CHARGE + "*Heaviside(v, 0)"},
            lines=(("shell", "x0", "0", "the shell, $v = 0$"),),
            points=(("mark", ("0", BV_TURN), "the sphere $r = q^2/2m$, where the energy density of the shell "
                                             "changes sign"),),
            marked=(("event", "1/1000", 1, "the event horizon"),), singular_runs=True, input=BV_INPUT["ingoing"]),
    Diagram("bonnor_vaidya", "eddington_finkelstein_outgoing", "shell", "a charged shell leaving", ("u", "r"),
            (0, 4, -2, 4), "$r/M$", "$(u + r)/M$", {}, EQUATOR, to_display=FINKELSTEIN_OUT,
            tau="u + r", orient="outgoing", areal=True, functions={"m": "Heaviside(-u, 0)", "q": BV_CHARGE + "*Heaviside(-u, 0)"},
            lines=(("shell", "x0", "0", "the shell, $u = 0$"),),
            points=(("mark", ("0", BV_TURN), "the sphere $r = q^2/2m$, where the energy density of the shell "
                                             "changes sign"),),
            marked=(("event", "-1/1000", 0, "the white hole's horizon"),), singular_runs=True,
            input=BV_INPUT["outgoing"]),
    # The homothetic collapse, m = mu v and q = (Q/M) mu v, in the chart that makes it static up to a
    # conformal factor: at Q = 24M/25 and mu = 1/20 the three rays that keep their R stand at
    # R = 0.67 M, 1.83 M and 7.50 M, the zeros of the published g^RR, and the dust's density changes
    # sign on R = Q^2/M.
    Diagram("bonnor_vaidya", "homothetic", "scaling", "the same at every scale", ("V", "R"), (0, 9, -4, 5),
            "$R/M$", "$V/M$", {"M": 1, "Q": BV_CHARGE, "mu": "1/20"}, EQUATOR, orient="ingoing", areal=True,
            lines=(("shell", "r", "576/625", "the sphere $R = Q^2/M$, where the energy density of the dust changes "
                                             "sign"),),
            input="$Q = 0.96\\,M$ and $\\mu = 1/20$: the mass grows as $m = v/20$ and the charge as $q = 0.96\\,m$."),
    # Robinson and Trautman's fronts, prolate at u = 0 and round by u = 3m, on the axis and on
    # the equator of the axisymmetric chart: the fronts are even about the equator, so d_theta H
    # vanishes on both and Gamma^theta_uu with it, and the null curves of each plane are null
    # geodesics. The equation runs toward the future only, so nothing is drawn before u = 0.
    Diagram("robinson_trautman", "axisymmetric", "axis", "the axis", ("u", "r"), (0, 5, 0, 9),
            "$r/m$", "$(cu + r)/m$", {}, {"theta": "0", "phi": "0"}, to_display=FINKELSTEIN_OUT, orient="outgoing",
            fronts=RT_FRONTS, input=RT_INPUT, singular_runs=True,
            lines=(("shell", "x0", "0", "the first front, $u = 0$"),)),
    Diagram("robinson_trautman", "axisymmetric", "equator", "the equator", ("u", "r"), (0, 5, 0, 9),
            "$r/m$", "$(cu + r)/m$", {}, EQUATOR, to_display=FINKELSTEIN_OUT, orient="outgoing",
            fronts=RT_FRONTS, input=RT_INPUT, singular_runs=True,
            lines=(("shell", "x0", "0", "the first front, $u = 0$"),)),
    # Kinnersley's photon rocket on the two halves of the axis it flies along, where sin(theta) = 0
    # kills g_u theta and every Gamma^theta of the plane, so the null curves are null geodesics.
    # The rectilinear chart measures theta in the rocket's rest frame from the direction opposite
    # to the acceleration, the Robinson-Trautman chart in the background frame from the direction
    # of motion, so the half behind the rocket is theta = 0 in the first and theta = pi in the second.
    *[Diagram("photon_rocket", system, view, label, ("u", "r"), (0, 12, -4, 22),
              "$r/m_0$", "$(cu + r)/m_0$", {}, {"theta": theta, "phi": "0"}, to_display=FINKELSTEIN_OUT,
              orient="outgoing", functions=functions, input=ROCKET_INPUT, singular_runs=True,
              lines=(("shell", "x0", "0", "the burn starts, $u = 0$"),
                     ("shell", "x0", "10", "the burn ends, $cu = 10\\,m_0$")))
      for system, functions, views in (
          ("rectilinear", ROCKET_BURN,
           (("behind", "behind the rocket, $\\theta = 0$", "0"), ("ahead", "ahead of the rocket, $\\theta = \\pi$", "pi"))),
          ("robinson_trautman", ROCKET_FLIGHT,
           (("ahead", "ahead of the rocket, $\\theta = 0$", "0"), ("behind", "behind the rocket, $\\theta = \\pi$", "pi"))))
      for view, label, theta in views],
    # The C-metric on the two halves of its axis, where sin(theta) = 0 kills Gamma^theta_tt and
    # Gamma^theta_rr and the null curves of the plane are null geodesics. Inside 2m and beyond
    # 1/alpha the time function is r, taken the way Griffiths, Krtous and Podolsky's extensions
    # take it: the black hole inside 2m and, beyond 1/alpha, the region to the future of the
    # acceleration horizon, which ends on null infinity. The Hong-Teo plane of the inner axis
    # runs on through y = 0, r = infinity, to null infinity at y = -1.
    Diagram("c_metric", "spherical", "inner", "the inner axis, $\\theta = 0$", ("t", "r"), (0, 9, -4.5, 4.5),
            "$r/m$", "$ct/m$", {"m": 1, "alpha": "1/6", "C": "3/4"}, {"theta": "0", "phi": "0"},
            tau="Piecewise((-r, r < 2), (t, r < 6), (r, True))"),
    Diagram("c_metric", "spherical", "outer", "the outer axis, $\\theta = \\pi$", ("t", "r"), (0, 6, -3, 3),
            "$r/m$", "$ct/m$", {"m": 1, "alpha": "1/6", "C": "3/4"}, {"theta": "pi", "phi": "0"},
            tau="Piecewise((-r, r < 2), (t, True))"),
    Diagram("c_metric", "hong_teo", "inner", "the inner axis, $x = 1$", ("\\tau", "y"), (-1, 5, -3, 3),
            "$y$", "$\\tau$", {"m": 1, "alpha": "1/6", "C": "3/4"}, {"x": "1", "phi": "0"},
            tau="Piecewise((-y, y < 1), (tau, y < 3), (y, True))",
            # y = 1/(alpha r) falls outward, so P, the family moving to smaller y, is the outgoing one.
            families=("outgoing", "ingoing")),
    # Bonnor's beam of light, in units of its radius. On the plane of the time and z at a fixed
    # place across the beam the rays moving with the beam keep ct - z and the null curves moving
    # against it keep ct + z + A (ct - z), A constant on the plane.
    Diagram("light_beam", "cartesian", "axis", "on the axis", ("t", "z"), LB_BOX, "$z/R$", "$ct/R$",
            {}, {"x": "0", "y": "0"}, families=LB_FAMILIES, functions={"A": LB_ONE}, input=LB_ONE_INPUT),
    Diagram("light_beam", "cartesian", "beside", "beside the beam, $x = 3R$", ("t", "z"), LB_BOX,
            "$z/R$", "$ct/R$", {}, {"x": "3", "y": "0"}, families=LB_FAMILIES, functions={"A": LB_ONE},
            input=LB_ONE_INPUT),
    Diagram("light_beam", "null_cartesian", "midway", "midway between two beams", ("u", "v"), LB_BOX,
            "$z/R$", "$ct/R$", {}, {"x": "0", "y": "0"}, to_display=LB_NULL_TO_TZ, tau="u + v",
            families=LB_FAMILIES, functions={"A": LB_TWO}, input=LB_TWO_INPUT),
    Diagram("light_beam", "null_cartesian", "one", "on the axis of one of two beams", ("u", "v"), LB_BOX,
            "$z/R$", "$ct/R$", {}, {"x": "2", "y": "0"}, to_display=LB_NULL_TO_TZ, tau="u + v",
            families=LB_FAMILIES, functions={"A": LB_TWO}, input=LB_TWO_INPUT),
    Diagram("light_beam", "null_cylindrical_interior", "edge", "at the edge of the beam, $\\rho = R$", ("u", "v"),
            LB_BOX, "$z/R$", "$ct/R$", LB, {"rho": "1", "phi": "0"}, to_display=LB_NULL_TO_TZ,
            tau="u + v", families=LB_FAMILIES),
    *[Diagram("light_beam", "null_cylindrical_exterior", view, f"$\\rho = {n}R$", ("u", "v"),
              LB_BOX, "$z/R$", "$ct/R$", LB, {"rho": str(n), "phi": "0"}, to_display=LB_NULL_TO_TZ,
              tau="u + v", families=LB_FAMILIES)
      for view, n in (("twice", 2), ("four", 4))],
    # Siklos's waves. The free charts on the plane of u and v at a fixed place on the wave front,
    # drawn against (v - u)/2 and (u + v)/2 with the declared pulse; Kaigorodov's spacetime on the
    # plane of its null coordinate v and its depth with the spacelike Killing direction d/du divided
    # out and y = 0, a totally geodesic slice, in each chart that holds the region x > 0, the
    # Poincare chart with d/dx divided out instead; and the stationary region on its plane of u and v.
    *[Diagram("siklos", "siklos", view, label, ("u", "v"), (-2, 2, -2, 2), "$(v - u)/2L$", "$(u + v)/2L$",
              SIKLOS, {"x": x, "y": "0"}, to_display=NULL_TO_TR, tau="u + v", families=SIDEWAYS,
              functions=SIKLOS_PULSE, lines=SIKLOS_CREST, input=SIKLOS_INPUT)
      for view, label, x in (("near", "$x = L$", "1"), ("far", "$x = 3L/2$", "3/2"))],
    Diagram("siklos", "ozsvath_robinson_rozga", "centre", "$\\xi = \\eta = 0$", ("u", "v"), (-2, 2, -2, 2),
            "$(v - u)/2L$", "$(u + v)/2L$", SIKLOS, {"xi": "0", "eta": "0"}, to_display=NULL_TO_TR, tau="u + v",
            families=SIDEWAYS, functions=SIKLOS_DISC_PULSE, lines=SIKLOS_CREST, input=SIKLOS_DISC_INPUT),
    Diagram("siklos", "kaigorodov", "depth", "$v$ and $x$", ("v", "x"), (0, 2, -1.5, 1.5), "$x/L$", "$v/L$",
            SIKLOS, {"y": "0"}, tau="v", families=SIDEWAYS, quotient="u", cone=SIKLOS_CONE.format("u")),
    Diagram("siklos", "kaigorodov_poincare", "depth", "$t$ and $z$", ("t", "z"), (0, 2, -1.5, 1.5), "$z/L$",
            "$ct/L$", SIKLOS, {"y": "0"}, families=SIDEWAYS, quotient="x", cone=SIKLOS_CONE.format("x")),
    Diagram("siklos", "kaigorodov_horospheric", "depth", "$v$ and $\\rho$", ("v", "\\rho"), (-0.75, 1.5, -1.5, 1.5),
            "$\\rho/L$", "$v/L$", SIKLOS, {"y": "0"}, tau="v", families=SIDEWAYS, quotient="u",
            cone=SIKLOS_CONE.format("u")),
    Diagram("siklos", "kaigorodov_stationary", "plane", "$u$ and $v$", ("u", "v"), (-2, 2, -2, 2),
            "$(v - u)/2L$", "$(u + v)/2L$", SIKLOS, {"rho": "0", "y": "0"}, to_display=NULL_TO_TR, tau="u + v",
            families=SIDEWAYS),
    Diagram("siklos", "kaigorodov_homogeneous", "depth", "$U$ and $Z$", ("U", "Z"), (-0.5, 0.35, -1.5, 1.5),
            "$Z$", "$-U/L$", {"L": 1, "k": "1/10"}, {"y": "0"}, to_display=((0, 1), (-1, 0)), tau="-U*exp(5*Z)",
            families=SIDEWAYS, quotient="X",
            cone=SIKLOS_CONE.format("X")),
    Diagram("siklos", "kaigorodov_kundt", "depth", "$V$ and $x$", ("V", "x"), (0.25, 2, -1.5, 1.5),
            "$x$", "$V/L$", SIKLOS, {"y": "0"}, tau="V*x**2", families=SIDEWAYS, quotient="U",
            cone=SIKLOS_CONE.format("U")),
]


# ---------------------------------------------------------------- the captions

def _as_caption(n):
    """The caption of the Aichelburg-Sexl view at rho = rho_0/n."""
    return [
        f"The plane of $u$ and $v$ ($x = \\rho_0/{n}$, $y = 0$), at distance $\\rho = \\rho_0/{n}$ from the axis "
        "the source moves along, drawn with $z = (v - u)/2$ and $ct = (u + v)/2$. Off the shock $u = 0$ the "
        "metric on this plane is $-du\\,dv$, flat, and light runs at 45°. The rays moving right keep their "
        "$u$ and run beside the shock without crossing it. Each ray moving left crosses it and comes out "
        f"moved along it by $\\Delta v = -(8GE/c^4)\\ln(\\rho/\\rho_0) = (8GE/c^4)\\ln {n}$, a delay of "
        f"$(4GE/c^4)\\ln {n}$ in $ct$.",
        "The Christoffel symbol $\\Gamma^x{}_{uu}$ turns every light ray crossing the shock toward the axis, out "
        "of this plane, so the curves drawn are null curves, and null geodesics everywhere off the shock. "
        "The jump grows by $(8GE/c^4)\\ln 4$ each time $\\rho$ is divided by $4$, and a change of "
        "$\\rho_0$ moves every ray behind the shock by the same amount.",
    ]


def _zipoy_voorhees_captions(system):
    """The four captions of one chart of Zipoy and Voorhees's metric: its axis and its equatorial plane,
    oblate and prolate, in the chart's own radius and parameter."""
    if system == "spherical":
        r, f, edge, neck = "r", "$f = 1 - 2m/r$", "r = 2m", "r = 3m"
        axis, equator = "$\\theta = 0$, $\\phi = 0$", "$\\theta = \\pi/2$, $\\phi = 0$"
        values = ("$q = 1$, $m = 1$", "$q = -1/2$, $m = 1$")
        h, m2, dr = "$h = (r - m)^2/(r(r - 2m))$", "", "dr"
        star = ("r_* = r + 4m\\ln(r/2m - 1) - 4m^2/(r - 2m)",
                "r_* = \\sqrt{r(r - 2m)} + 2m\\ln\\left(\\sqrt{r} + \\sqrt{r - 2m}\\right)")
        slope = ("r^{7/2}/(\\sqrt{r - 2m}\\,(r - m)^3)", "f^{-7/8}(1 - m/r)^{3/4}")
        K = ("192m^2(r - 3m)^2/r^8", "12m^2(r - 3m/2)^2/(r^5(r - 2m)^3)", "(r - 2m)^{-6}", "(r - 2m)^{-3/2}")
        at_neck, at_edge = "r = 3m", "3/(4m^4)"
    else:
        r, f, edge, neck = "x", "$f = (x - 1)/(x + 1)$", "x = 1", "x = 2"
        axis, equator = "$y = 1$, $\\phi = 0$", "$y = 0$, $\\phi = 0$"
        values = ("$\\delta = 2$, $m = 1$", "$\\delta = 1/2$, $m = 1$")
        h, m2, dr = "$h = x^2/(x^2 - 1)$", "m^2", "dx"
        star = ("r_* = m\\left(x + 4\\ln(x - 1) - 4/(x - 1)\\right)",
                "r_* = m\\left(\\sqrt{x^2 - 1} + 2\\ln\\left(\\sqrt{x + 1} + \\sqrt{x - 1}\\right)\\right)")
        slope = ("m(x + 1)^{7/2}/(\\sqrt{x - 1}\\,x^3)", "m\\,f^{-7/8}(x/(x + 1))^{3/4}")
        K = ("192(x - 2)^2/(m^4(x + 1)^8)", "12(x - 1/2)^2/(m^4(x + 1)^5(x - 1)^3)", "(x - 1)^{-6}", "(x - 1)^{-3/2}")
        at_neck, at_edge = "x = 2", "3/(4m^4)"
    geodesics = "and no Christoffel symbol turns them out of the plane, so they are null geodesics."
    one = "and no Christoffel symbol turns it out of the plane, so the rays are null geodesics."
    return {
        "axis_oblate": [
            f"The plane of $t$ and ${r}$ ({axis}) of the Zipoy-Voorhees metric ({values[0]}). The metric on it is "
            f"$-f^{{2}}c^2dt^2 + {m2}f^{{-2}}{dr}^2$ with {f}, so the rays are $ct = \\pm r_* + $ const with "
            f"${star[0]}$, {geodesics}",
            f"The cones close as $f^{{2}}$ toward ${edge}$, which a ray reaches only as $t \\to \\pm\\infty$, though after "
            f"a finite affine distance, since ${r}$ is an affine parameter along it. The Kretschmann scalar on the axis, "
            f"${K[0]}$, vanishes at ${at_neck}$ and is ${at_edge}$ at ${edge}$.",
        ],
        "equator_oblate": [
            f"The plane of $t$ and ${r}$ ({equator}) of the Zipoy-Voorhees metric ({values[0]}). The metric on it is "
            f"$-f^{{2}}c^2dt^2 + {m2}f^{{-2}}h^{{-3}}{dr}^2$ with {f} and {h}, so a ray has "
            f"$c\\,dt/{dr} = \\pm {slope[0]}$, {one}",
            f"The cones close toward ${edge}$, the ring, which every ingoing ray reaches in a finite time $t$ and where "
            f"the Kretschmann scalar diverges as ${K[2]}$. The circle about the axis is narrowest at ${neck}$.",
        ],
        "axis_prolate": [
            f"The plane of $t$ and ${r}$ ({axis}) of the Zipoy-Voorhees metric ({values[1]}). The metric on it is "
            f"$-f^{{1/2}}c^2dt^2 + {m2}f^{{-1/2}}{dr}^2$ with {f}, so the rays are $ct = \\pm r_* + $ const with "
            f"${star[1]}$, {geodesics}",
            f"The cones close as $f^{{1/2}}$ toward ${edge}$, which every ingoing ray reaches in a finite time $t$ and "
            f"where the Kretschmann scalar on the axis, ${K[1]}$, diverges.",
        ],
        "equator_prolate": [
            f"The plane of $t$ and ${r}$ ({equator}) of the Zipoy-Voorhees metric ({values[1]}). The metric on it is "
            f"$-f^{{1/2}}c^2dt^2 + {m2}f^{{-1/2}}h^{{3/4}}{dr}^2$ with {f} and {h}, so a ray has "
            f"$c\\,dt/{dr} = \\pm {slope[1]}$, {one}",
            f"The cones close toward ${edge}$, which every ingoing ray reaches in a finite time $t$ and where the "
            f"Kretschmann scalar diverges as ${K[3]}$. The circles about the axis shrink to zero there.",
        ],
    }


def _lb_drawn():
    return "drawn with $z = (v - u)/\\sqrt{2}$ and $ct = (u + v)/\\sqrt{2}$"


def _lb_outside_caption(n, last):
    """The caption of the view outside Bonnor's beam at rho = n R."""
    return [
        f"The plane of $u$ and $v$ outside the beam ($\\rho = {n}R$, $\\phi = 0$), {_lb_drawn()}, where "
        f"$g_{{uu}} = -(1 + 2\\ln {n})/4 = {-(1 + 2 * math.log(n)) / 4:.2f}$. The rays moving with the beam keep "
        "their $u$, run at 45°, and are null geodesics. A null curve moving against the beam keeps "
        "$v - g_{uu}u/2$.",
        last,
    ]


LB_CAPTIONS = {
    ("light_beam", "cartesian", "axis"): [
        "The plane of $t$ and $z$ on the axis of a uniform beam ($x = y = 0$), in units of its radius $R$. The "
        "profile $A$ and its gradient both vanish on the axis, so the metric on this plane is "
        "$-c^2dt^2 + dz^2$, flat, and light runs at 45° with the beam and against it.",
        "No Christoffel symbol is left on the axis, so every ray drawn is a null geodesic. A ray sent against the "
        "beam anywhere else is pulled toward the axis, since "
        "$\\ddot{x} = -\\tfrac{1}{2}\\partial_xA\\,(\\dot{t} - \\dot{z})^2$, and swings back and forth across it.",
    ],
    ("light_beam", "cartesian", "beside"): [
        "The plane of $t$ and $z$ beside a uniform beam ($x = 3R$, $y = 0$), where $A = (1 + 2\\ln 3)/8 = 0.40$. "
        "A ray moving with the beam keeps its $ct - z$, runs at 45°, and is a null geodesic whatever the profile. "
        "A null curve moving against the beam keeps $ct + z + A(ct - z)$, so it covers $(1 + A)/(1 - A) = 2.3$ "
        "units of $z$ for each unit of $ct$.",
        "That rate belongs to the coordinates: a constant added to $A$ amounts to the change of coordinates "
        "$ct + z \\to ct + z + \\text{const}\\,(ct - z)$, and only the difference in $A$ between two distances from "
        "the axis is free of that choice. A light ray launched against the beam along one of these curves is "
        "turned toward the axis by $\\Gamma^x{}_{tt}$, $\\Gamma^x{}_{tz}$, and $\\Gamma^x{}_{zz}$ and leaves the "
        "plane.",
    ],
    ("light_beam", "null_cartesian", "midway"): [
        "The plane of $u$ and $v$ midway between two parallel uniform beams ($x = y = 0$), each of radius $R$ with "
        f"its axis at $x = \\pm 2R$, {_lb_drawn()}. The two profiles add, $A = (1 + 2\\ln 2)/4 = 0.60$ here, and "
        "their gradients cancel, so no Christoffel symbol is left on this plane and every ray drawn is a null "
        "geodesic.",
        "The rays moving with the beams keep their $u$. The rays moving against them keep $v + Au$ and stay "
        "midway between the beams, pulled equally toward each.",
    ],
    ("light_beam", "null_cartesian", "one"): [
        "The plane of $u$ and $v$ on the axis of one of two parallel uniform beams ($x = 2R$, $y = 0$), the other "
        f"beam's axis at $x = -2R$, {_lb_drawn()}. A beam's own profile vanishes on its axis, so "
        "$A = (1 + 2\\ln 4)/8 = 0.47$ is the other beam's alone.",
        "The rays moving with the beams keep their $u$ and are null geodesics: light shining the same way as the "
        "two beams feels neither, and the beams themselves stay parallel. A light ray launched against them "
        "along a curve of constant $v + Au$ is turned toward the other beam by $\\Gamma^x{}_{uu} = \\partial_xA$ "
        "and leaves the plane.",
    ],
    ("light_beam", "null_cylindrical_interior", "edge"): [
        f"The plane of $u$ and $v$ at the edge of the beam ($\\rho = R$, $\\phi = 0$), {_lb_drawn()}, where "
        "$g_{uu} = -8\\pi G\\epsilon R^2/c^4 = -1/4$. The rays moving with the beam keep their $u$, run at 45°, "
        "and are null geodesics, the rays of the beam itself.",
        "A null curve moving against the beam keeps $v + u/8$. A light ray launched along one is turned toward "
        "the axis by $\\Gamma^\\rho{}_{uu} = 8\\pi G\\epsilon\\rho/c^4$, a pull in proportion to $\\rho$, so inside "
        "the beam every such ray swings about the axis as a pendulum does, with the period $4\\pi R$ in $u$.",
    ],
    ("light_beam", "null_cylindrical_exterior", "twice"): _lb_outside_caption(
        2, "Outside the beam spacetime is empty and $\\Gamma^\\rho{}_{uu} = 8\\pi G\\epsilon R^2/(c^4\\rho)$ falls "
           "off as $1/\\rho$, as the pull of a line of matter does in Newton's theory, so a light ray launched "
           "against the beam along one of these curves is turned toward the axis and leaves the plane."),
    ("light_beam", "null_cylindrical_exterior", "four"): _lb_outside_caption(
        4, "The profile grows as $\\ln\\rho$ without bound, so no ray sent against the beam escapes to "
           "infinity. One launched straight outward from here, with no angular momentum about the axis, turns back "
           "at a finite distance and falls into the beam."),
}

def _roberts_captions():
    """Roberts's fifteen views: a first paragraph for the chart, the same for each value of p, and a
    second for what that value of p does, said in the chart's own coordinates."""
    shown = {"disperses": "9/10", "threshold": "1", "collapses": "2"}
    chart = {
        "double_null": lambda p: (
            f"The plane of $u$ and $v$ ($\\theta = \\pi/2$, $\\phi = 0$), drawn for $p = {p}$ against $(v - u)/2$ and "
            "$(u + v)/2$ in units of any length $\\ell$, each point in the plane a 2-sphere of radius $R$. Only "
            "$g_{uv}$ is nonzero on it, so the light rays are the lines $u = $ const and $v = $ const, at 45°."),
        "advanced": lambda p: (
            f"The plane of Roberts's $v$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$), drawn for $p = {p}$ against $r$ and "
            "$(1 + p)v - r$ in units of any length $\\ell$. Ingoing rays are the lines $v = $ const, outgoing rays have "
            "$dr/dv = (1 + p)/2$, and in these axes both run at 45°."),
        "areal": lambda p: (
            f"The plane of $v$ and the areal radius $R$ ($\\theta = \\pi/2$, $\\phi = 0$), drawn for $p = {p}$ against "
            "$R$ and $(1 + p/2)v - R$ in units of any length $\\ell$. Ingoing rays are the lines $v = $ const, and "
            "outgoing rays have $dR/dv = (2\\lambda - p^2v)/4R$."),
        "diagonal": lambda p: (
            f"The plane of $t$ and $\\rho$ ($\\theta = \\pi/2$, $\\phi = 0$), drawn for $p = {p}$ in units of any length "
            "$\\ell$. The metric on it is $-c^2dt^2 + d\\rho^2$, so every ray is at 45°, and the collapse is all in the "
            "areal radius $R = \\sqrt{\\rho(\\rho - pct)}$."),
        "scaling": lambda p: (
            f"The plane of $\\tau$ and $x$ ($\\theta = \\pi/2$, $\\phi = 0$), drawn for $p = {p}$ against $x$ and "
            "$\\tau + x$. Outgoing rays are the lines $\\tau = $ const, and ingoing rays have "
            "$d\\tau/dx = 2/(1 - e^{-2x})$. A step of one in $\\tau$ at fixed $x$ shrinks every length by the factor $e$ "
            "and changes nothing else."),
    }
    outcome = {
        ("double_null", "disperses"):
            "The field arrives on the ray $v = 0$, with flat space before it. For $p < 1$ every outgoing ray gets "
            "away to infinity, the mass inside each sphere, $-p^2uv/8R$, is back to zero on the ray $u = 0$, and after "
            "that ray the spacetime is flat again.",
        ("double_null", "threshold"):
            "At $p = 1$ the areal radius is $R = \\tfrac{1}{2}\\sqrt{u(u - 2v)}$, which vanishes on the ray $u = 0$: the "
            "singularity is a light ray, on which the Kretschmann scalar diverges as $12/u^2v^2$. Every outgoing ray "
            "before it still reaches infinity, so no horizon forms.",
        ("double_null", "collapses"):
            "For $p = 2$ the singularity is the spacelike line $u = (1 - p)v$, here $u = -v$, and the apparent horizon is "
            "$u = (1 - p^2)v$, here $u = -3v$, where $\\partial_vR = 0$. Between the two $R$ shrinks along both families "
            "of rays, so every sphere there is trapped. The horizon's radius, $p\\sqrt{p^2 - 1}\\,v/2$, grows with $v$ "
            "for as long as the field keeps arriving.",
        ("advanced", "disperses"):
            "The field arrives on $v = 0$ and its last ray is $r = (1 + p)v/2$, with flat space before the one and after "
            "the other. The spheres have radius $R = \\sqrt{r(r - pv)}$, which stays positive between them.",
        ("advanced", "threshold"):
            "At $p = 1$ the last ray $r = (1 + p)v/2$ is $r = pv$, where $R = \\sqrt{r(r - pv)}$ vanishes and the "
            "Kretschmann scalar diverges: the singularity is a light ray.",
        ("advanced", "collapses"):
            "For $p = 2$ the apparent horizon is $r = p(1 + p)v/2$, here $r = 3v$, which is level in these axes, and "
            "the singularity $r = pv$ comes after it. Between the two both families of rays run to smaller $R$.",
        ("areal", "disperses"):
            "The field arrives on $v = 0$ and its last ray is $R = \\sqrt{1 - p^2}\\,v/2$, with flat space before the "
            "one and after the other. Outgoing rays gain $R$ everywhere, since $2\\lambda > p^2v$ for $p < 1$.",
        ("areal", "threshold"):
            "At $p = 1$ outgoing rays have $dR/dv = (2\\lambda - v)/4R$, which is positive, so each of them gains $R$ "
            "and no sphere is trapped. The singularity is $R = 0$ for $v > 0$, which every ingoing ray reaches.",
        ("areal", "collapses"):
            "For $p = 2$ the outgoing rays stand still in $R$ on $R = p\\sqrt{p^2 - 1}\\,v/2$, where $g^{RR} = 0$: the "
            "apparent horizon. Inside it both families of rays fall to the singularity $R = 0$.",
        ("diagonal", "disperses"):
            "The field fills the wedge $\\rho > c|t|$ between its first ray $ct = -\\rho$ and its last $ct = \\rho$, "
            "with flat space before the one and after the other. The areal radius is $\\rho$ at $t = 0$, larger than "
            "$\\rho$ before and smaller after.",
        ("diagonal", "threshold"):
            "At $p = 1$ the areal radius vanishes on the last ray $ct = \\rho$, where the Kretschmann scalar diverges: "
            "the singularity is a light ray.",
        ("diagonal", "collapses"):
            "For $p = 2$ the singularity is $ct = \\rho/p$ and the apparent horizon is $ct = (2 - p)\\rho/p$, which at "
            "this $p$ is the moment $t = 0$. After that moment $R$ shrinks along both families of rays.",
        ("scaling", "disperses"):
            "The field's first ray $v = 0$ is $x = 0$, and its last ray $u = 0$ is $\\tau \\to \\infty$. The spheres have "
            "radius $R = \\ell e^{x - \\tau}\\sqrt{\\cosh^2x - p^2\\sinh^2x}$, which is positive at every $x$ for $p < 1$.",
        ("scaling", "threshold"):
            "At $p = 1$ the areal radius is $R = \\ell e^{x - \\tau}$ and the field is $\\varphi = -x$. The singularity "
            "$u = 0$ lies at $\\tau \\to \\infty$ and $x \\to \\infty$ together, with $2x - \\tau$ finite, where "
            "$R \\to 0$.",
        ("scaling", "collapses"):
            "For $p = 2$ the singularity is the line $x = \\tanh^{-1}(1/p)$, here $x = 0.55$, and the apparent horizon "
            "is the line $x = \\tanh^{-1}(1/p^2)$, here $x = 0.26$. Both stand at the same $x$ for every $\\tau$, since "
            "the collapse is the same at every scale.",
    }
    return {("roberts", system, case): [chart[system](shown[case]), text] for (system, case), text in outcome.items()}


ROBERTS_CAPTIONS = _roberts_captions()

CAPTIONS = {
    **LB_CAPTIONS,
    **{("aichelburg_sexl", "null_cartesian", view): _as_caption(rho[2:]) for view, rho in AS_RHO.items()},
    ("schwarzschild", "spherical", "radial"): [
        "The plane of $t$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$), the same at every fixed angle "
        "by spherical symmetry. Outside $r_s$ the cones narrow toward "
        "the vertical as $r \\to r_s$, because $dt/dr = \\pm(1 - r_s/r)^{-1}$ diverges there. The "
        "ingoing family piles up against $r_s$ toward $t \\to +\\infty$, and the outgoing family "
        "peels away from it out of $t \\to -\\infty$; the two Eddington-Finkelstein charts carry "
        "them across.",
        "The domain of this chart stops at $r_s$. Evaluated inside it, the same components give "
        "cones lying on their side, because there $r$ is the time. The components alone do not fix "
        "whether that region is the black hole or the white hole; we take the future from the "
        "ingoing Eddington-Finkelstein chart, which runs smoothly across $r_s$, and this makes it "
        "the black hole, where every cone points to $r = 0$. The Kretschmann scalar $12r_s^2/r^6$ "
        "is finite at $r_s$ and diverges only at $r = 0$.",
    ],
    ("schwarzschild", "eddington_finkelstein_ingoing", "finkelstein"): [
        "The plane of $v$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$), drawn with $v - "
        "r$ as the vertical axis so that the ingoing rays, $v = $ const, run at 45°. The chart "
        "crosses the horizon. The outgoing family has $dv/dr = 2(1 - r_s/r)^{-1}$, so it stands "
        "exactly vertical at $r_s$: the horizon is itself an outgoing ray that stays where it is.",
        "The cones cross $r_s$ smoothly and keep tipping. Inside, both edges of every future cone "
        "point to smaller $r$, so every future directed ray ends at $r = 0$.",
    ],
    ("schwarzschild", "eddington_finkelstein_ingoing", "chart"): [
        "The same plane of $v$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$), drawn "
        "against the chart's own coordinates. The ingoing family is $v = $ const and runs "
        "horizontally here, since $v$ is itself a null coordinate. The outgoing family turns "
        "vertical at $r_s$ and leans back toward smaller $r$ inside it.",
    ],
    ("schwarzschild", "eddington_finkelstein_outgoing", "finkelstein"): [
        "The plane of $u$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$), drawn with $u + "
        "r$ as the vertical axis so that the outgoing rays, $u = $ const, run at 45°. The retarded "
        "chart crosses the other horizon. Inside $r_s$ both edges of every future cone point to "
        "larger $r$: this is the white hole, which nothing from outside can enter.",
        "The chart is the time reverse of the ingoing one: its $g_{ur}$ is $-1$ where the ingoing "
        "chart's $g_{vr}$ is $+1$.",
    ],
    ("schwarzschild", "eddington_finkelstein_outgoing", "chart"): [
        "The same plane of $u$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$), drawn "
        "against the chart's own coordinates. The outgoing family is $u = $ const and runs "
        "horizontally here, since $u$ is itself a null coordinate. The ingoing family turns "
        "vertical at $r_s$ and leans toward larger $r$ inside it.",
    ],
    ("schwarzschild_de_sitter", "static", "radial"): [
        "The plane of $t$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$), drawn for $\\Lambda = 0.2/r_s^2$, "
        "the same at every other angle by spherical symmetry. There $g^{rr} = 1 - r_s/r - \\Lambda r^2/3$ "
        "vanishes at the black hole horizon $r_h = 1.085\\,r_s$ and at the cosmological horizon "
        "$r_c = 3.215\\,r_s$, and the cones close at both, since $dt/dr = \\pm(1 - r_s/r - \\Lambda r^2/3)^{-1}$ "
        "diverges there. Between them the cones are widest at $r = (3r_s/2\\Lambda)^{1/3} = 1.957\\,r_s$, "
        "where $g^{rr}$ is greatest.",
        "Inside $r_h$ and beyond $r_c$, $t$ is a spacelike coordinate, and the components alone do not fix "
        "which way is future. We take it from the ingoing Eddington-Finkelstein chart inside $r_h$, which "
        "makes that region the black hole, where every cone points to $r = 0$, and from the outgoing one "
        "beyond $r_c$, which makes that region the expanding universe, where every cone points to larger $r$. "
        "The Kretschmann scalar $12r_s^2/r^6 + 8\\Lambda^2/3$ diverges only at $r = 0$.",
    ],
    ("schwarzschild_de_sitter", "eddington_finkelstein_ingoing", "finkelstein"): [
        "The plane of $v$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$), drawn for $\\Lambda = 0.2/r_s^2$ with "
        "$v - r$ as the vertical axis, so that the ingoing rays, $v = $ const, run at 45°. The outgoing family "
        "has $dv/dr = 2(1 - r_s/r - \\Lambda r^2/3)^{-1}$, so it stands vertical at both horizons: each horizon "
        "is an outgoing ray that stays where it is.",
        "The chart crosses the black hole horizon $r_h = 1.085\\,r_s$ into the black hole, where both edges of "
        "every future cone point to smaller $r$. It crosses the cosmological horizon $r_c = 3.215\\,r_s$ into "
        "the contracting region in the past of the static one, where both edges point to smaller $r$ as well.",
    ],
    ("schwarzschild_de_sitter", "eddington_finkelstein_ingoing", "chart"): [
        "The same plane of $v$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$), drawn against the chart's own "
        "coordinates. The ingoing family is $v = $ const and runs horizontally here, since $v$ is itself a "
        "null coordinate. The outgoing family turns vertical at both horizons and leans back toward smaller "
        "$r$ inside $r_h$ and beyond $r_c$.",
    ],
    ("schwarzschild_de_sitter", "eddington_finkelstein_outgoing", "finkelstein"): [
        "The plane of $u$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$), drawn for $\\Lambda = 0.2/r_s^2$ with "
        "$u + r$ as the vertical axis, so that the outgoing rays, $u = $ const, run at 45°. The ingoing family "
        "stands vertical at both horizons. Inside $r_h$ both edges of every future cone point to larger $r$: "
        "this is the white hole, which nothing from outside can enter. Beyond $r_c$ they point to larger $r$ "
        "as well, into the expanding region that the static observers' light goes on to reach.",
        "The chart is the time reverse of the ingoing one: its $g_{ur}$ is $-1$ where the ingoing chart's "
        "$g_{vr}$ is $+1$.",
    ],
    ("schwarzschild_de_sitter", "eddington_finkelstein_outgoing", "chart"): [
        "The same plane of $u$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$), drawn against the chart's own "
        "coordinates. The outgoing family is $u = $ const and runs horizontally here, since $u$ is itself a "
        "null coordinate. The ingoing family turns vertical at both horizons and leans toward larger $r$ "
        "inside $r_h$ and beyond $r_c$.",
    ],
    ("kiselev", "static", "radial"): [
        "The plane of $t$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$), drawn for $w = -2/3$ and $r_q = 8\\,r_s$, "
        "the same at every other angle by spherical symmetry. There $g^{rr} = 1 - r_s/r - r/r_q$ "
        "vanishes at the black hole horizon $r_- = 1.172\\,r_s$ and at the outer horizon "
        "$r_+ = 6.828\\,r_s$, and the cones close at both, since $dt/dr = \\pm(1 - r_s/r - r/r_q)^{-1}$ "
        "diverges there. Between them the cones are widest at $r = \\sqrt{r_sr_q} = 2.828\\,r_s$, "
        "where $g^{rr}$ is greatest.",
        "Inside $r_-$ and beyond $r_+$, $t$ is a spacelike coordinate, and the components alone do not fix "
        "which way is future. We take it from the ingoing Eddington-Finkelstein chart inside $r_-$, which "
        "makes that region the black hole, where every cone points to $r = 0$, and from the outgoing one "
        "beyond $r_+$, where every cone points to larger $r$. "
        "The Kretschmann scalar $12r_s^2/r^6 + 8/(r^2r_q^2)$ diverges only at $r = 0$.",
    ],
    ("kiselev", "eddington_finkelstein_ingoing", "finkelstein"): [
        "The plane of $v$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$), drawn for $w = -2/3$ and $r_q = 8\\,r_s$ with "
        "$v - r$ as the vertical axis, so that the ingoing rays, $v = $ const, run at 45°. The outgoing family "
        "has $dv/dr = 2(1 - r_s/r - r/r_q)^{-1}$, so it stands vertical at both horizons: each horizon "
        "is an outgoing ray that stays where it is.",
        "The chart crosses the black hole horizon $r_- = 1.172\\,r_s$ into the black hole, where both edges of "
        "every future cone point to smaller $r$. It crosses the outer horizon $r_+ = 6.828\\,r_s$ into "
        "the region in the past of the static one, where both edges point to smaller $r$ as well.",
    ],
    ("kiselev", "eddington_finkelstein_ingoing", "chart"): [
        "The same plane of $v$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$), drawn against the chart's own "
        "coordinates. The ingoing family is $v = $ const and runs horizontally here, since $v$ is itself a "
        "null coordinate. The outgoing family turns vertical at both horizons and leans back toward smaller "
        "$r$ inside $r_-$ and beyond $r_+$.",
    ],
    ("kiselev", "eddington_finkelstein_outgoing", "finkelstein"): [
        "The plane of $u$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$), drawn for $w = -2/3$ and $r_q = 8\\,r_s$ with "
        "$u + r$ as the vertical axis, so that the outgoing rays, $u = $ const, run at 45°. The ingoing family "
        "stands vertical at both horizons. Inside $r_-$ both edges of every future cone point to larger $r$: "
        "this is the white hole, which nothing from outside can enter. Beyond $r_+$ they point to larger $r$ "
        "as well, into the region that the static observers' light goes on to reach.",
        "The chart is the time reverse of the ingoing one: its $g_{ur}$ is $-1$ where the ingoing chart's "
        "$g_{vr}$ is $+1$.",
    ],
    ("kiselev", "eddington_finkelstein_outgoing", "chart"): [
        "The same plane of $u$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$), drawn against the chart's own "
        "coordinates. The outgoing family is $u = $ const and runs horizontally here, since $u$ is itself a "
        "null coordinate. The ingoing family turns vertical at both horizons and leans toward larger $r$ "
        "inside $r_-$ and beyond $r_+$.",
    ],
    ("kiselev", "linear", "radial"): [
        "The plane of $t$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$), drawn for $r_q = 8\\,r_s$, "
        "the same at every other angle by spherical symmetry. There $g^{rr} = 1 - r_s/r - r/r_q$ "
        "vanishes at the black hole horizon $r_- = 1.172\\,r_s$ and at the outer horizon "
        "$r_+ = 6.828\\,r_s$, and the cones close at both, since $dt/dr = \\pm(1 - r_s/r - r/r_q)^{-1}$ "
        "diverges there. Between them the cones are widest at $r = \\sqrt{r_sr_q} = 2.828\\,r_s$, "
        "where $g^{rr}$ is greatest.",
        "Inside $r_-$ and beyond $r_+$, $t$ is a spacelike coordinate, and the components alone do not fix "
        "which way is future. We take it from the ingoing Eddington-Finkelstein chart inside $r_-$, which "
        "makes that region the black hole, where every cone points to $r = 0$, and from the outgoing one "
        "beyond $r_+$, where every cone points to larger $r$. "
        "The Kretschmann scalar $12r_s^2/r^6 + 8/(r^2r_q^2)$ diverges only at $r = 0$.",
    ],
    ("kiselev", "hyperbolic", "radial"): [
        "The plane of $\\eta$ and $\\chi$ ($\\theta = \\pi/2$, $\\phi = 0$), with no black hole and $w = -2/3$. "
        "The plane is conformally flat, so every radial null ray is a straight line at 45°, "
        "$\\eta \\mp \\chi = $ const. The centre $\\chi = 0$ is a curvature singularity, where the Kretschmann "
        "scalar $2e^{2\\chi}/(r_q^4\\sinh^2\\chi)$ diverges, and the areal radius $r = r_q(1 - e^{-2\\chi})$ "
        "reaches the horizon $r = r_q$ only as $\\chi \\to \\infty$.",
    ],
    ("kiselev", "conformally_flat", "radial"): [
        "The plane of $\\tau$ and $\\rho$ ($\\theta = \\pi/2$, $\\phi = 0$), with no black hole and $w = -2/3$. "
        "The metric is $4r_q^2/(\\tau + \\rho)^2$ times Minkowski's, so every radial null ray is a straight line "
        "at 45°. The static region is $\\tau > \\rho$, and the horizon $r = r_q$ is the ray $\\tau = \\rho$ "
        "from the origin.",
        "Below the horizon the areal radius $r = 2r_q\\rho/(\\tau + \\rho)$ exceeds $r_q$, and it grows without "
        "bound toward $\\tau + \\rho = 0$, which is past null infinity. The centre $\\rho = 0$ is a curvature "
        "singularity, where the Kretschmann scalar $2(\\tau + \\rho)^2/(r_q^4\\rho^2)$ diverges.",
    ],
    ("reissner_nordstrom_de_sitter", "static", "radial"): [
        "The plane of $t$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$), drawn for the lukewarm hole ($r_q = r_s/2$, $\\Lambda = 27/(64\\,r_s^2)$), "
        "the same at every other angle by spherical symmetry. There $g^{rr} = 1 - r_s/r + r_q^2/r^2 - \\Lambda r^2/3$ "
        "vanishes at the inner horizon $r_- = 0.431\\,r_s$, at the black hole horizon $r_+ = 2r_s/3$, and at the "
        "cosmological horizon $r_c = 2\\,r_s$, and the cones close at each of the three, since $dt/dr = \\pm 1/g^{rr}$ "
        "diverges there. Between $r_+$ and $r_c$ the cones are widest at $r = 1.298\\,r_s$, where $g^{rr}$ is greatest.",
        "Between $r_-$ and $r_+$ and beyond $r_c$, $t$ is a spacelike coordinate, and the components alone do not fix "
        "which way is future. We take it from the ingoing Eddington-Finkelstein chart inside $r_+$, which "
        "makes that region the black hole, where every cone points to smaller $r$, and from the outgoing one "
        "beyond $r_c$, which makes that region the expanding universe, where every cone points to larger $r$. "
        "Inside $r_-$ the coordinate $t$ is a time again, and the singularity $r = 0$, where the Kretschmann scalar "
        "diverges, is timelike.",
    ],
    ("reissner_nordstrom_de_sitter", "eddington_finkelstein_ingoing", "finkelstein"): [
        "The plane of $v$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$), drawn for the lukewarm hole ($r_q = r_s/2$, $\\Lambda = 27/(64\\,r_s^2)$) with "
        "$v - r$ as the vertical axis, so that the ingoing rays, $v = $ const, run at 45°. The outgoing family "
        "has $dv/dr = 2/g^{rr}$ with $g^{rr} = 1 - r_s/r + r_q^2/r^2 - \\Lambda r^2/3$, so it stands vertical at each of the "
        "three horizons: each horizon is an outgoing ray that stays where it is.",
        "The chart crosses the black hole horizon $r_+ = 2r_s/3$ into the black hole, where both edges of "
        "every future cone point to smaller $r$, and the inner horizon $r_- = 0.431\\,r_s$ into the region about the "
        "singularity, where the outgoing edge points to larger $r$ again. It crosses the cosmological horizon "
        "$r_c = 2\\,r_s$ into the contracting region in the past of the static one, where both edges point to "
        "smaller $r$ as well.",
    ],
    ("reissner_nordstrom_de_sitter", "eddington_finkelstein_ingoing", "chart"): [
        "The same plane of $v$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$), drawn against the chart's own "
        "coordinates. The ingoing family is $v = $ const and runs horizontally here, since $v$ is itself a "
        "null coordinate. The outgoing family turns vertical at each of the three horizons and leans back toward smaller "
        "$r$ between $r_-$ and $r_+$ and beyond $r_c$.",
    ],
    ("reissner_nordstrom_de_sitter", "eddington_finkelstein_outgoing", "finkelstein"): [
        "The plane of $u$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$), drawn for the lukewarm hole ($r_q = r_s/2$, $\\Lambda = 27/(64\\,r_s^2)$) with "
        "$u + r$ as the vertical axis, so that the outgoing rays, $u = $ const, run at 45°. The ingoing family "
        "stands vertical at each of the three horizons. Between $r_-$ and $r_+$ both edges of every future cone point to "
        "larger $r$: this is the white hole, which nothing from outside can enter. Beyond $r_c$ they point to "
        "larger $r$ as well, into the expanding region that the static observers' light goes on to reach.",
        "The chart is the time reverse of the ingoing one: its $g_{ur}$ is $-1$ where the ingoing chart's "
        "$g_{vr}$ is $+1$.",
    ],
    ("reissner_nordstrom_de_sitter", "eddington_finkelstein_outgoing", "chart"): [
        "The same plane of $u$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$), drawn against the chart's own "
        "coordinates. The outgoing family is $u = $ const and runs horizontally here, since $u$ is itself a "
        "null coordinate. The ingoing family turns vertical at each of the three horizons and leans toward larger $r$ "
        "between $r_-$ and $r_+$ and beyond $r_c$.",
    ],
    ("reissner_nordstrom_de_sitter", "cosmological", "plane"): [
        "The plane of $\\tau$ and $\\rho$ ($\\theta = \\pi/2$, $\\phi = 0$), drawn for $H = 3c/(8\\,r_s)$, the same at "
        "every other angle by spherical symmetry. The rays obey $d\\rho/d(c\\tau) = \\pm(H\\tau + r_s/2\\rho)^{-2}$, so "
        "$\\tau$ is a time everywhere and no cone closes. The marked curves, where $|\\nabla r|^2$ vanishes for the "
        "areal radius $r = H\\tau\\rho + r_s/2$, are the three horizons, the hyperbolas $c\\tau\\rho = 4\\,r_s^2$ of $r_c$, "
        "$4\\,r_s^2/9$ of $r_+$, and $-0.185\\,r_s^2$ of $r_-$.",
        "The chart crosses the cosmological horizon into the expanding region and the black hole horizon into "
        "the white hole, as the outgoing Eddington-Finkelstein chart does. The line $\\tau = 0$ is the sphere "
        "$r = r_s/2$ inside the white hole, and below it the chart runs on through the inner horizon to the "
        "singularity $r = 0$ on $H\\tau\\rho = -r_s/2$, where the Kretschmann scalar diverges.",
    ],
    ("reissner_nordstrom_ads", "static", "radial"): [
        "The plane of $t$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$), drawn for $r_s = 27L/8$ and $r_q^2 = 11L^2/8$, "
        "the same at every other angle by spherical symmetry. There $g^{rr} = 1 - r_s/r + r_q^2/r^2 + r^2/L^2$ "
        "vanishes at the inner horizon $r_- = L/2$ and at the event horizon $r_+ = L$, and the cones close at "
        "both, since $dt/dr = \\pm 1/g^{rr}$ diverges there. Far outside, $g^{rr}$ grows as $r^2/L^2$ and the "
        "cones open toward the horizontal: a light ray runs from any radius to $r \\to \\infty$ in a finite time "
        "$t$, as in anti-de Sitter space.",
        "Between $r_-$ and $r_+$, $t$ is a spacelike coordinate, and the components alone do not fix which way is "
        "future. We take it from the ingoing Eddington-Finkelstein chart, which makes that region the black hole, "
        "where every cone points to smaller $r$. Inside $r_-$ the coordinate $t$ is a time again, the cones open "
        "wide as $g^{rr}$ grows as $r_q^2/r^2$, and the singularity $r = 0$, where the Kretschmann scalar diverges, "
        "is timelike.",
    ],
    ("reissner_nordstrom_ads", "eddington_finkelstein_ingoing", "finkelstein"): [
        "The plane of $v$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$), drawn for $r_s = 27L/8$ and $r_q^2 = 11L^2/8$ "
        "with $v - r$ as the vertical axis, so that the ingoing rays, $v = $ const, run at 45°. The outgoing family "
        "has $dv/dr = 2/g^{rr}$ with $g^{rr} = 1 - r_s/r + r_q^2/r^2 + r^2/L^2$, so it stands vertical at both "
        "horizons: each horizon is an outgoing ray that stays where it is.",
        "The chart crosses the event horizon $r_+ = L$ into the black hole, where both edges of every future cone "
        "point to smaller $r$, and the inner horizon $r_- = L/2$ into the region about the singularity, where the "
        "outgoing edge points to larger $r$ again. Far outside, and near $r = 0$, the outgoing edge leans toward "
        "the ingoing one, since $dv/dr$ falls as $2L^2/r^2$ and as $2r^2/r_q^2$.",
    ],
    ("reissner_nordstrom_ads", "eddington_finkelstein_ingoing", "chart"): [
        "The same plane of $v$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$), drawn against the chart's own "
        "coordinates. The ingoing family is $v = $ const and runs horizontally here, since $v$ is itself a "
        "null coordinate. The outgoing family turns vertical at both horizons and leans back toward smaller "
        "$r$ between them.",
    ],
    ("reissner_nordstrom_ads", "eddington_finkelstein_outgoing", "finkelstein"): [
        "The plane of $u$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$), drawn for $r_s = 27L/8$ and $r_q^2 = 11L^2/8$ "
        "with $u + r$ as the vertical axis, so that the outgoing rays, $u = $ const, run at 45°. The ingoing family "
        "stands vertical at both horizons. Between $r_-$ and $r_+$ both edges of every future cone point to "
        "larger $r$: this is the white hole, which nothing from outside can enter.",
        "The chart is the time reverse of the ingoing one: its $g_{ur}$ is $-1$ where the ingoing chart's "
        "$g_{vr}$ is $+1$.",
    ],
    ("reissner_nordstrom_ads", "eddington_finkelstein_outgoing", "chart"): [
        "The same plane of $u$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$), drawn against the chart's own "
        "coordinates. The outgoing family is $u = $ const and runs horizontally here, since $u$ is itself a "
        "null coordinate. The ingoing family turns vertical at both horizons and leans toward larger $r$ "
        "between them.",
    ],
    ("tangherlini", "spherical", "radial"): [
        "The plane of $t$ and $r$ ($\\psi = \\theta = \\pi/2$, $\\phi = 0$) in five dimensions, the same at every fixed angle by hyperspherical "
        "symmetry. Outside $r_h$ the cones narrow toward the vertical as $r \\to r_h$, because "
        "$dt/dr = \\pm(1 - r_h^2/r^2)^{-1}$ diverges there. Away from the horizon they open faster than "
        "Schwarzschild's, since $r_h^2/r^2$ falls faster than Schwarzschild's $r_s/r$ with $r_s$ its own horizon radius.",
        "Inside $r_h$ the same components make $r$ the time. We take the future from the ingoing "
        "Eddington-Finkelstein chart, which runs smoothly across $r_h$, and this makes that region the black "
        "hole, where every cone points to $r = 0$. The Kretschmann scalar $72r_h^4/r^8$ is finite at $r_h$ and "
        "diverges only at $r = 0$.",
    ],
    ("tangherlini", "eddington_finkelstein_ingoing", "finkelstein"): [
        "The plane of $v$ and $r$ ($\\psi = \\theta = \\pi/2$, $\\phi = 0$) in five dimensions, drawn with $v - r$ as the vertical axis so that the "
        "ingoing rays, $v = $ const, run at 45°. The outgoing family has $dv/dr = 2(1 - r_h^2/r^2)^{-1}$, so it "
        "stands exactly vertical at $r_h$: the horizon is itself an outgoing ray that stays where it is.",
        "The cones cross $r_h$ smoothly and keep tipping. Inside, both edges of every future cone point to "
        "smaller $r$, so every future directed ray ends at $r = 0$.",
    ],
    ("tangherlini", "eddington_finkelstein_ingoing", "chart"): [
        "The same plane of $v$ and $r$ ($\\psi = \\theta = \\pi/2$, $\\phi = 0$), drawn against the chart's own coordinates. The ingoing family is "
        "$v = $ const and runs horizontally here, since $v$ is itself a null coordinate. The outgoing family "
        "turns vertical at $r_h$ and leans back toward smaller $r$ inside it.",
    ],
    ("tangherlini", "eddington_finkelstein_outgoing", "finkelstein"): [
        "The plane of $u$ and $r$ ($\\psi = \\theta = \\pi/2$, $\\phi = 0$) in five dimensions, drawn with $u + r$ as the vertical axis so that the "
        "outgoing rays, $u = $ const, run at 45°. The ingoing family stands vertical at $r_h$. Inside $r_h$ both "
        "edges of every future cone point to larger $r$: this is the white hole, which nothing from outside can "
        "enter.",
        "The chart is the time reverse of the ingoing one: its $g_{ur}$ is $-1$ where the ingoing chart's "
        "$g_{vr}$ is $+1$.",
    ],
    ("tangherlini", "eddington_finkelstein_outgoing", "chart"): [
        "The same plane of $u$ and $r$ ($\\psi = \\theta = \\pi/2$, $\\phi = 0$), drawn against the chart's own coordinates. The outgoing family is "
        "$u = $ const and runs horizontally here, since $u$ is itself a null coordinate. The ingoing family "
        "turns vertical at $r_h$ and leans toward larger $r$ inside it.",
    ],
    ("boulware_deser", "spherical", "radial"): [
        "The plane of $t$ and $r$ ($\\psi = \\theta = \\pi/2$, $\\phi = 0$), drawn for $r_0 = 13r_h/12$ and "
        "$\\ell = 5r_h/12$, the same at every fixed angle by hyperspherical symmetry. Outside $r_h$ the cones narrow "
        "toward the vertical as $r \\to r_h$, because $dt/dr = \\pm 1/f$ with $f = 1 - 2r_0^2/(r^2 + W)$ diverges "
        "there. At $r = 2r_h$ it is $\\pm 1.41$, where Tangherlini's at twice its own horizon radius is $\\pm 4/3$.",
        "Inside $r_h$ the same components make $r$ the time, and we take the future as the ingoing rays carry it "
        "across the horizon, which makes that region the black hole, where every cone points to $r = 0$. There "
        "$f$ ends at the finite value $1 - r_0/\\ell = -8/5$, so the cones reach the singularity with their edges "
        "at $c\\,dt/dr = \\pm 5/8$, while the Kretschmann scalar grows as $12r_0^2/(\\ell^2r^4)$.",
    ],
    ("boulware_deser", "eddington_finkelstein_ingoing", "finkelstein"): [
        "The plane of $v$ and $r$ ($\\psi = \\theta = \\pi/2$, $\\phi = 0$), drawn for $r_0 = 13r_h/12$ and "
        "$\\ell = 5r_h/12$ with $v - r$ as the vertical axis, so that the ingoing rays, $v = $ const, run at 45°. The "
        "outgoing family has $dv/dr = 2/f$ with $f = 1 - 2r_0^2/(r^2 + W)$, so it stands exactly vertical at $r_h$: "
        "the horizon is an outgoing ray that stays where it is.",
        "The cones cross $r_h$ smoothly and keep tipping. Inside, both edges of every future cone point to "
        "smaller $r$, so every future directed ray on this plane ends at $r = 0$. The outgoing edge arrives there "
        "at $dv/dr = 2/(1 - r_0/\\ell) = -5/4$, with the cone still open, where Tangherlini's cone has closed onto "
        "its ingoing edge.",
    ],
    ("boulware_deser", "eddington_finkelstein_ingoing", "chart"): [
        "The same plane of $v$ and $r$ ($\\psi = \\theta = \\pi/2$, $\\phi = 0$), drawn against the chart's own "
        "coordinates. The ingoing family is $v = $ const and runs horizontally here, since $v$ is itself a null "
        "coordinate. The outgoing family turns vertical at $r_h$ and leans back toward smaller $r$ inside it.",
    ],
    ("boulware_deser", "eddington_finkelstein_outgoing", "finkelstein"): [
        "The plane of $u$ and $r$ ($\\psi = \\theta = \\pi/2$, $\\phi = 0$), drawn for $r_0 = 13r_h/12$ and "
        "$\\ell = 5r_h/12$ with $u + r$ as the vertical axis, so that the outgoing rays, $u = $ const, run at 45°. "
        "The ingoing family stands vertical at $r_h$. Inside $r_h$ both edges of every future cone point to larger "
        "$r$: this is the white hole, which nothing from outside can enter.",
        "The chart is the time reverse of the ingoing one: its $g_{ur}$ is $-1$ where the ingoing chart's "
        "$g_{vr}$ is $+1$.",
    ],
    ("boulware_deser", "eddington_finkelstein_outgoing", "chart"): [
        "The same plane of $u$ and $r$ ($\\psi = \\theta = \\pi/2$, $\\phi = 0$), drawn against the chart's own "
        "coordinates. The outgoing family is $u = $ const and runs horizontally here, since $u$ is itself a null "
        "coordinate. The ingoing family turns vertical at $r_h$ and leans toward larger $r$ inside it.",
    ],
    ("boulware_deser", "spherical_plus", "radial"): [
        "The plane of $t$ and $r$ ($\\psi = \\theta = \\pi/2$, $\\phi = 0$) of the anti-de Sitter branch, drawn "
        "for $r_0 = \\ell$. The rays run at $dr/d(ct) = \\pm f_+$ with $f_+ = 1 + 2r_0^2/(W - r^2)$, which is "
        "$1 + r_0/\\ell = 2$ at the centre and grows as $r^2/\\ell^2$ far away, so the cones open toward the "
        "horizontal and a light ray runs from any radius to $r \\to \\infty$ in a finite time $t$, as in anti-de "
        "Sitter space.",
        "No horizon stands anywhere on the plane: $f_+ > 1$ at every $r$, and $t$ is a time throughout. The "
        "singularity $r = 0$, where the Kretschmann scalar grows as $12r_0^2/(\\ell^2r^4)$, is timelike, and a light "
        "ray from far away reaches it and one from it reaches every radius.",
    ],
    ("black_string", "static", "radial"): [
        "The plane of $t$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$, $z = 0$), the same at every angle and at every "
        "$z$ along the string. It is Schwarzschild's plane: outside $r_s$ the cones narrow toward the vertical as "
        "$r \\to r_s$, because $dt/dr = \\pm(1 - r_s/r)^{-1}$ diverges there. The Gregory-Laflamme instability "
        "is a ripple along $z$, which is fixed on this plane.",
        "Inside $r_s$ the same components make $r$ the time. We take the future from the ingoing "
        "Eddington-Finkelstein chart, which runs smoothly across $r_s$, and this makes that region the black "
        "hole, where every cone points to $r = 0$, a singular line along $z$. The Kretschmann scalar "
        "$12r_s^2/r^6$ is finite at $r_s$ and diverges only at $r = 0$.",
    ],
    ("black_string", "eddington_finkelstein_ingoing", "finkelstein"): [
        "The plane of $v$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$, $z = 0$), drawn with $v - r$ as the vertical "
        "axis so that the ingoing rays, $v = $ const, run at 45°. The outgoing family has "
        "$dv/dr = 2(1 - r_s/r)^{-1}$, so it stands exactly vertical at $r_s$: each generator of the horizon is "
        "an outgoing ray that stays where it is.",
        "The cones cross $r_s$ smoothly and keep tipping. Inside, both edges of every future cone point to "
        "smaller $r$, so every future directed ray on this plane ends at $r = 0$. The apparent horizon of a "
        "rippled string, where the outgoing rays stop diverging, lies at $r_s$ plus a term that varies as a cosine along $z$.",
    ],
    ("black_string", "eddington_finkelstein_ingoing", "chart"): [
        "The same plane of $v$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$, $z = 0$), drawn against the chart's own "
        "coordinates. The ingoing family is $v = $ const and runs horizontally here, since $v$ is itself a null "
        "coordinate. The outgoing family turns vertical at $r_s$ and leans back toward smaller $r$ inside it.",
    ],
    ("black_string", "kerr_schild", "radial"): [
        "The plane of $T$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$, $z = 0$). The ingoing rays are the straight "
        "lines $cT + r = $ const, at 45° as in flat space, and the outgoing family has "
        "$c\\,dT/dr = (r + r_s)/(r - r_s)$, vertical at $r_s$. The surfaces of constant $T$ are spacelike at every "
        "$r$ and cross the horizon, and Choptuik and his collaborators set the initial data of their evolutions "
        "on one.",
        "Inside $r_s$ both edges of every future cone point to smaller $r$, and every future directed ray on "
        "this plane ends at $r = 0$.",
    ],
    ("black_string", "static_six", "radial"): [
        "The plane of $t$ and $r$ ($\\psi = \\theta = \\pi/2$, $\\phi = 0$, $z = 0$) of the black string in six "
        "dimensions, the same at every angle and at every $z$. It is the plane of the Schwarzschild-Tangherlini "
        "black hole of five dimensions: the cones narrow toward the vertical as $r \\to r_h$, because "
        "$dt/dr = \\pm(1 - r_h^2/r^2)^{-1}$ diverges there.",
        "Inside $r_h$ the same components make $r$ the time, and we take the future as the ingoing rays carry it "
        "across the horizon, which makes that region the black hole, where every cone points to $r = 0$. The "
        "Kretschmann scalar $72r_h^4/r^8$ is finite at $r_h$ and diverges only at $r = 0$.",
    ],
    ("tangherlini", "spherical_six", "radial"): [
        "The plane of $t$ and $r$ ($\\chi = \\psi = \\theta = \\pi/2$, $\\phi = 0$) in six dimensions, the same "
        "at every fixed angle by hyperspherical symmetry. Outside $r_h$ the cones narrow toward the vertical as "
        "$r \\to r_h$, because $dt/dr = \\pm(1 - r_h^3/r^3)^{-1}$ diverges there. At $r = 2r_h$ it is $\\pm 8/7$, "
        "where Schwarzschild's at twice its own horizon radius is $\\pm 2$.",
        "Inside $r_h$ the same components make $r$ the time, and we take the future as the ingoing rays carry it "
        "across the horizon, which makes that region the black hole, where every cone points to $r = 0$. The "
        "Kretschmann scalar $240r_h^6/r^{10}$ is finite at $r_h$ and diverges only at $r = 0$.",
    ],
    ("myers_perry", "boyer_lindquist", "transverse"): [
        "The plane of $t$ and $r$ ($\\theta = 0$, $\\psi = 0$) in the plane transverse to the rotation, drawn for "
        "$a = 0.6\\sqrt{\\mu}$, each point in the plane a circle of $\\psi$ of radius $r$. The metric on it is "
        "$-f\\,c^2dt^2 + dr^2/f$ with $f = (r^2 + a^2 - \\mu)/(r^2 + a^2)$, so the rays run at "
        "$dr/d(ct) = \\pm f$ and the cones close at the one horizon, $r_+ = \\sqrt{\\mu - a^2} = 0.8\\sqrt{\\mu}$. "
        "No Christoffel symbol turns the rays out of the plane, so they are null geodesics.",
        "Inside $r_+$, $r$ is the time, and we take the future from the ingoing chart, which makes that region "
        "the black hole. Every ray there ends at $r = 0$, where the circle of $\\psi$ has shrunk to a point. The "
        "Kretschmann scalar is $72\\mu^2/a^8$ there, finite, and the surface $r = 0$ is the locus of a conical "
        "singularity whose rim, in the plane of rotation, is where the curvature diverges.",
    ],
    ("myers_perry", "boyer_lindquist", "rotation"): [
        "The plane of rotation ($\\theta = \\pi/2$) drawn in $t$ and $r$ with $\\phi$ divided out, for "
        "$a = 0.6\\sqrt{\\mu}$: the metric orthogonal to the circles of $\\phi$. Its null curves are the shadows "
        "on $t$ and $r$ of the null geodesics of zero angular momentum, each turning in $\\phi$ at "
        "$d\\phi/d(ct) = \\mu a/(r^4 + a^2r^2 + \\mu a^2)$, and each cone is the future cone of the directions of "
        "zero angular momentum. They run at $dr/d(ct) = \\pm(r^2 + a^2 - \\mu)/\\sqrt{r^4 + a^2r^2 + \\mu a^2}$, "
        "and the cones close at the horizon $r_+ = 0.8\\sqrt{\\mu}$.",
        "The dotted line is the ergosurface, $g_{tt} = 0$ at $r = \\sqrt{\\mu}$, and between it and $r_+$ no "
        "observer keeps $\\phi$ fixed. Inside $r_+$ every ray ends at $r = 0$, the ring where the Kretschmann "
        "scalar $72\\mu^2/r^8$ diverges.",
    ],
    ("myers_perry", "ingoing_kerr", "transverse"): [
        "The plane of $v$ and $r$ ($\\theta = 0$, $\\psi = 0$) in the plane transverse to the rotation, drawn "
        "for $a = 0.6\\sqrt{\\mu}$ with $v - r$ as the vertical axis so that the ingoing rays, $v = $ const, run "
        "at 45°. The outgoing family has $dv/dr = 2(r^2 + a^2)/(r^2 + a^2 - \\mu)$, so it stands vertical at "
        "$r_+ = 0.8\\sqrt{\\mu}$: the horizon is an outgoing ray that stays where it is.",
        "The cones cross $r_+$ smoothly and keep tipping. Inside it both edges of every future cone point to "
        "smaller $r$, so every future directed ray ends at $r = 0$.",
    ],
    ("myers_perry", "ingoing_kerr", "rotation"): [
        "The plane of rotation ($\\theta = \\pi/2$) drawn in $v$ and $r$ with $\\tilde\\phi$ divided out, for "
        "$a = 0.6\\sqrt{\\mu}$, with $v - r$ as the vertical axis. Its null curves are the shadows of the null "
        "geodesics of zero angular momentum, and each cone is the future cone of those directions. The lines of "
        "constant $v$ and $\\tilde\\phi$ are the ingoing principal rays, which carry angular momentum, so the "
        "ingoing family drawn here crosses them.",
        "The outgoing family stands vertical at the horizon $r_+ = 0.8\\sqrt{\\mu}$, and the cones cross it "
        "smoothly. The dotted line is the ergosurface, $g_{vv} = 0$ at $r = \\sqrt{\\mu}$. Inside $r_+$ both "
        "edges of every cone point to smaller $r$, and every ray ends on the ring $r = 0$.",
    ],
    ("myers_perry", "equal_spins", "radial"): [
        "The plane of $t$ and $\\rho$ of the hole with equal spins ($a = 0.4\\sqrt{\\mu}$) with $\\psi$ divided "
        "out, each point in the plane a squashed 3-sphere of $\\theta$, $\\phi$, and $\\psi$. The metric depends "
        "on $\\rho$ alone, so the plane is the same at every $\\theta$ and $\\phi$. Its null curves are the "
        "shadows of the null geodesics with no angular momentum, which run at $d\\rho/d(ct) = "
        "\\pm(\\rho^4 - \\mu\\rho^2 + \\mu a^2)/(\\rho^2\\sqrt{\\rho^4 + \\mu a^2})$.",
        "The cones close at the horizons $\\rho_+ = 0.894\\sqrt{\\mu}$ and $\\rho_- = 0.447\\sqrt{\\mu}$, and "
        "the dotted line is the ergosurface, $g_{tt} = 0$ at $\\rho = \\sqrt{\\mu}$, a round value at every "
        "angle. Inside $\\rho_-$ the lines of constant $\\rho$ are timelike again, down to the singularity "
        "$\\rho = 0$, where the Kretschmann scalar $384\\mu^2a^4/\\rho^{12}$ diverges.",
    ],
    ("myers_perry", "boyer_lindquist_six", "transverse"): [
        "The plane of $t$ and $r$ ($\\theta = 0$, $\\chi = \\pi/2$, $\\psi = 0$) in the space transverse to the "
        "rotation in six dimensions, drawn for $a = 1.5\\,\\mu^{1/3}$, each point in the plane a 2-sphere of "
        "radius $r$. The metric on it is $-f\\,c^2dt^2 + dr^2/f$ with $f = 1 - \\mu/(r(r^2 + a^2))$, and the "
        "cones close at the one horizon, the root $r_+ = 0.413\\,\\mu^{1/3}$ of $r^3 + a^2r = \\mu$. No "
        "Christoffel symbol turns the rays out of the plane, so they are null geodesics.",
        "In five dimensions a hole of this mass parameter loses its horizon at $a = \\sqrt{\\mu}$. In six the "
        "cubic has a positive root for every $a$, and as the spin grows the horizon only moves in, toward "
        "$r_+ = \\mu/a^2$. Inside it every ray ends at the singularity $r = 0$, where the Kretschmann scalar "
        "diverges.",
    ],
    ("myers_perry", "boyer_lindquist_six", "rotation"): [
        "The plane of rotation ($\\theta = \\pi/2$) in six dimensions drawn in $t$ and $r$ with $\\phi$ divided "
        "out, for $a = 1.5\\,\\mu^{1/3}$. Its null curves are the shadows of the null geodesics of zero angular "
        "momentum, and each cone is the future cone of those directions. The cones close at the horizon "
        "$r_+ = 0.413\\,\\mu^{1/3}$.",
        "The dotted line is the ergosurface, $g_{tt} = 0$ at $r = \\mu^{1/3}$. In this plane the horizon is a "
        "circle of circumference radius $(r_+^2 + a^2)/r_+ = 5.86\\,\\mu^{1/3}$, and across the plane it is a "
        "2-sphere of radius $r_+$, fourteen times smaller: the pancake of Emparan and Myers.",
    ],
    ("schwarzschild_ads", "static", "radial"): [
        "The plane of $t$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$), drawn for $r_s = 2L$, the same at every "
        "other angle by spherical symmetry. There $g^{rr} = 1 - r_s/r + r^2/L^2$ vanishes at the horizon "
        "$r_h = L$, and the cones close on it, since $dt/dr = \\pm(1 - r_s/r + r^2/L^2)^{-1}$ diverges there. "
        "Far outside, $g^{rr}$ grows as $r^2/L^2$ and the cones open toward the horizontal: a light ray runs "
        "from any radius to $r \\to \\infty$ in a finite time $t$, as in anti-de Sitter space.",
        "Inside $r_h$, $t$ is a spacelike coordinate, and the components alone do not fix which way is future. "
        "We take it from the ingoing Eddington-Finkelstein chart, which makes that region the black hole, where "
        "every cone points to $r = 0$. The Kretschmann scalar $12r_s^2/r^6 + 24/L^4$ is finite at $r_h$ and "
        "diverges only at $r = 0$.",
    ],
    ("schwarzschild_ads", "eddington_finkelstein_ingoing", "finkelstein"): [
        "The plane of $v$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$), drawn for $r_s = 2L$ with $v - r$ as the "
        "vertical axis, so that the ingoing rays, $v = $ const, run at 45°. The outgoing family has "
        "$dv/dr = 2(1 - r_s/r + r^2/L^2)^{-1}$, so it stands vertical at $r_h = L$: the horizon is an outgoing "
        "ray that stays where it is.",
        "The cones cross $r_h$ smoothly and keep tipping. Inside, both edges of every future cone point to "
        "smaller $r$, so every future directed ray ends at $r = 0$. Far outside, the outgoing edge leans toward "
        "the ingoing one, since $dv/dr$ falls as $2L^2/r^2$.",
    ],
    ("schwarzschild_ads", "eddington_finkelstein_ingoing", "chart"): [
        "The same plane of $v$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$), drawn against the chart's own "
        "coordinates. The ingoing family is $v = $ const and runs horizontally here, since $v$ is itself a "
        "null coordinate. The outgoing family turns vertical at $r_h$ and leans back toward smaller $r$ inside it.",
    ],
    ("schwarzschild_ads", "eddington_finkelstein_outgoing", "finkelstein"): [
        "The plane of $u$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$), drawn for $r_s = 2L$ with $u + r$ as the "
        "vertical axis, so that the outgoing rays, $u = $ const, run at 45°. The ingoing family stands vertical "
        "at $r_h = L$. Inside $r_h$ both edges of every future cone point to larger $r$: this is the white hole, "
        "which nothing from outside can enter.",
        "The chart is the time reverse of the ingoing one: its $g_{ur}$ is $-1$ where the ingoing chart's "
        "$g_{vr}$ is $+1$.",
    ],
    ("schwarzschild_ads", "eddington_finkelstein_outgoing", "chart"): [
        "The same plane of $u$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$), drawn against the chart's own "
        "coordinates. The outgoing family is $u = $ const and runs horizontally here, since $u$ is itself a "
        "null coordinate. The ingoing family turns vertical at $r_h$ and leans toward larger $r$ inside it.",
    ],
    ("bardeen", "static", "radial"): [
        "The plane of $t$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$), drawn for $g = r_s/3$, the same at every "
        "other angle by spherical symmetry. There $g^{rr} = 1 - r_sr^2/(r^2 + g^2)^{3/2}$ vanishes twice, at "
        "$r_+ = 0.775\\,r_s$ and $r_- = 0.301\\,r_s$, and the cones close at both. Between them $r$ is the time "
        "and the cones point to smaller $r$. Inside $r_-$, $t$ is a time again, and the cones open toward 45° "
        "as $r \\to 0$, where $g^{rr} \\approx 1 - r_sr^2/g^3$, as in de Sitter space.",
        "The chart alone does not fix which way is future in the two inner regions. We take it from the ingoing "
        "Eddington-Finkelstein chart, which makes the region between the horizons the black hole. The left "
        "edge, $r = 0$, is a regular centre, where the Kretschmann scalar is $24r_s^2/g^6$.",
    ],
    ("bardeen", "eddington_finkelstein_ingoing", "finkelstein"): [
        "The plane of $v$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$), drawn for $g = r_s/3$ with $v - r$ as the "
        "vertical axis, so that the ingoing rays, $v = $ const, run at 45°. The outgoing family has "
        "$dv/dr = 2(1 - r_sr^2/(r^2 + g^2)^{3/2})^{-1}$, so it stands vertical at both horizons, "
        "$r_+ = 0.775\\,r_s$ and $r_- = 0.301\\,r_s$: each is an outgoing ray that stays where it is.",
        "The cones cross both horizons smoothly. Between them both edges of every future cone point to smaller "
        "$r$, and inside $r_-$ the outgoing edge points outward again, so an outgoing ray there climbs toward "
        "$r_-$ and never reaches it. An ingoing ray reaches the regular centre $r = 0$ at a finite $v$ and "
        "passes through it, to continue as an outgoing ray.",
    ],
    ("bardeen", "eddington_finkelstein_ingoing", "chart"): [
        "The same plane of $v$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$), drawn against the chart's own "
        "coordinates. The ingoing family is $v = $ const and runs horizontally here, since $v$ is itself a "
        "null coordinate. The outgoing family turns vertical at $r_+$, leans back toward smaller $r$ between the "
        "horizons, and turns vertical again at $r_-$.",
    ],
    ("bardeen", "eddington_finkelstein_outgoing", "finkelstein"): [
        "The plane of $u$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$), drawn for $g = r_s/3$ with $u + r$ as the "
        "vertical axis, so that the outgoing rays, $u = $ const, run at 45°. The ingoing family stands vertical "
        "at both horizons, $r_+ = 0.775\\,r_s$ and $r_- = 0.301\\,r_s$. Between them both edges of every "
        "future cone point to larger $r$: this is the white hole, which nothing from outside can enter.",
        "The chart is the time reverse of the ingoing one: its $g_{ur}$ is $-1$ where the ingoing chart's "
        "$g_{vr}$ is $+1$. An outgoing ray leaves the regular centre $r = 0$ and crosses both horizons on its "
        "way to infinity.",
    ],
    ("bardeen", "eddington_finkelstein_outgoing", "chart"): [
        "The same plane of $u$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$), drawn against the chart's own "
        "coordinates. The outgoing family is $u = $ const and runs horizontally here, since $u$ is itself a "
        "null coordinate. The ingoing family turns vertical at $r_-$ and at $r_+$, and leans toward larger $r$ "
        "between them.",
    ],
    ("kaluza_klein_monopole", "gross_perry", "radial"): [
        "The plane of $t$ and $r$ ($\\theta = 0$, $x_5 = 0$) of the Kaluza-Klein monopole, each point in the "
        "plane a squashed 3-sphere of $\\theta$, $\\phi$, and $x_5$. The metric on it is "
        "$-c^2dt^2 + (1 + 4m/r)\\,dr^2$, so the rays run at $dr/d(ct) = \\pm\\sqrt{r/(r + 4m)}$: at 45° far "
        "away, and slower in $r$ toward the nut $r = 0$, where the cones close. No Christoffel symbol turns "
        "them out of the plane, so they are null geodesics.",
        "The cones close because $r$ is a poor ruler near the nut: the proper distance from it is "
        "$\\int\\sqrt{1 + 4m/r}\\,dr \\approx 4\\sqrt{mr}$, and against that distance every ray runs at the "
        "speed of light. A ray reaches $r = 0$ in a finite time, and the Kretschmann scalar "
        "$384m^2/(r + 4m)^6$ is finite there. The monopole has no horizon.",
    ],
    ("kaluza_klein_monopole", "gross_perry", "through"): [
        "The line through the nut along the axis: $x = r$ on the right is $\\theta = 0$ and $x = -r$ on the left "
        "is $\\theta = \\pi$, and the monopole's spherical symmetry makes the two halves mirror images. Rays "
        "cross the nut smoothly, since $r = 0$ is a regular point of the geometry in five dimensions, where the "
        "Kretschmann scalar is $3/(32m^4)$.",
    ],
    ("kaluza_klein_monopole", "hopf", "radial"): [
        "The plane of $t$ and $r$ ($\\theta = 0$, $\\psi = 0$) of the Kaluza-Klein monopole in the chart of the "
        "Hopf angle, each point in the plane a squashed 3-sphere of the Euler angles $\\theta$, $\\phi$, and "
        "$\\psi$. The metric on it is $-c^2dt^2 + (1 + 4m/r)\\,dr^2$, as in Gross and Perry's chart, so the rays "
        "run at $dr/d(ct) = \\pm\\sqrt{r/(r + 4m)}$ and close toward the nut $r = 0$, which they reach in a finite "
        "time. No Christoffel symbol turns them out of the plane, so they are null geodesics.",
    ],
    ("kaluza_klein_monopole", "taub_nut", "radial"): [
        "The plane of $t$ and $\\rho$ ($\\theta = 0$, $\\psi = 0$) of the Kaluza-Klein monopole with the radius "
        "of the Taub-NUT line element, each point in the plane a squashed 3-sphere of the Euler angles. The "
        "metric on it is $-c^2dt^2 + (\\rho + 2m)/(\\rho - 2m)\\,d\\rho^2$, so the rays run at "
        "$d\\rho/d(ct) = \\pm\\sqrt{(\\rho - 2m)/(\\rho + 2m)}$ and the cones close at $\\rho = 2m$.",
        "The edge $\\rho = 2m$ is the nut, a single point of space where the 3-spheres have shrunk away, and the "
        "chart has no points with $\\rho < 2m$. A ray reaches it in a finite time, and the Kretschmann scalar "
        "$384m^2/(\\rho + 2m)^6$ is finite there.",
    ],
    ("kaluza_klein_black_hole", "electric", "radial"): [
        "The plane of $t$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$) of the electrically charged hole with the "
        "circle of $y$ divided out, drawn for $q = 2\\,r_s$: the metric orthogonal to the circles. Its null "
        "curves are the shadows on $t$ and $r$ of the null geodesics with no momentum along $y$, the light "
        "that carries no electric charge in four dimensions, and each cone is the future cone of those "
        "directions. They run at $c\\,dt/dr = \\pm\\sqrt{1 + (q - r_s)/r}\\,(1 - r_s/r)^{-1}$, and the cones "
        "close at the horizon $r_s$.",
        "The dotted line is the ergosurface, $g_{tt} = 0$ at $r = q$. The black string moves along $y$, and "
        "between $r = q$ and $r_s$ no observer keeps $y$ fixed. Inside $r_s$ we take the future from the "
        "ingoing chart, and every ray ends at $r = 0$, where the Kretschmann scalar $12r_s^2/r^6$ diverges.",
    ],
    ("kaluza_klein_black_hole", "eddington_finkelstein_ingoing", "finkelstein"): [
        "The plane of $v$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$) with the circle of $w$ divided out, drawn "
        "for $q = 2\\,r_s$ with $v - \\sqrt{q/r_s}\\,r$ as the vertical axis. The rays with no momentum along "
        "the circle run at $dv/dr = (\\sqrt{q/r_s} \\pm \\sqrt{1 + (q - r_s)/r})(1 - r_s/r)^{-1}$: the ingoing "
        "family crosses the horizon $r_s$ with a finite slope, and the outgoing family stands vertical there, "
        "an outgoing ray that stays where it is.",
        "Inside $r_s$ both edges of every future cone point to smaller $r$, so every such ray ends at $r = 0$. "
        "The dotted line is the ergosurface $r = q$. The curves of constant $v$ and $w$ are the ingoing rays "
        "of the string's rest frame, and in this chart they carry momentum along the circle.",
    ],
    ("kaluza_klein_black_hole", "magnetic", "radial"): [
        "The plane of $t$ and $r$ ($\\theta = 0$, $y = 0$) of the magnetically charged hole, drawn for "
        "$p = 2\\,r_s$ on the half axis where the potential vanishes, each point in the plane a squashed "
        "3-sphere of $\\theta$, $\\phi$, and $y$. The metric on it is "
        "$-(1 - r_s/r)\\,c^2dt^2 + (1 + (p - r_s)/r)(1 - r_s/r)^{-1}dr^2$, so the rays run at "
        "$c\\,dt/dr = \\pm\\sqrt{1 + (p - r_s)/r}\\,(1 - r_s/r)^{-1}$ and the cones close at the horizon $r_s$. "
        "No Christoffel symbol turns them out of the plane, so they are null geodesics.",
        "Inside $r_s$ we take the future from the ingoing rays, and every ray ends at $r = 0$, where the "
        "circle of $y$ has shrunk to a point and the Kretschmann scalar diverges. At $r_s = 0$ that point is "
        "the nut of the Kaluza-Klein monopole, a regular point.",
    ],
    ("kaluza_klein_black_hole", "dyonic", "radial"): [
        "The plane of $t$ and $\\rho$ ($\\theta = \\pi/2$, $\\phi = 0$) of the hole of equal charges with the "
        "circle of $y$ divided out, drawn for $p = 3\\,r_s$: the metric orthogonal to the circles, which is "
        "Reissner and Nordström's. $g^{\\rho\\rho}$ vanishes twice, at $\\rho_\\pm = (p \\pm r_s)/2$, which is "
        "$2\\,r_s$ and $r_s$, and the cones close at both. Between them $\\rho$ is the time and the cones "
        "point to smaller $\\rho$. Inside $\\rho_-$, $t$ is a time again.",
        "We take the future in the two inner regions from an ingoing chart, which makes the region between "
        "the horizons the black hole. The Kretschmann scalar diverges at $\\rho = 0$, a timelike singularity.",
    ],
    ("kaluza_klein_black_hole", "einstein", "radial"): [
        "The plane of $t$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$) of the Einstein metric of four dimensions, "
        "drawn for $q = 2\\,r_s$, each point in the plane a 2-sphere of area $4\\pi r\\sqrt{r(r + q - r_s)}$. "
        "The rays run at $c\\,dt/dr = \\pm\\sqrt{r(r + q - r_s)}/(r - r_s)$, as the rays with no momentum "
        "along the circle do in five dimensions, and the cones close at the horizon $r_s$.",
        "Inside $r_s$ we take the future from the ingoing Eddington-Finkelstein chart, which makes that region "
        "the black hole, where every cone points to the singularity $r = 0$. The spheres have zero area there "
        "and the Kretschmann scalar diverges.",
    ],
    ("kaluza_klein_black_hole", "einstein_eddington_finkelstein", "finkelstein"): [
        "The plane of $v$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$) of the Einstein metric of four dimensions, "
        "drawn for $q = 2\\,r_s$ with $v - r$ as the vertical axis, so that the ingoing rays, $v = $ const, run "
        "at 45°. The outgoing family has $dv/dr = 2\\sqrt{r(r + q - r_s)}/(r - r_s)$, so it stands vertical at "
        "the horizon $r_s$, an outgoing ray that stays where it is.",
        "The cones cross $r_s$ smoothly and keep tipping. Inside, both edges of every future cone point to "
        "smaller $r$, so every future directed ray ends at the singularity $r = 0$, where the area of the "
        "spheres vanishes.",
    ],
    ("hayward", "static", "radial"): [
        "The plane of $t$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$), drawn for $\\ell = 12m/7\\sqrt{7} = 0.648\\,m$, "
        "the same at every other angle by spherical symmetry. There $g^{rr} = 1 - 2mr^2/(r^3 + 2m\\ell^2)$ "
        "vanishes at $r_- = 6m/7$ and at $r_+ = 12m/7$, and the cones close on both horizons, since "
        "$c\\,dt/dr = \\pm 1/g^{rr}$ diverges there. Inside $r_-$ the cones open again, and at the centre, where "
        "$g^{rr} = 1$, they stand at 45°.",
        "Between the horizons $t$ is a spacelike coordinate, and the components alone do not fix which way is "
        "future. We take it from the ingoing Eddington-Finkelstein chart, which makes that region the black "
        "hole, where every cone points to smaller $r$. The Kretschmann scalar is finite at every radius, and at "
        "$r = 0$ it is $24/\\ell^4$, its value in de Sitter space of radius $\\ell$.",
    ],
    ("hayward", "eddington_finkelstein_ingoing", "finkelstein"): [
        "The plane of $v$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$), drawn for $\\ell = 0.648\\,m$ with $v - r$ as "
        "the vertical axis, so that the ingoing rays, $v = $ const, run at 45°. The outgoing family has "
        "$dv/dr = 2/g^{rr}$, so it stands vertical at $r_- = 6m/7$ and at $r_+ = 12m/7$: each horizon is an "
        "outgoing ray that stays where it is.",
        "The cones cross $r_+$ smoothly and keep tipping. Between the horizons both edges of every future cone "
        "point to smaller $r$. Inside $r_-$ the outgoing edge points to larger $r$ again, and an outgoing ray "
        "there climbs toward $r_-$, which it approaches as $v \\to \\infty$. An ingoing ray reaches the centre at "
        "a finite $v$, where the curvature is finite.",
    ],
    ("hayward", "eddington_finkelstein_ingoing", "chart"): [
        "The same plane of $v$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$), drawn against the chart's own "
        "coordinates. The ingoing family is $v = $ const and runs horizontally here, since $v$ is itself a "
        "null coordinate. The outgoing family turns vertical at $r_+$, leans toward smaller $r$ between the "
        "horizons, and turns vertical again at $r_-$.",
    ],
    ("hayward", "eddington_finkelstein_outgoing", "finkelstein"): [
        "The plane of $u$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$), drawn for $\\ell = 0.648\\,m$ with $u + r$ as "
        "the vertical axis, so that the outgoing rays, $u = $ const, run at 45°. The ingoing family stands "
        "vertical at $r_- = 6m/7$ and at $r_+ = 12m/7$. Between the horizons both edges of every future cone "
        "point to larger $r$: this is the white hole, which nothing from outside $r_+$ can enter.",
        "The chart is the time reverse of the ingoing one: its $g_{ur}$ is $-1$ where the ingoing chart's "
        "$g_{vr}$ is $+1$.",
    ],
    ("hayward", "eddington_finkelstein_outgoing", "chart"): [
        "The same plane of $u$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$), drawn against the chart's own "
        "coordinates. The outgoing family is $u = $ const and runs horizontally here, since $u$ is itself a "
        "null coordinate. The ingoing family turns vertical at $r_+$, leans toward larger $r$ between the "
        "horizons, and turns vertical again at $r_-$.",
    ],
    ("hayward", "evaporating", "history"): [
        "The plane of $v$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$) of a black hole that forms and evaporates "
        "($\\ell = 0.648\\,m_0$), drawn with $v - r$ as the vertical axis, each point in the plane a 2-sphere "
        "of area $4\\pi r^2$. Before $v = 0$ and after $v = 8\\,m_0$ the mass vanishes and the plane is "
        "Minkowski's. The curve $g^{rr} = 0$ is closed. It opens at $r = \\sqrt{3}\\,\\ell = 1.12\\,m_0$ when the "
        "mass passes $3\\sqrt{3}\\,\\ell/4 = 0.842\\,m_0$, at $v = 1.48\\,m_0$, reaches $r_- = 6m_0/7$ and "
        "$r_+ = 12m_0/7$ while the mass is $m_0$, and closes at $v = 5.04\\,m_0$.",
        "Inside the curve both edges of every future cone point to smaller $r$, so every sphere there is "
        "trapped. Its outer part is the outer trapping horizon and its inner part the inner one. An outgoing "
        "ray inside the curve loses $r$, and gains it again once the curve has closed or its inner part has "
        "swept past the ray. Each such ray reaches infinity, so this spacetime has no event horizon.",
    ],
    ("hiscock", "ingoing", "history"): [
        "The plane of $v$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$) of a black hole that forms and evaporates, drawn "
        "with $v - r$ as the vertical axis so that the ingoing rays, $v$ constant, run at 45°, each point in the "
        "plane a 2-sphere of area $4\\pi r^2$. A shell of light of mass $m_0$ falls in along $v = 0$ and makes the "
        "hole, and from $v_1 = 2\\,m_0$ to $v_0 = 8\\,m_0$ negative energy falls in and the mass drops to zero. "
        "The chart ends on the surface of pair creation, $R = 3m$, which shrinks with the mass onto the centre at "
        "$v_0$.",
        "On the apparent horizon $r = 2m$, where $g^{rr}$ vanishes, the outgoing rays stand still, and inside it "
        "both edges of every future cone point to smaller $r$. While the mass falls it moves inward as a timelike "
        "curve. The event horizon is the outgoing ray that reaches the centre just as the mass runs out: it leaves "
        "the centre at $v = -3.41\\,m_0$, before the shell arrives, and stays inside the apparent horizon to "
        "the end, so light that starts between the two escapes. The singularity $r = 0$, where the Kretschmann "
        "scalar $48m^2/r^6$ diverges, lasts from $v = 0$ to $v_0$.",
    ],
    ("hiscock", "ingoing", "shells"): [
        "The plane of $v$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$) of the simplest of the models, drawn with "
        "$v - r$ as the vertical axis. A shell of light of mass $m_0$ falls in along $v = 0$ and a shell of mass "
        "$-m_0$ follows along $v_1 = 6\\,m_0$, with flat space before the first, Schwarzschild's metric between "
        "them, and flat space again after the second.",
        "The event horizon is the outgoing ray that reaches the centre as the second shell does, and the "
        "singularity ends there. Inside $r = 2m_0$ an outgoing ray loses $r$, and one that is still above the "
        "event horizon when the second shell passes escapes into the flat space beyond. The outgoing ray from "
        "the last point of the singularity is the Cauchy horizon, along which Hiscock found the energy flux of a "
        "quantum field diverging.",
    ],
    ("hiscock", "outgoing", "history"): [
        "The plane of $u$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$) outside the surface of pair creation, drawn "
        "with $u + r$ as the vertical axis so that the outgoing rays, $u$ constant, run at 45°. Each outgoing ray "
        "carries the mass $m(u)$ it left the surface with, $m_0$ until $u_1 = -16\\,m_0$ and zero from "
        "$u_0 = 8\\,m_0$, and the radiation between them has positive energy. The chart begins on the collapsing "
        "shell before $u_1$ and on the surface $R = 3m$ after it.",
        "The ray $u_0$ is the Cauchy horizon, the first light from the last point of the singularity, and beyond "
        "it the plane is Minkowski's, with a regular centre. Everywhere on the plane $g^{rr} = 1 - 2m/r$ is "
        "positive, so no sphere here is trapped. A distant observer reads the mass off the rays as they arrive "
        "and sees it fall to zero with $dm/du \\to 0$.",
    ],
    ("hiscock", "flat", "after"): [
        "The plane of $u$ and $v$ ($\\theta = \\pi/2$, $\\phi = 0$) after the hole has gone, drawn with "
        "$r = (v - u)/2$ across and $(u + v)/2$ up, so that every ray runs at 45°. The chart begins on the Cauchy "
        "horizon $u_0 = 8\\,m_0$, the outgoing ray from the last point of the singularity, and its left edge is "
        "the centre $v = u$, a regular centre once more.",
    ],
    **{("topological_black_hole", system, "flat"): [
        "The plane of $t$ and $r$ (" + at + ") of the flat hole ($k = 0$, $\\mu = L$), the same at every "
        "other point of the horizon by symmetry. The cones close at the horizon $r_h = (\\mu L^2)^{1/3} = L$, "
        "where $g^{rr} = r^2/L^2 - \\mu/r$ vanishes and $c\\,dt/dr = \\pm(r^2/L^2 - \\mu/r)^{-1}$ diverges. Far "
        "out the rays flatten, and a ray reaches $r \\to \\infty$, the conformal boundary, in a finite time, "
        "as in anti-de Sitter space.",
        "Inside $r_h$, $r$ is the time. We take the future from the ingoing Eddington-Finkelstein chart, "
        "which makes that region the black hole, every cone pointing to $r = 0$. The Kretschmann scalar "
        "$12\\mu^2/r^6 + 24/L^4$ is finite at $r_h$ and diverges only at $r = 0$.",
    ] for system, at in (("static", "$\\rho = 1/2$, $\\phi = 0$"),)},
    **{("topological_black_hole", system, "massless"): [
        "The plane of $t$ and $r$ (" + at + ") of the hyperbolic hole without mass ($k = -1$, $\\mu = 0$). "
        "The cones close at the horizon $r_h = L$, where $g^{rr} = r^2/L^2 - 1$ vanishes, and the rays are "
        "those of the BTZ hole at $M = 1$, $c\\,dt/dr = \\pm L^2/(r^2 - L^2)$.",
        "Inside $r_h$, $r$ is the time, and we take the future from the ingoing Eddington-Finkelstein "
        "chart, every cone pointing to $r = 0$. The Kretschmann scalar is $24/L^4$ at every point, that of "
        "anti-de Sitter space, and at $r = 0$ the surfaces of constant $t$ and $r$ shrink to zero area.",
    ] for system, at in (("static", "$\\rho = 1/2$, $\\phi = 0$"), ("hyperbolic", "$\\theta = 1$, $\\phi = 0$"))},
    **{("topological_black_hole", system, "negative"): [
        "The plane of $t$ and $r$ (" + at + ") of the hyperbolic hole of negative mass ($k = -1$, "
        "$\\mu = -120L/343$), where $g^{rr} = (r - 5L/7)(r - 3L/7)(r + 8L/7)/(L^2r)$. The cones close at the "
        "event horizon $r_+ = 5L/7$ and again at the inner horizon $r_- = 3L/7$.",
        "Between the horizons $r$ is the time, and we take the future from the ingoing Eddington-Finkelstein "
        "chart, every cone pointing to smaller $r$. Inside $r_-$, $t$ is a time again, so the singularity "
        "$r = 0$, where the Kretschmann scalar $12\\mu^2/r^6 + 24/L^4$ diverges, is timelike, as in the "
        "Reissner-Nordström hole.",
    ] for system, at in (("static", "$\\rho = 1/2$, $\\phi = 0$"), ("hyperbolic", "$\\theta = 1$, $\\phi = 0$"))},
    ("topological_black_hole", "eddington_finkelstein_ingoing", "flat"): [
        "The plane of $v$ and $r$ ($\\rho = 1/2$, $\\phi = 0$) of the flat hole ($k = 0$, $\\mu = L$), with "
        "$v - r$ as the vertical axis, so that the ingoing rays, $v = $ const, run at 45°. The outgoing family "
        "has $dv/dr = 2(r^2/L^2 - \\mu/r)^{-1}$, so it stands vertical at $r_h = L$: the horizon is an "
        "outgoing ray that stays where it is.",
        "The cones cross $r_h$ smoothly and keep tipping. Inside, both edges of every future cone point to "
        "smaller $r$, so every future directed ray ends at $r = 0$. Far outside, the outgoing edge leans toward "
        "the ingoing one, since $dv/dr$ falls as $2L^2/r^2$.",
    ],
    ("topological_black_hole", "eddington_finkelstein_ingoing", "massless"): [
        "The plane of $v$ and $r$ ($\\rho = 1/2$, $\\phi = 0$) of the hyperbolic hole without mass ($k = -1$, "
        "$\\mu = 0$), with $v - r$ as the vertical axis, so that the ingoing rays, $v = $ const, run at 45°. "
        "The outgoing family, $dv/dr = 2L^2/(r^2 - L^2)$, stands vertical at the horizon $r_h = L$.",
        "The cones cross $r_h$ smoothly, and inside it both edges of every future cone point to smaller $r$. "
        "The curvature is that of anti-de Sitter space at every point of the plane, $r = 0$ included.",
    ],
    ("topological_black_hole", "eddington_finkelstein_ingoing", "negative"): [
        "The plane of $v$ and $r$ ($\\rho = 1/2$, $\\phi = 0$) of the hyperbolic hole of negative mass "
        "($k = -1$, $\\mu = -120L/343$), with $v - r$ as the vertical axis, so that the ingoing rays, "
        "$v = $ const, run at 45°. The outgoing family stands vertical at both horizons, $r_+ = 5L/7$ and "
        "$r_- = 3L/7$.",
        "Between the horizons both edges of every future cone point to smaller $r$. Inside $r_-$ the "
        "outgoing edge points to larger $r$ again, so a ray there need not reach the singularity $r = 0$.",
    ],
    ("topological_black_hole", "eddington_finkelstein_outgoing", "flat"): [
        "The plane of $u$ and $r$ ($\\rho = 1/2$, $\\phi = 0$) of the flat hole ($k = 0$, $\\mu = L$), with "
        "$u + r$ as the vertical axis, so that the outgoing rays, $u = $ const, run at 45°. The ingoing family "
        "stands vertical at $r_h = L$. Inside $r_h$ both edges of every future cone point to larger $r$: this "
        "is the white hole, which nothing from outside can enter.",
        "The chart is the time reverse of the ingoing one: its $g_{ur}$ is $-1$ where the ingoing chart's "
        "$g_{vr}$ is $+1$.",
    ],
    ("topological_black_hole", "eddington_finkelstein_outgoing", "massless"): [
        "The plane of $u$ and $r$ ($\\rho = 1/2$, $\\phi = 0$) of the hyperbolic hole without mass ($k = -1$, "
        "$\\mu = 0$), with $u + r$ as the vertical axis, so that the outgoing rays, $u = $ const, run at 45°. "
        "The ingoing family, $du/dr = -2L^2/(r^2 - L^2)$, stands vertical at the horizon $r_h = L$, and inside "
        "it both edges of every future cone point to larger $r$.",
    ],
    ("topological_black_hole", "eddington_finkelstein_outgoing", "negative"): [
        "The plane of $u$ and $r$ ($\\rho = 1/2$, $\\phi = 0$) of the hyperbolic hole of negative mass "
        "($k = -1$, $\\mu = -120L/343$), with $u + r$ as the vertical axis, so that the outgoing rays, "
        "$u = $ const, run at 45°. The ingoing family stands vertical at both horizons, $r_+ = 5L/7$ and "
        "$r_- = 3L/7$.",
        "Between the horizons both edges of every future cone point to larger $r$, the white hole. Inside "
        "$r_-$ the ingoing edge points to smaller $r$ again, toward the timelike singularity $r = 0$.",
    ],
    ("topological_black_hole", "black_string", "radial"): [
        "The plane of $t$ and $r$ ($\\phi = 0$, $z = 0$) of the black string ($\\mu = L$), the same at every "
        "$\\phi$ and $z$ by the symmetry of the cylinder. The cones close at the horizon "
        "$r_h = (\\mu L^2)^{1/3} = L$, where $g^{rr} = r^2/L^2 - \\mu/r$ vanishes, and far out a ray reaches "
        "$r \\to \\infty$ in a finite time $t$.",
        "Inside $r_h$, $r$ is the time. We take the future from the ingoing Eddington-Finkelstein chart, "
        "which makes that region the black hole, every cone pointing to the axis $r = 0$, where the "
        "Kretschmann scalar $12\\mu^2/r^6 + 24/L^4$ diverges.",
    ],
    ("topological_black_hole", "brane", "tz"): [
        "The plane of $t$ and $z$ ($x = 0$, $y = 0$) of the black brane ($z_h = L$), the same at every $x$ "
        "and $y$. The left edge $z = 0$ is the conformal boundary, where the rays run at 45° as in the "
        "Poincaré patch of anti-de Sitter space. The cones close at the horizon $z_h$, where "
        "$g^{zz} = (z^2/L^2)(1 - z^3/z_h^3)$ vanishes and $c\\,dt/dz = \\pm(1 - z^3/z_h^3)^{-1}$ diverges.",
        "Beyond $z_h$, $z$ is the time. We take the future toward larger $z$, which makes that region the "
        "black hole. The Kretschmann scalar $(12/L^4)(2 + z^6/z_h^6)$ is finite at $z_h$ and grows without "
        "bound as $z \\to \\infty$.",
    ],
    ("global_monopole", "static", "radial"): [
        "The plane of $t$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$) of Letelier's black hole in a cloud of strings, "
        "drawn for $\\Delta = 0.19$, each point in the plane a 2-sphere of area $4\\pi r^2$. The cones close at the "
        "horizon $r_h = r_s/(1 - \\Delta) = 1.235\\,r_s$, where $c\\,dt/dr = \\pm(1 - \\Delta - r_s/r)^{-1}$ diverges, "
        "and far from it their edges tend to $c\\,dt/dr = \\pm 1/(1 - \\Delta) = \\pm 1.235$, steeper than 45° at "
        "every radius.",
        "Inside $r_h$, $t$ is a spacelike coordinate, and we take the future from the ingoing Eddington-Finkelstein "
        "chart, which makes that region the black hole, where every cone points to $r = 0$. The Kretschmann scalar "
        "$(12r_s^2 + 8\\Delta\\,r_s\\,r + 4\\Delta^2r^2)/r^6$ is finite at $r_h$ and diverges only at $r = 0$.",
    ],
    ("global_monopole", "conical", "radial"): [
        "The plane of $t$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$) around a global monopole with no mass at its "
        "centre, drawn for $\\Delta = 0.19$, each point in the plane a 2-sphere of area $4\\pi(1 - \\Delta)r^2$. The "
        "metric on the plane is $-c^2dt^2 + dr^2$, so every ray runs at 45°, and $\\Delta$ enters only "
        "$g_{\\theta\\theta}$ and $g_{\\phi\\phi}$. The Kretschmann scalar "
        "$4\\Delta^2/\\left((1 - \\Delta)^2r^4\\right)$ diverges at the centre, $r = 0$.",
    ],
    ("global_monopole", "eddington_finkelstein_ingoing", "finkelstein"): [
        "The plane of $v$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$), drawn for $\\Delta = 0.19$ with $v - r$ as "
        "the vertical axis, so that the ingoing rays, $v = $ const, run at 45°. The outgoing family has "
        "$dv/dr = 2(1 - \\Delta - r_s/r)^{-1}$, so it stands vertical at the horizon $r_h = 1.235\\,r_s$, "
        "an outgoing ray that stays where it is.",
        "The cones cross $r_h$ smoothly and keep tipping. Inside, both edges of every future cone point to "
        "smaller $r$, so every future directed ray ends at $r = 0$.",
    ],
    ("global_monopole", "eddington_finkelstein_ingoing", "chart"): [
        "The same plane of $v$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$), drawn against the chart's own "
        "coordinates. The ingoing family is $v = $ const and runs horizontally here, since $v$ is itself a null "
        "coordinate. The outgoing family turns vertical at $r_h$ and leans back toward smaller $r$ inside it.",
    ],
    ("global_monopole", "eddington_finkelstein_outgoing", "finkelstein"): [
        "The plane of $u$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$), drawn for $\\Delta = 0.19$ with $u + r$ as "
        "the vertical axis, so that the outgoing rays, $u = $ const, run at 45°. The retarded chart crosses the "
        "other horizon. Inside $r_h$ both edges of every future cone point to larger $r$: this is the white hole, "
        "which nothing from outside can enter.",
        "The chart is the time reverse of the ingoing one: its $g_{ur}$ is $-1$ where the ingoing chart's "
        "$g_{vr}$ is $+1$.",
    ],
    ("global_monopole", "eddington_finkelstein_outgoing", "chart"): [
        "The same plane of $u$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$), drawn against the chart's own "
        "coordinates. The outgoing family is $u = $ const and runs horizontally here, since $u$ is itself a null "
        "coordinate. The ingoing family turns vertical at $r_h$ and leans toward larger $r$ inside it.",
    ],
    ("string_black_hole", "static", "radial"): [
        "The plane of $t$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$) of a black hole threaded by a cosmic string, "
        "drawn for $b = 0.9$, each point in the plane a sphere of radius $r$ with a wedge of $36°$ missing around "
        "the string. The string enters $g_{\\phi\\phi}$ alone, so the rays are Schwarzschild's: the cones narrow "
        "toward the vertical as $r \\to r_s$, where $c\\,dt/dr = \\pm(1 - r_s/r)^{-1}$ diverges.",
        "Inside $r_s$, $t$ is a spacelike coordinate, and we take the future from the ingoing Eddington-Finkelstein "
        "chart, which makes that region the black hole, where every cone points to $r = 0$. The Kretschmann scalar "
        "$12r_s^2/r^6$ is finite at $r_s$ and diverges only at $r = 0$.",
    ],
    ("string_black_hole", "wedge", "radial"): [
        "The plane of $t$ and $r$ ($\\theta = \\pi/2$, $\\tilde\\phi = 0$) of a black hole threaded by a cosmic "
        "string, drawn for $b = 0.9$, each point in the plane a sphere of radius $r$ whose angle $\\tilde\\phi$ runs "
        "over $2\\pi b$, $324°$. The line element is Schwarzschild's, and so are the rays: the cones narrow toward "
        "the vertical as $r \\to r_s$, where $c\\,dt/dr = \\pm(1 - r_s/r)^{-1}$ diverges.",
        "Inside $r_s$, $t$ is a spacelike coordinate, and we take the future from the ingoing Eddington-Finkelstein "
        "chart, which makes that region the black hole, where every cone points to $r = 0$.",
    ],
    ("string_black_hole", "eddington_finkelstein_ingoing", "finkelstein"): [
        "The plane of $v$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$), drawn for $b = 0.9$ with $v - r$ as the "
        "vertical axis, so that the ingoing rays, $v = $ const, run at 45°. The outgoing family has "
        "$dv/dr = 2(1 - r_s/r)^{-1}$, so it stands vertical at the horizon $r_s$, an outgoing ray that stays where "
        "it is.",
        "The cones cross $r_s$ smoothly and keep tipping. Inside, both edges of every future cone point to "
        "smaller $r$, so every future directed ray ends at $r = 0$, and the string ends there with them.",
    ],
    ("string_black_hole", "eddington_finkelstein_ingoing", "chart"): [
        "The same plane of $v$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$), drawn against the chart's own "
        "coordinates. The ingoing family is $v = $ const and runs horizontally here, since $v$ is itself a null "
        "coordinate. The outgoing family turns vertical at $r_s$ and leans back toward smaller $r$ inside it.",
    ],
    ("string_black_hole", "eddington_finkelstein_outgoing", "finkelstein"): [
        "The plane of $u$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$), drawn for $b = 0.9$ with $u + r$ as the "
        "vertical axis, so that the outgoing rays, $u = $ const, run at 45°. The retarded chart crosses the "
        "other horizon. Inside $r_s$ both edges of every future cone point to larger $r$: this is the white hole, "
        "which nothing from outside can enter.",
        "The chart is the time reverse of the ingoing one: its $g_{ur}$ is $-1$ where the ingoing chart's "
        "$g_{vr}$ is $+1$.",
    ],
    ("string_black_hole", "eddington_finkelstein_outgoing", "chart"): [
        "The same plane of $u$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$), drawn against the chart's own "
        "coordinates. The outgoing family is $u = $ const and runs horizontally here, since $u$ is itself a null "
        "coordinate. The ingoing family turns vertical at $r_s$ and leans toward larger $r$ inside it.",
    ],
    ("mcvittie", "isotropic", "radial"): [
        "The plane of $t$ and the comoving $r$ ($\\theta = \\pi/2$, $\\phi = 0$), the same at every other "
        "angle by spherical symmetry, for a universe of dust and a cosmological constant with "
        "$\\Lambda r_s^2 = 1/5$. The rays obey $dr/d(ct) = \\pm(1 - \\mu)/(a(1 + \\mu)^3)$ with "
        "$\\mu = r_s/4ar$, so far from the mass the cones are those of the spatially flat "
        "Friedmann-Lemaître-Robertson-Walker universe, and they close on the curve $r = r_s/4a$, where "
        "$\\mu = 1$.",
        "That curve is the sphere of areal radius $R = r_s$, a curvature singularity, and it falls toward "
        "$r = 0$ as $a$ grows. The dotted curve, where $|\\nabla R|^2 = 1 - r_s/R - H^2R^2/c^2$ vanishes "
        "for the areal radius $R = ar(1 + \\mu)^2$, appears at $ct = 2.10\\,r_s$, and its two branches "
        "approach the spheres $R = 1.085\\,r_s$ and $R = 3.215\\,r_s$.",
    ],
    ("mcvittie", "areal", "radial"): [
        "The plane of $t$ and the areal radius $R$ ($\\theta = \\pi/2$, $\\phi = 0$), the same at every "
        "other angle by spherical symmetry, for a universe of dust and a cosmological constant with "
        "$\\Lambda r_s^2 = 1/5$. The rays obey "
        "$dR/d(ct) = \\sqrt{1 - r_s/R}\\,(HR/c \\pm \\sqrt{1 - r_s/R})$, so every outgoing ray gains $R$, "
        "and an ingoing ray loses $R$ only where $g^{RR} = 1 - r_s/R - H^2R^2/c^2$ is positive. The curve "
        "where $g^{RR}$ vanishes appears at $ct = 2.10\\,r_s$ on $R = 3r_s/2$, and its two branches run "
        "toward $1.085\\,r_s$ and $3.215\\,r_s$, the horizons of the Schwarzschild-de Sitter black hole "
        "with the same $r_s$ and $\\Lambda$.",
        "On $R = r_s$ the Ricci scalar $12H^2/c^2 + 6\\dot{H}/(c\\sqrt{1 - r_s/R})$ diverges, and both "
        "families of rays leave that sphere, which lies in the past of every event of the plane. The "
        "curvature also diverges toward $t = 0$, where $H$ is infinite. An ingoing ray that starts between "
        "the two branches runs down to the inner one as $t \\to \\infty$ and reaches it at a finite affine "
        "parameter, as Nemanja Kaloper, Matthew Kleban, and Damien Martin showed, so the black hole "
        "horizon is the surface $R = 1.085\\,r_s$, $t = \\infty$.",
    ],
    ("frw", "comoving_spherical", "radial"): [
        "The plane of $t$ and the comoving $r$ ($\\theta = \\pi/2$, $\\phi = 0$), the same at "
        "every other angle by spherical symmetry. The line element leaves "
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
        "The same universe along a line through the observer at $r = 0$, in the plane "
        "$\\theta = \\pi/2$: $x = r$ on the right is $\\phi = 0$ and $x = -r$ on the left is "
        "$\\phi = \\pi$, and spherical symmetry makes the two halves mirror images. The past light "
        "cone of an event on the observer's world line $x = 0$ flares out as it runs back toward "
        "the bang, and the Hubble sphere lies at the same distance on either side.",
    ],
    ("frw", "conformal_spherical", "radial"): [
        "The plane of the conformal time $\\eta$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$). For $k "
        "= 0$ the metric on it is $a^2(-d\\eta^2 + dr^2)$; the scale factor "
        "multiplies both terms and drops out of the null condition, so the rays are straight 45° "
        "lines whatever $a(\\eta)$ is.",
        "The Hubble sphere and the Kretschmann scalar do depend on $a(\\eta)$, and for dust the "
        "Hubble sphere sits at $r = \\eta/2$, half the comoving radius of the particle horizon.",
    ],
    ("ellis_bronnikov", "spherical", "radial"): [
        "The plane of $t$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$), with $r$ running "
        "from one mouth of the wormhole, $r < 0$, through the throat at $r = 0$ to the other "
        "mouth, $r > 0$. The metric on this plane is $-dt^2 + dr^2$, exactly flat, so the rays are "
        "straight 45° lines that pass through the throat without bending.",
        "The throat is where the spheres are smallest. The angular part of the metric is "
        "$g_{\\theta\\theta} = r^2 + \\ell^2$, so the areal radius $R = \\sqrt{r^2 + \\ell^2}$ "
        "takes its least value, $\\ell$, at $r = 0$, where $\\partial_r R = 0$. The faint vertical "
        "lines are the spheres $R = 1.5\\ell$, $2\\ell$, and $3\\ell$, one of each size on either "
        "side of the throat. The Kretschmann scalar $12\\ell^4/(r^2 + \\ell^2)^4$ is finite "
        "everywhere.",
    ],
    ("thin_shell_wormhole", "throat", "radial"): [
        "The plane of $t$ and $\\ell$ ($\\theta = \\pi/2$, $\\phi = 0$), drawn for $a = 1.25\\,r_s$, with "
        "$\\ell < 0$ on one side of the throat and $\\ell > 0$ on the other. The edges of the cones are "
        "$d\\ell/d(ct) = \\pm\\left(1 - r_s/(a + |\\ell|)\\right)$, so $ct \\mp \\ell_*$ is constant along a ray, "
        "with $\\ell_* = \\ell + r_s\\ln\\left(1 + |\\ell|/(a - r_s)\\right)\\mathrm{sgn}(\\ell)$. The cones "
        "are narrowest at the throat, where $d\\ell/d(ct) = \\pm(1 - r_s/a) = \\pm 0.2$, and every ray crosses it "
        "in a finite time.",
        "The throat carries the shell. There the slope of a ray is continuous and its rate of change jumps "
        "sign, since $r = a + |\\ell|$ has a kink at $\\ell = 0$. The faint vertical lines are the spheres of areal "
        "radius $1.5\\,r_s$, $2\\,r_s$, and $3\\,r_s$, one of each on either side. Off the throat the "
        "Kretschmann scalar is $12r_s^2/(a + |\\ell|)^6$, at most $12r_s^2/a^6$.",
    ],
    ("thin_shell_wormhole", "spherical", "radial"): [
        "The plane of $t$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$) on one side of the throat, drawn for "
        "$a = 1.25\\,r_s$, the same on the other side and at every other angle. The edges of the cones are "
        "$dr/d(ct) = \\pm(1 - r_s/r)$, so $ct \\mp r_*$ is constant along a ray, with "
        "$r_* = r + r_s\\ln(r/r_s - 1)$, as outside a Schwarzschild black hole.",
        "The chart stops at the throat $r = a$, before the cones close at $r_s$, so there is no horizon. A ray "
        "moving in reaches the throat at a finite $t$ and goes on into the other side, where $r$ grows again. "
        "The Kretschmann scalar $12r_s^2/r^6$ is at most $12r_s^2/a^6$.",
    ],
    ("gravastar", "interior", "radial"): [
        "The plane of $t$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$) inside the shell, drawn for $R = 1.25\\,r_s$ and "
        "$L = 2\\,r_s$, where $C = 64/195$. The edges of the cones are $dr/d(ct) = \\pm\\sqrt{C}\\left(1 - r^2/L^2\\right)$, "
        "so $ct \\mp x$ is constant along a ray, with $x = (L/\\sqrt{C})\\,\\mathrm{artanh}(r/L)$ the tortoise "
        "coordinate. The cones are widest at the centre, where $dr/d(ct) = \\pm 0.57$, and narrowest on the shell, "
        "where it is $\\pm 0.35$.",
        "The chart stops at the shell $r = R$, before the cones close at the de Sitter horizon $r = L$, and $t$ "
        "is the time of the exterior chart. A ray from the centre reaches the shell after $ct = 2.56\\,r_s$ and "
        "goes on into Schwarzschild's vacuum. The Kretschmann scalar is $24/L^4$ at every point.",
    ],
    ("gravastar", "interior", "through"): [
        "The line through the centre of the gravastar in the plane $\\theta = \\pi/2$, from the shell to the "
        "shell ($R = 1.25\\,r_s$, $L = 2\\,r_s$): $x = r$ on the right is $\\phi = 0$ and $x = -r$ on the left is "
        "$\\phi = \\pi$. Rays cross the centre smoothly, and the cones are widest there.",
    ],
    ("gravastar", "interior_tortoise", "radial"): [
        "The plane of $t$ and $x$ ($\\theta = \\pi/2$, $\\phi = 0$) inside the shell, drawn for $R = 1.25\\,r_s$ and "
        "$L = 2\\,r_s$, where the shell is at $x = 2.56\\,r_s$. On this plane $ds^2 = C\\left(-c^2dt^2 + dx^2\\right)/"
        "\\cosh^2(\\sqrt{C}\\,x/L)$, so every ray runs at 45° and every cone is the same. The faint vertical lines are "
        "the spheres of areal radius $0.5\\,r_s$ and $r_s$, which crowd toward the shell since "
        "$r = L\\tanh(\\sqrt{C}\\,x/L)$.",
    ],
    ("gravastar", "exterior", "radial"): [
        "The plane of $t$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$) outside the shell, drawn for $R = 1.25\\,r_s$. "
        "The edges of the cones are $dr/d(ct) = \\pm(1 - r_s/r)$, so $ct \\mp r_*$ is constant along a ray, with "
        "$r_* = r + r_s\\ln(r/r_s - 1)$, as outside a Schwarzschild black hole.",
        "The chart stops at the shell $r = R$, before the cones close at $r_s$, so there is no horizon. A ray "
        "moving in reaches the shell at a finite $t$, where $dr/d(ct) = -0.2$, and crosses into the de Sitter "
        "interior. Light leaving the shell arrives far away with its frequency lowered by the factor "
        "$\\sqrt{1 - r_s/R} = 0.45$. The Kretschmann scalar $12r_s^2/r^6$ is at most $12r_s^2/R^6$.",
    ],
    ("teo_wormhole", "spherical", "axis"): [
        "The plane of $t$ and $r$ on the axis of rotation ($\\theta = 0$) on one side of the throat, drawn for "
        "$a = 1/4$, where $N = 1 + b_0/r$. The dragging of frames carries a factor $\\sin^2\\theta$ and "
        "vanishes on the axis, so these rays are null geodesics that stay on it, with edges "
        "$dr/d(ct) = \\pm N^{-1}\\sqrt{1 - b_0/r}$.",
        "$ct \\mp r_*$ is constant along a ray, with $r_* = \\sqrt{r(r - b_0)} - b_0\\ln\\left(\\sqrt{r/b_0} + "
        "\\sqrt{r/b_0 - 1}\\right) + \\sqrt{2}\\,b_0\\,\\mathrm{artanh}\\sqrt{(r - b_0)/2r}$. The cones close "
        "toward the throat $r = b_0$ because $g_{rr}$ diverges there, while $N = 2$ stays finite, so there is no "
        "horizon: a ray moving in reaches the throat at a finite $t$ and goes on into the other side, where $r$ "
        "grows again.",
    ],
    ("teo_wormhole", "spherical", "equator"): [
        "The plane of $t$ and $r$ on the equator ($\\theta = \\pi/2$) with $\\phi$ divided out, "
        "$-c^2dt^2 + dr^2/(1 - b_0/r)$, the metric orthogonal to the circles of $\\phi$, drawn for $a = 1$. Its "
        "null curves are the shadows on $t$ and $r$ of the null geodesics of zero angular momentum, each "
        "turning in $\\phi$ at $d\\phi/d(ct) = 2ab_0^2/r^3$, and each cone is the future cone of the "
        "directions of zero angular momentum.",
        "On the equator $N = 1$, so $ct \\mp l$ is constant along a ray, with $l = \\sqrt{r(r - b_0)} + "
        "b_0\\ln\\left(\\sqrt{r/b_0} + \\sqrt{r/b_0 - 1}\\right)$ the proper distance from the throat, and "
        "the rays are the same for every spin. The dotted line is the ergosurface, $g_{tt} = 0$ at "
        "$r = \\sqrt{2a}\\,b_0 = 1.41\\,b_0$. Between the throat and the ergosurface no observer keeps "
        "$\\phi$ fixed, and the region exists for $|a| > 1/2$.",
    ],
    ("teo_wormhole", "proper_radial", "axis"): [
        "The plane of $t$ and $l$ on the axis of rotation ($\\theta = 0$), drawn for $a = 1/4$, with $l < 0$ "
        "on one side of the throat and $l > 0$ on the other. The dragging of frames vanishes on the axis, so "
        "these rays are null geodesics that stay on it, with edges $dl/d(ct) = \\pm 1/N$ and $N = 1 + b_0/r$.",
        "The cones are narrowest at the throat, where $N = 2$ and $dl/d(ct) = \\pm 1/2$, and they open to "
        "45° far from it on both sides. Every ray crosses the throat in a finite time, and $ct \\mp l_*$ is "
        "constant along it, with $l_* = \\pm\\left(\\sqrt{r(r - b_0)} - b_0\\ln\\left(\\sqrt{r/b_0} + "
        "\\sqrt{r/b_0 - 1}\\right) + \\sqrt{2}\\,b_0\\,\\mathrm{artanh}\\sqrt{(r - b_0)/2r}\\right)$.",
    ],
    ("teo_wormhole", "proper_radial", "equator"): [
        "The plane of $t$ and $l$ on the equator ($\\theta = \\pi/2$) with $\\phi$ divided out, "
        "$-c^2dt^2 + dl^2$, the metric orthogonal to the circles of $\\phi$, drawn for $a = 1$. It is exactly "
        "flat, so the shadows on $t$ and $l$ of the null geodesics of zero angular momentum are straight 45° "
        "lines through the throat, each turning in $\\phi$ at $d\\phi/d(ct) = 2ab_0^2/r^3$, and each cone "
        "is the future cone of the directions of zero angular momentum.",
        "The dotted lines are the ergosurface, $g_{tt} = 0$ at $r = \\sqrt{2a}\\,b_0$, which is "
        "$l = \\pm 1.37\\,b_0$. Between them, through the throat, no observer keeps $\\phi$ fixed. The "
        "ergoregion is a tube round the equator of the throat and reaches neither pole.",
    ],
    ("wormhole_time_machine", "wormhole", "speeding"): [
        "The plane of $t$ and $l$ on the axis of the acceleration ($\\theta = 0$), where "
        "$g_{tt} = -(1 + glF/c^2)^2$, for three units of $r_0/c$ about the moment the right mouth speeds up "
        "hardest, $g = 0.2\\,c^2/r_0$. The edges of the cones are $dl/d(ct) = \\pm(1 + glF/c^2)$. On the left "
        "half, $l \\le 0$, $F = 0$ and every ray is at 45°.",
        "To the right of the throat the cones open as $F$ rises, and the rays there cross a given stretch of "
        "$l$ in less of the time $t$. Clocks to the right of the throat run fast against $t$ by the factor "
        "$1 + glF/c^2$, as clocks higher up in the accelerated frame of the right mouth.",
    ],
    ("wormhole_time_machine", "wormhole", "slowing"): [
        "The plane of $t$ and $l$ on the axis of the acceleration ($\\theta = 0$), for three units of $r_0/c$ "
        "about the moment the right mouth slows hardest on its way out, $g = -0.2\\,c^2/r_0$. The edges of the "
        "cones are $dl/d(ct) = \\pm(1 + glF/c^2)$, and to the right of the throat the cones close as $F$ rises.",
        "Along this axis the acceleration now points toward the throat, and clocks to the right of the throat "
        "run slow against $t$. On the left half, $l \\le 0$, $F = 0$ and every ray is at 45°, as in the static "
        "wormhole.",
    ],
    ("wormhole_time_machine", "short_throat", "radial"): [
        "The plane of $t$ and $l$ ($\\theta = \\pi/2$, $\\phi = 0$) in the rest frame of a mouth of radius $b$, "
        "with $l < 0$ outside one mouth and $l > 0$ outside the other. Here $g_{tt} = -1$ and $g_{ll} = 1$, so "
        "every ray is at 45° and crosses the throat as a straight line.",
        "The throat $l = 0$ carries the curvature, a delta function, since $r = b + |l|$ has a kink there. "
        "The faint vertical lines are the spheres of areal radius $1.5\\,b$, $2\\,b$, and $3\\,b$, one of each "
        "on either side.",
    ],
    ("simpson_visser", "spherical", "bounce"): [
        "The plane of $t$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$), drawn for $a = r_s/2$, the black bounce, "
        "the same at every other angle by spherical symmetry. The edges of the cones are "
        "$dr/d(ct) = \\pm(1 - r_s/\\rho)$ with $\\rho = \\sqrt{r^2 + a^2}$, and they close on the two horizons "
        "$r = \\pm\\sqrt{r_s^2 - a^2}$, which is $\\pm 0.87\\,r_s$ here.",
        "Between the horizons $r$ is the time and the cones lie on their side. We take the future from the "
        "ingoing Eddington-Finkelstein chart, which runs smoothly across both horizons, so every cone there "
        "points to smaller $r$. The areal radius shrinks from $r_s$ to $a$ at $r = 0$, where Schwarzschild's "
        "singularity would stand, and grows back to $r_s$ at the second horizon, beyond which $t$ is a time "
        "again in a second asymptotically flat region. The Kretschmann scalar at $r = 0$ is "
        "$(9r_s^2 - 16ar_s + 12a^2)/a^6$, which is $256/r_s^4$ here.",
    ],
    ("simpson_visser", "spherical", "null"): [
        "The plane of $t$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$), drawn for $a = r_s$, the one way wormhole. "
        "The cones close at $r = 0$ alone, where $1 - r_s/\\rho$ vanishes as $r^2/2r_s^2$, so the sphere of "
        "least area is an extremal horizon. Along a ray $ct \\mp r_*$ is constant, with "
        "$r_* = r + r_s\\,\\mathrm{arsinh}(r/r_s) - r_s(r_s + \\rho)/r$ and $\\rho = \\sqrt{r^2 + r_s^2}$, which "
        "diverges as $-2r_s^2/r$, so in this chart a ray reaches $r = 0$ only as $t \\to \\pm\\infty$.",
        "On both sides $t$ is a time, and the region $r < 0$ is the mirror image of the region $r > 0$. Each "
        "Eddington-Finkelstein chart carries one family of rays across $r = 0$. The Kretschmann scalar at "
        "$r = 0$ is $5/r_s^4$.",
    ],
    ("simpson_visser", "spherical", "wormhole"): [
        "The plane of $t$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$), drawn for $a = 2\\,r_s$, the traversable "
        "wormhole. No cone closes: $dr/d(ct) = \\pm(1 - r_s/\\rho)$ is least at the throat $r = 0$, where it is "
        "$\\pm(1 - r_s/a) = \\pm 1/2$, and every ray crosses the throat. Along a ray $ct \\mp r_*$ is constant, with "
        "$r_* = r + r_s\\,\\mathrm{arsinh}(r/a) + (r_s^2/k)\\left(\\arctan(r/k) + \\arctan(r_sr/k\\rho)\\right)$, "
        "$\\rho = \\sqrt{r^2 + a^2}$, and $k = \\sqrt{a^2 - r_s^2}$.",
        "A clock at rest in the throat ticks at $\\sqrt{1 - r_s/a} = 0.71$ of the rate of a distant one, where "
        "the Ellis-Bronnikov wormhole, $r_s = 0$, keeps one rate everywhere. The Kretschmann scalar at the "
        "throat is $(9r_s^2 - 16ar_s + 12a^2)/a^6$, which is $0.39/r_s^4$ here.",
    ],
    ("simpson_visser", "areal", "bounce"): [
        "The plane of $t$ and the areal radius $\\rho$ ($\\theta = \\pi/2$, $\\phi = 0$) on the side $r > 0$ of "
        "the sphere of least area, drawn for $a = r_s/2$. Here $g_{tt}$ and the spheres are Schwarzschild's, and "
        "the edges of the cones are $d\\rho/d(ct) = \\pm(1 - r_s/\\rho)\\sqrt{1 - a^2/\\rho^2}$, so they close on "
        "the horizon $\\rho = r_s$ and again on $\\rho = a$, the left edge of the chart.",
        "Inside the horizon $\\rho$ is the time, and we take the future toward smaller $\\rho$, as the ingoing "
        "Eddington-Finkelstein chart does. A ray reaches $\\rho = a$ at a finite $t$ with $d\\rho/d(ct) = 0$: "
        "the areal radius has stopped shrinking, and the ray goes on into the half $r < 0$, where $\\rho$ grows "
        "again. The Kretschmann scalar on $\\rho = a$ is $256/r_s^4$ here.",
    ],
    ("simpson_visser", "areal", "null"): [
        "The plane of $t$ and the areal radius $\\rho$ ($\\theta = \\pi/2$, $\\phi = 0$) on the side $r > 0$, "
        "drawn for $a = r_s$. The horizon and the sphere of least area are one sphere, $\\rho = r_s$, the left "
        "edge of the chart, where $d\\rho/d(ct) = \\pm(1 - r_s/\\rho)\\sqrt{1 - r_s^2/\\rho^2}$ vanishes as "
        "$(\\rho - r_s)^{3/2}$, and a ray moving in reaches it only as $t \\to \\infty$.",
        "Outside it $g_{tt}$ and the spheres are Schwarzschild's, so light circles at $\\rho = 3r_s/2$ and the "
        "innermost stable circular orbit is at $\\rho = 3r_s$, as for Schwarzschild's black hole. The "
        "Kretschmann scalar on $\\rho = r_s$ is $5/r_s^4$.",
    ],
    ("simpson_visser", "areal", "wormhole"): [
        "The plane of $t$ and the areal radius $\\rho$ ($\\theta = \\pi/2$, $\\phi = 0$) on one side of the "
        "throat, drawn for $a = 2\\,r_s$. The cones close toward the throat $\\rho = a$ because "
        "$g_{\\rho\\rho} = (1 - r_s/\\rho)^{-1}(1 - a^2/\\rho^2)^{-1}$ diverges there. On the throat "
        "$g_{tt} = -(1 - r_s/a) = -1/2$, so $\\partial_t$ is timelike right up to it, and $\\rho = a$ is the edge "
        "of this chart with no horizon on it.",
        "A ray moving in from $\\rho$ reaches the throat after the finite time $r_*/c$, which is $5.7\\,r_s/c$ "
        "from $\\rho = 4\\,r_s$ here, and goes on into the other side. The throat lies beyond Schwarzschild's "
        "circular orbit of light, $\\rho = 3r_s/2$, and the one circular orbit of light is on the throat itself. The "
        "Kretschmann scalar at the throat is $0.39/r_s^4$.",
    ],
    ("simpson_visser", "eddington_finkelstein_ingoing", "bounce"): [
        "The plane of $v$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$), drawn for $a = r_s/2$ with $v - r$ as the "
        "vertical axis, so that the rays moving left, $v = $ const, run at 45°. The rays moving right have "
        "$dv/dr = 2(1 - r_s/\\rho)^{-1}$ with $\\rho = \\sqrt{r^2 + a^2}$, and stand vertical on both horizons, "
        "$r = \\pm 0.87\\,r_s$.",
        "The chart runs from the region $r > 0.87\\,r_s$ through the black hole, across the sphere of least "
        "area $r = 0$, and through the second horizon into the region $r < -0.87\\,r_s$. Between the horizons "
        "both edges of every future cone point to smaller $r$, so whatever crosses the first horizon crosses "
        "$r = 0$ and comes out in the second region. Seen from there the horizon $r = -0.87\\,r_s$ is a white "
        "hole's: rays leave it and none go back in.",
    ],
    ("simpson_visser", "eddington_finkelstein_ingoing", "null"): [
        "The plane of $v$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$), drawn for $a = r_s$ with $v - r$ as the "
        "vertical axis, so that the rays moving left, $v = $ const, run at 45°. The rays moving right have "
        "$dv/dr = 2(1 - r_s/\\rho)^{-1}$ and stand vertical at $r = 0$ alone: the extremal horizon is itself a "
        "ray of that family, staying where it is.",
        "A ray moving left crosses $r = 0$ at a finite $v$ and goes on toward $r \\to -\\infty$, and no future "
        "directed curve returns across that horizon, which is why Simpson and Visser call this geometry a one "
        "way wormhole.",
    ],
    ("simpson_visser", "eddington_finkelstein_ingoing", "wormhole"): [
        "The plane of $v$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$), drawn for $a = 2\\,r_s$ with $v - r$ as the "
        "vertical axis, so that the rays moving left, $v = $ const, run at 45°. The rays moving right have "
        "$dv/dr = 2(1 - r_s/\\rho)^{-1}$, which is greatest at the throat $r = 0$, where it is $4$, so both "
        "families cross the throat.",
        "Here $v = ct + r_*$ with $r_*$ finite at every $r$, so this chart covers the same region as Simpson "
        "and Visser's own, the whole spacetime.",
    ],
    ("simpson_visser", "eddington_finkelstein_outgoing", "bounce"): [
        "The plane of $u$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$), drawn for $a = r_s/2$ with $u + r$ as the "
        "vertical axis, so that the rays moving right, $u = $ const, run at 45°. The rays moving left have "
        "$du/dr = -2(1 - r_s/\\rho)^{-1}$ with $\\rho = \\sqrt{r^2 + a^2}$, and stand vertical on both horizons, "
        "$r = \\pm 0.87\\,r_s$.",
        "The chart is the time reverse of the ingoing one. Between the horizons both edges of every future "
        "cone point to larger $r$: whatever enters from the region $r < -0.87\\,r_s$ crosses the sphere of "
        "least area $r = 0$ and comes out through $r = 0.87\\,r_s$, which is a white hole's horizon for the "
        "region beyond it.",
    ],
    ("simpson_visser", "eddington_finkelstein_outgoing", "null"): [
        "The plane of $u$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$), drawn for $a = r_s$ with $u + r$ as the "
        "vertical axis, so that the rays moving right, $u = $ const, run at 45°. The rays moving left have "
        "$du/dr = -2(1 - r_s/\\rho)^{-1}$ and stand vertical at $r = 0$ alone, the extremal horizon.",
        "A ray moving right crosses $r = 0$ at a finite $u$, coming from $r \\to -\\infty$. This is the "
        "horizon through which a region $r < 0$ empties into the region $r > 0$ to its future, the crossing "
        "that the ingoing chart leaves out.",
    ],
    ("simpson_visser", "eddington_finkelstein_outgoing", "wormhole"): [
        "The plane of $u$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$), drawn for $a = 2\\,r_s$ with $u + r$ as the "
        "vertical axis, so that the rays moving right, $u = $ const, run at 45°. The rays moving left have "
        "$du/dr = -2(1 - r_s/\\rho)^{-1}$, steepest at the throat $r = 0$, where it is $-4$, and both families "
        "cross the throat.",
        "Here $u = ct - r_*$ with $r_*$ finite at every $r$, so this chart, the ingoing one, and Simpson and "
        "Visser's own all cover the whole spacetime.",
    ],
    ("damour_solodukhin", "spherical", "radial"): [
        "The plane of $t$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$) on one side of the throat, drawn for "
        "$\\lambda = 0.2$, the same on the other side and at every other angle. The edges of the cones are "
        "$dr/d(ct) = \\pm\\sqrt{(1 - r_s/r)(1 - r_s/r + \\lambda^2)}$, so $ct \\mp r_*$ is constant along a ray, "
        "with $r_* = \\sqrt{(r - r_s)(ar - r_s)}/a + \\left((1 + a)r_s/2a^{3/2}\\right)"
        "\\ln\\left(\\left(2ar - (1 + a)r_s + 2\\sqrt{a(r - r_s)(ar - r_s)}\\right)/\\lambda^2r_s\\right)$ "
        "and $a = 1 + \\lambda^2$, which vanishes at the throat.",
        "The cones close toward the throat $r = r_s$ because $g_{rr} = (1 - r_s/r)^{-1}$ diverges there. "
        "There $g_{tt} = -\\lambda^2$, so $\\partial_t$ is timelike right up to the throat, and $r = r_s$ is "
        "the edge of this chart with no horizon on it. A ray moving in from $r$ reaches the throat after the "
        "finite time $r_*(r)/c$, which is $5.5\\,r_s/c$ from $r = 2\\,r_s$ here and grows as "
        "$(r_s/c)\\ln(1/\\lambda^2)$ as $\\lambda \\to 0$, and it goes on into the other side. The Kretschmann "
        "scalar at the throat is $(1 + 24\\lambda^4)/(4\\lambda^4r_s^4)$.",
    ],
    ("damour_solodukhin", "rescaled", "radial"): [
        "The plane of $t$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$) on one side of the throat, drawn for "
        "$\\lambda = 0.2$, with $t$ the proper time of a clock at rest far away. The edges of the cones are "
        "$dr/d(ct) = \\pm\\sqrt{(1 - r_s/r)(1 - r_s/ar)}$ with $a = 1 + \\lambda^2$, so they open to 45° far "
        "from the throat, and $ct \\mp \\sqrt{a}\\,r_*$ is constant along a ray, with $r_*$ the tortoise "
        "coordinate of Damour and Solodukhin's chart.",
        "The cones close toward the throat $r = r_s$ because $g_{rr} = (1 - r_s/r)^{-1}$ diverges there, and "
        "$g_{tt} = -\\lambda^2/a$ on it. A ray moving in from $r = 2\\,r_s$ reaches the throat after "
        "$5.6\\,r_s/c$ of this time and goes on into the other side. For small $\\lambda$ the time from $r$ "
        "is $\\left(r + r_s\\ln(r/r_s - 1) - r_s + r_s\\ln(4/\\lambda^2)\\right)/c$, Schwarzschild's tortoise "
        "coordinate and a constant that grows as $\\ln(1/\\lambda^2)$.",
    ],
    ("damour_solodukhin", "throat", "radial"): [
        "The plane of $t$ and $\\rho$ ($\\theta = \\pi/2$, $\\phi = 0$), drawn for $\\lambda = 0.2$, with "
        "$\\rho < 0$ on one side of the throat and $\\rho > 0$ on the other. The metric on it is "
        "$(1 - r_s/ar)\\left(-c^2dt^2 + r^2d\\rho^2\\right)$ with $a = 1 + \\lambda^2$, so the edges of the cones "
        "are $d\\rho/d(ct) = \\pm 1/r$ and $ct \\mp x$ is constant along a ray, with "
        "$x = r_s\\left((2 + \\lambda^2)\\rho + \\lambda^2\\sinh\\rho\\right)/2a$.",
        "The cones are widest at the throat $\\rho = 0$, where $d\\rho/d(ct) = \\pm 1/r_s$, and every ray "
        "crosses it. The coordinate $\\rho$ stretches the neighbourhood of the throat: the faint vertical "
        "lines are the spheres of areal radius $1.5\\,r_s$, $2\\,r_s$, and $3\\,r_s$, at $|\\rho| = 4.0$, "
        "$4.7$, and $5.3$, and between the two spheres of radius $2\\,r_s$ a ray spends $11.2\\,r_s/c$. "
        "The Kretschmann scalar is $(1 + 24\\lambda^4)/(4\\lambda^4r_s^4)$ at the throat.",
    ],
    ("damour_solodukhin", "isotropic", "radial"): [
        "The plane of $t$ and the isotropic radius $r$ ($\\theta = \\pi/2$, $\\phi = 0$), drawn for "
        "$\\lambda = 0.2$. The areal radius is $R = r(1 + r_s/4r)^2$, least at the throat $r = r_s/4$, and the "
        "whole of the other side lies between the throat and $r = 0$, which is its far end. The edges of the "
        "cones are $dr/d(ct) = \\pm\\sqrt{(4r - r_s)^2 + \\lambda^2(4r + r_s)^2}\\,16r^2/(4r + r_s)^3$, so "
        "$ct \\mp r_*\\,\\mathrm{sgn}(4r - r_s)$ is constant along a ray, with $r_*$ the tortoise coordinate of "
        "the sphere of areal radius $R$.",
        "Every component of the metric is finite at the throat, where $dr/d(ct) = \\pm\\lambda/2$, and every ray "
        "crosses it. A ray moving in slows as $r^2$ toward $r = 0$ and reaches it only as $t \\to \\infty$. The "
        "faint vertical lines are the spheres of areal radius $1.5\\,r_s$ and $2\\,r_s$, at $r = 0.93\\,r_s$ and "
        "$1.46\\,r_s$ on this side and at $r_s^2/16r = 0.067\\,r_s$ and $0.043\\,r_s$ on the other. The "
        "Kretschmann scalar is $(1 + 24\\lambda^4)/(4\\lambda^4r_s^4)$ at the throat.",
    ],
    ("damour_solodukhin", "einstein_rosen", "radial"): [
        "The plane of $t$ and $u$ ($\\theta = \\pi/2$, $\\phi = 0$), drawn for $\\lambda = 0.2$, with "
        "$r = r_s + u^2$, $u < 0$ on one side of the throat and $u > 0$ on the other. The edges of the cones are "
        "$du/d(ct) = \\pm\\sqrt{(1 + \\lambda^2)u^2 + \\lambda^2r_s}/2(u^2 + r_s)$, so $ct \\mp r_*\\,\\mathrm{sgn}(u)$ "
        "is constant along a ray, with $r_*$ the tortoise coordinate of the sphere $r = r_s + u^2$. Every "
        "component of the metric is finite at the throat $u = 0$, where the cones are narrowest, "
        "$du/d(ct) = \\pm\\lambda/2\\sqrt{r_s}$, and every ray crosses it.",
        "A ray takes the time $2r_*(r)/c$ to pass from the sphere of radius $r$ on one side to the sphere of "
        "the same radius on the other, $11\\,r_s/c$ for $r = 2\\,r_s$ here, of which the logarithm "
        "$(r_s/c)\\ln(1/\\lambda^2)$ on each side is spent beside the throat. The faint vertical lines are the "
        "spheres of areal radius $1.5\\,r_s$, $2\\,r_s$, and $3\\,r_s$, one of each on either side. The "
        "Kretschmann scalar is $(1 + 24\\lambda^4)/(4\\lambda^4r_s^4)$ at the throat.",
    ],
    ("fisher_jnw", "spherical", "radial"): [
        "The plane of $t$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$), drawn for $\\gamma = 1/2$, each point in the "
        "plane a 2-sphere of area $4\\pi r^2f^{1-\\gamma}$ with $f = 1 - b/r$. The edges of the cones are "
        "$dr/d(ct) = \\pm f^{\\gamma}$, so $ct \\mp r_*$ is constant along a ray, with "
        "$r_* = \\sqrt{r(r - b)} + b\\ln\\left((\\sqrt{r} + \\sqrt{r - b})/\\sqrt{b}\\right)$ at this $\\gamma$, "
        "which vanishes at $r = b$.",
        "The cones narrow toward $r = b$ and close only on it, where the spheres have zero area and the "
        "Kretschmann scalar diverges as $(r - b)^{2\\gamma - 4}$. A ray moving in from $r$ reaches that "
        "singularity after the finite time $r_*(r)/c$, which is $2.3\\,b/c$ from $r = 2b$ here, and a ray "
        "leaves it for infinity at every moment, so no horizon hides it. At $\\gamma = 1$ the time is "
        "infinite and $r = b$ is Schwarzschild's horizon.",
    ],
    ("fisher_jnw", "jnw", "radial"): [
        "The plane of $t$ and $R$ ($\\theta = \\pi/2$, $\\phi = 0$), drawn for $\\gamma = 1/2$, where "
        "$R = r - 3b/4$ and the singularity is $R = b/4$. The edges of the cones are "
        "$dR/d(ct) = \\pm f^{\\gamma}$ with $f = (4R - b)/(4R + 3b)$, the cones of the spherical chart moved "
        "over by $3b/4$, and $ct \\mp r_*$ is constant along a ray.",
        "In this radius $\\Gamma^R{}_{\\theta\\theta} = -R$ for every $\\gamma$, as in flat space. Each point "
        "in the plane is a 2-sphere of area $4\\pi\\left(R - b/4\\right)^{1/2}\\left(R + 3b/4\\right)^{3/2}$, which "
        "vanishes on the singularity, where the Kretschmann scalar diverges. A ray moving in from "
        "$R = 5b/4$ reaches it after $2.3\\,b/c$.",
    ],
    ("fisher_jnw", "isotropic", "radial"): [
        "The plane of $t$ and the isotropic radius $\\rho$ ($\\theta = \\pi/2$, $\\phi = 0$), drawn for "
        "$\\gamma = 1/2$, with $r = \\rho\\left(1 + b/4\\rho\\right)^2$ and the singularity at $\\rho = b/4$. The "
        "edges of the cones are $d\\rho/d(ct) = \\pm h^{2\\gamma - 1}\\left(1 + b/4\\rho\\right)^{-2}$ with "
        "$h = (4\\rho - b)/(4\\rho + b)$, and $ct \\mp r_*$ is constant along a ray.",
        "For $\\gamma > 1/2$ the cones close on the singularity and for $\\gamma < 1/2$ they open there without "
        "bound. At $\\gamma = 1/2$ the power of $h$ drops out and they keep the width "
        "$d\\rho/d(ct) = \\pm 1/4$ on it, so every ray drawn meets the singularity at that slope. The "
        "spheres have zero area there and the Kretschmann scalar diverges.",
    ],
    ("fisher_jnw", "harmonic", "radial"): [
        "The plane of $t$ and $u$ ($\\theta = \\pi/2$, $\\phi = 0$), drawn for $m = k/2$, which is "
        "$\\gamma = 1/2$. Spatial infinity is $u = 0$, the singularity is $u \\to \\infty$, and the scalar field "
        "grows in proportion to $u$. The edges of the cones are $du/d(ct) = \\pm e^{-2mu}\\sinh^2(ku)/k^2$, "
        "and $ct \\pm r_*$ is constant along a ray, with $r = 2k/(1 - e^{-2ku})$.",
        "The cones close as $u^2$ toward $u = 0$, which a ray reaches only as $t \\to \\pm\\infty$, and open as "
        "$e^{2(k - m)u}/4k^2$ at large $u$. A ray moving toward larger $u$ runs through all of it in a finite "
        "time and reaches the singularity $1.6\\,k/c$ after passing $ku = 1$. The faint vertical lines are "
        "the spheres of areal radius $4k$, $2k$, $k$, and $k/2$, at $ku = 0.28$, $0.64$, $1.5$, and $2.8$.",
    ],
    ("hartle_thorne", "hartle_thorne", "axis"): [
        "The plane of $t$ and $r$ on the axis of rotation ($\\theta = 0$) from the surface of the star out "
        "($R = 6\\,m$, $a = m/4$, $q = 4a^2$). The dragging of frames carries a factor $\\sin^2\\theta$ and vanishes "
        "on the axis, so these rays are null geodesics that stay on it, with edges "
        "$dr/d(ct) = \\pm F\\sqrt{(1 + 2h_2)/(1 - 2j_2)}$.",
        "The cones are narrowest at the surface, where $dr/d(ct) = \\pm 0.6691$, and open to 45° far from the "
        "star. Schwarzschild's edges at the same radius are $\\pm 0.6667$, so at these values the spin and the "
        "quadrupole moment move each edge by less than half a per cent, and the chart ends at $r = R$, where the "
        "star begins.",
    ],
    ("hartle_thorne", "hartle_thorne", "equator"): [
        "The plane of $t$ and $r$ on the equator ($\\theta = \\pi/2$) with $\\phi$ divided out, "
        "$-F(1 - h_2)c^2dt^2 + (1 + j_2)dr^2/F$, the metric orthogonal to the circles of $\\phi$, from the "
        "surface of the star out ($R = 6\\,m$, $a = m/4$, $q = 4a^2$). Its null curves are the shadows on $t$ and "
        "$r$ of the null geodesics of zero angular momentum, each turning in $\\phi$ at $d\\phi/d(ct) = 2am/r^3$, "
        "and each cone is the future cone of the directions of zero angular momentum.",
        "At the surface such a ray is carried around the axis at $2amc/R^2 = 0.014\\,c$, the dragging of inertial "
        "frames, which falls off as $1/r^2$. The edges of the cones there are $dr/d(ct) = \\pm 0.6656$, "
        "a little inside Schwarzschild's $\\pm 0.6667$, where on the axis they lie a little outside.",
    ],
    ("hartle_thorne", "painleve_gullstrand", "axis"): [
        "The plane of $t$ and $r$ on the axis of rotation ($\\theta = 0$) from the surface of the star out "
        "($R = 6\\,m$, $a = m/4$), with $t$ the proper time of an observer falling from rest far away. The edges of "
        "the cones are $dr/d(ct) = \\pm 1 - \\sqrt{2m/r}$: each is a 45° cone tipped inward by the speed of that "
        "observer, $0.58\\,c$ at the surface.",
        "$ct + r - 2\\sqrt{2mr} + 4m\\ln\\left(\\sqrt{r/2m} + 1\\right)$ is constant along an ingoing ray and "
        "$ct - r - 2\\sqrt{2mr} - 4m\\ln\\left(\\sqrt{r/2m} - 1\\right)$ along an outgoing one. The observer's "
        "world line, $dr/d(ct) = -\\sqrt{2m/r}$, runs midway between the two edges of every cone, and the "
        "chart ends at $r = R$, where the star begins.",
    ],
    ("hartle_thorne", "painleve_gullstrand", "equator"): [
        "The plane of $t$ and $r$ on the equator ($\\theta = \\pi/2$) with $\\phi$ divided out, "
        "$-c^2dt^2 + \\left(dr + \\sqrt{2m/r}\\,c\\,dt\\right)^2$, the metric orthogonal to the circles of $\\phi$, "
        "from the surface of the star out ($R = 6\\,m$, $a = m/4$). It is the metric of the axis for every spin, so "
        "the shadows on $t$ and $r$ of the null geodesics of zero angular momentum are the rays of the axis, each "
        "turning in $\\phi$ at $d\\phi/d(ct) = 2am/r^3$, and each cone is the future cone of the directions of "
        "zero angular momentum.",
        "An observer falling from rest far away with no angular momentum moves inward at $\\sqrt{2m/r}\\,c$, "
        "midway between the two edges of every cone, and turns about the axis at the same $2amc/r^3$, "
        "which at the surface is a speed of $0.014\\,c$ around it.",
    ],
    ("witten_black_hole", "witten", "radial"): [
        "The plane of $t$ and $r$, each point in the diagram a single event. The edges of the cones are "
        "$dr/d(ct) = \\pm\\tanh(\\lambda r)$, so $ct \\mp \\sigma$ is constant along a ray, with "
        "$\\sigma = \\ln\\sinh(\\lambda r)/\\lambda$.",
        "Far from the horizon the cones open to 45°, as in flat space, and they close on $r = 0$, the "
        "horizon, which a ray reaches only as $t \\to \\pm\\infty$. Since $r$ is the proper distance, the "
        "horizon lies a finite distance from every event, and the Ricci scalar there is $4\\lambda^2$.",
    ],
    ("witten_black_hole", "schwarzschild_gauge", "radial"): [
        "The plane of $t$ and $x$, drawn for $m = 1$, each point in the diagram a single event. The edges "
        "of the cones are $dx/d(ct) = \\pm\\left(1 - m\\,e^{-2\\lambda x}\\right)$, and the horizon is "
        "$x = \\ln(m)/2\\lambda$, here $x = 0$. Outside it the cones narrow toward the vertical as the "
        "horizon is approached, and far away they open to 45°.",
        "Inside the horizon $x$ is the time. We take the future from the ingoing Eddington-Finkelstein "
        "chart, which makes the region the black hole, every cone pointing to smaller $x$. The Ricci scalar "
        "$4\\lambda^2m\\,e^{-2\\lambda x}$ grows without bound as $x \\to -\\infty$, the singularity, and no "
        "observer spends a proper time longer than $\\pi/2\\lambda c$ between the horizon and it.",
    ],
    ("witten_black_hole", "dilaton", "radial"): [
        "The plane of $t$ and $w$, drawn for $m = 1$, each point in the diagram a single event where the "
        "dilaton is $\\Phi = -\\tfrac{1}{2}\\ln w$. The edges of the cones are "
        "$dw/d(ct) = \\pm 2\\lambda\\left(w - m\\right)$, so along a ray $w - m$ grows or decays as "
        "$e^{\\pm 2\\lambda ct}$.",
        "The cones close on the horizon $w = m$ and open again inside it, where $w$ is the time and we take "
        "the future toward smaller $w$, the black hole. They meet the singularity $w = 0$ at the finite "
        "slope $2\\lambda m$. There $e^{-2\\Phi}$ vanishes and the Ricci scalar $4\\lambda^2m/w$ diverges.",
    ],
    ("witten_black_hole", "conformal", "radial"): [
        "The plane of $t$ and $\\sigma$, drawn for $m = 1$, each point in the diagram a single event. The "
        "metric is $\\left(1 + m\\,e^{-2\\lambda\\sigma}\\right)^{-1}$ times Minkowski's, so the light rays are "
        "the lines of constant $ct \\pm \\sigma$, at 45°.",
        "The chart covers the outside of the horizon. The horizon is $\\sigma \\to -\\infty$, where the "
        "conformal factor vanishes as $e^{2\\lambda\\sigma}/m$, and toward $\\sigma \\to \\infty$ the factor "
        "tends to $1$, the plane is flat, and the dilaton is linear in $\\sigma$.",
    ],
    ("witten_black_hole", "kruskal", "plane"): [
        "The plane of $U$ and $V$, drawn with $(U + V)/2$ up and $(V - U)/2$ across, each point in the "
        "diagram a single event. Only $g_{UV}$ is nonzero, so the light rays are the lines $U = $ const and "
        "$V = $ const at 45°.",
        "The horizons are the lines $UV = 0$, and the singularity is the hyperbola $UV = 1$, with one branch "
        "in the future and one in the past, where the Ricci scalar $4\\lambda^2/(1 - UV)$ diverges. The "
        "quadrant $U < 0 < V$ is the outside of the black hole, with $V = \\sinh(\\lambda r)\\,e^{\\lambda ct}$ "
        "and $U = -\\sinh(\\lambda r)\\,e^{-\\lambda ct}$ in Witten's $t$ and $r$, and the opposite quadrant is "
        "a second copy of it. The four regions and their causal relations are those of Kruskal's plane for "
        "the Schwarzschild black hole.",
    ],
    ("witten_black_hole", "eddington_finkelstein_ingoing", "finkelstein"): [
        "The plane of $v$ and $x$, drawn for $m = 1$ with $v - x$ as the vertical axis so that the ingoing "
        "rays, $v = $ const, run at 45°. The outgoing family has "
        "$dv/dx = 2\\left(1 - m\\,e^{-2\\lambda x}\\right)^{-1}$, so it stands vertical at the horizon "
        "$x = \\ln(m)/2\\lambda$: the horizon is itself an outgoing ray that stays where it is.",
        "The cones cross the horizon smoothly and keep tipping. Inside, both edges of every future cone "
        "point to smaller $x$, toward the singularity $x \\to -\\infty$.",
    ],
    ("witten_black_hole", "eddington_finkelstein_outgoing", "finkelstein"): [
        "The plane of $u$ and $x$, drawn for $m = 1$ with $u + x$ as the vertical axis so that the outgoing "
        "rays, $u = $ const, run at 45°. The retarded chart crosses the other horizon. Inside "
        "$x = \\ln(m)/2\\lambda$ both edges of every future cone point to larger $x$: this is the white "
        "hole, which nothing from outside can enter.",
        "The chart is the time reverse of the ingoing one: its $g_{ux}$ is $-1$ where the ingoing chart's "
        "$g_{vx}$ is $+1$.",
    ],
    **ROBERTS_CAPTIONS,
    ("einstein_rosen_bridge", "bridge", "radial"): [
        "The plane of $t$ and $u$ ($\\theta = \\pi/2$, $\\phi = 0$), with $r = r_s + u^2$, $u < 0$ on one sheet and "
        "$u > 0$ on the other. The edges of the cones are $du/d(ct) = \\pm|u|/2(u^2 + r_s)^{3/2}$, so "
        "$ct \\mp r_*\\,\\mathrm{sgn}(u)$ is constant along a ray, with $r_* = r + r_s\\ln(r/r_s - 1)$ the tortoise "
        "coordinate of the sphere of radius $r$.",
        "Every component of the metric is finite on the bridge $u = 0$, and $g_{tt}$ and the determinant vanish "
        "there, so the cones close: $r_* \\to -\\infty$ on the bridge, and a ray sent toward it from either sheet "
        "arrives only as $t \\to \\infty$. The bridge is the horizon of the Schwarzschild spacetime. The faint "
        "vertical lines are the spheres of areal radius $1.5\\,r_s$, $2\\,r_s$, and $3\\,r_s$, one of each on "
        "either sheet, and the Kretschmann scalar $12r_s^2/(u^2 + r_s)^6$ is $12/r_s^4$ on the bridge.",
    ],
    ("einstein_rosen_bridge", "spherical", "radial"): [
        "The plane of $t$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$) on one sheet, the same on the other and at "
        "every other angle. The edges of the cones are $dr/d(ct) = \\pm(1 - r_s/r)$, so $ct \\mp r_*$ is constant "
        "along a ray, with $r_* = r + r_s\\ln(r/r_s - 1)$.",
        "The cones close toward the bridge $r = r_s$, where $g^{rr} = 0$, and a ray moving in reaches it only as "
        "$t \\to \\infty$. Einstein and Rosen's space has no sphere smaller than the bridge: the second sheet "
        "begins there, with $r$ growing again.",
    ],
    ("einstein_rosen_bridge", "isotropic", "radial"): [
        "The plane of $t$ and the isotropic radius $r$ ($\\theta = \\pi/2$, $\\phi = 0$). The areal radius is "
        "$R = r(1 + r_s/4r)^2$, least on the bridge $r = r_s/4$, and the whole of the other sheet lies between the "
        "bridge and $r = 0$, which is its far end. The edges of the cones are "
        "$dr/d(ct) = \\pm 16r^2|4r - r_s|/(4r + r_s)^3$, so $ct \\mp r_*\\,\\mathrm{sgn}(4r - r_s)$ is constant along "
        "a ray, with $r_* = R + r_s\\ln(R/r_s - 1)$.",
        "The cones close on the bridge, where $g_{tt} = 0$, and again toward $r = 0$, where a ray slows as $r^2$ "
        "and arrives only as $t \\to \\infty$. The faint vertical lines are the spheres of areal radius "
        "$1.5\\,r_s$ and $2\\,r_s$, at $r = 0.93\\,r_s$ and $1.46\\,r_s$ on one sheet and at "
        "$r_s^2/16r = 0.067\\,r_s$ and $0.043\\,r_s$ on the other.",
    ],
    ("einstein_rosen_bridge", "charged_spherical", "radial"): [
        "The plane of $t$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$) on one sheet of the charged bridge with a mass, "
        "drawn for $r_q = \\sqrt{3}\\,r_s/2$. The edges of the cones are $dr/d(ct) = \\pm(1 - r_s/r - r_q^2/r^2)$, "
        "so $ct \\mp r_*$ is constant along a ray, with "
        "$r_* = r + \\tfrac{9}{8}r_s\\ln(r/r_+ - 1) - \\tfrac{1}{8}r_s\\ln(2r/r_s + 1)$.",
        "The cones close toward the bridge $r_+ = 3r_s/2$, where $g^{rr} = 0$, and a ray moving in reaches it "
        "only as $t \\to \\infty$. With the charge entering as Einstein and Rosen took it, $g^{rr}$ has this one "
        "positive zero whatever the mass, and the Kretschmann scalar "
        "$(12r_s^2r^2 + 48r_sr_q^2r + 56r_q^4)/r^8$ is finite on it.",
    ],
    ("einstein_rosen_bridge", "charged_bridge", "radial"): [
        "The plane of $t$ and $u$ ($\\theta = \\pi/2$, $\\phi = 0$) of the charged bridge with no mass, with "
        "$r = \\sqrt{u^2 + r_q^2}$, $u < 0$ on one sheet and $u > 0$ on the other. The edges of the cones are "
        "$du/d(ct) = \\pm|u|/\\sqrt{u^2 + r_q^2}$, so $ct \\mp r_*\\,\\mathrm{sgn}(u)$ is constant along a ray, with "
        "$r_* = r + \\tfrac{1}{2}r_q\\ln\\left((r - r_q)/(r + r_q)\\right)$.",
        "The distance from the bridge along a radius is $|u|$, and far from the bridge the cones open to 45°. "
        "They close on the bridge $u = 0$, where $g_{tt} = 0$, and a ray sent toward it arrives only as "
        "$t \\to \\infty$: the bridge is a horizon, as the neutral one is. The faint vertical lines are the spheres "
        "of areal radius $1.5\\,r_q$, $2\\,r_q$, and $3\\,r_q$, one of each on either sheet, and the Kretschmann "
        "scalar $56r_q^4/(u^2 + r_q^2)^4$ is $56/r_q^4$ on the bridge.",
    ],
    ("morris_thorne", "spherical", "radial"): [
        "The plane of $t$ and the areal radius $r$ ($\\theta = \\pi/2$, $\\phi = 0$). The metric leaves $\\Phi(r)$ "
        "and $b(r)$ free. With $\\Phi = 0$ and $b = b_0^2/r$ it is the Ellis-Bronnikov wormhole of throat radius "
        "$\\ell = b_0$, with $r^2 = r_{\\rm EB}^2 + \\ell^2$, and these rays conserve $t \\mp \\sqrt{r^2 - b_0^2} "
        "= t \\mp r_{\\rm EB}$: they are the same rays as in the Ellis-Bronnikov chart.",
        "In this areal chart the cones close toward the throat at $r = b_0$, as they would at a "
        "horizon, because $g_{rr} = (1 - b_0^2/r^2)^{-1}$ diverges there. But $g_{tt} = -1$ stays "
        "finite, so $\\partial_t$ is timelike right up to the throat, and $r = b_0$ is only the "
        "edge of this chart. The rays reach the throat in finite $t$ and pass into the other "
        "mouth, which the Ellis-Bronnikov chart, running through the throat, covers in full. Below "
        "$b_0$ the formula gives a metric on this plane with no null directions.",
    ],
    ("misner", "misner", "plane"): [
        "The plane of $T$ and $\\psi$ ($y = z = 0$), a cylinder drawn unrolled, its edges $\\psi = 0$ and "
        "$\\psi = \\psi_0$ one line. One family of light rays runs straight up at constant $\\psi$ and "
        "crosses $T = 0$, and the other follows $dT/d\\psi = -T/2$, so the cones tip over as $T$ climbs.",
        "Below $T = 0$ the circles of constant $T$ are spacelike, and the rays of the second family wind round "
        "the cylinder toward $T = 0$ without reaching it. The circle $T = 0$ is itself a light ray, the "
        "chronology horizon, and above it every circle of constant $T$ is a closed timelike curve.",
    ],
    ("misner", "milne", "plane"): [
        "The plane of $t$ and $\\chi$ ($y = z = 0$) in the region $T < 0$, a cylinder drawn unrolled, its edges "
        "$\\chi = 0$ and $\\chi = \\psi_0/2$ one line, with $g_{\\chi\\chi} = c^2t^2$. The circles of "
        "constant $t$ shrink as $t$ climbs toward $0$, and the cones close up against them.",
        "Both families of light rays wind round the cylinder, $c\\,dt = \\pm ct\\,d\\chi$, and reach "
        "$t = 0$ only as $\\chi \\to \\pm\\infty$. One family is the rays of constant $\\psi$ in "
        "Misner's coordinates, which cross the chronology horizon.",
    ],
    ("misner", "rindler", "plane"): [
        "The plane of $\\eta$ and $\\xi$ ($y = z = 0$) in the region $T > 0$, its edges $\\eta = 0$ and "
        "$\\eta = \\psi_0/2$ one line, with $g_{\\eta\\eta} = -\\xi^2$. Every vertical line is a closed "
        "timelike curve, of proper length $\\xi\\psi_0/2$.",
        "The light rays run as $\\xi = \\xi_0e^{\\pm\\eta}$ for any number $\\xi_0$, and the cones close toward $\\xi = 0$, the "
        "chronology horizon, where the closed curves turn into closed null geodesics.",
    ],
    **{("string_wave", "null_conical", view): [
        f"The plane of $u$ and $v$ ($r = \\ell$, $\\phi = {phi}$), on the side of the string {side}, drawn with "
        "$z = (v - u)/2$ and $ct = (u + v)/2$. The rays moving right keep their $u$ and travel with the wave. A "
        f"curve moving left crosses the wave with $dv/du = F$, and here $F = {sign}\\tfrac{{1}}{{2}}\\ell A''$: it "
        f"keeps $v {other} \\tfrac{{1}}{{2}}\\ell A'$, so it is moved along $v$ while the string accelerates and "
        "comes out of the pulse where it would have been without it.",
        "The rays moving right are null geodesics. $\\Gamma^r{}_{uu} = -\\partial_rF/2$ turns light crossing the "
        "wave out of this plane, so a curve moving left is a null curve, and a null geodesic only outside the pulse.",
    ] for view, phi, side, sign, other in (
        ("toward", "0", "its displacement $A$ points to", "-", "+"),
        ("away", "\\pi", "its displacement $A$ points away from", "", "-"))},
    ("string_wave", "isotropic", "beside"): [
        "The plane of $u$ and $v$ ($x = \\ell/4$, $y = 0$), at a fixed isotropic distance from the string, drawn "
        "with $z = (v - u)/2$ and $ct = (u + v)/2$. The chart moves with the string, so the wave shows on this "
        "plane as $g_{uu} = -2xA''$ alone. The rays moving right keep their $u$ and travel with the wave, and a "
        "curve moving left keeps $v + 2xA'$.",
        "The rays moving right are null geodesics. $\\Gamma^x{}_{uu} = (\\rho/\\ell)^{2 - 2b}A''$ turns light "
        "crossing the wave out of this plane, so a curve moving left is a null curve, and a null geodesic only "
        "outside the pulse.",
    ],
    **{("string_wave", "moving_string", view): [
        f"The plane of $u$ and $V$ ($X = {X}$, $Y = 0$), a line parallel to the string's resting place, drawn with "
        f"$z = (V - u)/2$ and $ct = (u + V)/2$. The pulse carries the string from $X = 0$ out to $\\ell/2$ and "
        f"back, {passing}. On this plane $g_{{uu}} = ((\\rho/\\ell)^{{2b - 2}} - 1)A'^2$, positive inside "
        "$\\rho = \\ell$, so the cones lean toward $+z$ on either side of the crest, where the string is moving, "
        "and stand at 45° on the crest and away from the pulse.",
        "The rays moving right keep their $u$, travel with the wave, and are null geodesics. A curve moving left "
        f"is moved along $V$ by $\\int g_{{uu}}\\,du = {shift}\\,\\ell$ in crossing the pulse. Light crossing the "
        "wave is turned out of this plane, so that curve is a null curve, and a null geodesic only outside the pulse.",
    ] for view, X, passing, shift in (
        ("behind", "-\\ell/4", "between $\\ell/4$ and $3\\ell/4$ from this line", "0.65"),
        ("ahead", "3\\ell/4", "coming within $\\ell/4$ of this line", "0.78"))},
    ("spinning_string", "proper_radius", "inside"): [
        "The cylinder of $t$ and $\\phi$ ($r = r_c/2$, $z = 0$), opened along the line $\\phi = \\pm\\pi$ and drawn with $br\\phi$ across, so that its left and right edges are that one line. The metric on it is $-(c\\,dt + a\\,d\\phi)^2 + b^2r^2d\\phi^2$, the same at every point, so its null curves are straight: $c\\,dt = (br - a)\\,d\\phi$ and $c\\,dt = -(br + a)\\,d\\phi$. Inside $r_c = a/b$ both go down in $t$ toward $+\\phi$, so every horizontal line, run toward $+\\phi$, points into the future cones, and the circle of constant $t$, $r$, and $z$ is a closed timelike curve.",
        "Neither curve is a null geodesic: the spacetime is flat, and light launched along either one leaves the cylinder for larger $r$, turned by $\\Gamma^r{}_{\\phi\\phi} = -b^2r$.",
    ],
    ("spinning_string", "proper_radius", "outside"): [
        "The cylinder of $t$ and $\\phi$ ($r = 3r_c/2$, $z = 0$), opened along the line $\\phi = \\pm\\pi$ and drawn with $br\\phi$ across, so that its left and right edges are that one line. The metric on it is $-(c\\,dt + a\\,d\\phi)^2 + b^2r^2d\\phi^2$, the same at every point, so its null curves are straight: $c\\,dt = (br - a)\\,d\\phi$ and $c\\,dt = -(br + a)\\,d\\phi$. Outside $r_c = a/b$ the first climbs in $t$ toward $+\\phi$ and the second toward $-\\phi$, five times as steeply, so the cones lean toward $+\\phi$, and the horizontal lines, circles of constant $t$, lie outside every cone and are spacelike.",
        "Neither curve is a null geodesic: the spacetime is flat, and light launched along either one leaves the cylinder for larger $r$, turned by $\\Gamma^r{}_{\\phi\\phi} = -b^2r$.",
    ],
    ("spinning_string", "rescaled_radius", "inside"): [
        "The cylinder of $t$ and $\\phi$ ($\\rho = a/2$, $z = 0$), opened along the line $\\phi = \\pm\\pi$ and drawn with $\\rho\\phi$ across, so that its left and right edges are that one line. The metric on it is $-(c\\,dt + a\\,d\\phi)^2 + \\rho^2d\\phi^2$, the same at every point, so its null curves are straight: $c\\,dt = (\\rho - a)\\,d\\phi$ and $c\\,dt = -(\\rho + a)\\,d\\phi$. Inside $\\rho = a$ both go down in $t$ toward $+\\phi$, so every horizontal line, run toward $+\\phi$, points into the future cones, and the circle of constant $t$, $\\rho$, and $z$ is a closed timelike curve.",
        "Neither curve is a null geodesic: the spacetime is flat, and light launched along either one leaves the cylinder for larger $\\rho$, turned by $\\Gamma^\\rho{}_{\\phi\\phi} = -b^2\\rho$.",
    ],
    ("spinning_string", "rescaled_radius", "outside"): [
        "The cylinder of $t$ and $\\phi$ ($\\rho = 3a/2$, $z = 0$), opened along the line $\\phi = \\pm\\pi$ and drawn with $\\rho\\phi$ across, so that its left and right edges are that one line. The metric on it is $-(c\\,dt + a\\,d\\phi)^2 + \\rho^2d\\phi^2$, the same at every point, so its null curves are straight: $c\\,dt = (\\rho - a)\\,d\\phi$ and $c\\,dt = -(\\rho + a)\\,d\\phi$. Outside $\\rho = a$ the first climbs in $t$ toward $+\\phi$ and the second toward $-\\phi$, five times as steeply, so the cones lean toward $+\\phi$, and the horizontal lines, circles of constant $t$, lie outside every cone and are spacelike.",
        "Neither curve is a null geodesic: the spacetime is flat, and light launched along either one leaves the cylinder for larger $\\rho$, turned by $\\Gamma^\\rho{}_{\\phi\\phi} = -b^2\\rho$.",
    ],
    ("spinning_string", "circumference_radius", "outside"): [
        "The cylinder of $t$ and $\\phi$ ($R = a$, $z = 0$), opened along the line $\\phi = \\pm\\pi$ and drawn "
        "with $R\\phi$ across, so that its left and right edges are that one line. The metric on it is "
        "$-c^2dt^2 - 2ac\\,dt\\,d\\phi + R^2d\\phi^2$, the same at every point, so its null curves are straight: "
        "$c\\,dt = (\\sqrt{R^2 + a^2} - a)\\,d\\phi$ and $c\\,dt = -(\\sqrt{R^2 + a^2} + a)\\,d\\phi$. "
        "At every $R > 0$ the first climbs in $t$ toward $+\\phi$ and the second toward $-\\phi$, so the "
        "horizontal lines, circles of circumference $2\\pi R$, are spacelike. The cones lean farther toward "
        "$+\\phi$ as $R$ shrinks, and at $R = 0$, where the chart ends, the circle is null.",
        "Neither curve is a null geodesic: the spacetime is flat, and light launched along either one leaves the "
        "cylinder for larger $R$.",
    ],
    ("spinning_string", "helical", "inside"): [
        "The cylinder of $\\tau$ and $\\tilde\\phi$ ($r = r_c/2$, $z = 0$), opened along $\\tilde\\phi = 0$ and "
        "drawn with $r\\tilde\\phi$ across. The line element is Minkowski's, so the null curves run at 45°. "
        "The right edge, $\\tilde\\phi = 2\\pi b$, is the left edge moved up by $2\\pi a/c$ in $\\tau$, so a "
        "circuit of the string at constant $t$ is the straight line $c\\tau = a\\tilde\\phi/b$ from an event on "
        "the left edge to the same event on the right.",
        "Inside $r_c = a/b$ that line climbs by $2\\pi a$ while it crosses $2\\pi br < 2\\pi a$, steeper than "
        "the null curves, so it lies inside the cones: a closed timelike curve through each of its events.",
    ],
    ("spinning_string", "helical", "outside"): [
        "The cylinder of $\\tau$ and $\\tilde\\phi$ ($r = 3r_c/2$, $z = 0$), opened along $\\tilde\\phi = 0$ "
        "and drawn with $r\\tilde\\phi$ across. The line element is Minkowski's, so the null curves run at 45°. "
        "The right edge, $\\tilde\\phi = 2\\pi b$, is the left edge moved up by $2\\pi a/c$ in $\\tau$.",
        "Outside $r_c = a/b$ the circuit of the string at constant $t$, the straight line "
        "$c\\tau = a\\tilde\\phi/b$, climbs by $2\\pi a$ while it crosses $2\\pi br > 2\\pi a$, so it lies "
        "outside the cones and is spacelike. It ends on the right edge at the event it left on the left edge, "
        "one turn of a helix in $\\tau$.",
    ],
    ("gott_time_machine", "grant_rindler", "plane"): [
        "The plane of $\\eta$ and $\\xi$ ($Y = z = 0$) in the region of closed timelike curves, with "
        "$g_{\\eta\\eta} = -\\xi^2$. The edge $\\eta = a$ is the edge $\\eta = 0$ moved by $b$ along $Y$. "
        "The light rays run as $\\xi = \\xi_0e^{\\pm\\eta}$ for any number $\\xi_0$, and the cones close toward $\\xi = 0$, the "
        "chronology horizon.",
        "An event at $\\xi$ and its $n$th image lie $n^2b^2 - 4\\xi^2\\sinh^2(na/2)$ apart in squared "
        "interval, so a timelike line joins them beyond the $n$th polarised hypersurface, "
        "$\\xi = nb/(2\\sinh(na/2))$. Those hypersurfaces crowd toward the horizon as $n$ grows, and a "
        "closed timelike curve passes through every event with $\\xi > 0$.",
    ],
    ("gott_time_machine", "grant_milne", "plane"): [
        "The plane of $\\tau$ and $\\chi$ ($Y = z = 0$) to the past of the chronology horizon, with "
        "$g_{\\chi\\chi} = c^2\\tau^2$. The edge $\\chi = a$ is the edge $\\chi = 0$ moved by $b$ along "
        "$Y$. Both families of light rays wind toward $\\tau = 0$, $c\\,d\\tau = \\pm c\\tau\\,d\\chi$, "
        "and reach it only as $\\chi \\to \\pm\\infty$.",
        "An event and its $n$th image lie $n^2b^2 + 4c^2\\tau^2\\sinh^2(na/2)$ apart in squared interval, "
        "which is positive, so no closed timelike curve passes through this region. The closed curve of "
        "constant $\\tau$ through an event has length $\\sqrt{a^2c^2\\tau^2 + b^2}$, which shrinks to $b$ "
        "at the horizon.",
    ],
    ("ori_time_machine", "vacuum_core", "centre"): [
        "The plane of $T$ and $z$ ($x = y = 0$), a cylinder drawn unrolled, its edges $z = 0$ and $z = L$ one "
        "line. On it $f = 0$ and the metric is Misner's, $-2\\,dz\\,dT - T\\,dz^2$. One family of light rays "
        "runs straight up at constant $z$, and the other follows $dT/dz = -T/2$, so the cones tip over as $T$ "
        "climbs.",
        "The symmetry of $f$ under $x \\to -x$ and under $y \\to -y$ keeps every ray of this plane a null "
        "geodesic. The circle $T = 0$ is $N$, the one closed null geodesic of the core, and above it every "
        "circle of constant $T$ is a closed timelike curve.",
    ],
    ("ori_time_machine", "vacuum_core", "off_centre"): [
        "The cylinder of $T$ and $z$ ($x = 4\\,\\ell$, $y = 0$), drawn unrolled, its edges $z = 0$ and $z = L$ "
        "one line. On it $f = ax^2/2 = \\ell^2/2$ and the metric is $-2\\,dz\\,dT + (f - T)\\,dz^2$, so the "
        "circle of constant $T$ is spacelike below $T = f$, null there, and a closed timelike curve above.",
        "The vertical lines are light rays. The curves that tip over, $dT/dz = (f - T)/2$, are null, and light "
        "launched along one leaves the cylinder for larger $x$, turned by $\\Gamma^x{}_{zz} = -ax/2$. Each "
        "circle of the core turns null at its own $T = f(x, y)$, sooner along $y$, where $f < 0$, and later "
        "along $x$.",
    ],
    ("ori_time_machine", "foliation", "centre"): [
        "The plane of $t$ and $z$ ($x = y = 0$), a cylinder drawn unrolled, its edges $z = 0$ and $z = L$ one "
        "line. On it $t = T$ and the metric is $-2\\,dz\\,dt - t\\,dz^2$. One family of light rays runs "
        "straight up at constant $z$, and the other follows $dt/dz = -t/2$, winding round the cylinder toward "
        "$t = 0$ from below.",
        "Every surface of constant $t < 0$ is spacelike, and the surface $t = 0$ is spacelike everywhere but "
        "on this central circle, the closed null geodesic $N$. Above it the circles of constant $t$ are closed "
        "timelike curves.",
    ],
    ("ori_time_machine", "foliation", "off_centre"): [
        "The cylinder of $t$ and $z$ ($x = 4\\,\\ell$, $y = 0$), drawn unrolled, its edges $z = 0$ and $z = L$ "
        "one line. On it the metric is $-2\\,dz\\,dt + (ex^2 - t)\\,dz^2$ with $ex^2 = 2\\,\\ell^2$, so the "
        "circle of constant $t$ is spacelike below $t = 2\\,\\ell^2$, null there, and a closed timelike curve "
        "above.",
        "At $t = 0$, when the central circle is already null, this one still has circumference "
        "$L\\sqrt{2}\\,\\ell$. The vertical lines are light rays, and the curves that tip over, "
        "$dt/dz = (ex^2 - t)/2$, are null: light launched along one leaves the cylinder for larger $x$, turned "
        "by $\\Gamma^x{}_{zz} = -ax/2$.",
    ],
    ("ori_time_machine", "brinkmann", "plane"): [
        "The plane of $u$ and $v$ ($x = y = 0$), drawn against $(v - u)/2$ and $(u + v)/2$, where the metric "
        "is $-2\\,du\\,dv$ and every light ray runs at 45°. The coordinates cover $u < 0$, and a circuit of $z$ "
        "carries the line $u = -2$ onto $u = -2e^{-L/2}$ with $v$ stretched by $e^{L/2}$, so the strip between "
        "them is one copy of the core's central plane.",
        "Each hyperbola $uv = -2T$ is a circle of constant $T$, closed by that boost: spacelike where $v < 0$ "
        "and timelike where $v > 0$. The ray $v = 0$ is the closed null geodesic $N$. Each circuit of it is "
        "shorter in $u$, an affine parameter, by the factor $e^{-L/2}$, so $N$ runs round without end in a "
        "finite affine length and is incomplete to the future.",
    ],
    ("minkowski", "spherical", "radial"): [
        "The plane of $t$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$) in flat spacetime. "
        "The metric on it is $-c^2dt^2 + dr^2$, and every ray is at 45°. Each ingoing ray meets "
        "an outgoing one on the axis $r = 0$.",
    ],
    ("minkowski", "spherical_null", "radial"): [
        "The same plane of flat spacetime in the double null chart, $u$ and $v$ ($\\theta = "
        "\\pi/2$, $\\phi = 0$), drawn against $(v - u)/2$ and $(u + v)/2$. Only "
        "$g_{uv}$ is nonzero on it, so a null direction has $du\\,dv = 0$, and the light rays are "
        "the coordinate lines $u = $ const and $v = $ const themselves.",
    ],
    ("minkowski", "cartesian", "tx"): [
        "The plane of $t$ and $x$ ($y = z = 0$) in flat spacetime. The metric on it is "
        "$-c^2dt^2 + dx^2$, and every ray is at 45°.",
    ],
    ("minkowski", "rindler", "tx"): [
        "The plane of $T$ and $X$ ($Y = Z = 0$) in the Rindler chart, which covers the "
        "wedge $X > 0$ seen by observers of constant proper acceleration $a$, with $g_{TT} = "
        "-a^2X^2/c^4$. The cones close toward $X = 0$, where $g_{TT}$ vanishes, and a ray takes "
        "infinite $T$ to get there, $cT = \\pm(c^2/a)\\ln X + $ const.",
        "The edge of the wedge is the Rindler horizon, which only the accelerated observers have. "
        "The spacetime is flat, its Kretschmann scalar is zero, and the rays carry on across $X = "
        "0$ into the rest of Minkowski spacetime, which the Cartesian chart covers whole.",
    ],
    ("de_sitter", "static_spherical", "radial"): [
        "The plane of $t$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$) in the static chart, the same "
        "at every other angle by spherical symmetry. The cones close at "
        "the cosmological horizon $r = \\sqrt{3/\\Lambda}$, where $g^{rr} = 1 - \\Lambda r^2/3$ "
        "vanishes, on the far side from the observer at $r = 0$.",
        "The static chart covers only the observer's side of the horizon. Beyond it $t$ is "
        "spacelike, and the cones point along the outgoing rays to larger $r$, into the region "
        "the observer's own light goes on to reach.",
    ],
    ("de_sitter", "static_spherical", "through"): [
        "The static chart along a line through the observer in the plane $\\theta = "
        "\\pi/2$: $x = r$ on the right is $\\phi = 0$ and $x = -r$ on the left is $\\phi = \\pi$, "
        "and spherical symmetry makes the two halves mirror images. The horizon crosses the line "
        "on both sides, at $x = \\pm\\sqrt{3/\\Lambda}$, and beyond it the cones point away from "
        "the observer.",
    ],
    ("de_sitter", "flat_slicing", "tx"): [
        "The plane of $t$ and $x$ ($y = z = 0$) in the flat slicing ($ds^2 = -c^2dt^2 + "
        "e^{2Ht}dx^2$), the cones narrowing as $e^{-Ht}$ toward the future and opening out "
        "toward the past.",
        "A ray covers only a finite comoving distance however long it runs, $x = \\pm(c/H)e^{-Ht} "
        "+ $ const, so an observer at $x = 0$ has an event horizon.",
    ],
    ("einstein_static", "hyperspherical", "radial"): [
        "The plane of $t$ and $\\chi$ ($\\theta = \\pi/2$, $\\phi = 0$), the same at every other angle by "
        "spherical symmetry about the pole. The metric on it is $-c^2dt^2 + R^2d\\chi^2$, exactly flat, so "
        "every ray is a straight 45° line, $ct = \\pm R\\chi + $ const, and every cone is the same.",
        "The line $\\chi = 0$ is the pole and $\\chi = \\pi$ its antipode, where the spheres of constant "
        "$\\chi$, of area $4\\pi R^2\\sin^2\\chi$, shrink to points; a ray that reaches either passes "
        "through it and comes back as a ray of the other family. The Kretschmann scalar $12/R^4$ is the "
        "same everywhere.",
    ],
    ("einstein_static", "hyperspherical", "through"): [
        "The Einstein static universe along a great circle through the pole, in the plane $\\theta = "
        "\\pi/2$: $x = \\chi$ on the right is $\\phi = 0$ and $x = -\\chi$ on the left is $\\phi = \\pi$, "
        "and the two edges $x = \\pm\\pi$ are one line, the world line of the antipode.",
        "Light sent out from the pole at $t = 0$ reaches the antipode after $\\pi R/c$, where the rays "
        "sent out in every direction meet again, and is back at the pole after $2\\pi R/c$.",
    ],
    ("einstein_static", "static_areal", "radial"): [
        "The plane of $t$ and the areal radius $r = R\\sin\\chi$ ($\\theta = \\pi/2$, $\\phi = 0$). The "
        "cones narrow toward $r = R$ because $g_{rr} = R^2/(R^2 - r^2)$ diverges there, and a ray "
        "reaches it in the finite time $\\pi R/2c$, as $ct = \\pm R\\arcsin(r/R) + $ const.",
        "The line $r = R$ is the equator of the three sphere, its largest sphere, where the chart ends. "
        "The rays run on through it into the far hemisphere, which the hyperspherical chart covers, and "
        "the Kretschmann scalar $12/R^4$ is the same on both sides.",
    ],
    ("milne", "comoving_hyperbolic", "through"): [
        "The Milne universe along a line through the comoving particle at $\\chi = 0$, in the plane $\\theta = "
        "\\pi/2$: $x = \\chi$ on the right is $\\phi = 0$ and $x = -\\chi$ on the left is $\\phi = \\pi$. The "
        "edges of the cones are $d\\chi/d(ct) = \\pm 1/ct$, so a ray that crosses $\\chi = 0$ at the time $t_1$ "
        "runs along $\\chi = \\pm\\ln(t/t_1)$, and the cones open out toward $t = 0$.",
        "The whole line $t = 0$ is one event, $T = R = 0$ in the inertial chart, where every comoving "
        "particle starts, since $g_{\\chi\\chi} = c^2t^2$ vanishes there. The Kretschmann scalar is zero "
        "everywhere, and the past light cone of every event reaches every $\\chi$, so the Milne universe "
        "has no particle horizon.",
    ],
    ("milne", "comoving_spherical", "radial"): [
        "The plane of $t$ and the comoving radius $r = \\sinh\\chi$ ($\\theta = \\pi/2$, $\\phi = 0$), "
        "the same at every other angle by spherical symmetry. The edges of the cones are $dr/d(ct) = "
        "\\pm\\sqrt{1 + r^2}/ct$, so $\\ln t \\pm \\chi$ is constant along a ray, as in the hyperbolic "
        "chart, and the cones widen toward $t = 0$ and toward large $r$.",
        "The areal radius is $ctr$, and the square of its gradient is $1$ everywhere, so no sphere is "
        "trapped. The Kretschmann scalar is zero everywhere.",
    ],
    ("milne", "logarithmic_time", "radial"): [
        "The plane of the logarithmic time $\\tau = t_0\\ln(t/t_0)$ and $\\chi$ ($\\theta = \\pi/2$, "
        "$\\phi = 0$). The metric on it is $e^{2\\tau/t_0}\\left(-c^2d\\tau^2 + c^2t_0^2d\\chi^2\\right)$, flat "
        "up to the factor $e^{2\\tau/t_0}$, so every ray is a straight 45° line, $\\tau = \\pm t_0\\chi + $ "
        "const, and every cone is the same.",
        "The moment $t = t_0$ is $\\tau = 0$, and the event $t = 0$ lies at $\\tau \\to -\\infty$. The "
        "Kretschmann scalar is zero everywhere.",
    ],
    ("milne", "inertial", "through"): [
        "The Milne universe in the inertial chart of the comoving particle at $\\chi = 0$, along a line "
        "through it in the plane $\\theta = \\pi/2$: $x = R$ on the right is $\\phi = 0$ and $x = -R$ on "
        "the left is $\\phi = \\pi$. The metric is Minkowski's, so every ray is a straight 45° line, and "
        "the chart covers the inside $R < cT$ of the future light cone of the event $T = R = 0$.",
        "Each comoving particle moves along the straight line $R = cT\\tanh\\chi$ from that event, and "
        "each moment of constant $t$ is the hyperbola $c^2T^2 - R^2 = c^2t^2$.",
    ],
    ("coleman_de_luccia", "wall", "into_flat"): [
        "The plane of $\\psi$ and $\\xi$ ($\\theta = \\pi/2$, $\\phi = 0$) outside the light cone of the bubble's "
        "centre, for a decay into flat space, in units of the de Sitter radius $\\ell$. The metric on it is "
        "$d\\xi^2 - \\rho^2d\\psi^2$, so the edges of the cones are $d\\xi/d\\psi = \\pm\\rho$: they close toward "
        "the light cone $\\xi = 0$, which a ray reaches only as $\\psi \\to \\pm\\infty$, and they are widest "
        "where $\\rho$ is greatest, on de Sitter's equator.",
        "The wall is the line $\\xi = 4\\ell/5$, at rest in these coordinates, with flat space to its left, "
        "$\\rho = \\xi$, and de Sitter space to its right, $\\rho = \\ell\\sin((\\xi - \\xi_0)/\\ell)$ with "
        "$\\xi_0 = 4\\ell/5 - \\ell\\arcsin(4/5)$. A ray crosses the wall with a kink, since $\\rho$ is continuous "
        "and its slope is not. On the far side $\\rho$ falls to zero again at $\\xi = 3.01\\,\\ell$, the light "
        "cone of the event opposite the bubble's centre.",
    ],
    ("coleman_de_luccia", "wall", "out_of_flat"): [
        "The plane of $\\psi$ and $\\xi$ ($\\theta = \\pi/2$, $\\phi = 0$) outside the light cone of the bubble's "
        "centre, for a decay of flat space, in units of the anti-de Sitter radius $\\ell$. The metric on it is "
        "$d\\xi^2 - \\rho^2d\\psi^2$, so the edges of the cones are $d\\xi/d\\psi = \\pm\\rho$, and they close "
        "toward the light cone $\\xi = 0$, which a ray reaches only as $\\psi \\to \\pm\\infty$.",
        "The wall is the line $\\xi = \\ell\\ln 3$, at rest in these coordinates, with anti-de Sitter space to its "
        "left, $\\rho = \\ell\\sinh(\\xi/\\ell)$, and flat space to its right, where $\\rho$ grows by one for each "
        "unit of $\\xi$. To the right the plane is Rindler's: an observer at fixed $\\xi$ accelerates away from "
        "the bubble at $c^2/\\rho$ and the wall, which accelerates at $c^2/\\bar\\rho$ on the flat side, never "
        "reaches that observer.",
    ],
    ("coleman_de_luccia", "open", "zero"): [
        "The plane of $\\tau$ and $\\chi$ ($\\theta = \\pi/2$, $\\phi = 0$) inside a bubble whose new vacuum has "
        "zero energy, with $\\ell$ any length. The edges of the cones are $d\\chi/d(c\\tau) = \\pm 1/a$ with "
        "$a = c\\tau$, so $\\ln\\tau \\pm \\chi$ is constant along a ray and the cones open out toward $\\tau = 0$.",
        "The whole line $\\tau = 0$ is the light cone of the bubble's centre, where $a$ vanishes and the curvature "
        "stays finite. This universe is Milne's, flat space in hyperbolic slices, and the wall lies beyond the "
        "light cone, outside the plane.",
    ],
    ("coleman_de_luccia", "open", "positive"): [
        "The plane of $\\tau$ and $\\chi$ ($\\theta = \\pi/2$, $\\phi = 0$) inside a bubble whose new vacuum has "
        "positive energy, in units of its de Sitter radius $\\ell$. The edges of the cones are "
        "$d\\chi/d(c\\tau) = \\pm 1/a$ with $a = \\ell\\sinh(c\\tau/\\ell)$, so $\\ln\\tanh(c\\tau/2\\ell) \\pm \\chi$ "
        "is constant along a ray.",
        "The cones open out toward the light cone $\\tau = 0$ and close as the universe inflates, and a ray that "
        "leaves $\\chi = 0$ at the time $\\tau_1$ never passes $\\chi = -\\ln\\tanh(c\\tau_1/2\\ell)$, its event "
        "horizon. The Kretschmann scalar is $24/\\ell^4$ everywhere.",
    ],
    ("coleman_de_luccia", "open", "negative"): [
        "The plane of $\\tau$ and $\\chi$ ($\\theta = \\pi/2$, $\\phi = 0$) inside a bubble whose new vacuum has "
        "negative energy, in units of its anti-de Sitter radius $\\ell$. The edges of the cones are "
        "$d\\chi/d(c\\tau) = \\pm 1/a$ with $a = \\ell\\sin(c\\tau/\\ell)$, so $\\ln\\tan(c\\tau/2\\ell) \\pm \\chi$ "
        "is constant along a ray, and the cones are narrowest at $c\\tau = \\pi\\ell/2$, where the universe is "
        "largest.",
        "The scale factor returns to zero at $c\\tau = \\pi\\ell$, the top edge. In this exact solution the "
        "Kretschmann scalar is $24/\\ell^4$ there as everywhere, a coordinate singularity like the one at "
        "$\\tau = 0$, and Coleman and De Luccia argued that any small disturbance of the field makes it a true "
        "singularity, a crunch.",
    ],
    ("coleman_de_luccia", "open_conformal", "negative"): [
        "The plane of the conformal time $\\eta$ and $\\chi$ ($\\theta = \\pi/2$, $\\phi = 0$) inside the bubble. "
        "The metric on it is $a^2(-d\\eta^2 + d\\chi^2)$; the scale factor multiplies both terms and drops out of "
        "the null condition, so the rays are straight 45° lines whatever $a(\\eta)$ is.",
        "For a new vacuum of negative energy $\\eta = \\ln\\tan(c\\tau/2\\ell)$ and $a = \\ell/\\cosh\\eta$: the "
        "light cone of the bubble's centre is $\\eta \\to -\\infty$, the universe is largest at $\\eta = 0$, and "
        "the return of $a$ to zero is $\\eta \\to \\infty$.",
    ],
    ("coleman_de_luccia", "static_inside", "into_flat"): [
        "The plane of $t$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$) inside the wall for a decay into flat space, "
        "from the moment $t = 0$ the bubble appears, in units of the de Sitter radius $\\ell$ of the vacuum "
        "outside. The metric is Minkowski's, so every ray is a straight 45° line, and the chart ends at the wall, "
        "the hyperbola $r_w^2 = \\bar\\rho^2 + c^2t^2$ with $\\bar\\rho = 4\\ell/5$.",
        "The wall starts at rest and approaches the light cone $r = ct$ of the bubble's centre, with the "
        "acceleration $c^2/\\bar\\rho$ seen from this side. An outgoing ray that leaves the centre after $t = 0$ "
        "never catches it.",
    ],
    ("coleman_de_luccia", "static_inside", "out_of_flat"): [
        "The plane of $t$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$) inside the wall for a decay of flat space, "
        "from the moment $t = 0$ the bubble appears, in the static coordinates of the anti-de Sitter space "
        "inside, in units of its radius $\\ell$. The edges of the cones are $dr/d(ct) = \\pm(1 + r^2/\\ell^2)$, so "
        "$ct \\pm \\ell\\arctan(r/\\ell)$ is constant along a ray.",
        "The chart ends at the wall, $r_w^2 = \\bar\\rho^2 + (\\ell^2 + \\bar\\rho^2)\\tan^2(ct/\\ell)$ with "
        "$\\bar\\rho = 4\\ell/3$, which runs out to infinite $r$ at $ct = \\pi\\ell/2$. The open universe inside "
        "the bubble is the region above the ray $ct = \\ell\\arctan(r/\\ell)$ from the centre, and its moment "
        "$\\tau$ is the curve $\\sqrt{\\ell^2 + r^2}\\cos(ct/\\ell) = \\ell\\cos(c\\tau/\\ell)$.",
    ],
    ("coleman_de_luccia", "static_outside", "into_flat"): [
        "The plane of $t$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$) outside the wall for a decay into flat space, "
        "from the moment $t = 0$ the bubble appears, in the static coordinates of the de Sitter space outside, "
        "in units of its radius $\\ell$. The cones close at the cosmological horizon $r = \\ell$ of an observer "
        "at the bubble's centre, where $g^{rr}$ vanishes.",
        "The chart begins at the wall, $r_w^2 = \\bar\\rho^2 + (\\ell^2 - \\bar\\rho^2)\\tanh^2(ct/\\ell)$ with "
        "$\\bar\\rho = 4\\ell/5$, which starts at rest and reaches the horizon only as $t \\to \\infty$. Beyond "
        "the horizon $t$ is spacelike, and the cones point along the outgoing rays to larger $r$.",
    ],
    ("coleman_de_luccia", "static_outside", "out_of_flat"): [
        "The plane of $t$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$) outside the wall for a decay of flat space, "
        "from the moment $t = 0$ the bubble appears, in units of the anti-de Sitter radius $\\ell$ of the vacuum "
        "inside. The metric is Minkowski's, so every ray is a straight 45° line, and the chart begins at the "
        "wall, the hyperbola $r_w^2 = \\bar\\rho^2 + c^2t^2$ with $\\bar\\rho = 4\\ell/3$.",
        "An observer at rest at the radius $R$ meets the wall at $ct = \\sqrt{R^2 - \\bar\\rho^2}$, and a light ray sent "
        "inward from the wall at the moment it appears reaches that observer at $ct = R - \\bar\\rho$, so for a "
        "distant observer the warning comes a time $\\bar\\rho^2/2Rc$ ahead of the wall.",
    ],
    ("domain_wall", "planar", "tz"): [
        "The plane of $t$ and $z$ ($x = y = 0$), the same at every $x$ and $y$, in units of $1/k$. The "
        "metric on it is $-(1 - k|z|)^2c^2dt^2 + dz^2$, Rindler's on either side of the wall at $z = 0$, so "
        "the edges of the cones are $dz/d(ct) = \\pm(1 - k|z|)$, and a ray crosses the wall with the "
        "slope it had.",
        "The cones close toward the horizons $z = \\pm 1/k$, where $g_{tt}$ vanishes, and a ray takes an "
        "infinite time $t$ to reach either, running along $kct \\mp \\mathrm{sgn}(z)\\ln(1 - k|z|) = $ "
        "const. No Christoffel symbol turns a ray out of the plane, so every curve drawn is a null "
        "geodesic, and the only curvature is on the wall, $R_{tztz} = -2k\\,\\delta(z)$.",
    ],
    ("domain_wall", "inertial", "through"): [
        "One side of the wall in its inertial chart, along a line through the centre in the plane "
        "$\\theta = \\pi/2$: $x = R$ on the right is $\\phi = 0$ and $x = -R$ on the left is $\\phi = \\pi$, "
        "in units of $1/k$. The metric is Minkowski's, so every ray is a straight 45° line, and the chart "
        "covers the inside of the wall, the hyperbola $R^2 - c^2T^2 = 1/k^2$, which falls in, stops at "
        "$R = 1/k$ at $T = 0$, and recedes.",
        "Across the wall lies the same region of a second copy of Minkowski space, and a ray that reaches "
        "the wall runs on into it. The rays $R = c|T|$ through the centre at $T = 0$ bound the region "
        "the planar and global charts cover, between them and the wall, and are the horizons $|z| = 1/k$ "
        "of those charts.",
    ],
    ("randall_sundrum", "proper_distance", "ty"): [
        "The plane of $t$ and $y$ ($x_1 = x_2 = x_3 = 0$), the same at every place along the wall, in units of "
        "$1/k$. The metric on it is $-e^{-2k|y|}c^2dt^2 + dy^2$, so the edges of the cones are "
        "$dy/d(ct) = \\pm e^{-k|y|}$, widest on the wall at $y = 0$, and a ray crosses the wall with the slope it "
        "had.",
        "The cones close toward the horizons $y \\to \\pm\\infty$, where $g_{tt}$ vanishes, and a ray takes an "
        "infinite time $t$ and a finite affine parameter to reach either, running along "
        "$kct \\mp \\mathrm{sgn}(y)\\,(e^{k|y|} - 1) = $ const. No Christoffel symbol turns a ray out of the plane, "
        "so every curve drawn is a null geodesic.",
    ],
    ("randall_sundrum", "conformal", "tw"): [
        "The plane of $t$ and $w$ ($x_1 = x_2 = x_3 = 0$), in units of $1/k$. The metric on it is "
        "$(-c^2dt^2 + dw^2)/(1 + k|w|)^2$, conformal to Minkowski's plane, so every ray is a straight 45° line, "
        "$ct \\pm w = $ const, through the wall at $w = 0$ and on.",
        "The horizons lie at $w \\to \\pm\\infty$, off the drawing in every direction. No Christoffel symbol turns "
        "a ray out of the plane, so every curve drawn is a null geodesic.",
    ],
    ("randall_sundrum", "poincare", "tz"): [
        "The plane of $t$ and $z$ ($x_1 = x_2 = x_3 = 0$) on one side of the wall, in units of $1/k$. The metric "
        "on it is $(-c^2dt^2 + dz^2)/k^2z^2$, anti-de Sitter space in the Poincaré chart, so every ray is a "
        "straight 45° line, $ct \\pm z = $ const.",
        "The left edge is the wall, $z = 1/k$, where a ray passes into the mirror image of this side, and the "
        "conformal boundary of anti-de Sitter space, $z = 0$, lies in the part cut away. No Christoffel symbol "
        "turns a ray out of the plane, so every curve drawn is a null geodesic.",
    ],
    ("randall_sundrum", "two_walls", "tphi"): [
        "The plane of $t$ and $\\phi$ ($x_1 = x_2 = x_3 = 0$) between two walls with $kr_c = 1/2$, in units of "
        "$1/k$. The metric on it is $-e^{-2kr_c|\\phi|}c^2dt^2 + r_c^2d\\phi^2$, so the edges of the cones are "
        "$r_c\\,d\\phi/d(ct) = \\pm e^{-kr_c|\\phi|}$, and the cones at the edges are $e^{-\\pi/2} = 0.21$ as wide "
        "as on the wall of positive tension at $\\phi = 0$.",
        "The edges $\\phi = \\pm\\pi$ are one line, the wall of negative tension, and $\\phi$ and $-\\phi$ are one "
        "point, so the two halves of the drawing are two copies of the space between the walls and a ray that "
        "reaches either wall turns back from it. No Christoffel symbol turns a ray out of the plane, so every "
        "curve drawn is a null geodesic.",
    ],
    ("anti_de_sitter", "static_global", "radial"): [
        "The plane of $t$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$) in the global "
        "chart. The cones stay open everywhere, since $g^{rr} = 1 + r^2/L^2$ never vanishes. But "
        "$dt/dr = \\pm(1 + r^2/L^2)^{-1}$ falls off fast enough that a ray reaches $r \\to "
        "\\infty$, the conformal boundary, in the finite time $\\pi L/2$, and the rays flatten as "
        "they near it.",
    ],
    ("anti_de_sitter", "static_global", "through"): [
        "The global chart along a line through the centre in the plane $\\theta = \\pi/2$: "
        "$x = r$ on the right is $\\phi = 0$ and $x = -r$ on the left is $\\phi = \\pi$. A ray "
        "from the centre reaches the boundary on either side in the finite time $\\pi L/2$, and "
        "the rays flatten as they near it.",
    ],
    ("anti_de_sitter", "poincare", "tx"): [
        "The plane of $t$ and $x$ ($y = 0$, $z = L$) in the Poincaré chart. The metric "
        "on it is $(L^2/z^2)(-c^2dt^2 + dx^2)$, conformal to flat, so the rays are exact 45° "
        "lines, as they are at every $z$. The conformal boundary, which a ray reaches in finite "
        "time, lies off this plane, at $z \\to 0$.",
    ],
    ("dilaton_black_hole", "static", "radial"): [
        "The plane of $t$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$) of the Einstein metric, drawn for "
        "$r_d = r_s/2$, each point in the plane a 2-sphere of area $4\\pi r(r - r_d)$. On the plane the metric is "
        "Schwarzschild's for every charge, so the cones narrow toward the vertical as $r \\to r_s$, where "
        "$c\\,dt/dr = \\pm(1 - r_s/r)^{-1}$ diverges.",
        "Inside $r_s$, $t$ is a spacelike coordinate, and we take the future from the ingoing Eddington-Finkelstein "
        "chart, which makes that region the black hole, where every cone points to the singularity $r = r_d$. "
        "The spheres have zero area there and the Kretschmann scalar diverges, and at $r_s$ it is finite.",
    ],
    ("dilaton_black_hole", "eddington_finkelstein_ingoing", "finkelstein"): [
        "The plane of $v$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$) of the Einstein metric, drawn for "
        "$r_d = r_s/2$ with $v - r$ as the vertical axis, so that the ingoing rays, $v = $ const, run at 45°. "
        "The outgoing family has $dv/dr = 2(1 - r_s/r)^{-1}$, so it stands vertical at the horizon $r_s$, an "
        "outgoing ray that stays where it is.",
        "The cones cross $r_s$ smoothly and keep tipping. Inside, both edges of every future cone point to "
        "smaller $r$, so every future directed ray ends at the singularity $r = r_d$, where the area of the "
        "spheres vanishes.",
    ],
    ("dilaton_black_hole", "eddington_finkelstein_ingoing", "chart"): [
        "The same plane of $v$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$), drawn against the chart's own "
        "coordinates. The ingoing family is $v = $ const and runs horizontally here, since $v$ is itself a null "
        "coordinate. The outgoing family turns vertical at $r_s$ and leans back toward smaller $r$ inside it.",
    ],
    ("dilaton_black_hole", "eddington_finkelstein_outgoing", "finkelstein"): [
        "The plane of $u$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$) of the Einstein metric, drawn for "
        "$r_d = r_s/2$ with $u + r$ as the vertical axis, so that the outgoing rays, $u = $ const, run at 45°. "
        "The retarded chart crosses the other horizon. Inside $r_s$ both edges of every future cone point to "
        "larger $r$: this is the white hole, which nothing from outside can enter.",
        "The chart is the time reverse of the ingoing one: its $g_{ur}$ is $-1$ where the ingoing chart's "
        "$g_{vr}$ is $+1$.",
    ],
    ("dilaton_black_hole", "eddington_finkelstein_outgoing", "chart"): [
        "The same plane of $u$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$), drawn against the chart's own "
        "coordinates. The outgoing family is $u = $ const and runs horizontally here, since $u$ is itself a null "
        "coordinate. The ingoing family turns vertical at $r_s$ and leans toward larger $r$ inside it.",
    ],
    ("dilaton_black_hole", "string_magnetic", "radial"): [
        "The plane of $t$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$) of the string metric of the magnetically "
        "charged hole, drawn for $r_d = r_s/2$, each point in the plane a 2-sphere of area $4\\pi r^2$. The string "
        "metric is $e^{2\\varphi}$ times the Einstein metric, so its null rays on this plane are the Einstein "
        "metric's, $c\\,dt/dr = \\pm(1 - r_s/r)^{-1}$, and the cones close at $r_s$.",
        "The spheres keep the area $4\\pi r_d^2$ at the singularity $r = r_d$, where the Kretschmann scalar "
        "diverges. At the extremal charge $r_d = r_s$ the metric on the plane is "
        "$-c^2dt^2 + dr^2/(1 - r_s/r)^2$, and a ray takes an infinite time $t$ and an infinite affine parameter "
        "to reach $r_s$.",
    ],
    ("dilaton_black_hole", "string_electric", "radial"): [
        "The plane of $t$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$) of the string metric of the electrically "
        "charged hole, drawn for $r_d = r_s/2$, each point in the plane a 2-sphere of area $4\\pi(r - r_d)^2$. The "
        "string metric is $e^{2\\varphi}$ times the Einstein metric, so its null rays on this plane are the "
        "Einstein metric's and the cones close at $r_s$.",
        "With the areal radius $\\rho = r - r_d$ a moment of constant $t$ has the metric "
        "$d\\rho^2/(1 - (r_s - r_d)/\\rho) + \\rho^2d\\Omega^2$, Schwarzschild's with $r_s - r_d$ for $r_s$. The "
        "Kretschmann scalar diverges at $r = r_d$, where $g_{tt}$ vanishes with the area of the spheres.",
    ],
    ("mass_inflation", "ingoing", "tail"): [
        "The plane of $v$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$), drawn with $v - r$ as the vertical axis so "
        "that the ingoing rays, $v$ constant, run at 45°, for $r_q = 0.96\\,m_0$. Until $v_0$ the hole is "
        "Reissner-Nordström's with the mass $0.98\\,m_0$, and from $v_0$ on the tail falls in and the mass climbs to "
        "$m_0$. The curves $g^{rr} = 0$ are the apparent horizons, which move from $1.18\\,m_0$ and $0.78\\,m_0$ to "
        "$r_+ = 1.28\\,m_0$ and $r_- = 0.72\\,m_0$.",
        "The event horizon is the outgoing ray that ends on $r_+$, and it stands at $1.27\\,m_0$ when the tail "
        "begins, outside the apparent horizon. Every outgoing ray inside it closes on the inner apparent horizon as "
        "$v \\to \\infty$, which is the Cauchy horizon, and the ray marked is the one Ori's shell runs along, "
        "through $r = m_0$ at $v_0$. With this stream alone $m$ only settles to $m_0$, so the Kretschmann scalar "
        "$8(6m^2r^2 - 12mr_q^2r + 7r_q^4)/r^8$ stays finite toward the Cauchy horizon and diverges at $r = 0$.",
    ],
    ("mass_inflation", "ingoing", "behind"): [
        "The plane of $v$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$) behind Ori's shell, drawn with $v - r$ as the "
        "vertical axis, for $r_q = 0.96\\,m_0$. Here $v$ is the advanced time of this side, which reaches the Cauchy "
        "horizon at $v = 0$, and the shell is the edge of the chart, falling from $r = m_0$ where the tail begins to "
        "$r_- = 0.72\\,m_0$ on the Cauchy horizon. The mass function is $m_0$ until the tail arrives at "
        "$v = -3.75\\,m_0$ and $1.04\\,m_0$ soon after.",
        "It passes $2\\,m_0$ only $5 \\times 10^{-11}\\,m_0$ before the Cauchy horizon and $100\\,m_0$ at "
        "$4 \\times 10^{-14}\\,m_0$, and grows without limit as $1/(|v|\\ln^{12}(m_0/|v|))$. That growth is "
        "integrable in $v$, so every outgoing ray arrives on the Cauchy horizon at a radius above zero: the ray that "
        "starts at $0.9\\,m_0$ where the tail begins arrives at $0.67\\,m_0$. The inner apparent horizon, "
        "$g^{rr} = 0$, shrinks from $0.72\\,m_0$ toward $r = 0$ as the mass grows, and inside it the cones point to "
        "larger $r$, away from the singularity at $r = 0$.",
    ],
    ("rn_metric", "spherical", "radial"): [
        "The plane of $t$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$), drawn for $r_q = 0.48\\,r_s$, "
        "the same at every other angle by spherical symmetry. There "
        "$g^{rr}$ vanishes twice, at $r_\\pm = (r_s \\pm \\sqrt{r_s^2 - 4r_q^2})/2$, which is "
        "$0.64\\,r_s$ and $0.36\\,r_s$, and the cones close at both. Between them $r$ is the time "
        "and the cones point to smaller $r$. Inside $r_-$, $t$ is a time again.",
        "The chart alone does not fix which way is future in the two inner regions. An ingoing "
        "chart runs smoothly through both horizons, and we take the future from it, which makes "
        "the region between the horizons the black hole. The Kretschmann scalar diverges at $r = "
        "0$.",
    ],
    ("majumdar_papapetrou", "cartesian", "tz"): [
        "The plane of $t$ and $z$ on the axis through both holes ($x = y = 0$), which light launched along "
        "the axis never leaves, since $U$ is symmetric about it. Its rays are null geodesics with "
        "$dz/dt = \\pm c/U^2$. Near a hole $U \\approx m/|z \\mp 2m|$, so a ray slows as $(z \\mp 2m)^2/m^2$ and "
        "reaches the horizon only as $t \\to \\pm\\infty$; each horizon is a sphere of area $4\\pi m^2$ at the "
        "single coordinate point $z = \\pm 2m$, where $g^{zz} = 1/U^2$ vanishes.",
        "Between the holes $U \\ge 2$, with its least value at $z = 0$, so light crosses the gap no faster "
        "than $c/4$ in $t$. The Kretschmann scalar stays finite on the whole axis, the horizons included.",
    ],
    ("majumdar_papapetrou", "cartesian", "tx"): [
        "The plane of $t$ and $x$ midway between the holes ($y = z = 0$), which light launched in it never "
        "leaves, by the reflection $z \\to -z$ and the rotation about the axis. There "
        "$U = 1 + 2m/\\sqrt{x^2 + 4m^2}$ is greatest on the axis, where it is $2$, so the cones are narrowest "
        "there, with $dx/dt = \\pm c/4$, and open toward $45°$ far out.",
    ],
    ("majumdar_papapetrou", "cylindrical", "radial"): [
        "The plane of $t$ and $\\rho$ midway between the holes ($z = 0$, $\\phi = 0$), the same at every "
        "$\\phi$ by the symmetry about the axis. There $U = 1 + 2m/\\sqrt{\\rho^2 + 4m^2}$ and "
        "$d\\rho/dt = \\pm c/U^2$, which is $\\pm c/4$ on the axis and tends to $\\pm c$ far out.",
    ],
    ("majumdar_papapetrou", "isotropic", "radial"): [
        "The plane of $t$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$) about a single hole, the same at every "
        "other angle by spherical symmetry. There $dt/dr = \\pm(1 + m/r)^2/c$ diverges at $r = 0$, where "
        "$g^{rr} = r^2/(r + m)^2$ vanishes, so the cones close there and an ingoing ray reaches $r = 0$ only "
        "as $t \\to +\\infty$.",
        "$r = 0$ is the horizon, a sphere of areal radius $m$, since $g_{\\theta\\theta} = (r + m)^2$. Beyond it "
        "lies the interior of the extremal Reissner-Nordström black hole, $0 < R < m$ in the areal radius "
        "$R = r + m$. The Kretschmann scalar $8m^2(6r^2 + m^2)/(r + m)^8$ is finite at the horizon.",
    ],
    ("israel_wilson_perjes", "cylindrical", "midplane"): [
        "The plane of $t$ and $\\rho$ midway between the two sources ($z = 0$), with $\\phi$ divided out. Its "
        "rays are the null geodesics of no angular momentum, which the reflection $z \\to -z$ keeps in the "
        "plane, and each turns about the axis at $d\\phi/dt = c\\,\\omega/(W^4\\rho^2 - \\omega^2)$ as it "
        "moves in $\\rho$ at $d\\rho/dt = \\pm c\\,\\rho/\\sqrt{W^4\\rho^2 - \\omega^2}$.",
        "There $W = 1 + 2m/s$ and $\\omega = 4m^2/s + 4m^3/s^2$ with $s = \\sqrt{\\rho^2 + 4m^2}$. Far out "
        "$\\omega$ falls as $4m^2/\\rho$, the dragging of an angular momentum $l$ times the separation $4m$ "
        "of the sources, and the cones open toward $45°$. The drawing begins at $\\rho = m$, and inside "
        "$\\rho = 0.73\\,m$, where $W^4\\rho^2 = \\omega^2$, the circles of $\\phi$ are timelike.",
    ],
    ("israel_wilson_perjes", "spheroidal", "axis"): [
        "The plane of $t$ and $r$ on the axis of rotation ($\\theta = 0$), which passes through the middle of "
        "the ring, drawn for $a = m$. Its rays are null geodesics with "
        "$dr/dt = \\pm c\\,(r^2 + a^2)/((r + m)^2 + a^2)$, which never vanishes, so no cone closes and every "
        "ray runs the whole length of the axis: the source has no horizon.",
        "Through the ring the axis runs on to negative $r$, a second asymptotically flat region, where the "
        "mass is $-m$. The Kretschmann scalar is finite all along the axis, since the ring singularity, "
        "$r = -m$, lies on the equator.",
    ],
    ("israel_wilson_perjes", "spheroidal", "principal"): [
        "The equatorial plane ($\\theta = \\pi/2$) drawn in $t$ and $r$, with $\\phi$ left out, for "
        "$a = m$. Its rays are the principal null congruence, along the two repeated principal null "
        "directions of the Weyl tensor, and they are null geodesics. Each turns at "
        "$d\\phi/dr = \\pm a/(r^2 + a^2)$ as it goes, and the curves drawn are the rays' projections, "
        "$d(ct)/dr = \\pm((r + m)^2 + a^2)/(r^2 + a^2)$, the same at every $\\theta$.",
        "No cone closes anywhere. The ingoing rays cross $r = 0$, the ring of radius $a$ in the flat "
        "background where $U$ is infinite and $\\partial_t$ is null, and run on to the ring singularity at "
        "$r = -m$, which no horizon hides. The drawing begins at $r = -m/2$, and inside $r = -0.60\\,m$ the "
        "circles of $\\phi$ are timelike.",
    ],
    ("israel_wilson_perjes", "spherical", "radial"): [
        "The plane of $t$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$), drawn for $l = m/2$. There "
        "$dr/dt = \\pm c\\,(r - m)^2/(r^2 + l^2)$, so the cones close at $r = m$, where $g^{rr}$ has a double "
        "zero, and an ingoing ray reaches it only as $t \\to +\\infty$.",
        "$r = m$ is a degenerate horizon of area $4\\pi(m^2 + l^2)$. The twist sits in $g_{t\\phi}$, which "
        "drops out of the metric on this plane, and no Christoffel symbol turns these null curves out of "
        "it, so they are null geodesics, the paths light takes.",
    ],
    ("kastor_traschen", "cartesian", "tz"): [
        "The plane of $\\tau$ and $z$ on the axis through both holes ($x = y = 0$), which light launched along the "
        "axis never leaves, since $V$ is symmetric about it. Its rays are null geodesics with "
        "$dz/d\\tau = \\pm c/U^2$, $U = H\\tau + V$, and $\\tau$ is a time everywhere, so no cone closes. Each hole "
        "is the single coordinate point $z = \\pm 2m$, where $g^{zz} = 1/U^2$ vanishes, an infinitely long throat "
        "whose spheres close on the area $4\\pi m^2$.",
        "Since $H < 0$, $U$ falls at every place as $\\tau$ grows, and it reaches zero on the marked curve "
        "$c\\tau = (32m/3)V$, where the Kretschmann scalar diverges. That singularity comes in from $z = \\pm\\infty$ "
        "after $\\tau = 0$ and reaches the midpoint of the pair at $c\\tau = 32m/3$. Outside the pair, the marked ray "
        "on each side is the last along the axis to reach infinity, the event horizon there: it leaves its hole "
        "in the far past and runs out along $c\\tau z = -7.11\\,m^2$, where the areal radius "
        "$H\\tau z + 2m$ is $8m/3$, the horizon of the lukewarm hole of mass parameter $2m$.",
    ],
    ("kastor_traschen", "cartesian", "tx"): [
        "The plane of $\\tau$ and $x$ midway between the holes ($y = z = 0$), which light launched in it never "
        "leaves, by the reflection $z \\to -z$ and the rotation about the axis. There "
        "$U = H\\tau + 2m/\\sqrt{x^2 + 4m^2}$ and $dx/d\\tau = \\pm c/U^2$, so the cones open as $\\tau$ grows and "
        "$U$ falls, up to the singularity $U = 0$ on the marked curve $c\\tau = (64m/3)/\\sqrt{x^2/m^2 + 4}$.",
        "The two marked rays are the last in the plane to reach infinity, and they cross on the axis at "
        "$c\\tau = -6.20\\,m$. Before that event a ray from any point of the plane escapes, so the plane lies "
        "outside both holes. From it on, the event horizon cuts the plane in the circle the two rays trace, "
        "which is where the two horizons have joined into one, and the circle grows along "
        "$c\\tau|x| = -7.11\\,m^2$, the horizon of the lukewarm hole of mass parameter $2m$.",
    ],
    ("kastor_traschen", "cylindrical", "radial"): [
        "The plane of $\\tau$ and $\\rho$ midway between the holes ($z = 0$, $\\phi = 0$), the same at every "
        "$\\phi$ by the symmetry about the axis. There $U = H\\tau + 2m/\\sqrt{\\rho^2 + 4m^2}$ and "
        "$d\\rho/d\\tau = \\pm c/U^2$, and the singularity $U = 0$ is the marked curve "
        "$c\\tau = (64m/3)/\\sqrt{\\rho^2/m^2 + 4}$.",
        "The marked ray is the last from the axis to reach infinity. It leaves $\\rho = 0$ at $c\\tau = -6.20\\,m$, "
        "the event at which the horizons of the two holes join, and from then on it is the circle in which the "
        "one event horizon cuts the plane.",
    ],
    ("kastor_traschen", "isotropic", "radial"): [
        "The plane of $\\tau$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$) about a single hole, drawn for $H = -3c/(16m)$, "
        "the same at every other angle by spherical symmetry. The rays obey $dr/d\\tau = \\pm c/U^2$ with "
        "$U = H\\tau + m/r$, so $\\tau$ is a time everywhere and no cone closes. The marked curves, where "
        "$|\\nabla R|^2$ vanishes for the areal radius $R = rU = H\\tau r + m$, are the three horizons of the lukewarm "
        "hole, the hyperbolas $c\\tau r = -16\\,m^2$ of the cosmological horizon $R = 4m$, $-16\\,m^2/9$ of the "
        "black hole horizon $R = 4m/3$, and $0.74\\,m^2$ of the inner horizon $R = 0.86\\,m$.",
        "With $H < 0$ the chart crosses the cosmological horizon inward from the contracting region and the "
        "black hole horizon into the black hole, as the ingoing Eddington-Finkelstein chart does. The line "
        "$\\tau = 0$ is the sphere $R = m$ inside the black hole, and above it the chart runs on through the inner "
        "horizon to the singularity $R = 0$ on $c\\tau r = 16\\,m^2/3$, where the Kretschmann scalar diverges.",
    ],
    ("kastor_traschen", "comoving", "tz"): [
        "The plane of $t$ and $z$ on the axis through both holes ($x = y = 0$), whose rays are null geodesics with "
        "$dz/dt = \\pm c/(a\\Omega^2)$. The chart is the part $\\tau < 0$ of the Cartesian one, $H\\tau = e^{Ht}$, so "
        "it ends before the singularity, and $t \\to \\infty$ is the moment $\\tau = 0$. Each hole is the single "
        "coordinate point $z = \\pm 2m$.",
        "With $H < 0$ the scale factor $a = e^{Ht}$ shrinks, and every proper distance between comoving points "
        "shrinks with it. The marked ray on each side is the last along the axis to reach infinity, the event "
        "horizon there, which runs out along $|z| = (2m/3)e^{-Ht}$.",
    ],
    ("kastor_traschen", "comoving", "tx"): [
        "The plane of $t$ and $x$ midway between the holes ($y = z = 0$), which light launched in it never leaves. "
        "There $\\Omega = 1 + 2m/(a\\sqrt{x^2 + 4m^2})$ and $dx/dt = \\pm c/(a\\Omega^2)$.",
        "The two marked rays are the last in the plane to reach infinity. They cross on the axis at "
        "$ct = 5.79\\,m$, the event at which the horizons of the two holes join, and from then on they trace the "
        "circle in which the one event horizon cuts the plane, which grows as $|x| = (2m/3)e^{-Ht}$.",
    ],
    ("btz", "stationary", "static"): [
        "The plane of $t$ and $r$ ($\\phi = 0$) of the hole without rotation ($M = 1$, $J = 0$), the "
        "same at every $\\phi$ by circular symmetry. The cones close at the horizon $r_+ = \\sqrt{M}\\,\\ell$, "
        "where $g^{rr} = r^2/\\ell^2 - M$ vanishes and $c\\,dt/dr = \\pm\\ell^2/(r^2 - M\\ell^2)$ diverges. Far "
        "out the rays flatten, and a ray reaches $r \\to \\infty$, the conformal boundary, in a finite time, "
        "as in anti-de Sitter space.",
        "Inside $r_+$, $r$ is the time. We take the future from the ingoing Eddington-Finkelstein chart, "
        "which makes that region the black hole, every cone pointing to $r = 0$. The Kretschmann scalar is "
        "$12/\\ell^4$ at every point, and at $r = 0$ the circles of $\\phi$ shrink to zero length: $r = 0$ is "
        "a singularity in the causal structure, past which the circles would be closed timelike curves.",
    ],
    ("btz", "stationary", "rotating"): [
        "The plane of $t$ and $r$ of the rotating hole ($M = 1$, $J = 4\\ell/5$) with $\\phi$ divided out, "
        "$-N^2c^2dt^2 + dr^2/N^2$, the metric orthogonal to the circles of $\\phi$. Its null curves are the "
        "shadows on $t$ and $r$ of the null geodesics of zero angular momentum, each turning in $\\phi$ at "
        "$d\\phi/d(ct) = J/(2r^2)$, and each cone is the future cone of the directions of zero angular "
        "momentum.",
        "The cones close at both zeros of $N^2$, the horizons $r_+ = 2\\ell/\\sqrt{5}$ and $r_- = "
        "\\ell/\\sqrt{5}$. Between them $r$ is the time, and the future taken from the ingoing chart makes "
        "that region the black hole; inside $r_-$ the lines of constant $r$ are timelike again. The dotted "
        "line is the ergosurface, $g_{tt} = 0$ at $r = \\sqrt{M}\\,\\ell$, and between it and $r_+$ no "
        "observer keeps $\\phi$ fixed.",
    ],
    ("btz", "eddington_finkelstein_ingoing", "static"): [
        "The plane of $v$ and $r$ ($\\tilde\\phi = 0$) of the hole without rotation ($M = 1$, $J = 0$), "
        "drawn with $v - r$ as the vertical axis so that the ingoing rays, $v = $ const, run at 45°. The "
        "outgoing family has $dv/dr = 2\\ell^2/(r^2 - M\\ell^2)$, so it stands vertical at $r_+ = "
        "\\sqrt{M}\\,\\ell$: the horizon is an outgoing ray that stays where it is.",
        "The cones cross $r_+$ smoothly and keep tipping. Inside it both edges of every future cone point "
        "to smaller $r$, so every future directed ray ends at $r = 0$.",
    ],
    ("btz", "eddington_finkelstein_ingoing", "rotating"): [
        "The plane of $v$ and $r$ of the rotating hole ($M = 1$, $J = 4\\ell/5$) with $\\tilde\\phi$ divided "
        "out, $-N^2dv^2 + 2\\,dv\\,dr$, drawn with $v - r$ as the vertical axis. Its null curves are the "
        "shadows on $v$ and $r$ of the null geodesics of zero angular momentum. The ingoing family, $v = $ "
        "const, runs through both horizons to $r = 0$, and the outgoing family stands vertical at $r_+ = "
        "2\\ell/\\sqrt{5}$ and at $r_- = \\ell/\\sqrt{5}$.",
        "Between the horizons both edges of every future cone point to smaller $r$. Inside $r_-$ the "
        "outgoing edge turns back toward larger $r$, and the outgoing rays pile up against $r_-$, the inner "
        "horizon. The dotted line is the ergosurface, $g_{vv} = 0$ at $r = \\sqrt{M}\\,\\ell$.",
    ],
    ("btz", "eddington_finkelstein_outgoing", "static"): [
        "The plane of $u$ and $r$ ($\\tilde\\phi = 0$) of the hole without rotation ($M = 1$, $J = 0$), "
        "drawn with $u + r$ as the vertical axis so that the outgoing rays, $u = $ const, run at 45°. The "
        "retarded chart crosses the other horizon. Inside $r_+$ both edges of every future cone point to "
        "larger $r$, so that region is the white hole, which nothing from outside can enter.",
        "The chart is the time reverse of the ingoing one: its $g_{ur}$ is $-1$ where the ingoing chart's "
        "$g_{vr}$ is $+1$.",
    ],
    ("btz", "eddington_finkelstein_outgoing", "rotating"): [
        "The plane of $u$ and $r$ of the rotating hole ($M = 1$, $J = 4\\ell/5$) with $\\tilde\\phi$ divided "
        "out, $-N^2du^2 - 2\\,du\\,dr$, drawn with $u + r$ as the vertical axis. Its null curves are the "
        "shadows on $u$ and $r$ of the null geodesics of zero angular momentum. The outgoing family, $u = $ "
        "const, runs out from $r = 0$ through both horizons, and the ingoing family stands vertical at $r_+ "
        "= 2\\ell/\\sqrt{5}$ and at $r_- = \\ell/\\sqrt{5}$.",
        "Between the horizons both edges of every future cone point to larger $r$, the white hole. The "
        "dotted line is the ergosurface, $g_{uu} = 0$ at $r = \\sqrt{M}\\,\\ell$.",
    ],
    ("taub_nut", "spherical", "radial"): [
        "The plane of $t$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$), drawn for $l = "
        "m/2$. There $g^{rr}$ vanishes at $r_+ = m + \\sqrt{m^2 + l^2}$, about $2.118\\,m$, and "
        "inside it is the Taub region, where $r$ is the time and the cones, following the ingoing "
        "family, point to smaller $r$.",
        "The NUT parameter's twist sits in the cross term $g_{t\\phi}$, which drops out of the "
        "metric on this plane. No Christoffel symbol turns these null curves out of the plane "
        "either, so they are null geodesics, the paths light takes, even though the solution is "
        "not spherically symmetric.",
    ],
    ("bertotti_robinson", "static", "radial"): [
        "The plane of $t$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$), the AdS₂ factor of the "
        "product, $-(r^2/b^2)\\,c^2dt^2 + (b^2/r^2)\\,dr^2$. With $dt/dr = \\pm b^2/r^2$ the rays "
        "take infinite $t$ to reach $r = 0$, where $g^{rr} = r^2/b^2$ vanishes. That is a Poincaré "
        "horizon, the kind of edge that bounds the Poincaré chart of anti-de Sitter space. The "
        "other factor is the 2-sphere of radius $b$, the same everywhere, each point in the plane "
        "one such 2-sphere.",
    ],
    ("bertotti_robinson", "poincare", "tx"): [
        "The plane of $t$ and $x$ ($\\theta = \\pi/2$, $\\phi = 0$) in the Poincaré chart, where "
        "the spacetime is the product of $(b^2/x^2)(-c^2dt^2 + dx^2)$ with a 2-sphere of radius "
        "$b$. The first factor is conformal to flat, so the rays are at 45°. "
        "Here $x$ is a coordinate on the AdS₂ factor, with the boundary at $x \\to 0$ and the "
        "Poincaré horizon at $x \\to \\infty$.",
    ],
    ("near_horizon_extreme_kerr", "poincare", "equator"): [
        "The plane of $t$ and $r$ on the equator ($\\theta = \\pi/2$) with $\\phi$ divided out, "
        "$\\tfrac{1}{2}\\left(-(r^2/r_0^2)\\,c^2dt^2 + (r_0^2/r^2)\\,dr^2\\right)$, the metric orthogonal to the "
        "circles of $\\phi$ and half the Poincaré chart of a two dimensional anti-de Sitter space. Its null "
        "curves are the shadows on $t$ and $r$ of the null geodesics of zero angular momentum, each turning in "
        "$\\phi$ at $d\\phi/d(ct) = -r/r_0^2$, and each cone is the future cone of the directions of zero "
        "angular momentum.",
        "$ct \\pm r_0^2/r$ is constant along a ray, so a ray takes infinite $t$ to reach $r = 0$, a Poincaré "
        "horizon as in the Bertotti-Robinson throat, and reaches $r \\to \\infty$, the timelike boundary, in a "
        "finite one. The same rays stand over the plane of $t$ and $r$ at every $\\theta$: they are the two "
        "principal null directions of the throat. On the equator $g_{tt} = 3r^2/2r_0^2$ is positive at every "
        "$r$, so no observer there keeps $\\phi$ fixed.",
    ],
    ("near_horizon_extreme_kerr", "inverse_radius", "equator"): [
        "The plane of $t$ and $x$ on the equator ($\\theta = \\pi/2$) with $\\phi$ divided out, "
        "$(r_0^2/2x^2)\\left(-c^2dt^2 + dx^2\\right)$, the metric orthogonal to the circles of $\\phi$. It is "
        "conformal to flat, so the shadows on $t$ and $x$ of the null geodesics of zero angular momentum run "
        "at 45°, each turning in $\\phi$ at $d\\phi/d(ct) = -1/x$, and each cone is the future cone of the "
        "directions of zero angular momentum.",
        "The timelike boundary of the throat is $x \\to 0$, which a ray reaches in a finite $t$, and the "
        "Poincaré horizon is $x \\to \\infty$. The boundary is where the throat joined the rest of the extreme "
        "Kerr black hole before the limit was taken.",
    ],
    ("near_horizon_extreme_kerr", "global", "equator"): [
        "The plane of $\\tau$ and $y$ on the equator ($\\theta = \\pi/2$) with $\\phi$ divided out, "
        "$\\tfrac{1}{2}r_0^2\\left(-(1 + y^2)\\,d\\tau^2 + dy^2/(1 + y^2)\\right)$, the metric orthogonal to the "
        "circles of $\\phi$ and half the global chart of a two dimensional anti-de Sitter space. Its null "
        "curves are the shadows on $\\tau$ and $y$ of the null geodesics of zero angular momentum, each turning "
        "in $\\phi$ at $d\\phi/d\\tau = -y$, and each cone is the future cone of the directions of zero angular "
        "momentum.",
        "$\\tau \\mp \\arctan y$ is constant along a ray, so light crosses the whole throat, from the boundary "
        "$y \\to -\\infty$ to the boundary $y \\to \\infty$, in $\\Delta\\tau = \\pi$. The cones open as "
        "$1 + y^2$ and close at no $y$: this chart has no horizon, and it covers the whole spacetime. The "
        "dotted lines are $y = \\pm 1/\\sqrt{3}$, where $g_{\\tau\\tau} = \\tfrac{1}{2}r_0^2(3y^2 - 1)$ "
        "vanishes on the equator. Beyond them no observer keeps $\\phi$ fixed, the throat's counterpart of "
        "Kerr's ergoregion.",
    ],
    ("near_horizon_extreme_kerr", "near_nhek", "equator"): [
        "The plane of $t$ and $r$ on the equator ($\\theta = \\pi/2$) with $\\phi$ divided out, drawn for "
        "$k = r_0/2$: $\\tfrac{1}{2}\\left(-(r(r - 2k)/r_0^2)\\,c^2dt^2 + r_0^2\\,dr^2/r(r - 2k)\\right)$, the "
        "metric orthogonal to the circles of $\\phi$. Its null curves are the shadows on $t$ and $r$ of the "
        "null geodesics of zero angular momentum, each turning in $\\phi$ at $d\\phi/d(ct) = -(r - k)/r_0^2$, "
        "and each cone is the future cone of the directions of zero angular momentum.",
        "$ct \\mp r_*$ is constant along a ray, with $r_* = (r_0^2/2k)\\ln(1 - 2k/r)$. The cones close at "
        "$r = 2k$, a horizon of surface gravity $c^2k/r_0^2$ for the observers at fixed $r$, which a ray "
        "reaches only as $t \\to \\pm\\infty$. The patch is a part of the same throat, and its horizon is an "
        "acceleration horizon of the two dimensional anti-de Sitter space, as Rindler's is of flat space.",
    ],
    ("nariai", "static", "patch"): [
        "The plane of $t$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$) in the static chart, the two dimensional "
        "de Sitter factor of the product, $-(1 - \\Lambda r^2)\\,c^2dt^2 + dr^2/(1 - \\Lambda r^2)$, each point in "
        "the plane a 2-sphere of radius $1/\\sqrt{\\Lambda}$. The cones close at both edges, the horizons $r = "
        "\\pm 1/\\sqrt{\\Lambda}$, where $g^{rr} = 1 - \\Lambda r^2$ vanishes, and a ray takes infinite $t$ to "
        "reach either, $ct\\sqrt{\\Lambda} = \\pm\\mathrm{artanh}(r\\sqrt{\\Lambda}) + $ const.",
        "The patch is the same under $r \\to -r$. The Nariai universe is the limit of the Schwarzschild-de Sitter "
        "black hole as its black hole and cosmological horizons reach one radius, $1/\\sqrt{\\Lambda}$, and its two "
        "horizons are those two, at that radius and the same temperature.",
    ],
    ("nariai", "global", "circle"): [
        "The plane of $t$ and $\\chi$ ($\\theta = \\pi/2$, $\\phi = 0$) in the global chart, $-c^2dt^2 + "
        "\\cosh^2(\\sqrt{\\Lambda}\\,ct)\\,d\\chi^2/\\Lambda$, each point in the plane a 2-sphere of radius "
        "$1/\\sqrt{\\Lambda}$, with $\\chi = 0$ and $2\\pi$ one line. The cones narrow as "
        "$1/\\cosh(\\sqrt{\\Lambda}\\,ct)$ toward the future and the past, and along a ray $\\chi \\pm "
        "\\arctan(\\sinh(\\sqrt{\\Lambda}\\,ct))$ is constant.",
        "A ray therefore crosses half the circle, $\\Delta\\chi = \\pi$, between the infinite past and the "
        "infinite future, so observers at opposite points of the circle can never exchange a signal. The static "
        "patch is the region $\\sin\\chi > |\\tanh(\\sqrt{\\Lambda}\\,ct)|$, bounded by the four rays that leave $\\chi = 0$ and "
        "$\\chi = \\pi$ at $t = 0$.",
    ],
    ("interior_schwarzschild", "spherical", "radial"): [
        "The plane of $t$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$) over the whole domain of the "
        "chart ($r \\in [0, R]$), for a star with $R = 1.5\\,r_s$. The cones are "
        "narrowest at the centre, where $|g_{tt}|$ is least and the redshift greatest, and they "
        "stay open. They would close at the centre exactly when $3\\sqrt{1 - r_s/R} = 1$, which is "
        "Buchdahl's $R = 9r_s/8$. Beyond $R$ the spacetime is Schwarzschild's exterior, and the rays "
        "go on into it as they do in Schwarzschild's own chart.",
    ],
    ("interior_schwarzschild", "spherical", "through"): [
        "The line through the centre of the star in the plane $\\theta = \\pi/2$: $x = r$ "
        "on the right is $\\phi = 0$ and $x = -r$ on the left is $\\phi = \\pi$. Rays cross the "
        "centre smoothly, and the cones are narrowest there.",
    ],
    ("tolman_vii", "spherical", "radial"): [
        "The plane of $t$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$) over the whole domain of the chart "
        "($r \\in [0, R]$), for a star with $R = 2\\,r_s$ ($\\beta = 1/4$). The cones are narrowest at the "
        "centre, where $|g_{tt}| = 0.199$ and the redshift is greatest, and they stay open. They would close at "
        "the centre at $\\beta \\approx 0.3862$, where $\\psi$ reaches $\\pi/2$ there and the central pressure "
        "is infinite. Beyond $R$ the spacetime is Schwarzschild's exterior, and the rays go on into it as they do "
        "in Schwarzschild's own chart.",
    ],
    ("tolman_vii", "spherical", "through"): [
        "The line through the centre of the star in the plane $\\theta = \\pi/2$ ($R = 2\\,r_s$): $x = r$ "
        "on the right is $\\phi = 0$ and $x = -r$ on the left is $\\phi = \\pi$. Rays cross the "
        "centre smoothly, and the cones are narrowest there.",
    ],
    ("tolman_vii", "tolman", "radial"): [
        "The plane of $t$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$) for the same star in Tolman's constants "
        "($R^2 = 16r_s^2/5$, $A^4 = 256r_s^4/3$, $B^2 = 7/12$, $C \\approx 0.0782$), with its surface at "
        "$r_b = 2\\,r_s$. The phase $\\psi$ here is $\\pi/2$ less Lattimer and Prakash's, so $B^2\\sin^2\\psi$ "
        "is their $(1 - 5\\beta/3)\\cos^2\\psi$ and the rays are the same curves.",
    ],
    ("tolman_vii", "tolman", "through"): [
        "The line through the centre of the same star in Tolman's constants: $x = r$ on the right is "
        "$\\phi = 0$ and $x = -r$ on the left is $\\phi = \\pi$. Rays cross the centre smoothly, and the cones "
        "are narrowest there.",
    ],
    ("kerr", "boyer_lindquist", "radial"): [
        "The plane of $t$ and $r$ on the rotation axis ($\\theta = 0$), drawn for $a = 0.9\\,GM/c^2$. The curves "
        "drawn are null, and on the axis they are also null geodesics, the paths light takes. Off the axis a light "
        "ray launched along a curve of fixed $\\theta$ and $\\phi$ is turned out of the plane by "
        "$\\Gamma^\\theta{}_{tt}$, $\\Gamma^\\theta{}_{rr}$, and $\\Gamma^\\phi{}_{tr}$. On the axis $g^{rr} = "
        "\\Delta/\\Sigma$, with $\\Delta = r^2 - 2GMr/c^2 + a^2$ and $\\Sigma = r^2 + a^2\\cos^2\\theta$, vanishes "
        "at both roots of $\\Delta$, $r_\\pm = GM/c^2 \\pm \\sqrt{(GM/c^2)^2 - a^2}$, and the cones close at both; "
        "between them they point to smaller $r$, following the ingoing family.",
        "The domain of the chart begins at $r_+$, the outer zero of $g^{rr}$. On the axis the "
        "ergosurface touches the horizon, so the metric on the plane stays Lorentzian. The "
        "Kretschmann scalar stays finite at $r = 0$ on the axis, because the ring singularity "
        "lies in the equatorial plane.",
    ],
    ("kerr", "boyer_lindquist", "principal"): [
        "The equatorial plane ($\\theta = \\pi/2$) drawn in $t$ and $r$, with $\\phi$ left "
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
        "The equatorial plane ($\\theta = \\pi/2$) seen from above, along the axis from "
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
    ("kerr_de_sitter", "boyer_lindquist", "axis"): [
        "The plane of $t$ and $r$ on the rotation axis ($\\theta = 0$), drawn for $a = 0.45\\,r_s$ and $\\Lambda = 0.2/r_s^2$. On the axis the curves drawn are null geodesics, the paths light takes, and $g^{rr} = \\Delta_r/(r^2 + a^2)$ vanishes at $r_- = 0.279$, $r_+ = 0.785$, and $r_c = 3.232\\,r_s$. The cones close at each of them, since $d(ct)/dr = \\pm(r^2 + a^2)/\\Delta_r$ diverges there.",
        "The domain of the chart is the region between $r_+$ and $r_c$, where $t$ is a time. Inside $r_+$ and beyond $r_c$ the components alone do not fix which way is future. We take it from the ingoing Kerr chart inside $r_+$, which makes the region between $r_-$ and $r_+$ the black hole, and from the outgoing one beyond $r_c$, which makes that region the expanding universe. The Kretschmann scalar stays finite at $r = 0$ on the axis, because the ring singularity lies in the equatorial plane.",
    ],
    ("kerr_de_sitter", "boyer_lindquist", "principal"): [
        "The equatorial plane ($\\theta = \\pi/2$) drawn in $t$ and $r$, with $\\phi$ left out, for $a = 0.45\\,r_s$ and $\\Lambda = 0.2/r_s^2$. Its rays are the principal null congruence, the light rays that run straight in and straight out: the two null directions of the plane of $\\partial_r$ and $(r^2 + a^2)\\,\\partial_t + a\\Xi\\,\\partial_\\phi$, which are the repeated principal null directions of the Weyl tensor. Every ray turns as it goes, at $d\\phi/dr = \\pm a\\Xi/\\Delta_r$, and the curves drawn are the rays' projections, $d(ct)/dr = \\pm(r^2 + a^2)/\\Delta_r$, the same at every $\\theta$.",
        "The rays are null geodesics, and each cone is the future cone of the principal plane. The cones close at $r_- = 0.279$, $r_+ = 0.785$, and $r_c = 3.232\\,r_s$. The dotted lines are the two ergosurfaces, $g_{tt} = 0$, at $r = 1.105$ and $3.173\\,r_s$ on the equator. Between $r_+$ and the first nothing, light included, can keep $\\phi$ fixed against the turning of the black hole, and between the second and $r_c$ nothing can keep it fixed against the turning of the cosmological horizon. The ingoing rays end on the ring singularity at $r = 0$, which lies in this plane.",
    ],
    ("kerr_de_sitter", "boyer_lindquist", "above"): [
        "The equatorial plane ($\\theta = \\pi/2$) between the horizons seen from above, with $r$ and $\\phi$ drawn as polar coordinates and $t$ left out, for $a = 0.45\\,r_s$ and $\\Lambda = 0.2/r_s^2$. Its rays are the principal null congruence, each turning at $d\\phi/dr = \\pm a\\Xi/\\Delta_r$, so both families wind counterclockwise, the way the hole turns: the ingoing rays as they fall and the outgoing rays as they climb.",
        "At $r_+$ and at $r_c$ the angle $\\phi$ runs to infinity along every ray, as $t$ does, so the rays wind without end onto both horizons. The winding is in the coordinate $\\phi$ alone: along an ingoing ray $\\tilde\\phi = \\phi + \\int a\\Xi\\,dr/\\Delta_r$ stays fixed, and in it the ray crosses each horizon at a finite angle. The dotted circles are the ergosurfaces, $r = 1.105$ and $3.173\\,r_s$ on the equator.",
    ],
    ("kerr_de_sitter", "boyer_lindquist", "axis_ads"): [
        "The plane of $t$ and $r$ on the rotation axis ($\\theta = 0$), drawn for $r_s = 2\\ell$ and $a = \\ell/2$, with $\\Lambda = -3/\\ell^2$. On the axis the curves drawn are null geodesics, the paths light takes, and $g^{rr} = \\Delta_r/(r^2 + a^2)$ vanishes at $r_- = 0.137\\,\\ell$ and $r_+ = 0.859\\,\\ell$, where the cones close. Far out $d(ct)/dr \\to \\pm\\ell^2/r^2$, so a ray reaches $r \\to \\infty$ in a finite $t$, as in anti-de Sitter space.",
        "The domain of the chart begins at $r_+$. Inside it we take the future from the ingoing Kerr chart, which makes the region between the horizons the black hole, where every cone points to smaller $r$. The Kretschmann scalar stays finite at $r = 0$ on the axis, because the ring singularity lies in the equatorial plane.",
    ],
    ("kerr_de_sitter", "boyer_lindquist", "principal_ads"): [
        "The equatorial plane ($\\theta = \\pi/2$) drawn in $t$ and $r$, with $\\phi$ left out, for $r_s = 2\\ell$ and $a = \\ell/2$, with $\\Lambda = -3/\\ell^2$. Its rays are the principal null congruence, each turning at $d\\phi/dr = \\pm a\\Xi/\\Delta_r$, and the curves drawn are the rays' projections, $d(ct)/dr = \\pm(r^2 + a^2)/\\Delta_r$, the same at every $\\theta$.",
        "The rays are null geodesics, and each cone is the future cone of the principal plane. The cones close at $r_+ = 0.859\\,\\ell$ and point to smaller $r$ inside it. The dotted line is the ergosurface, $g_{tt} = 0$, at $r = 0.939\\,\\ell$ on the equator, and between it and $r_+$ nothing, light included, can keep $\\phi$ fixed.",
    ],
    ("kerr_de_sitter", "nonrotating", "above"): [
        "The equatorial plane ($\\theta = \\pi/2$) between the horizons seen from above, with $r$ and $\\Phi$ drawn as polar coordinates and $t$ left out, for $a = 0.45\\,r_s$ and $\\Lambda = 0.2/r_s^2$. In the frame that does not rotate as $r \\to \\infty$ the principal null rays turn at $d\\Phi/dr = \\pm a\\left(1 - \\Lambda r^2/3\\right)/\\Delta_r$, less than in Carter's coordinates by the factor $(1 - \\Lambda r^2/3)/\\Xi$, which falls from $0.946$ at $r_+$ to $0.299$ at $r_c$. The rays still wind without end onto both horizons, where $\\Delta_r = 0$.",
    ],
    ("kerr_de_sitter", "nonrotating", "above_ads"): [
        "The equatorial plane ($\\theta = \\pi/2$) outside the horizon seen from above, with $r$ and $\\Phi$ drawn as polar coordinates and $t$ left out, for $r_s = 2\\ell$ and $a = \\ell/2$, with $\\Lambda = -3/\\ell^2$. The angle $\\Phi$ is that of the frame that does not rotate at infinity, in which the conformal boundary is an Einstein universe at rest. The principal null rays turn at $d\\Phi/dr = \\pm a\\left(1 + r^2/\\ell^2\\right)/\\Delta_r$, which dies away as $a/r^2$ far out, and they wind without end onto the horizon $r_+ = 0.859\\,\\ell$. In this frame the horizon turns at $\\Omega = ac\\left(1 + r_+^2/\\ell^2\\right)/(r_+^2 + a^2) = 0.879\\,c/\\ell$, slower than light at the boundary.",
    ],
    ("kerr_de_sitter", "kerr_ingoing", "axis"): [
        "The plane of $v$ and $r$ on the rotation axis ($\\theta = 0$), drawn for $a = 0.45\\,r_s$ and $\\Lambda = 0.2/r_s^2$ with $v - r$ as the vertical axis, so that the ingoing rays, $v = $ const, run at 45°. The outgoing family has $dv/dr = 2(r^2 + a^2)/\\Delta_r$, so it stands vertical at $r_- = 0.279$, $r_+ = 0.785$, and $r_c = 3.232\\,r_s$: each horizon is an outgoing ray that stays where it is.",
        "The chart crosses $r_c$ from the contracting region in the past of the static one, and $r_+$ into the black hole, where both edges of every future cone point to smaller $r$. Inside $r_-$ the cones open again, and an ingoing ray reaches $r = 0$, the centre of the ring's disc, in a finite $v$.",
    ],
    ("kerr_de_sitter", "kerr_ingoing", "axis_ads"): [
        "The plane of $v$ and $r$ on the rotation axis ($\\theta = 0$), drawn for $r_s = 2\\ell$ and $a = \\ell/2$, with $\\Lambda = -3/\\ell^2$ with $v - r$ as the vertical axis, so that the ingoing rays, $v = $ const, run at 45°. The outgoing family has $dv/dr = 2(r^2 + a^2)/\\Delta_r$ and stands vertical at $r_- = 0.137\\,\\ell$ and $r_+ = 0.859\\,\\ell$. Between the horizons both edges of every future cone point to smaller $r$, and inside $r_-$ the cones open again.",
    ],
    ("kerr_de_sitter", "kerr_outgoing", "axis"): [
        "The plane of $u$ and $r$ on the rotation axis ($\\theta = 0$), drawn for $a = 0.45\\,r_s$ and $\\Lambda = 0.2/r_s^2$ with $u + r$ as the vertical axis, so that the outgoing rays, $u = $ const, run at 45°. The ingoing family stands vertical at $r_-$, $r_+$, and $r_c$. Between $r_-$ and $r_+$ both edges of every future cone point to larger $r$: this is the white hole, which nothing from outside can enter. Beyond $r_c$ they point to larger $r$ as well, into the expanding region.",
    ],
    ("kerr_de_sitter", "kerr_outgoing", "axis_ads"): [
        "The plane of $u$ and $r$ on the rotation axis ($\\theta = 0$), drawn for $r_s = 2\\ell$ and $a = \\ell/2$, with $\\Lambda = -3/\\ell^2$ with $u + r$ as the vertical axis, so that the outgoing rays, $u = $ const, run at 45°. The ingoing family stands vertical at $r_- = 0.137\\,\\ell$ and $r_+ = 0.859\\,\\ell$, and between them both edges of every future cone point to larger $r$, out of the white hole.",
    ],
    ("kerr_de_sitter", "kerr_schild", "axis"): [
        "The plane of $\\tau$ and $r$ on the rotation axis ($\\theta = 0$), drawn for $a = 0.45\\,r_s$ and $\\Lambda = 0.2/r_s^2$. On the axis the ingoing rays have $d(c\\tau)/dr = -1/(1 - \\Lambda r^2/3)$, as in de Sitter space, whatever the mass, since the null vector of the Kerr-Schild form is tangent to them. The outgoing rays stand vertical at $r_-$, $r_+$, and $r_c$, and the chart runs through all three.",
        "The chart ends at $r = \\sqrt{3/\\Lambda} = 3.873\\,r_s$, the horizon of the de Sitter background, which lies beyond $r_c = 3.232\\,r_s$.",
    ],
    ("kerr_de_sitter", "kerr_schild", "axis_ads"): [
        "The plane of $\\tau$ and $r$ on the rotation axis ($\\theta = 0$), drawn for $r_s = 2\\ell$ and $a = \\ell/2$, with $\\Lambda = -3/\\ell^2$. On the axis the ingoing rays have $d(c\\tau)/dr = -1/(1 + r^2/\\ell^2)$, as in anti-de Sitter space, whatever the mass, since the null vector of the Kerr-Schild form is tangent to them. The outgoing rays stand vertical at $r_- = 0.137\\,\\ell$ and $r_+ = 0.859\\,\\ell$, and the chart runs through both.",
    ],
    ("kerr_newman", "boyer_lindquist", "radial"): [
        "The plane of $t$ and $r$ on the rotation axis ($\\theta = 0$), drawn for $a = "
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
        "The equatorial plane ($\\theta = \\pi/2$) drawn in $t$ and $r$, with $\\phi$ left "
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
        "The equatorial plane ($\\theta = \\pi/2$) seen from above, along the axis from "
        "$\\theta = 0$, with $r$ and $\\phi$ drawn as polar coordinates and $t$ left out, for "
        "$a = 0.6\\,GM/c^2$ and $r_Q = 0.5\\,GM/c^2$. Its rays are the principal null congruence, "
        "each turning at $d\\phi/dr = \\pm a/\\Delta$ with $\\Delta = r^2 - 2GMr/c^2 + a^2 + r_Q^2$, "
        "so both families wind counterclockwise, the way the hole turns, and far out the turning "
        "dies away as $a/r^2$. The charge enters only through $\\Delta$: it pulls the horizon in "
        "to $r_+ = 1.624\\,GM/c^2$, and the rays wind onto it without end, as $\\phi$ runs to "
        "infinity there.",
        "The winding is in the coordinate $\\phi$ alone: along an ingoing ray "
        "$\\tilde\\phi = \\phi + \\int a\\,dr/\\Delta$ stays fixed, and in it the ray crosses the "
        "horizon at a finite angle. The dotted circles are the ergosurfaces, "
        "$r = GM/c^2 \\pm \\sqrt{(GM/c^2)^2 - r_Q^2}$ on the equator, and between the outer one "
        "and $r_+$ nothing, light included, can keep $\\phi$ fixed.",
    ],
    ("kasner", "cartesian", "tx"): [
        "The plane of $t$ and $x$ ($y = z = 0$). Along $x$ the scale factor $t^{-2/7}$ "
        "grows toward the singularity, so the cones close up as $t \\to 0$: $dx/dt = \\pm "
        "t^{2/7}$. The Kretschmann scalar $-16p_1p_2p_3/t^4$ diverges at $t = 0$, the "
        "singularity.",
    ],
    ("kasner", "cartesian", "tz"): [
        "The plane of $t$ and $z$ ($x = y = 0$). Along $z$ the scale factor $t^{6/7}$ "
        "goes to zero at the singularity, and the cones open out flat: $dz/dt = \\pm t^{-6/7}$. "
        "In the plane of $t$ and $x$ the same singularity closes the cones, since the "
        "scale factor there, $t^{-2/7}$, grows as $t \\to 0$.",
    ],
    ("kantowski_sachs", "comoving", "tr"): [
        "The plane of $t$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$), each point in the diagram a 2-sphere of "
        "radius $b(t)$. The line element leaves $a(t)$ and $b(t)$ free, and here they are those of dust at rest "
        "in the chart, solved from $G^r{}_r = G^\\theta{}_\\theta = 0$ with both at rest at the dashed line, "
        "where the spheres are largest. The universe lasts $\\pi b_0/c$ from one singularity to the other, "
        "and $t$ is counted from the first.",
        "The edges of the cones are $dr/d(ct) = \\pm 1/a$. At both singularities $b \\to 0$ while $a$ grows "
        "without bound, so the cones close up along $r$, and a ray crosses only $2.44\\,b_0$ of $r$ in the "
        "whole life of the universe.",
    ],
    ("kantowski_sachs", "dust", "etar"): [
        "The plane of $\\eta$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$) of the dust universe symmetric in time "
        "($\\kappa = 0$), each point in the diagram a 2-sphere of radius $b_0\\cos^2\\eta$. The edges of the "
        "cones are $dr/d\\eta = \\pm 2b_0\\cos^2\\eta/(1 + \\eta\\tan\\eta)$, widest at $\\eta = 0$, "
        "where the spheres are largest.",
        "The cones close toward $\\eta = \\pm\\pi/2$, where the spheres shrink to nothing, the lengths along "
        "$r$ grow without bound, and the Kretschmann scalar diverges. No Christoffel symbol turns a ray out of "
        "the plane, so every curve drawn is a null geodesic.",
    ],
    ("kantowski_sachs", "schwarzschild_interior", "Tr"): [
        "The plane of $T$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$) inside the horizon of a Schwarzschild "
        "black hole, each point in the diagram a 2-sphere of radius $T$. The future lies toward smaller $T$, so "
        "every cone points down, and its edges are $dr/dT = \\pm T/(r_s - T)$, so that "
        "$r \\pm \\left(T + r_s\\ln(1 - T/r_s)\\right)$ is constant along a ray.",
        "The cones lie flat at the horizon, $T = r_s$, which a ray leaves at any $r$, and close toward the "
        "singularity $T = 0$, where the Kretschmann scalar $12r_s^2/T^6$ diverges. From $T = r_s$ to $T = 0$ a "
        "ray's $r$ changes without bound near the horizon and by less and less near the singularity, so two "
        "observers at different $r$ lose sight of each other before the end.",
    ],
    ("bianchi", "type_i_cartesian", "tx"): [
        "The plane of $t$ and $x$ ($y = z = 0$). The singularity comes about $0.378/\\bar "
        "H$ before the dashed line, and $t$ is counted from it. Near its singularity the dust "
        "universe is a Kasner spacetime, its scale factors running as powers of $t$ whose "
        "exponents lie on the Kasner circle.",
        "Along $x$ the cones therefore close toward $t = 0$, as in Kasner's contracting direction, "
        "and open again later as $a_1$ turns round.",
    ],
    ("godel", "cartesian", "tx"): [
        "The plane of $t$ and $x$ ($y = z = 0$). The Gödel universe has no time function: "
        "its $g^{tt} = 2\\omega^2$ is positive, so the surfaces $t = $ const are not "
        "spacelike, and the cones are oriented by $\\partial_t$, which is timelike "
        "everywhere.",
        "The metric on this plane is $(-dt^2 + dx^2)/2\\omega^2$, so the curves drawn are null and "
        "run at 45°. They are not null geodesics, though. Since $\\Gamma^y{}_{tx} = -e^{-x}$ is not "
        "zero, a light ray launched along one of these curves is turned out of the plane into $y$. "
        "The closed timelike curves of the Gödel universe circle each world line of the dust beyond "
        "a critical radius, through $y$ as well as $x$, so they cross this plane.",
    ],
    ("alcubierre", "cartesian", "tx"): [
        "The plane of $t$ and $x$ on the bubble's axis of motion ($y = z = 0$), for a "
        "bubble moving at twice the speed of light along $x = 2ct$. There $\\partial_y f = "
        "\\partial_z f = 0$, so the null curves drawn are null geodesics, the paths light takes. "
        "Inside the bubble and far outside, the cones are Minkowski's; in the walls they tilt with "
        "the shift $v_s f$.",
        "The dashed lines are $g^{xx} = 1 - v_s^2 f^2 = 0$, at $v_s f = 1$. For $v_s = 2$ that is "
        "also where a forward ray keeps pace with the bubble, $f = 1 - 1/v_s$, so here they are "
        "the two horizons: forward rays from inside stall at the front wall, and forward rays from "
        "behind stall at the back wall. At other speeds the two places differ.",
    ],
    ("van_den_broeck", "cartesian", "tx"): [
        "The plane of $t$ and $x$ on the bubble's axis of motion ($y = z = 0$), for a bubble moving at twice "
        "the speed of light along $x = 2ct$. On the axis $\\partial_y$ and $\\partial_z$ of $f$ and $B$ "
        "vanish, so the null curves drawn are null geodesics, the paths light takes. Their slopes are "
        "$dx/d(ct) = v_s f \\pm 1/B$: in the wall the cones tilt with the shift $v_s f$, as Alcubierre's do, "
        "and through the neck they narrow about the world line of the ship, since a unit of $x$ there is "
        "$B$ units of distance.",
        "The dashed lines are $g^{xx} = (1 - v_s^2 f^2B^2)/B^2 = 0$. Outside the neck $B = 1$, so they "
        "stand where $v_s f = 1$, and for $v_s = 2$ a forward ray keeps pace with the bubble there. "
        "The pocket itself is the sliver $|x - 2ct| < 0.004\\,R$ about the centre, $3R$ across for "
        "whoever is inside it.",
    ],
    ("van_den_broeck", "pocket", "radial"): [
        "The plane of $t$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$) inside the bubble, in the chart that "
        "rides with the ship, which ends where the wall of the bubble begins, at $r = R$. The radial null "
        "geodesics are $c\\,dt = \\pm B\\,dr$. From the wall in to "
        "$r = R/2$ they run at 45°, and through the neck the cones close about the vertical as $B$ climbs "
        "to $380$.",
        "A ray takes $4R/c$ from the middle of the neck to the centre, a proper distance of $4R$, and "
        "nearly all of it is spent inside $r = 0.02\\,R$, where the lid, the rim, and the floor of the "
        "pocket lie.",
    ],
    ("van_den_broeck", "proper_radial", "radial"): [
        "The plane of $t$ and the proper distance $l$ ($\\theta = \\pi/2$, $\\phi = 0$), from the centre "
        "of the pocket at $l = -4R$ out to the wall of the bubble. Here $g_{tt} = -1$ and $g_{ll} = 1$, "
        "so every radial null geodesic runs at 45° whatever $r(l)$ is, and the pocket shows in "
        "$g_{\\theta\\theta} = r^2$ alone.",
        "The two lines marked are the spheres where $dr/dl = 0$: the widest, $r = 7R/4$ at $l = -2R$ on "
        "the rim of the pocket, and the narrowest, $r = R/4$ at $l = 0$ in the neck. Between them the "
        "spheres shrink as $l$ grows.",
    ],
    ("natario", "cartesian_flow", "tx"): [
        "The plane of $t$ and $x$ on the axis of motion ($y = z = 0$). There the field "
        "reduces to $u = 2nv_s = v_s f$, which is Alcubierre's shift, so on this plane the two "
        "metrics are the same, and so are their light rays.",
        "The drives differ only off the axis, where Natário's flow slides space sideways with no "
        "compression.",
    ],
    ("krasnikov", "cylindrical", "tx"): [
        "The plane of $t$ and $x$ along the axis of the tube ($r = 0$). Outside the tube $k "
        "= 1$ and the cones are Minkowski's. Inside, $k$ comes close to $\\delta - 1$ and the edge "
        "moving left tips below the horizontal, so a ray going back toward $x = 0$ loses about "
        "$0.8$ in $ct$ for every unit of $x$ it covers. The edge moving right, $c\\,dt = dx$, is "
        "the same everywhere, so it points to the future inside the tube as it does outside, where "
        "$t$ is a time.",
        "The dash dot line is $g^{tt} = 0$, where $k = 0$ and the surfaces $t = $ const stop being "
        "spacelike. On the axis the expression for the Kretschmann scalar is $0/0$; the tube's "
        "curvature is concentrated in its thin walls.",
    ],
    ("khan_penrose", "double_null", "plane"): [
        "The plane of $u$ and $v$ ($x = y = 0$) where both waves have passed, drawn with $u + v$ up and $v - u$ "
        "across, each point in the diagram a single event. Only $g_{uv}$ is nonzero on it, so the light rays are "
        "the lines $u = $ const and $v = $ const at 45°, and no Christoffel symbol turns a ray along them out of the "
        "plane, so each is a null geodesic.",
        "The region is bounded below by the fronts of the two impulsive waves, $u = 0$ and $v = 0$, which met at "
        "$u = v = 0$ and on which the Riemann tensor has a delta singularity, and above by the spacelike curvature "
        "singularity $u^2 + v^2 = 1$, where the Kretschmann scalar diverges. Every ray and every observer in the "
        "region ends on it. Its two ends, $u = 1$ on $v = 0$ and $v = 1$ on $u = 0$, are where the fold "
        "singularities behind each wave alone meet it.",
    ],
    ("khan_penrose", "cosmological", "plane"): [
        "The plane of $\\tau$ and $\\sigma$ ($x = y = 0$) where both waves have passed, each point in the diagram "
        "a single event. The metric on it is $L^2(\\cos\\tau)^{3/2}(-d\\tau^2 + d\\sigma^2)/(2\\sqrt{\\cos\\sigma})$, "
        "a multiple of Minkowski's, so the light rays are the lines of constant $\\tau + \\sigma = 2\\arcsin u$ "
        "and constant $\\tau - \\sigma = 2\\arcsin v$, at 45°, and each is a null geodesic.",
        "The fronts of the two waves are the lines $\\sigma = \\pm\\tau$, which leave the collision at "
        "$\\tau = \\sigma = 0$, and the curvature singularity $u^2 + v^2 = 1$ is the line $\\tau = \\pi/2$, where "
        "$g_{xx}$ grows as $1/\\cos\\tau$ and $g_{yy}$ falls as $\\cos^3\\tau$. Each surface of constant "
        "$\\tau$ is spacelike, and every observer in the region reaches $\\tau = \\pi/2$.",
    ],
    ("bell_szekeres", "double_null", "plane"): [
        "The plane of $u$ and $v$ ($x = y = 0$) where both waves have passed, drawn with $u + v$ up and $v - u$ "
        "across at $a = b$, each point in the diagram a single event. Only $g_{uv}$ is nonzero on it, so the light "
        "rays are the lines $u = $ const and $v = $ const at 45°, and no Christoffel symbol turns a ray along them "
        "out of the plane, so each is a null geodesic.",
        "The region is bounded below by the fronts of the two electromagnetic shock waves, $u = 0$ and $v = 0$, "
        "which met at $u = v = 0$ and which carry the impulsive gravitational waves made by the collision, and "
        "above by the Killing-Cauchy horizon $au + bv = \\pi/2$, where $g_{yy} = \\cos^2(au + bv)$ vanishes and "
        "the Kretschmann scalar is $32a^2b^2$ as everywhere else. Every ray in the region reaches the horizon. "
        "Its two ends, $au = \\pi/2$ on $v = 0$ and $bv = \\pi/2$ on $u = 0$, are where the fold singularities "
        "behind each wave alone meet it.",
    ],
    ("bell_szekeres", "time_space", "plane"): [
        "The plane of $\\xi$ and $\\eta$ ($x = y = 0$) where both waves have passed, each point in the diagram a "
        "single event. The metric on it is $(-d\\xi^2 + d\\eta^2)/(2ab)$, flat, so the light rays are the lines "
        "of constant $\\xi - \\eta = 2au$ and constant $\\xi + \\eta = 2bv$, at 45°, and each is a null "
        "geodesic.",
        "The fronts of the two waves are the lines $\\eta = \\pm\\xi$, which leave the collision at "
        "$\\xi = \\eta = 0$, and the Killing-Cauchy horizon is the line $\\xi = \\pi/2$, where "
        "$g_{yy} = \\cos^2\\xi$ vanishes. Each surface of constant $\\xi$ is spacelike, and an observer at rest "
        "in $\\eta$, $x$, and $y$ reaches the horizon after the proper time $\\pi/(2\\sqrt{2ab})$.",
    ],
    ("bell_szekeres", "regular", "plane"): [
        "The plane of $T$ and $Z$ ($X = 0.99$, $Y = 0$) at $a = b$, each point in the diagram a single event. "
        "The surfaces of constant $\\xi = au + bv$ are the hyperbolas $T^2 - Z^2 = \\cos^2\\xi$, and the "
        "lines through the origin are the surfaces of constant $y$. No Christoffel symbol turns a ray out of "
        "the plane, so each ray is a null geodesic.",
        "The collision $\\xi = 0$ is the hyperbola $T^2 - Z^2 = 1$, and on this plane the two wave fronts lie "
        "within $0.01$ of it, on $T^2 - Z^2 = X^2$. The Killing-Cauchy horizon is the pair of lines $T = -|Z|$, "
        "and the event $T = Z = 0$ where they cross is every point of the horizon at finite $y$. The metric is "
        "regular there, and the rays run on through the horizon into $T > -|Z|$.",
    ],
    ("bell_szekeres", "global", "plane"): [
        "The plane of $\\chi$ and $\\rho$ ($\\theta = \\pi/2$, $\\phi = 0$) at $a = b$, each point in the "
        "diagram a single event. The metric on it is $(-\\cosh^2\\rho\\,d\\chi^2 + d\\rho^2)/(2ab)$, the anti-de "
        "Sitter space of two dimensions in its global chart, and the light rays are the curves of constant "
        "$\\chi \\pm \\arctan(\\sinh\\rho)$, each a null geodesic.",
        "The collision is the line $\\chi = -\\pi/2$, and the Killing-Cauchy horizon is the pair of rays "
        "$\\cos\\chi\\cosh\\rho = 1$ that climb from its two ends, $\\rho \\to \\pm\\infty$, to meet at "
        "$\\chi = \\rho = 0$. Where both waves have passed is the region between them. The chart runs on above "
        "the horizon, where Chris Clarke and Sean Hayward continued the spacetime.",
    ],
    ("bell_szekeres", "kruskal_szekeres", "plane"): [
        "The plane of $U$ and $V$ ($\\eta = 0$, $x = 0$), drawn with $U + V$ up and $V - U$ across at $a = b$, "
        "each point in the diagram a single event. Only $g_{UV}$ is nonzero on it, so the light rays are the "
        "lines $U = $ const and $V = $ const at 45°, each a null geodesic.",
        "The collision is the hyperbola $2abUV = 1$ in the quadrant $U < 0$, $V < 0$, and the Killing-Cauchy "
        "horizon is the pair of lines $U = 0$ and $V = 0$, on which the metric is regular. The surfaces of "
        "constant $\\xi = au + bv$ are the hyperbolas $2abUV = (1 - \\sin\\xi)/(1 + \\sin\\xi)$, and the "
        "lines through the origin are the surfaces of constant $y$, so the translation along $y$ acts on the "
        "plane as a boost and is null on the horizon.",
    ],
    ("bell_szekeres", "bertotti_robinson", "plane"): [
        "The plane of $t$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$), each point in the diagram a single event. "
        "The metric on it is $(-dt^2 + dr^2)/(2ab\\,r^2)$, a multiple of Minkowski's, so the light rays are the "
        "lines of constant $t \\pm r$, at 45°, and each is a null geodesic.",
        "The collision is the line $t = 0$, and where both waves have passed is the wedge $0 \\le t < r$ above "
        "it. The surfaces of constant $\\xi = au + bv$ are the lines $t = r\\sin\\xi$ through the origin, and "
        "the Killing-Cauchy horizon is the ray $t = r$, reached in this chart only as $y \\to -\\infty$; the rest "
        "of the horizon lies at $t = \\infty$.",
    ],
    ("pp_wave", "exact_plane_wave", "tz"): [
        "The plane the wave travels in, on its axis ($x = y = 0$), drawn with $u = t - z$ and $v = (t + z)/2$ for "
        "a time $t$ and a distance $z$ along the wave, so that the axes are $t$ and $z$; the chart's own $u$ and "
        "$v$ are both null. On the axis the profile $A(x^2 - y^2) + 2Bxy$ vanishes whatever $A$ and $B$ are, so "
        "the metric on this plane is flat and the rays are at 45°.",
        "Off the axis a light ray is pushed out of the plane, since $\\ddot x = (Ax + By)\\dot "
        "u^2$ and $\\ddot y = -(Ay - Bx)\\dot u^2$: the wave squeezes a beam toward the axis in "
        "one transverse direction and stretches it in the other, across this plane in $x$ and "
        "$y$. On the axis both accelerations vanish, so the null curves drawn there are null "
        "geodesics.",
    ],
    ("stockum_dust", "cylindrical", "inside"): [
        "The cylinder of $t$ and $\\phi$ ($r = R/2$, $z = 0$), opened along the line "
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
        "The cylinder of $t$ and $\\phi$ ($r = 3R/2$, $z = 0$), opened along $\\phi = "
        "\\pm\\pi$ in the same way. Beyond $r = R$ the coefficient $g_{\\phi\\phi} = r^2(1 - r^2/R^2)$ is "
        "negative, and the cones have tipped over past the horizontal: the null curve moving to $+\\phi$, "
        "$dt = r(1 - r/R)\\,d\\phi$, goes down in $t$, while the one moving to $-\\phi$, $dt = -r(1 + "
        "r/R)\\,d\\phi$, climbs steeply. Every horizontal line, run toward $+\\phi$, points into the "
        "future cones, so the circle of constant $t$, $r$, and $z$ is a closed timelike curve.",
        "The curve moving to $+\\phi$ comes round to its own $\\phi$ at a $t$ earlier by $3\\pi R/2$ after "
        "each turn. None of the curves drawn here is a null geodesic: light launched along one moving "
        "to $+\\phi$ is turned toward the axis by $\\Gamma^r{}_{t\\phi}$ and $\\Gamma^r{}_{\\phi\\phi}$, and "
        "light launched along one moving to $-\\phi$ is turned away from it.",
    ],
    ("godel", "cylindrical", "inside"): [
        "The cylinder of $t$ and $\\phi$ ($r = r_c/2$, $z = 0$, $\\sinh r_c = 1$) about the axis $r = 0$, the "
        "world line of one particle of the dust. It is opened along the "
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
        "The cylinder of $t$ and $\\phi$ ($r = 3r_c/2$, $z = 0$), opened along $\\phi = \\pm\\pi$ "
        "in the same way. Beyond $r_c$ the coefficient $g_{\\phi\\phi} = 2\\sinh^2 r\\,(1 - \\sinh^2 "
        "r)/\\omega^2$ is negative, and the cones have tipped over past the horizontal: the null curve moving "
        "to $+\\phi$, $dt = -\\tfrac{1}{2}(3 - \\sqrt{2})\\,d\\phi$, goes down in $t$, while the one moving to "
        "$-\\phi$, $dt = -\\tfrac{1}{2}(17 - \\sqrt{2})\\,d\\phi$, climbs steeply. Every horizontal line, run "
        "toward $+\\phi$, points into the future cones, so the circle of constant $t$, $r$, and $z$ is a closed "
        "timelike curve.",
        "The curve moving to $+\\phi$ comes round to its own $\\phi$ at a $t$ earlier by "
        "$(3 - \\sqrt{2})\\pi$ after each turn. None of the curves drawn here is a null geodesic: light "
        "launched along one moving to $+\\phi$ is turned toward the axis by $\\Gamma^r{}_{t\\phi}$ and "
        "$\\Gamma^r{}_{\\phi\\phi}$, and light launched along one moving to $-\\phi$ is turned away from it.",
    ],
    ("som_raychaudhuri", "cylindrical", "inside"): [
        "The cylinder of $t$ and $\\phi$ ($r = r_c/2$, $z = 0$, $r_c = c/\\Omega$) about the axis $r = 0$, the "
        "world line of one particle of the dust. It is opened along the line $\\phi = \\pm\\pi$ and drawn "
        "with $r\\phi$ across, so that its left and right edges are that one line. The metric on it is the "
        "same at every point, so its null curves are straight: $c\\,dt = r(1 - r/r_c)\\,d\\phi$ moving to "
        "$+\\phi$, which is $c\\,dt = \\tfrac{1}{4}r_c\\,d\\phi$ here, and $c\\,dt = -r(1 + r/r_c)\\,d\\phi$ "
        "moving to $-\\phi$, which is $c\\,dt = -\\tfrac{3}{4}r_c\\,d\\phi$. The cross term tilts every cone "
        "toward $+\\phi$, and a curve moving that way covers three times the $\\phi$ in a given $t$ that one "
        "moving the other way does.",
        "At this radius the curve moving to $+\\phi$ is a null geodesic: $\\Gamma^r{}_{t\\phi}$ and "
        "$\\Gamma^r{}_{\\phi\\phi}$ cancel along it, and light sent that way circles the axis at "
        "$r = r_c/2$, a Larmor orbit of radius $c/2\\Omega$. The curve moving to $-\\phi$ is not a geodesic, "
        "and light launched along it is turned away from the axis. The horizontal lines, circles of constant "
        "$t$, lie outside every cone and are spacelike.",
    ],
    ("som_raychaudhuri", "cylindrical", "beyond"): [
        "The cylinder of $t$ and $\\phi$ ($r = 3r_c/2$, $z = 0$, $r_c = c/\\Omega$), opened along "
        "$\\phi = \\pm\\pi$ in the same way. Beyond $r_c$ the coefficient $g_{\\phi\\phi} = r^2(1 - r^2/r_c^2)$ "
        "is negative, and the cones have tipped over past the horizontal: the null curve moving to $+\\phi$, "
        "$c\\,dt = -\\tfrac{3}{4}r_c\\,d\\phi$, goes down in $t$, while the one moving to $-\\phi$, "
        "$c\\,dt = -\\tfrac{15}{4}r_c\\,d\\phi$, climbs steeply. Every horizontal line, run toward $+\\phi$, "
        "points into the future cones, so the circle of constant $t$, $r$, and $z$ is a closed timelike curve.",
        "The curve moving to $+\\phi$ comes round to its own $\\phi$ at a $t$ earlier by $3\\pi r_c/2c$ after "
        "each turn. None of the curves drawn here is a null geodesic: light launched along one moving to "
        "$+\\phi$ is turned toward the axis by $\\Gamma^r{}_{t\\phi}$ and $\\Gamma^r{}_{\\phi\\phi}$, and light "
        "launched along one moving to $-\\phi$ is turned away from it.",
    ],
    ("som_raychaudhuri", "cartesian", "tx"): [
        "The plane of $t$ and $x$ ($y = z = 0$), in units of $r_c = c/\\Omega$. The surfaces $t = $ const are "
        "spacelike only within $r_c$ of the axis, where $g^{tt} = (x^2 + y^2)/r_c^2 - 1$ is negative, so the "
        "cones are oriented by $\\partial_t$, the direction of the dust's world lines, which is timelike everywhere.",
        "The metric on this plane is $-c^2dt^2 + dx^2$, so the curves drawn are null and run at 45°. Light "
        "launched along one is turned out of the plane into $y$, since $\\Gamma^y{}_{tx} = -\\Omega/c$ is not "
        "zero, and its path projects onto the plane of $x$ and $y$ as a circle of radius $r_c/2$. The closed "
        "timelike curves circle each world line of the dust beyond $r_c$, through $y$ as well as $x$, so they "
        "cross this plane.",
    ],
    ("tov", "spherical", "radial"): [
        "The plane of $t$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$) through a star of fluid with a polytrope "
        "for its equation of state, the same at every other angle by spherical symmetry. Its mass and redshift "
        "functions come from $G^t{}_t$ and $G^r{}_r$ and its pressure "
        "from $\\partial_r p = -(\\rho c^2 + p)\\,\\partial_r\\Phi$, and beyond its surface, where the "
        "pressure falls to zero, the same coordinates carry on as Schwarzschild's exterior. The "
        "rays obey $c\\,dt = \\pm e^{-\\Phi}(1 - 2m/r)^{-1/2}\\,dr$.",
        "The cones are narrowest at the centre, where a clock runs at $e^\\Phi = 0.67$ of the rate "
        "$t$ counts, and they open toward the surface, where it runs at $0.84$, and on toward 45° far "
        "away. They stay open everywhere, since at $2GM/c^2R = 0.29$ the star is well inside "
        "Buchdahl's bound of $8/9$ and has no horizon.",
    ],
    ("tov", "spherical", "through"): [
        "The line through the centre of the star in the plane $\\theta = \\pi/2$: $x = r$ on the "
        "right is $\\phi = 0$ and $x = -r$ on the left is $\\phi = \\pi$, and spherical symmetry makes the "
        "two halves mirror images. Rays cross the centre smoothly, where the cones are narrowest and "
        "the Kretschmann scalar is finite, and the surface crosses the line on both sides.",
    ],
    ("einstein_cluster", "areal", "radial"): [
        "The plane of $t$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$) through a cluster with no surface "
        "($b = 2\\,r_s$), the same at every other angle by spherical symmetry. The rays obey "
        "$c\\,dt = \\pm e^{-\\Phi}(1 - 2m/r)^{-1/2}\\,dr$.",
        "The cones are narrowest at the centre, where a clock runs at $e^\\Phi = 0.71$ of the rate $t$ counts, "
        "and they open toward 45° far away, where $m$ approaches $r_s/2$. They stay open everywhere, since "
        "$3m < r$ throughout the cluster.",
    ],
    ("einstein_cluster", "areal", "through"): [
        "The line through the centre of the cluster in the plane $\\theta = \\pi/2$: $x = r$ on the right is "
        "$\\phi = 0$ and $x = -r$ on the left is $\\phi = \\pi$, and spherical symmetry makes the two halves "
        "mirror images. Rays cross the centre smoothly, where the cones are narrowest and the Kretschmann scalar "
        "is finite.",
    ],
    ("einstein_cluster", "constant_speed", "radial"): [
        "The plane of $t$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$) over the whole cluster ($r \\in (0, R]$, "
        "$V = 1/2$, $R = 3\\,r_s$). The rays obey $c\\,dt = \\pm(1 + 2V^2)(r/R)^{-V^2}\\,dr$, so "
        "$ct \\mp (1 + 2V^2)\\,r\\,(r/R)^{-V^2}/(1 - V^2)$ is constant along each.",
        "The cones narrow toward the centre and close on it, where the density grows as $1/r^2$ and a clock at "
        "rest stops, and a ray still reaches the centre in a finite time $t$. Beyond $R$ the spacetime is "
        "Schwarzschild's exterior.",
    ],
    ("einstein_cluster", "isotropic", "radial"): [
        "The plane of $t$ and $\\rho$ ($\\theta = \\pi/2$, $\\phi = 0$) over the whole cluster in Einstein's "
        "coordinates ($\\sigma = 5 - 2\\sqrt{6}$, the cluster with $V = 1/2$). The rays obey "
        "$c\\,dt = \\pm\\sqrt{a/b}\\,d\\rho$.",
        "The surface is $\\rho_0 = 2.47\\,r_s$, where the areal radius is $3\\,r_s$, and the cones close on the "
        "centre as they do in the areal radius.",
    ],
    ("einstein_cluster", "uniform", "radial"): [
        "The plane of $t$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$) over the whole cluster ($r \\in [0, R]$, "
        "$R = 3\\,r_s$). The rays obey $c\\,dt = \\pm(1 - r_s/R)^{-3/4}(1 - r_sr^2/R^3)^{-1/4}\\,dr$.",
        "A clock at the centre runs at $0.74$ of the rate $t$ counts and a clock at the surface at $0.82$, and "
        "the cones stay open. Beyond $R$ the spacetime is Schwarzschild's exterior.",
    ],
    ("einstein_cluster", "uniform", "through"): [
        "The line through the centre of the cluster in the plane $\\theta = \\pi/2$: $x = r$ on the right is "
        "$\\phi = 0$ and $x = -r$ on the left is $\\phi = \\pi$. Rays cross the centre smoothly, where the "
        "Kretschmann scalar is finite.",
    ],
    ("einstein_cluster", "hyperspherical", "radial"): [
        "The plane of $t$ and $\\chi$ ($\\theta = \\pi/2$, $\\phi = 0$) over the whole cluster "
        "($\\chi \\in [0, \\chi_0]$, $\\sin\\chi_0 = 1/\\sqrt{3}$, the cluster with $R = 3\\,r_s$ and "
        "$a = \\sqrt{27}\\,r_s$). The rays obey $c\\,dt = \\pm a\\sqrt{\\cos\\chi/\\cos^3\\chi_0}\\,d\\chi$.",
        "Past the surface the same metric runs on to the equator of the three sphere, $\\chi = \\pi/2$, where "
        "the cones open without limit and the particles of a cluster would move at the speed of light.",
    ],
    ("einstein_rosen_waves", "cylindrical", "radial"): [
        "The plane of $t$ and $\\rho$ ($\\phi = 0$, $z = 0$) through a pulse of Weber, Wheeler, and Bonnor. "
        "The metric on it is $e^{2(\\gamma - \\psi)}(-c^2dt^2 + d\\rho^2)$, and $e^{2(\\gamma - \\psi)}$ drops out "
        "of the null condition, so for every wave the rays are at 45°, each ingoing ray meeting an outgoing one "
        "on the axis $\\rho = 0$. No Christoffel symbol turns them out of the plane, so they are null geodesics.",
        "The pulse comes in from the past and goes out again with its crest just outside the rays "
        "$\\rho = |ct|$, and is greatest on the axis at $t = 0$, where $\\psi = 2C/a$.",
    ],
    ("einstein_rosen_waves", "null", "radial"): [
        "The same plane in the null chart, $u$ and $v$ ($\\phi = 0$, $z = 0$), drawn against $(v - u)/2 = \\rho$ "
        "and $(u + v)/2 = ct$. Only $g_{uv}$ is nonzero on it, so a null direction has $du\\,dv = 0$, and the "
        "light rays are the coordinate lines $u = $ const and $v = $ const themselves, for every wave.",
        "The pulse comes in with its crest just outside the ray $v = 0$, is greatest on the axis at $u = v = 0$, "
        "and goes out with its crest just outside $u = 0$.",
    ],
    ("gowdy", "areal", "plane"): [
        "The plane of $t$ and $\\theta$ ($\\sigma = 0$, $\\delta = 0$) of a polarised wave on the torus, the "
        "edges $\\theta = 0$ and $\\theta = 2\\pi$ one line. The metric on it is "
        "$L^2t^{-1/2}e^{\\lambda/2}(-dt^2 + d\\theta^2)$, and the factor drops out of the null condition, so for "
        "every wave the rays are at 45°. No Christoffel symbol turns them out of the plane, so they are null "
        "geodesics.",
        "A ray leaving the singularity $t = 0$ has covered the angle $t$ by the time $t$, so two events at "
        "the time $t$ whose $\\theta$ differ by more than $2t$ have no event in both their pasts, and a ray "
        "first comes back to its own $\\theta$ at $t = 2\\pi$.",
    ],
    ("gowdy", "logarithmic", "plane"): [
        "The plane of $\\tau$ and $\\theta$ ($\\sigma = 0$, $\\delta = 0$) of the same wave, drawn against "
        "$-\\tau$ so that the future is up and the singularity, $\\tau \\to \\infty$, lies below the drawing. "
        "The rays are $\\theta \\pm e^{-\\tau} = $ const, so the cones close up toward the singularity as "
        "$d\\theta/d\\tau = \\pm e^{-\\tau}$.",
        "A ray covers the angle $e^{-\\tau}$ between the singularity and $\\tau$, so the pasts of two events "
        "at $\\tau$ whose $\\theta$ differ by more than $2e^{-\\tau}$ do not meet, and each $\\theta$ reaches "
        "the singularity with its own asymptotic velocity, $v(\\theta) = \\tfrac{1}{2}\\cos\\theta$.",
    ],
    ("gowdy", "sphere", "plane"): [
        "The plane of $t$ and $\\theta$ ($\\sigma = 0$, $\\delta = 0$) of the inside of Schwarzschild's horizon "
        "as a Gowdy universe on $S^2 \\times S^1$, the edges $\\theta = 0$ and $\\theta = \\pi$ the poles of the "
        "sphere. The metric on it is $L^2e^{2a}(-dt^2 + d\\theta^2)$, so for every wave the rays are at 45°. "
        "They are null geodesics, light running along a meridian of the sphere at one value of "
        "Schwarzschild's time.",
        "The universe begins at the singularity $t = 0$, where $r = 0$, and ends at $t = \\pi$, the horizon "
        "$r = r_s$, which makes it the white hole's side of the horizon. A ray takes the whole life of the "
        "universe to run from one pole to the other. On the diagonals $t = \\theta$ and $t + \\theta = \\pi$ "
        "the gradient of the orbit area, $4\\pi^2L^2\\sin t\\sin\\theta$, is null.",
    ],
    ("curzon_chazy", "weyl", "axis"): [
        "The plane of $t$ and $z$ ($\\rho = 0$, $\\phi = 0$) of the Curzon-Chazy particle ($m = 1$). The metric on it is "
        "$-e^{-2m/z}c^2dt^2 + e^{2m/z}dz^2$, so the rays are $ct = \\pm z_* + $ const with "
        "$z_* = z\\,e^{2m/z} - 2m\\,\\mathrm{Ei}(2m/z)$, and no Christoffel symbol turns them out of the plane, so they "
        "are null geodesics.",
        "The cones close as $e^{-2m/z}$ toward $z = 0$, which a ray reaches only as $t \\to \\pm\\infty$, though after "
        "a finite affine distance, since $z$ is an affine parameter along it. The Kretschmann scalar on the axis, "
        "$48m^2(z - m)^2e^{-4m/z}/z^8$, vanishes at $z = m$ and goes to zero at $z = 0$.",
    ],
    ("curzon_chazy", "weyl", "equator"): [
        "The plane of $t$ and $\\rho$ ($\\phi = 0$, $z = 0$) of the Curzon-Chazy particle ($m = 1$). The metric on it is "
        "$-e^{-2m/\\rho}c^2dt^2 + e^{2m/\\rho - m^2/\\rho^2}d\\rho^2$, so a ray has $c\\,dt/d\\rho = \\pm e^{2m/\\rho - m^2/2\\rho^2}$, and no "
        "Christoffel symbol turns it out of the plane, so the rays are null geodesics.",
        "The cones are narrowest at $\\rho = m/2$ and open without bound toward $\\rho = 0$, the ring, which every ingoing ray "
        "reaches in a finite time $t$ and where the Kretschmann scalar diverges as $e^{2m^2/\\rho^2}$.",
    ],
    ("curzon_chazy", "spherical", "axis"): [
        "The plane of $t$ and $r$ ($\\theta = 0$, $\\phi = 0$) of the Curzon-Chazy particle ($m = 1$). The metric on it is "
        "$-e^{-2m/r}c^2dt^2 + e^{2m/r}dr^2$, so the rays are $ct = \\pm r_* + $ const with "
        "$r_* = r\\,e^{2m/r} - 2m\\,\\mathrm{Ei}(2m/r)$, and no Christoffel symbol turns them out of the plane, so they "
        "are null geodesics.",
        "The cones close as $e^{-2m/r}$ toward $r = 0$, which a ray reaches only as $t \\to \\pm\\infty$, though after "
        "a finite affine distance, since $r$ is an affine parameter along it. The Kretschmann scalar on the axis, "
        "$48m^2(r - m)^2e^{-4m/r}/r^8$, vanishes at $r = m$ and goes to zero at $r = 0$.",
    ],
    ("curzon_chazy", "spherical", "equator"): [
        "The plane of $t$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$) of the Curzon-Chazy particle ($m = 1$). The metric on it is "
        "$-e^{-2m/r}c^2dt^2 + e^{2m/r - m^2/r^2}dr^2$, so a ray has $c\\,dt/dr = \\pm e^{2m/r - m^2/2r^2}$, and no "
        "Christoffel symbol turns it out of the plane, so the rays are null geodesics.",
        "The cones are narrowest at $r = m/2$ and open without bound toward $r = 0$, the ring, which every ingoing ray "
        "reaches in a finite time $t$ and where the Kretschmann scalar diverges as $e^{2m^2/r^2}$.",
    ],
    **{("zipoy_voorhees", system, view): text for system in ("spherical", "prolate_spheroidal")
       for view, text in _zipoy_voorhees_captions(system).items()},
    ("senovilla", "cylindrical", "radial"): [
        "The plane of $t$ and $\\rho$ ($\\phi = 0$, $z = 0$) of Senovilla's universe ($a = 1$). The metric on it is "
        "$\\cosh^4(act)\\cosh^2(3a\\rho)(-c^2dt^2 + d\\rho^2)$, and the factor drops out of the null condition, so the "
        "rays are at 45°, each ingoing ray meeting an outgoing one on the axis $\\rho = 0$. No Christoffel symbol "
        "turns them out of the plane, so they are null geodesics, and the Kretschmann scalar is finite everywhere "
        "on it.",
        "The fluid's world lines are the vertical lines. It contracts before $t = 0$ and expands after it, and its "
        "density, $15a^2c^4/(8\\pi G\\cosh^4(act)\\cosh^4(3a\\rho))$, is greatest at the marked event and falls away "
        "from it in every direction.",
    ],
    ("melvin", "cylindrical", "radial"): [
        "The plane of $t$ and $\\rho$ ($\\phi = 0$, $z = 0$) of Melvin's universe ($B = 1$). The metric on it is "
        "$(1 + B^2\\rho^2/4)^2(-c^2dt^2 + d\\rho^2)$, and the factor drops out of the null condition, so the rays are "
        "at 45°, each ingoing ray meeting an outgoing one on the axis $\\rho = 0$. No Christoffel symbol turns them "
        "out of the plane, so they are null geodesics, and the Kretschmann scalar is finite everywhere on it.",
    ],
    ("melvin", "ernst", "radial"): [
        "The plane of $t$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$) of Ernst's black hole ($B = 1/(2r_s)$). The metric on "
        "it is $(1 + B^2r^2/4)^2$ times Schwarzschild's, $-(1 - r_s/r)c^2dt^2 + dr^2/(1 - r_s/r)$, and the factor drops "
        "out of the null condition, so the rays are Schwarzschild's curves $ct = \\pm(r + r_s\\ln|r/r_s - 1|) + $ const, "
        "and they are the same curves in every plane of constant $\\theta$.",
        "On the equator and on the axis no Christoffel symbol turns them out of the plane, so there they are null "
        "geodesics. The cones close at the horizon $r = r_s$, and inside it every cone points to $r = 0$, where the "
        "Kretschmann scalar diverges.",
    ],
    ("point_particle_2plus1", "conical", "radial"): [
        "The plane of $t$ and $r$ ($\\phi = 0$) around a point particle ($\\alpha = 3/4$), each point in the plane "
        "a circle of circumference $2\\pi\\alpha r$ about the particle. The metric on the plane is "
        "$-c^2dt^2 + dr^2$, so every ray runs at 45°, a null geodesic of the spacetime, and $\\alpha$ enters only "
        "$g_{\\phi\\phi}$. The left edge, $r = 0$, is the particle's world line.",
    ],
    ("point_particle_2plus1", "wedge", "radial"): [
        "The plane of $t$ and $r$ ($\\theta = 0$), one face of the missing wedge ($\\alpha = 3/4$), with the "
        "particle's world line on its left edge. The line element is Minkowski's, so every ray runs at 45°, a null "
        "geodesic of the spacetime, and the particle shows only in the range of $\\theta$, which closes at "
        "$2\\pi\\alpha$.",
    ],
    ("point_particle_2plus1", "circumference", "radial"): [
        "The plane of $t$ and $R$ ($\\phi = 0$) around a point particle ($\\alpha = 3/4$), each point in the plane "
        "a circle of circumference $2\\pi R$ about the particle. The metric on the plane is "
        "$-c^2dt^2 + dR^2/\\alpha^2$, so the rays are the straight lines $ct = \\pm R/\\alpha + $ const, null "
        "geodesics of the spacetime, steeper than 45° because the circle at $R$ lies $R/\\alpha$ from the "
        "particle.",
    ],
    ("point_particle_2plus1", "planet", "radial"): [
        "The plane of $t$ and $\\chi$ ($\\phi = 0$) inside Gott and Alpert's planet, from its centre to its edge "
        "at $\\chi_0$ ($\\cos\\chi_0 = 3/4$). The metric on it is $-c^2dt^2 + a^2d\\chi^2$, so every ray runs at "
        "45° against $a\\chi$, a null geodesic of the spacetime, and the dust of the planet stays where it is: "
        "$\\Gamma^\\chi{}_{tt} = 0$, so nothing at rest falls toward the centre.",
    ],
    ("point_particle_2plus1", "isotropic", "radial"): [
        "The plane of $t$ and $\\rho$ ($\\phi = 0$) around a point particle ($\\alpha = 3/4$). The metric on it is "
        "$-c^2dt^2 + (\\rho/\\ell)^{2\\alpha - 2}d\\rho^2$, so the rays are the curves "
        "$ct = \\pm(\\ell/\\alpha)(\\rho/\\ell)^\\alpha + $ const, null geodesics of the spacetime, with cones "
        "that close toward the particle and open wider than 45° beyond $\\rho = \\ell$.",
        "The particle, $\\rho = 0$, lies at a finite proper distance, $(\\ell/\\alpha)(\\rho/\\ell)^\\alpha$ from "
        "the circle at $\\rho$, and a ray reaches it in a finite time.",
    ],
    ("point_particle_2plus1", "two_bodies", "between"): [
        "The plane of $t$ and $x$ ($y = 0$) between two particles at rest at $x = \\pm d$ "
        "($\\alpha_1 = \\alpha_2 = 3/4$, $\\ell = d$), one on each edge. The metric on it is "
        "$-c^2dt^2 + \\Omega\\,dx^2$ with $\\Omega = |x^2/d^2 - 1|^{-1/2}$, so the cones close toward each "
        "particle, which a ray reaches in a finite time: the two are $2.40\\,d$ apart in proper distance. The "
        "plane is one of symmetry, so no Christoffel symbol turns the rays out of it and they are null geodesics.",
    ],
    ("point_particle_2plus1", "two_bodies", "beyond"): [
        "The plane of $t$ and $x$ ($y = 0$) beyond the particle at $x = d$, which is the left edge, the other "
        "particle behind it at $x = -d$ ($\\alpha_1 = \\alpha_2 = 3/4$, $\\ell = d$). The metric on it is "
        "$-c^2dt^2 + \\Omega\\,dx^2$ with $\\Omega = |x^2/d^2 - 1|^{-1/2}$, and the rays are null geodesics, with "
        "cones that close toward the particle and open past 45° beyond $x = \\sqrt{2}\\,d$.",
        "Far from both particles $\\Omega$ falls as $d/|x|$, the conformal factor of one particle with "
        "$\\alpha = \\alpha_1 + \\alpha_2 - 1 = 1/2$, so the two masses add.",
    ],
    ("levi_civita", "weyl", "radial"): [
        "The plane of $t$ and $\\rho$ ($\\phi = 0$, $z = 0$) of Levi-Civita's cylinder ($\\sigma = 1/4$). The metric "
        "on it is $\\rho^{4\\sigma}(-c^2dt^2 + \\rho^{8\\sigma(\\sigma - 1)}d\\rho^2)$, so the rays are the curves "
        "$ct = \\pm 4\\rho^{1/4} + $ const, with cones that close toward the axis and open wider than 45° beyond "
        "$\\rho = 1$. No Christoffel symbol turns them out of the plane, so they are null geodesics.",
        "The Kretschmann scalar $3/(4\\rho^3)$ diverges on the axis $\\rho = 0$, and a ray that leaves the axis "
        "reaches every $\\rho$ in a finite time, so the singularity is naked.",
    ],
    ("levi_civita", "kasner", "radial"): [
        "The plane of $t$ and $r$ ($\\phi = 0$, $z = 0$) of Levi-Civita's cylinder in its Kasner form "
        "($p_0 = p_2 = 2/3$, $p_3 = -1/3$), with $r$ the proper distance from the axis. The metric on it is "
        "$-r^{2p_0}c^2dt^2 + dr^2$, so the rays are the curves $ct = \\pm 3r^{1/3} + $ const, null geodesics of the "
        "spacetime, and the Kretschmann scalar $-16p_0p_2p_3/r^4$ diverges on the axis $r = 0$.",
    ],
    ("malament_hogarth", "cartesian", "tx"): [
        "The plane of $t$ and $x$ ($y = z = 0$) through the removed event at the origin. "
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
        "The line through the centre of the collapsing star in the plane $\\theta = \\pi/2$, in "
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
        "The plane of $t$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$) outside the collapsing "
        "star, where the metric is Schwarzschild's. The surface falls freely from rest at $r = 2r_s$ at "
        "$t = 0$, along its radial geodesic, and inside it, $r < R(t)$, lies the star, which these "
        "coordinates do not cover.",
        "The surface reaches $r_s$ only as $t \\to \\infty$, though its own clock reads a finite time "
        "there, and each outgoing ray it sends takes longer than the last to climb away. Outside the "
        "star the horizon and the black hole behind it lie beyond these coordinates; inside it the "
        "comoving coordinates carry on across both to the crunch.",
    ],
    ("szekeres", "axisymmetric", "north"): [
        "The plane of $t$ and the comoving $r$ along the axis of symmetry ($\\theta = 0$) through a cloud of dust whose shells "
        "fall as those of the Tolman-Bondi cloud do, with vacuum outside its surface $r_b$. Every shell reaches $R = 0$ at "
        "its own time, the centre first, at $ct = 0.60\\,r_b$, and the surface at $0.94\\,r_b$, and the singularity, "
        "where the Kretschmann scalar diverges, is that curve. The rays obey "
        "$c\\,dt = \\pm\\left(\\partial_r R + R\\,S'/S\\right)dr$, where $R\\,S'/S$ is the rate at which the centres of "
        "the spheres move along the axis toward $\\theta = 0$, so on this half of the axis neighbouring shells stand "
        "farther apart than on the other and the cones are narrower in $r$.",
        "The marked ray is the last along this half of the axis to reach infinity. It leaves the centre at "
        "$ct = -0.81\\,r_b$ and crosses the surface at $0.61\\,r_b$, and outside the cloud, where $S$ is constant, it runs "
        "along the horizon of Schwarzschild's exterior in Georges Lemaître's coordinates.",
    ],
    ("szekeres", "axisymmetric", "south"): [
        "The plane of $t$ and the comoving $r$ along the axis of symmetry ($\\theta = \\pi$) through the same cloud, "
        "the other half of the axis. The rays obey $c\\,dt = \\pm\\left(\\partial_r R - R\\,S'/S\\right)dr$, so "
        "on this half neighbouring shells stand closer together and the cones are wider in $r$, while every shell "
        "reaches $R = 0$ at the same time as on the other half, since $R$ depends on $t$ and $r$ alone.",
        "The last ray along this half of the axis to reach infinity leaves the centre at $ct = -0.12\\,r_b$, later than "
        "its counterpart toward $\\theta = 0$ by $0.69\\,r_b$, and crosses the surface at the same $0.61\\,r_b$.",
    ],
    ("tolman_bondi", "comoving_synchronous", "collapse"): [
        "The plane of $t$ and the comoving $r$ ($\\theta = \\pi/2$, $\\phi = 0$) through a cloud of dust whose density falls from its centre to "
        "zero at its surface $r_b$, with vacuum outside. Every shell falls on its own clock, $R^{3/2} = r^{3/2} - \\tfrac{3}{2}\\sqrt{2GM(r)/c^2}\\,ct$, and "
        "reaches $R = 0$ at its own time: the centre first, at $ct = 0.60\\,r_b$, and the surface at $0.94\\,r_b$. The singularity, where the Kretschmann "
        "scalar diverges, is that curve, and beyond it there is no spacetime. The rays obey $c\\,dt = \\pm\\partial_r R\\,dr$, and outside the cloud the same "
        "coordinates are Georges Lemaître's for Schwarzschild's exterior, carried by observers who fall freely from rest at infinity.",
        "The dotted curve, where $|\\nabla R|^2 = 1 - 2GM(r)/c^2R$ vanishes, bounds the trapped "
        "spheres. It meets the surface as the surface crosses $r_s = r_b/2$ and runs inward, dipping "
        "below the moment the centre is crushed, so trapped spheres form in the body of the cloud "
        "before its centre becomes singular. The event horizon leaves the centre at "
        "$ct = -0.42\\,r_b$, well before any of this, and outside the cloud it runs along the dotted "
        "curve, the sphere $R = 2GM/c^2$.",
    ],
    ("vaidya", "eddington_finkelstein_ingoing", "shell"): [
        "The plane of $v$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$), drawn with $cv - r$ as "
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
        "The plane of $u$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$), drawn with $cu + r$ as "
        "the vertical axis so that the outgoing rays, $u$ constant, run at 45°. It is the imploding shell run backward "
        "in time: a shell of null dust carries the whole mass $M$ out along $u = 0$, with "
        "Schwarzschild's metric in outgoing coordinates before it passes and flat space after.",
        "Before the shell, the region inside $r_s$ is a white hole: every future cone there points "
        "to larger $r$, and $r = 0$, where the Kretschmann scalar diverges, lies in its past. Its "
        "horizon is the ingoing ray that stays at $r_s$ until the shell leaves and then crosses flat "
        "space to the centre, arriving at $cu = 2r_s$. Every ingoing ray that reaches the centre "
        "before that moment came out of the white hole, and every one after came in from far away.",
    ],
    ("bonnor_vaidya", "eddington_finkelstein_ingoing", "shell"): [
        "The plane of $v$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$), drawn with $v - r$ as the vertical axis so "
        "that the ingoing rays, $v$ constant, run at 45°. A shell of charged null dust of mass $M$ and charge "
        "$q = 0.96\\,M$ falls in along $v = 0$: before it the metric is flat, and after it the metric is "
        "Reissner-Nordström's in ingoing coordinates, with $g^{rr} = 1 - 2M/r + q^2/r^2$ vanishing on "
        "$r_+ = 1.28\\,M$ and $r_- = 0.72\\,M$. The event horizon is the outgoing ray that reaches $r_+$ just as "
        "the shell does. It leaves the centre at $v = -2.56\\,M$ and crosses flat space.",
        "Between $r_+$ and $r_-$ every future cone points to smaller $r$, and inside $r_-$ the cones open "
        "again, so the outgoing rays there climb toward $r_-$ and $r = 0$, where the Kretschmann scalar "
        "diverges, is a timelike line. Across the shell $G_{vv}$ jumps by $2(Mr - q^2/2)/r^3$ times a delta "
        "function of $v$, so the energy density of the dust is positive outside the marked sphere "
        "$r = q^2/2M = 0.46\\,M$ and negative inside it. Amos Ori showed in 1991 that the Lorentz force turns "
        "the charged dust round on that sphere, and the stretch of the shell below the mark is drawn as the "
        "metric stands.",
    ],
    ("bonnor_vaidya", "eddington_finkelstein_outgoing", "shell"): [
        "The plane of $u$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$), drawn with $u + r$ as the vertical axis so "
        "that the outgoing rays, $u$ constant, run at 45°. It is the falling shell run backward in time, in "
        "Bonnor and Vaidya's own coordinates: a shell of charged null dust carries the whole mass $M$ and the "
        "whole charge $q = 0.96\\,M$ out along $u = 0$, with Reissner-Nordström's metric in outgoing "
        "coordinates before it passes and flat space after.",
        "Before the shell, the region between $r_- = 0.72\\,M$ and $r_+ = 1.28\\,M$ is a white hole, where every "
        "future cone points to larger $r$. Its horizon is the ingoing ray that stays at $r_+$ until the shell "
        "leaves and then crosses flat space to the centre, arriving at $u = 2.56\\,M$. The energy density of "
        "the shell is positive outside the marked sphere $r = q^2/2M = 0.46\\,M$ and negative inside it.",
    ],
    ("bonnor_vaidya", "homothetic", "scaling"): [
        "The plane of $V$ and $R$ ($\\theta = \\pi/2$, $\\phi = 0$) for a mass and a charge that grow in proportion "
        "to the advanced time, $m = \\mu v$ and $q = Qm/M$. A line of constant $R$ is a sphere whose areal radius "
        "grows with the mass, $r = Rm/M$, and a step up in $V$ multiplies every length by one factor, so the cones "
        "depend on $R$ alone. The ingoing rays are level, $V$ constant.",
        "The component $g^{RR}$ vanishes on $R = 0.67\\,M$, $1.83\\,M$, and $7.50\\,M$, and on each of those lines an "
        "outgoing ray keeps its $R$: all three come out of the point $v = 0$, $r = 0$, where the collapse begins. "
        "The marginally trapped spheres, $r^2 - 2mr + q^2 = 0$, are $R = 0.72\\,M$ and $1.28\\,M$, and the "
        "Kretschmann scalar diverges on $R = 0$. The energy density of the dust is positive outside the marked "
        "sphere $R = Q^2/M = 0.92\\,M$ and negative inside it.",
    ],
    ("robinson_trautman", "axisymmetric", "axis"): [
        "The plane of $u$ and $r$ on the axis of symmetry ($\\theta = 0$), drawn with $cu + r$ as "
        "the vertical axis so that the outgoing rays, $u$ constant, run at 45°. Each outgoing ray is one "
        "point of every wave front it crosses, and $r$ is the affine distance along it. The first front "
        "is $u = 0$, and the Robinson-Trautman equation carries it toward the future only. The fronts are "
        "even about the equator, so $\\partial_\\theta H$ vanishes on the axis and the null curves drawn "
        "are null geodesics.",
        "The ingoing rays obey $dr/d(cu) = -H$, with $2H = K - 2r\\,\\partial_u\\ln f - 2m/r$. On the "
        "axis the first front has $K = 2.25$ and $\\partial_u\\ln f = 0.80/m$, so $H < 0$ at every $r$, "
        "and until $cu = 1.18\\,m$ an ingoing ray far out on the axis gains $r$: the ray that leaves the first "
        "front at $r = 3.5\\,m$ reaches $4.07\\,m$ before it turns. The curve $g^{rr} = 2H = 0$ comes in from "
        "large $r$ at $cu = 1.18\\,m$ and closes on $r = 2m$, where it stays once the fronts are round, and "
        "from then on the plane is Schwarzschild's in outgoing coordinates. Inside it every future cone "
        "points to larger $r$, a white hole, and $r = 0$, where the Kretschmann scalar $48m^2/r^6$ "
        "diverges, lies in its past.",
    ],
    ("robinson_trautman", "axisymmetric", "equator"): [
        "The plane of $u$ and $r$ on the equator ($\\theta = \\pi/2$, $\\phi = 0$), drawn with $cu + r$ as "
        "the vertical axis so that the outgoing rays, $u$ constant, run at 45°. The fronts are even about "
        "the equator, so $\\partial_\\theta H$ vanishes here and the null curves drawn are null geodesics. "
        "On the equator the first front has $K = 0.49$ and $\\partial_u\\ln f = -0.14/m$: $f$ falls here while "
        "it rises at the poles.",
        "The curve $g^{rr} = 2H = 0$ meets the first front at $r = 1.85\\,m$, falls to $1.70\\,m$ at "
        "$cu = 0.2\\,m$, and rises to $r = 2m$ as the fronts grow round. Outside it the ingoing rays lose $r$ "
        "from the first front on, and inside it they gain $r$ and every future cone points to larger $r$, a "
        "white hole with $r = 0$ in its past. By $cu = 3m$ the fronts are round to within $0.2\\%$, and "
        "this plane and the axis's are both Schwarzschild's in outgoing coordinates.",
    ],
    ("photon_rocket", "rectilinear", "behind"): [
        "The plane of $u$ and $r$ on the axis behind the rocket ($\\theta = 0$), drawn with $cu + r$ as the vertical axis so that the outgoing rays, $u$ constant, run at 45°. Each outgoing ray left the rocket at the retarded time $u$, and $r$ is the affine distance along it. On the axis $\\sin\\theta = 0$ removes $g_{u\\theta}$ and every $\\Gamma^\\theta$ of the plane, so the null curves drawn are null geodesics.",
        "The ingoing rays obey $dr/d(cu) = -g^{rr}/2$, with $g^{rr} = 1 - 2m/r - 2\\alpha r$. Before the burn the plane is Schwarzschild's in outgoing coordinates, with $g^{rr} = 0$ on $r = 2m_0$: inside it every future cone points to larger $r$, a white hole, and $r = 0$, where the Kretschmann scalar $48m^2/r^6$ diverges, lies in its past. During the burn a second zero of $g^{rr}$ comes in from large $r$, reaches $r = 4.73\\,m_0$ at $cu = 4.59\\,m_0$, and goes out again. Beyond it a ray sent after the rocket gains $r$, since the rocket accelerates away from it. The inner zero swells to $2.17\\,m_0$ and then shrinks with the mass to $0.60\\,m_0$, where it stays once the burn is over and the plane is Schwarzschild's again.",
    ],
    ("photon_rocket", "rectilinear", "ahead"): [
        "The plane of $u$ and $r$ on the axis ahead of the rocket ($\\theta = \\pi$), drawn with $cu + r$ as the vertical axis so that the outgoing rays, $u$ constant, run at 45°. Each outgoing ray left the rocket at the retarded time $u$, and $r$ is the affine distance along it. On the axis $\\sin\\theta = 0$ removes $g_{u\\theta}$ and every $\\Gamma^\\theta$ of the plane, so the null curves drawn are null geodesics.",
        "The ingoing rays obey $dr/d(cu) = -g^{rr}/2$, with $g^{rr} = 1 - 2m/r + 2\\alpha r$, which has one zero. Before the burn it is $r = 2m_0$, the horizon of a white hole in Schwarzschild's outgoing coordinates, with $r = 0$, where the Kretschmann scalar $48m^2/r^6$ diverges, in its past. Through the burn it lies inside $2m$, at $0.85$ of it where $\\alpha m$ is greatest, and it ends on $r = 0.60\\,m_0$, the horizon of the lighter mass. The rocket accelerates toward the light coming at it, so an ingoing ray loses $r$ faster during the burn than before or after it.",
    ],
    ("photon_rocket", "robinson_trautman", "ahead"): [
        "The plane of $u$ and $r$ on the axis ahead of the rocket ($\\theta = 0$), drawn with $cu + r$ as the vertical axis so that the outgoing rays, $u$ constant, run at 45°. Each outgoing ray left the rocket at the retarded time $u$, and $r$ is the affine distance along it. On the axis $\\partial_\\theta p = 0$ removes every $\\Gamma^\\theta$ of the plane, so the null curves drawn are null geodesics, and there $\\partial_u p/p = -\\alpha$, the rocket's acceleration.",
        "The ingoing rays obey $dr/d(cu) = -g^{rr}/2$, with $g^{rr} = 1 - 2m/r - 2r\\,\\partial_u p/p$, which has one zero. Before the burn it is $r = 2m_0$, the horizon of a white hole in Schwarzschild's outgoing coordinates, with $r = 0$, where the Kretschmann scalar $48m^2/r^6$ diverges, in its past. Through the burn it lies inside $2m$, at $0.85$ of it where $\\alpha m$ is greatest, and it ends on $r = 0.60\\,m_0$, the horizon of the lighter mass. The rocket accelerates toward the light coming at it, so an ingoing ray loses $r$ faster during the burn than before or after it.",
    ],
    ("photon_rocket", "robinson_trautman", "behind"): [
        "The plane of $u$ and $r$ on the axis behind the rocket ($\\theta = \\pi$), drawn with $cu + r$ as the vertical axis so that the outgoing rays, $u$ constant, run at 45°. Each outgoing ray left the rocket at the retarded time $u$, and $r$ is the affine distance along it. On the axis $\\partial_\\theta p = 0$ removes every $\\Gamma^\\theta$ of the plane, so the null curves drawn are null geodesics, and there $\\partial_u p/p = \\alpha$, the rocket's acceleration.",
        "The ingoing rays obey $dr/d(cu) = -g^{rr}/2$, with $g^{rr} = 1 - 2m/r - 2r\\,\\partial_u p/p$. Before the burn the plane is Schwarzschild's in outgoing coordinates, with $g^{rr} = 0$ on $r = 2m_0$: inside it every future cone points to larger $r$, a white hole, and $r = 0$, where the Kretschmann scalar $48m^2/r^6$ diverges, lies in its past. During the burn a second zero of $g^{rr}$ comes in from large $r$, reaches $r = 4.73\\,m_0$ at $cu = 4.59\\,m_0$, and goes out again. Beyond it a ray sent after the rocket gains $r$, since the rocket accelerates away from it. The inner zero swells to $2.17\\,m_0$ and then shrinks with the mass to $0.60\\,m_0$, where it stays once the burn is over and the plane is Schwarzschild's again.",
    ],
    ("c_metric", "spherical", "inner"): [
        "The plane of $t$ and $r$ on the half axis between the black holes ($\\theta = 0$, $\\phi = 0$), "
        "drawn for $\\alpha m = 1/6$. The curves drawn are null, and on the axis they are also null "
        "geodesics, the paths light takes, since $\\Gamma^\\theta{}_{tt}$ and $\\Gamma^\\theta{}_{rr}$ carry a "
        "factor $\\sin\\theta$. There $g^{rr} = (1 - 2m/r)(1 - \\alpha^2r^2)(1 + \\alpha r)^2$ vanishes at the "
        "black hole horizon $r = 2m$ and at the acceleration horizon $r = 1/\\alpha = 6m$, and the cones close "
        "at both.",
        "Inside $2m$ and beyond $6m$ the coordinate $r$ is the time. We take the future from the extensions of "
        "Jerry Griffiths, Pavel Krtouš, and Jiří Podolský across both horizons, which make the region inside "
        "$2m$ the black hole, where the cones point to $r = 0$, and the region beyond $6m$ the one to the "
        "future of the acceleration horizon, where the cones point to larger $r$ and null infinity lies ahead, "
        "at $r = -1/\\alpha$ past $r = \\infty$. The Kretschmann scalar $48m^2(1 + \\alpha r)^6/r^6$ diverges "
        "at $r = 0$.",
    ],
    ("c_metric", "spherical", "outer"): [
        "The plane of $t$ and $r$ on the half axis beyond the black hole ($\\theta = \\pi$, $\\phi = 0$), drawn "
        "for $\\alpha m = 1/6$, with the null curves null geodesics as on the inner axis. There the factor "
        "$1 + \\alpha r\\cos\\theta = 1 - \\alpha r$ vanishes at $r = 1/\\alpha = 6m$, so on this half of the "
        "axis the acceleration horizon lies at null infinity, where the plane ends. The cones close at "
        "$r = 2m$, and inside it they point to $r = 0$, in the black hole.",
        "With $C = 1/(1 + 2\\alpha m)$ this half of the axis carries the cosmic string, whose deficit angle "
        "$8\\pi\\alpha m/(1 + 2\\alpha m)$ lies in the angle about the axis and leaves this plane unchanged. The "
        "Kretschmann scalar $48m^2(1 - \\alpha r)^6/r^6$ diverges at $r = 0$ and vanishes at null infinity.",
    ],
    ("c_metric", "hong_teo", "inner"): [
        "The plane of $\\tau$ and $y$ on the half axis between the black holes ($x = 1$, $\\phi = 0$), drawn "
        "for $\\alpha m = 1/6$, with $y = 1/(\\alpha r)$ and $\\tau = \\alpha ct$. The coordinate $y$ runs on "
        "through $y = 0$, where $r$ is infinite, to null infinity at $y = -1$, where $x + y$ vanishes. The cones "
        "close at the acceleration horizon $y = 1$ and at the black hole horizon $y = 1/(2\\alpha m) = 3$.",
        "Between $y = -1$ and $y = 1$ the coordinate $y$ is the time, and we take the future as Griffiths, "
        "Krtouš, and Podolský do beyond the acceleration horizon, so that the cones point to smaller $y$ and "
        "every future directed ray reaches null infinity. Beyond $y = 3$ they point to larger $y$, into the "
        "black hole, whose singularity lies at $y = \\infty$.",
    ],
    **{("siklos", "siklos", view): [
        f"The plane of $u$ and $v$ (${x}$, $y = 0$), at one place on the wave front, drawn with $(v - u)/2$ "
        "across and $(u + v)/2$ up. On it the metric is $(L^2/x^2)(-2\\,du\\,dv + H\\,du^2)$, so away from the "
        "pulse the cones stand at 45°, as in anti-de Sitter space, and inside it they lean the way the wave "
        "travels. The rays moving right keep their $u$ and travel with the wave, and a curve moving left is "
        f"carried along $v$ by $\\tfrac{{1}}{{2}}\\int H\\,du = {shift}\\,L$ in crossing the pulse.",
        "The rays moving right are null geodesics, and away from the pulse so are the curves moving left. "
        "Inside it light crossing the wave is pushed toward larger $x$, "
        "$\\ddot{x} = \\tfrac{1}{2}\\partial_xH\\,\\dot{u}^2$, out of this plane, so there a curve moving left "
        "is a null curve of the plane.",
    ] for view, x, shift in (("near", "x = L", "0.44"), ("far", "x = 3L/2", "1.50"))},
    ("siklos", "ozsvath_robinson_rozga", "centre"): [
        "The plane of $u$ and $v$ at the centre of the disc ($\\xi = \\eta = 0$), which is $x = L$, $y = 0$ of "
        "Siklos's chart, drawn with $(v - u)/2$ across and $(u + v)/2$ up. There $p = q = 1$ and the metric "
        "of the plane is $-2\\,du\\,dv + h\\,du^2$, that of a pp-wave. The rays moving right keep their $u$ and "
        "travel with the wave, and a curve moving left is carried along $v$ by "
        "$\\tfrac{1}{2}\\int h\\,du = 0.44\\,L$ in crossing the pulse.",
        "The rays moving right are null geodesics, and away from the pulse so are the curves moving left. "
        "Inside it light crossing the wave is pushed toward negative $\\xi$, "
        "$\\ddot{\\xi} = -(3h/2L)\\,\\dot{u}^2$, out of this plane, so there a curve moving left is a null "
        "curve of the plane.",
    ],
    ("siklos", "kaigorodov", "depth"): [
        "The plane of $v$ and $x$ ($y = 0$) with $u$ divided out, $-(L^5/x^5)\\,dv^2 + (L^2/x^2)\\,dx^2$, the "
        "metric orthogonal to the lines of $u$. Its null curves are the shadows on $v$ and $x$ of the null "
        "geodesics with no momentum along $u$, each moving in $u$ at $du/dv = L^3/x^3$, and each cone is the "
        "future cone of those directions.",
        "The rays keep $v \\pm \\tfrac{2}{5}x^{5/2}/L^{3/2}$. They reach the conformal boundary $x = 0$ at a "
        "finite $v$, as light reaches anti-de Sitter's. The cones close toward larger $x$: the singularity "
        "$x \\to \\infty$ lies at $v = \\pm\\infty$, a finite affine distance away, and no curvature scalar "
        "diverges there, only the tidal forces on whatever falls in.",
    ],
    ("siklos", "kaigorodov_poincare", "depth"): [
        "The plane of $t$ and $z$ ($y = 0$) with $x$ divided out, "
        "$(L^2/z^2)\\left(-c^2dt^2/(1 + z^3/2L^3) + dz^2\\right)$, the metric orthogonal to the lines of $x$. "
        "Its null curves are the shadows on $t$ and $z$ of the null geodesics with no momentum along $x$, each "
        "moving in $x$ at $dx/d(ct) = z^3/(z^3 + 2L^3)$, and each cone is the future cone of those directions.",
        "Near the conformal boundary $z = 0$ the cones stand at 45°, as on the Poincaré patch of anti-de Sitter "
        "space, and a ray reaches the boundary at a finite $t$. Deeper in the wave drags the light along $x$ "
        "and the cones close, $c\\,dt = \\pm\\sqrt{1 + z^3/2L^3}\\,dz$, so the singularity $z \\to \\infty$ "
        "lies at $t = \\pm\\infty$.",
    ],
    ("siklos", "kaigorodov_horospheric", "depth"): [
        "The plane of $v$ and $\\rho$ ($y = 0$) with $u$ divided out, $-e^{5\\rho/L}dv^2 + d\\rho^2$, the metric "
        "orthogonal to the lines of $u$. Its null curves are the shadows on $v$ and $\\rho$ of the null "
        "geodesics with no momentum along $u$, each moving in $u$ at $du/dv = e^{3\\rho/L}$, and each cone is "
        "the future cone of those directions.",
        "The rays keep $v \\pm \\tfrac{2}{5}Le^{-5\\rho/2L}$. Toward the conformal boundary, $\\rho \\to \\infty$, "
        "the cones open without limit, and a ray arrives there at a finite $v$. Toward the singularity, "
        "$\\rho \\to -\\infty$, they close, and the singularity lies at $v = \\pm\\infty$.",
    ],
    ("siklos", "kaigorodov_stationary", "plane"): [
        "The plane of $u$ and $v$ ($\\rho = 0$, $y = 0$) in the region $x < 0$ of Siklos's chart, drawn with "
        "$(v - u)/2$ across and $(u + v)/2$ up. On it the metric is $-du^2 - 2\\,du\\,dv$, the same at every "
        "event, so every cone is the same cone, opened wide toward the left. A line of constant $v$ is "
        "timelike here, the world line of an observer at rest in this stationary region.",
        "The rays moving right keep their $u$ and are null geodesics. A curve moving left keeps $v + u/2$, and "
        "light sent along it is turned toward larger $\\rho$, $\\ddot{\\rho} = (3/2L)\\,\\dot{u}^2$, out of "
        "this plane, so it is a null curve of the plane.",
    ],
    ("siklos", "kaigorodov_homogeneous", "depth"): [
        "The plane of $U$ and $Z$ ($y = 0$) with $X$ divided out, drawn at $k = 1/10$. Its null curves are the "
        "shadows on $U$ and $Z$ of the null geodesics with no momentum along $X$, the rays of Siklos's plane "
        "of $v$ and $x$ with $x = Le^{2Z}$ and $v = -Ue^{5Z}$, and each cone is the future cone of those "
        "directions.",
        "The rays keep $(U \\pm \\tfrac{2}{5}L)e^{5Z}$, so the lines $U = \\pm 2L/5$ are two of them, and the "
        "future is the direction of decreasing $Ue^{5Z}$. The conformal boundary is $Z \\to -\\infty$ and the "
        "singularity $Z \\to \\infty$.",
    ],
    ("siklos", "kaigorodov_kundt", "depth"): [
        "The plane of $V$ and $x$ ($y = 0$) with $U$ divided out. Its null curves are the shadows on $V$ and "
        "$x$ of the null geodesics with no momentum along $U$, the rays of Siklos's plane of $v$ and $x$ with "
        "$v = Vx^2/\\sqrt{2}$ and $Lx$ for his $x$, and each cone is the future cone of those directions.",
        "The rays keep $Vx^2/\\sqrt{2} \\pm \\tfrac{2}{5}Lx^{5/2}$, and the future is the direction of "
        "increasing $Vx^2$. The conformal boundary is $x = 0$ and the singularity $x \\to \\infty$.",
    ],
}


# ---------------------------------------------------------------- reading a chart

def strip_lhs(latex):
    """'K = ...' -> '...'; a bare value is returned unchanged, and a condition such as the
    domain wall's "\\;\\text{for}\\; z \\neq 0" is left off, since every ray is traced off the
    zero it names."""
    head, eq, tail = latex.partition("=")
    value = tail if eq and len(head.strip()) <= 3 else latex
    return vm.scalar_parts(value)[0]


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

class FrontSolver:
    """Robinson and Trautman's wave fronts for declared axisymmetric initial data, G = c = m = 1.

    The fronts' metric is r^2 f^-2 (dtheta^2 + sin^2 theta dphi^2), and in vacuum
    d_u f = -f^3 Lap K / 12, with Lap the Laplacian of the unit sphere and
    K = f^2 (1 + Lap ln f) = f^2 + f Lap f - (1 - x^2)(d_x f)^2 the curvature of the front at
    r = 1, x = cos theta. f is a series of Legendre polynomials in x, as Macedo and Saa take
    it, and the equation is projected onto them by Gauss-Legendre quadrature, exact for the
    polynomials it meets, and integrated by an implicit Runge-Kutta method from the data
    f(0, theta)^2 = f_0^2 (1 - epsilon^2 cos^2 theta), f_0^2 = ln((1 + epsilon)/(1 - epsilon))/(2 epsilon),
    which gives every front the area 4 pi r^2. The equation has no solution toward the past,
    so there is no front before u = 0, and from `end` on the fronts are round to rounding.
    The vacuum 2H = K - 2r d_u ln f - 2/r follows. null_rays draws light with it, embedding.py
    the fronts and conformal.py the plane of its axis.
    """

    funcs = ()
    end = 16.0
    OFF_AXIS = 1e-3     # where a published value that divides by sin(theta) is read beside the axis

    def __init__(self, epsilon, theta=0.0, modes=40):
        from numpy.polynomial import legendre as L
        self.epsilon, self.theta, self.n = float(epsilon), float(theta), modes
        x, w = L.leggauss(4 * modes)
        self.x = x
        eye = np.eye(2 * modes + 1)
        self.V = L.legvander(x, 2 * modes)
        self.V1 = np.stack([L.legval(x, L.legder(row, 1)) for row in eye], axis=1)
        self.V2 = np.stack([L.legval(x, L.legder(row, 2)) for row in eye], axis=1)
        self.project = (self.V * w[:, None]).T * ((2 * np.arange(2 * modes + 1) + 1) / 2)[:, None]
        self.ll = np.arange(2 * modes + 1) * (np.arange(2 * modes + 1) + 1.0)
        self.f0sq = math.log((1 + self.epsilon) / (1 - self.epsilon)) / (2 * self.epsilon)
        b0 = (self.project @ np.sqrt(self.f0sq * (1 - self.epsilon ** 2 * x ** 2)))[:modes + 1]
        self.ivp = solve_ivp(lambda _, b: self.rate(b), (0.0, self.end), b0, method="Radau",
                             jac=lambda _, b: self.jacobian(b), rtol=1e-10, atol=1e-12, dense_output=True)
        if not self.ivp.success:
            raise SystemExit("FrontSolver: the Robinson-Trautman equation was not integrated")
        self._series = {}

    def on_nodes(self, b):
        """f, K and Lap K at the quadrature nodes, and the Legendre coefficients of K."""
        n = self.n + 1
        x = self.x
        f, fx, fxx = self.V[:, :n] @ b, self.V1[:, :n] @ b, self.V2[:, :n] @ b
        K = f ** 2 + f * ((1 - x ** 2) * fxx - 2 * x * fx) - (1 - x ** 2) * fx ** 2
        k = self.project @ K
        return f, K, -(self.V @ (self.ll * k)), k

    def rate(self, b):
        f, _, lap, _ = self.on_nodes(b)
        return (self.project @ (-f ** 3 * lap / 12))[:self.n + 1]

    def jacobian(self, b):
        """The derivative of rate(b) in b, written out, since differences of a rate that carries
        the fourth derivative of f lose to rounding the digits Newton's iteration needs."""
        n = self.n + 1
        x = self.x[:, None]
        Vn, V1n, V2n = self.V[:, :n], self.V1[:, :n], self.V2[:, :n]
        f, fx, fxx = (Vn @ b)[:, None], (V1n @ b)[:, None], (V2n @ b)[:, None]
        dK = (2 * f + (1 - x ** 2) * fxx - 2 * x * fx) * Vn + f * ((1 - x ** 2) * V2n - 2 * x * V1n) \
            - 2 * (1 - x ** 2) * fx * V1n
        lap = self.on_nodes(b)[2][:, None]
        dlap = -(self.V @ (self.ll[:, None] * (self.project @ dK)))
        return (self.project @ (-(3 * f ** 2 * lap * Vn + f ** 3 * dlap) / 12))[:n]

    def modes_at(self, u):
        return self.ivp.sol(min(max(float(u), 0.0), self.end))

    def series(self, u):
        """The Legendre coefficients of f, d_u f, d_u^2 f, K and d_u K at retarded time u, d_u b
        being the rate and d_u^2 b the Jacobian applied to it; kept, since a front is read at
        many angles."""
        u = min(max(float(u), 0.0), self.end)
        if u not in self._series:
            b = self.modes_at(u)
            late = u >= self.end
            db = 0 * b if late else self.rate(b)
            ddb = 0 * b if late else self.jacobian(b) @ db
            n = self.n + 1
            fn, f1, f2 = (self.V[:, :n] @ b)[:, None], (self.V1[:, :n] @ b)[:, None], (self.V2[:, :n] @ b)[:, None]
            xn = self.x[:, None]
            dK = ((2 * fn + (1 - xn ** 2) * f2 - 2 * xn * f1) * self.V[:, :n]
                  + fn * ((1 - xn ** 2) * self.V2[:, :n] - 2 * xn * self.V1[:, :n]) - 2 * (1 - xn ** 2) * f1 * self.V1[:, :n])
            if len(self._series) > 20000:
                self._series.clear()
            self._series[u] = (b, db, ddb, self.on_nodes(b)[3], self.project @ (dK @ db))
        return self._series[u]

    def jet(self, u, theta):
        """Everything the published tensors read at retarded time u and angles theta: f and K with
        their derivatives in theta, d_u f with its derivatives in theta, d_u^2 f and d_u K, each
        from the series itself."""
        from numpy.polynomial import legendre as L
        theta = np.asarray(theta, dtype=float)
        x, s = np.cos(theta), np.sin(theta)
        b, db, ddb, k, dk = self.series(u)

        def along(c):
            g, gx, gxx = L.legval(x, c), L.legval(x, L.legder(c)), L.legval(x, L.legder(c, 2))
            return g, -s * gx, s ** 2 * gxx - x * gx
        f, ft, ftt = along(b)
        fu, fut, futt = along(db)
        K, Kt, Ktt = along(k)
        return {"f": f, "ft": ft, "ftt": ftt, "fu": fu, "fut": fut, "futt": futt, "fuu": L.legval(x, ddb),
                "K": K, "Kt": Kt, "Ktt": Ktt, "Ku": L.legval(x, dk)}

    def shape(self, u):
        """f and d_theta f at retarded time u as functions of theta, for a front's surface."""
        from numpy.polynomial import legendre as L
        b = self.series(u)[0]
        db = L.legder(b)
        return (lambda th: L.legval(np.cos(th), b), lambda th: -np.sin(th) * L.legval(np.cos(th), db))

    JET = ("f", "ft", "ftt", "fu", "fut", "futt", "fuu", "K", "Kt", "Ktt", "Ku")

    def table(self):
        """The jet on the solver's own theta as cubic splines in u, on a grid fine where the
        short modes of the first front die out, so that a plane's worth of points is cheap."""
        if not hasattr(self, "_table"):
            from scipy.interpolate import CubicSpline
            grid = np.concatenate([[0.0], np.geomspace(1e-7, self.end, 6000)])
            rows = [self.jet(u, self.theta) for u in grid]
            self._table = {name: CubicSpline(grid, np.array([float(row[name]) for row in rows])) for name in self.JET}
        return self._table

    def values(self, x0, r, exact=False):
        """Every derivative of f and H a published component can name, at chart points (x^0, r)
        on the solver's theta, keyed by the function and its orders in (u, theta) or
        (u, r, theta); NaN before the first front. H = (K - 2r d_u ln f - 2/r)/2. With `exact`
        each point's jet is taken from the series itself and the splines are left out, which
        is what a check of a few points wants."""
        shape = np.broadcast(np.asarray(x0), np.asarray(r)).shape
        u = np.broadcast_to(np.asarray(x0, dtype=float), shape)
        r = np.broadcast_to(np.asarray(r, dtype=float), shape)
        if exact:
            jets = [self.jet(v, self.theta) for v in u.ravel()]
            j = {name: np.array([float(jet[name]) for jet in jets]).reshape(shape) for name in self.JET}
        else:
            table = self.table()
            at = np.clip(u, 0.0, self.end)
            j = {name: np.where(u >= 0, table[name](at), np.nan) for name in self.JET}
        f = j["f"]
        ln_u = j["fu"] / f
        ln_ut = j["fut"] / f - j["fu"] * j["ft"] / f ** 2
        ln_utt = (j["futt"] / f - 2 * j["fut"] * j["ft"] / f ** 2 - j["fu"] * j["ftt"] / f ** 2
                  + 2 * j["fu"] * j["ft"] ** 2 / f ** 3)
        ln_uu = j["fuu"] / f - j["fu"] ** 2 / f ** 2
        with np.errstate(all="ignore"):
            return {
                ("f", (0, 0)): f, ("f", (0, 1)): j["ft"], ("f", (0, 2)): j["ftt"], ("f", (1, 0)): j["fu"],
                ("f", (1, 1)): j["fut"], ("f", (1, 2)): j["futt"], ("f", (2, 0)): j["fuu"],
                ("H", (0, 0, 0)): (j["K"] - 2 * r * ln_u - 2 / r) / 2,
                ("H", (0, 1, 0)): -ln_u + 1 / r ** 2, ("H", (0, 2, 0)): -2 / r ** 3,
                ("H", (0, 0, 1)): (j["Kt"] - 2 * r * ln_ut) / 2, ("H", (0, 0, 2)): (j["Ktt"] - 2 * r * ln_utt) / 2,
                ("H", (0, 1, 1)): -ln_ut, ("H", (1, 0, 0)): (j["Ku"] - 2 * r * ln_uu) / 2,
            }

    @staticmethod
    def placeholders(reader, expr):
        """expr with H, f and every derivative of either replaced by a plain symbol, and the
        symbols with the key of values() each stands for."""
        found = {}

        def symbol(fn, orders):
            key = (fn.func.__name__, tuple(orders))
            if key not in found:
                found[key] = sp.Symbol("_" + key[0] + "_" + "".join(map(str, key[1])))
            return found[key]
        declared = [reader.parameters[name] for name in ("H", "f")]
        for d in sorted(expr.atoms(sp.Derivative), key=lambda d: -len(d.variables)):
            if d.expr in declared:
                expr = expr.xreplace({d: symbol(d.expr, [d.variables.count(a) for a in d.expr.args])})
        for fn in declared:
            expr = expr.xreplace({fn: symbol(fn, [0] * len(fn.args))})
        return expr, list(found.items())

    LATE = 10.0     # the fronts are round to a part in 10^9 from here on

    def kruskal_v(self, u, r):
        """Kruskal's V of the ingoing ray through (u, r) on the solver's theta, for u >= 0: the
        ray is carried by dr/du = -H, every ray at once and each over its own stretch of u, to the
        retarded time LATE, where the plane is Schwarzschild's and V = (r/2 - 1) e^{(u + 2r)/4}."""
        u, r = np.atleast_1d(np.asarray(u, dtype=float)), np.atleast_1d(np.asarray(r, dtype=float))
        shape = np.broadcast(u, r).shape
        u, r = np.broadcast_to(u, shape).ravel(), np.broadcast_to(r, shape).ravel()
        span = np.maximum(self.LATE - u, 0.0)

        def rate(s, y):
            return -span * self.values(u + s * span, y)[("H", (0, 0, 0))]
        end = solve_ivp(rate, (0.0, 1.0), r, method="DOP853", rtol=1e-12, atol=1e-13).y[:, -1]
        at = np.maximum(u, self.LATE)
        return ((end / 2 - 1) * np.exp((at + 2 * end) / 4)).reshape(shape)

    def area(self, u):
        """The area of the front r = 1 over 4 pi, which the equation conserves."""
        f = self.on_nodes(self.modes_at(u))[0]
        x, w = np.polynomial.legendre.leggauss(4 * self.n)
        return float(np.sum(w / f ** 2) / 2)

    def bondi_mass(self, u):
        """Singleton's Bondi mass over m, the mean of f^-3 over the sphere."""
        f = self.on_nodes(self.modes_at(u))[0]
        x, w = np.polynomial.legendre.leggauss(4 * self.n)
        return float(np.sum(w / f ** 3) / 2)

    def radiated(self):
        """Macedo and Saa's closed form for the fraction of the first Bondi mass radiated."""
        e = self.epsilon
        return 1 - math.sqrt((1 - e ** 2) * math.log((1 + e) / (1 - e)) ** 3 / (8 * e ** 3))


def front_solver(fronts, theta=0.0):
    key = ("fronts", json.dumps(fronts, sort_keys=True), float(theta))
    if key not in _SOLVERS:
        _SOLVERS[key] = FrontSolver(float(number(fronts["epsilon"])), float(theta))
    return _SOLVERS[key]



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

def smoothed(expr, pulse):
    """expr with every Dirac delta, and every derivative of one, replaced by the declared pulse of
    the same argument and its derivatives. A pp-wave's profile may be any function of u, so the
    pulse is another exact solution, and the null curves of a plane that cross it are moved
    along it by the integral of the pulse, which is 1, exactly as across the delta."""
    s = sp.Symbol("s", real=True)
    body = sp.sympify(pulse, locals={"s": s})
    return expr.replace(lambda e: isinstance(e, sp.DiracDelta),
                        lambda e: sp.diff(body, s, int(e.args[1]) if len(e.args) > 1 else 0).subs(s, e.args[0]))


def numeric_modules(expr):
    """What lambdify evaluates an expression with: numpy, and scipy as well where a declared
    function holds a Bessel function, as Gowdy's wave does, which numpy does not have."""
    held = any(sp.sympify(e).has(sp.besselj, sp.bessely) for e in (expr if isinstance(expr, (list, tuple)) else [expr]))
    return ["scipy", "numpy"] if held else "numpy"


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
                       else star_solver(spec.metric, spec.system, spec.star) if spec.star
                       else front_solver(spec.fronts, float(number(spec.fixed["theta"]))) if spec.fronts else None)

        def prep(expr):
            expr = sp.sympify(expr)
            if reader.held:
                # A name the reader holds as a function, as the photon rocket's p, is written out.
                expr = expr.subs(reader.held).doit()
            for name, rep in (spec.functions or {}).items():
                expr = expr.replace(reader.parameters[name].func, _as_lambda(reader, name, rep)).doit()
            if spec.delta:
                expr = smoothed(expr, spec.delta)
            return expr.subs(subs)

        g = published_matrix(reader, entry, "metric_components")
        gi = published_matrix(reader, entry, "inverse_metric_components")
        self.fn = {}
        # A principal view's null directions come from the principal plane instead of the
        # coordinate plane; its markers still come from the published inverse metric.
        drawn = (() if spec.principal else (("g00", g[a, a]), ("g0r", g[a, b]), ("grr", g[b, b])))
        self.quotient = Quotient(self, g, prep, a, b) if spec.quotient else None
        if self.quotient:
            drawn = self.quotient.drawn
        for name, expr in drawn + (("gi00", gi[a, a]), ("gi0r", gi[a, b]), ("girr", gi[b, b])):
            self.fn[name] = self.lambdify(prep(expr))
        if spec.mark_gtt:
            self.fn["gtt"] = self.lambdify(prep(g[0, 0]))
        self.principal = PrincipalPlane(self, g, prep, a, b) if spec.principal else None
        K = prep(reader(strip_lhs(entry["kretschmann"]))) if spec.kretschmann else sp.Integer(0)
        # Exponentials are gathered into one, so that e^(-4m/R) e^(2m^2 rho^2/R^4) overflows to
        # infinity near R = 0 and never to zero times infinity.
        self.fn["K"] = self.lambdify(sp.powsimp(K, combine="exp"))
        if spec.singular_zero:
            self.K_expr = K
            self.zero_expr = prep(sp.sympify(spec.singular_zero, locals={**reader.local, **by_plain}))
            self.fn["szero"] = self.lambdify(self.zero_expr)
        if spec.where:
            self.fn["where"] = self.lambdify(prep(sp.sympify(spec.where, locals={**reader.local, **DECLARED_FUNCTIONS,
                                                                                **by_plain})))
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
        if self.solver and not isinstance(self.solver, FrontSolver):
            for name in self.solver.funcs:
                fn = self.reader.parameters[name]
                var = fn.args[0]
                A0, A1, A2 = sp.symbols(f"{name}_0 {name}_1 {name}_2")
                expr = expr.subs(sp.Derivative(fn, (var, 2)), A2).subs(sp.Derivative(fn, var), A1).subs(fn, A0)
                extra += [A0, A1, A2]
        fronts = isinstance(self.solver, FrontSolver)
        solver, fixed_vals = self.solver, self.fixed_vals
        if fronts:
            # A published scalar that divides by sin(theta), as the Kretschmann scalar does, is
            # smooth across the axis and is read OFF_AXIS beside it, with the fronts' own values
            # there; the metric on the plane has no such division and is read on the axis itself.
            theta = self.reader.symbol["\\theta"]
            here = self.fixed_vals[self.fixed_syms.index(theta)]
            beside = min(max(here, FrontSolver.OFF_AXIS), math.pi - FrontSolver.OFF_AXIS)
            if beside != here and expr.has(theta):
                solver = front_solver(self.spec.fronts, beside)
                fixed_vals = [beside if sym == theta else v for sym, v in zip(self.fixed_syms, self.fixed_vals)]
            expr, held = FrontSolver.placeholders(self.reader, expr)
            extra = [symbol for _, symbol in held]
        f = sp.lambdify((self.x0, self.xr, *self.fixed_syms, *extra), expr, numeric_modules(expr))

        def call(x0, r):
            x0 = np.asarray(x0, dtype=float)
            r = np.asarray(r, dtype=float)
            shape = np.broadcast(x0, r).shape
            args = [np.full(shape, v) for v in fixed_vals]
            if fronts:
                values = solver.values(x0, r)
                args += [values[key] for key, _ in held]
            elif self.solver:
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
        """Where the chart's published domain leaves off inside a star's surface, or where the
        expression a row declares in `where` is not positive."""
        gone = np.zeros(np.broadcast(np.asarray(x0), np.asarray(r)).shape, dtype=bool)
        with np.errstate(invalid="ignore"):
            if self.surface is not None:
                gone = gone | (np.asarray(r) < self.surface(np.asarray(x0, dtype=float)))
            if self.spec.where:
                gone = gone | ~(self.fn["where"](x0, r) > 0)
        return gone

    def null_dirs(self, x0, r):
        """The directions P and M in (dx^0, dr), and D. NaN where there are none, and inside a
        star's surface, which the chart does not cover."""
        if self.surface is not None or self.spec.where:
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
        # A view drawn against a time that runs down the chart, as Gowdy's tau does, mirrors the
        # plane, and there the family drawn moving left is the chart's M.
        show = self.spec.to_display
        if show != POLAR and show[0][0] * show[1][1] - show[0][1] * show[1][0] > 0:
            return M, P, D
        return P, M, D

    def orient(self, x0, r, P, M):
        """Signs making P and M future directed at one point, or (0, 0) where none do."""
        at = (np.array([x0]), np.array([r]))
        g00, g0r, grr = (value[0] for value in self.block(*at))

        def dot(p, q):
            return g00 * p[0] * q[0] + g0r * (p[0] * q[1] + p[1] * q[0]) + grr * p[1] * q[1]

        mode = self.spec.orient
        if mode == "split":
            mode = "ingoing" if r < self.spec.split else "outgoing"
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
    body = sp.sympify(rep, locals={**reader.local, **DECLARED_FUNCTIONS, **{str(arg): arg for arg in fn.args}})
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

class Quotient:
    """A plane of x^0 and r with the coordinate k, which the metric does not depend on, divided
    out: its metric h_ab = g_ab - g_ak g_bk/g_kk is the one the header's "Rays of no angular
    momentum" gives, and a null direction of h, lifted by dx^k = -(g_ka dx^a + g_kb dx^b)/g_kk,
    is null in the whole spacetime with no momentum along k."""

    def __init__(self, chart, g, prep, a, b):
        spec, entry, reader = chart.spec, chart.entry, chart.reader
        coords = entry["coords"]
        plain = [reader._plain(c) for c in coords]
        self.name, self.a, self.b = key(spec), a, b
        drawn = [reader._plain(spec.plane[0]), reader._plain(spec.plane[1])]
        if sorted(drawn + list(spec.fixed) + [spec.quotient]) != sorted(plain) or spec.dust or spec.principal:
            raise SystemExit(f"{self.name}: a view that divides out a coordinate needs every coordinate "
                             "drawn, held fixed or divided out exactly once, and no dust")
        self.k = k = plain.index(spec.quotient)
        n = len(coords)
        gamma = entry["christoffel"]["variants"]["ull"]["nonzero"]
        self.gamma_index = [tuple(coords.index(x) for x in comp["indices"]) for comp in gamma]
        exprs = [prep(g[i, j]) for i in range(n) for j in range(n)]
        gammas = [prep(reader(comp["value"])) for comp in gamma]
        along = reader.symbol[coords[k]]
        for expr in exprs + gammas:
            if along in expr.free_symbols:
                raise SystemExit(f"{self.name}: the published metric depends on {spec.quotient}, "
                                 "which the view divides out")
        self.drawn = tuple((name, g[i, j] - g[i, k] * g[j, k] / g[k, k])
                           for name, (i, j) in (("g00", (a, a)), ("g0r", (a, b)), ("grr", (b, b))))
        args = (chart.x0, chart.xr, *chart.fixed_syms)
        self.f = sp.lambdify(args, exprs, "numpy")
        self.f_gamma = sp.lambdify(args, gammas, "numpy")
        self.fixed_vals, self.n = chart.fixed_vals, n

    def _values(self, f, x0, r):
        x0, r = np.broadcast_arrays(np.asarray(x0, dtype=float), np.asarray(r, dtype=float))
        args = [np.full(x0.size, v) for v in self.fixed_vals]
        with np.errstate(all="ignore"):
            out = f(x0.ravel(), r.ravel(), *args)
        return np.array([np.broadcast_to(np.asarray(v, dtype=float), (x0.size,)) for v in out])

    def metric(self, x0, r):
        """The published g_ab at the points, stacked along the first axis."""
        v = self._values(self.f, x0, r)
        return v.T.reshape(v.shape[1], self.n, self.n)

    def christoffel(self, x0, r):
        """The published Gamma^a_bc at the points, stacked along the first axis."""
        v = self._values(self.f_gamma, x0, r)
        G = np.zeros((v.shape[1], self.n, self.n, self.n))
        for (i, j, k), value in zip(self.gamma_index, v):
            G[:, i, j, k] = value
        return G

    def lift(self, direction, g):
        """A direction (dx^0, dr) of the plane as a vector of the spacetime with no momentum
        along k, the coordinates held fixed not moving."""
        K = np.zeros(direction.shape[:1] + (self.n,))
        K[:, self.a], K[:, self.b] = direction[:, 0], direction[:, 1]
        K[:, self.k] = -(g[:, self.k, self.a] * K[:, self.a] + g[:, self.k, self.b] * K[:, self.b]) / g[:, self.k, self.k]
        return K


def quotient_checks(chart, n=241):
    """A view that divides out a coordinate, against what the drawing does not use.

    At n points along the drawn r, away from where g^rr vanishes, each family's direction is
    lifted into the spacetime, scaled to k^r = 1, and checked null against the published metric
    and to satisfy k^b nabla_b k^a = lambda k^a with the published Christoffel symbols,
    derivatives along the two drawn coordinates taken by central differences and none along
    the coordinate divided out, which nothing depends on; the orbits of that coordinate are
    checked spacelike, g_kk > 0. Returns the worst of each and refuses the view past GEODESIC.
    """
    spec, q = chart.spec, chart.quotient
    r = np.linspace(spec.box[0], spec.box[1], n)[1:-1]
    x0 = np.full_like(r, 0.5)
    with np.errstate(all="ignore"):
        keep = np.abs(chart.fn["girr"](x0, r)) > 1e-2
    x0, r = x0[keep], r[keep]
    step = 1e-5 * np.maximum(1.0, np.abs(r))

    def tangents(t, rr):
        g = q.metric(t, rr)
        P, M, _ = chart.null_dirs(t, rr)
        return [K / K[:, q.b:q.b + 1] for K in (q.lift(d, g) for d in (P, M))]

    here = tangents(x0, r)
    along = [[(p - m) / (2 * step[:, None]) for p, m in zip(tangents(*plus), tangents(*minus))]
             for plus, minus in (((x0 + step, r), (x0 - step, r)), ((x0, r + step), (x0, r - step)))]
    g, G = q.metric(x0, r), q.christoffel(x0, r)
    if not np.all(g[:, q.k, q.k] > 0):
        raise SystemExit(f"{key(spec)}: the orbits of {spec.quotient} are not spacelike everywhere on the view")
    geodesic, null, finite = [], [], np.ones(r.size, dtype=bool)
    for f, K in enumerate(here):
        dK = K[:, q.a:q.a + 1] * along[0][f] + K[:, q.b:q.b + 1] * along[1][f]
        GKK = np.einsum("nabc,nb,nc->na", G, K, K)
        acc = dK + GKK
        across = acc - (np.sum(acc * K, 1) / np.sum(K * K, 1))[:, None] * K
        # A family whose tangent is exactly parallel, as v = const of an Eddington-Finkelstein chart
        # is, has both terms zero, and misses by nothing; the smallest float keeps that 0/0 a zero.
        # Where both terms vanish at one point only, as at the throat of Teo's wormhole, where
        # dr/dl = 0, they are rounding there, so the miss is measured against no less than 1e-9
        # of k.k, a hundred times the error of the central differences.
        size = (np.linalg.norm(dK, axis=1) + np.linalg.norm(GKK, axis=1) + 1e-9 * np.sum(K * K, 1)
                + np.finfo(float).tiny)
        geodesic.append(np.linalg.norm(across, axis=1) / size)
        null.append(np.abs(np.einsum("nab,na,nb->n", g, K, K)) / (np.abs(g).max((1, 2)) * np.sum(K * K, 1)))
        finite &= np.isfinite(geodesic[-1]) & np.isfinite(null[-1])
    if finite.sum() < 0.8 * r.size:
        raise SystemExit(f"{key(spec)}: the checks of the divided out view are finite at {finite.sum()} of {r.size} points")
    geodesic = float(max(np.max(v[finite]) for v in geodesic))
    null = float(max(np.max(v[finite]) for v in null))
    if not (geodesic <= GEODESIC and null <= 1e-12):
        raise SystemExit(f"{key(spec)}: the lifted rays miss the geodesic equation by {geodesic:.1e} "
                         f"or the null condition by {null:.1e}")
    return {"points": int(finite.sum()), "geodesic": geodesic, "null": null}


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
        if callable(lo) or callable(hi):
            raise SystemExit(f"{key(self.c.spec)}: inside=True needs a domain of fixed ends")
        margin = 1e-3 * float(np.max(self.span))
        return (None if lo is None else lo + margin, None if hi is None else hi - margin)

    def beyond(self, u):
        if self.c.spec.ends_on_singular:
            x0, r = self.to_chart(self.from_unit(u))
            if float(self.c.fn["szero"](np.array([x0]), np.array([r]))[0]) <= 0:
                return True
        if self.edge is None:
            return False
        _, r = self.to_chart(self.from_unit(u))
        lo, hi = self.edge
        return bool((lo is not None and r < lo) or (hi is not None and r > hi))

    # ---- rays

    def trace(self, u0, family, sign):
        # Every ray may run as far as 5000 steps of the usual 0.0025 take it.
        h = self.c.spec.step
        steps = round(5000 * 0.0025 / h)
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
        if self.c.surface is not None or self.c.spec.where:
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
        """Edges along which the published Kretschmann scalar diverges. A row with
        singular_where_claimed judges each edge only on its part inside the published domains,
        for a chart whose formula runs on past its domain into no part of the spacetime, as
        Khan and Penrose's cosmological chart does beyond the wave fronts, where it diverges on
        sigma = +-pi/2; every other row judges the whole edge, since Schwarzschild's r = 0 is
        a singularity of the spacetime although the spherical chart claims only r > r_s."""
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
            lorentzian = True if self.c.principal else self.lorentzian_near(at)
            claimed = (self.claimed(*self.to_chart(self.from_unit(at(1e-4)))) if self.c.spec.singular_where_claimed
                       else np.ones(t.shape, dtype=bool))
            with np.errstate(all="ignore"):
                # A scalar that grows as e^(2m^2/rho^2), as the Curzon-Chazy particle's does in its plane
                # z = 0, is past the largest double this close to the edge, and reads as infinite.
                growing = (near / far > 50) | np.isposinf(near)
                # An edge the domains touch at one corner alone, as the tip of Roberts's wedge, is no edge
                # of the spacetime: as in singular_runs, a stretch is more than two points long.
                if claimed.sum() > 2 and ((near > 1e8) & growing & lorentzian)[claimed].mean() > 0.5:
                    out.append(name)
        return out

    def lorentzian_near(self, at):
        """Whether the plane's metric is Lorentzian beside an edge, at 1e-4 of the drawing from it,
        or, where a component is past the range of a double there, as the Curzon-Chazy particle's
        e^(2m/rho - m^2/rho^2) is beside its ring, where it reads as zero, at the nearest of 1e-3 and
        1e-2 where it is not."""
        with np.errstate(all="ignore"):
            signs = [self.c.null_dirs(*self.to_chart(self.from_unit(at(e))))[2] for e in (1e-4, 1e-3, 1e-2)]
        out = signs[0]
        for further in signs[1:]:
            out = np.where(np.isfinite(out) & (out != 0), out, further)
        return out > 0

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

    def weak_singularity(self):
        """The zero set a row declares singular in `singular_zero`, as lines of the unit square,
        checked to be a curvature singularity at a dozen points of each line. McVittie's
        Kretschmann scalar diverges on R = r_s only as (dH/dt)^2/(1 - r_s/R), and dH/dt falls
        off exponentially, so at 1e-5 of the drawing it is nowhere near 1e8. The test is the
        header's own, the scalar past 1e8 and growing fiftyfold, taken at 1e-20 and 1e-30 of the
        chart's unit from the curve, on the side the published domains claim, in 60 digits."""
        import mpmath
        lines = self.zero_set("szero")
        if not lines:
            raise SystemExit(f"{key(self.c.spec)}: {self.c.spec.singular_zero} = 0 nowhere in the drawing")
        fixed = dict(zip(self.c.fixed_syms, self.c.fixed_vals))
        K = sp.lambdify((self.c.x0, self.c.xr), self.c.K_expr.subs(fixed), "mpmath")
        Z = sp.lambdify((self.c.x0, self.c.xr), self.c.zero_expr.subs(fixed), "mpmath")
        with mpmath.workdps(60):
            for line in lines:
                for u in np.asarray(line)[1:-1:max(1, len(line) // 12)]:
                    x0, r = (float(v) for v in self.to_chart(self.from_unit(u)))
                    # The curve is met along r, or along x^0 where it runs level in r there, as
                    # Kastor and Traschen's does above the midpoint of their two holes.
                    try:
                        root = mpmath.findroot(lambda q: Z(mpmath.mpf(x0), q), mpmath.mpf(r))

                        def off(d, root=root, x0=x0):
                            return mpmath.mpf(x0), root + d
                    except ValueError:
                        root = mpmath.findroot(lambda q: Z(q, mpmath.mpf(r)), mpmath.mpf(x0))

                        def off(d, root=root, r=r):
                            return root + d, mpmath.mpf(r)
                    step = 1e-3 * float(np.max(self.span))
                    side = next((k for k in (1, -1) if self.claimed(*(float(v) for v in off(k * step)))), None)
                    if side is None:
                        raise SystemExit(f"{key(self.c.spec)}: neither side of the declared singularity at "
                                         f"{x0, r} lies in the published domains")
                    near, far = (abs(K(*off(side * mpmath.mpf(d)))) for d in ("1e-30", "1e-20"))
                    if not (near > 1e8 and near / far > 50):
                        raise SystemExit(f"{key(self.c.spec)}: the Kretschmann scalar does not diverge on the "
                                         f"declared singularity at {x0, r}: {float(far):.1e}, {float(near):.1e}")
        return lines

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
            closest = self.c.spec.singular_near
            near, far = (np.abs(self.c.fn["K"](*self.to_chart(self.from_unit(at(e))))) for e in (closest, 10 * closest))
            with np.errstate(all="ignore"):
                hit = (near > 1e8) & (near / far > 50)
            if self.c.spec.singular_where_claimed:
                # A stretch of an edge is singular only where the published domains reach it.
                hit &= self.claimed(*self.to_chart(self.from_unit(at(1e-4))))
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
        ways or, with toward="past" or "future", only that way. The point is an x^0, meaning the
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
        way = {"past": -1, "future": 1}[toward]
        return [self.trace(seed, f, way * signs[f]) for f in families]

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

    def claimed(self, x0, r, domains=None):
        """Where chart points lie inside the entry's published domains of the plane's two
        coordinates."""
        domains = parse_domains(self.c) if domains is None else domains
        x0, r = np.broadcast_arrays(np.asarray(x0, dtype=float), np.asarray(r, dtype=float))
        outside = np.zeros(x0.shape, dtype=bool)
        with np.errstate(invalid="ignore"):
            for name, values, others in ((self.c.spec.plane[0], x0, r), (self.c.spec.plane[1], r, x0)):
                if name in domains and name not in self.c.spec.periodic:
                    lo, hi, lo_open, hi_open = domains[name]
                    lo, hi = at(lo, others), at(hi, others)
                    if lo is not None:
                        outside |= (values < lo) | ((values == lo) & lo_open)
                    if hi is not None:
                        outside |= (values > hi) | ((values == hi) & hi_open)
            if self.c.spec.where:
                outside |= ~(self.c.fn["where"](x0, r) > 0)
        return ~outside

    def hatch(self):
        """Where a chart point lies outside the entry's published domains, as polygons."""
        domains = parse_domains(self.c)
        if not domains and self.c.surface is None and not self.c.spec.crunch and not self.c.spec.where:
            return []
        UU, VV, x0, r = self.grid(161)
        outside = ~self.claimed(x0, r, domains)
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
        # Beyond a declared singular curve there is no spacetime, so nothing is marked there.
        here = (lambda x0, r: self.claimed(x0, r)) if spec.singular_zero else (lambda x0, r: True)
        keep = here if spec.singular_zero else (lambda x0, r: self.claimed(x0, r)) if spec.where_is_infinity else None
        lines = [] if spec.any_factor else self.zero_set("girr", keep=keep)
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
            throat = self.zero_set("dRr", keep=lambda x0, r: (fn["R"](x0, r) > 1e-6) & here(x0, r), drop_edge=True)
            if throat and spec.null_radius:
                out.append({"kind": "apparent", "lines": throat})
            elif throat and not spec.no_throat:
                out.append({"kind": "throat", "lines": throat})
            if not self.c.same_as_grr and not spec.null_radius:
                apparent = self.zero_set("grad2", keep=lambda x0, r: (np.abs(fn["dRr"](x0, r)) > 1e-6) & here(x0, r))
                apparent = [l for l in apparent if not same_line(l, throat)]
                if apparent:
                    out.append({"kind": "apparent", "lines": apparent})
        for kind, which, at, legend, *span in spec.lines:
            at = float(number(at))
            ends = [float(number(v)) for v in span[0]] if span else (-1e6, 1e6)
            pairs = [(at, v) for v in ends] if which == "x0" else [(v, at) for v in ends]
            line = clip_unit(np.array([self.to_unit(self.to_display(*pair)) for pair in pairs]))
            if line is not None:
                # A line of a marker the page names itself, as g^rr = 0, carries no legend.
                out.append({"kind": kind, "lines": [rounded(line)], **({"legend": legend} if legend else {})})
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
        elif spec.singular_zero:
            out.append({"kind": "singular", "edges": self.singular_edges(), "lines": self.weak_singularity()})
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


def outer_root(chart, x0=0.0, back=0):
    """The outermost zero of g^rr along the drawn radial range at x^0, where r_+ sits, or with
    `back` the zero that many inward of it, as r_+ lies one inside a cosmological horizon r_c."""
    lo, hi = chart.spec.box[0], chart.spec.box[1]
    if chart.spec.to_display == POLAR:
        lo, hi = 0.0, float(np.hypot(np.max(np.abs(chart.spec.box[:2])), np.max(np.abs(chart.spec.box[2:]))))
    X = np.linspace(lo, hi, 20001)
    f = chart.fn["girr"](np.full_like(X, x0), X)
    crossings = np.where(np.sign(f[:-1]) * np.sign(f[1:]) < 0)[0]
    zeros = np.flatnonzero(f == 0)
    roots = sorted([float(X[i] - f[i] * (X[i + 1] - X[i]) / (f[i + 1] - f[i])) for i in crossings]
                   + [float(X[i]) for i in zeros])
    return roots[-1 - back] if len(roots) > back else None


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
                # Inside a cosmological horizon r_c, the outermost zero, r_+ is the next one in.
                ends.append(outer_root(chart, back=1 if body[cut + 1:].strip() == "r_c" else 0))
            elif s == "r_c":
                ends.append(outer_root(chart))
            else:
                subs = {chart.reader.c: 1}
                subs.update({chart.reader.parameters[k]: number(v) for k, v in chart.spec.params.items()})
                try:
                    ends.append(float(sp.N(chart.reader(s).subs(subs))))
                except (vm.LatexError, TypeError, ValueError):
                    ends.append(bound_along(chart, name, s, subs))
        out[name] = (ends[0], ends[1], interval[0] == "(", interval[-1] == ")")
    return out


def bound_along(chart, name, text, subs):
    """An end of a domain that is a function of the plane's other coordinate, as Milne's
    R \\in [0, cT) inside the light cone, as a function of that coordinate's values; None
    where the end names anything else, and it is then not hatched."""
    other = chart.x0 if name == chart.spec.plane[1] else chart.xr
    try:
        # A declared function in the end, as McVittie's r_s/(4a), takes the row's expression.
        expr = chart.prep(chart.reader(text)).subs(subs)
    except (vm.LatexError, TypeError, ValueError):
        return None
    if expr.free_symbols != {other}:
        # An end that names a coordinate the plane holds fixed, as V/H does on a plane of
        # Kastor and Traschen's two holes, is taken at the plane's value of it.
        expr = expr.subs(dict(zip(chart.fixed_syms, chart.fixed_vals)))
    if expr.free_symbols != {other}:
        return None
    f = sp.lambdify(other, expr, "numpy")

    def end(values):
        # A bound such as sqrt(1 - u^2) has no value past u = 1, where the other bound hatches.
        with np.errstate(invalid="ignore"):
            return np.asarray(f(np.asarray(values, dtype=float)), dtype=float)
    return end


def at(end, values):
    """A domain's end at the plane's other coordinate: a number, or a bound_along."""
    return end(values) if callable(end) else end


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
QUOTIENT_FIELDS = ["christoffel"]
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
        # The scale is G^0_0 where the chart publishes one; a pp-wave's only component is G^v_u, and
        # it is then measured against 1.
        scale_text = next((c["value"] for c in ul if c["indices"] == [chart.entry["coords"][0]] * 2), "0")
        size = chart.lambdify(chart.prep(chart.reader(scale_text)))
        with np.errstate(all="ignore"):
            v, scale = value(x0, r), np.abs(size(x0, r))
        ok = np.isfinite(v) & np.isfinite(scale)
        miss = float(np.max(np.abs(v[ok]) / (scale[ok] + 1)))
        if ok.sum() < 500 or not miss < 1e-9:
            raise SystemExit(f"{key(spec)}: the declared functions miss G^{index[0]}_{index[1]} = 0 by {miss:.1e} "
                             f"at {ok.sum()} finite points")


def front_checks(where, fronts, report=None):
    """Declared Robinson-Trautman fronts against what their construction does not use: every
    published component of the axisymmetric chart's Ricci tensor vanishes on them and its
    published Kretschmann scalar is Schwarzschild's 48/r^6, at points all over u, r and theta; the
    area of a front is conserved; the Bondi mass falls from Macedo and Saa's closed form to m;
    and the last mode to die falls as e^{-2u/m}, Foster and Newman's rate."""
    kept = ("front_checks", json.dumps(fronts, sort_keys=True))
    if kept not in _SOLVERS:
        _SOLVERS[kept] = _front_misses(fronts)
    out = _SOLVERS[kept]
    limits = {"ricci": 1e-8, "kretschmann": 1e-9, "area": 1e-10, "bondi": 1e-10, "settles": 1e-10,
              "falls": 1e-12, "rate": 1e-4}
    if report:
        for name, miss in out.items():
            report(f"{where}: the fronts, {name}", miss, limits[name])
    wrong = {k: v for k, v in out.items() if not v <= limits[k]}
    if wrong:
        raise SystemExit(f"{where}: the declared fronts miss {wrong}")
    return out


def _front_misses(fronts):
    _, entry, reader = load("robinson_trautman", "axisymmetric")
    rng = np.random.default_rng(11)
    out = {"ricci": 0.0, "kretschmann": 0.0}
    names = [reader.symbol[c] for c in entry["coords"]]
    for theta in rng.uniform(0.15, math.pi - 0.15, 12):
        solver = FrontSolver(float(number(fronts["epsilon"])), float(theta))
        u, r = rng.uniform(0.0, 3.0, 40), rng.uniform(0.3, 6.0, 40)
        values = solver.values(u, r, exact=True)

        def at(text):
            expr, held = FrontSolver.placeholders(reader, reader(text).subs(reader.c, 1))
            f = sp.lambdify((*names, *[symbol for _, symbol in held]), expr, "numpy")
            return np.asarray(f(u, r, theta, 0.0, *[values[k] for k, _ in held]), dtype=float) * np.ones_like(u)

        def size(text):
            # The same value with every term taken positive, which is what a miss is measured against.
            expr, held = FrontSolver.placeholders(reader, reader(text).subs(reader.c, 1))
            terms = sp.Add.make_args(sp.expand(sp.numer(sp.together(expr)))), sp.denom(sp.together(expr))
            f = sp.lambdify((*names, *[symbol for _, symbol in held]), [list(terms[0]), terms[1]], "numpy")
            top, bottom = f(u, r, theta, 0.0, *[values[k] for k, _ in held])
            return sum(np.abs(np.asarray(t, dtype=float) * np.ones_like(u)) for t in top) / np.abs(bottom)
        for c in entry["ricci_tensor"]["variants"]["ll"]["nonzero"]:
            out["ricci"] = max(out["ricci"], float(np.max(np.abs(at(c["value"])) / size(c["value"]))))
        K = at(strip_lhs(entry["kretschmann"]))
        out["kretschmann"] = max(out["kretschmann"], float(np.max(np.abs(K * r ** 6 / 48 - 1))))
    solver = front_solver(fronts)
    times = np.concatenate([[0.0], np.geomspace(1e-4, FrontSolver.end, 60)])
    mass = np.array([solver.bondi_mass(t) for t in times])
    out["area"] = float(max(abs(solver.area(t) - 1) for t in times))
    out["bondi"] = abs(mass[0] - 1 / (1 - solver.radiated()))
    out["settles"] = abs(mass[-1] - 1)
    out["falls"] = float(max(0.0, np.max(np.diff(mass))))
    out["rate"] = abs(math.log(abs(solver.modes_at(5.0)[2] / solver.modes_at(6.0)[2])) - 2.0)
    return {k: float(v) for k, v in out.items()}


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
              + (["ricci_tensor"] if spec.fronts else [])
              + (PRINCIPAL_FIELDS if spec.principal else []) + (QUOTIENT_FIELDS if spec.quotient else []))
    if spec.principal:
        principal_checks(chart)
    if spec.quotient:
        quotient_checks(chart)
    if spec.star:
        star_checks(spec, chart.solver)
    if spec.fronts:
        front_checks(key(spec), spec.fronts)
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
    marks = view_slices(spec, view["markers"])
    if marks:
        view["slices"] = marks
    return view


def view_slices(spec, markers=()):
    """The moments of the spacetime's embedding diagram that the view shows, which
    slices.flat declares in the chart's (x^0, r), carried into the unit square by the view's
    own map and box: lines cut where they leave the box and thinned as the rays are, points
    inside it, and regions cut to it, each a list of rings filled by the even odd rule. A
    mirrored view carries the half at r >= 0, as its rays do, and the page draws both."""
    import slices
    X0, X1, Y0, Y1 = spec.box
    lo, span = np.array([X0, Y0], float), np.array([X1 - X0, Y1 - Y0], float)

    def unit(P):
        P = np.atleast_2d(np.asarray(P, float))
        if spec.to_display == POLAR:
            Q = np.column_stack([P[:, 1] * np.cos(P[:, 0]), P[:, 1] * np.sin(P[:, 0])])
        else:
            Q = P @ np.array(spec.to_display, float).T
        return (Q - lo) / span

    out, drawn = [], []
    for mark in slices.flat(spec):
        lines = [rounded(thin(run, 0.0006)) for line in mark.lines for run in slices.clip_runs(unit(line))]
        points = [rounded(u) for p in mark.points for u in unit(p) if np.all((u >= 0) & (u <= 1))]
        fills = []
        for rings in mark.fills:
            cut = [slices.clip_ring(unit(ring)) for ring in rings]
            cut = [rounded(thin(np.vstack([ring, ring[:1]]), 0.0006)[:-1]) for ring in cut if ring is not None]
            if cut:
                fills.append(cut)
        if not (lines or points or fills):
            raise SystemExit(f"{key(spec)}: the moment {mark.label} of the embedding lies outside the box")
        out.append({**mark.moment.json(), "label": mark.label, "lines": lines, "points": points, "fills": fills})
        # A mirrored view draws each line and its mirror, so either end may take the label.
        both = lines + ([[[-u[0], u[1]] for u in line] for line in lines] if spec.mirror else [])
        drawn.append(slices.Mark(mark.moment, both, points, fills, mark.label))
    # The labels are placed on the plot in the units the page draws it in, 520 wide, a mirrored
    # view's unit square its right half, at the size the page sets them at on a desktop at the
    # usual text, the caption's 21px on a plot about 458px wide; where the page draws them
    # larger against the plot, as on a phone, its fitSliceLabels() moves them from there. The
    # reference line's label stands at the plot's right edge above its line.
    W = 520.0
    H = W * (Y1 - Y0) / ((X1 - X0) * (2 if spec.mirror else 1))
    size = 24.0

    def px(u):
        x = (0.5 + 0.5 * u[0]) * W if spec.mirror else u[0] * W
        return np.array([x, (1 - u[1]) * H])
    others = []
    for marker in markers:
        if marker["kind"] == "reference":
            y = (1 - marker["y"]) * H
            w, h = (v * size for v in slices.label_size(marker["label"]))
            others.append((W - 0.45 * size - w, y - 0.35 * size - h, W - 0.45 * size, y - 0.35 * size))
    for entry, place in zip(out, slices.place(drawn, px, (W, H), size, others)):
        if place:
            entry["place"] = place
    return out


def rewrite_slices(metric_ids=None):
    """Only the slices of every diagram file, from the embedding files as they stand, without
    tracing a ray again: each flat view's slices are drawn anew with its row's map and box,
    and each figure, which draws in seconds, is drawn anew whole."""
    import projections as pj
    for metric_id, (specs, figures) in by_metric().items():
        if metric_ids and metric_id not in metric_ids:
            continue
        path = DIAGRAMS_DIR / f"{metric_id}.json"
        data = json.loads(path.read_text(encoding="utf-8"))
        for spec in specs:
            view = next(v for v in data["systems"][spec.system] if v["id"] == spec.view)
            if (list(view["box"]) != list(spec.box) or view["mirror"] != spec.mirror
                    or view["settings"] != settings(spec, load(spec.metric, spec.system)[1])):
                raise SystemExit(f"{key(spec)}: the file was drawn with another box or other values; "
                                 "redraw it whole")
            view.pop("slices", None)
            marks = view_slices(spec, view["markers"])
            if marks:
                view["slices"] = marks
        for spec in figures:
            views = data["projections"][spec.system]
            i = next(i for i, v in enumerate(views) if v["id"] == spec.view)
            views[i] = pj.draw(spec)
        path.write_text(json.dumps(data, ensure_ascii=False, separators=(",", ":"), allow_nan=False) + "\n",
                        encoding="utf-8")
        print(f"wrote the slices of {path.relative_to(build.ROOT)}")


def by_metric():
    """Every spacetime with a diagram file: its flat views and its figures, in table order."""
    import projections as pj
    out = {}
    for spec in DIAGRAMS:
        out.setdefault(spec.metric, ([], []))[0].append(spec)
    for spec in pj.FIGURES:
        out.setdefault(spec.metric, ([], []))[1].append(spec)
    return out


def write(metric_ids=None):
    import projections as pj
    if not DIAGRAMS_DIR.exists():
        DIAGRAMS_DIR.mkdir(parents=True)
    drawn = by_metric()
    if not metric_ids:
        # A spacetime nothing is drawn of has no file, so a full redraw removes one left behind.
        for path in sorted(DIAGRAMS_DIR.glob("*.json")):
            if path.stem not in drawn and not build.CONFLICT_COPY.search(path.stem):
                path.unlink()
                print(f"removed {path.relative_to(build.ROOT)}")
    for metric_id, (specs, figures) in drawn.items():
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
        path = DIAGRAMS_DIR / f"{metric_id}.json"
        path.write_text(json.dumps(data, ensure_ascii=False, separators=(",", ":"), allow_nan=False) + "\n",
                        encoding="utf-8")
        print(f"wrote {path.relative_to(build.ROOT)}")


# ---------------------------------------------------------------- closed forms

def _ks_tau(eta):
    """The conformal time of the Kantowski-Sachs dust at kappa = 0 and b_0 = 1, the integral of
    2 cos^2(s)/(1 + s tan(s)) from 0 to eta, along which r -+ tau is constant on a ray."""
    return np.array([quad(lambda s: 2 * np.cos(s) ** 3 / (np.cos(s) + s * np.sin(s)), 0, e)[0]
                     for e in np.atleast_1d(eta)]).reshape(np.shape(eta))


def _teo_proper(r):
    """Teo's proper radial distance from the throat at b_0 = 1, his eq. (28)."""
    r = np.maximum(np.asarray(r, float), 1.0)
    return np.sqrt(r * (r - 1)) + np.log(np.sqrt(r) + np.sqrt(r - 1))


def _teo_axis(r):
    """r_* on the axis of Teo's example at b_0 = 1 and a = 1/4, zero at the throat:
    dr_*/dr = r^(3/2)/((r + 1) sqrt(r - 1))."""
    r = np.maximum(np.asarray(r, float), 1.0)
    return (np.sqrt(r * (r - 1)) - np.log(np.sqrt(r) + np.sqrt(r - 1))
            + np.sqrt(2) * np.arctanh(np.sqrt((r - 1) / (2 * r))))


def _gravastar_x(r):
    """The tortoise coordinate inside the gravastar the diagrams draw, dx/dr = 1/(sqrt(C)(1 - r^2/L^2))."""
    return 2 / math.sqrt(GRAVASTAR_C) * np.arctanh(np.asarray(r, float) / 2)


def _rstar(r, horizons):
    """The tortoise coordinate of f = prod(1 - r_i/r) with simple roots r_i, up to a constant."""
    out = np.asarray(r, float).copy()
    for i, ri in enumerate(horizons):
        others = [rj for j, rj in enumerate(horizons) if j != i]
        coefficient = ri ** len(horizons) / np.prod([ri - rj for rj in others]) if others else ri
        out = out + coefficient * np.log(np.abs(r - ri))
    return out


def _curzon_axis(z):
    """z_* on the axis of the Curzon-Chazy particle at m = 1, dz_*/dz = e^(2/z)."""
    z = np.asarray(z, float)
    return z * np.exp(2 / z) - 2 * scipy_expi(2 / z)


def _curzon_plane(rho):
    """rho_* in the plane z = 0 of the Curzon-Chazy particle at m = 1, the integral from 0 of
    e^(2/s - 1/(2 s^2))."""
    def one(r):
        return quad(lambda s: math.exp(2 / s - 0.5 / s ** 2) if s > 0 else 0.0, 0, r, epsabs=1e-13, epsrel=1e-13,
                    limit=200)[0]
    return np.vectorize(one, otypes=[float])(np.asarray(rho, float))


def _zv_axis(r, oblate):
    """r_* on the axis of Zipoy and Voorhees's metric at m = 1, dr_*/dr = f^(-1 - q): at q = 1,
    r + 4 ln(r - 2) - 4/(r - 2), and at q = -1/2, sqrt(r (r - 2)) + 2 ln(sqrt r + sqrt(r - 2))."""
    r = np.asarray(r, float)
    if oblate:
        return r + 4 * np.log(r - 2) - 4 / (r - 2)
    return np.sqrt(r * (r - 2)) + 2 * np.log(np.sqrt(r) + np.sqrt(r - 2))


def _zv_equator(r, oblate):
    """r_* in the equatorial plane at m = 1, the integral from 2 of f^(-1 - q) h^(-q (2 + q)/2):
    s^(7/2)/(sqrt(s - 2) (s - 1)^3) at q = 1 and (1 - 2/s)^(-7/8) (1 - 1/s)^(3/4) at q = -1/2."""
    alpha, g = ((-0.5, lambda s: s ** 3.5 / (s - 1) ** 3) if oblate
                else (-0.875, lambda s: s ** 0.875 * (1 - 1 / s) ** 0.75))

    def one(x):
        return quad(g, 2, x, weight="alg", wvar=(alpha, 0), epsabs=1e-13, epsrel=1e-13, limit=200)[0]
    return np.vectorize(one, otypes=[float])(np.asarray(r, float))


def _zv_forms(star, oblate, shift, edge):
    """t + r_* and t - r_* of one view, the prolate spheroidal chart's x = r/m - 1 moved by `shift`."""
    return (lambda t, r: t + star(r + shift, oblate), lambda t, r: t - star(r + shift, oblate),
            lambda t, r: r + shift > 2 + edge)


def _btz_rstar(r, roots):
    """The integral of 1/N^2 for the BTZ hole at l = 1, N^2 = (r^2 - a)(r^2 - b)/r^2 with the
    squared horizons a > b, or r^2 - a where b = 0, up to a constant."""
    a, b = roots
    one = lambda s: np.sqrt(s) / 2 * np.log(np.abs((r - np.sqrt(s)) / (r + np.sqrt(s))))
    return (one(a) - (one(b) if b else 0)) / (a - b)


def _tangherlini_rstar(r, D=5):
    """Tangherlini's tortoise coordinate at r_h = 1, the integral of 1/(1 - r^(3-D)), up to a constant:
    r + ln|(r - 1)/(r + 1)|/2 in five dimensions, and in six
    r + ln|r - 1|/3 - ln(r^2 + r + 1)/6 - arctan((2r + 1)/sqrt 3)/sqrt 3."""
    r = np.asarray(r, float)
    if D == 5:
        return r + 0.5 * np.log(np.abs((r - 1) / (r + 1)))
    return (r + np.log(np.abs(r - 1)) / 3 - np.log(r * r + r + 1) / 6
            - np.arctan((2 * r + 1) / np.sqrt(3)) / np.sqrt(3))


def _kk_rstar(r):
    """The Kaluza-Klein monopole's tortoise coordinate at m = 1, the integral of sqrt(1 + 4/r) from the
    nut: sqrt(r(r + 4)) + 4 arsinh(sqrt(r)/2)."""
    r = np.maximum(np.asarray(r, float), 0.0)
    return np.sqrt(r * (r + 4)) + 4 * np.arcsinh(np.sqrt(r) / 2)


def _kkbh_rstar(r, a=1.0):
    """The tortoise coordinate of a Kaluza-Klein black hole of one charge at r_s = 1, with
    a = q - r_s or p - r_s: the integral of sqrt(r(r + a))/(r - 1) from r = 0, slices.kkbh_rstar."""
    return slices.kkbh_rstar(r, a)


def _myers_rstar(r, D=5):
    """The tortoise coordinate on the plane transverse to the rotation of Myers and Perry's black hole
    with one spin. In five dimensions at mu = 1 and a = 3/5, 1/f = 1 + mu/(r^2 - r_+^2) with r_+ = 4/5:
    r + (5/8) ln|(r - 4/5)/(r + 4/5)|. In six at mu = 1 and a = 3/2, 1/f = 1 + mu/(r^3 + a^2 r - mu):
    r + Re sum_i A_i ln(1 - r/r_i) over the three roots of the cubic, A_i = mu/(3 r_i^2 + a^2)."""
    r = np.asarray(r, float)
    if D == 5:
        return r + 0.625 * np.log(np.abs((r - 0.8) / (r + 0.8)))
    roots = np.roots([1.0, 0.0, 2.25, -1.0]).astype(complex)
    return r + np.real(sum(np.log(1 - r.astype(complex) / ri) / (3 * ri ** 2 + 2.25) for ri in roots))


def _away(*radii):
    return lambda x0, r: np.all([np.abs(r - h) > 0.05 for h in radii], axis=0)


BTZ_STATIC, BTZ_ROTATING = (1.0, 0.0), (0.8, 0.2)


def _sds_rstar(r):
    """Kottler's tortoise coordinate at r_s = 1 and Lambda = 1/5, sum_i ln|r - r_i|/f'(r_i) over the
    three roots of Lambda r^3 - 3r + 3r_s, the negative one included, since 1/f has no polynomial part."""
    roots = np.roots([0.2, 0, -3, 3]).real
    return sum(np.log(np.abs(r - ri)) / (1 / ri ** 2 - 0.4 * ri / 3) for ri in roots)


def _rnds_rstar(r):
    """The lukewarm hole's tortoise coordinate at r_s = 1, r_q = 1/2 and Lambda = 27/64,
    sum_i ln|r - r_i|/f'(r_i) over the four roots of 9r^4 - 64r^2 + 64r - 16, the negative one
    included, since 1/f = -64r^2/(9 prod(r - r_i)) has no polynomial part."""
    return sum(np.log(np.abs(r - a)) * (-64 * a * a / (9 * np.prod([a - b for b in RNDS_ROOTS if b != a])))
               for a in RNDS_ROOTS)


def _rnds_away(x, r):
    return np.all([np.abs(r - a) > 0.05 for a in RNDS_ROOTS[:3]], axis=0) & (r > 0.05)


def _rnds_cosmic(tau, rho, sign):
    """What each family keeps on the cosmological plane at r_s = 1 and H = 3/8: with the areal radius
    r = H tau rho + 1/2, the static time is cT = ln|H tau|/H + F(r), F' = H r^2/((r - 1/2) f), a sum
    of logarithms over 1/2 and the four roots of f, and the rays keep cT -+ r_*."""
    r = 0.375 * tau * rho + 0.5
    poles = (0.5,) + RNDS_ROOTS
    F = sum(np.log(np.abs(r - a)) * (-(8 / 3) * a ** 4 / np.prod([a - b for b in poles if b != a])) for a in poles)
    return np.log(np.abs(0.375 * tau)) / 0.375 + F + sign * _rnds_rstar(r)


def _rnds_cosmic_away(tau, rho):
    r = 0.375 * tau * rho + 0.5
    return _rnds_away(tau, r) & (np.abs(r - 0.5) > 0.05) & (np.abs(tau) > 0.05)


def _kt_cosmic(tau, r, sign):
    """What each family keeps on one hole's plane at m = 1 and H = -3/16: the static time of
    slices.kt_static_t at the areal radius R = H tau r + 1, plus or minus the tortoise coordinate."""
    return slices.kt_static_t(tau, r) + sign * slices.kt_rstar(slices.KT_H_ONE * tau * r + 1)


def _kt_cosmic_away(tau, r):
    R = slices.KT_H_ONE * tau * r + 1
    return (np.all([np.abs(R - a) > 0.05 for a in slices.KT_ONE_ROOTS[:3]], axis=0) & (R > 0.05)
            & (np.abs(R - 1) > 0.05) & (np.abs(tau) > 0.05))


def _kt_label(plane, direction, comoving=False):
    """What a ray of a plane of the two holes keeps, there being no closed form: where the ray
    stood at c tau = -12 m, as 12 arctan(x/12), found by following it into the past with scipy's integrator, which
    the tracing never uses. direction is +1 for the family moving right, dx/d(c tau) = 1/U^2,
    and -1 for the one moving left. Into the past U only grows, so every ray is followed there
    without meeting the singularity, and a ray that left a hole closes on it."""
    H = slices.KT_H_TWO

    def one(tau, x):
        ray = solve_ivp(lambda q, y: [direction / (H * q + float(slices.kt_potential(plane, y[0]))) ** 2],
                        [tau, -12.0], [x], rtol=1e-11, atol=1e-12)
        # A ray that reaches the edge of the drawing near tau = 0 stood hundreds of m away then, so
        # the place is named by its arctangent, which keeps the comparison at the drawing's own scale.
        return 12 * math.atan(ray.y[0, -1] / 12)

    def label(time, x):
        time, x = np.asarray(time, dtype=float), np.asarray(x, dtype=float)
        tau = np.exp(H * time) / H if comoving else time
        return np.array([one(a, b) for a, b in zip(tau, x)])
    return label


def _kt_keep(tau, x):
    """A ray of the two holes is compared at every twelfth point, since each costs an integration,
    and before tau = 0: after it U is small, neighbouring rays followed into the past part as
    exp(int 2|dU/dx|/U^3 d(c tau)), and no integrator names the place to a part in 1e5."""
    return (np.arange(len(x)) % 12 == 0) & (tau < 0)


def _kt_every(time, x):
    return np.arange(len(x)) % 12 == 0


def _rnads_away(x, r):
    return (np.abs(r - 1) > 0.05) & (np.abs(r - 0.5) > 0.05) & (r > 0.05)


def _sads_rstar(r):
    """Schwarzschild-anti-de Sitter's tortoise coordinate at r_s = 2 and L = 1, where 1/f =
    r/((r - 1)(r^2 + r + 2)): (1/4) ln|r - 1| - (1/8) ln(r^2 + r + 2) + (5/(4 sqrt 7)) arctan((2r + 1)/sqrt 7),
    up to a constant."""
    w = np.sqrt(7.0)
    return 0.25 * np.log(np.abs(r - 1)) - 0.125 * np.log(r * r + r + 2) + 5 / (4 * w) * np.arctan((2 * r + 1) / w)


def _bardeen_rstar(r):
    """Bardeen's tortoise coordinate at r_s = 1 and g = 1/3, which the slices share."""
    import slices
    return slices.bardeen_rstar(r)


def _bardeen_away(x, r):
    return (np.abs(r - 0.30096) > 0.05) & (np.abs(r - 0.77542) > 0.05)


def _hayward_rstar(r):
    """Hayward's tortoise coordinate at m = 1 and ell = 12/(7 sqrt 7), where 1/F = 1 + 2r^2/((r - 6/7)
    (r - 12/7)(r + 4/7)): r - (6/5) ln|r - 6/7| + 3 ln|r - 12/7| + (1/5) ln(r + 4/7), up to a constant."""
    return r - 1.2 * np.log(np.abs(r - 6 / 7)) + 3 * np.log(np.abs(r - 12 / 7)) + 0.2 * np.log(r + 4 / 7)


def _hayward_away(x, r):
    return (np.abs(r - 6 / 7) > 0.05) & (np.abs(r - 12 / 7) > 0.05)


def _tbh_rstar(case):
    """The topological black holes' tortoise coordinates in units of L, up to a constant. The flat hole at
    mu = 1, 1/f = r/((r - 1)(r^2 + r + 1)): (1/3) ln|r - 1| - (1/6) ln(r^2 + r + 1) + arctan((2r + 1)/sqrt 3)/
    sqrt 3. The hyperbolic hole without mass, 1/f = 1/(r^2 - 1): (1/2) ln|(r - 1)/(r + 1)|. The hyperbolic
    hole at mu = -120/343, the sum of ln|r - r_i|/f'(r_i) over the roots 5/7, 3/7 and -8/7:
    (35/26) ln|r - 5/7| - (21/22) ln|r - 3/7| - (56/143) ln(r + 8/7)."""
    w = np.sqrt(3.0)
    return {
        "flat": lambda r: (np.log(np.abs(r - 1)) / 3 - np.log(r * r + r + 1) / 6
                           + np.arctan((2 * r + 1) / w) / w),
        "massless": lambda r: 0.5 * np.log(np.abs((r - 1) / (r + 1))),
        "negative": lambda r: (35 / 26 * np.log(np.abs(r - 5 / 7)) - 21 / 22 * np.log(np.abs(r - 3 / 7))
                               - 56 / 143 * np.log(r + 8 / 7)),
    }[case]


def _tbh_zstar(z):
    """The black brane's tortoise coordinate at z_h = 1, the integral of 1/(1 - z^3):
    -(1/3) ln|1 - z| + (1/6) ln(z^2 + z + 1) + arctan((2z + 1)/sqrt 3)/sqrt 3."""
    w = np.sqrt(3.0)
    return -np.log(np.abs(1 - z)) / 3 + np.log(z * z + z + 1) / 6 + np.arctan((2 * z + 1) / w) / w


def _tbh_closed():
    """The closed forms of every view of the topological black holes: c t +- r_*, v and v - 2 r_*,
    u + 2 r_* and u, and for the brane c t +- z_*."""
    away = {"flat": _away(1.0), "massless": _away(1.0), "negative": _away(5 / 7, 3 / 7)}
    forms = {}
    for case in away:
        rs = _tbh_rstar(case)
        static = (lambda t, r, rs=rs: t + rs(r), lambda t, r, rs=rs: t - rs(r), away[case])
        forms[("topological_black_hole", "static", case)] = static
        if case != "flat":
            forms[("topological_black_hole", "hyperbolic", case)] = static
        forms[("topological_black_hole", "eddington_finkelstein_ingoing", case)] = (
            lambda v, r: v, lambda v, r, rs=rs: v - 2 * rs(r), away[case])
        forms[("topological_black_hole", "eddington_finkelstein_outgoing", case)] = (
            lambda u, r, rs=rs: u + 2 * rs(r), lambda u, r: u, away[case])
    forms[("topological_black_hole", "black_string", "radial")] = forms[("topological_black_hole", "static", "flat")]
    forms[("topological_black_hole", "brane", "tz")] = (
        lambda t, z: t + _tbh_zstar(z), lambda t, z: t - _tbh_zstar(z), _away(1.0))
    return forms


def _sds_away(x, r):
    return (np.abs(r - 1.0852) > 0.05) & (np.abs(r - 3.2146) > 0.05)


def _kiselev_rstar(r):
    """Kiselev's tortoise coordinate at w = -2/3, r_s = 1 and r_q = 8, where 1/f = -8r/((r - a)(r - b))
    with a, b = 4 -+ 2 sqrt 2: (4 sqrt 2 - 4) ln|r - a| - (4 sqrt 2 + 4) ln|r - b|, up to a constant."""
    a, b = KISELEV_ROOTS
    return (4 * 2 ** 0.5 - 4) * np.log(np.abs(r - a)) - (4 * 2 ** 0.5 + 4) * np.log(np.abs(r - b))


def _kiselev_away(x, r):
    return (np.abs(r - KISELEV_ROOTS[0]) > 0.05) & (np.abs(r - KISELEV_ROOTS[1]) > 0.05)


def _mp_axis(z):
    """The integral of U^2 along the axis for two holes of m = 1 at z = +-2, U = 1 + 1/|z - 2| +
    1/|z + 2|, up to a constant in each of the three intervals the holes cut the axis into:
    z - 1/(z - 2) - 1/(z + 2) + 2 sgn(z - 2) ln|z - 2| + 2 sgn(z + 2) ln|z + 2|, plus
    (1/2) ln|(z + 2)/(2 - z)| between the holes and (1/2) ln|(z - 2)/(z + 2)| outside them."""
    z = np.asarray(z, dtype=float)
    a, b = z - 2, z + 2
    with np.errstate(divide="ignore", invalid="ignore"):
        cross = np.where(np.abs(z) < 2, 0.5 * np.log(np.abs(b / a)), 0.5 * np.log(np.abs(a / b)))
        return z - 1 / a - 1 / b + 2 * np.sign(a) * np.log(np.abs(a)) + 2 * np.sign(b) * np.log(np.abs(b)) + cross


def _iwp_rstar(r):
    """The tortoise coordinate of Israel, Wilson and Perjes's spinning source at m = a = 1, on its axis
    and along its principal null rays: dr_*/dr = ((r + m)^2 + a^2)/(r^2 + a^2), so
    r_* = r + m ln(r^2 + a^2) + (m^2/a) arctan(r/a)."""
    r = np.asarray(r, dtype=float)
    return r + np.log(r * r + 1) + np.arctan(r)


def _iwp_midplane(rho):
    """The quadrature of sqrt(W^4 rho^2 - omega^2)/rho from rho = 1 across the midplane of the two declared
    sources, W = 1 + 2/s and omega = 4/s + 4/s^2 with s = sqrt(rho^2 + 4): what ct keeps, up to its sign,
    along a null geodesic of no angular momentum."""
    from scipy.integrate import quad

    def integrand(x):
        s = math.sqrt(x * x + 4)
        return math.sqrt((1 + 2 / s) ** 4 * x * x - (4 / s + 4 / (s * s)) ** 2) / x

    return np.vectorize(lambda x: quad(integrand, 1.0, x, epsabs=1e-12, epsrel=1e-12)[0])(np.asarray(rho, dtype=float))


def _mp_midplane(x):
    """The integral of U^2 across the midplane of the same two holes, U = 1 + 2/sqrt(x^2 + 4):
    x + 4 arcsinh(x/2) + 2 arctan(x/2)."""
    x = np.asarray(x, dtype=float)
    return x + 4 * np.arcsinh(x / 2) + 2 * np.arctan(x / 2)


# (metric, system, view): (what P conserves, what M conserves, where to compare). None
# where a family has no closed form. P moves toward smaller r or x, M toward larger.
def _sv_rstar(r, a):
    """The tortoise coordinate of Simpson and Visser's black bounce at r_s = 1, dr_*/dr = rho/(rho - r_s)
    with rho = sqrt(r^2 + a^2): r + r_s arsinh(r/a) and a third term, a logarithm of the two horizons
    r = +-h, h = sqrt(r_s^2 - a^2), where a < r_s, -r_s(r_s + rho)/r where a = r_s, and two arctangents
    with k = sqrt(a^2 - r_s^2) where a > r_s. simpson_visser.md derives each."""
    r = np.asarray(r, float)
    rho = np.sqrt(r * r + a * a)
    out = r + np.arcsinh(r / a)
    if a < 1:
        h = math.sqrt(1 - a * a)
        return out + np.log(np.abs((r - h) * (r - h * rho) / ((r + h) * (r + h * rho)))) / (2 * h)
    if a == 1:
        return out - (1 + rho) / r
    k = math.sqrt(a * a - 1)
    return out + (np.arctan(r / k) + np.arctan(r / (k * rho))) / k


def _sv_forms(a, chart):
    """What each family keeps in a chart of the black bounce, and where the check is made: away from
    the horizons, where r_* diverges."""
    h = math.sqrt(max(1 - a * a, 0.0))

    def clear(x0, r):
        return (np.abs(r - h) > 0.05) & (np.abs(r + h) > 0.05) if a <= 1 else np.ones_like(np.asarray(r, float), bool)

    if chart == "spherical":
        return (lambda t, r: t + _sv_rstar(r, a), lambda t, r: t - _sv_rstar(r, a), clear)
    if chart == "areal":
        def star(rho):
            return _sv_rstar(np.sqrt(np.maximum(np.asarray(rho, float) ** 2 - a * a, 0)), a)
        return (lambda t, rho: t + star(rho), lambda t, rho: t - star(rho),
                lambda t, rho: (np.abs(rho - 1) > 0.05) & (rho > a + 0.002))
    if chart == "eddington_finkelstein_ingoing":
        return (lambda v, r: v, lambda v, r: v - 2 * _sv_rstar(r, a), clear)
    return (lambda u, r: u + 2 * _sv_rstar(r, a), lambda u, r: u, clear)


def _ds_rstar(r):
    """The tortoise coordinate of Damour and Solodukhin's wormhole at r_s = 1 and lambda = 1/5, zero at
    the throat: dr_*/dr = r/sqrt((r - r_s)(a r - r_s)) with a = 1 + lambda^2."""
    r = np.asarray(r, float)
    a = 1.04
    root = np.sqrt(np.maximum((r - 1) * (a * r - 1), 0))
    return root / a + (1 + a) / (2 * a ** 1.5) * np.log((2 * a * r - (1 + a) + 2 * math.sqrt(a) * root) / (a - 1))


def _witten_xstar(x):
    """The tortoise coordinate of Witten's black hole in the Schwarzschild gauge at lambda = m = 1:
    dx_*/dx = 1/(1 - e^(-2x)), so x_* = ln|e^(2x) - 1|/2."""
    return np.log(np.abs(np.expm1(2 * np.asarray(x, float)))) / 2


def _fjnw_rstar(r):
    """The tortoise coordinate of Fisher, Janis, Newman and Winicour's metric at gamma = 1/2 and b = 1,
    zero at the singularity: dr_*/dr = (1 - 1/r)^(-1/2)."""
    r = np.asarray(r, float)
    root = np.sqrt(np.maximum(r * (r - 1), 0))
    return root + np.log(np.sqrt(r) + np.sqrt(np.maximum(r - 1, 0)))


def _ds_xstar(rho):
    """The tortoise coordinate of the rescaled time in Bueno and his collaborators' rho, at r_s = 1
    and lambda = 1/5: r_s((2 + lambda^2) rho + lambda^2 sinh rho)/(2(1 + lambda^2)), odd in rho."""
    rho = np.asarray(rho, float)
    return (2.04 * rho + 0.04 * np.sinh(rho)) / 2.08


def _ds_areal(r):
    """The areal radius at the isotropic radius r, r(1 + r_s/4r)^2 at r_s = 1, never below the throat's 1."""
    r = np.asarray(r, float)
    return np.maximum(r * (1 + 1 / (4 * r)) ** 2, 1.0)


def _ds_isotropic(R):
    """The isotropic radius outside the throat of the sphere of areal radius R, at r_s = 1."""
    return (R - 0.5 + math.sqrt(R * (R - 1))) / 2


def _ds_radius(rho):
    """The areal radius at rho, r_s(2 + lambda^2(1 + cosh rho))/(2(1 + lambda^2)), at r_s = 1 and lambda = 1/5."""
    return (2.04 + 0.04 * np.cosh(np.asarray(rho, float))) / 2.08


def cdl_conformal_distance(xi, bubble):
    """X(xi), the integral of d xi/rho across the thin wall of each declared bubble of Coleman
    and De Luccia, l = 1, with psi -+ X constant along the rays of the wall chart. Decay into
    flat space: ln xi inside the wall at 4/5 and ln tan((xi - xi_0)/2) + ln(8/5) outside it, since
    tan(arcsin(4/5)/2) = 1/2. Decay of flat space: ln tanh(xi/2) inside the wall at ln 3, where
    it is ln(1/2), and ln(xi - ln 3 + 4/3) + ln(3/8) outside it."""
    xi = np.asarray(xi, dtype=float)
    with np.errstate(all="ignore"):
        if bubble == "into_flat":
            angle = xi - 0.8 + math.asin(0.8)
            return np.where(xi < 0.8, np.log(xi), np.log(np.tan(angle / 2)) + math.log(1.6))
        return np.where(xi < math.log(3), np.log(np.tanh(xi / 2)),
                        np.log(xi - math.log(3) + 4 / 3) + math.log(3 / 8))


def _tsw_lstar(l):
    """The tortoise coordinate of the thin shell wormhole's chart through the throat, r_s = 1 and
    a = 5/4, zero at the throat: dl_*/dl = (a + |l|)/(a + |l| - r_s)."""
    return l + np.sign(l) * np.log(1 + 4 * np.abs(l))


def _two_particles_distance(x):
    """The proper distance from the midpoint along the line through two particles of alpha = 3/4 at
    x = +-1, the integral of |s^2 - 1|^(-1/4): x F(1/4, 1/2; 3/2; x^2) between them, and beyond
    either one the same at 1 plus the integral of (4/3)(2 + w^(4/3))^(-1/4) dw up to
    w = (|x| - 1)^(3/4), which is s = 1 + w^(4/3) and has a smooth integrand."""
    from scipy.integrate import quad
    from scipy.special import hyp2f1
    x = np.asarray(x, float)
    a = np.abs(x)
    between = np.minimum(a, 1.0)
    out = between * hyp2f1(0.25, 0.5, 1.5, between ** 2)
    beyond = np.vectorize(lambda w: quad(lambda v: 4 / 3 * (2 + v ** (4 / 3)) ** -0.25, 0, w, epsabs=1e-13, epsrel=1e-13)[0])
    return np.sign(x) * (out + beyond(np.maximum(a - 1.0, 0.0) ** 0.75))


CLOSED_FORMS = {
    # Bonnor's beam: a ray with the beam keeps ct - z, or u, and a null curve against it keeps
    # ct + z + A (ct - z), or v + A u, with A the declared profile's value on the plane.
    ("light_beam", "cartesian", "axis"): (lambda t, z: t + z, lambda t, z: t - z, None),
    ("light_beam", "cartesian", "beside"):
        (lambda t, z: t + z + (1 + 2 * math.log(3)) / 8 * (t - z), lambda t, z: t - z, None),
    ("light_beam", "null_cartesian", "midway"): (lambda u, v: v + (1 + 2 * math.log(2)) / 4 * u, lambda u, v: u, None),
    ("light_beam", "null_cartesian", "one"): (lambda u, v: v + (1 + 2 * math.log(4)) / 8 * u, lambda u, v: u, None),
    ("light_beam", "null_cylindrical_interior", "edge"): (lambda u, v: v + u / 8, lambda u, v: u, None),
    **{("light_beam", "null_cylindrical_exterior", view):
       (lambda u, v, n=n: v + (1 + 2 * math.log(n)) / 8 * u, lambda u, v: u, None)
       for view, n in (("twice", 2), ("four", 4))},
    ("btz", "stationary", "static"):
        (lambda t, r: t + _btz_rstar(r, BTZ_STATIC), lambda t, r: t - _btz_rstar(r, BTZ_STATIC), _away(1.0)),
    ("btz", "stationary", "rotating"):
        (lambda t, r: t + _btz_rstar(r, BTZ_ROTATING), lambda t, r: t - _btz_rstar(r, BTZ_ROTATING),
         _away(np.sqrt(0.8), np.sqrt(0.2))),
    ("btz", "eddington_finkelstein_ingoing", "static"):
        (lambda v, r: v, lambda v, r: v - 2 * _btz_rstar(r, BTZ_STATIC), _away(1.0)),
    ("btz", "eddington_finkelstein_ingoing", "rotating"):
        (lambda v, r: v, lambda v, r: v - 2 * _btz_rstar(r, BTZ_ROTATING), _away(np.sqrt(0.8), np.sqrt(0.2))),
    ("btz", "eddington_finkelstein_outgoing", "static"):
        (lambda u, r: u + 2 * _btz_rstar(r, BTZ_STATIC), lambda u, r: u, _away(1.0)),
    ("btz", "eddington_finkelstein_outgoing", "rotating"):
        (lambda u, r: u + 2 * _btz_rstar(r, BTZ_ROTATING), lambda u, r: u, _away(np.sqrt(0.8), np.sqrt(0.2))),
    ("schwarzschild", "spherical", "radial"):
        (lambda t, r: t + _rstar(r, [1]), lambda t, r: t - _rstar(r, [1]), lambda t, r: np.abs(r - 1) > 0.05),
    ("schwarzschild", "eddington_finkelstein_ingoing", "finkelstein"):
        (lambda v, r: v, lambda v, r: v - 2 * _rstar(r, [1]), lambda v, r: np.abs(r - 1) > 0.05),
    ("schwarzschild", "eddington_finkelstein_outgoing", "finkelstein"):
        (lambda u, r: u + 2 * _rstar(r, [1]), lambda u, r: u, lambda u, r: np.abs(r - 1) > 0.05),
    ("myers_perry", "boyer_lindquist", "transverse"):
        (lambda t, r: t + _myers_rstar(r), lambda t, r: t - _myers_rstar(r), _away(0.8)),
    ("myers_perry", "ingoing_kerr", "transverse"):
        (lambda v, r: v, lambda v, r: v - 2 * _myers_rstar(r), _away(0.8)),
    ("myers_perry", "boyer_lindquist_six", "transverse"):
        (lambda t, r: t + _myers_rstar(r, 6), lambda t, r: t - _myers_rstar(r, 6), _away(0.4131104345566232)),
    ("kaluza_klein_monopole", "gross_perry", "radial"):
        (lambda t, r: t + _kk_rstar(r), lambda t, r: t - _kk_rstar(r), lambda t, r: r > 0.05),
    ("kaluza_klein_monopole", "hopf", "radial"):
        (lambda t, r: t + _kk_rstar(r), lambda t, r: t - _kk_rstar(r), lambda t, r: r > 0.05),
    ("kaluza_klein_monopole", "taub_nut", "radial"):
        (lambda t, rho: t + _kk_rstar(rho - 2), lambda t, rho: t - _kk_rstar(rho - 2), lambda t, rho: rho > 2.05),
    **{("kaluza_klein_black_hole", system, "radial"):
       (lambda t, r: t + _kkbh_rstar(r), lambda t, r: t - _kkbh_rstar(r), _away(1.0))
       for system in ("electric", "magnetic", "einstein")},
    ("kaluza_klein_black_hole", "eddington_finkelstein_ingoing", "finkelstein"):
        (lambda v, r: v - math.sqrt(2) * _rstar(r, [1]) + _kkbh_rstar(r),
         lambda v, r: v - math.sqrt(2) * _rstar(r, [1]) - _kkbh_rstar(r), _away(1.0)),
    ("kaluza_klein_black_hole", "dyonic", "radial"):
        (lambda t, rho: t + _rstar(rho, [2, 1]), lambda t, rho: t - _rstar(rho, [2, 1]), _away(2.0, 1.0)),
    ("kaluza_klein_black_hole", "einstein_eddington_finkelstein", "finkelstein"):
        (lambda v, r: v, lambda v, r: v - 2 * _kkbh_rstar(r), _away(1.0)),
    ("tangherlini", "spherical", "radial"):
        (lambda t, r: t + _tangherlini_rstar(r), lambda t, r: t - _tangherlini_rstar(r), _away(1.0)),
    ("tangherlini", "eddington_finkelstein_ingoing", "finkelstein"):
        (lambda v, r: v, lambda v, r: v - 2 * _tangherlini_rstar(r), _away(1.0)),
    ("tangherlini", "eddington_finkelstein_outgoing", "finkelstein"):
        (lambda u, r: u + 2 * _tangherlini_rstar(r), lambda u, r: u, _away(1.0)),
    ("tangherlini", "spherical_six", "radial"):
        (lambda t, r: t + _tangherlini_rstar(r, 6), lambda t, r: t - _tangherlini_rstar(r, 6), _away(1.0)),
    ("boulware_deser", "spherical", "radial"):
        (lambda t, r: t + slices.boulware_deser_rstar(r), lambda t, r: t - slices.boulware_deser_rstar(r), _away(1.0)),
    ("boulware_deser", "eddington_finkelstein_ingoing", "finkelstein"):
        (lambda v, r: v, lambda v, r: v - 2 * slices.boulware_deser_rstar(r), _away(1.0)),
    ("boulware_deser", "eddington_finkelstein_outgoing", "finkelstein"):
        (lambda u, r: u + 2 * slices.boulware_deser_rstar(r), lambda u, r: u, _away(1.0)),
    ("boulware_deser", "spherical_plus", "radial"):
        (lambda t, r: t + slices.boulware_deser_plus_rstar(r), lambda t, r: t - slices.boulware_deser_plus_rstar(r),
         lambda t, r: r > 0.02),
    ("black_string", "static", "radial"):
        (lambda t, r: t + _rstar(r, [1]), lambda t, r: t - _rstar(r, [1]), _away(1.0)),
    ("black_string", "eddington_finkelstein_ingoing", "finkelstein"):
        (lambda v, r: v, lambda v, r: v - 2 * _rstar(r, [1]), _away(1.0)),
    ("black_string", "kerr_schild", "radial"):
        (lambda T, r: T + r, lambda T, r: T + r - 2 * _rstar(r, [1]), _away(1.0)),
    ("black_string", "static_six", "radial"):
        (lambda t, r: t + _tangherlini_rstar(r), lambda t, r: t - _tangherlini_rstar(r), _away(1.0)),
    ("global_monopole", "static", "radial"):
        (lambda t, r: t + _rstar(r, [GM_RH]) / 0.81, lambda t, r: t - _rstar(r, [GM_RH]) / 0.81,
         lambda t, r: np.abs(r - GM_RH) > 0.05),
    ("global_monopole", "conical", "radial"): (lambda t, r: t + r, lambda t, r: t - r, None),
    ("global_monopole", "eddington_finkelstein_ingoing", "finkelstein"):
        (lambda v, r: v, lambda v, r: v - 2 * _rstar(r, [GM_RH]) / 0.81, lambda v, r: np.abs(r - GM_RH) > 0.05),
    ("global_monopole", "eddington_finkelstein_outgoing", "finkelstein"):
        (lambda u, r: u + 2 * _rstar(r, [GM_RH]) / 0.81, lambda u, r: u, lambda u, r: np.abs(r - GM_RH) > 0.05),
    ("string_black_hole", "static", "radial"):
        (lambda t, r: t + _rstar(r, [1]), lambda t, r: t - _rstar(r, [1]), lambda t, r: np.abs(r - 1) > 0.05),
    ("string_black_hole", "wedge", "radial"):
        (lambda t, r: t + _rstar(r, [1]), lambda t, r: t - _rstar(r, [1]), lambda t, r: np.abs(r - 1) > 0.05),
    ("string_black_hole", "eddington_finkelstein_ingoing", "finkelstein"):
        (lambda v, r: v, lambda v, r: v - 2 * _rstar(r, [1]), lambda v, r: np.abs(r - 1) > 0.05),
    ("string_black_hole", "eddington_finkelstein_outgoing", "finkelstein"):
        (lambda u, r: u + 2 * _rstar(r, [1]), lambda u, r: u, lambda u, r: np.abs(r - 1) > 0.05),
    **{("dilaton_black_hole", system, "radial"):
       (lambda t, r: t + _rstar(r, [1]), lambda t, r: t - _rstar(r, [1]), lambda t, r: np.abs(r - 1) > 0.05)
       for system in ("static", "string_magnetic", "string_electric")},
    ("dilaton_black_hole", "eddington_finkelstein_ingoing", "finkelstein"):
        (lambda v, r: v, lambda v, r: v - 2 * _rstar(r, [1]), lambda v, r: np.abs(r - 1) > 0.05),
    ("dilaton_black_hole", "eddington_finkelstein_outgoing", "finkelstein"):
        (lambda u, r: u + 2 * _rstar(r, [1]), lambda u, r: u, lambda u, r: np.abs(r - 1) > 0.05),
    ("schwarzschild_ads", "static", "radial"):
        (lambda t, r: t + _sads_rstar(r), lambda t, r: t - _sads_rstar(r), lambda t, r: np.abs(r - 1) > 0.05),
    ("schwarzschild_ads", "eddington_finkelstein_ingoing", "finkelstein"):
        (lambda v, r: v, lambda v, r: v - 2 * _sads_rstar(r), lambda v, r: np.abs(r - 1) > 0.05),
    ("schwarzschild_ads", "eddington_finkelstein_outgoing", "finkelstein"):
        (lambda u, r: u + 2 * _sads_rstar(r), lambda u, r: u, lambda u, r: np.abs(r - 1) > 0.05),
    ("reissner_nordstrom_de_sitter", "static", "radial"):
        (lambda t, r: t + _rnds_rstar(r), lambda t, r: t - _rnds_rstar(r), _rnds_away),
    ("reissner_nordstrom_de_sitter", "eddington_finkelstein_ingoing", "finkelstein"):
        (lambda v, r: v, lambda v, r: v - 2 * _rnds_rstar(r), _rnds_away),
    ("reissner_nordstrom_de_sitter", "eddington_finkelstein_outgoing", "finkelstein"):
        (lambda u, r: u + 2 * _rnds_rstar(r), lambda u, r: u, _rnds_away),
    ("reissner_nordstrom_de_sitter", "cosmological", "plane"):
        (lambda tau, rho: _rnds_cosmic(tau, rho, 1), lambda tau, rho: _rnds_cosmic(tau, rho, -1), _rnds_cosmic_away),
    ("bardeen", "static", "radial"):
        (lambda t, r: t + _bardeen_rstar(r), lambda t, r: t - _bardeen_rstar(r), _bardeen_away),
    ("bardeen", "eddington_finkelstein_ingoing", "finkelstein"):
        (lambda v, r: v, lambda v, r: v - 2 * _bardeen_rstar(r), _bardeen_away),
    ("bardeen", "eddington_finkelstein_outgoing", "finkelstein"):
        (lambda u, r: u + 2 * _bardeen_rstar(r), lambda u, r: u, _bardeen_away),
    ("hayward", "static", "radial"):
        (lambda t, r: t + _hayward_rstar(r), lambda t, r: t - _hayward_rstar(r), _hayward_away),
    ("hayward", "eddington_finkelstein_ingoing", "finkelstein"):
        (lambda v, r: v, lambda v, r: v - 2 * _hayward_rstar(r), _hayward_away),
    ("hayward", "eddington_finkelstein_outgoing", "finkelstein"):
        (lambda u, r: u + 2 * _hayward_rstar(r), lambda u, r: u, _hayward_away),
    **_tbh_closed(),
    ("reissner_nordstrom_ads", "static", "radial"):
        (lambda t, r: t + slices.rnads_rstar(r), lambda t, r: t - slices.rnads_rstar(r), _rnads_away),
    ("reissner_nordstrom_ads", "eddington_finkelstein_ingoing", "finkelstein"):
        (lambda v, r: v, lambda v, r: v - 2 * slices.rnads_rstar(r), _rnads_away),
    ("reissner_nordstrom_ads", "eddington_finkelstein_outgoing", "finkelstein"):
        (lambda u, r: u + 2 * slices.rnads_rstar(r), lambda u, r: u, _rnads_away),
    ("schwarzschild_de_sitter", "static", "radial"):
        (lambda t, r: t + _sds_rstar(r), lambda t, r: t - _sds_rstar(r), _sds_away),
    ("schwarzschild_de_sitter", "eddington_finkelstein_ingoing", "finkelstein"):
        (lambda v, r: v, lambda v, r: v - 2 * _sds_rstar(r), _sds_away),
    ("schwarzschild_de_sitter", "eddington_finkelstein_outgoing", "finkelstein"):
        (lambda u, r: u + 2 * _sds_rstar(r), lambda u, r: u, _sds_away),
    ("kiselev", "static", "radial"):
        (lambda t, r: t + _kiselev_rstar(r), lambda t, r: t - _kiselev_rstar(r), _kiselev_away),
    ("kiselev", "linear", "radial"):
        (lambda t, r: t + _kiselev_rstar(r), lambda t, r: t - _kiselev_rstar(r), _kiselev_away),
    ("kiselev", "eddington_finkelstein_ingoing", "finkelstein"):
        (lambda v, r: v, lambda v, r: v - 2 * _kiselev_rstar(r), _kiselev_away),
    ("kiselev", "eddington_finkelstein_outgoing", "finkelstein"):
        (lambda u, r: u + 2 * _kiselev_rstar(r), lambda u, r: u, _kiselev_away),
    # Kiselev's two charts of the matter alone are conformal to flat planes, so the rays keep the sums.
    ("kiselev", "hyperbolic", "radial"): (lambda eta, chi: eta + chi, lambda eta, chi: eta - chi, None),
    ("kiselev", "conformally_flat", "radial"): (lambda tau, rho: tau + rho, lambda tau, rho: tau - rho, None),
    ("rn_metric", "spherical", "radial"):
        (lambda t, r: t + _rstar(r, [0.64, 0.36]), lambda t, r: t - _rstar(r, [0.64, 0.36]),
         lambda t, r: (np.abs(r - 0.64) > 0.05) & (np.abs(r - 0.36) > 0.05)),
    ("kastor_traschen", "cartesian", "tz"): (_kt_label("axis", -1), _kt_label("axis", 1), _kt_keep),
    ("kastor_traschen", "cartesian", "tx"): (_kt_label("midplane", -1), _kt_label("midplane", 1), _kt_keep),
    ("kastor_traschen", "cylindrical", "radial"): (_kt_label("midplane", -1), _kt_label("midplane", 1), _kt_keep),
    ("kastor_traschen", "isotropic", "radial"):
        (lambda tau, r: _kt_cosmic(tau, r, 1), lambda tau, r: _kt_cosmic(tau, r, -1), _kt_cosmic_away),
    ("kastor_traschen", "comoving", "tz"):
        (_kt_label("axis", -1, True), _kt_label("axis", 1, True), _kt_every),
    ("kastor_traschen", "comoving", "tx"):
        (_kt_label("midplane", -1, True), _kt_label("midplane", 1, True), _kt_every),
    ("israel_wilson_perjes", "cylindrical", "midplane"):
        (lambda t, r: t + _iwp_midplane(r), lambda t, r: t - _iwp_midplane(r), None),
    ("israel_wilson_perjes", "spheroidal", "axis"):
        (lambda t, r: t + _iwp_rstar(r), lambda t, r: t - _iwp_rstar(r), None),
    ("israel_wilson_perjes", "spheroidal", "principal"):
        (lambda t, r: t + _iwp_rstar(r), lambda t, r: t - _iwp_rstar(r), None),
    ("israel_wilson_perjes", "spherical", "radial"):
        (lambda t, r: t + r + 2 * np.log(r - 1) - 1.25 / (r - 1), lambda t, r: t - r - 2 * np.log(r - 1) + 1.25 / (r - 1),
         lambda t, r: r > 1.05),
    ("majumdar_papapetrou", "cartesian", "tz"):
        (lambda t, z: t + _mp_axis(z), lambda t, z: t - _mp_axis(z), lambda t, z: np.abs(np.abs(z) - 2) > 0.05),
    ("majumdar_papapetrou", "cartesian", "tx"):
        (lambda t, x: t + _mp_midplane(x), lambda t, x: t - _mp_midplane(x), None),
    ("majumdar_papapetrou", "cylindrical", "radial"):
        (lambda t, r: t + _mp_midplane(r), lambda t, r: t - _mp_midplane(r), None),
    ("majumdar_papapetrou", "isotropic", "radial"):
        (lambda t, r: t + r + 2 * np.log(r) - 1 / r, lambda t, r: t - r - 2 * np.log(r) + 1 / r,
         lambda t, r: r > 0.05),
    ("de_sitter", "static_spherical", "radial"):
        (lambda t, r: t + 0.5 * np.log(np.abs((1 + r) / (1 - r))),
         lambda t, r: t - 0.5 * np.log(np.abs((1 + r) / (1 - r))), lambda t, r: np.abs(r - 1) > 0.05),
    ("einstein_static", "hyperspherical", "radial"): (lambda t, c: t + c, lambda t, c: t - c, None),
    ("einstein_static", "hyperspherical", "through"): (lambda t, c: t + c, lambda t, c: t - c, None),
    ("einstein_static", "static_areal", "radial"):
        (lambda t, r: t + np.arcsin(r), lambda t, r: t - np.arcsin(r), lambda t, r: r < 0.999),
    ("milne", "comoving_hyperbolic", "through"):
        (lambda t, c: np.log(t) + c, lambda t, c: np.log(t) - c, lambda t, c: t > 0.02),
    ("milne", "comoving_spherical", "radial"):
        (lambda t, r: np.log(t) + np.arcsinh(r), lambda t, r: np.log(t) - np.arcsinh(r), lambda t, r: t > 0.02),
    ("milne", "logarithmic_time", "radial"): (lambda tau, c: tau + c, lambda tau, c: tau - c, None),
    ("milne", "inertial", "through"): (lambda T, R: T + R, lambda T, R: T - R, None),
    ("anti_de_sitter", "static_global", "radial"):
        (lambda t, r: t + np.arctan(r), lambda t, r: t - np.arctan(r), None),
    # Randall and Sundrum's wall at k = 1: the rays keep ct -+ sgn(y)(e^|y| - 1), which is ct -+ w, and
    # ct -+ z on one side; between two walls at r_c = 1/2, ct -+ sgn(phi)(e^(|phi|/2) - 1).
    ("randall_sundrum", "proper_distance", "ty"):
        (lambda t, y: t + np.sign(y) * np.expm1(np.abs(y)), lambda t, y: t - np.sign(y) * np.expm1(np.abs(y)), None),
    ("randall_sundrum", "conformal", "tw"): (lambda t, w: t + w, lambda t, w: t - w, None),
    ("randall_sundrum", "poincare", "tz"): (lambda t, z: t + z, lambda t, z: t - z, None),
    ("randall_sundrum", "two_walls", "tphi"):
        (lambda t, f: t + np.sign(f) * np.expm1(np.abs(f) / 2), lambda t, f: t - np.sign(f) * np.expm1(np.abs(f) / 2), None),
    ("coleman_de_luccia", "wall", "into_flat"):
        (lambda psi, xi: psi + cdl_conformal_distance(xi, "into_flat"),
         lambda psi, xi: psi - cdl_conformal_distance(xi, "into_flat"), lambda psi, xi: (xi > 0.05) & (xi < 2.95)),
    ("coleman_de_luccia", "wall", "out_of_flat"):
        (lambda psi, xi: psi + cdl_conformal_distance(xi, "out_of_flat"),
         lambda psi, xi: psi - cdl_conformal_distance(xi, "out_of_flat"), lambda psi, xi: xi > 0.05),
    ("coleman_de_luccia", "open", "zero"):
        (lambda tau, c: np.log(tau) + c, lambda tau, c: np.log(tau) - c, lambda tau, c: tau > 0.02),
    ("coleman_de_luccia", "open", "positive"):
        (lambda tau, c: np.log(np.tanh(tau / 2)) + c, lambda tau, c: np.log(np.tanh(tau / 2)) - c,
         lambda tau, c: tau > 0.02),
    ("coleman_de_luccia", "open", "negative"):
        (lambda tau, c: np.log(np.tan(tau / 2)) + c, lambda tau, c: np.log(np.tan(tau / 2)) - c,
         lambda tau, c: (tau > 0.02) & (tau < math.pi - 0.02)),
    ("coleman_de_luccia", "open_conformal", "negative"): (lambda eta, c: eta + c, lambda eta, c: eta - c, None),
    ("coleman_de_luccia", "static_inside", "into_flat"): (lambda t, r: t + r, lambda t, r: t - r, None),
    ("coleman_de_luccia", "static_inside", "out_of_flat"):
        (lambda t, r: t + np.arctan(r), lambda t, r: t - np.arctan(r), None),
    ("coleman_de_luccia", "static_outside", "into_flat"):
        (lambda t, r: t + 0.5 * np.log(np.abs((1 + r) / (1 - r))),
         lambda t, r: t - 0.5 * np.log(np.abs((1 + r) / (1 - r))), lambda t, r: np.abs(r - 1) > 0.05),
    ("coleman_de_luccia", "static_outside", "out_of_flat"): (lambda t, r: t + r, lambda t, r: t - r, None),
    ("domain_wall", "planar", "tz"):
        (lambda t, z: t - np.sign(z) * np.log(1 - np.abs(z)), lambda t, z: t + np.sign(z) * np.log(1 - np.abs(z)),
         lambda t, z: np.abs(z) < 0.95),
    ("domain_wall", "inertial", "through"): (lambda T, R: T + R, lambda T, R: T - R, None),
    ("ellis_bronnikov", "spherical", "radial"): (lambda t, r: t + r, lambda t, r: t - r, None),
    ("wormhole_time_machine", "short_throat", "radial"): (lambda t, l: t + l, lambda t, l: t - l, None),
    ("morris_thorne", "spherical", "radial"):
        (lambda t, r: t + np.sqrt(r ** 2 - 1), lambda t, r: t - np.sqrt(r ** 2 - 1), lambda t, r: r > 1.0005),
    # Teo's example at b_0 = 1: on the equator the rays keep ct -+ l, and on the axis at a = 1/4,
    # where N = 1 + 1/r, they keep ct -+ r_* with dr_*/dr = 1/(N sqrt(1 - 1/r)).
    ("teo_wormhole", "spherical", "axis"):
        (lambda t, r: t + _teo_axis(r), lambda t, r: t - _teo_axis(r), lambda t, r: r > 1.0005),
    ("teo_wormhole", "spherical", "equator"):
        (lambda t, r: t + _teo_proper(r), lambda t, r: t - _teo_proper(r), lambda t, r: r > 1.0005),
    ("teo_wormhole", "proper_radial", "axis"):
        (lambda t, l: t + np.sign(l) * _teo_axis(teo_inverse(l)[0]),
         lambda t, l: t - np.sign(l) * _teo_axis(teo_inverse(l)[0]), None),
    ("teo_wormhole", "proper_radial", "equator"): (lambda t, l: t + l, lambda t, l: t - l, None),
    # With r_s = 1 and a = 5/4, l_* = l + ln(1 + 4|l|) sgn l on the throat's chart and r_* = r + ln(r - 1).
    ("thin_shell_wormhole", "throat", "radial"):
        (lambda t, l: t + _tsw_lstar(l), lambda t, l: t - _tsw_lstar(l), None),
    ("thin_shell_wormhole", "spherical", "radial"):
        (lambda t, r: t + _rstar(r, [1]), lambda t, r: t - _rstar(r, [1]), None),
    # The gravastar at r_s = 1, R = 5/4 and L = 2: ct -+ x inside, x = (L/sqrt(C)) artanh(r/L), and
    # Schwarzschild's ct -+ r_* outside.
    ("gravastar", "interior", "radial"):
        (lambda t, r: t + _gravastar_x(r), lambda t, r: t - _gravastar_x(r), None),
    ("gravastar", "interior", "through"):
        (lambda t, r: t + _gravastar_x(r), lambda t, r: t - _gravastar_x(r), None),
    ("gravastar", "interior_tortoise", "radial"): (lambda t, x: t + x, lambda t, x: t - x, None),
    ("gravastar", "exterior", "radial"):
        (lambda t, r: t + _rstar(r, [1]), lambda t, r: t - _rstar(r, [1]), None),
    # Mass inflation: an ingoing ray keeps its v, and an outgoing one the radius at which it crosses
    # one ingoing ray, carried there by ori_shell.py's own integrators, which the tracing does not use.
    # Inside the inner apparent horizon a ray comes from r = 0 and may never have crossed that
    # ingoing ray, so the check is made outside it.
    ("mass_inflation", "ingoing", "tail"):
        (lambda v, r: v, lambda v, r: ori_shell.label_before(v, r, 12.0),
         lambda v, r: r > ori_shell.horizons(ori_shell.mass_before(v))[0] + 0.03),
    ("mass_inflation", "ingoing", "behind"):
        # in tenths of m_0: behind the shell the rays span 0.2 m_0 of radius, and the check tells the two
        # families apart by a spread above one.
        (lambda v, r: v, lambda v, r: 10 * ori_shell.label_behind(np.minimum(v, -1e-3), r, -1.0),
         lambda v, r: (v < -1e-3) & (r < ori_shell.shell().radius_behind(np.minimum(v, -1e-3)) - 0.01)
         & (r > ori_shell.horizons(ori_shell.shell().mass_behind(np.minimum(v, -1e-3)))[0] + 0.03)),
    ("gowdy", "areal", "plane"): (lambda t, th: t + th, lambda t, th: t - th, lambda t, th: t > 0.02),
    ("gowdy", "sphere", "plane"): (lambda t, th: t + th, lambda t, th: t - th, None),
    ("gowdy", "logarithmic", "plane"): (lambda tau, th: np.exp(-tau) + th, lambda tau, th: np.exp(-tau) - th, None),
    # With r_s = 1 and a = 1/2, 1 and 2: ct -+ r_* in the two charts of t, and v - 2r_* and u + 2r_*.
    **{("simpson_visser", chart, case): _sv_forms(float(Fraction(a)), chart)
       for chart in ("spherical", "areal", "eddington_finkelstein_ingoing", "eddington_finkelstein_outgoing")
       for case, (_, a) in SV_CASES.items()},
    # With r_s = 1 and lambda = 1/5, and r = 1 + u^2 through the throat, where r_* changes sign with u.
    ("damour_solodukhin", "spherical", "radial"):
        (lambda t, r: t + _ds_rstar(r), lambda t, r: t - _ds_rstar(r), lambda t, r: r > 1.0005),
    # The rescaled time is sqrt(1 + lambda^2) t, and in rho the tortoise coordinate of that time is
    # r_s((2 + lambda^2) rho + lambda^2 sinh rho)/(2(1 + lambda^2)).
    ("damour_solodukhin", "rescaled", "radial"):
        (lambda t, r: t + math.sqrt(1.04) * _ds_rstar(r), lambda t, r: t - math.sqrt(1.04) * _ds_rstar(r),
         lambda t, r: r > 1.0005),
    ("damour_solodukhin", "throat", "radial"):
        (lambda t, rho: t + _ds_xstar(rho), lambda t, rho: t - _ds_xstar(rho), None),
    # The isotropic radius has the areal radius r(1 + 1/4r)^2, and r_* changes sign at the throat r = 1/4.
    ("damour_solodukhin", "isotropic", "radial"):
        (lambda t, r: t + np.sign(4 * r - 1) * _ds_rstar(_ds_areal(r)),
         lambda t, r: t - np.sign(4 * r - 1) * _ds_rstar(_ds_areal(r)), lambda t, r: r > 0.005),
    ("damour_solodukhin", "einstein_rosen", "radial"):
        (lambda t, u: t + np.sign(u) * _ds_rstar(1 + u ** 2), lambda t, u: t - np.sign(u) * _ds_rstar(1 + u ** 2), None),
    # With b = 1 and gamma = 1/2: r = R + 3/4, r = rho (1 + 1/(4 rho))^2, and in the harmonic chart, at
    # k = 1 and b = 2, r = 2/(1 - e^(-2u)) falls as u grows and the tortoise coordinate is twice that of b = 1.
    ("fisher_jnw", "spherical", "radial"):
        (lambda t, r: t + _fjnw_rstar(r), lambda t, r: t - _fjnw_rstar(r), lambda t, r: r > 1.0005),
    ("fisher_jnw", "jnw", "radial"):
        (lambda t, R: t + _fjnw_rstar(R + 0.75), lambda t, R: t - _fjnw_rstar(R + 0.75), lambda t, R: R > 0.2505),
    ("fisher_jnw", "isotropic", "radial"):
        (lambda t, rho: t + _fjnw_rstar(rho * (1 + 1 / (4 * rho)) ** 2),
         lambda t, rho: t - _fjnw_rstar(rho * (1 + 1 / (4 * rho)) ** 2), lambda t, rho: rho > 0.2505),
    ("fisher_jnw", "harmonic", "radial"):
        (lambda t, u: t - 2 * _fjnw_rstar(1 / (1 - np.exp(-2 * u))), lambda t, u: t + 2 * _fjnw_rstar(1 / (1 - np.exp(-2 * u))),
         lambda t, u: u > 0.05),
    # Witten's black hole at lambda = m = 1: the tortoise coordinate is ln sinh r in his own chart,
    # ln|e^(2x) - 1|/2 in the charts of x and ln|w - 1|/2 in the dilaton chart, and it is sigma itself.
    ("witten_black_hole", "witten", "radial"):
        (lambda t, r: t + np.log(np.sinh(r)), lambda t, r: t - np.log(np.sinh(r)), lambda t, r: r > 0.02),
    ("witten_black_hole", "schwarzschild_gauge", "radial"):
        (lambda t, x: t + _witten_xstar(x), lambda t, x: t - _witten_xstar(x), lambda t, x: np.abs(x) > 0.02),
    ("witten_black_hole", "dilaton", "radial"):
        (lambda t, w: t + np.log(np.abs(w - 1)) / 2, lambda t, w: t - np.log(np.abs(w - 1)) / 2,
         lambda t, w: np.abs(w - 1) > 0.02),
    ("witten_black_hole", "conformal", "radial"): (lambda t, s: t + s, lambda t, s: t - s, None),
    ("witten_black_hole", "kruskal", "plane"): (lambda U, V: V, lambda U, V: U, None),
    ("witten_black_hole", "eddington_finkelstein_ingoing", "finkelstein"):
        (lambda v, x: v, lambda v, x: v - 2 * _witten_xstar(x), lambda v, x: np.abs(x) > 0.02),
    ("witten_black_hole", "eddington_finkelstein_outgoing", "finkelstein"):
        (lambda u, x: u + 2 * _witten_xstar(x), lambda u, x: u, lambda u, x: np.abs(x) > 0.02),
    # Roberts's collapse: the advanced and the retarded null coordinate v and u, in each chart's own coordinates.
    **{("roberts", "double_null", case): (lambda u, v: v, lambda u, v: u, None) for case in ROBERTS_CASES},
    # Hiscock's evaporating hole: an ingoing ray keeps v and an outgoing one u, whatever the mass
    # does; the other family of each Vaidya chart has no closed form for the declared mass.
    ("hiscock", "ingoing", "history"): (lambda v, r: v, None, None),
    ("hiscock", "ingoing", "shells"): (lambda v, r: v, None, None),
    ("hiscock", "outgoing", "history"): (None, lambda u, r: u, None),
    ("hiscock", "flat", "after"): (lambda u, v: v, lambda u, v: u, None),
    **{("roberts", "advanced", case): (lambda v, r: v, lambda v, r, p=float(Fraction(p)): (1 + p) * v - 2 * r, None)
       for case, (_, p) in ROBERTS_CASES.items()},
    **{("roberts", "areal", case):
       (lambda v, R: v, lambda v, R, p=float(Fraction(p)): v - 2 * np.sqrt(R ** 2 + p ** 2 * v ** 2 / 4),
        lambda v, R: R > 0.05)
       for case, (_, p) in ROBERTS_CASES.items()},
    **{("roberts", "diagonal", case): (lambda t, rho: t + rho, lambda t, rho: t - rho, None) for case in ROBERTS_CASES},
    **{("roberts", "scaling", case): (lambda tau, x: np.exp(-tau) * (np.exp(2 * x) - 1), lambda tau, x: tau,
                                      lambda tau, x: x > 0.01) for case in ROBERTS_CASES},
    # The neutral bridge at r_s = 1, r_* = r + ln(r - 1), which changes sign with the sheet: r = 1 + u^2 in
    # Einstein and Rosen's u and r(1 + 1/4r)^2 in the isotropic radius. The charged bridge with a mass has
    # the roots 3/2 and -1/2, and with no mass, at r_q = 1, r = sqrt(u^2 + 1) and the roots 1 and -1.
    ("einstein_rosen_bridge", "bridge", "radial"):
        (lambda t, u: t + np.sign(u) * _rstar(1 + u ** 2, [1]), lambda t, u: t - np.sign(u) * _rstar(1 + u ** 2, [1]),
         lambda t, u: np.abs(u) > 0.2),
    ("einstein_rosen_bridge", "spherical", "radial"):
        (lambda t, r: t + _rstar(r, [1]), lambda t, r: t - _rstar(r, [1]), lambda t, r: r > 1.05),
    ("einstein_rosen_bridge", "isotropic", "radial"):
        (lambda t, r: t + np.sign(4 * r - 1) * _rstar(r * (1 + 1 / (4 * r)) ** 2, [1]),
         lambda t, r: t - np.sign(4 * r - 1) * _rstar(r * (1 + 1 / (4 * r)) ** 2, [1]),
         lambda t, r: (np.abs(4 * r - 1) > 0.2) & (r > 0.02)),
    ("einstein_rosen_bridge", "charged_spherical", "radial"):
        (lambda t, r: t + _rstar(r, [1.5, -0.5]), lambda t, r: t - _rstar(r, [1.5, -0.5]), lambda t, r: r > 1.55),
    ("einstein_rosen_bridge", "charged_bridge", "radial"):
        (lambda t, u: t + np.sign(u) * _rstar(np.sqrt(1 + u ** 2), [1, -1]),
         lambda t, u: t - np.sign(u) * _rstar(np.sqrt(1 + u ** 2), [1, -1]), lambda t, u: np.abs(u) > 0.2),
    ("einstein_rosen_waves", "cylindrical", "radial"): (lambda t, r: t + r, lambda t, r: t - r, None),
    ("einstein_rosen_waves", "null", "radial"): (lambda u, v: v, lambda u, v: u, None),
    # With m = 1, z_* = z e^(2/z) - 2 Ei(2/z) on the axis and rho_* the integral of e^(2/s - 1/(2 s^2)) in the plane.
    **{("curzon_chazy", system, "axis"): (lambda t, z: t + _curzon_axis(z), lambda t, z: t - _curzon_axis(z),
                                          lambda t, z: z > 0.35) for system in ("weyl", "spherical")},
    **{("curzon_chazy", system, "equator"): (lambda t, r: t + _curzon_plane(r), lambda t, r: t - _curzon_plane(r),
                                             lambda t, r: r > 0.02) for system in ("weyl", "spherical")},
    **{("zipoy_voorhees", system, f"{plane}_{shape}"): _zv_forms(star, shape == "oblate", shift, edge)
       for system, shift in (("spherical", 0), ("prolate_spheroidal", 1))
       for shape in ("oblate", "prolate")
       for plane, star, edge in (("axis", _zv_axis, 0.1), ("equator", _zv_equator, 0.02))},
    ("melvin", "cylindrical", "radial"): (lambda t, r: t + r, lambda t, r: t - r, None),
    ("senovilla", "cylindrical", "radial"): (lambda t, r: t + r, lambda t, r: t - r, None),
    ("melvin", "ernst", "radial"):
        (lambda t, r: t + _rstar(r, [1]), lambda t, r: t - _rstar(r, [1]), lambda t, r: np.abs(r - 1) > 0.05),
    ("levi_civita", "weyl", "radial"):
        (lambda t, r: t + 4 * r ** 0.25, lambda t, r: t - 4 * r ** 0.25, lambda t, r: r > 0.01),
    ("levi_civita", "kasner", "radial"):
        (lambda t, r: t + 3 * np.cbrt(r), lambda t, r: t - 3 * np.cbrt(r), lambda t, r: r > 0.01),
    ("point_particle_2plus1", "conical", "radial"): (lambda t, r: t + r, lambda t, r: t - r, None),
    ("point_particle_2plus1", "wedge", "radial"): (lambda t, r: t + r, lambda t, r: t - r, None),
    ("point_particle_2plus1", "circumference", "radial"): (lambda t, R: t + R / 0.75, lambda t, R: t - R / 0.75, None),
    ("point_particle_2plus1", "planet", "radial"): (lambda t, chi: t + chi, lambda t, chi: t - chi, None),
    ("point_particle_2plus1", "isotropic", "radial"):
        (lambda t, r: t + r ** 0.75 / 0.75, lambda t, r: t - r ** 0.75 / 0.75, lambda t, r: r > 0.01),
    **{("point_particle_2plus1", "two_bodies", view):
        (lambda t, x: t + _two_particles_distance(x), lambda t, x: t - _two_particles_distance(x),
         lambda t, x: np.abs(np.abs(x) - 1) > 0.01) for view in ("between", "beyond")},
    ("minkowski", "spherical", "radial"): (lambda t, r: t + r, lambda t, r: t - r, None),
    ("minkowski", "cartesian", "tx"): (lambda t, x: t + x, lambda t, x: t - x, None),
    ("minkowski", "rindler", "tx"):
        (lambda T, X: T + np.log(X), lambda T, X: T - np.log(X), lambda T, X: X > 1e-3),
    ("frw", "conformal_spherical", "radial"): (lambda e, r: e + r, lambda e, r: e - r, None),
    ("de_sitter", "flat_slicing", "tx"): (lambda t, x: x - np.exp(-t), lambda t, x: x + np.exp(-t), None),
    ("nariai", "static", "patch"):
        (lambda t, r: t + np.arctanh(r), lambda t, r: t - np.arctanh(r), lambda t, r: np.abs(r) < 0.95),
    ("nariai", "global", "circle"):
        (lambda t, c: c + np.arctan(np.sinh(t)), lambda t, c: c - np.arctan(np.sinh(t)), None),
    ("anti_de_sitter", "poincare", "tx"): (lambda t, x: t + x, lambda t, x: t - x, None),
    ("bertotti_robinson", "poincare", "tx"): (lambda t, x: t + x, lambda t, x: t - x, None),
    # The throat of extreme Kerr at r_0 = 1: each plane is conformal to a chart of AdS2, whose rays
    # keep ct -+ r_* with r_* = -1/r, x, arctan y, and ln(1 - 1/r) at k = 1/2.
    ("near_horizon_extreme_kerr", "poincare", "equator"): (lambda t, r: t - 1 / r, lambda t, r: t + 1 / r, None),
    ("near_horizon_extreme_kerr", "inverse_radius", "equator"): (lambda t, x: t + x, lambda t, x: t - x, None),
    ("near_horizon_extreme_kerr", "global", "equator"):
        (lambda t, y: t + np.arctan(y), lambda t, y: t - np.arctan(y), None),
    ("near_horizon_extreme_kerr", "near_nhek", "equator"):
        (lambda t, r: t + np.log(1 - 1 / r), lambda t, r: t - np.log(1 - 1 / r), lambda t, r: r > 1.02),
    ("godel", "cartesian", "tx"): (lambda t, x: t + x, lambda t, x: t - x, None),
    ("som_raychaudhuri", "cartesian", "tx"): (lambda t, x: t + x, lambda t, x: t - x, None),
    # Inside Schwarzschild's horizon dr/dT = +-T/(1 - T), so r +- (T + ln(1 - T)) is constant.
    ("kantowski_sachs", "schwarzschild_interior", "Tr"):
        (lambda T, r: r - T - np.log(1 - T), lambda T, r: r + T + np.log(1 - T), lambda T, r: T < 0.95),
    ("kantowski_sachs", "dust", "etar"):
        (lambda e, r: r + _ks_tau(e), lambda e, r: r - _ks_tau(e), lambda e, r: np.abs(e) < 1.5),
    ("kasner", "cartesian", "tx"):
        (lambda t, x: x + t ** (9 / 7) * 7 / 9, lambda t, x: x - t ** (9 / 7) * 7 / 9, lambda t, x: t > 1e-3),
    ("kasner", "cartesian", "tz"):
        (lambda t, z: z + 7 * t ** (1 / 7), lambda t, z: z - 7 * t ** (1 / 7), lambda t, z: t > 1e-3),
    ("pp_wave", "exact_plane_wave", "tz"): (lambda u, v: v, lambda u, v: u, None),
    # Across the declared pulse a ray moving left gains ln(rho_0/rho) times the pulse's integral.
    **{("aichelburg_sexl", "null_cartesian", view):
       (lambda u, v, n=float(sp.Rational(rho)): v + math.log(n) * (1 + scipy_erf(20 * u)) / 2, lambda u, v: u, None)
       for view, rho in AS_RHO.items()},
    ("khan_penrose", "double_null", "plane"): (lambda u, v: v, lambda u, v: u, None),
    ("khan_penrose", "cosmological", "plane"): (lambda tau, s: tau + s, lambda tau, s: tau - s, None),
    ("krasnikov", "cylindrical", "tx"): (None, lambda t, x: t - x, None),
    ("bell_szekeres", "double_null", "plane"): (lambda u, v: v, lambda u, v: u, None),
    ("bell_szekeres", "time_space", "plane"): (lambda xi, eta: xi + eta, lambda xi, eta: xi - eta, None),
    # cos(xi) e^(-+ky) = -T -+ Z, so (-T -+ Z)/(1 + sin(xi)) with sin(xi) = sqrt(1 - T^2 + Z^2) are null.
    ("bell_szekeres", "regular", "plane"):
        (lambda T, Z: -(T + Z) / (1 + np.sqrt(1 - T ** 2 + Z ** 2)), lambda T, Z: (Z - T) / (1 + np.sqrt(1 - T ** 2 + Z ** 2)),
         lambda T, Z: 1 - T ** 2 + Z ** 2 > 0.01),
    ("bell_szekeres", "global", "plane"):
        (lambda chi, rho: chi + np.arctan(np.sinh(rho)), lambda chi, rho: chi - np.arctan(np.sinh(rho)), None),
    ("bell_szekeres", "kruskal_szekeres", "plane"): (lambda U, V: V, lambda U, V: U, None),
    ("bell_szekeres", "bertotti_robinson", "plane"): (lambda t, r: t + r, lambda t, r: t - r, None),
    ("van_den_broeck", "pocket", "radial"):
        (lambda t, r: t + _vdb_distance(r), lambda t, r: t - _vdb_distance(r), None),
    ("van_den_broeck", "proper_radial", "radial"): (lambda t, l: t + l, lambda t, l: t - l, None),
}


def _c_metric_rstar(y):
    """alpha r* of the C-metric at alpha m = 1/6, in y = 1/(alpha r), up to a constant: Griffiths,
    Krtous and Podolsky's eq. (10) with k_c = 3/8, k_a = -3/4 and k_o = 3/8."""
    y = np.asarray(y, float)
    return 3 / 8 * np.log(np.abs(1 + y)) - 3 / 4 * np.log(np.abs(1 - y)) + 3 / 8 * np.log(np.abs(1 - y / 3))


CLOSED_FORMS.update({
    ("c_metric", "spherical", "inner"):
        (lambda t, r: t + 6 * _c_metric_rstar(6 / r), lambda t, r: t - 6 * _c_metric_rstar(6 / r),
         lambda t, r: (np.abs(r - 2) > 0.05) & (np.abs(r - 6) > 0.05)),
    ("c_metric", "spherical", "outer"):
        (lambda t, r: t + 6 * _c_metric_rstar(6 / r), lambda t, r: t - 6 * _c_metric_rstar(6 / r),
         lambda t, r: (np.abs(r - 2) > 0.05) & (np.abs(r - 6) > 0.05)),
    # P is outgoing here, and conserves the retarded tau - alpha r*.
    ("c_metric", "hong_teo", "inner"):
        (lambda tau, y: tau - _c_metric_rstar(y), lambda tau, y: tau + _c_metric_rstar(y),
         lambda tau, y: (np.abs(y - 1) > 0.02) & (np.abs(y - 3) > 0.02) & (y > -0.98)),
})
# The cylinders of t and phi, where each family runs straight, dt = k dphi, and conserves
# t - k phi: van Stockum's k = r(1 - r/R) moving right and -r(1 + r/R) moving left, and
# Godel's sinh r (cosh r - sqrt 2 sinh r) and -sinh r (cosh r + sqrt 2 sinh r), which are the
# numbers the captions state.
CYLINDERS = {
    ("stockum_dust", "cylindrical", "inside"): (-0.75, 0.25),
    ("stockum_dust", "cylindrical", "beyond"): (-3.75, -0.75),
    ("godel", "cylindrical", "inside"): (-(3 - math.sqrt(2)) / 2, (math.sqrt(2) - 1) / 2),
    ("godel", "cylindrical", "beyond"): (-(17 - math.sqrt(2)) / 2, -(3 - math.sqrt(2)) / 2),
    # Som and Raychaudhuri's block of t and phi is van Stockum's, with R = c/Omega.
    ("som_raychaudhuri", "cylindrical", "inside"): (-0.75, 0.25),
    ("som_raychaudhuri", "cylindrical", "beyond"): (-3.75, -0.75),
}
# The spinning string's cylinders: k = -(b r + a) moving left and b r - a moving right, with
# b r = 0.45 inside and 1.35 outside, in the proper radius and rescaled radius charts, and
# sqrt(R^2 + a^2) = 0.9 sqrt 2 at R = a in the circumference radius chart; in the helical chart
# the plane is Minkowski's and k = -r and r.
CYLINDERS.update({
    **{("spinning_string", system, view): (-(br + 0.9), br - 0.9)
       for system in ("proper_radius", "rescaled_radius") for view, br in (("inside", 0.45), ("outside", 1.35))},
    ("spinning_string", "circumference_radius", "outside"): (-0.9 * (math.sqrt(2) + 1), 0.9 * (math.sqrt(2) - 1)),
    ("spinning_string", "helical", "inside"): (-0.5, 0.5),
    ("spinning_string", "helical", "outside"): (-1.5, 1.5),
})


def _string_wave_slope(u):
    """A'(u) of the declared pulse A = exp(-4u^2)/2."""
    return -4 * u * np.exp(-4 * u ** 2)


def _string_wave_shift(X):
    """The integral of g_uu = (1/rho - 1) A'^2 from far before the pulse to u, on the line X of the
    moving string chart at b = 1/2, where rho = |X - A|, by a quadrature the tracing does not use."""
    def g_uu(s):
        return (1 / abs(X - math.exp(-4 * s * s) / 2) - 1) * (4 * s * math.exp(-4 * s * s)) ** 2
    return lambda u: np.array([quad(g_uu, -6.0, x, limit=200, epsabs=1e-12, epsrel=1e-12)[0] for x in np.atleast_1d(u)])


# The travelling wave on a string: a curve moving left keeps v minus the integral of g_uu, which
# is v + 2xA' where g_uu = -2xA'', at x = ell/4 or -ell/4, and a quadrature in the moving string
# chart; a ray moving right keeps u.
CLOSED_FORMS.update({
    ("string_wave", "null_conical", "toward"): (lambda u, v: v + _string_wave_slope(u) / 2, lambda u, v: u, None),
    ("string_wave", "null_conical", "away"): (lambda u, v: v - _string_wave_slope(u) / 2, lambda u, v: u, None),
    ("string_wave", "isotropic", "beside"): (lambda u, v: v + _string_wave_slope(u) / 2, lambda u, v: u, None),
    **{("string_wave", "moving_string", view): (lambda u, V, shift=_string_wave_shift(X): V - shift(u),
                                                lambda u, V: u, None)
       for view, X in (("behind", -0.25), ("ahead", 0.75))},
})
def _siklos_zstar(z):
    """The tortoise coordinate of Kaigorodov's Poincare plane with x divided out, the integral of
    sqrt(1 + s^3/2) from the boundary, by a quadrature the tracing does not use."""
    return np.array([quad(lambda s: math.sqrt(1 + s ** 3 / 2), 0.0, x, epsabs=1e-13, epsrel=1e-13)[0]
                     for x in np.atleast_1d(z)]).reshape(np.shape(z))


# Siklos's waves at L = 1. On the plane of u and v a curve moving left keeps v minus half the
# integral of H, which for the pulse is x^3 sqrt(pi) erf(2u)/8, and a ray moving right keeps u. With
# d/du divided out Kaigorodov's plane is -x^-5 dv^2 + x^-2 dx^2, whose rays keep v -+ (2/5) x^(5/2),
# written in each chart's own coordinates, and with d/dx divided out the Poincare plane's keep
# t -+ z_*. The stationary plane at rho = 0 is -du^2 - 2 du dv, whose other family keeps v + u/2.
CLOSED_FORMS.update({
    **{("siklos", "siklos", view): (lambda u, v, x=x: v - x ** 3 * math.sqrt(math.pi) * scipy_erf(2 * u) / 8,
                                     lambda u, v: u, None)
       for view, x in (("near", 1.0), ("far", 1.5))},
    ("siklos", "ozsvath_robinson_rozga", "centre"): (
        lambda u, v: v - math.sqrt(math.pi) * scipy_erf(2 * u) / 8, lambda u, v: u, None),
    ("siklos", "kaigorodov", "depth"): (lambda v, x: v + 0.4 * x ** 2.5, lambda v, x: v - 0.4 * x ** 2.5,
                                        lambda v, x: x > 0.02),
    ("siklos", "kaigorodov_poincare", "depth"): (lambda t, z: t + _siklos_zstar(z), lambda t, z: t - _siklos_zstar(z),
                                                 lambda t, z: z > 0.02),
    ("siklos", "kaigorodov_horospheric", "depth"): (lambda v, rho: v - 0.4 * np.exp(-2.5 * rho),
                                                    lambda v, rho: v + 0.4 * np.exp(-2.5 * rho), None),
    ("siklos", "kaigorodov_stationary", "plane"): (lambda u, v: v + u / 2, lambda u, v: u, None),
    ("siklos", "kaigorodov_homogeneous", "depth"): (lambda U, Z: (0.4 - U) * np.exp(5 * Z),
                                                    lambda U, Z: (-0.4 - U) * np.exp(5 * Z), None),
    ("siklos", "kaigorodov_kundt", "depth"): (lambda V, x: V * x ** 2 / math.sqrt(2) + 0.4 * x ** 2.5,
                                              lambda V, x: V * x ** 2 / math.sqrt(2) - 0.4 * x ** 2.5, None),
})
CLOSED_FORMS.update({where: (lambda t, phi, k=left: t - k * phi, lambda t, phi, k=right: t - k * phi, None)
                     for where, (left, right) in CYLINDERS.items()})
# Misner space: in Misner's plane one family keeps psi and the other T e^(psi/2), read through
# arcsinh so that the drift is measured on a scale the winding does not blow up; in the Milne and
# Rindler planes each family keeps ln|t| -+ chi or ln xi +- eta, the logarithms of the covering
# plane's null coordinates.
CLOSED_FORMS.update({
    ("misner", "misner", "plane"): (lambda T, psi: psi, lambda T, psi: np.arcsinh(T * np.exp(psi / 2)), None),
    ("misner", "milne", "plane"): (lambda t, chi: np.log(-t) - chi, lambda t, chi: np.log(-t) + chi,
                                   lambda t, chi: t < -0.02),
    ("misner", "rindler", "plane"): (lambda eta, xi: np.log(xi) + eta, lambda eta, xi: np.log(xi) - eta,
                                     lambda eta, xi: xi > 0.02),
    # Grant's charts of Gott's spacetime are Misner's Rindler and Milne planes.
    ("gott_time_machine", "grant_rindler", "plane"): (lambda eta, xi: np.log(xi) + eta,
                                                      lambda eta, xi: np.log(xi) - eta,
                                                      lambda eta, xi: xi > 0.02),
    ("gott_time_machine", "grant_milne", "plane"): (lambda t, chi: np.log(-t) - chi, lambda t, chi: np.log(-t) + chi,
                                                    lambda t, chi: t < -0.02),
})
# Ori's core: on a cylinder of the time and z where g_zz = k - (time), one family keeps z and the
# other (time - k) e^(z/2), with k = f or e x^2 there, read through arcsinh as Misner's is; on the
# Brinkmann plane each family keeps u or v.
CLOSED_FORMS.update({
    **{("ori_time_machine", system, view): (lambda T, z: z, lambda T, z, k=k: np.arcsinh((T - k) * np.exp(z / 2)), None)
       for system, view, k in (("vacuum_core", "centre", 0.0), ("vacuum_core", "off_centre", 0.5),
                               ("foliation", "centre", 0.0), ("foliation", "off_centre", 2.0))},
    ("ori_time_machine", "brinkmann", "plane"): (lambda u, v: v, lambda u, v: u, None),
})
# What light launched along each family of a cylinder does, as its caption says: stays on the
# cylinder as a null geodesic, or is turned toward or away from the axis.
TURNING = {
    ("stockum_dust", "cylindrical", "inside"): ("away", "geodesic"),
    ("stockum_dust", "cylindrical", "beyond"): ("away", "toward"),
    ("godel", "cylindrical", "inside"): ("away", "geodesic"),
    ("godel", "cylindrical", "beyond"): ("away", "toward"),
    ("som_raychaudhuri", "cylindrical", "inside"): ("away", "geodesic"),
    ("som_raychaudhuri", "cylindrical", "beyond"): ("away", "toward"),
    # The spinning string is flat, so light launched along a circle about it leaves for larger r.
    **{("spinning_string", system, view): ("away", "away")
       for system in ("proper_radius", "helical") for view in ("inside", "outside")},
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


def _kerr_de_sitter_forms():
    """Kerr-de Sitter's axis in each chart, for each sign of Lambda: ct -+ r_* in Carter's chart,
    v and v - 2r_* in the ingoing Kerr chart, u + 2r_* and u in the outgoing one, and in the
    Kerr-Schild chart c tau - s -+ r_* with s = c(tau - t), whose ingoing form is c tau plus the
    primitive of 1/(1 - Lambda r^2/3), de Sitter's own. r_* and s are slices.py's, zero at r = 0."""
    import slices
    out = {}
    for sign, suffix, edges in (("de_sitter", "", (0.2787501957398, 0.7847845522526, 3.2322636072011, 15 ** 0.5)),
                                ("anti_de_sitter", "_ads", (0.1368868045598, 0.8594395618858))):
        star, shift = slices.kds_rstar(sign), slices.kds_schild(sign)

        def away(x, r, edges=edges):
            return np.all([np.abs(r - z) > 0.05 for z in edges], axis=0)
        out[("kerr_de_sitter", "boyer_lindquist", "axis" + suffix)] = (
            lambda t, r, f=star: t + f(r), lambda t, r, f=star: t - f(r), away)
        out[("kerr_de_sitter", "kerr_ingoing", "axis" + suffix)] = (
            lambda v, r: v, lambda v, r, f=star: v - 2 * f(r), away)
        out[("kerr_de_sitter", "kerr_outgoing", "axis" + suffix)] = (
            lambda u, r, f=star: u + 2 * f(r), lambda u, r: u, away)
        out[("kerr_de_sitter", "kerr_schild", "axis" + suffix)] = (
            lambda T, r, f=star, g=shift: T - g(r) + f(r), lambda T, r, f=star, g=shift: T - g(r) - f(r), away)
    return out


CLOSED_FORMS.update(_kerr_de_sitter_forms())

# Hartle and Thorne's star: the tortoise coordinate of ht_star on the axis and on the equator, and
# the two quadratures of the Painleve-Gullstrand plane, which is Schwarzschild's on both.
CLOSED_FORMS.update({
    ("hartle_thorne", "hartle_thorne", "axis"):
        (lambda t, r: t + ht_star(r, 1.0), lambda t, r: t - ht_star(r, 1.0), None),
    ("hartle_thorne", "hartle_thorne", "equator"):
        (lambda t, r: t + ht_star(r, -0.5), lambda t, r: t - ht_star(r, -0.5), None),
    ("hartle_thorne", "painleve_gullstrand", "axis"): (lambda t, r: t + pg_in(r), lambda t, r: t - pg_out(r), None),
    ("hartle_thorne", "painleve_gullstrand", "equator"): (lambda t, r: t + pg_in(r), lambda t, r: t - pg_out(r), None),
})


def _bv_rstar(r):
    """Reissner-Nordstrom's tortoise coordinate at M = 1 and q = 24/25, r_+ = 32/25 and r_- = 18/25:
    r_* = r + (r_+^2 ln|r/r_+ - 1| - r_-^2 ln|r/r_- - 1|)/(r_+ - r_-)."""
    rp, rm = 1.28, 0.72
    return r + (rp ** 2 * np.log(np.abs(r / rp - 1)) - rm ** 2 * np.log(np.abs(r / rm - 1))) / (rp - rm)


def _bv_away(side):
    """Off the shell on the side where the metric is Reissner-Nordstrom's, and off both horizons and r = 0."""
    return lambda w, r: (side * w > 0.02) & (np.abs(r - 1.28) > 0.05) & (np.abs(r - 0.72) > 0.05) & (r > 0.05)


# Bonnor and Vaidya's charged shell: on the Reissner-Nordstrom side an ingoing ray keeps v and an
# outgoing one v - 2 r_*, and in the outgoing chart u + 2 r_* and u.
CLOSED_FORMS.update({
    ("bonnor_vaidya", "eddington_finkelstein_ingoing", "shell"):
        (lambda v, r: v, lambda v, r: v - 2 * _bv_rstar(r), _bv_away(1)),
    ("bonnor_vaidya", "eddington_finkelstein_outgoing", "shell"):
        (lambda u, r: u + 2 * _bv_rstar(r), lambda u, r: u, _bv_away(-1)),
    # In the homothetic chart an ingoing ray keeps V.
    ("bonnor_vaidya", "homothetic", "scaling"): (lambda V, R: V, None, None),
})


def verify(metrics=()):
    """Trace rays as the page does and measure how far each family's closed form drifts.
    Given metric ids, check only their closed forms and dust, skipping the checks of turning
    and of the principal null rays. Returns the number of failures."""
    failures = 0
    specs = {(s.metric, s.system, s.view): s for s in DIAGRAMS}

    def wanted(metric_id):
        return not metrics or metric_id in metrics
    forms = {where: form for where, form in CLOSED_FORMS.items() if wanted(where[0])}
    if wanted("frw"):
        solver = Chart(specs[("frw", "comoving_spherical", "radial")]).solver

        def eta(t):
            return np.array([quad(lambda s: 1 / solver.values("a", np.array([s]))[0][0], 1e-12, x, limit=400)[0]
                             for x in np.atleast_1d(t)])
        forms[("frw", "comoving_spherical", "radial")] = (lambda t, r: eta(t) + r, lambda t, r: eta(t) - r,
                                                          lambda t, r: t > 0.02)
    for view in ("axis", "equator"):
        where = ("robinson_trautman", "axisymmetric", view)
        if wanted(where[0]):
            # An ingoing ray keeps the Kruskal V it has once the fronts are round, found by an
            # integrator the tracing does not use, and an outgoing one its u.
            fronts = Chart(specs[where]).solver
            forms[where] = (lambda u, r, fronts=fronts: np.arctan(fronts.kruskal_v(u, r)), lambda u, r: u,
                            lambda u, r: (u > 0.01) & (r > 0.05))
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
    if wanted("frw"):
        t_sing = solver.t_sing
        t = np.linspace(0.01, 2.0, 400)
        eds = float(np.max(np.abs(solver.values("a", t)[0] / (t / -t_sing) ** (2 / 3) - 1)))
        ok = abs(t_sing + 2 / 3) < 1e-8 and eds < 1e-7
        failures += not ok
        print(f"FRW dust: the bang {t_sing:.9f} from a = 1, Einstein-de Sitter -2/3; "
              f"a(t) against (t/t0)^(2/3) to {eds:.1e}  {'ok' if ok else 'FAILED'}")
    if wanted("bianchi"):
        bianchi = Chart(specs[("bianchi", "type_i_cartesian", "tx")]).solver
        y = bianchi.state(np.array([1e-4, 2e-4]))
        p = [float(np.log(y[2 * i][1] / y[2 * i][0]) / np.log(2)) for i in range(3)]
        ok = abs(sum(p) - 1) < 1e-3 and abs(sum(q * q for q in p) - 1) < 1e-3
        failures += not ok
        print(f"Bianchi I dust: exponents at the singularity {', '.join(f'{q:.4f}' for q in p)}, "
              f"on the Kasner circle  {'ok' if ok else 'FAILED'}")
    if wanted("kantowski_sachs"):
        ks = Chart(specs[("kantowski_sachs", "comoving", "tr")]).solver
        eta = np.linspace(-1.5, 1.5, 301)
        t = math.pi / 2 + eta + np.sin(eta) * np.cos(eta)
        miss = max(float(np.max(np.abs(ks.values("a", t)[0] / (1 + eta * np.tan(eta)) - 1))),
                   float(np.max(np.abs(ks.values("b", t)[0] / np.cos(eta) ** 2 - 1))))
        ok = abs(ks.t_sing + math.pi / 2) < 1e-6 and miss < 1e-6
        failures += not ok
        print(f"Kantowski-Sachs dust: the first singularity {ks.t_sing:.9f} from the widest moment, -pi/2; a and b "
              f"against 1 + eta tan(eta) and cos^2(eta) to {miss:.1e}  {'ok' if ok else 'FAILED'}")
    if wanted("alcubierre") or wanted("natario"):
        alcubierre = Chart(specs[("alcubierre", "cartesian", "tx")])
        natario = Chart(specs[("natario", "cartesian_flow", "tx")])
        T, X = np.meshgrid(np.linspace(-2, 2, 201), np.linspace(-3, 3, 301))
        gap = max(float(np.nanmax(np.abs(alcubierre.fn[k](T, X) - natario.fn[k](T, X)))) for k in ("g00", "g0r", "grr"))
        ok = gap < 1e-12
        failures += not ok
        print(f"Natario against Alcubierre on the axis: the metrics on the plane differ by {gap:.1e}  {'ok' if ok else 'FAILED'}")
    if metrics:
        return failures
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
    """Every view and figure has a caption, and no two share a place."""
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


def main(argv=None):
    check_table()
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--metric", action="append", default=[], help="redraw only this spacetime, repeatable")
    parser.add_argument("--verify", action="store_true",
                        help="check the rays against closed forms instead of writing")
    parser.add_argument("--slices", action="store_true",
                        help="rewrite only the slices of the embedding diagrams, tracing no ray")
    args = parser.parse_args(argv)
    unknown = set(args.metric) - set(by_metric())
    if unknown:
        parser.error(f"no diagram is drawn for {sorted(unknown)}")
    if args.verify:
        return 1 if verify(set(args.metric)) else 0
    if args.slices:
        rewrite_slices(set(args.metric))
        return 0
    write(set(args.metric))
    return 0


if __name__ == "__main__":
    sys.exit(main())
