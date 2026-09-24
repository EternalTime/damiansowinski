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
held fixed, the metric is the 2 by 2 block h = [[g_00, g_0r], [g_0r, g_rr]], and a null
direction (dx^0, dr) solves

    g_00 (dx^0)^2 + 2 g_0r dx^0 dr + g_rr dr^2 = 0,

which has two real roots exactly where D = g_0r^2 - g_00 g_rr > 0. They are written as

    P ~ (-g_0r + sqrt(D), g_00) ~ (g_rr, -g_0r - sqrt(D))
    M ~ (-g_0r - sqrt(D), g_00) ~ (g_rr, -g_0r + sqrt(D))

taking at each point whichever form does not vanish, which is what lets one formula
serve a diagonal block, the Eddington-Finkelstein charts' g_rr = 0 and a double null
chart's g_uu = g_vv = 0. The block is divided by its largest entry first, which changes
no direction and keeps the numbers finite next to a horizon. P and M are continuous line
fields: P is the family that conserves an advanced coordinate, drawn as the ingoing or
left moving family, and M the retarded one. --verify confirms that labelling chart by
chart. Their integral curves are traced both ways by fourth order Runge-Kutta in the
drawing's own unit square, from seeds spaced evenly along its four edges.

A fixed coordinate null curve is a light ray only if nothing accelerates it out of the
plane, Gamma^A_ab k^a k^b = 0 for every fixed coordinate A. The planes in DIAGRAMS pass
that test, except Godel's, whose caption says its curves are not light rays.


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
  throat      d_r R = 0, with R^2 = g_theta_theta the areal radius.
  apparent    |grad R|^2 = 0 where grad R is not zero: marginally trapped spheres.
  singular    an edge along which the published Kretschmann scalar exceeds 1e8 and grows
              at least fiftyfold between 1e-4 and 1e-5 of the drawing, where the block is
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
"""

import argparse
import json
import math
import sys
from dataclasses import dataclass, field
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
    to_display: tuple = ((0, 1), (1, 0))            # rows: X = a.(x^0, r), Y = b.(x^0, r)
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
    input: str = None               # the declared input, in prose, printed beside the diagram


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

FRW_DUST = {"funcs": ["a"], "eqs": [["r", "r"]], "rates": [1.0], "params": {"k": 0}}
BIANCHI_DUST = {"funcs": ["a_1", "a_2", "a_3"], "eqs": [["x", "x"], ["y", "y"], ["z", "z"]],
                "rates": [-0.5, 1.5, 2.0]}

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
            input="Dust: $a(t)$ solved from the entry's own published $G^r{}_r = 0$ with $k = 0$, "
                  "starting from $a = 1$ and $\\dot a = H_0$ at the dashed line."),
    Diagram("frw", "comoving_spherical", "through", "through the observer", ("t", "r"),
            (0, 1.5, 0, 3), "$x\\;[c/H_0]$", "$ct\\;[c/H_0]$", {"k": 0}, EQUATOR, mirror=True,
            families=SIDEWAYS, cones=(4, 8), areal=True, dust=FRW_DUST, reference="$a = 1$",
            input="Dust: $a(t)$ solved from the entry's own published $G^r{}_r = 0$ with $k = 0$, "
                  "starting from $a = 1$ and $\\dot a = H_0$ at the dashed line."),
    Diagram("frw", "conformal_spherical", "radial", "$\\eta$ and $r$", ("\\eta", "r"), (0, 3, 0, 3),
            "$r\\;[c/H_0]$", "$\\eta\\;[c/H_0]$", {"k": 0}, EQUATOR, tau="eta", areal=True,
            dust=FRW_DUST, reference="$a = 1$",
            input="Dust, for the markers only: $a(\\eta)$ solved from this chart's own published "
                  "$G^r{}_r = 0$ with $k = 0$, starting from $a = 1$ and $a' = 1$ at the dashed line."),
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
    Diagram("kerr_newman", "boyer_lindquist", "radial", "$t$ and $r$ on the axis", ("t", "r"),
            (0, 4, -2, 2), "$r/(GM/c^2)$", "$ct/(GM/c^2)$",
            {"G": 1, "M": 1, "a": "3/5", "r_Q": "1/2"}, {"theta": "0", "phi": "0"},
            orient="ingoing"),
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
            input="Dust: the three scale factors solved from the entry's own published "
                  "$G^x{}_x = G^y{}_y = G^z{}_z = 0$, starting from $a_i = 1$ with rates "
                  "$(-0.5, 1.5, 2.0)\\,\\bar H$ at the dashed line, $\\bar H$ their mean."),
    Diagram("godel", "cartesian", "tx", "$t$ and $x$", ("t", "x"), (-2, 2, -2, 2), "$x$", "$t$",
            {"omega": 1}, {"y": "0", "z": "0"}, orient="vector", families=SIDEWAYS),
    Diagram("alcubierre", "cartesian", "tx", "$t$ and $x$ on the axis", ("t", "x"), (-3, 3, -2, 2),
            "$x/R$", "$ct/R$", {}, {"y": "0", "z": "0"}, families=SIDEWAYS, cones=(8, 7),
            functions={"v_s": "2", "f": _alcubierre_profile()},
            input="$v_s = 2$, and Alcubierre's own profile, which the entry names, "
                  "$f = [\\tanh\\sigma(r_s + R) - \\tanh\\sigma(r_s - R)]/(2\\tanh\\sigma R)$ "
                  "with $R = 1$ and $\\sigma = 4$."),
    Diagram("natario", "cartesian_flow", "tx", "$t$ and $x$ on the axis", ("t", "x"), (-3, 3, -2, 2),
            "$x/R$", "$ct/R$", {}, {"y": "0", "z": "0"}, families=SIDEWAYS, cones=(8, 7),
            functions=_natario_field(),
            input="$v_s = 2$, $n = f/2$ with Alcubierre's profile, and the zero expansion field "
                  "$X = v_s[(2n + \\rho n')\\,e_x - n'\\,x_r\\,(x_r, y, z)/\\rho]$, "
                  "$x_r = x - v_s t$, checked to have no divergence."),
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
]


# ---------------------------------------------------------------- the captions

CAPTIONS = {
    ("schwarzschild", "spherical", "radial"): [
        "Outside $r_s$ the cones narrow toward the vertical as $r \\to r_s$, because "
        "$dt/dr = \\pm(1 - r_s/r)^{-1}$ diverges there. The ingoing family piles up against $r_s$ "
        "toward $t \\to +\\infty$, and the outgoing family peels away from it out of $t \\to -\\infty$; "
        "neither crosses it in this chart.",
        "The published domain stops at $r_s$, and the region inside is hatched. Read there, the same "
        "components give cones lying on their side, because $r$ is the time. This chart alone cannot "
        "say whether that region is the black hole or the white hole; the cones follow the ingoing "
        "chart, which makes it the black hole, and point to $r = 0$. The Kretschmann scalar "
        "$12r_s^2/r^6$ is finite at $r_s$ and diverges only at $r = 0$.",
    ],
    ("schwarzschild", "eddington_finkelstein_ingoing", "finkelstein"): [
        "This chart crosses the horizon. The ingoing family is $v = $ const, drawn at 45° with "
        "$v - r$ as the vertical axis. The outgoing family has $dv/dr = 2(1 - r_s/r)^{-1}$, so it "
        "stands exactly vertical at $r_s$: the horizon is itself an outgoing ray that stays where it is.",
        "The cones cross $r_s$ smoothly and keep tipping. Inside, both future edges point to smaller "
        "$r$, so every future directed ray ends at $r = 0$.",
    ],
    ("schwarzschild", "eddington_finkelstein_ingoing", "chart"): [
        "The same rays against the chart's own coordinates. The ingoing family is $v = $ const and "
        "runs horizontally here, since $v$ is itself a null coordinate. The outgoing family turns "
        "vertical at $r_s$ and leans back toward smaller $r$ inside it.",
    ],
    ("schwarzschild", "eddington_finkelstein_outgoing", "finkelstein"): [
        "The retarded chart crosses the other horizon. The outgoing family is $u = $ const, drawn at "
        "45° with $u + r$ as the vertical axis. Inside $r_s$ both future edges point to larger $r$: "
        "this is the white hole, which nothing from outside can enter.",
        "It is the time reverse of the ingoing chart, as the published $g_{ur} = -1$ against "
        "$g_{vr} = +1$ says it must be.",
    ],
    ("schwarzschild", "eddington_finkelstein_outgoing", "chart"): [
        "The same rays against the chart's own coordinates. The outgoing family is $u = $ const and "
        "runs horizontally here, since $u$ is itself a null coordinate. The ingoing family turns "
        "vertical at $r_s$ and leans toward larger $r$ inside it.",
    ],
    ("frw", "comoving_spherical", "radial"): [
        "The entry leaves $a(t)$ free, so it is solved from the entry's own published "
        "$G^r{}_r = -2\\ddot a/a - \\dot a^2/a^2 - k/a^2$, set to zero, which is pressureless dust. "
        "Integrated back from the dashed line it reaches $a = 0$ at $H_0t = 2/3$ before it, the "
        "Einstein-de Sitter age, and that instant is $t = 0$ here.",
        "The cones are $dt/dr = \\pm a(t)$ and open out toward the bang, where the curvature "
        "diverges, so every ray leaves $t = 0$ almost flat and reaches only a finite comoving "
        "distance: the particle horizon. The dotted curve is where $|\\nabla R|^2 = 0$ for the areal "
        "radius $R = ar$, the Hubble sphere $R = c/H$, the apparent horizon of an observer at $r = 0$.",
    ],
    ("frw", "comoving_spherical", "through"): [
        "The same universe along a line through the observer at $r = 0$, drawn as $x = r$ on the "
        "right and $x = -r$ on the left, which spherical symmetry makes exact. The past light cone "
        "of an event on the axis flares out as it runs back toward the bang, and the Hubble sphere "
        "shows on both sides.",
    ],
    ("frw", "conformal_spherical", "radial"): [
        "In conformal time the radial block is $a^2(-d\\eta^2 + dr^2)$ for $k = 0$. The scale factor "
        "multiplies both terms and drops out of the null condition, so the rays are straight 45° "
        "lines whatever $a(\\eta)$ is.",
        "The scale factor is still solved, because the apparent horizon and the Kretschmann scalar "
        "need it. The Hubble sphere then sits at $r = \\eta/2$, half the comoving radius of the "
        "particle horizon.",
    ],
    ("ellis_bronnikov", "spherical", "radial"): [
        "The radial block is $-dt^2 + dr^2$, exactly flat, so the rays are straight 45° lines that "
        "pass through $r = 0$ without feeling it. The throat is not in this block at all.",
        "It is in $g_{\\theta\\theta} = r^2 + \\ell^2$. Its square root, the areal radius $R$, is least "
        "at $r = 0$, where $\\partial_r R = 0$ is marked. The faint vertical lines are the spheres "
        "$R = 1.5\\ell$, $2\\ell$ and $3\\ell$, the same sizes on both sides. The Kretschmann scalar "
        "$12\\ell^4/(r^2 + \\ell^2)^4$ is finite everywhere.",
    ],
    ("morris_thorne", "spherical", "radial"): [
        "The entry leaves $\\Phi(r)$ and $b(r)$ free. With $\\Phi = 0$ and $b = b_0^2/r$ it is the "
        "Ellis-Bronnikov wormhole, with $r^2 = r_{\\rm EB}^2 + \\ell^2$ and $b_0 = \\ell$, and these rays "
        "conserve $t \\mp \\sqrt{r^2 - b_0^2} = t \\mp r_{\\rm EB}$: they are the same rays as that "
        "entry's.",
        "In this areal chart the cones close toward $r = b_0$, as they would at a horizon, because "
        "$g_{rr} = (1 - b_0^2/r^2)^{-1}$ diverges there. But $g_{tt} = -1$ stays finite, so "
        "$\\partial_t$ is timelike right up to the throat, and what ends at $b_0$ is the chart, not "
        "the spacetime. The rays reach it in finite $t$ and carry on into the other mouth, which this "
        "chart does not cover. Below $b_0$ the formula gives a block with no null directions at all.",
    ],
    ("minkowski", "spherical", "radial"): [
        "A control. The block is $-c^2dt^2 + dr^2$ and every ray is at 45°. Each ingoing ray meets "
        "an outgoing one on the axis $r = 0$.",
    ],
    ("minkowski", "spherical_null", "radial"): [
        "A control for a double null chart, where $g_{uu} = g_{vv} = 0$ and only $g_{uv}$ is "
        "published. The null condition still gives both directions, which are the coordinate lines "
        "themselves, drawn here against $(v - u)/2$ and $(u + v)/2$.",
    ],
    ("minkowski", "cartesian", "tx"): [
        "A control. The block is $-c^2dt^2 + dx^2$ and every ray is at 45°.",
    ],
    ("minkowski", "rindler", "tx"): [
        "The Rindler chart covers the wedge $X > 0$ seen by observers of constant proper "
        "acceleration $a$, with $g_{TT} = -a^2X^2/c^4$. The cones close toward $X = 0$, where "
        "$g_{TT}$ vanishes, and a ray takes infinite $T$ to get there, "
        "$cT = \\pm(c^2/a)\\ln X + $ const.",
        "That is the Rindler horizon, a horizon of the accelerated observers only: the spacetime is "
        "flat, and its Kretschmann scalar is zero.",
    ],
    ("de_sitter", "static_spherical", "radial"): [
        "The cones close at the cosmological horizon $r = \\sqrt{3/\\Lambda}$, where "
        "$g^{rr} = 1 - \\Lambda r^2/3$ vanishes, on the far side from the observer at $r = 0$.",
        "The published domain is the observer's side. Beyond it $t$ is spacelike, and the cones "
        "follow the outgoing family toward larger $r$: the region the observer's own light reaches.",
    ],
    ("de_sitter", "static_spherical", "through"): [
        "The static chart along a line through the observer, $x = r$ on the right and $x = -r$ on "
        "the left, which spherical symmetry makes exact. The horizon shows on both sides at "
        "$x = \\pm\\sqrt{3/\\Lambda}$, and beyond it the cones point away from the observer.",
    ],
    ("de_sitter", "flat_slicing", "tx"): [
        "On this plane $ds^2 = -c^2dt^2 + e^{2Ht}dx^2$, so the cones narrow as $e^{-Ht}$ toward the "
        "future and open out toward the past.",
        "A ray covers only a finite comoving distance however long it runs, "
        "$x = \\pm(c/H)e^{-Ht} + $ const, so an observer at $x = 0$ has an event horizon.",
    ],
    ("anti_de_sitter", "static_global", "radial"): [
        "The cones never close, since $g^{rr} = 1 + r^2/L^2$ never vanishes. But "
        "$dt/dr = \\pm(1 + r^2/L^2)^{-1}$ falls off fast enough that a ray reaches $r \\to \\infty$ in "
        "the finite time $\\pi L/2$, and the rays flatten toward the right, where the boundary is.",
    ],
    ("anti_de_sitter", "static_global", "through"): [
        "The global chart along a line through the centre, reflected by spherical symmetry. A ray "
        "crosses the whole space and reaches the boundary on either side in finite time, $\\pi L/2$ "
        "from the centre, so the rays flatten toward both edges.",
    ],
    ("anti_de_sitter", "poincare", "tx"): [
        "At fixed $z$ the block is $(L^2/z^2)(-c^2dt^2 + dx^2)$, conformal to flat, so the rays are "
        "exact 45° lines at every $z$. The null structure of anti-de Sitter shows in its global chart "
        "instead.",
    ],
    ("rn_metric", "spherical", "radial"): [
        "$g^{rr}$ vanishes twice, at $r_\\pm = (r_s \\pm \\sqrt{r_s^2 - 4r_q^2})/2$, which is $0.8\\,r_s$ "
        "and $0.2\\,r_s$ for $r_q = 0.4\\,r_s$, and the cones close at both. Between them $r$ is the "
        "time and the cones point to smaller $r$. Inside $r_-$, $t$ is a time again.",
        "The chart alone does not orient the two inner regions; the cones there follow the ingoing "
        "family, as an ingoing chart carries them through both horizons. The Kretschmann scalar "
        "diverges at $r = 0$.",
    ],
    ("taub_nut", "spherical", "radial"): [
        "$g^{rr}$ vanishes at $r_+ = m + \\sqrt{m^2 + l^2}$, about $2.118\\,m$ for $l = m/2$, and inside "
        "it is the Taub region, where $r$ is the time and the cones, following the ingoing family, "
        "point to smaller $r$.",
        "The cross term $g_{t\\phi}$ does not enter this block, and nothing accelerates these curves "
        "out of the plane, so they are light rays although the solution is not spherically symmetric.",
    ],
    ("bertotti_robinson", "static", "radial"): [
        "The AdS₂ factor, as $-(r^2/b^2)\\,c^2dt^2 + (b^2/r^2)\\,dr^2$. With $dt/dr = \\pm b^2/r^2$ the "
        "rays take infinite $t$ to reach $r = 0$, where $g^{rr} = r^2/b^2$ vanishes: a horizon of "
        "Poincaré type, on the left edge. The two sphere keeps the radius $b$ throughout.",
    ],
    ("bertotti_robinson", "poincare", "tx"): [
        "$(b^2/x^2)(-c^2dt^2 + dx^2)$ times a two sphere of radius $b$: conformal to flat, so the rays "
        "are at 45°. Here $x$ runs along the AdS₂ factor rather than across space; the boundary is "
        "$x \\to 0$ and the horizon of Poincaré type is $x \\to \\infty$.",
    ],
    ("interior_schwarzschild", "spherical", "radial"): [
        "The whole of the published domain $r \\in [0, R]$ of a star with $R = 1.5\\,r_s$. The cones are "
        "narrowest at the centre, where $|g_{tt}|$ is least and the redshift greatest, and they never "
        "close. They would close at the centre exactly when $3\\sqrt{1 - r_s/R} = 1$, which is "
        "Buchdahl's $R = 9r_s/8$.",
    ],
    ("interior_schwarzschild", "spherical", "through"): [
        "The line through the centre of the star, reflected by spherical symmetry. Rays cross the "
        "centre smoothly, and the cones are narrowest there.",
    ],
    ("kerr", "boyer_lindquist", "radial"): [
        "Only on the axis is the plane of $t$ and $r$ in Boyer-Lindquist coordinates a plane of light "
        "rays: off it, $\\Gamma^\\theta$ and $\\Gamma^\\phi$ turn every fixed angle null curve out of "
        "the plane. $g^{rr} = \\Delta/\\Sigma$ vanishes at both roots of $\\Delta$, "
        "$r_\\pm = GM/c^2 \\pm \\sqrt{(GM/c^2)^2 - a^2}$, and the cones close at both; between them "
        "they point to smaller $r$, following the ingoing family.",
        "The published domain begins at $r_+$, read at the outer zero of $g^{rr}$, and below it the "
        "drawing is hatched. On the axis the ergosurface touches the horizon, so the plane stays "
        "Lorentzian, and the Kretschmann scalar stays finite at $r = 0$, because the ring singularity "
        "lies in the equatorial plane.",
    ],
    ("kerr_newman", "boyer_lindquist", "radial"): [
        "The same on the axis with charge: $g^{rr}$ vanishes at "
        "$r_\\pm = GM/c^2 \\pm \\sqrt{(GM/c^2)^2 - a^2 - r_Q^2}$, the cones close at both, and between "
        "them they point to smaller $r$.",
        "The published domain begins at $r_+$, read at the outer zero of $g^{rr}$. The Kretschmann "
        "scalar again stays finite at $r = 0$ on the axis.",
    ],
    ("kasner", "cartesian", "tx"): [
        "Along $x$ the scale factor $t^{-2/7}$ grows toward the singularity, so the cones close up as "
        "$t \\to 0$: $dx/dt = \\pm t^{2/7}$. The Kretschmann scalar $-16p_1p_2p_3/t^4$ diverges on "
        "the bottom edge.",
    ],
    ("kasner", "cartesian", "tz"): [
        "Along $z$ the scale factor $t^{6/7}$ goes to zero at the singularity, and the cones open out "
        "flat: $dz/dt = \\pm t^{-6/7}$. The same singularity closes the cones along $x$ and opens "
        "them along $z$.",
    ],
    ("bianchi", "type_i_cartesian", "tx"): [
        "The singularity comes about $0.378/\\bar H$ before the dashed line and is placed at $t = 0$. "
        "Next to it the scale factors run as powers of $t$ whose exponents lie on the Kasner circle: "
        "the dust model is Kasner at its singularity.",
        "So along $x$ the cones close toward $t = 0$, as in Kasner's contracting direction, and open "
        "again later as $a_1$ turns round.",
    ],
    ("godel", "cartesian", "tx"): [
        "Gödel has no time function: its published $g^{tt} = 2\\omega^2$ is positive, so the surfaces "
        "$t = $ const are not spacelike, and the cones are oriented instead by $\\partial_t$, which is "
        "timelike everywhere.",
        "The block is flat, but these are null curves, not light rays: $\\Gamma^y{}_{tx} = -e^{-x}$ "
        "turns every one of them out of the plane.",
    ],
    ("alcubierre", "cartesian", "tx"): [
        "A bubble moving at twice the speed of light along $x = 2ct$, drawn on its axis of motion, "
        "where these are light rays because $\\partial_y f = \\partial_z f = 0$ there. Inside the "
        "bubble and far outside, the cones are Minkowski's; in the walls they tilt with the shift "
        "$v_s f$.",
        "The dashed lines are $g^{xx} = 1 - v_s^2 f^2 = 0$, at $v_s f = 1$. For $v_s = 2$ that is also "
        "where a forward ray keeps pace with the bubble, $f = 1 - 1/v_s$, so here they are the two "
        "horizons: forward rays from inside stall at the front wall, and forward rays from behind "
        "stall at the back wall. At other speeds the two places differ.",
    ],
    ("natario", "cartesian_flow", "tx"): [
        "On the axis of motion the field reduces to $u = 2nv_s = v_s f$, which is Alcubierre's shift "
        "exactly, so this is Alcubierre's diagram, and the two blocks agree there to rounding.",
        "The drives differ only off the axis, where Natário's flow slides space sideways instead of "
        "compressing it.",
    ],
    ("krasnikov", "cylindrical", "tx"): [
        "Outside the tube $k = 1$ and the cones are Minkowski's. Inside, $k$ comes close to "
        "$\\delta - 1$ and the edge moving left tips below the horizontal, so a ray going back toward "
        "$x = 0$ loses about $0.8$ in $ct$ for every unit of $x$ it covers. The edge moving right, "
        "$c\\,dt = dx$, is the same everywhere, and the cones follow it, since $t$ orients the cones "
        "outside the tube.",
        "The dash dot line is $g^{tt} = 0$, where $k = 0$ and the surfaces $t = $ const stop being "
        "spacelike. The entry's Kretschmann scalar is written as $0/0$ on the axis, so this plane is "
        "not tested for singularities.",
    ],
    ("pp_wave", "exact_plane_wave", "tz"): [
        "A plane wave has no $t$ in the collection's charts; $u$ and $v$ are null. The nearest pair to "
        "$t$ and $x$ is the plane the wave travels in, drawn with $u = t - z$ and $v = (t + z)/2$, on "
        "the wave's axis $x = y = 0$. There the profile $A(x^2 - y^2) + 2Bxy$ vanishes whatever $A$ "
        "and $B$ are, so the block is flat and the rays are at 45°.",
        "Off the axis $\\Gamma^x{}_{uu} = -(xA + yB)$ pulls every ray toward it. The wave's effect is "
        "focusing, which no fixed plane can show.",
    ],
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

    def __init__(self, metric_id, system_id, time_name, funcs, eqs, rates, params=None):
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
        y0 = [v for rate in rates for v in (1.0, float(rate))]
        self.back = solve_ivp(rhs, (0, -20), y0, events=singular, rtol=1e-11, atol=1e-13, dense_output=True)
        self.fwd = solve_ivp(rhs, (0, 20), y0, events=singular, rtol=1e-11, atol=1e-13, dense_output=True)
        self.t_sing = self.back.t_events[0][0]
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

    def values(self, name, t):
        """(a, a dot, a double dot) of one scale factor at chart times t."""
        i = self.funcs.index(name)
        y = self.state(t)
        acc = self.f(*y)
        return y[2 * i], y[2 * i + 1], np.asarray(acc[i], dtype=float) * np.ones_like(y[0])


_SOLVERS = {}


def dust_solver(metric_id, system_id, time_name, dust):
    key = (metric_id, system_id, json.dumps(dust, sort_keys=True))
    if key not in _SOLVERS:
        _SOLVERS[key] = DustSolver(metric_id, system_id, time_name, dust["funcs"], dust["eqs"],
                                   dust["rates"], dust.get("params"))
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
        self.solver = (dust_solver(spec.metric, spec.system, spec.plane[0], spec.dust)
                       if spec.dust else None)

        def prep(expr):
            expr = sp.sympify(expr)
            for name, rep in (spec.functions or {}).items():
                expr = expr.replace(reader.parameters[name].func, _as_lambda(reader, name, rep)).doit()
            return expr.subs(subs)

        g = published_matrix(reader, entry, "metric_components")
        gi = published_matrix(reader, entry, "inverse_metric_components")
        self.fn = {}
        for name, expr in (("g00", g[a, a]), ("g0r", g[a, b]), ("grr", g[b, b]),
                           ("gi00", gi[a, a]), ("gi0r", gi[a, b]), ("girr", gi[b, b])):
            self.fn[name] = self.lambdify(prep(expr))
        K = prep(reader(strip_lhs(entry["kretschmann"]))) if spec.kretschmann else sp.Integer(0)
        self.fn["K"] = self.lambdify(K)
        tau = sp.sympify(spec.tau, locals={str(self.x0): self.x0, str(self.xr): self.xr})
        self.fn["dtau0"] = self.lambdify(sp.diff(tau, self.x0))
        self.fn["dtaur"] = self.lambdify(sp.diff(tau, self.xr))

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
                    args += list(self.solver.values(name, x0))
            with np.errstate(all="ignore"):
                out = np.asarray(f(x0, r, *args))
            if np.iscomplexobj(out):
                real = np.abs(out.imag) <= 1e-9 * np.abs(out.real) + 1e-12
                out = np.where(real, out.real, np.nan)
            return np.broadcast_to(out.astype(float), shape)

        return call

    def null_dirs(self, x0, r):
        """The directions P and M in (dx^0, dr), and D. NaN where there are none."""
        g00, g0r, grr = self.fn["g00"](x0, r), self.fn["g0r"](x0, r), self.fn["grr"](x0, r)
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
        g00, g0r, grr = (self.fn[k](*at)[0] for k in ("g00", "g0r", "grr"))

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


def _as_lambda(reader, name, rep):
    fn = reader.parameters[name]
    body = sp.sympify(rep, locals={**reader.local, **{str(arg): arg for arg in fn.args}})
    return sp.Lambda(fn.args, body)


# ---------------------------------------------------------------- one view of a chart

class Plot:
    """A Diagram drawn: rays, cones and markers in the unit square of its axes."""

    def __init__(self, chart):
        spec = chart.spec
        self.c = chart
        self.A = np.array(spec.to_display, float)
        self.Ai = np.linalg.inv(self.A)
        X0, X1, Y0, Y1 = spec.box
        self.lo = np.array([X0, Y0], float)
        self.span = np.array([X1 - X0, Y1 - Y0], float)

    def to_chart(self, q):
        """Drawn axes (X, Y) to the chart's (x^0, r): the map is linear, with no offset."""
        p = np.asarray(q) @ self.Ai.T
        return p[..., 0], p[..., 1]

    def from_unit(self, u):
        return np.asarray(u) * self.span + self.lo

    def dirs_unit(self, u):
        """P and M at unit square points u, as unit vectors of the square, with D, P and M."""
        x0, r = self.to_chart(self.from_unit(u))
        P, M, D = self.c.null_dirs(x0, r)
        out = []
        for k in (P, M):
            d = (k @ self.A.T) / self.span
            with np.errstate(all="ignore"):
                out.append(d / np.linalg.norm(d, axis=-1, keepdims=True))
        return out[0], out[1], D, P, M

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
        where the block has no null direction at all is recorded as such."""
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
            lorentzian = self.c.null_dirs(*self.to_chart(self.from_unit(at(1e-4))))[2] > 0
            with np.errstate(all="ignore"):
                if ((near > 1e8) & (near / far > 50) & lorentzian).mean() > 0.5:
                    out.append(name)
        return out

    def hatch(self):
        """Where a chart point lies outside the entry's published domains, as polygons."""
        domains = parse_domains(self.c)
        if not domains:
            return []
        UU, VV, x0, r = self.grid(161)
        outside = np.zeros_like(UU, dtype=bool)
        for name, values in ((self.c.spec.plane[0], x0), (self.c.spec.plane[1], r)):
            if name in domains:
                lo, hi, lo_open, hi_open = domains[name]
                if lo is not None:
                    outside |= (values < lo) | ((values == lo) & lo_open)
                if hi is not None:
                    outside |= (values > hi) | ((values == hi) & hi_open)
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
        lines = self.zero_set("girr")
        if lines:
            out.append({"kind": "grr", "lines": lines})
        if spec.areal:
            throat = self.zero_set("dRr", keep=lambda x0, r: fn["R"](x0, r) > 1e-6, drop_edge=True)
            if throat:
                out.append({"kind": "throat", "lines": throat})
            if not self.c.same_as_grr:
                apparent = self.zero_set("grad2", keep=lambda x0, r: np.abs(fn["dRr"](x0, r)) > 1e-6)
                apparent = [l for l in apparent if not same_line(l, throat)]
                if apparent:
                    out.append({"kind": "apparent", "lines": apparent})
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


def outer_root(chart):
    """The outermost zero of g^rr along the drawn radial range, where r_+ sits."""
    X = np.linspace(chart.spec.box[0], chart.spec.box[1], 20001)
    f = chart.fn["girr"](np.zeros_like(X), X)
    crossings = np.where(np.sign(f[:-1]) * np.sign(f[1:]) < 0)[0]
    if not crossings.size:
        return None
    i = crossings[-1]
    return float(X[i] - f[i] * (X[i + 1] - X[i]) / (f[i + 1] - f[i]))


def parse_domains(chart):
    """{coordinate: (low, high, low open, high open)} for the two coordinates of the plane."""
    out = {}
    for text in chart.entry.get("domains", []):
        t = text.replace("\\left", "").replace("\\right", "").replace("\\,", "")
        name, found, interval = t.partition("\\in")
        name, interval = name.strip(), interval.strip()
        if not found or name not in chart.spec.plane or len(interval) < 2:
            continue
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
        shown = sp.latex(number(value))
        parts.append(f"${names[plain]} = {shown}$")
    return ", ".join(parts)


def draw(spec):
    chart = Chart(spec)
    plot = Plot(chart)
    families = plot.rays()
    fields = BASE_FIELDS + (["einstein_tensor"] if spec.dust else [])
    view = {
        "id": spec.view, "label": spec.label, "plane": list(spec.plane), "families": list(spec.families),
        "xlabel": spec.xlabel, "ylabel": spec.ylabel, "ticks": axes(spec), "box": list(spec.box),
        "to_display": [list(map(float, row)) for row in spec.to_display], "mirror": spec.mirror,
        "rays": {name: [rounded(thin(l, 0.0006)) for l in families[i]] for i, name in ((0, "P"), (1, "M"))},
        "cones": plot.cones(), "markers": plot.markers(), "hatch": plot.hatch(),
        "settings": settings(spec, chart.entry), "input": spec.input,
        "caption": CAPTIONS[(spec.metric, spec.system, spec.view)],
        "source": {"fields": fields, "version": build.diagram_source_version(chart.entry, fields)},
    }
    if spec.areal_contours:
        view["areal"] = plot.level_sets("R", spec.areal_contours)
    return view


def write(metric_ids=None):
    if not DIAGRAMS_DIR.exists():
        DIAGRAMS_DIR.mkdir(parents=True)
    by_metric = {}
    for spec in DIAGRAMS:
        by_metric.setdefault(spec.metric, []).append(spec)
    for metric_id, specs in by_metric.items():
        if metric_ids and metric_id not in metric_ids:
            continue
        systems = {}
        for spec in specs:
            systems.setdefault(spec.system, []).append(draw(spec))
            print(f"  {key(spec)}", flush=True)
        path = DIAGRAMS_DIR / f"{metric_id}.json"
        path.write_text(json.dumps({"metric": metric_id, "systems": systems},
                                   ensure_ascii=False, separators=(",", ":")) + "\n", encoding="utf-8")
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
    for where, (own_P, own_M, keep) in forms.items():
        plot = Plot(Chart(specs[where]))
        drift, spread = {0: 0.0, 1: 0.0}, {0: 0.0, 1: 0.0}
        for s in np.linspace(0.05, 0.95, 7):
            for seed in ((s, 0.001), (0.999, s), (s, 0.999), (0.001, s)):
                for family, own, other in ((0, own_P, own_M), (1, own_M, own_P)):
                    if own is None:
                        continue
                    line = plot.ray_through(np.array(seed), family)
                    line = line[(line >= 0).all(1) & (line <= 1).all(1)]
                    x0, r = plot.to_chart(plot.from_unit(line))
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
    print(f"Natario against Alcubierre on the axis: the blocks differ by {gap:.1e}  {'ok' if ok else 'FAILED'}")
    return failures


def check_table():
    """Every view has a caption, and no two views share a place."""
    places = [(spec.metric, spec.system, spec.view) for spec in DIAGRAMS]
    if len(set(places)) != len(places):
        raise SystemExit("two views in DIAGRAMS share a metric, system and view id")
    missing = [key(spec) for spec in DIAGRAMS if (spec.metric, spec.system, spec.view) not in CAPTIONS]
    stray = set(CAPTIONS) - set(places)
    if missing or stray:
        raise SystemExit(f"views without a caption: {missing}; captions without a view: {sorted(stray)}")


def main(argv=None):
    check_table()
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--metric", action="append", default=[], help="redraw only this spacetime, repeatable")
    parser.add_argument("--verify", action="store_true",
                        help="check the rays against closed forms instead of writing")
    args = parser.parse_args(argv)
    unknown = set(args.metric) - {spec.metric for spec in DIAGRAMS}
    if unknown:
        parser.error(f"no diagram is drawn for {sorted(unknown)}")
    if args.verify:
        return 1 if verify() else 0
    write(set(args.metric))
    return 0


if __name__ == "__main__":
    sys.exit(main())
