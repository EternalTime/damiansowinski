#!/usr/bin/env python3
"""Embedding diagrams of the spacetimes, drawn from the published metrics.

An embedding diagram draws a two dimensional slice of a spacetime, the equatorial plane at
one moment, as a surface in ordinary flat three dimensional space, so that distances
measured along the surface are the metric distances. Every slice drawn here is
the surface of one spatial coordinate x and the angle phi of one coordinate system, every
other coordinate held fixed, read from the system's published metric_components through
the Reader of verify_metrics.py beside this file, in the x^0 = cT chart the collection
writes, with c = 1, by the load and published_matrix the null ray diagrams use. Its
metric is

    g_xx(x) dx^2 + g_phiphi(x) dphi^2,

with no cross term, which is checked, so turning phi carries the slice onto itself and it
is a surface of revolution about an axis z. A circle of constant x has circumference
2 pi sqrt(g_phiphi), so it is drawn at the radius rho = sqrt(g_phiphi) from the axis, and
the distance sqrt(g_xx) dx between neighbouring circles is the hypotenuse of drho and dz:

    dz/dx = sqrt(g_xx - (drho/dx)^2),

which for an areal radius, rho = r, is dz/dr = sqrt(g_rr - 1). The surface exists exactly
where g_xx >= (drho/dx)^2. Where a circle grows faster than the distance across to it, no
surface of revolution in flat space carries the slice, and the file says where the
construction stops and why rather than drawing past it; where every circle does, as on the
hyperbolic plane of anti-de Sitter space, the slice is drawn in three dimensional Minkowski
space instead, dX^2 + dY^2 - dZ^2, where it climbs at dZ/dx = sqrt((drho/dx)^2 - g_xx).
A flat slice is drawn as the plane it is, with what the spacetime does marked on it as
curves and points that are no circles: a ring of free particles. Three spacetimes whose slices
are flat draw a quantity as a height over a plane instead, a GridPiece: Alcubierre's drive the
expansion of the observers who ride its slices, Natario's the energy density they measure,
with the lines of its flow on it, and the Krasnikov tube 1 - k, how far it tips the light
cone. The height is sampled on a grid and drawn as flat triangles, each within GRID_SAG of the
drawing's size of it in space.

    python3 -m venv /tmp/mfs-venv && /tmp/mfs-venv/bin/pip install sympy numpy scipy contourpy
    /tmp/mfs-venv/bin/python _tools/derivations/embedding.py
    python3 _tools/build_mfs_data.py

writes MFS/assets/data/embedding/<metric_id>.json, one file per spacetime that has a
diagram, and the second command stamps each file's version into the index. Pass
--metric <metric_id> to redraw one spacetime, repeatable, and --verify to run every check
and print it without writing anything.


The construction
----------------

drho/dx is taken in sympy from the published g_phiphi, and z by adaptive quadrature of
sqrt(g_xx - (drho/dx)^2) between neighbouring points of the profile; scipy's quad never
evaluates an end, so it integrates through the integrable divergence of g_xx at a throat.
An interval of x is halved until the profile strays less than SAG of the drawing's size
from its chord and the chord is shorter than STEP of it, and every circle the file marks
is a point of the profile. No closed form is used to draw anything.


What is checked
---------------

Every surface is measured as an application will draw it, from the rounded numbers the
file holds, against the published metric:

  along     each chord of the profile, and the profile end to end, against the proper
            distance between the same two values of x, the quadrature of sqrt(g_xx);
  across    the straight line in space from each point to the next one DPHI further
            round, against the metric length of the straight coordinate line
            between them;
  around    the circumference 2 pi rho at every point against 2 pi sqrt(g_phiphi);
  joins     where two pieces meet, as a star's surface meets the exterior, they meet at
            one point with one tangent, which says g_xx agrees on both sides;
  forms     the closed form each surface is known by, as Flamm's paraboloid, the interior
            Schwarzschild cap, the catenoid, the cone, Gott's cap, the sphere, the
            cylinder, the plane and the hyperboloid, the circumference of a rotating hole's
            horizon, and the ellipse a ring of particles is stretched into;
  marks     every curve and point marked on a surface lies on its piece;
  heights   every node of a grid against the quantity it stands for, and every triangle within
            GRID_SAG of the drawing's size of it, measured in space;
  fields    a declared star, scale factor or dust cloud against the published Einstein
            tensor it is meant to solve;
  stops     where the file says a slice cannot be drawn, g_xx - (drho/dx)^2 is negative
            there, or g_phiphi is, where the circles are timelike, and where it says a
            slice is flat, that it is.

A chord is shorter than the arc it cuts by a part in (h kappa)^2/24, h its length and
kappa the profile's curvature, so ALONG and ACROSS allow that and a little rounding; a
wrong surface misses by far more. write() refuses to write while any check fails, and
--verify prints them all with the worst error of each.


Output
------

_tools/README.md, "Embedding diagrams", is the definition of the file, since the
application builds against it: the surface itself, as profiles of surfaces of
revolution with their units, marked circles and the kind of each end, and the page's
drawing of it, projected here from a fixed camera into the form a figure in three
dimensions takes in the diagram files, with what a client needs to draw it again from
another camera, as the page does with MFS/assets/turn.js when a reader turns it.
"""

import argparse
import json
import math
import sys
import time
import warnings
from pathlib import Path

import contourpy
import numpy as np
import sympy as sp
from scipy.integrate import IntegrationWarning, quad
from sympy.utilities.lambdify import implemented_function

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))
import build_mfs_data as build  # noqa: E402
import null_rays as nr  # noqa: E402
import conformal  # noqa: E402
from conformal import Sources  # noqa: E402
from projections import Camera  # noqa: E402

EMBEDDING_DIR = build.EMBEDDING_DIR
FIELDS = ("coords", "parameters", "metric_components")
EQUATOR = {"theta": "pi/2"}

SAG = 2e-5          # how far the profile may stray from a chord, as a part of the drawing's size
BEND = 0.004        # and as a part of the chord: the chord then turns through at most 0.032
STEP = 1 / 90       # the longest chord of the profile, as a part of the drawing's size
DIGITS = 1e-7       # rho and z are rounded to below this part of their piece's extent
LORENTZ_DIGITS = 1e-9   # and a surface in Minkowski space, whose chords near the light cone are short
DPHI = 0.02         # the turn of the "across" chords, in radians
ALONG = 2e-4        # how far a chord may miss the proper distance, as a part of it
ACROSS = 2e-4       # the same for the chords across
AROUND = 1e-6       # how far rho may miss sqrt(g_phiphi), as a part of the drawing's size
FORM = 2e-5         # how far a point may lie from its closed form, as a part of the drawing's size
JOIN = 1e-9         # how far two pieces may miss each other in place and in tangent
OUTLINE = 1e-5      # how far off the surface a point of its outline is judged, as a part of the drawing's size

CAMERA = Camera(-90, 22)
RING = 720          # points round a circle of the drawing
LAB = {"lab": 21, "small": 21}  # label sizes, in units of a figure 628 wide, the most the page sets them at
GRID_SAG = 5e-4     # how far a facet of a height over a plane may lie from the height, in space, as a part of the drawing's size
HEIGHT_CAMERA = Camera(-90, 30)  # a height over a plane seen from low enough that its relief shows


# ---------------------------------------------------------------- reading a slice

class Slice:
    """The published metric of one coordinate system on the surface of x and phi, every
    other coordinate held fixed, as sympy and as numpy functions of x.

    fixed      every other coordinate, by plain name, and its value;
    params     parameter -> value; functions, a declared function -> an expression for it,
               in the plain names of the coordinates and parameters, substituted first;
    numeric    a declared function solved numerically -> (f, df), numpy functions of x giving
               it and its derivative along x on the slice, as a star's mass m(r);
    along      a coordinate that moves with x on the slice -> an expression for it in x, as
               v = T + r on Vaidya's slice of constant v - r: the metric is pulled back along
               the slice, g_xx picking up the cross terms and the moving coordinate's own;
    swept      a coordinate that moves with phi -> an expression for it in phi, as psi = -phi on
               the Mixmaster great sphere: g_phiphi picks up the moving coordinate's terms;
    turn       where the angle is no coordinate of the chart: the chart coordinate that turning
               carries x into, about the chart's origin, as y on the Malament-Hogarth plane of
               x and y; phi is then None, and the published metric is pulled back to the polar
               coordinates x cos(phi), x sin(phi) of the plane, which must leave no phi in it;
    space      "flat" for a surface drawn in flat space, or "minkowski" for one drawn in three
               dimensional Minkowski space, dX^2 + dY^2 - dZ^2, where a profile whose circles
               grow faster than the distance out to them climbs at dZ/dx = sqrt((drho/dx)^2 -
               g_xx), as the hyperbolic plane does on its hyperboloid;
    rewrite    a form to write g_phiphi and g_xx - (drho/dx)^2 in before they are evaluated, where
               sympy's own loses digits, as half_angles() for the Mixmaster great sphere.
    """

    def __init__(self, sources, metric_id, system_id, x, phi, fixed, params=None, functions=None, numeric=None,
                 along=None, swept=None, turn=None, space="flat", rewrite=None):
        _, entry, reader = nr.load(metric_id, system_id)
        sources.note(metric_id, system_id, FIELDS)
        self.metric_id, self.system_id, self.coordinate = metric_id, system_id, x
        self.lorentz = {"flat": False, "minkowski": True}[space]
        R = reader
        names = {R._plain(n): s for n, s in R.symbol.items()}
        names.update(R.parameters)
        along, swept = along or {}, swept or {}
        plane = (x, turn) if turn else (x, phi)
        held = sorted(R._plain(c) for c in entry["coords"] if c not in plane)
        if held != sorted([*fixed, *along, *swept]):
            raise SystemExit(f"{metric_id}/{system_id}: the slice of {x} and {turn or phi} holds {held} fixed, "
                             f"and the table fixes {sorted(fixed)} and moves {sorted([*along, *swept])}")
        self.x = R.symbol[x]
        self.phi = sp.Symbol("varphi", real=True) if turn else R.symbol[phi]
        subs = {R.c: 1}
        subs.update({R.parameters[k]: sp.sympify(v) for k, v in (params or {}).items()})
        held_at = {names[k]: sp.sympify(v) for k, v in fixed.items()}
        moved = {names[k]: sp.sympify(v, locals=names) for k, v in along.items()}
        held_at.update(moved)
        funcs = {R.parameters[k]: sp.sympify(v, locals=names) for k, v in (functions or {}).items()}
        self.numeric = bool(numeric)
        funcs.update({R.parameters[k]: numeric_function(k, f, df)(R.symbol[x]) for k, (f, df) in (numeric or {}).items()})

        def prep(e):
            for fn, rep in funcs.items():
                e = e.subs(fn, rep).doit()
            return sp.simplify(e.subs(subs).subs(held_at))
        g = nr.published_matrix(R, entry, "metric_components")
        coords = [names[R._plain(c)] for c in entry["coords"]]
        n = len(coords)
        if turn or swept:
            # The surface as a map from (x, phi) into the chart, pulled back whole: the tangents
            # along x and along phi are the derivatives of the map, and the metric is taken at
            # the image of each point.
            image = dict(held_at)
            image.update({names[k]: sp.sympify(v, locals=names) for k, v in swept.items()})
            if turn:
                y = R.symbol[turn]
                image.update({self.x: self.x * sp.cos(self.phi), y: self.x * sp.sin(self.phi)})
            else:
                image.setdefault(self.x, self.x)
                image.setdefault(self.phi, self.phi)
            at = [image.get(c, c) for c in coords]
            along_x = [sp.diff(c, self.x) for c in at]
            along_phi = [sp.diff(c, self.phi) for c in at]

            # The metric at the image first, and only then contracted with the tangents, whose
            # components are already functions of x and phi.
            def placed(e):
                for fn, rep in funcs.items():
                    e = e.subs(fn, rep).doit()
                return e.subs(subs).subs(dict(zip(coords, at)), simultaneous=True)
            gi = g.applyfunc(placed)

            def pulled(u, v):
                return sum(u[a] * v[b] * gi[a, b] for a in range(n) for b in range(n))
            full = [pulled(along_x, along_x), pulled(along_x, along_phi), pulled(along_phi, along_phi)]
            if turn:
                # sympy seldom clears cos^2 + sin^2 inside a declared function, so the metric is
                # read where phi = 0, on the profile, and the turn is checked to carry it onto
                # itself at a thousand points of the plane, to a part in 1e12.
                points = np.random.default_rng(0).uniform([0.01, 0], [4, 2 * math.pi], (1000, 2))
                for e in full:
                    f = sp.lambdify((self.x, self.phi), e, "numpy")
                    with np.errstate(all="ignore"):
                        turned = np.array([complex(f(r, a)) for r, a in points])
                        base = np.array([complex(f(r, 0.0)) for r, _ in points])
                    if np.nanmax(np.abs(turned - base) / np.maximum(1, np.abs(base))) > 1e-12:
                        raise SystemExit(f"{metric_id}/{system_id}: turning {x} into {turn} does not carry the "
                                         "metric onto itself, so the plane is no surface of revolution")
                full = [e.subs(self.phi, 0) for e in full]
            self.gxx, self.gxp, self.gpp = (sp.simplify(e) for e in full)
        else:
            i, j = entry["coords"].index(x), entry["coords"].index(phi)
            # The slice's tangent along x in the chart: 1 along x and d(expression)/dx along a moving
            # coordinate, which is how g_xx picks up g_xv and g_vv on a slice of constant v - r.
            tangent = [sp.Integer(0)] * n
            tangent[i] = sp.Integer(1)
            for name, expr in moved.items():
                tangent[coords.index(name)] = sp.diff(expr, self.x)
            gxx = sum(tangent[a] * tangent[b] * g[a, b] for a in range(n) for b in range(n))
            gxp = sum(tangent[a] * g[a, j] for a in range(n))
            self.gxx, self.gxp, self.gpp = prep(gxx), prep(gxp), prep(g[j, j])
        self._surface(f"{metric_id}/{system_id} on the slice of {x} and {turn or phi}", rewrite)

    def _surface(self, where, rewrite=None):
        """Check the slice is a surface of revolution, and make the functions of x it is drawn by."""
        if self.gxp != 0:
            raise SystemExit(f"{where}: g_x phi = {self.gxp}, so the slice is not a surface of revolution")
        stray = (self.gxx.free_symbols | self.gpp.free_symbols) - {self.x}
        if stray:
            raise SystemExit(f"{where}: the metric still depends on {sorted(map(str, stray))}")
        # drho/dx as g_phiphi'/(2 sqrt(g_phiphi)), so that no absolute value is differentiated.
        self.rho = sp.sqrt(self.gpp)
        self.drho = sp.diff(self.gpp, self.x) / (2 * self.rho)
        self.defect = sp.simplify(self.gxx - sp.diff(self.gpp, self.x) ** 2 / (4 * self.gpp))
        if rewrite:
            self.gpp, self.defect = rewrite(self.gpp, self.x), rewrite(self.defect, self.x)
            self.rho = sp.sqrt(self.gpp)
            self.drho = sp.diff(self.gpp, self.x) / (2 * self.rho)
        self._gxx = sp.lambdify(self.x, self.gxx, "numpy")
        self._gpp = sp.lambdify(self.x, self.gpp, "numpy")
        self._rho = sp.lambdify(self.x, self.rho, "numpy")
        self._drho = sp.lambdify(self.x, self.drho, "numpy")
        self._defect = sp.lambdify(self.x, self.defect, "numpy")
        self.exact = {}             # a float the drawing uses -> the exact number it stands for
        self._near = {}

    def horizons(self):
        """The real roots of 1/g_xx where it is a rational function, largest first, as floats,
        each remembered exactly, so that slope() takes its limit at the root itself."""
        numerator = sp.numer(sp.together(1 / self.gxx))
        solved = sp.solve(numerator, self.x)
        roots = [z for z in solved if z.is_real]
        if any(z.is_real is None for z in solved):
            # An irreducible cubic with three real roots, as Kottler's f, comes back from solve
            # in radicals of complex numbers whose reality sympy cannot decide; its real roots
            # are then taken as exact algebraic numbers.
            roots = sp.Poly(numerator, self.x).real_roots()
        roots = sorted(roots, key=float, reverse=True)
        for z in roots:
            self.exact[float(z)] = z
        return [float(z) for z in roots]

    @staticmethod
    def _at(f, x):
        # numpy evaluates both branches of a Piecewise, as E and R of the dust clouds, and the
        # one not taken may divide by zero; every value that is used is measured by the checks.
        x = np.asarray(x, dtype=float)
        with np.errstate(all="ignore"):
            return np.broadcast_to(np.asarray(f(x), dtype=float), x.shape) * 1.0

    def gxx_at(self, x):
        return self._at(self._gxx, x)

    def gpp_at(self, x):
        return self._at(self._gpp, x)

    def rho_at(self, x):
        return self._at(self._rho, x)

    def defect_at(self, x):
        """g_xx - (drho/dx)^2: the surface exists where this is not negative."""
        return self._at(self._defect, x)

    def slope(self, x, side):
        """The unit tangent (drho, dz)/sqrt(g_xx) of the profile, z rising with x, as x tends
        to `x` from above (side '+') or below ('-'), taken in sympy so that it is exact at a
        throat, where g_xx diverges, or from the numbers where a function is numerical, at a
        point where g_xx is finite."""
        if x in self.exact:
            # A horizon: 1/g_xx vanishes there while drho/dx stays finite, which is checked,
            # so drho/sqrt(g_xx) -> 0 and the unit tangent is vertical, whichever side.
            d = sp.simplify(self.drho.subs(self.x, self.exact[x]))
            if not d.is_finite:
                raise AssertionError(f"{self.metric_id}: drho/dx is {d} at the horizon {x}")
            return np.array([0.0, 1.0])
        sign = -1 if self.lorentz else 1
        if self.numeric:
            g = float(self.gxx_at(x))
            return np.array([float(self._at(self._drho, x)) / math.sqrt(g),
                             math.sqrt(max(sign * float(self.defect_at(x)), 0.0) / g)])
        x0 = self.exact.get(x, sp.nsimplify(x))
        return np.array([float(sp.limit(e, self.x, x0, side)) for e in
                         (self.drho / sp.sqrt(self.gxx), sp.sqrt(sign * self.defect / self.gxx))])

    def _local(self, root):
        """g_xx, g_phiphi and g_xx - (drho/dx)^2 as numpy functions of u = x - root, expanded in
        sympy with the root exact, so that what cancels at the root cancels exactly and they keep
        their precision where u is far below the root's own rounding, as next to a throat."""
        if root not in self._near:
            u = sp.Symbol("u")

            def near(e):
                exact = self.exact[root]
                if exact.has(sp.CRootOf):
                    # A root with no form in radicals: each coefficient in u is reduced modulo
                    # the root's minimal polynomial, so that the terms which cancel at the root
                    # are exactly zero, and only then evaluated, to forty digits.
                    a = sp.Symbol("a")
                    lowest = sp.minimal_polynomial(exact, a)
                    parts = sp.fraction(sp.together(e.subs(self.x, a + u)))
                    coeffs = [[float(sp.rem(c, lowest, a).subs(a, exact).evalf(40))
                               for c in sp.Poly(sp.expand(part), u).all_coeffs()] for part in parts]
                    return lambda v: float(np.polyval(coeffs[0], v) / np.polyval(coeffs[1], v))
                n, d = sp.fraction(sp.together(e.subs(self.x, exact + u)))
                f = sp.lambdify(u, sp.expand(n) / sp.expand(d), "numpy")
                return lambda v: float(Slice._at(f, v))
            self._near[root] = [near(e) for e in (self.gxx, self.gpp, self.defect)]
        return self._near[root]

    def _between(self, a, b):
        """(g_xx, g_phiphi, defect, lo, hi): the three as functions of s from lo to hi along the
        interval from x = a to b, s being x - root from an end at a horizon and x otherwise."""
        for end in (a, b):
            if end in self.exact:
                return (*self._local(end), a - end, b - end)
        return (lambda v: float(self.gxx_at(v)), lambda v: float(self.gpp_at(v)),
                lambda v: float(self.defect_at(v)), a, b)

    def proper(self, a, b):
        """The proper distance along the slice from x = a to x = b at fixed phi."""
        g, _, _, lo, hi = self._between(a, b)
        return integrate(lambda s: math.sqrt(g(s)), lo, hi)

    def across(self, a, b, turn):
        """The metric length of the line from x = a to b that turns steadily through
        `turn` radians per unit of proper distance."""
        g, gpp, _, lo, hi = self._between(a, b)
        return integrate(lambda s: math.sqrt(g(s) * (1 + gpp(s) * turn ** 2)), lo, hi)

    def rise(self, a, b):
        """z(b) - z(a) of the surface of revolution, the quadrature of sqrt(g_xx - (drho/dx)^2),
        or in Minkowski space of sqrt((drho/dx)^2 - g_xx)."""
        g, _, defect, lo, hi = self._between(a, b)
        sign = -1.0 if self.lorentz else 1.0

        def f(s):
            d = sign * defect(s)
            if d < -1e-12 * max(1.0, g(s)):
                space = "Minkowski space" if self.lorentz else "flat space"
                raise SystemExit(f"{self.metric_id}/{self.system_id}: g_xx - (drho/dx)^2 = {defect(s)} at "
                                 f"{self.coordinate} = {s} from {a}, where no surface of revolution in {space} "
                                 "carries the slice")
            return math.sqrt(max(d, 0.0))
        return integrate(f, lo, hi)


class FlatPlane(Slice):
    """A plane of two coordinates x and y of a chart, every other coordinate held fixed, on which
    the published metric has constant coefficients and no cross term, g_xx dx^2 + g_yy dy^2, which
    is checked. In X = sqrt(g_xx) x and Y = sqrt(g_yy) y it is the Euclidean plane, so it is drawn
    as a flat disc about the chart's origin, a surface of revolution about any of its points,
    with its profile along x: the circle through the point x of the profile has the proper radius
    sqrt(g_xx) x, so g_phiphi = g_xx x^2 with phi the angle of the Euclidean plane, which need be
    no angle of the chart, since in Kasner's plane of x and z a circle of the chart is an ellipse.
    `draw` carries chart points (x, y) to the drawing's (X, Y, 0)."""

    def __init__(self, sources, metric_id, system_id, x, y, fixed, params=None, functions=None):
        _, entry, R = nr.load(metric_id, system_id)
        sources.note(metric_id, system_id, FIELDS)
        self.metric_id, self.system_id, self.coordinate = metric_id, system_id, x
        self.lorentz, self.numeric = False, False
        names = {R._plain(n): s for n, s in R.symbol.items()}
        names.update(R.parameters)
        held = sorted(R._plain(c) for c in entry["coords"] if c not in (x, y))
        if held != sorted(fixed):
            raise SystemExit(f"{metric_id}/{system_id}: the plane of {x} and {y} holds {held} fixed, "
                             f"and the table fixes {sorted(fixed)}")
        subs = {R.c: 1}
        subs.update({R.parameters[k]: sp.sympify(v) for k, v in (params or {}).items()})
        held_at = {names[k]: sp.sympify(v) for k, v in fixed.items()}
        funcs = {R.parameters[k]: sp.sympify(v, locals=names) for k, v in (functions or {}).items()}

        def prep(e):
            for fn, rep in funcs.items():
                e = e.subs(fn, rep).doit()
            return sp.simplify(e.subs(subs).subs(held_at))
        g = nr.published_matrix(R, entry, "metric_components")
        i, j = entry["coords"].index(x), entry["coords"].index(y)
        gxx, gyy, gxy = prep(g[i, i]), prep(g[j, j]), prep(g[i, j])
        where = f"{metric_id}/{system_id} on the plane of {x} and {y}"
        if gxy != 0 or gxx.free_symbols or gyy.free_symbols or not (gxx > 0 and gyy > 0):
            raise SystemExit(f"{where}: g_xx = {gxx}, g_yy = {gyy} and g_xy = {gxy}, so the plane is not flat "
                             "in these coordinates")
        self.x, self.phi = R.symbol[x], sp.Symbol("varphi", real=True)
        self.scale = np.array([math.sqrt(float(gxx)), math.sqrt(float(gyy))])
        self.gxx, self.gxp, self.gpp = gxx, sp.Integer(0), gxx * self.x ** 2
        self._surface(where)

    def draw(self, cx, cy):
        """The drawing's (X, Y, 0) of the chart points (cx, cy) of the plane."""
        cx, cy = np.broadcast_arrays(np.asarray(cx, dtype=float), np.asarray(cy, dtype=float))
        return np.column_stack([self.scale[0] * cx.ravel(), self.scale[1] * cy.ravel(), np.zeros(cx.size)])


def half_angles(e, theta):
    """e, a function of cos(theta) and sin(theta)^2 and of multiples of theta, written as a ratio
    of polynomials in sin^2(theta/2), factored: 1 - cos(theta) = 2 sin^2(theta/2) then keeps every
    digit at the pole theta = 0, where the great sphere of the Mixmaster slice and its defect
    vanish as theta^2 and sympy's form of either subtracts numbers that agree there."""
    u = sp.Symbol("u", positive=True)
    e = sp.expand_trig(e).subs(sp.sin(theta) ** 2, 1 - sp.cos(theta) ** 2).subs(sp.cos(theta), 1 - 2 * u)
    if e.has(theta):
        raise AssertionError(f"{e} is not a function of cos(theta) and sin(theta)^2")
    # Written with float coefficients once the form is set, since exact ones may be too large
    # for numpy, and in this form nothing cancels.
    return sp.N(sp.factor(sp.cancel(e)).subs(u, sp.sin(theta / 2) ** 2), 20)


def numeric_function(name, f, df):
    """A function of one variable for sympy, which numpy evaluates as f and whose derivative is
    df, so that a declared function solved numerically enters a slice's metric and drho/dx."""
    prime = implemented_function(sp.Function(f"{name}_prime"), df)
    return type(name, (sp.Function,), {"_imp_": staticmethod(f), "fdiff": lambda self, i=1: prime(self.args[0])})


def integrate(f, a, b):
    """The integral of f from a to b, where f may diverge integrably, as 1/sqrt, at either end:
    such an end is moved to the end of a new variable s, x = end +- (b - a) s^2, which makes
    the integrand finite there."""
    if b == a:
        return 0.0
    def finite(x):
        with np.errstate(all="ignore"):
            try:
                return math.isfinite(f(x))
            except (ValueError, ZeroDivisionError, SystemExit):
                return False
    left, right = not finite(a), not finite(b)
    if left and right:
        m = 0.5 * (a + b)
        return integrate(f, a, m) + integrate(f, m, b)
    w = b - a
    # Two units in the last place from the end, so that no point the rule takes rounds onto it.
    near_a = np.nextafter(np.nextafter(a, b), b)
    near_b = np.nextafter(np.nextafter(b, a), a)
    if left:
        g = lambda s: f(max(a + w * s * s, near_a)) * 2 * w * s  # noqa: E731
    elif right:
        g = lambda s: f(min(b - w * s * s, near_b)) * 2 * w * s  # noqa: E731
    else:
        g = lambda s: f(a + w * s) * w  # noqa: E731
    # Next to a divergence the integrand is rounded, which quad reports as a warning while it
    # still returns the integral to a part in 1e-12 or better; the checks measure the result.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", IntegrationWarning)
        val, _ = quad(g, 0, 1, epsabs=1e-14, epsrel=1e-11, limit=400)
    return val


# ---------------------------------------------------------------- the surface

class Piece:
    """One run of a profile: x from `lo` to `hi` on one slice, the surface rising with x
    (sense +1) or falling (sense -1) from z0 at lo.

    cls        what it is, which names its tint and the style of its circles;
    ends       (start, end), each (kind, text): axis, apex, join, throat, edge or stops;
    marks      values of x whose circles the file marks, each (x, class, TeX label or None);
    reference  a surface drawn for comparison that is not part of the slice.
    """

    def __init__(self, pid, cls, sl, lo, hi, z0=0.0, sense=1, ends=(("edge", None), ("edge", None)),
                 marks=(), size=1.0, reference=False, legend=None, digits=None):
        self.id, self.cls, self.sl, self.size, self.digits = pid, cls, sl, size, digits
        self.lo, self.hi, self.sense = lo, hi, sense
        self.ends, self.marks, self.reference, self.legend = ends, list(marks), reference, legend
        knots = sorted({lo, hi} | {m[0] for m in self.marks if lo <= m[0] <= hi})
        xs, zs = [knots[0]], [z0]
        for a, b in zip(knots, knots[1:]):
            self._refine(a, zs[-1], b, zs[-1] + sense * sl.rise(a, b), size, xs, zs)
        self.x = np.array(xs)
        self.z = np.array(zs)
        self.rho = sl.rho_at(self.x)

    def _refine(self, a, za, b, zb, size, xs, zs, depth=0):
        """Append the points after a up to b, halving until the chord follows the profile."""
        m = 0.5 * (a + b)
        zm = za + self.sense * self.sl.rise(a, m)
        pa = np.array([float(self.sl.rho_at(a)), za])
        pb = np.array([float(self.sl.rho_at(b)), zb])
        pm = np.array([float(self.sl.rho_at(m)), zm])
        chord = pb - pa
        length = float(np.hypot(*chord))
        off = pm - pa
        stray = (abs(float(chord[0] * off[1] - chord[1] * off[0])) / length if length > 0
                 else float(np.hypot(*off)))
        if depth < 40 and (stray > SAG * size or stray > BEND * length or length > STEP * size):
            self._refine(a, za, m, zm, size, xs, zs, depth + 1)
            self._refine(m, zm, b, zb, size, xs, zs, depth + 1)
            return
        xs.append(b)
        zs.append(zb)

    @property
    def decimals(self):
        return decimals(max(float(np.ptp(self.rho)), float(np.ptp(self.z)), 1e-9),
                        self.digits or (LORENTZ_DIGITS if self.sl.lorentz else DIGITS))

    def at(self, x):
        """(rho, z) of the circle at x, which must be a point of the profile."""
        i = int(np.argmin(np.abs(self.x - x)))
        if abs(self.x[i] - x) > 1e-12:
            raise AssertionError(f"piece {self.id}: {x} is not a point of the profile")
        return float(self.rho[i]), float(self.z[i])

    def inward(self, x):
        """The unit tangent of the profile at the end x, pointing into the piece."""
        if x == self.lo:
            return self.sl.slope(x, "+") * [1, self.sense]
        if x == self.hi:
            return -self.sl.slope(x, "-") * [1, self.sense]
        raise AssertionError(f"piece {self.id}: {x} is not an end")

    def data(self):
        out = {"id": self.id, "class": self.cls, "metric": self.sl.metric_id, "system": self.sl.system_id,
               "coordinate": self.sl.coordinate,
               "points": [[significant(x), fixed(r, self.decimals), fixed(z, self.decimals)]
                          for x, r, z in zip(self.x, self.rho, self.z)],
               "start": end_data(self.ends[0]), "end": end_data(self.ends[1])}
        if self.sl.lorentz:
            out["space"] = "minkowski"
        if self.reference:
            out["reference"] = True
        return out


class FormPiece(Piece):
    """A reference piece given by a closed form rather than built from a slice's metric, as the
    light cone anti-de Sitter's hyperboloid nears, which is no slice of the spacetime: its points
    are rho_of(x) and z_of(x) at the values xs of the coordinate of `sl` it is set beside."""

    def __init__(self, pid, sl, xs, rho_of, z_of, ends, size=1.0):
        self.id, self.cls, self.sl, self.size, self.digits = pid, "reference", sl, size, None
        self.lo, self.hi, self.sense = float(xs[0]), float(xs[-1]), 1
        self.ends, self.marks, self.reference, self.legend = ends, [], True, None
        self.x = np.asarray(xs, dtype=float)
        self.rho, self.z = np.asarray(rho_of(self.x), dtype=float), np.asarray(z_of(self.x), dtype=float)


# Adding 0.0 turns a -0.0 that rounding leaves into 0.0, which is how every number is written.
def fixed(v, digits):
    return round(float(v), digits) + 0.0


def decimals(extent, digits=DIGITS):
    """The decimals written for rho and z of a piece whose rho or z spans `extent`, at least
    six, and more for a small piece, as a collapsing star's shrinking cap, whose chords are
    short, so the rounding is below `digits` of it."""
    return max(6, math.ceil(-math.log10(digits * extent)))


# x is written as the double itself, the shortest decimal that reads back as it: next to a
# throat g_xx diverges, and an irrational horizon, as Kerr's, rounded to fewer figures would
# fall inside it.
def significant(x):
    return float(x) + 0.0


def end_data(end):
    kind, text = end
    return {"kind": kind, "text": text} if text else {"kind": kind}


class Curve:
    """A curve marked on a piece of a surface that is no circle about its axis, as a ring of free
    particles stretched into an ellipse or a line of flow: points (X, Y, Z) in the surface's own
    frame, the axis along Z, closed or open."""

    def __init__(self, piece, cls, points, closed=False, label=None, time=None):
        self.piece, self.cls, self.closed = piece, cls, closed
        self.points = np.asarray(points, dtype=float)
        self.label, self.time = label, time

    def data(self):
        d = self.piece.decimals
        out = {"piece": self.piece.id, "class": self.cls}
        if self.label:
            out["label"] = self.label
            out["time"] = fixed(self.time, 6)
        if self.closed:
            out["closed"] = True
        out["points"] = [[fixed(v, d) for v in P] for P in self.points]
        return out


class Surface:
    """A slice as a surface of revolution about the axis z: pieces and marked circles, and the
    curves and points marked on it that are no circles. In a sequence, one surface per moment,
    with its `label` and `time`."""

    def __init__(self, pieces, label=None, time=None, curves=(), dots=(), axis=None):
        self.pieces, self.label, self.time = pieces, label, time
        self.curves = list(curves)
        self.dots = list(dots)          # (piece, class, (X, Y, Z))
        self.axis = axis                # (z from, z to): the axis of time a stack of moments stands on

    def rings(self):
        out = []
        for p in self.pieces:
            for x, cls, text in p.marks:
                rho, z = p.at(x)
                d = p.decimals
                ring = {"piece": p.id, "class": cls, "x": significant(x), "rho": fixed(rho, d), "z": fixed(z, d)}
                if text:
                    ring["label"] = text
                out.append(ring)
        return out

    def data(self):
        out = {}
        if self.label:
            out["label"] = self.label
            if self.time is not None:
                out["time"] = fixed(self.time, 6)
        out["pieces"] = [p.data() for p in self.pieces]
        out["rings"] = self.rings()
        if self.curves:
            out["curves"] = [c.data() for c in self.curves]
        if self.dots:
            out["dots"] = [{"piece": p.id, "class": cls, "at": [fixed(v, p.decimals) for v in P]}
                           for p, cls, P in self.dots]
        if self.axis:
            d = max(p.decimals for p in self.pieces)
            out["axis"] = {"from": fixed(self.axis[0], d), "to": fixed(self.axis[1], d)}
        return out


# ---------------------------------------------------------------- a height over a plane

class Facets:
    """Flat triangles, corners A, B and C of shape (T, 3), seen along one direction: for points of
    the page, which triangles lie over them and how near the viewer, found exactly. Seen along
    `toward`, a ray from a point P toward the viewer meets a triangle exactly where the triangle
    covers the point of the page P falls on, and there it meets it after the difference of their
    depths, so each point is tested against the triangles whose box on the page holds it, found
    through a grid of bins laid over the page."""

    def __init__(self, A, B, C):
        self.A, self.B, self.C = (np.asarray(M, dtype=float) for M in (A, B, C))
        self._seen = {}

    def _frame(self, right, up, toward):
        key = tuple(np.round(np.concatenate([right, up, toward]), 15))
        if key in self._seen:
            return self._seen[key]
        page = [np.stack([M @ right, M @ up], -1) for M in (self.A, self.B, self.C)]
        depth = [M @ toward for M in (self.A, self.B, self.C)]
        (ax, ay), (bx, by), (cx, cy) = ((p[:, 0], p[:, 1]) for p in page)
        area = (bx - ax) * (cy - ay) - (cx - ax) * (by - ay)
        # A triangle seen edge on covers no area of the page and hides nothing.
        keep = np.abs(area) > 1e-15 * max(1.0, float(np.max(np.abs(area))))
        lo = np.minimum(np.minimum(page[0], page[1]), page[2])
        hi = np.maximum(np.maximum(page[0], page[1]), page[2])
        idx = np.flatnonzero(keep)
        n = max(1, int(math.sqrt(max(len(idx), 1))))
        g0 = lo[idx].min(0) if len(idx) else np.zeros(2)
        g1 = hi[idx].max(0) if len(idx) else np.ones(2)
        step = np.maximum((g1 - g0) / n, 1e-12)
        i0 = np.clip(((lo[idx] - g0) // step).astype(int), 0, n - 1)
        i1 = np.clip(((hi[idx] - g0) // step).astype(int), 0, n - 1)
        nx, ny = i1[:, 0] - i0[:, 0] + 1, i1[:, 1] - i0[:, 1] + 1
        count = nx * ny
        tri = np.repeat(idx, count)
        k = np.arange(int(count.sum())) - np.repeat(np.cumsum(count) - count, count)
        bx_ = np.repeat(i0[:, 0], count) + k % np.repeat(nx, count)
        by_ = np.repeat(i0[:, 1], count) + k // np.repeat(nx, count)
        bins = by_ * n + bx_
        order = np.argsort(bins, kind="stable")
        tri, bins = tri[order], bins[order]
        start = np.searchsorted(bins, np.arange(n * n))
        end = np.searchsorted(bins, np.arange(n * n), side="right")
        frame = {"page": page, "depth": depth, "area": area, "n": n, "g0": g0, "g1": g1, "step": step,
                 "tri": tri, "start": start, "end": end}
        self._seen[key] = frame
        return frame

    def over(self, S, right, up, toward):
        """For points S of the page, (point, triangle, depth of the triangle there) for every
        triangle that covers the point, edges included."""
        F = self._frame(np.asarray(right, dtype=float), np.asarray(up, dtype=float), np.asarray(toward, dtype=float))
        S = np.asarray(S, dtype=float).reshape(-1, 2)
        n, g0, g1, step = F["n"], F["g0"], F["g1"], F["step"]
        # A point on the far edge of the bins falls in the last of them, as a triangle's box does.
        cell = np.clip(np.floor((S - g0) / step).astype(int), 0, n - 1)
        inside = np.all((S >= g0) & (S <= g1), axis=1)
        b = np.where(inside, cell[:, 1] * n + cell[:, 0], 0)
        count = np.where(inside, F["end"][b] - F["start"][b], 0)
        q = np.repeat(np.arange(len(S)), count)
        k = np.arange(int(count.sum())) - np.repeat(np.cumsum(count) - count, count)
        t = F["tri"][np.repeat(F["start"][b], count) + k]
        (ax, ay), (bx, by), (cx, cy) = ((p[t, 0], p[t, 1]) for p in F["page"])
        px, py = S[q, 0], S[q, 1]
        area = F["area"][t]
        wa = ((bx - px) * (cy - py) - (cx - px) * (by - py)) / area
        wb = ((cx - px) * (ay - py) - (ax - px) * (cy - py)) / area
        wc = 1.0 - wa - wb
        hit = (wa >= -1e-12) & (wb >= -1e-12) & (wc >= -1e-12)
        d = wa * F["depth"][0][t] + wb * F["depth"][1][t] + wc * F["depth"][2][t]
        return q[hit], t[hit], d[hit]

    def reach(self, P, camera, eps):
        """For rays from P toward the camera, the largest t > eps at which each meets a triangle,
        or -inf."""
        P = np.asarray(P, dtype=float)
        q, _, d = self.over(camera.screen(P), camera.right, camera.up, camera.toward)
        t = d - P[q] @ camera.toward
        out = np.full(len(P), -np.inf)
        ok = t > eps
        np.maximum.at(out, q[ok], t[ok])
        return out

    def height(self, X, Y):
        """The height of the triangles over the points (X, Y), seen from straight above, or nan
        where none lies."""
        X, Y = np.broadcast_arrays(np.asarray(X, dtype=float), np.asarray(Y, dtype=float))
        q, _, d = self.over(np.column_stack([X.ravel(), Y.ravel()]), [1.0, 0.0, 0.0], [0.0, 1.0, 0.0],
                            [0.0, 0.0, 1.0])
        out = np.full(X.size, -np.inf)
        np.maximum.at(out, q, d)
        return np.where(np.isfinite(out), out, np.nan).reshape(X.shape)


class GridPiece:
    """A quantity drawn as a height over a plane, z(X, Y) in the surface's own frame, sampled on a
    grid and drawn as flat triangles. It is no slice of the spacetime, and nothing is measured
    along it: the height stands for the quantity, at the scale the view states.

    frame      "polar": u is the distance from the axis and v the angle from X toward Y, running
               round; "cartesian": u is X and v is Y;
    u, v       the values the grid must hold, as the lines it draws and the curves it marks along
               them; u is halved, and for a Cartesian grid v as well, until every triangle lies
               within GRID_SAG of the drawing's size of the height, measured in space;
    height     z(X, Y), vectorised;
    edge       the TeX text for what the surface does beyond its rim;
    lines      the grid lines the figure draws, {"u": [values], "v": [values]}.

    Each cell between neighbouring u_i, u_(i+1) and v_j, v_(j+1) is cut into two triangles
    along its diagonal from (i, j) to (i+1, j+1). The heights are rounded as they are written,
    and every triangle and check is taken from the rounded heights, so the surface drawn is the
    surface the file holds."""

    def __init__(self, pid, cls, frame, u, v, height, metric_id, system_id, edge, size, lines):
        self.id, self.cls, self.frame, self.size = pid, cls, frame, size
        self.metric_id, self.system_id, self.edge = metric_id, system_id, edge
        self.reference, self.marks = False, []
        self.wrap = frame == "polar"
        self.fun = height
        self.u = np.array(sorted(set(float(x) for x in u)))
        self.v = np.array(sorted(set(float(x) for x in v)))
        while True:
            err = self._errors()
            bad_u = np.flatnonzero(err.max(axis=1) > GRID_SAG * size)
            bad_v = np.flatnonzero(err.max(axis=0) > GRID_SAG * size) if not self.wrap else np.array([], int)
            if not len(bad_u) and not len(bad_v):
                break
            # Halving u cannot bring a polar grid's cells nearer the height than its angles allow,
            # and the halving then runs on without end.
            if len(self.u) * len(self.v) > 200000:
                raise AssertionError(f"piece {pid}: the grid cannot hold the height within GRID_SAG by halving; "
                                     "a polar grid needs more angles")
            self.u = np.sort(np.concatenate([self.u, 0.5 * (self.u[bad_u] + self.u[bad_u + 1])]))
            if len(bad_v):
                self.v = np.sort(np.concatenate([self.v, 0.5 * (self.v[bad_v] + self.v[bad_v + 1])]))
        X, Y = self.plane()
        Z = height(X, Y)
        self.decimals = decimals(max(float(np.ptp(X)), float(np.ptp(Y)), float(np.ptp(Z)), 1e-9))
        self.Z = np.round(Z, self.decimals) + 0.0
        self.lines = {k: [self.index(k, x) for x in lines.get(k, [])] for k in ("u", "v")}
        self.facets = Facets(*self.triangles())
        # The width of the grid across the plane, which the tests of what hides a point on it
        # are measured in, so that the page, which knows the grid and not the drawing's size,
        # judges them as the figure is written.
        self.extent = self.facets.extent = max(float(np.ptp(X)), float(np.ptp(Y)))

    def plane(self, u=None, v=None):
        u = self.u if u is None else u
        v = self.v if v is None else v
        if self.wrap:
            return np.outer(u, np.cos(v)), np.outer(u, np.sin(v))
        return np.meshgrid(u, v, indexing="ij")

    def index(self, which, value):
        grid = self.u if which == "u" else self.v
        i = int(np.argmin(np.abs(grid - value)))
        if abs(grid[i] - value) > 1e-12:
            raise AssertionError(f"piece {self.id}: {value} is not a value of {which} on the grid")
        return i

    def nodes(self):
        X, Y = self.plane()
        return np.stack([X, Y, self.Z], -1)

    def _corners(self, P):
        if self.wrap:
            P = np.concatenate([P, P[:, :1]], axis=1)
        a, b, c, d = P[:-1, :-1], P[1:, :-1], P[1:, 1:], P[:-1, 1:]
        return a, b, c, d

    def triangles(self):
        a, b, c, d = (M.reshape(-1, 3) for M in self._corners(self.nodes()))
        return np.vstack([a, a]), np.vstack([b, c]), np.vstack([c, d])

    def _errors(self, Z=None):
        """The distance in space from each cell's triangles to the height, the worst of samples
        over each, with the heights Z at the corners, unrounded unless given."""
        X, Y = self.plane()
        P = np.stack([X, Y, self.fun(X, Y) if Z is None else Z], -1)
        a, b, c, d = self._corners(P)
        g = np.linspace(0, 1, 9)
        s, t = (w.ravel() for w in np.meshgrid(g, g))
        keep = s + t <= 1 + 1e-12
        s, t = s[keep], t[keep]
        worst = np.zeros(a.shape[:2])
        h = 1e-6
        for A, B, C in ((a, b, c), (a, c, d)):
            Q = A[..., None, :] + s[:, None] * (B - A)[..., None, :] + t[:, None] * (C - A)[..., None, :]
            x, y = Q[..., 0], Q[..., 1]
            gx = (self.fun(x + h, y) - self.fun(x - h, y)) / (2 * h)
            gy = (self.fun(x, y + h) - self.fun(x, y - h)) / (2 * h)
            miss = np.abs(Q[..., 2] - self.fun(x, y)) / np.sqrt(1 + gx * gx + gy * gy)
            worst = np.maximum(worst, miss.max(axis=-1))
        return worst

    def sag(self):
        """The worst distance in space from a triangle, as the file holds it, to the height."""
        return float(self._errors(self.Z).max())

    def rim(self):
        """The edge of the grid, its nodes in order round it, to be closed back to the first."""
        P = self.nodes()
        if self.wrap:
            return P[-1]
        return np.vstack([P[0, :-1], P[:-1, -1], P[-1, :0:-1], P[:0:-1, 0]])

    def data(self):
        d = self.decimals
        return {"id": self.id, "class": self.cls, "metric": self.metric_id, "system": self.system_id,
                "grid": {"frame": self.frame, "u": [significant(x) for x in self.u],
                         "v": [significant(x) for x in self.v],
                         "z": [[fixed(z, d) for z in row] for row in self.Z]},
                "edge": end_data(("edge", self.edge))}


class EllipsesPiece(GridPiece):
    """A surface drawn as a stack of ellipses about a vertical axis, the row at u the ellipse
    (a(u) cos v, b(u) sin v) at the height z(u), v the angle of its columns, sampled on a grid
    and drawn as flat triangles as a height over a plane is. It carries two kinds of surface:

      a tube, the ellipses of one ring of free particles at every moment between its first and
      its last, each at the height of its time, so that the height is time and nothing is
      measured along it but round each row, one moment of a flat plane;
      a cone cut along one meridian and opened, a = b = r sin(alpha) and z = r cos(alpha) with r
      the distance from its apex, running round only part of the axis, `open`, as a cone unrolls
      onto the plane, every distance along it the distance on the cone.

    rows(u) gives (a, b, z) for an array of u; u is halved until every triangle lies within
    GRID_SAG of the drawing's size of the surface rows() describes, measured in space, which is
    the parametric miss at points spread over each triangle, at least the distance to the
    surface. `lines` names the lines drawn on it, [(class, "u" or "v", values)], each a value
    of the grid. The heights and axes are rounded as they are written, and every triangle and
    check is taken from the rounded numbers."""

    def __init__(self, pid, cls, u, v, rows, metric_id, system_id, edge, size, lines, open_=False):
        self.id, self.cls, self.frame, self.size = pid, cls, "ellipses", size
        self.metric_id, self.system_id, self.edge = metric_id, system_id, edge
        self.reference, self.marks = False, []
        self.open = open_
        self.wrap = not open_
        self.rows = rows
        self.u = np.array(sorted(set(float(x) for x in u)))
        self.v = np.array(sorted(set(float(x) for x in v)))
        while True:
            err = self._errors()
            bad = np.flatnonzero(err.max(axis=1) > GRID_SAG * size)
            if not len(bad):
                break
            if len(self.u) > 20000:
                raise AssertionError(f"piece {pid}: the rows cannot hold the surface within GRID_SAG by halving")
            self.u = np.sort(np.concatenate([self.u, 0.5 * (self.u[bad] + self.u[bad + 1])]))
        a, b, z = (np.asarray(w, dtype=float) for w in rows(self.u))
        self.decimals = decimals(max(2 * float(np.max(a)), 2 * float(np.max(b)), float(np.ptp(z)), 1e-9))
        self.A, self.B, self.Zr = (np.round(w, self.decimals) + 0.0 for w in (a, b, z))
        self.Z = np.repeat(self.Zr[:, None], len(self.v), axis=1)
        self.drawn = [(c, which, [self.index(which, x) for x in values]) for c, which, values in lines]
        self.lines = {"u": [], "v": []}
        self.facets = Facets(*self.triangles())
        X, Y = self.plane()
        self.extent = self.facets.extent = max(float(np.ptp(X)), float(np.ptp(Y)))

    def plane(self, u=None, v=None):
        if u is None and v is None:
            return np.outer(self.A, np.cos(self.v)), np.outer(self.B, np.sin(self.v))
        a, b, _ = self.rows(self.u if u is None else u)
        v = self.v if v is None else v
        return np.outer(a, np.cos(v)), np.outer(b, np.sin(v))

    def _errors(self, Z=None):
        """The distance in space from each cell's triangles to the surface at the same place of
        the grid, the worst of samples over each, from the unrounded rows."""
        u, v = self.u, self.v
        vv = np.append(v, v[0] + 2 * math.pi) if self.wrap else v
        a, b, z = (np.asarray(w, dtype=float) for w in self.rows(u))
        P = np.stack([np.outer(a, np.cos(vv)), np.outer(b, np.sin(vv)), np.repeat(z[:, None], len(vv), axis=1)], -1)
        g = np.linspace(0, 1, 9)
        s, t = (w.ravel() for w in np.meshgrid(g, g))
        keep = s + t <= 1 + 1e-12
        s, t = s[keep], t[keep]
        du, dv = np.diff(u), np.diff(vv)
        worst = np.zeros((len(u) - 1, len(vv) - 1))
        A, B, C, D = P[:-1, :-1], P[1:, :-1], P[1:, 1:], P[:-1, 1:]
        for (Pa, Pb, Pc), (fu, fv) in (((A, B, C), (s + t, t)), ((A, C, D), (s, s + t))):
            Q = Pa[..., None, :] + s[:, None] * (Pb - Pa)[..., None, :] + t[:, None] * (Pc - Pa)[..., None, :]
            uu = u[:-1, None, None] + du[:, None, None] * fu[None, None, :]
            ww = np.broadcast_to(vv[None, :-1, None] + dv[None, :, None] * fv[None, None, :],
                                 (len(u) - 1, len(vv) - 1, len(fv)))
            ra, rb, rz = (np.asarray(w, dtype=float).reshape(uu.shape[0], 1, -1) for w in self.rows(uu[:, 0, :].ravel()))
            T = np.stack([ra * np.cos(ww), rb * np.sin(ww), np.broadcast_to(rz, ww.shape)], -1)
            worst = np.maximum(worst, np.linalg.norm(Q - T, axis=-1).max(axis=-1))
        return worst

    def sag(self):
        return float(self._errors().max())

    def nodes(self):
        X, Y = self.plane()
        return np.stack([X, Y, self.Z], -1)

    def normals(self):
        """The unit normal at each node, pointing away from the axis: the sum of the normals of
        the triangles that meet there, each as long as twice its area, taken from the grid's own
        order, rows up and columns round, which turns every triangle the same way."""
        P = self.nodes()
        m, n = P.shape[:2]
        Q = np.concatenate([P, P[:, :1]], axis=1) if self.wrap else P
        a, b, c, d = Q[:-1, :-1], Q[1:, :-1], Q[1:, 1:], Q[:-1, 1:]
        N = np.zeros((m, Q.shape[1], 3))
        for A, B, C, corners in ((a, b, c, ((0, 0), (1, 0), (1, 1))), (a, c, d, ((0, 0), (1, 1), (0, 1)))):
            f = -np.cross(B - A, C - A)
            for di, dj in corners:
                N[di:di + m - 1, dj:dj + Q.shape[1] - 1] += f
        if self.wrap:
            N[:, 0] += N[:, -1]
            N = N[:, :-1]
        # A node that no triangle with any area meets, as a cone's apex at the end of its cut,
        # takes the normal of the node before it in its row, or after it at the row's start.
        L = np.linalg.norm(N, axis=-1)
        for i, j in zip(*np.nonzero(L == 0)):
            N[i, j] = N[i, j - 1] if j > 0 else N[i, j + 1]
        return N / np.linalg.norm(N, axis=-1, keepdims=True)

    def rim(self):
        raise AssertionError(f"piece {self.id}: a stack of ellipses draws only the lines it names")

    def data(self):
        d = self.decimals
        grid = {"frame": "ellipses", "u": [significant(x) for x in self.u], "v": [significant(x) for x in self.v],
                "a": [fixed(x, d) for x in self.A], "b": [fixed(x, d) for x in self.B], "z": [fixed(x, d) for x in self.Zr]}
        if self.open:
            grid["open"] = True
        grid["lines"] = [{"class": c, "u": idx if which == "u" else [], "v": idx if which == "v" else []}
                         for c, which, idx in self.drawn]
        return {"id": self.id, "class": self.cls, "metric": self.metric_id, "system": self.system_id,
                "grid": grid, "edge": end_data(("edge", self.edge))}


def sides(piece, P, size):
    """The direction each point of a line on a grid piece is judged along, as long as the piece's
    extent is a part of the drawing's size: up, for a height over a plane, whose two sides are
    above and below it; and straight out from the axis for a stack of ellipses, which that
    crosses wherever it is not flat, or up on the axis itself."""
    P = np.asarray(P, dtype=float)
    k = piece.extent / size
    if not isinstance(piece, EllipsesPiece):
        return np.tile([0.0, 0.0, k], (len(P), 1))
    r = np.hypot(P[:, 0], P[:, 1])
    out = np.zeros((len(P), 3))
    on = r > 0
    out[on, 0], out[on, 1] = P[on, 0] / r[on], P[on, 1] / r[on]
    out[~on, 2] = 1.0
    return out * k


def level_curves(piece, level):
    """The curves on which the triangles of a grid piece stand at the height `level`, each a list of
    points on the edges of the triangles, closed where it closes, a corner at the level counted
    as above it."""
    A, B, C = piece.triangles()
    m, n = len(piece.u), len(piece.v)
    cols = n if piece.wrap else n - 1
    # The node index of each corner, to name the edges two triangles share.
    ii, jj = np.meshgrid(np.arange(m - 1), np.arange(cols), indexing="ij")
    na = (ii * n + jj).ravel()
    nb = ((ii + 1) * n + jj).ravel()
    nc = ((ii + 1) * n + (jj + 1) % n).ravel()
    nd = (ii * n + (jj + 1) % n).ravel()
    corners = [np.concatenate([na, na]), np.concatenate([nb, nc]), np.concatenate([nc, nd])]
    pts = [A, B, C]
    links, at = {}, {}
    for k in range(len(A)):
        above = [pts[e][k][2] >= level for e in range(3)]
        if all(above) or not any(above):
            continue
        ends = []
        for e0, e1 in ((0, 1), (1, 2), (2, 0)):
            if above[e0] != above[e1]:
                key = tuple(sorted((int(corners[e0][k]), int(corners[e1][k]))))
                if key[0] == key[1]:
                    continue
                if key not in at:
                    P, Q = pts[e0][k], pts[e1][k]
                    w = (level - P[2]) / (Q[2] - P[2])
                    at[key] = P + w * (Q - P)
                ends.append(key)
        if len(ends) == 2:
            links.setdefault(ends[0], []).append(ends[1])
            links.setdefault(ends[1], []).append(ends[0])
    out, done = [], set()

    def walk(start):
        line, prev, cur = [start], None, start
        done.add(start)
        while True:
            nxt = [e for e in links[cur] if e != prev and e not in done]
            if not nxt:
                break
            prev, cur = cur, nxt[0]
            line.append(cur)
            done.add(cur)
        closed = len(line) > 2 and start in links[cur]
        return np.array([at[e] for e in line]), closed
    for key, nb_ in links.items():
        if len(nb_) == 1 and key not in done:
            out.append(walk(key))
    for key in links:
        if key not in done:
            out.append(walk(key))
    return out


def grid_isolines(F, ni, nj):
    """The polylines where F = 0 on a grid of ni rows and nj columns, F[i * nj + j], by marching
    squares with a zero counted as positive and a saddle settled by the mean of its corners, each
    a list of (i, j) in grid units, a closed one repeating its first point at its end: the
    isolines() of MFS/assets/turn.js, step for step, so the page draws the outline of a
    grid as the figure is written."""
    H = ni * nj
    at, links = {}, {}

    def place(key, i0, j0, i1, j1, p, q):
        if key not in at:
            w = p / (p - q)
            at[key] = (i0 + w * (i1 - i0), j0 + w * (j1 - j0))
        return key

    def link(a, b):
        links.setdefault(a, []).append(b)
        links.setdefault(b, []).append(a)
    for i in range(ni - 1):
        for j in range(nj - 1):
            a, b, c, d = F[i * nj + j], F[i * nj + j + 1], F[(i + 1) * nj + j + 1], F[(i + 1) * nj + j]
            ia, ib, ic, id_ = a >= 0, b >= 0, c >= 0, d >= 0
            if ia == ib == ic == id_:
                continue
            T = place(i * nj + j, i, j, i, j + 1, a, b) if ia != ib else -1
            R = place(H + i * nj + j + 1, i, j + 1, i + 1, j + 1, b, c) if ib != ic else -1
            B = place((i + 1) * nj + j, i + 1, j, i + 1, j + 1, d, c) if id_ != ic else -1
            L = place(H + i * nj + j, i, j, i + 1, j, a, d) if ia != id_ else -1
            if T >= 0 and R >= 0 and B >= 0 and L >= 0:
                if (a + b + c + d >= 0) == ia:
                    link(T, R)
                    link(B, L)
                else:
                    link(T, L)
                    link(R, B)
            else:
                ends = [e for e in (T, R, B, L) if e >= 0]
                link(ends[0], ends[1])
    out, done = [], set()

    def walk(start):
        line, prev, cur = [start], -1, start
        done.add(start)
        while True:
            nxt = [e for e in links[cur] if e != prev and e not in done]
            if not nxt:
                break
            line.append(nxt[0])
            done.add(nxt[0])
            prev, cur = cur, nxt[0]
        if len(line) > 2 and start in links[cur]:
            line.append(start)
        return [at[k] for k in line]
    for key, nb_ in links.items():
        if len(nb_) == 1 and key not in done:
            out.append(walk(key))
    for key in links:
        if key not in done:
            out.append(walk(key))
    return out


def grid_normals(piece):
    """The unit normal at each node of a grid piece, pointing up: the sum of the normals of the
    triangles that meet there, each as long as twice its area."""
    P = piece.nodes()
    m, n = P.shape[:2]
    Q = np.concatenate([P, P[:, :1]], axis=1) if piece.wrap else P
    a, b, c, d = Q[:-1, :-1], Q[1:, :-1], Q[1:, 1:], Q[:-1, 1:]
    N = np.zeros((m, Q.shape[1], 3))
    for A, B, C, corners in ((a, b, c, ((0, 0), (1, 0), (1, 1))), (a, c, d, ((0, 0), (1, 1), (0, 1)))):
        f = np.cross(B - A, C - A)
        f *= np.where(f[..., 2:3] < 0, -1.0, 1.0)
        for di, dj in corners:
            N[di:di + m - 1, dj:dj + Q.shape[1] - 1] += f
    if piece.wrap:
        N[:, 0] += N[:, -1]
        N = N[:, :-1]
    return N / np.linalg.norm(N, axis=-1, keepdims=True)


def grid_outline(piece, cam):
    """Where a grid piece turns edge on to the camera: the lines where the normal at the nodes,
    carried linearly along each edge of the grid, is square to the line of sight, found by
    marching squares over the grid, each point on an edge of the grid and so on the triangles."""
    P = piece.nodes()
    N = piece.normals() if isinstance(piece, EllipsesPiece) else grid_normals(piece)
    F = N @ cam.toward
    if piece.wrap:
        P = np.concatenate([P, P[:, :1]], axis=1)
        F = np.concatenate([F, F[:, :1]], axis=1)
    ni, nj = F.shape
    runs = []
    for line in grid_isolines(F.ravel(), ni, nj):
        pts = []
        for gi, gj in line:
            i0, j0 = int(math.floor(gi)), int(math.floor(gj))
            i1, j1 = min(i0 + 1, ni - 1), min(j0 + 1, nj - 1)
            wi, wj = gi - i0, gj - j0
            pts.append((1 - wi) * (1 - wj) * P[i0, j0] + wi * (1 - wj) * P[i1, j0]
                       + (1 - wi) * wj * P[i0, j1] + wi * wj * P[i1, j1])
        runs.append(np.array(pts))
    return runs


# ---------------------------------------------------------------- the checks

class Checks:
    """Every surface measured against the metric it was drawn from, as it is written."""

    def __init__(self):
        self.items = []

    def add(self, name, error, tol):
        self.items.append({"name": name, "error": float(error), "tol": tol, "ok": bool(error <= tol)})

    def isometry(self, where, piece):
        """along, across and around, from the rounded numbers the file holds."""
        self.size = piece.size
        pts = np.array(piece.data()["points"])
        x, rho, z = pts[:, 0], pts[:, 1], pts[:, 2]
        sl = piece.sl
        # In Minkowski space a chord's length is sqrt(drho^2 - dz^2), and every chord of a
        # surface drawn there is spacelike.
        chords = (np.sqrt(np.diff(rho) ** 2 - np.diff(z) ** 2) if sl.lorentz
                  else np.hypot(np.diff(rho), np.diff(z)))
        proper = np.array([sl.proper(a, b) for a, b in zip(x, x[1:])])
        self.add(f"{where}: along each chord", float(np.max(np.abs(chords - proper) / proper)), ALONG)
        self.add(f"{where}: along the whole profile", abs(chords.sum() - proper.sum()) / proper.sum(), ALONG)
        # Across: the line that runs out at a steady proper distance while it turns steadily
        # through DPHI, whose length is the quadrature of sqrt(g_xx (1 + g_phiphi (DPHI/L)^2)),
        # L the proper distance out, against the straight line in space between its ends.
        worst = 0.0
        for i in range(len(x) - 1):
            a, b, L = x[i], x[i + 1], proper[i]
            P = np.array([rho[i], 0.0, z[i]])
            Q = np.array([rho[i + 1] * math.cos(DPHI), rho[i + 1] * math.sin(DPHI), z[i + 1]])
            d = Q - P
            space = float(math.sqrt(d[0] ** 2 + d[1] ** 2 - d[2] ** 2) if sl.lorentz else np.linalg.norm(d))
            length = sl.across(a, b, DPHI / L)
            worst = max(worst, abs(space - length) / length)
        self.add(f"{where}: across, {DPHI} round the axis", worst, ACROSS)
        want = np.sqrt(np.maximum(sl.gpp_at(x), 0.0))
        self.add(f"{where}: around, rho = sqrt(g_phiphi)", float(np.max(np.abs(rho - want))) / self.size, AROUND)

    def join(self, where, a, x_a, b, x_b):
        """Piece a at x_a meets piece b at x_b, each at one of its ends, in one point, and the
        tangent pointing into a is the opposite of the tangent pointing into b."""
        pa, pb = np.array(a.at(x_a)), np.array(b.at(x_b))
        self.add(f"{where}: one point", float(np.max(np.abs(pa - pb))), JOIN)
        self.add(f"{where}: one tangent", float(np.max(np.abs(a.inward(x_a) + b.inward(x_b)))), JOIN)

    def form(self, where, piece, z_of, size):
        """The profile against a closed form z(x) derived by hand, as a part of the size."""
        pts = np.array(piece.data()["points"])
        self.add(f"{where}: the closed form", float(np.max(np.abs(pts[:, 2] - z_of(pts[:, 0])))) / size, FORM)

    def radius(self, where, piece, rho_of, size):
        pts = np.array(piece.data()["points"])
        self.add(f"{where}: the closed form of rho", float(np.max(np.abs(pts[:, 1] - rho_of(pts[:, 0])))) / size, FORM)

    def on_piece(self, where, piece, points):
        """Every point (X, Y, Z) of a curve or a dot lies on the piece: its distance from the axis
        within the piece's circles, and its height the profile's there, between the two points of
        the profile it falls between, whose chord follows the profile to SAG of the size."""
        P = np.asarray(points, dtype=float).reshape(-1, 3)
        if isinstance(piece, EllipsesPiece):
            # Each point is a node of the row at its height.
            N = piece.nodes()
            worst = 0.0
            for q in P:
                i = int(np.argmin(np.abs(piece.Zr - q[2])))
                worst = max(worst, abs(piece.Zr[i] - q[2]), float(np.min(np.linalg.norm(N[i] - q, axis=1))))
            self.add(f"{where}: on the nodes of its row", worst / piece.size, FORM)
            return
        if isinstance(piece, GridPiece):
            off = np.abs(P[:, 2] - piece.facets.height(P[:, 0], P[:, 1]))
            self.add(f"{where}: on the triangles of the grid", float(np.max(np.where(np.isnan(off), np.inf, off))) / piece.size,
                     FORM)
            return
        pts = np.array(piece.data()["points"])
        rho, z = pts[:, 1], pts[:, 2]
        if np.any(np.diff(rho) <= 0) and np.any(np.diff(rho) >= 0):
            raise AssertionError(f"{where}: the profile of {piece.id} turns, so a height is not one of its radius")
        order = np.argsort(rho)
        r = np.hypot(P[:, 0], P[:, 1])
        outside = float(max(0.0, np.max(rho.min() - r), np.max(r - rho.max())))
        off = float(np.max(np.abs(P[:, 2] - np.interp(r, rho[order], z[order]))))
        self.add(f"{where}: on the surface", max(outside, off) / piece.size, FORM)

    def stops(self, where, sl, xs):
        """g_xx - (drho/dx)^2 < 0 at every sample: no surface of revolution in flat space."""
        d = sl.defect_at(np.asarray(xs, dtype=float))
        self.add(f"{where}: g_xx - (drho/dx)^2 < 0 throughout", float(max(0.0, np.max(d))), 0.0)
        if not np.all(d < 0):
            self.items[-1]["ok"] = False

    def plane(self, where, sl, xs):
        """g_xx - (drho/dx)^2 = 0 at every sample: the slice is a plane."""
        d = sl.defect_at(np.asarray(xs, dtype=float))
        self.add(f"{where}: g_xx - (drho/dx)^2 = 0 throughout", float(np.max(np.abs(d))), 1e-12)

    def exact(self, where, ok):
        self.add(where, 0.0 if ok else 1.0, 0.5)

    def failures(self):
        return [c["name"] for c in self.items if not c["ok"]]

    def worst(self, word):
        errors = [c["error"] for c in self.items if word in c["name"]]
        return max(errors) if errors else 0.0

    def report(self):
        for c in self.items:
            print(f"{'ok    ' if c['ok'] else 'FAILED'} {c['name']}  ({c['error']:.2e}, allowed {c['tol']:.0e})")


# ---------------------------------------------------------------- the drawing

class Scene:
    """Surfaces of revolution about vertical axes, each placed at an offset in the drawing's
    (X, Y, Z), with Z up, and seen from one camera. What hides what is found by casting rays
    through every truncated cone between two neighbouring circles of every profile, which
    is the surface a client draws from the file's points."""

    def __init__(self, camera=CAMERA, size=1.0):
        self.camera, self.size = camera, size
        self.solids = []            # (offset, rho, z, fill class or None)
        self.grids = []             # (offset, Facets, fill class or None), after the solids

    def add(self, piece, offset=(0.0, 0.0, 0.0), fill=None):
        if isinstance(piece, GridPiece):
            self.grids.append((np.asarray(offset, dtype=float), piece.facets, fill))
            return
        self.solids.append((np.asarray(offset, dtype=float), piece.rho.astype(float), piece.z.astype(float), fill))

    def fills(self):
        """The fill class of each solid and then of each grid, in the order reach() gives them."""
        return [fill for *_, fill in self.solids] + [fill for *_, fill in self.grids]

    def reach(self, P, v, chunk=1500):
        """For rays P + t v, t > 0, the largest t at which each ray meets each solid, or -inf:
        an array of shape (len(P), number of solids)."""
        P = np.asarray(P, dtype=float)
        A = v[0] ** 2 + v[1] ** 2
        eps = 1e-7 * self.size
        out = np.full((len(P), len(self.solids)), -np.inf)
        for k, (off, rho, z, _) in enumerate(self.solids):
            r1, dr = rho[:-1], np.diff(rho)
            z1, dz = z[:-1], np.diff(z)
            beta = dz / v[2]
            qa = A * beta ** 2 - dr ** 2
            for s in range(0, len(P), chunk):
                Q = P[s:s + chunk] - off
                Bp = 2 * (Q[:, :1] * v[0] + Q[:, 1:2] * v[1])
                C = Q[:, :1] ** 2 + Q[:, 1:2] ** 2
                alpha = (z1[None, :] - Q[:, 2:3]) / v[2]
                # |P + t v| off the axis equals rho1 + u drho where z1 + u dz = P_z + t v_z.
                qb = Bp * beta + 2 * A * alpha * beta - 2 * r1 * dr
                qc = C + Bp * alpha + A * alpha ** 2 - r1 ** 2
                disc = qb ** 2 - 4 * qa * qc
                ok = disc >= 0
                root = np.sqrt(np.where(ok, disc, 0.0))
                with np.errstate(divide="ignore", invalid="ignore"):
                    lin = np.abs(qa) < 1e-14 * np.maximum(1.0, np.abs(qb))
                    u1 = np.where(lin, -qc / qb, (-qb - root) / (2 * qa))
                    u2 = np.where(lin, -qc / qb, (-qb + root) / (2 * qa))
                best = out[s:s + chunk, k]
                for u in (u1, u2):
                    t = alpha + beta * u
                    hit = ok & np.isfinite(u) & (u >= 0) & (u <= 1) & (t > eps)
                    best = np.maximum(best, np.max(np.where(hit, t, -np.inf), axis=1))
                out[s:s + chunk, k] = best
        if self.grids:
            if not np.allclose(v, self.camera.toward):
                raise AssertionError("a grid is only cast through toward the camera")
            out = np.hstack([out] + [facets.reach(P - off, self.camera, 1e-7 * facets.extent)[:, None]
                                     for off, facets, _ in self.grids])
        return out

    def hidden(self, P):
        """Whether each point, on a surface, is hidden from the camera by any surface."""
        return np.any(np.isfinite(self.reach(P, self.camera.toward)), axis=1)

    def front(self, S):
        """The solid nearest the camera at each point S of the page, or -1 where there is none."""
        S = np.asarray(S, dtype=float)
        cam = self.camera
        P = S[:, :1] * cam.right + S[:, 1:2] * cam.up - 4 * self.size * cam.toward
        t = self.reach(P, cam.toward)
        return np.where(np.isfinite(t).any(axis=1), np.argmax(t, axis=1), -1)


class Figure:
    """The page's drawing: layers painted in the order given, TeX labels at points, a legend,
    in the plane of the page, the form a figure in three dimensions takes. It also records
    what a client needs to draw the same figure from another camera, as the page does when a
    reader turns it: where each surface stands, how many meridians it has, which pieces are
    tinted, the lines marked on a surface besides its circles, the circle each label names,
    and which layers lie flat in the plane of the page and never turn."""

    def __init__(self, scene):
        self.scene, self.camera = scene, scene.camera
        self.fills, self.lines, self.dots, self.labels, self.legend_items = [], [], [], [], []
        self.extent = []
        self.surfaces, self.meridians, self.tint, self.marks = [], None, {}, []
        self.grids = []             # the lines of each grid piece a client draws again, by index
        self.shades = []            # the parts of a piece shaded while a conformal view is shown

    def screen(self, P):
        return self.camera.screen(np.asarray(P, dtype=float))

    def line(self, cls, P, closed=False, normals=None):
        """A line on the surfaces, split into the parts seen, `cls`, and the parts hidden,
        `cls`-far, cut halfway between neighbouring points of opposite kinds. A line along the
        outline comes with the surface's unit normal at each point, and each point is judged by
        two points OUTLINE of the drawing's size off the surface on either side along it: the
        line of sight only grazes the surface there, and between two circles of the profile the
        cone that carries it dips across that line by less than OUTLINE where the profile is
        concave, as on Flamm's paraboloid near its throat, which would hide the outline from
        itself. On the side the surface folds toward the fold hides the point, and on the other
        only what truly lies in front does, so the point is hidden only if both are."""
        P = np.asarray(P, dtype=float)
        if closed:
            P = np.vstack([P, P[:1]])
        if normals is None:
            hid = self.scene.hidden(P)
        else:
            d = OUTLINE * self.scene.size * np.asarray(normals, dtype=float)
            hid = self.scene.hidden(P + d) & self.scene.hidden(P - d)
        S = self.screen(P)
        self.extent.append(S)
        cuts = np.flatnonzero(hid[1:] != hid[:-1])
        start = 0
        for c in list(cuts) + [len(P) - 1]:
            run = S[start:c + 1]
            if c < len(P) - 1:
                run = np.vstack([run, 0.5 * (S[c] + S[c + 1])])
            if start > 0:
                run = np.vstack([0.5 * (S[start - 1] + S[start]), run])
            if len(run) > 1:
                self.lines.append((cls + ("-far" if hid[start] else ""), nr.thin(run, 4e-4 * self.scene.size), False))
            start = c + 1

    def mark(self, cls, k, piece, phi):
        """The meridian at `phi` of `piece` on the figure's surface k, marked in a class of its
        own, as the cut the cone is laid flat along."""
        off = np.asarray(self.surfaces[k][1], dtype=float)
        self.line(cls, off + densify(np.column_stack([piece.rho * math.cos(phi), piece.rho * math.sin(phi),
                                                      piece.z]), 4))
        self.marks.append({"class": cls, "surface": k, "piece": piece.id, "phi": phi})

    def flat(self, cls, S):
        """A line already in the plane of the page, never hidden."""
        S = np.asarray(S, dtype=float)
        self.extent.append(S)
        self.lines.append((cls, S, True))

    def fills_seen(self, n=420, scene=None):
        """Where each fill class is the surface nearest the camera, as polygons of a grid n
        points across the drawing, painted under every line: of the figure's own scene, or of
        `scene`, the same surfaces tinted otherwise, as a shade tints them, on the same grid."""
        scene = scene or self.scene
        fills = scene.fills()
        classes = sorted({fill for fill in fills if fill})
        out = []
        if not classes:
            return out
        E = np.vstack(self.extent)
        lo, hi = E.min(0), E.max(0)
        step = float(np.max(hi - lo)) / n
        X = np.arange(lo[0] - 2 * step, hi[0] + 3 * step, step)
        Y = np.arange(lo[1] - 2 * step, hi[1] + 3 * step, step)
        GX, GY = np.meshgrid(X, Y)
        near = scene.front(np.column_stack([GX.ravel(), GY.ravel()])).reshape(GX.shape)
        for cls in classes:
            which = [k for k, fill in enumerate(fills) if fill == cls]
            mask = np.isin(near, which).astype(float)
            gen = contourpy.contour_generator(X, Y, mask, fill_type=contourpy.FillType.OuterOffset)
            points, offsets = gen.filled(0.5, 1.5)
            for pts, offs in zip(points, offsets):
                rings = [nr.thin(pts[a:b], 5e-4 * self.scene.size) for a, b in zip(offs[:-1], offs[1:])]
                layer = {"kind": "fill", "class": cls, "points": rings[0]}
                if len(rings) > 1:
                    layer["holes"] = rings[1:]
                out.append(layer)
        return out

    def shade(self, view, k, piece, lo, hi, T, reach, legend, caption):
        """Tint the circles of `piece` on the figure's surface k from x = lo to x = hi, both
        points of its profile, in the class shade in place of the piece's own tint, the rest of
        the piece clear, while the view `view` of the spacetime's conformal diagram is shown: the
        part of the moment that the region the view draws covers at its conformal time T.
        `reach` is that part at other conformal times, each [T, lo, hi], or [T, None, None]
        where the region covers none of the moment; `legend` names the tint and `caption` is a
        paragraph read after the view's caption while the shade is shown."""
        if isinstance(piece, GridPiece) or piece.reference:
            raise AssertionError(f"piece {piece.id}: only a profile of the slice is shaded")
        for x in (lo, hi):
            piece.at(x)
        self.shades.append({"view": view, "surface": k, "piece": piece, "from": lo, "to": hi, "T": T,
                            "reach": reach, "legend": legend, "caption": caption})

    def shaded(self, shade):
        """The figure's scene with `shade` in place of its piece's tint: the piece cut at the
        circles lo and hi into cones tinted shade between them and left clear outside."""
        scene = Scene(self.camera, self.scene.size)
        for k, (surface, off) in enumerate(self.surfaces):
            for p in surface.pieces:
                if p.reference:
                    continue
                if k != shade["surface"] or p is not shade["piece"]:
                    scene.add(p, off, self.tint.get(p.cls))
                    continue
                i, j = sorted(int(np.argmin(np.abs(p.x - x))) for x in (shade["from"], shade["to"]))
                for a, b, fill in ((0, i, None), (i, j, "shade"), (j, len(p.x) - 1, None)):
                    if b > a:
                        scene.solids.append((np.asarray(off, dtype=float), p.rho[a:b + 1].astype(float),
                                             p.z[a:b + 1].astype(float), fill))
        return scene

    def plain_fill(self, cls, S):
        """A polygon already in the plane of the page, painted under every line."""
        S = np.asarray(S, dtype=float)
        self.extent.append(S)
        self.fills.append({"kind": "fill", "class": cls, "points": S, "flat": True})

    def dot(self, cls, P):
        S = self.screen(np.asarray(P, dtype=float)[None, :])[0]
        self.dots.append((cls, S))

    def label(self, S, text, anchor="l", cls="lab", dx=0, dy=0, ring=None, clear=None, curve=None, axis=None,
              frame=False):
        """TeX at a point of the page, with an anchor and an offset in units of a figure 628 wide;
        `ring` names the circle it stands beside, and `clear` how far above and below the
        circle's end it looks for the outline to stand past; `curve` names the curve it stands
        beside the end of, `axis` the surface whose axis of time it names at its top, and
        `frame` says it names the frame of a movie shown."""
        L = {"at": S, "text": text, "anchor": anchor, "class": cls, "dx": dx, "dy": dy}
        if ring:
            L["ring"] = ring
        if clear:
            L["clear"] = fixed(clear, 4)
        if curve:
            L["curve"] = curve
        if axis:
            L["axis"] = axis
        if frame:
            L["frame"] = True
        self.labels.append(L)

    def legend(self, kind, cls, text):
        self.legend_items.append([kind, cls, text])

    def done(self, pad=0.04):
        """The figure as it is written, boxed so that every label fits inside at its size, and
        refused if two labels overlap, since each is set on its own ground over the lines."""
        flat_fills = list(self.fills)
        self.fills += self.fills_seen()
        P = np.vstack(self.extent)
        lo, hi = P.min(0), P.max(0)
        for _ in range(3):
            width = (hi[0] - lo[0]) * (1 + 2 * pad)
            unit = width / 560          # data units per unit of a figure 628 wide, of which 560 drawn
            box_lo, box_hi = lo - pad * width, hi + pad * width
            for L in self.labels:
                x0, y0, x1, y1 = label_box(L, unit)
                box_lo = np.minimum(box_lo, [x0 - pad * width / 4, y0 - pad * width / 4])
                box_hi = np.maximum(box_hi, [x1 + pad * width / 4, y1 + pad * width / 4])
            lo, hi = box_lo + pad * width, box_hi - pad * width
        box = [box_lo[0], box_hi[0], box_lo[1], box_hi[1]]
        boxes = [label_box(L, unit) for L in self.labels]
        for i, a in enumerate(boxes):
            for j, b in enumerate(boxes[:i]):
                if a[0] < b[2] and b[0] < a[2] and a[1] < b[3] and b[1] < a[3]:
                    raise AssertionError(f"the labels {self.labels[j]['text']} and {self.labels[i]['text']} overlap")
        def filled(fills):
            return [dict(layer, points=rounded(layer["points"]), **({"holes": [rounded(h) for h in layer["holes"]]}
                                                                    if "holes" in layer else {}))
                    for layer in fills]
        layers = filled(self.fills)
        order = {cls: i for i, (_, cls, _) in enumerate(self.legend_items)}
        # Hidden lines first, then the lines seen, each in the order of the legend.
        for far in (True, False):
            for cls, S, flat in sorted(self.lines, key=lambda item: order.get(item[0].removesuffix("-far"), 99)):
                if cls.endswith("-far") == far:
                    layers.append({"kind": "line", "class": cls, "points": rounded(S), **({"flat": True} if flat else {})})
        for cls, S in self.dots:
            layers.append({"kind": "point", "class": cls, "at": rounded(S)})
        drawn = {layer["class"].removesuffix("-far") for layer in layers} | getattr(self, "later", set())
        for _, cls, _ in self.legend_items:
            if cls not in drawn:
                raise AssertionError(f"the legend names {cls}, which is not drawn")
        labels = [dict(L, at=rounded(L["at"])) for L in self.labels]
        turn = {"origins": [rounded(self.screen(off)) for _, off in self.surfaces], "meridians": self.meridians,
                "tint": self.tint, "marks": self.marks}
        if self.grids:
            turn["grid"] = self.grids
        out = {"box": [fixed(b, 4) for b in box],
               "camera": {"azimuth": self.camera.azimuth, "elevation": self.camera.elevation},
               "layers": layers, "labels": labels, "legend": self.legend_items, "turn": turn}
        if self.shades:
            # The fills painted in place of the figure's own while each shade is shown; view()
            # sets the shades beside the figure.
            out["shades"] = [{"view": sh["view"], "surface": sh["surface"], "piece": sh["piece"].id,
                              "from": significant(sh["from"]), "to": significant(sh["to"]), "T": sh["T"],
                              "reach": [[T] + ([None, None] if lh is None else [fixed(lh[0], 7), fixed(lh[1], 7)])
                                        for T, lh in sh["reach"]],
                              "layers": filled(flat_fills + self.fills_seen(scene=self.shaded(sh))),
                              "legend": sh["legend"], "caption": sh["caption"]}
                             for sh in self.shades]
        return out


def rounded(points):
    return (np.round(np.asarray(points, dtype=float), 4) + 0.0).tolist()


def label_width(text):
    """A generous width of a TeX label in ems: what MathJax sets, less the markup."""
    plain = text.replace("$", "").replace("\\,", " ")
    for cmd in ("\\chi", "\\ell", "\\pi", "\\phi", "\\eta", "\\theta", "\\sqrt", "\\frac", "\\infty", "\\mu",
                "\\delta"):
        plain = plain.replace(cmd, "x")
    plain = plain.replace("{", "").replace("}", "").replace("_", "").replace("^", "")
    return 0.55 * len(plain) + 0.3


ANCHOR = {"l": (0, -0.5), "r": (-1, -0.5), "t": (-0.5, 0), "b": (-0.5, -1), "c": (-0.5, -0.5),
          "tl": (0, 0), "tr": (-1, 0), "bl": (0, -1), "br": (-1, -1)}


def label_box(L, unit):
    """The box a label takes in the drawing, x0, y0, x1, y1, with `unit` the drawing's length
    per unit of a figure 628 wide: generous, as label_width() is, and 1.25 of its size tall,
    which holds the ground the page sets it on."""
    size = LAB[L["class"]] * unit
    w, h = label_width(L["text"]) * size, 1.25 * size
    ax, ay = ANCHOR[L["anchor"]]
    x0 = L["at"][0] + L["dx"] * unit + ax * w
    # The page's y runs down, the drawing's up.
    y0 = L["at"][1] - L["dy"] * unit - (ay + 1) * h
    return x0, y0, x0 + w, y0 + h


def draw_surface(fig, surface, offset=(0.0, 0.0, 0.0), meridians=24):
    """Every piece of a surface on the figure: its outline where it turns edge on to the
    camera, its meridians, its marked circles and the circles at its ends. A reference piece
    is drawn in dashes, with every other meridian. A grid piece is drawn as draw_grid() draws it."""
    off = np.asarray(offset, dtype=float)
    for p in surface.pieces:
        if isinstance(p, GridPiece):
            draw_grid(fig, surface, p, off)
            continue
        style = "reference" if p.reference else "meridian"
        for k in range(0, meridians, 2 if p.reference else 1):
            phi = 2 * math.pi * k / meridians
            P = off + np.column_stack([p.rho * math.cos(phi), p.rho * math.sin(phi), p.z])
            fig.line(style, densify(P, 4))
        for run, normals in outline(p, fig.camera):
            fig.line("reference" if p.reference else "outline", off + run, normals=normals)
        for (kind, _), x in zip(p.ends, (p.lo, p.hi)):
            if kind in ("edge", "stops") and not p.reference:
                fig.line("outline", off + circle(*p.at(x)))
        for x, cls, _ in p.marks:
            rho, z = p.at(x)
            if rho > 0:
                fig.line(cls, off + circle(rho, z))
            else:
                fig.dot(cls, off + [0, 0, z])
    for c in surface.curves:
        if isinstance(c.piece, GridPiece):
            Q = np.vstack([c.points, c.points[:1]]) if c.closed else c.points
            fig.line(c.cls, off + Q, normals=sides(c.piece, Q, fig.scene.size))
        else:
            fig.line(c.cls, off + c.points, closed=c.closed)
    for _, cls, P in surface.dots:
        fig.dot(cls, off + np.asarray(P, dtype=float))
    if surface.axis:
        # The axis of time, straight up the surface's axis, hidden where the surface lies in front.
        fig.line("axis", off + finely(np.array([[0.0, 0.0, surface.axis[0]], [0.0, 0.0, surface.axis[1]]]),
                                      fig.scene.size / 360))


def draw_grid(fig, surface, piece, off):
    """A grid piece on the figure: the lines of its grid that `lines` names, in the class grid, its
    rim and its outline, each densified as a meridian is. Every point of every line on it is
    judged by the two points OUTLINE of the grid's extent above and below it, and hidden only if
    both are: a line lies on the triangles, and above and below them are the two sides of a
    height, as the two sides of an outline are. A stack of ellipses draws the lines it names in
    their own classes and no rim, each point judged straight out from the axis and back, the two
    sides of a tube or a cone."""
    k = next(i for i, (s, _) in enumerate(fig.surfaces) if s is surface)
    P = piece.nodes()

    def line(cls, Q, closed=False):
        Q = densify(np.vstack([Q, Q[:1]]) if closed else Q, 4)
        fig.line(cls, Q + off, normals=sides(piece, Q, fig.scene.size))
    if isinstance(piece, EllipsesPiece):
        for cls, which, idx in piece.drawn:
            for i in idx:
                if which == "u":
                    line(cls, P[i], closed=piece.wrap)
                else:
                    line(cls, P[:, i])
    else:
        for i in piece.lines["u"]:
            line("grid", P[i], closed=piece.wrap)
        for j in piece.lines["v"]:
            line("grid", P[:, j])
        line("outline", piece.rim(), closed=True)
    for run in grid_outline(piece, fig.camera):
        fig.line("outline", run + off, normals=sides(piece, run, fig.scene.size))
    if not isinstance(piece, EllipsesPiece):
        fig.grids.append({"class": "grid", "surface": k, "piece": piece.id, "u": piece.lines["u"], "v": piece.lines["v"]})


def densify(P, n):
    """n points on each segment of a polyline, so that a line changes from seen to hidden
    close to where it does."""
    P = np.asarray(P, dtype=float)
    s = np.linspace(0, 1, n + 1)[:-1]
    out = (P[:-1, None, :] + s[None, :, None] * (P[1:] - P[:-1])[:, None, :]).reshape(-1, 3)
    return np.vstack([out, P[-1:]])


def finely(P, most):
    """A polyline with each segment cut into as few equal pieces as keep every piece shorter than
    `most`, as a meridian's chords, a quarter of at most 1/90 of the drawing, are shorter than
    1/360 of it, so that a line on a grid changes from seen to hidden close to where it does."""
    P = np.asarray(P, dtype=float)
    out = [P[:1]]
    for a, b in zip(P[:-1], P[1:]):
        n = max(1, int(math.ceil(float(np.linalg.norm(b - a)) / most)))
        s = np.linspace(0, 1, n + 1)[1:]
        out.append(a + s[:, None] * (b - a))
    return np.vstack(out)


def circle(rho, z, n=RING):
    phi = np.linspace(0, 2 * math.pi, n + 1)
    return np.column_stack([rho * np.cos(phi), rho * np.sin(phi), np.full_like(phi, z)])


def outline(p, cam):
    """Where a piece turns edge on to the camera: on the circle of each point of the profile,
    the angles at which the normal (dz cos phi, dz sin phi, -drho) is square to the line of
    sight, cos(phi - azimuth) = (drho/dz) tan(elevation), joined from point to point. Each run
    is its points and the unit normal at each, which Figure.line() judges them by."""
    rho, z = p.rho, p.z
    drho, dz = np.gradient(rho), np.gradient(z)
    length = np.hypot(drho, dz)
    tan_e = math.tan(math.radians(cam.elevation))
    a = math.radians(cam.azimuth)
    with np.errstate(divide="ignore", invalid="ignore"):
        c = drho * tan_e / dz
    ok = np.isfinite(c) & (np.abs(c) <= 1) & (rho > 0)
    runs, i = [], 0
    while i < len(rho):
        if not ok[i]:
            i += 1
            continue
        j = i
        while j < len(rho) and ok[j]:
            j += 1
        ang = np.arccos(np.clip(c[i:j], -1, 1))
        nr, nz = dz[i:j] / length[i:j], -drho[i:j] / length[i:j]
        side = [(np.column_stack([rho[i:j] * np.cos(a + s * ang), rho[i:j] * np.sin(a + s * ang), z[i:j]]),
                 np.column_stack([nr * np.cos(a + s * ang), nr * np.sin(a + s * ang), nz])) for s in (1, -1)]
        # The two sides meet where the outline turns, |c| = 1; otherwise each runs to an end.
        if j - i > 1:
            if i > 0 and j < len(rho):
                runs.append((np.vstack([side[0][0][::-1], side[1][0]]), np.vstack([side[0][1][::-1], side[1][1]])))
            else:
                runs += side
        i = j
    return runs


def ring_label(fig, off, rho, z, text, side=1, cls="small", dx=8, dy=0, clear=False):
    """A label beside the right (side 1) or left end of a circle, on the page. With `clear`
    it stands past the surface's outline instead, where the outline runs outside the circle's
    end within the label's height, as a cone's sides do below its rim. The circle must be
    exactly one of the circles the figure's surfaces mark, which the label names, so that a
    client turning the figure keeps the label beside it."""
    found = [(k, i) for k, (s, o) in enumerate(fig.surfaces) if np.array_equal(np.asarray(o, dtype=float), off)
             for i, (p, x) in enumerate((p, x) for p in s.pieces for x, _, _ in p.marks) if p.at(x) == (rho, z)]
    if len(found) != 1:
        raise AssertionError(f"the label {text} stands beside {len(found)} marked circles, not one")
    S = fig.screen(np.asarray(off, dtype=float) + [side * rho, 0, z])
    band = 0.02 * fig.scene.size
    if clear:
        for cls_, P, _ in fig.lines:
            if cls_ != "outline":
                continue
            for y in (S[1] - band, S[1], S[1] + band):
                a, b = P[:-1], P[1:]
                cross = ((a[:, 1] - y) * (b[:, 1] - y) <= 0) & (a[:, 1] != b[:, 1])
                x = a[cross, 0] + (y - a[cross, 1]) * (b[cross, 0] - a[cross, 0]) / (b[cross, 1] - a[cross, 1])
                if x.size:
                    S[0] = max(S[0], x.max()) if side > 0 else min(S[0], x.min())
    (k, i), = found
    fig.label(S, text, "l" if side > 0 else "r", cls, dx=side * dx, dy=dy,
              ring={"surface": k, "ring": i, "side": side}, clear=band if clear else None)


def curve_label(fig, k, i, text, side=1, cls="small", dx=8, dy=0):
    """A label beside the end of the curve i of the figure's surface k, its point farthest right
    (side 1) or left on the page, the first of them where two tie, as a client turning the
    figure finds it again from the curve's points as the file holds them."""
    surface, off = fig.surfaces[k]
    P = np.array(surface.curves[i].data()["points"], dtype=float)
    S = fig.screen(np.asarray(off, dtype=float) + P)
    j = int(np.argmax(side * S[:, 0]))
    fig.label(S[j], text, "l" if side > 0 else "r", cls, dx=side * dx, dy=dy,
              curve={"surface": k, "curve": i, "side": side})


def axis_label(fig, k, text):
    """The name of the axis of time of the figure's surface k, at its top, as a figure of light
    cones names its axis."""
    surface, off = fig.surfaces[k]
    top = np.asarray(off, dtype=float) + [0.0, 0.0, surface.data()["axis"]["to"]]
    fig.label(fig.screen(top), text, "b", "lab", dy=-4, axis={"surface": k})


def movie_values(keys, step):
    """The values a movie's frames stand at: each interval between neighbouring key values cut
    into as many equal steps as come nearest `step`, so that every key value is a frame. Returns
    the values and the index of each key among them."""
    out, at = [keys[0]], [0]
    for a, b in zip(keys, keys[1:]):
        n = max(1, round((b - a) / step))
        out += [a + (b - a) * k / n for k in range(1, n)] + [b]
        at.append(len(out) - 1)
    return out, at


def movie_figure(frames, fills, size, camera=CAMERA, meridians=12):
    """A movie's drawing: its first frame, on its own axis at the origin, in a box that holds
    every frame, since every frame is drawn in turn at one scale about one place, with a label
    below the lowest point any frame reaches that names the frame shown."""
    fig = figure_of(frames[:1], fills, size, camera, meridians=meridians)
    phi = np.linspace(0, 2 * math.pi, 73)
    for s in frames:
        for p in s.pieces:
            if isinstance(p, GridPiece):
                fig.extent.append(fig.screen(p.nodes().reshape(-1, 3)))
            elif not p.reference:
                fig.extent.append(fig.screen(np.column_stack([np.outer(p.rho, np.cos(phi)).ravel(),
                                                              np.outer(p.rho, np.sin(phi)).ravel(),
                                                              np.repeat(p.z, len(phi))])))
    E = np.vstack(fig.extent)
    fig.label(np.array([float(fig.screen(np.zeros(3))[0]), float(E[:, 1].min())]), frames[0].label, "t", "small",
              dy=8, frame=True)
    # What later frames mark, which the legend may name though the first frame lacks it.
    fig.later = {cls for s in frames for p in s.pieces for _, cls, _ in p.marks} | {c.cls for s in frames for c in s.curves}
    return fig


def movie(frames, variable, values, seconds=5, turns=True, loop=None):
    """What a view carries to play its frames in turn, each with its label and its value of the
    movie's `variable`, one pass taking `seconds`, and played forward, back and forward again
    where `loop` is "pingpong"."""
    out = {"variable": variable, "seconds": seconds, "frames": []}
    if loop:
        out["loop"] = loop
    for s, v in zip(frames, values):
        d = s.data()
        d.pop("time", None)
        out["frames"].append({"label": d.pop("label"), "value": fixed(v, 6), **d})
    if not turns:
        out["turns"] = False
    return out


# ---------------------------------------------------------------- the spacetimes

def view(vid, label, unit, surfaces, figure, **fields):
    shades = figure.pop("shades", None)
    out = {"id": vid, "label": label, "unit": unit, "surfaces": [s.data() for s in surfaces],
           "figure": figure}
    out.update({k: v for k, v in fields.items() if v is not None})
    if shades:
        out["shades"] = shades
    return out


def figure_of(surfaces, fills, size, camera=CAMERA, offsets=None, meridians=24):
    """A figure of surfaces, each piece tinted by `fills` where it is nearest the camera. A
    reference piece is not part of the slice, so it hides nothing."""
    offsets = offsets or [(0.0, 0.0, 0.0)] * len(surfaces)
    scene = Scene(camera, size)
    for s, off in zip(surfaces, offsets):
        for p in s.pieces:
            if not p.reference:
                scene.add(p, off, fills.get(p.cls))
    fig = Figure(scene)
    fig.surfaces, fig.meridians, fig.tint = list(zip(surfaces, offsets)), meridians, dict(fills)
    for s, off in zip(surfaces, offsets):
        draw_surface(fig, s, off, meridians)
    return fig


def sequence_figure(surfaces, fills, size, columns, camera=CAMERA, meridians=12, gap=0.15):
    """A sequence of surfaces in rows of `columns`, read left to right and down, each on its own
    axis and centred in its column, the tops of a row level, and each moment's label set below
    its row. `gap` is the space between columns as a part of the widest surface."""
    phi = np.linspace(0, 2 * math.pi, 73)
    boxes = []
    for s in surfaces:
        P = np.vstack([np.column_stack([np.outer(p.rho, np.cos(phi)).ravel(), np.outer(p.rho, np.sin(phi)).ravel(),
                                        np.repeat(p.z, len(phi))]) for p in s.pieces if not p.reference])
        S = camera.screen(P)
        boxes.append((S.min(0), S.max(0)))
    width = max(hi[0] - lo[0] for lo, hi in boxes)
    height = max(hi[1] - lo[1] for lo, hi in boxes)
    unit = columns * width * (1 + gap) / 560
    row = height + (8 + 1.25 * LAB["small"] + 14) * unit
    lift = float(camera.screen([0.0, 0.0, 1.0])[1] - camera.screen([0.0, 0.0, 0.0])[1])
    offsets = []
    for k, (lo, hi) in enumerate(boxes):
        top = -(k // columns) * row
        offsets.append(((k % columns) * width * (1 + gap) - 0.5 * (lo[0] + hi[0]), 0.0, (top - hi[1]) / lift))
    fig = figure_of(surfaces, fills, size, camera, offsets, meridians)
    for k, (s, off) in enumerate(zip(surfaces, offsets)):
        at = np.array([float(fig.screen(np.asarray(off))[0]), -(k // columns) * row - height])
        fig.label(at, s.label, "t", "small", dy=8)
    return fig


def schwarzschild(ck, src):
    """Flamm's paraboloid. The slice of constant Schwarzschild t through the equator has
    g_rr = 1/(1 - r_s/r) and g_phiphi = r^2, so dz/dr = sqrt(r_s/(r - r_s)) and z^2 =
    4 r_s (r - r_s): every slice of constant t passes through the bifurcation sphere r = r_s,
    and the slice of Kruskal's T = 0 runs on to the other exterior, the same paraboloid
    upside down. Drawn at r_s = 1 to r = 6 on both sheets."""
    sl = Slice(src, "schwarzschild", "spherical", "r", "\\phi", {"t": 0, **EQUATOR}, {"r_s": 1})
    top, radii = 6.0, (1.5, 2, 3, 4, 5)
    size = 2 * top
    near = Piece("exterior", "sheet", sl, 1.0, top, 0.0, 1,
                 (("throat", "the throat $r = r_s$, the bifurcation sphere, where the other exterior begins"),
                  ("edge", "the paraboloid runs on to $r \\to \\infty$")),
                 [(1.0, "horizon", "$r = r_s$")] + [(r, "r", None) for r in radii] + [(top, "r", "$r = 6\\,r_s$")],
                 size)
    far = Piece("other_exterior", "sheet2", sl, 1.0, top, 0.0, -1,
                (("throat", "the throat $r = r_s$"), ("edge", "the paraboloid runs on to $r \\to \\infty$")),
                [(r, "r2", None) for r in radii] + [(top, "r2", None)], size)
    surface = Surface([near, far])
    for p in (near, far):
        ck.isometry(f"Schwarzschild, {p.id}", p)
        ck.form(f"Schwarzschild, {p.id}, Flamm's z = 2 sqrt(r_s (r - r_s))", p,
                lambda r, s=p.sense: s * 2 * np.sqrt(np.maximum(r - 1, 0)), size)
    ck.join("Schwarzschild, the two sheets at the throat", near, 1.0, far, 1.0)

    fig = figure_of([surface], {"sheet": "cover"}, size)
    ring_label(fig, [0, 0, 0], 1.0, 0.0, "$r = r_s$", dx=14)
    ring_label(fig, [0, 0, 0], *near.at(3.0), "$3\\,r_s$")
    ring_label(fig, [0, 0, 0], *near.at(top), "$6\\,r_s$")
    fig.legend("fill", "cover", "the exterior $r > r_s$ that $t$ and $r$ cover")
    fig.legend("line", "r", "$r$ constant, at $1.5$, $2$, $3$, $4$, $5$ and $6\\,r_s$")
    fig.legend("line", "r2", "the same radii on the other exterior")
    fig.legend("line", "horizon", "the throat $r = r_s$, where the slice crosses the horizon")
    fig.legend("line", "meridian", "$\\phi$ constant, every $15°$")
    return [view("flamm", "Flamm's paraboloid", "$r_s$", [surface], fig.done(),
                 settings="$r_s = 1$, the unit of every length.")]


def interior_schwarzschild(ck, src):
    """A star of uniform density, R = 1.5 r_s as the conformal diagram draws it. Inside,
    g_rr = 1/(1 - r^2 r_s/R^3), the slice is a cap of a sphere of radius a = sqrt(R^3/r_s):
    z = a - sqrt(a^2 - r^2) from the centre. Outside it is Flamm's paraboloid from the
    schwarzschild entry, and at r = R both give g_rr = 1/(1 - r_s/R), so the cap meets the
    paraboloid with one tangent. The vacuum paraboloid is drawn on under the cap, down to
    the throat the star does not have."""
    R, top = 1.5, 4.0
    size = 2 * top
    inner = Slice(src, "interior_schwarzschild", "spherical", "r", "\\phi", {"t": 0, **EQUATOR},
                  {"r_s": 1, "R": "3/2"})
    outer = Slice(src, "schwarzschild", "spherical", "r", "\\phi", {"t": 0, **EQUATOR}, {"r_s": 1})
    vacuum = Piece("vacuum", "reference", outer, 1.0, R, 0.0, 1,
                   (("throat", "the throat $r = r_s$ of the vacuum, which the star replaces"), ("join", None)),
                   [(1.0, "reference", None)], size, reference=True)
    zR = vacuum.at(R)[1]
    ext = Piece("exterior", "sheet", outer, R, top, zR, 1,
                (("join", "the surface of the star, $r = R$"), ("edge", "the paraboloid runs on to $r \\to \\infty$")),
                [(R, "surface", "$r = R$")] + [(r, "r", None) for r in (2, 3)] + [(top, "r", None)], size)
    star = Piece("star", "star", inner, 0.0, R, 0.0, 1,
                 (("axis", "the centre $r = 0$, where the cap is smooth"), ("join", "the surface of the star, $r = R$")),
                 [(r, "r", None) for r in (0.5, 1.0)], size)
    # The cap is built from its centre and moved up to meet the exterior at R.
    star.z = star.z + (zR - star.z[-1])
    surface = Surface([star, ext, vacuum])
    ck.isometry("interior Schwarzschild, the star", star)
    ck.isometry("interior Schwarzschild, the exterior", ext)
    ck.join("interior Schwarzschild, the star meets the exterior at r = R", star, R, ext, R)
    a = math.sqrt(R ** 3)
    ck.form("interior Schwarzschild, the cap z = a - sqrt(a^2 - r^2)", star,
            lambda r: star.z[0] + a - np.sqrt(a * a - r * r), size)
    ck.form("interior Schwarzschild, the exterior is Flamm's", ext, lambda r: 2 * np.sqrt(r - 1), size)

    fig = figure_of([surface], {"star": "star", "sheet": "cover"}, size, Camera(-90, 32))
    ring_label(fig, [0, 0, 0], *ext.at(R), "$r = R$", dx=10)
    ring_label(fig, [0, 0, 0], *ext.at(3.0), "$3\\,r_s$")
    ring_label(fig, [0, 0, 0], *ext.at(top), "$4\\,r_s$")
    ring_label(fig, [0, 0, 0], *vacuum.at(1.0), "$r_s$", side=-1, dx=8)
    fig.legend("fill", "star", "the star, $r \\le R$, a cap of a sphere of radius $\\sqrt{R^3/r_s}$")
    fig.legend("fill", "cover", "the exterior, Flamm's paraboloid")
    fig.legend("line", "r", "$r$ constant, at $0.5$ and $1\\,r_s$ inside and $2$, $3$ and $4\\,r_s$ outside")
    fig.legend("line", "surface", "the surface of the star, $r = R$")
    fig.legend("line", "reference", "the vacuum paraboloid inside $R$, down to its throat at $r_s$")
    fig.legend("line", "meridian", "$\\phi$ constant, every $15°$")
    return [view("star", "The star and its exterior", "$r_s$", [surface], fig.done(),
                 settings="$r_s = 1$, the unit of every length, and $R = 1.5\\,r_s$.")]


def tov(ck, src):
    """The declared neutron star: the polytrope p = K rho_0^2 at K = 100 and central rho_0 =
    1.28e-3, G = c = M_sun = 1, which null_rays.StarSolver solves from this spacetime's own
    Einstein tensor and the conformal diagram draws. On the slice g_rr = r/(r - 2m(r)), with the
    solver's m, so dz/dr = sqrt(2m/(r - 2m)). Its Gaussian curvature is m'/r^2 - m/r^3 =
    4 pi (rho - rho_mean/3): positive at the centre, negative in the outer layers, where the
    density falls below a third of the mean inside, and -M/r^3 outside, where m = M and the
    slice is Flamm's paraboloid of that mass, which is checked. The vacuum paraboloid is drawn
    on under the star down to its throat at 2M, as under Schwarzschild's star."""
    solver = nr.StarSolver("tov", "spherical", 100.0, 1.28e-3)
    src.note("tov", "spherical", ["einstein_tensor"])
    M, R = solver.M, solver.R
    ck.add("TOV: the declared star solves the published G^theta_theta = 8 pi p",
           float(np.max(np.abs(solver.theta_theta(np.linspace(0.05 * R, 0.95 * R, 200))))), 1e-7)

    def mass(k):
        return lambda x: solver.values(x)["m"][k].reshape(np.shape(x))
    top = 3 * R
    size = 2 * top
    inner = Slice(src, "tov", "spherical", "r", "\\phi", {"t": 0, **EQUATOR}, numeric={"m": (mass(0), mass(1))})
    outer = Slice(src, "schwarzschild", "spherical", "r", "\\phi", {"t": 0, **EQUATOR}, {"r_s": repr(2 * M)})
    vacuum = Piece("vacuum", "reference", outer, 2 * M, R, 0.0, 1,
                   (("throat", "the throat $r = 2GM/c^2$ of the vacuum, which the star replaces"), ("join", None)),
                   [(2 * M, "reference", None)], size, reference=True)
    zR = vacuum.at(R)[1]
    ext = Piece("exterior", "sheet", inner, R, top, zR, 1,
                (("join", "the surface of the star, $r = R$"), ("edge", "the paraboloid runs on to $r \\to \\infty$")),
                [(R, "surface", "$r = R$"), (2 * R, "r", None), (top, "r", None)], size)
    star = Piece("star", "star", inner, 0.0, R, 0.0, 1,
                 (("axis", "the centre $r = 0$, where the surface is flat"), ("join", "the surface of the star, $r = R$")),
                 [(R / 3, "r", None), (2 * R / 3, "r", None)], size)
    star.z = star.z + (zR - star.z[-1])
    surface = Surface([star, ext, vacuum])
    ck.isometry("TOV, the star", star)
    ck.isometry("TOV, the exterior", ext)
    ck.join("TOV, the star meets the exterior at r = R", star, R, ext, R)
    ck.form("TOV, the exterior is Flamm's of mass M", ext, lambda r: 2 * np.sqrt(2 * M * (r - 2 * M)), size)

    fig = figure_of([surface], {"star": "star", "sheet": "cover"}, size, Camera(-90, 32))
    ring_label(fig, [0, 0, 0], *ext.at(R), "$r = R$", dx=10)
    ring_label(fig, [0, 0, 0], *ext.at(2 * R), "$2R$")
    ring_label(fig, [0, 0, 0], *ext.at(top), "$3R$")
    ring_label(fig, [0, 0, 0], *vacuum.at(2 * M), "$r_s$", side=-1, dx=8)
    fig.legend("fill", "star", "the star, where $p > 0$")
    fig.legend("fill", "cover", "the exterior, Flamm's paraboloid of the star's mass")
    fig.legend("line", "r", "$r$ constant, at $R/3$ and $2R/3$ inside and $2R$ and $3R$ outside")
    fig.legend("line", "surface", "the surface of the star, $r = R$, where the pressure falls to zero")
    fig.legend("line", "reference", "the vacuum paraboloid inside $R$, down to its throat at $r_s = 2GM/c^2$")
    fig.legend("line", "meridian", "$\\phi$ constant, every $15°$")
    km = 1.4766250614  # GM_sun/c^2 in km
    return [view("star", "The star and its exterior", "$GM_\\odot/c^2$", [surface], fig.done(),
                 settings=f"$G = c = M_\\odot = 1$, so that the unit of every length is $GM_\\odot/c^2 = {km:.2f}$ km.",
                 input="A polytrope, $p = K\\rho_0^2$ with rest mass density $\\rho_0$ and energy density "
                       "$\\rho c^2 = \\rho_0c^2 + p$, at $K = 100$ and a central $\\rho_0 = 1.28\\times10^{-3}$, "
                       f"solved from this spacetime's own $G^t{{}}_t$ and $G^r{{}}_r$: a star of $M = {M:.2f}\\,M_\\odot$ "
                       f"and $R = {R * km:.1f}$ km, the star on which numerical relativists test their codes.")]


def morris_thorne(ck, src):
    """The Ellis-Bronnikov member, Phi = 0 and b = b_0^2/r, as the conformal diagram and the
    null rays draw it; the embedding reads b alone. g_rr = r^2/(r^2 - b_0^2) gives dz/dr =
    b_0/sqrt(r^2 - b_0^2), the catenoid z = b_0 arccosh(r/b_0), on both sides of the throat.
    The proper radial chart, with r(l) = sqrt(l^2 + b_0^2), is checked to give the same
    surface."""
    sl = Slice(src, "morris_thorne", "spherical", "r", "\\phi", {"t": 0, **EQUATOR}, {"b_0": 1},
               {"Phi": "0", "b": "b_0**2/r"})
    top, radii = 5.0, (1.5, 2, 3, 4)
    size = 2 * top
    near = Piece("near", "sheet", sl, 1.0, top, 0.0, 1,
                 (("throat", "the throat $r = b_0$, where the two sides join"),
                  ("edge", "the side runs on, flattening, to $r \\to \\infty$")),
                 [(1.0, "throat", "$r = b_0$")] + [(r, "r", None) for r in radii] + [(top, "r", None)], size)
    far = Piece("far", "sheet2", sl, 1.0, top, 0.0, -1,
                (("throat", "the throat $r = b_0$"), ("edge", "the other side runs on, flattening, to $r \\to \\infty$")),
                [(r, "r2", None) for r in radii] + [(top, "r2", None)], size)
    for p in (near, far):
        ck.isometry(f"Morris-Thorne, {p.id} side", p)
        ck.form(f"Morris-Thorne, {p.id} side, the catenoid z = b_0 arccosh(r/b_0)", p,
                lambda r, s=p.sense: s * np.arccosh(np.maximum(r, 1)), size)
    ck.join("Morris-Thorne, the two sides at the throat", near, 1.0, far, 1.0)
    proper = Slice(src, "morris_thorne", "proper_radial", "l", "\\phi", {"t": 0, **EQUATOR}, {"b_0": 1},
                   {"Phi": "0", "r": "sqrt(l**2 + b_0**2)"})
    whole = Piece("proper", "sheet", proper, 0.0, math.sqrt(top ** 2 - 1), 0.0, 1, size=size)
    ck.isometry("Morris-Thorne, the proper radial chart", whole)
    # The same surface: at each of its points, the areal chart's height at r = rho.
    heights = np.array([near.sl.rise(1.0, r) for r in whole.rho])
    ck.add("Morris-Thorne, the proper radial chart gives the same surface",
           float(np.max(np.abs(heights - whole.z))) / size, FORM)
    surface = Surface([near, far])

    fig = figure_of([surface], {"sheet": "cover"}, size)
    ring_label(fig, [0, 0, 0], 1.0, 0.0, "$r = b_0$", dx=14)
    ring_label(fig, [0, 0, 0], *near.at(3.0), "$3\\,b_0$")
    ring_label(fig, [0, 0, 0], *far.at(3.0), "$3\\,b_0$")
    fig.legend("fill", "cover", "the side $l > 0$ that $t$ and $r$ cover")
    fig.legend("line", "r", "$r$ constant, at $1.5$, $2$, $3$, $4$ and $5\\,b_0$")
    fig.legend("line", "r2", "the same radii on the other side, $l < 0$")
    fig.legend("line", "throat", "the throat $r = b_0$, the smallest circle")
    fig.legend("line", "meridian", "$\\phi$ constant, every $15°$")
    return [view("wormhole", "The wormhole", "$b_0$", [surface], fig.done(),
                 settings="$b_0 = 1$, the unit of every length.",
                 input="$\\Phi = 0$ and $b = b_0^2/r$, the member of the family that is the Ellis-Bronnikov "
                       "wormhole; $\\Phi$ does not enter the surface.")]


def thin_shell_wormhole(ck, src):
    """Flamm's paraboloid cut at the throat r = a and joined to its mirror image, drawn at r_s = 1
    and a = 5/4. On either side dz/dr = sqrt(r_s/(r - r_s)), so z = +-2(sqrt(r_s(r - r_s)) -
    sqrt(r_s(a - r_s))) from the throat, where both sides climb at sqrt(r_s/(a - r_s)) = 2 and
    so meet at a crease: the tangent turns through 2 arctan(sqrt((a - r_s)/r_s)) there, the shell.
    The chart through the throat, r = a + |l|, is checked to give the same surface."""
    params = {"r_s": 1, "a": "5/4"}
    a, top, radii = 1.25, 5.0, (1.5, 2, 3, 4)
    size = 2 * top
    sl = Slice(src, "thin_shell_wormhole", "spherical", "r", "\\phi", {"t": 0, **EQUATOR}, params)
    near = Piece("near", "sheet", sl, a, top, 0.0, 1,
                 (("crease", "the throat $r = a$, where the two sides meet at the shell"),
                  ("edge", "the side runs on, flattening, to $r \\to \\infty$")),
                 [(a, "throat", "$r = a$")] + [(r, "r", None) for r in radii] + [(top, "r", None)], size)
    far = Piece("far", "sheet2", sl, a, top, 0.0, -1,
                (("crease", "the throat $r = a$"), ("edge", "the other side runs on, flattening, to $r \\to \\infty$")),
                [(r, "r2", None) for r in radii] + [(top, "r2", None)], size)
    for p in (near, far):
        ck.isometry(f"thin shell wormhole, {p.id} side", p)
        ck.form(f"thin shell wormhole, {p.id} side, Flamm's z = 2 sqrt(r_s (r - r_s)) less its value at a", p,
                lambda r, s=p.sense: s * 2 * (np.sqrt(np.maximum(r - 1, 0)) - 0.5), size)
    # The two sides meet in one circle, and each leaves it climbing away at dz/drho = 2, so the
    # tangents pointing into them are mirror images across the plane of the throat.
    pa, pb = np.array(near.at(a)), np.array(far.at(a))
    ck.add("thin shell wormhole, the two sides at the throat: one circle", float(np.max(np.abs(pa - pb))), JOIN)
    ta, tb = near.inward(a), far.inward(a)
    ck.add("thin shell wormhole, the crease: each side leaves the throat at dz/drho = 2",
           float(max(abs(ta[1] / ta[0] - 2), abs(tb[1] / tb[0] + 2))), JOIN)
    through = Slice(src, "thin_shell_wormhole", "throat", "\\ell", "\\phi", {"t": 0, **EQUATOR}, params)
    whole = Piece("through", "sheet", through, 1e-9, top - a, 0.0, 1, size=size)
    ck.isometry("thin shell wormhole, the chart through the throat", whole)
    heights = np.array([near.sl.rise(a, r) for r in whole.rho])
    ck.add("thin shell wormhole, the chart through the throat gives the same surface",
           float(np.max(np.abs(heights - whole.z))) / size, FORM)
    surface = Surface([near, far])

    fig = figure_of([surface], {"sheet": "cover"}, size)
    ring_label(fig, [0, 0, 0], *near.at(a), "$r = a$", dx=14)
    ring_label(fig, [0, 0, 0], *near.at(3.0), "$3\\,r_s$")
    ring_label(fig, [0, 0, 0], *far.at(3.0), "$3\\,r_s$")
    fig.legend("fill", "cover", "the side $\\ell > 0$ that $t$ and $r$ cover")
    fig.legend("line", "r", "$r$ constant, at $1.5$, $2$, $3$, $4$ and $5\\,r_s$")
    fig.legend("line", "r2", "the same radii on the other side, $\\ell < 0$")
    fig.legend("line", "throat", "the throat $r = a$, the smallest circle, where the shell sits")
    fig.legend("line", "meridian", "$\\phi$ constant, every $15°$")
    return [view("wormhole", "The wormhole", "$r_s$", [surface], fig.done(),
                 settings="$r_s = 1$, the unit of every length, and $a = 1.25$.")]


def two_sheets(ck, name, sl, throat, top, radii, size, near_marks=(), texts=("", "")):
    """A slice of constant t through a bifurcation sphere, as Schwarzschild's: the exterior from
    the throat out to `top`, tinted, and the same surface turned over on the other side."""
    near = Piece("exterior", "sheet", sl, throat, top, 0.0, 1,
                 (("throat", f"the throat $r = r_+$, the bifurcation sphere{texts[0]}, where the other exterior begins"),
                  ("edge", "the surface runs on to $r \\to \\infty$")),
                 list(near_marks) + [(r, "r", None) for r in radii] + [(top, "r", None)], size)
    far = Piece("other_exterior", "sheet2", sl, throat, top, 0.0, -1,
                (("throat", "the throat $r = r_+$"), ("edge", "the surface runs on to $r \\to \\infty$")),
                [(r, "r2", None) for r in radii] + [(top, "r2", None)], size)
    for p in (near, far):
        ck.isometry(f"{name}, {p.id}", p)
    ck.join(f"{name}, the two sheets at the throat", near, throat, far, throat)
    return near, far


def rn_metric(ck, src):
    """Reissner-Nordstrom at r_q = 0.48 r_s, as the conformal diagram draws it, so that r+ =
    0.64 and r- = 0.36 r_s. g_rr = r^2/(r^2 - r_s r + r_q^2): outside r+ the slice of constant t
    runs through the outer bifurcation sphere into a second exterior; between the horizons g_rr
    < 0 and it is not a moment of space; inside r- it runs through the inner bifurcation sphere,
    the widest circle there, into a second region inside r-, and g_rr - 1 = (r_s r - r_q^2)/
    (...) falls to zero at r = r_q^2/r_s, where the surface lies level, and is negative nearer
    the singularity, where no surface of revolution in flat space carries the slice."""
    sl = Slice(src, "rn_metric", "spherical", "r", "\\phi", {"t": 0, **EQUATOR}, {"r_s": 1, "r_q": "12/25"})
    rp, rm = sl.horizons()
    level = 0.48 ** 2
    ck.add("Reissner-Nordstrom: the horizons are at 0.64 and 0.36 r_s", abs(rp - 0.64) + abs(rm - 0.36), 1e-12)
    ck.stops("Reissner-Nordstrom, between the horizons", sl, np.linspace(rm, rp, 402)[1:-1])
    ck.stops("Reissner-Nordstrom, nearer the singularity than r_q^2/r_s", sl, np.linspace(0, level, 402)[1:-1])
    top, radii = 6.0, (1.0, 2.0, 3.0, 4.0, 5.0)
    size = 2 * top
    near, far = two_sheets(ck, "Reissner-Nordstrom outside", sl, rp, top, radii, size,
                           [(rp, "horizon", "$r = r_+$")], (" of the outer horizon", ""))
    outside = Surface([near, far])
    fig = figure_of([outside], {"sheet": "cover"}, size)
    ring_label(fig, [0, 0, 0], rp, 0.0, "$r = r_+$", dx=14)
    ring_label(fig, [0, 0, 0], *near.at(3.0), "$3\\,r_s$")
    ring_label(fig, [0, 0, 0], *near.at(top), "$6\\,r_s$")
    fig.legend("fill", "cover", "the exterior $r > r_+$ that $t$ and $r$ cover")
    fig.legend("line", "r", "$r$ constant, at $1$, $2$, $3$, $4$, $5$ and $6\\,r_s$")
    fig.legend("line", "r2", "the same radii on the other exterior")
    fig.legend("line", "horizon", "the throat $r = r_+$, where the slice crosses the outer horizon")
    fig.legend("line", "meridian", "$\\phi$ constant, every $15°$")
    settings = "$r_s = 1$, the unit of every length, and $r_q = 0.48\\,r_s$, so that $r_+ = 0.64\\,r_s$ and $r_- = 0.36\\,r_s$."
    between = ("Between the horizons, $r_- < r < r_+$, $g_{rr} < 0$: $r$ is a time there, and a slice of "
               "constant $t$ is not a moment of space.")
    views = [view("outside", "Outside $r_+$", "$r_s$", [outside], fig.done(), settings=settings, stops=[between])]

    # Inside r-: each side from where the surface lies level up to the widest circle, r-.
    size = 2 * rm
    lo = Piece("inside", "sheet", sl, level, rm, 0.0, 1,
               (("stops", "at $r = r_q^2/r_s$ the surface lies level, and nearer the singularity the circles grow "
                          "faster than the distance out to them"),
                ("join", "the inner horizon $r = r_-$, the widest circle, where the slice runs on into the other "
                         "region inside $r_-$")),
               [(level, "chartedge", None), (0.3, "r", None), (rm, "horizon", "$r = r_-$")], size)
    lo.z = lo.z - lo.z[-1]
    hi = Piece("other_inside", "sheet2", sl, level, rm, -lo.z[0], -1,
               (("stops", "at $r = r_q^2/r_s$"), ("join", "the inner horizon $r = r_-$")),
               [(level, "chartedge", None), (0.3, "r2", None)], size)
    for p in (lo, hi):
        ck.isometry(f"Reissner-Nordstrom inside, {p.id}", p)
    ck.join("Reissner-Nordstrom inside, the two sides at r-", lo, rm, hi, rm)
    inside = Surface([lo, hi])
    fig = figure_of([inside], {"sheet": "cover"}, size, Camera(-90, 22))
    ring_label(fig, [0, 0, 0], rm, 0.0, "$r = r_-$", dx=10)
    ring_label(fig, [0, 0, 0], *lo.at(level), "$r_q^2/r_s$", side=-1)
    fig.legend("fill", "cover", "the region $r < r_-$ that $t$ and $r$ cover")
    fig.legend("line", "r", "$r$ constant, at $0.3\\,r_s$")
    fig.legend("line", "r2", "the same radius in the other region inside $r_-$")
    fig.legend("line", "horizon", "the widest circle $r = r_-$, where the slice crosses the inner horizon")
    fig.legend("line", "chartedge", "$r = r_q^2/r_s$, where the surface lies level and the drawing stops")
    fig.legend("line", "meridian", "$\\phi$ constant, every $15°$")
    views.append(view("inside", "Inside $r_-$", "$r_s$", [inside], fig.done(), settings=settings,
                      stops=["Nearer the singularity than $r = r_q^2/r_s$, $g_{rr} < 1$: the circles grow faster "
                             "than the distance out to them, and no surface in flat space carries that part of the "
                             "slice.", between]))
    return views


def kerr_family(ck, src, metric_id, name, params, ergo):
    """The equatorial slice of constant Boyer-Lindquist t outside r+, which has no cross term,
    g_tphi dropping out at constant t: rho = sqrt(g_phiphi), the circumference radius, and the
    two sheets through the bifurcation sphere, with the ergosphere's edge, where g_tt = 0 on
    the equator, marked."""
    sl = Slice(src, metric_id, "boyer_lindquist", "r", "\\phi", {"t": 0, **EQUATOR}, params)
    rp = sl.horizons()[0]
    top, radii = 8.0, (3.0, 4.0, 5.0, 6.0, 7.0)
    size = 2 * float(sl.rho_at(top))
    near, far = two_sheets(ck, name, sl, rp, top, radii, size,
                           [(rp, "horizon", "$r = r_+$"), (ergo, "ergo", None)])
    surface = Surface([near, far])
    fig = figure_of([surface], {"sheet": "cover"}, size)
    ring_label(fig, [0, 0, 0], *near.at(rp), "$r = r_+$", dx=14)
    ring_label(fig, [0, 0, 0], *near.at(ergo), "$r_E$")
    ring_label(fig, [0, 0, 0], *near.at(top), "$8\\,GM/c^2$")
    return sl, rp, near, surface, fig


def kerr(ck, src):
    """a = 0.9 GM/c^2, as the conformal diagram draws Kerr. On the equator the circumference
    radius at r+ is (r+^2 + a^2)/r+ = 2GM/c^2 whatever the spin, which is checked, and the
    ergosphere's edge is r = 2GM/c^2."""
    sl, rp, near, surface, fig = kerr_family(ck, src, "kerr", "Kerr", {"G": 1, "M": 1, "a": "9/10"}, 2.0)
    ck.add("Kerr: the throat's circumference radius is 2GM/c^2", abs(near.at(rp)[0] - 2.0), 1e-6)
    fig.legend("fill", "cover", "the exterior $r > r_+$ that $t$ and $r$ cover")
    fig.legend("line", "r", "$r$ constant, at $3$, $4$, $5$, $6$, $7$ and $8\\,GM/c^2$")
    fig.legend("line", "r2", "the same radii on the other exterior")
    fig.legend("line", "horizon", "the throat $r = r_+$, where the slice crosses the horizon")
    fig.legend("line", "ergo", "the edge of the ergosphere, $r_E = 2GM/c^2$ on the equator")
    fig.legend("line", "meridian", "$\\phi$ constant, every $15°$")
    return [view("equator", "The equator", "$GM/c^2$", [surface], fig.done(),
                 settings="$G = c = M = 1$, so that $GM/c^2$ is the unit of every length, and $a = 0.9\\,GM/c^2$, "
                          f"so that $r_+ = {rp:.3f}\\,GM/c^2$.")]


def kerr_newman(ck, src):
    """a = 0.6 GM/c^2 and r_Q = 0.5 GM/c^2, as the conformal diagram draws Kerr-Newman; the
    ergosphere's edge on the equator is where g_tt = 0, r^2 - 2GMr/c^2 + r_Q^2 = 0."""
    ergo = 1 + math.sqrt(0.75)
    sl, rp, near, surface, fig = kerr_family(ck, src, "kerr_newman", "Kerr-Newman",
                                             {"G": 1, "M": 1, "a": "3/5", "r_Q": "1/2"}, ergo)
    ck.add("Kerr-Newman: the throat's circumference radius is 2GM/c^2 - r_Q^2/r+", abs(near.at(rp)[0] - (2 - 0.25 / rp)), 1e-6)
    fig.legend("fill", "cover", "the exterior $r > r_+$ that $t$ and $r$ cover")
    fig.legend("line", "r", "$r$ constant, at $3$, $4$, $5$, $6$, $7$ and $8\\,GM/c^2$")
    fig.legend("line", "r2", "the same radii on the other exterior")
    fig.legend("line", "horizon", "the throat $r = r_+$, where the slice crosses the outer horizon")
    fig.legend("line", "ergo", "the edge of the ergosphere, $r_E = GM/c^2 + \\sqrt{(GM/c^2)^2 - r_Q^2}$ on the equator")
    fig.legend("line", "meridian", "$\\phi$ constant, every $15°$")
    return [view("equator", "The equator", "$GM/c^2$", [surface], fig.done(),
                 settings="$G = c = M = 1$, so that $GM/c^2$ is the unit of every length, $a = 0.6\\,GM/c^2$ and "
                          f"$r_Q = 0.5\\,GM/c^2$, so that $r_+ = {rp:.3f}\\,GM/c^2$ and $r_E = {ergo:.3f}\\,GM/c^2$.")]


def de_sitter(ck, src):
    """The static slice t = 0 at Lambda = 3, so that l = sqrt(3/Lambda) = 1: g_rr = 1/(1 - r^2),
    the metric of a sphere of radius l, z = -sqrt(1 - r^2) on the hemisphere the static chart
    covers out to its horizon r = l, and the same turned over on the antipodal observer's patch,
    which the slice runs on into through the horizon's bifurcation sphere. It is the waist of
    the hyperboloid, the smallest slice of the closed slicing. The flat slicing's slices are
    flat, which is checked and stated."""
    sl = Slice(src, "de_sitter", "static_spherical", "r", "\\phi", {"t": 0, **EQUATOR}, {"Lambda": 3})
    horizon = sl.horizons()[0]
    flat_slices(ck, src, "de_sitter", "flat_slicing")
    size = 2.0
    near = Piece("near", "sheet", sl, 0.0, horizon, -1.0, 1,
                 (("axis", "the observer at $r = 0$, the pole of its hemisphere"),
                  ("join", "the horizon $r = \\ell$, the equator, where the antipodal observer's patch begins")),
                 [(0.5, "r", None), (math.sqrt(3) / 2, "r", None), (horizon, "horizon", "$r = \\ell$")], size)
    far = Piece("far", "sheet2", sl, 0.0, horizon, 1.0, -1,
                (("axis", "the antipodal observer"), ("join", "the horizon")),
                [(0.5, "r2", None), (math.sqrt(3) / 2, "r2", None)], size)
    for p in (near, far):
        ck.isometry(f"de Sitter, the {p.id} hemisphere", p)
        ck.radius(f"de Sitter, the {p.id} hemisphere, the sphere rho = r", p, lambda r: r, size)
        ck.form(f"de Sitter, the {p.id} hemisphere, the sphere z = -+sqrt(l^2 - r^2)", p,
                lambda r, s=p.sense: -s * np.sqrt(np.maximum(1 - r * r, 0)), size)
    ck.join("de Sitter, the hemispheres meet at the horizon", near, horizon, far, horizon)
    surface = Surface([near, far])
    fig = figure_of([surface], {"sheet": "cover"}, size)
    ring_label(fig, [0, 0, 0], *near.at(horizon), "$r = \\ell$", dx=10)
    fig.legend("fill", "cover", "the static patch $r < \\ell$ that $t$ and $r$ cover")
    fig.legend("line", "r", "$r$ constant, at $\\ell/2$ and $\\sqrt{3}\\,\\ell/2$")
    fig.legend("line", "r2", "the same radii in the antipodal observer's patch")
    fig.legend("line", "horizon", "the horizon $r = \\ell$, the equator, where the slice crosses from one patch to the other")
    fig.legend("line", "meridian", "$\\phi$ constant, every $15°$")
    return [view("static", "The static patch", "$\\ell$", [surface], fig.done(),
                 settings="$\\Lambda = 3$, so that $\\ell = \\sqrt{3/\\Lambda} = 1$, the unit of every length.",
                 stops=["Every slice of constant $t$ of the flat slicing is flat, $e^{2Ht}(dx^2 + dy^2 + dz^2)$ "
                        "being Euclidean space scaled by $e^{Ht}$, so its equator is a plane."])]


def einstein_static(ck, src):
    """The equator of one moment of the hyperspherical chart at R = 1: g_chichi = 1 and
    g_phiphi = sin^2 chi, the round sphere of radius R, rho = sin chi and z = -cos chi from the
    pole chi = 0 to the antipode chi = pi, the same at every moment. The areal chart, r = R sin chi,
    and Einstein's projection cover the hemisphere chi < pi/2 and end at its equator, which is
    marked. The areal chart's own slice, g_rr = 1/(1 - r^2), is checked to be the same hemisphere."""
    sl = Slice(src, "einstein_static", "hyperspherical", "\\chi", "\\phi", {"t": 0, **EQUATOR}, {"R": 1})
    size = 2.0
    sphere = Piece("sphere", "sheet", sl, 0.0, math.pi, -1.0, 1,
                   (("axis", "the pole $\\chi = 0$"), ("axis", "the antipode $\\chi = \\pi$")),
                   [(math.pi / 4, "r", None), (math.pi / 2, "chartedge", "$\\chi = \\pi/2$"),
                    (3 * math.pi / 4, "r", None)], size)
    ck.isometry("Einstein static, the sphere", sphere)
    ck.radius("Einstein static, the sphere rho = R sin chi", sphere, np.sin, size)
    ck.form("Einstein static, the sphere z = -R cos chi", sphere, lambda chi: -np.cos(chi), size)
    areal = Slice(src, "einstein_static", "static_areal", "r", "\\phi", {"t": 0, **EQUATOR}, {"R": 1})
    near = Piece("near", "sheet", areal, 0.0, areal.horizons()[0], -1.0, 1, size=size)
    ck.isometry("Einstein static, the areal hemisphere", near)
    ck.form("Einstein static, the areal hemisphere z = -sqrt(R^2 - r^2)", near,
            lambda r: -np.sqrt(np.maximum(1 - r * r, 0)), size)
    surface = Surface([sphere])
    fig = figure_of([surface], {"sheet": "cover"}, size)
    ring_label(fig, [0, 0, 0], *sphere.at(math.pi / 2), "$\\chi = \\pi/2$", dx=10)
    fig.legend("fill", "cover", "the whole moment, a sphere of radius $R$")
    fig.legend("line", "r", "$\\chi$ constant, at $\\pi/4$ and $3\\pi/4$")
    fig.legend("line", "chartedge", "the equator $\\chi = \\pi/2$, $r = R$, where the areal chart and Einstein's end")
    fig.legend("line", "meridian", "$\\phi$ constant, every $15°$")
    # The moment is t = 0, eta = ct/R = 0, the slice the conformal diagram draws. While the
    # conformal diagram shows Minkowski space, de Sitter space or anti-de Sitter space in its
    # strip, the sphere is shaded where that spacetime's own map covers it at eta = 0, and the
    # map's reach at other times is checked against Hawking and Ellis's regions: chi < pi - |eta|,
    # the whole sphere for |eta| < pi/2, and chi < pi/2 at every eta.
    hawking_ellis = {
        "minkowski": lambda eta: (0.0, math.pi - abs(eta)) if abs(eta) < math.pi else None,
        "de_sitter": lambda eta: (0.0, math.pi) if abs(eta) < math.pi / 2 else None,
        "anti_de_sitter": lambda eta: (0.0, math.pi / 2),
    }
    words = {
        "minkowski": ("Minkowski space at $t = 0$, the sphere but its antipode $i^0$",
                      "The part of the sphere Minkowski space covers (pink, $\\eta = 0$): all of it but the antipode, "
                      "spatial infinity $i^0$, and the cap $\\chi < \\pi - |\\eta|$ at any other $\\eta$."),
        "de_sitter": ("de Sitter space at its waist, the whole sphere",
                      "The part of the sphere de Sitter space covers (pink, $\\eta = 0$): the whole sphere, at every "
                      "$\\eta$ from $-\\pi/2$ to $\\pi/2$."),
        "anti_de_sitter": ("anti-de Sitter space at $t = 0$, the hemisphere $\\chi < \\pi/2$",
                           "The part of the sphere anti-de Sitter space covers (pink, $\\eta = 0$): the hemisphere "
                           "$\\chi < \\pi/2$, the same at every $\\eta$."),
    }
    for vid, (name, fmap, reach) in conformal.ESU_EMBEDDED.items():
        reaches = [(eta, conformal.covered(fmap, reach, eta)) for eta in conformal.ESU_ETAS]
        err = 0.0
        for eta, got in reaches:
            want = hawking_ellis[vid](eta)
            err = max(err, 0.0 if got is None and want is None else math.inf if got is None or want is None
                      else max(abs(got[0] - want[0]), abs(got[1] - want[1])))
        ck.add(f"Einstein static: the part of the sphere {name} covers is Hawking and Ellis's, at every eta checked",
               err, 1e-6)
        lo, hi = conformal.covered(fmap, reach, 0.0)
        # Each end is a circle of the profile: the pole, the equator or the antipode.
        lo, hi = (float(sphere.x[int(np.argmin(np.abs(sphere.x - x)))]) for x in (lo, hi))
        ck.add(f"Einstein static: the part of the sphere {name} covers at eta = 0 ends on circles of the profile",
               max(abs(a - b) for a, b in zip((lo, hi), conformal.covered(fmap, reach, 0.0))), 1e-6)
        fig.shade(vid, 0, sphere, lo, hi, 0.0, reaches, ["fill", "shade", words[vid][0]], words[vid][1])
    return [view("sphere", "The three sphere", "$R$", [surface], fig.done(),
                 settings="$R = 1$, the unit of every length.")]


def schwarzschild_de_sitter(ck, src):
    """Kottler's static slice t = 0 at r_s = 1 and Lambda = 1/5, so that the horizons, the two
    positive roots of the published g^rr, sit at r_h = 1.085 and r_c = 3.215. g_rr = 1/f with
    f = 1 - 1/r - r^2/15 < 1, so dz/dr = sqrt(1/f - 1) is real from r_h to r_c and diverges at
    both ends: the surface stands vertical at the throat r_h, its smallest circle, and at r_c,
    its widest. The slice of the maximal extension runs through the bifurcation sphere at r_c
    into the next static region, the same surface turned over, and through each throat into
    another, without end; one period, from a throat to the next, is drawn. z(r) is an elliptic
    integral and is drawn by quadrature, checked against the published metric."""
    sl = Slice(src, "schwarzschild_de_sitter", "static", "r", "\\phi", {"t": 0, **EQUATOR}, {"r_s": 1, "Lambda": "1/5"})
    rc, rh, rn = sl.horizons()
    ck.add("Schwarzschild-de Sitter: the horizons are the roots of 3r - 3 - r^3/5", 
           abs(rh - 1.0851996154371) + abs(rc - 3.2146274073952) + abs(rn + 4.2998270228323), 1e-12)
    size = 2 * rc
    radii = (1.5, 2.0, 2.5, 3.0)
    near = Piece("static", "sheet", sl, rh, rc, 0.0, 1,
                 (("throat", "the throat $r = r_h$, the bifurcation sphere of the black hole horizon, where the "
                             "slice runs on into another static region"),
                  ("join", "the widest circle $r = r_c$, the bifurcation sphere of the cosmological horizon, "
                           "where the slice runs on into the next static region")),
                 [(rh, "horizon", "$r = r_h$")] + [(r, "r", None) for r in radii] + [(rc, "horizon", "$r = r_c$")],
                 size)
    top = 2 * near.z[-1]
    far = Piece("next", "sheet2", sl, rh, rc, top, -1,
                (("throat", "the throat $r = r_h$ of the next black hole horizon"),
                 ("join", "the widest circle $r = r_c$")),
                [(rh, "horizon", None)] + [(r, "r2", None) for r in radii], size)
    for p in (near, far):
        ck.isometry(f"Schwarzschild-de Sitter, the {p.id} region", p)
        ck.radius(f"Schwarzschild-de Sitter, the {p.id} region, rho = r", p, lambda r: r, size)
    ck.join("Schwarzschild-de Sitter, the two static regions at r_c", near, rc, far, rc)
    surface = Surface([near, far])
    fig = figure_of([surface], {"sheet": "cover"}, size)
    ring_label(fig, [0, 0, 0], *near.at(rh), "$r = r_h$", dx=10)
    ring_label(fig, [0, 0, 0], *near.at(rc), "$r = r_c$", dx=10)
    fig.legend("fill", "cover", "the static region $r_h < r < r_c$ that $t$ and $r$ cover")
    fig.legend("line", "r", "$r$ constant, at $1.5$, $2$, $2.5$ and $3\\,r_s$")
    fig.legend("line", "r2", "the same radii in the next static region")
    fig.legend("line", "horizon", "the throats $r = r_h$ and the widest circle $r = r_c$, where the slice crosses the horizons")
    fig.legend("line", "meridian", "$\\phi$ constant, every $15°$")
    return [view("static", "Between the horizons", "$r_s$", [surface], fig.done(),
                 settings="$r_s = 1$, the unit of every length, and $\\Lambda = 0.2/r_s^2$, so that "
                          "$r_h = 1.085\\,r_s$ and $r_c = 3.215\\,r_s$.")]


def global_monopole(ck, src):
    """Two moments at Delta = 0.19, so that sqrt(1 - Delta) = 0.9. The monopole with no mass at
    its centre, the Barriola-Vilenkin chart's equator at t = 0: g_rr = 1 and g_phiphi =
    (1 - Delta) r^2, so rho = 0.9 r and dz/dr = sqrt(Delta), a cone of half angle arcsin 0.9
    whose apex is the centre, where the Kretschmann scalar diverges, the same cone as the cosmic
    string's at 4G mu/c^2 = 0.1. Letelier's black hole, the static chart's equator at t = 0 with
    r_s = 1: g_rr = 1/(1 - Delta - 1/r) and rho = r, so dz/dr = sqrt((Delta r + 1)/((1 - Delta) r - 1)),
    vertical at the throat r_h = 1/(1 - Delta) = 100/81 and tending to the cone's sqrt(Delta/(1 - Delta))
    far out; with w = sqrt((1 - Delta) r - 1) its closed form is
    z = (w sqrt(Delta w^2 + 1) + arsinh(sqrt(Delta) w)/sqrt(Delta))/(1 - Delta)^(3/2). The slice runs
    through the bifurcation sphere into the other exterior, the same surface turned over."""
    D, fold = 0.19, 0.9
    size, top = 6.0, 3.0
    cone = Slice(src, "global_monopole", "conical", "r", "\\phi", {"t": 0, **EQUATOR}, {"Delta": "19/100"})
    apex = Piece("cone", "sheet", cone, 0.0, top, 0.0, 1,
                 (("apex", "the centre $r = 0$, where the Kretschmann scalar diverges"),
                  ("edge", "the cone runs on to $r \\to \\infty$")),
                 [(r, "r", None) for r in (1.0, 2.0)] + [(top, "r", None)], size)
    ck.isometry("global monopole, the cone", apex)
    ck.form("global monopole, the cone z = sqrt(Delta) r", apex, lambda r: math.sqrt(D) * r, size)
    ck.radius("global monopole, the cone rho = sqrt(1 - Delta) r", apex, lambda r: fold * r, size)
    monopole = Surface([apex])
    fig = figure_of([monopole], {"sheet": "cover"}, size, Camera(-90, 20))
    ring_label(fig, [0, 0, 0], *apex.at(1.0), "$r = \\ell$", side=-1, clear=True)
    ring_label(fig, [0, 0, 0], *apex.at(top), "$3\\ell$", side=-1)
    fig.legend("fill", "cover", "the equatorial plane, a cone of half angle $\\arcsin\\sqrt{1 - \\Delta}$")
    fig.legend("line", "r", "$r$ constant, at $\\ell$, $2\\ell$ and $3\\ell$")
    fig.legend("line", "meridian", "$\\phi$ constant, every $15°$")
    first = view("monopole", "The monopole", "$\\ell$", [monopole], fig.done(),
                 settings="$\\Delta = 0.19$, so that $\\sqrt{1 - \\Delta} = 0.9$, and any length $\\ell$ as the "
                          "unit.")

    size, top = 12.0, 6.0
    sl = Slice(src, "global_monopole", "static", "r", "\\phi", {"t": 0, **EQUATOR}, {"Delta": "19/100", "r_s": 1})
    rh = sl.horizons()[0]
    ck.add("global monopole: the horizon is r_s/(1 - Delta) = 100/81", abs(rh - 100 / 81), 1e-12)
    radii = (2.0, 3.0, 4.0, 5.0)
    near = Piece("exterior", "sheet", sl, rh, top, 0.0, 1,
                 (("throat", "the throat $r = r_h$, the bifurcation sphere, where the other exterior begins"),
                  ("edge", "the surface runs on, opening into the cone of the monopole, to $r \\to \\infty$")),
                 [(rh, "horizon", "$r = r_h$")] + [(r, "r", None) for r in radii] + [(top, "r", None)], size)
    far = Piece("other_exterior", "sheet2", sl, rh, top, 0.0, -1,
                (("throat", "the throat $r = r_h$"),
                 ("edge", "the surface runs on, opening into the cone of the monopole, to $r \\to \\infty$")),
                [(r, "r2", None) for r in radii] + [(top, "r2", None)], size)

    def closed(r):
        w = np.sqrt(np.maximum((1 - D) * r - 1, 0))
        return (w * np.sqrt(D * w * w + 1) + np.arcsinh(math.sqrt(D) * w) / math.sqrt(D)) / (1 - D) ** 1.5
    for p in (near, far):
        ck.isometry(f"global monopole, Letelier's black hole, {p.id}", p)
        ck.form(f"global monopole, Letelier's black hole, {p.id}, the closed form of z", p,
                lambda r, s=p.sense: s * closed(r), size)
        ck.radius(f"global monopole, Letelier's black hole, {p.id}, rho = r", p, lambda r: r, size)
    ck.join("global monopole, Letelier's black hole, the two sheets at the throat", near, rh, far, rh)
    hole = Surface([near, far])
    fig = figure_of([hole], {"sheet": "cover"}, size)
    ring_label(fig, [0, 0, 0], *near.at(rh), "$r = r_h$", dx=14)
    ring_label(fig, [0, 0, 0], *near.at(3.0), "$3\\,r_s$")
    ring_label(fig, [0, 0, 0], *near.at(top), "$6\\,r_s$")
    fig.legend("fill", "cover", "the exterior $r > r_h$ that $t$ and $r$ cover")
    fig.legend("line", "r", "$r$ constant, at $2$, $3$, $4$, $5$ and $6\\,r_s$")
    fig.legend("line", "r2", "the same radii on the other exterior")
    fig.legend("line", "horizon", "the throat $r = r_h$, where the slice crosses the horizon")
    fig.legend("line", "meridian", "$\\phi$ constant, every $15°$")
    second = view("black_hole", "Letelier's black hole", "$r_s$", [hole], fig.done(),
                  settings="$r_s = 1$, the unit of every length, and $\\Delta = 0.19$, so that "
                           "$r_h = r_s/(1 - \\Delta) = 1.235\\,r_s$.")
    return [first, second]


def vaidya(ck, src):
    """The imploding shell of radiation the conformal diagram draws: the ingoing chart with m = 0
    for v < 0 and M for v > 0, r_s = 2GM/c^2 = 1. A slice of constant v is null, so the moments
    are slices of constant v - r = T, spacelike everywhere, inside the horizon too, on which the
    published metric pulls back to (1 + 2Gm/c^2r) dr^2 + r^2 dphi^2 with v = T + r. Inside the
    shell, at r < -T, m = 0 and the surface is a flat disc, which is checked; outside, dz/dr =
    sqrt(r_s/r), z = 2 sqrt(r_s r), Flamm's paraboloid moved in by r_s, which reaches the axis in
    a spike at the singularity once the shell has gone. The shell folds the surface where they
    meet. The event horizon, u = v - 2r = -2 r_s inside and r = r_s outside, is r = T + 2 on the
    disc while the shell is outside r_s and r_s after."""
    top = 4.0
    size = 2 * top

    def slice_at(T, m):
        return Slice(src, "vaidya", "eddington_finkelstein_ingoing", "r", "\\phi", {"theta": "pi/2"}, {"G": 1},
                     {"m": m}, along={"v": f"{T} + r"})
    rim = ("edge", "the surface runs on, as Flamm's paraboloid moved in by $r_s$, to $r \\to \\infty$")

    def moment(T):
        where = f"Vaidya, v - r = {T:g}"
        hole = slice_at(T, "1/2")
        rings = [(r, "r", None) for r in (2.0, 3.0, top)]
        if T < 0:
            R = -T
            flat = slice_at(T, "0")
            ck.plane(f"{where}, inside the shell", flat, np.linspace(1e-3, R, 200))
            inside = Piece("inside", "sheet", flat, 0.0, R, 0.0, 1,
                           (("axis", "the centre $r = 0$, where space is flat until the shell arrives"),
                            ("crease", "the shell of radiation, where the surface folds")),
                           [(R, "surface", None)] + ([(T + 2, "horizon", None)] if 0 < T + 2 < R else []), size)
            outside = Piece("outside", "sheet", hole, R, top, 0.0, 1, (("crease", "the shell"), rim),
                            ([(1.0, "horizon", None)] if R < 1 else []) + [m for m in rings if m[0] > R], size)
            pieces = [inside, outside]
            # The rim of the drawing, r = 4 r_s, stands at z = 0 at every moment.
            shift = -outside.z[-1]
            inside.z, outside.z = inside.z + shift, outside.z + shift
            ck.add(f"{where}, the two sides meet at the shell: one point",
                   float(np.max(np.abs(np.array(inside.at(R)) - outside.at(R)))), JOIN)
            ck.form(f"{where}, outside the shell z = 2 sqrt(r_s r) - 2 sqrt(r_s R)", outside,
                    lambda r, R=R, shift=shift: 2 * (np.sqrt(r) - np.sqrt(R)) + shift, size)
        else:
            whole = Piece("whole", "sheet", hole, 0.0, top, 0.0, 1,
                          (("apex", "the singularity $r = 0$, where the surface closes in a spike"), rim),
                          [(1.0, "horizon", None)] + rings, size)
            shift = -whole.z[-1]
            whole.z = whole.z + shift
            pieces = [whole]
            ck.form(f"{where}, z = 2 sqrt(r_s r)", whole, lambda r, shift=shift: 2 * np.sqrt(r) + shift, size)
        for p in pieces:
            ck.isometry(f"{where}, {p.id}", p)
        return Surface(pieces, label=f"$v - r = {T:g}\\,r_s$", time=T)

    # A frame every r_s/16 of v - r, which holds each moment of the flat views.
    frames = [moment(-3.0 + k / 16) for k in range(65)]
    surfaces = [frames[k] for k in (0, 24, 40, 64)]
    fig = movie_figure(frames, {"sheet": "cover"}, size)
    fig.legend("fill", "cover", "the slice of constant $v - r$, which $v$ and $r$ cover")
    fig.legend("line", "r", "$r$ constant, at $2$, $3$ and $4\\,r_s$")
    fig.legend("line", "surface", "the shell of radiation, where the surface folds")
    fig.legend("line", "horizon", "the event horizon, which forms at the centre at $v = -2\\,r_s$ and grows through flat "
                                  "space to meet the shell at $r_s$")
    fig.legend("line", "meridian", "$\\phi$ constant, every $30°$")
    return [view("shell", "The falling shell", "$r_s$", surfaces, fig.done(),
                 movie=movie(frames, "$v - r$", [f.time for f in frames]),
                 settings="$r_s = 2GM/c^2 = 1$, the unit of every length; each moment is a slice of constant $v - r$.",
                 input="An imploding shell of radiation, $m = 0$ for $v < 0$ and $m = M$ for $v > 0$, as in "
                       "the conformal diagram.")]


def einstein(src, metric_id, system_id, index):
    """A published G^i_i with c = 1, every declared function and derivative a plain symbol named
    as R, R_t, R_tt, R_tr or E_r, as a numpy function of those symbols by keyword."""
    _, entry, R = nr.load(metric_id, system_id)
    src.note(metric_id, system_id, ["einstein_tensor"])
    ul = entry["einstein_tensor"]["variants"]["ul"]["nonzero"]
    e = R(next(c["value"] for c in ul if c["indices"] == [index, index])).subs(R.c, 1)
    named = {}
    for d in sorted(e.atoms(sp.Derivative), key=lambda d: -len(d.variables)):
        name = f"{d.expr.func}_{''.join(str(v) for v in d.variables)}"
        named[name] = sp.Symbol(name)
        e = e.subs(d, named[name])
    for f in e.atoms(sp.core.function.AppliedUndef):
        named[str(f.func)] = sp.Symbol(str(f.func))
        e = e.subs(f, named[str(f.func)])
    order = sorted(named)
    fn = sp.lambdify([named[k] for k in order], e, "numpy")
    return lambda **v: fn(*[v[k] for k in order])


class RestCloud:
    """Dust released from rest at t = 0 with R(r, 0) = r, G = c = 1, described by k(r) = 2M(r)/r^3,
    M the mass inside the shell r: each shell falls on its own cycloid, R = (r/2)(1 + cos eta)
    and t = (eta + sin eta)/(2 sqrt k), with 1 + 2E = 1 - k r^2 and E = -M/r, the energy of a
    shell at rest at R = r. At fixed t, d eta/dr = (eta + sin eta) k'/(2k (1 + cos eta)), which
    gives dR/dr; dR/dt = -r sqrt(k) tan(eta/2) and d^2R/dt^2 = -M/R^2 follow from the
    parametrisation, and are checked against the published field equations rather than assumed.
    Outside the dust M is constant and the slices of constant t are Novikov's, the moment of
    clocks falling freely from rest, which run on through the horizon."""

    def __init__(self, k, dk):
        self.k, self.dk = k, dk

    def eta(self, r, t):
        """eta + sin eta = 2 t sqrt(k(r)) by bisection, which is exact to double precision."""
        target = 2 * t * np.sqrt(self.k(r))
        if np.any(target > math.pi):
            raise AssertionError(f"a shell has reached R = 0 by t = {t}")
        lo, hi = np.zeros_like(target), np.full_like(target, math.pi)
        for _ in range(64):
            mid = 0.5 * (lo + hi)
            below = mid + np.sin(mid) < target
            lo, hi = np.where(below, mid, lo), np.where(below, hi, mid)
        return 0.5 * (lo + hi)

    def fields(self, r, t):
        r = np.asarray(r, dtype=float)
        k, dk = self.k(r), self.dk(r)
        e = self.eta(r, t)
        c, sn = np.cos(e), np.sin(e)
        with np.errstate(invalid="ignore", divide="ignore"):
            de = np.where(e > 0, (e + sn) * dk / (2 * k * (1 + c)), 0.0)
        half = np.cos(e / 2)
        return {"R": r * (1 + c) / 2, "R_r": (1 + c) / 2 - r * sn * de / 2,
                "R_t": -r * np.sqrt(k) * np.tan(e / 2), "R_tt": -k * r / (2 * half ** 4),
                "R_tr": -(np.sqrt(k) + r * dk / (2 * np.sqrt(k))) * np.tan(e / 2) - r * np.sqrt(k) * de / (2 * half ** 2),
                "E": -k * r * r / 2, "E_r": -(dk * r * r + 2 * k * r) / 2,
                "M": k * r ** 3 / 2, "M_r": (dk * r ** 3 + 3 * k * r * r) / 2, "eta": e}

    def check(self, ck, src, name, rs, ts):
        """The published G^r_r vanishes, which is dust, and G^t_t = -2 M'/(R^2 R') is its density."""
        Grr = einstein(src, "tolman_bondi", "comoving_synchronous", "r")
        Gtt = einstein(src, "tolman_bondi", "comoving_synchronous", "t")
        worst_r = worst_t = 0.0
        for t in ts:
            f = self.fields(rs, t)
            v = {"R": f["R"], "R_t": f["R_t"], "R_tt": f["R_tt"], "R_r": f["R_r"], "R_tr": f["R_tr"],
                 "E": f["E"], "E_r": f["E_r"]}
            scale = f["M"] / f["R"] ** 3 + 1e-300
            worst_r = max(worst_r, float(np.max(np.abs(Grr(**v)) / scale)))
            density = -2 * f["M_r"] / (f["R"] ** 2 * f["R_r"])
            worst_t = max(worst_t, float(np.max(np.abs(Gtt(**v) - density) / scale)))
        ck.add(f"{name}: released from rest, the published G^r_r vanishes", worst_r, 1e-9)
        ck.add(f"{name}: the published G^t_t is -2M'/(R^2 R'), the declared density", worst_t, 1e-9)

    def slice(self, src, t, E):
        """Tolman-Bondi's comoving chart at time t, with R and dR/dr from the cycloids and E closed."""
        def value(key):
            return lambda x: self.fields(np.atleast_1d(x), t)[key].reshape(np.shape(x))
        return Slice(src, "tolman_bondi", "comoving_synchronous", "r", "\\phi", {"t": repr(t), **EQUATOR},
                     functions={"E": E}, numeric={"R": (value("R"), value("R_r"))})

    def horizon(self, t, lo, hi):
        """The outermost shell in (lo, hi) with R = 2M, the apparent horizon, or None."""
        r = np.linspace(lo, hi, 4001)[1:]
        f = self.fields(r, t)
        g = f["R"] - 2 * f["M"]
        cross = np.flatnonzero((g[:-1] < 0) & (g[1:] >= 0))
        if not cross.size:
            return None
        a, b = r[cross[-1]], r[cross[-1] + 1]
        for _ in range(64):
            m = 0.5 * (a + b)
            fm = self.fields(np.array([m]), t)
            a, b = (m, b) if fm["R"][0] - 2 * fm["M"][0] < 0 else (a, m)
        return 0.5 * (a + b)


def oppenheimer_snyder(ck, src):
    """Dust released from rest at R0 = 2 r_s, chi0 = pi/4 and a_m = 2 sqrt(2) r_s, as the
    conformal diagram draws the collapse, at four moments of its proper time tau. Inside, the
    published interior chart at a = (a_m/2)(1 + cos eta), tau = (a_m/2)(eta + sin eta), which is
    checked to make the published G^chi_chi vanish: a cap of a sphere of radius a out to chi0.
    Outside, the moment of constant tau carries on as Novikov's slice, the moment of clocks
    released from rest at t = 0 at every radius, read from Tolman-Bondi's comoving chart with
    no dust, k = r_s/r^3; a slice of constant Schwarzschild t would meet the dust at an angle
    after the release and cannot reach it once it is inside r_s. The two meet with one tangent,
    1 + 2E = 1 - r_s/R0 = cos^2 chi0 on both sides, which is checked; at the release the outside
    is Flamm's paraboloid."""
    R0, chi0 = 2.0, math.pi / 4
    am = R0 / math.sin(chi0)
    cloud = RestCloud(lambda r: np.where(r < R0, 1 / R0 ** 3, 1 / np.maximum(r, R0) ** 3),
                      lambda r: np.where(r < R0, 0.0, -3 / np.maximum(r, R0) ** 4))
    # The interior's own field equation, 2 a a'' + a'^2 + 1 = 0 in tau, along the cycloid.
    Gcc = einstein(src, "oppenheimer_snyder", "interior_comoving", "\\chi")
    etas = np.linspace(0.05, 0.95 * math.pi, 50)
    half = np.cos(etas / 2)
    ck.add("Oppenheimer-Snyder: a = (a_m/2)(1 + cos eta) makes the published G^chi_chi vanish",
           float(np.max(np.abs(Gcc(a=am * half ** 2, a_tau=-np.tan(etas / 2), a_tautau=-1 / (2 * am * half ** 4)))
                        * (am * half ** 2) ** 2)), 1e-9)
    cloud.check(ck, src, "Oppenheimer-Snyder outside", np.linspace(R0, 8, 60), [0.0, 1.0, 2.0, 3.0, 4.0])
    top = 4.0
    size = 2 * top
    from scipy.optimize import brentq

    def moment(eta):
        a = am * (1 + math.cos(eta)) / 2
        tau = am * (eta + math.sin(eta)) / 2
        where = f"Oppenheimer-Snyder, eta = {eta / math.pi:.4f} pi"
        inner = Slice(src, "oppenheimer_snyder", "interior_comoving", "\\chi", "\\phi", {"tau": "0", **EQUATOR},
                      {"chi_0": "pi/4", "a_m": repr(am)}, {"a": repr(a)})
        outer = cloud.slice(src, tau, "-1/(2*r)")
        edge = cloud.horizon(tau, 0.0, top)
        dust_marks = [(chi0 / 3, "r", None), (2 * chi0 / 3, "r", None), (chi0, "surface", None)]
        out_marks = [(r, "r", None) for r in (3.0, top)]
        if edge is not None and edge > R0:
            out_marks.append((edge, "horizon", None))
        elif edge is not None:
            dust_marks.append((math.asin(edge / am), "horizon", None))
        ext = Piece("exterior", "sheet", outer, R0, top, 0.0, 1,
                    (("join", "the surface of the dust"), ("edge", "the slice runs on to $r \\to \\infty$")), out_marks, size)
        dust = Piece("dust", "star", inner, 0.0, chi0, 0.0, 1,
                     (("axis", "the centre $\\chi = 0$, where the cap is smooth"), ("join", "the surface $\\chi = \\chi_0$")),
                     dust_marks, size)
        dust.z = dust.z + (ext.z[0] - dust.z[-1])
        # The rim of the drawing, the clocks released at r = 4 r_s, stands at z = 0 at every
        # moment, and the dust sinks below it.
        shift = -ext.z[-1]
        ext.z, dust.z = ext.z + shift, dust.z + shift
        ck.isometry(f"{where}, the dust", dust)
        ck.isometry(f"{where}, outside", ext)
        ck.join(f"{where}, the dust meets the outside", dust, chi0, ext, R0)
        ck.form(f"{where}, the dust is a cap of a sphere of radius a", dust,
                lambda c, a=a, z0=dust.z[0]: z0 + a * (1 - np.cos(c)), size)
        if eta == 0:
            ck.form(f"{where}, outside it is Flamm's paraboloid", ext, lambda r: 2 * np.sqrt(r - 1) - 2 + shift, size)
        return Surface([dust, ext], label=f"$c\\tau = {tau:.2f}\\,r_s$", time=tau)

    # The movie runs through the moments at a steady proper time of the dust, a frame about every
    # 0.07 r_s of c tau.
    etas = [f * math.pi for f in (0.0, 0.3, 0.6, 0.8)]
    taus, keys = movie_values([am * (e + math.sin(e)) / 2 for e in etas], 0.07)
    frames = [moment(etas[keys.index(i)] if i in keys else
                     brentq(lambda e, t=t: am * (e + math.sin(e)) / 2 - t, 0.0, math.pi, xtol=1e-15))
              for i, t in enumerate(taus)]
    surfaces = [frames[i] for i in keys]
    fig = movie_figure(frames, {"star": "star", "sheet": "cover"}, size)
    fig.legend("fill", "star", "the dust, a cap of a sphere of radius $a(\\tau)$, which $\\tau$ and $\\chi$ cover")
    fig.legend("fill", "cover", "outside it, the moment of clocks released from rest with the dust")
    fig.legend("line", "r", "$\\chi$ constant in the dust, at $\\chi_0/3$ and $2\\chi_0/3$, and outside the clocks "
                            "released at $3$ and $4\\,r_s$")
    fig.legend("line", "surface", "the surface of the dust, $\\chi = \\chi_0$")
    fig.legend("line", "horizon", "the apparent horizon, $R = 2GM/c^2$ for the mass inside it")
    fig.legend("line", "meridian", "$\\phi$ constant, every $30°$")
    return [view("collapse", "The collapse", "$r_s$", surfaces, fig.done(),
                 movie=movie(frames, "$c\\tau$", [f.time for f in frames]),
                 settings="$R_0 = 2\\,r_s$, so that $\\chi_0 = \\pi/4$ and $a_m = 2\\sqrt{2}\\,r_s$, with $r_s = 1$ the "
                          "unit of every length; the moments are the dust's proper time $\\tau$ since the release.",
                 input="Outside the dust, the slices of Tolman-Bondi's comoving chart with no dust in it, $1 + 2E = "
                       "1 - r_s/r$, each shell of clocks released from rest at $R = r$ when the dust is.")]


VERTICAL = 2   # the Tolman-Bondi cloud is drawn this many times taller than its embedding


def tolman_bondi(ck, src):
    """The cloud the spacetime diagram draws, density falling as 1 - r^2/r_b^2 to zero at r_b and
    2GM(r_b)/c^2 = r_b/2, so k = 2M/r^3 = (5 - 3r^2)/4 inside and 1/(2r^3) outside, r_b = 1,
    but released from rest, E = -M/r, rather than marginally bound: with E = 0 the slice has
    g_rr = (dR/dr)^2 and is a plane, which is checked for the spacetime diagram's own cloud and
    stated. Four moments of t until just before the centre is crushed, at t = pi/sqrt(5).
    Every surface is checked as the metric gives it and then drawn with its heights doubled,
    VERTICAL, which the caption states: the cloud rises by a third of its width and would read
    on the page as a flat disc."""
    cloud = RestCloud(lambda r: np.where(r < 1, (5 - 3 * r * r) / 4, 1 / (2 * np.maximum(r, 1.0) ** 3)),
                      lambda r: np.where(r < 1, -1.5 * r, -1.5 / np.maximum(r, 1.0) ** 4))
    cloud.check(ck, src, "Tolman-Bondi", np.linspace(0.02, 4, 80), [0.0, 0.5, 1.0, 1.3])
    E = "Piecewise((-(5 - 3*r**2)*r**2/8, r < 1), (-1/(4*r), True))"
    for t in (0.0, 0.3):
        flat = Slice(src, "tolman_bondi", "comoving_synchronous", "r", "\\phi", {"t": repr(t), **EQUATOR},
                     functions={"E": "0", "R": nr._TB_R})
        ck.plane(f"Tolman-Bondi, marginally bound, t = {t:g}", flat, np.linspace(0.02, 3, 300))
    top = 2.5
    size = 2 * top
    surfaces = []
    for t in (0.0, 0.6, 1.0, 1.3):
        sl = cloud.slice(src, t, E)
        where = f"Tolman-Bondi, t = {t:g}"
        edge = cloud.horizon(t, 0.0, top)
        horizon = [(edge, "horizon", None)] if edge is not None else []
        dust = Piece("cloud", "star", sl, 0.0, 1.0, 0.0, 1,
                     (("axis", "the centre $r = 0$, where the surface is smooth"), ("join", "the surface $r = r_b$")),
                     [(1 / 3, "r", None), (2 / 3, "r", None), (1.0, "surface", None)] + [h for h in horizon if h[0] < 1], size)
        ext = Piece("exterior", "sheet", sl, 1.0, top, dust.z[-1], 1,
                    (("join", "the surface $r = r_b$"), ("edge", "the slice runs on to $r \\to \\infty$")),
                    [(r, "r", None) for r in (1.5, 2.0, top)] + [h for h in horizon if h[0] > 1], size)
        ck.isometry(f"{where}, the cloud", dust)
        ck.isometry(f"{where}, outside", ext)
        ck.join(f"{where}, the cloud meets the outside", dust, 1.0, ext, 1.0)
        if t == 0:
            ck.form(f"{where}, outside it is Flamm's paraboloid", ext,
                    lambda r, z1=ext.z[0]: z1 + 2 * np.sqrt(0.5 * (r - 0.5)) - 2 * np.sqrt(0.25), size)
        for p in (dust, ext):
            p.z = VERTICAL * p.z
        surfaces.append(Surface([dust, ext], label=f"$ct = {t:g}\\,r_b$", time=t))
    fig = sequence_figure(surfaces, {"star": "star", "sheet": "cover"}, size, columns=2)
    fig.legend("fill", "star", "the cloud, $r < r_b$")
    fig.legend("fill", "cover", "outside it, the moment of clocks released from rest with the cloud")
    fig.legend("line", "r", "$r$ constant: the shells $r_b/3$ and $2r_b/3$ of the cloud, and outside it the clocks "
                            "released at $1.5$, $2$ and $2.5\\,r_b$")
    fig.legend("line", "surface", "the surface of the cloud, $r = r_b$")
    if any(ring["class"] == "horizon" for s in surfaces for ring in s.rings()):
        fig.legend("line", "horizon", "the apparent horizon, $R = 2GM(r)/c^2$ for the mass inside it")
    fig.legend("line", "meridian", "$\\phi$ constant, every $30°$")
    return [view("cloud", "The collapsing cloud", "$r_b$", surfaces, fig.done(), vertical=VERTICAL,
                 settings="$r_b = 1$, the unit of every length, and $2GM/c^2 = r_b/2$.",
                 input="The cloud of the spacetime diagram, its density falling as $1 - r^2/r_b^2$ to zero at $r_b$ "
                       "with $R(r, 0) = r$, but released from rest, $E = -GM(r)/c^2r$, each shell falling on its own "
                       "cycloid, checked to solve this spacetime's own $G^r{}_r = 0$ and to give its density.",
                 stops=["The spacetime diagram's cloud is marginally bound, $E = 0$, and then every slice of constant "
                        "$t$ is flat: $g_{rr} = (\\partial_rR)^2$ makes the distance between two shells the difference "
                        "of their areal radii, so the equator is a plane and the collapse lies in how the slices are "
                        "stacked."])]


def bertotti_robinson(ck, src):
    """AdS2 x S2 with one radius b = 1. The equator at one moment is a line of the first factor
    times a great circle of the second, b^2 dr^2/r^2 + b^2 dphi^2: rho = b, z = b ln r, a
    cylinder, flat, with r -> 0 and r -> infinity infinitely far. The sphere of theta and phi at
    one event of t and r is the second factor, drawn as the second view; the first factor is
    Lorentzian and has no surface in flat space."""
    sl = Slice(src, "bertotti_robinson", "static", "r", "\\phi", {"t": 0, **EQUATOR}, {"b": 1})
    size = 2.0
    tube = Piece("cylinder", "sheet", sl, math.exp(-1), math.exp(1), -1.0, 1,
                 (("edge", "the cylinder runs on for ever toward $r \\to 0$"),
                  ("edge", "the cylinder runs on for ever toward $r \\to \\infty$")),
                 [(math.exp(k / 2), "r", None) for k in (-2, -1, 0, 1, 2)], size)
    ck.isometry("Bertotti-Robinson, the equator", tube)
    ck.form("Bertotti-Robinson, the cylinder z = b ln r", tube, np.log, size)
    ck.radius("Bertotti-Robinson, the cylinder rho = b", tube, lambda r: np.ones_like(r), size)
    equator = Surface([tube])
    fig = figure_of([equator], {"sheet": "cover"}, size)
    ring_label(fig, [0, 0, 0], *tube.at(1.0), "$r = b$")
    ring_label(fig, [0, 0, 0], *tube.at(math.exp(1)), "$eb$")
    ring_label(fig, [0, 0, 0], *tube.at(math.exp(-1)), "$b/e$")
    fig.legend("fill", "cover", "the equator at one moment, which $t$ and $r$ cover")
    fig.legend("line", "r", "$r$ constant, at $e^{-1}$, $e^{-1/2}$, $1$, $e^{1/2}$ and $e$ times $b$, a step $b/2$ apart along the cylinder")
    fig.legend("line", "meridian", "$\\phi$ constant, every $15°$")
    settings = "$b = 1$, the unit of every length."
    views = [view("equator", "The equator", "$b$", [equator], fig.done(), settings=settings,
                  stops=["The factor of $t$ and $r$ is a two dimensional anti-de Sitter space, which is Lorentzian and "
                         "has no surface in flat space; its moment of constant $t$ is the line along the cylinder."])]
    sphere_slice = Slice(src, "bertotti_robinson", "static", "\\theta", "\\phi", {"t": 0, "r": "1"}, {"b": 1})
    size = 2.0
    ball = Piece("sphere", "sheet", sphere_slice, 0.0, math.pi, 0.0, 1,
                 (("axis", "the pole $\\theta = 0$"), ("axis", "the pole $\\theta = \\pi$")),
                 [(math.pi / 4, "r", None), (math.pi / 2, "r", None), (3 * math.pi / 4, "r", None)], size)
    ck.isometry("Bertotti-Robinson, the sphere", ball)
    ck.form("Bertotti-Robinson, the sphere z = b(1 - cos theta)", ball, lambda c: 1 - np.cos(c), size)
    ck.radius("Bertotti-Robinson, the sphere rho = b sin theta", ball, np.sin, size)
    sphere = Surface([ball])
    fig = figure_of([sphere], {"sheet": "cover"}, size)
    ring_label(fig, [0, 0, 0], *ball.at(math.pi / 2), "$\\theta = \\pi/2$")
    fig.legend("fill", "cover", "the sphere, which $\\theta$ and $\\phi$ cover but for its poles")
    fig.legend("line", "r", "$\\theta$ constant, at $\\pi/4$, $\\pi/2$ and $3\\pi/4$")
    fig.legend("line", "meridian", "$\\phi$ constant, every $15°$")
    views.append(view("sphere", "The sphere", "$b$", [sphere], fig.done(), settings=settings))
    return views


# Two of Majumdar and Papapetrou's holes, each of mass parameter m, the unit, at z = +-2m on the
# axis, as their spacetime diagrams declare.
MP_TWO_CYLINDRICAL = "1 + 1/sqrt(rho**2 + (z - 2)**2) + 1/sqrt(rho**2 + (z + 2)**2)"


def mp_throat(r):
    """z(r) of one hole's equator, m = 1, from r = 1/50: the quadrature of
    sqrt((1 + 1/r)^2 - 1) = sqrt(2r + 1)/r, which is 2w + ln((w - 1)/(w + 1)) with w = sqrt(2r + 1)."""
    def F(x):
        w = np.sqrt(2 * np.asarray(x, dtype=float) + 1)
        return 2 * w + np.log((w - 1) / (w + 1))
    return F(r) - F(MP_THROAT_LO)


MP_THROAT_LO = 1 / 50


def majumdar_papapetrou(ck, src):
    """Two moments t = 0, in two views. One hole, m = 1, in the isotropic chart: the equator's
    g_rr = (1 + 1/r)^2 and circumference radius rho = r + 1, so dz/dr = sqrt(2r + 1)/r, which
    grows as 1/r toward the horizon r = 0, where z falls as ln r without end while rho closes on
    m: an infinitely long throat, drawn from r = 1/50 out to 4m. Two holes at z = +-2m, in the
    cylindrical chart: the plane z = 0 midway between them, U = 1 + 2/sqrt(rho^2 + 4), a surface of
    revolution about the axis with rho_emb = rho U and dz/drho = sqrt(U^2 - (U + rho U')^2), which
    is real since U + rho U' = 1 + 8/(rho^2 + 4)^(3/2) lies between 0 and U; it is flat at the
    axis, where U' = 0, and drawn out to rho = 4m."""
    sl = Slice(src, "majumdar_papapetrou", "isotropic", "r", "\\phi", {"t": 0, **EQUATOR}, {"m": 1})
    lo, top = MP_THROAT_LO, 4.0
    size = 2 * (top + 1)
    radii = (0.1, 0.25, 0.5, 1.0, 2.0, 3.0)
    throat = Piece("one_hole", "sheet", sl, lo, top, 0.0, 1,
                   (("edge", "the throat runs on without end toward the horizon $r = 0$, its circles closing on "
                             "the circumference radius $m$"),
                    ("edge", "the surface runs on to $r \\to \\infty$")),
                   [(lo, "r", None)] + [(r, "r", None) for r in radii] + [(top, "r", None)], size)
    ck.isometry("Majumdar-Papapetrou, one hole", throat)
    ck.radius("Majumdar-Papapetrou, one hole, rho = r + m", throat, lambda r: r + 1, size)
    ck.form("Majumdar-Papapetrou, one hole, z = 2w + ln((w - 1)/(w + 1)), w = sqrt(2r/m + 1)", throat, mp_throat, size)
    one = Surface([throat])
    fig = figure_of([one], {"sheet": "cover"}, size)
    ring_label(fig, [0, 0, 0], *throat.at(lo), "$r = m/50$")
    ring_label(fig, [0, 0, 0], *throat.at(1.0), "$r = m$")
    ring_label(fig, [0, 0, 0], *throat.at(top), "$4m$")
    fig.legend("fill", "cover", "the equator at one moment, which $t$ and $r > 0$ cover")
    fig.legend("line", "r", "$r$ constant, at $1/50$, $1/10$, $1/4$, $1/2$, $1$, $2$, $3$ and $4$ times $m$")
    fig.legend("line", "meridian", "$\\phi$ constant, every $15°$")
    views = [view("one_hole", "One hole", "$m$", [one], fig.done(),
                  settings="$m = 1$, the unit of every length, and $U = 1 + m/r$.")]

    sl = Slice(src, "majumdar_papapetrou", "cylindrical", "\\rho", "\\phi", {"t": 0, "z": 0},
               functions={"U": MP_TWO_CYLINDRICAL})
    top = 4.0
    size = 2 * float(sl.rho_at(top))
    radii = (1.0, 2.0, 3.0)
    plane = Piece("two_holes", "sheet", sl, 0.0, top, 0.0, 1,
                  (("axis", "the axis $\\rho = 0$, midway between the holes"),
                   ("edge", "the surface runs on to $\\rho \\to \\infty$")),
                  [(r, "r", None) for r in radii] + [(top, "r", None)], size)
    ck.isometry("Majumdar-Papapetrou, two holes, the midplane", plane)
    ck.radius("Majumdar-Papapetrou, two holes, rho = rho U", plane,
              lambda r: r * (1 + 2 / np.sqrt(r * r + 4)), size)
    two = Surface([plane])
    fig = figure_of([two], {"sheet": "cover"}, size)
    ring_label(fig, [0, 0, 0], *plane.at(2.0), "$\\rho = 2m$")
    ring_label(fig, [0, 0, 0], *plane.at(top), "$4m$")
    fig.legend("fill", "cover", "the plane $z = 0$ at one moment, which $t$, $\\rho$ and $\\phi$ cover")
    fig.legend("line", "r", "$\\rho$ constant, at $1$, $2$, $3$ and $4$ times $m$")
    fig.legend("line", "meridian", "$\\phi$ constant, every $15°$")
    views.append(view("two_holes", "Two holes", "$m$", [two], fig.done(),
                      settings="$m = 1$, the unit of every length, and two holes of mass parameter $m$ on the axis at "
                               "$z = \\pm 2m$, so that $U = 1 + 2m/\\sqrt{\\rho^2 + 4m^2}$ on the plane $z = 0$."))
    return views


def stockum_dust(ck, src):
    """The plane z = 0 at one moment, R = 1: g_rr = e^(-r^2), g_phiphi = r^2 (1 - r^2). The
    circles grow out to r = R/sqrt(2) and shrink after, so the surface curls back toward the
    axis; g_rr - (drho/dr)^2 falls to zero at r = 0.83 R, found here, beyond which the circles
    shrink faster than the distance out to them; at r = R they are null and beyond it timelike,
    the closed timelike curves, both checked."""
    sl = Slice(src, "stockum_dust", "cylindrical", "r", "\\phi", {"t": 0, "z": 0}, {"R": 1})
    stop = float(sp.nsolve(sl.defect, sl.x, 0.83))
    ck.stops("van Stockum, between the last surface and r = R", sl, np.linspace(stop, 1, 202)[1:-1])
    beyond = sl.gpp_at(np.linspace(1, 3, 201)[1:])
    ck.add("van Stockum: beyond r = R the circles are timelike, g_phiphi < 0", float(max(0.0, np.max(beyond))), 0.0)
    size = 2 * 0.72
    widest = 1 / math.sqrt(2)
    dust = Piece("dust", "star", sl, 0.0, stop, 0.0, 1,
                 (("axis", "the axis $r = 0$"),
                  ("stops", "the circles shrink faster than the distance out to them, and nothing in flat space carries the slice on")),
                 [(0.25, "r", None), (0.5, "r", None), (widest, "r", None), (stop, "chartedge", None)], size)
    ck.isometry("van Stockum, the dust", dust)
    ck.add("van Stockum: the widest circle is at r = R/sqrt(2)", abs(float(np.max(dust.rho)) - 0.5), 1e-9)
    surface = Surface([dust])
    fig = figure_of([surface], {"star": "star"}, size, Camera(-90, 32))
    ring_label(fig, [0, 0, 0], *dust.at(widest), "$R/\\sqrt{2}$")
    ring_label(fig, [0, 0, 0], *dust.at(stop), f"${stop:.2f}\\,R$", side=-1)
    fig.legend("fill", "star", "the rotating dust, which $t$, $r$ and $\\phi$ cover")
    fig.legend("line", "r", "$r$ constant, at $R/4$, $R/2$ and $R/\\sqrt{2}$, the widest circle")
    fig.legend("line", "chartedge", f"$r = {stop:.2f}\\,R$, where the drawing stops")
    fig.legend("line", "meridian", "$\\phi$ constant, every $15°$")
    return [view("dust", "The rotating dust", "$R$", [surface], fig.done(),
                 settings="$R = 1$, the unit of every length.",
                 stops=[f"From $r = {stop:.2f}\\,R$ the circles shrink faster than the distance out to them, "
                        "$g_{rr} < (\\partial_r\\sqrt{g_{\\phi\\phi}})^2$, and no surface of revolution in flat space "
                        "carries the slice on.",
                        "At $r = R$ the circles about the axis are null, and beyond it they are closed timelike "
                        "curves, so a surface of constant $t$ is not a moment of space there."])]


def taub_nut(ck, src):
    """The equator of a slice of constant t at m = 1 and l = 1/2, as the spacetime diagram draws
    Taub-NUT: g_tphi carries cos(theta) and vanishes there, so the slice has g_rr = (r^2 + l^2)/
    (r^2 - 2mr - l^2) and g_phiphi = r^2 + l^2, a surface of revolution out of the horizon r+ = m
    + sqrt(m^2 + l^2), where it stands vertical with circumference radius sqrt(r+^2 + l^2) =
    sqrt(2(m r+ + l^2)), more than Schwarzschild's 2m, which is checked with r+ from its formula.
    Across the horizon lies Taub's cosmology, where r is a time."""
    sl = Slice(src, "taub_nut", "spherical", "r", "\\phi", {"t": 0, **EQUATOR}, {"m": 1, "l": "1/2"})
    rp = sl.horizons()[0]
    top = 8.0
    size = 2 * float(sl.rho_at(top))
    near = Piece("exterior", "sheet", sl, rp, top, 0.0, 1,
                 (("throat", "the horizon $r = r_+$, where the NUT region meets Taub's cosmology"),
                  ("edge", "the surface runs on to $r \\to \\infty$")),
                 [(rp, "horizon", "$r = r_+$")] + [(r, "r", None) for r in (3.0, 4.0, 5.0, 6.0, 7.0, top)], size)
    ck.isometry("Taub-NUT, the equator", near)
    formula = 1 + math.sqrt(1.25)
    ck.add("Taub-NUT: the horizon is r+ = m + sqrt(m^2 + l^2)", abs(rp - formula), 1e-12)
    ck.add("Taub-NUT: the horizon's circumference radius is sqrt(2(m r+ + l^2))",
           abs(near.at(rp)[0] - math.sqrt(2 * (formula + 0.25))), 1e-6)
    ck.stops("Taub-NUT, across the horizon", sl, np.linspace(1 - math.sqrt(1.25), rp, 202)[1:-1])
    surface = Surface([near])
    fig = figure_of([surface], {"sheet": "cover"}, size)
    ring_label(fig, [0, 0, 0], *near.at(rp), "$r = r_+$", dx=14)
    ring_label(fig, [0, 0, 0], *near.at(top), "$8\\,m$")
    fig.legend("fill", "cover", "the NUT region $r > r_+$, which $t$ and $r$ cover")
    fig.legend("line", "r", "$r$ constant, at $3$, $4$, $5$, $6$, $7$ and $8\\,m$")
    fig.legend("line", "horizon", "the horizon $r = r_+$, where the surface stands vertical")
    fig.legend("line", "meridian", "$\\phi$ constant, every $15°$")
    return [view("equator", "The equator", "$m$", [surface], fig.done(),
                 settings=f"$m = GM/c^2 = 1$, the unit of every length, and $l = m/2$, so that $r_+ = {rp:.3f}\\,m$.",
                 stops=["Across the horizon lies Taub's cosmology, where $r$ is a time and a slice of constant $t$ is "
                        "not a moment of space, so the surface ends at $r_+$."])]


def godel(ck, src):
    """The plane z = 0 about one world line of the dust, at one moment of the cylindrical chart's
    t, omega = 1: g_rr = 2 and g_phiphi = 2 (1 - sinh^2 r) sinh^2 r, so g_rr - (drho/dr)^2 =
    2 (1 - (1 + s)(1 - 2s)^2/(1 - s)) with s = sinh^2 r, which vanishes where 4 s^3 = 2 s: the
    construction stops at sinh^2 r = 1/sqrt(2), exactly, which is checked, before the circles
    turn null at sinh r = 1, r_c, and timelike beyond, both checked."""
    sl = Slice(src, "godel", "cylindrical", "r", "\\phi", {"t": 0, "z": 0}, {"omega": 1})
    stop = math.asinh(2 ** -0.25)
    rc = math.asinh(1)
    ck.add("Godel: the surface stops at sinh^2 r = 1/sqrt(2), where g_rr = (drho/dr)^2",
           abs(float(sl.defect_at(stop))), 1e-12)
    ck.stops("Godel, between the last surface and r_c", sl, np.linspace(stop, rc, 202)[1:-1])
    beyond = sl.gpp_at(np.linspace(rc, 3, 201)[1:])
    ck.add("Godel: beyond r_c the circles are timelike, g_phiphi < 0", float(max(0.0, np.max(beyond))), 0.0)
    size = 2 * 0.72
    widest = math.asinh(math.sqrt(0.5))
    dust = Piece("dust", "star", sl, 0.0, stop, 0.0, 1,
                 (("axis", "the world line $r = 0$ of the dust"),
                  ("stops", "the circles shrink faster than the distance out to them, and nothing in flat space carries the slice on")),
                 [(rc / 4, "r", None), (rc / 2, "r", None), (widest, "r", None), (stop, "chartedge", None)], size)
    ck.isometry("Godel, about one world line", dust)
    surface = Surface([dust])
    fig = figure_of([surface], {"star": "star"}, size, Camera(-90, 32))
    ring_label(fig, [0, 0, 0], *dust.at(widest), "$\\sinh^2 r = 1/2$")
    fig.legend("fill", "star", "the dust about one of its world lines, which $t$, $r$ and $\\phi$ cover")
    fig.legend("line", "r", "$r$ constant, at $r_c/4$ and $r_c/2$, and the widest circle, $\\sinh^2 r = 1/2$")
    fig.legend("line", "chartedge", "$\\sinh^2 r = 1/\\sqrt{2}$, where the drawing stops")
    fig.legend("line", "meridian", "$\\phi$ constant, every $15°$")
    return [view("dust", "About one world line", "$1/\\omega$", [surface], fig.done(),
                 settings="$\\omega = 1$, so that $c/\\omega$ is the unit of every length and the proper distance out "
                          "from the axis is $\\sqrt{2}\\,r$.",
                 stops=["From $\\sinh^2 r = 1/\\sqrt{2}$ the circles shrink faster than the distance out to them and no "
                        "surface of revolution in flat space carries the slice on.",
                        "At $\\sinh r = 1$, $r_c$, the circles about the axis are null, and beyond it they are closed "
                        "timelike curves, so a surface of constant $t$ is not a moment of space there. Gödel's universe "
                        "has no moment of time that is space everywhere, and every world line of its dust is an axis "
                        "like this one."])]


def ellis_bronnikov(ck, src):
    """In its own chart r is the proper distance from the throat, g_rr = 1 and g_phiphi = r^2 +
    l^2, so dz/dr = l/sqrt(r^2 + l^2) and z = l arcsinh(r/l): the catenoid rho = l cosh(z/l),
    the surface Morris-Thorne draws, here as one piece through the throat, since the chart
    runs from one side to the other."""
    sl = Slice(src, "ellis_bronnikov", "spherical", "r", "\\phi", {"t": 0, **EQUATOR}, {"ell": 1})
    top = 5.0
    size = 2 * math.sqrt(top ** 2 + 1)
    steps = (1.0, 2.0, 3.0, 4.0)
    whole = Piece("whole", "sheet", sl, -top, top, -math.asinh(top), 1,
                  (("edge", "the side $r < 0$ runs on, flattening, to $r \\to -\\infty$"),
                   ("edge", "the side $r > 0$ runs on, flattening, to $r \\to \\infty$")),
                  [(-top, "r", None)] + [(-r, "r", None) for r in reversed(steps)] + [(0.0, "throat", "$r = 0$")]
                  + [(r, "r", None) for r in steps] + [(top, "r", None)], size)
    ck.isometry("Ellis-Bronnikov, through the throat", whole)
    ck.form("Ellis-Bronnikov, the catenoid z = l arcsinh(r/l)", whole, np.arcsinh, size)
    ck.radius("Ellis-Bronnikov, the catenoid rho = sqrt(r^2 + l^2)", whole, lambda r: np.sqrt(r * r + 1), size)
    surface = Surface([whole])

    fig = figure_of([surface], {"sheet": "cover"}, size)
    ring_label(fig, [0, 0, 0], *whole.at(0.0), "$r = 0$", dx=14)
    ring_label(fig, [0, 0, 0], *whole.at(2.0), "$2\\ell$")
    ring_label(fig, [0, 0, 0], *whole.at(-2.0), "$-2\\ell$")
    fig.legend("fill", "cover", "the whole slice, which $t$ and $r$ cover from one side to the other")
    fig.legend("line", "r", "$r$ constant, at $\\pm\\ell$, $\\pm 2\\ell$, $\\pm 3\\ell$, $\\pm 4\\ell$ and $\\pm 5\\ell$")
    fig.legend("line", "throat", "the throat $r = 0$, the smallest circle, of radius $\\ell$")
    fig.legend("line", "meridian", "$\\phi$ constant, every $15°$")
    return [view("wormhole", "The wormhole", "$\\ell$", [surface], fig.done(),
                 settings="$\\ell = 1$, the unit of every length.")]


def cosmic_string(ck, src):
    """The conical exterior: g_rr = 1 and g_phiphi = (1 - 4G mu/c^2)^2 r^2, so rho =
    (1 - 4G mu/c^2) r and dz/dr = sqrt(1 - (1 - 4G mu/c^2)^2): a cone of half angle
    arcsin(1 - 4G mu/c^2), its apex the string. Cut along a line from the apex and laid
    flat it is a plane missing the wedge 8 pi G mu/c^2, which the figure draws beside it.
    Gott's core, l^2(dchi^2 + sin^2 chi dphi^2), is a cap of a sphere of radius l, and meets
    the cone at chi_0 with one tangent, cos chi_0 = 1 - 4G mu/c^2, where the exterior's
    r = l tan chi_0. Drawn at 4G mu/c^2 = 0.1 with l = 1."""
    fold = 0.9                          # 1 - 4G mu/c^2
    chi0 = math.acos(fold)
    join = math.tan(chi0)
    top = 3.0
    size = 2 * top
    cone = Slice(src, "cosmic_string", "conical", "r", "\\phi", {"t": 0, "z": 0}, {"mu": "1/40", "G": 1})
    cap = Slice(src, "cosmic_string", "interior_cap", "\\chi", "\\phi", {"t": 0, "z": 0}, {"ell": 1})
    ideal = Piece("apex", "reference", cone, 0.0, join, 0.0, 1,
                  (("apex", "the string itself, a conical singularity at $r = 0$"), ("join", None)),
                  [], size, reference=True)
    zj = ideal.at(join)[1]
    ext = Piece("exterior", "sheet", cone, join, top, zj, 1,
                (("join", "the edge of Gott's core, $r = \\ell\\tan\\chi_0$"), ("edge", "the cone runs on to $r \\to \\infty$")),
                [(join, "surface", None)] + [(r, "r", None) for r in (1, 2)] + [(top, "r", None)], size)
    core = Piece("core", "star", cap, 0.0, chi0, 0.0, 1,
                 (("axis", "the axis $\\chi = 0$, where the cap is smooth"), ("join", "the edge of the core, $\\chi = \\chi_0$")),
                 [], size)
    core.z = core.z + (zj - core.z[-1])
    ck.isometry("cosmic string, the cone", ext)
    ck.isometry("cosmic string, the cone to its apex", ideal)
    ck.isometry("cosmic string, Gott's core", core)
    ck.join("cosmic string, Gott's core meets the cone at chi_0", core, chi0, ext, join)
    slope = math.sqrt(1 - fold ** 2)
    ck.form("cosmic string, the cone z = sqrt(1 - (1 - 4G mu/c^2)^2) r", ext, lambda r: slope * r, size)
    ck.radius("cosmic string, the cone rho = (1 - 4G mu/c^2) r", ext, lambda r: fold * r, size)
    ck.form("cosmic string, Gott's cap z = l (1 - cos chi)", core, lambda c: core.z[0] + 1 - np.cos(c), size)
    ck.radius("cosmic string, Gott's cap rho = l sin chi", core, np.sin, size)
    surface = Surface([core, ext, ideal])

    # Seen from 20 degrees up, below the 26 degree slope of the cone's wall, so that its sides
    # stand out as a cone while its inside, with the core at the bottom, still shows.
    camera = Camera(-90, 20)
    fig = figure_of([surface], {"star": "star", "sheet": "cover"}, size, camera)
    # The cut, the meridian at phi = 0, which the development opens along.
    fig.mark("cut", 0, ext, 0.0)
    ring_label(fig, [0, 0, 0], *ext.at(1.0), "$r = \\ell$", side=-1, clear=True)
    ring_label(fig, [0, 0, 0], *ext.at(top), "$3\\ell$", side=-1)

    # The development: the cone laid flat beside it, a disc of radius r missing the wedge
    # delta = 2 pi (1 - fold), in the plane of the page. A circle of radius r on the cone
    # has length 2 pi fold r, so it is an arc of angle 2 pi fold; meridian phi lies at the
    # angle fold phi from the cut's first edge. It is drawn at half the cone's scale,
    # which its label says, so that the cone, the surface itself, is the larger drawing.
    half = 0.5
    E = np.vstack(fig.extent)
    centre = np.array([E[:, 0].max() + 0.5 + half * top, 0.5 * (E[:, 1].min() + E[:, 1].max())])
    gap = 2 * math.pi * (1 - fold)
    first = gap / 2                         # the wedge opens to the right, as the cut lies on the cone

    def flat(r, psi):
        return centre + half * np.column_stack([r * np.cos(psi), r * np.sin(psi)])
    arc = np.linspace(first, first + 2 * math.pi * fold, 400)
    fig.plain_fill("cover", np.vstack([flat(top, arc), flat(join, arc[::-1])]))
    fig.plain_fill("wedge", np.vstack([centre, flat(top, np.linspace(first - gap, first, 60))]))
    for k in range(1, 24):
        psi = first + fold * 2 * math.pi * k / 24
        fig.flat("meridian", flat(np.array([join, top]), np.array([psi, psi])))
        if k % 2 == 0:
            fig.flat("reference", flat(np.array([0, join]), np.array([psi, psi])))
    for psi in (first, first + 2 * math.pi * fold):
        fig.flat("cut", flat(np.array([0, top]), np.array([psi, psi])))
    for r in (1.0, 2.0):
        fig.flat("r", flat(r, arc))
    fig.flat("outline", flat(top, arc))
    fig.flat("surface", flat(join, arc))
    fig.label(centre + [half * top * 0.62, 0], "$\\delta$", "c", "lab")
    fig.label(centre + [0, half * top], "laid flat, at half the scale", "b", "small", dy=-6)
    fig.legend("fill", "star", "Gott's core, $\\chi \\le \\chi_0$, a cap of a sphere of radius $\\ell$")
    fig.legend("fill", "cover", "the conical exterior, a cone of half angle $\\arcsin(1 - 4G\\mu/c^2)$")
    fig.legend("fill", "wedge", "the wedge $\\delta = 8\\pi G\\mu/c^2$ the cone lacks, laid flat")
    fig.legend("line", "cut", "the line the cone is cut along to lay it flat")
    fig.legend("line", "r", "$r$ constant, at $\\ell$, $2\\ell$ and $3\\ell$")
    fig.legend("line", "surface", "the edge of the core, $\\chi = \\chi_0$ and $r = \\ell\\tan\\chi_0$")
    fig.legend("line", "reference", "the cone inside the core, down to its apex, where an ideal string lies")
    fig.legend("line", "meridian", "$\\phi$ constant, every $15°$")
    settings = ("$4G\\mu/c^2 = 0.1$, so that $\\cos\\chi_0 = 0.9$ and $\\delta = 36°$, and $\\ell = 1$, the unit "
                "of every length.")
    return [unrolling(ck, cone, fold, top, size, settings),
            view("cone", "Gott's core", "$\\ell$", [surface], fig.done(), settings=settings)]


def unrolling(ck, cone, fold, top, size, settings):
    """The cone of an ideal string, cut along the meridian phi = 0 and unrolled onto the plane: a
    movie of cones of half angle alpha, sin(alpha) = 2 pi fold/(2 pi - d), each carrying the whole
    cone, every circle of it running round the axis through 2 pi - d, so that the cut's two
    edges stand d apart, from d = 0, the cone itself, to d = delta = 2 pi (1 - fold) = 8 pi G mu/c^2,
    where alpha = 90 degrees and the cone lies flat, a disc missing the wedge delta. The point at r
    and chart angle phi stands at the distance r from the apex and the angle d/2 + phi (2 pi -
    d)/2 pi round the axis, so each frame is the cone moved without stretching: along each
    straight line from the apex every distance is r, and each circle is as long as the published
    g_phiphi makes it, which is checked on every frame, as is the wedge of the last."""
    phi = np.linspace(0, 2 * math.pi, 361)
    delta = 2 * math.pi * (1 - fold)
    ck.add("cosmic string: the deficit angle 2 pi (1 - 4G mu/c^2) is 8 pi G mu/c^2 at mu = 1/40, G = 1",
           abs(delta - 8 * math.pi / 40), 1e-15)
    whole = Piece("cone", "sheet", cone, 0.0, top, 0.0, 1,
                  (("apex", "the string itself, a conical singularity at $r = 0$"), ("edge", "the cone runs on to $r \\to \\infty$")),
                  [(1.0, "r", None), (2.0, "r", None)], size)
    ck.isometry("cosmic string, the ideal string's cone", whole)
    rs = np.array([0.0, 1.0, 2.0, top])
    frames, openings = [], [math.radians(k) for k in range(37)]
    openings[-1] = delta
    for d in openings:
        sin_a = 1.0 if d == delta else 2 * math.pi * fold / (2 * math.pi - d)
        cos_a = math.sqrt(1 - sin_a * sin_a)
        v = d / 2 + phi * (2 * math.pi - d) / (2 * math.pi)
        v[0], v[-1] = d / 2, 2 * math.pi - d / 2

        def rows(u, s=sin_a, c=cos_a):
            u = np.asarray(u, dtype=float)
            return s * u, s * u, c * u
        piece = EllipsesPiece("cone", "sheet", rs, v, rows, "cosmic_string", "conical",
                              "the cone runs on to $r \\to \\infty$", size,
                              [("cut", "v", [v[0], v[-1]]), ("meridian", "v", list(v[15:-1:15])),
                               ("r", "u", [1.0, 2.0]), ("outline", "u", [top])], open_=True)
        where = f"cosmic string, unrolled by {math.degrees(d):.2f} degrees"
        N = piece.nodes()
        along = np.linalg.norm(np.diff(N, axis=0), axis=-1)
        ck.add(f"{where}: along each line from the apex, every distance is r", float(np.max(np.abs(along - np.diff(piece.u)[:, None]))),
               1e-6)
        arc = piece.A[:, None] * np.diff(piece.v)[None, :]
        want = np.array([float(cone.rho_at(r)) for r in piece.u])[:, None] * np.diff(phi)[None, :]
        ck.add(f"{where}: round each circle, as long as sqrt(g_phiphi) makes it", float(np.max(np.abs(arc - want))), 1e-6)
        ck.add(f"{where}: the cut's edges stand d apart", abs(2 * math.pi - (piece.v[-1] - piece.v[0]) - d), 1e-12)
        label = f"$\\Delta\\phi = {math.degrees(d):.0f}°" + (" = \\delta$" if d == delta else "$")
        frames.append(Surface([piece], label=label, time=d))
    first, last = frames[0].pieces[0], frames[-1].pieces[0]
    ck.add("cosmic string, unrolled: the first frame is the cone, rho = (1 - 4G mu/c^2) r and z = sqrt(1 - (1 - 4G mu/c^2)^2) r",
           float(max(np.max(np.abs(first.A - fold * first.u)), np.max(np.abs(first.Zr - math.sqrt(1 - fold ** 2) * first.u)))),
           1e-6)
    ck.add("cosmic string, unrolled: the last frame lies flat", float(np.max(np.abs(last.Zr))), 0.0)
    ck.add("cosmic string, unrolled: the wedge the last frame lacks is the deficit angle 8 pi G mu/c^2",
           abs(2 * math.pi - (last.v[-1] - last.v[0]) - 8 * math.pi / 40), 1e-12)
    fig = movie_figure(frames, {"sheet": "cover"}, size, Camera(-90, 30))
    fig.legend("fill", "cover", "the cone of an ideal string, cut along $\\phi = 0$ and unrolling onto the plane")
    fig.legend("line", "cut", "the two edges of the cut, which part by $\\Delta\\phi$, to the deficit angle $\\delta = "
                              "8\\pi G\\mu/c^2$ where the cone lies flat")
    fig.legend("line", "r", "$r$ constant, at $\\ell$ and $2\\ell$")
    fig.legend("line", "meridian", "$\\phi$ constant, every $15°$")
    return view("unroll", "The cone unrolled", "$\\ell$", [Surface([whole])], fig.done(),
                movie=movie(frames, "$\\Delta\\phi$", [math.degrees(d) for d in openings], loop="pingpong"),
                settings=settings)


def frw(ck, src):
    """The closed universe of dust, a = 1 - cos eta and ct = eta - sin eta in units of
    1/sqrt(k), the scale factor the conformal diagram declares, checked here to make the
    published G^r_r vanish. At one moment the equator of the three sphere, g_rr = a^2/(1 -
    r^2) and g_phiphi = a^2 r^2, is a sphere of radius a: z = a(1 - sqrt(1 - r^2)) on the
    hemisphere r = sin chi covers, and the same again on the far hemisphere. Drawn at five
    moments from the bang to the crunch. For k = 0 the slice is a plane, and for k = -1 it
    has g_rr - (drho/dr)^2 = -a^2 r^2/(1 + r^2) < 0 at every r > 0, so no part of it about
    its centre is a surface of revolution in flat space; both are checked at the scale
    factors of dust, eta^2 and cosh eta - 1."""
    src.note("frw", "conformal_spherical", ["einstein_tensor"])
    _, entry, reader = nr.load("frw", "conformal_spherical")
    ul = entry["einstein_tensor"]["variants"]["ul"]["nonzero"]
    G = reader(next(c["value"] for c in ul if c["indices"] == ["r", "r"]))
    eta = reader.symbol["\\eta"]
    G = G.subs(reader.parameters["a"], 1 - sp.cos(eta)).doit().subs({reader.parameters["k"]: 1, reader.c: 1})
    ck.exact("FRW: a = 1 - cos eta makes the published G^r_r vanish, k = 1", sp.simplify(G) == 0)

    moments = [math.pi * f for f in (1 / 3, 2 / 3, 1, 4 / 3, 5 / 3)]
    size = 4.0
    from scipy.optimize import brentq

    def moment(e):
        a = 1 - math.cos(e)
        sl = Slice(src, "frw", "comoving_spherical", "r", "\\phi", {"t": 0, **EQUATOR}, {"k": 1}, {"a": repr(a)})
        near = Piece("near", "sheet", sl, 0.0, 1.0, -a, 1,
                     (("axis", "the centre $\\chi = 0$"), ("join", "the equator $r = 1$, where the far hemisphere begins")),
                     [(0.5, "r", None), (math.sqrt(3) / 2, "r", None), (1.0, "chartedge", None)], size)
        far = Piece("far", "sheet2", sl, 0.0, 1.0, a, -1,
                    (("axis", "the antipode $\\chi = \\pi$"), ("join", "the equator")),
                    [(0.5, "r2", None), (math.sqrt(3) / 2, "r2", None)], size)
        where = f"FRW closed, eta = {e / math.pi:.4f} pi"
        ck.isometry(f"{where}, near hemisphere", near)
        ck.isometry(f"{where}, far hemisphere", far)
        ck.join(f"{where}, the hemispheres meet at the equator", near, 1.0, far, 1.0)
        ck.form(f"{where}, the sphere z = a(1 - sqrt(1 - r^2)) - a", near,
                lambda r, a=a: a * (1 - np.sqrt(np.maximum(1 - r * r, 0))) - a, size)
        ck.radius(f"{where}, the sphere rho = a r", near, lambda r, a=a: a * r, size)
        t = e - math.sin(e)
        return Surface([near, far], label=f"$ct = {t:.2f}$", time=t)

    # The movie runs through the moments at a steady cosmic time, a frame about every 0.1 of ct.
    times, keys = movie_values([e - math.sin(e) for e in moments], 0.1)
    etas = [moments[keys.index(i)] if i in keys else brentq(lambda e, t=t: e - math.sin(e) - t, 1e-9, 2 * math.pi, xtol=1e-15)
            for i, t in enumerate(times)]
    frames = [moment(e) for e in etas]
    surfaces = [frames[i] for i in keys]
    for k, name, a_of in ((0, "flat", lambda e: e * e), (-1, "open", lambda e: math.cosh(e) - 1)):
        for e in moments:
            sl = Slice(src, "frw", "comoving_spherical", "r", "\\phi", {"t": 0, **EQUATOR}, {"k": k},
                       {"a": repr(a_of(e))})
            xs = np.linspace(1e-3, 20, 400)
            if k == 0:
                ck.plane(f"FRW flat, eta = {e / math.pi:.3f} pi", sl, xs)
            else:
                ck.stops(f"FRW open, eta = {e / math.pi:.3f} pi", sl, xs)

    fig = movie_figure(frames, {"sheet": "cover"}, size, meridians=12)
    fig.legend("fill", "cover", "the hemisphere $\\chi < \\pi/2$ that $r = \\sin\\chi$ covers")
    fig.legend("line", "r", "$\\chi$ constant, at $\\pi/6$ and $\\pi/3$")
    fig.legend("line", "r2", "$\\chi$ constant on the far hemisphere, at $2\\pi/3$ and $5\\pi/6$")
    fig.legend("line", "chartedge", "the equator $\\chi = \\pi/2$, where $r = 1$")
    fig.legend("line", "meridian", "$\\phi$ constant, every $30°$")
    return [view("closed", "Closed dust, $k = +1$", "$1/\\sqrt{k}$", surfaces, fig.done(),
                 movie=movie(frames, "$ct$", [f.time for f in frames], turns=False),
                 settings="$k = 1$, with $1/\\sqrt{k}$ the unit of every length.",
                 input="Dust: $a = 1 - \\cos\\eta$ with $ct = \\eta - \\sin\\eta$, solved from this "
                       "spacetime's own $G^r{}_r = 0$.",
                 stops=["For $k = 0$ every slice of constant $t$ is flat, and its equator is a plane.",
                        "For $k = -1$ the circles about any point grow faster than the distance out to "
                        "them, $g_{rr} < (\\partial_r\\sqrt{g_{\\phi\\phi}})^2$ at every $r > 0$, so no "
                        "surface of revolution about a point carries the slice. A piece of it lies in "
                        "flat space on Eugenio Beltrami's pseudosphere, and David Hilbert proved in 1901 "
                        "that no surface in flat space carries the whole of it."])]



def misner(ck, src):
    """The contracting region T < 0 of Misner space in the Milne chart, at Li and Gott's
    psi_0 = 4 pi, where chi runs once round 2 pi. The slice z = 0 of one moment t has the metric
    dy^2 + c^2t^2 dchi^2: a flat cylinder of radius c|t| about the axis y, z = y. Drawn at four
    moments as the closed direction shrinks toward the chronology horizon t = 0, each over
    -2 <= y <= 2 in any length l, since flat space names none."""
    moments = (-2.0, -1.5, -1.0, -0.5)
    surfaces, size = [], 4.0
    for ct in moments:
        sl = Slice(src, "misner", "milne", "y", "\\chi", {"t": repr(ct), "z": 0}, {"psi_0": "4*pi"})
        tube = Piece("tube", "sheet", sl, -2.0, 2.0, -2.0, 1,
                     (("edge", "the cylinder runs on for ever toward $y \\to -\\infty$"),
                      ("edge", "the cylinder runs on for ever toward $y \\to \\infty$")),
                     [(-1.0, "r", None), (0.0, "r", None), (1.0, "r", None)], size)
        where = f"Misner, the cylinder at ct = {ct}"
        ck.isometry(where, tube)
        ck.radius(f"{where}, rho = c|t|", tube, lambda y, ct=ct: np.full_like(y, abs(ct)), size)
        ck.form(f"{where}, z = y", tube, lambda y: y, size)
        surfaces.append(Surface([tube], label=f"$ct = {ct:g}$", time=ct))
    # The cylinders stand far enough apart that their names under them, at the size the page
    # sets them at, stand clear of each other.
    offsets, x, gap = [], 0.0, 2.0
    for s in surfaces:
        rho = float(s.pieces[0].rho[0])
        offsets.append((x + rho, 0.0, 0.0))
        x += 2 * rho + gap
    fig = figure_of(surfaces, {"sheet": "cover"}, 8.0, offsets=offsets, meridians=12)
    base = min(fig.screen(np.asarray(off) + [0, 0, -2.0])[1] for off in offsets)
    for s, off in zip(surfaces, offsets):
        fig.label(np.array([fig.screen(np.asarray(off))[0], base]), s.label, "t", "small", dy=8)
    fig.legend("fill", "cover", "the slice $z = 0$ of a moment of constant $t$, which $\\chi$ and $y$ cover")
    fig.legend("line", "r", "$y$ constant, at $-\\ell$, $0$ and $\\ell$, each a circle of circumference $2\\pi c|t|$")
    fig.legend("line", "meridian", "$\\chi$ constant, every $30°$")
    return [view("cylinders", "The contracting region", "$\\ell$", surfaces, fig.done(),
                 settings="$\\psi_0 = 4\\pi$, and $\\ell$, any length, the unit of every length, since flat space "
                          "has none of its own.")]




# ---------------------------------------------------------------- the eleven drawn last

FLAT_CAMERA = Camera(-90, 55)   # a flat plane seen from well above it, so that what is marked on it shows


def published_christoffel(src, metric_id, system_id):
    """The published Gamma^a_bc of a chart, as {(a, b, c): sympy}, with c = 1."""
    src.note(metric_id, system_id, ["christoffel"])
    _, entry, R = nr.load(metric_id, system_id)
    ull = entry["christoffel"]["variants"]["ull"]["nonzero"]
    return {tuple(c["indices"]): R(c["value"]).subs(R.c, 1) for c in ull}, R


def disc(sl, pid, top, centre, rim, marks=(), size=1.0, cls="sheet"):
    """A flat piece from the centre of a plane out to the proper radius `top`, level at z = 0."""
    return Piece(pid, cls, sl, 0.0, top / math.sqrt(float(sl.gxx)), 0.0, 1, (("axis", centre), ("edge", rim)),
                 marks, size)


def particles(ck, where, sl, piece, cx, cy, every=30):
    """A ring of free particles at the chart points (cx, cy) of a flat plane, drawn as the curve
    through all of them and as a point at every `every` of them, each checked to lie on the
    plane, and the drawing's lengths checked against the plane's metric."""
    P = sl.draw(cx, cy)
    ck.on_piece(f"{where}, the ring", piece, P)
    chord = np.linalg.norm(np.diff(P, axis=0), axis=1)
    metric = np.sqrt(float(sl.gxx) * np.diff(cx) ** 2 + float(sl.scale[1] ** 2) * np.diff(cy) ** 2)
    ck.add(f"{where}, the ring's steps against the plane's metric", float(np.max(np.abs(chord - metric))), 1e-12)
    return Curve(piece, "particles", P, closed=True), [(piece, "particles", Q) for Q in P[:-1:every]]


def minkowski(ck, src):
    """The equator of a slice of constant t in the spherical chart: g_rr = 1 and g_phiphi = r^2,
    so rho = r and the surface is the plane, on which a circle of radius r has circumference
    2 pi r, the surface every other embedding diagram is measured against. Drawn out to r = 4, as
    its spacetime diagram runs, in any unit, since flat space has none of its own. Every slice of
    constant t of the Cartesian chart is the same flat space, which is checked."""
    flat_slices(ck, src, "minkowski", "cartesian")
    sl = Slice(src, "minkowski", "spherical", "r", "\\phi", {"t": 0, **EQUATOR})
    ck.plane("Minkowski, the equator of the spherical chart", sl, np.linspace(1e-3, 20, 400))
    top = 4.0
    size = 2 * top
    plane = Piece("plane", "sheet", sl, 0.0, top, 0.0, 1,
                  (("axis", "the centre $r = 0$"), ("edge", "the plane runs on to $r \\to \\infty$")),
                  [(r, "r", None) for r in (1.0, 2.0, 3.0, top)], size)
    ck.isometry("Minkowski, the plane", plane)
    ck.radius("Minkowski, the plane rho = r", plane, lambda r: r, size)
    ck.form("Minkowski, the plane z = 0", plane, np.zeros_like, size)
    surface = Surface([plane])
    fig = figure_of([surface], {"sheet": "cover"}, size, FLAT_CAMERA)
    ring_label(fig, [0, 0, 0], *plane.at(1.0), "$r = \\ell$")
    ring_label(fig, [0, 0, 0], *plane.at(top), "$4\\,\\ell$")
    fig.legend("fill", "cover", "the plane, which $t$ and $r$ cover")
    fig.legend("line", "r", "$r$ constant, at $\\ell$, $2\\ell$, $3\\ell$ and $4\\ell$, each of circumference $2\\pi r$")
    fig.legend("line", "meridian", "$\\phi$ constant, every $15°$")
    return [view("plane", "The plane", "$\\ell$", [surface], fig.done(),
                 settings="$\\ell$, any length, the unit of every length, since flat space has none of its own.")]


def krasnikov(ck, src):
    """How far the tube tips the light cone, 1 - k, drawn as a height over the plane of the tube's
    axis, phi = 0 and pi, at ct = 5, a unit after the ship reached the far end of the tube the
    spacetime diagram declares, along x from 0 to D = 4 behind a ship that left x = 0 at t = 0 at
    the speed of light, so that the whole tube stands. In the plane of t and x the null directions
    of the published metric are (1, 1) and (k, -1), so a light signal sent home over a length L
    takes k L of ct where flat space takes L, and 1 - k, twice the published g_tx, is what it
    gains per unit of length: 0 outside the tube, 1 on the curve k = 0, where the published g_xx
    vanishes and inside which the direction along the tube is a time, and 2 - delta deep inside.
    Over -1 <= x <= 5, as the spacetime diagram runs, and |y| <= 2, as the cross section was
    drawn, on a Cartesian grid with its lines every rho_0/2, the path along the drawing's X from
    x = 0 to D, and the level line k = 0 marked. The slice of constant t is no moment of space
    inside that line, so the height is a plot over the coordinates and no embedding.
    _tools/derivations/krasnikov_height.md is the derivation."""
    T, D = 5.0, 4.0
    src.note("krasnikov", "cylindrical", FIELDS)
    _, entry, R = nr.load("krasnikov", "cylindrical")
    g = nr.published_matrix(R, entry, "metric_components")
    k = R.parameters["k"]
    for name, L in (("outbound", [1, 1, 0, 0]), ("back", [k, -1, 0, 0])):
        L = sp.Matrix(L)
        ck.exact(f"Krasnikov: the {name} direction of the plane of t and x is null for every k",
                 sp.simplify((L.T * g * L)[0]) == 0)
    ck.exact("Krasnikov: 1 - k is twice the published g_tx", sp.simplify(2 * g[0, 1] - (1 - k)) == 0)
    names = {R._plain(n): s for n, s in R.symbol.items()}
    tube = sp.sympify(nr._KRASNIKOV_TUBE, locals=names)
    at = {R.symbol["t"]: T, R.c: 1}
    x, r = R.symbol["x"], R.symbol["r"]
    two_gtx = sp.lambdify((x, r), 2 * g[0, 1].subs(k, tube).subs(at), "numpy")
    gxx = sp.lambdify((x, r), g[1, 1].subs(k, tube).subs(at), "numpy")

    # The drawing's X runs along the axis from the middle of the tube, and Y = +-r across it.
    def height(X, Y):
        return two_gtx(np.asarray(X) + D / 2, np.abs(Y))
    top, half = 3.0, 2.0
    size = 2 * top
    along = [-top + 0.5 * n for n in range(13)]
    across = [-half + 0.5 * n for n in range(9)]
    plane = GridPiece("plane", "sheet", "cartesian", along, across, height, "krasnikov", "cylindrical",
                      "the plane runs on, flat, where $k = 1$", size, {"u": along[1:-1], "v": across[1:-1]})
    X, Y = plane.plane()
    ck.add("Krasnikov: the height at every node is 1 - k, twice the published g_tx",
           float(np.max(np.abs(plane.Z - height(X, Y)))), 1e-7)
    ck.add("Krasnikov: every triangle of the grid lies on the height, in space", plane.sag() / size, GRID_SAG)
    # The path, from x = 0 to D along the axis, through the nodes of the grid.
    j0 = plane.index("v", 0.0)
    i0, i1 = plane.index("u", -D / 2), plane.index("u", D / 2)
    path = finely(plane.nodes()[i0:i1 + 1, j0], size / 360)
    # The level line k = 0, one closed curve round the tube.
    levels = level_curves(plane, 1.0)
    ck.exact("Krasnikov: k = 0 is one closed curve", len(levels) == 1 and levels[0][1])
    wall = finely(np.vstack([levels[0][0], levels[0][0][:1]]), size / 360)[:-1]
    ck.add("Krasnikov: the published g_xx vanishes on the curve k = 0, to the triangles' height",
           float(np.max(np.abs(gxx(wall[:, 0] + D / 2, np.abs(wall[:, 1]))))), 5 * GRID_SAG * size)
    inside = np.linspace(0.02, D - 0.02, 400)
    ck.add("Krasnikov: inside it, on the axis, g_xx < 0, so the direction along the tube is a time",
           float(max(0.0, np.max(gxx(inside, 0.0)))), 0.0)
    ck.items[-1]["ok"] = bool(np.all(gxx(inside, 0.0) < 0))
    for where, P in (("the path", path), ("the curve k = 0", wall)):
        ck.on_piece(f"Krasnikov, {where}", plane, P)
    surface = Surface([plane], curves=[Curve(plane, "wall", wall, closed=True), Curve(plane, "path", path)])
    fig = figure_of([surface], {"sheet": "cover"}, size, HEIGHT_CAMERA)
    fig.legend("fill", "cover", "$1 - k$ as a height over the plane of the tube's axis at $ct = 5\\rho_0$")
    fig.legend("line", "grid", "$x$ constant and $r$ constant, every $\\rho_0/2$")
    fig.legend("line", "path", "the ship's path, from $x = 0$ to $x = D$")
    fig.legend("line", "wall", "$k = 0$, at the height $\\rho_0$: inside it the direction along the tube is a time")
    deep = float(height(np.array(0.0), np.array(0.0)))
    return [view("plane", "Along the tube", "$\\rho_0$", [surface], fig.done(),
                 settings="$ct = 5\\rho_0$, a unit after the ship reached the far end of the tube, with $\\rho_0 = 1$, "
                          "the unit of every length.",
                 input="A tube along $x$ from $0$ to $D = 4$, built by a ship that left $x = 0$ at $t = 0$ at the "
                       "speed of light: $k = 1 - (2 - \\delta)\\,S(\\tfrac{\\rho_0^2 - r^2}{2\\rho_0})\\,S(ct - "
                       "x)\\,S(x)\\,S(D - x)$ with $\\delta = 0.2$, $\\rho_0 = 1$, and $S$ a step of width $0.15$ "
                       "built from $\\tanh$, as in the spacetime diagram.",
                 height=f"$1 - k = 2g_{{tx}}$, a height of $\\rho_0$ for $1 - k = 1$, and ${deep:.4f}\\rho_0$ on the "
                        "axis in the middle of the tube.",
                 stops=["The slice of constant $t$ itself, which is no moment of space inside the curve $k = 0$, "
                        "where the direction along the tube is a time."])]


def expansion(src, metric_id, system_id, functions, t0):
    """The expansion of the observers who ride the slices of constant t, theta = div n with n_mu
    = -N dt, from the published metric with the declared functions, as a numpy function of the
    plane z = 0 at t = t0."""
    _, entry, R = nr.load(metric_id, system_id)
    names = {R._plain(n): s for n, s in R.symbol.items()}
    names.update(R.parameters)
    g = nr.published_matrix(R, entry, "metric_components").subs(R.c, 1)
    for k, v in functions.items():
        g = g.subs(R.parameters[k], sp.sympify(v, locals=names)).doit()
    X = [R.symbol[c] for c in entry["coords"]]
    gi, root = g.inv(), sp.sqrt(-g.det())
    lapse = 1 / sp.sqrt(-gi[0, 0])
    n = [-lapse * gi[a, 0] for a in range(4)]
    theta = sum(sp.diff(root * n[a], X[a]) for a in range(4)) / root
    return sp.lambdify((X[1], X[2]), theta.subs({X[0]: t0, X[3]: 0}), "numpy"), sp.simplify(root)


def alcubierre(ck, src):
    """The expansion of the observers who ride the slices, drawn as a height over the plane z = 0
    of the ship's path at t = 0, where the declared profile centres the bubble on x = 0, v_s = 2
    and Alcubierre's f with R = 1 and sigma = 4, as the spacetime diagram declares. The observers'
    n^mu = (1, v_s f, 0, 0) comes from the published inverse metric and sqrt(-g) = 1, so theta =
    div n = v_s df/dx per unit of ct, c v_s (x - x_s)/r_s df/dr_s per unit of time: negative ahead
    of the ship, where space contracts, and positive behind it, where it expands, as Miguel
    Alcubierre drew it in 1994. The height is theta R^2/4c over the disc of radius 3R about the
    ship, the path along the drawing's X with the ship heading toward +X, on a polar grid of 72
    angles with its circles every R/4 and its lines from the ship every 15 degrees; marked on it,
    the circle v_s f = 1, the zero of the published g_tt, and the level lines where theta is half
    its greatest value ahead of the ship and behind it. The slice itself is flat, which is checked.
    _tools/derivations/alcubierre_expansion.md is the derivation."""
    fns = {"v_s": "2", "f": nr._alcubierre_profile()}
    # The slice is flat, dx^2 + dy^2, which FlatPlane checks as it reads it.
    FlatPlane(src, "alcubierre", "cartesian", "x", "y", {"t": 0, "z": 0}, functions=fns)
    _, entry, R = nr.load("alcubierre", "cartesian")
    names = {R._plain(n): s for n, s in R.symbol.items()}
    names.update(R.parameters)
    g = nr.published_matrix(R, entry, "metric_components")
    gtt = g[0, 0]
    for k, v in fns.items():
        gtt = gtt.subs(R.parameters[k], sp.sympify(v, locals=names)).doit()
    gtt = gtt.subs({R.symbol["t"]: 0, R.symbol["y"]: 0, R.symbol["z"]: 0, R.c: 1})
    wall = float(sp.nsolve(gtt, R.symbol["x"], 1.0))
    ck.add("Alcubierre: the published g_tt vanishes on the circle v_s f = 1", abs(float(gtt.subs(R.symbol["x"], wall))), 1e-12)
    theta, root = expansion(src, "alcubierre", "cartesian", fns, 0)
    ck.exact("Alcubierre: sqrt(-g) = 1, so the expansion is the divergence of n", root == 1)
    rs = sp.Symbol("rs", positive=True)
    f_of = (sp.tanh(4 * (rs + 1)) - sp.tanh(4 * (rs - 1))) / (2 * sp.tanh(4))
    df = sp.lambdify(rs, sp.diff(f_of, rs), "numpy")
    pts = np.random.default_rng(3).uniform(-3, 3, (2000, 2))
    rr = np.hypot(pts[:, 0], pts[:, 1])
    ck.add("Alcubierre: the expansion is v_s (x/r_s) df/dr_s", float(np.max(np.abs(theta(pts[:, 0], pts[:, 1])
                                                                                 - 2 * pts[:, 0] / rr * df(rr)))), 1e-12)
    # The ship heads toward +x, so ahead of it space contracts and behind it expands.
    for (px, py), sign, where in (((1.0, 0.0), -1, "one R ahead of the ship"), ((-1.0, 0.0), 1, "one R behind it"),
                                  ((0.6, 0.8), -1, "ahead, off the path"), ((-0.6, -0.8), 1, "behind, off the path")):
        value = float(theta(px, py))
        ck.exact(f"Alcubierre: theta is {'negative' if sign < 0 else 'positive'} {where}, {value:.6f} c/R",
                 sign * value > 0)
    from scipy.optimize import minimize_scalar
    best = minimize_scalar(lambda r: -float(theta(-r, 0.0)), bounds=(0.5, 1.5), method="bounded",
                           options={"xatol": 1e-12})
    most = float(theta(-best.x, 0.0))

    def height(X, Y):
        with np.errstate(invalid="ignore", divide="ignore"):
            z = theta(X, Y) / 4
        # (x - x_s)/r_s df/dr_s runs to 0 at the centre of the bubble, where df/dr_s does.
        return np.where(np.hypot(X, Y) > 0, z, 0.0)
    top = 3.0
    size = 2 * top
    angles = 2 * math.pi * np.arange(72) / 72
    circles = [0.25 * n for n in range(1, 12)]
    plane = GridPiece("plane", "sheet", "polar", [0.0, *circles, top, wall], angles, height, "alcubierre",
                      "cartesian", "the plane runs on to $r_s \\to \\infty$, where $\\theta$ vanishes", size,
                      {"u": circles, "v": angles[::3]})
    X, Y = plane.plane()
    ck.add("Alcubierre: the height at every node is theta R^2/4c", float(np.max(np.abs(plane.Z - height(X, Y)))), 1e-7)
    ck.add("Alcubierre: every triangle of the grid lies on the height, in space", plane.sag() / size, GRID_SAG)
    # The path through the nodes on the lines from the ship behind it and ahead of it.
    P = plane.nodes()
    back, ahead = plane.index("v", math.pi), plane.index("v", 0.0)
    path = finely(np.vstack([P[::-1, back], P[1:, ahead]]), size / 360)
    ring = plane.index("u", wall)
    circle_ = finely(np.vstack([P[ring], P[ring, :1]]), size / 360)[:-1]
    ck.add("Alcubierre: the circle v_s f = 1 lies at the published zero of g_tt",
           float(np.max(np.abs(np.hypot(P[ring, :, 0], P[ring, :, 1]) - wall))), 1e-12)
    crescents = {}
    for cls, level in (("contract", -most / 8), ("expand", most / 8)):
        found = level_curves(plane, level)
        ck.exact(f"Alcubierre: theta = {4 * level:.4f} c/R is one closed curve", len(found) == 1 and found[0][1])
        C = finely(np.vstack([found[0][0], found[0][0][:1]]), size / 360)[:-1]
        ck.add(f"Alcubierre, the level line {cls}: theta there is half its greatest value, to the triangles' height",
               float(np.max(np.abs(height(C[:, 0], C[:, 1]) - level))), 5 * GRID_SAG * size)
        ck.exact(f"Alcubierre, the level line {cls}: it lies {'ahead of' if level < 0 else 'behind'} the ship",
                 bool(np.all(np.sign(C[:, 0]) == np.sign(-level))))
        crescents[cls] = C
    for where, Q in (("the path", path), ("the circle v_s f = 1", circle_), ("contract", crescents["contract"]),
                     ("expand", crescents["expand"])):
        ck.on_piece(f"Alcubierre, {where}", plane, Q)
    surface = Surface([plane], curves=[Curve(plane, "path", path), Curve(plane, "wall", circle_, closed=True),
                                       Curve(plane, "contract", crescents["contract"], closed=True),
                                       Curve(plane, "expand", crescents["expand"], closed=True)])
    fig = figure_of([surface], {"sheet": "cover"}, size, HEIGHT_CAMERA)
    fig.legend("fill", "cover", "$\\theta$ as a height over the plane $z = 0$ of the ship's path at $t = 0$")
    fig.legend("line", "grid", "circles about the ship every $R/4$, and lines from it every $15°$")
    fig.legend("line", "path", "the ship's path, along $x$, the ship heading toward the trough")
    fig.legend("line", "wall", "$v_sf = 1$, where $g_{tt} = 0$")
    fig.legend("line", "contract", f"$\\theta = -{most / 2:.2f}\\,c/R$, half the fastest contraction, ahead of the ship")
    fig.legend("line", "expand", f"$\\theta = {most / 2:.2f}\\,c/R$, half the fastest expansion, behind it")
    return [view("plane", "The plane of the path", "$R$", [surface], fig.done(),
                 settings="$t = 0$, when the bubble is centred on $x = 0$, with $R = 1$, the unit of every length.",
                 input="$v_s = 2$, and Alcubierre's own profile, $f = [\\tanh\\sigma(r_s + R) - \\tanh\\sigma(r_s - "
                       "R)]/(2\\tanh\\sigma R)$ with $R = 1$ and $\\sigma = 4$, as in the spacetime diagram.",
                 height=f"$\\theta$, a height of $R$ for $\\theta = 4c/R$, and ${most / 4:.4f}R$ where the expansion "
                        f"is fastest, $\\theta = {most:.4f}\\,c/R$.")]


def natario(ck, src):
    """The energy density of the observers who ride the slices, drawn as a height over the plane
    z = 0 of the ship's path at t = 0, with the zero expansion field the spacetime diagram
    declares: Alcubierre's profile with R = 1 and sigma = 4, n = f/2 and v_s = 2. It is read from
    the published G^tt, which for these observers, n_mu = (-1, 0, 0, 0), is 8 pi G/c^4 times the
    energy density epsilon they measure, and checked against -K_ij K^ij/2 from the published
    shift, K_ij = (d_i g_tj + d_j g_ti)/2 on the flat slice, whose trace, the expansion, is checked
    to vanish. So epsilon = -c^4 K_ij K^ij/16 pi G, negative wherever space shears and zero inside
    the bubble and far outside it. The height is G^tt R^2/16 over the disc of radius 3R about the
    ship, the path along the drawing's X with the ship heading toward +X, as Alcubierre's, on a
    polar grid of 72 angles with its circles every R/4 and its lines from the ship every 15
    degrees. Marked on it, lifted onto the triangles, the circle r_s = R, the middle of the wall,
    and three pairs of the lines along which the declared field carries space, found by following
    the field from where each crosses the ship's plane x = 0 at y = +-0.2, 0.4 and 0.6 R until it
    returns there. Inside the bubble they run along the path, and every one closes through the
    wall, since the field is divergence free and zero outside. The slice itself is flat, which is
    checked."""
    fns = nr._natario_field()
    # The slice is flat, dx^2 + dy^2, which FlatPlane checks as it reads it.
    FlatPlane(src, "natario", "cartesian_flow", "x", "y", {"t": 0, "z": 0}, functions=fns)
    _, entry, R = nr.load("natario", "cartesian_flow")
    src.note("natario", "cartesian_flow", ["einstein_tensor"])
    names = {R._plain(n): s for n, s in R.symbol.items()}
    names.update(R.parameters)
    tt, xx, yy, zz = (R.symbol[c] for c in "txyz")

    def declared(e):
        for k, v in fns.items():
            e = e.subs(R.parameters[k], sp.sympify(v, locals=names))
        return e.doit().subs(R.c, 1)
    on_plane = {tt: 0, zz: 0}
    uu = entry["einstein_tensor"]["variants"]["uu"]["nonzero"]
    gtt_up = declared(R(next(c["value"] for c in uu if c["indices"] == ["t", "t"])))
    G = sp.lambdify((xx, yy), gtt_up.subs(on_plane), "numpy")
    # K_ij from the published shift, the metric's g_ti, on the flat slice of unit lapse.
    g = nr.published_matrix(R, entry, "metric_components")
    ck.exact("Natario: the published g_tt is -1 + g_ti g_ti, a unit lapse, and g_ij = delta_ij",
             sp.simplify(g[0, 0] + 1 - sum(g[0, i] ** 2 for i in (1, 2, 3))) == 0
             and all(sp.simplify(g[i, j] - (1 if i == j else 0)) == 0 for i in (1, 2, 3) for j in (1, 2, 3)))
    X3 = (xx, yy, zz)
    shift = [declared(g[0, i]) for i in (1, 2, 3)]
    K = [[(sp.diff(shift[j], X3[i]) + sp.diff(shift[i], X3[j])) / 2 for j in range(3)] for i in range(3)]
    trace = sp.lambdify((xx, yy), sum(K[i][i] for i in range(3)).subs(on_plane), "numpy")
    KK = sp.lambdify((xx, yy), sum(K[i][j] ** 2 for i in range(3) for j in range(3)).subs(on_plane), "numpy")
    pts = np.random.default_rng(3).uniform(-3, 3, (2000, 2))
    ck.add("Natario: the expansion of the riding observers, the trace of K_ij, vanishes",
           float(np.max(np.abs(trace(pts[:, 0], pts[:, 1])))), 1e-12)
    ck.add("Natario: the published G^tt is -K_ij K^ij/2 from the published shift",
           float(np.max(np.abs(G(pts[:, 0], pts[:, 1]) + KK(pts[:, 0], pts[:, 1]) / 2))), 1e-10)
    ck.add("Natario: G^tt runs to 0 at the centre of the bubble, where the flow is uniform",
           abs(float(G(1e-6, 1e-6))), 1e-9)

    def height(X, Y):
        with np.errstate(invalid="ignore", divide="ignore"):
            z = G(X, Y) / 16
        # G^tt is 0/0 as written at the centre, where it runs to 0.
        return np.where(np.hypot(X, Y) > 0, z, 0.0)
    top = 3.0
    size = 2 * top
    angles = 2 * math.pi * np.arange(72) / 72
    circles = [0.25 * n for n in range(1, 12)]
    plane = GridPiece("plane", "sheet", "polar", [0.0, *circles, top], angles, height, "natario",
                      "cartesian_flow", "the plane runs on to $r_s \\to \\infty$, where $\\varepsilon$ vanishes", size,
                      {"u": circles, "v": angles[::3]})
    Xn, Yn = plane.plane()
    ck.add("Natario: the height at every node is G^tt R^2/16", float(np.max(np.abs(plane.Z - height(Xn, Yn)))), 1e-7)
    ck.add("Natario: every triangle of the grid lies on the height, in space", plane.sag() / size, GRID_SAG)
    ck.add("Natario: the energy density is nowhere positive", float(max(0.0, np.max(plane.Z))), 0.0)
    ck.items[-1]["ok"] = bool(np.all(plane.Z <= 0))
    # The wall is deepest beside the ship, across the path, and one R ahead and behind on it.
    from scipy.optimize import minimize_scalar
    side = minimize_scalar(lambda y: float(G(0.0, y)), bounds=(0.5, 1.5), method="bounded", options={"xatol": 1e-12})
    on_path = minimize_scalar(lambda x: float(G(x, 0.0)), bounds=(0.5, 1.5), method="bounded",
                              options={"xatol": 1e-12})
    least = float(side.fun)
    ck.exact(f"Natario: no node is deeper than the wall beside the ship, epsilon = {least / (8 * math.pi):.4f} c^4/GR^2 "
             f"at r_s = {side.x:.4f} R", bool(np.min(plane.Z) >= least / 16 - 1e-7))
    ck.add("Natario: on the path the wall is deepest one R ahead of the ship", abs(on_path.x - 1.0), 1e-6)
    for text, value in (("0.68", least), ("0.48", float(on_path.fun))):
        ck.exact(f"Natario: the caption's {text} c^4/GR^2 is {value / (8 * math.pi):.4f}",
                 f"{-value / (8 * math.pi):.2f}" == text)
    ck.exact(f"Natario: the caption's 0.88R is {side.x:.4f} R", f"{side.x:.2f}" == "0.88")
    U = sp.lambdify((xx, yy), sp.sympify(fns["u"], locals=names).subs(on_plane), "numpy")
    V = sp.lambdify((xx, yy), sp.sympify(fns["v"], locals=names).subs(on_plane), "numpy")
    from scipy.integrate import solve_ivp

    def flow(_, P):
        u, v = float(U(*P)), float(V(*P))
        speed = math.hypot(u, v)
        return [u / speed, v / speed]

    def back(_, P):
        return P[0]
    back.direction = 1
    lines = []
    for y0 in (0.2, 0.4, 0.6, -0.2, -0.4, -0.6):
        # Start just past x = 0 so that the first crossing of x = 0 met is the return.
        run = solve_ivp(flow, (0, 40), [1e-9, y0], rtol=1e-12, atol=1e-13, events=back, dense_output=True,
                        max_step=0.01)
        end = run.t_events[0][0]
        P = run.sol(np.linspace(0, end, 601)).T
        ck.add(f"Natario, the flow line through y = {y0}: it closes", float(np.hypot(*(P[-1] - P[0]))), 1e-8)
        # A field that is axisymmetric about the path and divergence free runs along the level
        # curves of Stokes's stream function, Psi = integral of y u dy from the path, which is
        # taken from the declared field and checked to hold one value along the line.
        psi = [integrate(lambda w, x=x: w * float(U(x, w)), 0.0, y) for x, y in P[::20]]
        ck.add(f"Natario, the flow line through y = {y0}: Stokes's stream function holds one value on it",
               float(np.ptp(psi) / abs(np.mean(psi))), 1e-8)
        lines.append(P[:-1])

    def lifted(P, closed=False):
        """Chart points (x, y) of the plane, cut finely and set on the triangles of the grid."""
        Q = np.column_stack([P, np.zeros(len(P))])
        Q = finely(np.vstack([Q, Q[:1]]) if closed else Q, size / 360)
        Q = Q[:-1] if closed else Q
        Q[:, 2] = plane.facets.height(Q[:, 0], Q[:, 1])
        return Q
    flows = [lifted(P, closed=True) for P in lines]
    P = plane.nodes()
    back_, ahead = plane.index("v", math.pi), plane.index("v", 0.0)
    path = finely(np.vstack([P[::-1, back_], P[1:, ahead]]), size / 360)
    ring = plane.index("u", 1.0)
    wall = finely(np.vstack([P[ring], P[ring, :1]]), size / 360)[:-1]
    for k, Q in enumerate(flows):
        ck.on_piece(f"Natario, flow line {k}", plane, Q)
    for where, Q in (("the path", path), ("the circle r_s = R", wall)):
        ck.on_piece(f"Natario, {where}", plane, Q)
    surface = Surface([plane], curves=[Curve(plane, "path", path), Curve(plane, "wall", wall, closed=True)]
                      + [Curve(plane, "flow", Q, closed=True) for Q in flows])
    fig = figure_of([surface], {"sheet": "cover"}, size, HEIGHT_CAMERA)
    fig.legend("fill", "cover", "$\\varepsilon$ as a height over the plane $z = 0$ of the ship's path at $t = 0$")
    fig.legend("line", "grid", "circles about the ship every $R/4$, and lines from it every $15°$")
    fig.legend("line", "path", "the ship's path, along $x$, the ship heading toward $+x$")
    fig.legend("line", "wall", "$r_s = R$, the middle of the wall")
    fig.legend("line", "flow", "lines of the flow, forward through the bubble and back round it through the wall")
    return [view("plane", "The plane of the path", "$R$", [surface], fig.done(),
                 settings="$t = 0$, when the bubble is centred on $x = 0$, with $R = 1$, the unit of every length.",
                 input="$v_s = 2$, $n = f/2$ with Alcubierre's profile, and the zero expansion field $X = v_s[(2n + "
                       "\\rho n')\\,e_x - n'\\,x_r\\,(x_r, y, z)/\\rho]$, $x_r = x - v_s t$, as in the spacetime "
                       "diagram.",
                 height=f"$\\varepsilon = -c^4K_{{ij}}K^{{ij}}/16\\pi G$, a height of $-R$ for $\\varepsilon = "
                        f"-2c^4/\\pi GR^2$, and ${least / 16:.4f}R$ where it is most negative, "
                        f"$\\varepsilon = {least / (8 * math.pi):.4f}\\,c^4/GR^2$.")]


def at_rest(ck, src, metric_id, system_id):
    """The published Christoffel symbols have no Gamma^i_tt, so a particle at rest in the chart
    stays at rest, and a ring of them keeps its chart positions."""
    gamma, _ = published_christoffel(src, metric_id, system_id)
    ck.exact(f"{metric_id}: no published Gamma^i_tt, so particles at rest in the chart stay there",
             not any(ix[1:] == ("t", "t") and ix[0] != "t" for ix in gamma))


def ring_sequence(ck, src, name, metric_id, system_id, moments, top, params=None, axes=("x", "z")):
    """A flat plane of the two axes, x and z unless named, at each moment, (label, time, fixed,
    functions), with a ring of particles at rest on the unit circle of the chart, which each
    moment stretches into an ellipse of semi-axes sqrt(g_xx) and sqrt(g_zz), checked."""
    size = 2 * top
    alpha = np.linspace(0, 2 * math.pi, 361)
    surfaces = []
    for label, time, fixed_at, functions in moments:
        sl = FlatPlane(src, metric_id, system_id, *axes, fixed_at, params, functions)
        where = f"{name}, {label}"
        plane = disc(sl, "plane", top, "the centre of the ring", "the plane runs on, flat, to infinity", [], size)
        ck.isometry(where, plane)
        ring, dots = particles(ck, where, sl, plane, np.cos(alpha), np.sin(alpha))
        surfaces.append(Surface([plane], label=label, time=time, curves=[ring], dots=dots))
    return surfaces


STACK_CAMERA = Camera(-60, 18)   # a stack of moments seen from low and across both axes of its ellipses


def stack(ck, name, surfaces, rows, lift, planes, size, edge):
    """The moments of a ring of particles on a flat plane stacked up the axis of time into its
    world tube: at every time u from the first moment to the last, the ellipse rows(u) = (a, b)
    of semi-axes along X and Y, at the height lift * u, as an EllipsesPiece whose columns are the
    ring's particles, one every degree, so that the tube's rows at the moments are the rings of
    `surfaces`, which is checked, and its columns their world lines. planes(u) gives the plane of
    each row as a FlatPlane, or the ring's own points, whose semi-axes are checked against the row
    it stands at. The axis of time runs from the first moment to a quarter of the tube's length
    past the last."""
    alpha = np.linspace(0, 2 * math.pi, 361)[:-1]
    times = [s.time for s in surfaces]

    def at(u):
        a, b = rows(np.asarray(u, dtype=float))
        return a, b, lift * np.asarray(u, dtype=float)
    tube = EllipsesPiece("tube", "sheet", times, alpha, at, surfaces[0].pieces[0].sl.metric_id,
                         surfaces[0].pieces[0].sl.system_id, edge, size, [("worldline", "v", alpha[::30])])
    ck.add(f"{name}, the tube's triangles against its ellipses, in space", tube.sag() / size, GRID_SAG)
    worst = 0.0
    for u, a, b in zip(tube.u, tube.A, tube.B):
        worst = max(worst, planes(u, a, b))
    # Within the rounding of the file, and for rings run as geodesics the 1e-9 they are run to.
    ck.add(f"{name}, every row of the tube is the ring at its time", worst, 0.5 * 10.0 ** -tube.decimals + 1e-9)
    N = tube.nodes()
    curves, dots = [], []
    for s in surfaces:
        i = tube.index("u", s.time)
        ring = s.curves[0].points[:360]
        ck.add(f"{name}, {s.label}: the tube's row is the ring of the flat moment, at the height of its time",
               float(max(np.max(np.abs(N[i][:, :2] - ring[:, :2])), abs(N[i][0, 2] - lift * s.time))) / size, FORM)
        curves.append(Curve(tube, "particles", N[i], closed=True, label=s.label, time=s.time))
        dots += [(tube, "particles", Q) for Q in N[i][::30]]
    surface = Surface([tube], curves=curves, dots=dots,
                      axis=(tube.Zr[0], tube.Zr[-1] + (tube.Zr[-1] - tube.Zr[0]) / 4))
    for c in curves:
        ck.on_piece(f"{name}, the ring at {c.label}", tube, c.points)
    ck.on_piece(f"{name}, the particles marked", tube, [Q for _, _, Q in dots])
    return surface


def stack_figure(surface, size, axis_name, legend):
    """A stack of moments on its axis of time, each moment named beside an end of its ring, the
    right and the left in turn from the first up, so that moments close in time stand apart, and
    the axis named at its top."""
    fig = figure_of([surface], {"sheet": "cover"}, size, STACK_CAMERA)
    for i, c in enumerate(surface.curves):
        curve_label(fig, 0, i, c.label, side=1 if i % 2 == 0 else -1)
    axis_label(fig, 0, axis_name)
    for item in legend:
        fig.legend(*item)
    return fig.done()


def kasner(ck, src):
    """The plane y = 0 at t = 1/4, 1/2, 1 and 2, at the exponents (-2/7, 3/7, 6/7) the spacetime
    diagrams declare, each flat, t^(2 p_1) dx^2 + t^(2 p_3) dz^2, with the ring of particles at
    rest on x^2 + z^2 = l^2, which the published Christoffel symbols keep at rest: the ellipse of
    semi-axes t^p_1 l and t^p_3 l, along the direction that contracts and the one that expands
    fastest."""
    p = (sp.Rational(-2, 7), sp.Rational(3, 7), sp.Rational(6, 7))
    ck.exact("Kasner: the exponents sum to 1, and so do their squares", sum(p) == 1 and sum(q * q for q in p) == 1)
    at_rest(ck, src, "kasner", "cartesian")
    params = {"p_1": "-2/7", "p_2": "3/7", "p_3": "6/7"}
    moments = [(f"$t = {s}$", float(sp.Rational(s)), {"t": s, "y": 0}, None) for s in ("1/4", "1/2", "1", "2")]
    surfaces = ring_sequence(ck, src, "Kasner", "kasner", "cartesian", moments, 2.0, params)
    for s, (_, time, _, _) in zip(surfaces, moments):
        P = s.curves[0].points
        a = np.linspace(0, 2 * math.pi, 361)
        want = np.column_stack([time ** (-2 / 7) * np.cos(a), time ** (6 / 7) * np.sin(a)])
        ck.add(f"Kasner, t = {time}: the ellipse of semi-axes t^p_1 and t^p_3", float(np.max(np.abs(P[:, :2] - want))), 1e-12)
    def kasner_plane(u, a, b):
        sl = FlatPlane(src, "kasner", "cartesian", "x", "z", {"t": repr(float(u)), "y": 0}, params)
        return float(np.max(np.abs(sl.scale - [a, b])))
    tube = stack(ck, "Kasner", surfaces, lambda u: (u ** (-2 / 7), u ** (6 / 7)), 1.5, kasner_plane, 5.0,
                 "the ring's world tube runs on before $t = 1/4$ and after $t = 2$")
    tube_fig = stack_figure(tube, 5.0, "$t$", [
        ("fill", "cover", "the ring at every moment from $t = 1/4$ to $t = 2$, each at the height of its time"),
        ("line", "particles", "the ring at the four moments of the flat view, twelve of its particles marked"),
        ("line", "worldline", "the world lines of the twelve particles, at rest in the chart"),
        ("line", "axis", "the axis of time, through the centre of the ring")])
    fig = sequence_figure(surfaces, {"sheet": "cover"}, 4.0, columns=2, camera=FLAT_CAMERA, meridians=12)
    fig.legend("fill", "cover", "the plane $y = 0$ at each moment, flat")
    fig.legend("line", "particles", "a ring of particles at rest in the chart on $x^2 + z^2 = \\ell^2$, with twelve of them "
                                    "marked: an ellipse reaching $t^{p_1}\\ell$ along $x$ and $t^{p_3}\\ell$ along $z$")
    fig.legend("line", "meridian", "straight lines from the centre, every $30°$")
    settings = ("$(p_1, p_2, p_3) = (-2/7, 3/7, 6/7)$ and $t$ in the unit of time in which the powers are evaluated, "
                "with $\\ell$ the ring's radius at $t = 1$, the unit of every length.")
    given = ("Exponents $(p_1, p_2, p_3) = (-2/7, 3/7, 6/7)$, a point on the Kasner circle, as in the spacetime "
             "diagrams.")
    return [view("tube", "The ring's world tube", "$\\ell$", [tube], tube_fig, settings=settings, input=given,
                 height="$t$, a height of $1.5\\,\\ell$ for each unit of $t$"),
            view("ring", "A ring of particles", "$\\ell$", surfaces, fig.done(),
                 settings="$(p_1, p_2, p_3) = (-2/7, 3/7, 6/7)$ and $t$ in the unit of time in which the powers are "
                          "evaluated, with $\\ell$ the ring's radius at $t = 1$, the unit of every length.",
                 input="Exponents $(p_1, p_2, p_3) = (-2/7, 3/7, 6/7)$, a point on the Kasner circle, as in the "
                       "spacetime diagrams.")]


def bianchi(ck, src):
    """The plane y = 0 at four moments of the dust the spacetime diagram declares, solved from
    this spacetime's own G^x_x = G^y_y = G^z_z = 0 by null_rays.DustSolver: flat at each,
    a_1^2 dx^2 + a_3^2 dz^2, with the ring of dust on x^2 + z^2 = l^2, at rest in the chart as the
    published Christoffel symbols keep it: the ellipse of semi-axes a_1 l and a_3 l."""
    solver = nr.DustSolver("bianchi", "type_i_cartesian", "t", **nr.BIANCHI_DUST)
    src.note("bianchi", "type_i_cartesian", ["einstein_tensor"])
    at_rest(ck, src, "bianchi", "type_i_cartesian")
    moments = []
    for t in (0.1, solver.t_ref, 1.0, 2.0):
        y = solver.state(np.array([t]))[:, 0]
        functions = {"a_1": repr(float(y[0])), "a_2": repr(float(y[2])), "a_3": repr(float(y[4]))}
        moments.append((f"$c\\bar Ht = {t:.2f}$", float(t), {"t": repr(float(t)), "y": 0}, functions))
    surfaces = ring_sequence(ck, src, "Bianchi I", "bianchi", "type_i_cartesian", moments, 3.5)
    for s, (_, t, _, fn) in zip(surfaces, moments):
        P = s.curves[0].points
        a = np.linspace(0, 2 * math.pi, 361)
        want = np.column_stack([float(fn["a_1"]) * np.cos(a), float(fn["a_3"]) * np.sin(a)])
        ck.add(f"Bianchi I, t = {t:.4f}: the ellipse of semi-axes a_1 and a_3", float(np.max(np.abs(P[:, :2] - want))), 1e-12)
    def bianchi_rows(u):
        y = solver.state(np.atleast_1d(u))
        return y[0], y[4]

    def bianchi_plane(u, a, b):
        y = solver.state(np.array([u]))[:, 0]
        sl = FlatPlane(src, "bianchi", "type_i_cartesian", "x", "z", {"t": repr(float(u)), "y": 0}, None,
                       {"a_1": repr(float(y[0])), "a_2": repr(float(y[2])), "a_3": repr(float(y[4]))})
        return float(np.max(np.abs(sl.scale - [a, b])))
    tube = stack(ck, "Bianchi I", surfaces, bianchi_rows, 2.5, bianchi_plane, 7.0,
                 "the world tube of the ring runs on before $c\\bar Ht = 0.10$ and after $c\\bar Ht = 2.00$")
    tube_fig = stack_figure(tube, 7.0, "$c\\bar Ht$", [
        ("fill", "cover", "the ring of dust at every moment from $c\\bar Ht = 0.10$ to $2.00$, each at the height of its "
                          "time"),
        ("line", "particles", "the ring at the four moments of the flat view, twelve of its grains marked"),
        ("line", "worldline", "the world lines of the twelve grains, at rest in the chart"),
        ("line", "axis", "the axis of time, through the centre of the ring")])
    fig = sequence_figure(surfaces, {"sheet": "cover"}, 7.0, columns=2, camera=FLAT_CAMERA, meridians=12)
    fig.legend("fill", "cover", "the plane $y = 0$ at each moment, flat")
    fig.legend("line", "particles", "a ring of the dust on $x^2 + z^2 = \\ell^2$, with twelve of its grains marked: an "
                                    "ellipse reaching $a_1\\ell$ along $x$ and $a_3\\ell$ along $z$")
    fig.legend("line", "meridian", "straight lines from the centre, every $30°$")
    settings = ("$t$ in units of $1/\\bar H$ from the singularity, with $\\ell$ the ring's radius where "
                "$a_1 = a_2 = a_3 = 1$, the unit of every length.")
    given = ("Dust: the three scale factors solved from this spacetime's own $G^x{}_x = G^y{}_y = G^z{}_z = 0$, "
             "starting from $a_i = 1$ with rates $(-0.5, 1.5, 2.0)\\,\\bar H$, $\\bar H$ their mean, as in the "
             "spacetime diagram.")
    return [view("tube", "The ring's world tube", "$\\ell$", [tube], tube_fig, settings=settings, input=given,
                 height="$c\\bar Ht$, a height of $2.5\\,\\ell$ for each $1/\\bar H$"),
            view("ring", "A ring of dust", "$\\ell$", surfaces, fig.done(), settings=settings, input=given)]


def pp_wave(ck, src):
    """The wave front at four values of u of the pulse A = exp(-u^2), B = 0 the spacetime diagram
    declares, in units of L: a surface of constant u has the metric dx^2 + dy^2 whatever v is on
    it, flat. The ring is 360 free particles at rest on x^2 + y^2 = L^2 before the pulse, each run
    from cu = -8 with the published Christoffel symbols, u being affine on every geodesic since no
    Gamma^u is published: x'' = -Gamma^x_uu and y'' = -Gamma^y_uu. The pulse stretches the ring
    along x and squeezes it along y until every particle reaches the x axis at once."""
    gamma, R = published_christoffel(src, "pp_wave", "exact_plane_wave")
    ck.exact("pp-wave: no published Gamma^u, so u is an affine parameter", not any(ix[0] == "u" for ix in gamma))
    ck.exact("pp-wave: the only published Gamma^x and Gamma^y are Gamma^x_uu and Gamma^y_uu",
             {ix for ix in gamma if ix[0] in ("x", "y")} == {("x", "u", "u"), ("y", "u", "u")})
    names = {R._plain(n): s for n, s in R.symbol.items()}
    decl = {R.parameters["A"]: sp.exp(-names["u"] ** 2), R.parameters["B"]: sp.Integer(0)}
    args = (names["u"], names["x"], names["y"])
    gx = sp.lambdify(args, gamma[("x", "u", "u")].subs(decl).doit(), "numpy")
    gy = sp.lambdify(args, gamma[("y", "u", "u")].subs(decl).doit(), "numpy")
    from scipy.integrate import solve_ivp
    from scipy.optimize import brentq
    alpha = np.linspace(0, 2 * math.pi, 361)[:-1]
    n = len(alpha)

    def rhs(u, w):
        x, y = w[:n], w[n:2 * n]
        return np.concatenate([w[2 * n:3 * n], w[3 * n:], -gx(u, x, y) * np.ones(n), -gy(u, x, y) * np.ones(n)])
    start = np.concatenate([np.cos(alpha), np.sin(alpha), np.zeros(2 * n)])
    run = solve_ivp(rhs, (-8, 2), start, rtol=1e-12, atol=1e-14, dense_output=True)
    focus = brentq(lambda u: run.sol(u)[n + 90], 0, 1.5, xtol=1e-15)
    ck.add("pp-wave: at the focus every particle is on the x axis", float(np.max(np.abs(run.sol(focus)[n:2 * n]))), 1e-9)
    top = 3.0
    surfaces = []
    for u in (-3.0, -0.5, 0.0, focus):
        sl = FlatPlane(src, "pp_wave", "exact_plane_wave", "x", "y", {"u": repr(u), "v": 0},
                       functions={"A": "exp(-u**2)", "B": "0"})
        where = f"pp-wave, cu = {u:.4f}"
        plane = disc(sl, "plane", top, "the centre of the ring", "the wave front runs on, flat, to infinity", [], 2 * top)
        ck.isometry(where, plane)
        w = run.sol(u)
        cx, cy = np.append(w[:n], w[0]), np.append(w[n:2 * n], w[n])
        X, Y = w[0], w[n + 90]
        ck.add(f"{where}: the ring is the ellipse X cos(alpha), Y sin(alpha)",
               float(np.max(np.abs(np.column_stack([cx, cy]) - np.column_stack([X * np.cos(np.append(alpha, 0)),
                                                                                  Y * np.sin(np.append(alpha, 0))])))), 1e-9)
        ring, dots = particles(ck, where, sl, plane, cx, cy)
        label = "$cu = " + (f"{u:g}" if u != focus else f"{u:.2f}") + "\\,L$"
        surfaces.append(Surface([plane], label=label, time=u, curves=[ring], dots=dots))
    def pp_rows(u):
        w = run.sol(np.atleast_1d(u))
        return w[0], w[n + 90]

    def pp_ring(u, a, b):
        # Every particle of the ring, run with the published Christoffel symbols, on the row.
        w = run.sol(u)
        return float(np.max(np.abs(np.column_stack([w[:n], w[n:2 * n]]) - np.column_stack([a * np.cos(alpha),
                                                                                        b * np.sin(alpha)]))))
    tube = stack(ck, "pp-wave", surfaces, pp_rows, 0.5, pp_ring, 2 * top,
                 "the world tube of the ring runs on before $cu = -3\\,L$, and past the focus the ring turns inside out")
    tube_fig = stack_figure(tube, 2 * top, "$cu$", [
        ("fill", "cover", "the ring at every wave front from $cu = -3\\,L$ to the focus, each at the height of its $u$"),
        ("line", "particles", "the ring at the four wave fronts of the flat view, twelve of its particles marked"),
        ("line", "worldline", "the world lines of the twelve particles"),
        ("line", "axis", "the axis of retarded time, through the centre of the ring")])
    fig = sequence_figure(surfaces, {"sheet": "cover"}, 2 * top, columns=2, camera=FLAT_CAMERA, meridians=12)
    fig.legend("fill", "cover", "the wave front at each moment, flat")
    fig.legend("line", "particles", "a ring of free particles at rest on $x^2 + y^2 = L^2$ before the pulse, with twelve of them "
                                    "marked")
    fig.legend("line", "meridian", "straight lines from the centre, every $30°$")
    settings = "$L = 1$, the unit of every length and of $cu$; each moment is the wave front of one $u$."
    given = "A pulse of the plus polarisation, $A = e^{-u^2}/L^2$ and $B = 0$, as in the spacetime diagram."
    return [view("tube", "The ring's world tube", "$L$", [tube], tube_fig, settings=settings, input=given,
                 height="$cu$, a height of $L/2$ for each $L$"),
            view("ring", "A ring of particles", "$L$", surfaces, fig.done(), settings=settings, input=given)]


def aichelburg_sexl(ck, src):
    """The wave front at four values of u of the null Cartesian chart, u = -1, 0.5, 1 and 1.5 in
    units of 8GE/c^4 = rho_0 = 1, with the delta of the shock drawn as the pulse its spacetime
    diagram declares: a surface of constant u has the metric dx^2 + dy^2 whatever v is on it, flat.
    The ring is 360 free particles at rest on x^2 + y^2 = 1 before the shock, the circle of radius
    8GE/c^4 about the axis the source moves along, each run from u = -1 with the published
    Christoffel symbols, u being affine on every geodesic since no Gamma^u is published:
    x'' = -Gamma^x_uu and y'' = -Gamma^y_uu. Across the delta each particle gains the velocity
    -(4GE/c^4) x/rho^2 toward the axis, so the ring stays a circle of radius 1 - u/2 behind the
    shock and closes on the axis at u = 2. The particles move while they cross a pulse of width w,
    which moves the focus by the order of w, so the delta is drawn here as a pulse a tenth as wide
    as the spacetime diagram's, whose width no drawing shows, and every check holds to twice it."""
    width = 0.005
    gamma, R = published_christoffel(src, "aichelburg_sexl", "null_cartesian")
    ck.exact("Aichelburg-Sexl: no published Gamma^u, so u is an affine parameter", not any(ix[0] == "u" for ix in gamma))
    ck.exact("Aichelburg-Sexl: the only published Gamma^x and Gamma^y are Gamma^x_uu and Gamma^y_uu",
             {ix for ix in gamma if ix[0] in ("x", "y")} == {("x", "u", "u"), ("y", "u", "u")})
    names = {R._plain(n): s for n, s in R.symbol.items()}
    decl = {R.parameters["G"]: 1, R.parameters["E"]: sp.Rational(1, 8), R.parameters["rho_0"]: 1}
    args = (names["u"], names["x"], names["y"])
    pulse = f"exp(-s**2/{width}**2)/({width}*sqrt(pi))"
    gx = sp.lambdify(args, nr.smoothed(gamma[("x", "u", "u")], pulse).subs(decl), "numpy")
    gy = sp.lambdify(args, nr.smoothed(gamma[("y", "u", "u")], pulse).subs(decl), "numpy")
    from scipy.integrate import solve_ivp
    from scipy.optimize import brentq
    alpha = np.linspace(0, 2 * math.pi, 361)[:-1]
    n = len(alpha)

    def rhs(u, w):
        x, y = w[:n], w[n:2 * n]
        return np.concatenate([w[2 * n:3 * n], w[3 * n:], -gx(u, x, y) * np.ones(n), -gy(u, x, y) * np.ones(n)])
    start = np.concatenate([np.cos(alpha), np.sin(alpha), np.zeros(2 * n)])
    run = solve_ivp(rhs, (-1, 2.5), start, rtol=1e-12, atol=1e-14, dense_output=True, max_step=width / 5)

    def radius(u):
        w = run.sol(np.atleast_1d(u))
        return np.hypot(w[0], w[n])
    # The particle that starts at (1, 0) moves along the x axis, and reaches the axis x = 0 at the focus.
    focus = brentq(lambda u: float(run.sol(u)[0]), 1.5, 2.2, xtol=1e-12)
    ck.add("Aichelburg-Sexl: the ring closes on the axis at u = 2, to twice the width of the pulse", abs(focus - 2.0), 2 * width)
    ck.add("Aichelburg-Sexl: behind the pulse the ring shrinks as 1 - u/2", float(np.max(np.abs(radius(np.linspace(0.3, 1.5, 50))
                                                                                                - (1 - np.linspace(0.3, 1.5, 50) / 2)))), 2 * width)
    top = 1.5
    surfaces = []
    for u in (-1.0, 0.5, 1.0, 1.5):
        sl = FlatPlane(src, "aichelburg_sexl", "null_cartesian", "x", "y", {"u": repr(u), "v": 0}, params={"G": 1, "E": "1/8", "rho_0": 1})
        where = f"Aichelburg-Sexl, u = {u:g}"
        plane = disc(sl, "plane", top, "the axis the source moves along", "the wave front runs on, flat, to infinity", [], 2 * top)
        ck.isometry(where, plane)
        w = run.sol(u)
        cx, cy = np.append(w[:n], w[0]), np.append(w[n:2 * n], w[n])
        a = float(np.hypot(w[0], w[n]))
        ck.add(f"{where}: the ring is the circle of radius {a:.6f} about the axis",
               float(np.max(np.abs(np.column_stack([cx, cy]) - a * np.column_stack([np.cos(np.append(alpha, 0)),
                                                                                    np.sin(np.append(alpha, 0))])))), 1e-9)
        ring, dots = particles(ck, where, sl, plane, cx, cy)
        surfaces.append(Surface([plane], label=f"$u = {u:g}$", time=u, curves=[ring], dots=dots))

    def as_rows(u):
        w = run.sol(np.atleast_1d(u))
        return w[0], w[n + 90]

    def as_ring(u, a, b):
        # Every particle of the ring, run with the published Christoffel symbols, on the row.
        w = run.sol(u)
        return float(np.max(np.abs(np.column_stack([w[:n], w[n:2 * n]]) - np.column_stack([a * np.cos(alpha),
                                                                                        b * np.sin(alpha)]))))
    tube = stack(ck, "Aichelburg-Sexl", surfaces, as_rows, 0.5, as_ring, 2 * top,
                 "the world tube of the ring runs on before $u = -1$, and closes on the axis at $u = 2$")
    tube_fig = stack_figure(tube, 2 * top, "$u$", [
        ("fill", "cover", "the ring at every wave front from $u = -1$ to $u = 1.5$, each at the height of its $u$"),
        ("line", "particles", "the ring at the four wave fronts of the flat view, twelve of its particles marked"),
        ("line", "worldline", "the world lines of the twelve particles"),
        ("line", "axis", "the axis the source moves along, crossed by the shock at $u = 0$")])
    fig = sequence_figure(surfaces, {"sheet": "cover"}, 2 * top, columns=2, camera=FLAT_CAMERA, meridians=12)
    fig.legend("fill", "cover", "the wave front at each moment, flat")
    fig.legend("line", "particles", "a ring of free particles at rest on $x^2 + y^2 = (8GE/c^4)^2$ before the shock, "
                                    "with twelve of them marked")
    fig.legend("line", "meridian", "straight lines from the axis, every $30°$")
    settings = ("$8GE/c^4 = \\rho_0 = 1$, the unit of every length and of $u$; each moment is the wave front of one "
                "$u$.")
    given = ("$\\delta(u)$ drawn as the pulse $e^{-u^2/w^2}/(w\\sqrt{\\pi})$ with $w = 0.005$, the field of a pulse of "
             "light of that length in $u$ along the axis, carrying the same energy $E$.")
    return [view("tube", "The ring's world tube", "$8GE/c^4$", [tube], tube_fig, settings=settings, input=given,
                 height="$u$, a height of $1/2$ for each unit of $u$"),
            view("ring", "A ring of particles", "$8GE/c^4$", surfaces, fig.done(), settings=settings, input=given)]


def khan_penrose(ck, src):
    """The wave front, the plane of x and y, at the four events tau = 0, 0.6, 1.0 and 1.3 on
    sigma = 0 of the cosmological chart, where u = v, at L = 1: each flat, g_xx dx^2 + g_yy dy^2
    with g_xx = (1 + sin tau)^2/cos tau and g_yy = (1 - sin tau)^2/cos tau, with the ring of free
    particles at rest on x^2 + y^2 = l^2, which the published Christoffel symbols keep at rest
    on sigma = 0: they have no Gamma^x_tautau or Gamma^y_tautau, and Gamma^sigma_tautau vanishes
    there. At tau = 0, the collision, the ring is the circle, and toward the singularity at
    tau = pi/2 it is drawn out along x and squeezed along y."""
    gamma, R = published_christoffel(src, "khan_penrose", "cosmological")
    ck.exact("Khan-Penrose: no published Gamma^x_tautau or Gamma^y_tautau",
             not any(ix[0] in ("x", "y") and ix[1:] == ("\\tau", "\\tau") for ix in gamma))
    sigma_tt = gamma.get(("\\sigma", "\\tau", "\\tau"), sp.Integer(0))
    ck.exact("Khan-Penrose: Gamma^sigma_tautau vanishes on sigma = 0", sigma_tt.subs(R.symbol["\\sigma"], 0) == 0)
    times = (0.0, 0.6, 1.0, 1.3)
    moments = [(f"$\\tau = {t:g}$", t, {"tau": repr(t), "sigma": 0}, None) for t in times]
    surfaces = ring_sequence(ck, src, "Khan-Penrose", "khan_penrose", "cosmological", moments, 4.0, {"L": 1},
                             axes=("x", "y"))

    def rows(tau):
        tau = np.asarray(tau, dtype=float)
        return (1 + np.sin(tau)) / np.sqrt(np.cos(tau)), (1 - np.sin(tau)) / np.sqrt(np.cos(tau))
    a = np.linspace(0, 2 * math.pi, 361)
    for s, t in zip(surfaces, times):
        A, B = rows(t)
        ck.add(f"Khan-Penrose, tau = {t}: the ellipse of semi-axes (1 +- sin tau)/sqrt(cos tau)",
               float(np.max(np.abs(s.curves[0].points[:, :2] - np.column_stack([A * np.cos(a), B * np.sin(a)])))), 1e-12)

    def plane(tau, A, B):
        sl = FlatPlane(src, "khan_penrose", "cosmological", "x", "y", {"tau": repr(float(tau)), "sigma": 0}, {"L": 1})
        return float(np.max(np.abs(sl.scale - [A, B])))
    tube = stack(ck, "Khan-Penrose", surfaces, rows, 3.0, plane, 8.0,
                 "the world tube of the ring runs on before the collision, a cylinder, and after $\\tau = 1.3$ to the "
                 "singularity")
    tube_fig = stack_figure(tube, 8.0, "$\\tau$", [
        ("fill", "cover", "the ring at every moment from the collision, $\\tau = 0$, to $\\tau = 1.3$, each at the height "
                          "of its $\\tau$"),
        ("line", "particles", "the ring at the four moments of the flat view, twelve of its particles marked"),
        ("line", "worldline", "the world lines of the twelve particles, at rest in the chart"),
        ("line", "axis", "the axis of $\\tau$, through the centre of the ring")])
    fig = sequence_figure(surfaces, {"sheet": "cover"}, 8.0, columns=2, camera=FLAT_CAMERA, meridians=12)
    fig.legend("fill", "cover", "the wave front at each moment, flat")
    fig.legend("line", "particles", "a ring of free particles at rest on $x^2 + y^2 = \\ell^2$, with twelve of them marked: an "
                                    "ellipse reaching $(1 + \\sin\\tau)\\ell/\\sqrt{\\cos\\tau}$ along $x$ and "
                                    "$(1 - \\sin\\tau)\\ell/\\sqrt{\\cos\\tau}$ along $y$")
    fig.legend("line", "meridian", "straight lines from the centre, every $30°$")
    settings = ("$L = 1$, each moment the wave front at one $\\tau$ on $\\sigma = 0$, where $u = v$, with $\\ell$ the "
                "ring's radius before the collision, the unit of every length.")
    return [view("tube", "The ring's world tube", "$\\ell$", [tube], tube_fig, settings=settings,
                 height="$\\tau$, a height of $3\\,\\ell$ for each unit of $\\tau$"),
            view("ring", "A ring of particles", "$\\ell$", surfaces, fig.done(), settings=settings)]


def malament_hogarth(ck, src):
    """The plane z = 0 about the removed event at ct = -0.7, -0.3, -0.1 and 0, with the conformal
    factor the spacetime diagram declares, which depends only on c^2t^2 + x^2 + y^2 + z^2: turned
    about the origin the plane is a surface of revolution with g_ss = Omega^2 and rho = s Omega,
    and -s Omega'(2 Omega + s Omega') >= 0 everywhere, so it is drawn whole: a well inside the
    unit ball, where Omega > 1, flat outside. At ct = 0 the well has no bottom: rho -> 1 while
    the distance down grows as ln(1/s), a tube, drawn to s = 0.03, and its length from the rim to
    s is the proper time of the computer on the axis from ct = -1 to -s, which is checked."""
    fn = {"Omega": nr._MH_FACTOR}
    top = 1.5
    size = 2 * top

    def moment(T):
        t = float(sp.Rational(T))
        sl = Slice(src, "malament_hogarth", "cartesian", "x", None, {"t": T, "z": 0}, functions=fn, turn="y")
        # The edge of the region where Omega > 1, sqrt(1 - c^2t^2) correctly rounded, moved on by
        # a unit in the last place while the declared Omega's two branches disagree there, since
        # within one unit of the edge the float c^2t^2 + s^2 can fall below 1 while s^2 - (1 -
        # c^2t^2) rounds above 0, and the inside branch then divides by a number of the wrong sign.
        rim = float(sp.sqrt(1 - sp.Rational(T) ** 2))
        while not float(sl.gxx_at(rim)) == 1.0:
            rim = float(np.nextafter(rim, 2.0))
        where = f"Malament-Hogarth, ct = {T}"
        lo = 0.03 if t == 0 else 0.0
        start = (("edge", "the tube runs on for ever toward the removed event") if t == 0
                 else ("axis", "the centre, straight below the removed event in time"))
        marks = [(s, "r", None) for s in (0.1, 0.3, 0.5) if lo < s < rim] + [(rim, "surface", None)]
        well = Piece("well", "star", sl, lo, rim, 0.0, 1, (start, ("join", "the edge of the region where $\\Omega > 1$")),
                     marks, size)
        # The flat plane beyond the unit ball stands at z = 0 at every moment, and the well sinks below it.
        well.z = well.z - well.z[-1]
        flat = Piece("flat", "sheet", sl, rim, top, 0.0, 1,
                     (("join", "the edge of the region where $\\Omega > 1$"), ("edge", "the plane runs on, flat")),
                     [(1.0, "r", None)] if rim < 1.0 - 1e-9 else [], size)
        ck.isometry(f"{where}, the well", well)
        ck.isometry(f"{where}, the flat plane", flat)
        ck.join(f"{where}, the well meets the flat plane", well, rim, flat, rim)
        ck.plane(f"{where}, beyond the region where Omega > 1", sl, np.linspace(rim + 1e-6, 3, 300))
        if t == 0:
            _, entry, R = nr.load("malament_hogarth", "cartesian")
            g = nr.published_matrix(R, entry, "metric_components")
            names = {R._plain(n): s for n, s in R.symbol.items()}
            omega = sp.sympify(nr._MH_FACTOR, locals=names)
            gtt = g[0, 0].subs(R.parameters["Omega"], omega).subs({names["x"]: 0, names["y"]: 0, names["z"]: 0, R.c: 1})
            clock = sp.lambdify(names["t"], sp.sqrt(-gtt), "numpy")
            for s in (0.03, 0.1, 0.3):
                length = sl.proper(s, 1.0)
                ticks = integrate(lambda v: float(clock(v)), -1.0, -s)
                ck.add(f"Malament-Hogarth: the tube from its rim down to s = {s} is as long as the computer's clock "
                       f"runs from ct = -1 to -{s}", abs(length - ticks) / ticks, 1e-10)
        return Surface([well, flat], label=f"$ct = {t:g}$", time=t)

    # A frame every 1/80 of ct from -0.7, which holds each moment of the flat views, up to the last
    # whose well is no deeper than the tube is drawn at ct = 0: the wells just before it reach
    # deeper still, while the bottomless tube is drawn only to s = 0.03.
    last = moment("0")
    depth = -last.pieces[0].z[0]
    frames = []
    for k in range(-56, 0):
        f = moment(str(sp.Rational(k, 80)))
        if -f.pieces[0].z[0] > depth:
            break
        frames.append(f)
    frames.append(last)
    surfaces = [f for f in frames if f.time in (-0.7, -0.3, -0.1, 0.0)]
    ck.exact("Malament-Hogarth: the movie holds every moment of the flat views", len(surfaces) == 4)
    fig = movie_figure(frames, {"star": "star", "sheet": "cover"}, size, Camera(-90, 22), meridians=12)
    fig.legend("fill", "star", "inside the unit ball about the removed event, where $\\Omega > 1$")
    fig.legend("fill", "cover", "outside it, where $\\Omega = 1$ and the plane is flat")
    fig.legend("line", "surface", "the edge of the region where $\\Omega > 1$")
    fig.legend("line", "r", "$s$ constant, at $0.1$, $0.3$, $0.5$ and $1$")
    fig.legend("line", "meridian", "$\\phi$ constant, every $30°$")
    return [view("plane", "Toward the removed event", "$1$", surfaces, fig.done(),
                 movie=movie(frames, "$ct$", [f.time for f in frames]),
                 settings="$c\\,t$ and every length in the unit in which the declared $\\Omega$ is written, the radius "
                          "of the region where $\\Omega > 1$; each moment is the plane $z = 0$ of one $t$.",
                 input="$\\Omega = 1 + e^{1 - 1/(1 - \\varrho^2)}/\\varrho$ for $\\varrho^2 = c^2t^2 + x^2 + y^2 + z^2 < 1$ and "
                       "$\\Omega = 1$ beyond, as in the spacetime diagram.",
                 stops=["At $ct = 0$ the tube runs on without end toward the removed event, and is drawn down to "
                        "$s = 0.03$."])]


def anti_de_sitter(ck, src):
    """The static slice's equator, g_rr = 1/(1 + r^2/L^2) and g_phiphi = r^2, at L = 1: its circles
    grow faster than the distance out to them at every r > 0, which is checked, so no surface of
    revolution in flat space carries it. In three dimensional Minkowski space, dX^2 + dY^2 - dZ^2,
    it climbs at dZ/dr = sqrt((drho/dr)^2 - g_rr) = r/sqrt(1 + r^2): the hyperboloid Z = sqrt(1 +
    r^2) - 1, on which the whole hyperbolic plane lies, drawn to r = 4 with the light cone it nears
    as a reference."""
    sl = Slice(src, "anti_de_sitter", "static_global", "r", "\\phi", {"t": 0, **EQUATOR}, {"L": 1}, space="minkowski")
    ck.stops("anti-de Sitter, the equator of the static chart in flat space", sl, np.linspace(1e-3, 20, 400))
    top = 4.0
    size = 2 * top
    sheet = Piece("sheet", "sheet", sl, 0.0, top, 0.0, 1,
                  (("axis", "the centre $r = 0$"), ("edge", "the sheet runs on toward the light cone, to $r \\to \\infty$")),
                  [(r, "r", None) for r in (1.0, 2.0, 3.0, top)], size)
    cone = FormPiece("cone", sl, np.linspace(0.0, top, 81), lambda r: r, lambda r: r - 1.0,
                     (("apex", "the apex of the light cone, a distance $L$ below the centre"),
                      ("edge", "the cone runs on")), size)
    ck.isometry("anti-de Sitter, the sheet", sheet)
    ck.form("anti-de Sitter, the hyperboloid Z = sqrt(L^2 + r^2) - L", sheet, lambda r: np.sqrt(1 + r * r) - 1, size)
    ck.radius("anti-de Sitter, rho = r", sheet, lambda r: r, size)
    surface = Surface([sheet, cone])
    fig = figure_of([surface], {"sheet": "cover"}, size)
    ring_label(fig, [0, 0, 0], *sheet.at(1.0), "$r = L$")
    ring_label(fig, [0, 0, 0], *sheet.at(top), "$4L$")
    fig.legend("fill", "cover", "the equator of the static slice at $t = 0$, which $t$ and $r$ cover whole")
    fig.legend("line", "r", "$r$ constant, at $L$, $2L$, $3L$ and $4L$, a proper distance $0.88$, $0.56$, $0.37$ and "
                            "$0.28\\,L$ apart")
    fig.legend("line", "reference", "the light cone of the Minkowski space it is drawn in, which the sheet nears as "
                                    "$r \\to \\infty$")
    fig.legend("line", "meridian", "$\\phi$ constant, every $15°$")
    return [view("hyperboloid", "In Minkowski space", "$L$", [surface], fig.done(),
                 settings="$L = 1$, the unit of every length, and every length along the sheet measured with "
                          "$dX^2 + dY^2 - dZ^2$.",
                 stops=["At every $r > 0$ the circles grow faster than the distance out to them, $g_{rr} < "
                        "(\\partial_r\\sqrt{g_{\\phi\\phi}})^2$, and no surface of revolution in flat space carries the "
                        "slice; Minkowski space carries it."])]


def btz(ck, src):
    """The moment t = 0 of the hole without rotation, M = 1 and J = 0 at l = 1, as its conformal
    diagram's square draws it, through both exteriors as one surface. g_rr = 1/N^2 with N^2 =
    r^2 - 1 and g_phiphi = r^2, so in flat space dz/dr = sqrt(1/N^2 - 1) = sqrt((2 - r^2)/(r^2 -
    1)): the slice runs through the bifurcation circle r_+ = 1, its throat, into the second
    exterior, and lies level at r = sqrt 2, where N^2 = 1. Farther out 1/N^2 < 1, the circles grow
    faster than the distance out to them, which is checked, and the slice is drawn on in three
    dimensional Minkowski space, climbing at dZ/dr = sqrt(1 - 1/N^2) from level at sqrt 2 toward a
    light cone of that space, as anti-de Sitter space's own static slice does; between r_+ and
    sqrt 2 no surface of Minkowski space carries it, since there 1/N^2 > 1, which is checked as
    well. Each exterior is a piece in flat space from the throat to sqrt 2 and a piece in
    Minkowski space from there to r = 3, started at the height where the first ends; both lie
    level at sqrt 2, so they meet in one circle with one tangent, which is checked from the
    numbers and from the file, and the two exteriors meet at the throat, which is checked too.
    The drawing stops at r = 3, where dZ/dr = sqrt(7/8) and the sheet runs close to its cone,
    so that the throat keeps a third of the width."""
    params = {"ell": 1, "M": 1, "J": 0}
    sl = Slice(src, "btz", "stationary", "r", "\\phi", {"t": 0}, params)
    msl = Slice(src, "btz", "stationary", "r", "\\phi", {"t": 0}, params, space="minkowski")
    rp = sl.horizons()[0]
    level = math.sqrt(2.0)
    top = 3.0
    size = 2 * top
    ck.add("BTZ: the horizon is at r_+ = l", abs(rp - 1.0), 1e-12)
    ck.add("BTZ: the surface lies level where N^2 = 1, at r = sqrt(1 + M) l",
           abs(float(sl.defect_at(np.array([level]))[0])), 1e-12)
    ck.stops("BTZ, beyond r = sqrt(1 + M) l in flat space", sl, np.linspace(level, 40, 402)[1:])
    inward = np.linspace(rp, level, 402)[1:-1]
    ck.add("BTZ: between r_+ and sqrt(1 + M) l no surface in Minkowski space carries the slice, (drho/dr)^2 - g_rr < 0",
           float(max(0.0, np.max(-msl.defect_at(inward)))), 0.0)
    if not np.all(msl.defect_at(inward) > 0):
        ck.items[-1]["ok"] = False

    join = "at $r = \\sqrt{1 + M}\\,\\ell$ the surface lies level, in flat space nearer the horizon and in Minkowski space beyond"
    edge = "the sheet runs on toward a light cone, to $r \\to \\infty$"
    near = Piece("exterior", "sheet", sl, rp, level, 0.0, 1,
                 (("throat", "the throat $r = r_+$, the bifurcation circle, where the other exterior begins"), ("join", join)),
                 [(rp, "horizon", "$r = r_+$"), (1.2, "r", None), (level, "space", None)], size)
    out = Piece("exterior_minkowski", "sheet", msl, level, top, near.at(level)[1], 1, (("join", join), ("edge", edge)),
                [(2.0, "r", None), (top, "r", None)], size)
    far = Piece("other_exterior", "sheet2", sl, rp, level, 0.0, -1, (("throat", "the throat $r = r_+$"), ("join", join)),
                [(1.2, "r2", None), (level, "space", None)], size)
    far_out = Piece("other_exterior_minkowski", "sheet2", msl, level, top, far.at(level)[1], -1,
                    (("join", join), ("edge", edge)), [(2.0, "r2", None), (top, "r2", None)], size)
    for p in (near, out, far, far_out):
        space = "in Minkowski space" if p.sl.lorentz else "in flat space"
        ck.isometry(f"BTZ, {p.id} {space}", p)
        ck.radius(f"BTZ, {p.id}, rho = r {space}", p, lambda r: r, size)
    ck.join("BTZ, the two sheets at the throat", near, rp, far, rp)
    for a, b in ((near, out), (far, far_out)):
        ck.join(f"BTZ, {a.id} in flat space and {b.id} in Minkowski space at r = sqrt(1 + M) l", a, level, b, level)
        # And as the file holds them: the last point of the one is the first of the other, to
        # the rounding of the coarser of the two.
        pa, pb = a.data()["points"][-1], b.data()["points"][0]
        ck.add(f"BTZ, {a.id} and {b.id} as written: one point at r = sqrt(1 + M) l",
               max(abs(pa[0] - pb[0]), abs(pa[1] - pb[1]), abs(pa[2] - pb[2])), 10.0 ** -min(a.decimals, b.decimals))

    surface = Surface([near, out, far, far_out])
    fig = figure_of([surface], {"sheet": "cover"}, size)
    ring_label(fig, [0, 0, 0], rp, 0.0, "$r = r_+$", dx=14)
    ring_label(fig, [0, 0, 0], *near.at(level), "$\\sqrt{2}\\,\\ell$")
    ring_label(fig, [0, 0, 0], *out.at(top), "$3\\ell$")
    fig.legend("fill", "cover", "the exterior $r > r_+$ that $t$ and $r$ cover")
    fig.legend("line", "r", "$r$ constant, at $1.2\\,\\ell$, $2\\ell$ and $3\\ell$")
    fig.legend("line", "r2", "the same radii on the other exterior")
    fig.legend("line", "horizon", "the throat $r = r_+$, where the slice crosses the horizon")
    fig.legend("line", "space", "$r = \\sqrt{2}\\,\\ell$, where $N^2 = 1$: flat space inside, Minkowski space beyond")
    fig.legend("line", "meridian", "$\\phi$ constant, every $15°$")
    return [view("throat", "Through the throat", "$\\ell$", [surface], fig.done(),
                 settings="$\\ell = 1$, the unit of every length, $M = 1$ and $J = 0$, so that $r_+ = \\ell$. Every "
                          "length along the surface beyond $r = \\sqrt{2}\\,\\ell$ is measured with $dX^2 + dY^2 - dZ^2$.")]


MILNE = (0.5, 1.0, 2.0, 3.0)    # the moments of ct drawn, in any unit of length l


def milne(ck, src):
    """The equator of a moment of the comoving hyperbolic chart, g_chichi = c^2t^2 and g_phiphi =
    c^2t^2 sinh^2 chi, at ct = 0.5, 1, 2 and 3 in any unit of length and as a movie with a frame
    every 0.05 of ct. At every chi > 0 its circles grow faster than the distance out to them,
    which is checked, so it is drawn in three dimensional Minkowski space, dX^2 + dY^2 - dZ^2,
    where it climbs at dZ/dchi = sqrt((drho/dchi)^2 - g_chichi) = ct sinh chi from Z = ct on the
    axis: rho = ct sinh chi and Z = ct cosh chi, the sheet Z^2 - rho^2 = c^2t^2 of a hyperboloid.
    With Z = cT and rho = R that is the moment itself on the plane theta = pi/2 of the inertial
    chart, T = t cosh chi and R = ct sinh chi, so every moment is drawn where it lies in spacetime
    and all of them nest inside the light cone Z = rho of the event T = R = 0, drawn in every frame
    as a reference. Each runs out to rho = 4, where chi = arcsinh(4/ct)."""
    top = 4.0
    size = 2 * top
    rings = (0.5, 1.0, 1.5, 2.0, 2.5)

    def moment(t, check):
        sl = Slice(src, "milne", "comoving_hyperbolic", "\\chi", "\\phi", {"t": repr(t), **EQUATOR}, space="minkowski")
        hi = math.asinh(top / t)
        sheet = Piece("sheet", "sheet", sl, 0.0, hi, t, 1,
                      (("axis", "the comoving particle at $\\chi = 0$"),
                       ("edge", "the sheet runs on toward the light cone, to $\\chi \\to \\infty$")),
                      [(x, "r", None) for x in rings if x < hi], size)
        cone = FormPiece("cone", sl, np.linspace(0.0, hi, 81), lambda x: t * np.sinh(x), lambda x: t * np.sinh(x),
                         (("apex", "the event $T = R = 0$, where every comoving particle starts"),
                          ("edge", "the cone runs on")), size)
        where = f"Milne, ct = {t:g}"
        if check:
            flat = Slice(src, "milne", "comoving_hyperbolic", "\\chi", "\\phi", {"t": repr(t), **EQUATOR})
            ck.stops(f"{where}, the equator in flat space", flat, np.linspace(1e-3, 20, 400))
            ck.isometry(f"{where}, the sheet", sheet)
            ck.form(f"{where}, the hyperboloid Z = ct cosh chi", sheet, lambda x: t * np.cosh(x), size)
            ck.radius(f"{where}, rho = ct sinh chi", sheet, lambda x: t * np.sinh(x), size)
        return Surface([sheet, cone], label=f"$ct = {t:g}$", time=t)

    times, keys = movie_values(list(MILNE), 0.05)
    frames = [moment(t, True) for t in times]
    surfaces = [frames[i] for i in keys]
    fig = movie_figure(frames, {"sheet": "cover"}, size, meridians=12)
    fig.legend("fill", "cover", "the equator of the moment, which $\\chi$ and $\\phi$ cover whole")
    fig.legend("line", "r", "$\\chi$ constant, every $1/2$, the circle of comoving particles at that $\\chi$")
    fig.legend("line", "reference", "the light cone of the event $T = R = 0$, which every moment nears as "
                                    "$\\chi \\to \\infty$")
    fig.legend("line", "meridian", "$\\phi$ constant, every $30°$")
    return [view("hyperboloids", "In Minkowski space", "$\\ell$", surfaces, fig.done(),
                 movie=movie(frames, "$ct$", times, turns=False),
                 settings="$\\ell$, any length, the unit of every length and of $ct$, and every length along a "
                          "sheet measured with $dX^2 + dY^2 - dZ^2$.",
                 stops=["At every $\\chi > 0$ the circles grow faster than the distance out to them, "
                        "$g_{\\chi\\chi} < (\\partial_\\chi\\sqrt{g_{\\phi\\phi}})^2$, and no surface of revolution "
                        "in flat space carries the slice; Minkowski space carries it."])]


# At ct = a the pulse bends the surface by the axis on a scale of a/4, against a drawing 20a
# across, so the chords there are a few thousandths of a long and are rounded a decade finer.
ER_DIGITS = 1e-8


def einstein_rosen_waves(ck, src):
    """The pulse of Weber, Wheeler, and Bonnor at C = a = 1, going out from the axis, at the moments
    ct = a to 8a of the cylindrical chart, on the plane z = 0: g_rhorho = e^{2(gamma - psi)} and
    g_phiphi = rho^2 e^{-2psi}, so the circle at rho has radius rho e^{-psi} and the surface climbs
    at dz/drho = e^{-psi} sqrt(e^{2gamma} - (1 - rho d_rho psi)^2). psi and gamma enter as numbers,
    with d_rho psi taken in sympy, and the pulse is checked to make every published Einstein
    component vanish. On the axis gamma = 0 and d_rho psi = 0, so the surface starts level; far
    out psi -> 0 and gamma -> C^2/a^2, a cone of deficit 2 pi (1 - e^{-1}). The pulse is even in t,
    so the moments before it reaches the axis are these in reverse; at t = 0 no surface of
    revolution in flat space carries the moment near the axis, where e^{2gamma} < (1 - rho
    d_rho psi)^2, which is checked, and every moment drawn is checked to embed whole."""
    tt, rr = sp.symbols("t rho", real=True)
    pulse = {k: sp.sympify(v, locals={"t": tt, "rho": rr}) for k, v in nr._er_pulse("t", "rho").items()}
    as_num = {k: (sp.lambdify((tt, rr), e, "numpy"), sp.lambdify((tt, rr), sp.diff(e, rr), "numpy"))
              for k, e in pulse.items()}

    # The published Einstein tensor, mixed, with the pulse, at points all over the plane.
    src.note("einstein_rosen_waves", "cylindrical", ["einstein_tensor"])
    _, entry, R = nr.load("einstein_rosen_waves", "cylindrical")
    points = np.random.default_rng(7).uniform([-8, 0.05], [8, 12], (400, 2))
    worst = 0.0
    for c in entry["einstein_tensor"]["variants"]["ul"]["nonzero"]:
        e = R(c["value"]).subs(R.c, 1)
        for name, text in nr._er_pulse("t", "rho").items():
            e = e.replace(R.parameters[name].func, nr._as_lambda(R, name, text)).doit()
        f = sp.lambdify((R.symbol["t"], R.symbol["\\rho"]), e, "numpy")
        worst = max(worst, float(np.max(np.abs(f(points[:, 0], points[:, 1])))))
    ck.add("Einstein-Rosen: the pulse makes every published G^mu_nu vanish", worst, 1e-10)

    def slice_at(T):
        numeric = {k: (lambda x, f=f, T=T: f(T, x), lambda x, df=df, T=T: df(T, x)) for k, (f, df) in as_num.items()}
        return Slice(src, "einstein_rosen_waves", "cylindrical", "\\rho", "\\phi", {"t": T, "z": 0}, numeric=numeric)

    from scipy.optimize import brentq
    at0 = slice_at(0.0)
    edge0 = brentq(lambda x: float(at0.defect_at(np.array([x]))[0]), 1.0, 3.0, xtol=1e-12)
    ck.stops("Einstein-Rosen, t = 0, out to the edge of the band", at0, np.linspace(1e-3, edge0, 402)[1:-1])
    ck.add("Einstein-Rosen, t = 0: the band ends at rho = 1.72 a", abs(round(edge0, 2) - 1.72), 0.0)

    top = 10.0
    size = 2 * top
    far = as_num["gamma"][0](8.0, 1e6), as_num["psi"][0](8.0, 1e6)
    ck.add("Einstein-Rosen: far out gamma -> C^2/a^2 and psi -> 0", max(abs(far[0] - 1), abs(far[1])), 1e-5)

    def moment(T):
        where = f"Einstein-Rosen, ct = {T:g} a"
        sl = slice_at(T)
        whole = Piece("whole", "sheet", sl, 0.0, top, 0.0, 1,
                      (("axis", "the axis $\\rho = 0$, where the surface starts level"),
                       ("edge", "the surface runs on toward the cone of deficit angle $2\\pi(1 - e^{-C^2/a^2})$, "
                                "to $\\rho \\to \\infty$")),
                      [(r, "r", None) for r in (2.0, 4.0, 6.0, 8.0)] + [(top, "chartedge", None)], size,
                      digits=ER_DIGITS)
        ck.isometry(where, whole)
        ck.radius(f"{where}, rho e^(-psi)", whole, lambda x, T=T: x * np.exp(-as_num["psi"][0](T, x)), size)
        return Surface([whole], label=f"$ct = {T:g}\\,a$", time=T)

    times, keys = movie_values([1.0, 2.0, 4.0, 8.0], 0.1)
    frames = [moment(T) for T in times]
    surfaces = [frames[i] for i in keys]
    fig = movie_figure(frames, {"sheet": "cover"}, size)
    fig.legend("fill", "cover", "the plane $z = 0$ of the moment, which $\\rho$ and $\\phi$ cover")
    fig.legend("line", "r", "$\\rho$ constant, at $2a$, $4a$, $6a$ and $8a$")
    fig.legend("line", "chartedge", "$\\rho = 10\\,a$, where the drawing stops")
    fig.legend("line", "meridian", "$\\phi$ constant, every $30°$")
    return [view("pulse", "The outgoing pulse", "$a$", surfaces, fig.done(),
                 movie=movie(frames, "$ct$", [f.time for f in frames]),
                 settings="$C = a$, with $a$ the unit of every length; each moment is a slice of constant $t$ and $z$.",
                 input=nr.ER_INPUT,
                 stops=[f"At $t = 0$ the circles about the axis out to $\\rho = {edge0:.2f}\\,a$ grow faster than the "
                        "distance out to them, and no surface of revolution in flat space carries that part of the "
                        "moment."])]


MELVIN_ERNST = {"r_s": 1, "B": sp.Rational(1, 2)}   # Ernst's hole at B r_s = 1/2, its widest circle at r = 4 r_s


def melvin(ck, src):
    """Melvin's plane z = 0 at t = 0 in the cylindrical chart at B = 1, and the equator of Ernst's hole
    at t = 0 with r_s = 1 and B = 1/2.

    Melvin's slice has g_rhorho = Lambda^2 and g_phiphi = rho^2/Lambda^2, Lambda = 1 + B^2 rho^2/4, so
    its circles sit at rho/Lambda, which grows to 1/B at the Melvin radius rho = 2/B and falls as 4/(B^2 rho)
    beyond, while the distance out to them grows as rho + B^2 rho^3/12: drho/dx = (1 - B^2 rho^2/4)/Lambda^2
    stays below Lambda, so the slice embeds everywhere, a vase that widens to radius 1/B and narrows into a
    spike without end, drawn out to rho = 4/B. Ernst's equator has g_rr = Lambda^2/(1 - r_s/r) and circles
    of radius r/Lambda, Lambda = 1 + B^2 r^2/4: Flamm's throat at r_s, where 1/g_rr has its root, opening to
    the widest circle at r = 2/B = 4 r_s and closing into the same kind of spike, drawn out to r = 5 r_s on
    both sheets of the slice through the bifurcation sphere."""
    sl = Slice(src, "melvin", "cylindrical", "\\rho", "\\phi", {"t": 0, "z": 0}, {"B": 1})
    top, size = 4.0, 9.0
    vase = Piece("universe", "sheet", sl, 0.0, top, 0.0, 1,
                 (("axis", "the axis $\\rho = 0$"), ("edge", "the surface runs on, narrowing, to $\\rho \\to \\infty$")),
                 [(1.0, "r", None), (2.0, "surface", "$\\rho = 2/B$"), (3.0, "r", None), (top, "r", "$4/B$")], size)
    ck.isometry("Melvin, the plane z = 0", vase)
    ck.radius("Melvin, the plane z = 0, rho = rho/(1 + B^2 rho^2/4)", vase, lambda r: r / (1 + r * r / 4), size)
    surface = Surface([vase])
    fig = figure_of([surface], {"sheet": "cover"}, size)
    ring_label(fig, [0, 0, 0], *vase.at(2.0), "$\\rho = 2/B$", dx=10)
    ring_label(fig, [0, 0, 0], *vase.at(top), "$4/B$")
    fig.legend("fill", "cover", "the plane $z = 0$ at $t = 0$, which $\\rho$ and $\\phi$ cover")
    fig.legend("line", "r", "$\\rho$ constant, at $1$, $3$ and $4/B$")
    fig.legend("line", "surface", "the widest circle, of radius $1/B$, at the Melvin radius $\\rho = 2/B$")
    fig.legend("line", "meridian", "$\\phi$ constant, every $15°$")
    views = [view("universe", "Melvin's universe", "$1/B$", [surface], fig.done(),
                  settings="$B = 1$, so that $1/B$ is the unit of every length.")]

    sl = Slice(src, "melvin", "ernst", "r", "\\phi", {"t": 0, **EQUATOR}, MELVIN_ERNST)
    throat, = [r for r in sl.horizons() if r > 0]
    ck.add("Ernst: the throat is the root r = r_s of 1/g_rr", abs(throat - 1.0), 1e-12)
    top, size = 5.0, 18.0
    radii = (2.0, 3.0)
    near = Piece("exterior", "sheet", sl, throat, top, 0.0, 1,
                 (("throat", "the throat $r = r_s$, the bifurcation sphere, where the other exterior begins"),
                  ("edge", "the surface runs on, narrowing, to $r \\to \\infty$")),
                 [(throat, "horizon", "$r = r_s$")] + [(r, "r", None) for r in radii]
                 + [(4.0, "surface", "$r = 2/B$"), (top, "r", None)], size)
    far = Piece("other_exterior", "sheet2", sl, throat, top, 0.0, -1,
                (("throat", "the throat $r = r_s$"), ("edge", "the surface runs on, narrowing, to $r \\to \\infty$")),
                [(r, "r2", None) for r in radii] + [(4.0, "surface", None), (top, "r2", None)], size)
    ck.join("Ernst, the two sheets at the throat", near, throat, far, throat)
    for p in (near, far):
        ck.isometry(f"Ernst, {p.id}", p)
        ck.radius(f"Ernst, {p.id}, rho = r/(1 + B^2 r^2/4)", p, lambda r: r / (1 + r * r / 16), size)
    surface = Surface([near, far])
    fig = figure_of([surface], {"sheet": "cover"}, size)
    ring_label(fig, [0, 0, 0], *near.at(throat), "$r = r_s$", dx=10)
    ring_label(fig, [0, 0, 0], *near.at(4.0), "$r = 2/B$", dx=10)
    fig.legend("fill", "cover", "the exterior $r > r_s$ that $t$ and $r$ cover")
    fig.legend("line", "r", "$r$ constant, at $2$, $3$ and $5\\,r_s$")
    fig.legend("line", "r2", "the same radii on the other exterior")
    fig.legend("line", "horizon", "the throat $r = r_s$, where the slice crosses the horizon")
    fig.legend("line", "surface", "the widest circle, of radius $1/B$, at $r = 2/B = 4\\,r_s$")
    fig.legend("line", "meridian", "$\\phi$ constant, every $15°$")
    views.append(view("ernst", "Ernst's black hole", "$r_s$", [surface], fig.done(),
                      settings="$r_s = 1$, the unit of every length, and $B = 1/(2r_s)$, so that the widest circle "
                               "of the equator sits at $r = 4\\,r_s$."))
    return views


def levi_civita(ck, src):
    """Levi-Civita's plane z = 0 at t = 0 in Weyl's coordinates at sigma = 1/4 and C = 1.

    The slice has g_rhorho = rho^{4 sigma (2 sigma - 1)} and circles of radius rho^{1 - 2 sigma}/C,
    which grow at (1 - 2 sigma) rho^{-4 sigma^2}/C times the distance out to them: faster than the
    distance inside rho_0 = ((1 - 2 sigma)/C)^{1/(4 sigma^2)}, where no surface of revolution in flat
    space carries the slice, and ever more slowly beyond it. At sigma = 1/4 and C = 1, rho_0 = 1/16,
    the circles have radius sqrt(rho), and dz/drho = sqrt(4 sqrt(rho) - 1)/(2 sqrt(rho)), so the
    surface is z = (4 sqrt(rho) - 1)^{3/2}/6: level on the circle rho_0, of radius 1/4, and a horn
    that steepens without end above it, drawn out to rho = 4. The curvature of the profile diverges
    on the level circle, as (rho - rho_0)^{-1/2}, so its chords there are short and the piece is
    written to the digits a surface in Minkowski space takes."""
    sl = Slice(src, "levi_civita", "weyl", "\\rho", "\\phi", {"t": 0, "z": 0}, nr.LC_WEYL)
    stop, top, size = 1 / 16, 4.0, 4.0
    ck.add("Levi-Civita: g_rhorho = (drho/dx)^2 at rho = 1/16", abs(float(sl.defect_at(np.array([stop]))[0])), 1e-12)
    ck.stops("Levi-Civita, inside rho = 1/16", sl, np.geomspace(1e-6, stop, 202)[:-1])
    horn = Piece("field", "sheet", sl, stop, top, 0.0, 1,
                 (("stops", "the circles grow faster than the distance out to them, and nothing in flat space "
                            "carries the slice on toward the axis"),
                  ("edge", "the surface runs on, steepening, to $\\rho \\to \\infty$")),
                 [(stop, "chartedge", None), (1.0, "r", None), (2.0, "r", None), (3.0, "r", None), (top, "r", None)], size,
                 digits=LORENTZ_DIGITS)
    ck.isometry("Levi-Civita, the plane z = 0", horn)
    ck.radius("Levi-Civita, the plane z = 0, the circle's radius sqrt(rho)", horn, np.sqrt, size)
    ck.form("Levi-Civita, the plane z = 0, z = (4 sqrt(rho) - 1)^(3/2)/6", horn,
            lambda r: np.maximum(4 * np.sqrt(r) - 1, 0) ** 1.5 / 6, size)
    surface = Surface([horn])
    fig = figure_of([surface], {"sheet": "cover"}, size)
    ring_label(fig, [0, 0, 0], *horn.at(stop), "$\\rho = 1/16$", dx=10)
    ring_label(fig, [0, 0, 0], *horn.at(1.0), "$\\rho = 1$", dx=10)
    ring_label(fig, [0, 0, 0], *horn.at(top), "$4$")
    fig.legend("fill", "cover", "the plane $z = 0$ at $t = 0$, which $\\rho$ and $\\phi$ cover")
    fig.legend("line", "r", "$\\rho$ constant, at $1$, $2$, $3$ and $4$")
    fig.legend("line", "chartedge", "$\\rho = 1/16$, where the drawing stops")
    fig.legend("line", "meridian", "$\\phi$ constant, every $15°$")
    return [view("field", "The plane $z = 0$", "$\\rho_0$", [surface], fig.done(),
                 settings="$\\sigma = 1/4$ and $C = 1$, with $\\rho_0$ the unit of length in which the powers of "
                          "$\\rho$ are evaluated.",
                 stops=["Inside $\\rho = 1/16$ the circles grow faster than the distance out to them, "
                        "$g_{\\rho\\rho} < (\\partial_\\rho\\sqrt{g_{\\phi\\phi}})^2$, and no surface of revolution in "
                        "flat space carries the slice on to the axis."])]


def curzon_chazy(ck, src):
    """The plane z = 0 of the Curzon-Chazy particle at t = 0 in Weyl's chart, m = 1.

    On it g_rhorho = e^(2m/rho - m^2/rho^2) and g_phiphi = rho^2 e^(2m/rho), so the circle of Weyl's
    radius rho has radius rho e^(m/rho) on the surface, least, m e, at rho = m and growing without
    bound both as rho grows and as rho -> 0. The slice has a surface of revolution in flat space where
    e^(-m^2/rho^2) >= (1 - m/rho)^2, which is rho >= 0.7226 m, found here: a funnel that narrows to a
    neck at rho = m and widens again below it, until it lies level at rho = 0.7226 m, inside which the
    circles grow faster than the distance in to them, which is checked. Where it lies level the profile
    rises as the 3/2 power of the distance from that circle, so its chords there are short, and it is
    written to nine decimals, as a surface in Minkowski space is. The spherical chart's equator, where r
    is rho, is checked to give the same surface."""
    sl = Slice(src, "curzon_chazy", "weyl", "\\rho", "\\phi", {"t": 0, "z": 0}, {"m": 1})
    stop = float(sp.nsolve(sl.defect, sl.x, 0.72))
    ck.add("Curzon-Chazy: the surface starts where e^(-m^2/rho^2) = (1 - m/rho)^2",
           abs(math.exp(-1 / stop ** 2) - (1 - 1 / stop) ** 2), 1e-12)
    ck.stops("Curzon-Chazy, inside the first surface", sl, np.linspace(0.05, stop, 202)[:-1])
    top = 5.0
    size = 2 * float(sl.rho_at(top))
    plane = Piece("plane", "sheet", sl, stop, top, 0.0, 1,
                  (("stops", "the circles grow faster than the distance in to them, and nothing in flat space carries the slice on"),
                   ("edge", "the surface runs on to $\\rho \\to \\infty$")),
                  [(stop, "chartedge", None), (1.0, "surface", "$\\rho = m$"), (2.0, "r", None), (3.0, "r", None),
                   (4.0, "r", None), (top, "r", "$5\\,m$")], size, digits=LORENTZ_DIGITS)
    ck.isometry("Curzon-Chazy, the plane z = 0", plane)
    ck.radius("Curzon-Chazy, the plane z = 0, rho e^(m/rho)", plane, lambda r: r * np.exp(1 / r), size)
    ck.add("Curzon-Chazy: the narrowest circle is at rho = m, of radius m e",
           abs(float(np.min(plane.rho)) - math.e) + abs(plane.at(1.0)[0] - math.e), 1e-9)
    other = Slice(src, "curzon_chazy", "spherical", "r", "\\phi", {"t": 0, **EQUATOR}, {"m": 1})
    xs = plane.x[1:]
    ck.add("Curzon-Chazy, the spherical chart's equator gives the same surface",
           float(max(np.max(np.abs(other.rho_at(xs) - sl.rho_at(xs))), np.max(np.abs(other.defect_at(xs) - sl.defect_at(xs))))),
           1e-9)
    surface = Surface([plane])
    fig = figure_of([surface], {"sheet": "cover"}, size)
    ring_label(fig, [0, 0, 0], *plane.at(1.0), "$\\rho = m$", dx=10)
    ring_label(fig, [0, 0, 0], *plane.at(stop), f"${stop:.2f}\\,m$", side=-1)
    ring_label(fig, [0, 0, 0], *plane.at(top), "$5\\,m$")
    fig.legend("fill", "cover", "the plane $z = 0$ at $t = 0$, which $\\rho$ and $\\phi$ cover")
    fig.legend("line", "r", "$\\rho$ constant, at $2$, $3$, $4$ and $5\\,m$")
    fig.legend("line", "surface", "the narrowest circle, of radius $m\\,e$, at $\\rho = m$")
    fig.legend("line", "chartedge", f"$\\rho = {stop:.2f}\\,m$, where the drawing stops")
    fig.legend("line", "meridian", "$\\phi$ constant, every $15°$")
    return [view("equator", "The plane $z = 0$", "$m$", [surface], fig.done(),
                 settings="$m = GM/c^2 = 1$, the unit of every length.",
                 stops=[f"Inside $\\rho = {stop:.4f}\\,m$ the circles grow faster than the distance in to them, "
                        "$g_{\\rho\\rho} < (\\partial_\\rho\\sqrt{g_{\\phi\\phi}})^2$, and no surface of revolution in flat "
                        "space carries the slice on."])]


TAUB = (1, sp.Rational(1, 2))   # m and l of Taub's universe, as Taub-NUT's spacetime diagram declares


def taub(T):
    """a_1 = a_2 = a and a_3 = b of Taub's universe at Taub's time T, exact in sympy."""
    m, l = TAUB
    U = (-T ** 2 + 2 * m * T + l ** 2) / (T ** 2 + l ** 2)
    return sp.sqrt(T ** 2 + l ** 2), 2 * l * sp.sqrt(U), U


def mixmaster(ck, src):
    """The great two sphere of the slice at five moments of Abraham Taub's universe, the member of
    the family with a_1 = a_2, at m = 1 and l = 1/2 as Taub-NUT declares, checked to make every
    published Einstein component vanish. Through the identity the great sphere is psi + phi = 0 mod
    2 pi: the hemispheres psi = -phi and psi = 2 pi - phi of theta and phi, on which, with a_1 = a_2
    = a and a_3 = b, the published metric pulls back to a^2 dtheta^2 + (a^2 sin^2 theta + b^2 (1 -
    cos theta)^2) dphi^2, a surface of revolution about the axis through the identity, meeting its
    twin at the equator theta = pi, a fibre of psi of circumference 4 pi b. Near the pole g_thth -
    (drho/dtheta)^2 = (1 - 3b^2/4a^2) a^2 theta^2, so where b > 2a/sqrt(3) the drawing is the band
    about the equator out to where that vanishes."""
    m, l = TAUB
    T = sp.Symbol("T", real=True)
    a, b, U = taub(T)
    _, entry, R = nr.load("mixmaster", "euler_angles")
    src.note("mixmaster", "euler_angles", ["einstein_tensor"])
    t = R.symbol["t"]
    D = lambda e: sp.sqrt(U) * sp.diff(e, T)  # noqa: E731, d/d(c tau) along Taub's time
    values = {"a_1": a, "a_2": a, "a_3": b}
    zero = True
    for c in entry["einstein_tensor"]["variants"]["ul"]["nonzero"]:
        e = R(c["value"]).subs(R.c, 1)
        for name, v in values.items():
            fn = R.parameters[name]
            e = e.subs(sp.Derivative(fn, (t, 2)), D(D(v))).subs(sp.Derivative(fn, t), D(v)).subs(fn, v)
        zero = zero and sp.simplify(e) == 0
    ck.exact("Mixmaster: Taub's universe makes every published Einstein component vanish", zero)
    lo = float(m - sp.sqrt(m * m + l * l))
    Uf = sp.lambdify(T, U, "numpy")
    from scipy.optimize import brentq
    round_T = brentq(lambda s: float(sp.N((a - b).subs(T, s))), 0.9, 0.95, xtol=1e-15)
    moments = [sp.Rational(-1, 10), sp.Rational(1, 5), sp.Rational(repr(round_T)), sp.Rational(3, 2), sp.Integer(2)]
    size = 8.0

    def proper(Tm):
        return integrate(lambda s: 1 / math.sqrt(Uf(s)), lo, float(Tm))

    def moment(Tm):
        A, B = a.subs(T, Tm), b.subs(T, Tm)
        af, bf = float(A), float(B)
        c_tau = proper(Tm)
        halves = [Slice(src, "mixmaster", "euler_angles", "\\theta", "\\phi", {"t": 0},
                        functions={"a_1": str(A), "a_2": str(A), "a_3": str(B)}, swept={"psi": sweep}, rewrite=half_angles)
                  for sweep in ("-phi", "2*pi - phi")]
        where = f"Mixmaster, Taub's T = {float(Tm):.6f}"
        if bf > 2 * af / math.sqrt(3):
            begin = brentq(lambda th: float(halves[0].defect_at(th)), 1e-3, math.pi - 1e-9, xtol=1e-14)
            ck.stops(f"{where}, about the poles", halves[0], np.linspace(0, begin, 202)[1:-1])
            ends = ("stops", "the circles about the pole grow faster than the distance out to them, and no surface of "
                             "revolution in flat space carries the sphere there")
        else:
            begin = 0.0
            ends = ("axis", "the pole, the identity" )
        H = halves[0].rise(begin, math.pi)
        rings = [(th, "r", None) for th in (math.pi / 3, 2 * math.pi / 3) if th > begin]
        near = Piece("near", "sheet", halves[0], begin, math.pi, -H, 1,
                     (ends, ("join", "the equator $\\theta = \\pi$, a fibre of $\\psi$")),
                     rings + [(math.pi, "chartedge", None)], size)
        far = Piece("far", "sheet", halves[1], begin, math.pi, H, -1,
                    ((ends[0], ends[1].replace("the identity", "its antipode")), ("join", "the equator")), rings, size)
        for p in (near, far):
            ck.isometry(f"{where}, the {p.id} hemisphere", p)
        ck.join(f"{where}, the hemispheres meet at the equator", near, math.pi, far, math.pi)
        ck.add(f"{where}: the equator has radius 2b", abs(near.at(math.pi)[0] - 2 * bf), 1e-9)
        if begin == 0.0:
            ck.add(f"{where}: pole to pole is pi a long", abs(halves[0].proper(0, math.pi) - math.pi * af) / af, 1e-10)
        if abs(af - bf) < 1e-12:
            ck.form(f"{where}, the round sphere z = -2a cos(theta/2)", near, lambda th: -2 * af * np.cos(th / 2), size)
            ck.radius(f"{where}, the round sphere rho = 2a sin(theta/2)", near, lambda th: 2 * af * np.sin(th / 2), size)
        return Surface([near, far], label=f"$c\\tau = {c_tau:.2f}\\,m$", time=c_tau)

    # The movie runs through the moments at a steady proper time, a frame about every 0.05 of c tau,
    # each at the T whose proper time it is.
    taus, keys = movie_values([proper(Tm) for Tm in moments], 0.05)
    Ts = [moments[keys.index(i)] if i in keys else
          sp.Rational(repr(brentq(lambda x, c=c: proper(x) - c, float(moments[0]), float(moments[-1]), xtol=1e-14)))
          for i, c in enumerate(taus)]
    frames = [moment(Tm) for Tm in Ts]
    surfaces = [frames[i] for i in keys]
    fig = movie_figure(frames, {"sheet": "cover"}, size, meridians=12)
    fig.legend("fill", "cover", "the great sphere through the identity, which the Euler angles cover but for its poles and "
                                "equator")
    fig.legend("line", "r", "$\\theta$ constant, at $\\pi/3$ and $2\\pi/3$ on each hemisphere")
    fig.legend("line", "chartedge", "the equator $\\theta = \\pi$, a fibre of $\\psi$ of circumference $4\\pi a_3$")
    fig.legend("line", "meridian", "$\\phi$ constant, every $30°$, running on through the equator into the other hemisphere")
    return [view("sphere", "The great sphere", "$m$", surfaces, fig.done(),
                 movie=movie(frames, "$c\\tau$", [f.time for f in frames]),
                 settings="$m = 1$, the unit of every length, and $l = m/2$; each moment is labelled by the proper time "
                          "$c\\tau$ from Taub's first horizon.",
                 input="Taub's universe, $a_1 = a_2 = \\sqrt{T^2 + l^2}$ and $a_3 = 2l\\sqrt{U}$ with $U = (-T^2 + 2mT + "
                       "l^2)/(T^2 + l^2)$ and $c\\,d\\tau = dT/\\sqrt{U}$, checked to make every Einstein component of "
                       "this spacetime vanish.",
                 stops=["When the three scale factors differ, as in every other Mixmaster universe, the great sphere's "
                        "metric depends on $\\phi$ as well as $\\theta$, and no surface of revolution carries it.",
                        "While $a_3$ is more than $2/\\sqrt{3}$ times $a_1$, the curvature about the poles "
                        "is negative, and only the band about the equator has a surface of revolution in flat space."])]


# ---------------------------------------------------------------- the spacetimes with nothing to draw

def c_metric(ck, src):
    """The C-metric at alpha m = 1/6 and C = 1/(1 + 2 alpha m) = 3/4, as its other diagrams
    declare, m = 1, in its spherical chart. Two views.

    The equator of the moment t = 0: there 1 + alpha r cos(theta) and 1 + 2 alpha m cos(theta)
    are both 1, so the slice is dr^2/Q + C^2 r^2 dphi^2 with Q = (1 - alpha^2 r^2)(1 - 2m/r),
    and rho = C r. Q < 1 between the horizons, so the profile climbs from the black hole
    horizon r = 2m, where it is vertical, to the acceleration horizon r = 1/alpha, where it is
    vertical again at its widest circle. Through the acceleration horizon's bifurcation the
    moment runs on into the second black hole's exterior, the same surface turned over, and the
    two meet there with one tangent.

    The black hole horizon itself, r = 2m at one moment, its bifurcation sphere: there
    1 + alpha r cos(theta) = 1 + 2 alpha m cos(theta) = P, so the metric is
    4m^2 dtheta^2/P^3 + 4m^2 C^2 sin^2(theta) dphi^2/P and rho = 2mC sin(theta)/sqrt(P). The
    pole theta = 0 is smooth, since drho/ds = C(1 + 2 alpha m) = 1 there, and theta = pi is the
    apex of a cone, drho/ds = C(1 - 2 alpha m) = 1/2, where the cosmic string meets the horizon,
    the deficit angle 8 pi alpha m/(1 + 2 alpha m) being 2 pi times what that rate falls short
    of 1. Its area is checked against Griffiths, Krtous and Podolsky's eq. (21),
    16 pi C m^2/(1 - 4 alpha^2 m^2) = 13.5 pi."""
    params = {"m": 1, "alpha": "1/6", "C": "3/4"}
    C, alpha = 0.75, 1 / 6
    sl = Slice(src, "c_metric", "spherical", "r", "\\phi", {"t": 0, **EQUATOR}, params)
    found = sorted(h for h in sl.horizons() if h > 0)
    ck.add("C-metric: the horizons on the equator are at r = 2m and 1/alpha = 6m",
           abs(found[0] - 2) + abs(found[-1] - 6) + abs(len(found) - 2), 1e-12)
    inner, outer = 2.0, 6.0
    size = 2 * C * outer
    near = Piece("near", "sheet", sl, inner, outer, 0.0, 1,
                 (("throat", "the black hole horizon $r = 2m$, its bifurcation sphere"),
                  ("join", "the acceleration horizon $r = 1/\\alpha$, the widest circle, where the moment runs on "
                           "into the second black hole's exterior")),
                 [(inner, "horizon", "$r = 2m$")] + [(r, "r", None) for r in (3.0, 4.0, 5.0)]
                 + [(outer, "horizon", "$r = 1/\\alpha$")], size)
    top = float(near.z[-1])
    far = Piece("far", "sheet2", sl, inner, outer, 2 * top, -1,
                (("throat", "the second black hole's horizon, $r = 2m$"), ("join", "the acceleration horizon")),
                [(inner, "horizon", None)] + [(r, "r2", None) for r in (3.0, 4.0, 5.0)], size)
    for p in (near, far):
        ck.isometry(f"C-metric equator, {p.id}", p)
        ck.radius(f"C-metric equator, {p.id}, rho = C r", p, lambda r: C * r, size)
    ck.join("C-metric equator, the two exteriors at the acceleration horizon", near, outer, far, outer)
    surface = Surface([near, far])
    fig = figure_of([surface], {"sheet": "cover"}, size)
    ring_label(fig, [0, 0, 0], *near.at(inner), "$r = 2m$", dx=10)
    ring_label(fig, [0, 0, 0], *near.at(outer), "$r = 1/\\alpha = 6m$")
    fig.legend("fill", "cover", "the exterior $2m < r < 1/\\alpha$ that $t$ and $r$ cover")
    fig.legend("line", "r", "$r$ constant, at $3m$, $4m$ and $5m$")
    fig.legend("line", "r2", "the same radii in the second black hole's exterior")
    fig.legend("line", "horizon", "the horizons: the throats $r = 2m$ and the widest circle $r = 1/\\alpha$")
    fig.legend("line", "meridian", "$\\phi$ constant, every $15°$")
    settings = ("$m = 1$, the unit of every length, $\\alpha = 1/(6m)$, and $C = 1/(1 + 2\\alpha m) = 3/4$, "
                "which leaves the half axis $\\theta = 0$ regular.")
    views = [view("equator", "The equator", "$m$", [surface], fig.done(), settings=settings)]

    hz = Slice(src, "c_metric", "spherical", "\\theta", "\\phi", {"t": 0, "r": 2}, params)
    hsize = 4 * C
    horizon = Piece("horizon", "sheet", hz, 0.0, math.pi, 0.0, -1,
                    (("axis", "the pole $\\theta = 0$, on the half axis between the black holes, where the "
                              "horizon is smooth"),
                     ("apex", "the pole $\\theta = \\pi$, where the cosmic string meets the horizon")),
                    [(math.pi / 4, "r", None), (math.pi / 2, "r", "$\\theta = \\pi/2$"), (3 * math.pi / 4, "r", None)],
                    hsize)
    ck.isometry("C-metric horizon", horizon)
    ck.radius("C-metric horizon, rho = 2mC sin(theta)/sqrt(1 + 2 alpha m cos(theta))", horizon,
              lambda th: 2 * C * np.sin(th) / np.sqrt(1 + np.cos(th) / 3), hsize)
    for pole, want in ((0.0, 1.0), (math.pi, 0.5)):
        h = 1e-6 if pole == 0 else -1e-6
        rate = abs(float(hz.rho_at(pole + h)) - float(hz.rho_at(pole))) / hz.proper(min(pole, pole + h), max(pole, pole + h))
        ck.add(f"C-metric horizon: drho/ds at theta = {pole:.4f} is C(1 {'+' if pole == 0 else '-'} 2 alpha m)",
               abs(rate - want), 1e-5)
    pts = np.array(horizon.data()["points"])
    ds = np.hypot(np.diff(pts[:, 1]), np.diff(pts[:, 2]))
    area = float(np.sum(2 * math.pi * 0.5 * (pts[1:, 1] + pts[:-1, 1]) * ds))
    ck.add("C-metric horizon: the area is 16 pi C m^2/(1 - 4 alpha^2 m^2), Griffiths, Krtous and Podolsky's (21)",
           abs(area - 16 * math.pi * C / (1 - 4 * alpha ** 2)) / (13.5 * math.pi), 1e-4)
    hsurface = Surface([horizon])
    fig = figure_of([hsurface], {"sheet": "cover"}, hsize)
    ring_label(fig, [0, 0, 0], *horizon.at(math.pi / 2), "$\\theta = \\pi/2$")
    fig.legend("fill", "cover", "the black hole horizon $r = 2m$ at one moment")
    fig.legend("line", "r", "$\\theta$ constant, at $\\pi/4$, $\\pi/2$ and $3\\pi/4$")
    fig.legend("line", "meridian", "$\\phi$ constant, every $15°$")
    views.append(view("horizon", "The black hole horizon", "$m$", [hsurface], fig.done(), settings=settings))
    return views


NARIAI_MOMENTS = (-1.5, -0.75, 0.0, 0.75, 1.5)     # ct in units of 1/sqrt(Lambda)


def nariai(ck, src):
    """dS2 x S2 with both radii a = 1/sqrt(Lambda) = 1. A moment ct of the global chart is a
    circle of chi of radius a cosh(ct/a) times the round sphere of radius a, and its surface of
    chi and a great circle of the sphere, theta at phi = 0 down from the north pole and back up
    at phi = pi, has g_thetatheta = a^2 and g_chichi = a^2 cosh^2(ct/a): a cylinder of radius
    a cosh(ct/a) turned about by chi, z = a theta on the first half and a(2 pi - theta) on the
    second, its top and bottom edges one circle, the north pole times the circle of chi. Drawn at
    five moments from ct = -1.5 a to 1.5 a, the circle contracting to a at ct = 0 and expanding
    after, played as a movie with a frame every 0.1 a of ct. The sphere of theta and phi at one
    event is the second view, the same at every event."""
    size = 2 * math.pi
    top = 2 * math.pi

    def moment(t):
        near_sl = Slice(src, "nariai", "global", "\\theta", "\\chi", {"t": repr(t), "phi": "0"}, {"Lambda": 1})
        far_sl = Slice(src, "nariai", "global", "\\theta", "\\chi", {"t": repr(t), "phi": "pi"}, {"Lambda": 1})
        near = Piece("near", "sheet", near_sl, 0.0, math.pi, 0.0, 1,
                     (("edge", "the north pole, $\\theta = 0$, one circle with the top edge"),
                      ("join", "the south pole, $\\theta = \\pi$, where the great circle turns to $\\phi = \\pi$")),
                     [(math.pi / 4, "r", None), (math.pi / 2, "r", None), (3 * math.pi / 4, "r", None)], size)
        far = Piece("far", "sheet2", far_sl, 0.0, math.pi, top, -1,
                    (("edge", "the north pole, $\\theta = 0$, one circle with the bottom edge"),
                     ("join", "the south pole")),
                    [(math.pi / 4, "r2", None), (math.pi / 2, "r2", None), (3 * math.pi / 4, "r2", None)], size)
        where = f"Nariai, ct = {t:+.2f} a"
        radius = math.cosh(t)
        for p, z_of in ((near, lambda th: th), (far, lambda th: top - th)):
            ck.isometry(f"{where}, the half great circle at phi = {'0' if p is near else 'pi'}", p)
            ck.radius(f"{where}, {p.id}: the cylinder rho = a cosh(ct/a)", p, lambda th, k=radius: np.full_like(th, k), size)
            ck.form(f"{where}, {p.id}: z = a theta along the great circle", p, z_of, size)
        ck.join(f"{where}, the halves meet at the south pole", near, math.pi, far, math.pi)
        return Surface([near, far], label=f"$ct = {t:.2f}$", time=t)

    values, keys = movie_values(list(NARIAI_MOMENTS), 0.1)
    frames = [moment(t) for t in values]
    surfaces = [frames[i] for i in keys]
    fig = movie_figure(frames, {"sheet": "cover"}, size, meridians=12)
    fig.legend("fill", "cover", "the half great circle $\\phi = 0$ times the circle of $\\chi$")
    fig.legend("line", "r", "$\\theta$ constant at $\\phi = 0$, at $\\pi/4$, $\\pi/2$ and $3\\pi/4$")
    fig.legend("line", "r2", "the same at $\\phi = \\pi$")
    fig.legend("line", "meridian", "$\\chi$ constant, every $30°$")
    settings = "$\\Lambda = 1$, so that $1/\\sqrt{\\Lambda}$ is the unit of every length."
    views = [view("universe", "The circle and a great circle", "$1/\\sqrt{\\Lambda}$", surfaces, fig.done(),
                  movie=movie(frames, "$ct$", [f.time for f in frames]), settings=settings)]

    sphere_slice = Slice(src, "nariai", "static", "\\theta", "\\phi", {"t": 0, "r": 0}, {"Lambda": 1})
    ball = Piece("sphere", "sheet", sphere_slice, 0.0, math.pi, 0.0, 1,
                 (("axis", "the pole $\\theta = 0$"), ("axis", "the pole $\\theta = \\pi$")),
                 [(math.pi / 4, "r", None), (math.pi / 2, "r", None), (3 * math.pi / 4, "r", None)], 2.0)
    ck.isometry("Nariai, the sphere", ball)
    ck.form("Nariai, the sphere z = a(1 - cos theta)", ball, lambda c: 1 - np.cos(c), 2.0)
    ck.radius("Nariai, the sphere rho = a sin theta", ball, np.sin, 2.0)
    sphere = Surface([ball])
    fig = figure_of([sphere], {"sheet": "cover"}, 2.0)
    ring_label(fig, [0, 0, 0], *ball.at(math.pi / 2), "$\\theta = \\pi/2$")
    fig.legend("fill", "cover", "the sphere, which $\\theta$ and $\\phi$ cover but for its poles")
    fig.legend("line", "r", "$\\theta$ constant, at $\\pi/4$, $\\pi/2$ and $3\\pi/4$")
    fig.legend("line", "meridian", "$\\phi$ constant, every $15°$")
    views.append(view("sphere", "The sphere", "$1/\\sqrt{\\Lambda}$", [sphere], fig.done(), settings=settings))
    return views


# ct in units of 1/k. Beyond |kct| = 1.1 the cones tip so near the light cone, |drho/dz| = cosh(kct),
# that a chord turned 0.02 round the axis beside the apex misses the surface by more than ACROSS.
WALL_MOMENTS = (-1.0, -0.5, 0.0, 0.5, 1.0)


def domain_wall(ck, src):
    """The equator of a moment of the global chart, g_zz = 1 and g_phiphi = (1 - k|z|)^2
    cosh^2(kct)/k^2, at k = 1, from the centre of one side, z = -1/k, through the wall at z = 0 to
    the centre of the other, z = 1/k. Its circles grow at |drho/dz| = cosh(kct), faster than the
    distance out to them at every moment but ct = 0, which is checked, so it is drawn in three
    dimensional Minkowski space, dX^2 + dY^2 - dZ^2, where it climbs at dZ/dz = sqrt(cosh^2(kct) - 1)
    = |sinh(kct)|: rho = (1 - |z|) cosh(kct) and Z = z sinh(kct), two cones joined rim to rim at the
    wall. The side z < 0 is then its own moment where it lies in its own inertial chart, Z = cT and
    rho = R, the cone cT = R tanh(kct) from the centre at T = 0 out to the wall, moved so that the
    wall stands at Z = 0; the side z > 0, the same cone of the second copy, is its mirror image in
    the wall, which is an isometry of Minkowski space. At ct = 0 both cones are the flat disc of
    radius 1/k, and the moment is the disc taken twice, joined at its rim. Drawn at five moments
    from ct = -1/k to 1/k and played as a movie with a frame every 0.05/k of ct."""
    size = 4.0
    rings = (0.25, 0.5, 0.75)

    def moment(t):
        sl = Slice(src, "domain_wall", "global", "z", "\\phi", {"t": repr(t), **EQUATOR}, {"k": 1}, space="minkowski")
        sense = -1 if t < 0 else 1
        where = f"domain wall, ct = {t:+.2f}/k"
        near = Piece("near", "sheet", sl, -1.0, 0.0, 0.0, sense,
                     (("apex", "the centre of the side $z < 0$ at $T = 0$, where the horizon meets it"),
                      ("crease", "the wall, where the surface folds")),
                     [(-x, "r", None) for x in rings] + [(0.0, "surface", None)], size)
        far = Piece("far", "sheet2", sl, 0.0, 1.0, float(near.z[-1]), sense,
                    (("crease", "the wall"), ("apex", "the centre of the side $z > 0$ at $T = 0$")),
                    [(x, "r2", None) for x in rings], size)
        shift = -float(near.z[-1])
        near.z, far.z = near.z + shift, far.z + shift
        ck.add(f"{where}, the two sides meet at the wall: one point",
               float(np.max(np.abs(np.array(near.at(0.0)) - far.at(0.0)))), JOIN)
        flat = Slice(src, "domain_wall", "global", "z", "\\phi", {"t": repr(t), **EQUATOR}, {"k": 1})
        xs = np.concatenate([np.linspace(-0.999, -1e-3, 200), np.linspace(1e-3, 0.999, 200)])
        if t == 0:
            ck.plane(f"{where}, the equator is a plane", flat, xs)
        else:
            ck.stops(f"{where}, the equator in flat space", flat, xs)
        for piece in (near, far):
            ck.isometry(f"{where}, {piece.id}", piece)
            ck.radius(f"{where}, {piece.id}: rho = (1 - k|z|) cosh(kct)/k", piece,
                      lambda z: (1 - np.abs(z)) * math.cosh(t), size)
            ck.form(f"{where}, {piece.id}: Z = z sinh(kct)", piece, lambda z: z * math.sinh(t), size)
        return Surface([near, far], label=f"$kct = {t:.2f}$", time=t)

    values, keys = movie_values(list(WALL_MOMENTS), 0.05)
    frames = [moment(t) for t in values]
    surfaces = [frames[i] for i in keys]
    fig = movie_figure(frames, {"sheet": "cover"}, size, meridians=12)
    fig.legend("fill", "cover", "the side $z < 0$ of the equator of the moment, which $z$ and $\\phi$ cover")
    fig.legend("line", "r", "$z$ constant on the side $z < 0$, at $-1/4k$, $-1/2k$ and $-3/4k$")
    fig.legend("line", "r2", "the same on the side $z > 0$, at $1/4k$, $1/2k$ and $3/4k$")
    fig.legend("line", "surface", "the wall, $z = 0$, where the surface folds")
    fig.legend("line", "meridian", "$\\phi$ constant, every $30°$")
    return [view("moments", "The closed moments", "$1/k$", surfaces, fig.done(),
                 movie=movie(frames, "$kct$", values),
                 settings="$k = 1$, so that $1/k$, the radius of the wall when it stops, is the unit of every length "
                          "and of $ct$, and every length along the surface measured with $dX^2 + dY^2 - dZ^2$.",
                 stops=["At every moment but $ct = 0$ the circles grow faster than the distance out to them, "
                        "$g_{zz} < (\\partial_z\\sqrt{g_{\\phi\\phi}})^2$, and no surface of revolution in flat "
                        "space carries the slice; Minkowski space carries it."])]


KS_DUST_MOMENTS = (0.0, 0.3, 0.6, 0.9, 1.2)         # eta, the dust chart's parametric time, from the widest moment on
KS_VACUUM_MOMENTS = (0.9, 0.7, 0.5, 0.3, 0.1)       # T in units of r_s, toward the singularity


def kantowski_sachs(ck, src):
    """The equator of a moment of a Kantowski-Sachs universe has the metric a^2 dr^2 + b^2 dphi^2
    with a and b constant on it: a flat cylinder of radius b about the axis r, on which the
    stretch |r| <= 1 is 2a long, rho = b and z = a r.

    The first view is the dust universe symmetric in time, kappa = 0 of the dust chart at
    b_0 = 1, where a = 1 + eta tan(eta) and b = cos^2(eta): widest and shortest at eta = 0, and
    narrowing while it lengthens toward either singularity, eta -> +-pi/2. Drawn at five moments
    of its collapse, from eta = 0 to 1.2, and played as a movie with a frame every 0.03 of eta;
    the expansion before eta = 0 is the same movie run backward.

    The second is the vacuum member, the inside of Schwarzschild's horizon at r_s = 1, where
    a = sqrt(1/T - 1) and b = T, at five moments from T = 0.9 down to 0.1, set side by side in
    the order of the proper time of an observer at fixed r since the horizon,
    c tau = eta + sin(eta) cos(eta) with T = cos^2(eta), since the future lies toward smaller T."""
    reach = 1.0
    marks = [(-0.5, "r", None), (0.0, "r", None), (0.5, "r", None)]
    ends = (("edge", "the cylinder runs on for ever toward $r \\to -\\infty$"),
            ("edge", "the cylinder runs on for ever toward $r \\to \\infty$"))

    def tube(sl, a, b, where, size):
        piece = Piece("tube", "sheet", sl, -reach, reach, -a * reach, 1, ends, marks, size)
        ck.isometry(where, piece)
        ck.radius(f"{where}, rho = b", piece, lambda r: np.full_like(r, b), size)
        ck.form(f"{where}, z = a r", piece, lambda r: a * r, size)
        return piece

    size = 2 * (1 + 1.2 * math.tan(1.2)) * reach

    def moment(eta):
        sl = Slice(src, "kantowski_sachs", "dust", "r", "\\phi", {"eta": repr(eta), **EQUATOR}, {"b_0": 1, "kappa": 0})
        piece = tube(sl, 1 + eta * math.tan(eta), math.cos(eta) ** 2, f"Kantowski-Sachs dust, eta = {eta:+.2f}", size)
        return Surface([piece], label=f"$\\eta = {eta:.2f}$", time=eta)

    values, keys = movie_values(list(KS_DUST_MOMENTS), 0.03)
    frames = [moment(eta) for eta in values]
    surfaces = [frames[i] for i in keys]
    fig = movie_figure(frames, {"sheet": "cover"}, size, meridians=12)
    fig.legend("fill", "cover", "the stretch $|r| \\le b_0$ of the equator of a moment, which $r$ and $\\phi$ cover")
    fig.legend("line", "r", "$r$ constant, at $-b_0/2$, $0$ and $b_0/2$, each a circle of circumference $2\\pi b$")
    fig.legend("line", "meridian", "$\\phi$ constant, every $30°$")
    views = [view("dust", "Dust", "$b_0$", surfaces, fig.done(), movie=movie(frames, "$\\eta$", values),
                  settings="$b_0 = 1$, the unit of every length, and $\\kappa = 0$, with the moments named by the "
                           "dust chart's time $\\eta$.")]

    vacuum = []
    for T in KS_VACUUM_MOMENTS:
        sl = Slice(src, "kantowski_sachs", "schwarzschild_interior", "r", "\\phi", {"T": repr(T), **EQUATOR}, {"r_s": 1})
        piece = tube(sl, math.sqrt(1 / T - 1), T, f"Kantowski-Sachs vacuum, T = {T:g} r_s", 6.0)
        eta = math.acos(math.sqrt(T))
        vacuum.append(Surface([piece], label=f"$T = {T:g}$", time=eta + math.sin(eta) * math.cos(eta)))
    fig = sequence_figure(vacuum, {"sheet": "cover"}, 6.0, columns=5)
    fig.legend("fill", "cover", "the stretch $|r| \\le r_s$ of the equator of a moment, which $r$ and $\\phi$ cover")
    fig.legend("line", "r", "$r$ constant, at $-r_s/2$, $0$ and $r_s/2$, each a circle of circumference $2\\pi T$")
    fig.legend("line", "meridian", "$\\phi$ constant, every $30°$")
    views.append(view("vacuum", "Vacuum", "$r_s$", vacuum, fig.done(), system="schwarzschild_interior",
                      settings="$r_s = 1$, the unit of every length, with the moments in the order of the proper time "
                               "$c\\tau = r_s\\left(\\eta + \\sin\\eta\\cos\\eta\\right)$ of an observer at fixed $r$ since "
                               "the horizon, where $T = r_s\\cos^2\\eta$."))
    return views


def flat_slices(ck, src, metric_id, system_id, time="t"):
    """Check that every slice of constant `time` of a coordinate system is flat: its spatial
    metric has no cross term and no component that depends on a spatial coordinate, so at
    each moment it is Euclidean space with its axes scaled, and its equator is a plane."""
    _, entry, R = nr.load(metric_id, system_id)
    src.note(metric_id, system_id, FIELDS)
    g = nr.published_matrix(R, entry, "metric_components")
    space = [i for i, c in enumerate(entry["coords"]) if c != time]
    where = {R.symbol[entry["coords"][i]] for i in space}
    ok = (all(sp.simplify(g[i, j]) == 0 for i in space for j in space if i != j)
          and all(g[i, i] != 0 and not (sp.sympify(g[i, i]).free_symbols & where) for i in space))
    ck.exact(f"{metric_id}/{system_id}: every slice of constant {time} is flat", ok)


# The spacetimes with no surface to draw that say why, each function checking what it states
# from the published metric and returning the sentences.
STATED = {}


# ---------------------------------------------------------------- the tables

DRAWN = {
    "aichelburg_sexl": aichelburg_sexl,
    "schwarzschild": schwarzschild,
    "misner": misner,
    "interior_schwarzschild": interior_schwarzschild,
    "tov": tov,
    "morris_thorne": morris_thorne,
    "ellis_bronnikov": ellis_bronnikov,
    "rn_metric": rn_metric,
    "de_sitter": de_sitter,
    "einstein_static": einstein_static,
    "schwarzschild_de_sitter": schwarzschild_de_sitter,
    "vaidya": vaidya,
    "oppenheimer_snyder": oppenheimer_snyder,
    "tolman_bondi": tolman_bondi,
    "bertotti_robinson": bertotti_robinson,
    "stockum_dust": stockum_dust,
    "taub_nut": taub_nut,
    "godel": godel,
    "kerr": kerr,
    "kerr_newman": kerr_newman,
    "cosmic_string": cosmic_string,
    "frw": frw,
    "milne": milne,
    "domain_wall": domain_wall,
    "minkowski": minkowski,
    "anti_de_sitter": anti_de_sitter,
    "malament_hogarth": malament_hogarth,
    "mixmaster": mixmaster,
    "kasner": kasner,
    "bianchi": bianchi,
    "pp_wave": pp_wave, "khan_penrose": khan_penrose,
    "krasnikov": krasnikov,
    "alcubierre": alcubierre,
    "natario": natario,
    "btz": btz,
    "c_metric": c_metric,
    "einstein_rosen_waves": einstein_rosen_waves,
    "nariai": nariai,
    "global_monopole": global_monopole,
    "majumdar_papapetrou": majumdar_papapetrou,
    "melvin": melvin,
    "thin_shell_wormhole": thin_shell_wormhole,
    "levi_civita": levi_civita,
    "kantowski_sachs": kantowski_sachs,
    "curzon_chazy": curzon_chazy,
}

# The spacetimes with no embedding diagram, for which nothing is written. Every slice of constant
# t in Lentz's class is flat, h_ij = delta_ij, for every potential phi, and his soliton exists only
# as a numerical integral over rhomboid sources whose sizes, charges and places his paper does not
# give, so the plane of its path could carry nothing of the soliton.
NOT_DRAWN = {"lentz"}

CAPTIONS = {
    ("btz", "throat"): [
        "The moment $t = 0$ of the hole without rotation ($M = 1$, $J = 0$) through both of its exteriors, "
        "joined at the bifurcation circle $r = r_+$, in flat space out to the circle $r = \\sqrt{1 + M}\\,\\ell$ "
        "and in three dimensional Minkowski space ($dX^2 + dY^2 - dZ^2$) beyond it, where the slice's radial "
        "stretching $g_{rr} = 1/N^2$ drops below that of flat space, every distance along the surface the "
        "metric distance.",
        "On the slice the metric is $dr^2/N^2 + r^2d\\phi^2$ with $N^2 = r^2/\\ell^2 - M$. Inside the circle the "
        "surface climbs at $dz/dr = \\sqrt{1/N^2 - 1}$, from vertical at the throat, as Schwarzschild's does at "
        "$r_s$, to level where $N^2 = 1$. Beyond it the circles grow faster than the distance out to them, as on "
        "the static slice of anti-de Sitter space, and the surface climbs at $dZ/dr = \\sqrt{1 - 1/N^2}$ from level "
        "toward a light cone of Minkowski space, as the hyperboloid of anti-de Sitter space does. Both parts lie "
        "level at the circle, so they meet there with one tangent plane.",
    ],
    ("c_metric", "equator"): [
        "The equatorial plane ($\\theta = \\pi/2$) of the C-metric at one moment of $t$ ($\\alpha m = 1/6$, "
        "$C = 3/4$), drawn as a surface in flat space with every distance along it the metric distance. On the "
        "equator $1 + \\alpha r\\cos\\theta$ and $1 + 2\\alpha m\\cos\\theta$ are both 1, so the slice is "
        "$dr^2/Q + C^2r^2\\,d\\phi^2$ with $Q = (1 - \\alpha^2r^2)(1 - 2m/r)$, and the circle of coordinate $r$ "
        "has radius $Cr$. Between the horizons $Q < 1$, so the surface climbs from the black hole horizon at "
        "$r = 2m$, where it stands vertical, out to the acceleration horizon at $r = 1/\\alpha = 6m$, where it "
        "stands vertical again at its widest circle.",
        "Through the bifurcation of the acceleration horizon the moment runs on into the exterior of the second "
        "black hole, the same surface turned over, and closes in on that black hole's horizon. Each throat is "
        "the bifurcation sphere of a black hole horizon, through which the extension of Jerry Griffiths, Pavel "
        "Krtouš, and Jiří Podolský continues into a further exterior.",
    ],
    ("c_metric", "horizon"): [
        "The black hole horizon $r = 2m$ at one moment of $t$, its bifurcation sphere ($\\alpha m = 1/6$, "
        "$C = 3/4$), drawn as a surface in flat space with every distance along it the metric distance. There "
        "$1 + \\alpha r\\cos\\theta$ equals $P = 1 + 2\\alpha m\\cos\\theta$, so the metric is "
        "$4m^2\\,d\\theta^2/P^3 + 4m^2C^2\\sin^2\\theta\\,d\\phi^2/P$, and the circle of latitude $\\theta$ has "
        "radius $2mC\\sin\\theta/\\sqrt{P}$.",
        "With $C = 1/(1 + 2\\alpha m)$ the surface is smooth at the pole $\\theta = 0$, on the half axis between "
        "the black holes, and comes to the apex of a cone at $\\theta = \\pi$, where the cosmic string meets the "
        "horizon. There the radius grows with distance from the pole at the rate $(1 - 2\\alpha m)/(1 + 2\\alpha m)$, "
        "one half here, and the rate falls short of 1 by the deficit angle $8\\pi\\alpha m/(1 + 2\\alpha m)$ "
        "divided by $2\\pi$. The area is $16\\pi Cm^2/(1 - 4\\alpha^2m^2) = 13.5\\pi m^2$.",
    ],
    ("nariai", "universe"): [
        "The circle of $\\chi$ times a great circle of the sphere ($\\phi = 0$ and $\\pi$) at moments of the "
        "global chart from $ct = -1.5/\\sqrt{\\Lambda}$ to $1.5/\\sqrt{\\Lambda}$, each drawn as a surface in flat "
        "space with every distance along it the metric distance. On it the metric is "
        "$\\left(\\cosh^2(\\sqrt{\\Lambda}\\,ct)\\,d\\chi^2 + d\\theta^2\\right)/\\Lambda$, a cylinder of radius "
        "$\\cosh(\\sqrt{\\Lambda}\\,ct)/\\sqrt{\\Lambda}$ and height $2\\pi/\\sqrt{\\Lambda}$, its top and bottom "
        "edges one circle, the north pole of the sphere times the circle of $\\chi$.",
        "The circle of $\\chi$ contracts to its least radius, $1/\\sqrt{\\Lambda}$, at $ct = 0$ and expands after, "
        "and every sphere keeps the radius $1/\\sqrt{\\Lambda}$, so the cylinder widens at a fixed height. At "
        "$ct = 0$ the static chart covers the half $0 < \\chi < \\pi$, between its two horizons, and a second "
        "static patch covers the other half.",
    ],
    ("nariai", "sphere"): [
        "The sphere of $\\theta$ and $\\phi$ of the Nariai universe at one event, drawn as a surface in flat space "
        "with every distance along it the metric distance. It is the product's second factor, a sphere of radius "
        "$1/\\sqrt{\\Lambda}$, the same at every event.",
    ],
    ("schwarzschild", "flamm"): [
        "The equatorial plane ($\\theta = \\pi/2$) of the Schwarzschild spacetime at one moment of $t$, "
        "drawn as a surface in flat space with every distance along it the metric distance. On "
        "it the metric is $dr^2/(1 - r_s/r) + r^2d\\phi^2$: the circle of radius $r$ has circumference "
        "$2\\pi r$, while the distance out to the next circle, $dr/\\sqrt{1 - r_s/r}$, is longer than "
        "$dr$. Both hold on the paraboloid $z^2 = 4r_s(r - r_s)$, a surface of revolution that Ludwig "
        "Flamm found in 1916.",
        "Every slice of constant $t$ passes through the bifurcation sphere $r = r_s$, where the circles "
        "are smallest, and runs on through it into a second exterior, the same paraboloid turned over. "
        "Albert Einstein and Nathan Rosen took this bridge between the two sheets as a model of a "
        "particle in 1935. Robert Fuller and John Wheeler showed in 1962 that its throat closes before "
        "light can cross it, so nothing passes from one exterior to the other.",
    ],
    ("interior_schwarzschild", "star"): [
        "The equatorial plane ($\\theta = \\pi/2$) of a static star of uniform density at one moment "
        "of $t$, drawn as a surface in flat space with every distance along it the "
        "metric distance. Inside, $g_{rr} = 1/(1 - r^2r_s/R^3)$ is the metric of a sphere of "
        "radius $\\sqrt{R^3/r_s}$, so the slice is a cap of that sphere, curved alike at every point "
        "because the density is the same everywhere. Outside it is Flamm's paraboloid.",
        "At $r = R$, $g_{rr} = 1/(1 - r_s/R)$ on both sides, so the cap meets the paraboloid in one "
        "circle with one tangent plane, and the surface is smooth across the surface of the star. The "
        "vacuum paraboloid would run on down to a throat at $r_s$. The star, at $R = 1.5\\,r_s$, ends it "
        "above there, and its circles shrink to a point at the centre.",
    ],
    ("tov", "star"): [
        "The equatorial plane ($\\theta = \\pi/2$) of a neutron star at one moment of $t$, drawn as a "
        "surface in flat space with every distance along it the metric distance. On it "
        "$g_{rr} = r/(r - 2m)$, with $m = GM(r)/c^2$ for the mass $M(r)$ inside $r$, so the surface climbs "
        "at $dz/dr = \\sqrt{2m/(r - 2m)}$: level at the centre, where $m$ grows as $r^3$, and steeper "
        "outward as the mass inside grows. Outside the star $m$ no longer grows, and the surface is Flamm's "
        "paraboloid for the star's mass.",
        "The curvature of the surface at radius $r$ is $4\\pi G(\\rho - \\bar\\rho/3)/c^2$, with "
        "$\\bar\\rho$ the mean density inside $r$. At the centre the two are equal and the surface curves "
        "like a cap. In the outer layers the density falls below a third of the mean and the surface curves "
        "like a saddle, as Flamm's paraboloid does everywhere, while Karl Schwarzschild's star of uniform "
        "density is a cap of a sphere all the way out. This star, at $2GM/c^2R = 0.29$, ends far above the "
        "throat the vacuum paraboloid would have at $r_s = 2GM/c^2$, drawn dashed below it.",
    ],
    ("morris_thorne", "wormhole"): [
        "The equatorial plane ($\\theta = \\pi/2$) of the Morris-Thorne wormhole at one moment of $t$, "
        "drawn as a surface in flat space with every distance along it the metric "
        "distance. The slice has $g_{rr} = 1/(1 - b/r)$, so $dz/dr = \\pm 1/\\sqrt{r/b - 1}$, and "
        "the shape function $b(r)$ alone fixes the surface, which is why Michael Morris and Kip Thorne "
        "gave it that name in 1988; the redshift function $\\Phi$ does not enter. At the throat "
        "$r = b_0$ the surface stands vertical and joins a second side, flat far away as the first is.",
        "The throat flares out, $d^2r/dz^2 > 0$, exactly when $\\partial_r b(b_0) < 1$, and through the "
        "Einstein equation flaring out demands a radial tension at the throat greater than its energy "
        "density, which observers passing through it fast measure as a negative energy density. With "
        "$b = b_0^2/r$ the surface is the catenoid $r = b_0\\cosh(z/b_0)$, the shape of a soap film "
        "stretched between two rings.",
    ],
    ("thin_shell_wormhole", "wormhole"): [
        "The equatorial plane ($\\theta = \\pi/2$) of Visser's thin shell wormhole at one moment of $t$ "
        "($a = 1.25\\,r_s$), drawn as a surface in flat space with every distance along it the metric distance. "
        "On either side of the throat it is Flamm's paraboloid, $z = \\pm 2\\left(\\sqrt{r_s(r - r_s)} - "
        "\\sqrt{r_s(a - r_s)}\\right)$, cut at $r = a$ before it reaches the horizon.",
        "The two halves meet at the throat with a crease: each climbs away from it at "
        "$dz/dr = \\sqrt{r_s/(a - r_s)} = 2$, so the surface turns through a finite angle there. The crease "
        "is the shell, whose surface energy density $-(c^4/2\\pi Ga)\\sqrt{1 - r_s/a}$ is negative, since the "
        "throat flares out on both sides.",
    ],
    ("rn_metric", "outside"): [
        "The equatorial plane ($\\theta = \\pi/2$) of a charged black hole at one moment of $t$ outside its "
        "outer horizon, drawn as a surface in flat space with every distance along it the "
        "metric distance. On it $g_{rr} = r^2/(r^2 - r_sr + r_q^2)$, so $dz/dr = \\sqrt{(r_sr - "
        "r_q^2)/((r - r_+)(r - r_-))}$, and the slice passes through the outer horizon's bifurcation sphere "
        "$r = r_+$, its throat, into a second exterior, as Schwarzschild's does through $r_s$.",
        "The charge pulls the throat in from $r_s$ to $r_+ = (r_s + \\sqrt{r_s^2 - 4r_q^2})/2$, and far out the "
        "surface rises as Flamm's paraboloid of the same mass does, $dz/dr \\to \\sqrt{r_s/r}$. Between the "
        "horizons $r$ is a time, and no slice of constant $t$ enters there.",
    ],
    ("rn_metric", "inside"): [
        "The equatorial plane ($\\theta = \\pi/2$) of the same black hole at one moment of $t$ inside its "
        "inner horizon, where $r$ is again a distance and $t$ a time, drawn as a surface in flat space with "
        "every distance along it the metric distance. The slice runs through the inner "
        "horizon's bifurcation sphere $r = r_-$, its widest circle, into a second region inside $r_-$, the "
        "same surface turned over.",
        "Moving in from $r_-$, $g_{rr} = r^2/((r_+ - r)(r_- - r))$ falls to $1$ at $r = r_q^2/r_s$, where the "
        "surface lies level. Nearer the singularity at $r = 0$, $g_{rr} < 1$: the circles grow faster than "
        "the distance out to them, and no surface in flat space carries that part of the slice.",
    ],
    ("kerr", "equator"): [
        "The equatorial plane ($\\theta = \\pi/2$) of a rotating black hole at one moment of "
        "Boyer-Lindquist $t$, drawn as a surface in flat space with every distance along it the "
        "metric distance. The rotation's $g_{t\\phi}$ drops out at constant $t$, and the spin "
        "enters through the circles, whose circumference is $2\\pi\\sqrt{r^2 + a^2 + 2GMa^2/c^2r}$, so the "
        "drawing's distance from the axis is this radius, which exceeds $r$. As Schwarzschild's does, the "
        "slice passes through the bifurcation sphere at $r_+$, its throat, into a second exterior.",
        "On the equator the throat's circumference is $4\\pi GM/c^2$ whatever the spin, since $r_+^2 + a^2 = "
        "2GMr_+/c^2$ there. The dotted circle is the edge of the ergosphere, $r = 2GM/c^2$ on the equator, "
        "inside which nothing can stand still against the rotation. The ergosphere leaves the shape of the "
        "slice unmarked and lies in how the slices are stacked, the rotation dragging each one round past the "
        "next.",
    ],
    ("kerr_newman", "equator"): [
        "The equatorial plane ($\\theta = \\pi/2$) of a charged rotating black hole at one moment of "
        "Boyer-Lindquist $t$, drawn as a surface in flat space with every distance along it the "
        "metric distance. The rotation's $g_{t\\phi}$ drops out at constant $t$, and the circles "
        "have circumference $2\\pi\\sqrt{r^2 + a^2 + a^2(2GMr/c^2 - r_Q^2)/r^2}$, the drawing's distance from "
        "the axis. As Schwarzschild's does, the slice passes through the bifurcation sphere at $r_+$, its "
        "throat, into a second exterior.",
        "At the horizon $r_+^2 + a^2 = 2GMr_+/c^2 - r_Q^2$, so the throat's circumference radius is $2GM/c^2 "
        "- r_Q^2/r_+$, which the charge pulls in below Kerr's $2GM/c^2$. The dotted circle is the edge of the "
        "ergosphere, where $g_{tt} = 0$ on the equator, inside which nothing can stand still against the "
        "rotation.",
    ],
    ("de_sitter", "static"): [
        "The equatorial plane ($\\theta = \\pi/2$) of de Sitter space at the moment $t = 0$ of its static "
        "chart, drawn as a surface in flat space with every distance along it the metric "
        "distance. On it $g_{rr} = 1/(1 - r^2/\\ell^2)$, with $\\ell = \\sqrt{3/\\Lambda}$, the metric of a "
        "sphere of radius $\\ell$. The static chart covers the hemisphere about its observer out to the "
        "horizon $r = \\ell$, the equator, where the surface stands vertical, and the slice runs on through "
        "the horizon's bifurcation sphere into the static patch of an observer at the antipode.",
        "The whole slice is the three sphere of radius $\\ell$ at the waist of de Sitter's hyperboloid, the "
        "smallest moment of the closed slicing, in which space is a three sphere of radius "
        "$\\ell\\cosh(ct/\\ell)$ that contracts to this waist and expands after it. The two observers can "
        "never exchange light: each hemisphere lies outside the other observer's past and future alike.",
    ],
    ("einstein_static", "sphere"): [
        "The equatorial plane ($\\theta = \\pi/2$) of the Einstein static universe at one moment of $t$, "
        "drawn as a surface in flat space with every distance along it the metric distance. On it the "
        "metric is $R^2\\left(d\\chi^2 + \\sin^2\\chi\\,d\\phi^2\\right)$, the round sphere of radius $R$, "
        "the same at every moment.",
        "The areal radius $r = R\\sin\\chi$ and Einstein's projected coordinates cover the hemisphere "
        "$\\chi < \\pi/2$ and end at the equator $r = R$, where the surface stands vertical.",
    ],
    ("schwarzschild_de_sitter", "static"): [
        "The equatorial plane ($\\theta = \\pi/2$) of Kottler's spacetime at the moment $t = 0$ of its static chart, "
        "drawn as a surface in flat space with every distance along it the metric distance. On it "
        "$g_{rr} = (1 - r_s/r - \\Lambda r^2/3)^{-1}$, larger than 1 everywhere between the horizons, so the "
        "surface climbs at $dz/dr = \\sqrt{g_{rr} - 1}$ from the throat $r_h$, its smallest circle, to $r_c$, its "
        "widest, and stands vertical at both.",
        "The slice runs through the bifurcation sphere at $r_c$ into the next static region, the same surface "
        "turned over, and through its throat at $r_h$ into another, a chain of throats and widest circles "
        "without end, one period of which is drawn. Identifying the two throats closes the slice into a space "
        "of topology $S^1 \\times S^2$.",
    ],
    ("curzon_chazy", "equator"): [
        "The plane $z = 0$ of the Curzon-Chazy particle at one moment ($m = 1$), drawn as a surface in flat space "
        "with every distance along it the metric distance. On it $g_{\\rho\\rho} = e^{2m/\\rho - m^2/\\rho^2}$ and the "
        "circle of Weyl's radius $\\rho$ has circumference $2\\pi\\rho\\,e^{m/\\rho}$, least, $2\\pi m\\,e$, at "
        "$\\rho = m$, and growing without bound as $\\rho \\to 0$.",
        "The surface narrows from far out to a neck at $\\rho = m$ and widens again toward the singularity, the ring "
        "of infinite circumference, until it lies level at $\\rho = 0.72\\,m$, where the drawing stops; the same "
        "plane at every moment is the same surface.",
    ],
    ("melvin", "universe"): [
        "The plane $z = 0$ of Melvin's universe at one moment ($B = 1$), drawn as a surface in flat space with "
        "every distance along it the metric distance. On it $g_{\\rho\\rho} = (1 + B^2\\rho^2/4)^2$ and the circle "
        "about the axis has radius $\\rho/(1 + B^2\\rho^2/4)$, which grows to $1/B$ at the Melvin radius $\\rho = 2/B$ "
        "and shrinks as $4/(B^2\\rho)$ beyond it, while the distance out to it grows as $\\rho + B^2\\rho^3/12$.",
        "The surface widens from the axis like a vase to its widest circle and then narrows into a spike of "
        "unbounded length, drawn out to $\\rho = 4/B$; every plane of constant $z$ at every moment is the same.",
    ],
    ("melvin", "ernst"): [
        "The equatorial plane ($\\theta = \\pi/2$) of Ernst's black hole at the moment $t = 0$ ($B = 1/(2r_s)$), drawn "
        "as a surface in flat space. On it $g_{rr} = (1 + B^2r^2/4)^2/(1 - r_s/r)$ and the circle of constant $r$ "
        "has radius $r/(1 + B^2r^2/4)$, so the surface stands vertical at Flamm's throat $r = r_s$, opens to its "
        "widest circle, of radius $1/B$, at $r = 2/B$, and closes beyond it into a spike, as Melvin's plane does.",
        "The slice runs through the bifurcation sphere at $r_s$ into the other exterior, the same surface turned "
        "over. Where $B r_s \\ge 2$ the widest circle is the throat itself.",
    ],
    ("levi_civita", "field"): [
        "The plane $z = 0$ around Levi-Civita's line of mass at one moment ($\\sigma = 1/4$, $C = 1$), drawn as a "
        "surface in flat space with every distance along it the metric distance. On it "
        "$g_{\\rho\\rho} = \\rho^{4\\sigma(2\\sigma - 1)}$ and the circle about the axis has radius "
        "$\\rho^{1 - 2\\sigma}/C = \\sqrt{\\rho}$, which grows as the power $2/3$ of the distance out to it, so the "
        "surface is the horn $z = (4\\sqrt{\\rho} - 1)^{3/2}/6$, level on the circle $\\rho = 1/16$ and steeper at "
        "every circle beyond it.",
        "Every plane of constant $z$ at every moment is the same. With $\\sigma = 0$ the surface is the cone of a "
        "cosmic string, or the plane where $C = 1$.",
    ],
    ("einstein_rosen_waves", "pulse"): [
        "The plane $z = 0$ of space around a pulse of Weber, Wheeler, and Bonnor ($C = a$) going out from the axis, "
        "from $ct = a$ to $8a$, each moment drawn as a surface in flat space with every distance along it the metric "
        "distance. The metric on it is $e^{2(\\gamma - \\psi)}d\\rho^2 + \\rho^2e^{-2\\psi}d\\phi^2$, so the circle at "
        "$\\rho$ has radius $\\rho e^{-\\psi}$, and on the axis $\\gamma = 0$ and the surface starts level.",
        "As the pulse goes out it leaves the space behind it flat about the axis, where $\\psi$ and $\\gamma$ fall "
        "toward zero. Far ahead of it $\\psi \\to 0$ and $\\gamma \\to C^2/a^2$, and the surface becomes a cone of "
        "deficit angle $2\\pi(1 - e^{-C^2/a^2})$, $228°$ at $C = a$.",
    ],
    ("vaidya", "shell"): [
        "The equatorial plane ($\\theta = \\pi/2$) of space around a shell of radiation falling inward, from "
        "$v - r = -3\\,r_s$ to $r_s$, each moment drawn as a surface in flat space with every distance along it "
        "the metric distance. A slice of constant $v$ is a light cone, so the moments are slices "
        "of constant $v - r$, which are spacelike everywhere, inside the horizon as well, and carry the metric "
        "$(1 + 2Gm/c^2r)\\,dr^2 + r^2d\\phi^2$. Inside the shell $m = 0$ and the surface is a flat disc. "
        "Outside it $m = M$ and the surface is $z^2 = 4r_sr$, Flamm's paraboloid moved in by $r_s$, and the "
        "energy of the shell folds the surface where the two meet.",
        "The event horizon forms at the centre at $v = -2\\,r_s$, before the shell arrives, and grows "
        "through the flat interior, where nothing yet marks it, to meet the shell at $r_s$, where it stays. "
        "Once the shell has reached the centre the whole slice is Schwarzschild's, and the paraboloid runs "
        "on through the horizon to close in a spike at the singularity $r = 0$.",
    ],
    ("oppenheimer_snyder", "collapse"): [
        "The equatorial plane ($\\theta = \\pi/2$) of a star of dust collapsing from rest, as the dust's own time "
        "$\\tau$ runs from the release to $c\\tau = 4.39\\,r_s$, each moment drawn as a surface in flat space with "
        "every distance along it the "
        "metric distance. The dust is a piece of a closed universe, and its slice is "
        "a cap of a sphere of radius $a(\\tau)$ out to $\\chi = \\chi_0$, which shrinks as the dust falls while "
        "keeping its angle $\\chi_0$. Outside, the moment carries on as the moment of clocks released from rest "
        "at every radius when the dust was, Igor Novikov's slicing of Schwarzschild's exterior, and the two "
        "meet with one tangent, since $1 + 2E = \\cos^2\\chi_0$ on both sides of the surface.",
        "At the release the outside is Flamm's paraboloid and the dust sits in it as Schwarzschild's star does. "
        "As the dust falls, the cap shrinks and the outside follows it down, meeting it at the same angle at "
        "every moment. Once the surface is inside $r_s$ the slice runs through the horizon, drawn where "
        "$R = 2GM/c^2$, down to the dust, which no slice of constant Schwarzschild $t$ reaches. J. Robert "
        "Oppenheimer and Hartland Snyder worked out this collapse in 1939.",
    ],
    ("tolman_bondi", "cloud"): [
        "The equatorial plane ($\\theta = \\pi/2$) of a cloud of dust collapsing from rest, densest at its "
        "centre, at four moments of $t$, each drawn as its embedding in flat space (vertical scale $\\times 2$). "
        "On it $g_{rr} = (\\partial_rR)^2/(1 + 2E)$, with $R(r, t)$ "
        "the areal radius of the shell $r$, so in the areal radius the surface climbs at "
        "$dz/dR = \\sqrt{-2E/(1 + 2E)}$, set by the energy $E = -GM(r)/c^2r$ of the shell there alone. Every "
        "shell falls on its own clock, the centre first, and the surface follows the shells as they go.",
        "Richard Tolman found these solutions in 1934, and Hermann Bondi took them up in 1947. In the "
        "spacetime diagram the cloud is marginally bound, $E = 0$, falling from rest at infinity, and then "
        "every slice of constant $t$ is flat; released from rest from the same density at $t = 0$, as it is "
        "here, its slices curve.",
    ],
    ("majumdar_papapetrou", "one_hole"): [
        "The equatorial plane ($\\theta = \\pi/2$) of a single hole ($U = 1 + m/r$) at one moment of $t$, drawn as a "
        "surface in flat space with every distance along it the metric distance. On it $g_{rr} = (1 + m/r)^2$ and "
        "the circumference radius is $r + m$, so $dz/dr = \\sqrt{2m/r + m^2/r^2}$, which grows as $m/r$ toward the "
        "horizon $r = 0$.",
        "The surface falls as $m\\ln r$ without end while its circles close on the radius $m$: the throat of the "
        "extremal Reissner-Nordström black hole is infinitely long, and it is drawn down to $r = m/50$. Down the "
        "throat the equator tends to Bertotti-Robinson's cylinder of radius $m$, $m^2dr^2/r^2 + m^2d\\phi^2$.",
    ],
    ("majumdar_papapetrou", "two_holes"): [
        "The plane $z = 0$ midway between two holes of mass parameter $m$ at $z = \\pm 2m$, at one moment of $t$, "
        "drawn as a surface in flat space with every distance along it the metric distance. On it "
        "$U = 1 + 2m/\\sqrt{\\rho^2 + 4m^2}$, the circle of coordinate radius $\\rho$ has the circumference radius "
        "$\\rho U$, and $dz/d\\rho = \\sqrt{U^2 - (U + \\rho\\,\\partial_\\rho U)^2}$ is real, since "
        "$U + \\rho\\,\\partial_\\rho U = 1 + 8m^3/(\\rho^2 + 4m^2)^{3/2}$ lies between $0$ and $U$.",
        "The surface is flat where it crosses the axis, halfway between the horizons, and far out it rises as "
        "Flamm's paraboloid of the total mass $2m$ does, $dz/d\\rho \\to \\sqrt{4m/\\rho}$.",
    ],
    ("bertotti_robinson", "equator"): [
        "The equatorial plane ($\\theta = \\pi/2$) of the Bertotti-Robinson universe at one moment of $t$, drawn "
        "as a surface in flat space with every distance along it the metric "
        "distance. The spacetime is the product of a two dimensional anti-de Sitter space, which carries $t$ and "
        "$r$, and a sphere of radius $b$, which carries $\\theta$ and $\\phi$, so the equator at one moment is "
        "a line of the one times a great circle of the other, $b^2dr^2/r^2 + b^2d\\phi^2$: a cylinder of "
        "radius $b$, on which $z = b\\ln r$ puts $r \\to 0$ and $r \\to \\infty$ both infinitely far away.",
        "Every circle on it has the circumference $2\\pi b$, and the cylinder is flat, as a rolled sheet of "
        "paper is. It is the throat of an extreme Reissner-Nordström black hole, whose funnel narrows into "
        "this cylinder, infinitely long, as its charge reaches its mass. Bruno Bertotti and Ivor Robinson "
        "found the solution independently in 1959, a uniform electromagnetic field whose energy holds both "
        "factors at the one radius $b$.",
    ],
    ("bertotti_robinson", "sphere"): [
        "The sphere of $\\theta$ and $\\phi$ of the Bertotti-Robinson universe at one moment of $t$ and one "
        "$r$, drawn as a surface in flat space with every distance along it the metric "
        "distance. It is the product's second factor, a sphere of radius $b$, the same at every $t$ and "
        "$r$, so a slice of constant $t$ is the line along the cylinder of the other view times this sphere.",
    ],
    ("stockum_dust", "dust"): [
        "The plane $z = 0$ across Cornelius Lanczos's cylinder of rotating dust at one moment of $t$, drawn "
        "about its axis as a surface in flat space with every distance along it the metric "
        "distance. On it $g_{rr} = e^{-r^2/R^2}$, and the circle of radius $r$ has circumference "
        "$2\\pi r\\sqrt{1 - r^2/R^2}$, which grows only out to $r = R/\\sqrt{2}$ and then shrinks, so the "
        "surface curls back toward the axis.",
        "At $r = 0.83\\,R$ the circles shrink faster than the distance out to them and the drawing stops. At "
        "$r = R$ they are null, and beyond it they are closed timelike curves, which Willem Jacob van Stockum "
        "found in 1937, more than a decade before Gödel's universe.",
    ],
    ("taub_nut", "equator"): [
        "The equatorial plane ($\\theta = \\pi/2$) of Taub-NUT space at one moment of $t$, drawn as a surface in "
        "flat space with every distance along it the metric distance. The NUT "
        "parameter enters $g_{t\\phi}$ through $\\cos\\theta$, which vanishes on the equator, so the slice there "
        "is a surface of revolution: $g_{rr} = (r^2 + l^2)/(r^2 - 2mr - l^2)$, with circles of circumference "
        "$2\\pi\\sqrt{r^2 + l^2}$, standing vertical at the horizon $r_+ = m + \\sqrt{m^2 + l^2}$ as Flamm's "
        "paraboloid does at $r_s$.",
        "The horizon's circumference radius is $\\sqrt{2(mr_+ + l^2)}$, larger than Schwarzschild's $2m$ for "
        "the same mass. Across the horizon lies Taub's cosmology, where $r$ is a time, so the slice ends there. The "
        "equator stays clear of the Misner string, the singular axis $\\theta = 0$ and $\\pi$, which a periodic "
        "$t$ removes only at the price of closed timelike curves through every point.",
    ],
    ("godel", "dust"): [
        "The plane $z = 0$ about one world line of the dust in Gödel's universe, at one moment of the cylindrical "
        "chart's $t$, drawn as a surface in flat space with every distance along it the "
        "metric distance. On it the circle $r$ has circumference $2\\pi\\sqrt{2}\\,\\sinh r"
        "\\sqrt{1 - \\sinh^2 r}/\\omega$, which grows out to $\\sinh^2 r = 1/2$ and then shrinks, so the surface "
        "curls back toward the axis, and at $\\sinh^2 r = 1/\\sqrt{2}$ the circles shrink faster than the "
        "distance out to them and the drawing stops.",
        "At $\\sinh r = 1$ the circles are null, and beyond it they are closed timelike curves, which Kurt "
        "Gödel found in 1949. The universe is homogeneous, so every world line of the dust is an axis like "
        "this one, and it has no moment of time that is space everywhere: this surface is a moment only near "
        "its axis.",
    ],
    ("ellis_bronnikov", "wormhole"): [
        "The equatorial plane ($\\theta = \\pi/2$) of the Ellis-Bronnikov wormhole at one moment of $t$, drawn "
        "as a surface in flat space with every distance along it the metric distance. Its $r$ is "
        "the proper distance from the throat, running from $-\\infty$ on one side to $\\infty$ on the other, so "
        "the circles of constant $r$ stand at equal steps along the surface, and the circle at $r$ has "
        "circumference $2\\pi\\sqrt{r^2 + \\ell^2}$. Both hold on the catenoid $\\sqrt{r^2 + \\ell^2} = "
        "\\ell\\cosh(z/\\ell)$, one piece through the throat at $r = 0$, where the circles are smallest and the "
        "surface stands vertical.",
        "The same surface, with the same metric, is the Morris-Thorne wormhole of Michael Morris and Kip "
        "Thorne, who set it out in 1988 as the simplest traversable wormhole, with the shape function $b = "
        "\\ell^2/R$ of the areal radius $R = \\sqrt{r^2 + \\ell^2}$, unaware that Homer Ellis and Kirill "
        "Bronnikov had each found it in 1973. The areal radius turns back at the throat, so it covers one side "
        "at a time, while the proper $r$ runs straight through.",
    ],
    ("cosmic_string", "unroll"): [
        "The cone of an ideal cosmic string ($4G\\mu/c^2 = 0.1$) cut along $\\phi = 0$ and unrolled onto the plane "
        "without stretching, every distance from the apex and every circle keeping its length. Laid flat it is a "
        "plane missing the wedge between the cut's two edges, the deficit angle $\\delta = 8\\pi G\\mu/c^2 = 36°$.",
    ],
    ("global_monopole", "monopole"): [
        "The equatorial plane ($\\theta = \\pi/2$) around a global monopole with no mass at its centre, at one "
        "moment of $t$, drawn as a surface in flat space with every distance along it the metric distance. On it "
        "the metric is $dr^2 + (1 - \\Delta)r^2d\\phi^2$, so the circle at the proper distance $r$ from the centre "
        "has circumference $2\\pi\\sqrt{1 - \\Delta}\\,r$ and the surface is a cone of half angle "
        "$\\arcsin\\sqrt{1 - \\Delta}$. Cut along a line from its apex and laid flat, it is a plane missing a wedge "
        "of angle $2\\pi(1 - \\sqrt{1 - \\Delta})$, $36°$ at $\\Delta = 0.19$.",
        "Every plane through the centre is the same cone. A light ray in one runs straight across it laid flat, so "
        "it turns through $\\pi(1/\\sqrt{1 - \\Delta} - 1)$, $20°$ at $\\Delta = 0.19$, whatever its impact "
        "parameter. The apex is the centre, where the Kretschmann scalar $4\\Delta^2/\\left((1 - \\Delta)^2r^4\\right)$ "
        "diverges.",
    ],
    ("global_monopole", "black_hole"): [
        "The equatorial plane ($\\theta = \\pi/2$) of Letelier's black hole in a cloud of strings, at one moment "
        "of $t$, drawn as a surface in flat space with every distance along it the metric distance. On it the "
        "metric is $dr^2/(1 - \\Delta - r_s/r) + r^2d\\phi^2$, so $dz/dr = \\sqrt{(\\Delta r + r_s)/((1 - \\Delta)r - "
        "r_s)}$, vertical at the horizon $r_h = r_s/(1 - \\Delta)$ and tending to the slope "
        "$\\sqrt{\\Delta/(1 - \\Delta)}$ of the monopole's cone far out.",
        "The slice passes through the bifurcation sphere $r = r_h$, its smallest circle, and runs on into a second "
        "exterior, the same surface turned over. At $\\Delta = 0$ it is Flamm's paraboloid of the Schwarzschild "
        "spacetime.",
    ],
    ("cosmic_string", "cone"): [
        "The plane $z = 0$ across a straight cosmic string at one moment of $t$, drawn as a surface in flat "
        "space with every distance along it the metric distance. The "
        "circle of radius $r$ about the string has circumference $2\\pi(1 - 4G\\mu/c^2)\\,r$, short of "
        "$2\\pi r$, so the surface is a cone of half angle $\\arcsin(1 - 4G\\mu/c^2)$, flat everywhere but "
        "at its apex. Cut along a line from the apex and laid flat, it is a plane with a wedge of angle "
        "$\\delta = 8\\pi G\\mu/c^2$ missing, the deficit angle.",
        "Richard Gott's core, a cylinder of uniform density, rounds the apex off. Its slice is a cap of "
        "a sphere of radius $\\ell$, and it meets the cone where their tangents agree, at $\\cos\\chi_0 = "
        "1 - 4G\\mu/c^2$. The cone of an ideal string runs on below the cap to its apex.",
    ],
    ("kantowski_sachs", "dust"): [
        "The equatorial plane ($\\theta = \\pi/2$) of the dust universe symmetric in time ($\\kappa = 0$) as $\\eta$ "
        "runs from $0$ to $1.2$, each moment drawn as a surface in flat space with every distance along it the "
        "metric distance. On it the metric is $a^2dr^2 + b^2d\\phi^2$ with $a = 1 + \\eta\\tan\\eta$ and "
        "$b = b_0\\cos^2\\eta$, so each moment is a flat cylinder of radius $b$ on which the stretch "
        "$|r| \\le b_0$ is $2ab_0$ long.",
        "The cylinder is widest and shortest at $\\eta = 0$, where $b = b_0$ and $a = 1$. Toward the last "
        "singularity, $\\eta \\to \\pi/2$, its radius goes to zero while the length between any two circles of "
        "constant $r$ grows without bound, and the expansion from the first singularity is the same collapse "
        "reversed in time.",
    ],
    ("kantowski_sachs", "vacuum"): [
        "The equatorial plane ($\\theta = \\pi/2$) inside the horizon of a Schwarzschild black hole at five moments "
        "of $T$, from $0.9\\,r_s$ to $0.1\\,r_s$, each drawn as a surface in flat space with every distance along it "
        "the metric distance. On it the metric is $(r_s/T - 1)\\,dr^2 + T^2d\\phi^2$, so each moment is a flat "
        "cylinder of radius $T$ on which the stretch $|r| \\le r_s$ is $2r_s\\sqrt{r_s/T - 1}$ long.",
        "At the horizon, $T = r_s$, the cylinder has its greatest radius and no length. As $T$ falls toward the "
        "singularity its circles shrink and the same stretch of $r$ lengthens without bound, so a body falling in "
        "is squeezed around the axis and stretched along it.",
    ],
    ("misner", "cylinders"): [
        "The slice $z = 0$ of Misner space at four moments of the Milne time $t$ in the region $T < 0$, each "
        "drawn as a surface in flat space with every distance along it the metric distance. On it the metric "
        "is $dy^2 + c^2t^2d\\chi^2$ with $\\chi$ periodic in $\\psi_0/2$, so each moment is a flat cylinder of "
        "circumference $\\psi_0c|t|/2$, here $2\\pi c|t|$.",
        "The cylinders narrow as $t$ climbs toward $0$, where the circles become the closed null geodesics of "
        "the chronology horizon. Beyond it, where $T > 0$, the same circles are closed timelike curves.",
    ],
    ("minkowski", "plane"): [
        "The equatorial plane ($\\theta = \\pi/2$) of Minkowski space at one moment of $t$, drawn as a "
        "surface in flat space with every distance along it the metric distance. On it $g_{rr} = 1$ "
        "and the circle of radius $r$ has circumference $2\\pi r$, so the surface is the flat plane itself.",
        "Every other slice is measured against this flat plane. Where a circle's circumference falls short of "
        "$2\\pi$ times the distance out to it, as around a star, the plane curves into a bowl; where it exceeds it, "
        "as in anti-de Sitter space, no surface of revolution in flat space carries the plane. Hermann Minkowski set "
        "out in 1908 the geometry in which space at one moment of any inertial observer is this flat space of "
        "Euclid.",
    ],
    ("anti_de_sitter", "hyperboloid"): [
        "The equatorial plane ($\\theta = \\pi/2$) of anti-de Sitter space at the moment $t = 0$ of its static chart, "
        "drawn as a surface in three dimensional Minkowski space with every distance along it, measured with $dX^2 + "
        "dY^2 - dZ^2$, the metric distance. On it $g_{rr} = 1/(1 + r^2/L^2)$ while the "
        "circle of radius $r$ has circumference $2\\pi r$, so every circle grows faster than the distance out to it, "
        "which no surface of revolution in flat space allows. In Minkowski space the plane is one sheet of the "
        "hyperboloid $(Z + L)^2 - X^2 - Y^2 = L^2$, and the whole hyperbolic plane of curvature $-1/L^2$ lies on it.",
        "Where the sheet is steep a step along it is shorter than it looks: the circles at $L$, $2L$, $3L$, and $4L$ "
        "stand $0.88$, $0.56$, $0.37$, and $0.28\\,L$ apart. The sheet nears the light cone of the space it is drawn "
        "in, dashed, without ever reaching it, and the conformal boundary of anti-de Sitter space lies along that "
        "cone at infinity. Wilhelm Killing in 1880 and Henri Poincaré in 1881 each described the hyperbolic plane as "
        "this sheet, and David Hilbert proved in 1901 that no surface in flat space carries the whole of it.",
    ],
    ("malament_hogarth", "plane"): [
        "The plane $z = 0$ about the removed event of a Malament-Hogarth spacetime as $t$ runs from $ct = -0.7$ up to "
        "the removed event, each moment drawn as a surface in flat space with every distance along it the metric "
        "distance. "
        "The metric is $\\Omega^2$ times Minkowski's, so the circle of radius $s$ has circumference $2\\pi s\\Omega$ "
        "and every distance is $\\Omega$ times its flat value: where $\\Omega$ grows toward the removed event the "
        "plane sinks into a well, flat again beyond the unit ball where $\\Omega = 1$.",
        "As $t$ runs up to the moment of the removed event the well deepens without limit, and at $ct = 0$ it has no "
        "bottom: its circles close in on the radius $1$ while the distance down to them grows as $\\ln(1/s)$, a tube "
        "that runs on for ever. The computer of John Earman and John Norton's toy of 1993 rides the axis into the "
        "removed event and ages $\\int\\Omega\\,c\\,dt$, which diverges. Because $\\Omega$ depends only on $c^2t^2 + "
        "x^2 + y^2 + z^2$, the tube from its rim down to the circle $s$ is exactly as long as the computer's clock "
        "runs from $ct = -1$ to $ct = -s$, so the infinite time the computer spends on its way is the infinite length "
        "of the tube.",
    ],
    ("mixmaster", "sphere"): [
        "The great two sphere of the Mixmaster universe's three sphere as its proper time runs from $c\\tau = "
        "0.09\\,m$ to $2.83\\,m$, each moment drawn as a surface in flat space with every distance along it the "
        "metric distance. Every "
        "great sphere of a Mixmaster slice is congruent to every other, and the three great circles in which it meets "
        "its planes of symmetry have circumferences $4\\pi a_1$, $4\\pi a_2$, and $4\\pi a_3$, so the sphere carries "
        "all three scale factors. When all three agree it is a round sphere of radius $2a$, the equator of a round "
        "three sphere.",
        "The moments are those of Abraham Taub's universe of 1951, the vacuum member of the family with $a_1 = a_2$, "
        "which is the region across the horizon of Taub-NUT space where its $r$ is a time. Its fibres, the circles "
        "of $\\psi$, open from nothing at its first horizon and close again at its last, so its great sphere, a "
        "surface of revolution about the axis through the identity with a fibre as its equator, runs from a sphere "
        "squeezed about its equator, through a wide band and a round sphere, to two lobes joined at a narrow waist. "
        "While $a_3$ is more than $2/\\sqrt{3}$ times $a_1$, as at $c\\tau = 0.39\\,m$, where it is $2.69$ times, the "
        "curvature about the poles is negative, and only the "
        "band about the equator has a surface of revolution in flat space. Charles Misner let all three scale "
        "factors differ in 1969 and found them oscillating without end toward the singularity; the great sphere of "
        "such a moment has no axis.",
    ],
    ("kasner", "tube"): [
        "The world tube of a ring of particles at rest in the chart ($(p_1, p_2, p_3) = (-2/7, 3/7, 6/7)$), each "
        "moment from $t = 1/4$ to $t = 2$ an ellipse at the height of its time. Its cross section reaches "
        "$t^{p_1}\\ell$ along $x$ and $t^{p_3}\\ell$ along $z$, so the tube narrows along $x$ as it widens along $z$.",
    ],
    ("bianchi", "tube"): [
        "The world tube of a ring of dust at rest in the chart (rates $(-0.5, 1.5, 2.0)\\,\\bar H$ where "
        "$a_1 = a_2 = a_3 = 1$), each moment from $c\\bar Ht = 0.10$ to $2.00$ an ellipse at the height of its "
        "time. Its cross section reaches $a_1\\ell$ along $x$ and $a_3\\ell$ along $z$: squeezed along $z$ near the "
        "singularity, narrowest along $x$ at $c\\bar Ht = 1.18$, and widening along both afterwards.",
    ],
    ("aichelburg_sexl", "tube"): [
        "The world tube of a ring of free particles at rest about the axis the source moves along before the shock, "
        "each wave front from $u = -1$ to $u = 1.5$ a circle at the height of its $u$. Behind the shock the tube "
        "narrows as a cone and closes on the axis at $u = 2$, in units of $8GE/c^4$.",
    ],
    ("aichelburg_sexl", "ring"): [
        "The wave front of the Aichelburg-Sexl shock at four values of $u = ct - z$, each drawn as a surface in flat "
        "space with every distance along it the metric distance. A surface of constant $u$ has the metric "
        "$dx^2 + dy^2$ whatever $v$ is on it, so the drawing is a flat disc, and the shock shows in a ring of free "
        "particles at rest on the circle of radius $8GE/c^4$ about the axis before it arrives.",
        "Crossing the shock, each particle gains the velocity $d\\rho/du = -4GE/(c^4\\rho)$ toward the axis, so the "
        "ring stays a circle of radius $\\rho - (4GE/c^4)\\,u/\\rho$ and closes on the axis at $u = "
        "c^4\\rho^2/4GE$, twice $8GE/c^4$ for this ring. A ring of radius $\\rho$ focuses at a $u$ growing as "
        "$\\rho^2$, so the shock focuses like a lens whose focal length grows with the square of the distance from "
        "its axis.",
    ],
    ("pp_wave", "tube"): [
        "The world tube of a ring of free particles at rest before the pulse ($A = e^{-u^2}/L^2$), each wave front "
        "from $cu = -3\\,L$ to the focus at $cu = 0.66\\,L$ an ellipse at the height of its $u$. The pulse stretches "
        "the tube along $x$ and squeezes it along $y$ until it closes on a segment of the $x$ axis at the focus.",
    ],
    ("khan_penrose", "tube"): [
        "The world tube of a ring of free particles at rest on the wave front where both waves have passed ($L = 1$), "
        "each moment on $\\sigma = 0$ from the collision at $\\tau = 0$ to $\\tau = 1.3$ an ellipse at the height of "
        "its $\\tau$. The two waves stretch the tube along $x$ and squeeze it along $y$, and toward the singularity "
        "at $\\tau = \\pi/2$ it flattens without bound onto the plane of $x$ and $\\tau$.",
    ],
    ("khan_penrose", "ring"): [
        "The wave front, the plane of $x$ and $y$, at four moments on $\\sigma = 0$ where both waves have passed, "
        "each drawn as a surface in flat space with every distance along it the metric distance. At each moment the "
        "front has the metric $g_{xx}dx^2 + g_{yy}dy^2$ with constant coefficients, so the drawing is a flat disc, "
        "and the waves show in a ring of free particles at rest on the circle $x^2 + y^2 = \\ell^2$ before they "
        "arrive.",
        "At the collision, $\\tau = 0$, the ring is still that circle, and afterwards it is the ellipse reaching "
        "$(1 + \\sin\\tau)\\ell/\\sqrt{\\cos\\tau}$ along $x$ and $(1 - \\sin\\tau)\\ell/\\sqrt{\\cos\\tau}$ "
        "along $y$. The area it encloses falls as $\\pi\\ell^2\\cos\\tau$ and vanishes at the curvature "
        "singularity $\\tau = \\pi/2$, where the ring is drawn out without bound along $x$ and closes along $y$. "
        "K. A. Khan and Roger Penrose found this spacetime in 1971, the exact vacuum solution for two impulsive plane "
        "waves colliding head on.",
    ],
    ("kasner", "ring"): [
        "The plane $y = 0$ of Kasner's universe at four moments of $t$, each drawn as a surface in flat space with "
        "every distance along it the metric distance. At every moment the plane is flat, $t^{2p_1}dx^2 + "
        "t^{2p_3}dz^2$ being Euclid's plane with its axes scaled, so the drawing is a flat disc, and the uneven "
        "expansion shows in a ring of particles at rest in the chart, which stay at rest because the metric has no "
        "$\\Gamma^i{}_{tt}$.",
        "The ring is the circle $x^2 + z^2 = \\ell^2$ at $t = 1$, and at time $t$ the ellipse reaching $t^{p_1}\\ell$ "
        "along $x$ and $t^{p_3}\\ell$ along $z$. With $(p_1, p_2, p_3) = (-2/7, 3/7, 6/7)$ the direction $x$ "
        "contracts while $y$ and $z$ expand, so toward the singularity at $t = 0$ every sphere of particles is drawn "
        "out into a needle along $x$. Edward Kasner found the solution in 1921: its exponents sum to $1$, so volumes "
        "grow as $t$, and the vacuum field equations require their squares to sum to $1$ as well.",
    ],
    ("bianchi", "ring"): [
        "The plane $y = 0$ of a Bianchi type I universe of dust at four moments of cosmic time, each drawn as a "
        "surface in flat space with every distance along it the metric distance. At every moment the plane is "
        "flat, $a_1^2dx^2 + a_3^2dz^2$ being Euclid's plane with its axes scaled, so the drawing is a flat disc, and "
        "the uneven expansion shows in a ring of the dust itself, whose grains stay at rest in the chart.",
        "The ring is the circle $x^2 + z^2 = \\ell^2$ at the moment all three scale factors are $1$, and at every "
        "other moment the ellipse reaching $a_1\\ell$ along $x$ and $a_3\\ell$ along $z$. Near the singularity the dust "
        "behaves as Kasner's vacuum does, drawn out along $x$ and flattened along $z$, and then the contraction along "
        "$x$ turns round: $a_1$ reaches its least, $0.89$, at $c\\bar Ht = 1.18$, and afterwards every "
        "direction expands. Luigi Bianchi sorted the homogeneous geometries of three dimensions into his nine types "
        "in 1898, and type I is the one whose slices are flat.",
    ],
    ("pp_wave", "ring"): [
        "The wave front of a plane gravitational wave at four values of its retarded time $u$, each drawn as a surface "
        "in flat space with every distance along it the metric distance. A surface of constant $u$ has "
        "the metric $dx^2 + dy^2$ whatever $v$ is on it, so the drawing is a flat disc, and the wave shows in a ring "
        "of free particles at rest on the circle $x^2 + y^2 = L^2$ before the pulse arrives.",
        "The pulse, $A = e^{-u^2}/L^2$ of the plus polarisation, pulls the particles as $d^2x/d(cu)^2 = Ax$ and "
        "$d^2y/d(cu)^2 = -Ay$, stretching the ring along $x$ and squeezing it along $y$, and the area it encloses "
        "falls although the spacetime is a vacuum: the Weyl curvature shears the ring, and the shear alone focuses "
        "it. At $cu = 0.66\\,L$ every particle reaches the $x$ axis at once, and past it the ring turns inside out. "
        "Roger Penrose showed in 1965 that this focusing keeps every plane wave spacetime from being globally "
        "hyperbolic.",
    ],
    ("krasnikov", "plane"): [
        "The tilt of the light cone in the Krasnikov tube, $1 - k$, drawn as a height over the plane of the tube's "
        "axis at the moment $ct = 5\\rho_0$, the height the tilt alone ($\\rho_0$ for $1 - k = 1$). Along the back "
        "edge of the light cone $c\\,dt = -k\\,dx$, so a light signal sent home over a "
        "length $L$ arrives $(1 - k)L/c$ sooner than in flat space: $1 - k$ is $0$ outside the tube, $1$ on the curve "
        "$k = 0$, where the signal arrives at the moment it left, and close to $1.8$ deep inside, where it arrives "
        "$0.8L/c$ before it left.",
        "The ship left $x = 0$ at $t = 0$ at the speed of light and reached the far end, $x = D = 4\\rho_0$, a unit "
        "of $ct$ before this moment, so the whole tube stands behind it as a ridge along the path; at any earlier "
        "moment the ridge ends at the ship, since nothing is built ahead of it. Inside the curve $k = 0$ the "
        "direction along the tube is a time, so a slice of constant $t$ is no moment of space there, and a round "
        "trip over the tube, out at the speed of light and home along the tipped cones, takes $(1 + k)D$ of $ct$, "
        "which for Krasnikov's $k = -1 + \\delta$ is $\\delta D = 0.8\\rho_0$ where flat space takes "
        "$2D = 8\\rho_0$. Serguei Krasnikov proposed the tube in 1995, and in 1997 Allen Everett and Thomas Roman "
        "built it in four dimensions and showed that its wall needs negative energy.",
    ],
    ("alcubierre", "plane"): [
        "The expansion $\\theta$ of the observers who ride the slices of Alcubierre's warp drive, drawn as a height "
        "over the plane $z = 0$ of the ship's path at the moment $t = 0$, the height $\\theta$ alone ($R$ for an "
        "expansion of $4c/R$) and the slices themselves flat. The observers are "
        "carried along $x$ at $v_sf$ times the speed of light, faster than light inside the circle $v_sf = 1$, and a "
        "small volume of them changes at the rate $\\theta = c\\,v_s\\,\\partial_xf = c\\,v_s\\,\\frac{x - "
        "x_s}{r_s}\\frac{df}{dr_s}$.",
        "Ahead of the ship $f$ falls toward the front, so the observers behind catch up with those ahead and space "
        "contracts, a trough; behind it they draw apart and space expands, a crest. Each is deepest or highest on "
        "the path one $R$ from the ship, $4.00\\,c/R$, and the level lines mark half of that. Miguel Alcubierre "
        "drew this surface in 1994 over $x$ and $\\rho = \\sqrt{y^2 + z^2}$. The energy density the same "
        "observers measure is $-\\frac{c^4v_s^2}{32\\pi G}\\frac{y^2 + z^2}{r_s^2}\\left(\\frac{df}{dr_s}\\right)^2$, "
        "negative in the wall except on the path, where it vanishes.",
    ],
    ("natario", "plane"): [
        "The energy density $\\varepsilon$ measured by the observers who ride the slices of Natário's warp bubble, "
        "drawn as a height over the plane $z = 0$ of the ship's path at the moment $t = 0$, the height $\\varepsilon$ "
        "alone ($-R$ for $\\varepsilon = -2c^4/\\pi GR^2$) and the slices themselves flat. The observers are carried "
        "by the flow of space $X$, whose divergence vanishes, so no volume of them grows or shrinks anywhere, and "
        "$\\varepsilon = -c^4K_{ij}K^{ij}/16\\pi G$, negative wherever the flow shears space.",
        "Inside the bubble space moves forward at $v_s$ with the ship and outside it is at rest, so $\\varepsilon$ "
        "vanishes in both, and the wall between them is a moat, deepest beside the ship, $-0.68\\,c^4/GR^2$ at "
        "$r_s = 0.88R$ across the path, and $-0.48\\,c^4/GR^2$ on the path one $R$ ahead and behind. The lines "
        "marked are lines of the flow, each closing back through the wall, as the flow of a fluid that cannot be "
        "compressed closes round an obstacle. José Natário built the drive in 2002, carrying the ship with no "
        "contraction of space ahead of it and no expansion behind.",
    ],
    ("domain_wall", "moments"): [
        "The equatorial plane ($\\theta = \\pi/2$) of the global chart as $kct$ runs from $-1$ to $1$, from the "
        "centre of one side through the wall to the centre of the other, each moment drawn as a surface in three "
        "dimensional Minkowski space with every distance along it, measured with $dX^2 + dY^2 - dZ^2$, the metric "
        "distance. On it $g_{zz} = 1$ while the circle at $z$ has radius $(1 - k|z|)\\cosh(kct)/k$, so away from "
        "$ct = 0$ every circle grows faster than the distance out to it, and each side is a cone, "
        "$Z = z\\sinh(kct)$.",
        "The two cones meet rim to rim in a fold along the wall, the circle of radius $\\cosh(kct)/k$, and each "
        "closes at its apex, the centre of its side at $T = 0$, where the horizon meets it. At $ct = 0$ both are the "
        "flat disc of radius $1/k$, the moment the wall stops, and the whole equator is that disc taken twice, "
        "joined at its rim.",
    ],
    ("milne", "hyperboloids"): [
        "The equatorial plane ($\\theta = \\pi/2$) of the Milne universe as cosmic time runs from $ct = 0.5$ to $3$, "
        "each moment drawn as a surface in three dimensional Minkowski space with every distance along it, measured "
        "with $dX^2 + dY^2 - dZ^2$, the metric distance. On it $g_{\\chi\\chi} = c^2t^2$ while the circle at $\\chi$ "
        "has circumference $2\\pi ct\\sinh\\chi$, so every circle grows faster than the distance out to it, and the "
        "moment is the hyperbolic plane of curvature $-1/c^2t^2$. It lies whole on one sheet of the hyperboloid "
        "$Z^2 - X^2 - Y^2 = c^2t^2$.",
        "That Minkowski space is the plane $\\theta = \\pi/2$ of the inertial chart, with $Z = cT$, so each sheet "
        "stands where its moment lies in spacetime. The sheets nest inside the light cone of the event $T = R = 0$, "
        "dashed, and each nears it as $\\chi$ grows. The circle of comoving particles at each $\\chi$ moves out along "
        "the straight line $Z = \\rho\\coth\\chi$ through the apex of the cone, at the speed $c\\tanh\\chi$.",
    ],
    ("frw", "closed"): [
        "The equator ($\\theta = \\pi/2$) of space in a closed universe of dust as cosmic time runs from $ct = 0.18$ "
        "to $6.10$, each moment drawn as a surface in flat space with every distance along it the "
        "metric distance. Each slice of constant $t$ is a three sphere, and its equator is a "
        "sphere of radius $a(t)/\\sqrt{k}$. The comoving circles of constant $\\chi$ keep their places on "
        "it while every distance between them grows and shrinks with $a$, from zero at the bang to the "
        "largest at $ct = \\pi/\\sqrt{k}$ and back to zero at the crunch at $ct = 2\\pi/\\sqrt{k}$.",
        "The radius $r = \\sin\\chi$ covers one hemisphere and ends at the equator $r = 1$, where "
        "$g_{rr} = a^2/(1 - kr^2)$ diverges while the distance across it stays finite. The other "
        "hemisphere is the same cap turned over.",
    ],
}


def draw(metric_id, ck):
    src = Sources()
    if metric_id in STATED:
        stops = STATED[metric_id](ck, src)
        return {"metric": metric_id, "source": src.stamps(), "views": [], "stops": stops}
    views = DRAWN[metric_id](ck, src)
    for v in views:
        v["caption"] = CAPTIONS[(metric_id, v["id"])]
    return {"metric": metric_id, "source": src.stamps(), "views": views}


def check_table():
    """Every view has a caption, and every spacetime is drawn, stated or named as not drawn,
    in exactly one of the three tables."""
    ids = {p.stem for p in build.METRICS_DIR.glob("*.json") if not build.CONFLICT_COPY.search(p.stem)}
    tables = (set(DRAWN), set(STATED), set(NOT_DRAWN))
    twice = {m for m in ids if sum(m in t for t in tables) > 1}
    neither = ids - set().union(*tables)
    if twice or neither:
        raise SystemExit(f"in more than one table: {sorted(twice)}; in none: {sorted(neither)}")
    stray = {metric for metric, _ in CAPTIONS} - set(DRAWN)
    if stray:
        raise SystemExit(f"captions for spacetimes that are not drawn: {sorted(stray)}")


def main(argv=None):
    check_table()
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--metric", action="append", default=[], help="redraw only this spacetime, repeatable")
    parser.add_argument("--verify", action="store_true", help="run and print every check instead of writing")
    args = parser.parse_args(argv)
    unknown = set(args.metric) - set(DRAWN) - set(STATED)
    if unknown:
        parser.error(f"no embedding diagram is drawn for {sorted(unknown)}")
    wanted = [m for m in [*DRAWN, *STATED] if not args.metric or m in args.metric]
    start = time.time()
    ck = Checks()
    files = {metric_id: draw(metric_id, ck) for metric_id in wanted}
    if args.verify:
        ck.report()
    failed = ck.failures()
    print(f"{len(ck.items)} checks in {time.time() - start:.1f} s, {len(failed)} failed; worst along "
          f"{ck.worst('along'):.1e}, across {ck.worst('across'):.1e}, around {ck.worst('around'):.1e}",
          file=sys.stderr if failed else sys.stdout)
    for name in failed:
        print(f"FAILED {name}", file=sys.stderr)
    if failed or args.verify:
        return 1 if failed else 0
    EMBEDDING_DIR.mkdir(parents=True, exist_ok=True)
    if not args.metric:
        for path in sorted(EMBEDDING_DIR.glob("*.json")):
            if path.stem not in DRAWN and path.stem not in STATED and not build.CONFLICT_COPY.search(path.stem):
                path.unlink()
                print(f"removed {path.relative_to(build.ROOT)}")
    for metric_id, data in files.items():
        path = EMBEDDING_DIR / f"{metric_id}.json"
        path.write_text(json.dumps(data, ensure_ascii=False, separators=(",", ":")) + "\n", encoding="utf-8")
        print(f"wrote {path.relative_to(build.ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
