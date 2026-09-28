#!/usr/bin/env python3
"""Embedding diagrams of the spacetimes, drawn from the published metrics.

An embedding diagram draws a two dimensional slice of a spacetime, the equatorial plane at
one moment, as a surface in ordinary flat three dimensional space, so that distances
measured along the surface are the distances the metric gives. Every slice drawn here is
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
construction stops and why rather than drawing past it.

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
            round, against the length the metric gives the straight coordinate line
            between them;
  around    the circumference 2 pi rho at every point against 2 pi sqrt(g_phiphi);
  joins     where two pieces meet, as a star's surface meets the exterior, they meet at
            one point with one tangent, which says g_xx agrees on both sides;
  forms     the closed form each surface is known by, as Flamm's paraboloid, the interior
            Schwarzschild cap, the catenoid, the cone, Gott's cap, the sphere and the
            cylinder, and the circumference of a rotating hole's horizon;
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
dimensions takes in the diagram files.
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
from conformal import Sources  # noqa: E402
from projections import Camera  # noqa: E402

EMBEDDING_DIR = build.EMBEDDING_DIR
FIELDS = ("coords", "parameters", "metric_components")
EQUATOR = {"theta": "pi/2"}

SAG = 2e-5          # how far the profile may stray from a chord, as a part of the drawing's size
BEND = 0.004        # and as a part of the chord: the chord then turns through at most 0.032
STEP = 1 / 90       # the longest chord of the profile, as a part of the drawing's size
DIGITS = 1e-7       # rho and z are rounded to below this part of their piece's extent
DPHI = 0.02         # the turn of the "across" chords, in radians
ALONG = 2e-4        # how far a chord may miss the proper distance, as a part of it
ACROSS = 2e-4       # the same for the chords across
AROUND = 1e-6       # how far rho may miss sqrt(g_phiphi), as a part of the drawing's size
FORM = 2e-5         # how far a point may lie from its closed form, as a part of the drawing's size
JOIN = 1e-9         # how far two pieces may miss each other in place and in tangent

CAMERA = Camera(-90, 22)
RING = 720          # points round a circle of the drawing
LAB = {"lab": 15, "small": 13}  # label sizes, in units of a figure 628 wide, as on the page


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
               the slice, g_xx picking up the cross terms and the moving coordinate's own.
    """

    def __init__(self, sources, metric_id, system_id, x, phi, fixed, params=None, functions=None, numeric=None,
                 along=None):
        _, entry, reader = nr.load(metric_id, system_id)
        sources.note(metric_id, system_id, FIELDS)
        self.metric_id, self.system_id, self.coordinate = metric_id, system_id, x
        R = reader
        names = {R._plain(n): s for n, s in R.symbol.items()}
        names.update(R.parameters)
        along = along or {}
        held = sorted(R._plain(c) for c in entry["coords"] if c not in (x, phi))
        if held != sorted([*fixed, *along]):
            raise SystemExit(f"{metric_id}/{system_id}: the slice of {x} and {phi} holds {held} fixed, "
                             f"and the table fixes {sorted(fixed)} and moves {sorted(along)}")
        self.x, self.phi = R.symbol[x], R.symbol[phi]
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
        i, j = entry["coords"].index(x), entry["coords"].index(phi)
        # The slice's tangent along x in the chart: 1 along x and d(expression)/dx along a moving
        # coordinate, which is how g_xx picks up g_xv and g_vv on a slice of constant v - r.
        tangent = [sp.Integer(0)] * len(entry["coords"])
        tangent[i] = sp.Integer(1)
        for name, expr in moved.items():
            tangent[[names[R._plain(c)] for c in entry["coords"]].index(name)] = sp.diff(expr, self.x)
        n = len(tangent)
        gxx = sum(tangent[a] * tangent[b] * g[a, b] for a in range(n) for b in range(n))
        gxp = sum(tangent[a] * g[a, j] for a in range(n))
        self.gxx, self.gxp, self.gpp = prep(gxx), prep(gxp), prep(g[j, j])
        where = f"{metric_id}/{system_id} on the slice of {x} and {phi}"
        if self.gxp != 0:
            raise SystemExit(f"{where}: g_x phi = {self.gxp}, so the slice is not a surface of revolution")
        stray = (self.gxx.free_symbols | self.gpp.free_symbols) - {self.x}
        if stray:
            raise SystemExit(f"{where}: the metric still depends on {sorted(map(str, stray))}")
        # drho/dx as g_phiphi'/(2 sqrt(g_phiphi)), so that no absolute value is differentiated.
        self.rho = sp.sqrt(self.gpp)
        self.drho = sp.diff(self.gpp, self.x) / (2 * self.rho)
        self.defect = sp.simplify(self.gxx - sp.diff(self.gpp, self.x) ** 2 / (4 * self.gpp))
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
        roots = [z for z in sp.solve(sp.numer(sp.together(1 / self.gxx)), self.x) if z.is_real]
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
        if self.numeric:
            g = float(self.gxx_at(x))
            return np.array([float(self._at(self._drho, x)) / math.sqrt(g), math.sqrt(max(float(self.defect_at(x)), 0.0) / g)])
        x0 = self.exact.get(x, sp.nsimplify(x))
        return np.array([float(sp.limit(e, self.x, x0, side)) for e in
                         (self.drho / sp.sqrt(self.gxx), sp.sqrt(self.defect / self.gxx))])

    def _local(self, root):
        """g_xx, g_phiphi and g_xx - (drho/dx)^2 as numpy functions of u = x - root, expanded in
        sympy with the root exact, so that what cancels at the root cancels exactly and they keep
        their precision where u is far below the root's own rounding, as next to a throat."""
        if root not in self._near:
            u = sp.Symbol("u")

            def near(e):
                n, d = sp.fraction(sp.together(e.subs(self.x, self.exact[root] + u)))
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
        """The length the metric gives the line from x = a to b that turns steadily through
        `turn` radians per unit of proper distance."""
        g, gpp, _, lo, hi = self._between(a, b)
        return integrate(lambda s: math.sqrt(g(s) * (1 + gpp(s) * turn ** 2)), lo, hi)

    def rise(self, a, b):
        """z(b) - z(a) of the surface of revolution, the quadrature of sqrt(g_xx - (drho/dx)^2)."""
        g, _, defect, lo, hi = self._between(a, b)

        def f(s):
            d = defect(s)
            if d < -1e-12 * max(1.0, g(s)):
                raise SystemExit(f"{self.metric_id}/{self.system_id}: g_xx - (drho/dx)^2 = {d} < 0 at "
                                 f"{self.coordinate} = {s} from {a}, where the slice has no surface of revolution")
            return math.sqrt(max(d, 0.0))
        return integrate(f, lo, hi)


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
                 marks=(), size=1.0, reference=False, legend=None):
        self.id, self.cls, self.sl, self.size = pid, cls, sl, size
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
        return decimals(max(float(np.ptp(self.rho)), float(np.ptp(self.z)), 1e-9))

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
        if self.reference:
            out["reference"] = True
        return out


# Adding 0.0 turns a -0.0 that rounding leaves into 0.0, which is how every number is written.
def fixed(v, digits):
    return round(float(v), digits) + 0.0


def decimals(extent):
    """The decimals written for rho and z of a piece whose rho or z spans `extent`, at least
    six, and more for a small piece, as a collapsing star's shrinking cap, whose chords are
    short, so the rounding is below DIGITS of it."""
    return max(6, math.ceil(-math.log10(DIGITS * extent)))


# x is written as the double itself, the shortest decimal that reads back as it: next to a
# throat g_xx diverges, and an irrational horizon, as Kerr's, rounded to fewer figures would
# fall inside it.
def significant(x):
    return float(x) + 0.0


def end_data(end):
    kind, text = end
    return {"kind": kind, "text": text} if text else {"kind": kind}


class Surface:
    """A slice as a surface of revolution about the axis z: pieces and marked circles. In a
    sequence, one surface per moment, with its `label` and `time`."""

    def __init__(self, pieces, label=None, time=None):
        self.pieces, self.label, self.time = pieces, label, time

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
            out["time"] = fixed(self.time, 6)
        out["pieces"] = [p.data() for p in self.pieces]
        out["rings"] = self.rings()
        return out


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
        chords = np.hypot(np.diff(rho), np.diff(z))
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
            space = float(np.linalg.norm(Q - P))
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

    def add(self, piece, offset=(0.0, 0.0, 0.0), fill=None):
        self.solids.append((np.asarray(offset, dtype=float), piece.rho.astype(float), piece.z.astype(float), fill))

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
    in the plane of the page, the form a figure in three dimensions takes."""

    def __init__(self, scene):
        self.scene, self.camera = scene, scene.camera
        self.fills, self.lines, self.dots, self.labels, self.legend_items = [], [], [], [], []
        self.extent = []

    def screen(self, P):
        return self.camera.screen(np.asarray(P, dtype=float))

    def line(self, cls, P, closed=False):
        """A line on the surfaces, split into the parts seen, `cls`, and the parts hidden,
        `cls`-far, cut halfway between neighbouring points of opposite kinds."""
        P = np.asarray(P, dtype=float)
        if closed:
            P = np.vstack([P, P[:1]])
        hid = self.scene.hidden(P)
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
                self.lines.append((cls + ("-far" if hid[start] else ""), nr.thin(run, 4e-4 * self.scene.size)))
            start = c + 1

    def flat(self, cls, S):
        """A line already in the plane of the page, never hidden."""
        S = np.asarray(S, dtype=float)
        self.extent.append(S)
        self.lines.append((cls, S))

    def fills_seen(self, n=420):
        """Where each fill class is the surface nearest the camera, as polygons of a grid n
        points across the drawing, painted under every line."""
        classes = sorted({fill for *_, fill in self.scene.solids if fill})
        if not classes:
            return
        E = np.vstack(self.extent)
        lo, hi = E.min(0), E.max(0)
        step = float(np.max(hi - lo)) / n
        X = np.arange(lo[0] - 2 * step, hi[0] + 3 * step, step)
        Y = np.arange(lo[1] - 2 * step, hi[1] + 3 * step, step)
        GX, GY = np.meshgrid(X, Y)
        near = self.scene.front(np.column_stack([GX.ravel(), GY.ravel()])).reshape(GX.shape)
        for cls in classes:
            which = [k for k, (*_, fill) in enumerate(self.scene.solids) if fill == cls]
            mask = np.isin(near, which).astype(float)
            gen = contourpy.contour_generator(X, Y, mask, fill_type=contourpy.FillType.OuterOffset)
            points, offsets = gen.filled(0.5, 1.5)
            for pts, offs in zip(points, offsets):
                rings = [nr.thin(pts[a:b], 5e-4 * self.scene.size) for a, b in zip(offs[:-1], offs[1:])]
                layer = {"kind": "fill", "class": cls, "points": rings[0]}
                if len(rings) > 1:
                    layer["holes"] = rings[1:]
                self.fills.append(layer)

    def plain_fill(self, cls, S):
        """A polygon already in the plane of the page, painted under every line."""
        S = np.asarray(S, dtype=float)
        self.extent.append(S)
        self.fills.append({"kind": "fill", "class": cls, "points": S})

    def dot(self, cls, P):
        S = self.screen(np.asarray(P, dtype=float)[None, :])[0]
        self.dots.append((cls, S))

    def label(self, S, text, anchor="l", cls="lab", dx=0, dy=0):
        """TeX at a point of the page, with an anchor and an offset in units of a figure 628 wide."""
        self.labels.append({"at": S, "text": text, "anchor": anchor, "class": cls, "dx": dx, "dy": dy})

    def legend(self, kind, cls, text):
        self.legend_items.append([kind, cls, text])

    def done(self, pad=0.04):
        """The figure as it is written, boxed so that every label fits inside at its size, and
        refused if two labels overlap, since each is set on its own ground over the lines."""
        self.fills_seen()
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
        layers = [dict(layer, points=rounded(layer["points"]), **({"holes": [rounded(h) for h in layer["holes"]]}
                                                                    if "holes" in layer else {}))
                  for layer in self.fills]
        order = {cls: i for i, (_, cls, _) in enumerate(self.legend_items)}
        # Hidden lines first, then the lines seen, each in the order of the legend.
        for far in (True, False):
            for cls, S in sorted(self.lines, key=lambda item: order.get(item[0].removesuffix("-far"), 99)):
                if cls.endswith("-far") == far:
                    layers.append({"kind": "line", "class": cls, "points": rounded(S)})
        for cls, S in self.dots:
            layers.append({"kind": "point", "class": cls, "at": rounded(S)})
        drawn = {layer["class"].removesuffix("-far") for layer in layers}
        for _, cls, _ in self.legend_items:
            if cls not in drawn:
                raise AssertionError(f"the legend names {cls}, which is not drawn")
        labels = [dict(L, at=rounded(L["at"])) for L in self.labels]
        return {"box": [fixed(b, 4) for b in box],
                "camera": {"azimuth": self.camera.azimuth, "elevation": self.camera.elevation},
                "layers": layers, "labels": labels, "legend": self.legend_items}


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
    is drawn in dashes, with every other meridian."""
    off = np.asarray(offset, dtype=float)
    for p in surface.pieces:
        style = "reference" if p.reference else "meridian"
        for k in range(0, meridians, 2 if p.reference else 1):
            phi = 2 * math.pi * k / meridians
            P = off + np.column_stack([p.rho * math.cos(phi), p.rho * math.sin(phi), p.z])
            fig.line(style, densify(P, 4))
        for run in outline(p, fig.camera):
            fig.line("reference" if p.reference else "outline", off + run)
        for (kind, _), x in zip(p.ends, (p.lo, p.hi)):
            if kind in ("edge", "stops") and not p.reference:
                fig.line("outline", off + circle(*p.at(x)))
        for x, cls, _ in p.marks:
            rho, z = p.at(x)
            if rho > 0:
                fig.line(cls, off + circle(rho, z))
            else:
                fig.dot(cls, off + [0, 0, z])


def densify(P, n):
    """n points on each segment of a polyline, so that a line changes from seen to hidden
    close to where it does."""
    P = np.asarray(P, dtype=float)
    s = np.linspace(0, 1, n + 1)[:-1]
    out = (P[:-1, None, :] + s[None, :, None] * (P[1:] - P[:-1])[:, None, :]).reshape(-1, 3)
    return np.vstack([out, P[-1:]])


def circle(rho, z, n=RING):
    phi = np.linspace(0, 2 * math.pi, n + 1)
    return np.column_stack([rho * np.cos(phi), rho * np.sin(phi), np.full_like(phi, z)])


def outline(p, cam):
    """Where a piece turns edge on to the camera: on the circle of each point of the profile,
    the angles at which the normal (dz cos phi, dz sin phi, -drho) is square to the line of
    sight, cos(phi - azimuth) = (drho/dz) tan(elevation), joined from point to point."""
    rho, z = p.rho, p.z
    drho, dz = np.gradient(rho), np.gradient(z)
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
        side = [np.column_stack([rho[i:j] * np.cos(a + s * ang), rho[i:j] * np.sin(a + s * ang), z[i:j]])
                for s in (1, -1)]
        # The two sides meet where the outline turns, |c| = 1; otherwise each runs to an end.
        if j - i > 1:
            if i > 0 and j < len(rho):
                runs.append(np.vstack([side[0][::-1], side[1]]))
            else:
                runs += side
        i = j
    return runs


def ring_label(fig, off, rho, z, text, side=1, cls="small", dx=8, dy=0, clear=False):
    """A label beside the right (side 1) or left end of a circle, on the page. With `clear`
    it stands past the surface's outline instead, where the outline runs outside the circle's
    end within the label's height, as a cone's sides do below its rim."""
    S = fig.screen(np.asarray(off, dtype=float) + [side * rho, 0, z])
    if clear:
        band = 0.02 * fig.scene.size
        for cls_, P in fig.lines:
            if cls_ != "outline":
                continue
            for y in (S[1] - band, S[1], S[1] + band):
                a, b = P[:-1], P[1:]
                cross = ((a[:, 1] - y) * (b[:, 1] - y) <= 0) & (a[:, 1] != b[:, 1])
                x = a[cross, 0] + (y - a[cross, 1]) * (b[cross, 0] - a[cross, 0]) / (b[cross, 1] - a[cross, 1])
                if x.size:
                    S[0] = max(S[0], x.max()) if side > 0 else min(S[0], x.min())
    fig.label(S, text, "l" if side > 0 else "r", cls, dx=side * dx, dy=dy)


# ---------------------------------------------------------------- the spacetimes

def view(vid, label, unit, surfaces, figure, **fields):
    out = {"id": vid, "label": label, "unit": unit, "surfaces": [s.data() for s in surfaces],
           "figure": figure}
    out.update({k: v for k, v in fields.items() if v is not None})
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
                       f"and $R = {R * km:.1f}$ km, the one numerical relativity tests its codes on.")]


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
    between = ("Between the horizons, $r_- < r < r_+$, $g_{rr} < 0$: $r$ is a time there and a slice of "
               "constant $t$ is not a moment of space, so nothing is drawn.")
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
    surfaces = []
    rim = ("edge", "the surface runs on, as Flamm's paraboloid moved in by $r_s$, to $r \\to \\infty$")
    for T in (-3.0, -1.5, -0.5, 1.0):
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
            ck.add(f"{where}, the two sides meet at the shell: one point",
                   float(np.max(np.abs(np.array(inside.at(R)) - outside.at(R)))), JOIN)
            ck.form(f"{where}, outside the shell z = 2 sqrt(r_s r) - 2 sqrt(r_s R)", outside,
                    lambda r, R=R: 2 * (np.sqrt(r) - np.sqrt(R)), size)
        else:
            whole = Piece("whole", "sheet", hole, 0.0, top, 0.0, 1,
                          (("apex", "the singularity $r = 0$, where the surface closes in a spike"), rim),
                          [(1.0, "horizon", None)] + rings, size)
            pieces = [whole]
            ck.form(f"{where}, z = 2 sqrt(r_s r)", whole, lambda r: 2 * np.sqrt(r), size)
        for p in pieces:
            ck.isometry(f"{where}, {p.id}", p)
        surfaces.append(Surface(pieces, label=f"$v - r = {T:g}\\,r_s$", time=T))
    fig = sequence_figure(surfaces, {"sheet": "cover"}, size, columns=2)
    fig.legend("fill", "cover", "the slice of constant $v - r$, which $v$ and $r$ cover")
    fig.legend("line", "r", "$r$ constant, at $2$, $3$ and $4\\,r_s$")
    fig.legend("line", "surface", "the shell of radiation, where the surface folds")
    fig.legend("line", "horizon", "the event horizon, which forms at the centre at $v = -2\\,r_s$ and grows through flat "
                                  "space to meet the shell at $r_s$")
    fig.legend("line", "meridian", "$\\phi$ constant, every $30°$")
    return [view("shell", "The falling shell", "$r_s$", surfaces, fig.done(),
                 settings="$r_s = 2GM/c^2 = 1$, the unit of every length; each moment is a slice of constant $v - r$.",
                 input="An imploding shell of radiation, $m = 0$ for $v < 0$ and $m = M$ for $v > 0$, as the "
                       "conformal diagram draws it.")]


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
    top = 6.0
    size = 2 * top
    surfaces = []
    for f in (0.0, 0.3, 0.6, 0.8):
        eta = f * math.pi
        a = am * (1 + math.cos(eta)) / 2
        tau = am * (eta + math.sin(eta)) / 2
        where = f"Oppenheimer-Snyder, eta = {f:g} pi"
        inner = Slice(src, "oppenheimer_snyder", "interior_comoving", "\\chi", "\\phi", {"tau": "0", **EQUATOR},
                      {"chi_0": "pi/4", "a_m": repr(am)}, {"a": repr(a)})
        outer = cloud.slice(src, tau, "-1/(2*r)")
        edge = cloud.horizon(tau, 0.0, top)
        dust_marks = [(chi0 / 3, "r2", None), (2 * chi0 / 3, "r2", None), (chi0, "surface", None)]
        out_marks = [(r, "r", None) for r in (3.0, 4.0, 5.0, top)]
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
        ck.isometry(f"{where}, the dust", dust)
        ck.isometry(f"{where}, outside", ext)
        ck.join(f"{where}, the dust meets the outside", dust, chi0, ext, R0)
        ck.form(f"{where}, the dust is a cap of a sphere of radius a", dust,
                lambda c, a=a, z0=dust.z[0]: z0 + a * (1 - np.cos(c)), size)
        if f == 0:
            ck.form(f"{where}, outside it is Flamm's paraboloid", ext, lambda r: 2 * np.sqrt(r - 1) - 2, size)
        surfaces.append(Surface([dust, ext], label=f"$c\\tau = {tau:.2f}\\,r_s$", time=tau))
    fig = sequence_figure(surfaces, {"star": "star", "sheet": "cover"}, size, columns=2)
    fig.legend("fill", "star", "the dust, a cap of a sphere of radius $a(\\tau)$, which $\\tau$ and $\\chi$ cover")
    fig.legend("fill", "cover", "outside it, the moment of clocks released from rest with the dust")
    fig.legend("line", "r2", "$\\chi$ constant in the dust, at $\\chi_0/3$ and $2\\chi_0/3$")
    fig.legend("line", "r", "the clocks released at $3$, $4$, $5$ and $6\\,r_s$")
    fig.legend("line", "surface", "the surface of the dust, $\\chi = \\chi_0$")
    fig.legend("line", "horizon", "the apparent horizon, $R = 2GM/c^2$ for the mass inside it")
    fig.legend("line", "meridian", "$\\phi$ constant, every $30°$")
    return [view("collapse", "The collapse", "$r_s$", surfaces, fig.done(),
                 settings="$R_0 = 2\\,r_s$, so that $\\chi_0 = \\pi/4$ and $a_m = 2\\sqrt{2}\\,r_s$, with $r_s = 1$ the "
                          "unit of every length; the moments are the dust's proper time $\\tau$ since the release.",
                 input="Outside the dust, the slices of Tolman-Bondi's comoving chart with no dust in it, $1 + 2E = "
                       "1 - r_s/r$, each shell of clocks released from rest at $R = r$ when the dust is.")]


def tolman_bondi(ck, src):
    """The cloud the spacetime diagram draws, density falling as 1 - r^2/r_b^2 to zero at r_b and
    2GM(r_b)/c^2 = r_b/2, so k = 2M/r^3 = (5 - 3r^2)/4 inside and 1/(2r^3) outside, r_b = 1,
    but released from rest, E = -M/r, rather than marginally bound: with E = 0 the slice has
    g_rr = (dR/dr)^2 and is a plane, which is checked for the spacetime diagram's own cloud and
    stated. Four moments of t until just before the centre is crushed, at t = pi/sqrt(5)."""
    cloud = RestCloud(lambda r: np.where(r < 1, (5 - 3 * r * r) / 4, 1 / (2 * np.maximum(r, 1.0) ** 3)),
                      lambda r: np.where(r < 1, -1.5 * r, -1.5 / np.maximum(r, 1.0) ** 4))
    cloud.check(ck, src, "Tolman-Bondi", np.linspace(0.02, 4, 80), [0.0, 0.5, 1.0, 1.3])
    E = "Piecewise((-(5 - 3*r**2)*r**2/8, r < 1), (-1/(4*r), True))"
    for t in (0.0, 0.3):
        flat = Slice(src, "tolman_bondi", "comoving_synchronous", "r", "\\phi", {"t": repr(t), **EQUATOR},
                     functions={"E": "0", "R": nr._TB_R})
        ck.plane(f"Tolman-Bondi, marginally bound, t = {t:g}", flat, np.linspace(0.02, 3, 300))
    top = 4.0
    size = 2 * top
    surfaces = []
    for t in (0.0, 0.6, 1.0, 1.3):
        sl = cloud.slice(src, t, E)
        where = f"Tolman-Bondi, t = {t:g}"
        edge = cloud.horizon(t, 0.0, top)
        horizon = [(edge, "horizon", None)] if edge is not None else []
        dust = Piece("cloud", "star", sl, 0.0, 1.0, 0.0, 1,
                     (("axis", "the centre $r = 0$, where the surface is smooth"), ("join", "the surface $r = r_b$")),
                     [(1 / 3, "r2", None), (2 / 3, "r2", None), (1.0, "surface", None)] + [h for h in horizon if h[0] < 1], size)
        ext = Piece("exterior", "sheet", sl, 1.0, top, dust.z[-1], 1,
                    (("join", "the surface $r = r_b$"), ("edge", "the slice runs on to $r \\to \\infty$")),
                    [(r, "r", None) for r in (2.0, 3.0, top)] + [h for h in horizon if h[0] > 1], size)
        ck.isometry(f"{where}, the cloud", dust)
        ck.isometry(f"{where}, outside", ext)
        ck.join(f"{where}, the cloud meets the outside", dust, 1.0, ext, 1.0)
        if t == 0:
            ck.form(f"{where}, outside it is Flamm's paraboloid", ext,
                    lambda r, z1=ext.z[0]: z1 + 2 * np.sqrt(0.5 * (r - 0.5)) - 2 * np.sqrt(0.25), size)
        surfaces.append(Surface([dust, ext], label=f"$ct = {t:g}\\,r_b$", time=t))
    fig = sequence_figure(surfaces, {"star": "star", "sheet": "cover"}, size, columns=2)
    fig.legend("fill", "star", "the cloud, $r < r_b$")
    fig.legend("fill", "cover", "outside it, the moment of clocks released from rest with the cloud")
    fig.legend("line", "r2", "the shells $r = r_b/3$ and $2r_b/3$ of the cloud")
    fig.legend("line", "r", "the clocks released at $2$, $3$ and $4\\,r_b$")
    fig.legend("line", "surface", "the surface of the cloud, $r = r_b$")
    if any(ring["class"] == "horizon" for s in surfaces for ring in s.rings()):
        fig.legend("line", "horizon", "the apparent horizon, $R = 2GM(r)/c^2$ for the mass inside it")
    fig.legend("line", "meridian", "$\\phi$ constant, every $30°$")
    return [view("cloud", "The collapsing cloud", "$r_b$", surfaces, fig.done(),
                 settings="$r_b = 1$, the unit of every length, and $2GM/c^2 = r_b/2$.",
                 input="The cloud the spacetime diagram draws, its density falling as $1 - r^2/r_b^2$ to zero at $r_b$ "
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
    size = 4.0
    tube = Piece("cylinder", "sheet", sl, math.exp(-2), math.exp(2), -2.0, 1,
                 (("edge", "the cylinder runs on for ever toward $r \\to 0$"),
                  ("edge", "the cylinder runs on for ever toward $r \\to \\infty$")),
                 [(math.exp(k), "r", None) for k in (-2, -1, 0, 1, 2)], size)
    ck.isometry("Bertotti-Robinson, the equator", tube)
    ck.form("Bertotti-Robinson, the cylinder z = b ln r", tube, np.log, size)
    ck.radius("Bertotti-Robinson, the cylinder rho = b", tube, lambda r: np.ones_like(r), size)
    equator = Surface([tube])
    fig = figure_of([equator], {"sheet": "cover"}, size)
    ring_label(fig, [0, 0, 0], *tube.at(1.0), "$r = b$")
    ring_label(fig, [0, 0, 0], *tube.at(math.exp(2)), "$e^2b$")
    ring_label(fig, [0, 0, 0], *tube.at(math.exp(-2)), "$e^{-2}b$")
    fig.legend("fill", "cover", "the equator at one moment, which $t$ and $r$ cover")
    fig.legend("line", "r", "$r$ constant, at $e^{-2}$, $e^{-1}$, $1$, $e$ and $e^2$ times $b$, a step $b$ apart along the cylinder")
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
    cut = 0.0
    fig.line("cut", densify(np.column_stack([ext.rho * math.cos(cut), ext.rho * math.sin(cut), ext.z]), 4))
    ring_label(fig, [0, 0, 0], *ext.at(1.0), "$r = \\ell$", side=-1, clear=True)
    ring_label(fig, [0, 0, 0], *ext.at(top), "$3\\ell$", side=-1)

    # The development: the cone laid flat beside it, a disc of radius r missing the wedge
    # delta = 2 pi (1 - fold), in the plane of the page. A circle of radius r on the cone
    # has length 2 pi fold r, so it is an arc of angle 2 pi fold; meridian phi lies at the
    # angle fold (phi - cut) from the cut's first edge. It is drawn at half the cone's scale,
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
    return [view("cone", "The cone", "$\\ell$", [surface], fig.done(),
                 settings="$4G\\mu/c^2 = 0.1$, so that $\\cos\\chi_0 = 0.9$ and $\\delta = 36°$, and "
                          "$\\ell = 1$, the unit of every length.")]


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
    surfaces, size = [], 4.0
    for e in moments:
        a = 1 - math.cos(e)
        sl = Slice(src, "frw", "comoving_spherical", "r", "\\phi", {"t": 0, **EQUATOR}, {"k": 1}, {"a": repr(a)})
        near = Piece("near", "sheet", sl, 0.0, 1.0, -a, 1,
                     (("axis", "the centre $\\chi = 0$"), ("join", "the equator $r = 1$, where the far hemisphere begins")),
                     [(0.5, "r", None), (math.sqrt(3) / 2, "r", None), (1.0, "chartedge", None)], size)
        far = Piece("far", "sheet2", sl, 0.0, 1.0, a, -1,
                    (("axis", "the antipode $\\chi = \\pi$"), ("join", "the equator")),
                    [(0.5, "r2", None), (math.sqrt(3) / 2, "r2", None)], size)
        where = f"FRW closed, eta = {e / math.pi:.3f} pi"
        ck.isometry(f"{where}, near hemisphere", near)
        ck.isometry(f"{where}, far hemisphere", far)
        ck.join(f"{where}, the hemispheres meet at the equator", near, 1.0, far, 1.0)
        ck.form(f"{where}, the sphere z = a(1 - sqrt(1 - r^2)) - a", near,
                lambda r, a=a: a * (1 - np.sqrt(np.maximum(1 - r * r, 0))) - a, size)
        ck.radius(f"{where}, the sphere rho = a r", near, lambda r, a=a: a * r, size)
        t = e - math.sin(e)
        surfaces.append(Surface([near, far], label=f"$ct = {t:.2f}$", time=t))
    for k, name, a_of in ((0, "flat", lambda e: e * e), (-1, "open", lambda e: math.cosh(e) - 1)):
        for e in moments:
            sl = Slice(src, "frw", "comoving_spherical", "r", "\\phi", {"t": 0, **EQUATOR}, {"k": k},
                       {"a": repr(a_of(e))})
            xs = np.linspace(1e-3, 20, 400)
            if k == 0:
                ck.plane(f"FRW flat, eta = {e / math.pi:.3f} pi", sl, xs)
            else:
                ck.stops(f"FRW open, eta = {e / math.pi:.3f} pi", sl, xs)

    offsets, x, gap = [], 0.0, 0.6
    for s in surfaces:
        a = -s.pieces[0].z[0]
        offsets.append((x + a, 0.0, 0.0))
        x += 2 * a + gap
    fig = figure_of(surfaces, {"sheet": "cover"}, 12.0, offsets=offsets, meridians=12)
    base = min(fig.screen(np.asarray(off) + [0, 0, 0])[1] - (-s.pieces[0].z[0]) for s, off in zip(surfaces, offsets))
    for s, off in zip(surfaces, offsets):
        fig.label(np.array([fig.screen(np.asarray(off))[0], base]), s.label, "t", "small", dy=8)
    fig.legend("fill", "cover", "the hemisphere $\\chi < \\pi/2$ that $r = \\sin\\chi$ covers")
    fig.legend("line", "r", "$\\chi$ constant, at $\\pi/6$ and $\\pi/3$")
    fig.legend("line", "r2", "$\\chi$ constant on the far hemisphere, at $2\\pi/3$ and $5\\pi/6$")
    fig.legend("line", "chartedge", "the equator $\\chi = \\pi/2$, where $r = 1$")
    fig.legend("line", "meridian", "$\\phi$ constant, every $30°$")
    return [view("closed", "Closed dust, $k = +1$", "$1/\\sqrt{k}$", surfaces, fig.done(),
                 settings="$k = 1$, with $1/\\sqrt{k}$ the unit of every length.",
                 input="Dust: $a = 1 - \\cos\\eta$ with $ct = \\eta - \\sin\\eta$, solved from this "
                       "spacetime's own $G^r{}_r = 0$.",
                 stops=["For $k = 0$ every slice of constant $t$ is flat, and its equator is a plane.",
                        "For $k = -1$ the circles about any point grow faster than the distance out to "
                        "them, $g_{rr} < (\\partial_r\\sqrt{g_{\\phi\\phi}})^2$ at every $r > 0$, so no "
                        "surface of revolution about a point carries the slice. A piece of it lies in "
                        "flat space on Eugenio Beltrami's pseudosphere, and David Hilbert proved in 1901 "
                        "that no surface in flat space carries the whole of it."])]




# ---------------------------------------------------------------- the spacetimes with nothing to draw

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


def warp_drive(metric_id, shift, systems):
    def stated(ck, src):
        for system_id in systems:
            flat_slices(ck, src, metric_id, system_id)
        return [f"Every slice of constant $t$ is flat: the metric on it is $dx^2 + dy^2 + dz^2$, so its equator "
                f"is a plane. The warp drive is in how the slices are "
                f"stacked, the shift {shift} carrying each one past the next, and not in the shape of any one of them."]
    return stated


def krasnikov(ck, src):
    """At constant t the published metric is k dx^2 + dr^2 + r^2 dphi^2: flat where k = 1,
    outside the tube, and not spacelike where k < 0, deep inside it."""
    _, entry, R = nr.load("krasnikov", "cylindrical")
    src.note("krasnikov", "cylindrical", FIELDS)
    g = nr.published_matrix(R, entry, "metric_components")
    k, r = R.parameters["k"], R.symbol["r"]
    want = sp.diag(k, 1, r ** 2)
    ck.exact("Krasnikov: at constant t the metric is k dx^2 + dr^2 + r^2 dphi^2",
             sp.simplify(g[1:, 1:] - want) == sp.zeros(3, 3))
    return ["Outside the tube $k = 1$ and every slice of constant $t$ is flat, its equator a plane. Deep "
            "inside it $k = -1 + \\delta$ is negative, so the direction along the tube at constant $t$ is "
            "timelike, and a surface of constant $t$ is not a moment of space there."]


def kasner(ck, src):
    flat_slices(ck, src, "kasner", "cartesian")
    return ["Every slice of constant $t$ is flat, $t^{2p_1}dx^2 + t^{2p_2}dy^2 + t^{2p_3}dz^2$ being "
            "Euclidean space with its three axes scaled, so its equator is a plane. The cosmology is in how "
            "the scales change from one moment to the next, one axis contracting while the other two expand."]


def bianchi(ck, src):
    flat_slices(ck, src, "bianchi", "type_i_cartesian")
    return ["Every slice of constant $t$ of type I is flat, $a_1^2dx^2 + a_2^2dy^2 + a_3^2dz^2$ being "
            "Euclidean space with its three axes scaled, so its equator is a plane. The cosmology is in how "
            "the three scale factors change from one moment to the next."]


def minkowski(ck, src):
    flat_slices(ck, src, "minkowski", "cartesian")
    sl = Slice(src, "minkowski", "spherical", "r", "\\phi", {"t": 0, **EQUATOR})
    ck.plane("Minkowski, the equator of the spherical chart", sl, np.linspace(1e-3, 20, 400))
    return ["Every slice of constant $t$ is flat and its equator is the plane, on which $g_{rr} = 1$ and a "
            "circle of radius $r$ has circumference $2\\pi r$: the surface every other embedding diagram is "
            "measured against."]


def anti_de_sitter(ck, src):
    """The static slice has g_rr = 1/(1 + r^2/L^2) < 1 = (drho/dr)^2 at every r > 0: with r =
    L sinh s it is L^2 (ds^2 + sinh^2 s dphi^2), the hyperbolic plane of curvature -1/L^2."""
    sl = Slice(src, "anti_de_sitter", "static_global", "r", "\\phi", {"t": 0, **EQUATOR}, {"L": 1})
    ck.stops("anti-de Sitter, the equator of the static chart", sl, np.linspace(1e-3, 20, 400))
    return ["At every $r > 0$ the circles grow faster than the distance out to them, $g_{rr} = 1/(1 + r^2/L^2) "
            "< 1$, so no surface of revolution in flat space carries the equator of a slice of constant $t$. It "
            "is the hyperbolic plane of curvature $-1/L^2$: a piece of it lies in flat space on Eugenio "
            "Beltrami's pseudosphere, and David Hilbert proved in 1901 that no surface in flat space carries "
            "the whole of it."]


# The spacetimes with no surface to draw that say why, each function checking what it states
# from the published metric and returning the sentences.
STATED = {
    "alcubierre": warp_drive("alcubierre", "$v_sf$", ["cartesian"]),
    "natario": warp_drive("natario", "$(u, v, w)$", ["cartesian_flow", "plane_flow"]),
    "lentz": warp_drive("lentz", "$\\partial_i\\phi$", ["cartesian"]),
    "krasnikov": krasnikov,
    "kasner": kasner,
    "bianchi": bianchi,
    "minkowski": minkowski,
    "anti_de_sitter": anti_de_sitter,
}


# ---------------------------------------------------------------- the tables

DRAWN = {
    "schwarzschild": schwarzschild,
    "interior_schwarzschild": interior_schwarzschild,
    "tov": tov,
    "morris_thorne": morris_thorne,
    "ellis_bronnikov": ellis_bronnikov,
    "rn_metric": rn_metric,
    "de_sitter": de_sitter,
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
}

# The spacetimes with no embedding diagram yet, for which nothing is written.
# Mixmaster's slices are squashed three spheres and a pp-wave's spacelike slices carry its
# profile, so neither is flat and neither has one surface that says anything; the
# Malament-Hogarth slices take whatever shape an arbitrary conformal factor gives them.
NOT_DRAWN = {"mixmaster", "pp_wave", "malament_hogarth"}

CAPTIONS = {
    ("schwarzschild", "flamm"): [
        "This is the equatorial plane $\\theta = \\pi/2$ of the Schwarzschild spacetime at one moment "
        "of $t$, drawn as a surface in flat space so that every distance along it is the distance the "
        "metric gives. On it the metric is $dr^2/(1 - r_s/r) + r^2d\\phi^2$: the circle of radius $r$ "
        "has circumference $2\\pi r$, while the distance out to the next circle, $dr/\\sqrt{1 - r_s/r}$, "
        "is longer than $dr$. The surface of revolution that carries both is the paraboloid $z^2 = "
        "4r_s(r - r_s)$, which Ludwig Flamm found in 1916.",
        "Every slice of constant $t$ passes through the bifurcation sphere $r = r_s$, where the circles "
        "are smallest, and runs on through it into a second exterior, the same paraboloid turned over. "
        "Einstein and Rosen took this bridge between the two sheets as a model of a particle in 1935. "
        "Robert Fuller and John Wheeler showed in 1962 that its throat closes before light can cross "
        "it, so nothing passes from one exterior to the other.",
    ],
    ("interior_schwarzschild", "star"): [
        "This is the equatorial plane $\\theta = \\pi/2$ of a static star of uniform density at one "
        "moment of $t$, drawn as a surface in flat space so that every distance along it is the "
        "distance the metric gives. Inside, $g_{rr} = 1/(1 - r^2r_s/R^3)$ is the metric of a sphere of "
        "radius $\\sqrt{R^3/r_s}$, so the slice is a cap of that sphere, curved alike at every point "
        "because the density is the same everywhere. Outside it is Flamm's paraboloid.",
        "At $r = R$ both sides give $g_{rr} = 1/(1 - r_s/R)$, so the cap meets the paraboloid in one "
        "circle with one tangent plane, and the surface is smooth across the surface of the star. The "
        "vacuum paraboloid would run on down to a throat at $r_s$. The star, at $R = 1.5\\,r_s$, ends it "
        "above there, and its circles shrink to a point at the centre instead.",
    ],
    ("tov", "star"): [
        "This is the equatorial plane $\\theta = \\pi/2$ of a neutron star at one moment of $t$, drawn as "
        "a surface in flat space so that every distance along it is the distance the metric gives. On it "
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
        "This is the equatorial plane $\\theta = \\pi/2$ of the Morris-Thorne wormhole at one moment of "
        "$t$, drawn as a surface in flat space so that every distance along it is the distance the "
        "metric gives. The slice has $g_{rr} = 1/(1 - b/r)$, so $dz/dr = \\pm 1/\\sqrt{r/b - 1}$, and "
        "the shape function $b(r)$ alone fixes the surface, which is why Michael Morris and Kip Thorne "
        "gave it that name in 1988; the redshift function $\\Phi$ does not enter. At the throat "
        "$r = b_0$ the surface stands vertical and joins a second side, flat far away as the first is.",
        "The throat flares out, $d^2r/dz^2 > 0$, exactly when $\\partial_r b(b_0) < 1$, and through the "
        "Einstein equation flaring out demands a radial tension at the throat greater than its energy "
        "density, which observers passing through it fast measure as a negative energy density. With "
        "$b = b_0^2/r$ the surface is the catenoid $r = b_0\\cosh(z/b_0)$, the shape of a soap film "
        "stretched between two rings.",
    ],
    ("rn_metric", "outside"): [
        "This is the equatorial plane $\\theta = \\pi/2$ of a charged black hole at one moment of $t$ "
        "outside its outer horizon, drawn as a surface in flat space so that every distance along it is the "
        "distance the metric gives. On it $g_{rr} = r^2/(r^2 - r_sr + r_q^2)$, so $dz/dr = \\sqrt{(r_sr - "
        "r_q^2)/((r - r_+)(r - r_-))}$, and the slice passes through the outer horizon's bifurcation sphere "
        "$r = r_+$, its throat, into a second exterior, as Schwarzschild's does through $r_s$.",
        "The charge pulls the throat in from $r_s$ to $r_+ = (r_s + \\sqrt{r_s^2 - 4r_q^2})/2$, and far out the "
        "surface rises as Flamm's paraboloid of the same mass does, $dz/dr \\to \\sqrt{r_s/r}$. Between the "
        "horizons $r$ is a time, and no slice of constant $t$ enters there.",
    ],
    ("rn_metric", "inside"): [
        "This is the equatorial plane $\\theta = \\pi/2$ of the same black hole at one moment of $t$ inside "
        "its inner horizon, where $r$ is again a distance and $t$ a time, drawn as a surface in flat space so "
        "that every distance along it is the distance the metric gives. The slice runs through the inner "
        "horizon's bifurcation sphere $r = r_-$, its widest circle, into a second region inside $r_-$, the "
        "same surface turned over.",
        "Moving in from $r_-$, $g_{rr} = r^2/((r_+ - r)(r_- - r))$ falls to $1$ at $r = r_q^2/r_s$, where the "
        "surface lies level. Nearer the singularity at $r = 0$, $g_{rr} < 1$: the circles grow faster than "
        "the distance out to them, and no surface in flat space carries that part of the slice.",
    ],
    ("kerr", "equator"): [
        "This is the equatorial plane $\\theta = \\pi/2$ of a rotating black hole at one moment of "
        "Boyer-Lindquist $t$, drawn as a surface in flat space so that every distance along it is the "
        "distance the metric gives. The rotation's $g_{t\\phi}$ drops out at constant $t$, and the spin "
        "enters through the circles, whose circumference is $2\\pi\\sqrt{r^2 + a^2 + 2GMa^2/c^2r}$, so the "
        "drawing's distance from the axis is this radius rather than $r$. As Schwarzschild's does, the "
        "slice passes through the bifurcation sphere at $r_+$, its throat, into a second exterior.",
        "On the equator the throat's circumference is $4\\pi GM/c^2$ whatever the spin, since $r_+^2 + a^2 = "
        "2GMr_+/c^2$ there. The dotted circle is the edge of the ergosphere, $r = 2GM/c^2$ on the equator, "
        "inside which nothing can stand still against the rotation. Nothing in the shape of the slice marks "
        "it: the ergosphere lies in how the slices are stacked, the rotation dragging each one round past "
        "the next.",
    ],
    ("kerr_newman", "equator"): [
        "This is the equatorial plane $\\theta = \\pi/2$ of a charged rotating black hole at one moment of "
        "Boyer-Lindquist $t$, drawn as a surface in flat space so that every distance along it is the "
        "distance the metric gives. The rotation's $g_{t\\phi}$ drops out at constant $t$, and the circles "
        "have circumference $2\\pi\\sqrt{r^2 + a^2 + a^2(2GMr/c^2 - r_Q^2)/r^2}$, the drawing's distance from "
        "the axis. As Schwarzschild's does, the slice passes through the bifurcation sphere at $r_+$, its "
        "throat, into a second exterior.",
        "At the horizon $r_+^2 + a^2 = 2GMr_+/c^2 - r_Q^2$, so the throat's circumference radius is $2GM/c^2 "
        "- r_Q^2/r_+$, which the charge pulls in below Kerr's $2GM/c^2$. The dotted circle is the edge of the "
        "ergosphere, where $g_{tt} = 0$ on the equator, inside which nothing can stand still against the "
        "rotation.",
    ],
    ("de_sitter", "static"): [
        "This is the equatorial plane $\\theta = \\pi/2$ of de Sitter space at the moment $t = 0$ of its "
        "static chart, drawn as a surface in flat space so that every distance along it is the distance the "
        "metric gives. On it $g_{rr} = 1/(1 - r^2/\\ell^2)$, with $\\ell = \\sqrt{3/\\Lambda}$, the metric of a "
        "sphere of radius $\\ell$. The static chart covers the hemisphere about its observer out to the "
        "horizon $r = \\ell$, the equator, where the surface stands vertical, and the slice runs on through "
        "the horizon's bifurcation sphere into the static patch of an observer at the antipode.",
        "The whole slice is the three sphere of radius $\\ell$ at the waist of de Sitter's hyperboloid, the "
        "smallest moment of the closed slicing, in which space is a three sphere of radius "
        "$\\ell\\cosh(ct/\\ell)$ that contracts to this waist and expands after it. The two observers can "
        "never exchange light: each hemisphere lies outside the other observer's past and future alike.",
    ],
    ("vaidya", "shell"): [
        "This is the equatorial plane $\\theta = \\pi/2$ of space around a shell of radiation falling "
        "inward, at four moments, each drawn as a surface in flat space so that every distance along it is "
        "the distance the metric gives. A slice of constant $v$ is a light cone, so the moments are slices "
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
        "This is the equatorial plane $\\theta = \\pi/2$ of a star of dust collapsing from rest, at four "
        "moments of the dust's own time $\\tau$, each drawn as a surface in flat space so that every distance "
        "along it is the distance the metric gives. The dust is a piece of a closed universe, and its slice is "
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
        "This is the equatorial plane $\\theta = \\pi/2$ of a cloud of dust collapsing from rest, densest at "
        "its centre, at four moments of $t$, each drawn as a surface in flat space so that every distance "
        "along it is the distance the metric gives. On it $g_{rr} = (\\partial_rR)^2/(1 + 2E)$, with $R(r, t)$ "
        "the areal radius of the shell $r$, so in the areal radius the surface climbs at "
        "$dz/dR = \\sqrt{-2E/(1 + 2E)}$, set by the energy $E = -GM(r)/c^2r$ of the shell there alone. Every "
        "shell falls on its own clock, the centre first, and the surface follows the shells as they go.",
        "Richard Tolman found these solutions in 1934 and Hermann Bondi took them up in 1947. The spacetime "
        "diagram draws this cloud marginally bound, $E = 0$, falling from rest at infinity, and then every "
        "slice of constant $t$ is flat; drawn here released from rest from the same density at $t = 0$, its "
        "slices curve.",
    ],
    ("bertotti_robinson", "equator"): [
        "This is the equatorial plane $\\theta = \\pi/2$ of the Bertotti-Robinson universe at one moment of "
        "$t$, drawn as a surface in flat space so that every distance along it is the distance the metric "
        "gives. The spacetime is the product of a two dimensional anti-de Sitter space, which carries $t$ and "
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
        "This is the sphere of $\\theta$ and $\\phi$ of the Bertotti-Robinson universe at one moment of $t$ "
        "and one $r$, drawn as a surface in flat space so that every distance along it is the distance the "
        "metric gives. It is the product's second factor, a sphere of radius $b$, the same at every $t$ and "
        "$r$, so a slice of constant $t$ is the line along the cylinder of the other view times this sphere.",
    ],
    ("stockum_dust", "dust"): [
        "This is the plane $z = 0$ across Cornelius Lanczos's cylinder of rotating dust at one moment of $t$, "
        "drawn about its axis as a surface in flat space so that every distance along it is the distance the "
        "metric gives. On it $g_{rr} = e^{-r^2/R^2}$, and the circle of radius $r$ has circumference "
        "$2\\pi r\\sqrt{1 - r^2/R^2}$, which grows only out to $r = R/\\sqrt{2}$ and then shrinks, so the "
        "surface curls back toward the axis.",
        "At $r = 0.83\\,R$ the circles shrink faster than the distance out to them and the drawing stops. At "
        "$r = R$ they are null, and beyond it they are closed timelike curves, which Willem Jacob van Stockum "
        "found in 1937, more than a decade before Gödel's universe.",
    ],
    ("taub_nut", "equator"): [
        "This is the equatorial plane $\\theta = \\pi/2$ of Taub-NUT space at one moment of $t$, drawn as a "
        "surface in flat space so that every distance along it is the distance the metric gives. The NUT "
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
        "This is the plane $z = 0$ about one world line of the dust in Gödel's universe, at one moment of the "
        "cylindrical chart's $t$, drawn as a surface in flat space so that every distance along it is the "
        "distance the metric gives. On it the circle $r$ has circumference $2\\pi\\sqrt{2}\\,\\sinh r"
        "\\sqrt{1 - \\sinh^2 r}/\\omega$, which grows out to $\\sinh^2 r = 1/2$ and then shrinks, so the surface "
        "curls back toward the axis, and at $\\sinh^2 r = 1/\\sqrt{2}$ the circles shrink faster than the "
        "distance out to them and the drawing stops.",
        "At $\\sinh r = 1$ the circles are null, and beyond it they are closed timelike curves, which Kurt "
        "Gödel found in 1949. The universe is homogeneous, so every world line of the dust is an axis like "
        "this one, and it has no moment of time that is space everywhere: this surface is a moment only near "
        "its axis.",
    ],
    ("ellis_bronnikov", "wormhole"): [
        "This is the equatorial plane $\\theta = \\pi/2$ of the Ellis-Bronnikov wormhole at one moment of "
        "$t$, drawn as a surface in flat space so that every distance along it is the distance the metric "
        "gives. Its $r$ is the proper distance from the throat, running from $-\\infty$ on one side to "
        "$\\infty$ on the other, so the circles of constant $r$ stand at equal steps along the surface, and "
        "the circle at $r$ has circumference $2\\pi\\sqrt{r^2 + \\ell^2}$. The surface that carries both is "
        "the catenoid $\\sqrt{r^2 + \\ell^2} = \\ell\\cosh(z/\\ell)$, one piece through the throat at $r = 0$, "
        "where the circles are smallest and the surface stands vertical.",
        "It is the surface the Morris-Thorne wormhole draws, and the same metric: Michael Morris and Kip "
        "Thorne set it out in 1988 as the simplest traversable wormhole, with the shape function $b = "
        "\\ell^2/R$ of the areal radius $R = \\sqrt{r^2 + \\ell^2}$, unaware that Homer Ellis and Kirill "
        "Bronnikov had each found it in 1973. The areal radius turns back at the throat, so it covers one side "
        "at a time; the proper $r$ runs straight through.",
    ],
    ("cosmic_string", "cone"): [
        "This is the plane $z = 0$ across a straight cosmic string at one moment of $t$, drawn as a "
        "surface in flat space so that every distance along it is the distance the metric gives. The "
        "circle of radius $r$ about the string has circumference $2\\pi(1 - 4G\\mu/c^2)\\,r$, short of "
        "$2\\pi r$, so the surface is a cone of half angle $\\arcsin(1 - 4G\\mu/c^2)$, flat everywhere but "
        "at its apex. Cut along a line from the apex and laid flat, it is a plane with a wedge of angle "
        "$\\delta = 8\\pi G\\mu/c^2$ missing, the deficit angle.",
        "Richard Gott's core, a cylinder of uniform density, rounds the apex off. Its slice is a cap of "
        "a sphere of radius $\\ell$, and it meets the cone where their tangents agree, at $\\cos\\chi_0 = "
        "1 - 4G\\mu/c^2$. The cone of an ideal string runs on below the cap to its apex.",
    ],
    ("frw", "closed"): [
        "This is the equator $\\theta = \\pi/2$ of space in a closed universe of dust at five moments "
        "of cosmic time, each drawn as a surface in flat space so that every distance along it is the "
        "distance the metric gives. Each slice of constant $t$ is a three sphere, and its equator is a "
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
