#!/usr/bin/env python3
"""Conformal diagrams of the spacetimes, drawn from the published metrics.

Every diagram is a map of a whole spacetime, or of a totally geodesic surface in it where
nothing else is faithful, into a finite region of a plane, with each null coordinate sent
to one of the drawing's null coordinates. The metric on the surface is read from a
coordinate system's published metric_components through the Reader of verify_metrics.py
beside this file, in the x^0 = cT chart the collection writes, with c = 1, by the same
load and published_matrix the null ray diagrams use.

    python3 -m venv /tmp/mfs-venv && /tmp/mfs-venv/bin/pip install sympy numpy scipy contourpy
    /tmp/mfs-venv/bin/python _tools/derivations/conformal.py
    python3 _tools/build_mfs_data.py

writes MFS/assets/data/conformal/<metric_id>.json, one file per spacetime that has a
diagram, and the second command stamps each file's version into the index. Pass
--metric <metric_id> to redraw one spacetime, repeatable, and --verify to run every check
and print it without writing anything.


The drawing
-----------

The drawing's null coordinates are p and q, each a function of one null coordinate of the
spacetime, and every view is drawn with T = p + q up and X = q - p across at one scale on
both axes. A line of constant p or q is then a light ray exactly when the gradient of p or
q is null, and it runs at 45 degrees because the data puts it there. Nothing is hand
fitted: every line, fill and boundary is the image of a coordinate line or a limit under
the map a construction derives, and each construction's docstring carries its derivation.


What is checked
---------------

Every map drawn is checked where it is drawn, by the construction that draws it, at 4000
random points of the region it covers:

  inverse   the published inverse metric on the surface is the inverse of the published
            metric there, which also says the surface is orthogonal to the coordinates
            held fixed;
  null      g^ij n_i n_j = 0 for the unit gradients n of p and of q, measured against the
            largest entry of g^ij, with the published inverse: lines of constant p and q
            are light rays;
  distinct  the two gradients are not parallel: they are two families of rays;
  future    a future directed timelike vector of the coordinates raises T.

Points where the arctangent has flattened below what a central difference resolves are
counted out. Each construction adds the limits that place its horizons, infinities and
singularities, and the published Kretschmann scalar is evaluated where a line is drawn as
a curvature singularity, which must diverge, and at every centre or throat drawn as
regular, which must stay finite. write() runs every check and refuses to write if one
fails; --verify prints them all.


Which spacetimes
----------------

DRAWN lists the seventeen spacetimes that have a diagram. NOT_DRAWN gives each of the
others the reason it has none, as prose written into its file in place of views, which the
page prints where the diagram would be, so no spacetime is left with a silent gap.


Output
------

A file carries the fields each view was drawn from, as the list `source` of
{metric, system, fields, version}, with version from build_mfs_data.diagram_source_version,
the function the null ray diagrams use; build_mfs_data.py --check fails when any of those
fields changes without the file being redrawn. A view is

  id, label       its place and the name of its button, TeX in $...$;
  system          the coordinate system whose region the view tints, if it has one;
  box             [Xmin, Xmax, Tmin, Tmax], in the drawing's units;
  layers          fills, lines, zigzags and points, each with a class, in (X, T);
  labels          TeX at a point, with an anchor and an offset (dx, dy) in units of a
                  figure 628 wide, the size the page draws;
  legend          [kind, class, text] for every class the view uses;
  caption         paragraphs of prose with TeX in $...$;
  restriction     on a view of a surface that is not the whole spacetime, what that
                  surface is, printed on the diagram itself;
  settings, input the parameter values and any declared function, as TeX prose;
  fade            {top, bottom}: how far the drawing fades out where it continues.

The file of a spacetime with no diagram carries `none` in place of `views`: the paragraphs
of NOT_DRAWN, stamped like views with the coords, parameters and metric components of each
of its coordinate systems, so a changed metric asks for the reason again.
"""

import argparse
import json
import math
import sys
import time
from pathlib import Path

import numpy as np
import sympy as sp
from scipy.integrate import cumulative_trapezoid, solve_ivp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))
import build_mfs_data as build  # noqa: E402
import null_rays as nr  # noqa: E402

CONFORMAL_DIR = build.CONFORMAL_DIR
BASE_FIELDS = ("coords", "parameters", "metric_components", "inverse_metric_components")

PI = np.pi
HALF = PI / 2
Q4 = PI / 4
EQUATOR = {"theta": "pi/2", "phi": "0"}

# Why each spacetime without a diagram has none, as prose the page prints in the diagram's
# place, so that no spacetime is left with a silent gap.
NOT_DRAWN = {
    "godel": [
        "The Gödel universe has no conformal diagram. A closed timelike curve passes through every "
        "event, so every event lies in the chronological future of every other, and there is no "
        "boundary between where light can reach and where it cannot."],
    "stockum_dust": [
        "Van Stockum's cylinder has no conformal diagram, for the reason the Gödel universe has none. "
        "Beyond $r = R$ the circles of constant $t$, $r$ and $z$ are closed timelike curves, and a "
        "timelike curve that runs out to them and back again joins any two events, so every event "
        "lies in the chronological future of every other."],
    "taub_nut": [
        "Taub-NUT has no conformal diagram. Its time is wound into the angles through "
        "$c\\,dt + 2l\\cos\\theta\\,d\\phi$, so the axis is singular unless $t$ is made periodic, and "
        "made periodic it runs in closed timelike curves outside the Taub region, wherever "
        "$\\partial_t$ is timelike, and continues across $r_+$ in more than one way."],
    "kasner": [
        "The Kasner spacetime has no conformal diagram that tells its directions apart. Each plane "
        "of $t$ and one of $x$, $y$ and $z$ is totally geodesic, and on it the metric "
        "$-c^2dt^2 + t^{2p_i}dx_i^2$ is flat in the time $c\\,t^{1 - p_i}/(1 - p_i)$, which runs from "
        "zero to infinity for every exponent below one. So all three planes give the same half of "
        "Minkowski's diamond, the singularity along its lower edge, and they cannot be told apart."],
    "bianchi": [
        "Bianchi I has no conformal diagram of its own. Each plane of $t$ and one of $x$, $y$ and $z$ "
        "is flat in the time $\\int c\\,dt/a_i$, and where that time begins and ends is up to the "
        "free scale factors. For dust each one begins at a finite value at the singularity and runs "
        "on without end, so each plane gives the half of Minkowski's diamond that Kasner's planes "
        "give, the three alike."],
    "tolman_bondi": [
        "Tolman-Bondi has no single conformal diagram. Its free functions make it anything from a "
        "Friedmann universe to the Oppenheimer-Snyder collapse, with shell crossings, black holes "
        "and naked singularities between, and each choice of them has a diagram of its own."],
    "alcubierre": [
        "The Alcubierre warp drive has no conformal diagram of its own. Where light can reach "
        "depends on the whole history of the bubble, $v_s(t)$ and $f$ from its start to its end, "
        "which the metric leaves free."],
    "natario": [
        "The Natário warp drive has no conformal diagram of its own, for the reason the Alcubierre "
        "drive has none: its horizons are set by the whole history of the flow, three free "
        "functions of all four coordinates."],
    "krasnikov": [
        "The Krasnikov tube has no conformal diagram of its own. Light returning through the tube "
        "arrives when the shape $k(t, x, r)$ lets it, so the null coordinate a diagram would be "
        "built on depends on the whole history of the tube's construction."],
    "pp_wave": [
        "The plane wave has no conformal diagram. On its axis the plane of $u$ and $v$ is flat, "
        "and its diamond would be Minkowski's, but the wave focuses light across that plane, in "
        "$x$ and $y$, so strongly that the spacetime has no Cauchy surface, as Penrose showed in "
        "1965, and its causal structure is not Minkowski's."],
    "mixmaster": [
        "The Mixmaster universe has no conformal diagram. Its slices are three spheres squashed "
        "differently along three directions, with no spherical symmetry to reduce them to a line, "
        "so no surface of two dimensions carries the causal structure of the whole. Toward the "
        "singularity its vacuum solutions pass through an endless sequence of Kasner epochs, and "
        "the directions light can cross the universe in change from one epoch to the next."],
    "lentz": [
        "Lentz's soliton has no conformal diagram. Where light can reach depends on the whole "
        "history of the soliton, as it does for every warp drive, and the potential is known only "
        "as a numerical integral over the rhomboid sources Lentz laid out, so no soliton of the "
        "class can be written down to draw."],
}
NONE_FIELDS = ["coords", "parameters", "metric_components"]


# ---------------------------------------------------------------- the drawing

def xt(p, q):
    """The drawing's axes: X = q - p across, T = p + q up."""
    p, q = np.asarray(p, dtype=float), np.asarray(q, dtype=float)
    return q - p, q + p


def point(p, q):
    X, T = xt(p, q)
    return [float(X), float(T)]


def spread(lo, hi, n=500, ends=12.0):
    """n parameter values on (lo, hi), crowding toward an infinite end or toward both ends."""
    if np.isinf(lo) and np.isinf(hi):
        return np.sinh(np.linspace(-ends, ends, n))
    if np.isinf(hi):
        return lo + np.exp(np.linspace(-ends, ends, n))
    if np.isinf(lo):
        return hi - np.exp(np.linspace(ends, -ends, n))
    s = np.linspace(-ends, ends, n)
    return lo + (hi - lo) / (1 + np.exp(-s))


def runs(p, q, tol=0.002):
    """A curve given by (p, q) as polylines in (X, T), split where a value is not finite and
    thinned to the fewest points within tol of it."""
    X, T = xt(p, q)
    X, T = np.ravel(X), np.ravel(T)
    good = np.isfinite(X) & np.isfinite(T)
    out, start = [], None
    for i in range(len(X) + 1):
        if i < len(X) and good[i]:
            start = i if start is None else start
            continue
        if start is not None and i - start > 1:
            pts = np.column_stack([X[start:i], T[start:i]])
            out.append(np.round(nr.thin(pts, tol), 4).tolist())
        start = None
    return out


def rounded(points):
    return np.round(np.asarray(points, dtype=float), 4).tolist()


# The order a legend lists its classes in, the same in every view: what is tinted, the
# coordinate lines, what is marked inside, and then the edges of the spacetime.
LEGEND_ORDER = ["cover", "cover2", "star", "past", "r", "r2", "t", "t2", "null", "world", "cone", "mark",
                "removed", "surface", "throat", "centre", "horizon", "event", "apparent", "chartedge",
                "boundary", "singular", "scri"]


class View:
    """One view of a spacetime: what it draws, its words, and the coordinates it tints."""

    def __init__(self, vid, label, box, system=None):
        self.d = {"id": vid, "label": label, "box": [round(float(b), 4) for b in box]}
        if system:
            self.d["system"] = system
        self.layers, self.labels, self.legend_items = [], [], []

    def fill(self, cls, pts):
        self.layers.append({"kind": "fill", "class": cls, "points": rounded(pts)})

    def line(self, cls, polylines, zig=False):
        for pts in polylines:
            if len(pts) > 1:
                self.layers.append({"kind": "zig" if zig else "line", "class": cls, "points": rounded(pts)})

    def segment(self, cls, a, b, zig=False):
        """The straight line from a to b, both given as (p, q)."""
        self.line(cls, [[point(*a), point(*b)]], zig)

    def curve(self, cls, p, q, zig=False, tol=0.002):
        self.line(cls, runs(p, q, tol), zig)

    def point(self, cls, pq):
        self.layers.append({"kind": "point", "class": cls, "at": rounded(point(*pq))})

    def label(self, pq, text, anchor="c", cls="lab", dx=0, dy=0):
        self.labels.append({"at": rounded(point(*pq)), "text": text, "anchor": anchor,
                            "class": cls, "dx": dx, "dy": dy})

    def label_xt(self, at, text, anchor="c", cls="lab", dx=0, dy=0):
        self.labels.append({"at": rounded(at), "text": text, "anchor": anchor,
                            "class": cls, "dx": dx, "dy": dy})

    def legend(self, cls, text):
        if not any(item[1] == cls and item[2] == text for item in self.legend_items):
            self.legend_items.append([None, cls, text])

    def set(self, **fields):
        self.d.update(fields)
        return self

    def done(self):
        """The view as it is written; every class the legend names must be drawn."""
        kinds = {}
        for layer in self.layers:
            kinds.setdefault(layer["class"], layer["kind"])
        for item in self.legend_items:
            if item[1] not in kinds:
                raise AssertionError(f"view {self.d['id']}: the legend names {item[1]}, which is not drawn")
            item[0] = kinds[item[1]]
        self.legend_items.sort(key=lambda item: LEGEND_ORDER.index(item[1]))
        order = {"fill": 0, "line": 1, "zig": 2, "point": 3}
        self.d["layers"] = sorted(self.layers, key=lambda layer: order[layer["kind"]])
        self.d["labels"] = self.labels
        self.d["legend"] = self.legend_items
        return self.d


# ---------------------------------------------------------------- reading a surface

class Sources:
    """The published fields every construction read, per coordinate system."""

    def __init__(self):
        self.read = {}

    def note(self, metric_id, system_id, fields):
        self.read.setdefault((metric_id, system_id), set()).update(fields)

    def stamps(self):
        out = []
        for (metric_id, system_id), fields in sorted(self.read.items()):
            metric = json.loads((build.METRICS_DIR / f"{metric_id}.json").read_text(encoding="utf-8"))
            system = next(s for s in metric["coordinates"] if s["id"] == system_id)
            names = sorted(fields)
            out.append({"metric": metric_id, "system": system_id, "fields": names,
                        "version": build.diagram_source_version(system, names)})
        return out


class Plane:
    """The published metric of one coordinate system on the surface of two of its
    coordinates, every other coordinate held fixed, as sympy and as numbers.

    functions  a declared function -> an expression for it, in the plain names of the
               coordinates and parameters, substituted before anything else;
    numeric    declared functions given as numbers instead: each, with its first two
               derivatives, becomes a plain symbol, and every evaluation is handed a
               callable (x0, x1) -> {name: (value, first, second)};
    axis       (coordinate, value) reached as a limit, for a surface where the coordinates
               degenerate, as the axis of Kerr does at theta = 0.
    """

    def __init__(self, sources, metric_id, system_id, plane, fixed=None, params=None,
                 functions=None, numeric=(), axis=None):
        metric, entry, reader = nr.load(metric_id, system_id)
        sources.note(metric_id, system_id, BASE_FIELDS)
        self.sources, self.metric_id, self.system_id = sources, metric_id, system_id
        self.entry, self.reader = entry, reader
        R = reader
        self.names = {R._plain(n): s for n, s in R.symbol.items()}
        self.names.update(R.parameters)
        self.x0, self.x1 = R.symbol[plane[0]], R.symbol[plane[1]]
        self.subs = {R.c: 1}
        self.subs.update({R.parameters[k]: sp.sympify(v) for k, v in (params or {}).items()})
        self.fixed = {self.names[k]: sp.sympify(v) for k, v in (fixed or {}).items()}
        self.functions = {name: sp.sympify(text, locals=self.names) for name, text in (functions or {}).items()}
        self.axis = (self.names[axis[0]], sp.sympify(axis[1])) if axis else None
        self.numeric = {}
        for name in numeric:
            fn = R.parameters[name]
            self.numeric[name] = (fn, fn.args[0], sp.symbols(f"_{name}0 _{name}1 _{name}2"))
        coords = entry["coords"]
        i, j = coords.index(plane[0]), coords.index(plane[1])
        g = nr.published_matrix(R, entry, "metric_components")
        gi = nr.published_matrix(R, entry, "inverse_metric_components")
        self.g = sp.Matrix([[self.prep(g[i, i]), self.prep(g[i, j])], [self.prep(g[j, i]), self.prep(g[j, j])]])
        self.gi = sp.Matrix([[self.prep(gi[i, i]), self.prep(gi[i, j])], [self.prep(gi[j, i]), self.prep(gi[j, j])]])
        self._metric = self.lambdify([self.g[0, 0], self.g[0, 1], self.g[1, 1],
                                      self.gi[0, 0], self.gi[0, 1], self.gi[1, 1]])
        self._K = None

    def prep(self, expr, limit=True):
        e = sp.sympify(expr)
        for name, rep in self.functions.items():
            e = e.subs(self.reader.parameters[name], rep).doit()
        for fn, var, (F0, F1, F2) in self.numeric.values():
            e = e.subs(sp.Derivative(fn, (var, 2)), F2).subs(sp.Derivative(fn, var), F1).subs(fn, F0)
        e = e.subs(self.subs)
        if self.axis:
            e = sp.limit(e, *self.axis) if limit else e.subs(*self.axis)
        return sp.simplify(e.subs(self.fixed))

    def lambdify(self, exprs):
        args = [self.x0, self.x1] + [F for _, _, Fs in self.numeric.values() for F in Fs]
        return sp.lambdify(args, exprs, "numpy")

    def call(self, f, x0, x1, fvals=None):
        x0, x1 = np.asarray(x0, dtype=float), np.asarray(x1, dtype=float)
        extra = []
        if self.numeric:
            values = fvals(x0, x1)
            extra = [np.asarray(v, dtype=float) for name in self.numeric for v in values[name]]
        out = f(x0, x1, *extra)
        return [np.broadcast_to(np.asarray(o, dtype=float), np.broadcast(x0, x1).shape) for o in out]

    def metric(self, x0, x1, fvals=None):
        """g00, g01, g11 and the published g^00, g^01, g^11 on the surface."""
        return self.call(self._metric, x0, x1, fvals)

    def kretschmann(self, x0, x1, fvals=None):
        """The published Kretschmann scalar at points of the surface."""
        if self._K is None:
            self.sources.note(self.metric_id, self.system_id, ["kretschmann"])
            K = self.prep(self.reader(nr.strip_lhs(self.entry["kretschmann"])), limit=False)
            self._K = self.lambdify([K])
        return self.call(self._K, x0, x1, fvals)[0]


def mink_pq(t, x, ell=1.0):
    """Minkowski's compactification: p = arctan((t - x)/l), q = arctan((t + x)/l)."""
    t, x = np.asarray(t, dtype=float), np.asarray(x, dtype=float)
    return np.arctan((t - x) / ell), np.arctan((t + x) / ell)


def atan_exp(x):
    """arctan(exp(x)), without overflow."""
    x = np.asarray(x, dtype=float)
    with np.errstate(over="ignore"):
        return np.where(x < 0, np.arctan(np.exp(np.minimum(x, 0))), HALF - np.arctan(np.exp(-np.maximum(x, 0))))


# ---------------------------------------------------------------- the checks

class Checks:
    """Every chart map and limit a construction draws with, checked as it is drawn."""

    def __init__(self):
        self.charts, self.limits = [], []
        self.rng = np.random.default_rng(20260924)

    def uniform(self, lo, hi, n=4000):
        return self.rng.uniform(lo, hi, n)

    def chart(self, name, plane, fmap, x0, x1, future, fvals=None):
        """fmap: (x0, x1) -> (p, q); future: (x0, x1) -> a future directed timelike vector."""
        x0, x1 = np.asarray(x0, dtype=float), np.asarray(x1, dtype=float)
        g00, g01, g11, h00, h01, h11 = plane.metric(x0, x1, fvals)
        det = g00 * g11 - g01 ** 2
        lorentzian = bool(np.all(det < 0))
        scale = np.maximum.reduce([np.abs(h00), np.abs(h01), np.abs(h11)])
        inverse = float(np.max(np.maximum.reduce([np.abs(h00 - g11 / det), np.abs(h01 + g01 / det),
                                                  np.abs(h11 - g00 / det)]) / scale))
        hx0, hx1 = 1e-6 * np.maximum(1, np.abs(x0)), 1e-6 * np.maximum(1, np.abs(x1))
        pa, qa = fmap(x0 + hx0, x1)
        pb, qb = fmap(x0 - hx0, x1)
        pc, qc = fmap(x0, x1 + hx1)
        pd, qd = fmap(x0, x1 - hx1)
        grads = [np.stack([(pa - pb) / (2 * hx0), (pc - pd) / (2 * hx1)]),
                 np.stack([(qa - qb) / (2 * hx0), (qc - qd) / (2 * hx1)])]
        units, sizes, residual = [], [], []
        for gr in grads:
            size = np.linalg.norm(gr, axis=0)
            n = gr / size
            units.append(n)
            sizes.append(size)
            residual.append(np.abs(h00 * n[0] ** 2 + 2 * h01 * n[0] * n[1] + h11 * n[1] ** 2) / scale)
        # Where the arctangent has flattened, a gradient below 1e-6 is within a factor 1e4 of
        # the rounding a central difference makes, and is left out.
        resolved = (sizes[0] > 1e-6) & (sizes[1] > 1e-6)
        null = float(np.max(np.maximum(*residual)[resolved])) if resolved.any() else math.inf
        sine = np.abs(units[0][0] * units[1][1] - units[0][1] * units[1][0])
        distinct = bool(np.all(sine[resolved] > 1e-9))
        v0, v1 = future(x0, x1)
        v0, v1 = np.broadcast_to(v0, x0.shape), np.broadcast_to(v1, x0.shape)
        timelike = bool(np.all(g00 * v0 * v0 + 2 * g01 * v0 * v1 + g11 * v1 * v1 < 0))
        rise = (grads[0][0] + grads[1][0]) * v0 + (grads[0][1] + grads[1][1]) * v1
        up = timelike and bool(np.all(rise[resolved] > 0))
        ok = (lorentzian and inverse < 1e-8 and null < 1e-3 and distinct and up
              and resolved.sum() > 0.25 * x0.size)
        self.charts.append({"name": name, "samples": int(x0.size), "resolved": int(resolved.sum()),
                            "inverse": inverse, "null": null, "distinct": distinct, "future": up, "ok": ok})

    def limit(self, name, got, want, tol=1e-6):
        got, want = np.atleast_1d(np.asarray(got, dtype=float)), np.atleast_1d(np.asarray(want, dtype=float))
        err = float(np.max(np.abs(got - want)))
        self.limits.append({"name": name, "error": err, "ok": bool(err < tol)})

    def diverges(self, name, K_near, K_nearer):
        """The Kretschmann scalar at two points approaching a line drawn as a singularity, the
        second ten times closer: it is large and still growing fast."""
        K1, K2 = np.abs(np.asarray(K_near, dtype=float)), np.abs(np.asarray(K_nearer, dtype=float))
        ok = bool(np.all(K2 > 1e8) and np.all(K2 > 50 * K1))
        self.limits.append({"name": name, "error": float(np.min(K2)), "ok": ok, "kind": "diverges"})

    def finite(self, name, K):
        """The Kretschmann scalar at points on or approaching a line drawn as regular."""
        K = np.abs(np.asarray(K, dtype=float))
        ok = bool(np.all(np.isfinite(K)) and np.all(K < 1e4))
        self.limits.append({"name": name, "error": float(np.max(K)), "ok": ok, "kind": "finite"})

    def failures(self):
        return [c["name"] for c in self.charts if not c["ok"]] + [x["name"] for x in self.limits if not x["ok"]]

    def report(self):
        print(f"{'chart':52s} {'resolved':>9s} {'inverse':>9s} {'null':>9s}  distinct  future")
        for c in self.charts:
            print(f"{c['name']:52s} {c['resolved']:5d}/{c['samples']:<4d}{c['inverse']:9.1e} {c['null']:9.1e}"
                  f"  {'yes' if c['distinct'] else 'NO':8s}  {'yes' if c['future'] else 'NO':6s}"
                  f"  {'ok' if c['ok'] else 'FAILED'}")
        print()
        for x in self.limits:
            shown = {"diverges": "K", "finite": "K"}.get(x.get("kind"), "err")
            print(f"{'ok    ' if x['ok'] else 'FAILED'} {x['name']}  ({shown} = {x['error']:.2e})")


# ---------------------------------------------------------------- shared pictures

TRIANGLE = [[0, -PI], [PI, 0], [0, PI]]
DIAMOND = [[PI, 0], [0, PI], [-PI, 0], [0, -PI]]


def triangle_edges(v, centre="$r = 0$", centre_class="centre"):
    """The half diamond of a static, asymptotically flat spacetime, each point a sphere,
    drawn with p, q = arctan((t -+ r*)/l): the centre on X = 0, i0 at (pi, 0), i+- at
    (0, +-pi)."""
    v.line(centre_class, [[[0, -PI], [0, PI]]])
    v.line("scri", [[[0, PI], [PI, 0]], [[PI, 0], [0, -PI]]])
    for at in ((PI, 0), (0, PI), (0, -PI)):
        v.label_xt(at, {(PI, 0): "$i^0$", (0, PI): "$i^+$", (0, -PI): "$i^-$"}[at],
                   {(PI, 0): "l", (0, PI): "b", (0, -PI): "t"}[at],
                   dx=6 if at == (PI, 0) else 0, dy={(0, PI): -6, (0, -PI): 6}.get(at, 0))
        v.layers.append({"kind": "point", "class": "infinity", "at": [round(at[0], 4), round(at[1], 4)]})
    v.label_xt([HALF, HALF], "$\\mathscr{I}^+$", "bl", dx=5, dy=-3)
    v.label_xt([HALF, -HALF], "$\\mathscr{I}^-$", "tl", dx=5, dy=3)
    if centre:
        v.label_xt([0, 0.25], centre, "r", dx=-6)
    v.legend("scri", "null infinity $\\mathscr{I}^\\pm$")


def diamond_edges(v):
    """The full diamond of a flat plane, or of a wormhole's two sides."""
    v.line("scri", [[[PI, 0], [0, PI]], [[0, PI], [-PI, 0]], [[-PI, 0], [0, -PI]], [[0, -PI], [PI, 0]]])
    for at, text, anchor, dx, dy in (((PI, 0), "$i^0$", "l", 6, 0), ((-PI, 0), "$i^0$", "r", -6, 0),
                                     ((0, PI), "$i^+$", "b", 0, -6), ((0, -PI), "$i^-$", "t", 0, 6)):
        v.layers.append({"kind": "point", "class": "infinity", "at": [round(at[0], 4), round(at[1], 4)]})
        v.label_xt(at, text, anchor, dx=dx, dy=dy)
    for sx in (1, -1):
        v.label_xt([sx * HALF, HALF], "$\\mathscr{I}^+$", "bl" if sx > 0 else "br", dx=5 * sx, dy=-3)
        v.label_xt([sx * HALF, -HALF], "$\\mathscr{I}^-$", "tl" if sx > 0 else "tr", dx=5 * sx, dy=3)
    v.legend("scri", "null infinity $\\mathscr{I}^\\pm$")


def grid(v, cls, fmap, constants, s, first=True):
    """Lines of one coordinate held at each constant, the other running over s:
    fmap(constant, s) if first, fmap(s, constant) otherwise."""
    for c in constants:
        cs = np.full_like(s, c)
        v.curve(cls, *(fmap(cs, s) if first else fmap(s, cs)))


def label_on(v, pq, text, anchor="b", cls="coord", dx=0, dy=-3):
    v.label((float(pq[0]), float(pq[1])), text, anchor, cls, dx, dy)


S_ALL = spread(-np.inf, np.inf, 500, 9)
S_POS = spread(0, np.inf, 500, 12)


# ---------------------------------------------------------------- Minkowski

def minkowski(ck, src):
    """Minkowski's triangle and diamond.

    With u, v = ct -+ r the metric on the plane of t and r is -du dv, and
    p = arctan(u/l), q = arctan(v/l) give -du dv = -l^2 sec^2 p sec^2 q dp dq, conformal to
    -dT^2 + dX^2 with T = p + q, X = q - p, for any length l, drawn here as 1. r >= 0 is the
    triangle X >= 0 below the lines X + |T| = pi, and the Cartesian plane, x over the whole
    line, is the full diamond. Rindler's ct = X sinh(aT/c), x = X cosh(aT/c) make
    ct -+ x = -+X exp(-+aT/c), so the wedge x > c|t| is p < 0 < q.
    """
    views = []
    tri_box = [-0.35, PI + 0.35, -PI - 0.25, PI + 0.25]
    dia_box = [-PI - 0.35, PI + 0.35, -PI - 0.25, PI + 0.25]
    TS, RS, XS = (-4, -2, -1, 0, 1, 2, 4), (0.5, 1, 2, 4), (-4, -2, -1, 0, 1, 2, 4)

    sph = Plane(src, "minkowski", "spherical", ("t", "r"), EQUATOR)
    ck.chart("Minkowski spherical", sph, mink_pq, ck.uniform(-20, 20), ck.uniform(0.01, 20), lambda t, r: (1, 0))
    ck.finite("Minkowski: r = 0 is a regular centre", sph.kretschmann(ck.uniform(-5, 5, 50), np.full(50, 1e-6)))
    v = View("spherical", "Spherical", tri_box, "spherical")
    v.fill("region", TRIANGLE)
    v.fill("cover", TRIANGLE)
    grid(v, "r", lambda r, t: mink_pq(t, r), RS, S_ALL)
    grid(v, "t", mink_pq, TS, S_POS)
    triangle_edges(v)
    label_on(v, mink_pq(0, 1), "$r = 1$")
    label_on(v, mink_pq(0, 4), "$4$")
    v.legend("cover", "the whole spacetime, which $t$ and $r$ cover")
    v.legend("r", "$r$ constant, in units of $\\ell$")
    v.legend("t", "$ct$ constant")
    v.legend("centre", "$r = 0$, a regular centre")
    views.append(v)

    null = Plane(src, "minkowski", "spherical_null", ("u", "v"), EQUATOR)

    def null_map(u, w):
        return np.arctan(np.asarray(u, dtype=float)), np.arctan(np.asarray(w, dtype=float))
    u = ck.uniform(-20, 20)
    ck.chart("Minkowski spherical null", null, null_map, u, u + ck.uniform(0.01, 30), lambda u, w: (1, 1))
    v = View("spherical_null", "Spherical null", tri_box, "spherical_null")
    v.fill("region", TRIANGLE)
    v.fill("cover", TRIANGLE)
    for c in TS:
        s = np.linspace(c, 60, 400)
        v.curve("null", *null_map(np.full_like(s, c), s))
        s = np.linspace(-60, c, 400)
        v.curve("null", *null_map(s, np.full_like(s, c)))
    triangle_edges(v)
    v.legend("cover", "the whole spacetime, which $u$ and $v$ cover")
    v.legend("null", "$u$ constant and $v$ constant, every one a light ray")
    v.legend("centre", "$r = 0$, where $u = v$")
    views.append(v)

    cart = Plane(src, "minkowski", "cartesian", ("t", "x"), {"y": "0", "z": "0"})
    ck.chart("Minkowski Cartesian", cart, mink_pq, ck.uniform(-20, 20), ck.uniform(-20, 20), lambda t, x: (1, 0))
    v = View("cartesian", "Cartesian", dia_box, "cartesian")
    v.fill("region", DIAMOND)
    v.fill("cover", DIAMOND)
    grid(v, "r", lambda x, t: mink_pq(t, x), XS, S_ALL)
    grid(v, "t", mink_pq, TS, S_ALL)
    diamond_edges(v)
    v.legend("r", "$x$ constant")
    v.legend("t", "$ct$ constant")
    v.set(restriction="The plane $y = z = 0$ only, which is totally geodesic; each point of the "
                      "diagram is a single event rather than a sphere of them.")
    views.append(v)

    dn = Plane(src, "minkowski", "double_null", ("u", "v"), {"y": "0", "z": "0"})
    ck.chart("Minkowski double null", dn, null_map, ck.uniform(-20, 20), ck.uniform(-20, 20), lambda u, w: (1, 1))
    v = View("double_null", "Double null", dia_box, "double_null")
    v.fill("region", DIAMOND)
    v.fill("cover", DIAMOND)
    s = np.linspace(-80, 80, 600)
    for c in XS:
        v.curve("null", *null_map(np.full_like(s, c), s))
        v.curve("null", *null_map(s, np.full_like(s, c)))
    diamond_edges(v)
    v.legend("null", "$u$ constant and $v$ constant, every one a light ray")
    v.set(restriction="The plane $y = z = 0$ only, which is totally geodesic; each point of the "
                      "diagram is a single event rather than a sphere of them.")
    views.append(v)

    rind = Plane(src, "minkowski", "rindler", ("T", "X"), {"Y": "0", "Z": "0"}, {"a": 1})

    def rindler(T, X):
        T, X = np.asarray(T, dtype=float), np.asarray(X, dtype=float)
        return np.arctan(-X * np.exp(-T)), np.arctan(X * np.exp(T))
    ck.chart("Minkowski Rindler", rind, rindler, ck.uniform(-5, 5), ck.uniform(0.01, 10), lambda T, X: (1, 0))
    p, q = rindler(np.array([-2.0, 2.0]), np.array([1e-12, 1e-12]))
    ck.limit("Rindler: X -> 0 is the pair of null lines p = 0 and q = 0", [q[0], p[1]], [0, 0], 1e-9)
    v = View("rindler", "Rindler", dia_box, "rindler")
    v.fill("region", DIAMOND)
    v.fill("cover", [[0, 0], [HALF, -HALF], [PI, 0], [HALF, HALF]])
    grid(v, "r", lambda X, T: rindler(T, X), (0.25, 0.5, 1, 2, 4), S_ALL)
    grid(v, "t", rindler, (-2, -1, -0.5, 0, 0.5, 1, 2), S_POS)
    v.line("horizon", [[[0, 0], [HALF, HALF]], [[0, 0], [HALF, -HALF]],
                       [[0, 0], [-HALF, HALF]], [[0, 0], [-HALF, -HALF]]])
    diamond_edges(v)
    v.label_xt([Q4 - 0.05, Q4 + 0.05], "$X = 0$", "br", "small", dx=-3, dy=-2)
    v.legend("cover", "the wedge $x > c|t|$, which $T$ and $X$ cover")
    v.legend("r", "$X$ constant, a uniformly accelerated observer, in units of $c^2/a$")
    v.legend("t", "$T$ constant")
    v.legend("horizon", "the horizon $X = 0$ and the null lines that continue it")
    v.set(restriction="The plane $Y = Z = 0$ only, which is totally geodesic; each point of the "
                      "diagram is a single event rather than a sphere of them.")
    views.append(v)
    return views


# ---------------------------------------------------------------- the tower

class Tower:
    """Static surfaces whose metric is -f dt^2 + dr^2/f with f having simple roots.

    The tortoise coordinate is r* = int dr/f, split into partial fractions,

        r* = r + sum_i A_i ln|r/r_i - 1|,      A_i = 1/f'(r_i),

    every logarithm's argument made dimensionless by its own root, so that r*(0) = 0; both
    identities are checked in sympy on construction, and |A_i| = 1/2k_i for the surface
    gravity k_i = |f'(r_i)|/2 of each horizon. With u = t - r* and v = t + r*, the Kruskal
    coordinate of the outer horizon is |U+| = exp(-k+ u), and every cell of the drawing is
    written in it,

        G(u) = arctan exp(-k+ u) = arctan |U+|,

    which makes the drawing smooth across every r+ horizon and continuous across every r-
    horizon, where it cannot also be smooth: the rays that reach null infinity late are the
    rays that pile up at the Cauchy horizon, so one function of the ray serves both, and the
    exterior is kept exactly as Kruskal draws it. Cells, as diamonds of side pi/2 in (p, q):

        I    exterior, right       p = -G(u),   q = G(-v)
        I'   exterior, left        p = G(u),    q = -G(-v)
        II   future interior       p = G(u),    q = G(-v)
        IV   past interior         p = -G(u),   q = -G(-v)
        III  inside r-, right      p = G(u),    q = pi - G(-v)
        III' inside r-, left       the mirror of III, p <-> q

    with t the static coordinate of each cell. The cells above the inner horizon's
    bifurcation are the reflection (p, q) -> (pi - q, pi - p), the t -> -t isometry of III,
    and the whole repeats with period pi in both p and q. With one root the same cells I, II,
    IV and I' are Kruskal and Szekeres's extension of Schwarzschild, with U = tan p and
    V = tan q.
    """

    def __init__(self, f, r, roots):
        self.f_sym = sp.simplify(f)
        self.roots = [sp.nsimplify(x) for x in roots]
        inv = sp.apart(sp.together(1 / self.f_sym), r)
        self.A = [sp.simplify(sp.limit((r - ri) / self.f_sym, r, ri)) for ri in self.roots]
        rest = sp.simplify(inv - sum(Ai / (r - ri) for Ai, ri in zip(self.A, self.roots)))
        assert rest == 1, f"1/f is not 1 plus simple poles: the remainder is {rest}"
        rstar = r + sum(Ai * sp.log(r / ri - 1) for Ai, ri in zip(self.A, self.roots))
        assert sp.simplify(sp.diff(rstar, r) - 1 / self.f_sym) == 0, "dr*/dr is not 1/f"
        fp = sp.diff(self.f_sym, r)
        self.kappa = [sp.simplify(sp.Abs(fp.subs(r, ri)) / 2) for ri in self.roots]
        for Ai, ki in zip(self.A, self.kappa):
            assert sp.simplify(sp.Abs(Ai) - 1 / (2 * ki)) == 0, "a residue is not 1/2k"
        self.kp = float(self.kappa[0])
        self.Af = [float(a) for a in self.A]
        self.rf = [float(x) for x in self.roots]

    def G(self, u):
        return atan_exp(-self.kp * np.asarray(u, dtype=float))

    def rstar(self, r):
        r = np.asarray(r, dtype=float)
        with np.errstate(divide="ignore"):
            return r + sum(A * np.log(np.abs(r / ri - 1)) for A, ri in zip(self.Af, self.rf))

    def pq(self, cell, t, r):
        """(p, q) of points (t, r) of a cell, t the static coordinate of its region."""
        t = np.asarray(t, dtype=float)
        rs = self.rstar(r)
        u, v = t - rs, t + rs
        G = self.G
        if cell == "I":
            return -G(u), G(-v)
        if cell == "II":
            return G(u), G(-v)
        if cell == "IV":
            return -G(u), -G(-v)
        if cell == "I'":
            return G(u), -G(-v)
        if cell == "III":
            return G(u), np.pi - G(-v)
        if cell == "III'":
            q, p = self.pq("III", t, r)
            return p, q
        raise KeyError(cell)


CELL_POLYGON = {
    "I": [(-HALF, 0), (-HALF, HALF), (0, HALF), (0, 0)],
    "I'": [(0, -HALF), (HALF, -HALF), (HALF, 0), (0, 0)],
    "II": [(0, 0), (HALF, 0), (HALF, HALF), (0, HALF)],
    "IV": [(-HALF, -HALF), (0, -HALF), (0, 0), (-HALF, 0)],
    "III": [(0, HALF), (HALF, HALF), (HALF, PI), (0, PI)],
    "III'": [(HALF, 0), (PI, 0), (PI, HALF), (HALF, HALF)],
}
CELL_RANGE = {"I": "out", "I'": "out", "II": "mid", "IV": "mid", "III": "in", "III'": "in"}
# One period of the tower and a little more: the cells as (cell, reflected above r-).
TOWER = [("IV", False), ("I", False), ("I'", False), ("II", False), ("III", False),
         ("III'", False), ("II", True), ("I", True), ("I'", True), ("IV", True)]


def reflect(p, q, up):
    p, q = np.asarray(p, dtype=float), np.asarray(q, dtype=float)
    return (PI - q, PI - p) if up else (p, q)


class TowerDrawing:
    """The cells of a Tower, drawn with their grids and edges."""

    def __init__(self, T, singular, rmin):
        self.T, self.singular, self.rmin = T, singular, rmin
        self.rp, self.rm = T.rf[0], T.rf[1]

    def r_range(self, cell):
        return {"out": (self.rp, np.inf), "mid": (self.rm, self.rp), "in": (self.rmin, self.rm)}[CELL_RANGE[cell]]

    def singularity(self, n=600):
        """r = 0 in cell III, which r*(0) = 0 puts on u = -v, the vertical line X = pi/2."""
        s = spread(-np.inf, np.inf, n, 9)
        return self.T.G(s), np.pi - self.T.G(-s)

    def polygon(self, cell, up=False):
        if cell in ("III", "III'") and self.singular:
            ps, qs = self.singularity()
            if cell == "III":
                pts = [(0, HALF), (HALF, HALF), (HALF, PI)] + list(zip(ps, qs))
            else:
                pts = [(HALF, 0), (HALF, HALF), (PI, HALF)] + list(zip(qs, ps))
        else:
            pts = CELL_POLYGON[cell]
        return [point(*reflect(p, q, up)) for p, q in pts]

    def curve(self, v, cls, cell, t, r, up=False):
        v.curve(cls, *reflect(*self.T.pq(cell, t, r), up))

    def grid(self, v, cell, radii, times, up=False):
        t = spread(-np.inf, np.inf, 500, 10)
        for r in radii:
            self.curve(v, "r", cell, t, np.full_like(t, r), up)
        lo, hi = self.r_range(cell)
        rr = spread(lo, hi, 600, 16)
        for tt in times:
            self.curve(v, "t", cell, np.full_like(rr, tt), rr, up)

    def edges(self, v, cell, up=False):
        """The horizons between cells, null infinity, and r = 0."""
        def seg(cls, a, b, zig=False):
            v.segment(cls, reflect(*a, up), reflect(*b, up), zig)
        if cell == "I":
            seg("scri", (-HALF, 0), (-HALF, HALF))
            seg("scri", (-HALF, HALF), (0, HALF))
            seg("horizon", (0, 0), (0, HALF))
            seg("horizon", (-HALF, 0), (0, 0))
        elif cell == "I'":
            seg("scri", (0, -HALF), (HALF, -HALF))
            seg("scri", (HALF, -HALF), (HALF, 0))
            seg("horizon", (0, 0), (HALF, 0))
            seg("horizon", (0, -HALF), (0, 0))
        elif cell in ("II", "IV"):
            sign = 1 if cell == "II" else -1
            for a, b in (((0, 0), (HALF, 0)), ((HALF, 0), (HALF, HALF)), ((HALF, HALF), (0, HALF)), ((0, HALF), (0, 0))):
                seg("horizon", (sign * a[0], sign * a[1]), (sign * b[0], sign * b[1]))
        elif self.singular:
            ps, qs = self.singularity()
            if cell == "III'":
                ps, qs = qs, ps
            v.curve("singular", *reflect(ps, qs, up), zig=True, tol=0.01)
        elif cell == "III":
            seg("scri", (0, HALF), (0, PI))
            seg("scri", (0, PI), (HALF, PI))
        else:
            seg("scri", (HALF, 0), (PI, 0))
            seg("scri", (PI, 0), (PI, HALF))

    def draw(self, v, grids, cover=()):
        for cell, up in TOWER:
            v.fill("region", self.polygon(cell, up))
        for cell, up in cover:
            v.fill("cover", self.polygon(cell, up))
        for cell, up in TOWER:
            key = cell.rstrip("'")
            if key in grids:
                self.grid(v, cell, *grids[key], up)
        for cell, up in TOWER:
            self.edges(v, cell, up)

    def labels(self, v, inner="$r < r_-$"):
        for sx in (1, -1):
            for base in (0, 2 * PI):
                v.layers.append({"kind": "point", "class": "infinity", "at": [round(sx * PI, 4), round(base, 4)]})
                v.label_xt([sx * PI, base], "$i^0$", "l" if sx > 0 else "r", dx=6 * sx)
                for dT, text, anchor, dy in ((HALF, "$i^+$", "b", -2), (-HALF, "$i^-$", "t", 2)):
                    v.layers.append({"kind": "point", "class": "infinity", "at": [round(sx * HALF, 4), round(base + dT, 4)]})
                    v.label_xt([sx * HALF, base + dT], text, anchor + ("l" if sx > 0 else "r"), dx=5 * sx, dy=dy)
                v.label_xt([sx * 3 * Q4, base + Q4], "$\\mathscr{I}^+$", "bl" if sx > 0 else "br", dx=4 * sx, dy=-3)
                v.label_xt([sx * 3 * Q4, base - Q4], "$\\mathscr{I}^-$", "tl" if sx > 0 else "tr", dx=4 * sx, dy=3)
                v.label_xt([sx * HALF, base], "exterior", cls="region")
        v.label_xt([0, HALF + 0.1], "black hole", cls="region")
        v.label_xt([0, 1.5 * PI - 0.1], "white hole", cls="region")
        v.label_xt([0, -HALF], "white hole", cls="region")
        v.label_xt([0, 2.5 * PI], "black hole", cls="region")
        v.label_xt([Q4, Q4], "$r_+$", "tl", "small", dx=5, dy=1)
        v.label_xt([Q4, 3 * Q4], "$r_-$", "bl", "small", dx=5, dy=-1)
        for sx in (1, -1):
            if self.singular:
                v.label_xt([sx * HALF, PI + 0.35], "$r = 0$", "l" if sx > 0 else "r", dx=8 * sx)
                v.label_xt([sx * 1.05, PI - 0.2], inner, cls="region")
            else:
                v.layers.append({"kind": "point", "class": "infinity", "at": [round(sx * PI, 4), round(PI, 4)]})
                v.label_xt([sx * PI, PI], "$i^0$", "l" if sx > 0 else "r", dx=6 * sx)
                v.label_xt([sx * 3 * Q4, PI + Q4], "$\\mathscr{I}^+$", "bl" if sx > 0 else "br", dx=4 * sx, dy=-3)
                v.label_xt([sx * 3 * Q4, PI - Q4], "$\\mathscr{I}^-$", "tl" if sx > 0 else "tr", dx=4 * sx, dy=3)
                v.label_xt([sx * 2.2, PI], inner, cls="region")


def nice(x, keep=()):
    """The fewest significant figures, two at least, that keep x clear of the values in keep."""
    if x == 0:
        return 0.0
    for sig in range(2, 16):
        y = round(x, -(int(np.floor(np.log10(abs(x)))) - sig + 1))
        gap = min([abs(x - w) for w in keep] or [abs(x)])
        if abs(y - x) < 0.2 * gap:
            return y
    return x


def nice_all(xs, bounds):
    return [nice(x, list(bounds) + [w for j, w in enumerate(xs) if j != i]) for i, x in enumerate(xs)]


def invert_rstar(T, target, lo, hi):
    """r in (lo, hi) with r*(r) = target, by bisection, r* being monotone in each cell."""
    a, b = lo, hi
    if np.isinf(b):
        b = max(lo, 1.0) * 2
        while T.rstar(b) < target:
            b *= 2
    if np.isinf(a):
        a = -1.0
        while T.rstar(a) > target:
            a *= 2
    fa = T.rstar(a + 1e-15 * max(1, abs(a))) - target
    for _ in range(200):
        m = (a + b) / 2
        fm = T.rstar(m) - target
        if (fm > 0) == (fa > 0):
            a, fa = m, fm
        else:
            b = m
    return (a + b) / 2


def even_radii(T, cell, n, lo, hi, xmax=PI):
    """Radii whose lines cross the middle of a cell evenly: in I, X at T = 0 is
    2 arctan e^(k r*) from 0 to pi; in II, T at X = 0 is the same; in III, X at T = pi is
    pi - 2 arctan e^(k r*), across (0, xmax)."""
    out = []
    for i in range(1, n + 1):
        f = i / (n + 1)
        if cell in ("I", "II"):
            target = np.log(np.tan(f * HALF)) / T.kp
        else:
            target = np.log(np.tan((PI - xmax * (1 - f)) / 2)) / T.kp
        out.append(invert_rstar(T, target, lo, hi))
    return out


def listed(xs):
    return ", ".join("%.15g" % x for x in xs)


def tower_checks(ck, name, plane, T, rmin, span):
    """The tower's cells against the published metric, with t the static coordinate of each."""
    rp, rm = T.rf
    ck.chart(f"{name}, exterior", plane, lambda t, r: T.pq("I", t, r),
             ck.uniform(-span, span), ck.uniform(rp + 1e-3, 40), lambda t, r: (1, 0))
    ck.chart(f"{name}, black hole", plane, lambda t, r: T.pq("II", t, r),
             ck.uniform(-span, span), ck.uniform(rm + 1e-3, rp - 1e-3), lambda t, r: (0, -1))
    ck.chart(f"{name}, inside r-", plane, lambda t, r: T.pq("III", t, r),
             ck.uniform(-span, span), ck.uniform(rmin, rm - 1e-3), lambda t, r: (-1, 0))


# ---------------------------------------------------------------- Schwarzschild

def schwarzschild(ck, src):
    """Kruskal and Szekeres's extension, compactified by p = arctan U, q = arctan V.

    U = -exp(-u/2r_s) and V = exp(v/2r_s), u, v = ct -+ r*, r* = r + r_s ln|r/r_s - 1|, give
    UV = (1 - r/r_s) exp(r/r_s), analytic through r_s: the tower of one root, cells I, II,
    IV and I'. Since tan(p + q) = (U + V)/(1 - UV), the singularity UV = 1 is exactly the
    pair of straight lines T = +-pi/2. The ingoing coordinates are V = exp(v/2r_s),
    U = (1 - r/r_s) exp(r/r_s)/V, one formula for every r > 0 covering I and II; the
    outgoing ones their time reverse, U = -exp(-u/2r_s), V = (r/r_s - 1) exp(r/r_s)/(-U),
    covering I and IV.
    """
    sph = Plane(src, "schwarzschild", "spherical", ("t", "r"), EQUATOR, {"r_s": 1})
    assert sph.g[0, 1] == 0 and sp.simplify(sph.g[0, 0] * sph.g[1, 1] + 1) == 0
    T = Tower(-sph.g[0, 0], sph.x1, [1])
    ck.chart("Schwarzschild spherical, exterior", sph, lambda t, r: T.pq("I", t, r),
             ck.uniform(-15, 15), ck.uniform(1.001, 30), lambda t, r: (1, 0))
    ck.chart("Schwarzschild spherical, black hole", sph, lambda t, r: T.pq("II", t, r),
             ck.uniform(-15, 15), ck.uniform(0.01, 0.999), lambda t, r: (0, -1))
    ck.chart("Schwarzschild spherical, white hole", sph, lambda t, r: T.pq("IV", t, r),
             ck.uniform(-15, 15), ck.uniform(0.01, 0.999), lambda t, r: (0, 1))
    ck.chart("Schwarzschild spherical, other exterior", sph, lambda t, r: T.pq("I'", t, r),
             ck.uniform(-15, 15), ck.uniform(1.001, 30), lambda t, r: (-1, 0))

    def ingoing(w, r):
        w, r = np.asarray(w, dtype=float), np.asarray(r, dtype=float)
        return np.arctan((1 - r) * np.exp(r - w / 2)), atan_exp(w / 2)

    def outgoing(u, r):
        u, r = np.asarray(u, dtype=float), np.asarray(r, dtype=float)
        return -atan_exp(-u / 2), np.arctan((r - 1) * np.exp(r + u / 2))
    ein = Plane(src, "schwarzschild", "eddington_finkelstein_ingoing", ("v", "r"), EQUATOR, {"r_s": 1})
    ck.chart("Schwarzschild ingoing Eddington-Finkelstein", ein, ingoing,
             ck.uniform(-15, 15), ck.uniform(0.01, 30), lambda w, r: (1, -60))
    eout = Plane(src, "schwarzschild", "eddington_finkelstein_outgoing", ("u", "r"), EQUATOR, {"r_s": 1})
    ck.chart("Schwarzschild outgoing Eddington-Finkelstein", eout, outgoing,
             ck.uniform(-15, 15), ck.uniform(0.01, 30), lambda u, r: (1, 60))

    p, q = T.pq("II", np.array([-5.0, 0, 5]), np.full(3, 1e-9))
    ck.limit("Schwarzschild: r -> 0 in the black hole lands on T = pi/2", p + q, [HALF] * 3)
    p, q = T.pq("I", np.array([0.0]), np.array([1e8]))
    ck.limit("Schwarzschild: r -> infinity at t = 0 lands on i0, (X, T) = (pi, 0)", point(p[0], q[0]), [PI, 0], 1e-3)
    p, q = T.pq("I", np.array([3.0]), np.array([1 + 1e-12]))
    ck.limit("Schwarzschild: r -> r_s at fixed t lands on the bifurcation sphere", point(p[0], q[0]), [0, 0], 1e-4)
    rr = np.linspace(0.05, 5, 50)
    for cell, sel in (("I", rr > 1.001), ("II", rr < 0.999)):
        pp, qq = T.pq(cell, 0.3 + 0 * rr[sel], rr[sel])
        ck.limit(f"Schwarzschild {cell}: tan p tan q is Kruskal's UV = (1 - r/r_s) e^(r/r_s)",
                 np.tan(pp) * np.tan(qq), (1 - rr[sel]) * np.exp(rr[sel]), 1e-8)
    ck.limit("Schwarzschild: the ingoing and spherical coordinates put one event at one point",
             ingoing(2.0 + 3 + np.log(2), 3.0), T.pq("I", 2.0, 3.0), 1e-12)
    K = sph.kretschmann
    ck.diverges("Schwarzschild: the Kretschmann scalar diverges at r = 0", K(0, 1e-2), K(0, 1e-3))
    ck.finite("Schwarzschild: the Kretschmann scalar is finite at r = r_s", K(np.zeros(3), np.array([0.999, 1, 1.001])))

    box = [-PI - 0.25, PI + 0.25, -HALF - 0.25, HALF + 0.25]
    hexagon = [[PI, 0], [HALF, HALF], [-HALF, HALF], [-PI, 0], [-HALF, -HALF], [HALF, -HALF]]
    exterior = [[0, 0], [HALF, -HALF], [PI, 0], [HALF, HALF]]
    R_OUT, R_IN, TS = (1.05, 1.25, 1.5, 2, 3), (0.5, 0.75, 0.9), (-4, -2, -1, 0, 1, 2, 4)

    def edges(v):
        v.line("scri", [[[PI, 0], [HALF, HALF]], [[PI, 0], [HALF, -HALF]],
                        [[-PI, 0], [-HALF, HALF]], [[-PI, 0], [-HALF, -HALF]]])
        v.line("horizon", [[[-HALF, -HALF], [HALF, HALF]], [[HALF, -HALF], [-HALF, HALF]]])
        v.line("singular", [[[-HALF, HALF], [HALF, HALF]], [[-HALF, -HALF], [HALF, -HALF]]], zig=True)
        for at in ((PI, 0), (-PI, 0), (HALF, HALF), (HALF, -HALF), (-HALF, HALF), (-HALF, -HALF)):
            v.layers.append({"kind": "point", "class": "infinity", "at": rounded(at)})
        v.label_xt([PI, 0], "$i^0$", "l", dx=6)
        v.label_xt([-PI, 0], "$i^0$", "r", dx=-6)
        for sx in (1, -1):
            v.label_xt([sx * HALF, HALF], "$i^+$", "b", dy=-6)
            v.label_xt([sx * HALF, -HALF], "$i^-$", "t", dy=6)
            v.label_xt([sx * 3 * Q4, Q4], "$\\mathscr{I}^+$", "bl" if sx > 0 else "br", dx=4 * sx, dy=-4)
            v.label_xt([sx * 3 * Q4, -Q4], "$\\mathscr{I}^-$", "tl" if sx > 0 else "tr", dx=4 * sx, dy=4)
        v.label_xt([0, HALF], "$r = 0$", "b", dy=-8)
        v.label_xt([0, -HALF], "$r = 0$", "t", dy=8)
        v.label_xt([-Q4, Q4], "$r = r_s$", "tr", "small", dx=-6, dy=2)
        v.label_xt([HALF, -0.95], "exterior", cls="region")
        v.label_xt([-HALF, 0], "exterior", cls="region")
        v.label_xt([0, 1.15], "black hole", cls="region")
        v.label_xt([0, -1.15], "white hole", cls="region")
        v.legend("horizon", "the horizon $r = r_s$")
        v.legend("singular", "$r = 0$, where the Kretschmann scalar diverges")
        v.legend("scri", "null infinity $\\mathscr{I}^\\pm$")

    views = []
    t = spread(-np.inf, np.inf, 500, 9)
    v = View("spherical", "Spherical", box, "spherical")
    v.fill("region", hexagon)
    v.fill("cover", exterior)
    for r in R_OUT:
        v.curve("r", *T.pq("I", t, np.full_like(t, r)))
    rr = spread(1, np.inf, 500, 14)
    for tt in TS:
        v.curve("t", *T.pq("I", np.full_like(rr, tt), rr))
    edges(v)
    for r, text in ((1.25, "$1.25\\,r_s$"), (2, "$2\\,r_s$")):
        label_on(v, T.pq("I", 0.0, r), text)
    v.legend("cover", "the region that $t$ and $r > r_s$ cover")
    v.legend("r", "$r$ constant")
    v.legend("t", "$ct$ constant, in units of $r_s$")
    views.append(v)

    v = View("ingoing", "Ingoing Eddington-Finkelstein", box, "eddington_finkelstein_ingoing")
    v.fill("region", hexagon)
    v.fill("cover", [[0, 0], [HALF, -HALF], [PI, 0], [HALF, HALF], [-HALF, HALF]])
    for r in R_OUT + R_IN:
        v.curve("r", *ingoing(t, np.full_like(t, r)))
    rr = spread(0, np.inf, 600, 14)
    for w in (-6, -4, -2, 0, 2, 4, 6):
        v.curve("null", *ingoing(np.full_like(rr, w), rr))
    edges(v)
    v.legend("cover", "the region that $v$ and $r > 0$ cover")
    v.legend("r", "$r$ constant")
    v.legend("null", "$v$ constant, an ingoing light ray")
    views.append(v)

    v = View("outgoing", "Outgoing Eddington-Finkelstein", box, "eddington_finkelstein_outgoing")
    v.fill("region", hexagon)
    v.fill("cover", [[0, 0], [-HALF, -HALF], [HALF, -HALF], [PI, 0], [HALF, HALF]])
    for r in R_OUT + R_IN:
        v.curve("r", *outgoing(t, np.full_like(t, r)))
    for u in (-6, -4, -2, 0, 2, 4, 6):
        v.curve("null", *outgoing(np.full_like(rr, u), rr))
    edges(v)
    v.legend("cover", "the region that $u$ and $r > 0$ cover")
    v.legend("r", "$r$ constant")
    v.legend("null", "$u$ constant, an outgoing light ray")
    views.append(v)
    return views


# ---------------------------------------------------------------- Reissner-Nordstrom, and the axis of Kerr

def reissner_nordstrom(ck, src):
    """The tower of Reissner-Nordstrom, horizons at the roots of the published g^rr.

    Drawn at r_q = 0.48 r_s, so that k-/k+ = 3.2: at r_q = 0.4 r_s the ratio is 16 and a
    single Kruskal coordinate squeezes everything inside r- to within 1e-6 of r = 0. Every
    event beyond the Cauchy horizon has the whole exterior I in its causal past, since
    p <= p_e and q <= q_e hold at every point of it and both grow along every future causal
    curve: the Malament-Hogarth property on a member of the family that has it.
    """
    pl = Plane(src, "rn_metric", "spherical", ("t", "r"), EQUATOR, {"r_s": 1, "r_q": "12/25"})
    assert pl.g[0, 1] == 0 and sp.simplify(pl.g[0, 0] * pl.g[1, 1] + 1) == 0
    roots = sorted(sp.solve(sp.numer(sp.together(pl.gi[1, 1])), pl.x1), reverse=True)
    T = Tower(-pl.g[0, 0], pl.x1, roots)
    rp, rm = T.rf
    tower_checks(ck, "Reissner-Nordstrom", pl, T, 0.01, 15)
    p, q = T.pq("III", np.array([-6.0, 0, 6]), np.full(3, 1e-12))
    ck.limit("Reissner-Nordstrom: r -> 0 lands on the vertical line X = pi/2", q - p, [HALF] * 3, 1e-9)
    ck.limit("Reissner-Nordstrom: the horizons are the roots of the published g^rr", [rp, rm], [0.64, 0.36], 1e-12)
    ck.diverges("Reissner-Nordstrom: the Kretschmann scalar diverges at r = 0",
                pl.kretschmann(0, 1e-2), pl.kretschmann(0, 1e-3))
    ck.finite("Reissner-Nordstrom: the Kretschmann scalar is finite at both horizons",
              pl.kretschmann(np.zeros(2), np.array([rp, rm])))

    D = TowerDrawing(T, True, 0.0)
    box = [-PI - 0.45, PI + 0.45, -PI - 0.1, 3 * PI + 0.1]
    times = [c / T.kp for c in (-1.6, -0.6, 0, 0.6, 1.6)]
    rI = nice_all(even_radii(T, "I", 4, rp, np.inf), [rp])
    rII = nice_all(even_radii(T, "II", 4, rm, rp), [rm, rp])
    rIII = [0.3] + nice_all(even_radii(T, "III", 2, 0, rm, xmax=HALF), [0, rm, 0.3])
    grids = {"I": (rI, times), "II": (rII, times), "IV": (rII, times), "III": (rIII, times)}
    views = []
    v = View("tower", "Spherical", box, "spherical")
    D.draw(v, grids, cover=[("I", False), ("II", False), ("III", False)])
    D.labels(v)
    v.set(fade={"top": 0.9, "bottom": 0.9})
    v.legend("cover", "one exterior, one region between the horizons and one inside $r_-$, which "
                      "$t$ and $r > 0$ cover")
    v.legend("r", f"$r$ constant: {listed(rI)} outside, {listed(rII)} between, {listed(rIII)} inside $r_-$, "
                  "in units of $r_s$")
    v.legend("t", "$t$ constant")
    v.legend("horizon", f"the horizons $r_+ = {rp:g}\\,r_s$ and $r_- = {rm:g}\\,r_s$")
    v.legend("singular", "$r = 0$, a timelike singularity, where the Kretschmann scalar diverges")
    v.legend("scri", "null infinity $\\mathscr{I}^\\pm$")
    views.append(v)

    v = View("malament_hogarth", "Malament-Hogarth", box)
    for cell, up in TOWER:
        v.fill("region", D.polygon(cell, up))
    v.fill("past", D.polygon("I"))
    for cell, up in TOWER:
        D.edges(v, cell, up)
    t = spread(-np.inf, np.inf, 500, 10)
    D.curve(v, "world", "I", t, np.full_like(t, 1.2))
    pe, qe = 0.55, HALF + 0.3
    ps, qs = D.singularity()
    k = np.argmin(np.abs(qs - qe))
    assert qe - pe < qs[k] - ps[k], "the event is not inside r = 0"
    v.segment("cone", (pe, qe), (pe, -HALF))
    v.segment("cone", (pe, qe), (ps[k], qe))
    v.point("mark", (pe, qe))
    D.labels(v)
    v.label((pe, qe), "an event beyond $r_-$", "l", "small", dx=7)
    v.label(T.pq("I", 0.0, 1.2), "observer at $r = 1.2\\,r_s$", "r", "small", dx=-6)
    v.set(fade={"top": 0.9, "bottom": 0.9})
    v.legend("past", "the whole exterior, in the event's causal past")
    v.legend("world", "a static observer, whose proper time is infinite")
    v.legend("cone", "the edges of the event's past light cone")
    v.legend("mark", "an event beyond the Cauchy horizon $r_-$")
    v.legend("horizon", f"the horizons $r_+ = {rp:g}\\,r_s$ and $r_- = {rm:g}\\,r_s$")
    v.legend("singular", "$r = 0$, where the Kretschmann scalar diverges")
    views.append(v)
    settings = ("$r_q = 0.48\\,r_s$, so that $r_+ = 0.64\\,r_s$, $r_- = 0.36\\,r_s$ and "
                "$\\kappa_-/\\kappa_+ = 3.2$; at $r_q = 0.4\\,r_s$ the ratio is 16 and every line "
                "inside $r_-$ would lie within $10^{-6}$ of the singularity.")
    for view in views:
        view.set(settings=settings)
    return views


def kerr_axis(ck, src, metric_id, params, name):
    """The symmetry axis theta = 0, where the published metric is, in the limit,
    -Delta/(r^2 + a^2) c^2dt^2 + (r^2 + a^2)/Delta dr^2, with g_tt g_rr = -1 checked: a tower
    with two simple roots and no singularity, r running through the ring's disc to a second
    asymptotically flat end at r -> -infinity. That is Carter's diagram of 1966.
    """
    pl = Plane(src, metric_id, "boyer_lindquist", ("t", "r"), {"phi": "0"}, params, axis=("theta", "0"))
    assert pl.g[0, 1] == 0 and sp.simplify(pl.g[0, 0] * pl.g[1, 1] + 1) == 0
    roots = sorted([x for x in sp.solve(sp.numer(sp.together(pl.gi[1, 1])), pl.x1) if x.is_real], reverse=True)
    T = Tower(-pl.g[0, 0], pl.x1, roots)
    rp, rm = T.rf
    tower_checks(ck, f"{name} axis", pl, T, -40, 30)
    p, q = T.pq("III", np.array([0.0]), np.array([-1e9]))
    ck.limit(f"{name} axis: r -> -infinity at t = 0 lands on the far i0, (X, T) = (pi, pi)",
             point(p[0], q[0]), [PI, PI], 1e-3)
    ck.finite(f"{name} axis: the Kretschmann scalar is finite at r = 0 on the axis",
              pl.kretschmann(np.zeros(3), np.array([-1e-3, 0.0, 1e-3])))

    D = TowerDrawing(T, False, -np.inf)
    box = [-PI - 0.45, PI + 0.45, -PI - 0.1, 3 * PI + 0.1]
    times = [c / T.kp for c in (-1.6, -0.6, 0, 0.6, 1.6)]
    rI = nice_all(even_radii(T, "I", 4, rp, np.inf), [rp])
    rII = nice_all(even_radii(T, "II", 4, rm, rp), [rm, rp])
    rIII = [x for x in nice_all(even_radii(T, "III", 5, -np.inf, rm), [0, rm]) if abs(x) > 1e-9]
    v = View("axis", "Symmetry axis", box, "boyer_lindquist")
    D.draw(v, {"I": (rI, times), "II": (rII, times), "IV": (rII, times), "III": (rIII, times)},
           cover=[("I", False)])
    t = spread(-np.inf, np.inf, 500, 10)
    for cell in ("III", "III'"):
        D.curve(v, "centre", cell, t, np.zeros_like(t))
    D.labels(v)
    v.label_xt([HALF + 0.08, PI - 0.55], "$r = 0$", "l", "small", dx=4)
    v.set(fade={"top": 0.9, "bottom": 0.9},
          restriction="The symmetry axis $\\theta = 0$ only, a totally geodesic surface. The ring "
                      "singularity is at $r = 0$ in the equatorial plane $\\theta = \\pi/2$, off this "
                      "surface; on the axis $r$ runs through the centre of the ring's disc to $r < 0$.")
    v.legend("cover", "the exterior $r > r_+$, which $t$ and $r$ cover")
    v.legend("r", f"$r$ constant: {listed(rI)} outside, {listed(rII)} between, {listed(rIII)} inside $r_-$, "
                  "in units of $GM/c^2$")
    v.legend("t", "$ct$ constant")
    v.legend("horizon", f"the horizons on the axis, $r_+ = {rp:.3f}$ and $r_- = {rm:.3f}\\,GM/c^2$")
    v.legend("centre", "$r = 0$ on the axis, the centre of the ring's disc, where the curvature is finite")
    v.legend("scri", "null infinity, of $r \\to +\\infty$ and of $r \\to -\\infty$")
    return [v]


def kerr(ck, src):
    views = kerr_axis(ck, src, "kerr", {"G": 1, "M": 1, "a": "9/10"}, "Kerr")
    views[0].set(settings="$a = 0.9\\,GM/c^2$.")
    return views


def kerr_newman(ck, src):
    views = kerr_axis(ck, src, "kerr_newman", {"G": 1, "M": 1, "a": "3/5", "r_Q": "1/2"}, "Kerr-Newman")
    views[0].set(settings="$a = 0.6\\,GM/c^2$ and $r_Q = 0.5\\,GM/c^2$.")
    return views


# ---------------------------------------------------------------- de Sitter

def de_sitter(ck, src):
    """The global square, from the hyperboloid -X0^2 + X1^2 + ... + X4^2 = L^2, L = sqrt(3/Lambda),
    drawn at L = 1.

    Global coordinates X0 = L tan T, X4 = L cos(chi)/cos T and |X_1..3| = L sin(chi)/cos T give
    L^2/cos^2 T (-dT^2 + dchi^2 + sin^2 chi dOmega^2) on |T| < pi/2, 0 <= chi <= pi, so
    X = chi. The static coordinates, X0 = sqrt(L^2 - r^2) sinh(ct/L),
    X4 = sqrt(L^2 - r^2) cosh(ct/L), |X_1..3| = r, enter as tan p = tanh(u/2L),
    tan q = tanh(v/2L) with u, v = ct -+ L artanh(r/L). The flat slicing, conformal time
    eta = -exp(-Ht)/H and a = exp(Ht), enters as p = pi/4 + arctan(H eta - H rho/c),
    q = pi/4 + arctan(H eta + H rho/c). Both closed forms are checked against the embedding
    at random points, as well as against the published metric.
    """
    st = Plane(src, "de_sitter", "static_spherical", ("t", "r"), EQUATOR, {"Lambda": 3})
    fl = Plane(src, "de_sitter", "flat_slicing", ("t", "x"), {"y": "0", "z": "0"}, {"H": 1})

    def static(t, r):
        t, r = np.asarray(t, dtype=float), np.asarray(r, dtype=float)
        rs = np.arctanh(r)
        return np.arctan(np.tanh((t - rs) / 2)), np.arctan(np.tanh((t + rs) / 2))

    def flat(t, rho):
        eta = -np.exp(-np.asarray(t, dtype=float))
        return Q4 + np.arctan(eta - rho), Q4 + np.arctan(eta + rho)
    ck.chart("de Sitter static", st, static, ck.uniform(-10, 10), ck.uniform(0.001, 0.999), lambda t, r: (1, 0))
    ck.chart("de Sitter flat slicing", fl, flat, ck.uniform(-6, 6), ck.uniform(0.001, 20), lambda t, x: (1, 0))

    def embedding(p, q):
        T, chi = p + q, q - p
        return np.tan(T), np.cos(chi) / np.cos(T), np.sin(chi) / np.cos(T)
    t, r = ck.uniform(-6, 6, 2000), ck.uniform(0, 0.999, 2000)
    X0, X4, R = embedding(*static(t, r))
    scale = 1 + np.abs(X0) + np.abs(X4)
    ck.limit("de Sitter: the static coordinates land where the hyperboloid puts them",
             np.concatenate([(X0 - np.sqrt(1 - r ** 2) * np.sinh(t)) / scale,
                             (X4 - np.sqrt(1 - r ** 2) * np.cosh(t)) / scale, (R - r) / scale]), 0, 1e-10)
    t, rho = ck.uniform(-4, 4, 2000), ck.uniform(0, 8, 2000)
    X0, X4, R = embedding(*flat(t, rho))
    scale = 1 + np.abs(X0) + np.abs(X4) + np.abs(R)
    a = np.exp(t)
    ck.limit("de Sitter: the flat slicing lands where the hyperboloid puts it",
             np.concatenate([(X0 - np.sinh(t) - rho ** 2 * a / 2) / scale,
                             (X4 - np.cosh(t) + rho ** 2 * a / 2) / scale, (R - a * rho) / scale]), 0, 1e-10)
    p, q = static(np.array([0.0, 3, -3]), np.full(3, 1 - 1e-14))
    ck.limit("de Sitter: r -> sqrt(3/Lambda) lands on the horizons q = pi/4 (t > 0) and p = -pi/4 (t < 0)",
             [q[0], q[1], p[2]], [Q4, Q4, -Q4], 1e-6)
    p, q = flat(np.array([50.0]), np.array([1.0]))
    ck.limit("de Sitter: t -> infinity in the flat slicing lands on T = pi/2", p + q, [HALF], 1e-6)
    ck.finite("de Sitter: r = 0 is a regular centre", st.kretschmann(ck.uniform(-5, 5, 50), np.full(50, 1e-6)))

    box = [-0.35, PI + 0.35, -HALF - 0.3, HALF + 0.3]
    square = [[0, -HALF], [PI, -HALF], [PI, HALF], [0, HALF]]

    def frame(v):
        v.fill("region", square)
        v.line("scri", [[[0, HALF], [PI, HALF]], [[0, -HALF], [PI, -HALF]]])
        v.line("centre", [[[0, -HALF], [0, HALF]], [[PI, -HALF], [PI, HALF]]])
        v.line("horizon", [[[0, -HALF], [PI, HALF]], [[0, HALF], [PI, -HALF]]])
        v.label_xt([HALF, HALF], "$\\mathscr{I}^+$", "b", dy=-5)
        v.label_xt([HALF, -HALF], "$\\mathscr{I}^-$", "t", dy=5)
        v.label_xt([0, 0.55], "$r = 0$", "r", dx=-6)
        v.label_xt([PI, 0.55], "$r = 0$", "l", dx=6)
        v.label_xt([0, -1.25], "observer", "r", "coord", dx=-6)
        v.label_xt([PI, -1.25], "antipode", "l", "coord", dx=6)
        v.legend("scri", "future and past infinity $\\mathscr{I}^\\pm$, spacelike")
        v.legend("centre", "$r = 0$ of the observer and of its antipode")
        v.legend("horizon", "the cosmological horizons $r = \\sqrt{3/\\Lambda}$")

    views = []
    v = View("static", "Static spherical", box, "static_spherical")
    frame(v)
    v.fill("cover", [[0, -HALF], [HALF, 0], [0, HALF]])
    grid(v, "r", lambda r, t: static(t, r), (0.3, 0.6, 0.85, 0.97), S_ALL)
    grid(v, "t", static, (-3, -1.5, -0.5, 0, 0.5, 1.5, 3), spread(0, 1, 500, 14))
    label_on(v, static(0, 0.6), "$0.6$")
    label_on(v, static(0, 0.3), "$r = 0.3$")
    v.label_xt([PI * 0.8, -0.28], "the antipode's static patch", cls="region")
    v.label_xt([HALF, 1.05], "expanding", cls="region")
    v.label_xt([HALF, -1.05], "contracting", cls="region")
    v.legend("cover", "the static patch, which $t$ and $r$ cover")
    v.legend("r", "$r$ constant, in units of $\\sqrt{3/\\Lambda}$")
    v.legend("t", "$ct$ constant")
    views.append(v)

    v = View("flat", "Flat slicing", box, "flat_slicing")
    frame(v)
    v.fill("cover", [[0, -HALF], [PI, HALF], [0, HALF]])
    grid(v, "r", lambda rho, t: flat(t, rho), (0.25, 0.5, 1, 2, 4), S_ALL)
    grid(v, "t", flat, (-2, -1, 0, 1, 2), S_POS)
    v.line("chartedge", [[[0, -HALF], [PI, HALF]]])
    v.label_xt([PI * 0.62, 0.27], "$t \\to -\\infty$", "tl", "small", dx=4, dy=2)
    label_on(v, flat(0, 1), "$\\rho = c/H$")
    v.legend("cover", "the half that the flat slicing covers")
    v.legend("r", "$\\rho = \\sqrt{x^2 + y^2 + z^2}$ constant, in units of $c/H$")
    v.legend("t", "$t$ constant")
    v.legend("chartedge", "$t \\to -\\infty$, the past edge of the flat slicing")
    views.append(v)
    return views


# ---------------------------------------------------------------- anti-de Sitter and Bertotti-Robinson

def poincare_pq(t, z):
    """The Poincare patch in the strip of AdS2: from X_0 = Lt/z, X_-1 + X_1 = L^2/z and
    X_-1 - X_1 = z - t^2/z with L = 1, p = -pi/4 + arctan(t + z), q = pi/4 + arctan(t - z)."""
    t, z = np.asarray(t, dtype=float), np.asarray(z, dtype=float)
    return -Q4 + np.arctan(t + z), Q4 + np.arctan(t - z)


def strip(v, full, T0, T1, name="$\\mathscr{I}$"):
    X0 = -HALF if full else 0
    v.fill("region", [[X0, T0], [HALF, T0], [HALF, T1], [X0, T1]])
    v.line("boundary", [[[HALF, T0], [HALF, T1]]])
    if full:
        v.line("boundary", [[[-HALF, T0], [-HALF, T1]]])
        v.label_xt([-HALF, (T0 + T1) / 2 + 0.6], name, "r", dx=-6)
    else:
        v.line("centre", [[[0, T0], [0, T1]]])
    v.label_xt([HALF, (T0 + T1) / 2 + 0.6], name, "l", dx=6)


def anti_de_sitter(ck, src):
    """The strip. With sigma = arctan(r/L) the static global metric on the plane of t and r is
    (-c^2dt^2 + L^2 dsigma^2)/cos^2 sigma, conformal to the strip 0 <= sigma < pi/2 as it
    stands, p, q = (ct/L -+ sigma)/2. The Poincare plane x = y = 0 is a totally geodesic AdS2,
    drawn in its own strip by poincare_pq, checked against the embedding.
    """
    gl = Plane(src, "anti_de_sitter", "static_global", ("t", "r"), EQUATOR, {"L": 1})
    po = Plane(src, "anti_de_sitter", "poincare", ("t", "z"), {"x": "0", "y": "0"}, {"L": 1})

    def global_pq(t, r):
        t, r = np.asarray(t, dtype=float), np.asarray(r, dtype=float)
        return (t - np.arctan(r)) / 2, (t + np.arctan(r)) / 2
    ck.chart("anti-de Sitter global", gl, global_pq, ck.uniform(-10, 10), ck.uniform(0.001, 50), lambda t, r: (1, 0))
    ck.chart("anti-de Sitter Poincare", po, poincare_pq, ck.uniform(-10, 10), ck.uniform(0.001, 20), lambda t, z: (1, 0))
    t, z = ck.uniform(-6, 6, 2000), ck.uniform(1e-3, 8, 2000)
    p, q = poincare_pq(t, z)
    tau, s = p + q, q - p
    Xm, X0, X1 = np.cos(tau) / np.cos(s), np.sin(tau) / np.cos(s), np.tan(s)
    scale = 1 + np.abs(Xm) + np.abs(X0) + np.abs(X1)
    ck.limit("anti-de Sitter: the Poincare plane lands where the embedding puts it",
             np.concatenate([(X0 - t / z) / scale, (Xm + X1 - 1 / z) / scale, (Xm - X1 - z + t ** 2 / z) / scale]),
             0, 1e-10)
    p, q = poincare_pq(np.array([0.0, 3.0]), np.array([1e-12, 1e-12]))
    ck.limit("anti-de Sitter: z -> 0 lands on the boundary X = pi/2", q - p, [HALF, HALF], 1e-9)
    # sin(sigma) = k sin(ct/L) is a radial timelike geodesic: along it the energy
    # E = -g_tt dt/dtau, with dtau from the published metric, does not change.
    tt = np.linspace(0.05, PI - 0.05, 400)
    for k in (0.5, 0.9):
        sig = np.arcsin(k * np.sin(tt))
        r = np.tan(sig)
        drdt = k * np.cos(tt) / np.sqrt(1 - (k * np.sin(tt)) ** 2) / np.cos(sig) ** 2
        g00, _, g11, _, _, _ = gl.metric(tt, r)
        dtau = np.sqrt(-(g00 + g11 * drdt ** 2))
        E = -g00 / dtau
        ck.limit(f"anti-de Sitter: sin(sigma) = {k} sin(ct/L) keeps its energy, a timelike geodesic",
                 np.ptp(E) / np.mean(E), 0, 1e-10)
    ck.finite("anti-de Sitter: r = 0 is a regular centre", gl.kretschmann(ck.uniform(-5, 5, 50), np.full(50, 1e-6)))

    views = []
    T0, T1 = -0.3 * PI, 1.3 * PI
    v = View("global", "Static global", [-0.95, HALF + 0.55, T0, T1], "static_global")
    strip(v, False, T0, T1)
    v.fill("cover", [[0, T0], [HALF, T0], [HALF, T1], [0, T1]])
    for r in (0.25, 0.5, 1, 2, 4):
        X = float(np.arctan(r))
        v.line("r", [[[X, T0], [X, T1]]])
    for k in range(-1, 6):
        v.line("t", [[[0, k * Q4], [HALF, k * Q4]]])
    v.segment("null", (0, 0), (0, HALF))
    v.segment("null", (0, HALF), (HALF, HALF))
    tau = np.linspace(0, PI, 300)
    for k in (0.5, 0.9):
        sig = np.arcsin(k * np.sin(tau))
        v.curve("world", (tau - sig) / 2, (tau + sig) / 2)
    v.label_xt([0, 0], "$t = 0$", "r", "coord", dx=-6)
    v.label_xt([0, HALF], "$\\pi L/2c$", "r", "coord", dx=-6)
    v.label_xt([0, PI], "$\\pi L/c$", "r", "coord", dx=-6)
    v.label_xt([float(np.arctan(1)), -0.55], "$r = L$", "b", "coord", dy=-2)
    v.set(fade={"top": 0.7, "bottom": 0.7})
    v.legend("cover", "the whole universal cover, which $t$ and $r$ cover")
    v.legend("r", "$r$ constant, at $L/4$, $L/2$, $L$, $2L$ and $4L$")
    v.legend("t", "$ct$ constant, every $\\pi L/4$")
    v.legend("boundary", "the conformal boundary, timelike")
    v.legend("centre", "$r = 0$, a regular centre")
    v.legend("null", "a radial light ray from the centre")
    v.legend("world", "radial timelike geodesics from the centre, $\\sin\\sigma = k\\sin(ct/L)$ at $k = 0.5$ and $0.9$")
    views.append(v)

    T0, T1 = -1.1 * PI, 1.1 * PI
    v = View("poincare", "Poincaré patch", [-HALF - 0.55, HALF + 0.55, T0, T1], "poincare")
    strip(v, True, T0, T1)
    v.fill("cover", [[-HALF, 0], [HALF, -PI], [HALF, PI]])
    grid(v, "r", lambda z, t: poincare_pq(t, z), (0.25, 0.5, 1, 2, 4), S_ALL)
    grid(v, "t", poincare_pq, (-4, -2, -1, 0, 1, 2, 4), S_POS)
    v.line("chartedge", [[[-HALF, 0], [HALF, PI]], [[-HALF, 0], [HALF, -PI]]])
    label_on(v, poincare_pq(0, 1), "$z = L$")
    v.label_xt([0.05, 1.75], "Poincaré horizon", "br", "small", dx=-4, dy=-2)
    v.set(fade={"top": 0.7, "bottom": 0.7},
          restriction="The plane $x = y = 0$ of the Poincaré patch only, a totally geodesic "
                      "anti-de Sitter space of two dimensions; each point of the diagram is a "
                      "single event rather than a sphere of them.")
    v.legend("cover", "the wedge that $t$ and $z$ cover")
    v.legend("r", "$z$ constant")
    v.legend("t", "$ct$ constant")
    v.legend("boundary", "the conformal boundary, timelike")
    v.legend("chartedge", "$z \\to \\infty$, the Poincaré horizon, where $t$ and $z$ end")
    views.append(v)
    return views


def bertotti_robinson(ck, src):
    """AdS2 of radius b times a sphere of radius b: the strip of AdS2 is the whole diagram.
    Both coordinate systems are the Poincare patch of it, the throat's r becoming the
    Poincare x = b^2/r, so the static map is poincare_pq(t, 1/r) at b = 1."""
    st = Plane(src, "bertotti_robinson", "static", ("t", "r"), EQUATOR, {"b": 1})
    po = Plane(src, "bertotti_robinson", "poincare", ("t", "x"), EQUATOR, {"b": 1})
    ck.chart("Bertotti-Robinson static", st, lambda t, r: poincare_pq(t, 1 / np.asarray(r, dtype=float)),
             ck.uniform(-10, 10), ck.uniform(0.01, 20), lambda t, r: (1, 0))
    ck.chart("Bertotti-Robinson Poincare", po, poincare_pq, ck.uniform(-10, 10), ck.uniform(0.01, 20),
             lambda t, x: (1, 0))
    ck.finite("Bertotti-Robinson: the curvature is the same everywhere",
              st.kretschmann(ck.uniform(-5, 5, 50), ck.uniform(1e-3, 50, 50)))

    views = []
    T0, T1 = -1.1 * PI, 1.1 * PI
    box = [-HALF - 0.55, HALF + 0.55, T0, T1]
    for vid, label, system, z_of in (("static", "Static throat", "static", lambda r: 1.0 / r),
                                     ("poincare", "Poincaré", "poincare", lambda x: x)):
        v = View(vid, label, box, system)
        strip(v, True, T0, T1)
        v.fill("cover", [[-HALF, 0], [HALF, -PI], [HALF, PI]])
        for c in (0.25, 0.5, 1, 2, 4):
            v.curve("r", *poincare_pq(S_ALL, np.full_like(S_ALL, z_of(c))))
        grid(v, "t", poincare_pq, (-4, -2, -1, 0, 1, 2, 4), S_POS)
        v.line("chartedge", [[[-HALF, 0], [HALF, PI]], [[-HALF, 0], [HALF, -PI]]])
        if vid == "static":
            label_on(v, poincare_pq(0, 1.0), "$r = b$")
            v.legend("cover", "the wedge that $t$ and $r$ cover")
            v.legend("r", "$r$ constant, from $b/4$ to $4b$")
            v.legend("chartedge", "$r = 0$, the Poincaré horizon, where $t$ and $r$ end")
        else:
            label_on(v, poincare_pq(0, 1.0), "$x = b$")
            v.legend("cover", "the wedge that $t$ and $x$ cover, the same as the throat's")
            v.legend("r", "$x$ constant, from $b/4$ to $4b$")
            v.legend("chartedge", "$x \\to \\infty$, the Poincaré horizon, where $t$ and $x$ end")
        v.legend("t", "$ct$ constant")
        v.legend("boundary", "the conformal boundary, timelike")
        v.set(fade={"top": 0.7, "bottom": 0.7})
        views.append(v)
    return views


# ---------------------------------------------------------------- wormholes

def ellis_bronnikov(ck, src):
    """The metric on the plane of t and r is -c^2dt^2 + dr^2 with r over the whole line:
    the full diamond, p, q = arctan((ct -+ r)/l), drawn at l = 1."""
    pl = Plane(src, "ellis_bronnikov", "spherical", ("t", "r"), EQUATOR, {"ell": 1})
    ck.chart("Ellis-Bronnikov", pl, mink_pq, ck.uniform(-20, 20), ck.uniform(-20, 20), lambda t, r: (1, 0))
    ck.finite("Ellis-Bronnikov: the curvature is finite at the throat r = 0",
              pl.kretschmann(ck.uniform(-5, 5, 50), ck.uniform(-0.01, 0.01, 50)))
    v = View("spherical", "Spherical", [-PI - 0.35, PI + 0.35, -PI - 0.25, PI + 0.25], "spherical")
    v.fill("region", DIAMOND)
    v.fill("cover", DIAMOND)
    grid(v, "r", lambda r, t: mink_pq(t, r), (-4, -2, -1, -0.5, 0.5, 1, 2, 4), S_ALL)
    grid(v, "t", mink_pq, (-4, -2, -1, 0, 1, 2, 4), S_ALL)
    v.line("throat", [[[0, -PI], [0, PI]]])
    diamond_edges(v)
    label_on(v, mink_pq(0, 1), "$r = \\ell$")
    label_on(v, mink_pq(0, -1), "$-\\ell$")
    v.label_xt([0, 0.3], "throat", "l", "small", dx=6)
    v.legend("cover", "the whole spacetime, which $t$ and $r$ cover")
    v.legend("r", "$r$ constant, a sphere of area $4\\pi(r^2 + \\ell^2)$")
    v.legend("t", "$ct$ constant")
    v.legend("throat", "the throat $r = 0$, where the spheres are smallest")
    return [v]


def morris_thorne(ck, src):
    """Drawn with Phi = 0 and b = b_0^2/r, the member of the family that is Ellis-Bronnikov:
    the proper distance is l = +-sqrt(r^2 - b_0^2), and p, q = arctan((ct -+ l)/b_0) give the
    full diamond, drawn at b_0 = 1. The areal coordinates cover the side l > 0."""
    areal = Plane(src, "morris_thorne", "spherical", ("t", "r"), EQUATOR, {"b_0": 1},
                  functions={"Phi": "0", "b": "b_0**2/r"})
    proper = Plane(src, "morris_thorne", "proper_radial", ("t", "l"), EQUATOR, {"b_0": 1}, functions={"Phi": "0"})
    ck.chart("Morris-Thorne areal", areal, lambda t, r: mink_pq(t, np.sqrt(np.asarray(r) ** 2 - 1)),
             ck.uniform(-20, 20), ck.uniform(1.001, 20), lambda t, r: (1, 0))
    ck.chart("Morris-Thorne proper radial", proper, mink_pq, ck.uniform(-20, 20), ck.uniform(-20, 20),
             lambda t, l: (1, 0))
    ck.finite("Morris-Thorne: the curvature is finite at the throat r = b_0",
              areal.kretschmann(ck.uniform(-5, 5, 50), 1 + ck.uniform(1e-6, 1e-3, 50)))
    box = [-PI - 0.35, PI + 0.35, -PI - 0.25, PI + 0.25]
    views = []
    for vid, label, system in (("spherical", "Areal radius", "spherical"),
                               ("proper_radial", "Proper radial distance", "proper_radial")):
        v = View(vid, label, box, system)
        v.fill("region", DIAMOND)
        if vid == "spherical":
            v.fill("cover", TRIANGLE)
            for r in (1.1, 1.5, 2.5, 5):
                ell = np.sqrt(r * r - 1)
                v.curve("r", *mink_pq(S_ALL, ell))
                v.curve("r2", *mink_pq(S_ALL, -ell))
            label_on(v, mink_pq(0, np.sqrt(1.5 ** 2 - 1)), "$r = 1.5\\,b_0$")
            v.legend("cover", "the side $l > 0$, which $t$ and $r$ cover")
            v.legend("r", "$r$ constant, the areal radius")
            v.legend("r2", "the same radii on the other side, $l < 0$")
        else:
            v.fill("cover", DIAMOND)
            grid(v, "r", lambda ell, t: mink_pq(t, ell), (-4, -2, -1, -0.5, 0.5, 1, 2, 4), S_ALL)
            label_on(v, mink_pq(0, 1), "$l = b_0$")
            v.legend("cover", "the whole spacetime, which $t$ and $l$ cover")
            v.legend("r", "$l$ constant, in units of $b_0$")
        grid(v, "t", mink_pq, (-4, -2, -1, 0, 1, 2, 4), S_ALL)
        v.line("throat", [[[0, -PI], [0, PI]]])
        diamond_edges(v)
        v.label_xt([0, 0.3], "throat", "l", "small", dx=6)
        v.legend("t", "$ct$ constant")
        v.legend("throat", "the throat, $r = b_0$ and $l = 0$")
        v.set(input="$\\Phi = 0$ and $b = b_0^2/r$, the member of the family that is the "
                    "Ellis-Bronnikov wormhole.")
        views.append(v)
    return views


# ---------------------------------------------------------------- the cosmic string

def cosmic_string(ck, src):
    """The half plane of fixed phi and z, flat and totally geodesic: -c^2dt^2 + dr^2 outside,
    and -c^2dt^2 + l^2 dchi^2 in Gott's core. With rho the proper distance from the axis,
    rho = l chi inside and l chi_0 + r - l tan chi_0 outside, the core's edge r = l tan chi_0
    being where the circumferences agree, 2 pi l sin chi_0 = 2 pi (1 - 4G mu/c^2) r with
    cos chi_0 = 1 - 4G mu/c^2, the whole half plane is -c^2dt^2 + drho^2: Minkowski's half
    diamond, p, q = arctan((ct -+ rho)/l), drawn at l = 1 and 4G mu/c^2 = 0.1."""
    chi0 = float(np.arccos(0.9))
    shift = chi0 - np.tan(chi0)
    con = Plane(src, "cosmic_string", "conical", ("t", "r"), {"phi": "0", "z": "0"})
    cap = Plane(src, "cosmic_string", "interior_cap", ("t", "\\chi"), {"phi": "0", "z": "0"}, {"ell": 1})
    ck.chart("cosmic string, conical exterior", con, lambda t, r: mink_pq(t, np.asarray(r) + shift),
             ck.uniform(-20, 20), ck.uniform(np.tan(chi0), 20), lambda t, r: (1, 0))
    ck.chart("cosmic string, Gott's core", cap, mink_pq, ck.uniform(-20, 20), ck.uniform(0.001, chi0),
             lambda t, c: (1, 0))
    ck.finite("cosmic string: the curvature of Gott's core is finite on its axis",
              cap.kretschmann(ck.uniform(-5, 5, 50), np.full(50, 1e-6)))
    box = [-0.35, PI + 0.35, -PI - 0.25, PI + 0.25]
    TS = (-4, -2, -1, 0, 1, 2, 4)
    views = []

    v = View("conical", "Conical exterior", box, "conical")
    v.fill("region", TRIANGLE)
    v.fill("cover", TRIANGLE)
    grid(v, "r", lambda r, t: mink_pq(t, r), (0.5, 1, 2, 4), S_ALL)
    grid(v, "t", mink_pq, TS, S_POS)
    triangle_edges(v, centre_class="surface")
    v.label_xt([0, -0.4], "the string", "r", "small", dx=-6)
    v.legend("cover", "the region that $t$ and $r$ cover")
    v.legend("r", "$r$ constant, the proper distance from the string")
    v.legend("t", "$ct$ constant")
    v.legend("surface", "the string, a conical singularity at $r = 0$")
    v.set(restriction="The half plane of fixed $\\phi$ and $z$ only, which is totally geodesic. Each "
                      "point of the diagram stands for a circle around the string times a line "
                      "along it, and the string's deficit angle $\\delta = 8\\pi G\\mu/c^2$ shows in "
                      "those circles.")
    views.append(v)

    v = View("gott", "Gott's core", box, "interior_cap")
    v.fill("region", TRIANGLE)
    pc, qc = mink_pq(S_ALL, chi0)
    edge = [point(p, q) for p, q in zip(pc, qc)]
    v.fill("cover2", TRIANGLE)
    v.fill("star", [[0, -PI]] + edge + [[0, PI]])
    for c in (0.15, 0.3):
        v.curve("r2", *mink_pq(S_ALL, c))
    for r in (1, 2, 4):
        v.curve("r", *mink_pq(S_ALL, r + shift))
    grid(v, "t", mink_pq, TS, S_POS)
    v.curve("surface", pc, qc)
    triangle_edges(v, centre="$\\chi = 0$")
    v.legend("star", "the core, $\\chi \\le \\chi_0$, which $t$ and $\\chi$ cover")
    v.legend("cover2", "the conical exterior")
    v.legend("surface", "the edge of the core, $\\chi = \\chi_0$")
    v.legend("r2", "$\\chi$ constant")
    v.legend("r", "$r$ constant, at $\\ell$, $2\\ell$ and $4\\ell$")
    v.legend("t", "$ct$ constant")
    v.legend("centre", "$\\chi = 0$, a regular axis")
    v.set(settings="$4G\\mu/c^2 = 0.1$, so that $\\cos\\chi_0 = 0.9$.",
          restriction="The half plane of fixed $\\phi$ and $z$ only, in the proper distance $\\rho$ "
                      "from the axis; each point of the diagram stands for a circle around the "
                      "axis times a line along it.")
    views.append(v)
    return views


# ---------------------------------------------------------------- the interior Schwarzschild star

def interior_schwarzschild(ck, src):
    """The published interior, r <= R, joined at R = 1.5 r_s to the Schwarzschild exterior of
    the schwarzschild entry. g_tt is continuous at R, checked, so t is one coordinate, and
    r* = int sqrt(g_rr/(-g_tt)) dr runs from the centre through the surface; p, q =
    arctan((t -+ r*)/R) give Minkowski's triangle with the star a timelike tube."""
    inner = Plane(src, "interior_schwarzschild", "spherical", ("t", "r"), EQUATOR, {"r_s": 1, "R": "3/2"})
    outer = Plane(src, "schwarzschild", "spherical", ("t", "r"), EQUATOR, {"r_s": 1})
    R = 1.5
    ck.limit("interior Schwarzschild: g_tt is continuous at r = R",
             [float(inner.g[0, 0].subs(inner.x1, R))], [float(outer.g[0, 0].subs(outer.x1, R))], 1e-12)
    speed = sp.lambdify(inner.x1, sp.sqrt(inner.g[1, 1] / (-inner.g[0, 0])), "numpy")
    rq = np.linspace(0, R, 300001)
    rsq = cumulative_trapezoid(speed(rq), rq, initial=0)

    def rstar(r):
        r = np.asarray(r, dtype=float)
        with np.errstate(invalid="ignore", divide="ignore"):
            outside = rsq[-1] + (r + np.log(np.abs(r - 1))) - (R + np.log(R - 1))
        return np.where(r <= R, np.interp(r, rq, rsq), outside)

    def star(t, r):
        return mink_pq(t, rstar(r), R)
    ck.chart("interior Schwarzschild, the star", inner, star, ck.uniform(-20, 20), ck.uniform(0.01, 1.49),
             lambda t, r: (1, 0))
    ck.chart("interior Schwarzschild, the exterior", outer, star, ck.uniform(-20, 20), ck.uniform(1.51, 30),
             lambda t, r: (1, 0))
    ck.finite("interior Schwarzschild: r = 0 is a regular centre",
              inner.kretschmann(ck.uniform(-5, 5, 50), np.full(50, 1e-6)))

    v = View("spherical", "Spherical", [-0.35, PI + 0.35, -PI - 0.25, PI + 0.25], "spherical")
    v.fill("region", TRIANGLE)
    ps, qs = mink_pq(S_ALL, rstar(R), R)
    v.fill("cover", [[0, -PI]] + [point(p, q) for p, q in zip(ps, qs)] + [[0, PI]])
    for r in (0.5, 1.0):
        v.curve("r", *star(S_ALL, np.full_like(S_ALL, r)))
    for r in (2.0, 3.0, 6.0):
        v.curve("r2", *star(S_ALL, np.full_like(S_ALL, r)))
    rr = np.concatenate([np.linspace(0, R, 60)[:-1], R + np.exp(np.linspace(-6, 8, 200)) - np.exp(-6)])
    for t in (-6, -3, -1.5, 0, 1.5, 3, 6):
        v.curve("t", *star(np.full_like(rr, t), rr))
    v.curve("surface", ps, qs)
    triangle_edges(v)
    label_on(v, mink_pq(0, rstar(R), R), "$r = R$")
    v.label_xt([0.35, 0.0], "star", cls="region")
    v.legend("cover", "the star, $r \\le R$, which the interior solution covers")
    v.legend("r", "$r$ constant inside, at $0.5$ and $1\\,r_s$")
    v.legend("r2", "$r$ constant outside, at $2$, $3$ and $6\\,r_s$")
    v.legend("t", "$t$ constant, one $t$ on both sides")
    v.legend("surface", "the surface $r = R$")
    v.legend("centre", "$r = 0$, a regular centre")
    v.set(settings="$R = 1.5\\,r_s$.")
    return [v]


# ---------------------------------------------------------------- FRW

def frw(ck, src):
    """Dust universes, the scale factor solved from the published G^r_r = 0 of the conformal
    coordinates, lengths in units of 1/sqrt|k| where k is not zero:

        k = 0    a = eta^2,          conformal to eta > 0 of Minkowski space: the triangle
        k = +1   a = 1 - cos eta,    r = sin chi, conformal to the Einstein static universe
                                     as it stands: the rectangle 0 < eta < 2 pi, 0 <= chi <= pi
        k = -1   a = cosh eta - 1,   r = sinh chi, and tan((T +- X)/2) = tanh((eta +- chi)/2)

    each scale factor checked in sympy to make the published G^r_r vanish. The comoving
    coordinates are the same map through t = int a d eta.
    """
    k0 = Plane(src, "frw", "conformal_spherical", ("\\eta", "r"), EQUATOR, {"k": 0}, functions={"a": "eta**2"})
    src.note("frw", "conformal_spherical", ["einstein_tensor"])
    for k, a in ((0, "eta**2"), (1, "1 - cos(eta)"), (-1, "cosh(eta) - 1")):
        pl = Plane(src, "frw", "conformal_spherical", ("\\eta", "r"), EQUATOR, {"k": k}, functions={"a": a})
        ul = pl.entry["einstein_tensor"]["variants"]["ul"]["nonzero"]
        G = pl.prep(pl.reader(next(c["value"] for c in ul if c["indices"] == ["r", "r"])))
        ck.limit(f"FRW: a = {a} makes the published G^r_r vanish, k = {k}", [float(sp.simplify(G) != 0)], [0], 0.5)

    def closed(e, r):
        return (np.asarray(e) - np.arcsin(r)) / 2, (np.asarray(e) + np.arcsin(r)) / 2

    def open_(e, r):
        chi = np.arcsinh(r)
        return np.arctan(np.tanh((np.asarray(e) - chi) / 2)), np.arctan(np.tanh((np.asarray(e) + chi) / 2))
    k1 = Plane(src, "frw", "conformal_spherical", ("\\eta", "r"), EQUATOR, {"k": 1}, functions={"a": "1 - cos(eta)"})
    km = Plane(src, "frw", "conformal_spherical", ("\\eta", "r"), EQUATOR, {"k": -1}, functions={"a": "cosh(eta) - 1"})
    ck.chart("FRW conformal, dust, k = 0", k0, mink_pq, ck.uniform(0.01, 20), ck.uniform(0.001, 20), lambda e, r: (1, 0))
    ck.chart("FRW conformal, dust, k = 1", k1, closed, ck.uniform(0.01, 2 * PI - 0.01), ck.uniform(0.001, 0.999),
             lambda e, r: (1, 0))
    ck.chart("FRW conformal, dust, k = -1", km, open_, ck.uniform(0.01, 20), ck.uniform(0.001, 20), lambda e, r: (1, 0))
    # The comoving coordinates: t = int a d eta, eta^3/3, eta - sin eta and sinh eta - eta.
    etas = np.linspace(1e-6, 2 * PI - 1e-6, 200001)
    for k, tau, a_of, pq, top in ((0, etas ** 3 / 3, lambda e: e ** 2, mink_pq, 20.0),
                                  (1, etas - np.sin(etas), lambda e: 1 - np.cos(e), closed, 2 * PI - 0.05),
                                  (-1, np.sinh(etas) - etas, lambda e: np.cosh(e) - 1, open_, 5.5)):
        def eta_of(t, tau=tau):
            return np.interp(t, tau, etas)
        pl = Plane(src, "frw", "comoving_spherical", ("t", "r"), EQUATOR, {"k": k}, numeric=["a"])
        fvals = (lambda a_of, eta_of: lambda t, r: {"a": (a_of(eta_of(t)), 0 * t, 0 * t)})(a_of, eta_of)
        hi = float(np.interp(top, etas, tau))
        ck.chart(f"FRW comoving, dust, k = {k}", pl, lambda t, r, pq=pq, eta_of=eta_of: pq(eta_of(t), r),
                 ck.uniform(0.01, hi), ck.uniform(0.001, 0.999 if k == 1 else 20), lambda t, r: (1, 0), fvals)
    for k, pl, a in ((0, k0, lambda e: e ** 2), (1, k1, lambda e: 1 - np.cos(e)), (-1, km, lambda e: np.cosh(e) - 1)):
        ck.diverges(f"FRW: the Kretschmann scalar diverges at the bang, k = {k}",
                    pl.kretschmann(1e-2, 0.5), pl.kretschmann(1e-3, 0.5))
    ck.diverges("FRW: the Kretschmann scalar diverges at the crunch, k = 1",
                k1.kretschmann(2 * PI - 1e-2, 0.5), k1.kretschmann(2 * PI - 1e-3, 0.5))

    views = []
    tri = [[0, 0], [PI, 0], [0, PI]]
    v = View("flat", "Flat dust, $k = 0$", [-0.35, PI + 0.35, -0.25, PI + 0.25])
    v.fill("region", tri)
    v.fill("cover", tri)
    rr = spread(0, np.inf, 500, 12)
    for eta in (0.25, 0.5, 1, 2, 4):
        v.curve("t", *mink_pq(np.full_like(rr, eta), rr))
    for r in (0.25, 0.5, 1, 2, 4):
        v.curve("r", *mink_pq(rr, np.full_like(rr, r)))
    v.line("singular", [[[0, 0], [PI, 0]]], zig=True)
    v.line("scri", [[[0, PI], [PI, 0]]])
    v.line("centre", [[[0, 0], [0, PI]]])
    now = mink_pq(1, 0)
    v.segment("null", now, mink_pq(0, 1))
    v.point("mark", now)
    v.layers.append({"kind": "point", "class": "infinity", "at": [0, round(PI, 4)]})
    v.layers.append({"kind": "point", "class": "infinity", "at": [round(PI, 4), 0]})
    v.label_xt([0, PI], "$i^+$", "b", dy=-6)
    v.label_xt([PI, 0], "$i^0$", "l", dx=6)
    v.label_xt([HALF, HALF], "$\\mathscr{I}^+$", "bl", dx=5, dy=-3)
    v.label_xt([HALF, 0], "the big bang, $t = 0$", "t", dy=8)
    v.label(now, "here, now", "r", "small", dx=-6)
    label_on(v, mink_pq(1, 1.6), "$t_0$")
    label_on(v, mink_pq(2, 2.4), "$8t_0$")
    v.legend("t", "cosmic time $t$ constant: $t_0/64$, $t_0/8$, $t_0$, $8t_0$ and $64t_0$, with $t \\propto \\eta^3$")
    v.legend("r", "comoving $r$ constant, in units of $\\eta_0$")
    v.legend("null", "our past light cone, which meets the bang at the particle horizon $r = \\eta_0$")
    v.legend("singular", "the big bang, where the Kretschmann scalar diverges")
    v.legend("centre", "$r = 0$, our world line")
    v.legend("cover", "the whole spacetime, which $r$ covers with $t$ or with $\\eta$")
    views.append(v)

    v = View("closed", "Closed dust, $k = +1$", [-0.35, PI + 0.35, -0.25, 2 * PI + 0.25])
    v.fill("region", [[0, 0], [PI, 0], [PI, 2 * PI], [0, 2 * PI]])
    v.fill("cover", [[0, 0], [HALF, 0], [HALF, 2 * PI], [0, 2 * PI]])
    for eta in np.arange(1, 7) * PI / 3.5:
        v.line("t", [[[0, eta], [PI, eta]]])
    for chi in (PI / 6, PI / 3):
        v.line("r", [[[chi, 0], [chi, 2 * PI]]])
    for chi in (2 * PI / 3, 5 * PI / 6):
        v.line("r2", [[[chi, 0], [chi, 2 * PI]]])
    v.line("chartedge", [[[HALF, 0], [HALF, 2 * PI]]])
    v.line("singular", [[[0, 0], [PI, 0]], [[0, 2 * PI], [PI, 2 * PI]]], zig=True)
    v.line("centre", [[[0, 0], [0, 2 * PI]], [[PI, 0], [PI, 2 * PI]]])
    v.line("null", [[[0, 0], [PI, PI]]])
    v.label_xt([HALF, 0], "the big bang", "t", dy=8)
    v.label_xt([HALF, 2 * PI], "the big crunch", "b", dy=-8)
    v.label_xt([0, PI], "$\\chi = 0$", "r", dx=-6)
    v.label_xt([PI, PI], "$\\chi = \\pi$", "l", dx=6)
    v.label_xt([HALF, 2 * PI - 0.35], "$r = 1$", "l", "small", dx=4)
    v.legend("t", "conformal time $\\eta$ constant, every $\\pi/3.5$")
    v.legend("r", "$\\chi$ constant, with $r = \\sin\\chi$")
    v.legend("r2", "$\\chi$ constant on the far hemisphere, $\\chi > \\pi/2$")
    v.legend("chartedge", "$\\chi = \\pi/2$, the equator, where $1 - kr^2$ vanishes at $r = 1$")
    v.legend("null", "the first light to cross the whole universe, which takes until maximum expansion")
    v.legend("singular", "the bang and the crunch, where the Kretschmann scalar diverges")
    v.legend("centre", "$\\chi = 0$ and its antipode $\\chi = \\pi$")
    v.legend("cover", "the hemisphere $\\chi < \\pi/2$, which $r$ covers, with $t$ or with $\\eta$")
    views.append(v)

    v = View("open", "Open dust, $k = -1$", [-0.35, HALF + 0.35, -0.2, HALF + 0.2])
    tri = [[0, 0], [HALF, 0], [0, HALF]]
    v.fill("region", tri)
    v.fill("cover", tri)
    cc = spread(0, np.inf, 500, 5)
    for eta in (0.5, 1, 2, 3, 4):
        v.curve("t", *open_(np.full_like(cc, eta), np.sinh(cc)))
    for chi in (0.5, 1, 2, 3):
        v.curve("r", *open_(cc, np.full_like(cc, np.sinh(chi))))
    v.line("singular", [[[0, 0], [HALF, 0]]], zig=True)
    v.line("scri", [[[0, HALF], [HALF, 0]]])
    v.line("centre", [[[0, 0], [0, HALF]]])
    v.layers.append({"kind": "point", "class": "infinity", "at": [0, round(HALF, 4)]})
    v.layers.append({"kind": "point", "class": "infinity", "at": [round(HALF, 4), 0]})
    v.label_xt([0, HALF], "$i^+$", "b", dy=-6)
    v.label_xt([HALF, 0], "$i^0$", "l", dx=6)
    v.label_xt([Q4, Q4], "$\\mathscr{I}^+$", "bl", dx=5, dy=-3)
    v.label_xt([Q4, 0], "the big bang", "t", dy=8)
    v.legend("t", "conformal time $\\eta$ constant, at $0.5$, $1$, $2$, $3$ and $4$")
    v.legend("r", "$\\chi$ constant, with $r = \\sinh\\chi$")
    v.legend("singular", "the big bang, where the Kretschmann scalar diverges")
    v.legend("scri", "null infinity $\\mathscr{I}^+$")
    v.legend("centre", "$\\chi = 0$")
    v.legend("cover", "the whole spacetime, which $r$ covers with $t$ or with $\\eta$")
    views.append(v)
    for view in views:
        view.set(input="Dust: $a \\propto \\eta^2$, $1 - \\cos\\eta$ and $\\cosh\\eta - 1$ for $k = 0$, $+1$ "
                       "and $-1$, each solved from this spacetime's own $G^r{}_r = 0$, with lengths in "
                       "units of $1/\\sqrt{|k|}$ where $k$ is not zero.")
    return views


# ---------------------------------------------------------------- Oppenheimer-Snyder

class Collapse:
    """Dust released from rest at R0, r_s = 1, drawn with the dust in its own conformal time
    and the exterior fitted to it.

    Inside, the dust is closed FRW, a^2(-d eta^2 + d chi^2 + ...) with
    a = (a_m/2)(1 + cos eta) and c tau = (a_m/2)(eta + sin eta), so its plane of eta and chi
    is drawn as it stands: p = (eta - chi)/2, q = (eta + chi)/2, a rectangle. Outside, the
    drawing is P(U) and Q(V) of the Kruskal coordinates, fixed by three conditions:

      1. on the surface chi = chi0 the two sides agree: P(U_s(eta)) = (eta - chi0)/2 and
         Q(V_s(eta)) = (eta + chi0)/2, for the U and V the surface reaches;
      2. the moment of rest is the line T = 0: P(U) = -Q(-U), the exterior's t -> -t
         symmetry about it, U <-> -V;
      3. r = 0 is the line T = pi: Q(V) = pi - P(1/V), since UV = 1 there.

    These leave nothing free: 1 fixes P and Q where the surface runs, 3 extends Q to the
    ingoing rays that arrive after the star has gone, and 2 extends P to the outgoing rays
    that start outside the star. The surface is the exterior's radial geodesic of energy
    cos chi0 carried in v, dv/d eta = a/(E + sqrt(E^2 - f)), which is regular at r_s; the
    origin of Schwarzschild time drops out, since P and Q are defined through the surface.
    Near the crunch U and V stall as (pi - eta)^4, below double precision, so the stalled
    points are dropped and the end point set on UV = 1.
    """

    def __init__(self, R0=2.0):
        self.R0 = R0
        self.chi0 = float(np.arcsin(np.sqrt(1 / R0)))
        self.am = R0 / np.sin(self.chi0)
        E, am = np.cos(self.chi0), self.am

        def rhs(eta, y):
            R = (R0 / 2) * (1 + np.cos(eta))
            a = (am / 2) * (1 + np.cos(eta))
            f = 1 - 1 / R if R > 0 else -np.inf
            return [a / (E + np.sqrt(max(E * E - f, 0.0)))]
        eta = np.concatenate([np.linspace(0, np.pi - 0.2, 20000), np.pi - 0.2 * np.logspace(0, -6, 20000)[1:], [np.pi]])
        sol = solve_ivp(rhs, (0, np.pi), [R0 + np.log(R0 - 1)], t_eval=eta, rtol=1e-13, atol=1e-13, method="DOP853")
        self.eta = sol.t
        R = (R0 / 2) * (1 + np.cos(self.eta))
        self.V = np.exp(sol.y[0] / 2)
        self.U = (1 - R) * np.exp(R) / self.V
        self.R = R
        keep = np.flatnonzero(np.concatenate([[True], (np.diff(self.U) > 0) & (np.diff(self.V) > 0)]))
        keep = keep[np.concatenate([[True], np.diff(self.U[keep]) > 0])]
        self.eta, self.U, self.V, self.R = (x[keep] for x in (self.eta, self.U, self.V, self.R))
        if self.R[-1] == 0:
            self.V[-1] = 1 / self.U[-1]
        assert np.all(np.diff(self.U) > 0) and np.all(np.diff(self.V) > 0), "the surface is not monotone"

    def P(self, U):
        U = np.asarray(U, dtype=float)
        inside = (np.interp(U, self.U, self.eta) - self.chi0) / 2
        outside = -self.Q(np.maximum(-U, self.V[0]))
        return np.where(U >= self.U[0], inside, outside)

    def Q(self, V):
        V = np.asarray(V, dtype=float)
        on_surface = (np.interp(V, self.V, self.eta) + self.chi0) / 2
        with np.errstate(divide="ignore"):
            late = np.pi - (np.interp(1 / np.maximum(V, 1e-300), self.U, self.eta) - self.chi0) / 2
        return np.where(V <= self.V[-1], on_surface, late)


def oppenheimer_snyder(ck, src):
    o = Collapse(2.0)
    c = o.chi0
    ext = Plane(src, "oppenheimer_snyder", "exterior_schwarzschild", ("t", "r"), EQUATOR, {"r_s": 1})
    T = Tower(-ext.g[0, 0], ext.x1, [1])
    inner = Plane(src, "oppenheimer_snyder", "interior_comoving", ("\\tau", "\\chi"), EQUATOR,
                  {"chi_0": sp.nsimplify(c), "a_m": sp.nsimplify(o.am)}, numeric=["a"])
    etas = np.linspace(0, np.pi, 200001)
    taus = (o.am / 2) * (etas + np.sin(etas))

    def eta_of(tau):
        return np.interp(tau, taus, etas)

    def scale(tau, chi):
        e = eta_of(tau)
        return {"a": ((o.am / 2) * (1 + np.cos(e)), -np.sin(e) / (1 + np.cos(e)),
                      -1 / ((o.am / 2) * (1 + np.cos(e)) ** 2))}

    def inside(tau, chi):
        e = eta_of(tau)
        return (e - chi) / 2, (e + chi) / 2

    def outside(t, r, cell="I"):
        p, q = T.pq(cell, t, r)
        return o.P(np.tan(p)), o.Q(np.tan(q))
    ck.chart("Oppenheimer-Snyder, the dust", inner, inside, ck.uniform(0.01, taus[-1] * 0.98),
             ck.uniform(0.001, c), lambda tau, chi: (1, 0), scale)
    ck.chart("Oppenheimer-Snyder, the exterior", ext, outside, ck.uniform(0.2, 12), ck.uniform(2.05, 30),
             lambda t, r: (1, 0))
    ck.limit("Oppenheimer-Snyder: the surface reaches r = 0 on UV = 1", [o.U[-1] * o.V[-1]], [1.0], 1e-9)
    ck.limit("Oppenheimer-Snyder: the moment of rest is one line, P(-V) + Q(V) = 0",
             [o.P(-V) + o.Q(V) for V in (3.0, 10.0, 1e4)], [0, 0, 0], 1e-9)
    ck.limit("Oppenheimer-Snyder: r = 0 is one line, P(1/V) + Q(V) = pi",
             [o.P(1 / V) + o.Q(V) for V in (30.0, 100.0, 1e5)], [PI] * 3, 1e-9)
    ck.limit("Oppenheimer-Snyder: the event horizon U = 0 is P = (pi - 3 chi0)/2", [o.P(0.0)], [(PI - 3 * c) / 2], 1e-7)
    ck.limit("Oppenheimer-Snyder: the surface crosses r = r_s at eta = pi - 2 chi0",
             [o.eta[np.argmin(np.abs(o.R - 1))]], [PI - 2 * c], 1e-3)
    # Marginally trapped spheres: the areal radius a sin(chi) has a null gradient on
    # eta = pi - 2 chi, by the published interior metric.
    chi = np.linspace(0.02, c, 50)
    e = PI - 2 * chi
    tau = (o.am / 2) * (e + np.sin(e))
    a, adot, _ = scale(tau, chi)["a"]
    _, _, _, h00, h01, h11 = inner.metric(tau, chi, scale)
    dR = [adot * np.sin(chi), a * np.cos(chi)]
    ck.limit("Oppenheimer-Snyder: the areal radius has a null gradient on eta = pi - 2 chi",
             (h00 * dR[0] ** 2 + 2 * h01 * dR[0] * dR[1] + h11 * dR[1] ** 2) / (dR[1] ** 2 * np.abs(h11)), 0, 1e-6)
    near, nearer = taus[-1] * (1 - 1e-2), taus[-1] * (1 - 1e-3)
    ck.diverges("Oppenheimer-Snyder: the Kretschmann scalar diverges at the crunch",
                inner.kretschmann(near, 0.3, scale), inner.kretschmann(nearer, 0.3, scale))
    ck.diverges("Oppenheimer-Snyder: the Kretschmann scalar diverges at r = 0 outside",
                ext.kretschmann(0, 1e-2), ext.kretschmann(0, 1e-3))
    ck.finite("Oppenheimer-Snyder: the centre chi = 0 is regular before the crunch",
              inner.kretschmann(taus[-1] * np.linspace(0, 0.9, 10), np.full(10, 1e-6), scale))

    Xmax = PI + 3 * c
    iplus = [3 * c, PI]
    v = View("collapse", "Interior and exterior", [-0.35, Xmax + 0.35, -0.25, PI + 0.3])
    outer_region = [[c, 0], [c, PI], iplus, [Xmax, 0]]
    dust = [[0, 0], [c, 0], [c, PI], [0, PI]]
    v.fill("region", outer_region)
    v.fill("region", dust)
    v.fill("star", dust)
    v.fill("cover", outer_region)
    for eta in (0.5, 1.0, 1.5, 2.0, 2.5, 2.9):
        v.line("t2", [[[0, eta], [c, eta]]])
    for x in (c / 3, 2 * c / 3):
        v.line("r2", [[[x, 0], [x, PI]]])

    def ext_curve(cls, t, r, cell):
        P, Q = outside(t, r, cell)
        keep = (Q - P > c + 1e-9) & (P + Q > -1e-9)
        v.curve(cls, np.where(keep, P, np.nan), np.where(keep, Q, np.nan))
    t = spread(-np.inf, np.inf, 2000, 11)
    for r in (0.3, 0.6, 1.25, 1.6, 2.5, 5, 10):
        ext_curve("r", t, np.full_like(t, r), "I" if r > 1 else "II")
    rr = spread(1, np.inf, 2000, 14)
    for tt in (1, 2.5, 5, 10):
        ext_curve("t", np.full_like(rr, tt), rr, "I")
    v.line("event", [[[0, PI - 3 * c], iplus]])
    v.line("apparent", [[[c, PI - 2 * c], [0, PI]]])
    v.line("surface", [[[c, 0], [c, PI]]])
    v.line("centre", [[[0, 0], [0, PI]]])
    v.line("singular", [[[0, PI], iplus]], zig=True)
    v.line("scri", [[iplus, [Xmax, 0]]])
    v.line("chartedge", [[[0, 0], [Xmax, 0]]])
    v.layers.append({"kind": "point", "class": "infinity", "at": rounded(iplus)})
    v.layers.append({"kind": "point", "class": "infinity", "at": rounded([Xmax, 0])})
    v.label_xt(iplus, "$i^+$", "bl", dx=4, dy=-3)
    v.label_xt([Xmax, 0], "$i^0$", "l", dx=6)
    v.label_xt([(3 * c + Xmax) / 2, HALF], "$\\mathscr{I}^+$", "bl", dx=5, dy=-3)
    v.label_xt([1.5 * c, PI], "$r = 0$", "b", dy=-8)
    v.label_xt([Xmax / 2, 0], "released from rest: $\\tau = 0$ inside, $t = 0$ outside", "t", "coord", dy=7)
    v.label_xt([c / 2, 0.25], "dust", cls="region")
    v.label_xt([c, 0.9], "$\\chi = \\chi_0$", "l", "small", dx=5)
    v.legend("star", "the dust, in its conformal time $\\eta$ and $\\chi$, which $\\tau$ and $\\chi$ cover")
    v.legend("cover", "the Schwarzschild exterior, $r \\ge R(t)$, which $t$ and $r$ cover")
    v.legend("surface", "the surface $\\chi = \\chi_0$, a radial geodesic of the exterior")
    v.legend("t2", "$\\eta$ constant inside")
    v.legend("r2", "$\\chi$ constant inside, the world lines of the dust")
    v.legend("r", "$r$ constant outside: $0.3$, $0.6$, $1.25$, $1.6$, $2.5$, $5$ and $10\\,r_s$")
    v.legend("t", "$ct$ constant outside")
    v.legend("event", "the event horizon, from the centre at $\\eta = \\pi - 3\\chi_0$ to $i^+$")
    v.legend("apparent", "marginally trapped spheres inside the dust, $\\eta = \\pi - 2\\chi$")
    v.legend("singular", "$r = 0$: the crunch inside and the singularity outside, where the "
                         "Kretschmann scalar diverges")
    v.legend("centre", "$\\chi = 0$, the centre")
    v.legend("chartedge", "the moment of rest")
    v.set(settings="$R_0 = 2\\,r_s$, so that $\\chi_0 = \\pi/4$ and $a_m = 2\\sqrt{2}\\,r_s$; the collapse "
                   "starts at the moment of rest, $\\tau = 0$.")
    return [v]


# ---------------------------------------------------------------- Vaidya

def vaidya(ck, src):
    """An imploding null shell: the ingoing coordinates with m = 0 for v < 0 and M for v > 0,
    r_s = 2GM/c^2 = 1. Outside the shell, Kruskal's p = arctan U, q = arctan V. Inside, flat
    space with retarded time u = v - 2r, and each outgoing ray drawn where it crosses the
    shell: there U = (1 - r) e^r and r = -u/2, so p = F(u) = arctan((1 + u/2) e^(-u/2)), and
    q = F(v), which puts the centre u = v on the straight line X = 0. The event horizon
    U = 0 is u = -2 inside, reaching the centre at v = -2."""
    flat = Plane(src, "vaidya", "eddington_finkelstein_ingoing", ("v", "r"), EQUATOR, {"G": 1}, functions={"m": "0"})
    hole = Plane(src, "vaidya", "eddington_finkelstein_ingoing", ("v", "r"), EQUATOR, {"G": 1}, functions={"m": "1/2"})

    def F(w):
        w = np.asarray(w, dtype=float)
        return np.arctan((1 + w / 2) * np.exp(-w / 2))

    def inside(w, r):
        return F(np.asarray(w) - 2 * np.asarray(r)), F(w)

    def outside(w, r):
        w, r = np.asarray(w, dtype=float), np.asarray(r, dtype=float)
        return np.arctan((1 - r) * np.exp(r - w / 2)), atan_exp(w / 2)
    ck.chart("Vaidya, flat inside the shell", flat, inside, ck.uniform(-15, -0.01), ck.uniform(0.01, 20),
             lambda w, r: (1, -60))
    ck.chart("Vaidya, Schwarzschild outside the shell", hole, outside, ck.uniform(0.01, 15), ck.uniform(0.01, 20),
             lambda w, r: (1, -60))
    r = ck.uniform(0.01, 10, 200)
    ck.limit("Vaidya: the two sides put the shell v = 0 at one place", inside(np.zeros_like(r), r),
             outside(np.zeros_like(r), r), 1e-12)
    ck.limit("Vaidya: the event horizon reaches the centre at v = -2r_s", [F(-2.0)], [0], 1e-12)
    ck.diverges("Vaidya: the Kretschmann scalar diverges at r = 0 after the shell",
                hole.kretschmann(1, 1e-2), hole.kretschmann(1, 1e-3))
    ck.finite("Vaidya: r = 0 before the shell is a regular centre", flat.kretschmann(ck.uniform(-9, -1, 20), np.full(20, 1e-6)))

    v = View("shell", "Imploding null shell", [-0.35, PI + 0.35, -PI - 0.25, HALF + 0.35], "eddington_finkelstein_ingoing")
    out_region = [point(-HALF, Q4), point(Q4, Q4), point(0, HALF), point(-HALF, HALF)]
    in_region = [point(-HALF, -HALF), point(Q4, Q4), point(-HALF, Q4)]
    v.fill("region", out_region)
    v.fill("region", in_region)
    v.fill("cover", out_region)
    v.fill("cover2", in_region)
    w = -spread(0, np.inf, 600, 10)[::-1]
    for r in (0.5, 1, 2, 4):
        v.curve("r2", *inside(w, np.full_like(w, r)))
    for tt in (-6, -4, -2, -1):
        rr = np.linspace(0, -tt, 400)
        v.curve("t2", *inside(tt + rr, rr))
    w = spread(0, np.inf, 600, 10)
    for r in (0.5, 0.8, 1.25, 1.6, 2.5, 5):
        v.curve("r", *outside(w, np.full_like(w, r)))
    v.segment("surface", (-HALF, Q4), (Q4, Q4))
    v.segment("event", (0, 0), (0, HALF))
    v.segment("centre", (-HALF, -HALF), (Q4, Q4))
    v.segment("singular", (Q4, Q4), (0, HALF), zig=True)
    v.segment("scri", (0, HALF), (-HALF, HALF))
    v.segment("scri", (-HALF, HALF), (-HALF, -HALF))
    for pq, text, anchor, dx, dy in (((-HALF, -HALF), "$i^-$", "t", 0, 6), ((-HALF, HALF), "$i^0$", "l", 6, 0),
                                     ((0, HALF), "$i^+$", "bl", 4, -3)):
        v.point("infinity", pq)
        v.label(pq, text, anchor, dx=dx, dy=dy)
    v.label((-0.9, HALF), "$\\mathscr{I}^+$", "bl", dx=4, dy=-3)
    v.label((-HALF, -0.2), "$\\mathscr{I}^-$", "tl", dx=5, dy=3)
    v.label_xt([Q4, HALF], "$r = 0$", "b", dy=-8)
    v.label((-0.95, Q4), "the shell, $v = 0$", "bl", "small", dx=4, dy=-3)
    v.label((0, 0.45), "event horizon", "l", "small", dx=5)
    v.legend("cover", "outside the shell, $m = M$: Schwarzschild")
    v.legend("cover2", "inside the shell, $m = 0$: flat")
    v.legend("surface", "the shell $v = 0$, a light ray")
    v.legend("r", "$r$ constant outside: $0.5$, $0.8$, $1.25$, $1.6$, $2.5$ and $5\\,r_s$")
    v.legend("r2", "$r$ constant inside: $0.5$, $1$, $2$ and $4\\,r_s$")
    v.legend("t2", "$ct = cv - r$ constant inside")
    v.legend("event", "the event horizon, from the centre at $cv = -2r_s$ to $i^+$")
    v.legend("singular", "$r = 0$ after the shell arrives, where the Kretschmann scalar diverges")
    v.legend("centre", "$r = 0$ before it, a regular centre")
    v.legend("scri", "null infinity $\\mathscr{I}^\\pm$")
    v.set(input="$m(v) = 0$ for $v < 0$ and $M$ for $v > 0$, with $r_s = 2GM/c^2$.")
    return [v]


# ---------------------------------------------------------------- Tolman-Oppenheimer-Volkoff

def tov(ck, src):
    """A static star of fluid, its redshift and mass functions solved for a declared equation of
    state from this spacetime's own Einstein tensor, G = c = M_sun = 1.

    G^t_t = -8 pi rho and G^r_r = 8 pi p, read from the published components with Phi and m
    as plain symbols, are linear in m' and Phi', and give them; the pressure follows from
    the conservation law p' = -(rho + p) Phi', which the Bianchi identity makes the same
    statement as G^theta_theta = G^r_r, so the published G^theta_theta is checked along the
    solution rather than used. The equation of state is the polytrope p = K rho_0^2 with
    energy density rho = rho_0 + p, at K = 100 and central rho_0 = 1.28e-3, the star
    numerical relativity tests its codes on. Outside the surface, where p = 0, m = M and
    e^(2 Phi) = 1 - 2M/r, the Schwarzschild exterior in the same coordinates, and Phi inside
    is shifted to meet it. Then r* = int e^(-Phi) (1 - 2m/r)^(-1/2) dr runs from 0 at the
    centre to infinity and p, q = arctan((t -+ r*)/R) give Minkowski's triangle.
    """
    kappa, rho_c = 100.0, 1.28e-3
    pl = Plane(src, "tov", "spherical", ("t", "r"), EQUATOR, numeric=["Phi", "m"])
    src.note("tov", "spherical", ["einstein_tensor"])
    ul = pl.entry["einstein_tensor"]["variants"]["ul"]["nonzero"]

    def published(index):
        return pl.prep(pl.reader(next(c["value"] for c in ul if c["indices"] == [index, index])))
    Gtt, Grr, Gthth = published("t"), published("r"), published("\\theta")
    (_, F1, F2), (M0, M1, _) = pl.numeric["Phi"][2], pl.numeric["m"][2]
    rho, pres, r = sp.Symbol("rho"), sp.Symbol("p"), pl.x1
    solved = sp.solve([Gtt + 8 * sp.pi * rho, Grr - 8 * sp.pi * pres], [M1, F1], dict=True)[0]
    dm = sp.lambdify((r, M0, rho), solved[M1], "numpy")
    dPhi = sp.lambdify((r, M0, pres), solved[F1], "numpy")
    theta_theta = sp.lambdify((r, F1, F2, M0, M1), Gthth, "numpy")

    def energy(p):
        rho0 = np.sqrt(np.maximum(p, 0) / kappa)
        return rho0 + p, rho0

    def rhs(x, y):
        m, _, p = y
        e, _ = energy(p)
        f = dPhi(x, m, p)
        return [dm(x, m, e), f, -(e + max(p, 0)) * f]

    def surface(x, y):
        return y[2]
    surface.terminal = True
    p_c = kappa * rho_c ** 2
    r0 = 1e-6
    ivp = solve_ivp(rhs, (r0, 100), [4 * PI / 3 * (rho_c + p_c) * r0 ** 3, 0.0, p_c], events=surface,
                    rtol=1e-12, atol=1e-16, dense_output=True, method="DOP853")
    R = float(ivp.t_events[0][0])
    M = float(ivp.sol(R)[0])
    shift = 0.5 * np.log(1 - 2 * M / R) - float(ivp.sol(R)[1])

    def values(x):
        """(value, first, second derivative) of Phi and of m at radii x."""
        x = np.atleast_1d(np.asarray(x, dtype=float))
        inside = x < R
        out = {"Phi": [np.empty_like(x) for _ in range(3)], "m": [np.empty_like(x) for _ in range(3)]}
        if inside.any():
            xi = np.maximum(x[inside], r0)
            m, phi, p = ivp.sol(xi)
            p = np.maximum(p, 0)
            e, rho0 = energy(p)
            f = dPhi(xi, m, p)
            me = dm(xi, m, e)
            dp = -(e + p) * f
            de = -(1 + 2 * kappa * rho0) * f / (2 * kappa) + dp
            num, den = m + 4 * PI * xi ** 3 * p, xi * (xi - 2 * m)
            dnum = me + 12 * PI * xi ** 2 * p + 4 * PI * xi ** 3 * dp
            dden = 2 * xi - 2 * m - 2 * me * xi
            for store, v in zip(out["Phi"], (phi + shift, f, (dnum * den - num * dden) / den ** 2)):
                store[inside] = v
            for store, v in zip(out["m"], (m, me, 8 * PI * xi * e + 4 * PI * xi ** 2 * de)):
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

    def fvals(t, x):
        v = values(x)
        return {"Phi": v["Phi"], "m": v["m"]}
    grid_r = np.linspace(0, R, 200001)
    v = values(grid_r)
    speed = np.exp(-v["Phi"][0]) / np.sqrt(1 - 2 * v["m"][0] / np.maximum(grid_r, r0))
    inner_rs = cumulative_trapezoid(speed, grid_r, initial=0)
    outer_c = inner_rs[-1] - (R + 2 * M * np.log(R / (2 * M) - 1))

    def rstar(x):
        x = np.asarray(x, dtype=float)
        with np.errstate(invalid="ignore", divide="ignore"):
            outside = x + 2 * M * np.log(x / (2 * M) - 1) + outer_c
        return np.where(x <= R, np.interp(x, grid_r, inner_rs), outside)

    def star(t, x):
        return mink_pq(t, rstar(x), R)
    ck.chart("Tolman-Oppenheimer-Volkoff, the declared star", pl, star, ck.uniform(-60, 60),
             ck.uniform(0.01, R), lambda t, x: (1, 0), fvals)
    ck.chart("Tolman-Oppenheimer-Volkoff, outside it", pl, star, ck.uniform(-60, 60),
             ck.uniform(R, 80), lambda t, x: (1, 0), fvals)
    x = np.linspace(0.05 * R, 0.95 * R, 200)
    vi = values(x)
    p = np.maximum(ivp.sol(x)[2], 0)
    ck.limit("Tolman-Oppenheimer-Volkoff: the declared star solves the published G^theta_theta = 8 pi p",
             (theta_theta(x, vi["Phi"][1], vi["Phi"][2], vi["m"][0], vi["m"][1]) - 8 * PI * p) / (8 * PI * p_c),
             0, 1e-7)
    ck.limit("Tolman-Oppenheimer-Volkoff: the surface clears Buchdahl's bound, 2M/R < 8/9",
             [float(2 * M / R < 8 / 9)], [1], 0.5)
    ck.limit("Tolman-Oppenheimer-Volkoff: Phi is continuous at the surface",
             [values(R * (1 - 1e-9))["Phi"][0][0]], [values(R * (1 + 1e-9))["Phi"][0][0]], 1e-8)
    ck.finite("Tolman-Oppenheimer-Volkoff: the centre is regular",
              pl.kretschmann(np.zeros(5), np.array([1e-4, 1e-3, 1e-2, 0.5, 1.0]), fvals))

    km = 1.4766250614  # GM_sun/c^2 in km
    v = View("spherical", "Spherical", [-0.35, PI + 0.35, -PI - 0.25, PI + 0.25], "spherical")
    v.fill("region", TRIANGLE)
    v.fill("cover", TRIANGLE)
    ps, qs = mink_pq(S_ALL, rstar(R), R)
    v.fill("star", [[0, -PI]] + [point(a, b) for a, b in zip(ps, qs)] + [[0, PI]])
    inner = [round(R / 3, 1), round(2 * R / 3, 1)]
    outer = [2 * round(R, 1), 4 * round(R, 1)]
    for x in inner:
        v.curve("r", *star(S_ALL, np.full_like(S_ALL, x)))
    for x in outer:
        v.curve("r2", *star(S_ALL, np.full_like(S_ALL, x)))
    rr = np.concatenate([np.linspace(0, R, 60)[:-1], R + np.exp(np.linspace(-6, 9, 200)) - np.exp(-6)])
    times = (-4, -2, -1, 0, 1, 2, 4)
    for tt in times:
        v.curve("t", *star(np.full_like(rr, tt * R), rr))
    v.curve("surface", ps, qs)
    triangle_edges(v)
    label_on(v, mink_pq(0, rstar(R), R), "$r = R$")
    v.label_xt([0.33, 0.0], "star", cls="region")
    v.legend("star", "the star, where $p > 0$")
    v.legend("cover", "the whole spacetime, which $t$ and $r$ cover")
    v.legend("r", f"$r$ constant inside, at {inner[0] * km:.1f} and {inner[1] * km:.1f} km")
    v.legend("r2", f"$r$ constant outside, at {outer[0] * km:.1f} and {outer[1] * km:.1f} km")
    v.legend("t", "$ct$ constant, at $0$, $\\pm R$, $\\pm 2R$ and $\\pm 4R$")
    v.legend("surface", "the surface $r = R$, where the pressure falls to zero")
    v.legend("centre", "$r = 0$, a regular centre")
    v.set(input="A polytrope, $p = K\\rho_0^2$ with rest mass density $\\rho_0$ and energy density "
                "$\\rho c^2 = \\rho_0c^2 + p$, at $K = 100$ and a central $\\rho_0 = 1.28\\times10^{-3}$ in "
                "units where $G = c = M_\\odot = 1$, solved from this spacetime's own $G^t{}_t$ and "
                f"$G^r{{}}_r$: a star of $M = {M:.2f}\\,M_\\odot$ and $R = {R * km:.1f}$ km, "
                "the one numerical relativity tests its codes on.")
    return [v]


# ---------------------------------------------------------------- the Malament-Hogarth toy

def malament_hogarth(ck, src):
    """Minkowski space less its origin, times Omega^2. A conformal factor changes no null
    direction, so for every Omega the causal structure is Minkowski's less that event, which
    is symmetric about the t axis through it: the plane of t and x at y = z = 0 with x > 0 is
    every half plane through the axis, and p, q = arctan((ct -+ r)/l), with r the distance
    from the axis, give Minkowski's triangle less one point of its axis. The check hands the
    metric a different random Omega at every sample point, which is how little it matters.
    """
    pl = Plane(src, "malament_hogarth", "cartesian", ("t", "x"), {"y": "0", "z": "0"}, numeric=["Omega"])

    def any_omega(t, x):
        w = ck.rng.uniform(0.2, 50, np.shape(t))
        return {"Omega": (w, 0 * w, 0 * w)}
    ck.chart("Malament-Hogarth, for any conformal factor", pl, mink_pq, ck.uniform(-20, 20),
             ck.uniform(0.01, 20), lambda t, x: (1, 0), any_omega)
    ck.limit("Malament-Hogarth: the removed event lands on the axis, (X, T) = (0, 0)", point(*mink_pq(0, 0)), [0, 0])

    v = View("cartesian", "Cartesian", [-0.35, PI + 0.35, -PI - 0.25, PI + 0.25], "cartesian")
    v.fill("region", TRIANGLE)
    v.fill("cover", TRIANGLE)
    event = mink_pq(1, 0)
    past = [point(*event), point(-HALF, float(event[1])), [0, -PI]]
    v.fill("past", past)
    grid(v, "r", lambda x, t: mink_pq(t, x), (0.5, 1, 2, 4), S_ALL)
    grid(v, "t", mink_pq, (-4, -2, -1, 1, 2, 4), S_POS)
    v.line("world", [[[0, -PI], [0, 0]]])
    v.segment("cone", event, (-HALF, float(event[1])))
    v.line("scri", [[[0, PI], [PI, 0]], [[PI, 0], [0, -PI]]])
    v.line("centre", [[[0, 0], [0, PI]]])
    for at, text, anchor, dx, dy in (((PI, 0), "$i^0$", "l", 6, 0), ((0, PI), "$i^+$", "b", 0, -6),
                                     ((0, -PI), "$i^-$", "t", 0, 6)):
        v.layers.append({"kind": "point", "class": "infinity", "at": rounded(at)})
        v.label_xt(at, text, anchor, dx=dx, dy=dy)
    v.label_xt([HALF, HALF], "$\\mathscr{I}^+$", "bl", dx=5, dy=-3)
    v.label_xt([HALF, -HALF], "$\\mathscr{I}^-$", "tl", dx=5, dy=3)
    v.point("mark", event)
    v.point("removed", (0, 0))
    v.label(event, "$p$", "r", "small", dx=-6)
    v.label((0, 0), "the removed event", "l", "small", dx=8)
    v.label_xt([0, -1.6], "the computer", "l", "small", dx=8)
    v.legend("cover", "the region that $t$, $x$, $y$ and $z$ cover: all but the removed event")
    v.legend("past", "the causal past of $p$, which holds the whole world line")
    v.legend("r", "$r = \\sqrt{x^2 + y^2 + z^2}$ constant, in units of $\\ell$")
    v.legend("t", "$ct$ constant")
    v.legend("world", "the computer's world line, $x = y = z = 0$ with $t < 0$")
    v.legend("cone", "the edge of $p$'s past light cone")
    v.legend("mark", "the event $p$ at $(ct, r) = (\\ell, 0)$")
    v.legend("removed", "the removed event, the origin")
    v.legend("centre", "the axis $r = 0$ above it")
    v.legend("scri", "null infinity $\\mathscr{I}^\\pm$")
    return [v]


# ---------------------------------------------------------------- the table

DRAWN = {
    "minkowski": minkowski, "schwarzschild": schwarzschild, "rn_metric": reissner_nordstrom,
    "kerr": kerr, "kerr_newman": kerr_newman, "de_sitter": de_sitter, "anti_de_sitter": anti_de_sitter,
    "bertotti_robinson": bertotti_robinson, "ellis_bronnikov": ellis_bronnikov, "morris_thorne": morris_thorne,
    "cosmic_string": cosmic_string, "interior_schwarzschild": interior_schwarzschild, "frw": frw,
    "oppenheimer_snyder": oppenheimer_snyder, "vaidya": vaidya, "tov": tov,
    "malament_hogarth": malament_hogarth,
}

# ---------------------------------------------------------------- the captions

# Each view's caption, prose under the rules of _tools/README.md: no dashes but in a name,
# and every sentence about the spacetime, never about the page or the collection. A caption
# opens by naming what is drawn, the whole spacetime or the surface in it.
CAPTIONS = {
    ("minkowski", "spherical"): [
        "This is the whole of Minkowski spacetime, and each point of the diagram stands for a "
        "sphere of radius $r$. With $u = ct - r$ and $v = ct + r$ the metric on the plane of $t$ "
        "and $r$ is $-du\\,dv$, and for any length $\\ell$ the maps $p = \\arctan(u/\\ell)$ and "
        "$q = \\arctan(v/\\ell)$ bring all of it into a finite triangle, drawn with $T = p + q$ up "
        "and $X = q - p$ across. Lines of constant $p$ or $q$ are light rays, and they run at "
        "45°.",
        "Every line of constant $r$ runs from $i^-$ to $i^+$, and every line of constant $t$ from "
        "the centre to $i^0$. Ingoing light starts on $\\mathscr{I}^-$, passes through the centre "
        "and ends on $\\mathscr{I}^+$.",
    ],
    ("minkowski", "spherical_null"): [
        "This is the whole of Minkowski spacetime in its spherical null coordinates, "
        "$u = t - r/c$ and $v = t + r/c$, and each point of the diagram stands for a sphere of "
        "radius $c(v - u)/2$. The lines of constant $u$ and of constant $v$ are light rays, the "
        "45° lines of the triangle, with $p = \\arctan(cu/\\ell)$ and $q = \\arctan(cv/\\ell)$, and "
        "the centre is the line $u = v$.",
    ],
    ("minkowski", "cartesian"): [
        "This is the plane $y = z = 0$, which is flat and totally geodesic, brought by "
        "$p, q = \\arctan((ct \\mp x)/\\ell)$ into the whole diamond. It has two ends, "
        "$x \\to +\\infty$ and $x \\to -\\infty$, each with its own null infinity.",
        "Turned about the line $x = 0$, each half of the diamond sweeps out the spherical "
        "triangle, which is the whole spacetime.",
    ],
    ("minkowski", "double_null"): [
        "This is the plane $y = z = 0$ in the coordinates $u = t - x/c$ and $v = t + x/c$, in "
        "which the metric on the plane is $-c^2\\,du\\,dv$. The lines of constant $u$ and of "
        "constant $v$ are light rays, the 45° lines of the diamond, with $p = \\arctan(cu/\\ell)$ "
        "and $q = \\arctan(cv/\\ell)$.",
    ],
    ("minkowski", "rindler"): [
        "This is the plane $Y = Z = 0$ in Rindler's coordinates, $ct = X\\sinh(aT/c)$ and "
        "$x = X\\cosh(aT/c)$, so that $ct - x = -Xe^{-aT/c}$ and $ct + x = Xe^{aT/c}$. They cover "
        "the wedge $x > c|t|$, and $X = 0$, where $g_{TT}$ vanishes, is the pair of null lines "
        "through the origin.",
        "An observer at constant $X$ accelerates uniformly, at $c^2/X$, and the null line "
        "$ct = x$ is that observer's horizon: no event beyond it can send a signal into the "
        "wedge, as nothing inside $r_s$ can reach a static observer outside a black hole.",
    ],
    ("schwarzschild", "spherical"): [
        "This is the whole of the Schwarzschild spacetime, maximally extended, and each point of "
        "the diagram stands for a sphere of radius $r$. Kruskal and Szekeres's "
        "$U = -e^{-u/2r_s}$ and $V = e^{v/2r_s}$, with $u, v = ct \\mp r_*$ and "
        "$r_* = r + r_s\\ln|r/r_s - 1|$, make the metric regular through $r = r_s$, where "
        "$UV = (1 - r/r_s)e^{r/r_s}$ vanishes. With $p = \\arctan U$ and $q = \\arctan V$ the "
        "singularity $UV = 1$ lies exactly on the straight lines $T = \\pm\\pi/2$, since "
        "$\\tan(p + q) = (U + V)/(1 - UV)$ diverges there.",
        "The coordinates $t$ and $r > r_s$ cover the right exterior alone. The horizon is the "
        "pair of null lines $U = 0$ and $V = 0$, crossing at the bifurcation sphere. The black "
        "hole above it ends at $r = 0$ on $T = \\pi/2$, the white hole below it begins at "
        "$r = 0$ on $T = -\\pi/2$, and both singularities are spacelike.",
    ],
    ("schwarzschild", "ingoing"): [
        "This is the whole Schwarzschild spacetime with the ingoing Eddington-Finkelstein "
        "coordinates $v$ and $r$ on it. From $V = e^{v/2r_s}$ and $U = (1 - r/r_s)e^{r/r_s}/V$, "
        "one formula for every $r > 0$, they cover the exterior and the black hole together, "
        "and their lines of constant $v$ are ingoing light rays, which cross the horizon at "
        "45° and end at $r = 0$.",
    ],
    ("schwarzschild", "outgoing"): [
        "This is the whole Schwarzschild spacetime with the outgoing Eddington-Finkelstein "
        "coordinates $u$ and $r$ on it, the time reverse of the ingoing ones. From "
        "$U = -e^{-u/2r_s}$ and $V = (r/r_s - 1)e^{r/r_s}/(-U)$ they cover the exterior and the "
        "white hole, and their lines of constant $u$ are outgoing light rays, which leave $r = 0$ "
        "and cross the horizon outward.",
    ],
    ("rn_metric", "tower"): [
        "This is the Reissner-Nordström spacetime maximally extended, and each point of the "
        "diagram stands for a sphere of radius $r$. The extension is a tower of regions that "
        "repeats up and down without end. Its tortoise coordinate is "
        "$r_* = r + \\frac{1}{2\\kappa_+}\\ln|r/r_+ - 1| - \\frac{1}{2\\kappa_-}\\ln|r/r_- - 1|$ "
        "with $\\kappa_\\pm = (r_+ - r_-)/2r_\\pm^2$, and each logarithm is made dimensionless by "
        "its own root, so that $r_*(0) = 0$.",
        "Every region is placed by the Kruskal coordinate of the outer horizon, "
        "$p = \\pm\\arctan e^{-\\kappa_+ u}$ and $q = \\pm\\arctan e^{\\kappa_+ v}$ with "
        "$u, v = t \\mp r_*$, and the regions above the inner horizon are the reflection "
        "$(p, q) \\to (\\pi - q, \\pi - p)$ of those below. This map is smooth across $r_+$ and puts "
        "the singularity $r = 0$ exactly on the vertical lines $X = \\pm\\pi/2$, where it is "
        "timelike. Across $r_-$ it is continuous and cannot also be smooth, because the late "
        "light rays that reach $\\mathscr{I}^+$ are the rays that pile up at the Cauchy horizon "
        "$r_-$, and one function of the ray has to serve both.",
        "The coordinates $t$ and $r > 0$ cover one region of each kind: an exterior, a region "
        "between the horizons and a region inside $r_-$. Between the horizons their $t$ alone "
        "cannot tell the black hole from the white hole, and the region is taken to be the "
        "black hole an infalling observer enters.",
    ],
    ("rn_metric", "malament_hogarth"): [
        "This is the same tower with one event beyond the Cauchy horizon $r_-$ marked on it. "
        "Every point of the exterior below it has $p \\le p_e$ and $q \\le q_e$, so the whole "
        "exterior lies in the event's causal past. A static observer at $r = 1.2\\,r_s$ lives "
        "from $i^-$ to $i^+$ for an infinite proper time, and every moment of that life can send "
        "a signal that reaches the event, arriving at the Cauchy horizon infinitely "
        "blueshifted.",
    ],
    ("kerr", "axis"): [
        "This is the symmetry axis $\\theta = 0$ of the maximally extended Kerr spacetime, the "
        "surface the rotations leave fixed and so totally geodesic. On it the metric is "
        "$-\\frac{\\Delta}{r^2 + a^2}c^2dt^2 + \\frac{r^2 + a^2}{\\Delta}dr^2$ with "
        "$\\Delta = r^2 - 2GMr/c^2 + a^2$, and its two simple roots give it the tower of "
        "Reissner-Nordström. Carter extended the axis this way in 1966.",
        "Where Reissner-Nordström ends at $r = 0$, the axis runs on through the centre of the "
        "ring's disc, where the curvature is finite, into $r < 0$, a second asymptotically flat "
        "end with its own null infinity. The ring singularity itself is at $r = 0$ in the "
        "equatorial plane $\\theta = \\pi/2$. The coordinates $t$ and $r > r_+$ cover the "
        "exterior.",
    ],
    ("kerr_newman", "axis"): [
        "This is the symmetry axis $\\theta = 0$ of the maximally extended Kerr-Newman spacetime, "
        "totally geodesic as Kerr's is. On it the metric is again "
        "$-\\frac{\\Delta}{r^2 + a^2}c^2dt^2 + \\frac{r^2 + a^2}{\\Delta}dr^2$, now with "
        "$\\Delta = r^2 - 2GMr/c^2 + a^2 + r_Q^2$, and its two simple roots give the same tower.",
        "The axis runs through the centre of the ring's disc into $r < 0$ as Kerr's does, and "
        "the ring singularity is at $r = 0$ in the equatorial plane $\\theta = \\pi/2$. The "
        "coordinates $t$ and $r > r_+$ cover the exterior.",
    ],
    ("de_sitter", "static"): [
        "This is the whole of de Sitter spacetime, the hyperboloid "
        "$-X_0^2 + X_1^2 + \\dots + X_4^2 = L^2$ with $L = \\sqrt{3/\\Lambda}$, and each point of "
        "the diagram stands for a sphere. Its global coordinates give the metric "
        "$\\frac{L^2}{\\cos^2 T}(-dT^2 + d\\chi^2 + \\sin^2\\chi\\,d\\Omega^2)$ on the square "
        "$|T| < \\pi/2$, $0 \\le \\chi \\le \\pi$, with $X = \\chi$ across.",
        "The static coordinates enter the square as $\\tan p = \\tanh(u/2L)$ and "
        "$\\tan q = \\tanh(v/2L)$, with $u, v = ct \\mp L\\,\\mathrm{artanh}(r/L)$, and cover the "
        "triangle about the observer at $\\chi = 0$; their horizon $r = L$ is the pair of null "
        "lines through the centre of the square. Infinity is spacelike, past and future, so every "
        "observer has an event horizon: nothing beyond the line from $(\\chi, T) = (0, \\pi/2)$ to "
        "$(\\pi, -\\pi/2)$ ever reaches the observer at $\\chi = 0$.",
    ],
    ("de_sitter", "flat"): [
        "This is the whole of de Sitter spacetime with the flat slicing on it. Its conformal "
        "time $\\eta = -e^{-Ht}/H$ makes the metric "
        "$\\frac{1}{H^2\\eta^2}(-c^2d\\eta^2 + d\\rho^2 + \\rho^2d\\Omega^2)$, conformal to half of "
        "Minkowski space, which enters the square as $p = \\pi/4 + \\arctan(H\\eta - H\\rho/c)$ and "
        "$q = \\pi/4 + \\arctan(H\\eta + H\\rho/c)$, with $\\rho^2 = x^2 + y^2 + z^2$.",
        "The slicing covers the half above the observer's past horizon, so its $t \\to -\\infty$ "
        "is a null line where the coordinates end, and the spacetime goes on below it.",
    ],
    ("anti_de_sitter", "global"): [
        "This is the whole of anti-de Sitter spacetime, its universal cover, and each point of "
        "the diagram stands for a sphere. With $\\sigma = \\arctan(r/L)$ the metric on the plane "
        "of $t$ and $r$ is $\\frac{1}{\\cos^2\\sigma}(-c^2dt^2 + L^2d\\sigma^2)$, already conformal "
        "to the strip $0 \\le \\sigma < \\pi/2$, which is unbounded in $t$. Its edge "
        "$\\sigma = \\pi/2$ is a timelike boundary.",
        "A radial light ray from the centre reaches the boundary at $ct = \\pi L/2$ and is back "
        "at $ct = \\pi L$. The radial timelike geodesics, $\\sin\\sigma = k\\sin(ct/L)$ with "
        "$k < 1$, all return to the centre at the same $ct = \\pi L$, whatever their energy.",
    ],
    ("anti_de_sitter", "poincare"): [
        "This is the plane $x = y = 0$ of the Poincaré patch, which passes through the centre and "
        "is conformal to the strip $-\\pi/2 < \\sigma < \\pi/2$. The coordinates $t$ and $z$ "
        "enter the strip as "
        "$p = -\\pi/4 + \\arctan((ct + z)/L)$ and $q = \\pi/4 + \\arctan((ct - z)/L)$.",
        "They cover a wedge of it. Their $z \\to 0$ is the conformal boundary, and "
        "$z \\to \\infty$ is the Poincaré horizon, the pair of null lines from "
        "$(\\sigma, ct/L) = (-\\pi/2, 0)$. The curvature there is the same as everywhere else, "
        "and the global coordinates run smoothly across it.",
    ],
    ("bertotti_robinson", "static"): [
        "This is the whole of the Bertotti-Robinson spacetime, the product of an anti-de Sitter "
        "space of two dimensions and radius $b$ with a sphere of radius $b$. Its conformal "
        "diagram is the strip of the first factor, and since the sphere has the same radius "
        "everywhere, each point of the strip stands for a sphere of radius $b$.",
        "The throat coordinates $t$ and $r$ cover a Poincaré wedge of the strip. Under "
        "$x = b^2/r$ their $-\\frac{r^2}{b^2}c^2dt^2 + \\frac{b^2}{r^2}dr^2$ becomes "
        "$\\frac{b^2}{x^2}(-c^2dt^2 + dx^2)$, so the throat's $r = 0$ is the Poincaré horizon, a "
        "horizon of the coordinates alone, across which the spacetime continues.",
    ],
    ("bertotti_robinson", "poincare"): [
        "This is the whole of the Bertotti-Robinson spacetime with the Poincaré coordinates $t$ "
        "and $x$ on it. They cover the same wedge as the throat coordinates, with "
        "$x = b^2/r$.",
    ],
    ("ellis_bronnikov", "spherical"): [
        "This is the whole of the Ellis-Bronnikov wormhole, and each point of the diagram stands "
        "for a sphere of area $4\\pi(r^2 + \\ell^2)$. The metric on the plane of $t$ and $r$ is "
        "$-c^2dt^2 + dr^2$ with $r$ over the whole line, which is Minkowski space of two "
        "dimensions, and $p, q = \\arctan((ct \\mp r)/\\ell)$ bring it into the full diamond.",
        "The two ends, $r \\to +\\infty$ and $r \\to -\\infty$, are two asymptotically flat "
        "universes, each with its own $i^0$ and $\\mathscr{I}^\\pm$, joined at the throat $r = 0$, "
        "where the spheres are smallest. Light crosses the throat at 45°, as it does everywhere "
        "else, so the wormhole has no horizon.",
    ],
    ("morris_thorne", "spherical"): [
        "This is the whole of the Morris-Thorne wormhole with $\\Phi = 0$ and $b = b_0^2/r$, and "
        "each point of the diagram stands for a sphere of radius $r$. The proper radial distance "
        "is $l = \\pm\\sqrt{r^2 - b_0^2}$, and $p, q = \\arctan((ct \\mp l)/b_0)$ bring the plane "
        "of $t$ and $l$ into the full diamond, with the throat $r = b_0$ on its axis.",
        "Every $\\Phi$ and $b$ that give no horizon and two flat ends give the same diamond. With "
        "$\\Phi$ bounded and tending to a constant and $b/r \\to 0$ at both ends, "
        "$\\int e^{-\\Phi}\\,dl$ runs over the whole line, and another choice moves only the "
        "surfaces of constant $t$ and $r$ inside it. The areal coordinates $t$ and $r$ cover one "
        "side and end at the throat, where $g_{rr}$ diverges.",
    ],
    ("morris_thorne", "proper_radial"): [
        "This is the same wormhole in the proper distance coordinates $t$ and $l$, which cover "
        "both sides and run smoothly through the throat at $l = 0$. With $\\Phi = 0$ the metric "
        "on the plane is $-c^2dt^2 + dl^2$, and $p, q = \\arctan((ct \\mp l)/b_0)$ bring it into "
        "the diamond.",
    ],
    ("cosmic_string", "conical"): [
        "This is the half plane of $t$ and $r$ at fixed $\\phi$ and $z$, which is totally "
        "geodesic. The metric on it is $-c^2dt^2 + dr^2$, so "
        "$p, q = \\arctan((ct \\mp r)/\\ell)$ bring it into Minkowski's half diamond, with the "
        "string at $r = 0$ in place of a regular centre.",
        "The string's gravity is its deficit angle $\\delta = 8\\pi G\\mu/c^2$, which shows in "
        "the circles of constant $r$ around it: each has circumference $2\\pi(1 - 4G\\mu/c^2)\\,r$, "
        "short of $2\\pi r$ by $\\delta r$.",
    ],
    ("cosmic_string", "gott"): [
        "This is the same half plane with Gott's core, which makes the axis regular. Inside the "
        "core the proper distance from the axis is $\\rho = \\ell\\chi$, and outside it is "
        "$\\ell\\chi_0 + r - \\ell\\tan\\chi_0$; the edge of the core, $r = \\ell\\tan\\chi_0$, is "
        "where the circumferences inside and outside agree. In $\\rho$ the metric on the whole "
        "half plane is $-c^2dt^2 + d\\rho^2$.",
    ],
    ("interior_schwarzschild", "spherical"): [
        "This is the whole of a static star of uniform density, the interior solution for "
        "$r \\le R$ joined at $R = 1.5\\,r_s$ to the Schwarzschild exterior, and each point of the "
        "diagram stands for a sphere of radius $r$. That radius clears Buchdahl's bound, "
        "$R > \\frac{9}{8}r_s$, so there is no horizon.",
        "Both sides give $g_{tt} = -(1 - r_s/R)$ at the surface, so $t$ is one coordinate "
        "throughout. The tortoise coordinate $r_* = \\int\\sqrt{g_{rr}/(-g_{tt})}\\,dr$ runs from "
        "the centre through the surface, and $p, q = \\arctan((t \\mp r_*)/R)$ bring the "
        "spacetime into Minkowski's triangle, the causal structure of empty space, with the star "
        "a timelike tube from $i^-$ to $i^+$.",
    ],
    ("frw", "flat"): [
        "This is the whole of a flat universe of dust, and each point of the diagram stands for "
        "a sphere. With $k = 0$, $G^r{}_r = 0$ gives $a \\propto \\eta^2$, and the metric "
        "$a^2(-d\\eta^2 + dr^2 + r^2d\\Omega^2)$ is conformal to the half $\\eta > 0$ of Minkowski "
        "space, which $p, q = \\arctan((\\eta \\mp r)/\\eta_0)$ bring into a triangle, with "
        "$\\eta_0$ the conformal time today.",
        "The big bang is the straight line $T = 0$ and is spacelike, while future infinity is "
        "null, as Minkowski's is. Our past light cone meets the bang at the comoving radius "
        "$r = \\eta_0$, the particle horizon, and light from anything farther away has not "
        "reached us yet.",
    ],
    ("frw", "closed"): [
        "This is the whole of a closed universe of dust, and each point of the diagram stands "
        "for a sphere. With $k = +1$, dust gives $a \\propto 1 - \\cos\\eta$, with $\\eta$ from "
        "$0$ to $2\\pi$, and with $r = \\sin\\chi$ the metric is already conformal to the Einstein "
        "static universe, the rectangle $0 \\le \\chi \\le \\pi$ with the bang along its bottom "
        "and the crunch along its top.",
        "A light ray that leaves $\\chi = 0$ at the bang reaches the antipode $\\chi = \\pi$ at "
        "maximum expansion and is back at the crunch. The radius $r = \\sin\\chi$ covers one "
        "hemisphere, $\\chi < \\pi/2$, and ends at the equator $r = 1$, where $1 - kr^2$ "
        "vanishes.",
    ],
    ("frw", "open"): [
        "This is the whole of an open universe of dust, and each point of the diagram stands for "
        "a sphere. With $k = -1$, dust gives $a \\propto \\cosh\\eta - 1$, and with "
        "$r = \\sinh\\chi$ the map $\\tan((T \\pm X)/2) = \\tanh((\\eta \\pm \\chi)/2)$ sends it into "
        "the Einstein static universe. It has the causal structure of the flat universe, a "
        "triangle with the bang along its base and null infinity above, and differs from it only "
        "in where its surfaces of constant $\\eta$ and $\\chi$ lie.",
    ],
    ("oppenheimer_snyder", "collapse"): [
        "This is the whole of a ball of dust collapsing from rest at $R_0 = 2\\,r_s$, and each "
        "point of the diagram stands for a sphere. Inside, the dust is a closed universe, "
        "$a^2(-d\\eta^2 + d\\chi^2 + \\sin^2\\chi\\,d\\Omega^2)$ with $a = \\frac{a_m}{2}(1 + \\cos\\eta)$, "
        "already conformal to the Einstein static universe in $\\eta$ and $\\chi$. Outside, "
        "$p = P(U)$ and $q = Q(V)$ are functions of the Kruskal coordinates, and three conditions "
        "fix them completely: the two sides agree on the surface $\\chi = \\chi_0$, a radial "
        "geodesic of the exterior with energy $\\cos\\chi_0$; the moment of rest is the line "
        "$T = 0$, the exterior's symmetry $U \\leftrightarrow -V$; and $r = 0$ is the line "
        "$T = \\pi$, where $UV = 1$.",
        "The event horizon $U = 0$ enters the dust as the outgoing light ray from the centre at "
        "$\\eta = \\pi - 3\\chi_0$, before the surface crosses $r_s$ at $\\eta = \\pi - 2\\chi_0$. "
        "Marginally trapped spheres appear at the surface then and move inward along "
        "$\\eta = \\pi - 2\\chi$, a timelike curve, and the crunch of the dust and the singularity "
        "outside it are one spacelike line.",
    ],
    ("tov", "spherical"): [
        "This is the whole of a static star of fluid with a polytrope for its equation of state, "
        "and each point of the diagram stands for a sphere of radius $r$. Its mass and redshift "
        "functions come from $G^t{}_t = -8\\pi G\\rho/c^2$ and $G^r{}_r = 8\\pi Gp/c^4$ and its "
        "pressure from $\\partial_r p = -(\\rho c^2 + p)\\,\\partial_r\\Phi$. Where the pressure falls "
        "to zero, the same coordinates carry on as Schwarzschild's exterior, with $m = GM/c^2$ and "
        "$e^{2\\Phi} = 1 - 2m/r$.",
        "The tortoise coordinate $r_* = \\int e^{-\\Phi}(1 - 2m/r)^{-1/2}\\,dr$ runs from $0$ at the "
        "centre to infinity, and $p, q = \\arctan((ct \\mp r_*)/R)$ bring the spacetime into "
        "Minkowski's triangle, with the star a timelike tube from $i^-$ to $i^+$. Every such star "
        "has this causal structure: by Buchdahl's theorem a static ball of fluid whose density "
        "does not grow outward has $2GM/c^2R \\le 8/9$, so it has no horizon, and another "
        "equation of state moves only the surfaces of constant $t$ and $r$ inside the triangle.",
    ],
    ("malament_hogarth", "cartesian"): [
        "This is the whole of the Malament-Hogarth toy spacetime, Minkowski space with one event "
        "removed and its metric multiplied by $\\Omega^2$. A conformal factor changes no null "
        "direction, so for every $\\Omega$ the causal structure is Minkowski's less that event. "
        "It is symmetric about the $t$ axis through the event, so each point of the triangle is "
        "the sphere of events at one $t$ and one distance $r = \\sqrt{x^2 + y^2 + z^2}$ from the "
        "axis, and $p, q = \\arctan((ct \\mp r)/\\ell)$ bring it into Minkowski's triangle with one "
        "point of its axis removed.",
        "The computer's world line runs up the axis into the removed event, and every event above "
        "it, such as $p$ at $(ct, r) = (\\ell, 0)$, has the whole of that world line in its causal "
        "past, since a signal can pass round the missing point. Where $\\Omega$ grows at least as "
        "fast as $1/|t|$ along the axis, the world line's proper time $\\int\\Omega\\,dt$ is "
        "infinite while every light ray runs as it does in Minkowski space. The computer then "
        "runs for an infinite proper time before the removed event, and its signals reach $p$ "
        "blueshifted by $\\Omega$, without bound.",
    ],
    ("vaidya", "shell"): [
        "This is the whole of a spacetime into which a spherical shell of null dust of mass $M$ "
        "falls along $v = 0$, and each point of the diagram stands for a sphere. With $m = 0$ for "
        "$v < 0$ the metric inside the shell is flat, and with $m = M$ for $v > 0$ it is "
        "Schwarzschild's in ingoing coordinates, placed outside the shell by $p = \\arctan U$ and "
        "$q = \\arctan V$. Inside, each outgoing light ray keeps the $p$ it has where it crosses "
        "the shell, $p = \\arctan((1 + cu/2r_s)e^{-cu/2r_s})$ with $u = v - 2r/c$ the flat "
        "retarded time, and $q$ is the same function of $v$, which puts the centre on the "
        "straight line $X = 0$.",
        "The event horizon forms at the centre at $cv = -2r_s$, before the shell arrives, and "
        "grows through flat space to meet the shell at $r = r_s$. The lines inside crowd toward "
        "the shell because $q$, chosen to make the centre straight, has zero slope there.",
    ],
}


def draw(metric_id, ck):
    src = Sources()
    if metric_id in NOT_DRAWN:
        metric = json.loads((build.METRICS_DIR / f"{metric_id}.json").read_text(encoding="utf-8"))
        for system in metric["coordinates"]:
            src.note(metric_id, system["id"], NONE_FIELDS)
        return {"metric": metric_id, "source": src.stamps(), "none": NOT_DRAWN[metric_id]}
    views = DRAWN[metric_id](ck, src)
    out = []
    for v in views:
        v.set(caption=CAPTIONS[(metric_id, v.d["id"])])
        out.append(v.done())
    return {"metric": metric_id, "source": src.stamps(), "views": out}


def check_table():
    """Every view has a caption, and every spacetime is either drawn or says why not."""
    ids = {p.stem for p in build.METRICS_DIR.glob("*.json") if not build.CONFLICT_COPY.search(p.stem)}
    both = set(DRAWN) & set(NOT_DRAWN)
    neither = ids - set(DRAWN) - set(NOT_DRAWN)
    if both or neither:
        raise SystemExit(f"drawn and not drawn: {sorted(both)}; neither: {sorted(neither)}")
    stray = {metric for metric, _ in CAPTIONS} - set(DRAWN)
    if stray:
        raise SystemExit(f"captions for spacetimes that are not drawn: {sorted(stray)}")


def main(argv=None):
    check_table()
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--metric", action="append", default=[], help="redraw only this spacetime, repeatable")
    parser.add_argument("--verify", action="store_true", help="run and print every check instead of writing")
    args = parser.parse_args(argv)
    unknown = set(args.metric) - set(DRAWN) - set(NOT_DRAWN)
    if unknown:
        parser.error(f"no conformal diagram is drawn for {sorted(unknown)}")
    wanted = [m for m in list(DRAWN) + list(NOT_DRAWN) if not args.metric or m in args.metric]
    start = time.time()
    ck = Checks()
    # A compactification sends the ends of every coordinate line to infinity on purpose:
    # exp overflows to inf there and arctan(inf) is exactly the boundary, pi/2. Anything
    # wrong that such an overflow could hide is what the checks look for.
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        files = {metric_id: draw(metric_id, ck) for metric_id in wanted}
    if args.verify:
        ck.report()
    failed = ck.failures()
    print(f"{len(ck.charts)} charts and {len(ck.limits)} limits checked in {time.time() - start:.1f} s, "
          f"{len(failed)} failed", file=sys.stderr if failed else sys.stdout)
    for name in failed:
        print(f"FAILED {name}", file=sys.stderr)
    if failed or args.verify:
        return 1 if failed else 0
    CONFORMAL_DIR.mkdir(parents=True, exist_ok=True)
    for metric_id, data in files.items():
        path = CONFORMAL_DIR / f"{metric_id}.json"
        path.write_text(json.dumps(data, ensure_ascii=False, separators=(",", ":")) + "\n", encoding="utf-8")
        print(f"wrote {path.relative_to(build.ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
