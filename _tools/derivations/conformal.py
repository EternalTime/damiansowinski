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

DRAWN lists the thirty-one spacetimes that have a diagram and NOT_DRAWN the others, which
have no file: a full redraw removes one left behind. The script stops if a metric file is
in neither, so a new spacetime needs a decision.


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
  fade            {top, bottom}: how far the drawing fades out where it continues;
  slices          each moment of the spacetime's embedding diagram that the view shows, as a
                  spacetime diagram's view carries it, its lines and points in (X, T), drawn
                  by the view's own maps over the part of the moment the embedding reaches;
                  slices.py reads the moments and their reach.
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

import slices

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))
import build_mfs_data as build  # noqa: E402
import null_rays as nr  # noqa: E402
import verify_metrics as vm  # noqa: E402

CONFORMAL_DIR = build.CONFORMAL_DIR
BASE_FIELDS = ("coords", "parameters", "metric_components", "inverse_metric_components")

PI = np.pi
HALF = PI / 2
Q4 = PI / 4
EQUATOR = {"theta": "pi/2", "phi": "0"}

# The spacetimes with no conformal diagram, for which nothing is written.
NOT_DRAWN = {"godel", "stockum_dust", "taub_nut", "kasner", "bianchi", "tolman_bondi", "alcubierre",
             "natario", "krasnikov", "pp_wave", "mixmaster", "lentz"}


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


# The page's drawing, as _layouts/mfs.html draws it: CD_W wide in a margin CD_M, a label's
# size in units of that 628, the most the page sets it at, since the page sets every label at
# the caption's size and never draws the drawing so narrow that that is more, and where its
# anchor pins it, as a fraction of its box.
CD_W, CD_M = 560, 34
CD_LABEL_SIZE = {"lab": 21, "small": 21, "region": 21, "coord": 21}
CD_ANCHOR_SHIFT = {"c": (-0.5, -0.5), "l": (0, -0.5), "r": (-1, -0.5), "t": (-0.5, 0), "b": (-0.5, -1),
                   "tl": (0, 0), "tr": (-1, 0), "bl": (0, -1), "br": (-1, -1)}


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
        self.layers, self.labels, self.legend_items, self.slices = [], [], [], []

    def slice(self, moment, lines=(), points=(), xt=(), label=None):
        """A moment of the embedding diagram, drawn over the view's own lines: `lines` as
        pairs of arrays (p, q), `xt` as polylines already in (X, T), and `points` as (p, q)."""
        drawn = [run for p, q in lines for run in runs(p, q)] + [rounded(line) for line in xt]
        self.slices.append({**moment.json(), "label": label or moment.label, "lines": drawn,
                            "points": [rounded(point(*pq)) for pq in points], "fills": []})

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
        self.clear_labels()
        if self.slices:
            X0, X1, T0, T1 = self.d["box"]
            for mark in self.slices:
                for at in [p for line in mark["lines"] for p in line] + mark["points"]:
                    if not (X0 <= at[0] <= X1 and T0 <= at[1] <= T1):
                        raise AssertionError(f"view {self.d['id']}: the slice {mark['label']} leaves the box at {at}")
            self.place_slice_labels()
            self.d["slices"] = self.slices
            # A slice's label that finds no side clear of every other label is named in the
            # legend instead, as a slice that is a region alone is.
            placed = [self.label_box(L["at"], L["text"], L["anchor"], L["dx"], L["dy"]) for L in self.labels]
            for mark in self.slices:
                if not mark.get("place"):
                    continue
                P = mark["place"]
                box = self.label_box(P["at"], mark["label"], P["anchor"], P["dx"], P["dy"])
                if any(overlap(box, other) for other in placed):
                    del mark["place"]
                else:
                    placed.append(box)
        return self.d

    def label_box(self, at, text, anchor, dx, dy):
        """The box a label takes on the page's drawing, x0, y0, x1, y1 from the box's top left
        corner in units of the 628, y down, at the most the page sets it at, as slices.py's
        label_size() gives it, which is never smaller than MathJax sets it."""
        X0, X1, T0, T1 = self.d["box"]
        s = CD_W / (X1 - X0)
        w, h = (v * CD_LABEL_SIZE["lab"] for v in slices.label_size(text))
        ax, ay = CD_ANCHOR_SHIFT[anchor]
        x, y = (at[0] - X0) * s + dx, (T1 - at[1]) * s + dy
        return x + ax * w, y + ay * h, x + (ax + 1) * w, y + (ay + 1) * h

    def clear_labels(self):
        """No label overlaps another at the size the page sets them at. A label that would
        overlap one before it stands on the other side of its point, above for below or left
        for right, its offset turned with it; one that overlaps from every side, or a label
        centred on its point, which has no other side, stops the script, naming both."""
        placed = []
        for L in self.labels:
            anchor = L["anchor"]
            flips = [(anchor, 1, 1)]
            if anchor != "c":
                vertical = {"t": "b", "b": "t"}
                horizontal = {"l": "r", "r": "l"}
                v_flip = "".join(vertical.get(c, c) for c in anchor)
                h_flip = "".join(horizontal.get(c, c) for c in anchor)
                both = "".join(horizontal.get(c, vertical.get(c, c)) for c in anchor)
                flips += [(a, sx, sy) for a, sx, sy in ((v_flip, 1, -1), (h_flip, -1, 1), (both, -1, -1)) if a != anchor]
            for a, sx, sy in flips:
                box = self.label_box(L["at"], L["text"], a, sx * L["dx"], sy * L["dy"])
                clash = [text for other, text in placed if overlap(box, other)]
                if not clash:
                    break
            else:
                box = self.label_box(L["at"], L["text"], anchor, L["dx"], L["dy"])
                clash = [text for other, text in placed if overlap(box, other)]
                raise AssertionError(f"view {self.d['id']}: the label {L['text']} overlaps {clash[0]}")
            L["anchor"], L["dx"], L["dy"] = a, sx * L["dx"], sy * L["dy"]
            placed.append((box, L["text"]))

    def place_slice_labels(self):
        """Each slice's label, placed by slices.place on the drawing as the page draws it,
        628 units wide with its margin, at the size of the small labels, clear of every
        label the view already carries."""
        X0, X1, T0, T1 = self.d["box"]
        s = CD_W / (X1 - X0)

        def px(at):
            return np.array([CD_M + (at[0] - X0) * s, CD_M + (T1 - at[1]) * s])
        others = []
        for L in self.labels:
            size = CD_LABEL_SIZE.get(L["class"], 14)
            w, h = (v * size for v in slices.label_size(L["text"]))
            ax, ay = CD_ANCHOR_SHIFT[L["anchor"]]
            x, y = px(L["at"]) + [L["dx"], L["dy"]]
            others.append((x + ax * w, y + ay * h, x + (ax + 1) * w, y + (ay + 1) * h))
        marks = [slices.Mark(None, mark["lines"], mark["points"], mark["fills"], mark["label"]) for mark in self.slices]
        box = (CD_W + 2 * CD_M, (T1 - T0) * s + 2 * CD_M)
        for mark, where in zip(self.slices, slices.place(marks, px, box, CD_LABEL_SIZE["small"], others)):
            if where:
                mark["place"] = where


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
               degenerate, as the axis of Kerr does at theta = 0;
    quotient   a coordinate the metric does not depend on, divided out: the surface's metric is
               g_ab - g_ak g_bk/g_kk, orthogonal to its orbits, as null_rays.py's rays of no
               angular momentum take it, and the published inverse's block is its inverse;
    off_shock  every Dirac delta, and every derivative of one, read as zero, for a surface
               evaluated only away from the hypersurface the delta sits on.
    """

    def __init__(self, sources, metric_id, system_id, plane, fixed=None, params=None,
                 functions=None, numeric=(), axis=None, quotient=None, off_shock=False):
        metric, entry, reader = nr.load(metric_id, system_id)
        self.off_shock = off_shock
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
        if quotient:
            k = [R._plain(c) for c in coords].index(quotient)
            if any(R.symbol[coords[k]] in self.prep(value).free_symbols for value in g):
                raise SystemExit(f"{metric_id}/{system_id}: the published metric depends on {quotient}, "
                                 "which the surface divides out")
            g = sp.Matrix(len(coords), len(coords), lambda a, b: g[a, b] - g[a, k] * g[b, k] / g[k, k])
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
        if self.off_shock:
            e = e.replace(lambda x: isinstance(x, sp.DiracDelta), lambda x: sp.Integer(0))
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


def ds_static_pq(t, r):
    """de Sitter's static patch in its global square at L = 1: tan p = tanh(u/2) and
    tan q = tanh(v/2), with u, v = t -+ artanh(r)."""
    t, r = np.asarray(t, dtype=float), np.asarray(r, dtype=float)
    rs = np.arctanh(r)
    return np.arctan(np.tanh((t - rs) / 2)), np.arctan(np.tanh((t + rs) / 2))


def ds_hyperboloid(p, q):
    """The point (X0, X4, |X_1..3|) of de Sitter's hyperboloid at L = 1 at p, q of its global
    square: X0 = tan T, X4 = cos(chi)/cos T and |X_1..3| = sin(chi)/cos T, with T = p + q and
    chi = q - p."""
    T, chi = p + q, q - p
    return np.tan(T), np.cos(chi) / np.cos(T), np.sin(chi) / np.cos(T)


def ds_closed_pq(t, chi):
    """de Sitter's closed slicing at L = 1, X0 = sinh t, X4 = cosh t cos chi and
    |X_1..3| = cosh t sin chi, in its global square by ds_hyperboloid read backwards: tan T = X0,
    and chi the same angle. Its time runs over the whole line and chi from 0 to pi."""
    T = np.arctan(np.sinh(np.asarray(t, dtype=float)))
    chi = np.asarray(chi, dtype=float)
    return (T - chi) / 2, (T + chi) / 2


def ads_global_pq(t, r):
    """Anti-de Sitter's static global chart in its strip at L = 1: sigma = arctan r and
    p, q = (t -+ sigma)/2."""
    t, r = np.asarray(t, dtype=float), np.asarray(r, dtype=float)
    return (t - np.arctan(r)) / 2, (t + np.arctan(r)) / 2


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

    def settles(self, name, K_near, K_nearer):
        """The Kretschmann scalar at two points approaching a line drawn as regular, the second ten
        times closer: both finite and apart by less than a part in 10^4 of the largest, for a
        spacetime whose curvature there is large in the units it is drawn in."""
        K1, K2 = np.asarray(K_near, dtype=float), np.asarray(K_nearer, dtype=float)
        change = float(np.max(np.abs(K2 - K1)) / np.max(np.abs(K2)))
        ok = bool(np.all(np.isfinite(K1)) and np.all(np.isfinite(K2)) and change < 1e-4)
        self.limits.append({"name": name, "error": change, "ok": ok, "kind": "settles"})

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
            shown = {"diverges": "K", "finite": "K", "settles": "change"}.get(x.get("kind"), "err")
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


def overlap(a, b):
    """Whether two boxes x0, y0, x1, y1 overlap; boxes that touch do not."""
    return a[0] < b[2] and b[0] < a[2] and a[1] < b[3] and b[1] < a[3]


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
    # The embedding's moment t = 0, out to where it reaches, on every view: along r on the
    # triangles, and through the centre along x on the diamonds of the plane y = z = 0.
    plane = slices.moments("minkowski")[0]
    reach = plane.reach("spherical", "r")[1]
    along_r, along_x = np.linspace(0, reach, 2), np.linspace(-reach, reach, 3)
    v.slice(plane, [mink_pq(0 * along_r, along_r)])
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
    v.slice(plane, [null_map(-along_r, along_r)])
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
    v.set(restriction="The plane $y = z = 0$ only, totally geodesic, each point in the diagram a "
                      "single event.")
    v.slice(plane, [mink_pq(0 * along_x, along_x)])
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
    v.set(restriction="The plane $y = z = 0$ only, totally geodesic, each point in the diagram a "
                      "single event.")
    v.slice(plane, [null_map(-along_x, along_x)])
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
    v.set(restriction="The plane $Y = Z = 0$ only, totally geodesic, each point in the diagram a "
                      "single event.")
    # T = 0 in the wedge, where x = X, and its mirror x < 0 beyond the horizon, one line.
    wedge = np.linspace(1e-12, reach, 2)
    left = mink_pq(0 * wedge, -wedge[::-1])
    right = rindler(0 * wedge, wedge)
    v.slice(plane, [(np.concatenate([left[0], right[0]]), np.concatenate([left[1], right[1]]))])
    views.append(v)
    return views


# ---------------------------------------------------------------- Misner space

MISNER_PSI0 = 2.0       # a boost of rapidity 1, so that several copies fit the drawing


def misner(ck, src):
    """Misner space in the plane y = z = 0 of the Minkowski space that covers it, one view per
    chart, each tinting one copy of the spacetime between two lines the boost carries onto one
    another.

    In the covering plane, with l = 1, Levanony and Ori's ct - x = -2 e^(-psi/2) and
    ct + x = 2T e^(psi/2) give -d(ct - x) d(ct + x) = -2 dT dpsi - T dpsi^2, so Misner's
    coordinates cover the half x > ct and p = arctan(ct - x), q = arctan(ct + x) bring it to the
    half square p < 0. The Milne chart's ct cosh chi, ct sinh chi, t < 0, is the past light cone
    of the origin, p, q < 0, and the Rindler chart's xi sinh eta, xi cosh eta the wedge
    p < 0 < q. The boost of rapidity psi_0/2 is psi -> psi + psi_0, chi -> chi + psi_0/2 and
    eta -> eta + psi_0/2. At psi_0 = 4 pi, the value the other diagrams are drawn at, one copy
    stretches by e^(2 pi) along each light ray and fills the drawing to within a pixel, so these
    views are drawn at psi_0 = 2. The horizon T = 0 is q = 0, and psi -> infinity is p = 0."""
    views = []
    half_box = [-HALF - 0.35, PI + 0.35, -PI - 0.25, HALF + 0.25]
    past_box = [-HALF - 0.35, HALF + 0.35, -PI - 0.25, 0.25]
    wedge_box = [-0.35, PI + 0.35, -HALF - 0.25, HALF + 0.25]
    psi0 = MISNER_PSI0
    params = {"psi_0": repr(psi0)}
    settings = f"$\\psi_0 = {psi0:g}$, a boost of rapidity $1$, and $\\ell = 1$."
    restriction = ("The plane $y = z = 0$ of the covering Minkowski space only, each point in the diagram a "
                   "single event, each event of Misner space drawn once in every copy.")
    ks = range(-4, 5)

    def pq(u, v):
        return np.arctan(np.asarray(u, dtype=float)), np.arctan(np.asarray(v, dtype=float))

    def misner_pq(T, psi):
        T, psi = np.asarray(T, dtype=float), np.asarray(psi, dtype=float)
        return pq(-2 * np.exp(-psi / 2), 2 * T * np.exp(psi / 2))

    def milne_pq(t, chi):
        t, chi = np.asarray(t, dtype=float), np.asarray(chi, dtype=float)
        return pq(t * np.exp(-chi), t * np.exp(chi))

    def rindler_pq(eta, xi):
        eta, xi = np.asarray(eta, dtype=float), np.asarray(xi, dtype=float)
        return pq(-xi * np.exp(-eta), xi * np.exp(eta))

    def corners(*pqs):
        return [point(p, q) for p, q in pqs]

    moments = slices.moments("misner")

    def hyperbola(t):
        """The moment ct = t in the covering plane: (ct - x)(ct + x) = t^2 with both negative,
        every copy of its circle, out to where either null coordinate is 50."""
        s = np.geomspace(t * t / 50, 50, 801)
        return pq(-s, -t * t / s)

    past_null = [[point(-HALF, -HALF), point(0, -HALF)], [point(-HALF, -HALF), point(-HALF, 0)]]

    # Misner's own coordinates.
    mis = Plane(src, "misner", "misner", ("T", "\\psi"), {"y": "0", "z": "0"}, params)
    ck.chart("Misner", mis, misner_pq, ck.uniform(-5, 5), ck.uniform(-6, 6), lambda T, psi: (np.abs(T) / 2 + 1, 1))
    p, q = misner_pq(np.zeros(5), np.linspace(-4, 4, 5))
    ck.limit("Misner: T = 0 is the null line q = 0", q, np.zeros(5), 1e-12)
    p, q = misner_pq(np.linspace(-3, 3, 5), np.full(5, 200.0))
    ck.limit("Misner: psi -> infinity is the null line p = 0", p, np.zeros(5), 1e-12)
    ck.finite("Misner: the chronology horizon T = 0 is regular", mis.kretschmann(np.zeros(50), ck.uniform(-5, 5, 50)))
    v = View("misner", "Misner", half_box, "misner")
    v.fill("region", corners((-HALF, -HALF), (-HALF, HALF), (0, HALF), (0, -HALF)))
    lo, hi = math.atan(-2.0), math.atan(-2 * math.exp(-psi0 / 2))
    v.fill("cover", corners((lo, -HALF), (lo, HALF), (hi, HALF), (hi, -HALF)))
    grid(v, "t", misner_pq, (-2, -1, -0.5, 0.5, 1, 2), S_ALL)
    grid(v, "r", misner_pq, [k * psi0 for k in ks], S_ALL, first=False)
    v.line("horizon", [corners((-HALF, 0), (0, 0))])
    v.line("chartedge", [corners((0, -HALF), (0, HALF))])
    v.line("scri", [corners((-HALF, -HALF), (-HALF, HALF), (0, HALF))] + past_null[:1])
    for at, text, anchor, dx, dy in (((PI, 0), "$i^0$", "l", 6, 0), ((0, -PI), "$i^-$", "t", 0, 6)):
        v.layers.append({"kind": "point", "class": "infinity", "at": [round(at[0], 4), round(at[1], 4)]})
        v.label_xt(at, text, anchor, dx=dx, dy=dy)
    v.label_xt([0, -PI / 2 - 0.35], "$T < 0$", "c", "small")
    v.label_xt([PI / 2 + 0.45, 0], "$T > 0$", "c", "small")
    v.legend("cover", "one copy of Misner space, $0 \\le \\psi < \\psi_0$")
    v.legend("r", "$\\psi = k\\psi_0$ for integer $k$, each the same light ray of Misner space")
    v.legend("t", "$T$ constant, at $\\pm1/2$, $\\pm1$ and $\\pm2$ times $\\ell^2$")
    v.legend("horizon", "$T = 0$, the chronology horizon")
    v.legend("chartedge", "$\\psi \\to \\infty$, where the coordinates end")
    v.set(restriction=restriction, settings=settings)
    for m in moments:
        v.slice(m, [hyperbola(m.time)])
    views.append(v)

    # The Milne chart of the region T < 0.
    mil = Plane(src, "misner", "milne", ("t", "\\chi"), {"y": "0", "z": "0"}, params)
    ck.chart("Misner, Milne", mil, milne_pq, ck.uniform(-10, -0.01), ck.uniform(-3, 3), lambda t, chi: (1, 0))
    p, q = milne_pq(np.full(5, -1e-12), np.linspace(-2, 2, 5))
    ck.limit("Misner, Milne: t -> 0 is the null lines p = 0 and q = 0", np.concatenate([p, q]), np.zeros(10), 1e-11)
    v = View("milne", "Milne", past_box, "milne")
    v.fill("region", corners((-HALF, -HALF), (-HALF, 0), (0, 0), (0, -HALF)))
    s = -np.geomspace(1e-6, 1e6, 400)
    a, b = milne_pq(s, np.zeros_like(s))
    c, d = milne_pq(s[::-1], np.full_like(s, psi0 / 2))
    v.fill("cover", [point(x, y) for x, y in zip(np.concatenate([a, c]), np.concatenate([b, d]))])
    grid(v, "t", milne_pq, (-0.5, -1, -2, -4), S_ALL)
    grid(v, "r", milne_pq, [k * psi0 / 2 for k in ks], -np.geomspace(1e-6, 1e6, 400), first=False)
    v.line("horizon", [corners((-HALF, 0), (0, 0)), corners((0, -HALF), (0, 0))])
    v.line("scri", past_null)
    v.layers.append({"kind": "point", "class": "infinity", "at": [0.0, round(-PI, 4)]})
    v.label_xt([0, -PI], "$i^-$", "t", dy=6)
    v.legend("cover", "one copy of the region $T < 0$, $0 \\le \\chi < \\psi_0/2$")
    v.legend("r", "$\\chi = k\\psi_0/2$ for integer $k$, each the same line of Misner space")
    v.legend("t", "$ct$ constant, at $-1/2$, $-1$, $-2$ and $-4$ times $\\ell$")
    v.legend("horizon", "$t \\to 0$, the chronology horizon of each extension")
    v.set(restriction=restriction, settings=settings)
    for m in moments:
        v.slice(m, [hyperbola(m.time)])
    views.append(v)

    # The Rindler chart of the region T > 0.
    rin = Plane(src, "misner", "rindler", ("\\eta", "\\xi"), {"y": "0", "z": "0"}, params)
    ck.chart("Misner, Rindler", rin, rindler_pq, ck.uniform(-3, 3), ck.uniform(0.01, 10), lambda eta, xi: (1, 0))
    p, q = rindler_pq(np.linspace(-2, 2, 5), np.full(5, 1e-12))
    ck.limit("Misner, Rindler: xi -> 0 is the null lines p = 0 and q = 0", np.concatenate([p, q]), np.zeros(10), 1e-11)
    v = View("rindler", "Rindler", wedge_box, "rindler")
    v.fill("region", corners((-HALF, 0), (-HALF, HALF), (0, HALF), (0, 0)))
    s = np.geomspace(1e-6, 1e6, 400)
    a, b = rindler_pq(np.zeros_like(s), s)
    c, d = rindler_pq(np.full_like(s, psi0 / 2), s[::-1])
    v.fill("cover", [point(x, y) for x, y in zip(np.concatenate([a, c]), np.concatenate([b, d]))])
    grid(v, "r", lambda xi, eta: rindler_pq(eta, xi), (0.25, 0.5, 1, 2, 4), S_ALL)
    grid(v, "t", rindler_pq, [k * psi0 / 2 for k in ks], s)
    v.line("horizon", [corners((-HALF, 0), (0, 0)), corners((0, HALF), (0, 0))])
    v.line("scri", [corners((-HALF, 0), (-HALF, HALF), (0, HALF))])
    v.layers.append({"kind": "point", "class": "infinity", "at": [round(PI, 4), 0.0]})
    v.label_xt([PI, 0], "$i^0$", "l", dx=6)
    v.legend("cover", "one copy of the region $T > 0$, $0 \\le \\eta < \\psi_0/2$")
    v.legend("r", "$\\xi$ constant, each a closed timelike curve of Misner space, at $1/4$, $1/2$, $1$, $2$ and $4$ times $\\ell$")
    v.legend("t", "$\\eta = k\\psi_0/2$ for integer $k$, each the same line of Misner space")
    v.legend("horizon", "$\\xi \\to 0$, the chronology horizon of each extension")
    v.set(restriction=restriction, settings=settings)
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
                # Above the moment t = 0, which runs through the middle of the exterior.
                v.label_xt([sx * HALF, base + 0.3], "exterior", cls="region")
        v.label_xt([0, HALF + 0.1], "black hole", cls="region")
        v.label_xt([0, 1.5 * PI - 0.1], "white hole", cls="region")
        v.label_xt([0, -HALF], "white hole", cls="region")
        v.label_xt([0, 2.5 * PI], "black hole", cls="region")
        v.label_xt([Q4, Q4], "$r_+$", "tl", "small", dx=5, dy=1)
        v.label_xt([Q4, 3 * Q4], "$r_-$", "bl", "small", dx=5, dy=-1)
        for sx in (1, -1):
            if self.singular:
                v.label_xt([sx * HALF, PI + 0.35], "$r = 0$", "l" if sx > 0 else "r", dx=8 * sx)
                v.label_xt([sx * 1.05, PI + 0.3], inner, cls="region")
            else:
                v.layers.append({"kind": "point", "class": "infinity", "at": [round(sx * PI, 4), round(PI, 4)]})
                v.label_xt([sx * PI, PI], "$i^0$", "l" if sx > 0 else "r", dx=6 * sx)
                v.label_xt([sx * 3 * Q4, PI + Q4], "$\\mathscr{I}^+$", "bl" if sx > 0 else "br", dx=4 * sx, dy=-3)
                v.label_xt([sx * 3 * Q4, PI - Q4], "$\\mathscr{I}^-$", "tl" if sx > 0 else "tr", dx=4 * sx, dy=3)
                v.label_xt([sx * 2.2, PI], inner, cls="region")


def through_bifurcation(T, cells, far, near):
    """The moment t = 0 of a tower, from r = far in the first cell to the bifurcation point
    at the horizon r = near and on to r = far in the second, as one line of (p, q)."""
    r = np.linspace(far, near, 2)
    a, b = T.pq(cells[0], 0 * r, r), T.pq(cells[1], 0 * r, r[::-1])
    return (np.concatenate([a[0], b[0]]), np.concatenate([a[1], b[1]]))


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

    flamm = slices.moments("schwarzschild")[0]
    lo, hi = flamm.reach("spherical", "r")
    rr_moment = np.linspace(lo, hi, 2)
    moment = [T.pq("I'", 0 * rr_moment, rr_moment[::-1]), T.pq("I", 0 * rr_moment, rr_moment)]
    moment = [(np.concatenate([moment[0][0], moment[1][0]]), np.concatenate([moment[0][1], moment[1][1]]))]
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
    v.slice(flamm, moment)
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
    v.slice(flamm, moment)
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
    v.slice(flamm, moment)
    views.append(v)
    return views


def global_monopole(ck, src):
    """Letelier's black hole in a cloud of strings at Delta = 0.19 and r_s = 1, and the monopole
    with no mass at its centre in the Barriola-Vilenkin chart.

    On the plane of t and r the static metric is -f dt^2 + dr^2/f with f = (1 - Delta)(1 - r_h/r),
    r_h = r_s/(1 - Delta), which is 1/(1 - Delta) times -F dtau^2 + dr^2/F with F = 1 - r_h/r and
    tau = (1 - Delta)t: Schwarzschild's plane with r_h for r_s, a constant factor away. So the Tower
    of F draws it at tau, Kruskal's square, and the tortoise coordinate of the Eddington-Finkelstein
    charts, r_* = r/(1 - Delta) + r_s ln|(1 - Delta)r/r_s - 1|/(1 - Delta)^2, is the Tower's own
    r + r_h ln|r/r_h - 1| over 1 - Delta, so their v and u are the Tower's over 1 - Delta as well:
    V = exp((1 - Delta)v/2r_h), U = (1 - r/r_h) e^(r/r_h)/V ingoing, and the time reverse outgoing.

    The Barriola-Vilenkin plane is -c^2dt^2 + dr^2, Minkowski's, drawn as its triangle by
    p, q = arctan(ct -+ r); its centre r = 0 is a curvature singularity, where the published
    Kretschmann scalar 4 Delta^2/((1 - Delta)^2 r^4) diverges.
    """
    D, params = 0.19, {"Delta": "19/100", "r_s": 1}
    rh = 100 / 81
    sph = Plane(src, "global_monopole", "static", ("t", "r"), EQUATOR, params)
    assert sph.g[0, 1] == 0 and sp.simplify(sph.g[0, 0] * sph.g[1, 1] + 1) == 0
    T = Tower(sp.simplify(-sph.g[0, 0] * sp.Rational(100, 81)), sph.x1, [sp.Rational(100, 81)])

    def cell(name):
        return lambda t, r: T.pq(name, (1 - D) * np.asarray(t, dtype=float), r)
    ck.chart("global monopole static, exterior", sph, cell("I"),
             ck.uniform(-15, 15), ck.uniform(rh + 0.001, 30), lambda t, r: (1, 0))
    ck.chart("global monopole static, black hole", sph, cell("II"),
             ck.uniform(-15, 15), ck.uniform(0.01, rh - 0.001), lambda t, r: (0, -1))
    ck.chart("global monopole static, white hole", sph, cell("IV"),
             ck.uniform(-15, 15), ck.uniform(0.01, rh - 0.001), lambda t, r: (0, 1))
    ck.chart("global monopole static, other exterior", sph, cell("I'"),
             ck.uniform(-15, 15), ck.uniform(rh + 0.001, 30), lambda t, r: (-1, 0))

    def ingoing(w, r):
        w, r = (1 - D) * np.asarray(w, dtype=float), np.asarray(r, dtype=float) / rh
        return np.arctan((1 - r) * np.exp(r - w / (2 * rh))), atan_exp(w / (2 * rh))

    def outgoing(u, r):
        u, r = (1 - D) * np.asarray(u, dtype=float), np.asarray(r, dtype=float) / rh
        return -atan_exp(-u / (2 * rh)), np.arctan((r - 1) * np.exp(r + u / (2 * rh)))
    ein = Plane(src, "global_monopole", "eddington_finkelstein_ingoing", ("v", "r"), EQUATOR, params)
    ck.chart("global monopole ingoing Eddington-Finkelstein", ein, ingoing,
             ck.uniform(-15, 15), ck.uniform(0.01, 30), lambda w, r: (1, -60))
    eout = Plane(src, "global_monopole", "eddington_finkelstein_outgoing", ("u", "r"), EQUATOR, params)
    ck.chart("global monopole outgoing Eddington-Finkelstein", eout, outgoing,
             ck.uniform(-15, 15), ck.uniform(0.01, 30), lambda u, r: (1, 60))

    p, q = cell("II")(np.array([-5.0, 0, 5]), np.full(3, 1e-9))
    ck.limit("global monopole: r -> 0 in the black hole lands on T = pi/2", p + q, [HALF] * 3)
    p, q = cell("I")(np.array([0.0]), np.array([1e8]))
    ck.limit("global monopole: r -> infinity at t = 0 lands on i0, (X, T) = (pi, 0)", point(p[0], q[0]), [PI, 0], 1e-3)
    p, q = cell("I")(np.array([3.0]), np.array([rh * (1 + 1e-12)]))
    ck.limit("global monopole: r -> r_h at fixed t lands on the bifurcation sphere", point(p[0], q[0]), [0, 0], 1e-4)
    rstar = lambda r: T.rstar(r) / (1 - D)
    ck.limit("global monopole: r_* is the published r/(1 - Delta) + r_s ln|(1 - Delta)r/r_s - 1|/(1 - Delta)^2",
             rstar(np.array([0.5, 3.0])),
             np.array([0.5, 3.0]) / (1 - D) + np.log(np.abs((1 - D) * np.array([0.5, 3.0]) - 1)) / (1 - D) ** 2, 1e-12)
    ck.limit("global monopole: the ingoing and static coordinates put one event at one point",
             ingoing(2.0 + rstar(3.0), 3.0), cell("I")(2.0, 3.0), 1e-12)
    ck.limit("global monopole: the outgoing and static coordinates put one event at one point",
             outgoing(2.0 - rstar(3.0), 3.0), cell("I")(2.0, 3.0), 1e-12)
    K = sph.kretschmann
    ck.diverges("global monopole: the Kretschmann scalar diverges at r = 0", K(0, 1e-2), K(0, 1e-3))
    ck.finite("global monopole: the Kretschmann scalar is finite at r = r_h",
              K(np.zeros(3), np.array([rh - 0.001, rh, rh + 0.001])))

    box = [-PI - 0.25, PI + 0.25, -HALF - 0.25, HALF + 0.25]
    hexagon = [[PI, 0], [HALF, HALF], [-HALF, HALF], [-PI, 0], [-HALF, -HALF], [HALF, -HALF]]
    R_OUT, R_IN, TS = (1.3, 1.5, 2, 3, 4), (0.5, 0.9, 1.1), (-4, -2, -1, 0, 1, 2, 4)

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
        v.label_xt([-Q4, Q4], "$r = r_h$", "tr", "small", dx=-6, dy=2)
        v.label_xt([HALF, -0.95], "exterior", cls="region")
        v.label_xt([-HALF, 0], "exterior", cls="region")
        v.label_xt([0, 1.15], "black hole", cls="region")
        v.label_xt([0, -1.15], "white hole", cls="region")
        v.legend("horizon", "the horizon $r = r_h = r_s/(1 - \\Delta)$")
        v.legend("singular", "$r = 0$, where the Kretschmann scalar diverges")
        v.legend("scri", "null infinity $\\mathscr{I}^\\pm$")

    hole = slices.moments("global_monopole", "black_hole")[0]
    lo, hi = hole.reach("static", "r")
    rr_moment = np.linspace(lo, hi, 2)
    moment = [cell("I'")(0 * rr_moment, rr_moment[::-1]), cell("I")(0 * rr_moment, rr_moment)]
    moment = [(np.concatenate([moment[0][0], moment[1][0]]), np.concatenate([moment[0][1], moment[1][1]]))]
    views = []
    t = spread(-np.inf, np.inf, 500, 9)
    v = View("static", "Static Spherical", box, "static")
    v.fill("region", hexagon)
    v.fill("cover", [[0, 0], [HALF, -HALF], [PI, 0], [HALF, HALF]])
    for r in R_OUT:
        v.curve("r", *cell("I")(t, np.full_like(t, r)))
    rr = spread(rh, np.inf, 500, 14)
    for tt in TS:
        v.curve("t", *cell("I")(np.full_like(rr, tt), rr))
    edges(v)
    for r, text in ((1.5, "$1.5\\,r_s$"), (3, "$3\\,r_s$")):
        label_on(v, cell("I")(0.0, r), text)
    v.legend("cover", "the region that $t$ and $r > r_h$ cover")
    v.legend("r", "$r$ constant")
    v.legend("t", "$ct$ constant, in units of $r_s$")
    v.slice(hole, moment)
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
    v.slice(hole, moment)
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
    v.slice(hole, moment)
    views.append(v)

    cone = Plane(src, "global_monopole", "conical", ("t", "r"), EQUATOR, {"Delta": "19/100"})
    ck.chart("global monopole Barriola-Vilenkin", cone, mink_pq, ck.uniform(-20, 20), ck.uniform(0.01, 20),
             lambda t, r: (1, 0))
    K = cone.kretschmann
    ck.diverges("global monopole: the Barriola-Vilenkin Kretschmann scalar diverges at r = 0", K(0, 1e-2), K(0, 1e-3))
    tri_box = [-0.35, PI + 0.35, -PI - 0.25, PI + 0.25]
    v = View("conical", "Barriola-Vilenkin", tri_box, "conical")
    v.fill("region", TRIANGLE)
    v.fill("cover", TRIANGLE)
    grid(v, "r", lambda r, t: mink_pq(t, r), (0.5, 1, 2, 4), S_ALL)
    grid(v, "t", mink_pq, TS, S_POS)
    v.line("singular", [[[0, -PI], [0, PI]]], zig=True)
    v.line("scri", [[[0, PI], [PI, 0]], [[PI, 0], [0, -PI]]])
    for at, text, anchor, dx, dy in (((PI, 0), "$i^0$", "l", 6, 0), ((0, PI), "$i^+$", "b", 0, -6),
                                     ((0, -PI), "$i^-$", "t", 0, 6)):
        v.layers.append({"kind": "point", "class": "infinity", "at": rounded(at)})
        v.label_xt(at, text, anchor, dx=dx, dy=dy)
    v.label_xt([HALF, HALF], "$\\mathscr{I}^+$", "bl", dx=5, dy=-3)
    v.label_xt([HALF, -HALF], "$\\mathscr{I}^-$", "tl", dx=5, dy=3)
    v.label_xt([0, 0.25], "$r = 0$", "r", dx=-6)
    label_on(v, mink_pq(0, 1), "$r = \\ell$")
    label_on(v, mink_pq(0, 4), "$4\\ell$")
    v.legend("cover", "the whole spacetime, which $t$ and $r$ cover")
    v.legend("r", "$r$ constant, in units of $\\ell$")
    v.legend("t", "$ct$ constant")
    v.legend("singular", "$r = 0$, the monopole, where the Kretschmann scalar diverges")
    v.legend("scri", "null infinity $\\mathscr{I}^\\pm$")
    monopole = slices.moments("global_monopole", "monopole")[0]
    along_r = np.linspace(*monopole.reach("conical", "r"), 2)
    v.slice(monopole, [mink_pq(0 * along_r, along_r)])
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
    v = View("tower", "Maximal extension", box, "spherical")
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
    outside, inside = slices.moments("rn_metric", "outside")[0], slices.moments("rn_metric", "inside")[0]
    # Outside r_+ the moment runs in to the outer bifurcation sphere at r_+, and inside r_- out
    # to the inner one at r_-.
    moments = [(outside, through_bifurcation(T, ("I'", "I"), *outside.reach("spherical", "r")[::-1])),
               (inside, through_bifurcation(T, ("III'", "III"), *inside.reach("spherical", "r")))]
    for m, line in moments:
        v.slice(m, [line])
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
    # Low on the observer's world line, clear of the region's own name at its centre.
    v.label(T.pq("I", -0.9 / T.kp, 1.2), "observer at $r = 1.2\\,r_s$", "r", "small", dx=-6)
    v.set(fade={"top": 0.9, "bottom": 0.9})
    v.legend("past", "the whole exterior, in the event's causal past")
    v.legend("world", "a static observer, whose proper time is infinite")
    v.legend("cone", "the edges of the event's past light cone")
    v.legend("mark", "an event beyond the Cauchy horizon $r_-$")
    v.legend("horizon", f"the horizons $r_+ = {rp:g}\\,r_s$ and $r_- = {rm:g}\\,r_s$")
    v.legend("singular", "$r = 0$, where the Kretschmann scalar diverges")
    for m, line in moments:
        v.slice(m, [line])
    views.append(v)
    settings = ("$r_q = 0.48\\,r_s$, so that $r_+ = 0.64\\,r_s$, $r_- = 0.36\\,r_s$, and "
                "$\\kappa_-/\\kappa_+ = 3.2$; at $r_q = 0.4\\,r_s$ the ratio is 16, and every line "
                "inside $r_-$ would lie within $10^{-6}$ of the singularity.")
    for view in views:
        view.set(settings=settings)
    return views


# ---------------------------------------------------------------- Majumdar-Papapetrou

def clip_polygon(pts, T0, T1):
    """The part of a convex polygon, given in (X, T), between the lines T = T0 and T = T1."""
    def cut(poly, keep, T):
        out = []
        for a, b in zip(poly, poly[1:] + poly[:1]):
            ina, inb = keep(a[1]), keep(b[1])
            if ina:
                out.append(a)
            if ina != inb:
                s = (T - a[1]) / (b[1] - a[1])
                out.append([a[0] + s * (b[0] - a[0]), T])
        return out
    return cut(cut([list(p) for p in pts], lambda T: T >= T0, T0), lambda T: T <= T1, T1)


def majumdar_papapetrou(ck, src):
    """One hole alone, U = 1 + m/r at m = 1: the extremal Reissner-Nordstrom black hole in
    isotropic coordinates, R = r + m the areal radius. f = (1 - m/R)^2 has a double root at R = m,
    so the surface gravity vanishes and Kruskal's exponential, which needs one, is not available. The tortoise
    coordinate outside, in the published isotropic chart, is r_* = r + 2 ln r - 1/r, and inside,
    in Reissner-Nordstrom's published chart at r_s = 2m and r_q = m, R_* = R + 2 ln|R - 1| -
    1/(R - 1) - 1, which vanishes at R = 0; each is checked against the published g^rr. With
    u, v = t -+ r_*, every region is a copy of one of two cells, shifted by (k pi, k pi):

        E_k   exterior        p = k pi + arctan u,        q = k pi + arctan v
        B_k   inside R = m    p = (k + 1) pi + arctan u,  q = k pi + arctan v

    Each null coordinate runs on across a horizon, v from E_k into B_k and u from B_k into
    E_(k+1), so the map is continuous there; the horizon is at u -> +-infinity or v -> +infinity,
    where the arctangent reaches pi/2. R = 0 is u = v in B_k, the vertical line X = -pi. The
    moment t = 0 of the embedding diagram, p = -q = -arctan r_*, runs across E_0 at T = 0, from
    X = -pi at the horizon's end of the throat to i^0."""
    ext = Plane(src, "majumdar_papapetrou", "isotropic", ("t", "r"), EQUATOR, {"m": 1})
    inner = Plane(src, "rn_metric", "spherical", ("t", "r"), EQUATOR, {"r_s": 2, "r_q": 1})
    r, R = ext.x1, inner.x1
    for pl, rstar, name in ((ext, r + 2 * sp.log(r) - 1 / r, "r_*"),
                            (inner, R + 2 * sp.log(1 - R) - 1 / (R - 1) - 1, "R_*")):
        assert pl.g[0, 1] == 0 and sp.simplify(pl.g[0, 0] * pl.g[1, 1] + 1) == 0
        # g_tt g_rr = -1, so a radial light ray has dr_*/dr = g_rr.
        ck.limit(f"Majumdar-Papapetrou: d{name}/dr is the published g_rr",
                 0.0 if sp.simplify(sp.diff(rstar, pl.x1) - pl.g[1, 1]) == 0 else 1.0, 0.0, 1e-12)
    ck.limit("Majumdar-Papapetrou: the horizon is the double root r = 0 of the published isotropic g^rr",
             [float(x) for x in sp.roots(sp.numer(sp.together(ext.gi[1, 1])), r, multiple=True)], [0.0, 0.0], 1e-12)
    ck.limit("Majumdar-Papapetrou: inside, the horizon is the double root R = m of Reissner-Nordstrom's g^rr at r_q = r_s/2",
             [float(x) for x in sp.roots(sp.numer(sp.together(inner.gi[1, 1])), R, multiple=True)], [1.0, 1.0], 1e-12)

    def rs_out(x):
        x = np.asarray(x, dtype=float)
        return x + 2 * np.log(x) - 1 / x

    def rs_in(x):
        x = np.asarray(x, dtype=float)
        with np.errstate(divide="ignore"):
            return x + 2 * np.log(np.abs(1 - x)) - 1 / (x - 1) - 1

    def E(k):
        return lambda t, x: (k * PI + np.arctan(np.asarray(t, dtype=float) - rs_out(x)),
                             k * PI + np.arctan(np.asarray(t, dtype=float) + rs_out(x)))

    def B(k):
        return lambda t, x: ((k + 1) * PI + np.arctan(np.asarray(t, dtype=float) - rs_in(x)),
                             k * PI + np.arctan(np.asarray(t, dtype=float) + rs_in(x)))

    span = 20
    ck.chart("Majumdar-Papapetrou, exterior", ext, E(0), ck.uniform(-span, span), ck.uniform(1e-3, 40), lambda t, x: (1, 0))
    ck.chart("Majumdar-Papapetrou, inside the horizon", inner, B(0), ck.uniform(-span, span), ck.uniform(1e-3, 1 - 1e-3),
             lambda t, x: (1, 0))
    ts = np.array([-6.0, 0.0, 6.0])
    p, q = E(0)(ts, np.full(3, 1e-9))
    ck.limit("Majumdar-Papapetrou: r -> 0 lands on the horizon p = pi/2", p, [HALF] * 3, 1e-6)
    p, q = B(0)(ts, np.full(3, 1 - 1e-9))
    ck.limit("Majumdar-Papapetrou: R -> m from inside lands on the same horizon p = pi/2", p, [HALF] * 3, 1e-6)
    p, q = B(0)(ts, np.full(3, 1e-12))
    ck.limit("Majumdar-Papapetrou: R -> 0 lands on the vertical line X = -pi", q - p, [-PI] * 3, 1e-9)
    ck.diverges("Majumdar-Papapetrou: the Kretschmann scalar diverges at R = 0",
                inner.kretschmann(0, 1e-2), inner.kretschmann(0, 1e-3))
    ck.finite("Majumdar-Papapetrou: the Kretschmann scalar is finite at the horizon, from outside and inside",
              np.concatenate([ext.kretschmann(np.zeros(2), np.array([1e-6, 1e-3])),
                              inner.kretschmann(np.zeros(2), np.array([1 - 1e-6, 1 - 1e-3]))]))

    box = [-PI - 0.45, PI + 0.45, -PI - 0.1, 3 * PI + 0.1]
    T0, T1 = box[2], box[3]
    v = View("one_hole", "One hole", box, "isotropic")

    def cell_polygon(kind, k):
        base = 2 * k * PI
        if kind == "E":
            pts = [[0, base - PI], [-PI, base], [0, base + PI], [PI, base]]
        else:
            pts = [[-PI, base], [0, base + PI], [-PI, base + 2 * PI]]
        return clip_polygon(pts, T0, T1)

    def clipped(p, q):
        X, T = xt(p, q)
        bad = (T < T0) | (T > T1)
        return np.where(bad, np.nan, p), np.where(bad, np.nan, q)

    cells = [("B", -1), ("E", 0), ("B", 0), ("E", 1), ("B", 1)]
    for kind, k in cells:
        v.fill("region", cell_polygon(kind, k))
    v.fill("cover", cell_polygon("E", 0))
    times = (-3.0, -1.0, 0.0, 1.0, 3.0)
    r_out, r_in = (0.2, 0.5, 1.0, 2.0), (0.25, 0.5, 0.75)
    t = spread(-np.inf, np.inf, 500, 10)
    for kind, k in cells:
        fmap, radii, lo, hi = (E(k), r_out, 0.0, np.inf) if kind == "E" else (B(k), r_in, 0.0, 1.0)
        for x in radii:
            v.curve("r", *clipped(*fmap(t, np.full_like(t, x))))
        xs = spread(lo, hi, 600, 16)
        for tt in times:
            v.curve("t", *clipped(*fmap(np.full_like(xs, tt), xs)))
    # The horizons, null infinity and the singularity, each an edge of a cell.
    for k in (0, 1):
        base = 2 * k * PI
        v.line("horizon", [[[-PI, base], [0, base + PI]], [[0, base - PI], [-PI, base]]])
        v.line("scri", [[[0, base + PI], [PI, base]], [[PI, base], [0, base - PI]]])
    v.line("singular", [[[-PI, T0], [-PI, T1]]], zig=True)
    for k in (0, 1):
        base = 2 * k * PI
        v.layers.append({"kind": "point", "class": "infinity", "at": [round(PI, 4), round(base, 4)]})
        v.label_xt([PI, base], "$i^0$", "l", dx=6)
        v.label_xt([3 * Q4, base + Q4], "$\\mathscr{I}^+$", "bl", dx=4, dy=-3)
        v.label_xt([3 * Q4, base - Q4], "$\\mathscr{I}^-$", "tl", dx=4, dy=3)
        v.label_xt([HALF, base + 0.45], "exterior", cls="region")
        v.label_xt([-2.1, base + PI], "$R < m$", cls="region")
    for at, text, anchor in (((0, PI), "$i^\\pm$", "l"), ((0, -PI), "$i^-$", "l"), ((0, 3 * PI), "$i^+$", "l")):
        v.layers.append({"kind": "point", "class": "infinity", "at": [round(at[0], 4), round(at[1], 4)]})
        v.label_xt(list(at), text, anchor, dx=6)
    v.label_xt([-PI, 2 * PI], "$R = 0$", "r", dx=-6)
    v.set(fade={"top": 0.9, "bottom": 0.9})
    v.legend("cover", "the exterior $r > 0$, which $t$ and $r$ cover")
    v.legend("r", f"$r$ constant outside, at {listed(r_out)}, and $R$ constant inside, at {listed(r_in)}, "
                  "in units of $m$")
    v.legend("t", "$t$ constant")
    v.legend("horizon", "the horizon $r = 0$, where $R = m$")
    v.legend("singular", "$R = 0$, a timelike singularity, where the Kretschmann scalar diverges")
    v.legend("scri", "null infinity $\\mathscr{I}^\\pm$")
    moment = slices.moments("majumdar_papapetrou", "one_hole")[0]
    lo, hi = moment.reach("isotropic", "r")
    rr = np.geomspace(lo, hi, 400)
    v.slice(moment, [E(0)(np.zeros_like(rr), rr)])
    v.set(settings="$m = 1$, the unit of every length; inside the horizon, Reissner-Nordström's own chart at "
                   "$r_s = 2m$ and $r_q = m$, where its $r$ is the areal radius $R$.")
    return [v]


# ---------------------------------------------------------------- Schwarzschild-de Sitter

class Kottler:
    """Kottler's spacetime, f = 1 - r_s/r - Lambda r^2/3 with r_s = 1, whose horizons are the
    positive roots r_h < r_c of the published g^rr and whose third root r_n is negative.

    1/f has no polynomial part, so r* = sum_i A_i ln|1 - r/r_i| with A_i = 1/f'(r_i), which
    vanishes at r = 0 and tends to R = -sum_i A_i ln|r_i| as r -> infinity, since sum_i A_i = 0;
    A_h = 1/2k_h and A_c = -1/2k_c for the surface gravities k_h = f'(r_h)/2 and
    k_c = -f'(r_c)/2. With u, v = t -+ r* in each region's own static coordinates, every region
    is drawn through one function,

        g(x) = arctan W(x),   W(x) = exp(-a x - b sqrt(x^2 + 1)),   a, b = (k_h +- k_c)/2,

    decreasing from pi/2 to 0, with g ~ e^(-k_h x) as x -> +infinity and pi/2 - g ~ e^(k_c x)
    as x -> -infinity. So p and q are, up to a factor that tends to 1, the Kruskal coordinates
    of the black hole horizon next to it and of the cosmological horizon next to it, and the
    drawing is continuous with a continuous tangent through both kinds of horizon, where one
    exponential alone would be smooth at one kind and kinked at the other. The regions, as
    cells of side pi/2 in (p, q):

        S    static, r_h < r < r_c        p = -g(u),       q = g(-v)
        S'   the static region beyond r_h  p = g(u),        q = -g(-v)
        B    black hole, r < r_h          p = g(u),        q = g(-v)
        W    white hole, r < r_h          p = -g(u),       q = -g(-v)
        C+   expanding, r > r_c           p = -g(u),       q = pi - g(-v)
        C-   contracting, r > r_c         p = g(u) - pi,   q = g(-v)
        S''  the static region beyond r_c  p = g(u) - pi,   q = pi - g(-v)

    with t the static coordinate of each. The singularity r = 0 is (g(t), g(-t)) in B, and
    future infinity r = infinity is (-g(t - R), pi - g(-t - R)) in C+, each a spacelike curve;
    neither is straight, since g(x) + g(-x) = pi/2 only where k_h = k_c. The chain repeats
    with period pi in p and q, S'' being S' moved by (-pi, pi).
    """

    def __init__(self, plane):
        r = plane.x1
        numerator = sp.numer(sp.together(plane.gi[1, 1]))
        roots = sorted(float(z) for z in sp.Poly(numerator, r).real_roots())
        self.rn, self.rh, self.rc = roots
        self.f = sp.lambdify(r, plane.gi[1, 1], "numpy")
        fp = sp.lambdify(r, sp.diff(plane.gi[1, 1], r), "numpy")
        self.A = [1 / float(fp(ri)) for ri in roots]
        self.kh, self.kc = float(fp(self.rh)) / 2, -float(fp(self.rc)) / 2
        self.a, self.b = (self.kh + self.kc) / 2, (self.kh - self.kc) / 2
        self.R = -sum(A * math.log(abs(ri)) for A, ri in zip(self.A, roots))
        self.roots = roots

    def rstar(self, r):
        r = np.asarray(r, dtype=float)
        with np.errstate(divide="ignore"):
            return sum(A * np.log(np.abs(1 - r / ri)) for A, ri in zip(self.A, self.roots))

    def g(self, x):
        x = np.asarray(x, dtype=float)
        return atan_exp(-self.a * x - self.b * np.sqrt(x * x + 1))

    def pq(self, cell, t, r):
        t = np.asarray(t, dtype=float)
        rs = self.rstar(r)
        return self.null(cell, t - rs, t + rs)

    def null(self, cell, u, v):
        g = self.g
        return {"S": lambda: (-g(u), g(-v)), "S'": lambda: (g(u), -g(-v)), "B": lambda: (g(u), g(-v)),
                "W": lambda: (-g(u), -g(-v)), "C+": lambda: (-g(u), PI - g(-v)),
                "C-": lambda: (g(u) - PI, g(-v)), "S''": lambda: (g(u) - PI, PI - g(-v))}[cell]()

    def singularity(self, sign=1, n=700):
        t = spread(-np.inf, np.inf, n, 40)
        p, q = self.g(t), self.g(-t)
        return (p, q) if sign > 0 else (-q, -p)

    def infinity(self, sign=1, n=700):
        """Future infinity r -> infinity in C+, or with sign -1 past infinity in C-, its image
        under the time reflection T -> -T, (p, q) -> (-q, -p)."""
        t = spread(-np.inf, np.inf, n, 40)
        p, q = -self.g(t - self.R), PI - self.g(-t - self.R)
        return (p, q) if sign > 0 else (-q, -p)


def schwarzschild_de_sitter(ck, src):
    """The maximal extension of Kottler's spacetime, which Lake and Roeder described, drawn at
    r_s = 1 and Lambda = 1/5 through the map of Kottler above. The static chart covers S, the
    ingoing Eddington-Finkelstein chart B, S and C- through q = g(-v) and p from u = v - 2r*
    in each, and the outgoing one W, S and C+ through p = -g(u) and q from v = u + 2r*."""
    st = Plane(src, "schwarzschild_de_sitter", "static", ("t", "r"), EQUATOR, {"r_s": 1, "Lambda": "1/5"})
    assert st.g[0, 1] == 0 and sp.simplify(st.g[0, 0] * st.g[1, 1] + 1) == 0
    K = Kottler(st)
    rh, rc = K.rh, K.rc
    ck.limit("Schwarzschild-de Sitter: the horizons are the roots of the published g^rr",
             [rh, rc], [1.0851996154371, 3.2146274073952], 1e-12)
    ck.limit("Schwarzschild-de Sitter: dr*/dr = 1/f", 
             (K.rstar(np.array([0.5, 2.0, 5.0]) + 1e-6) - K.rstar(np.array([0.5, 2.0, 5.0]) - 1e-6)) / 2e-6
             * K.f(np.array([0.5, 2.0, 5.0])), [1, 1, 1], 1e-6)
    ck.limit("Schwarzschild-de Sitter: the residues are 1/2k_h and -1/2k_c", [K.A[1] * 2 * K.kh, K.A[2] * 2 * K.kc],
             [1, -1], 1e-12)
    span = 15
    for cell, lo, hi, future in (("S", rh + 1e-3, rc - 1e-3, (1, 0)), ("S'", rh + 1e-3, rc - 1e-3, (-1, 0)),
                                 ("S''", rh + 1e-3, rc - 1e-3, (-1, 0)), ("B", 0.01, rh - 1e-3, (0, -1)),
                                 ("W", 0.01, rh - 1e-3, (0, 1)), ("C+", rc + 1e-3, 60, (0, 1)),
                                 ("C-", rc + 1e-3, 60, (0, -1))):
        ck.chart(f"Schwarzschild-de Sitter static, {cell}", st, lambda t, r, c=cell: K.pq(c, t, r),
                 ck.uniform(-span, span), ck.uniform(lo, hi), lambda t, r, f=future: f)

    def ingoing(w, r):
        w, r = np.asarray(w, dtype=float), np.asarray(r, dtype=float)
        u = w - 2 * K.rstar(r)
        q = K.g(-w)
        p = np.where(r < rh, K.g(u), np.where(r < rc, -K.g(u), K.g(u) - PI))
        return p, q

    def outgoing(u, r):
        u, r = np.asarray(u, dtype=float), np.asarray(r, dtype=float)
        v = u + 2 * K.rstar(r)
        p = -K.g(u)
        q = np.where(r < rh, -K.g(-v), np.where(r < rc, K.g(-v), PI - K.g(-v)))
        return p, q
    ein = Plane(src, "schwarzschild_de_sitter", "eddington_finkelstein_ingoing", ("v", "r"), EQUATOR,
                {"r_s": 1, "Lambda": "1/5"})
    eout = Plane(src, "schwarzschild_de_sitter", "eddington_finkelstein_outgoing", ("u", "r"), EQUATOR,
                 {"r_s": 1, "Lambda": "1/5"})
    for name, plane, fmap, sign in (("ingoing", ein, ingoing, 1), ("outgoing", eout, outgoing, -1)):
        for lo, hi in ((0.01, rh - 1e-3), (rh + 1e-3, rc - 1e-3), (rc + 1e-3, 60)):
            # d_v - (1 + |f|) d_r for the ingoing chart and d_u + (1 + |f|) d_r for the outgoing
            # one are timelike wherever f < 2, which is everywhere, and raise T.
            ck.chart(f"Schwarzschild-de Sitter {name} Eddington-Finkelstein, {lo:.2f} < r < {hi:.2f}",
                     plane, fmap, ck.uniform(-span, span), ck.uniform(lo, hi),
                     lambda w, r, s=sign: (1, -s * (1 + np.abs(K.f(r)))))

    ck.limit("Schwarzschild-de Sitter: the ingoing and static coordinates put one event at one point",
             np.concatenate([ingoing(2.0 + K.rstar(r), r) for r in (0.5, 2.0, 5.0)]),
             np.concatenate([K.pq(c, 2.0, r) for c, r in (("B", 0.5), ("S", 2.0), ("C-", 5.0))]), 1e-12)
    ck.limit("Schwarzschild-de Sitter: the outgoing and static coordinates put one event at one point",
             np.concatenate([outgoing(2.0 - K.rstar(r), r) for r in (0.5, 2.0, 5.0)]),
             np.concatenate([K.pq(c, 2.0, r) for c, r in (("W", 0.5), ("S", 2.0), ("C+", 5.0))]), 1e-12)
    for t in (-3.0, 0.0, 3.0):
        ck.limit(f"Schwarzschild-de Sitter: r -> r_h at t = {t:g} lands on the black hole's bifurcation sphere",
                 point(*K.pq("S", t, rh * (1 + 1e-13))), [0, 0], 1e-4)
        ck.limit(f"Schwarzschild-de Sitter: r -> r_c at t = {t:g} lands on the cosmological bifurcation sphere",
                 point(*K.pq("S", t, rc * (1 - 1e-13))), [PI, 0], 1e-4)
    t = np.array([-5.0, 0, 5])
    ck.limit("Schwarzschild-de Sitter: r -> 0 in the black hole lands on the singularity (g(t), g(-t))",
             np.concatenate(K.pq("B", t, np.full(3, 1e-12))), np.concatenate([K.g(t), K.g(-t)]), 1e-9)
    ck.limit("Schwarzschild-de Sitter: r -> infinity in the expanding region lands on (-g(t - R), pi - g(-t - R))",
             np.concatenate(K.pq("C+", t, np.full(3, 1e9))), np.concatenate([-K.g(t - K.R), PI - K.g(-t - K.R)]), 1e-6)
    ck.diverges("Schwarzschild-de Sitter: the Kretschmann scalar diverges at r = 0",
                st.kretschmann(0, 1e-2), st.kretschmann(0, 1e-3))
    ck.finite("Schwarzschild-de Sitter: the Kretschmann scalar is finite at both horizons",
              st.kretschmann(np.zeros(2), np.array([rh, rc])))

    H = HALF
    sing, past_sing = K.singularity(1), K.singularity(-1)
    scri, past_scri = K.infinity(1), K.infinity(-1)

    def region(cell):
        if cell == "B":
            return [(0, 0), (0, H)] + list(zip(*sing))[::-1] + [(H, 0)]
        if cell == "W":
            return [(0, 0), (0, -H)] + list(zip(*past_sing)) + [(-H, 0)]
        if cell == "C+":
            return [(-H, H), (-H, PI)] + list(zip(*scri))[::-1] + [(0, H)]
        if cell == "C-":
            return [(-H, H), (-PI, H)] + list(zip(*past_scri)) + [(-H, 0)]
        corners = {"S": (-H, 0), "S'": (0, -H), "S''": (-PI, H)}[cell]
        p0, q0 = corners
        return [(p0, q0), (p0 + H, q0), (p0 + H, q0 + H), (p0, q0 + H)]
    CELLS = ["S'", "B", "W", "S", "C+", "C-", "S''"]
    t_all = spread(-np.inf, np.inf, 500, 9)

    def edges(v):
        v.line("horizon", [[point(0, -H), point(0, H)], [point(-H, 0), point(H, 0)]])
        v.line("horizon", [[point(-H, 0), point(-H, PI)], [point(-PI, H), point(0, H)]])
        v.line("horizon", [[point(H, -H), point(H, 0)], [point(0, -H), point(H, -H)]])
        v.line("horizon", [[point(-PI, H), point(-PI, PI)], [point(-PI, PI), point(-H, PI)]])
        v.curve("singular", *sing, zig=True, tol=0.004)
        v.curve("singular", *past_sing, zig=True, tol=0.004)
        v.curve("scri", *scri, tol=0.004)
        v.curve("scri", *past_scri, tol=0.004)
        for cx in (-H, H, 3 * H):
            for ct in (H, -H):
                v.layers.append({"kind": "point", "class": "infinity", "at": [round(cx, 4), round(ct, 4)]})
        for cx in (-H, H, 3 * H):
            v.label_xt([cx, H], "$i^+$", "b", dy=-6)
            v.label_xt([cx, -H], "$i^-$", "t", dy=6)
        v.label_xt([0, 1.2], "black hole", cls="region")
        v.label_xt([0, -1.2], "white hole", cls="region")
        v.label_xt([PI, 1.2], "expanding", cls="region")
        v.label_xt([PI, -1.2], "contracting", cls="region")
        for cx in (-H, H, 3 * H):
            v.label_xt([cx, -0.45], "static", cls="region")
        v.label_xt([PI, 1.62], "$\\mathscr{I}^+$", "b", dy=-3)
        v.label_xt([PI, -1.62], "$\\mathscr{I}^-$", "t", dy=3)
        v.label_xt([0, H + 0.15], "$r = 0$", "b", dy=-4)
        v.label_xt([0, -H - 0.15], "$r = 0$", "t", dy=4)
        # Each horizon's name stands on the outer side of its point, which leaves the middle
        # of the static region to the circle r = 2 r_s.
        v.label_xt([Q4, Q4], "$r_h$", "tr", "small", dx=-5, dy=1)
        v.label_xt([PI - Q4, Q4], "$r_c$", "tl", "small", dx=5, dy=1)
        v.legend("horizon", f"the black hole horizons $r_h = {rh:.4g}\\,r_s$ and the cosmological horizons "
                            f"$r_c = {rc:.4g}\\,r_s$")
        v.legend("singular", "$r = 0$, where the Kretschmann scalar diverges")
        v.legend("scri", "future and past infinity $\\mathscr{I}^\\pm$, spacelike")

    moment = slices.moments("schwarzschild_de_sitter")[0]
    lo, hi = moment.reach("static", "r")
    # The moment ends on the bifurcation spheres, where the map of the static chart is 0/0; the
    # limits above place them at (p, q) = (0, 0) and (-pi/2, pi/2), and the next static
    # region's black hole sphere at (-pi, pi).
    rr = np.linspace(lo, hi, 9)[1:-1]
    a, b = K.pq("S", 0 * rr, rr), K.pq("S''", 0 * rr, rr[::-1])
    moment_line = [(np.concatenate([[0], a[0], [-H], b[0], [-PI]]), np.concatenate([[0], a[1], [H], b[1], [PI]]))]
    box = [-PI - 0.3, 2 * PI + 0.3, -HALF - 0.45, HALF + 0.45]
    R_S, R_IN, R_OUT = (1.2, 1.5, 2.0, 2.5, 3.0), (0.4, 0.7, 0.95), (4.0, 6.0, 12.0)
    TS = (-8, -4, -2, 0, 2, 4, 8)
    views = []

    def frame(v, cover):
        for cell in CELLS:
            v.fill("region", [point(*pq) for pq in region(cell)])
        for cell in cover:
            v.fill("cover", [point(*pq) for pq in region(cell)])

    v = View("static", "Static", box, "static")
    frame(v, ["S"])
    for r in R_S:
        v.curve("r", *K.pq("S", t_all, np.full_like(t_all, r)))
    rr = spread(rh, rc, 600, 16)
    for tt in TS:
        v.curve("t", *K.pq("S", np.full_like(rr, tt), rr))
    edges(v)
    label_on(v, K.pq("S", 0.0, 2.0), "$2\\,r_s$", "br", dx=-3)
    v.legend("cover", "the static region $r_h < r < r_c$, which $t$ and $r$ cover")
    v.legend("r", "$r$ constant, at " + listed(R_S) + " in units of $r_s$")
    v.legend("t", "$ct$ constant, in units of $r_s$")
    v.slice(moment, moment_line)
    views.append(v)

    for vid, label, system, fmap, cells, null_text in (
            ("ingoing", "Ingoing Eddington-Finkelstein", "eddington_finkelstein_ingoing", ingoing, ["B", "S", "C-"],
             "$v$ constant, an ingoing light ray"),
            ("outgoing", "Outgoing Eddington-Finkelstein", "eddington_finkelstein_outgoing", outgoing, ["W", "S", "C+"],
             "$u$ constant, an outgoing light ray")):
        v = View(vid, label, box, system)
        frame(v, cells)
        for r in R_IN + R_S + R_OUT:
            v.curve("r", *fmap(t_all, np.full_like(t_all, r)))
        for lo_, hi_ in ((0, rh), (rh, rc), (rc, np.inf)):
            rr = spread(lo_, hi_, 600, 16)
            for w in (-8, -4, 0, 4, 8):
                v.curve("null", *fmap(np.full_like(rr, w), rr))
        edges(v)
        v.legend("cover", ("the black hole, the static region and the contracting region, which $v$ and $r > 0$ cover"
                           if vid == "ingoing" else
                           "the white hole, the static region and the expanding region, which $u$ and $r > 0$ cover"))
        v.legend("r", "$r$ constant, at " + listed(R_IN + R_S + R_OUT) + " in units of $r_s$")
        v.legend("null", null_text)
        v.slice(moment, moment_line)
        views.append(v)
    for view in views:
        view.set(settings=f"$\\Lambda = 0.2/r_s^2$, so that $r_h = {rh:.4g}\\,r_s$, $r_c = {rc:.4g}\\,r_s$, and "
                          f"$\\kappa_c/\\kappa_h = {K.kc / K.kh:.3g}$.")
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
    # The moment of the embedding's equator meets the axis through the outer bifurcation point.
    m = slices.moments(metric_id)[0]
    v.slice(m, [through_bifurcation(T, ("I'", "I"), *m.reach("boyer_lindquist", "r")[::-1])])
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

    def flat(t, rho):
        eta = -np.exp(-np.asarray(t, dtype=float))
        return Q4 + np.arctan(eta - rho), Q4 + np.arctan(eta + rho)
    ck.chart("de Sitter static", st, ds_static_pq, ck.uniform(-10, 10), ck.uniform(0.001, 0.999), lambda t, r: (1, 0))
    ck.chart("de Sitter flat slicing", fl, flat, ck.uniform(-6, 6), ck.uniform(0.001, 20), lambda t, x: (1, 0))

    t, r = ck.uniform(-6, 6, 2000), ck.uniform(0, 0.999, 2000)
    X0, X4, R = ds_hyperboloid(*ds_static_pq(t, r))
    scale = 1 + np.abs(X0) + np.abs(X4)
    ck.limit("de Sitter: the static coordinates land where the hyperboloid puts them",
             np.concatenate([(X0 - np.sqrt(1 - r ** 2) * np.sinh(t)) / scale,
                             (X4 - np.sqrt(1 - r ** 2) * np.cosh(t)) / scale, (R - r) / scale]), 0, 1e-10)
    t, rho = ck.uniform(-4, 4, 2000), ck.uniform(0, 8, 2000)
    X0, X4, R = ds_hyperboloid(*flat(t, rho))
    scale = 1 + np.abs(X0) + np.abs(X4) + np.abs(R)
    a = np.exp(t)
    ck.limit("de Sitter: the flat slicing lands where the hyperboloid puts it",
             np.concatenate([(X0 - np.sinh(t) - rho ** 2 * a / 2) / scale,
                             (X4 - np.cosh(t) + rho ** 2 * a / 2) / scale, (R - a * rho) / scale]), 0, 1e-10)
    p, q = ds_static_pq(np.array([0.0, 3, -3]), np.full(3, 1 - 1e-14))
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

    static_moment = slices.moments("de_sitter")[0]
    lo, hi = static_moment.reach("static_spherical", "r")
    r = np.linspace(lo, hi, 2)
    p, q = ds_static_pq(0 * r, r)
    X, T = xt(p, q)
    # The observer's hemisphere out to its horizon, then the antipode's, the static patch's
    # mirror X -> pi - X, back in to its centre.
    moment_xt = [np.column_stack([np.concatenate([X, PI - X[::-1]]), np.concatenate([T, T[::-1]])])]
    views = []
    v = View("static", "Static spherical", box, "static_spherical")
    frame(v)
    v.fill("cover", [[0, -HALF], [HALF, 0], [0, HALF]])
    grid(v, "r", lambda r, t: ds_static_pq(t, r), (0.3, 0.6, 0.85, 0.97), S_ALL)
    grid(v, "t", ds_static_pq, (-3, -1.5, -0.5, 0, 0.5, 1.5, 3), spread(0, 1, 500, 14))
    label_on(v, ds_static_pq(0, 0.6), "$0.6$")
    label_on(v, ds_static_pq(0, 0.3), "$r = 0.3$")
    v.label_xt([PI * 0.8, -0.28], "the antipode's static patch", cls="region")
    v.label_xt([HALF, 1.05], "expanding", cls="region")
    v.label_xt([HALF, -1.05], "contracting", cls="region")
    v.legend("cover", "the static patch, which $t$ and $r$ cover")
    v.legend("r", "$r$ constant, in units of $\\sqrt{3/\\Lambda}$")
    v.legend("t", "$ct$ constant")
    v.slice(static_moment, xt=moment_xt)
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
    v.slice(static_moment, xt=moment_xt, label="static $t = 0$")
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

    ck.chart("anti-de Sitter global", gl, ads_global_pq, ck.uniform(-10, 10), ck.uniform(0.001, 50), lambda t, r: (1, 0))
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
    hyperboloid = slices.moments("anti_de_sitter")[0]
    lo, hi = hyperboloid.reach("static_global", "r")
    r = np.linspace(lo, hi, 2)
    v.slice(hyperboloid, [ads_global_pq(0 * r, r)])
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
                      "anti-de Sitter space of two dimensions, each point in the diagram a "
                      "single event.")
    v.legend("cover", "the wedge that $t$ and $z$ cover")
    v.legend("r", "$z$ constant")
    v.legend("t", "$ct$ constant")
    v.legend("boundary", "the conformal boundary, timelike")
    v.legend("chartedge", "$z \\to \\infty$, the Poincaré horizon, where $t$ and $z$ end")
    # The static moment t = 0 is the Poincare moment t = 0, and on the plane x = y = 0 its
    # static radius is |1 - z^2|/2z, so the embedding's reach r <= 4L runs from
    # z = sqrt(17) - 4 through the centre z = L to z = sqrt(17) + 4.
    z = np.array([np.sqrt(hi * hi + 1) - hi, np.sqrt(hi * hi + 1) + hi])
    v.slice(hyperboloid, [poincare_pq(0 * z, z)])
    views.append(v)
    return views


# ---------------------------------------------------------------- the BTZ black hole

class BTZTower(Tower):
    """A Tower for f = N^2 = r^2/l^2 - M + J^2/(4r^2), which grows as r^2 and so is no Tower's
    1 plus simple poles: 1/f is a sum of simple poles alone, at every root of the numerator of
    f, the negative ones included,

        r* = sum_i A_i ln|r/r_i - 1|,      A_i = 1/f'(r_i),

    so that r*(0) = 0 as before and now also r* -> 0 as r -> infinity, since the terms of each
    pair of roots +-r_i cancel there. Both identities, and dr*/dr = 1/f, are checked in sympy
    on construction. The horizons are the positive roots, largest first, and every cell is
    written as a Tower's, in G(u) = arctan exp(-k+ u). The conformal boundary r -> infinity has
    u = v = t, which puts it on X = q - p = G(-t) + G(t) = pi/2 in cell I, and r = 0 has
    r* = 0 as well, which puts it on the same vertical line in cell III.
    """

    def __init__(self, f, r):
        self.f_sym = sp.simplify(f)
        poles = sp.solve(sp.numer(sp.together(self.f_sym)), r)
        self.roots = sorted((x for x in poles if x.is_positive), reverse=True)
        self.A = [sp.simplify(sp.limit((r - ri) / self.f_sym, r, ri)) for ri in poles]
        rest = sp.simplify(sp.apart(sp.together(1 / self.f_sym), r, full=True).doit()
                           - sum(Ai / (r - ri) for Ai, ri in zip(self.A, poles)))
        assert rest == 0, f"1/f is not a sum of simple poles: the remainder is {rest}"
        rstar = sum(Ai * sp.log(r / ri - 1) for Ai, ri in zip(self.A, poles))
        assert sp.simplify(sp.diff(rstar, r) - 1 / self.f_sym) == 0, "dr*/dr is not 1/f"
        fp = sp.diff(self.f_sym, r)
        self.kappa = [sp.simplify(sp.Abs(fp.subs(r, ri)) / 2) for ri in self.roots]
        self.kp = float(self.kappa[0])
        self.rf = [float(x) for x in self.roots]
        self.poles = [float(x) for x in poles]
        self.Ap = [float(a) for a in self.A]

    def rstar(self, r):
        r = np.asarray(r, dtype=float)
        with np.errstate(divide="ignore"):
            return sum(A * np.log(np.abs(r / ri - 1)) for A, ri in zip(self.Ap, self.poles))


# The cells of the BTZ tower in (p, q). An exterior is the half of a Tower's cell inside the
# boundary X = pi/2, and the region inside r_- the half inside r = 0, on the same line.
BTZ_CELL = {
    "I": [(0, 0), (-HALF, 0), (0, HALF)],
    "I'": [(0, 0), (0, -HALF), (HALF, 0)],
    "II": [(0, 0), (HALF, 0), (HALF, HALF), (0, HALF)],
    "IV": [(0, 0), (0, -HALF), (-HALF, -HALF), (-HALF, 0)],
    "III": [(0, HALF), (HALF, HALF), (HALF, PI)],
    "III'": [(HALF, 0), (HALF, HALF), (PI, HALF)],
}


def clip_in_t(v, T0, T1):
    """Cut every fill and line of a view to T0 <= T <= T1, where a drawing that continues
    up and down fades out: a fill against the two horizontal lines by Sutherland and Hodgman,
    since every fill cut so is convex, and a line into the runs that stay between them."""
    def cross(a, b, T):
        s = (T - a[1]) / (b[1] - a[1])
        return [a[0] + s * (b[0] - a[0]), T]

    def polygon(pts):
        for T, inside in ((T0, lambda P: P[1] >= T0), (T1, lambda P: P[1] <= T1)):
            out = []
            for i, b in enumerate(pts):
                a = pts[i - 1]
                if inside(b):
                    if not inside(a):
                        out.append(cross(a, b, T))
                    out.append(b)
                elif inside(a):
                    out.append(cross(a, b, T))
            pts = out
            if not pts:
                break
        return pts

    def polyline(pts):
        runs_, run = [], []
        for i, b in enumerate(pts):
            ok = T0 <= b[1] <= T1
            if i and (ok != (T0 <= pts[i - 1][1] <= T1)):
                a = pts[i - 1]
                edge = T0 if min(a[1], b[1]) < T0 else T1
                run.append(cross(a, b, edge))
                if not ok:
                    runs_.append(run)
                    run = []
            if ok:
                run.append(list(b))
        if run:
            runs_.append(run)
        return [r for r in runs_ if len(r) > 1]

    layers = []
    for layer in v.layers:
        if layer["kind"] == "fill":
            pts = polygon([list(x) for x in layer["points"]])
            if len(pts) > 2:
                layers.append({**layer, "points": rounded(pts)})
        elif layer["kind"] in ("line", "zig"):
            layers += [{**layer, "points": rounded(r)} for r in polyline(layer["points"])]
        elif T0 <= layer["at"][1] <= T1:
            layers.append(layer)
    v.layers = layers


def btz(ck, src):
    """The hole without rotation, J = 0, and the rotating hole, J = 4l/5, at M = 1 and l = 1.

    Without rotation f = r^2 - 1 and r* = (1/2) ln|(r - 1)/(r + 1)|, so Kruskal's U = -exp(-u),
    V = exp(v) of the exterior, with k = 1, give UV = (1 - r)/(1 + r): -1 on the conformal
    boundary, 0 on the horizon and 1 at r = 0. With p = arctan U and q = arctan V, tan(q - p) =
    (V - U)/(1 + UV) puts the boundary on X = +-pi/2 and tan(p + q) = (U + V)/(1 - UV) puts r = 0
    on T = +-pi/2: the square of Banados, Henneaux, Teitelboim and Zanelli. The ingoing chart's
    V = exp(v), U = (1 - r)/((1 + r)V) is one formula for every r > 0, and the outgoing chart
    is its time reverse.

    The rotating hole is drawn on its plane of t and r with phi divided out, the metric
    -N^2 dt^2 + dr^2/N^2 its rays of no angular momentum follow, with horizons at r_+^2 = 4/5
    and r_-^2 = 1/5 and k_-/k_+ = r_+/r_- = 2. Its cells are a Tower's, cut by the boundary and
    by r = 0 on the vertical lines X = +-pi/2, and they repeat up the strip between them.
    r = 0 is where the circles of phi shrink to nothing, and the curvature is finite there.
    """
    J0, JR = {"ell": 1, "M": 1, "J": 0}, {"ell": 1, "M": 1, "J": "4/5"}
    st = Plane(src, "btz", "stationary", ("t", "r"), {"phi": "0"}, J0)
    assert st.g[0, 1] == 0 and sp.simplify(st.g[0, 0] * st.g[1, 1] + 1) == 0
    T = BTZTower(st.gi[1, 1], st.x1)
    ck.limit("BTZ, J = 0: the horizon is the root of the published g^rr, at r_+ = l", T.rf, [1.0], 1e-12)
    ck.chart("BTZ stationary J = 0, exterior", st, lambda t, r: T.pq("I", t, r),
             ck.uniform(-15, 15), ck.uniform(1.001, 40), lambda t, r: (1, 0))
    ck.chart("BTZ stationary J = 0, black hole", st, lambda t, r: T.pq("II", t, r),
             ck.uniform(-15, 15), ck.uniform(0.01, 0.999), lambda t, r: (0, -1))
    ck.chart("BTZ stationary J = 0, white hole", st, lambda t, r: T.pq("IV", t, r),
             ck.uniform(-15, 15), ck.uniform(0.01, 0.999), lambda t, r: (0, 1))
    ck.chart("BTZ stationary J = 0, other exterior", st, lambda t, r: T.pq("I'", t, r),
             ck.uniform(-15, 15), ck.uniform(1.001, 40), lambda t, r: (-1, 0))

    def ingoing(w, r):
        w, r = np.asarray(w, dtype=float), np.asarray(r, dtype=float)
        return np.arctan((1 - r) / (1 + r) * np.exp(-w)), atan_exp(w)

    def outgoing(u, r):
        u, r = np.asarray(u, dtype=float), np.asarray(r, dtype=float)
        return -atan_exp(-u), np.arctan((r - 1) / (r + 1) * np.exp(u))
    ein = Plane(src, "btz", "eddington_finkelstein_ingoing", ("v", "r"), {"tildephi": "0"}, J0)
    ck.chart("BTZ ingoing Eddington-Finkelstein J = 0", ein, ingoing,
             ck.uniform(-15, 15), ck.uniform(0.01, 40), lambda w, r: (1, -60))
    eout = Plane(src, "btz", "eddington_finkelstein_outgoing", ("u", "r"), {"tildephi": "0"}, J0)
    ck.chart("BTZ outgoing Eddington-Finkelstein J = 0", eout, outgoing,
             ck.uniform(-15, 15), ck.uniform(0.01, 40), lambda u, r: (1, 60))

    p, q = T.pq("II", np.array([-5.0, 0, 5]), np.full(3, 1e-9))
    ck.limit("BTZ, J = 0: r -> 0 in the black hole lands on T = pi/2", p + q, [HALF] * 3)
    p, q = T.pq("I", np.array([-5.0, 0, 5]), np.full(3, 1e9))
    ck.limit("BTZ, J = 0: r -> infinity lands on the boundary X = pi/2", q - p, [HALF] * 3, 1e-8)
    p, q = T.pq("I", np.array([3.0]), np.array([1 + 1e-12]))
    ck.limit("BTZ, J = 0: r -> r_+ at fixed t lands on the bifurcation circle", point(p[0], q[0]), [0, 0], 1e-4)
    rr = np.linspace(0.05, 5, 50)
    for cell, sel in (("I", rr > 1.001), ("II", rr < 0.999)):
        pp, qq = T.pq(cell, 0.3 + 0 * rr[sel], rr[sel])
        ck.limit(f"BTZ, J = 0, {cell}: tan p tan q is Kruskal's UV = (1 - r)/(1 + r)",
                 np.tan(pp) * np.tan(qq), (1 - rr[sel]) / (1 + rr[sel]), 1e-8)
    ck.limit("BTZ, J = 0: the ingoing and stationary coordinates put one event at one point",
             ingoing(2.0 + T.rstar(3.0), 3.0), T.pq("I", 2.0, 3.0), 1e-12)
    ck.finite("BTZ, J = 0: the Kretschmann scalar is finite at r = 0 and at r_+",
              st.kretschmann(np.zeros(3), np.array([1e-9, 1.0, 5.0])))

    rot = Plane(src, "btz", "stationary", ("t", "r"), None, JR, quotient="phi")
    assert rot.g[0, 1] == 0 and sp.simplify(rot.g[0, 0] * rot.g[1, 1] + 1) == 0
    R = BTZTower(rot.gi[1, 1], rot.x1)
    rp, rm = R.rf
    tower_checks(ck, "BTZ stationary J = 4l/5", rot, R, 0.01, 15)
    ck.limit("BTZ, J = 4l/5: the horizons are the roots of the published g^rr",
             [rp * rp, rm * rm], [0.8, 0.2], 1e-12)
    ck.limit("BTZ, J = 4l/5: k_-/k_+ = r_+/r_- = 2", float(R.kappa[1] / R.kappa[0]), 2.0, 1e-12)
    p, q = R.pq("III", np.array([-6.0, 0, 6]), np.full(3, 1e-12))
    ck.limit("BTZ, J = 4l/5: r -> 0 lands on the vertical line X = pi/2", q - p, [HALF] * 3, 1e-9)
    p, q = R.pq("I", np.array([-6.0, 0, 6]), np.full(3, 1e9))
    ck.limit("BTZ, J = 4l/5: r -> infinity lands on the boundary X = pi/2", q - p, [HALF] * 3, 1e-8)
    ck.finite("BTZ, J = 4l/5: the Kretschmann scalar is finite at r = 0 and at both horizons",
              rot.kretschmann(np.zeros(3), np.array([1e-9, rm, rp])))

    views = []
    box = [-HALF - 0.55, HALF + 0.55, -HALF - 0.25, HALF + 0.25]
    square = [[HALF, -HALF], [HALF, HALF], [-HALF, HALF], [-HALF, -HALF]]
    R_OUT, R_IN, TS = (1.25, 1.5, 2, 3), (0.25, 0.5, 0.75), (-2, -1, 0, 1, 2)

    def edges(v):
        v.line("boundary", [[[HALF, -HALF], [HALF, HALF]], [[-HALF, -HALF], [-HALF, HALF]]])
        v.line("horizon", [[[-HALF, -HALF], [HALF, HALF]], [[HALF, -HALF], [-HALF, HALF]]])
        v.line("singular", [[[-HALF, HALF], [HALF, HALF]], [[-HALF, -HALF], [HALF, -HALF]]], zig=True)
        v.label_xt([0, HALF], "$r = 0$", "b", dy=-8)
        v.label_xt([0, -HALF], "$r = 0$", "t", dy=8)
        v.label_xt([HALF, 0.9], "$r \\to \\infty$", "l", "small", dx=6)
        v.label_xt([-HALF, 0.9], "$r \\to \\infty$", "r", "small", dx=-6)
        v.label_xt([-Q4, Q4], "$r_+$", "tr", "small", dx=-6, dy=2)
        v.label_xt([0.3, -0.62], "exterior", cls="region")
        v.label_xt([-0.3, 0.62], "exterior", cls="region")
        v.label_xt([0, 1.2], "black hole", cls="region")
        v.label_xt([0, -1.2], "white hole", cls="region")
        v.legend("horizon", "the horizon $r_+ = \\sqrt{M}\\,\\ell$")
        v.legend("boundary", "the conformal boundary, timelike")
        v.legend("singular", "$r = 0$, a singularity in the causal structure, where the curvature is finite")

    t = spread(-np.inf, np.inf, 500, 9)
    v = View("static", "$J = 0$", box, "stationary")
    v.fill("region", square)
    v.fill("cover", [[0, 0], [HALF, -HALF], [HALF, HALF]])
    for r in R_OUT:
        v.curve("r", *T.pq("I", t, np.full_like(t, r)))
    rr = spread(1, np.inf, 500, 14)
    for tt in TS:
        v.curve("t", *T.pq("I", np.full_like(rr, tt), rr))
    edges(v)
    for r, text in ((1.25, "$1.25\\,r_+$"), (2, "$2\\,r_+$")):
        label_on(v, T.pq("I", 0.0, r), text)
    v.legend("cover", "the region that $t$ and $r > r_+$ cover")
    v.legend("r", "$r$ constant")
    v.legend("t", "$ct$ constant, in units of $\\ell$")
    views.append(v)

    v = View("ingoing", "Ingoing Eddington-Finkelstein", box, "eddington_finkelstein_ingoing")
    v.fill("region", square)
    v.fill("cover", [[0, 0], [HALF, -HALF], [HALF, HALF], [-HALF, HALF]])
    for r in R_OUT + R_IN:
        v.curve("r", *ingoing(t, np.full_like(t, r)))
    rr = spread(0, np.inf, 600, 14)
    for w in (-3, -2, -1, 0, 1, 2, 3):
        v.curve("null", *ingoing(np.full_like(rr, w), rr))
    edges(v)
    v.legend("cover", "the region that $v$ and $r > 0$ cover")
    v.legend("r", "$r$ constant")
    v.legend("null", "$v$ constant, an ingoing light ray")
    views.append(v)

    v = View("outgoing", "Outgoing Eddington-Finkelstein", box, "eddington_finkelstein_outgoing")
    v.fill("region", square)
    v.fill("cover", [[0, 0], [HALF, HALF], [HALF, -HALF], [-HALF, -HALF]])
    for r in R_OUT + R_IN:
        v.curve("r", *outgoing(t, np.full_like(t, r)))
    for w in (-3, -2, -1, 0, 1, 2, 3):
        v.curve("null", *outgoing(np.full_like(rr, w), rr))
    edges(v)
    v.legend("cover", "the region that $u$ and $r > 0$ cover")
    v.legend("r", "$r$ constant")
    v.legend("null", "$u$ constant, an outgoing light ray")
    views.append(v)
    for view in views:
        view.set(settings="$M = 1$ and $J = 0$, so that $r_+ = \\ell$.")
    # The embedded moment t = 0: through the bifurcation circle into both exteriors, as far out
    # as the surface reaches.
    moment, = slices.moments("btz")
    lo, hi = moment.reach("stationary", "r")
    for view in views:
        view.slice(moment, [through_bifurcation(T, ("I'", "I"), hi, lo)])

    # The rotating hole: one period of the strip, from the exteriors at T = 0 to those at 2 pi.
    T0, T1 = -0.45, 2 * PI + 0.45
    v = View("rotating", "$J = 4\\ell/5$", [-HALF - 0.55, HALF + 0.55, T0, T1], "stationary")
    cells = [("I", False), ("I'", False), ("II", False), ("III", False), ("III'", False),
             ("II", True), ("I", True), ("I'", True), ("IV", False), ("IV", True)]
    for cell, up in cells:
        v.fill("region", [point(*reflect(p, q, up)) for p, q in BTZ_CELL[cell]])
    for cell, up in (("I", False), ("II", False), ("III", False)):
        v.fill("cover", [point(*reflect(p, q, up)) for p, q in BTZ_CELL[cell]])
    radii = {"I": (1.25, 2, 4), "II": (0.55, 0.7, 0.8), "III": (0.15, 0.3)}
    lo_hi = {"I": (rp, np.inf), "II": (rm, rp), "III": (0, rm)}
    times = [c / R.kp for c in (-1.6, -0.6, 0, 0.6, 1.6)]
    for cell, up in cells:
        key = cell.rstrip("'")
        if key == "IV":
            key = "II"
        for r in radii[key]:
            v.curve("r", *reflect(*R.pq(cell, t, np.full_like(t, r)), up))
        rr = spread(*lo_hi[key], 600, 16)
        for tt in times:
            v.curve("t", *reflect(*R.pq(cell, np.full_like(rr, tt), rr), up))
    for base, top in ((0, False), (2 * PI, True)):
        for sx in (1, -1):
            v.segment("boundary", reflect(*((0, HALF) if sx > 0 else (HALF, 0)), top),
                      reflect(*((-HALF, 0) if sx > 0 else (0, -HALF)), top))
    for sx in (1, -1):
        a, b = ((0, HALF), (HALF, PI)) if sx > 0 else ((HALF, 0), (PI, HALF))
        v.segment("singular", a, b, zig=True)
    for up in (False, True):
        for a, b in (((0, 0), (HALF, 0)), ((HALF, 0), (HALF, HALF)), ((HALF, HALF), (0, HALF)), ((0, HALF), (0, 0)),
                     ((0, 0), (-HALF, 0)), ((0, 0), (0, -HALF))):
            v.segment("horizon", reflect(*a, up), reflect(*b, up))
    # Each region's name stands where its wedge is wide enough for it, well out from the vertex.
    for base in (0, 2 * PI):
        for sx in (1, -1):
            v.label_xt([sx * 0.95, base], "exterior", cls="region")
    v.label_xt([0, HALF + 0.1], "black hole", cls="region")
    v.label_xt([0, 1.5 * PI - 0.1], "white hole", cls="region")
    for sx in (1, -1):
        v.label_xt([sx * 0.95, PI], "$r < r_-$", cls="region")
        v.label_xt([sx * HALF, PI], "$r = 0$", "l" if sx > 0 else "r", dx=8 * sx)
        for base in (0, 2 * PI):
            v.label_xt([sx * HALF, base + (0.25 if base == 0 else -0.25)], "$r \\to \\infty$",
                       "l" if sx > 0 else "r", "small", dx=6 * sx)
    v.label_xt([Q4, Q4], "$r_+$", "tl", "small", dx=5, dy=1)
    v.label_xt([Q4, 3 * Q4], "$r_-$", "bl", "small", dx=5, dy=-1)
    clip_in_t(v, T0, T1)
    v.set(fade={"top": 0.9, "bottom": 0.9},
          settings="$M = 1$ and $J = 4\\ell/5$, so that $r_+ = 2\\ell/\\sqrt{5}$, $r_- = \\ell/\\sqrt{5}$, "
                   "and $\\kappa_-/\\kappa_+ = 2$.")
    v.legend("cover", "one exterior, one region between the horizons and one inside $r_-$, which $t$ and "
                      "$r > 0$ cover")
    v.legend("r", f"$r$ constant: {listed(radii['I'])} outside, {listed(radii['II'])} between, "
                  f"{listed(radii['III'])} inside $r_-$, in units of $\\ell$")
    v.legend("t", "$t$ constant")
    v.legend("horizon", "the horizons $r_+ = 2\\ell/\\sqrt{5}$ and $r_- = \\ell/\\sqrt{5}$")
    v.legend("boundary", "the conformal boundary, timelike")
    v.legend("singular", "$r = 0$, a singularity in the causal structure, where the curvature is finite")
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
        # The equator's moment t = 0 over the stretch the embedding reaches, and the sphere at
        # the event r = b, x = b.
        equator = slices.moments("bertotti_robinson", "equator")[0]
        sphere = slices.moments("bertotti_robinson", "sphere")[0]
        lo, hi = equator.reach("static", "r")
        c = np.linspace(lo, hi, 2) if vid == "static" else 1 / np.linspace(hi, lo, 2)
        v.slice(equator, [poincare_pq(0 * c, z_of(c))])
        v.slice(sphere, points=[poincare_pq(0.0, 1.0)],
                label="$t = 0$, $r = b$" if vid == "static" else "$t = 0$, $x = b$")
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
    wormhole = slices.moments("ellis_bronnikov")[0]
    r = np.linspace(*wormhole.reach("spherical", "r"), 2)
    v.slice(wormhole, [mink_pq(0 * r, r)])
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
        wormhole = slices.moments("morris_thorne")[0]
        ell = np.sqrt(wormhole.reach("spherical", "r")[1] ** 2 - 1)
        ls = np.array([-ell, 0.0, ell])
        v.slice(wormhole, [mink_pq(0 * ls, ls)])
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

    # Wider to the left than Gott's core, so that the string's name beside it stays on the drawing.
    v = View("conical", "Conical exterior", [-0.5] + box[1:], "conical")
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
    v.set(restriction="The half plane of fixed $\\phi$ and $z$ only, totally geodesic, each point in "
                      "the diagram a circle around the string times a line along it, the string's "
                      "deficit angle $\\delta = 8\\pi G\\mu/c^2$ showing in those circles.")
    # The ideal string's moment is the embedding's reference cone and the sheet outside the
    # core together, out from the string; with the core, its proper distance from the axis
    # runs out to l chi_0 + r - l tan chi_0.
    # The ideal string's cone unrolled is the same moment, out from the string.
    for m in slices.moments("cosmic_string"):
        r = np.linspace(*m.reach("conical", "r", reference=True), 2)
        v.slice(m, [mink_pq(0 * r, r)])
    cone = slices.moments("cosmic_string", "cone")[0]
    reach = cone.reach("conical", "r", reference=True)
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
                      "from the axis, each point in the diagram a circle around the axis times a "
                      "line along it.")
    if abs(cone.reach("interior_cap", "\\chi")[1] - chi0) > 1e-12:
        raise SystemExit("cosmic string: the embedding's core is not the core drawn here")
    rho = np.array([0.0, chi0, reach[1] + shift])
    v.slice(cone, [mink_pq(0 * rho, rho)])
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
    star_moment = slices.moments("interior_schwarzschild")[0]
    lo, hi = star_moment.reach("spherical", "r")
    rr = np.concatenate([np.linspace(lo, R, 40), R + np.geomspace(1e-6, hi - R, 200)])
    v.slice(star_moment, [star(0 * rr, rr)])
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
    for m in slices.moments("frw", "closed"):
        eta = float(np.interp(m.time, etas - np.sin(etas), etas))
        if not (m.reach("comoving_spherical", "r") == (0, 1) and abs(eta - np.sin(eta) - m.time) < 1e-9):
            raise SystemExit(f"FRW: the moment {m.label} is not a whole moment of the closed universe")
        # X = chi across the near hemisphere and the far one, as the view draws its lines of eta.
        v.slice(m, xt=[[[0.0, eta], [PI, eta]]])
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
        view.set(input="Dust: $a \\propto \\eta^2$, $1 - \\cos\\eta$, and $\\cosh\\eta - 1$ for $k = 0$, "
                       "$+1$, and $-1$, each solved from this spacetime's own $G^r{}_r = 0$, with lengths in "
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
    v.label_xt([Xmax / 2, 0], "$\\tau = 0$ inside, $t = 0$ outside", "t", "coord", dy=7)
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
    # Each moment of the dust's proper time: inside, the line eta across the dust; outside,
    # Novikov's slice, the shells released from rest with it at every R the embedding reaches,
    # each at the same proper time, through Kruskal's U and V to P(U) and Q(V). The surface
    # is the shell from R_0, so the two meet there, at (chi_0, eta).
    ck.limit("Oppenheimer-Snyder: slices.novikov_kruskal puts the surface where the Collapse does",
             [float(o.P(slices.novikov_kruskal(np.array([o.R0]), (o.am / 2) * (e + np.sin(e)))[0])[0])
              for e in (0.5, 1.5, 2.5)], [(e - c) / 2 for e in (0.5, 1.5, 2.5)], 1e-6)
    for m in slices.moments("oppenheimer_snyder"):
        if abs(m.reach("interior_comoving", "\\chi")[1] - c) > 1e-12 or abs(slices.OS_AM - o.am) > 1e-12:
            raise SystemExit("Oppenheimer-Snyder: the embedding's star is not the star drawn here")
        e = float(eta_of(m.time))
        chi = np.array([0.0, c])
        lo, hi = m.reach("comoving_synchronous", "r")
        U, V, _ = slices.novikov_kruskal(np.linspace(lo, hi, 401), m.time)
        # Two lines, so that the corner where they meet on the surface is kept exactly.
        v.slice(m, [((e - chi) / 2, (e + chi) / 2), (o.P(U), o.Q(V))])
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
    # Each moment v - r = w through each side's own map, meeting on the shell at r = -w.
    for m in slices.moments("vaidya"):
        lo, hi = m.reach("eddington_finkelstein_ingoing", "r")
        w, cross = m.time, max(lo, -m.time)
        r_in = np.linspace(lo, cross, 200) if cross > lo else np.zeros(0)
        r_out = cross + (hi - cross) * np.linspace(0, 1, 400) ** 2
        p_in, q_in = inside(w + r_in, r_in)
        p_out, q_out = outside(w + r_out, r_out)
        v.slice(m, [(np.concatenate([p_in, p_out]), np.concatenate([q_in, q_out]))])
    return [v]


# ---------------------------------------------------------------- Tolman-Oppenheimer-Volkoff

def tov(ck, src):
    """A static star of fluid, its redshift and mass functions solved for a declared equation of
    state from this spacetime's own Einstein tensor, G = c = M_sun = 1, by null_rays.StarSolver,
    which the spacetime diagram draws the same star with.

    G^t_t = -8 pi rho and G^r_r = 8 pi p, read from the published components with Phi and m
    as plain symbols, are linear in m' and Phi', and give them; the pressure follows from
    the conservation law p' = -(rho + p) Phi', which the Bianchi identity makes the same
    statement as G^theta_theta = G^r_r, so the published G^theta_theta is checked along the
    solution rather than used. The equation of state is the polytrope p = K rho_0^2 with
    energy density rho = rho_0 + p, at K = 100 and central rho_0 = 1.28e-3, the star on
    which numerical relativists test their codes. Outside the surface, where p = 0, m = M and
    e^(2 Phi) = 1 - 2M/r, the Schwarzschild exterior in the same coordinates, and Phi inside
    is shifted to meet it. Then r* = int e^(-Phi) (1 - 2m/r)^(-1/2) dr runs from 0 at the
    centre to infinity and p, q = arctan((t -+ r*)/R) give Minkowski's triangle.
    """
    solver = nr.StarSolver("tov", "spherical", 100.0, 1.28e-3)
    pl = Plane(src, "tov", "spherical", ("t", "r"), EQUATOR, numeric=["Phi", "m"])
    src.note("tov", "spherical", ["einstein_tensor"])
    R, M, r0, values = solver.R, solver.M, solver.r0, solver.values

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
    ck.limit("Tolman-Oppenheimer-Volkoff: the declared star solves the published G^theta_theta = 8 pi p",
             solver.theta_theta(np.linspace(0.05 * R, 0.95 * R, 200)), 0, 1e-7)
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
    star_moment = slices.moments("tov")[0]
    lo, hi = star_moment.reach("spherical", "r")
    rr = np.concatenate([np.linspace(lo, R, 40), R + np.geomspace(1e-6, hi - R, 200)])
    v.slice(star_moment, [star(0 * rr, rr)])
    v.set(input="A polytrope, $p = K\\rho_0^2$ with rest mass density $\\rho_0$ and energy density "
                "$\\rho c^2 = \\rho_0c^2 + p$, at $K = 100$ and a central $\\rho_0 = 1.28\\times10^{-3}$ in "
                "units where $G = c = M_\\odot = 1$, solved from this spacetime's own $G^t{}_t$ and "
                f"$G^r{{}}_r$: a star of $M = {M:.2f}\\,M_\\odot$ and $R = {R * km:.1f}$ km, "
                "the star on which numerical relativists test their codes.")
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
    for m in slices.moments("malament_hogarth"):
        x = np.linspace(*m.reach("cartesian", "x"), 400)
        v.slice(m, [mink_pq(np.full_like(x, m.time), x)])
    return [v]


def einstein_rosen_waves(ck, src):
    """The half plane of fixed phi and z, totally geodesic, whose metric is
    e^{2(gamma - psi)}(-c^2dt^2 + drho^2) in the cylindrical chart and -e^{2(gamma - psi)} du dv in
    the null chart, u, v = ct -+ rho. A conformal factor changes no null direction, so for every
    wave the causal structure of the half plane is Minkowski's half diamond,
    p, q = arctan((ct -+ rho)/a) = arctan(u/a), arctan(v/a), with the axis rho = 0 on X = 0. The
    maps are checked with random values of psi and gamma at every sample, and with the pulse of
    Weber, Wheeler, and Bonnor at C = a = 1, whose published Kretschmann scalar is checked to settle to a
    finite value on the axis; its rays through the event on the axis where it is greatest, rho = |ct|, are
    the lines q = 0 below and p = 0 above."""
    numeric = ["psi", "gamma"]

    def any_wave(x0, x1):
        values = {}
        for name in numeric:
            w = ck.rng.uniform(-3, 3, np.shape(x0))
            values[name] = (w, 0 * w, 0 * w)
        return values
    cyl = Plane(src, "einstein_rosen_waves", "cylindrical", ("t", "\\rho"), {"phi": "0", "z": "0"}, numeric=numeric)
    ck.chart("Einstein-Rosen cylindrical, for any wave", cyl, mink_pq, ck.uniform(-20, 20), ck.uniform(0.01, 20),
             lambda t, r: (1, 0), any_wave)
    pulse = nr._er_pulse("t", "rho")
    wave = Plane(src, "einstein_rosen_waves", "cylindrical", ("t", "\\rho"), {"phi": "0", "z": "0"}, functions=pulse)
    ck.chart("Einstein-Rosen cylindrical, the pulse", wave, mink_pq, ck.uniform(-20, 20), ck.uniform(0.01, 20),
             lambda t, r: (1, 0))
    # The published Kretschmann scalar with the pulse, differentiated in sympy and not simplified.
    _, entry, reader = nr.load("einstein_rosen_waves", "cylindrical")
    src.note("einstein_rosen_waves", "cylindrical", ["kretschmann"])
    K = reader(nr.strip_lhs(entry["kretschmann"]))
    for name, text in pulse.items():
        K = K.replace(reader.parameters[name].func, nr._as_lambda(reader, name, text)).doit()
    K = sp.lambdify((reader.symbol["t"], reader.symbol["\\rho"]), K.subs({reader.c: 1, reader.symbol["\\phi"]: 0,
                                                                        reader.symbol["z"]: 0}), "numpy")
    # On the axis at t = 0 the pulse has psi = 2, and K carries e^{4(psi - gamma)}: it settles near
    # 5.7e5 as rho -> 0, large in units of a but finite.
    t_axis = ck.uniform(-5, 5, 50)
    ck.settles("Einstein-Rosen: the axis rho = 0 is regular for the pulse",
               K(t_axis, np.full(50, 1e-3)), K(t_axis, np.full(50, 1e-4)))

    def null_map(u, w):
        return np.arctan(np.asarray(u, dtype=float)), np.arctan(np.asarray(w, dtype=float))
    null = Plane(src, "einstein_rosen_waves", "null", ("u", "v"), {"phi": "0", "z": "0"}, numeric=numeric)
    u = ck.uniform(-20, 20)
    ck.chart("Einstein-Rosen null, for any wave", null, null_map, u, u + ck.uniform(0.01, 30), lambda u, w: (1, 1),
             any_wave)

    box = [-0.35, PI + 0.35, -PI - 0.25, PI + 0.25]
    TS, RS = (-4, -2, -1, 0, 1, 2, 4), (0.5, 1, 2, 4)
    restriction = ("The half plane of fixed $\\phi$ and $z$ only, totally geodesic, each point in the diagram a "
                   "circle around the axis times a line along it.")
    moments = slices.moments("einstein_rosen_waves")
    views = []
    for system in ("cylindrical", "null"):
        v = View(system, {"cylindrical": "Cylindrical", "null": "Null"}[system], box, system)
        v.fill("region", TRIANGLE)
        v.fill("cover", TRIANGLE)
        if system == "cylindrical":
            grid(v, "r", lambda r, t: mink_pq(t, r), RS, S_ALL)
            grid(v, "t", mink_pq, TS, S_POS)
            v.legend("cover", "the whole spacetime, which $t$ and $\\rho$ cover")
            v.legend("r", "$\\rho$ constant, in units of $a$")
            v.legend("t", "$ct$ constant")
        else:
            for c in TS:
                s = np.linspace(c, 60, 400)
                v.curve("null", *null_map(np.full_like(s, c), s))
                s = np.linspace(-60, c, 400)
                v.curve("null", *null_map(s, np.full_like(s, c)))
            v.legend("cover", "the whole spacetime, which $u$ and $v$ cover")
            v.legend("null", "$u$ constant and $v$ constant, every one a light ray")
        triangle_edges(v, centre="$\\rho = 0$" if system == "cylindrical" else "$u = v$")
        v.segment("cone", (-HALF, 0), (0, 0))
        v.segment("cone", (0, 0), (0, HALF))
        v.point("mark", (0, 0))
        v.label((0, 0), "$\\psi = 2C/a$", "l", "small", dx=8)
        rays = "$\\rho = |ct|$" if system == "cylindrical" else "$v = 0$ and $u = 0$"
        v.legend("cone", f"the rays {rays} through that event, just inside the crest of the pulse")
        v.legend("mark", "the event on the axis where the pulse is greatest")
        v.legend("centre", "the axis $\\rho = 0$, regular where $\\gamma = 0$" if system == "cylindrical"
                 else "the axis $u = v$, regular where $\\gamma = 0$")
        v.set(restriction=restriction, input=nr.ER_INPUT)
        for m in moments:
            r = np.linspace(*m.reach("cylindrical", "\\rho"), 2)
            v.slice(m, [mink_pq(np.full_like(r, m.time), r)])
        views.append(v)
    return views


def melvin(ck, src):
    """Melvin's half plane of fixed phi and z, and Ernst's equator.

    The metric on Melvin's half plane is Lambda^2(-c^2dt^2 + drho^2), Lambda = 1 + B^2 rho^2/4, and the
    factor changes no null direction, so p, q = arctan((ct -+ rho)B) bring it into Minkowski's half
    diamond with the axis on X = 0; the affine parameter along a ray grows as the integral of Lambda^2,
    without bound, so the far edges are the half plane's null infinity. Ernst's equator has
    Lambda^2(-(1 - r_s/r)c^2dt^2 + dr^2/(1 - r_s/r)), Lambda = 1 + B^2 r^2/4, conformal to Schwarzschild's
    plane of t and r, so Kruskal and Szekeres's cells, the tower of the one root r_s, draw it, checked
    against the published metric at B = 1/(2 r_s); the same null curves serve every plane of constant theta,
    and on the equator and the axis they are geodesics."""
    cyl = Plane(src, "melvin", "cylindrical", ("t", "\\rho"), {"phi": "0", "z": "0"}, {"B": 1})
    ck.chart("Melvin cylindrical", cyl, mink_pq, ck.uniform(-20, 20), ck.uniform(0.01, 20), lambda t, r: (1, 0))
    K = cyl.kretschmann
    ck.finite("Melvin: the Kretschmann scalar is finite on the axis rho = 0",
              K(ck.uniform(-5, 5, 50), np.full(50, 1e-6)))
    ck.limit("Melvin: the Kretschmann scalar on the axis is 20 B^4", K(np.zeros(1), np.full(1, 1e-9)), [20.0], 1e-6)

    ern = Plane(src, "melvin", "ernst", ("t", "r"), EQUATOR, {"r_s": 1, "B": "1/2"})
    T = Tower(1 - 1 / ern.x1, ern.x1, [1])
    ck.chart("Ernst equator, exterior", ern, lambda t, r: T.pq("I", t, r),
             ck.uniform(-15, 15), ck.uniform(1.001, 30), lambda t, r: (1, 0))
    ck.chart("Ernst equator, black hole", ern, lambda t, r: T.pq("II", t, r),
             ck.uniform(-15, 15), ck.uniform(0.01, 0.999), lambda t, r: (0, -1))
    ck.chart("Ernst equator, white hole", ern, lambda t, r: T.pq("IV", t, r),
             ck.uniform(-15, 15), ck.uniform(0.01, 0.999), lambda t, r: (0, 1))
    ck.chart("Ernst equator, other exterior", ern, lambda t, r: T.pq("I'", t, r),
             ck.uniform(-15, 15), ck.uniform(1.001, 30), lambda t, r: (-1, 0))
    p, q = T.pq("II", np.array([-5.0, 0, 5]), np.full(3, 1e-9))
    ck.limit("Ernst: r -> 0 in the black hole lands on T = pi/2", p + q, [HALF] * 3)
    KE = ern.kretschmann
    ck.diverges("Ernst: the Kretschmann scalar diverges at r = 0", KE(0, 1e-2), KE(0, 1e-3))
    ck.finite("Ernst: the Kretschmann scalar is finite at r = r_s", KE(np.zeros(3), np.array([0.999, 1, 1.001])))

    views = []
    box = [-0.35, PI + 0.35, -PI - 0.25, PI + 0.25]
    TS, RS = (-4, -2, -1, 0, 1, 2, 4), (0.5, 1, 4)
    universe = slices.moments("melvin", "universe")[0]
    v = View("cylindrical", "Cylindrical", box, "cylindrical")
    v.fill("region", TRIANGLE)
    v.fill("cover", TRIANGLE)
    grid(v, "r", lambda r, t: mink_pq(t, r), RS, S_ALL)
    grid(v, "surface", lambda r, t: mink_pq(t, r), (2,), S_ALL)
    grid(v, "t", mink_pq, TS, S_POS)
    triangle_edges(v, centre="$\\rho = 0$")
    v.legend("cover", "the whole spacetime, which $t$ and $\\rho$ cover")
    v.legend("r", "$\\rho$ constant, in units of $1/B$")
    v.legend("surface", "the Melvin radius $\\rho = 2/B$, where the circles about the axis are widest")
    v.legend("t", "$ct$ constant")
    v.legend("centre", "the axis $\\rho = 0$")
    r = np.linspace(*universe.reach("cylindrical", "\\rho"), 2)
    v.slice(universe, [mink_pq(0 * r, r)])
    v.set(restriction="The half plane of fixed $\\phi$ and $z$ only, totally geodesic, each point in the diagram a "
                      "circle around the axis times a line along it.",
          settings="$B = 1$, the scale of $p = \\arctan((ct - \\rho)B)$ and $q = \\arctan((ct + \\rho)B)$.")
    views.append(v)

    box = [-PI - 0.25, PI + 0.25, -HALF - 0.25, HALF + 0.25]
    hexagon = [[PI, 0], [HALF, HALF], [-HALF, HALF], [-PI, 0], [-HALF, -HALF], [HALF, -HALF]]
    exterior = [[0, 0], [HALF, -HALF], [PI, 0], [HALF, HALF]]
    t = spread(-np.inf, np.inf, 500, 9)
    v = View("ernst", "Ernst", box, "ernst")
    v.fill("region", hexagon)
    v.fill("cover", exterior)
    for rr in (1.25, 2, 3, 6):
        v.curve("r", *T.pq("I", t, np.full_like(t, rr)))
    v.curve("surface", *T.pq("I", t, np.full_like(t, 4.0)))
    rr = spread(1, np.inf, 500, 14)
    for tt in TS:
        v.curve("t", *T.pq("I", np.full_like(rr, tt), rr))
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
    label_on(v, T.pq("I", 0.0, 4.0), "$r = 2/B$")
    v.legend("cover", "the region that $t$ and $r > r_s$ cover")
    v.legend("r", "$r$ constant, at $1.25$, $2$, $3$ and $6\\,r_s$")
    v.legend("surface", "$r = 2/B$, the widest circle of the equator")
    v.legend("t", "$ct$ constant, in units of $r_s$")
    v.legend("horizon", "the horizon $r = r_s$")
    v.legend("singular", "$r = 0$, where the Kretschmann scalar diverges")
    v.legend("scri", "null infinity of the plane, $\\mathscr{I}^\\pm$")
    moment = slices.moments("melvin", "ernst")[0]
    lo, hi = moment.reach("ernst", "r")
    rm = np.linspace(lo, hi, 2)
    a, b = T.pq("I'", 0 * rm, rm[::-1]), T.pq("I", 0 * rm, rm)
    v.slice(moment, [(np.concatenate([a[0], b[0]]), np.concatenate([a[1], b[1]]))])
    v.set(restriction="The plane of $t$ and $r$ at $\\theta = \\pi/2$ and $\\phi = 0$ only, totally geodesic, and "
                      "the same drawing for every plane of constant $\\theta$.",
          settings="$r_s = 1$ and $B = 1/(2r_s)$; the drawing is the same for every $B$.")
    views.append(v)
    return views


def published_gthth(src, metric_id, system_id, params):
    """The published g_thetatheta of a spherical chart, c = 1, as a numpy function of (t, r) on
    the equator."""
    _, entry, reader = nr.load(metric_id, system_id)
    src.note(metric_id, system_id, BASE_FIELDS)
    g = nr.published_matrix(reader, entry, "metric_components")
    i = entry["coords"].index("\\theta")
    names = {reader._plain(n): s for n, s in reader.symbol.items()}
    e = g[i, i].subs({reader.c: 1, **{reader.parameters[k]: sp.sympify(v) for k, v in params.items()}})
    e = e.subs(names["theta"], sp.pi / 2)
    x0, x1 = (reader.symbol[c] for c in entry["coords"][:2])
    return sp.lambdify((x0, x1), e, "numpy")


def covered(fmap, reach, eta):
    """The values of chi that a spacetime drawn in the strip of the Einstein static universe covers
    at the conformal time eta, as (lowest, highest), or None where it covers none. fmap draws it
    from its own chart, whose time runs over the whole line and whose radial coordinate over
    `reach`. At each of 4001 values of the radial coordinate the time at which the map reaches eta
    is found by bisection, since eta = p + q grows with the time, and the values of chi = q - p
    that reach it are kept. An unbounded reach is sampled out to 10^8, where Minkowski space and
    anti-de Sitter space fall short of their edges in chi by less than 10^-7."""
    lo, hi = reach
    r = np.concatenate([[lo], lo + np.geomspace(1e-9, 1e8, 4000)]) if np.isinf(hi) else np.linspace(lo, hi, 4001)
    a, b = np.full_like(r, -1e10), np.full_like(r, 1e10)
    with np.errstate(over="ignore"):
        for _ in range(200):
            m = (a + b) / 2
            p, q = fmap(m, r)
            below = p + q < eta
            a, b = np.where(below, m, a), np.where(below, b, m)
        p, q = fmap((a + b) / 2, r)
    hit = np.abs(p + q - eta) < 1e-7
    if not hit.any():
        return None
    chi = (q - p)[hit]
    return float(chi.min()), float(chi.max())


def across(polygon, T):
    """The least and greatest X at which the line of constant T crosses a convex polygon of (X, T),
    or None where it misses it."""
    P = np.asarray(polygon, dtype=float)
    xs = []
    for (x0, t0), (x1, t1) in zip(P, np.roll(P, -1, axis=0)):
        if min(t0, t1) <= T <= max(t0, t1):
            xs += [x0, x1] if t0 == t1 else [x0 + (T - t0) * (x1 - x0) / (t1 - t0)]
    return (min(xs), max(xs)) if xs else None


# The spacetimes the Einstein static universe's conformal diagram draws inside its strip, by the
# id of the view that draws each: its name, the map from its own chart into the strip at R = 1,
# and that chart's reach in its radial coordinate, its time running over the whole line.
# Minkowski space and anti-de Sitter space enter by the charts the conformal check pulls back;
# de Sitter space by its closed slicing, since its static patch covers only the diamond
# |eta| + chi < pi/2 about the pole.
ESU_EMBEDDED = {
    "minkowski": ("Minkowski space", mink_pq, (0.0, np.inf)),
    "de_sitter": ("de Sitter space", ds_closed_pq, (0.0, PI)),
    "anti_de_sitter": ("anti-de Sitter space", ads_global_pq, (0.0, np.inf)),
}
# The conformal times at which each region is checked against Hawking and Ellis's, none on an
# edge of a region.
ESU_ETAS = tuple(float(e) for e in np.arange(-3, 3.01, 0.5))


def einstein_static(ck, src):
    """The strip. In the hyperspherical chart the metric on the plane of t and chi is
    -c^2dt^2 + R^2dchi^2, flat as it stands, so with eta = ct/R the map p, q = (eta -+ chi)/2 draws
    it with X = chi and T = eta: the strip 0 <= chi <= pi, the pole on X = 0 and its antipode on
    X = pi, each point a 2-sphere of radius R sin chi. The areal chart, r = R sin chi, and
    Einstein's projection, |x| = R sin chi on the plane y = z = 0, enter it through
    chi = arcsin(r/R) and cover the half chi < pi/2.

    Minkowski space, de Sitter space and anti-de Sitter space are drawn inside it by the maps
    their own conformal diagrams use, whose X and T are this chi and eta. For each, the published
    metric of the Einstein static universe pulled back through the map is checked to be Omega^2
    times the other spacetime's published metric on its plane of t and r, with the one Omega^2
    that also carries the sphere, sin^2 chi = Omega^2 g_thetatheta, so that the whole of the other
    spacetime is conformal to the region drawn.
    """
    es = Plane(src, "einstein_static", "hyperspherical", ("t", "\\chi"), EQUATOR, {"R": 1})
    ar = Plane(src, "einstein_static", "static_areal", ("t", "r"), EQUATOR, {"R": 1})
    ca = Plane(src, "einstein_static", "einstein_cartesian", ("t", "x"), {"y": "0", "z": "0"}, {"R": 1})

    def hyper(t, chi):
        t, chi = np.asarray(t, dtype=float), np.asarray(chi, dtype=float)
        return (t - chi) / 2, (t + chi) / 2

    def areal(t, r):
        return hyper(t, np.arcsin(np.abs(np.asarray(r, dtype=float))))
    ck.chart("Einstein static hyperspherical", es, hyper, ck.uniform(-10, 10), ck.uniform(0.001, PI - 0.001),
             lambda t, c: (1, 0))
    ck.chart("Einstein static areal", ar, areal, ck.uniform(-10, 10), ck.uniform(0.001, 0.999), lambda t, r: (1, 0))
    ck.chart("Einstein static, Einstein's projection", ca, areal, ck.uniform(-10, 10), ck.uniform(0.001, 0.999),
             lambda t, x: (1, 0))
    ck.finite("Einstein static: the pole chi = 0 is regular", es.kretschmann(ck.uniform(-5, 5, 50), np.full(50, 1e-6)))
    ck.finite("Einstein static: the antipode chi = pi is regular",
              es.kretschmann(ck.uniform(-5, 5, 50), np.full(50, PI - 1e-6)))
    ck.finite("Einstein static: the equator r = R is regular",
              ar.kretschmann(ck.uniform(-5, 5, 50), np.full(50, 1 - 1e-9)))

    def conformal(name, metric_id, system_id, params, fmap, r_lo, r_hi):
        other = Plane(src, metric_id, system_id, ("t", "r"), EQUATOR, params)
        ck.chart(f"{name} in the Einstein static universe", other, fmap, ck.uniform(-5, 5), ck.uniform(r_lo, r_hi),
                 lambda t, r: (1, 0))
        t, r = ck.uniform(-5, 5), ck.uniform(r_lo, r_hi)
        h = 1e-6

        def image(t, r):
            p, q = fmap(t, r)
            return p + q, q - p
        eta, chi = image(t, r)
        e0, c0 = [(a - b) / (2 * h) for a, b in zip(image(t + h, r), image(t - h, r))]
        e1, c1 = [(a - b) / (2 * h) for a, b in zip(image(t, r + h), image(t, r - h))]
        G00, G01, G11, *_ = es.metric(eta, chi)

        def pull(a0, a1, b0, b1):
            return G00 * a0 * b0 + G01 * (a0 * b1 + a1 * b0) + G11 * a1 * b1
        g00, g01, g11, *_ = other.metric(t, r)
        omega2 = np.sin(chi) ** 2 / published_gthth(src, metric_id, system_id, params)(t, r)
        scale = omega2 * np.maximum.reduce([np.abs(g00), np.abs(g01), np.abs(g11)])
        err = np.maximum.reduce([np.abs(pull(e0, c0, e0, c0) - omega2 * g00), np.abs(pull(e0, c0, e1, c1) - omega2 * g01),
                                 np.abs(pull(e1, c1, e1, c1) - omega2 * g11)]) / scale
        ck.limit(f"Einstein static: {name} is conformal to the region drawn, plane and sphere with one factor",
                 float(np.max(err)), 0, 1e-6)

    conformal("Minkowski space", "minkowski", "spherical", {}, mink_pq, 0.001, 20)
    conformal("de Sitter space", "de_sitter", "static_spherical", {"Lambda": 3}, ds_static_pq, 0.001, 0.999)
    conformal("anti-de Sitter space", "anti_de_sitter", "static_global", {"L": 1}, ads_global_pq, 0.001, 50)

    T0, T1 = -PI - 0.3, PI + 0.3
    box = [-0.35, PI + 0.35, T0, T1]
    strip_ = [[0, T0], [PI, T0], [PI, T1], [0, T1]]
    half = [[0, T0], [HALF, T0], [HALF, T1], [0, T1]]
    moment = slices.moments("einstein_static")[0]
    if moment.reach("hyperspherical", "\\chi") != (0, PI):
        raise SystemExit("Einstein static: the embedding is not the whole moment, pole to antipode")

    # Each region drawn in the strip is the one its map covers, at every conformal time checked,
    # which is the part of the sphere the embedding diagram shades while the region's view is shown.
    covers = {"minkowski": TRIANGLE, "de_sitter": [[0, -HALF], [PI, -HALF], [PI, HALF], [0, HALF]],
              "anti_de_sitter": half}
    for vid, (name, fmap, reach) in ESU_EMBEDDED.items():
        err = []
        for eta in ESU_ETAS:
            got, want = covered(fmap, reach, eta), across(covers[vid], eta)
            err.append(0.0 if got is None and want is None else np.inf if got is None or want is None
                       else max(abs(got[0] - want[0]), abs(got[1] - want[1])))
        ck.limit(f"Einstein static: the region drawn for {name} is the region its map covers, at every eta checked",
                 err, 0, 1e-6)

    def frame(v, times=True):
        v.fill("region", strip_)
        v.line("centre", [[[0, T0], [0, T1]], [[PI, T0], [PI, T1]]])
        if times:
            for k in (-2, -1, 1, 2):
                v.line("t", [[[0, k * HALF], [PI, k * HALF]]])
        v.slice(moment, xt=[[[0.0, 0.0], [PI, 0.0]]])
        v.set(fade={"top": 0.7, "bottom": 0.7})

    views = []
    v = View("hyperspherical", "Hyperspherical", box, "hyperspherical")
    frame(v)
    v.fill("cover", strip_)
    for chi in (Q4, HALF, 3 * Q4):
        v.line("r", [[[chi, T0], [chi, T1]]])
    v.line("null", [[[0, -PI], [PI, 0], [0, PI]]])
    v.label_xt([0, 2.3], "pole", "r", "coord", dx=-6)
    v.label_xt([PI, 2.3], "antipode", "l", "coord", dx=6)
    v.legend("cover", "the whole spacetime, which $t$ and $\\chi$ cover")
    v.legend("r", "$\\chi$ constant, at $\\pi/4$, $\\pi/2$ and $3\\pi/4$")
    v.legend("t", "$ct$ constant, every $\\pi R/2$")
    v.legend("null", "a light ray from the pole, at the antipode after $\\pi R/c$ and back after $2\\pi R/c$")
    v.legend("centre", "the pole $\\chi = 0$ and its antipode $\\chi = \\pi$")
    views.append(v)

    for vid, label, system, what, edge in (
            ("areal", "Areal", "static_areal", "$r$", "$r = R$, the equator, where the areal chart ends"),
            ("einstein_cartesian", "Einstein's Cartesian", "einstein_cartesian", "$\\sqrt{x^2 + y^2 + z^2}$",
             "$x^2 + y^2 + z^2 = R^2$, the equator, where Einstein's coordinates end")):
        v = View(vid, label, box, system)
        frame(v)
        v.fill("cover", half)
        for r in (0.5, math.sqrt(3) / 2):
            chi = float(np.arcsin(r))
            v.line("r", [[[chi, T0], [chi, T1]]])
            v.line("r2", [[[PI - chi, T0], [PI - chi, T1]]])
        v.line("chartedge", [[[HALF, T0], [HALF, T1]]])
        v.label_xt([0, 2.3], "pole", "r", "coord", dx=-6)
        v.label_xt([PI, 2.3], "antipode", "l", "coord", dx=6)
        v.legend("cover", "the hemisphere $\\chi < \\pi/2$, which " + what + " covers")
        v.legend("r", what + " constant, at $R/2$ and $\\sqrt{3}\\,R/2$")
        v.legend("r2", "the same spheres on the far hemisphere, $\\chi > \\pi/2$")
        v.legend("t", "$ct$ constant, every $\\pi R/2$")
        v.legend("chartedge", edge)
        v.legend("centre", "the pole $\\chi = 0$ and its antipode $\\chi = \\pi$")
        views.append(v)

    v = View("minkowski", "Minkowski space", box)
    frame(v, times=False)
    v.fill("cover", covers["minkowski"])
    rr = spread(0, np.inf, 500, 12)
    s = S_ALL
    for r in (0.5, 1, 2):
        v.curve("r", *mink_pq(s, np.full_like(s, r)))
    for tm in (-1, 1):
        v.curve("t", *mink_pq(np.full_like(rr, tm), rr))
    triangle_edges(v, centre=None)
    v.legend("cover", "Minkowski space, conformal to the triangle $|\\eta| + \\chi < \\pi$")
    v.legend("r", "Minkowski's $r$ constant, at $R/2$, $R$ and $2R$")
    v.legend("t", "Minkowski's $ct$ constant, at $\\pm R$")
    v.legend("centre", "the pole $\\chi = 0$, Minkowski's $r = 0$, and the antipode")
    views.append(v)

    v = View("de_sitter", "de Sitter space", box)
    frame(v, times=False)
    v.fill("cover", covers["de_sitter"])
    v.line("scri", [[[0, HALF], [PI, HALF]], [[0, -HALF], [PI, -HALF]]])
    v.label_xt([HALF, HALF], "$\\mathscr{I}^+$", "b", dy=-5)
    v.label_xt([HALF, -HALF], "$\\mathscr{I}^-$", "t", dy=5)
    v.legend("cover", "de Sitter space, conformal to $-\\pi/2 < \\eta < \\pi/2$")
    v.legend("centre", "the pole $\\chi = 0$ and its antipode $\\chi = \\pi$")
    v.legend("scri", "de Sitter's future and past infinity $\\mathscr{I}^\\pm$, spacelike")
    views.append(v)

    v = View("anti_de_sitter", "anti-de Sitter space", box)
    frame(v, times=False)
    v.fill("cover", covers["anti_de_sitter"])
    for r in (0.5, 1, 2):
        X = float(np.arctan(r))
        v.line("r", [[[X, T0], [X, T1]]])
    v.line("boundary", [[[HALF, T0], [HALF, T1]]])
    v.label_xt([HALF, 2.3], "$\\mathscr{I}$", "l", dx=6)
    v.legend("cover", "anti-de Sitter space, conformal to the half $\\chi < \\pi/2$")
    v.legend("r", "anti-de Sitter's $r$ constant, at $L/2$, $L$ and $2L$, with $L = R$")
    v.legend("boundary", "the conformal boundary of anti-de Sitter space, timelike")
    v.legend("centre", "the pole $\\chi = 0$ and its antipode $\\chi = \\pi$")
    views.append(v)
    return views


# ---------------------------------------------------------------- the table

# ---------------------------------------------------------------- the C-metric

C_M = {"m": 1, "alpha": "1/6", "C": "3/4"}
C_KAPPA = 1 / 9                 # alpha (1 - 2 alpha m), the surface gravity of the acceleration horizon


def c_rstar(r):
    """Griffiths, Krtous and Podolsky's tortoise coordinate, eq. (10), at m = 1 and alpha = 1/6:
    alpha r* = k_c ln|1 + alpha r| + k_a ln|1 - alpha r| + k_o ln|1 - r/2m| with k_c = 3/8,
    k_a = -3/4 and k_o = 3/8, so that r*(0) = 0."""
    r = np.asarray(r, dtype=float)
    with np.errstate(divide="ignore"):
        return 6 * (3 / 8 * np.log(np.abs(1 + r / 6)) - 3 / 4 * np.log(np.abs(1 - r / 6))
                    + 3 / 8 * np.log(np.abs(1 - r / 2)))


def c_rstar_y(y):
    """The same r* in y = 1/(alpha r), which runs on through y = 0, r = infinity, to null
    infinity at y = -1 on the inner axis: alpha r* = k_c ln|1 + y| + k_a ln|1 - y|
    + k_o ln|1 - 2 alpha m y| - k_o ln(2 alpha m)."""
    y = np.asarray(y, dtype=float)
    with np.errstate(divide="ignore"):
        return 6 * (3 / 8 * np.log(np.abs(1 + y)) - 3 / 4 * np.log(np.abs(1 - y))
                    + 3 / 8 * np.log(np.abs(1 - y / 3)) + 3 / 8 * np.log(3))


def c_cell(cell, t, rs):
    """(p, q) of the points of a cell at the static time t and tortoise coordinate rs.

    With u = t - r* and v = t + r* and A(w) = arctan exp(k w), k the surface gravity of the
    acceleration horizon, the drawing is Kruskal's across the acceleration horizon, U = tan q
    and V = tan p near it, and continuous across the black hole horizon:

        II   the static region of the first black hole    p = -A(-v),  q = A(u)
        I    beyond the acceleration horizon, its future  p = A(-v),   q = A(u)
        BH   inside the first black hole                  p = -A(-v),  q = pi - A(u)

    II has the acceleration horizon on its left, v -> +infinity and u -> -infinity, and the
    black hole horizon on its right. The singularity r* = 0 of BH is v = u, where
    p + q = pi - A(-u) - A(u) = pi/2: the straight line T = pi/2. Null infinity on the inner
    axis, r* -> -infinity with u -> +infinity or v -> -infinity in I, is its two upper edges.
    The other cells are these reflected: X -> -X, the second black hole's side; T -> -T, the
    past; and X -> 2 pi - X, the exterior beyond the first black hole's Einstein-Rosen bridge.
    """
    t, rs = np.asarray(t, dtype=float), np.asarray(rs, dtype=float)
    u, v = t - rs, t + rs

    def A(w):
        return atan_exp(C_KAPPA * np.asarray(w, dtype=float))
    if cell == "II":
        return -A(-v), A(u)
    if cell == "I":
        return A(-v), A(u)
    if cell == "BH":
        return -A(-v), PI - A(u)
    raise KeyError(cell)


def c_images(p, q, where):
    """The reflections of c_cell's cells: 'L' X -> -X, 'P' T -> -T, 'F' X -> 2 pi - X, in turn."""
    p, q = np.asarray(p, dtype=float), np.asarray(q, dtype=float)
    for w in where:
        if w == "L":
            p, q = q, p
        elif w == "P":
            p, q = -q, -p
        elif w == "F":
            p, q = q - PI, p + PI
    return p, q


C_POLYGONS = {
    "II": [(0, 0), (-HALF, 0), (-HALF, HALF), (0, HALF)],
    "I": [(0, 0), (HALF, 0), (HALF, HALF), (0, HALF)],
    "BH": [(-HALF, HALF), (-HALF, PI), (0, HALF)],
}


def c_metric(ck, src):
    """The C-metric on the two halves of its axis, at alpha m = 1/6 and m = 1, as Jerry
    Griffiths, Pavel Krtous and Jiri Podolsky draw it in their figure 2 (b) and (c). On the axis
    sin(theta) = 0, so the plane of t and r is totally geodesic, and its metric is
    (-Q dt^2 + dr^2/Q)/(1 + alpha r cos(theta))^2 with Q = (1 - alpha^2 r^2)(1 - 2m/r): the null
    directions are Q's alone, dt = +-dr/Q = +-dr*, and the conformal factor places null
    infinity where 1 + alpha r cos(theta) vanishes. On the inner axis theta = 0 that is
    y = 1/(alpha r) = -1, past r = infinity, reached in the Hong-Teo plane of tau and y; on the
    outer axis theta = pi it is r = 1/alpha, the acceleration horizon itself, so there the left
    edges of the static region are null infinity. c_cell's docstring gives the maps.
    """
    views = []
    rinf = float(c_rstar_y(0.0))

    def planes(theta):
        sph = Plane(src, "c_metric", "spherical", ("t", "r"), {"theta": theta, "phi": "0"}, C_M)
        assert sph.g[0, 1] == 0
        return sph
    inner, outer = planes("0"), planes("pi")
    ht = Plane(src, "c_metric", "hong_teo", ("\\tau", "y"), {"x": "1", "phi": "0"}, C_M)
    assert ht.g[0, 1] == 0
    span = 40
    for name, pl in (("inner axis", inner), ("outer axis", outer)):
        ck.chart(f"C-metric {name}, static region", pl, lambda t, r: c_cell("II", t, c_rstar(r)),
                 ck.uniform(-span, span), ck.uniform(2.001, 5.999), lambda t, r: (1, 0))
        ck.chart(f"C-metric {name}, black hole", pl, lambda t, r: c_cell("BH", t, c_rstar(r)),
                 ck.uniform(-span, span), ck.uniform(0.01, 1.999), lambda t, r: (0, -1))
    ck.chart("C-metric inner axis, beyond the acceleration horizon", inner, lambda t, r: c_cell("I", t, c_rstar(r)),
             ck.uniform(-span, span), ck.uniform(6.001, 300), lambda t, r: (0, 1))

    def ht_map(cell):
        return lambda tau, y: c_cell(cell, 6 * np.asarray(tau, dtype=float), c_rstar_y(y))
    ck.chart("C-metric Hong-Teo inner axis, beyond the acceleration horizon", ht, ht_map("I"),
             ck.uniform(-span / 6, span / 6), ck.uniform(-0.999, 0.999), lambda tau, y: (0, -1))
    ck.chart("C-metric Hong-Teo inner axis, static region", ht, ht_map("II"),
             ck.uniform(-span / 6, span / 6), ck.uniform(1.001, 2.999), lambda tau, y: (1, 0))
    ck.chart("C-metric Hong-Teo inner axis, black hole", ht, ht_map("BH"),
             ck.uniform(-span / 6, span / 6), ck.uniform(3.001, 200), lambda tau, y: (0, 1))
    ck.limit("C-metric: r* in y is r* in r, y = 1/(alpha r)", c_rstar_y(6 / np.array([0.5, 1.5, 3, 5, 9, 40])),
             c_rstar(np.array([0.5, 1.5, 3, 5, 9, 40])), 1e-9)
    p, q = c_cell("BH", np.array([-9.0, 0.0, 9.0]), c_rstar(np.full(3, 1e-12)))
    ck.limit("C-metric: r -> 0 in the black hole lands on T = pi/2", p + q, [HALF] * 3, 1e-9)
    p, q = c_cell("II", np.array([3.0]), c_rstar(np.array([6 - 1e-12])))
    ck.limit("C-metric: r -> 1/alpha at fixed t lands on the acceleration horizon's bifurcation, (0, 0)",
             point(p[0], q[0]), [0, 0], 1e-4)
    # r* -> -infinity as r -> 2m only as 2.25 ln(r - 2m), too slowly to follow in r itself.
    ck.limit("C-metric: r* -> -infinity as r -> 2m", float(c_rstar(np.array([2 + 1e-12]))[0] < -60), 1, 0.5)
    p, q = c_cell("II", np.array([3.0]), np.array([-1e4]))
    ck.limit("C-metric: r* -> -infinity at fixed t lands on the black hole's bifurcation, (pi, 0)",
             point(p[0], q[0]), [PI, 0], 1e-4)
    p, q = c_cell("I", np.array([1.0]), c_rstar_y(np.array([-1 + 1e-15])))
    ck.limit("C-metric inner axis: y -> -1 at fixed t lands on the top of the region beyond the "
             "acceleration horizon, (0, pi), where its two edges of null infinity meet", point(p[0], q[0]), [0, PI], 2e-3)
    uu = np.array([-20.0, 0.0, 20.0])
    ys = np.full(3, -1 + 1e-15)
    p, q = c_cell("I", uu + c_rstar_y(ys), c_rstar_y(ys))
    ck.limit("C-metric inner axis: y -> -1 along a line of constant u lands on the edge p = pi/2", p, [HALF] * 3, 1e-6)
    ck.diverges("C-metric: the Kretschmann scalar diverges at r = 0 on the inner axis",
                inner.kretschmann(0, 1e-2), inner.kretschmann(0, 1e-3))
    ck.diverges("C-metric: the Kretschmann scalar diverges at r = 0 on the outer axis",
                outer.kretschmann(0, 1e-2), outer.kretschmann(0, 1e-3))
    ck.finite("C-metric: the Kretschmann scalar is finite at both horizons on both halves of the axis",
              np.concatenate([pl.kretschmann(np.zeros(2), np.array([2.0, 6.0])) for pl in (inner, outer)]))
    ck.finite("C-metric: the Kretschmann scalar is finite at null infinity on the inner axis",
              ht.kretschmann(np.zeros(3), np.array([-0.999, -0.9999, -1.0])))
    ck.finite("C-metric: the Kretschmann scalar is finite at null infinity on the outer axis",
              outer.kretschmann(np.zeros(2), np.array([5.999, 6.0])))

    times = [c / C_KAPPA for c in (-1.6, -0.6, 0, 0.6, 1.6)]
    r_II, r_BH, r_I = (2.3, 3, 4, 5.3), (0.8, 1.5), (7, 10, 20)
    y_I = (0.6, 0.2, -0.2, -0.6)
    rr_II, rr_BH = spread(2, 6, 600, 16), spread(0, 2, 600, 16)
    tt = spread(-np.inf, np.inf, 500, 10)

    def draw_cell(v, cell, where, radii, rmap, rr_range):
        for r in radii:
            v.curve("r", *c_images(*c_cell(cell, tt, np.full_like(tt, rmap(r))), where))
        for t0 in times:
            v.curve("t", *c_images(*c_cell(cell, np.full_like(rr_range, t0), rmap(rr_range)), where))

    def region(v, cell, where):
        pts = [c_images(np.array([a]), np.array([b]), where) for a, b in C_POLYGONS[cell]]
        v.fill("region", [point(p[0], q[0]) for p, q in pts])

    def seg(v, cls, a, b, where, zig=False):
        pa = c_images(np.array([a[0]]), np.array([a[1]]), where)
        pb = c_images(np.array([b[0]]), np.array([b[1]]), where)
        v.segment(cls, (pa[0][0], pa[1][0]), (pb[0][0], pb[1][0]), zig)

    def statics(v, sides, scri_left):
        """The static regions and the black holes beside them; `sides` lists the reflections,
        and on the outer axis the left edges of the static region are null infinity."""
        for where in sides:
            # The black hole of an exterior reflected by F is the black hole it shares with the
            # exterior it is reflected from, so it is drawn once.
            hole = "F" not in where
            for past in ("", "P"):
                region(v, "II", where + past)
                draw_cell(v, "II", where + past, r_II, c_rstar, rr_II)
                if hole:
                    region(v, "BH", where + past)
                    draw_cell(v, "BH", where + past, r_BH, c_rstar, rr_BH)
                    seg(v, "singular", (-HALF, PI), (0, HALF), where + past, zig=True)
                seg(v, "horizon", (-HALF, 0), (-HALF, HALF), where + past)
                seg(v, "horizon", (-HALF, HALF), (0, HALF), where + past)
                seg(v, "scri" if scri_left else "horizon", (0, 0), (0, HALF), where + past)

    box_in = [-2 * PI - 0.3, 2 * PI + 0.3, -PI - 0.25, PI + 0.25]
    box_out = [-0.3, 2 * PI + 0.3, -HALF - 0.45, HALF + 0.45]
    equator = slices.moments("c_metric", "equator")[0]
    horizon = slices.moments("c_metric", "horizon", label="$t = 0$, $r = 2m$")[0]
    lo, hi = equator.reach("spherical", "r")
    rr = np.linspace(lo, hi, 2)

    for system in ("spherical", "hong_teo"):
        v = View(f"inner_{system}", "the inner axis", box_in, system)
        statics(v, ["", "F", "L", "FL"], False)
        for where in ("", "P"):
            region(v, "I", where)
            if system == "spherical":
                draw_cell(v, "I", where, r_I, c_rstar, spread(6, np.inf, 600, 16))
            else:
                draw_cell(v, "I", where, y_I, c_rstar_y, spread(-1, 1, 600, 16)[::-1])
            seg(v, "scri", (0, HALF), (HALF, HALF), where)
            seg(v, "scri", (HALF, 0), (HALF, HALF), where)
        v.fill("cover", [point(a, b) for a, b in C_POLYGONS["II"]])
        v.fill("cover", [point(a, b) for a, b in C_POLYGONS["BH"]])
        if system == "spherical":
            ts = spread(-np.inf, np.inf, 500, 10)
            p, q = c_cell("I", ts, np.full_like(ts, rinf))
            v.fill("cover", [[0.0, 0.0]] + [point(a, b) for a, b in zip(p, q)])
            v.curve("chartedge", p, q)
            v.legend("chartedge", "$r = \\infty$, where the coordinate $r$ ends and $y = 1/(\\alpha r)$ passes through 0")
        else:
            v.fill("cover", [point(a, b) for a, b in C_POLYGONS["I"]])
        for sx in (1, -1):
            v.label_xt([sx * Q4 * 0.9, 3 * Q4 + 0.05], "$\\mathscr{I}^+$", "bl" if sx > 0 else "br", dx=4 * sx, dy=-3)
            v.label_xt([sx * Q4 * 0.9, -3 * Q4 - 0.05], "$\\mathscr{I}^-$", "tl" if sx > 0 else "tr", dx=4 * sx, dy=3)
            v.label_xt([sx * (PI + 0.15), 1.15], "black hole", cls="region")
            v.label_xt([sx * (PI + 0.15), -1.15], "white hole", cls="region")
        # Each horizon's name stands on the outer side of its point, clear of the other's.
        v.label_xt([Q4, Q4], "$r = 1/\\alpha$", "tr", "small", dx=-5, dy=1)
        v.label_xt([PI - Q4, Q4], "$r = 2m$", "tl", "small", dx=5, dy=1)
        v.slice(equator, [c_images(*c_cell("II", 0 * rr, c_rstar(rr)), w) for w in ("", "L")])
        v.slice(horizon, points=[(-HALF, HALF)])
        v.legend("cover", "the regions that " + ("$t$ and $r > 0$" if system == "spherical" else "$\\tau$ and $y > -1$")
                 + " cover")
        v.legend("r", "$r$ constant, in units of $m$" if system == "spherical" else "$y = 1/(\\alpha r)$ constant")
        v.legend("t", "$ct$ constant" if system == "spherical" else "$\\tau$ constant")
        v.legend("horizon", "the black hole horizons $r = 2m$ and the acceleration horizons $r = 1/\\alpha$")
        v.legend("singular", "$r = 0$, where the Kretschmann scalar diverges")
        v.legend("scri", "null infinity, where $1 + \\alpha r\\cos\\theta$ vanishes")
        v.set(restriction="The half axis $\\theta = 0$ between the black holes only, a totally geodesic surface, "
                          "each point in the diagram a single event.")
        views.append(v)

        v = View(f"outer_{system}", "the outer axis", box_out, system)
        statics(v, ["", "F"], True)
        v.fill("cover", [point(a, b) for a, b in C_POLYGONS["II"]])
        v.fill("cover", [point(a, b) for a, b in C_POLYGONS["BH"]])
        v.label_xt([PI, 1.15], "black hole", cls="region")
        v.label_xt([PI, -1.15], "white hole", cls="region")
        v.label_xt([PI - Q4, Q4], "$r = 2m$", "tr", "small", dx=-5, dy=1)
        for at in ((0, 0), (2 * PI, 0)):
            v.layers.append({"kind": "point", "class": "infinity", "at": rounded(at)})
        for x0, sx in ((Q4, -1), (2 * PI - Q4, 1)):
            v.label_xt([x0, Q4], "$\\mathscr{I}^+$", "bl" if sx > 0 else "br", dx=4 * sx, dy=-3)
            v.label_xt([x0, -Q4], "$\\mathscr{I}^-$", "tl" if sx > 0 else "tr", dx=4 * sx, dy=3)
        v.slice(equator, [c_cell("II", 0 * rr, c_rstar(rr))])
        v.slice(horizon, points=[(-HALF, HALF)])
        v.legend("cover", "the static region and the black hole that " + ("$t$ and $r$" if system == "spherical" else "$\\tau$ and $y$") + " cover")
        v.legend("r", "$r$ constant, in units of $m$" if system == "spherical" else "$y = 1/(\\alpha r)$ constant")
        v.legend("t", "$ct$ constant" if system == "spherical" else "$\\tau$ constant")
        v.legend("horizon", "the black hole horizon $r = 2m$")
        v.legend("singular", "$r = 0$, where the Kretschmann scalar diverges")
        v.legend("scri", "null infinity, $r = 1/\\alpha$, where the acceleration horizon lies on this half of the axis")
        v.set(restriction="The half axis $\\theta = \\pi$ beyond the black hole only, a totally geodesic surface, "
                          "each point in the diagram a single event.")
        views.append(v)
    settings = "$m = 1$, $\\alpha = 1/(6m)$, so that the horizons are at $r = 2m$ and $r = 1/\\alpha = 6m$."
    for view in views:
        view.set(settings=settings)
    return views


def milne(ck, src):
    """The Milne universe as the wedge of Minkowski's triangle above the future light cone of the
    event T = R = 0. With T = t cosh chi and R = ct sinh chi, the null coordinates of Minkowski's
    own diagram are u = cT - R = ct e^-chi and v = cT + R = ct e^chi, so p = arctan(u/l) and
    q = arctan(v/l) draw the comoving hyperbolic chart with Minkowski's map for any length l,
    drawn as 1: u > 0 and v > 0 are p > 0 and q > 0, the wedge under the centre and I+ with
    corners (X, T) = (0, 0), (pi/2, pi/2) and (0, pi). The comoving spherical chart enters through
    chi = arcsinh r, the logarithmic time through t = t_0 e^(tau/t_0) at t_0 = 1, and the inertial
    chart is Minkowski's spherical chart inside R < cT. Every point of the line t = 0 goes to the
    event p = q = 0, and chi -> infinity at fixed t goes to I+, which the limits check."""
    views = []
    tri_box = [-0.35, PI + 0.35, -PI - 0.25, PI + 0.25]
    wedge = [[0, 0], [HALF, HALF], [0, PI]]

    def hyper(t, chi):
        t, chi = np.asarray(t, dtype=float), np.asarray(chi, dtype=float)
        return np.arctan(t * np.exp(-chi)), np.arctan(t * np.exp(chi))

    def spherical(t, r):
        return hyper(t, np.arcsinh(np.asarray(r, dtype=float)))

    def logarithmic(tau, chi):
        return hyper(np.exp(np.asarray(tau, dtype=float)), chi)

    def inertial(T, R):
        return mink_pq(T, R)

    hp = Plane(src, "milne", "comoving_hyperbolic", ("t", "\\chi"), EQUATOR)
    sp_ = Plane(src, "milne", "comoving_spherical", ("t", "r"), EQUATOR)
    lp = Plane(src, "milne", "logarithmic_time", ("\\tau", "\\chi"), EQUATOR, {"t_0": 1})
    ip = Plane(src, "milne", "inertial", ("T", "R"), EQUATOR)
    ck.chart("Milne comoving hyperbolic", hp, hyper, ck.uniform(0.01, 20), ck.uniform(0.001, 6), lambda t, c: (1, 0))
    ck.chart("Milne comoving spherical", sp_, spherical, ck.uniform(0.01, 20), ck.uniform(0.001, 50), lambda t, r: (1, 0))
    ck.chart("Milne logarithmic time", lp, logarithmic, ck.uniform(-4, 3), ck.uniform(0.001, 6), lambda tau, c: (1, 0))
    T = ck.uniform(0.01, 20)
    ck.chart("Milne inertial", ip, inertial, T, T * ck.uniform(0.001, 0.999), lambda T, R: (1, 0))
    chis = np.linspace(0, 8, 41)
    ck.limit("Milne: t -> 0 at every chi is the event p = q = 0", np.concatenate(hyper(np.full(41, 1e-15), chis)),
             np.zeros(82), 1e-9)
    p, q = hyper(np.linspace(0.1, 10, 41), np.full(41, 40.0))
    ck.limit("Milne: chi -> infinity at fixed t is I+, q = pi/2", q, np.full(41, HALF), 1e-9)
    p, q = inertial(np.linspace(0.1, 10, 41), np.linspace(0.1, 10, 41))
    ck.limit("Milne: the light cone R = cT is the edge p = 0 of the wedge", p, np.zeros(41), 1e-12)
    ck.finite("Milne: the centre chi = 0 is regular", hp.kretschmann(ck.uniform(0.1, 5, 50), np.full(50, 1e-6)))
    ck.finite("Milne: the curvature vanishes toward t = 0", hp.kretschmann(np.full(50, 1e-6), ck.uniform(0, 5, 50)))

    moments = slices.moments("milne")

    def frame(v, legend):
        v.fill("region", TRIANGLE)
        v.fill("cover", wedge)
        triangle_edges(v, centre="$R = 0$")
        v.line("chartedge", [[[0, 0], [HALF, HALF]]])
        v.point("mark", (0, 0))
        v.label_xt([Q4, Q4], "$cT = R$", "tl", "small", dx=6, dy=2)
        v.legend("cover", legend)
        v.legend("chartedge", "the future light cone $cT = R$ of the event $T = R = 0$, where the Milne universe ends")
        v.legend("mark", "the event $T = R = 0$, from which every comoving particle moves off")
        v.legend("centre", "$R = 0$, the world line of the comoving particle at $\\chi = 0$")
        for m in moments:
            hi = m.reach("comoving_hyperbolic", "\\chi")[1]
            chi = np.linspace(0, hi, 400)
            v.slice(m, [hyper(np.full_like(chi, m.time), chi)])

    S_T = spread(0, np.inf, 500, 12)
    S_CHI = np.linspace(0, 40, 2000)
    TIMES, CHIS = (0.25, 0.5, 1, 2, 4), (0.5, 1, 1.5, 2, 3)

    v = View("comoving_hyperbolic", "Comoving hyperbolic", tri_box, "comoving_hyperbolic")
    grid(v, "r", lambda c, t: hyper(t, c), CHIS, S_T)
    grid(v, "t", hyper, TIMES, S_CHI)
    frame(v, "the Milne universe, which $t$ and $\\chi$ cover")
    v.legend("r", "$\\chi$ constant, a comoving particle, at $1/2$, $1$, $3/2$, $2$, and $3$")
    v.legend("t", "$ct$ constant, at $1/4$, $1/2$, $1$, $2$, and $4$ in units of $\\ell$")
    views.append(v)

    v = View("comoving_spherical", "Comoving spherical", tri_box, "comoving_spherical")
    RS = tuple(float(np.sinh(c)) for c in CHIS)
    grid(v, "r", lambda r, t: spherical(t, r), RS, S_T)
    grid(v, "t", spherical, TIMES, np.sinh(S_CHI))
    frame(v, "the Milne universe, which $t$ and $r$ cover")
    v.legend("r", "$r = \\sinh\\chi$ constant, a comoving particle, at $\\chi = 1/2$, $1$, $3/2$, $2$, and $3$")
    v.legend("t", "$ct$ constant, at $1/4$, $1/2$, $1$, $2$, and $4$ in units of $\\ell$")
    views.append(v)

    v = View("logarithmic_time", "Logarithmic time", tri_box, "logarithmic_time")
    grid(v, "r", lambda c, tau: logarithmic(tau, c), CHIS, np.linspace(-30, 12, 2000))
    grid(v, "t", logarithmic, (-2, -1, 0, 1, 2), S_CHI)
    frame(v, "the Milne universe, which $\\tau$ and $\\chi$ cover")
    v.legend("r", "$\\chi$ constant, a comoving particle, at $1/2$, $1$, $3/2$, $2$, and $3$")
    v.legend("t", "$\\tau$ constant, at $-2t_0$, $-t_0$, $0$, $t_0$, and $2t_0$, with $ct_0 = \\ell$")
    views.append(v)

    v = View("inertial", "Inertial", tri_box, "inertial")
    for R in (0.5, 1, 2, 4):
        T = np.concatenate([[R], R + spread(0, np.inf, 500, 12)[1:]])
        v.curve("r", *inertial(T, np.full_like(T, R)))
    for T0 in (0.5, 1, 2, 4):
        R = np.linspace(0, T0, 400)
        v.curve("t", *inertial(np.full_like(R, T0), R))
    frame(v, "the Milne universe, $R < cT$, which $T$ and $R$ cover")
    v.legend("r", "$R$ constant, at $\\ell/2$, $\\ell$, $2\\ell$, and $4\\ell$")
    v.legend("t", "$cT$ constant, at $\\ell/2$, $\\ell$, $2\\ell$, and $4\\ell$")
    views.append(v)
    for v in views:
        v.set(settings="$\\ell$, any length, the scale of $p = \\arctan((cT - R)/\\ell)$ and $q = \\arctan((cT + R)/\\ell)$.")
    return views


def nariai(ck, src):
    """dS2 x S2, both radii a = 1/sqrt(Lambda), drawn at a = 1. The de Sitter factor is the
    hyperboloid -Z0^2 + Z1^2 + Z2^2 = a^2; with Z0 = a tan(eta), Z1 = a sin(chi)/cos(eta) and
    Z2 = -a cos(chi)/cos(eta) its metric is a^2(-deta^2 + dchi^2)/cos^2(eta), conformal to the strip
    |eta| < pi/2 with chi running round a circle, so p, q = (eta -+ chi)/2 draw it with X = chi and
    T = eta, the edges X = 0 and 2 pi one line. The conformal chart is that map as it stands, the
    global chart enters through tan(eta) = sinh(ct/a), and the static chart, with
    Z0 = sqrt(a^2 - r^2) sinh(ct/a), Z1 = sqrt(a^2 - r^2) cosh(ct/a) and Z2 = r, through
    tan(eta) = sqrt(1 - r^2) sinh(ct) and chi = atan2(sqrt(1 - r^2) cosh(ct), -r) at a = 1, which
    covers the diamond about (chi, eta) = (pi/2, 0) whose edges are its horizons r = -+a. Every point
    is a 2-sphere of radius a."""
    st = Plane(src, "nariai", "static", ("t", "r"), EQUATOR, {"Lambda": 1})
    gl = Plane(src, "nariai", "global", ("t", "\\chi"), EQUATOR, {"Lambda": 1})
    co = Plane(src, "nariai", "conformal", ("\\eta", "\\chi"), EQUATOR, {"Lambda": 1})

    def conformal_pq(eta, chi):
        eta, chi = np.asarray(eta, dtype=float), np.asarray(chi, dtype=float)
        return (eta - chi) / 2, (eta + chi) / 2

    def global_pq(t, chi):
        return conformal_pq(np.arctan(np.sinh(np.asarray(t, dtype=float))), chi)

    def static_pq(t, r):
        t, r = np.asarray(t, dtype=float), np.asarray(r, dtype=float)
        s = np.sqrt(1 - r * r)
        return conformal_pq(np.arctan(s * np.sinh(t)), np.arctan2(s * np.cosh(t), -r))
    ck.chart("Nariai conformal", co, conformal_pq, ck.uniform(-HALF + 1e-3, HALF - 1e-3), ck.uniform(0, 2 * PI),
             lambda e, c: (1, 0))
    ck.chart("Nariai global", gl, global_pq, ck.uniform(-8, 8), ck.uniform(0, 2 * PI), lambda t, c: (1, 0))
    ck.chart("Nariai static", st, static_pq, ck.uniform(-8, 8), ck.uniform(-0.999, 0.999), lambda t, r: (1, 0))

    t, r = ck.uniform(-4, 4, 2000), ck.uniform(-0.999, 0.999, 2000)
    p, q = static_pq(t, r)
    eta, chi = p + q, q - p
    scale = 1 + np.cosh(t)
    ck.limit("Nariai: the static chart lands where the hyperboloid puts it",
             np.concatenate([(np.tan(eta) - np.sqrt(1 - r * r) * np.sinh(t)) / scale,
                             (np.sin(chi) / np.cos(eta) - np.sqrt(1 - r * r) * np.cosh(t)) / scale,
                             (-np.cos(chi) / np.cos(eta) - r) / scale]), 0, 1e-10)
    p, q = static_pq(np.array([0.0, 2.0, -2.0]), np.array([1 - 1e-14, -1 + 1e-14, 1 - 1e-14]))
    ck.limit("Nariai: r -> +-1/sqrt(Lambda) at fixed t lands on the bifurcation points chi = pi and 0, eta = 0",
             np.concatenate([q - p, p + q]), [PI, 0, PI, 0, 0, 0], 1e-6)
    p, q = static_pq(np.array([60.0, -60.0]), np.array([0.3, 0.3]))
    ck.limit("Nariai: t -> +-infinity in the static chart lands on the corners (pi/2, +-pi/2)",
             np.concatenate([q - p, p + q]), [HALF, HALF, HALF, -HALF], 1e-6)
    p, q = global_pq(np.array([60.0, -60.0]), np.array([1.0, 1.0]))
    ck.limit("Nariai: t -> +-infinity in the global chart lands on eta = +-pi/2", p + q, [HALF, -HALF], 1e-6)
    ck.finite("Nariai: the curvature is the same everywhere",
              gl.kretschmann(ck.uniform(-5, 5, 50), ck.uniform(0, 2 * PI, 50)))

    box = [-0.35, 2 * PI + 0.35, -HALF - 0.3, HALF + 0.3]
    strip_ = [[0, -HALF], [2 * PI, -HALF], [2 * PI, HALF], [0, HALF]]
    diamond = [[0, 0], [HALF, HALF], [PI, 0], [HALF, -HALF]]
    universe = slices.moments("nariai", "universe")
    sphere = slices.moments("nariai", "sphere")[0]

    def frame(v):
        v.fill("region", strip_)
        v.line("scri", [[[0, HALF], [2 * PI, HALF]], [[0, -HALF], [2 * PI, -HALF]]])
        # The horizons of the observer at chi = pi/2 and of the antipode at 3 pi/2.
        v.line("horizon", [[[0, 0], [HALF, HALF]], [[HALF, HALF], [PI, 0]], [[PI, 0], [HALF, -HALF]],
                           [[HALF, -HALF], [0, 0]], [[PI, 0], [3 * HALF, HALF]], [[3 * HALF, HALF], [2 * PI, 0]],
                           [[2 * PI, 0], [3 * HALF, -HALF]], [[3 * HALF, -HALF], [PI, 0]]])
        v.label_xt([PI, HALF], "$\\mathscr{I}^+$", "b", dy=-5)
        v.label_xt([PI, -HALF], "$\\mathscr{I}^-$", "t", dy=5)
        v.label_xt([0, -1.25], "$\\chi = 0$", "r", "coord", dx=-6)
        v.label_xt([2 * PI, -1.25], "$\\chi = 2\\pi$", "l", "coord", dx=6)
        v.label_xt([3 * HALF, -0.28], "the antipode's static patch", cls="region")
        v.label_xt([PI, 1.05], "expanding", cls="region")
        v.label_xt([PI, -1.05], "contracting", cls="region")
        v.legend("scri", "future and past infinity $\\mathscr{I}^\\pm$, spacelike")
        v.legend("horizon", "the horizons $r = \\pm 1/\\sqrt{\\Lambda}$ of two antipodal static patches")
        for m in universe:
            eta = float(np.arctan(np.sinh(m.time)))
            v.slice(m, xt=[[[0.0, eta], [2 * PI, eta]]])
        v.slice(sphere, points=[conformal_pq(0.0, HALF)], label="$t = 0$, $\\chi = \\pi/2$")

    views = []
    v = View("static", "Static", box, "static")
    frame(v)
    v.fill("cover", diamond)
    grid(v, "r", lambda r, t: static_pq(t, r), (-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9), S_ALL)
    grid(v, "t", static_pq, (-2, -1, 0, 1, 2), np.tanh(np.linspace(-14, 14, 1001)))
    # Below its point, which leaves the room above it to the moment of the embedding.
    label_on(v, static_pq(0, 0), "$r = 0$", "t", dy=3)
    v.legend("cover", "the static patch, which $t$ and $r$ cover")
    v.legend("r", "$r$ constant, from $-0.9$ to $0.9$ in units of $1/\\sqrt{\\Lambda}$")
    v.legend("t", "$ct$ constant, every $1/\\sqrt{\\Lambda}$")
    views.append(v)

    for vid, label, system, time, name, values in (
            ("global", "Global", "global", "$t$", "$ct$ constant, every $1/\\sqrt{\\Lambda}$", (-2, -1, 0, 1, 2)),
            ("conformal", "Conformal", "conformal", "$\\eta$", "$\\eta$ constant, every $\\pi/8$",
             (-3 * PI / 8, -PI / 4, -PI / 8, 0, PI / 8, PI / 4, 3 * PI / 8))):
        v = View(vid, label, box, system)
        frame(v)
        v.fill("cover", strip_)
        for chi in (HALF, PI, 3 * HALF):
            v.line("r", [[[chi, -HALF], [chi, HALF]]])
        for c in values:
            eta = float(np.arctan(np.sinh(c))) if vid == "global" else c
            v.line("t", [[[0, eta], [2 * PI, eta]]])
        v.legend("cover", f"the whole spacetime, which {time} and $\\chi$ cover")
        v.legend("r", "$\\chi$ constant, at $\\pi/2$, $\\pi$ and $3\\pi/2$")
        v.legend("t", name)
        views.append(v)
    return views


def kp_pq(x):
    """Khan and Penrose's u or v into the drawing: itself where it is positive, where a wave has
    passed, and its arctangent where it is negative, ahead of the wave. The map is continuous
    with its first derivative across the front, and sends the past edge to -pi/2."""
    x = np.asarray(x, dtype=float)
    return np.where(x >= 0, x, np.arctan(x))


class KPFlatRegion:
    """The published metric on the plane of u and v with u, or v, or both, set to zero, which
    is the metric ahead of that wave, as the chart's convention states."""

    def __init__(self, plane, zero_u, zero_v):
        self.plane, self.zero_u, self.zero_v = plane, zero_u, zero_v

    def metric(self, u, v, fvals=None):
        u, v = np.asarray(u, dtype=float), np.asarray(v, dtype=float)
        return self.plane.metric(0 * u if self.zero_u else u, 0 * v if self.zero_v else v)


def khan_penrose(ck, src):
    """The plane x = y = 0, totally geodesic, since x -> -x and y -> -y are isometries that fix
    it, drawn in p, q = kp_pq(u), kp_pq(v). Where both waves have passed, the published double
    null and cosmological charts are checked, the second through p, q = sin((tau +- sigma)/2),
    and the published Kretschmann scalar is checked to diverge on u^2 + v^2 = 1. Ahead of one
    wave or both the metric is the published one with u or v set to zero, checked through the
    same maps, and its Riemann tensor, computed from those components, is checked to vanish:
    each of those regions is flat, up to the fold singularities u = 1 and v = 1, where the
    metric degenerates with the curvature still zero."""
    params = {"L": 1}
    fixed = {"x": "0", "y": "0"}
    null = Plane(src, "khan_penrose", "double_null", ("u", "v"), fixed, params)
    r, th = np.sqrt(ck.uniform(0, 0.999)), ck.uniform(0, HALF)
    u, w = r * np.cos(th), r * np.sin(th)
    ck.chart("Khan-Penrose double null, both waves passed", null, lambda a, b: (kp_pq(a), kp_pq(b)), u, w,
             lambda a, b: (1, 1))
    cos = Plane(src, "khan_penrose", "cosmological", ("\\tau", "\\sigma"), fixed, params)
    tau = ck.uniform(0.001, HALF - 0.001)
    sigma = tau * ck.uniform(-0.999, 0.999)
    ck.chart("Khan-Penrose cosmological", cos, lambda t, s: (np.sin((t + s) / 2), np.sin((t - s) / 2)), tau, sigma,
             lambda t, s: (1, 0))
    for name, zu, zv, a, b in (("ahead of both waves", True, True, ck.uniform(-30, 0), ck.uniform(-30, 0)),
                               ("behind the wave on u = 0 alone", False, True, ck.uniform(0, 0.999), ck.uniform(-30, 0)),
                               ("behind the wave on v = 0 alone", True, False, ck.uniform(-30, 0), ck.uniform(0, 0.999))):
        ck.chart(f"Khan-Penrose {name}", KPFlatRegion(null, zu, zv), lambda a, b: (kp_pq(a), kp_pq(b)), a, b,
                 lambda a, b: (1, 1))
    # The Riemann tensor of the published metric with v, or u, set to zero, from its components.
    _, entry, reader = nr.load("khan_penrose", "double_null")
    g = nr.published_matrix(reader, entry, "metric_components").subs({reader.parameters["L"]: 1})
    symbols = [reader.symbol[c] for c in entry["coords"]]
    for name, zero in (("u = 0", symbols[0]), ("v = 0", symbols[1])):
        geo = vm.Geometry(g.subs(zero, 0), symbols, 10 ** 6)
        riemann = geo.riemann_llll()
        nonzero = sum(1 for i in vm._indices(4, 4) if vm.norm(vm._at(riemann, i)) != 0)
        ck.limit(f"Khan-Penrose: flat behind one wave, the published metric at {name}", nonzero, 0, tol=0.5)
    th = ck.uniform(0.01, HALF - 0.01, 200)
    d = 1e-4
    near, nearer = ((1 - d) * np.cos(th), (1 - d) * np.sin(th)), ((1 - d / 10) * np.cos(th), (1 - d / 10) * np.sin(th))
    ck.diverges("Khan-Penrose: u^2 + v^2 = 1 is a curvature singularity", null.kretschmann(*near),
                null.kretschmann(*nearer))
    edge_sigma = ck.uniform(-HALF, HALF, 200)
    ck.limit("Khan-Penrose: tau = pi/2 is u^2 + v^2 = 1",
             np.sin((HALF + edge_sigma) / 2) ** 2 + np.sin((HALF - edge_sigma) / 2) ** 2, 1.0, tol=1e-12)

    edge = 1 + HALF
    box = [-edge, edge, -PI, math.sqrt(2)]
    arc = np.linspace(0, HALF, 181)
    s_neg = spread(-np.inf, 0, 400, 9)[:-1]
    whole = ([point(p, -HALF) for p in (-HALF, 1)] + [point(1, 0)]
             + [point(np.cos(a), np.sin(a)) for a in arc] + [point(-HALF, 1)])
    region_iv = [point(0, 0)] + [point(np.cos(a), np.sin(a)) for a in arc]
    restriction = "The plane $x = y = 0$ only, totally geodesic, each point in the diagram a single event."
    moments = slices.moments("khan_penrose")
    views = []
    for system in ("double_null", "cosmological"):
        v = View(system, {"double_null": "Double Null", "cosmological": "Cosmological"}[system], box, system)
        v.fill("region", whole)
        v.fill("cover", region_iv)
        if system == "double_null":
            for c in (0.2, 0.4, 0.6, 0.8):
                top = math.sqrt(1 - c * c)
                v.curve("null", np.full(2, c), np.array([0, top]))
                v.curve("null", np.array([0, top]), np.full(2, c))
            v.legend("cover", "the region where both waves have passed, which $u$ and $v$ cover")
            v.legend("null", "$u$ constant and $v$ constant, every one a light ray")
        else:
            for c in (0.3, 0.6, 0.9, 1.2):
                s = np.linspace(-c, c, 200)
                v.curve("t", np.sin((c + s) / 2), np.sin((c - s) / 2))
            for c in (-0.9, -0.45, 0.45, 0.9):
                t = np.linspace(abs(c), HALF, 200)
                v.curve("r", np.sin((t + c) / 2), np.sin((t - c) / 2))
            v.legend("cover", "the region where both waves have passed, which $\\tau$ and $\\sigma$ cover")
            v.legend("t", "$\\tau$ constant, spacelike")
            v.legend("r", "$\\sigma$ constant")
        v.curve("surface", np.zeros_like(s_neg), kp_pq(s_neg))
        v.curve("surface", kp_pq(s_neg), np.zeros_like(s_neg))
        v.segment("surface", (0, 0), (0, 1))
        v.segment("surface", (0, 0), (1, 0))
        v.segment("singular", (1, -HALF), (1, 0), zig=True)
        v.segment("singular", (-HALF, 1), (0, 1), zig=True)
        v.curve("singular", np.cos(arc), np.sin(arc), zig=True)
        v.segment("scri", (-HALF, -HALF), (1, -HALF))
        v.segment("scri", (-HALF, -HALF), (-HALF, 1))
        v.layers.append({"kind": "point", "class": "infinity", "at": rounded(point(-HALF, -HALF))})
        v.point("mark", (0, 0))
        v.label_xt(point(-HALF, -HALF), "$i^-$", "b", dy=-6)
        v.label((1, -HALF / 2), "$u = 1$", "l", "small", dx=6)
        v.label((-HALF / 2, 1), "$v = 1$", "r", "small", dx=-6)
        v.label((math.sqrt(0.5), math.sqrt(0.5)), "$u^2 + v^2 = 1$", "t", "small", dy=6)
        v.label((0.25 * -HALF, -HALF), "$\\mathscr{I}^-$", "bl", dx=5, dy=-3)
        v.label((-HALF, 0.25 * -HALF), "$\\mathscr{I}^-$", "br", dx=-5, dy=-3)
        v.legend("surface", "the fronts of the two waves, $u = 0$ and $v = 0$, where the Riemann tensor has a "
                            "delta singularity")
        v.legend("mark", "the collision, $u = v = 0$")
        v.legend("singular", "$u^2 + v^2 = 1$, a spacelike curvature singularity, where the Kretschmann scalar "
                             "diverges")
        v.legend("singular", "$u = 1$ and $v = 1$ behind one wave alone, fold singularities, where the curvature "
                             "is zero")
        v.legend("scri", "past null infinity $\\mathscr{I}^-$")
        v.set(restriction=restriction)
        for m in moments:
            c = math.sin(m.time / 2)
            v.slice(m, points=[(c, c)])
        views.append(v)
    return views


# ---------------------------------------------------------------- Aichelburg-Sexl

AS_JUMP = float(np.log(8.0))    # the jump of a ray moving left at rho = rho_0/8, in units of 8GE/c^4


def aichelburg_sexl(ck, src):
    """The plane of u and v at x = rho_0/8, y = 0 of the null Cartesian chart, in units of
    8GE/c^4 = rho_0 = 1. Off the shock u = 0 the metric on it is -du dv, flat, and a ray
    moving left crosses the shock with v raised by Delta v = -ln(rho/rho_0) = ln 8, which is
    the published geodesic equation's v'' = -(4GE/c^4) ln((x^2 + y^2)/rho_0^2) delta'(u) u'^2
    integrated across it at fixed x and y. So p = arctan u and q = arctan(v - Delta v theta(u))
    give each half the half of Minkowski's diamond on its side of p = 0, and carry every ray
    moving left across the shock as one line of constant q. The plane is totally geodesic off
    the shock, where every Christoffel symbol vanishes."""
    fixed = {"x": "1/8", "y": "0"}
    params = {"G": 1, "E": "1/8", "rho_0": 1}
    plane = Plane(src, "aichelburg_sexl", "null_cartesian", ("u", "v"), fixed, params, off_shock=True)

    def before(u, w):
        return np.arctan(np.asarray(u, dtype=float)), np.arctan(np.asarray(w, dtype=float))

    def after(u, w):
        return np.arctan(np.asarray(u, dtype=float)), np.arctan(np.asarray(w, dtype=float) - AS_JUMP)
    ck.chart("Aichelburg-Sexl, before the shock", plane, before, ck.uniform(-20, -0.01), ck.uniform(-20, 20),
             lambda u, w: (1, 1))
    ck.chart("Aichelburg-Sexl, behind the shock", plane, after, ck.uniform(0.01, 20), ck.uniform(-20, 20),
             lambda u, w: (1, 1))
    # The jump from the published geodesic equation: v'' = -A delta'(u) with A the coefficient of
    # delta'(u) in Gamma^v_uu, integrated twice across u = 0 at u' = 1, gives Delta v = -A.
    metric, entry, R = nr.load("aichelburg_sexl", "null_cartesian")
    src.note("aichelburg_sexl", "null_cartesian", ["christoffel"])
    gamma = next(c for c in entry["christoffel"]["variants"]["ull"]["nonzero"] if c["indices"] == ["v", "u", "u"])
    names = {R._plain(n): s for n, s in R.symbol.items()}
    value = R(gamma["value"]).subs({R.c: 1, R.parameters["G"]: 1, R.parameters["E"]: sp.Rational(1, 8),
                                    R.parameters["rho_0"]: 1, names["x"]: sp.Rational(1, 8), names["y"]: 0})
    A = value.coeff(sp.DiracDelta(names["u"], 1))
    ck.limit("Aichelburg-Sexl: -Gamma^v_uu integrated twice across the shock is the jump ln 8", [-float(A)], [AS_JUMP], 1e-12)
    w = ck.uniform(-20, 20, 200)
    ck.limit("Aichelburg-Sexl: a ray moving left keeps its q across the shock",
             before(np.full_like(w, -1e-12), w)[1], after(np.full_like(w, 1e-12), w + AS_JUMP)[1], 1e-9)

    box = [-PI - 0.35, PI + 0.35, -PI - 0.25, PI + 0.25]
    v = View("shock", "The shock at $\\rho = \\rho_0/8$", box, "null_cartesian")
    v.fill("region", DIAMOND)
    v.fill("cover", DIAMOND)
    s = spread(-np.inf, np.inf, 600, 10)
    for c in (-4, -2, -1, 1, 2, 4):
        v.curve("null", *before(np.full_like(s, c), s))
    for c in (-4, -2, -1, 0, 1, 2, 4):
        u = -spread(0, np.inf, 400, 10)[::-1]
        v.curve("null", *before(u, np.full_like(u, c)))
        u = spread(0, np.inf, 400, 10)
        v.curve("null", *after(u, np.full_like(u, c)))
    v.segment("surface", (0, -HALF), (0, HALF))
    diamond_edges(v)
    v.label((0, 0.35), "the shock, $u = 0$", "tl", "small", dx=4, dy=3)
    v.legend("null", "$u$ constant and $v$ constant at $0$, $\\pm 1$, $\\pm 2$ and $\\pm 4$, light rays off the shock")
    v.legend("surface", "the shock $u = 0$, a light ray")
    v.set(restriction="The plane $x = \\rho_0/8$, $y = 0$ only, totally geodesic off the shock, each point in the "
                      "diagram a single event.",
          settings="$8GE/c^4 = \\rho_0$, the unit of $u$ and $v$.")
    # Each wave front of the embedding is the null line u = u_k, every v, the same flat front.
    w = spread(-np.inf, np.inf, 400, 10)
    for m in slices.moments("aichelburg_sexl"):
        fmap = before if m.time < 0 else after
        v.slice(m, [fmap(np.full_like(w, m.time), w)])
    return [v]


def domain_wall(ck, src):
    """The domain wall at k = 1, so that 1/k is the unit. Each side is the inside of the
    hyperbola R^2 - c^2T^2 = 1 of Minkowski's plane of T and R, and Minkowski's own maps
    p = arctan(cT - R) and q = arctan(cT + R) send the hyperbola, tan p tan q = -1, to the
    vertical line X = q - p = pi/2, which runs from the middle of I- at T = -pi/2 to the middle
    of I+ at T = pi/2, so each side is the part 0 <= X <= pi/2 of Minkowski's triangle. The side
    z < 0 is drawn so, and the side z > 0 is its mirror image in the wall, X -> pi - X, which is
    p = arctan(cT + R) - pi/2, q = arctan(cT - R) + pi/2; the two meet along the wall with the
    same T, so a ray reaching the wall runs on as one line. In the planar, global and conformal
    charts the point of each side at (t, z) is at R = (1 - |z|) cosh(ct), cT = (1 - |z|) sinh(ct),
    with 1 - |z| = e^{-|w|} in the conformal chart: the region R > c|T| between the wall and the
    light cone of the centre at T = 0, whose edges R = c|T| are the horizons z = -+1."""
    def left(T, R):
        return mink_pq(T, R)

    def right(T, R):
        T, R = np.asarray(T, dtype=float), np.asarray(R, dtype=float)
        return np.arctan(T + R) - HALF, np.arctan(T - R) + HALF

    def sides(T, R, z):
        pl, ql = left(T, R)
        pr, qr = right(T, R)
        z = np.asarray(z, dtype=float)
        return np.where(z > 0, pr, pl), np.where(z > 0, qr, ql)

    def planar(t, z):
        t, z = np.asarray(t, dtype=float), np.asarray(z, dtype=float)
        zeta = 1 - np.abs(z)
        return sides(zeta * np.sinh(t), zeta * np.cosh(t), z)

    def conformal(t, w):
        w = np.asarray(w, dtype=float)
        return planar(t, np.sign(w) * (1 - np.exp(-np.abs(w))))

    pl = Plane(src, "domain_wall", "planar", ("t", "z"), {"x": "0", "y": "0"}, {"k": 1})
    gl = Plane(src, "domain_wall", "global", ("t", "z"), EQUATOR, {"k": 1})
    cf = Plane(src, "domain_wall", "conformal", ("t", "w"), EQUATOR, {"k": 1})
    ip = Plane(src, "domain_wall", "inertial", ("T", "R"), EQUATOR, {"k": 1})
    for name, plane, fmap, lo, hi, x in (("planar", pl, planar, 0.001, 0.999, "z"),
                                         ("global", gl, planar, 0.001, 0.999, "z"),
                                         ("conformal", cf, conformal, 0.001, 8, "w")):
        for side, sign in (("<", -1), (">", 1)):
            ck.chart(f"domain wall {name}, {x} {side} 0", plane, fmap, ck.uniform(-4, 4), sign * ck.uniform(lo, hi),
                     lambda t, z: (1, 0))
    T = ck.uniform(-10, 10)
    ck.chart("domain wall inertial", ip, left, T, np.sqrt(1 + T ** 2) * ck.uniform(0.001, 0.999), lambda T, R: (1, 0))
    ck.chart("domain wall inertial, mirrored", ip, right, T, np.sqrt(1 + T ** 2) * ck.uniform(0.001, 0.999),
             lambda T, R: (1, 0))
    T = np.linspace(-20, 20, 81)
    ck.limit("domain wall: the wall R^2 - c^2T^2 = 1 is the line X = pi/2 on both sides",
             np.concatenate([xt(*left(T, np.sqrt(1 + T ** 2)))[0], xt(*right(T, np.sqrt(1 + T ** 2)))[0]]),
             np.full(162, HALF), 1e-12)
    t = np.linspace(-4, 4, 81)
    ck.limit("domain wall: both sides reach the wall at one T for every t",
             np.concatenate(xt(*planar(t, np.full(81, -1e-15)))), np.concatenate(xt(*planar(t, np.full(81, 1e-15)))),
             1e-12)
    ck.limit("domain wall: the wall ends at the middle of I+ and I-, (X, T) = (pi/2, +-pi/2)",
             np.concatenate([xt(*left(np.array([1e9, -1e9]), np.sqrt(1 + 1e18)))[1]]), [HALF, -HALF], 1e-8)
    X, Tz = xt(*planar(t, np.full(81, -1 + 1e-12)))
    ck.limit("domain wall: the horizon z = -1 is the light cone |T| = X of the centre at T = 0",
             np.abs(Tz) - X, np.zeros(81), 1e-9)
    X, Tz = xt(*planar(t, np.full(81, 1 - 1e-12)))
    ck.limit("domain wall: the horizon z = 1 is the light cone |T| = pi - X of the other centre",
             np.abs(Tz) - (PI - X), np.zeros(81), 1e-9)
    ck.finite("domain wall: the centre R = 0 is regular", ip.kretschmann(ck.uniform(-5, 5, 50), np.full(50, 1e-6)))

    box = [-0.35, PI + 0.35, -PI - 0.25, PI + 0.25]
    region = [[0, -PI], [0, PI], [HALF, HALF], [PI, PI], [PI, -PI], [HALF, -HALF]]
    diamond = [[0, 0], [HALF, HALF], [PI, 0], [HALF, -HALF]]
    pentagon = [[0, -PI], [0, PI], [HALF, HALF], [HALF, -HALF]]
    S = spread(-np.inf, np.inf, 600, 10)

    def frame(v, cover, legend):
        v.fill("region", region)
        v.fill("cover", cover)
        v.line("centre", [[[0, -PI], [0, PI]], [[PI, -PI], [PI, PI]]])
        v.line("scri", [[[0, PI], [HALF, HALF]], [[HALF, HALF], [PI, PI]],
                        [[0, -PI], [HALF, -HALF]], [[HALF, -HALF], [PI, -PI]]])
        v.line("horizon", [[[0, 0], [HALF, HALF]], [[0, 0], [HALF, -HALF]],
                           [[PI, 0], [HALF, HALF]], [[PI, 0], [HALF, -HALF]]])
        v.line("surface", [[[HALF, -HALF], [HALF, HALF]]])
        for at, text, anchor, dy in (((0, PI), "$i^+$", "b", -6), ((PI, PI), "$i^+$", "b", -6),
                                     ((0, -PI), "$i^-$", "t", 6), ((PI, -PI), "$i^-$", "t", 6)):
            v.layers.append({"kind": "point", "class": "infinity", "at": [round(at[0], 4), round(at[1], 4)]})
            v.label_xt(at, text, anchor, dy=dy)
        v.label_xt([Q4, 3 * Q4], "$\\mathscr{I}^+$", "bl", dx=5, dy=-3)
        v.label_xt([3 * Q4, 3 * Q4], "$\\mathscr{I}^+$", "br", dx=-5, dy=-3)
        v.label_xt([Q4, -3 * Q4], "$\\mathscr{I}^-$", "tl", dx=5, dy=3)
        v.label_xt([3 * Q4, -3 * Q4], "$\\mathscr{I}^-$", "tr", dx=-5, dy=3)
        v.label_xt([HALF, 0.3], "the wall", "l", "small", dx=5)
        v.legend("cover", legend)
        v.legend("centre", "the centre of each side, where its spheres shrink to a point")
        v.legend("horizon", "the horizons $z = \\pm 1/k$, the light cones of the two centres at $T = 0$")
        v.legend("surface", "the wall, the hyperbola $R^2 - c^2T^2 = 1/k^2$ of each side")
        v.legend("scri", "null infinity $\\mathscr{I}^\\pm$")

    views = []
    for vid, label, plane_map, coord, values, ts in (
            ("planar", "Planar", planar, "z", (-0.75, -0.5, -0.25, 0.25, 0.5, 0.75), (-1, 0, 1)),
            ("global", "Global", planar, "z", (-0.75, -0.5, -0.25, 0.25, 0.5, 0.75), (-1, 0, 1)),
            ("conformal", "Conformal", conformal, "w", (-2, -1, -0.5, 0.5, 1, 2), (-1, 0, 1))):
        v = View(vid, label, box, vid)
        grid(v, "r", lambda c, s: plane_map(s, c), values, S)
        span = np.linspace(-0.999999, 0.999999, 801) if coord == "z" else np.sinh(np.linspace(-6, 6, 801))
        for t0 in ts:
            for half in (span[span < 0], span[span > 0]):
                v.curve("t", *plane_map(np.full_like(half, t0), half))
        frame(v, diamond, f"the region between the horizons, which $t$ and ${coord}$ cover")
        v.legend("r", f"${coord}$ constant, at " + {"z": "$\\pm 1/4k$, $\\pm 1/2k$, and $\\pm 3/4k$",
                                                   "w": "$\\pm 1/2k$, $\\pm 1/k$, and $\\pm 2/k$"}[coord])
        v.legend("t", "$ct$ constant, at $-1/k$, $0$, and $1/k$")
        views.append(v)

    v = View("inertial", "Inertial", box, "inertial")
    for R0 in (0.5, 1, 2):
        T = np.sinh(np.linspace(-8, 8, 1200))
        inside = np.where(R0 <= np.sqrt(1 + T ** 2), T, np.nan)
        v.curve("r", *left(inside, np.full_like(T, R0)))
    for T0 in (-1, 0, 1):
        R = np.linspace(0, math.sqrt(1 + T0 ** 2), 400)
        v.curve("t", *left(np.full_like(R, T0), R))
    frame(v, pentagon, "the side $z < 0$, inside the wall, which $T$ and $R$ cover")
    v.label_xt([0, 0.25], "$R = 0$", "r", dx=-6)
    v.legend("r", "$R$ constant, at $1/2k$, $1/k$, and $2/k$")
    v.legend("t", "$cT$ constant, at $-1/k$, $0$, and $1/k$")
    views.append(v)
    # Each moment of the embedding, kct = t_k from the centre of one side through the wall to the
    # centre of the other, drawn by the planar map, which is the global chart's on this plane.
    z = np.concatenate([np.linspace(-1 + 1e-9, -1e-12, 400), np.linspace(1e-12, 1 - 1e-9, 400)])
    for v in views:
        v.set(settings="$1/k$, the radius of the wall when it stops, the unit of every length and of $ct$.")
        for m in slices.moments("domain_wall"):
            v.slice(m, [planar(np.full_like(z, m.time), z)])
    return views


DRAWN = {
    "aichelburg_sexl": aichelburg_sexl,
    "domain_wall": domain_wall,
    "minkowski": minkowski, "schwarzschild": schwarzschild, "rn_metric": reissner_nordstrom,
    "kerr": kerr, "kerr_newman": kerr_newman, "de_sitter": de_sitter,
    "schwarzschild_de_sitter": schwarzschild_de_sitter, "global_monopole": global_monopole, "anti_de_sitter": anti_de_sitter,
    "bertotti_robinson": bertotti_robinson, "ellis_bronnikov": ellis_bronnikov, "morris_thorne": morris_thorne,
    "cosmic_string": cosmic_string, "interior_schwarzschild": interior_schwarzschild, "frw": frw,
    "oppenheimer_snyder": oppenheimer_snyder, "vaidya": vaidya, "tov": tov,
    "malament_hogarth": malament_hogarth, "einstein_static": einstein_static, "btz": btz, "c_metric": c_metric,
    "misner": misner, "milne": milne,
    "einstein_rosen_waves": einstein_rosen_waves,
    "nariai": nariai, "khan_penrose": khan_penrose,
    "majumdar_papapetrou": majumdar_papapetrou,
    "melvin": melvin,
}

# ---------------------------------------------------------------- the captions

# Each view's caption, prose under the rules of _tools/README.md: no dashes but in a name,
# and every sentence about the spacetime, never about the page or the collection. A caption
# opens by naming what is drawn, the whole spacetime or the surface in it.
CAPTIONS = {
    ("aichelburg_sexl", "shock"): [
        "The plane of $u$ and $v$ at distance $\\rho = \\rho_0/8$ from the axis the source moves along, brought by "
        "$p = \\arctan u$ and $q = \\arctan(v - \\Delta v\\,\\theta(u))$ into the whole diamond, with "
        "$\\Delta v = (8GE/c^4)\\ln 8$ the jump of a ray moving left across the shock and $\\theta$ the unit step. "
        "The shock $u = 0$ is a light ray from $\\mathscr{I}^-$ to $\\mathscr{I}^+$, with flat Minkowski space on "
        "either side of it.",
        "Every ray moving left keeps its $q$ across the shock, and every line of constant $v$ breaks there, its part "
        "behind the shock moved along it by $\\Delta v$. The rays moving right never meet the shock.",
    ],
    ("domain_wall", "planar"): [
        "The whole spacetime of the wall, each point in the diagram a 2-sphere, two copies of the inside of the "
        "hyperbola $R^2 - c^2T^2 = 1/k^2$ in Minkowski's triangle joined along it, the second drawn as the mirror "
        "image of the first. The wall is the vertical line in the middle, from the middle of $\\mathscr{I}^-$ to "
        "the middle of $\\mathscr{I}^+$, and each side has its own centre, $i^\\pm$, and null infinity, and no "
        "spatial infinity.",
        "The planar chart reaches the diamond between the wall and the horizons $z = \\pm 1/k$, the light cones of "
        "the two centres at $T = 0$, where $(1 - k|z|)\\cosh(kct)$ and $(1 - k|z|)\\sinh(kct)$ are $kR$ and $kcT$. "
        "A light ray crosses the wall as one straight line.",
    ],
    ("domain_wall", "global"): [
        "The whole spacetime of the wall, each point in the diagram a 2-sphere, two copies of the inside of the "
        "hyperbola $R^2 - c^2T^2 = 1/k^2$ in Minkowski's triangle joined along it, the second drawn as the mirror "
        "image of the first. The global chart reaches the diamond between the wall and the horizons $z = \\pm 1/k$, "
        "where the spheres of constant $t$ and $z$ have radius $(1 - k|z|)\\cosh(kct)/k$.",
        "Each line of constant $t$ runs from the centre of one side at $T = 0$ through the wall to the centre of the "
        "other, and together with its spheres it closes up into a 3-sphere, the whole of space at one moment.",
    ],
    ("domain_wall", "conformal"): [
        "The whole spacetime of the wall, each point in the diagram a 2-sphere, two copies of the inside of the "
        "hyperbola $R^2 - c^2T^2 = 1/k^2$ in Minkowski's triangle joined along it, the second drawn as the mirror "
        "image of the first. The conformal chart reaches the diamond between the wall and the horizons, which lie "
        "at $w \\to \\pm\\infty$, and on its plane of $t$ and $w$ the metric is $e^{-2k|w|}(-c^2dt^2 + dw^2)$, so "
        "its light rays are $ct \\pm w = $ const.",
    ],
    ("domain_wall", "inertial"): [
        "The whole spacetime of the wall, each point in the diagram a 2-sphere of radius $R$ on the side the inertial "
        "chart covers, two copies of the inside of the hyperbola $R^2 - c^2T^2 = 1/k^2$ in Minkowski's triangle "
        "joined along it, the second drawn as the mirror image of the first. The chart is Minkowski's spherical "
        "chart inside the wall, drawn with Minkowski's own maps $p = \\arctan(k(cT - R))$ and "
        "$q = \\arctan(k(cT + R))$, which send the hyperbola to the vertical line $q - p = \\pi/2$.",
        "Its lines of constant $R$ beyond $1/k$ start and end on the wall, which reaches them only for $c|T| > "
        "\\sqrt{R^2 - 1/k^2}$.",
    ],
    ("c_metric", "inner_spherical"): [
        "The half axis $\\theta = 0$ between the black holes of the maximally extended C-metric ($\\alpha m = 1/6$), "
        "totally geodesic, each point in the diagram a single event. On it the metric is "
        "$(-Q\\,c^2dt^2 + dr^2/Q)/(1 + \\alpha r)^2$ with $Q = (1 - \\alpha^2r^2)(1 - 2m/r)$, whose roots are the "
        "black hole horizon $r = 2m$ and the acceleration horizon $r = 1/\\alpha$, and maps $p$ and $q$ of "
        "$ct \\mp r_*$, with $dr_* = dr/Q$, bring each region into a finite diamond, drawn with $T = p + q$ up and "
        "$X = q - p$ across. The maps are Kruskal's across the acceleration horizon and continuous across the black "
        "hole horizon.",
        "The two static regions in the middle are the exteriors of the two black holes, causally separated by the "
        "acceleration horizon. Above and below them lie the regions beyond the acceleration horizon, whose outer "
        "edges are null infinity, where $1 + \\alpha r$ vanishes at $r = -1/\\alpha$, past $r = \\infty$. Beyond "
        "each black hole horizon lie the black hole and the white hole, which end on the singularity $r = 0$, and "
        "the exterior across the Einstein-Rosen bridge, from which the chain of regions repeats in both directions. "
        "Jerry Griffiths, Pavel Krtouš, and Jiří Podolský drew this diagram in 2006.",
        "The coordinates $t$ and $r > 0$ cover one static region, its black hole, and the part of the region "
        "beyond the acceleration horizon out to $r = \\infty$.",
    ],
    ("c_metric", "inner_hong_teo"): [
        "The half axis $\\theta = 0$ between the black holes of the maximally extended C-metric ($\\alpha m = 1/6$), "
        "totally geodesic, each point in the diagram a single event. On it the metric is "
        "$(-Q\\,c^2dt^2 + dr^2/Q)/(1 + \\alpha r)^2$ with $Q = (1 - \\alpha^2r^2)(1 - 2m/r)$, whose roots are the "
        "black hole horizon $r = 2m$ and the acceleration horizon $r = 1/\\alpha$, and maps $p$ and $q$ of "
        "$ct \\mp r_*$, with $dr_* = dr/Q$, bring each region into a finite diamond, drawn with $T = p + q$ up and "
        "$X = q - p$ across. The maps are Kruskal's across the acceleration horizon and continuous across the black "
        "hole horizon.",
        "The two static regions in the middle are the exteriors of the two black holes, causally separated by the "
        "acceleration horizon. Above and below them lie the regions beyond the acceleration horizon, whose outer "
        "edges are null infinity, where $1 + \\alpha r$ vanishes at $r = -1/\\alpha$, past $r = \\infty$. Beyond "
        "each black hole horizon lie the black hole and the white hole, which end on the singularity $r = 0$, and "
        "the exterior across the Einstein-Rosen bridge, from which the chain of regions repeats in both directions. "
        "Jerry Griffiths, Pavel Krtouš, and Jiří Podolský drew this diagram in 2006.",
        "The coordinates $\\tau$ and $y > -1$ cover one static region, its black hole, and the whole region "
        "beyond the acceleration horizon, out to null infinity at $y = -1$.",
    ],
    ("c_metric", "outer_spherical"): [
        "The half axis $\\theta = \\pi$ beyond the black hole of the maximally extended C-metric ($\\alpha m = 1/6$), "
        "totally geodesic, each point in the diagram a single event. There the metric is "
        "$(-Q\\,c^2dt^2 + dr^2/Q)/(1 - \\alpha r)^2$, whose conformal factor diverges at $r = 1/\\alpha$, so on "
        "this half of the axis the acceleration horizon lies at null infinity and bounds the static region on the "
        "left, as null infinity bounds the exterior of Schwarzschild's black hole.",
        "Beyond the black hole horizon $r = 2m$ lie the black hole and the white hole, which end on the singularity "
        "$r = 0$, and the exterior across the Einstein-Rosen bridge, whose own outer axis runs to null infinity on "
        "the right. With $C = 1/(1 + 2\\alpha m)$ this half of the axis carries the cosmic string, whose deficit "
        "angle lies in the angle about the axis and leaves this surface unchanged.",
        "The coordinates $t$ and $r$ cover the static region and its black hole.",
    ],
    ("c_metric", "outer_hong_teo"): [
        "The half axis $\\theta = \\pi$ beyond the black hole of the maximally extended C-metric ($\\alpha m = 1/6$), "
        "totally geodesic, each point in the diagram a single event. There the metric is "
        "$(-Q\\,c^2dt^2 + dr^2/Q)/(1 - \\alpha r)^2$, whose conformal factor diverges at $r = 1/\\alpha$, so on "
        "this half of the axis the acceleration horizon lies at null infinity and bounds the static region on the "
        "left, as null infinity bounds the exterior of Schwarzschild's black hole.",
        "Beyond the black hole horizon $r = 2m$ lie the black hole and the white hole, which end on the singularity "
        "$r = 0$, and the exterior across the Einstein-Rosen bridge, whose own outer axis runs to null infinity on "
        "the right. With $C = 1/(1 + 2\\alpha m)$ this half of the axis carries the cosmic string, whose deficit "
        "angle lies in the angle about the axis and leaves this surface unchanged.",
        "The coordinates $\\tau$ and $y > 1$ cover the static region and its black hole.",
    ],
    ("milne", "comoving_hyperbolic"): [
        "The Milne universe, each point in the diagram a 2-sphere of radius $ct\\sinh\\chi$. With $T = "
        "t\\cosh\\chi$ and $R = ct\\sinh\\chi$ it is the inside of the future light cone of the event $T = R = 0$ "
        "of Minkowski spacetime, and Minkowski's own maps $p = \\arctan(ct\\,e^{-\\chi}/\\ell)$ and $q = "
        "\\arctan(ct\\,e^{\\chi}/\\ell)$ draw it as the wedge above that cone, with $T = p + q$ up and $X = q - p$ "
        "across.",
        "Every comoving particle runs from the event at the foot of the wedge to $i^+$, and every moment of "
        "constant $t$ runs from the centre to the corner where the cone meets $\\mathscr{I}^+$. Light crosses the "
        "cone into the Milne universe from the rest of Minkowski spacetime and leaves it through $\\mathscr{I}^+$.",
    ],
    ("milne", "comoving_spherical"): [
        "The Milne universe in its comoving spherical chart, each point in the diagram a 2-sphere of radius $ctr$. "
        "With $r = \\sinh\\chi$ the chart is the hyperbolic one, and Minkowski's maps draw it as the same wedge "
        "above the future light cone of the event $T = R = 0$.",
    ],
    ("milne", "logarithmic_time"): [
        "The Milne universe in its logarithmic time, each point in the diagram a 2-sphere of radius "
        "$ct_0e^{\\tau/t_0}\\sinh\\chi$. The metric is $e^{2\\tau/t_0}$ times that of a static universe whose "
        "space is the hyperbolic space of radius $ct_0$, and the whole line $\\tau \\to -\\infty$ goes to the "
        "event at the foot of the wedge, where $t = t_0e^{\\tau/t_0}$ vanishes.",
    ],
    ("milne", "inertial"): [
        "The Milne universe in the inertial chart of the comoving particle at $\\chi = 0$, each point in the "
        "diagram a 2-sphere of radius $R$. The chart is Minkowski's spherical chart inside $R < cT$, drawn with "
        "Minkowski's own maps, and its lines of constant $R$ enter the Milne universe across the light cone $cT = R$ "
        "and run on to $i^+$.",
    ],
    ("minkowski", "spherical"): [
        "Minkowski spacetime, each point in the diagram a 2-sphere of radius $r$. With $u = ct - r$ and $v = ct + r$ the metric on the plane of $t$ "
        "and $r$ is $-du\\,dv$, and for any length $\\ell$ the maps $p = \\arctan(u/\\ell)$ and "
        "$q = \\arctan(v/\\ell)$ bring all of it into a finite triangle, drawn with $T = p + q$ up "
        "and $X = q - p$ across. Lines of constant $p$ or $q$ are light rays, and they run at "
        "45°.",
        "Every line of constant $r$ runs from $i^-$ to $i^+$, and every line of constant $t$ from "
        "the centre to $i^0$. Ingoing light starts on $\\mathscr{I}^-$, passes through the centre "
        "and ends on $\\mathscr{I}^+$.",
    ],
    ("minkowski", "spherical_null"): [
        "Minkowski spacetime in its spherical null coordinates ($u = t - r/c$, $v = t + r/c$), "
        "each point in the diagram a 2-sphere of radius $c(v - u)/2$. The lines of constant $u$ and of constant $v$ are light rays, the "
        "45° lines of the triangle, with $p = \\arctan(cu/\\ell)$ and $q = \\arctan(cv/\\ell)$, and "
        "the centre is the line $u = v$.",
    ],
    ("minkowski", "cartesian"): [
        "The plane $y = z = 0$, flat and totally geodesic, brought by "
        "$p, q = \\arctan((ct \\mp x)/\\ell)$ into the whole diamond. It has two ends, "
        "$x \\to +\\infty$ and $x \\to -\\infty$, each with its own null infinity.",
        "Turned about the line $x = 0$, each half of the diamond sweeps out the spherical "
        "triangle, which is the whole spacetime.",
    ],
    ("minkowski", "double_null"): [
        "The plane $y = z = 0$ in the coordinates $u = t - x/c$ and $v = t + x/c$, the metric "
        "on it $-c^2\\,du\\,dv$. The lines of constant $u$ and of "
        "constant $v$ are light rays, the 45° lines of the diamond, with $p = \\arctan(cu/\\ell)$ "
        "and $q = \\arctan(cv/\\ell)$.",
    ],
    ("minkowski", "rindler"): [
        "The plane $Y = Z = 0$ in Rindler's coordinates ($ct = X\\sinh(aT/c)$, "
        "$x = X\\cosh(aT/c)$), with $ct - x = -Xe^{-aT/c}$ and $ct + x = Xe^{aT/c}$. They cover "
        "the wedge $x > c|t|$, and $X = 0$, where $g_{TT}$ vanishes, is the pair of null lines "
        "through the origin.",
        "An observer at constant $X$ accelerates uniformly, at $c^2/X$, and the null line "
        "$ct = x$ is that observer's horizon: no event beyond it can send a signal into the "
        "wedge, as nothing inside $r_s$ can reach a static observer outside a black hole.",
    ],
    ("misner", "misner"): [
        "The plane $y = z = 0$ of the Minkowski space that covers Misner space, in Misner's coordinates, with "
        "$ct - x = -2\\ell e^{-\\psi/2}$ and $ct + x = 2Te^{\\psi/2}/\\ell$ for any length $\\ell$, brought into a "
        "finite drawing by $p = \\arctan((ct - x)/\\ell)$ and $q = \\arctan((ct + x)/\\ell)$. The coordinates cover "
        "the half $x > ct$: the past light cone of the origin, where $T < 0$, and the wedge $x > c|t|$, where "
        "$T > 0$.",
        "The boost of rapidity $\\psi_0/2$ carries each light ray $\\psi = k\\psi_0$ onto the next, and one copy of "
        "Misner space lies between two neighbouring rays, which are one ray of it. The null line $T = 0$ is the "
        "chronology horizon, and in the quotient each hyperbola of constant $T > 0$ is a closed timelike curve.",
    ],
    ("misner", "milne"): [
        "The past light cone of the origin in the plane $y = z = 0$ of the covering Minkowski space, in the "
        "Milne chart, $ct\\cosh\\chi$ and $ct\\sinh\\chi$ with $t < 0$. Lines of constant $\\chi$ run straight "
        "into the origin, and the boost carries each line $\\chi = k\\psi_0/2$ onto the next.",
        "Each hyperbola of constant $t$ is, in the quotient, a circle of circumference $\\psi_0c|t|/2$, and the "
        "circles shrink toward the null lines $t \\to 0$. Misner's coordinates continue the region across the "
        "line $ct + x = 0$ into the wedge $x > c|t|$, and the other extension across $ct - x = 0$.",
    ],
    ("misner", "rindler"): [
        "The wedge $x > c|t|$ in the plane $y = z = 0$ of the covering Minkowski space, in the Rindler chart, "
        "$\\xi\\sinh\\eta$ and $\\xi\\cosh\\eta$. The boost carries each line $\\eta = k\\psi_0/2$ onto the next, "
        "so each hyperbola of constant $\\xi$ is, in the quotient, a closed timelike curve of proper length "
        "$\\xi\\psi_0/2$.",
        "The closed curves shrink toward the null lines $\\xi \\to 0$, the chronology horizon, where they become "
        "closed null geodesics.",
    ],
    ("schwarzschild", "spherical"): [
        "The Schwarzschild spacetime, maximally extended, each point in the diagram a 2-sphere "
        "of radius $r$. The coordinates of Martin Kruskal and "
        "George Szekeres, $U = -e^{-u/2r_s}$ and $V = e^{v/2r_s}$, with $u, v = ct \\mp r_*$ and "
        "$r_* = r + r_s\\ln|r/r_s - 1|$, make the metric regular through $r = r_s$, where $UV = (1 "
        "- r/r_s)e^{r/r_s}$ vanishes. With $p = \\arctan U$ and $q = \\arctan V$ the singularity "
        "$UV = 1$ lies exactly on the straight lines $T = \\pm\\pi/2$, since $\\tan(p + q) = (U + "
        "V)/(1 - UV)$ diverges there.",
        "The coordinates $t$ and $r > r_s$ cover the right exterior alone. The horizon is the "
        "pair of null lines $U = 0$ and $V = 0$, crossing at the bifurcation sphere. The black "
        "hole above it ends at $r = 0$ on $T = \\pi/2$, the white hole below it begins at "
        "$r = 0$ on $T = -\\pi/2$, and both singularities are spacelike.",
    ],
    ("schwarzschild", "ingoing"): [
        "The whole Schwarzschild spacetime with the ingoing Eddington-Finkelstein "
        "coordinates $v$ and $r$ on it. From $V = e^{v/2r_s}$ and $U = (1 - r/r_s)e^{r/r_s}/V$, "
        "one formula for every $r > 0$, they cover the exterior and the black hole together, "
        "and their lines of constant $v$ are ingoing light rays, which cross the horizon at "
        "45° and end at $r = 0$.",
    ],
    ("schwarzschild", "outgoing"): [
        "The whole Schwarzschild spacetime with the outgoing Eddington-Finkelstein "
        "coordinates $u$ and $r$ on it, the time reverse of the ingoing ones. From "
        "$U = -e^{-u/2r_s}$ and $V = (r/r_s - 1)e^{r/r_s}/(-U)$ they cover the exterior and the "
        "white hole, and their lines of constant $u$ are outgoing light rays, which leave $r = 0$ "
        "and cross the horizon outward.",
    ],
    ("global_monopole", "static"): [
        "Letelier's black hole in a cloud of strings, maximally extended ($\\Delta = 0.19$), each point in the "
        "diagram a 2-sphere of area $4\\pi r^2$. On the plane of $t$ and $r$ the metric is $1/(1 - \\Delta)$ "
        "times Schwarzschild's with $r_h = r_s/(1 - \\Delta)$ in place of $r_s$ and $(1 - \\Delta)\\,ct$ in place "
        "of $ct$, so Kruskal and Szekeres's $U = -e^{-(1 - \\Delta)u/2r_h}$ and $V = e^{(1 - \\Delta)v/2r_h}$, with "
        "$u, v = ct \\mp r_*$, draw it as they draw Schwarzschild's, and $p = \\arctan U$ and $q = \\arctan V$ put "
        "the singularity $UV = 1$ on the straight lines $T = \\pm\\pi/2$.",
        "The coordinates $t$ and $r > r_h$ cover the right exterior alone. The horizon is the pair of null lines "
        "$U = 0$ and $V = 0$, crossing at the bifurcation sphere, and the black hole above it ends at $r = 0$ on "
        "$T = \\pi/2$. Far out the proper distance between neighbouring spheres is $dr/\\sqrt{1 - \\Delta}$, so a "
        "sphere of area $4\\pi r^2$ lacks the solid angle $4\\pi\\Delta$ of a Euclidean one.",
    ],
    ("global_monopole", "ingoing"): [
        "The whole of Letelier's black hole with the ingoing Eddington-Finkelstein coordinates $v$ and $r$ on it. "
        "From $V = e^{(1 - \\Delta)v/2r_h}$ and $U = (1 - r/r_h)e^{r/r_h}/V$, one formula for every $r > 0$, they "
        "cover the exterior and the black hole together, and their lines of constant $v$ are ingoing light rays, "
        "which cross the horizon at 45° and end at $r = 0$.",
    ],
    ("global_monopole", "outgoing"): [
        "The whole of Letelier's black hole with the outgoing Eddington-Finkelstein coordinates $u$ and $r$ on "
        "it, the time reverse of the ingoing ones. From $U = -e^{-(1 - \\Delta)u/2r_h}$ and "
        "$V = (r/r_h - 1)e^{r/r_h}/(-U)$ they cover the exterior and the white hole, and their lines of constant "
        "$u$ are outgoing light rays, which leave $r = 0$ and cross the horizon outward.",
    ],
    ("global_monopole", "conical"): [
        "The global monopole with no mass at its centre ($\\Delta = 0.19$), each point in the diagram a 2-sphere "
        "of area $4\\pi(1 - \\Delta)r^2$. On the plane of $t$ and $r$ the metric is $-c^2dt^2 + dr^2$, "
        "Minkowski's, and $p = \\arctan((ct - r)/\\ell)$ and $q = \\arctan((ct + r)/\\ell)$ bring it into "
        "Minkowski's triangle, drawn with $T = p + q$ up and $X = q - p$ across.",
        "The deficit enters only $g_{\\theta\\theta}$ and $g_{\\phi\\phi}$, so the plane and its triangle are "
        "Minkowski's. The edge $X = 0$ is the monopole, $r = 0$, where the Kretschmann scalar "
        "$4\\Delta^2/\\left((1 - \\Delta)^2r^4\\right)$ diverges, a timelike singularity.",
    ],
    ("btz", "static"): [
        "The black hole without rotation ($M = 1$, $J = 0$), maximally extended, each point in the "
        "diagram a circle of circumference $2\\pi r$. With $r_* = \\frac{\\ell}{2\\sqrt{M}}\\ln|(r - r_+)/(r + "
        "r_+)|$, which vanishes as $r \\to \\infty$, the Kruskal coordinates $U = -e^{-\\kappa u}$ and $V = "
        "e^{\\kappa v}$, with $u, v = ct \\mp r_*$ and $\\kappa = \\sqrt{M}/\\ell$, make the metric regular "
        "through $r_+$, where $UV = (r_+ - r)/(r_+ + r)$ vanishes. With $p = \\arctan U$ and $q = \\arctan V$ "
        "the conformal boundary, $UV = -1$, lies on the vertical lines $X = \\pm\\pi/2$, and $r = 0$, where "
        "$UV = 1$, on the horizontal lines $T = \\pm\\pi/2$: the square of Máximo Bañados, Marc Henneaux, "
        "Claudio Teitelboim, and Jorge Zanelli.",
        "The coordinates $t$ and $r > r_+$ cover the right exterior alone. The horizon is the pair of null "
        "lines $U = 0$ and $V = 0$, crossing at the bifurcation circle. The black hole above it ends at $r = "
        "0$ on $T = \\pi/2$ and the white hole below it begins at $r = 0$ on $T = -\\pi/2$, where the circles "
        "shrink to zero length while the curvature stays $R = -6/\\ell^2$; continued past $r = 0$, the "
        "circles would be closed timelike curves.",
    ],
    ("btz", "ingoing"): [
        "The black hole without rotation ($M = 1$, $J = 0$) with the ingoing Eddington-Finkelstein "
        "coordinates $v$ and $r$ on it. From $V = e^{\\kappa v}$ and $U = (r_+ - r)/((r_+ + r)V)$, one "
        "formula for every $r > 0$, they cover the exterior and the black hole together, and their lines of "
        "constant $v$ are ingoing light rays, which cross the horizon at 45° and end at $r = 0$.",
    ],
    ("btz", "outgoing"): [
        "The black hole without rotation ($M = 1$, $J = 0$) with the outgoing Eddington-Finkelstein "
        "coordinates $u$ and $r$ on it, the time reverse of the ingoing ones. From $U = -e^{-\\kappa u}$ and "
        "$V = (r - r_+)/((r + r_+)(-U))$ they cover the exterior and the white hole, and their lines of "
        "constant $u$ are outgoing light rays, which leave $r = 0$ and cross the horizon outward.",
    ],
    ("btz", "rotating"): [
        "The rotating black hole ($M = 1$, $J = 4\\ell/5$), maximally extended on its plane of $t$ and $r$ "
        "with $\\phi$ divided out, $-N^2c^2dt^2 + dr^2/N^2$, each point in the diagram a circle of "
        "circumference $2\\pi r$. The extension is a tower of regions that repeats up and down without end. "
        "Its tortoise coordinate is $r_* = \\frac{1}{2\\kappa_+}\\ln\\left|\\frac{r - r_+}{r + r_+}\\right| - "
        "\\frac{1}{2\\kappa_-}\\ln\\left|\\frac{r - r_-}{r + r_-}\\right|$ with $\\kappa_\\pm = (r_+^2 - "
        "r_-^2)/(\\ell^2r_\\pm)$, which vanishes both at $r = 0$ and as $r \\to \\infty$.",
        "We place every region by the Kruskal coordinate of the outer horizon, $p = \\pm\\arctan "
        "e^{-\\kappa_+u}$ and $q = \\pm\\arctan e^{\\kappa_+v}$ with $u, v = ct \\mp r_*$, and the regions above "
        "the inner horizon are the reflection $(p, q) \\to (\\pi - q, \\pi - p)$ of those below. Since $r_*$ "
        "vanishes at both ends, the conformal boundary beside each exterior and $r = 0$ beside each region "
        "inside $r_-$ lie on the same vertical lines $X = \\pm\\pi/2$, and both are timelike. Across $r_-$ "
        "the map is continuous and cannot also be smooth, because the late light rays that reach the "
        "boundary are the rays that pile up at the inner horizon, and one function of the ray has to serve "
        "both.",
        "The coordinates $t$ and $r > 0$ cover one region of each kind: an exterior, the black hole between "
        "the horizons, and a region inside $r_-$. At $r = 0$ the circles shrink to zero length while the "
        "curvature stays $R = -6/\\ell^2$, and continued past it they would be closed timelike curves.",
    ],
    ("schwarzschild_de_sitter", "static"): [
        "Kottler's spacetime, maximally extended, each point in the diagram a 2-sphere of radius $r$. "
        "Its tortoise coordinate is $r_* = \\sum_i \\ln|1 - r/r_i|/f'(r_i)$ over the three roots $r_i$ of "
        "$f = 1 - r_s/r - \\Lambda r^2/3$, the third of them negative, so that $r_*(0) = 0$, and $u, v = ct \\mp r_*$ "
        "in each region. We place every region by $p = \\pm\\arctan W(u)$ and $q = \\pm\\arctan W(-v)$, shifted by "
        "$\\pi$ where the region lies beyond a cosmological horizon, with $W(x) = \\exp(-ax - b\\sqrt{x^2 + r_s^2})$ "
        "and $a, b = (\\kappa_h \\pm \\kappa_c)/2$: toward each horizon $W$ is its Kruskal coordinate, $e^{-\\kappa_hx}$ "
        "toward $r_h$ and $e^{-\\kappa_cx}$ toward $r_c$, so every line crosses both kinds of horizon with a "
        "continuous tangent.",
        "The coordinates $t$ and $r_h < r < r_c$ cover one static region. Static regions alternate along the chain "
        "with the black hole above the white hole across $r_h$ and the expanding region above the contracting one "
        "across $r_c$, and the chain runs on past both ends of the drawing without end. The singularity $r = 0$ and "
        "future and past infinity are spacelike curves, and neither is straight, since $\\kappa_h \\neq \\kappa_c$.",
    ],
    ("schwarzschild_de_sitter", "ingoing"): [
        "Kottler's spacetime with the ingoing Eddington-Finkelstein coordinates $v$ and $r$ on it. With "
        "$q = \\arctan W(-v)$ and $u = v - 2r_*$, one chart covers the contracting region, the static region and "
        "the black hole together, and its lines of constant $v$ are ingoing light rays, which start on "
        "$\\mathscr{I}^-$, cross both horizons at 45° and end at $r = 0$.",
    ],
    ("schwarzschild_de_sitter", "outgoing"): [
        "Kottler's spacetime with the outgoing Eddington-Finkelstein coordinates $u$ and $r$ on it, the time "
        "reverse of the ingoing ones. With $p = -\\arctan W(u)$ and $v = u + 2r_*$ they cover the white hole, the "
        "static region and the expanding region, and their lines of constant $u$ are outgoing light rays, which "
        "leave $r = 0$, cross both horizons outward and end on $\\mathscr{I}^+$.",
    ],
    ("majumdar_papapetrou", "one_hole"): [
        "One hole alone ($U = 1 + m/r$), which is the extremal Reissner-Nordström black hole, maximally "
        "extended, each point in the diagram a 2-sphere of areal radius $R = r + m$. The extension is a "
        "tower of exteriors and interiors that repeats up and down without end. Outside the horizon the "
        "tortoise coordinate is $r_* = r + 2m\\ln(r/m) - m^2/r$, and inside it "
        "$R_* = R + 2m\\ln|R/m - 1| - m^2/(R - m) - m$, which vanishes at $R = 0$.",
        "We place each region by $p = \\arctan(u/m)$ and $q = \\arctan(v/m)$ with $u, v = ct \\mp r_*$, "
        "shifted by $\\pi$ from one region to the next. The horizon is a double root of $g^{rr}$, where the "
        "surface gravity vanishes, and the map is continuous across it. The singularity $R = 0$ is timelike "
        "and lies on the vertical line $X = -\\pi$.",
        "The coordinates $t$ and $r > 0$ cover one exterior. Its moment $t = 0$ runs from spatial infinity "
        "down the infinitely long throat toward the corner $X = -\\pi$, $T = 0$, where the past and future "
        "horizons meet at an infinite distance.",
    ],
    ("rn_metric", "tower"): [
        "The Reissner-Nordström spacetime, maximally extended, each point in the diagram a "
        "2-sphere of radius $r$. The extension is a tower of regions that "
        "repeats up and down without end. Its tortoise coordinate is "
        "$r_* = r + \\frac{1}{2\\kappa_+}\\ln|r/r_+ - 1| - \\frac{1}{2\\kappa_-}\\ln|r/r_- - 1|$ "
        "with $\\kappa_\\pm = (r_+ - r_-)/2r_\\pm^2$, and each logarithm is made dimensionless by "
        "its own root, so that $r_*(0) = 0$.",
        "We place every region by the Kruskal coordinate of the outer horizon, $p = \\pm\\arctan "
        "e^{-\\kappa_+ u}$ and $q = \\pm\\arctan e^{\\kappa_+ v}$ with $u, v = t \\mp r_*$, and the "
        "regions above the inner horizon are the reflection $(p, q) \\to (\\pi - q, \\pi - p)$ of "
        "those below. This map is smooth across $r_+$ and puts the singularity $r = 0$ exactly on the "
        "vertical lines $X = \\pm\\pi/2$, where it is timelike. Across $r_-$ it is continuous and "
        "cannot also be smooth, because the late light rays that reach $\\mathscr{I}^+$ are the rays "
        "that pile up at the Cauchy horizon $r_-$, and one function of the ray has to serve both.",
        "The coordinates $t$ and $r > 0$ cover one region of each kind: an exterior, a region "
        "between the horizons, and a region inside $r_-$. Between the horizons their $t$ alone "
        "cannot tell the black hole from the white hole, and we take the region to be the black "
        "hole an infalling observer enters.",
    ],
    ("rn_metric", "malament_hogarth"): [
        "The same tower with one event beyond the Cauchy horizon $r_-$ marked on it. "
        "Every point of the exterior below it has $p \\le p_e$ and $q \\le q_e$, so the whole "
        "exterior lies in the event's causal past. A static observer at $r = 1.2\\,r_s$ lives "
        "from $i^-$ to $i^+$ for an infinite proper time, and every moment of that life can send "
        "a signal that reaches the event, arriving at the Cauchy horizon infinitely "
        "blueshifted.",
    ],
    ("kerr", "axis"): [
        "The symmetry axis $\\theta = 0$ of the maximally extended Kerr spacetime, the "
        "surface the rotations leave fixed and so totally geodesic. On it the metric is "
        "$-\\frac{\\Delta}{r^2 + a^2}c^2dt^2 + \\frac{r^2 + a^2}{\\Delta}dr^2$ with $\\Delta = r^2 "
        "- 2GMr/c^2 + a^2$, and its two simple roots give it the tower of Reissner-Nordström. "
        "Brandon Carter extended the axis this way in 1966.",
        "Where Reissner-Nordström ends at $r = 0$, the axis runs on through the centre of the "
        "ring's disc, where the curvature is finite, into $r < 0$, a second asymptotically flat "
        "end with its own null infinity. The ring singularity itself is at $r = 0$ in the "
        "equatorial plane $\\theta = \\pi/2$. The coordinates $t$ and $r > r_+$ cover the "
        "exterior.",
    ],
    ("kerr_newman", "axis"): [
        "The symmetry axis $\\theta = 0$ of the maximally extended Kerr-Newman spacetime, "
        "totally geodesic as Kerr's is. On it the metric is again "
        "$-\\frac{\\Delta}{r^2 + a^2}c^2dt^2 + \\frac{r^2 + a^2}{\\Delta}dr^2$, now with "
        "$\\Delta = r^2 - 2GMr/c^2 + a^2 + r_Q^2$, and its two simple roots give the same tower.",
        "The axis runs through the centre of the ring's disc into $r < 0$ as Kerr's does, and "
        "the ring singularity is at $r = 0$ in the equatorial plane $\\theta = \\pi/2$. The "
        "coordinates $t$ and $r > r_+$ cover the exterior.",
    ],
    ("nariai", "static"): [
        "The Nariai universe, the product of the hyperboloid $-Z_0^2 + Z_1^2 + Z_2^2 = 1/\\Lambda$ with a sphere of "
        "radius $1/\\sqrt{\\Lambda}$, each point in the diagram a 2-sphere of that radius. With "
        "$\\tan\\eta = \\sqrt{\\Lambda}\\,Z_0$ the metric of the first factor is "
        "$(-d\\eta^2 + d\\chi^2)/(\\Lambda\\cos^2\\eta)$ on the strip $|\\eta| < \\pi/2$, with $X = \\chi$ across "
        "and its edges $\\chi = 0$ and $2\\pi$ one line.",
        "The static chart covers the diamond about $\\chi = \\pi/2$, and its horizons $r = \\pm 1/\\sqrt{\\Lambda}$ "
        "are the diamond's four edges, two for each, which meet the antipode's at $\\chi = 0$ and $\\pi$. Infinity is spacelike, "
        "past and future, and a light ray crosses half the circle between them, so the observers at "
        "$\\chi = \\pi/2$ and $3\\pi/2$ never exchange a signal.",
    ],
    ("nariai", "global"): [
        "The Nariai universe in its global chart, on the strip $|\\eta| < \\pi/2$ with $X = \\chi$ and "
        "$\\tan\\eta = \\sinh(\\sqrt{\\Lambda}\\,ct)$, each point in the diagram a 2-sphere of radius "
        "$1/\\sqrt{\\Lambda}$. The chart covers the whole spacetime, and on each line of constant $t$ the circle of "
        "$\\chi$ has the radius $\\cosh(\\sqrt{\\Lambda}\\,ct)/\\sqrt{\\Lambda}$.",
    ],
    ("nariai", "conformal"): [
        "The Nariai universe in its conformal chart, $(-d\\eta^2 + d\\chi^2)/(\\Lambda\\cos^2\\eta) + "
        "d\\Omega^2/\\Lambda$, drawn with $X = \\chi$ and $T = \\eta$ on the strip $|\\eta| < \\pi/2$, each point in the "
        "diagram a 2-sphere of radius $1/\\sqrt{\\Lambda}$. The light rays are the lines of constant "
        "$\\eta \\pm \\chi$, at 45°, and $\\eta = \\pm\\pi/2$ is the infinite future and past, where "
        "$1/\\cos^2\\eta$ diverges.",
    ],
    ("de_sitter", "static"): [
        "De Sitter spacetime, the hyperboloid $-X_0^2 + X_1^2 + \\dots + X_4^2 = L^2$ "
        "($L = \\sqrt{3/\\Lambda}$), each point in the diagram a 2-sphere. In its global coordinates the metric is $\\frac{L^2}{\\cos^2 T}(-dT^2 + d\\chi^2 "
        "+ \\sin^2\\chi\\,d\\Omega^2)$ on the square $|T| < \\pi/2$, $0 \\le \\chi \\le \\pi$, "
        "with $X = \\chi$ across.",
        "The static coordinates enter the square as $\\tan p = \\tanh(u/2L)$ and "
        "$\\tan q = \\tanh(v/2L)$, with $u, v = ct \\mp L\\,\\mathrm{artanh}(r/L)$, and cover the "
        "triangle about the observer at $\\chi = 0$; their horizon $r = L$ is the pair of null "
        "lines through the centre of the square. Infinity is spacelike, past and future, so every "
        "observer has an event horizon: nothing beyond the line from $(\\chi, T) = (0, \\pi/2)$ to "
        "$(\\pi, -\\pi/2)$ ever reaches the observer at $\\chi = 0$.",
    ],
    ("de_sitter", "flat"): [
        "The whole de Sitter spacetime with the flat slicing on it. Its conformal "
        "time $\\eta = -e^{-Ht}/H$ makes the metric "
        "$\\frac{1}{H^2\\eta^2}(-c^2d\\eta^2 + d\\rho^2 + \\rho^2d\\Omega^2)$, conformal to half of "
        "Minkowski space, which enters the square as $p = \\pi/4 + \\arctan(H\\eta - H\\rho/c)$ and "
        "$q = \\pi/4 + \\arctan(H\\eta + H\\rho/c)$, with $\\rho^2 = x^2 + y^2 + z^2$.",
        "The slicing covers the half above the observer's past horizon, so its $t \\to -\\infty$ "
        "is a null line where the coordinates end, and the spacetime goes on below it.",
    ],
    ("anti_de_sitter", "global"): [
        "Anti-de Sitter spacetime, its universal cover, each point in the diagram a 2-sphere. "
        "With $\\sigma = \\arctan(r/L)$ the metric on the plane "
        "of $t$ and $r$ is $\\frac{1}{\\cos^2\\sigma}(-c^2dt^2 + L^2d\\sigma^2)$, already conformal "
        "to the strip $0 \\le \\sigma < \\pi/2$, which is unbounded in $t$. Its edge "
        "$\\sigma = \\pi/2$ is a timelike boundary.",
        "A radial light ray from the centre reaches the boundary at $ct = \\pi L/2$ and is back "
        "at $ct = \\pi L$. The radial timelike geodesics, $\\sin\\sigma = k\\sin(ct/L)$ with "
        "$k < 1$, all return to the centre at the same $ct = \\pi L$, whatever their energy.",
    ],
    ("einstein_static", "hyperspherical"): [
        "The Einstein static universe, each point in the diagram a 2-sphere of radius $R\\sin\\chi$. The "
        "metric on the plane of $t$ and $\\chi$ is $-c^2dt^2 + R^2d\\chi^2$, flat as it stands, so with the "
        "conformal time $\\eta = ct/R$ the maps $p = (\\eta - \\chi)/2$ and $q = (\\eta + \\chi)/2$ draw it as "
        "the strip $0 \\le \\chi \\le \\pi$, infinite in both directions of time, with light at 45°.",
        "The strip is the Einstein cylinder, $\\mathbb{R} \\times S^3$, with each of its 2-spheres about the "
        "pole drawn as a point. A light ray from the pole reaches the antipode after $\\pi R/c$ and returns "
        "to the pole every $2\\pi R/c$.",
    ],
    ("einstein_static", "areal"): [
        "The Einstein static universe in the strip of its hyperspherical chart, each point in the diagram "
        "a 2-sphere of radius $r = R\\sin\\chi$. The areal chart covers the half $\\chi < \\pi/2$ and ends "
        "at the equator $r = R$, where $g_{rr} = R^2/(R^2 - r^2)$ diverges and the curvature is the same as "
        "everywhere else.",
        "Each value of $r$ below $R$ names two spheres, one on either hemisphere, and the areal chart "
        "takes the one about the pole.",
    ],
    ("einstein_static", "einstein_cartesian"): [
        "The Einstein static universe in the strip of its hyperspherical chart, each point in the diagram "
        "a 2-sphere of radius $R\\sin\\chi$. Albert Einstein's coordinates of 1917 project the hemisphere "
        "$\\chi < \\pi/2$ onto its equatorial plane, with $x^2 + y^2 + z^2 = R^2\\sin^2\\chi$, and end at the "
        "equator.",
    ],
    ("einstein_static", "minkowski"): [
        "Minkowski space inside the Einstein static universe, each point in the diagram a 2-sphere. With "
        "$ct \\pm r = R\\tan((\\eta \\pm \\chi)/2)$ the metric of Minkowski space is "
        "$\\Omega^{-2}$ times the metric of the Einstein static universe, with "
        "$\\Omega = 2\\cos((\\eta + \\chi)/2)\\cos((\\eta - \\chi)/2)$, on the triangle $|\\eta| + \\chi < \\pi$.",
        "Null infinity $\\mathscr{I}^\\pm$ is the pair of light cones from the antipode's point $i^0$, and "
        "$i^\\pm$ lie on the pole at $\\eta = \\pm\\pi$.",
    ],
    ("einstein_static", "de_sitter"): [
        "De Sitter space inside the Einstein static universe, each point in the diagram a 2-sphere. "
        "De Sitter space of radius $\\ell = R$ is $\\ell^2/\\cos^2\\eta$ times the metric of the Einstein "
        "static universe on the band $-\\pi/2 < \\eta < \\pi/2$, which fills the whole three sphere at every "
        "$\\eta$.",
        "Its future and past infinity $\\mathscr{I}^\\pm$ are the spacelike lines $\\eta = \\pm\\pi/2$, and a "
        "light ray from the pole at $\\mathscr{I}^-$ reaches the antipode just at $\\mathscr{I}^+$.",
    ],
    ("einstein_static", "anti_de_sitter"): [
        "Anti-de Sitter space inside the Einstein static universe, each point in the diagram a 2-sphere. "
        "With $r = L\\tan\\chi$ and $L = R$, anti-de Sitter space is $1/\\cos^2\\chi$ times the metric of the "
        "Einstein static universe on the half $\\chi < \\pi/2$, the hemisphere about the pole at every "
        "moment.",
        "Its conformal boundary $\\mathscr{I}$ is the equator $\\chi = \\pi/2$, a timelike line of the strip, "
        "which light from the pole reaches after $\\pi R/2c$.",
    ],
    ("anti_de_sitter", "poincare"): [
        "The plane $x = y = 0$ of the Poincaré patch, through the centre and conformal to the "
        "strip $-\\pi/2 < \\sigma < \\pi/2$. The coordinates $t$ and $z$ "
        "enter the strip as "
        "$p = -\\pi/4 + \\arctan((ct + z)/L)$ and $q = \\pi/4 + \\arctan((ct - z)/L)$.",
        "The Poincaré coordinates cover a wedge of the strip. Their $z \\to 0$ is the conformal "
        "boundary, and $z \\to \\infty$ is the Poincaré horizon, the pair of null lines from "
        "$(\\sigma, ct/L) = (-\\pi/2, 0)$. The curvature there is the same as everywhere else, and "
        "the global coordinates run smoothly across it.",
    ],
    ("bertotti_robinson", "static"): [
        "The Bertotti-Robinson spacetime, the product of an anti-de Sitter space of two "
        "dimensions and radius $b$ with a sphere of radius $b$. Its conformal diagram is the "
        "strip of the first factor, each point in the strip a 2-sphere of radius $b$, the same "
        "radius everywhere.",
        "The throat coordinates $t$ and $r$ cover a Poincaré wedge of the strip. Under "
        "$x = b^2/r$ their $-\\frac{r^2}{b^2}c^2dt^2 + \\frac{b^2}{r^2}dr^2$ becomes "
        "$\\frac{b^2}{x^2}(-c^2dt^2 + dx^2)$, so the throat's $r = 0$ is the Poincaré horizon, a "
        "horizon of the coordinates alone, across which the spacetime continues.",
    ],
    ("bertotti_robinson", "poincare"): [
        "The whole Bertotti-Robinson spacetime with the Poincaré coordinates $t$ "
        "and $x$ on it. They cover the same wedge as the throat coordinates, with "
        "$x = b^2/r$.",
    ],
    ("ellis_bronnikov", "spherical"): [
        "The Ellis-Bronnikov wormhole, each point in the diagram a 2-sphere of area "
        "$4\\pi(r^2 + \\ell^2)$. The metric on the plane of $t$ and $r$ is "
        "$-c^2dt^2 + dr^2$ with $r$ over the whole line, which is Minkowski space of two "
        "dimensions, and $p, q = \\arctan((ct \\mp r)/\\ell)$ bring it into the full diamond.",
        "The two ends, $r \\to +\\infty$ and $r \\to -\\infty$, are two asymptotically flat "
        "universes, each with its own $i^0$ and $\\mathscr{I}^\\pm$, joined at the throat $r = 0$, "
        "where the spheres are smallest. Light crosses the throat at 45°, as it does everywhere "
        "else, so the wormhole has no horizon.",
    ],
    ("morris_thorne", "spherical"): [
        "The Morris-Thorne wormhole ($\\Phi = 0$, $b = b_0^2/r$), each point in the diagram a "
        "2-sphere of radius $r$. The proper radial distance "
        "is $l = \\pm\\sqrt{r^2 - b_0^2}$, and $p, q = \\arctan((ct \\mp l)/b_0)$ bring the plane "
        "of $t$ and $l$ into the full diamond, with the throat $r = b_0$ on its axis.",
        "Every $\\Phi$ and $b$ that give no horizon and two flat ends give the same diamond. With "
        "$\\Phi$ bounded and tending to a constant and $b/r \\to 0$ at both ends, "
        "$\\int e^{-\\Phi}\\,dl$ runs over the whole line, and another choice moves only the "
        "surfaces of constant $t$ and $r$ inside it. The areal coordinates $t$ and $r$ cover one "
        "side and end at the throat, where $g_{rr}$ diverges.",
    ],
    ("morris_thorne", "proper_radial"): [
        "The same wormhole in the proper distance coordinates $t$ and $l$, which cover "
        "both sides and run smoothly through the throat at $l = 0$. With $\\Phi = 0$ the metric "
        "on the plane is $-c^2dt^2 + dl^2$, and $p, q = \\arctan((ct \\mp l)/b_0)$ bring it into "
        "the diamond.",
    ],
    ("cosmic_string", "conical"): [
        "The half plane of $t$ and $r$ at fixed $\\phi$ and $z$, totally geodesic. The metric on it is $-c^2dt^2 + dr^2$, so "
        "$p, q = \\arctan((ct \\mp r)/\\ell)$ bring it into Minkowski's half diamond, with the "
        "string at $r = 0$ in place of a regular centre.",
        "The string's gravity is its deficit angle $\\delta = 8\\pi G\\mu/c^2$, which shows in "
        "the circles of constant $r$ around it: each has circumference $2\\pi(1 - 4G\\mu/c^2)\\,r$, "
        "short of $2\\pi r$ by $\\delta r$.",
    ],
    ("cosmic_string", "gott"): [
        "The same half plane with Gott's core, which makes the axis regular. Inside the "
        "core the proper distance from the axis is $\\rho = \\ell\\chi$, and outside it is "
        "$\\ell\\chi_0 + r - \\ell\\tan\\chi_0$; the edge of the core, $r = \\ell\\tan\\chi_0$, is "
        "where the circumferences inside and outside agree. In $\\rho$ the metric on the whole "
        "half plane is $-c^2dt^2 + d\\rho^2$.",
    ],
    ("interior_schwarzschild", "spherical"): [
        "A static star of uniform density, the interior solution for $r \\le R$ joined at "
        "$R = 1.5\\,r_s$ to the Schwarzschild exterior, each point in the diagram a 2-sphere of "
        "radius $r$. That radius clears Buchdahl's bound, "
        "$R > \\frac{9}{8}r_s$, so there is no horizon.",
        "At the surface $g_{tt} = -(1 - r_s/R)$ on both sides, so $t$ is one coordinate "
        "throughout. The tortoise coordinate $r_* = \\int\\sqrt{g_{rr}/(-g_{tt})}\\,dr$ runs from "
        "the centre through the surface, and $p, q = \\arctan((t \\mp r_*)/R)$ bring the spacetime "
        "into Minkowski's triangle, the causal structure of empty space, with the star a timelike "
        "tube from $i^-$ to $i^+$.",
    ],
    ("frw", "flat"): [
        "A flat universe of dust, each point in the diagram a 2-sphere. With $k = 0$, $G^r{}_r = 0$ gives $a \\propto \\eta^2$, and the metric "
        "$a^2(-d\\eta^2 + dr^2 + r^2d\\Omega^2)$ is conformal to the half $\\eta > 0$ of Minkowski "
        "space, which $p, q = \\arctan((\\eta \\mp r)/\\eta_0)$ bring into a triangle, with "
        "$\\eta_0$ the conformal time today.",
        "The big bang is the straight line $T = 0$ and is spacelike, while future infinity is "
        "null, as Minkowski's is. Our past light cone meets the bang at the comoving radius "
        "$r = \\eta_0$, the particle horizon, and light from anything farther away has not "
        "reached us yet.",
    ],
    ("frw", "closed"): [
        "A closed universe of dust, each point in the diagram a 2-sphere. With $k = +1$, dust gives $a \\propto 1 - \\cos\\eta$, with $\\eta$ from "
        "$0$ to $2\\pi$, and with $r = \\sin\\chi$ the metric is already conformal to the Einstein "
        "static universe, the rectangle $0 \\le \\chi \\le \\pi$ with the bang along its bottom "
        "and the crunch along its top.",
        "A light ray that leaves $\\chi = 0$ at the bang reaches the antipode $\\chi = \\pi$ at "
        "maximum expansion and is back at the crunch. The radius $r = \\sin\\chi$ covers one "
        "hemisphere, $\\chi < \\pi/2$, and ends at the equator $r = 1$, where $1 - kr^2$ "
        "vanishes.",
    ],
    ("frw", "open"): [
        "An open universe of dust, each point in the diagram a 2-sphere. With $k = -1$, dust gives $a \\propto \\cosh\\eta - 1$, and with "
        "$r = \\sinh\\chi$ the map $\\tan((T \\pm X)/2) = \\tanh((\\eta \\pm \\chi)/2)$ sends it into "
        "the Einstein static universe. It has the causal structure of the flat universe, a "
        "triangle with the bang along its base and null infinity above, and differs from it only "
        "in where its surfaces of constant $\\eta$ and $\\chi$ lie.",
    ],
    ("oppenheimer_snyder", "collapse"): [
        "A spherically symmetric distribution of dust collapsing from rest ($R_0 = 2\\,r_s$), "
        "each point in the diagram a 2-sphere. Inside, the dust is a closed universe, "
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
        "A static star of fluid with a polytrope for its equation of state, each point in the "
        "diagram a 2-sphere of radius $r$. Its mass and redshift "
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
    ("khan_penrose", "double_null"): [
        "The plane $x = y = 0$ of the Khan-Penrose spacetime, totally geodesic, each point in the diagram a single "
        "event. Two impulsive plane waves travel toward each other through flat space and collide at $u = v = 0$, and "
        "on their fronts, $u = 0$ and $v = 0$, the Riemann tensor has a delta singularity. Maps $p$ and $q$ of $u$ and "
        "$v$, each the coordinate itself where it is positive and its arctangent where it is negative, bring the "
        "whole plane into a finite drawing with $T = p + q$ up and $X = q - p$ across, and light at 45°.",
        "Ahead of both waves the spacetime is flat, and its past edges are null infinity, $u \\to -\\infty$ and "
        "$v \\to -\\infty$. Behind one wave alone it is flat again, $-2L^2du\\,dv + (1 + u)^2dx^2 + (1 - u)^2dy^2$ "
        "behind the wave on $u = 0$, and ends at $u = 1$ in a fold singularity, where the curvature is zero, as the "
        "region behind the other wave ends at $v = 1$. Where both waves have passed each focuses the other, and the "
        "region ends on the spacelike curvature singularity $u^2 + v^2 = 1$, where the Kretschmann scalar diverges.",
        "The coordinates $u$ and $v$ cover the region where both waves have passed, $0 \\le v < \\sqrt{1 - u^2}$ "
        "with $0 \\le u < 1$, and ahead of either wave the metric is the same with that wave's coordinate set to "
        "zero.",
    ],
    ("khan_penrose", "cosmological"): [
        "The plane $x = y = 0$ of the Khan-Penrose spacetime, totally geodesic, each point in the diagram a single "
        "event. Two impulsive plane waves travel toward each other through flat space and collide at $u = v = 0$, and "
        "on their fronts, $u = 0$ and $v = 0$, the Riemann tensor has a delta singularity. Maps $p$ and $q$ of $u$ and "
        "$v$, each the coordinate itself where it is positive and its arctangent where it is negative, bring the "
        "whole plane into a finite drawing with $T = p + q$ up and $X = q - p$ across, and light at 45°.",
        "Ahead of both waves the spacetime is flat, and its past edges are null infinity, $u \\to -\\infty$ and "
        "$v \\to -\\infty$. Behind one wave alone it is flat again, $-2L^2du\\,dv + (1 + u)^2dx^2 + (1 - u)^2dy^2$ "
        "behind the wave on $u = 0$, and ends at $u = 1$ in a fold singularity, where the curvature is zero, as the "
        "region behind the other wave ends at $v = 1$. Where both waves have passed each focuses the other, and the "
        "region ends on the spacelike curvature singularity $u^2 + v^2 = 1$, where the Kretschmann scalar diverges.",
        "The coordinates $\\tau = \\arcsin u + \\arcsin v$ and $\\sigma = \\arcsin u - \\arcsin v$ cover the "
        "region where both waves have passed, $|\\sigma| \\le \\tau < \\pi/2$, each surface of constant $\\tau$ "
        "spacelike, and the curvature singularity is $\\tau = \\pi/2$.",
    ],
    ("einstein_rosen_waves", "cylindrical"): [
        "The half plane of $t$ and $\\rho$ at fixed $\\phi$ and $z$, totally geodesic. The metric on it is "
        "$e^{2(\\gamma - \\psi)}(-c^2dt^2 + d\\rho^2)$, and a conformal factor changes no null direction, so for "
        "every wave $p, q = \\arctan((ct \\mp \\rho)/a)$ bring it into Minkowski's half diamond, with the axis "
        "$\\rho = 0$ in place of a regular centre wherever $\\gamma = 0$ there.",
        "The pulse of Weber, Wheeler, and Bonnor comes in from $\\mathscr{I}^-$, is greatest on the axis at "
        "$t = 0$, and goes out to $\\mathscr{I}^+$, its crest just outside the rays $\\rho = |ct|$.",
    ],
    ("einstein_rosen_waves", "null"): [
        "The same half plane in the null chart ($u = ct - \\rho$, $v = ct + \\rho$), where the metric on it is "
        "$-e^{2(\\gamma - \\psi)}du\\,dv$. The lines of constant $u$ and of constant $v$ are light rays, the 45° "
        "lines of the triangle, with $p = \\arctan(u/a)$ and $q = \\arctan(v/a)$, and the axis is the line $u = v$.",
    ],
    ("melvin", "cylindrical"): [
        "The half plane of $t$ and $\\rho$ of Melvin's universe at fixed $\\phi$ and $z$, totally geodesic. The metric "
        "on it is $(1 + B^2\\rho^2/4)^2(-c^2dt^2 + d\\rho^2)$, and a conformal factor changes no null direction, so "
        "$p, q = \\arctan((ct \\mp \\rho)B)$ bring it into Minkowski's half diamond, with the regular axis $\\rho = 0$ "
        "on its left edge.",
        "A light ray reaches $\\rho \\to \\infty$ only at an infinite value of its affine parameter, which grows as "
        "$\\int(1 + B^2\\rho^2/4)^2\\,d\\rho$, so the far edges are the null infinity of the half plane, although the "
        "circles about the axis shrink to zero there.",
    ],
    ("melvin", "ernst"): [
        "The plane of $t$ and $r$ of Ernst's black hole at $\\theta = \\pi/2$, totally geodesic. The metric on it is "
        "$(1 + B^2r^2/4)^2$ times Schwarzschild's, so Kruskal and Szekeres's extension draws it for every $B$: two "
        "exteriors, the black hole above, the white hole below, and the curvature singularity $r = 0$ at the top "
        "and the bottom.",
        "The plane of $t$ and $r$ at every other constant $\\theta$ has the same null curves and the same drawing, "
        "on the axis with the factor equal to $1$.",
    ],
    ("malament_hogarth", "cartesian"): [
        "The Malament-Hogarth toy spacetime, Minkowski space with one event "
        "removed and its metric multiplied by $\\Omega^2$. A conformal factor changes no null "
        "direction, so for every $\\Omega$ the causal structure is Minkowski's less that event. "
        "It is symmetric about the $t$ axis through the event, so each point in the triangle is "
        "a 2-sphere of events at one $t$ and one distance $r = \\sqrt{x^2 + y^2 + z^2}$ from the "
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
        "A spacetime into which a spherical shell of null dust of mass $M$ falls along $v = 0$, "
        "each point in the diagram a 2-sphere. With $m = 0$ for "
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
    views = DRAWN[metric_id](ck, src)
    out = []
    for v in views:
        v.set(caption=CAPTIONS[(metric_id, v.d["id"])])
        out.append(v.done())
    return {"metric": metric_id, "source": src.stamps(), "views": out}


def check_table():
    """Every view has a caption, and every spacetime is either drawn or named as not drawn."""
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
    unknown = set(args.metric) - set(DRAWN)
    if unknown:
        parser.error(f"no conformal diagram is drawn for {sorted(unknown)}")
    wanted = [m for m in DRAWN if not args.metric or m in args.metric]
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
    if not args.metric:
        # A spacetime that is not drawn has no file, so a full redraw removes one left behind.
        for path in sorted(CONFORMAL_DIR.glob("*.json")):
            if path.stem not in DRAWN and not build.CONFLICT_COPY.search(path.stem):
                path.unlink()
                print(f"removed {path.relative_to(build.ROOT)}")
    for metric_id, data in files.items():
        path = CONFORMAL_DIR / f"{metric_id}.json"
        path.write_text(json.dumps(data, ensure_ascii=False, separators=(",", ":")) + "\n", encoding="utf-8")
        print(f"wrote {path.relative_to(build.ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
