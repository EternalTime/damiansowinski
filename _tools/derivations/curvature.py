#!/usr/bin/env python3
"""Curvature figures of the spacetimes, drawn from the published Kretschmann and Ricci scalars.

A curvature figure plots the two scalars a chart publishes, `kretschmann` and `ricci_scalar`,
as numbers, read from the published LaTeX through the Reader of verify_metrics.py beside this
file, in the x^0 = cT chart the collection writes, with c = 1, by the same `load` the null ray
diagrams use. Nothing in a figure is typed by hand: every value is the published formula
evaluated at the point it is drawn at, and every horizon marked is a zero of the published
inverse metric.

    _tools/setup-env.sh
    .venv.noindex/bin/python _tools/derivations/curvature.py
    python3 _tools/build_mfs_data.py

writes MFS/assets/data/curvature/<metric_id>.json for every spacetime in DRAWN, and the second
command stamps each file's version into the index. Pass --metric <metric_id> to redraw one
spacetime, repeatable, and --verify to run every check and print it without writing anything.
_tools/README.md, under "Curvature figures", defines the file.


The three shapes
----------------

line   one coordinate across, and one or more panels stacked over it, each a scalar on an axis
       of its own: the spacetime's own curve, any comparison curve read from another
       spacetime's published scalar, and vertical lines at the horizons, the zeros of the
       published g^xx along the drawn coordinate.
map    two coordinates, at one scale on both axes, and the scalar as a colour at the centre of
       each cell of an even grid over the box, with points and curves marked on it; in the frame
       "meridian" the box is the meridian plane of a radius and an angle from the axis.
movie  either of those stepped through a parameter or a time, a frame for each value, as an
       embedding diagram's movie is stepped.

Each axis of values, a panel's or a map's colour, is linear, logarithmic, or signed
logarithmic, the position of v along it being v, log10 v, or asinh(v/2L)/ln 10 for the
scale's `linear` width L, which is linear within |v| < L and log10|v/L| times the sign of v
far outside it. A value beyond the range is cut at the edge, and the view's `clip` says where.


What is checked
---------------

write() refuses to write while any of these fails:

  values    every value written is the published formula at its point: each is evaluated in
            floating point and a sample of them again in sixty digits through sympy, which
            must agree to a part in 10^9, and every value written is finite, and positive on a
            logarithmic axis;
  depends   the published scalar depends on nothing but the coordinates drawn and the
            parameters set, so the figure is the same at every value of the others;
  marks     every horizon marked is a sign change of the published g^xx, found by
            bisection, and every point marked is checked as its spacetime's function says;
  shape     every panel and the plot of every map is between 1:2 and 2:1, as the spacetime
            diagrams are, and a stack of panels is too;
  drawn     no polyline written strays from the curve sampled at 4001 points by more than
            TOLERANCE of the 520 units a plot is drawn across;
  physics   each spacetime's function checks what its caption states, as that function's
            docstring says.

_tools/test_curvature.py holds every value on disk to the published formula again, from the
numbers in the files and the metric files alone.


Which spacetimes
----------------

DRAWN lists the spacetimes that have a curvature figure, each with the function that draws
it, and a full redraw removes the file of a spacetime it does not list. A spacetime earns one
where its curvature carries physics the formula alone does not show at a glance, as the
research list of October 2026 sets out; the others have no file and the page no section.
"""

import argparse
import json
import math
import sys
import time
from fractions import Fraction
from pathlib import Path

import numpy as np
import sympy as sp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))

import build_mfs_data as build  # noqa: E402
import null_rays as nr  # noqa: E402

CURVATURE_DIR = build.CURVATURE_DIR

# The width a plot is drawn across on the page, in the units its shapes are measured in, and how
# far a polyline written may stray from the curve it draws, in those units.
PLOT_W = 520
TOLERANCE = 0.05
# How many points a curve is sampled at before it is thinned.
SAMPLES = 4001
# The shapes a panel or a map's plot may take, width over height, as the spacetime diagrams'.
WIDEST, TALLEST = build.ASPECT_WIDEST, build.ASPECT_TALLEST
# How closely a value written agrees with the published formula evaluated in sixty digits.
AGREE = 1e-9
SCALARS = {"kretschmann": "K", "ricci_scalar": "R"}


class Sources:
    """The published fields every figure read, per coordinate system."""

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


class Checks:
    """Every check made, by name, with whether it held and what was measured."""

    def __init__(self):
        self.made = []

    def __call__(self, name, ok, detail=""):
        self.made.append((name, bool(ok), detail))
        return bool(ok)

    def failures(self):
        return [f"{name}: {detail}" for name, ok, detail in self.made if not ok]

    def report(self):
        for name, ok, detail in self.made:
            print(f"{'ok    ' if ok else 'FAILED'} {name}  {detail}")


# ---------------------------------------------------------------- reading a published function

class Published:
    """A published value of one chart, a curvature scalar, a component of the metric or of its
    inverse, one of the Riemann tensor with its first index raised, or one of the Einstein tensor
    with its first index raised, as a numpy function of the coordinates drawn and then the
    parameters, in the order named.

    The value is read through the checker's Reader, every name the reader holds as a function
    written out, every function the chart leaves free that `functions` declares put in as its
    closed form, written as null_rays.py's rows write theirs, and c set to 1. A coordinate or a
    parameter the value depends on and the figure does not name stops the script, since the
    figure would then hold only at one value of it, and so does a free function left undeclared."""

    def __init__(self, sources, metric_id, system_id, field, coords, params, indices=None, functions=None):
        _, entry, reader = nr.load(metric_id, system_id)
        self.metric, self.system, self.field = metric_id, system_id, field
        self.functions = dict(functions or {})
        if field in ("metric_components", "inverse_metric_components", "riemann", "einstein_tensor"):
            variant = {"riemann": "ulll", "einstein_tensor": "ul"}.get(field)
            listed = entry[field]["variants"][variant]["nonzero"] if variant else entry[field]
            value = next(c["value"] for c in listed if c["indices"] == list(indices))
            expr = reader(value)
        else:
            expr = reader(nr.strip_lhs(entry[field]))
        sources.note(metric_id, system_id, ["coords", "parameters", field])
        if reader.held:
            expr = expr.subs(reader.held).doit()
        expr = declared(reader, expr, self.functions)
        expr = expr.subs(reader.c, 1)
        by_plain = {reader._plain(name): symbol for name, symbol in reader.symbol.items()}
        self.coords = [by_plain[name] for name in coords]
        self.params = [reader.parameters[name] for name in params]
        free = expr.atoms(sp.Function) & {f for f in reader.parameters.values() if isinstance(f, sp.Function)}
        if free:
            raise SystemExit(f"{metric_id}/{system_id}: the published {field} holds the free functions "
                             f"{sorted(map(str, free))}, which the figure does not declare")
        stray = expr.free_symbols - set(self.coords) - set(self.params)
        if stray:
            raise SystemExit(f"{metric_id}/{system_id}: the published {field} depends on "
                             f"{sorted(map(str, stray))}, which the figure neither draws nor sets")
        # Exponentials are gathered into one, so that e^(-4 gamma) e^(4 psi) overflows to infinity
        # beside a ring and never to zero times infinity.
        self.expr = sp.powsimp(expr, combine="exp")
        self.f = sp.lambdify(self.coords + self.params, self.expr, nr.numeric_modules(self.expr))

    def __call__(self, *values):
        with np.errstate(all="ignore"):
            out = self.f(*values)
        shape = np.broadcast(*[np.asarray(v, dtype=float) for v in values]).shape
        return np.array(np.broadcast_to(np.asarray(out, dtype=float), shape))

    def exact(self, values, digits=60):
        """The value at one point in `digits` digits, the coordinates and parameters taken as the
        doubles given, exactly."""
        point = {s: sp.Float(float(v), digits + 10) for s, v in zip(self.coords + self.params, values)}
        return self.expr.xreplace(point).evalf(digits)


class Held:
    """A published value of a chart in more coordinates than a curve draws, with the others held
    at the numbers `fixed` gives them, by the names the chart spells them with, as a function of
    the coordinate drawn and then the parameters, as a curve along the axis of a chart in r and
    theta is drawn at theta = 0. The published value is read with the coordinate drawn first and
    the held ones after it, in the order `fixed` names them."""

    def __init__(self, published, fixed):
        self.published, self.fixed = published, {name: float(v) for name, v in fixed.items()}
        self.metric, self.system, self.field = published.metric, published.system, published.field
        self.functions = published.functions

    def __call__(self, x, *params):
        return self.published(x, *self.fixed.values(), *params)

    def exact(self, values, digits=60):
        return self.published.exact([values[0], *self.fixed.values(), *values[1:]], digits)


def declared(reader, expr, functions):
    """A published expression with each free function of its chart that `functions` names put in
    as the closed form written there, a sympy expression in the chart's coordinates and parameters,
    as conformal.py's planes put in theirs."""
    names = {reader._plain(n): s for n, s in reader.symbol.items()}
    names.update(reader.parameters)
    for name, text in functions.items():
        if name not in reader.parameters or not isinstance(reader.parameters[name], sp.Function):
            raise SystemExit(f"{name!r} is no function the chart leaves free")
        rep = sp.sympify(text, locals={**nr.DECLARED_FUNCTIONS, **names})
        expr = expr.subs(reader.parameters[name], rep).doit()
    return expr


def agrees(ck, name, published, points, values):
    """Check `values`, written at `points`, against the published formula in sixty digits."""
    worst = 0.0
    for point, value in zip(points, values):
        exact = published.exact(point)
        scale = max(abs(float(exact)), 1e-300)
        worst = max(worst, abs(float(exact) - value) / scale)
    return ck(f"{name} agrees with the published {published.field} in sixty digits",
              worst < AGREE, f"worst {worst:.1e} at {len(points)} points")


# ---------------------------------------------------------------- scales and ticks

def position(scale, v):
    """Where v stands along an axis of values, in the axis's own units."""
    v = np.asarray(v, dtype=float)
    with np.errstate(all="ignore"):
        if scale["kind"] == "log":
            return np.where(v > 0, np.log10(np.where(v > 0, v, 1.0)), np.nan)
        if scale["kind"] == "symlog":
            return np.arcsinh(v / (2 * scale["linear"])) / math.log(10)
    return v


def span(scale):
    return float(position(scale, scale["from"])), float(position(scale, scale["to"]))


def decade(k, sign=1):
    return f"${'-' if sign < 0 else ''}10^{{{k}}}$"


def value_ticks(scale, target=6):
    """The ticks of an axis of values, each at a value, with its label as TeX: steps of 1, 2, 2.5
    or 5 times a power of ten on a linear axis, the decades on a logarithmic one, and on a signed
    one zero and the decades from ten times its linear width out, of both signs."""
    lo, hi = scale["from"], scale["to"]
    if scale["kind"] == "linear":
        return nr.ticks(lo, hi, target)
    if scale["kind"] == "log":
        return [{"at": float(10.0 ** k), "label": decade(k)}
                for k in range(math.ceil(math.log10(lo) - 1e-12), math.floor(math.log10(hi) + 1e-12) + 1)]
    first = math.ceil(math.log10(10 * scale["linear"]) - 1e-12)
    out = []
    if lo < 0:
        out += [{"at": -float(10.0 ** k), "label": decade(k, -1)}
                for k in range(math.floor(math.log10(-lo) + 1e-12), first - 1, -1)]
    if lo <= 0 <= hi:
        out.append({"at": 0.0, "label": "$0$"})
    if hi > 0:
        out += [{"at": float(10.0 ** k), "label": decade(k)}
                for k in range(first, math.floor(math.log10(hi) + 1e-12) + 1)]
    return out


def axis(name, lo, hi, target=6):
    return {"name": name, "from": lo, "to": hi, "ticks": nr.ticks(lo, hi, target)}


def tex_number(v, digits=3):
    """A number as the captions write it, to `digits` significant figures, without a trailing
    point."""
    text = f"{v:.{digits}g}"
    if "e" in text:
        mantissa, exponent = text.split("e")
        return f"{mantissa} \\times 10^{{{int(exponent)}}}"
    return text


# ---------------------------------------------------------------- sampling a curve

def thinned(X, Y, tolerance):
    """The indices Ramer, Douglas and Peucker keep of the polyline (X, Y): every point dropped
    lies within `tolerance` of the chord that replaces it."""
    keep = np.zeros(len(X), dtype=bool)
    keep[0] = keep[-1] = True
    stack = [(0, len(X) - 1)]
    while stack:
        i, j = stack.pop()
        if j <= i + 1:
            continue
        dx, dy = X[j] - X[i], Y[j] - Y[i]
        length = math.hypot(dx, dy)
        k = np.arange(i + 1, j)
        if length == 0:
            d = np.hypot(X[k] - X[i], Y[k] - Y[i])
        else:
            d = np.abs(dy * (X[k] - X[i]) - dx * (Y[k] - Y[i])) / length
        worst = int(np.argmax(d))
        if d[worst] > tolerance:
            keep[k[worst]] = True
            stack += [(i, k[worst]), (k[worst], j)]
    return np.flatnonzero(keep)


def runs(mask):
    """The runs of True in a boolean array, as (start, stop) index pairs, stop exclusive."""
    edges = np.flatnonzero(np.diff(np.concatenate([[0], mask.astype(int), [0]])))
    return list(zip(edges[::2], edges[1::2]))


def sampled(ck, name, published, x_axis, scale, height, params, within=None):
    """A published scalar along the drawn coordinate, as the polylines a panel draws.

    The scalar is sampled at SAMPLES points across the axis, or across the part of it `within`
    names where the scalar's own chart or spacetime covers only that part, as a star's own chart
    holds from its centre to its surface, and kept where it is finite and within a quarter of the
    panel's height of its box, each stretch carried one sample further each way where that sample
    is finite, so that the line drawn runs to the edge it leaves by.
    Each stretch is thinned in the panel's own units, PLOT_W across and `height` up, so that the
    line drawn lies within TOLERANCE of the samples."""
    lo, hi = x_axis["from"], x_axis["to"]
    xs = np.linspace(*(within or (lo, hi)), SAMPLES)
    vs = published(xs, *params)
    p0, p1 = span(scale)
    pos = position(scale, vs)
    X = PLOT_W * (xs - lo) / (hi - lo)
    Y = height * (1 - (pos - p0) / (p1 - p0))
    finite = np.isfinite(vs) & np.isfinite(Y)
    near = finite & (Y >= -0.25 * height) & (Y <= 1.25 * height)
    out = []
    for start, stop in runs(near):
        if start > 0 and finite[start - 1]:
            start -= 1
        if stop < len(xs) and finite[stop]:
            stop += 1
        if stop - start < 2:
            continue
        keep = start + thinned(X[start:stop], Y[start:stop], TOLERANCE)
        # Adding zero writes a zero as 0.0, never as -0.0.
        out.append([[float(x) + 0.0, float(v) + 0.0] for x, v in zip(xs[keep], vs[keep])])
    ck(f"{name} has a curve to draw", out, f"{len(out)} stretches")
    if scale["kind"] == "log":
        ck(f"{name} is positive where a logarithmic axis draws it",
           all(v > 0 for line in out for _, v in line))
    return out


def beyond(published, x_axis, scale, params, within=None):
    """Where a scalar lies above its axis and where below it, as the intervals of the drawn
    coordinate, from the same samples sampled() takes."""
    xs = np.linspace(*(within or (x_axis["from"], x_axis["to"])), SAMPLES)
    vs = published(xs, *params)
    pos, (p0, p1) = position(scale, vs), span(scale)
    above = ~np.isfinite(vs) | (pos > p1)
    below = np.isfinite(vs) & (pos < p0)
    return ([(xs[a], xs[b - 1]) for a, b in runs(above)], [(xs[a], xs[b - 1]) for a, b in runs(below)])


def zeros(ck, name, published, lo, hi, params, n=SAMPLES):
    """The zeros of a published function of the drawn coordinate on (lo, hi), each a change of
    sign between two samples found by bisection to the last digit, as the null ray diagrams find
    their horizons."""
    xs = np.linspace(lo, hi, n)
    vs = published(xs, *params)
    out = []
    for i in range(n - 1):
        a, b, fa, fb = xs[i], xs[i + 1], vs[i], vs[i + 1]
        if not (np.isfinite(fa) and np.isfinite(fb)) or fa == 0 or np.sign(fa) == np.sign(fb):
            continue
        for _ in range(200):
            m = 0.5 * (a + b)
            fm = float(published(m, *params))
            if m in (a, b) or fm == 0:
                break
            a, b, fa = (m, b, fm) if np.sign(fm) == np.sign(fa) else (a, m, fa)
        out.append(float(0.5 * (a + b)))
    ck(f"{name} has its zeros found", True, ", ".join(f"{z:.6f}" for z in out))
    return out


def meridian(X, Y):
    """The radius r and the angle theta from the axis of the point (X, Y) of a meridian plane,
    drawn across as r sin theta and up as r cos theta."""
    return np.hypot(X, Y), np.arctan2(X, Y)


def cells(published, x_axis, y_axis, nx, ny, params, frame=None):
    """A published scalar at the centres of an even grid of nx by ny cells over a map's box, a
    row for each cell up, from the bottom, each with a value for each cell across. In the frame
    "meridian" the box is the meridian plane of the chart's radius and polar angle, and each
    cell is read at the r and theta of its centre."""
    xs = x_axis["from"] + (np.arange(nx) + 0.5) * (x_axis["to"] - x_axis["from"]) / nx
    ys = y_axis["from"] + (np.arange(ny) + 0.5) * (y_axis["to"] - y_axis["from"]) / ny
    X, Y = np.meshgrid(xs, ys)
    if frame == "meridian":
        X, Y = meridian(X, Y)
    return xs, ys, published(X, Y, *params)


def rounded(v, digits=6):
    """A value written to `digits` significant figures, which the page and the checks read as
    the value at its point."""
    return float(f"{v:.{digits}g}") + 0.0


def shape_checks(ck, where, aspects):
    for k, a in enumerate(aspects):
        ck(f"{where} panel {k} is between 1:2 and 2:1", TALLEST <= a <= WIDEST, f"{a:.3f}")
    whole = 1 / sum(1 / a for a in aspects)
    ck(f"{where} is between 1:2 and 2:1 as a whole", TALLEST <= whole <= WIDEST, f"{whole:.3f}")


def curve_entry(panel, cls, published, coordinate, params, points):
    """A curve as the file holds it: what it is drawn from, so that every value can be evaluated
    again, the chart's other coordinates it holds fixed among it, and its points [x, value]."""
    entry = {"panel": panel, "class": cls, "metric": published.metric, "system": published.system,
             "scalar": published.field, "coordinate": coordinate}
    if isinstance(published, Held):
        entry["fixed"] = dict(published.fixed)
    entry["params"] = {name: float(value) for name, value in params.items()}
    if published.functions:
        entry["functions"] = dict(published.functions)
    entry["points"] = points
    return entry


def spread(points, n=24):
    """About n of a curve's points, evenly through it, its two ends among them."""
    flat = [p for line in points for p in line]
    picks = sorted({round(k * (len(flat) - 1) / (n - 1)) for k in range(n)})
    return [flat[k] for k in picks]


def line_curves(ck, name, panel, cls, published, coordinate, x_axis, scale, height, params, within=None):
    """The curves one published scalar draws on a panel, each checked in sixty digits."""
    values = [params[p] for p in params]
    lines = sampled(ck, f"{name} {cls}", published, x_axis, scale, height, values, within)
    picked = spread(lines)
    agrees(ck, f"{name} {cls}", published, [[x] + values for x, _ in picked], [v for _, v in picked])
    return [curve_entry(panel, cls, published, coordinate, params, line) for line in lines]


# ---------------------------------------------------------------- Bardeen's regular black hole

BARDEEN_G = Fraction(1, 3)        # g in units of r_s, as every drawing of Bardeen's black hole takes it
BARDEEN_STEP = Fraction(1, 60)    # the step of g from one frame of the movie to the next
BARDEEN_LAST = Fraction(1, 2)     # the last frame's g


def bardeen(ck, src):
    """Bardeen's regular black hole at r_s = 1, its Kretschmann and Ricci scalars against r.

    The static chart's published K and R depend on r alone, and every chart publishes the same
    two, so one figure serves all three charts. Two panels share the axis of r, from 0 to 2 r_s,
    the box of the spacetime diagrams: K on a logarithmic axis from 10^-1 to 10^5, which holds
    Bardeen's curve whole, its centre 24 r_s^2/g^6 = 17496/r_s^4 and its least value near
    0.14/r_s^4 at 2 r_s, and R on a signed one, linear within 0.01/r_s^2 of zero, from -10 to
    10^4, since R changes sign at r = 2g. Schwarzschild's two scalars of the same mass, read
    from schwarzschild.json, are the comparison, and the horizons are the zeros of the static
    chart's published g^rr. The first view is the black hole at g = r_s/3, the charge every
    other drawing of it takes; the second steps g from 0, Schwarzschild's black hole, to r_s/2,
    through the extremal 2 r_s/3 sqrt 3, where the two horizons meet.

    What the captions state is checked here: K at the centre is 24 r_s^2/g^6 and K falls from
    the centre outward; R at the centre is 12 r_s/g^3, vanishes at r = 2g, between the horizons,
    and is least just outside r_+, at the radius and value stated; the published scalars at
    g = 0 are Schwarzschild's; and past the extremal charge g^rr is positive at every r drawn."""
    x_axis = axis("$r/r_s$", 0, 2)
    K_scale = {"kind": "log", "from": 0.1, "to": 1e5}
    R_scale = {"kind": "symlog", "from": -10, "to": 1e4, "linear": 0.01}
    panels = [{"id": "K", "name": "$K\\,r_s^4$", "aspect": 2, "scale": K_scale, "ticks": value_ticks(K_scale)},
              {"id": "R", "name": "$R\\,r_s^2$", "aspect": 2, "scale": R_scale, "ticks": value_ticks(R_scale)}]
    shape_checks(ck, "bardeen", [p["aspect"] for p in panels])
    height = {p["id"]: PLOT_W / p["aspect"] for p in panels}
    K = nr_published(src, "bardeen", "static", "kretschmann", ["r"], ["r_s", "g"])
    R = nr_published(src, "bardeen", "static", "ricci_scalar", ["r"], ["r_s", "g"])
    grr = nr_published(src, "bardeen", "static", "inverse_metric_components", ["r"], ["r_s", "g"], ("r", "r"))
    K_s = nr_published(src, "schwarzschild", "spherical", "kretschmann", ["r"], ["r_s"])
    R_s = nr_published(src, "schwarzschild", "spherical", "ricci_scalar", ["r"], ["r_s"])

    def references():
        return (line_curves(ck, "bardeen Schwarzschild", "K", "reference", K_s, "r", x_axis, K_scale, height["K"],
                            {"r_s": 1})
                + line_curves(ck, "bardeen Schwarzschild", "R", "reference", R_s, "r", x_axis, R_scale, height["R"],
                              {"r_s": 1}))

    def own(g, name):
        params = {"r_s": 1, "g": float(g)}
        curves = (line_curves(ck, name, "K", "K", K, "r", x_axis, K_scale, height["K"], params)
                  + line_curves(ck, name, "R", "R", R, "r", x_axis, R_scale, height["R"], params))
        marks = [{"class": "horizon", "x": x} for x in zeros(ck, f"{name} horizons", grr, 0, 2, [1, float(g)])]
        return curves, marks

    g = float(BARDEEN_G)
    curves, marks = own(BARDEEN_G, "bardeen at g = r_s/3")
    ck("bardeen at g = r_s/3 has its two horizons", len(marks) == 2, f"{[m['x'] for m in marks]}")
    r_minus, r_plus = (m["x"] for m in marks)
    # The figures of the captions.
    centre_K, centre_R = float(K.exact([0, 1, g])), float(R.exact([0, 1, g]))
    ck("bardeen K at the centre is 24 r_s^2/g^6", abs(centre_K / (24 / g ** 6) - 1) < 1e-12, f"{centre_K}")
    ck("bardeen R at the centre is 12 r_s/g^3", abs(centre_R / (12 / g ** 3) - 1) < 1e-12, f"{centre_R}")
    rs = np.linspace(0, 2, SAMPLES)
    ck("bardeen K falls from the centre outward", np.all(np.diff(K(rs, 1, g)) < 0))
    sign_change = zeros(ck, "bardeen R", R, 0, 2, [1, g])
    ck("bardeen R vanishes at r = 2g alone, between the horizons",
       len(sign_change) == 1 and abs(sign_change[0] - 2 * g) < 1e-12 and r_minus < 2 * g < r_plus, f"{sign_change}")
    fine = np.linspace(r_plus, 2, 200001)
    least = int(np.argmin(R(fine, 1, g)))
    r_least, R_least = float(fine[least]), float(R(fine[least], 1, g))
    ck("bardeen R is least just outside r_+", r_plus < r_least < 1, f"r = {r_least:.4f}, R = {R_least:.4f}")
    ck("bardeen R is least for all r drawn there", R_least <= float(np.min(R(rs, 1, g))) + 1e-12)
    far = 1e3
    ck("bardeen K approaches Schwarzschild's far out", abs(float(K.exact([far, 1, g])) * far ** 6 / 12 - 1) < 1e-5)
    ck("bardeen R falls off as -3 r_s g^2/r^5 far out", abs(float(R.exact([far, 1, g])) * far ** 5 / (-3 * g ** 2) - 1) < 1e-5)
    still_refs = references()
    above_K, _ = beyond(K_s, x_axis, K_scale, [1])
    ck("Schwarzschild's K leaves the top of its panel once, at the centre",
       len(above_K) == 1 and above_K[0][0] == 0, f"{above_K}")
    for name, fn, scale in (("K", K, K_scale), ("R", R, R_scale)):
        up, down = beyond(fn, x_axis, scale, [1, g])
        ck(f"bardeen {name} at g = r_s/3 stays inside its panel", not up and not down, f"{up} {down}")
    cut = above_K[0][1]
    clip = ("We draw $K\\,r_s^4$ on a logarithmic scale from $10^{-1}$ to $10^{5}$ and $R\\,r_s^2$ on a signed "
            "logarithmic scale, linear within $0.01$ of zero, from $-10$ to $10^{4}$, and cut each curve where it "
            f"leaves its panel, as Schwarzschild's $K$ does inside $r = {cut:.3f}\\,r_s$.")
    legend = [["line", "K", "$K$, the Kretschmann scalar"], ["line", "R", "$R$, the Ricci scalar"],
              ["line", "reference", "Schwarzschild's black hole of the same mass"],
              ["line", "horizon", "a horizon, where $g^{rr} = 0$"]]
    still = {
        "id": "charge", "label": "$g = r_s/3$", "kind": "line", "x": x_axis, "panels": panels,
        "curves": still_refs + curves, "marks": marks, "legend": legend,
        "caption": [
            "The Kretschmann scalar $K$ and the Ricci scalar $R$ of Bardeen's black hole against $r$ "
            "($g = r_s/3$), the same at every $t$, $\\theta$, and $\\phi$, beside Schwarzschild's of the same mass. "
            "Toward the centre Schwarzschild's $K = 12r_s^2/r^6$ grows without bound, and Bardeen's rises to a "
            f"finite $24r_s^2/g^6 = {centre_K:.0f}/r_s^4$ at $r = 0$, the curvature of de Sitter space of radius "
            f"$\\sqrt{{g^3/r_s}} = {math.sqrt(g ** 3):.3f}\\,r_s$.",
            f"At the centre $R = 12r_s/g^3 = {centre_R:.0f}/r_s^2$. It changes sign at $r = 2g$, between the two "
            f"horizons, and is least, ${R_least:.3f}/r_s^2$, at $r = {r_least:.3f}\\,r_s$, just outside $r_+$. Far "
            "out $K$ approaches Schwarzschild's and $R$ falls off as $-3r_sg^2/r^5$.",
        ],
        "settings": f"$r_s = 1$, the unit of every length, and $g = r_s/3$, so that $r_+ = {r_plus:.3f}\\,r_s$ and "
                    f"$r_- = {r_minus:.3f}\\,r_s$.",
        "clip": clip,
    }

    # The movie in g, at fixed mass, from Schwarzschild's black hole to past the extremal charge.
    frames, horizons = [], []
    steps = int(BARDEEN_LAST / BARDEEN_STEP)
    for k in range(steps + 1):
        value = BARDEEN_STEP * k
        curves, marks = own(value, f"bardeen at g = {float(value):.4f}")
        frames.append({"label": f"$g = {float(value):.3f}\\,r_s$", "value": float(value),
                       "curves": curves, "marks": marks})
        horizons.append(len(marks))
    extremal = 2 / (3 * math.sqrt(3))
    gaps = [f["marks"][1]["x"] - f["marks"][0]["x"] for f in frames if len(f["marks"]) == 2]
    ck("bardeen's two horizons approach each other as g grows", all(b < a for a, b in zip(gaps, gaps[1:])))
    ck("bardeen has two horizons below the extremal charge and none above it",
       all(n == (2 if 0 < v["value"] < extremal else 1 if v["value"] == 0 else 0)
           for n, v in zip(horizons, frames)), f"{horizons}")
    past = np.linspace(1e-9, 2, SAMPLES)
    ck("bardeen g^rr is positive at every r drawn past the extremal charge",
       all(np.all(grr(past, 1, v["value"]) > 0) for v in frames if v["value"] > extremal))
    at_zero = np.linspace(0.05, 2, 41)
    ck("bardeen at g = 0 is Schwarzschild's", np.allclose(K(at_zero, 1, 0.0), K_s(at_zero, 1), rtol=1e-13)
       and np.allclose(R(at_zero, 1, 0.0), 0, atol=0))
    last = float(BARDEEN_LAST)
    centre_last = float(K.exact([0, 1, last]))
    ck("bardeen K at the centre at the last frame is 24 r_s^2/g^6", abs(centre_last * last ** 6 / 24 - 1) < 1e-12)

    def first(predicate):
        """The least g at which `predicate` holds, by bisection between frames, where it holds at
        the last frame and fails at the first that is not Schwarzschild's."""
        lo, hi = float(BARDEEN_STEP) / 4, last
        if predicate(lo) or not predicate(hi):
            raise SystemExit("bardeen: the clip of the movie does not change once as g grows")
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            lo, hi = (lo, mid) if predicate(mid) else (mid, hi)
        return hi

    def inside(fn, scale, top):
        def held(value):
            up, down = beyond(fn, x_axis, scale, [1, value])
            return not (up if top else down)
        return held
    K_top = first(inside(K, K_scale, True))
    R_top = first(inside(R, R_scale, True))
    R_bottom = first(inside(R, R_scale, False))
    ck("bardeen K at the centre leaves the top where 24 r_s^2/g^6 = 10^5",
       abs(K_top - (24e-5) ** (1 / 6)) < 1e-3, f"{K_top:.4f}")
    movie_clip = ("The scales are those of $g = r_s/3$, and we cut each curve where it leaves its panel: "
                  f"Bardeen's $K$ passes the top for $g < {K_top:.3f}\\,r_s$ and its $R$ for "
                  f"$0 < g < {R_top:.3f}\\,r_s$, and $R$ passes the bottom for $0 < g < {R_bottom:.3f}\\,r_s$.")
    movie_view = {
        "id": "growing", "label": f"$0 \\le g \\le r_s/{int(1 / BARDEEN_LAST)}$", "kind": "line", "x": x_axis,
        "panels": panels, "curves": references(), "marks": [], "legend": legend,
        "caption": [
            "The Kretschmann scalar $K$ and the Ricci scalar $R$ of Bardeen's black hole against $r$ as its charge "
            f"grows at fixed mass ($0 \\le g \\le r_s/{int(1 / BARDEEN_LAST)}$), beside Schwarzschild's. At $g = 0$ "
            "the metric is Schwarzschild's, and as $g$ grows the curvature at the centre, $K = 24r_s^2/g^6$, falls, "
            f"to ${centre_last:.0f}/r_s^4$ at $g = r_s/{int(1 / BARDEEN_LAST)}$.",
            "The two horizons approach each other as $g$ grows and meet at the extremal charge "
            f"$g = 2r_s/3\\sqrt{{3}} = {extremal:.3f}\\,r_s$. Past it $g^{{rr}}$ is positive at every $r$, and no "
            "horizon is left.",
        ],
        "settings": f"$r_s = 1$, the unit of every length, and $g$ in steps of $r_s/{int(1 / BARDEEN_STEP)}$.",
        "clip": movie_clip,
        "movie": {"variable": "$g$", "seconds": 5, "frames": frames},
    }
    return [still, movie_view]


def nr_published(src, metric_id, system_id, field, coords, params, indices=None, functions=None):
    """A published value read as Published reads it, its parameters named as the chart spells them."""
    return Published(src, metric_id, system_id, field, coords, [p.replace("\\", "") for p in params], indices,
                     functions)


# ---------------------------------------------------------------- Bach and Weyl's ring

BW_MASS = Fraction(1, 2)          # m in units of the ring's radius a, as every drawing of the ring takes it
BW_CELLS = (150, 200)             # the cells of the map across and up, each a fiftieth of a by a fiftieth
BW_RANGE = (1e-3, 1e3)            # the colours of K a^4 run from the first to the second


def bach_weyl_ring(ck, src):
    """Bach and Weyl's ring at m = a/2, its Kretschmann scalar over the half plane of Weyl's
    rho and z.

    Weyl's chart publishes K in rho, z, Weyl's two functions psi and gamma and the derivatives of
    psi, and the reader writes psi and gamma out in the complete elliptic integrals of their
    definitions, so K is a function of rho and z alone, the same at every t and phi. The ring is
    a vacuum and publishes R = 0, so the figure is K alone, a map over the box 0 <= rho <= 3a,
    |z| <= 2a, which holds the ring, the disc it spans and the axis, at the centres of cells a
    fiftieth of a across, coloured on a logarithmic scale from 10^-3 to 10^3 in units of 1/a^4.
    The ring is marked at rho = a, z = 0, and the two points of the axis where K vanishes, found
    by the golden section in sixty digits on the axis.

    What the captions state is checked here, in sixty digits where floating point would overflow:
    K is finite at every cell, on the axis and on the disc; approaching the ring, at distances
    l_1 of 10^-2, 10^-3 and 10^-4 in Weyl's coordinates, it grows faster than any power from
    outside, along the plane, as a power from directly above, and falls to zero from inside,
    along the disc; on the axis it is Semerak's closed form of 2016 and vanishes where his
    bracket does, at z = +-z_0, and every published component of the Riemann tensor vanishes
    there with it; and the clip stated is where K leaves the colours, on a grid four times as
    fine as the map's."""
    m, a = float(BW_MASS), 1.0
    params = {"m": m, "a": a}
    K = nr_published(src, "bach_weyl_ring", "weyl", "kretschmann", ["rho", "z"], ["m", "a"])
    x_axis, y_axis = axis("$\\rho/a$", 0, 3), axis("$z/a$", -2, 2)
    aspect = (x_axis["to"] - x_axis["from"]) / (y_axis["to"] - y_axis["from"])
    shape_checks(ck, "bach_weyl_ring", [aspect])
    scale = {"kind": "log", "from": BW_RANGE[0], "to": BW_RANGE[1]}
    nx, ny = BW_CELLS
    xs, ys, values = cells(K, x_axis, y_axis, nx, ny, [m, a])
    ck("bach_weyl_ring K is finite and positive at every cell", np.all(np.isfinite(values)) and np.all(values > 0),
       f"from {np.min(values):.3e} to {np.max(values):.3e}")
    # Sixty digits at forty cells spread over the map and at the four beside the ring.
    rng = np.random.default_rng(1922)
    picks = {(int(i), int(j)) for i, j in zip(rng.integers(0, nx, 40), rng.integers(0, ny, 40))}
    picks |= {(i, j) for i in (nx // 3 - 1, nx // 3) for j in (ny // 2 - 1, ny // 2)}
    picks = sorted(picks)
    agrees(ck, "bach_weyl_ring map", K, [[xs[i], ys[j], m, a] for i, j in picks], [values[j, i] for i, j in picks])

    def exact(rho, z):
        return K.exact([rho, z, m, a])

    # The ring approached along three directions at distances l_1 of 10^-2, 10^-3 and 10^-4.
    # The ring approached along three directions, at distances l_1 from 10^-6 to 10^-8 in Weyl's
    # coordinates, in eighty digits, since K reaches e^(10^7) there.
    lengths = (1e-6, 1e-7, 1e-8)

    def logs(angle):
        return [float(sp.log(K.exact([a + l * math.cos(angle), l * math.sin(angle), m, a], 80))) for l in lengths]
    rate = 4 * m ** 2 / (math.pi ** 2 * a)
    outside, inside, above = logs(0.0), logs(math.pi), logs(math.pi / 2)
    ck("bach_weyl_ring K grows as e^(4m^2/pi^2 a l_1) from outside the ring, along the plane",
       all(abs(l * v / rate - 1) < 1e-3 for l, v in zip(lengths, outside)), f"l_1 ln K {[l * v for l, v in zip(lengths, outside)]}")
    ck("bach_weyl_ring K falls as e^(-4m^2/pi^2 a l_1) from inside the ring, along the disc",
       all(abs(-l * v / rate - 1) < 1e-3 for l, v in zip(lengths, inside)), f"l_1 ln K {[l * v for l, v in zip(lengths, inside)]}")
    # Directly above, ln K is a quadratic in ln l_1 with the leading coefficient 2m^2/pi^2 a^2, so its
    # second difference over equal steps of ln l_1 is twice that times the step squared.
    step = math.log(10)
    curving = (above[2] - 2 * above[1] + above[0]) / step ** 2
    ck("bach_weyl_ring K grows as e^((2m^2/pi^2 a^2) ln^2 l_1) from directly above the ring",
       abs(curving / (4 * m ** 2 / (math.pi ** 2 * a ** 2)) - 1) < 1e-3, f"{curving:.6f}")
    ck("bach_weyl_ring the drawn mass gives the rates the caption writes", rate == a / math.pi ** 2
       and 2 * m ** 2 / (math.pi ** 2 * a ** 2) == 1 / (2 * math.pi ** 2))

    # The axis, where K is finite and vanishes at z = +-z_0.
    axis_rho = 1e-12
    on_axis = [float(exact(axis_rho, z)) for z in np.linspace(0, 2, 21)]
    ck("bach_weyl_ring K is finite on the axis", all(math.isfinite(v) for v in on_axis), f"{on_axis[0]:.4f} at the centre")
    lo, hi = 0.6, 1.2
    golden = (math.sqrt(5) - 1) / 2
    for _ in range(70):
        c, d = hi - golden * (hi - lo), lo + golden * (hi - lo)
        if exact(axis_rho, c) < exact(axis_rho, d):
            hi = d
        else:
            lo = c
    z0 = 0.5 * (lo + hi)
    least = float(exact(axis_rho, z0))
    ck("bach_weyl_ring K vanishes on the axis at z_0", least < 1e-20 * on_axis[0], f"z_0 = {z0:.6f}, K = {least:.1e}")
    ck("bach_weyl_ring K is the same at -z_0", float(exact(axis_rho, -z0)) < 1e-20 * on_axis[0])
    # Semerák's closed form of K on the axis, his equation of section V D (Physical Review D 94,
    # 104021, 2016), 12m^2((2z^2 - a^2)sqrt(z^2 + a^2) - 2mz^2)^2 / ((z^2 + a^2)^6 e^(4m/sqrt(z^2 + a^2))),
    # which has one zero on each half of the axis, the root of the bracket.
    def semerak(z):
        w = math.sqrt(z * z + a * a)
        return 12 * m ** 2 * ((2 * z * z - a * a) * w - 2 * m * z * z) ** 2 / (w ** 12 * math.exp(4 * m / w))
    zs = np.linspace(0, 2, 21)
    ck("bach_weyl_ring K on the axis is Semerak's closed form",
       all(abs(v - semerak(z)) <= 1e-12 * max(1.0, semerak(z)) for v, z in zip(on_axis, zs)), f"centre {semerak(0):.6f}")
    lo, hi = 0.5, 1.5
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        w = math.sqrt(mid * mid + a * a)
        lo, hi = (mid, hi) if (2 * mid * mid - a * a) * w - 2 * m * mid * mid < 0 else (lo, mid)
    ck("bach_weyl_ring K vanishes on the axis where Semerak's bracket does", abs(z0 - lo) < 1e-9, f"{lo:.9f}")
    entry = json.loads((build.METRICS_DIR / "bach_weyl_ring.json").read_text(encoding="utf-8"))
    weyl = next(s for s in entry["coordinates"] if s["id"] == "weyl")
    worst = 0.0
    for comp in weyl["riemann"]["variants"]["ulll"]["nonzero"]:
        riemann = nr_published(src, "bach_weyl_ring", "weyl", "riemann", ["rho", "z"], ["m", "a"], comp["indices"])
        worst = max(worst, abs(float(riemann.exact([axis_rho, z0, m, a]))))
    ck("bach_weyl_ring every published Riemann component vanishes at z_0 on the axis", worst < 1e-9, f"{worst:.1e}")

    # The circumference of the circle about the axis through a point, 2 pi sqrt(g_phiphi), grows
    # without bound toward the ring along every direction drawn.
    gpp = nr_published(src, "bach_weyl_ring", "weyl", "metric_components", ["rho", "z"], ["m", "a"], ("\\phi", "\\phi"))
    growth = [[math.sqrt(float(gpp.exact([a + l * math.cos(angle), l * math.sin(angle), m, a]))) for l in (1e-2, 1e-4, 1e-6)]
              for angle in (0.0, math.pi / 2, math.pi)]
    ck("bach_weyl_ring the circles about the axis grow without bound toward the ring",
       all(r[0] < r[1] < r[2] and r[2] > 10 * a for r in growth), f"{growth}")

    # Where K leaves the colours, on a grid four times as fine as the map's.
    fx, fy, fine = cells(K, x_axis, y_axis, 4 * nx, 4 * ny, [m, a])
    FX, FY = np.meshgrid(fx, fy)
    high, low = ~(fine <= BW_RANGE[1]), fine < BW_RANGE[0]
    from_ring = np.hypot(FX - a, FY)
    from_zero = np.minimum(np.hypot(FX, FY - z0), np.hypot(FX, FY + z0))
    d_high, d_low = float(from_ring[high].max()), float(from_zero[low].max())
    ck("bach_weyl_ring K passes the top of its colours only beside the ring", d_high < 0.2, f"within {d_high:.4f} a")
    ck("bach_weyl_ring K passes the bottom of its colours only beside the zeros on the axis", d_low < 0.2,
       f"within {d_low:.4f} a")
    rows = [[rounded(v) for v in row] for row in values]
    view = {
        "id": "meridian", "label": "the half plane of $\\rho$ and $z$", "kind": "map",
        "x": x_axis, "y": y_axis,
        "colour": {"name": "$K\\,a^4$", "scale": scale, "ticks": value_ticks(scale)},
        "map": {"metric": "bach_weyl_ring", "system": "weyl", "scalar": "kretschmann", "x": "\\rho", "y": "z",
                "params": params, "nx": nx, "ny": ny, "cells": rows},
        "marks": [{"class": "singular", "at": [a, 0.0]}, {"class": "zero", "at": [0.0, z0]},
                  {"class": "zero", "at": [0.0, -z0]}],
        "legend": [["point", "singular", "the ring, $\\rho = a$, $z = 0$"],
                   ["point", "zero", "where $K = 0$"]],
        "caption": [
            "The Kretschmann scalar $K$ of the Bach-Weyl ring over the half plane of Weyl's $\\rho$ and $z$ "
            "($m = a/2$), the same at every $t$ and $\\phi$. The ring is a vacuum, so its Ricci scalar vanishes, "
            "and $K$ is finite everywhere off the ring, on the axis and on the disc the ring spans.",
            "At the ring, $\\rho = a$, $z = 0$, $K$ diverges in a way that depends on the direction it is approached "
            "from, $l_1$ being the distance to the ring in Weyl's coordinates. From outside, along the plane, it "
            "grows as $e^{a/\\pi^2l_1}$, and from directly above and below as $e^{\\ln^2(l_1/a)/2\\pi^2}$, each "
            "faster than any power of $l_1$. From inside, along the disc, it falls to zero as $e^{-a/\\pi^2l_1}$.",
            f"On the axis $K$ vanishes at $z = \\pm {z0:.3f}\\,a$, and there every component of the Riemann tensor "
            "vanishes with it. Weyl's $\\rho$ and $z$ are the cylindrical coordinates of the flat space in which "
            "$c^2\\psi$ is Newton's potential of the ring, and the metric stretches them: the ring, the circle "
            "$\\rho = a$, has an infinite circumference.",
        ],
        "settings": "$m = GM/c^2 = a/2$, with the radius $a$ of the ring in Weyl's coordinates the unit of every "
                    "length.",
        "clip": ("We colour $K\\,a^4$ on a logarithmic scale from $10^{-3}$ to $10^{3}$ and give a value beyond it "
                 f"the colour of the end it passes: $K$ passes $10^{{3}}/a^4$ within ${d_high:.2f}\\,a$ of the ring "
                 f"and falls below $10^{{-3}}/a^4$ within ${d_low:.2f}\\,a$ of the two points where it vanishes."),
    }
    return [view]


# ---------------------------------------------------------------- Bonnor and Vaidya's radiating star

BV_BOX = ((0, 3), (-1, 3))        # r/M and u/M: the star's trapped spheres, and the stretch over which it radiates
BV_CELLS = (150, 200)             # the cells of the map across and up, each a fiftieth of M on a side
BV_RANGE = (1e-3, 1e5)            # the colours of K M^4 run from the first to the second
# The star as it radiates, at M = 1: its mass falls from M to M/2 and its charge from 24M/25 to
# 27M/50 by one tanh of width l = M/2 about u = M, where both are 3M/4.
BV_FUNCTIONS = {"m": "(3 - tanh(2*(u - 1)))/4",
                "q": "Rational(3, 4) - Rational(21, 100)*tanh(2*(u - 1))"}


def bonnor_vaidya(ck, src):
    """Bonnor and Vaidya's radiating star, its Kretschmann scalar over the plane of the outgoing
    chart's u and r, as the star sends out charged null dust and its mass and charge fall.

    The outgoing chart publishes K in r and the free functions m(u) and q(u) of the retarded time,
    the chart of a star that radiates, which the figure declares as one smooth stretch of radiation:
    m = M(3 - tanh x)/4 and q = (3/4 - (21/100) tanh x)M with x = (u - M)/l and l = M/2,
    so that the mass falls from M to M/2 and the charge from 0.96 M to 0.54 M, the mass faster, and
    the two are equal, 3M/4, at u = M. The published Ricci scalar vanishes, so the figure is K
    alone: a map over 0 <= r <= 3M and -M <= u <= 3M, at M = 1, in 150 cells across and 200 up,
    each a fiftieth of M on a side, coloured on a logarithmic scale from 10^-3 to 10^5 in units of
    1/M^4, with the marginally anti-trapped spheres, the zeros of the published g^rr, marked as the one
    curve they make, the inner and the outer meeting at u = M.

    What the caption states is checked here: the published Ricci scalar vanishes; K at every cell
    is the closed form 8(6m^2r^2 - 12mq^2r + 7q^4)/r^8 at that u's m and q, and grows toward r = 0
    as 56q^4/r^8; K falls as u grows at every r drawn; the declared m and q start at M and 0.96 M,
    end at M/2 and 0.54 M and are equal at u = M alone, with q < m before it and q > m after; the
    anti-trapped spheres are r = m +- sqrt(m^2 - q^2), start at 1.28 M and 0.72 M, meet at 3M/4 at
    u = M, and every sphere between them has g^rr < 0, so that every future light cone there
    points to larger r; after u = M g^rr is positive at every r drawn; and the clip stated is where
    K leaves the colours, on a grid four times as fine."""
    M = 1.0
    OUT = "eddington_finkelstein_outgoing"
    (r0, r1), (u0, u1) = BV_BOX
    x_axis, y_axis = axis("$r/M$", r0, r1), axis("$u/M$", u0, u1)
    aspect = (x_axis["to"] - x_axis["from"]) / (y_axis["to"] - y_axis["from"])
    shape_checks(ck, "bonnor_vaidya", [aspect])
    scale = {"kind": "log", "from": BV_RANGE[0], "to": BV_RANGE[1]}

    def read(field, indices=None):
        return nr_published(src, "bonnor_vaidya", OUT, field, ["r", "u"], [], indices, BV_FUNCTIONS)

    ricci = read("ricci_scalar")
    ck("bonnor_vaidya the published Ricci scalar vanishes", ricci.expr == 0, f"{ricci.expr}")
    K, grr = read("kretschmann"), read("inverse_metric_components", ("r", "r"))

    def mass(u):
        return (3 - np.tanh(2 * (np.asarray(u, dtype=float) - 1))) / 4

    def charge(u):
        return 0.75 - 0.21 * np.tanh(2 * (np.asarray(u, dtype=float) - 1))

    def closed(r, u):
        m, q = mass(u), charge(u)
        return 8 * (6 * m ** 2 * r ** 2 - 12 * m * q ** 2 * r + 7 * q ** 4) / r ** 8

    # The declared profile.
    us = np.linspace(u0, u1, SAMPLES)
    ck("bonnor_vaidya the mass falls from M to M/2 and the charge from 0.96 M to 0.54 M",
       np.all(np.diff(mass(us)) < 0) and np.all(np.diff(charge(us)) < 0)
       and abs(mass(-1e3) - 1) < 1e-15 and abs(mass(1e3) - 0.5) < 1e-15
       and abs(charge(-1e3) - 0.96) < 1e-15 and abs(charge(1e3) - 0.54) < 1e-15,
       f"m from {mass(u0):.4f} to {mass(u1):.4f}, q from {charge(u0):.4f} to {charge(u1):.4f} across the box")
    ck("bonnor_vaidya the mass and the charge are equal, 3M/4, at u = M alone, q < m before and q > m after",
       mass(1.0) == 0.75 and charge(1.0) == 0.75 and np.all((charge(us) < mass(us)) == (us < 1)))

    xs, ys, values = cells(K, x_axis, y_axis, *BV_CELLS, [])
    X, Y = np.meshgrid(xs, ys)
    ck("bonnor_vaidya K is finite and positive at every cell", np.all(np.isfinite(values)) and np.all(values > 0),
       f"from {np.min(values):.3e} to {np.max(values):.3e}")
    ck("bonnor_vaidya K at every cell is 8(6m^2r^2 - 12mq^2r + 7q^4)/r^8 at that u's m and q",
       np.allclose(values, closed(X, Y), rtol=1e-12, atol=0))
    rng = np.random.default_rng(1970)
    nx, ny = BV_CELLS
    picks = sorted({(int(i), int(j)) for i, j in zip(rng.integers(0, nx, 40), rng.integers(0, ny, 40))}
                   | {(i, j) for i in (0, 1) for j in (ny // 2 - 1, ny // 2)})
    agrees(ck, "bonnor_vaidya map", K, [[xs[i], ys[j]] for i, j in picks], [values[j, i] for i, j in picks])
    centre = [float(K.exact([R, U]) * R ** 8) / (56 * float(charge(U)) ** 4) for R in (1e-6, 1e-8) for U in (-1.0, 1.0, 3.0)]
    ck("bonnor_vaidya K grows toward r = 0 as 56q^4/r^8", all(abs(c - 1) < 1e-5 for c in centre), f"{centre[:3]}")
    fine_r = np.linspace(1e-2, r1, 600)
    FU, FR = np.meshgrid(us[::4], fine_r)
    ck("bonnor_vaidya K falls as u grows at every r drawn", np.all(np.diff(K(FR, FU), axis=1) < 0))

    # The marginally anti-trapped spheres: the zeros of the published g^rr at each u, from the inner
    # to the outer, closing at u = M, where the two meet at r = 3M/4.
    meet = 1.0
    ts = np.linspace(0, 1, 400)[:-1]
    track = u0 + (meet - u0) * (1 - (1 - ts) ** 2)
    inner, outer = [], []
    for u in track:
        z = sorted(_g_zeros(grr, u, r1))
        if len(z) != 2:
            ck(f"bonnor_vaidya two anti-trapped spheres at u = {u:.4f}", False, f"{z}")
            break
        inner.append([z[0], float(u)])
        outer.append([z[1], float(u)])
    m_t, q_t = mass(track), charge(track)
    ck("bonnor_vaidya the anti-trapped spheres are r = m +- sqrt(m^2 - q^2)",
       np.allclose([p[0] for p in inner], m_t - np.sqrt(m_t ** 2 - q_t ** 2), rtol=0, atol=1e-10)
       and np.allclose([p[0] for p in outer], m_t + np.sqrt(m_t ** 2 - q_t ** 2), rtol=0, atol=1e-10))
    start = (outer[0][0], inner[0][0])
    ck("bonnor_vaidya the anti-trapped spheres start at 1.28 M and 0.72 M", abs(start[0] - 1.28) < 5e-4 and abs(start[1] - 0.72) < 5e-4,
       f"{start[0]:.6f}, {start[1]:.6f} at u = -M")
    at_meet = float(grr(0.75, meet))
    slope = float((grr(0.75 + 1e-6, meet) - grr(0.75 - 1e-6, meet)) / 2e-6)
    ck("bonnor_vaidya the two meet at r = 3M/4 at u = M, a double zero of g^rr", abs(at_meet) < 1e-14 and abs(slope) < 1e-8,
       f"g^rr = {at_meet:.1e}, its slope {slope:.1e}")
    between = [float(grr(0.5 * (a[0] + b[0]), a[1])) for a, b in zip(inner, outer)]
    ck("bonnor_vaidya every sphere between the two has g^rr < 0, its future cone pointing to larger r",
       all(g < 0 for g in between), f"{max(between):.3e}")
    later = grr(FR[:, FU[0] > meet + 1e-6], FU[:, FU[0] > meet + 1e-6])
    ck("bonnor_vaidya after u = M g^rr is positive at every r drawn", np.all(later > 0), f"{np.min(later):.3e}")
    curve = inner + [[0.75, meet]] + outer[::-1]
    ck("bonnor_vaidya the anti-trapped spheres lie inside the box", all(r0 <= x <= r1 and u0 <= y <= u1 for x, y in curve))

    # The clip, on a grid four times as fine as the map's.
    fx, fy, fine = cells(K, x_axis, y_axis, 4 * nx, 4 * ny, [])
    FX = np.meshgrid(fx, fy)[0]
    d_high = float(FX[~(fine <= BV_RANGE[1])].max())
    lo, hi = 1e-3, 2.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if closed(mid, u0) > BV_RANGE[1] else (lo, mid)
    cut = lo
    ck("bonnor_vaidya K passes the top of its colours only inside the radius stated, reached at u = -M",
       abs(d_high - cut) < 2 * (r1 - r0) / (4 * nx), f"{cut:.6f} against {d_high:.6f}")
    least = float(closed(r1, u1))
    ck("bonnor_vaidya K is least at r = 3M, u = 3M, inside the colours", least <= float(np.min(fine)) and least > BV_RANGE[0],
       f"{least:.4e}")
    view = {
        "id": "radiating", "label": "the plane of $u$ and $r$", "system": OUT, "kind": "map", "x": x_axis, "y": y_axis,
        "colour": {"name": "$K\\,M^4$", "scale": scale, "ticks": value_ticks(scale)},
        "map": {"metric": "bonnor_vaidya", "system": OUT, "scalar": "kretschmann", "x": "r", "y": "u", "params": {},
                "functions": dict(BV_FUNCTIONS), "nx": nx, "ny": ny,
                "cells": [[rounded(v) for v in row] for row in values]},
        "marks": [{"class": "horizon", "points": curve}],
        "legend": [["line", "horizon", "a marginally anti-trapped sphere, $r^2 - 2mr + q^2 = 0$"]],
        "caption": [
            "The Kretschmann scalar $K$ of Bonnor and Vaidya's radiating star over the plane of $u$ and $r$ ($m$ from $M$ "
            "to $M/2$, $q$ from $0.96\\,M$ to $0.54\\,M$), the same at every $\\theta$ and $\\phi$. The star sends out "
            "charged null dust over a stretch of retarded time about $u = M$, and its mass and charge fall together, the "
            "mass faster, so that at each $u$ the curvature is $K = 8(6m^2r^2 - 12mq^2r + 7q^4)/r^8$ with that moment's "
            "$m$ and $q$. At every $r$ it falls as the star radiates, and toward $r = 0$ it grows without bound, as "
            "$56q^4/r^8$.",
            f"The marginally anti-trapped spheres, $r = m \\pm \\sqrt{{m^2 - q^2}}$, start at $r = {start[0]:.2f}\\,M$ and "
            f"$r = {start[1]:.2f}\\,M$, close in on each other as the star radiates, and meet at $r = 0.75\\,M$ at "
            "$u = M$, where $m = q = 0.75\\,M$. Between them every future light cone points to larger $r$, as in a white "
            "hole. After $u = M$ the charge exceeds the mass, and no sphere is anti-trapped. The Ricci scalar vanishes "
            "throughout, since the stress tensors of the null dust and of the electric field are both traceless.",
        ],
        "settings": "$M = 1$, the unit of every length, and $\\ell = M/2$.",
        "input": ("$m = M(3 - \\tanh x)/4$ and $q = (0.75 - 0.21\\tanh x)\\,M$, with $x = (u - M)/\\ell$: a star whose "
                  "mass falls from $M$ to $M/2$ and whose charge falls from $0.96\\,M$ to $0.54\\,M$ as it radiates, the "
                  "two equal at $u = M$."),
        "clip": ("We colour $K\\,M^4$ on a logarithmic scale from $10^{-3}$ to $10^{5}$ and give a value beyond it the "
                 f"colour of the end it passes: $K$ passes $10^{{5}}/M^4$ inside $r = {cut:.3f}\\,M$ and is least, "
                 f"${tex_number(least)}/M^4$, at $r = 3M$, $u = 3M$."),
    }
    return [view]


def _g_zeros(grr, u, r1, n=SAMPLES):
    """The zeros of the published g^rr along r at one u, each a change of sign between two samples
    found by bisection to the last digit."""
    rs = np.linspace(1e-9, r1, n)
    g = grr(rs, u)
    out = []
    for i in range(n - 1):
        a, b, fa, fb = rs[i], rs[i + 1], g[i], g[i + 1]
        if not (np.isfinite(fa) and np.isfinite(fb)) or fa == 0 or np.sign(fa) == np.sign(fb):
            continue
        for _ in range(200):
            m = 0.5 * (a + b)
            fm = float(grr(m, u))
            if m in (a, b) or fm == 0:
                break
            a, b, fa = (m, b, fm) if np.sign(fm) == np.sign(fa) else (a, m, fa)
        out.append(float(0.5 * (a + b)))
    return out


# ---------------------------------------------------------------- Boulware and Deser's black hole

BD_R0 = Fraction(13, 12)          # r_0 in units of the horizon radius r_h, as every drawing of the black hole takes it
BD_ELL = Fraction(5, 12)          # the Gauss-Bonnet length in units of r_h, as every drawing of the black hole takes it
BD_STEP = Fraction(1, 20)         # the step of l/r_0 from one frame of the movie to the next, and its first frame
BD_LAST = Fraction(3, 2)          # the last frame's l/r_0


def boulware_deser(ck, src):
    """Boulware and Deser's black hole and its anti-de Sitter branch, their Kretschmann and Ricci
    scalars against r.

    Every chart of the black hole publishes the same K and R, functions of r alone through
    W = sqrt(r^4 + 4 l^2 r_0^2), and the ingoing Eddington-Finkelstein chart covers every r > 0,
    so it is the one read; the horizon is the zero of its published g^rr. The first view is the
    black hole at r_0 = 13 r_h/12 and l = 5 r_h/12, the values every other drawing of it takes,
    with r from 0 to 6 r_h, the box of its spacetime diagrams: K on a logarithmic axis from 10^-5
    to 10^5 and R on a signed one, linear within 10^-3/r_h^2 of zero, from -1 to 10^4, beside
    Tangherlini's two scalars of the same mass, read from tangherlini.json at its r_h = r_0. The
    second is the anti-de Sitter branch, its own chart, at r_0 = l, the unit of length, with r
    from 0 to 3 l: K on a logarithmic axis from 10 to 10^5 and R, negative at every r, on a
    signed one from -10^5 to -10. The third steps l from r_0/20 to 3 r_0/2 at fixed mass, with
    r from 0 to 3 r_0, through l = r_0, where the horizon r_h = sqrt(r_0^2 - l^2) closes. The
    movie starts at r_0/20 and not at l = 0, Tangherlini's black hole, whose published R is
    0 times a sum that rounding in W - r^2 leaves at 10^-9 rather than at zero; Tangherlini's
    scalars are dashed in every frame instead.

    What the captions state is checked here, in sixty digits: near the centre K is
    12 r_0^2/l^2 r^4 on both branches and R is 6 r_0/l r^2 on the black hole and its negative on
    the branch; the black hole's metric function is 1 - r_0/l at the centre; its K falls from the
    centre outward, its R changes sign once, inside the horizon, at the radius stated, and is
    least, -(10 - 7 sqrt 2)/l^2, at r = sqrt(2 l r_0); far out its K is Tangherlini's and its R
    is -12 l^2 r_0^4/r^8, and the branch's K and R are anti-de Sitter space's 40/l^4 and -20/l^2
    in five dimensions; the branch's g^rr is positive at every r; in the movie the horizon moves
    in as l grows and is gone from l = r_0, and K near the centre and the least R fall and rise;
    and every cut the clips state is where the curve leaves its panel."""
    r0, ell = float(BD_R0), float(BD_ELL)
    hole_params = {"r_0": r0, "\\ell": ell}
    K = nr_published(src, "boulware_deser", "eddington_finkelstein_ingoing", "kretschmann", ["r"], ["r_0", "\\ell"])
    R = nr_published(src, "boulware_deser", "eddington_finkelstein_ingoing", "ricci_scalar", ["r"], ["r_0", "\\ell"])
    grr = nr_published(src, "boulware_deser", "eddington_finkelstein_ingoing", "inverse_metric_components", ["r"],
                       ["r_0", "\\ell"], ("r", "r"))
    Kb = nr_published(src, "boulware_deser", "spherical_plus", "kretschmann", ["r"], ["r_0", "\\ell"])
    Rb = nr_published(src, "boulware_deser", "spherical_plus", "ricci_scalar", ["r"], ["r_0", "\\ell"])
    grr_b = nr_published(src, "boulware_deser", "spherical_plus", "inverse_metric_components", ["r"], ["r_0", "\\ell"],
                         ("r", "r"))
    K_t = nr_published(src, "tangherlini", "spherical", "kretschmann", ["r"], ["r_h"])
    R_t = nr_published(src, "tangherlini", "spherical", "ricci_scalar", ["r"], ["r_h"])
    legend_K_R = [["line", "K", "$K$, the Kretschmann scalar"], ["line", "R", "$R$, the Ricci scalar"]]
    legend = legend_K_R + [["line", "reference", "Tangherlini's black hole of the same mass"],
                           ["line", "horizon", "a horizon, where $g^{rr} = 0$"]]

    def stack(K_name, R_name, K_scale, R_scale):
        panels = [{"id": "K", "name": K_name, "aspect": 2, "scale": K_scale, "ticks": value_ticks(K_scale)},
                  {"id": "R", "name": R_name, "aspect": 2, "scale": R_scale, "ticks": value_ticks(R_scale)}]
        return panels, {p["id"]: PLOT_W / p["aspect"] for p in panels}

    def cut(fn, x_axis, scale, values):
        """The radius inside which a curve leaves its panel, by its top or its bottom, where it
        leaves it about the centre alone, and None where it never leaves it."""
        up, down = beyond(fn, x_axis, scale, values)
        outside = sorted(up + down)
        step = (x_axis["to"] - x_axis["from"]) / (SAMPLES - 1)
        joined = all(abs(b[0] - a[1] - step) < 1e-9 for a, b in zip(outside, outside[1:]))
        if outside and (outside[0][0] != x_axis["from"] or not joined):
            raise SystemExit(f"boulware_deser: a curve leaves its panel other than about the centre: {up} {down}")
        return outside[-1][1] if outside else None

    def near(fn, values, law, where):
        """A published scalar against the law it follows near the centre, in sixty digits."""
        ratios = [float(fn.exact([r] + values)) / law(r) for r in (1e-6, 1e-7, 1e-8)]
        ck(where, all(abs(q - 1) < 1e-5 for q in ratios), f"{ratios}")

    def far(fn, values, law, where):
        ratios = [float(fn.exact([r] + values)) / law(r) for r in (1e4, 1e5)]
        ck(where, all(abs(q - 1) < 1e-6 for q in ratios), f"{ratios}")

    # The black hole at r_0 = 13 r_h/12 and l = 5 r_h/12.
    x_axis = axis("$r/r_h$", 0, 6)
    K_scale = {"kind": "log", "from": 1e-5, "to": 1e5}
    R_scale = {"kind": "symlog", "from": -1, "to": 1e4, "linear": 1e-3}
    panels, height = stack("$K\\,r_h^4$", "$R\\,r_h^2$", K_scale, R_scale)
    shape_checks(ck, "boulware_deser hole", [p["aspect"] for p in panels])
    reference = {"r_h": r0}
    curves = (line_curves(ck, "boulware_deser Tangherlini", "K", "reference", K_t, "r", x_axis, K_scale, height["K"],
                          reference)
              + line_curves(ck, "boulware_deser Tangherlini", "R", "reference", R_t, "r", x_axis, R_scale, height["R"],
                            reference)
              + line_curves(ck, "boulware_deser hole", "K", "K", K, "r", x_axis, K_scale, height["K"], hole_params)
              + line_curves(ck, "boulware_deser hole", "R", "R", R, "r", x_axis, R_scale, height["R"], hole_params))
    values = [r0, ell]
    horizons = zeros(ck, "boulware_deser hole horizons", grr, 0, 6, values)
    ck("boulware_deser has one horizon, at r_h = sqrt(r_0^2 - l^2) = 1",
       len(horizons) == 1 and abs(horizons[0] - 1) < 1e-12, f"{horizons}")
    marks = [{"class": "horizon", "x": x} for x in horizons]
    ck("boulware_deser Tangherlini's published K is 72 r_h^4/r^8",
       all(abs(float(K_t.exact([r, r0])) * r ** 8 / (72 * r0 ** 4) - 1) < 1e-12 for r in (0.1, 1.0, 5.0)))
    centre_f = float(grr.exact([1e-30] + values))
    ck("boulware_deser the metric function at the centre is 1 - r_0/l = -8/5",
       abs(centre_f - (1 - r0 / ell)) < 1e-12 and abs(centre_f + 1.6) < 1e-12, f"{centre_f}")
    near(K, values, lambda r: 12 * r0 ** 2 / (ell ** 2 * r ** 4), "boulware_deser K grows as 12 r_0^2/l^2 r^4 at the centre")
    near(R, values, lambda r: 6 * r0 / (ell * r ** 2), "boulware_deser R grows as 6 r_0/l r^2 at the centre")
    rs = np.linspace(0, 6, SAMPLES)[1:]
    ck("boulware_deser K falls from the centre outward", np.all(np.diff(K(rs, *values)) < 0))
    sign = zeros(ck, "boulware_deser R", R, 0, 6, values)
    # R changes sign where 2r^4 + 4r^2 W - 3W^2 vanishes, W = (2 + sqrt 10) r^2/3.
    r_sign = math.sqrt(2 * ell * r0) * (9 / (5 + 4 * math.sqrt(10))) ** 0.25
    ck("boulware_deser R changes sign once, inside the horizon, where W = (2 + sqrt 10) r^2/3",
       len(sign) == 1 and abs(sign[0] - r_sign) < 1e-12 and sign[0] < horizons[0], f"{sign} {r_sign}")
    r_least = math.sqrt(2 * ell * r0)
    R_least = float(R.exact([r_least] + values))
    fine = np.linspace(1e-3, 6, 600001)
    ck("boulware_deser R is least at r = sqrt(2 l r_0), just inside the horizon, at -(10 - 7 sqrt 2)/l^2",
       abs(R_least * ell ** 2 / -(10 - 7 * math.sqrt(2)) - 1) < 1e-12 and R_least <= float(np.min(R(fine, *values)))
       and abs(float(fine[np.argmin(R(fine, *values))]) - r_least) < 1e-4 and r_sign < r_least < horizons[0],
       f"r = {r_least:.6f}, R = {R_least:.6f}")
    far(K, values, lambda r: 72 * r0 ** 4 / r ** 8, "boulware_deser K approaches Tangherlini's far out")
    far(R, values, lambda r: -12 * ell ** 2 * r0 ** 4 / r ** 8, "boulware_deser R falls off as -12 l^2 r_0^4/r^8 far out")
    cut_K, cut_t, cut_R = (cut(K, x_axis, K_scale, values), cut(K_t, x_axis, K_scale, [r0]),
                           cut(R, x_axis, R_scale, values))
    ck("boulware_deser R of Tangherlini stays inside its panel", cut(R_t, x_axis, R_scale, [r0]) is None)
    hole = {
        "id": "hole", "label": "the black hole", "kind": "line", "x": x_axis, "panels": panels,
        "curves": curves, "marks": marks, "legend": legend,
        "caption": [
            "The Kretschmann scalar $K$ and the Ricci scalar $R$ of Boulware and Deser's black hole against $r$ "
            "($r_0 = 13r_h/12$, $\\ell = 5r_h/12$), the same at every $t$, $\\psi$, $\\theta$, and $\\phi$, beside "
            "Tangherlini's of the same mass. Toward the centre Tangherlini's $K = 72r_0^4/r^8$ grows as $r^{-8}$. "
            "Boulware and Deser's metric function stays finite there, $1 - r_0/\\ell = -8/5$ at $r = 0$, and its "
            "$K$ grows as $12r_0^2/\\ell^2r^4$.",
            "Near the centre $R$ grows as $6r_0/\\ell r^2$. It changes sign at "
            f"$r = {sign[0]:.3f}\\,r_h$, inside the horizon, and is least, "
            f"$-(10 - 7\\sqrt{{2}})/\\ell^2 = {R_least:.3f}/r_h^2$, at $r = \\sqrt{{2\\ell r_0}} = {r_least:.3f}\\,r_h$, "
            "just inside it. Far out $K$ approaches Tangherlini's and $R$ falls off as $-12\\ell^2r_0^4/r^8$.",
        ],
        "settings": "$r_0 = 13r_h/12$ and $\\ell = 5r_h/12$, with the horizon radius $r_h = \\sqrt{r_0^2 - \\ell^2}$ "
                    "the unit of every length.",
        "clip": ("We draw $K\\,r_h^4$ on a logarithmic scale from $10^{-5}$ to $10^{5}$ and $R\\,r_h^2$ on a signed "
                 "logarithmic scale, linear within $0.001$ of zero, from $-1$ to $10^{4}$, and cut each curve where it "
                 f"leaves its panel, as Boulware and Deser's $K$ does inside $r = {cut_K:.3f}\\,r_h$, Tangherlini's "
                 f"inside $r = {cut_t:.3f}\\,r_h$, and $R$ inside $r = {cut_R:.3f}\\,r_h$."),
    }

    # The anti-de Sitter branch at r_0 = l, the unit of length.
    x_b = axis("$r/\\ell$", 0, 3)
    Kb_scale = {"kind": "log", "from": 10, "to": 1e5}
    Rb_scale = {"kind": "symlog", "from": -1e5, "to": -10, "linear": 1}
    panels_b, height_b = stack("$K\\,\\ell^4$", "$R\\,\\ell^2$", Kb_scale, Rb_scale)
    shape_checks(ck, "boulware_deser branch", [p["aspect"] for p in panels_b])
    branch_params = {"r_0": 1.0, "\\ell": 1.0}
    curves_b = (line_curves(ck, "boulware_deser branch", "K", "K", Kb, "r", x_b, Kb_scale, height_b["K"], branch_params)
                + line_curves(ck, "boulware_deser branch", "R", "R", Rb, "r", x_b, Rb_scale, height_b["R"],
                              branch_params))
    unit = [1.0, 1.0]
    wide = np.concatenate([np.linspace(0, 3, SAMPLES)[1:], np.geomspace(3, 1e4, 2001)])
    ck("boulware_deser branch g^rr is positive at every r", np.all(grr_b(wide, *unit) > 0)
       and not zeros(ck, "boulware_deser branch g^rr", grr_b, 0, 3, unit))
    near(Kb, unit, lambda r: 12 / r ** 4, "boulware_deser branch K grows as 12 r_0^2/l^2 r^4 at the centre")
    near(Rb, unit, lambda r: -6 / r ** 2, "boulware_deser branch R grows as -6 r_0/l r^2 at the centre")
    far(Kb, unit, lambda r: 40.0, "boulware_deser branch K tends to 40/l^4, anti-de Sitter space's in five dimensions")
    far(Rb, unit, lambda r: -20.0, "boulware_deser branch R tends to -20/l^2, anti-de Sitter space's in five dimensions")
    ck("boulware_deser branch R is negative at every r", np.all(Rb(wide, *unit) < 0))
    cut_Kb, cut_Rb = cut(Kb, x_b, Kb_scale, unit), cut(Rb, x_b, Rb_scale, unit)
    branch = {
        "id": "branch", "label": "the anti-de Sitter branch", "system": "spherical_plus", "kind": "line", "x": x_b,
        "panels": panels_b, "curves": curves_b, "marks": [], "legend": legend_K_R,
        "caption": [
            "The Kretschmann scalar $K$ and the Ricci scalar $R$ of the anti-de Sitter branch against $r$ "
            "($r_0 = \\ell$), the same at every $t$, $\\psi$, $\\theta$, and $\\phi$. Toward the centre $K$ grows as "
            "$12r_0^2/\\ell^2r^4$, as the black hole's does, and $R$ as $-6r_0/\\ell r^2$, the black hole's with its "
            "sign changed. Here $g^{rr}$ is positive at every $r$, and the singularity at $r = 0$ has no horizon.",
            "Far out $K$ tends to $40/\\ell^4$ and $R$ to $-20/\\ell^2$, the curvature of anti-de Sitter space of "
            "radius $\\ell$ in five dimensions, set by the Gauss-Bonnet coupling.",
        ],
        "settings": "$\\ell = 1$, the unit of every length, and $r_0 = \\ell$.",
        "clip": ("We draw $K\\,\\ell^4$ on a logarithmic scale from $10$ to $10^{5}$ and $R\\,\\ell^2$ on a signed "
                 "logarithmic scale, linear within $1$ of zero, from $-10^{5}$ to $-10$, and cut each curve where it "
                 f"leaves its panel, as $K$ does inside $r = {cut_Kb:.3f}\\,\\ell$ and $R$ inside "
                 f"$r = {cut_Rb:.4f}\\,\\ell$."),
    }

    # The movie in l at fixed mass, r_0 = 1, through l = r_0, where the horizon closes.
    x_m = axis("$r/r_0$", 0, 3)
    Km_scale = {"kind": "log", "from": 1e-3, "to": 1e5}
    Rm_scale = {"kind": "symlog", "from": -100, "to": 1e4, "linear": 0.01}
    panels_m, height_m = stack("$K\\,r_0^4$", "$R\\,r_0^2$", Km_scale, Rm_scale)
    shape_checks(ck, "boulware_deser movie", [p["aspect"] for p in panels_m])
    mass = {"r_h": 1.0}
    references = (line_curves(ck, "boulware_deser movie Tangherlini", "K", "reference", K_t, "r", x_m, Km_scale,
                              height_m["K"], mass)
                  + line_curves(ck, "boulware_deser movie Tangherlini", "R", "reference", R_t, "r", x_m, Rm_scale,
                                height_m["R"], mass))
    frames, cuts_K, cuts_R, least = [], [], [], []
    for k in range(1, int(BD_LAST / BD_STEP) + 1):
        value = float(BD_STEP * k)
        name, params = f"boulware_deser at l = {value:.2f} r_0", {"r_0": 1.0, "\\ell": value}
        frame_curves = (line_curves(ck, name, "K", "K", K, "r", x_m, Km_scale, height_m["K"], params)
                        + line_curves(ck, name, "R", "R", R, "r", x_m, Rm_scale, height_m["R"], params))
        found = zeros(ck, f"{name} horizons", grr, 0, 3, [1.0, value])
        expected = [math.sqrt(1 - value ** 2)] if value < 1 else []
        ck(f"{name} has its horizon at sqrt(r_0^2 - l^2), and none from l = r_0",
           len(found) == len(expected) and all(abs(a - b) < 1e-12 for a, b in zip(found, expected)), f"{found}")
        frames.append({"label": f"$\\ell = {value:.2f}\\,r_0$", "value": value, "curves": frame_curves,
                       "marks": [{"class": "horizon", "x": x} for x in found]})
        cuts_K.append(cut(K, x_m, Km_scale, [1.0, value]))
        cuts_R.append(cut(R, x_m, Rm_scale, [1.0, value]))
        least.append(float(np.min(R(np.linspace(1e-3, 3, SAMPLES), 1.0, value))))
        if value > 1:
            ck(f"{name} g^rr is positive at every r", np.all(grr(np.linspace(0, 3, SAMPLES), 1.0, value) > 0))
    ck("boulware_deser the movie's horizon moves in as l grows",
       all(b["marks"][0]["x"] < a["marks"][0]["x"] for a, b in zip(frames, frames[1:]) if b["marks"]))
    ck("boulware_deser the movie's K and R leave their panels nearer the centre as l grows, never farther",
       all(c is not None for c in cuts_K + cuts_R)
       and all(b <= a for a, b in zip(cuts_K, cuts_K[1:])) and all(b <= a for a, b in zip(cuts_R, cuts_R[1:]))
       and cuts_K[-1] < cuts_K[0] and cuts_R[-1] < cuts_R[0],
       f"K {cuts_K[0]:.3f} to {cuts_K[-1]:.3f}, R {cuts_R[0]:.3f} to {cuts_R[-1]:.3f}")
    ck("boulware_deser the movie's least R rises toward zero as l grows", all(b > a for a, b in zip(least, least[1:])),
       f"{least[0]:.3f} to {least[-1]:.4f}")
    cut_tm = cut(K_t, x_m, Km_scale, [1.0])
    ck("boulware_deser the movie's Tangherlini R stays inside its panel", cut(R_t, x_m, Rm_scale, [1.0]) is None)
    first, last = float(BD_STEP), float(BD_LAST)
    movie = {
        "id": "length", "label": f"$r_0/{int(1 / BD_STEP)} \\le \\ell \\le 3r_0/2$", "kind": "line", "x": x_m,
        "panels": panels_m, "curves": references, "marks": [], "legend": legend,
        "caption": [
            "The Kretschmann scalar $K$ and the Ricci scalar $R$ of Boulware and Deser's black hole against $r$ as "
            f"the Gauss-Bonnet length grows at fixed mass ($r_0/{int(1 / BD_STEP)} \\le \\ell \\le 3r_0/2$), beside "
            "Tangherlini's of the same mass, the black hole at $\\ell = 0$. As $\\ell$ grows, $K$ near the centre, "
            "$12r_0^2/\\ell^2r^4$, falls, and the least value of $R$, $-(10 - 7\\sqrt{2})/\\ell^2$ at "
            "$r = \\sqrt{2\\ell r_0}$, rises toward zero.",
            "The horizon $r_h = \\sqrt{r_0^2 - \\ell^2}$ moves in as $\\ell$ grows and closes at $\\ell = r_0$. Past "
            "it $g^{rr}$ is positive at every $r$, and the singularity at $r = 0$ has no horizon.",
        ],
        "settings": f"$r_0 = 1$, the unit of every length, and $\\ell$ in steps of $r_0/{int(1 / BD_STEP)}$.",
        "clip": ("We draw $K\\,r_0^4$ on a logarithmic scale from $10^{-3}$ to $10^{5}$ and $R\\,r_0^2$ on a signed "
                 "logarithmic scale, linear within $0.01$ of zero, from $-100$ to $10^{4}$, and cut each curve where it "
                 f"leaves its panel: Boulware and Deser's $K$ passes the top inside $r = {cuts_K[0]:.3f}\\,r_0$ at "
                 f"$\\ell = {first:.2f}\\,r_0$ and inside $r = {cuts_K[-1]:.3f}\\,r_0$ at $\\ell = {last:.2f}\\,r_0$, "
                 f"its $R$ inside $r = {cuts_R[0]:.3f}\\,r_0$ and $r = {cuts_R[-1]:.3f}\\,r_0$, and Tangherlini's $K$ "
                 f"inside $r = {cut_tm:.3f}\\,r_0$."),
        "movie": {"variable": "$\\ell$", "seconds": 5, "frames": frames},
    }
    return [hole, branch, movie]


# ---------------------------------------------------------------- Brill's charged Taub-NUT

BRILL_HOLE = {"m": Fraction(1), "l": Fraction(3, 4), "r_q": Fraction(1)}         # in units of m, as every drawing of the black hole takes it
BRILL_WORMHOLE = {"m": Fraction(0), "l": Fraction(1), "r_q": Fraction(3, 2)}     # in units of l, as every drawing of the wormhole takes it
BRILL_SCALE = {"kind": "symlog", "from": -1e3, "to": 1e3, "linear": 1e-3}


def brill_charged_taub_nut(ck, src):
    """Brill's charged Taub-NUT, its Kretschmann scalar over the whole line of r, for the black
    hole and for the wormhole.

    The published K depends on r alone, Sigma = r^2 + l^2 never vanishes while l is not zero, and
    every chart in r publishes the same K, so one figure serves them all; R = 0 in every chart,
    since the stress tensor of the electromagnetic field is traceless, and the figure is K alone.
    Each view is one panel over the box of the spacetime diagrams, r from -3m to 6m for the black
    hole at l = 3m/4, r_q = m, and from -3l to 3l for the wormhole at m = 0, r_q = 3l/2, K on a
    signed logarithmic axis, linear within 10^-3 of zero, from -10^3 to 10^3, since K changes
    sign. Reissner and Nordstrom's K of the same mass and charge, read from rn_metric.json at
    r_s = 2m, is the comparison, drawn for r > 0 alone, where its spacetime lies, and the horizons
    are the zeros of the published g^rr.

    What the captions state is checked here: the published K at l = 0 is Reissner and
    Nordstrom's; Brill's K is finite at every r drawn and stays inside its panel; K at r = 0 is
    8(r_q^4 + 6(l^2 - r_q^2)^2 - 6m^2l^2)/l^8; the black hole's horizons are r_- = m/4 and
    r_+ = 7m/4, its greatest and least K lie inside r_-, at the radii and values stated, and K
    changes sign five times on the range drawn, the last time between the horizons; the wormhole has no horizon, its K is even in r,
    greatest at the throat, and changes sign twice, at the radius stated, and is least at the
    radius and value stated; far out K falls off as 48(m^2 - l^2)/r^6; and Reissner and
    Nordstrom's K, of an extremal black hole beside the black hole, falls off as 48m^2/r^6, is
    56r_q^4/r^8 beside the wormhole, and leaves the top of its panel once, inside the radius the
    clip states."""
    K = nr_published(src, "brill_charged_taub_nut", "spherical", "kretschmann", ["r"], ["m", "l", "r_q"])
    grr = nr_published(src, "brill_charged_taub_nut", "spherical", "inverse_metric_components", ["r"],
                       ["m", "l", "r_q"], ("r", "r"))
    K_rn = nr_published(src, "rn_metric", "spherical", "kretschmann", ["r"], ["r_s", "r_q"])
    scale, aspect = BRILL_SCALE, 1.5
    shape_checks(ck, "brill_charged_taub_nut", [aspect])
    height = PLOT_W / aspect
    near = np.linspace(0.05, 6, 41)
    ck("brill_charged_taub_nut K at l = 0 is Reissner and Nordstrom's",
       np.allclose(K(near, 1, 0.0, 1), K_rn(near, 2, 1), rtol=1e-13, atol=0)
       and np.allclose(K(near, 0, 0.0, 1.5), K_rn(near, 0, 1.5), rtol=1e-13, atol=0))

    def view(case, unit, lo, hi, name):
        params = {key: float(value) for key, value in case.items()}
        m, l, r_q = params["m"], params["l"], params["r_q"]
        x_axis = axis(f"$r/{unit}$", lo, hi)
        panel = {"id": "K", "name": f"$K\\,{unit}^4$", "aspect": aspect, "scale": scale, "ticks": value_ticks(scale)}
        own = line_curves(ck, name, "K", "K", K, "r", x_axis, scale, height, params)
        reference = line_curves(ck, f"{name} Reissner and Nordstrom", "K", "reference", K_rn, "r", x_axis, scale,
                                height, {"r_s": 2 * m, "r_q": r_q}, within=(0, hi))
        marks = [{"class": "horizon", "x": x} for x in zeros(ck, f"{name} horizons", grr, lo, hi, [m, l, r_q])]
        rs = np.linspace(lo, hi, SAMPLES)
        ck(f"{name} K is finite at every r drawn", np.all(np.isfinite(K(rs, m, l, r_q))))
        up, down = beyond(K, x_axis, scale, [m, l, r_q])
        ck(f"{name} K stays inside its panel", not up and not down, f"{up} {down}")
        above, below = beyond(K_rn, x_axis, scale, [2 * m, r_q], within=(0, hi))
        ck(f"{name} Reissner and Nordstrom's K leaves the top of its panel once, at r = 0",
           len(above) == 1 and above[0][0] == 0 and not below, f"{above} {below}")
        centre = float(K.exact([0, m, l, r_q]))
        closed = 8 * (r_q ** 4 + 6 * (l ** 2 - r_q ** 2) ** 2 - 6 * m ** 2 * l ** 2) / l ** 8
        ck(f"{name} K at r = 0 is 8(r_q^4 + 6(l^2 - r_q^2)^2 - 6m^2l^2)/l^8", abs(centre / closed - 1) < 1e-12,
           f"{centre}")
        far = 1e6
        ck(f"{name} K falls off as 48(m^2 - l^2)/r^6 far out",
           abs(float(K.exact([far, m, l, r_q])) * far ** 6 / (48 * (m ** 2 - l ** 2)) - 1) < 1e-4)
        fine = np.linspace(lo, hi, 900001)
        values = K(fine, m, l, r_q)
        top, bottom = int(np.argmax(values)), int(np.argmin(values))
        extremes = (float(fine[top]), float(values[top]), float(fine[bottom]), float(values[bottom]))
        ck(f"{name} K is greatest and least at the samples kept", values[top] >= np.max(K(rs, m, l, r_q))
           and values[bottom] <= np.min(K(rs, m, l, r_q)),
           "greatest {1:.4f} at r = {0:.4f}, least {3:.4f} at r = {2:.4f}".format(*extremes))
        sign = zeros(ck, f"{name} K", K, lo, hi, [m, l, r_q])
        cut = above[0][1]
        view = {
            "kind": "line", "x": x_axis, "panels": [panel], "curves": reference + own, "marks": marks,
            "clip": ("We draw $K\\," f"{unit}^4$ on a signed logarithmic scale, linear within $0.001$ of zero, from "
                     "$-10^{3}$ to $10^{3}$, and cut each curve where it leaves its panel, as Reissner and "
                     f"Nordström's $K$ does inside $r = {cut:.3f}\\,{unit}$."),
        }
        return view, centre, extremes, sign, marks

    hole, centre, (r_top, K_top, r_bottom, K_bottom), sign, marks = view(BRILL_HOLE, "m", -3, 6,
                                                                         "brill_charged_taub_nut black hole")
    ck("brill_charged_taub_nut black hole has its horizons at m/4 and 7m/4",
       len(marks) == 2 and abs(marks[0]["x"] - 0.25) < 1e-12 and abs(marks[1]["x"] - 1.75) < 1e-12,
       f"{[m['x'] for m in marks]}")
    ck("brill_charged_taub_nut black hole K is greatest and least inside r_-", max(r_top, r_bottom) < 0.25)
    ck("brill_charged_taub_nut black hole K changes sign five times", len(sign) == 5, f"{sign}")
    ck("brill_charged_taub_nut black hole K changes sign last between the horizons", 0.25 < sign[-1] < 1.75)
    m = float(BRILL_HOLE["m"])
    ck("brill_charged_taub_nut black hole Reissner and Nordstrom's of the same charge is extremal",
       BRILL_HOLE["r_q"] == BRILL_HOLE["m"])
    ck("brill_charged_taub_nut black hole Reissner and Nordstrom's K falls off as 48m^2/r^6 far out",
       abs(float(K_rn.exact([1e6, 2 * m, float(BRILL_HOLE["r_q"])])) * 1e36 / (48 * m ** 2) - 1) < 1e-4)
    hole.update({
        "id": "black_hole", "label": "The black hole",
        "legend": [["line", "K", "$K$, the Kretschmann scalar"],
                   ["line", "reference", "Reissner and Nordström's black hole of the same mass and charge"],
                   ["line", "horizon", "a horizon, where $g^{rr} = 0$"]],
        "caption": [
            "The Kretschmann scalar $K$ of Brill's black hole against $r$ ($l = 3m/4$, $r_q = m$), the same at "
            "every $t$, $\\theta$, and $\\phi$, beside Reissner and Nordström's of the same mass and charge, their "
            "extremal black hole. The stress tensor of the electromagnetic field is traceless, so the Ricci "
            "scalar vanishes. Toward $r = 0$ Reissner and Nordström's $K$ grows without bound. Brill's "
            "$\\Sigma = r^2 + l^2$ is never less than $l^2$, his $K$ is finite at every $r$, and inside $r_-$ the "
            "radius runs on through $r = 0$ to negative values, a second region far from the hole.",
            "At $r = 0$, where $\\Sigma = l^2$, $K$ is $8(r_q^4 + 6(l^2 - r_q^2)^2 - 6m^2l^2)/l^8$, which is "
            f"${centre:.1f}/m^4$ here. Inside $r_-$, $K$ is "
            f"greatest, ${K_top:.1f}/m^4$, at $r = {r_top:.3f}\\,m$ and least, ${K_bottom:.1f}/m^4$, at "
            f"$r = {r_bottom:.3f}\\,m$. It changes sign five times, the first at $r = {sign[0]:.3f}\\,m$ and the "
            f"last at $r = {sign[-1]:.3f}\\,m$, between the horizons. Far out it falls off as $48(m^2 - l^2)/r^6$, and Reissner and "
            "Nordström's as $48m^2/r^6$.",
        ],
        "settings": "$m = 1$, the unit of every length, $l = 3m/4$, and $r_q = m$, so that $r_+ = 7m/4$ and "
                    "$r_- = m/4$.",
    })

    worm, centre, (r_top, K_top, r_bottom, K_bottom), sign, marks = view(BRILL_WORMHOLE, "l", -3, 3,
                                                                         "brill_charged_taub_nut wormhole")
    rs = np.linspace(0, 3, SAMPLES)
    params = [float(v) for v in BRILL_WORMHOLE.values()]
    ck("brill_charged_taub_nut wormhole has no horizon", not marks and np.all(grr(rs, *params) > 0))
    r_q = float(BRILL_WORMHOLE["r_q"])
    ck("brill_charged_taub_nut wormhole Reissner and Nordstrom's K is 56r_q^4/r^8",
       np.allclose(K_rn(rs[1:], 0, r_q), 56 * r_q ** 4 / rs[1:] ** 8, rtol=1e-13, atol=0))
    ck("brill_charged_taub_nut wormhole K is even in r", np.allclose(K(-rs, *params), K(rs, *params), rtol=1e-13, atol=0))
    ck("brill_charged_taub_nut wormhole K is greatest at the throat", abs(r_top) < 1e-5 and abs(K_top / centre - 1) < 1e-9,
       f"{r_top} {K_top}")
    ck("brill_charged_taub_nut wormhole K changes sign twice, at r = +-r_0",
       len(sign) == 2 and abs(sign[0] + sign[1]) < 1e-12, f"{sign}")
    ck("brill_charged_taub_nut wormhole K is least where the caption says, on either side",
       abs(K(-r_bottom, *params) / K_bottom - 1) < 1e-9 and abs(r_bottom) > sign[1], f"{r_bottom} {K_bottom}")
    worm.update({
        "id": "wormhole", "label": "The wormhole",
        "legend": [["line", "K", "$K$, the Kretschmann scalar"],
                   ["line", "reference", "Reissner and Nordström's spacetime of the same mass and charge"]],
        "caption": [
            "The Kretschmann scalar $K$ of Brill's wormhole against $r$ ($m = 0$, $r_q = 3l/2$), the same at every "
            "$t$, $\\theta$, and $\\phi$, beside Reissner and Nordström's of the same mass and charge. The stress "
            "tensor of the electromagnetic field is traceless, so the Ricci scalar vanishes. With no mass the two "
            "sides of the throat $r = 0$ are mirror images, and $K$ is the same at $r$ and at $-r$.",
            "At the throat, where $\\Sigma = l^2$, $K$ is greatest, $8(r_q^4 + 6(l^2 - r_q^2)^2)/l^8$, which is "
            f"${centre:.1f}/l^4$ here, and there Reissner and Nordström's $56r_q^4/r^8$ grows without bound. "
            f"It changes sign at $r = \\pm {sign[1]:.3f}\\,l$ and "
            f"is least, ${K_bottom:.3f}/l^4$, at $r = \\pm {abs(r_bottom):.3f}\\,l$. Far out it falls off as "
            "$-48l^2/r^6$.",
        ],
        "settings": "$l = 1$, the unit of every length, $m = 0$, and $r_q = 3l/2$, so that $\\Delta = r^2 + 5l^2/4$ "
                    "has no root.",
    })
    order = ("id", "label", "kind", "x", "panels", "curves", "marks", "legend", "caption", "settings", "clip")
    return [{key: v[key] for key in order} for v in (hole, worm)]


# ---------------------------------------------------------------- Born and Infeld's point charge

BI_STEPS = 30                     # the steps of r_s from Hoffmann's particle to the black hole in the movie
BI_NUMBER = math.gamma(0.25) ** 2 / (6 * math.sqrt(math.pi))   # Born and Infeld's 1.2360


def born_infeld_charge(ck, src):
    """Born and Infeld's point charge at r_0 = 1 and r_q = r_0/2, its Kretschmann and Ricci
    scalars against r.

    The static chart's published K and R depend on r alone, and every chart publishes the same
    two, so one figure serves all three charts. Two panels share the axis of r, from 0 to 4 r_0,
    the box of the spacetime diagrams: K on a logarithmic axis from 10^-4 to 10^6, and R, which is
    positive at every r, on a logarithmic axis from 10^-6 to 10^4. Reissner and Nordstrom's K of
    the same mass and charge, read from rn_metric.json, is the comparison; their R vanishes, so it
    has no curve on a logarithmic axis. R depends on r_q and r_0 alone, and is read once with those
    two parameters, so that the reading stops if it ever depends on the mass; each curve of R
    still carries the frame's r_s, as every curve of a frame does. The horizons are the
    zeros of the static chart's published g^rr. The first view is Hoffmann's particle, whose whole
    mass is the energy of its field, r_s = 2 (1.2360) r_q^2/r_0, and the second the black hole at
    r_s = 2 r_0, the two cases every other drawing of it takes; the third steps r_s in equal steps
    from the first to the second.

    What the captions state is checked here: for the particle, g^rr tends to 1 - 2 r_q^2/r_0^2 at
    the centre, so m vanishes there as r_q^2 r/r_0^2, g^rr is positive at every r drawn, and K r^4
    tends to 16 r_q^4/r_0^4; for the black hole, K r^6 tends to 48 m_0^2 with m_0 = r_s/2 - 1.2360
    r_q^2/r_0, and g^rr vanishes once, at r_h; Reissner and Nordstrom's K r^8 tends to 56 r_q^4 and
    their R vanishes; R r^2 tends to 4 r_q^2/r_0^2 at the centre and R r^8 to r_q^2 r_0^4 far out;
    K falls from the centre outward and approaches Reissner and Nordstrom's far out; in the movie
    the horizon forms at the centre as soon as r_s passes the particle's and moves out to r_h,
    one in every frame after the first; and every radius a clip states is where its curve
    leaves the top of its panel."""
    rq, r0 = 0.5, 1.0
    particle = float(nr.BI_PARTICLE["r_s"])
    hole = float(nr.BI_HOLE["r_s"])
    x_axis = axis("$r/r_0$", 0, 4)
    K_scale = {"kind": "log", "from": 1e-4, "to": 1e6}
    R_scale = {"kind": "log", "from": 1e-6, "to": 1e4}
    panels = [{"id": "K", "name": "$K\\,r_0^4$", "aspect": 2, "scale": K_scale, "ticks": value_ticks(K_scale)},
              {"id": "R", "name": "$R\\,r_0^2$", "aspect": 2, "scale": R_scale, "ticks": value_ticks(R_scale)}]
    shape_checks(ck, "born_infeld_charge", [p["aspect"] for p in panels])
    height = {p["id"]: PLOT_W / p["aspect"] for p in panels}
    names = ["r_s", "r_q", "r_0"]
    K = nr_published(src, "born_infeld_charge", "static", "kretschmann", ["r"], names)
    R = nr_published(src, "born_infeld_charge", "static", "ricci_scalar", ["r"], names)
    grr = nr_published(src, "born_infeld_charge", "static", "inverse_metric_components", ["r"], names, ("r", "r"))
    K_rn = nr_published(src, "rn_metric", "spherical", "kretschmann", ["r"], ["r_s", "r_q"])
    R_rn = nr_published(src, "rn_metric", "spherical", "ricci_scalar", ["r"], ["r_s", "r_q"])
    # Read with the charge alone, R stops the script if it ever depends on the mass.
    nr_published(src, "born_infeld_charge", "static", "ricci_scalar", ["r"], ["r_q", "r_0"])
    ck("born_infeld_charge R depends on the charge and not on the mass", R.expr.free_symbols.isdisjoint(R.params[:1]),
       f"{R.expr}")
    ck("Reissner and Nordstrom's R vanishes", R_rn.expr == 0, f"{R_rn.expr}")

    def top(fn, params):
        """Where a curve that falls from the centre outward leaves the top of its panel, by
        bisection on the drawn coordinate."""
        scale = K_scale if fn in (K, K_rn) else R_scale
        up, _ = beyond(fn, x_axis, scale, params)
        ck(f"{fn.metric} {fn.field} at {params} leaves the top of its panel once, at the centre",
           len(up) == 1 and up[0][0] == 0, f"{up}")
        lo, hi = 0.0, float(x_axis["to"])
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if mid in (lo, hi):
                break
            lo, hi = (mid, hi) if float(fn(mid, *params)) > scale["to"] else (lo, mid)
        ck(f"{fn.metric} {fn.field} at {params} leaves its panel where the samples do",
           up[0][1] <= hi <= up[0][1] + (x_axis["to"] - x_axis["from"]) / (SAMPLES - 1), f"{hi:.6f}")
        return hi

    def own(rs, name):
        params = {"r_s": float(rs), "r_q": rq, "r_0": r0}
        curves = (line_curves(ck, name, "K", "reference", K_rn, "r", x_axis, K_scale, height["K"],
                              {"r_s": float(rs), "r_q": rq})
                  + line_curves(ck, name, "K", "K", K, "r", x_axis, K_scale, height["K"], params)
                  + line_curves(ck, name, "R", "R", R, "r", x_axis, R_scale, height["R"], params))
        marks = [{"class": "horizon", "x": x} for x in zeros(ck, f"{name} horizons", grr, 0, 4, [float(rs), rq, r0])]
        rs_grid = np.linspace(0, 4, SAMPLES)[1:]
        ck(f"{name} K falls from the centre outward", np.all(np.diff(K(rs_grid, float(rs), rq, r0)) < 0))
        ck(f"{name} K stays above the bottom of its panel", not beyond(K, x_axis, K_scale, [float(rs), rq, r0])[1]
           and not beyond(K_rn, x_axis, K_scale, [float(rs), rq])[1])
        far = 1e3
        ck(f"{name} K approaches Reissner and Nordstrom's far out",
           abs(float(K.exact([far, rs, rq, r0])) / float(K_rn.exact([far, rs, rq])) - 1) < 1e-5)
        return curves, marks

    # R, the same for every mass, checked at the particle's.
    ck("born_infeld_charge R is positive at every r drawn",
       np.all(R(np.linspace(0, 4, SAMPLES)[1:], particle, rq, r0) > 0))
    ck("born_infeld_charge R stays above the bottom of its panel", not beyond(R, x_axis, R_scale, [particle, rq, r0])[1])
    centre_R = float(R.exact([1e-5, particle, rq, r0])) * 1e-10
    ck("born_infeld_charge R r^2 tends to 4 r_q^2/r_0^2 at the centre", abs(centre_R / (4 * rq ** 2 / r0 ** 2) - 1) < 1e-6,
       f"{centre_R}")
    far_R = float(R.exact([1e3, particle, rq, r0])) * 1e24
    ck("born_infeld_charge R r^8 tends to r_q^2 r_0^4 far out", abs(far_R / (rq ** 2 * r0 ** 4) - 1) < 1e-5, f"{far_R}")
    R_cut = top(R, [particle, rq, r0])
    ck("born_infeld_charge 1.2360 is Born and Infeld's number", f"{BI_NUMBER:.4f}" == "1.2360", f"{BI_NUMBER}")
    ck("born_infeld_charge Hoffmann's particle has r_s = 2 (1.2360) r_q^2/r_0",
       abs(particle / (2 * BI_NUMBER * rq ** 2 / r0) - 1) < 1e-15, f"{particle}")
    K_rn_centre = [float(K_rn.exact([1e-4, rs, rq])) * 1e-32 / (56 * rq ** 4) for rs in (particle, hole)]
    ck("Reissner and Nordstrom's K r^8 tends to 56 r_q^4 at the centre", all(abs(v - 1) < 1e-3 for v in K_rn_centre),
       f"{K_rn_centre}")
    legend_lines = [["line", "K", "$K$, the Kretschmann scalar"], ["line", "R", "$R$, the Ricci scalar"],
                    ["line", "reference", "Reissner and Nordström's $K$ for the same mass and charge"]]
    horizon_legend = ["line", "horizon", "a horizon, where $g^{rr} = 0$"]

    # Hoffmann's particle.
    curves, marks = own(particle, "born_infeld_charge particle")
    ck("born_infeld_charge Hoffmann's particle has no horizon", not marks
       and np.all(grr(np.linspace(0, 4, SAMPLES)[1:], particle, rq, r0) > 0))
    cone = float(grr.exact([1e-8, particle, rq, r0]))
    ck("born_infeld_charge g^rr tends to 1 - 2 r_q^2/r_0^2 at the particle's centre, so m vanishes as r_q^2 r/r_0^2",
       abs(cone - (1 - 2 * rq ** 2 / r0 ** 2)) < 1e-6, f"{cone}")
    growth = [float(K.exact([r, particle, rq, r0])) * r ** 4 / (16 * rq ** 4 / r0 ** 4) for r in (1e-4, 1e-5)]
    ck("born_infeld_charge K r^4 tends to 16 r_q^4/r_0^4 at the particle's centre", all(abs(v - 1) < 1e-6 for v in growth),
       f"{growth}")
    K_cut, K_rn_cut = top(K, [particle, rq, r0]), top(K_rn, [particle, rq])
    still_particle = {
        "id": "particle", "label": "Hoffmann's particle", "kind": "line", "x": x_axis, "panels": panels,
        "curves": curves, "marks": marks, "legend": legend_lines,
        "caption": [
            "The Kretschmann scalar $K$ and the Ricci scalar $R$ of Hoffmann's particle against $r$ "
            "($r_q = r_0/2$, $r_s = 0.618\\,r_0$), the same at every $t$, $\\theta$, and $\\phi$, beside the $K$ of "
            "Reissner and Nordström's metric of the same mass and charge. The whole mass is the energy of the field, "
            "so $m$ vanishes at the centre, as $r_q^2r/r_0^2$, and $K$ grows toward it as $16r_q^4/(r_0^4r^4)$, "
            "where Reissner and Nordström's grows as $56r_q^4/r^8$.",
            "Reissner and Nordström's $R$ vanishes, since Maxwell's field has a traceless stress tensor. Born and "
            "Infeld's $R$ is positive at every $r$, grows toward the centre as $4r_q^2/(r_0^2r^2)$, and falls off "
            "as $r_q^2r_0^4/r^8$ far out, where $K$ approaches Reissner and Nordström's.",
        ],
        "settings": "$r_0 = 1$, the unit of every length, $r_q = r_0/2$, and $r_s = 0.618\\,r_0$, the mass of the "
                    "field alone.",
        "clip": ("We draw $K\\,r_0^4$ on a logarithmic scale from $10^{-4}$ to $10^{6}$ and $R\\,r_0^2$ on one from "
                 "$10^{-6}$ to $10^{4}$, and cut each curve where it leaves its panel, as the particle's $K$ does "
                 f"inside $r = {K_cut:.3f}\\,r_0$, its $R$ inside $r = {R_cut:.3f}\\,r_0$, and Reissner and "
                 f"Nordström's $K$ inside $r = {K_rn_cut:.3f}\\,r_0$."),
    }

    # The black hole.
    curves, marks = own(hole, "born_infeld_charge hole")
    ck("born_infeld_charge the black hole has one horizon, at r_h", len(marks) == 1
       and abs(marks[0]["x"] - nr.BI_HORIZON) < 1e-12, f"{[m['x'] for m in marks]}")
    r_h = marks[0]["x"]
    m0 = hole / 2 - BI_NUMBER * rq ** 2 / r0
    centre = [float(K.exact([r, hole, rq, r0])) * r ** 6 / (48 * m0 ** 2) for r in (1e-5, 1e-6)]
    ck("born_infeld_charge K r^6 tends to 48 m_0^2 at the black hole's centre", all(abs(v - 1) < 1e-4 for v in centre),
       f"m_0 = {m0:.6f}, {centre}")
    K_cut, K_rn_cut = top(K, [hole, rq, r0]), top(K_rn, [hole, rq])
    still_hole = {
        "id": "hole", "label": "The black hole", "kind": "line", "x": x_axis, "panels": panels,
        "curves": curves, "marks": marks, "legend": legend_lines + [horizon_legend],
        "caption": [
            "The Kretschmann scalar $K$ and the Ricci scalar $R$ of the black hole against $r$ "
            "($r_q = r_0/2$, $r_s = 2\\,r_0$), the same at every $t$, $\\theta$, and $\\phi$, beside the $K$ of "
            "Reissner and Nordström's metric of the same mass and charge. With more mass than the field holds, the "
            f"mass at the centre is $m_0 = r_s/2 - 1.2360\\,r_q^2/r_0 = {m0:.3f}\\,r_0$, and toward it $K$ grows as "
            "$48m_0^2/r^6$, as Schwarzschild's does for the mass $m_0$, where Reissner and Nordström's grows as "
            "$56r_q^4/r^8$.",
            f"The horizon is at $r_h = {r_h:.3f}\\,r_0$, where $g^{{rr}} = 1 - 2m/r$ vanishes. $R$ is a function "
            "of $r$, $r_q$, and $r_0$ alone, so it is the particle's, and far out $K$ approaches Reissner and "
            "Nordström's.",
        ],
        "settings": f"$r_0 = 1$, the unit of every length, $r_q = r_0/2$, and $r_s = 2\\,r_0$, so that the horizon "
                    f"is at $r_h = {r_h:.3f}\\,r_0$.",
        "clip": ("We draw $K\\,r_0^4$ on a logarithmic scale from $10^{-4}$ to $10^{6}$ and $R\\,r_0^2$ on one from "
                 "$10^{-6}$ to $10^{4}$, and cut each curve where it leaves its panel, as the black hole's $K$ does "
                 f"inside $r = {K_cut:.3f}\\,r_0$, its $R$ inside $r = {R_cut:.3f}\\,r_0$, and Reissner and "
                 f"Nordström's $K$ inside $r = {K_rn_cut:.3f}\\,r_0$."),
    }

    # The movie in r_s, at fixed charge, from Hoffmann's particle to the black hole.
    frames, horizons, K_cuts, K_rn_cuts = [], [], [], []
    for k in range(BI_STEPS + 1):
        rs = particle + (hole - particle) * k / BI_STEPS if k < BI_STEPS else hole
        curves, marks = own(rs, f"born_infeld_charge at r_s = {rs:.4f}")
        frames.append({"label": f"$r_s = {rs:.3f}\\,r_0$", "value": float(rs), "curves": curves, "marks": marks})
        horizons.append([m["x"] for m in marks])
        K_cuts.append(top(K, [rs, rq, r0]))
        K_rn_cuts.append(top(K_rn, [rs, rq]))
    ck("born_infeld_charge has no horizon in the first frame and one in every other",
       horizons[0] == [] and all(len(h) == 1 for h in horizons[1:]), f"{horizons}")
    outward = [h[0] for h in horizons[1:]]
    ck("born_infeld_charge the horizon moves outward as r_s grows", all(b > a for a, b in zip(outward, outward[1:])))
    ck("born_infeld_charge the last frame is the black hole", frames[-1]["value"] == hole and outward[-1] == r_h)
    just_past = particle * (1 + 1e-6)
    newborn = zeros(ck, "born_infeld_charge horizon just past the particle", grr, 1e-12, 0.01, [just_past, rq, r0])
    ck("born_infeld_charge the horizon forms at the centre as r_s passes the particle's",
       len(newborn) == 1 and newborn[0] < 1e-4, f"{newborn}")
    sixth = []
    for frame in frames[1:]:
        rs = frame["value"]
        m_c = rs / 2 - BI_NUMBER * rq ** 2 / r0
        sixth.append(float(K.exact([1e-6, rs, rq, r0])) * 1e-36 / (48 * m_c ** 2))
    ck("born_infeld_charge K grows as 48 m_0^2/r^6 at the centre in every frame after the first",
       all(abs(v - 1) < 1e-3 for v in sixth), f"worst {max(abs(v - 1) for v in sixth):.1e}")
    ck("born_infeld_charge the radius where K leaves the top grows from frame to frame",
       all(b > a for a, b in zip(K_cuts, K_cuts[1:])), f"{K_cuts[0]:.4f} to {K_cuts[-1]:.4f}")
    movie_view = {
        "id": "mass", "label": "$0.618\\,r_0 \\le r_s \\le 2\\,r_0$", "kind": "line", "x": x_axis,
        "panels": panels, "curves": [], "marks": [], "legend": legend_lines + [horizon_legend],
        "caption": [
            "The Kretschmann scalar $K$ and the Ricci scalar $R$ of Born and Infeld's point charge against $r$ as its "
            "mass grows at fixed charge ($r_q = r_0/2$, $0.618\\,r_0 \\le r_s \\le 2\\,r_0$), beside the $K$ of "
            "Reissner and Nordström's metric of the same mass and charge. At $r_s = 0.618\\,r_0$ the whole mass is "
            "the energy of the field, and $K$ grows toward the centre as $16r_q^4/(r_0^4r^4)$. With any more mass "
            "the mass at the centre, $m_0 = r_s/2 - 1.2360\\,r_q^2/r_0$, is positive, and $K$ grows there as "
            "$48m_0^2/r^6$.",
            "A horizon forms at the centre as soon as $r_s$ passes $0.618\\,r_0$, and moves out to "
            f"$r_h = {r_h:.3f}\\,r_0$ at $r_s = 2\\,r_0$. $R$ is a function of $r$, $r_q$, and $r_0$ alone, so "
            "it is the same in every frame.",
        ],
        "settings": f"$r_0 = 1$, the unit of every length, $r_q = r_0/2$, and $r_s$ in {BI_STEPS} equal steps from "
                    "$0.618\\,r_0$ to $2\\,r_0$.",
        "clip": ("The scales are those of Hoffmann's particle and the black hole, and we cut each curve where it "
                 "leaves its panel: Born and Infeld's $K$ passes the top inside a radius that grows from "
                 f"$r = {K_cuts[0]:.3f}\\,r_0$ in the first frame to ${K_cuts[-1]:.3f}\\,r_0$ in the last, "
                 f"Reissner and Nordström's inside a radius between ${min(K_rn_cuts):.3f}\\,r_0$ and "
                 f"${max(K_rn_cuts):.3f}\\,r_0$, "
                 f"and $R$ inside $r = {R_cut:.3f}\\,r_0$ in every frame."),
        "movie": {"variable": "$r_s$", "seconds": 5, "frames": frames},
    }
    return [still_particle, still_hole, movie_view]


# ---------------------------------------------------------------- the C-metric

CM_ALPHA = Fraction(1, 6)         # alpha in units of 1/m, as every drawing of the C-metric takes it
CM_REACH = 8                      # the box runs out to r sin theta = 8 m and up to r cos theta = 8 m
CM_CELLS = 20                     # the cells of the map to each m, across and up
CM_RANGE = (1e-6, 1e2)            # the colours of K m^4 run from the first to the second


def c_metric(ck, src):
    """The C-metric at m = 1, alpha = 1/(6m) and C = 1/(1 + 2 alpha m) = 3/4, its Kretschmann
    scalar over the meridian half plane of the spherical chart's r and theta.

    The spherical chart publishes K = 48 m^2 (1 + alpha r cos theta)^6/r^6, a function of r and
    theta alone, the same at every t and phi, and R = 0, so the figure is K alone, a map with r
    and theta drawn as polar coordinates about r = 0, r sin theta across and r cos theta up, over
    0 <= r sin theta <= 8 m and -1/alpha <= r cos theta <= 8 m. The chart ends below on the line
    r cos theta = -1/alpha, null infinity, where 1 + alpha r cos theta vanishes, and the box ends
    there with it; it reaches past the acceleration horizon r = 1/alpha = 6 m, out to the 8 m the
    spacetime diagram of the inner axis reaches. The cells are a twentieth of m on a side,
    coloured on a logarithmic scale from 10^-6 to 10^2 in units of 1/m^4. The horizons r = 2m and
    r = 1/alpha, the zeros of the published g^rr, are marked as semicircles, the singularity
    r = 0 as a point, and the cosmic string along the half axis theta = pi as a line.

    What the captions state is checked here: R vanishes; at alpha = 0, K is Schwarzschild's
    published scalar of the same mass; on the black hole horizon K falls from 1024/243 m^4 at
    theta = 0 to 16/243 m^4 at theta = pi, 64 times less; K is 48 alpha^6 m^2 (d/r)^6 with
    d = r cos theta + 1/alpha, so that every level is a conic of focus r = 0 and directrix d = 0,
    the level 48 alpha^6 m^2 is the parabola r = d, and K tends to that level far out along
    theta = 0; K diverges as 48 m^2/r^6 toward r = 0 from every direction; the Hong-Teo chart's
    published K is the same scalar; with C = 3/4 a small circle about the half axis theta = 0 has
    2 pi times its radius as its circumference and one about theta = pi half that, the cosmic
    string, along which K is finite; and the clip stated is where K leaves the colours, on a grid
    four times as fine as the map's, against where the conics of the two levels reach."""
    m, alpha = 1.0, float(CM_ALPHA)
    C = float(1 / (1 + 2 * CM_ALPHA))
    params = {"m": m, "\\alpha": alpha}
    K = nr_published(src, "c_metric", "spherical", "kretschmann", ["r", "theta"], ["m", "\\alpha"])
    R = nr_published(src, "c_metric", "spherical", "ricci_scalar", ["r", "theta"], ["m", "\\alpha"])
    grr = nr_published(src, "c_metric", "spherical", "inverse_metric_components", ["r", "theta"], ["m", "\\alpha"],
                       ("r", "r"))
    gthth = nr_published(src, "c_metric", "spherical", "metric_components", ["r", "theta"], ["m", "\\alpha"],
                         ("\\theta", "\\theta"))
    gphph = nr_published(src, "c_metric", "spherical", "metric_components", ["r", "theta"], ["m", "\\alpha", "C"],
                         ("\\phi", "\\phi"))
    K_ht = nr_published(src, "c_metric", "hong_teo", "kretschmann", ["x", "y"], ["m", "\\alpha"])
    K_s = nr_published(src, "schwarzschild", "spherical", "kretschmann", ["r"], ["r_s"])
    bottom = float(-1 / CM_ALPHA)
    # Ticks every 2 m on both axes, as the spacetime diagram of the inner axis has them.
    x_axis, y_axis = axis("$r\\sin\\theta/m$", 0, CM_REACH), axis("$r\\cos\\theta/m$", bottom, CM_REACH, 8)
    aspect = (x_axis["to"] - x_axis["from"]) / (y_axis["to"] - y_axis["from"])
    shape_checks(ck, "c_metric", [aspect])
    scale = {"kind": "log", "from": CM_RANGE[0], "to": CM_RANGE[1]}
    nx = round(CM_CELLS * (x_axis["to"] - x_axis["from"]))
    ny = round(CM_CELLS * (y_axis["to"] - y_axis["from"]))
    xs, ys, values = cells(K, x_axis, y_axis, nx, ny, [m, alpha], frame="meridian")
    ck("c_metric K is finite and positive at every cell", np.all(np.isfinite(values)) and np.all(values > 0),
       f"from {np.min(values):.3e} to {np.max(values):.3e}")
    # Sixty digits at forty cells spread over the map, at the two beside the singularity and at
    # four along null infinity.
    rng = np.random.default_rng(1970)
    picks = {(int(i), int(j)) for i, j in zip(rng.integers(0, nx, 40), rng.integers(0, ny, 40))}
    centre = round(-bottom * CM_CELLS)
    picks |= {(0, centre - 1), (0, centre)} | {(i, 0) for i in (0, nx // 3, 2 * nx // 3, nx - 1)}
    picks = sorted(picks)
    points = [[float(v) for v in meridian(xs[i], ys[j])] + [m, alpha] for i, j in picks]
    agrees(ck, "c_metric map", K, points, [values[j, i] for i, j in picks])

    r, th = K.coords
    mm, aa = K.params
    ck("c_metric R vanishes", sp.simplify(R.expr) == 0, f"R = {R.expr}")
    schwarzschild = K_s.expr.subs({K_s.params[0]: 2 * mm, K_s.coords[0]: r})
    ck("c_metric K at alpha = 0 is Schwarzschild's of the same mass", sp.simplify(K.expr.subs(aa, 0) - schwarzschild) == 0,
       f"{schwarzschild}")
    exact = {mm: 1, aa: sp.Rational(CM_ALPHA.numerator, CM_ALPHA.denominator), r: 2}
    toward, beyond_ = K.expr.subs({**exact, th: 0}), K.expr.subs({**exact, th: sp.pi})
    ck("c_metric K on the black hole horizon is 1024/243 m^4 at theta = 0 and 16/243 m^4 at theta = pi",
       toward == sp.Rational(1024, 243) and beyond_ == sp.Rational(16, 243), f"{toward}, {beyond_}")
    ck("c_metric K on the black hole horizon at theta = 0 is ((1 + 2 alpha m)/(1 - 2 alpha m))^6 = 64 times that at pi",
       toward / beyond_ == ((1 + 2 * exact[aa]) / (1 - 2 * exact[aa])) ** 6 == 64)
    around = np.linspace(0, math.pi, SAMPLES)
    ck("c_metric K falls along the black hole horizon from theta = 0 to theta = pi",
       np.all(np.diff(K(2.0, around, m, alpha)) < 0))
    focus = 48 * aa ** 6 * mm ** 2 * ((r * sp.cos(th) + 1 / aa) / r) ** 6
    ck("c_metric K is 48 alpha^6 m^2 (d/r)^6 with d = r cos theta + 1/alpha", sp.simplify(K.expr - focus) == 0)
    level = 48 * alpha ** 6 * m ** 2
    ck("c_metric 48 alpha^6 m^2 is 1/972 m^4", abs(level * 972 - 1) < 1e-15, tex_number(level))
    on_parabola = []
    for z in np.linspace(bottom / 2, CM_REACH, 9):
        radius = z - bottom
        on_parabola.append(float(K.exact([radius, math.acos(z / radius), m, alpha])) / level - 1)
    ck("c_metric K is 48 alpha^6 m^2 on the parabola r = r cos theta + 1/alpha",
       max(map(abs, on_parabola)) < 1e-12, f"worst {max(map(abs, on_parabola)):.1e}")
    far = [float(K.exact([radius, 0, m, alpha])) / level - 1 for radius in (1e4, 1e8, 1e12)]
    ck("c_metric K tends to 48 alpha^6 m^2 from above far out along theta = 0",
       all(v > 0 for v in far) and far[0] > far[1] > far[2] and far[2] < 1e-10, f"{far}")
    near = [float(K.exact([1e-8, angle, m, alpha])) * 1e-48 / 48 - 1 for angle in np.linspace(0, math.pi, 7)]
    ck("c_metric K diverges as 48 m^2/r^6 toward r = 0 from every direction", max(map(abs, near)) < 1e-6,
       f"worst {max(map(abs, near)):.1e}")
    x, y = K_ht.coords
    hong_teo = K_ht.expr.subs({x: sp.cos(th), y: 1 / (aa * r), K_ht.params[0]: mm, K_ht.params[1]: aa})
    ck("c_metric the Hong-Teo chart's K is the same scalar at x = cos theta, y = 1/(alpha r)",
       sp.simplify(hong_teo - K.expr) == 0)

    # The horizons, the zeros of g^rr along rays from r = 0 across the box, at 2m and 1/alpha on
    # every ray but theta = pi, where null infinity ends the ray at 1/alpha.
    rays = []
    for angle in np.linspace(0, math.pi, 13)[:-1]:
        reach = min(CM_REACH / max(math.sin(angle), 1e-300), CM_REACH / math.cos(angle) if math.cos(angle) > 0
                    else bottom / math.cos(angle) if math.cos(angle) < 0 else math.inf)
        # A ray that ends on null infinity stops just short of it, where g^rr has its double zero.
        rays.append(zeros(ck, f"c_metric g^rr at theta = {angle:.4f}", grr, 0.0, reach * (1 - 1e-9), [angle, m, alpha],
                          n=4000))
    radii = rays[0]
    ck("c_metric has its horizons at r = 2m and r = 1/alpha on every ray",
       all(len(z) == 2 and abs(z[0] - 2) < 1e-12 and abs(z[1] - 6) < 1e-12 for z in rays), f"{radii}")
    per_unit = PLOT_W / (x_axis["to"] - x_axis["from"])

    def semicircle(radius):
        """The circle r = radius from theta = 0 to pi, in chords that stray from it by less than
        TOLERANCE of the plot's units."""
        n = math.ceil(math.pi / (2 * math.acos(1 - TOLERANCE / (radius * per_unit)))) + 1
        angles = np.linspace(0, math.pi, n)
        sag = radius * per_unit * (1 - math.cos(math.pi / (n - 1) / 2))
        ck(f"c_metric the circle r = {radius:g} strays from its chords by less than the tolerance", sag < TOLERANCE,
           f"{sag:.4f} in {n} points")
        return [[float(radius * math.sin(a)) + 0.0, float(radius * math.cos(a)) + 0.0] for a in angles]
    horizons = [semicircle(radius) for radius in radii]
    flips = []
    for line in horizons:
        for px, py in line:
            radius, angle = (float(v) for v in meridian(px, py))
            if 1 + alpha * radius * (1 + 1e-7) * math.cos(angle) <= 0:
                continue
            inside, outside = grr(radius * (1 - 1e-7), angle, m, alpha), grr(radius * (1 + 1e-7), angle, m, alpha)
            flips.append(float(inside) * float(outside) < 0)
    ck("c_metric g^rr changes sign across every point of both horizons marked", all(flips), f"{len(flips)} points")

    # The string: a small circle about each half axis, its circumference over 2 pi times the
    # distance out to it along theta, from the published g_theta_theta and g_phi_phi.
    delta = 1e-8

    def turn(pole, toward_):
        angle = pole + toward_ * delta
        return float(sp.sqrt(gphph.exact([1.0, angle, m, alpha, C])) / (delta * sp.sqrt(gthth.exact([1.0, angle, m, alpha]))))
    regular, string = turn(0.0, 1), turn(math.pi, -1)
    ck("c_metric with C = 3/4 the half axis theta = 0 is regular", abs(regular - 1) < 1e-6, f"{regular:.9f}")
    ck("c_metric with C = 3/4 the half axis theta = pi has circles of half the circumference, a string",
       abs(string - 0.5) < 1e-6 and abs((1 - string) * 2 * math.pi - 8 * math.pi * alpha * m / (1 + 2 * alpha * m)) < 1e-5,
       f"{string:.9f}")
    along = K(np.linspace(1e-3, -bottom, SAMPLES)[1:-1], math.pi, m, alpha)
    ck("c_metric K is finite and positive along the string", np.all(np.isfinite(along)) and np.all(along > 0))

    # Where K leaves the colours, on a grid four times as fine as the map's, against the conics
    # of the two levels: the top level an ellipse reaching farthest from r = 0 on theta = 0, and
    # the bottom level a hyperbola reaching highest above null infinity on the right edge.
    fx, fy, fine = cells(K, x_axis, y_axis, 4 * nx, 4 * ny, [m, alpha], frame="meridian")
    FX, FY = np.meshgrid(fx, fy)
    high, low = ~(fine <= CM_RANGE[1]), fine < CM_RANGE[0]
    d_high = float(np.hypot(FX, FY)[high].max())
    d_low = float((FY - bottom)[low].max())
    q_high = (CM_RANGE[1] / level) ** (1 / 6)
    r_high = 1 / (alpha * (q_high - 1))
    q_low = (CM_RANGE[0] / level) ** (1 / 6)
    lo_, hi_ = 0.0, CM_REACH - bottom
    for _ in range(200):
        mid = 0.5 * (lo_ + hi_)
        lo_, hi_ = (mid, hi_) if mid < q_low * math.hypot(CM_REACH, mid + bottom) else (lo_, mid)
    cell = 1 / (4 * CM_CELLS)
    ck("c_metric K passes the top of its colours only near r = 0, inside the ellipse of that level",
       r_high - 2 * cell < d_high <= r_high, f"within r = {d_high:.4f}, the ellipse reaching {r_high:.4f}")
    ck("c_metric K passes the bottom of its colours only near null infinity, below the hyperbola of that level",
       lo_ - 2 * cell < d_low <= lo_, f"within {d_low:.4f}, the hyperbola reaching {lo_:.4f}")
    rows = [[rounded(v) for v in row] for row in values]
    view = {
        "id": "meridian", "label": "the half plane of $r$ and $\\theta$", "kind": "map",
        "x": x_axis, "y": y_axis,
        "colour": {"name": "$K\\,m^4$", "scale": scale, "ticks": value_ticks(scale)},
        "map": {"metric": "c_metric", "system": "spherical", "scalar": "kretschmann", "x": "r", "y": "\\theta",
                "frame": "meridian", "params": params, "nx": nx, "ny": ny, "cells": rows},
        "marks": [{"class": "horizon", "points": line} for line in horizons]
                 + [{"class": "string", "points": [[0.0, 0.0], [0.0, bottom]]},
                    {"class": "singular", "at": [0.0, 0.0]}],
        "legend": [["line", "horizon", "the horizons $r = 2m$ and $r = 1/\\alpha$"],
                   ["line", "string", "the cosmic string, $\\theta = \\pi$"],
                   ["point", "singular", "$r = 0$, where $K$ diverges"]],
        "caption": [
            "The Kretschmann scalar $K$ of the C-metric over the half plane of $r$ and $\\theta$, drawn as polar "
            "coordinates about the black hole ($\\alpha m = 1/6$), the same at every $t$ and $\\phi$. The spacetime "
            "is a vacuum, so its Ricci scalar vanishes. At $\\alpha = 0$ the scalar "
            "$K = 48m^2(1 + \\alpha r\\cos\\theta)^6/r^6$ is Schwarzschild's $48m^2/r^6$, the same in every "
            "direction, and for $\\alpha > 0$ it is larger on the half axis $\\theta = 0$, between the black "
            f"holes, than on the half axis $\\theta = \\pi$. On the black hole horizon it falls from "
            f"${tex_number(float(toward))}/m^4$ at $\\theta = 0$ to ${tex_number(float(beyond_))}/m^4$ at "
            "$\\theta = \\pi$, by the factor $((1 + 2\\alpha m)/(1 - 2\\alpha m))^6 = 64$.",
            "In these coordinates every level of $K$ is a conic with its focus at the singularity $r = 0$ and the "
            "line $r\\cos\\theta = -1/\\alpha$ as its directrix, since $K = 48\\alpha^6m^2(d/r)^6$ with "
            "$d = r\\cos\\theta + 1/\\alpha$ the height above that line. The level "
            f"$K = 48\\alpha^6m^2 = {level * 1e3:.2f} \\times 10^{{-3}}/m^4$, which $K$ approaches far out along $\\theta = 0$, "
            "is a parabola, the levels above it are ellipses, and those below it are hyperbolas. On the "
            "directrix $1 + \\alpha r\\cos\\theta$ vanishes, and $K$ with it: the line is null infinity, where "
            "the chart ends.",
            "The black hole horizon $r = 2m$ and the acceleration horizon $r = 1/\\alpha = 6m$, where $g^{rr}$ "
            "vanishes, are circles about $r = 0$, and on the half axis $\\theta = \\pi$ the acceleration horizon lies at null "
            "infinity. With $C = 1/(1 + 2\\alpha m)$ that half axis carries the cosmic string, a conical "
            "singularity along which $K$ stays finite.",
        ],
        "settings": "$m = GM/c^2 = 1$, the unit of every length, $\\alpha = 1/(6m)$, and "
                    "$C = 1/(1 + 2\\alpha m) = 3/4$, which leaves the half axis $\\theta = 0$ regular.",
        "clip": ("We colour $K\\,m^4$ on a logarithmic scale from $10^{-6}$ to $10^{2}$ and give a value beyond it "
                 f"the colour of the end it passes: $K$ passes $10^{{2}}/m^4$ inside $r = {d_high:.2f}\\,m$ "
                 f"and falls below $10^{{-6}}/m^4$ within ${d_low:.2f}\\,m$ of null infinity."),
    }
    return [view]


# ---------------------------------------------------------------- Brans and Dicke's static sphere

BDS_N = Fraction(1, 4)            # n at m = 0, as every drawing of the sphere takes it, which is omega = 6
BDS_STEP = Fraction(1, 100)       # the step of n from one frame of the movie to the next, and its first frame
BDS_CASSINI = 40000               # the least omega the Cassini spacecraft allows, as the page's convention says


def brans_dicke_sphere(ck, src):
    """Brans and Dicke's static sphere at m = 0 and r_0 = 1, its Kretschmann and Ricci scalars
    against r in the spherical chart.

    The spherical chart's published K and R depend on r alone. Two panels share the axis of r,
    from the singularity at 2 r_0 to 8 r_0, the range of the embedding diagram: K on a logarithmic
    axis from 10^-4 to 10^2 and R on one from 10^-6 to 10^2, since both are positive on the range
    and both diverge at r = 2 r_0, where the spheres shrink to a point. Schwarzschild's K of the
    same mass, GM/c^2 = (m + 1) r_0 = r_0, so r_s = 2 r_0, read from schwarzschild.json, is the
    comparison; Schwarzschild's R vanishes and has no place on a logarithmic axis. No horizon is
    marked, since the published g^rr = A^(1 - n) is positive at every r drawn. The first view is
    the sphere at n = 1/4, omega = 6, the values every other drawing of it takes; the second steps
    n from 1/100 to 1/4 at fixed mass, from omega = 198 down to 6.

    What the captions state is checked here: omega = 2(1 - n)/n at m = 0; Schwarzschild's K at
    2 r_0 is 3/4 r_0^4; the published g_thetatheta vanishes at 2 r_0; K and R are positive on the
    range and diverge there as (r - 2 r_0)^-(2 + 2n) and (r - 2 r_0)^-(1 + n), measured in sixty
    digits; K crosses Schwarzschild's once on the range, at the radius stated, and far out tends
    to (n^2 - 2n + 2)/2 of it, with R falling off as 2n(1 - n) r_0^2/r^4; in the movie, K stands
    above twice Schwarzschild's only inside the radii stated, which shrink as n does; and at the
    n of the least omega Cassini allows K is within the fraction stated of Schwarzschild's outside
    the radius stated."""
    x_axis = axis("$r/r_0$", 2, 8)
    K_scale = {"kind": "log", "from": 1e-4, "to": 1e2}
    R_scale = {"kind": "log", "from": 1e-6, "to": 1e2}
    # R's eight decades are ticked at every second one, which the page fits at every size, where
    # nine labels stand as close as it allows and are thinned in one view and not the other.
    R_ticks = [t for t in value_ticks(R_scale) if round(math.log10(t["at"])) % 2 == 0]
    panels = [{"id": "K", "name": "$K\\,r_0^4$", "aspect": 2, "scale": K_scale, "ticks": value_ticks(K_scale)},
              {"id": "R", "name": "$R\\,r_0^2$", "aspect": 2, "scale": R_scale, "ticks": R_ticks}]
    shape_checks(ck, "brans_dicke_sphere", [p["aspect"] for p in panels])
    height = {p["id"]: PLOT_W / p["aspect"] for p in panels}
    names = ["r_0", "m", "n"]
    K = nr_published(src, "brans_dicke_sphere", "spherical", "kretschmann", ["r"], names)
    R = nr_published(src, "brans_dicke_sphere", "spherical", "ricci_scalar", ["r"], names)
    grr = nr_published(src, "brans_dicke_sphere", "spherical", "inverse_metric_components", ["r"], names, ("r", "r"))
    g_thth = nr_published(src, "brans_dicke_sphere", "spherical", "metric_components", ["r"], names, ("\\theta", "\\theta"))
    K_s = nr_published(src, "schwarzschild", "spherical", "kretschmann", ["r"], ["r_s"])
    n = float(BDS_N)

    def omega(n):
        m = 0
        return -2 * (m * m + n * n + m * n + m - n) / (m + n) ** 2
    ck("brans_dicke_sphere omega is 2(1 - n)/n at m = 0, and 6 at n = 1/4",
       omega(BDS_N) == 6 and all(omega(Fraction(k, 100)) == 2 * (1 - Fraction(k, 100)) / Fraction(k, 100)
                                for k in range(1, 26)))

    def reference():
        return line_curves(ck, "brans_dicke_sphere Schwarzschild", "K", "reference", K_s, "r", x_axis, K_scale,
                           height["K"], {"r_s": 2})

    def own(value, name):
        params = {"r_0": 1, "m": 0, "n": float(value)}
        curves = (line_curves(ck, name, "K", "K", K, "r", x_axis, K_scale, height["K"], params)
                  + line_curves(ck, name, "R", "R", R, "r", x_axis, R_scale, height["R"], params))
        ck(f"{name} has no horizon on the range", not zeros(ck, f"{name} g^rr", grr, 2, 8, [1, 0, float(value)])
           and np.all(grr(np.linspace(2, 8, SAMPLES)[1:], 1, 0, float(value)) > 0))
        return curves

    def ratio(value):
        """K over Schwarzschild's K of the same mass, less one, as a function of r."""
        return lambda r, *_: K(r, 1, 0, value) / K_s(r, 2) - 1

    def crossing(value, level, hi=8):
        """The greatest r between 2 r_0 and `hi` at which K stands `level` times Schwarzschild's."""
        found = zeros(ck, f"brans_dicke_sphere at n = {value:.6f} K = {level} K_s",
                      lambda r, *_: K(r, 1, 0, value) / K_s(r, 2) - level, 2, hi, [])
        return max(found) if found else None

    def edge(fn, scale, value):
        """Where a curve leaves the top of its panel, the r at which it reaches the top, found by
        bisection in log(r - 2 r_0); the curve falls at every r drawn, so it leaves once."""
        up, down = beyond(fn, x_axis, scale, [1, 0, value])
        falls = np.all(np.diff(fn(rs, 1, 0, value)) < 0)
        ck(f"brans_dicke_sphere at n = {value:.4f} falls across the range and leaves its panel once, at the top",
           falls and len(up) == 1 and up[0][0] == 2 and not down, f"{up} {down}")
        lo, hi = -40.0, math.log(6.0)
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if float(fn(2 + math.exp(mid), 1, 0, value)) > scale["to"] else (lo, mid)
        return 2 + math.exp(hi)

    rs = np.linspace(2, 8, SAMPLES)[1:]
    curves = reference() + own(BDS_N, "brans_dicke_sphere at n = 1/4")
    # The figures of the captions.
    centre_s = float(K_s.exact([2, 2]))
    ck("Schwarzschild's K at r = 2 r_0 is 3/4 r_0^4", centre_s == 0.75, f"{centre_s}")
    ck("brans_dicke_sphere the spheres shrink to a point at r = 2 r_0",
       float(g_thth.exact([2 + 2.0 ** -40, 1, 0, n])) < 1e-2 and float(g_thth.exact([2, 1, 0, n])) == 0)
    ck("brans_dicke_sphere K and R are positive on the range",
       np.all(K(rs, 1, 0, n) > 0) and np.all(R(rs, 1, 0, n) > 0))
    ck("Schwarzschild's K stays inside its panel", not any(beyond(K_s, x_axis, K_scale, [2])))
    for fn, power, label in ((K, 2 + 2 * n, "K"), (R, 1 + n, "R")):
        near = [2 + 2.0 ** -k for k in (20, 27, 34)]
        logs = [float(sp.log(fn.exact([r, 1, 0, n]))) for r in near]
        slopes = [(b - a) / (math.log(2.0 ** -k2) - math.log(2.0 ** -k1))
                  for (a, k1), (b, k2) in zip(zip(logs, (20, 27)), zip(logs[1:], (27, 34)))]
        ck(f"brans_dicke_sphere {label} diverges as (r - 2 r_0)^-{power}",
           all(abs(-s / power - 1) < 1e-5 for s in slopes), f"{slopes}")
    ck("brans_dicke_sphere the powers the caption writes are 2 + 2n and 1 + n", 2 + 2 * n == 5 / 2 and 1 + n == 5 / 4)
    crosses = zeros(ck, "brans_dicke_sphere K = K_s", ratio(n), 2, 8, [])
    ck("brans_dicke_sphere K crosses Schwarzschild's once on the range", len(crosses) == 1, f"{crosses}")
    cross = crosses[0]
    ck("brans_dicke_sphere K stands above Schwarzschild's inside the crossing and below it outside",
       np.all(ratio(n)(rs[rs < cross]) > 0) and np.all(ratio(n)(rs[rs > cross]) < 0))
    far = 1e4
    tail = (n * n - 2 * n + 2) / 2
    ck("brans_dicke_sphere K tends to (n^2 - 2n + 2)/2 of Schwarzschild's far out, 25/32 at n = 1/4",
       abs(float(K.exact([far, 1, 0, n])) / float(K_s.exact([far, 2])) / tail - 1) < 1e-3 and tail == 25 / 32)
    ck("brans_dicke_sphere R falls off as 2n(1 - n) r_0^2/r^4, 3 r_0^2/8 r^4 at n = 1/4",
       abs(float(R.exact([far, 1, 0, n])) * far ** 4 / (2 * n * (1 - n)) - 1) < 1e-3 and 2 * n * (1 - n) == 3 / 8)
    K_edge, R_edge = edge(K, K_scale, n), edge(R, R_scale, n)
    legend = [["line", "K", "$K$, the Kretschmann scalar"], ["line", "R", "$R$, the Ricci scalar"],
              ["line", "reference", "$K$ of Schwarzschild's black hole of the same mass"]]
    still = {
        "id": "omega", "label": "$n = 1/4$", "kind": "line", "x": x_axis, "panels": panels,
        "curves": curves, "marks": [], "legend": legend,
        "caption": [
            "The Kretschmann scalar $K$ and the Ricci scalar $R$ of Brans and Dicke's static sphere against $r$ "
            "($m = 0$, $n = 1/4$, which is $\\omega = 6$), the same at every $t$, $\\theta$, and $\\phi$, beside "
            "Schwarzschild's black hole of the same mass. At $r = 2r_0$ Schwarzschild's $K = 48r_0^2/r^6$ is a finite "
            f"${centre_s:g}/r_0^4$ on its horizon. Brans and Dicke's spheres shrink to zero area there, with no horizon "
            "to hide them, and $K$ grows without bound as $(r - 2r_0)^{-5/2}$.",
            "Schwarzschild's $R$ vanishes, and Brans and Dicke's is positive at every $r$ and grows without bound as "
            f"$(r - 2r_0)^{{-5/4}}$. $K$ crosses Schwarzschild's at $r = {cross:.3f}\\,r_0$ and far out falls to "
            "$25/32$ of it, while $R$ falls off as $3r_0^2/8r^4$.",
        ],
        "settings": "$m = 0$, $n = 1/4$, and $r_0 = 1$, the unit of every length, so that $GM/c^2 = r_0$, and "
                    "Schwarzschild's $r_s = 2r_0$.",
        "clip": ("We draw $K\\,r_0^4$ on a logarithmic scale from $10^{-4}$ to $10^{2}$ and $R\\,r_0^2$ on one from "
                 "$10^{-6}$ to $10^{2}$, and cut each curve where it leaves its panel, as $K$ does inside "
                 f"$r = {K_edge:.3f}\\,r_0$ and $R$ inside $r = {R_edge:.3f}\\,r_0$."),
    }

    # The movie in n, at fixed mass, from omega = 198 to the 6 of the first view.
    frames = []
    steps = int(BDS_N / BDS_STEP)
    for k in range(1, steps + 1):
        value = BDS_STEP * k
        frames.append({"label": f"$n = {float(value):.2f}$", "value": float(value),
                       "curves": own(value, f"brans_dicke_sphere at n = {float(value):.2f}"), "marks": []})
    first, last = float(BDS_STEP), n
    doubled = [crossing(f["value"], 2) for f in frames]
    ck("brans_dicke_sphere K stands above twice Schwarzschild's in a layer that thins as n falls",
       all(d is not None for d in doubled) and all(a < b for a, b in zip(doubled, doubled[1:])), f"{doubled}")
    for value in (first, last):
        r2 = crossing(value, 2)
        ck(f"brans_dicke_sphere at n = {value} K is above twice Schwarzschild's inside r = {r2:.4f} alone",
           np.all(K(rs[rs < r2], 1, 0, value) > 2 * K_s(rs[rs < r2], 2))
           and np.all(K(rs[rs > r2], 1, 0, value) < 2 * K_s(rs[rs > r2], 2)))
    values = [f["value"] for f in frames]
    ck("brans_dicke_sphere R shrinks toward zero at every r as n falls",
       all(np.all(R(rs, 1, 0, a) < R(rs, 1, 0, b)) for a, b in zip(values, values[1:])))
    outer = rs[rs >= 3]
    gaps = [float(np.max(np.abs(ratio(v)(outer)))) for v in values]
    ck("brans_dicke_sphere K approaches Schwarzschild's outside r = 3 r_0 as n falls",
       all(a < b for a, b in zip(gaps, gaps[1:])), f"{gaps[0]:.4f} to {gaps[-1]:.4f}")
    edges = [(edge(K, K_scale, f["value"]), edge(R, R_scale, f["value"])) for f in frames]
    reach = max(max(e) for e in edges) - 2
    stated = math.ceil(reach * 100) / 100
    ck("brans_dicke_sphere every curve of the movie leaves its panel nearer 2 r_0 as n falls, within the reach stated",
       all(a[0] < b[0] and a[1] < b[1] for a, b in zip(edges, edges[1:])) and reach <= stated, f"{reach:.4f}")
    # The least omega Cassini allows, and how near Schwarzschild's K the sphere's is there.
    n_cassini = 2 / (BDS_CASSINI + 2)
    ck("brans_dicke_sphere the least omega Cassini allows is n = 2/(omega + 2) at m = 0, below 5 x 10^-5",
       omega(Fraction(2, BDS_CASSINI + 2)) == BDS_CASSINI and n_cassini < 5e-5)
    fine = np.linspace(2, 8, 600001)[1:]
    layer = []
    for value in (n_cassini, n_cassini / 10, n_cassini / 100):
        apart = np.abs(K(fine, 1, 0, value) / K_s(fine, 2) - 1)
        layer.append(crossing(value, 1.01, 2.01))
        wide = np.logspace(math.log10(8), 8, 2001)
        ck(f"brans_dicke_sphere at n = {value:.3e} K is within a percent of Schwarzschild's outside r = {layer[-1]:.6f}",
           np.all(apart[fine > layer[-1]] <= 0.01) and np.all(apart[fine < layer[-1] - 1e-5] > 0.01)
           and np.all(np.abs(K(wide, 1, 0, value) / K_s(wide, 2) - 1) <= 0.01))
    r_percent = layer[0]
    ck("brans_dicke_sphere the layer of a percent thins further above Cassini's omega",
       layer[0] > layer[1] > layer[2], f"{layer}")
    movie_view = {
        "id": "falling", "label": "$1/100 \\le n \\le 1/4$", "kind": "line", "x": x_axis, "panels": panels,
        "curves": reference(), "marks": [], "legend": legend,
        "caption": [
            "The Kretschmann scalar $K$ and the Ricci scalar $R$ of Brans and Dicke's static sphere against $r$ as "
            "$n$ grows at fixed mass ($m = 0$, $1/100 \\le n \\le 1/4$, which is $198 \\ge \\omega \\ge 6$), beside "
            "Schwarzschild's black hole. As $n$ falls, $\\omega = 2(1 - n)/n$ grows, $R$ shrinks toward zero, and $K$ "
            "approaches Schwarzschild's except in a layer beside $r = 2r_0$ that thins with $n$ and where it still "
            f"grows without bound: $K$ is above twice Schwarzschild's inside $r = {crossing(last, 2):.3f}\\,r_0$ at "
            f"$n = 1/4$ and inside $r = {crossing(first, 2):.3f}\\,r_0$ at $n = 1/100$.",
            "In 1961 Brans and Dicke took $\\omega \\ge 6$ from Mercury's perihelion, and the last frame is that least "
            f"value. The Cassini spacecraft puts $\\omega$ above ${BDS_CASSINI // 1000}{{,}}000$, where "
            f"$n < {tex_number(n_cassini, 2)}$ and $K$ is within a percent of Schwarzschild's outside "
            f"$r = {r_percent:.4f}\\,r_0$.",
        ],
        "settings": "$m = 0$ and $r_0 = 1$, the unit of every length, so that $GM/c^2 = r_0$, and $n$ in steps of "
                    f"$1/{int(1 / BDS_STEP)}$.",
        "clip": ("The scales are those of $n = 1/4$, and we cut each curve where it leaves the top of its panel, "
                 f"within ${stated:g}\\,r_0$ of $r = 2r_0$ in every frame; no curve passes the bottom."),
        "movie": {"variable": "$n$", "seconds": 5, "frames": frames},
    }
    return [still, movie_view]


# ---------------------------------------------------------------- Bowers and Liang's star

BL_Q = Fraction(1, 4)             # Q, as every drawing of the star takes it
BL_RADIUS = Fraction(9, 8)        # R in units of r_s, as every drawing of the star takes it
BL_REACH = 3                      # the plot runs from the centre to r = 3 r_s
BL_STEP = Fraction(1, 36)         # the step of each movie, which lands on 8/9 and 1/4 and on no zero of R
BL_COMPACT = [BL_STEP * k for k in range(15, 36)]       # r_s/R of the frames of the first movie
BL_QS = [BL_STEP * k for k in range(1, 19)]             # Q of the frames of the second movie


def bowers_liang(ck, src):
    """Bowers and Liang's star at r_s = 1, its Kretschmann and Ricci scalars against r, inside the
    star and in Schwarzschild's vacuum outside it.

    The areal chart publishes K and R in r alone, through Z, N and P, and holds from the centre to
    the surface r = R alone, so its scalars are drawn from 0 to R and Schwarzschild's, read from
    schwarzschild.json, from R to the edge of the plot at 3 r_s, each a stretch of the same class.
    Two panels share the axis of r: K on a logarithmic axis from 10^-2 to 10^4 and R on a signed
    one, linear within 0.1/r_s^2 of zero, from -10^2 to 10. The radius R is a parameter of the
    chart, so the Ricci scalar is written R^mu_mu. The first view is the star at Q = 1/4 and
    R = 9 r_s/8, as every other drawing of it takes; the second steps r_s/R from 5/12 to 35/36 at
    Q = 1/4 and fixed mass, through 8/9, toward the limit 1 - 3^(-1/Q) = 80/81 of a finite
    central pressure; the third steps Q from 1/36 to 1/2 at R = 9 r_s/8, through 1/4, to
    Schwarzschild's star of equal pressures, whose central pressure is infinite at that radius.

    What the captions state is checked here: g^rr is positive on every stretch drawn, so no
    horizon is marked; the density and the two pressures are read from the published Einstein
    tensor, G^t_t = -8 pi G rho/c^2 and G^r_r and G^theta_theta 8 pi G/c^4 times p_r and p_perp,
    and the published R is minus its trace; at the centre K = 3 r_s^2(P^2 + 4)/R^6 and
    R = 3 r_s(2 - P)/R^3, and for the first view P = 1 + sqrt 3 and both pressures are
    rho c^2/sqrt 3 there; K is least inside the star at the radius and value stated and greatest
    at the surface, where Schwarzschild's is 12 r_s^2/R^6; R is negative throughout the star and
    -3 r_s/R^3 at its surface, where p_r = 0 and p_perp = rho c^2; along the first movie K at the
    centre grows, without bound toward 80/81, where N vanishes at the centre, and R is positive
    throughout the star below r_s/R = 4/5, where it vanishes at the surface, and negative
    throughout above 65/81, where it vanishes at the centre; along the second K at the centre
    grows with Q, the star at Q = 1/2 has equal pressures and is interior_schwarzschild.json's,
    N vanishes at its centre, and R is positive about the centre alone below
    Q = ln(3/2)/ln 9, negative throughout up to Q = 3/8 and positive about the surface alone
    above it; and the clip stated is where each curve leaves its panel."""
    x_axis = axis("$r/r_s$", 0, BL_REACH)
    K_scale = {"kind": "log", "from": 1e-2, "to": 1e4}
    R_scale = {"kind": "symlog", "from": -100, "to": 10, "linear": 0.1}
    panels = [{"id": "K", "name": "$K\\,r_s^4$", "aspect": 2, "scale": K_scale, "ticks": value_ticks(K_scale)},
              {"id": "R", "name": "$R^\\mu{}_\\mu\\,r_s^2$", "aspect": 2, "scale": R_scale,
               "ticks": value_ticks(R_scale)}]
    shape_checks(ck, "bowers_liang", [p["aspect"] for p in panels])
    height = {p["id"]: PLOT_W / p["aspect"] for p in panels}
    names = ["R", "r_s", "Q"]
    K = nr_published(src, "bowers_liang", "areal", "kretschmann", ["r"], names)
    Ric = nr_published(src, "bowers_liang", "areal", "ricci_scalar", ["r"], names)
    grr = nr_published(src, "bowers_liang", "areal", "inverse_metric_components", ["r"], names, ("r", "r"))
    G = {i: nr_published(src, "bowers_liang", "areal", "einstein_tensor", ["r"], names, (i, i))
         for i in ("t", "r", "\\theta", "\\phi")}
    K_s = nr_published(src, "schwarzschild", "spherical", "kretschmann", ["r"], ["r_s"])
    R_s = nr_published(src, "schwarzschild", "spherical", "ricci_scalar", ["r"], ["r_s"])
    grr_s = nr_published(src, "schwarzschild", "spherical", "inverse_metric_components", ["r"], ["r_s"], ("r", "r"))

    def inside(radius):
        return (0.0, min(float(radius), BL_REACH))

    def own(radius, q, name):
        """The star's own scalars from its centre to its surface, the stretch its chart holds."""
        params = {"R": float(radius), "r_s": 1.0, "Q": float(q)}
        rs = np.linspace(*inside(radius), SAMPLES)
        ck(f"{name} has no horizon inside the star", np.all(grr(rs, float(radius), 1, float(q)) > 0))
        return (line_curves(ck, name, "K", "K", K, "r", x_axis, K_scale, height["K"], params, inside(radius))
                + line_curves(ck, name, "R", "R", Ric, "r", x_axis, R_scale, height["R"], params, inside(radius)))

    def vacuum(radius, name):
        """Schwarzschild's scalars of the same mass from the surface to the edge of the plot."""
        outside = (float(radius), float(BL_REACH))
        rs = np.linspace(*outside, SAMPLES)
        ck(f"{name} has no horizon outside the star", np.all(grr_s(rs, 1) > 0))
        return (line_curves(ck, f"{name} outside", "K", "K", K_s, "r", x_axis, K_scale, height["K"], {"r_s": 1},
                            outside)
                + line_curves(ck, f"{name} outside", "R", "R", R_s, "r", x_axis, R_scale, height["R"], {"r_s": 1},
                              outside))

    def centre_P(radius, q):
        """P at the centre, where Z = 1, exactly."""
        return 2 / (3 * (1 - 1 / sp.Rational(radius)) ** sp.sympify(q) - 1)

    def cut(fn, params, within):
        """Where a curve leaves the top and the bottom of each panel, as beyond() finds it."""
        return {name: beyond(fn_, x_axis, scale, params, within)
                for name, fn_, scale in (("K", fn[0], K_scale), ("R", fn[1], R_scale))}

    def signs(radius, q):
        """How R changes sign through the star: how many times, and whether it is positive at the
        centre and at the surface."""
        values = Ric(np.linspace(0, float(radius), SAMPLES), float(radius), 1, float(q))
        return int(np.sum(np.sign(values[:-1]) * np.sign(values[1:]) < 0)), values[0] > 0, values[-1] > 0

    legend = [["line", "K", "$K$, the Kretschmann scalar"], ["line", "R", "$R^\\mu{}_\\mu$, the Ricci scalar"]]
    clip = ("We draw $K\\,r_s^4$ on a logarithmic scale from $10^{-2}$ to $10^{4}$ and $R^\\mu{}_\\mu\\,r_s^2$ on a "
            "signed logarithmic scale, linear within $0.1$ of zero, from $-10^{2}$ to $10$")

    # The star every other drawing takes.
    R0, q0 = float(BL_RADIUS), float(BL_Q)
    name = "bowers_liang at Q = 1/4, R = 9 r_s/8"
    curves = own(BL_RADIUS, BL_Q, name) + vacuum(BL_RADIUS, name)
    at = [R0, 1, q0]
    P0 = centre_P(BL_RADIUS, BL_Q)
    ck(f"{name} P at the centre is 1 + sqrt 3", abs(float(P0 - 1 - sp.sqrt(3))) < 1e-50, f"{float(P0)}")
    centre_K, centre_R = float(K.exact([0] + at)), float(Ric.exact([0] + at))
    ck(f"{name} K at the centre is 3 r_s^2(P^2 + 4)/R^6",
       abs(centre_K / float(3 * (P0 ** 2 + 4) / sp.Rational(BL_RADIUS) ** 6) - 1) < 1e-12, f"{centre_K}")
    ck(f"{name} R at the centre is 3 r_s(2 - P)/R^3",
       abs(centre_R / float(3 * (2 - P0) / sp.Rational(BL_RADIUS) ** 3) - 1) < 1e-12, f"{centre_R}")
    rho = -float(G["t"].exact([0] + at))       # 8 pi G rho/c^2, the same at every r
    rs = np.linspace(0, R0, SAMPLES)
    ck(f"{name} the density is uniform", np.allclose(-G["t"](rs, *at), rho, rtol=1e-13, atol=0), f"{rho}")
    trace = [sum(float(G[i].exact([r] + at)) for i in G) for r in (0, 0.3, 0.7, 1.0, R0)]
    published_R = [float(Ric.exact([r] + at)) for r in (0, 0.3, 0.7, 1.0, R0)]
    ck(f"{name} R is minus the trace of the Einstein tensor",
       all(abs(a + b) < 1e-12 * abs(b) for a, b in zip(published_R, trace)))
    ck(f"{name} the two pressures at the centre are rho c^2/sqrt 3",
       all(abs(float(G[i].exact([0] + at)) * math.sqrt(3) / rho - 1) < 1e-12 for i in ("r", "\\theta")))
    ck(f"{name} at the surface p_r = 0 and p_perp = rho c^2",
       abs(float(G["r"].exact([R0] + at))) < 1e-12 * rho and abs(float(G["\\theta"].exact([R0] + at)) / rho - 1) < 1e-12)
    fine = np.linspace(0, R0, 200001)
    Ks = K(fine, *at)
    least = int(np.argmin(Ks))
    r_least, K_least = float(fine[least]), float(Ks[least])
    ck(f"{name} K falls from the centre to its least and rises from there to the surface",
       0 < least < len(fine) - 1 and np.all(np.diff(Ks[:least + 1]) < 0) and np.all(np.diff(Ks[least:]) > 0),
       f"r = {r_least:.4f}, K = {K_least:.4f}")
    surface_K, outside_K = float(K.exact([R0] + at)), float(K_s.exact([R0, 1]))
    ck(f"{name} Schwarzschild's K at the surface is 12 r_s^2/R^6", abs(outside_K * R0 ** 6 / 12 - 1) < 1e-12)
    surface_R = float(Ric.exact([R0] + at))
    ck(f"{name} R is negative throughout the star", np.all(Ric(rs, *at) < 0))
    ck(f"{name} R at the surface is -3 r_s/R^3", abs(surface_R * R0 ** 3 / -3 - 1) < 1e-12, f"{surface_R}")
    ck(f"{name} R vanishes outside", np.all(R_s(np.linspace(R0, BL_REACH, 101), 1) == 0))
    for k, (up, down) in cut((K, Ric), at, inside(BL_RADIUS)).items():
        ck(f"{name} {k} stays inside its panel", not up and not down, f"{up} {down}")
    for k, (up, down) in cut((K_s, R_s), [1], (R0, float(BL_REACH))).items():
        ck(f"{name} {k} outside stays inside its panel", not up and not down, f"{up} {down}")
    still = {
        "id": "star", "label": "$Q = 1/4$, $R = 9r_s/8$", "kind": "line", "x": x_axis, "panels": panels,
        "curves": curves, "marks": [], "legend": legend,
        "caption": [
            "The Kretschmann scalar $K$ and the Ricci scalar $R^\\mu{}_\\mu$ of Bowers and Liang's star against $r$ "
            "($Q = 1/4$, $R = 9r_s/8$), the same at every $t$, $\\theta$, and $\\phi$, inside the star and in "
            f"Schwarzschild's vacuum outside it. Inside, $K$ falls from ${centre_K:.2f}/r_s^4$ at the centre to "
            f"${K_least:.2f}/r_s^4$ at $r = {r_least:.3f}\\,r_s$ and rises to ${surface_K:.2f}/r_s^4$ at the surface. "
            f"Outside it is Schwarzschild's $12r_s^2/r^6$, which is ${outside_K:.3f}/r_s^4$ at $r = R$.",
            "The Ricci scalar is $R^\\mu{}_\\mu = 8\\pi G(\\rho c^2 - p_r - 2p_\\perp)/c^4$, negative throughout the "
            f"star. At the centre both pressures are $\\rho c^2/\\sqrt{{3}}$ and $R^\\mu{{}}_\\mu = {centre_R:.3f}/r_s^2$, "
            "and at the surface $p_r = 0$, $p_\\perp = \\rho c^2$, and "
            f"$R^\\mu{{}}_\\mu = -3r_s/R^3 = {surface_R:.3f}/r_s^2$. Outside it vanishes. Both scalars jump at the "
            "surface, where the density and the pressure across the radius fall to zero.",
        ],
        "settings": "$r_s = 1$, the unit of every length, $Q = 1/4$, and $R = 9r_s/8$.",
        "clip": clip + ", and every value drawn lies within them.",
    }

    # The star of fixed mass made more compact, toward the limit of a finite central pressure.
    frames, centres = [], []
    for value in BL_COMPACT:
        radius = 1 / value
        name = f"bowers_liang at r_s/R = {float(value):.4f}"
        curves = own(radius, BL_Q, name) + vacuum(radius, name)
        params = [float(radius), 1, q0]
        for k, (up, down) in cut((K, Ric), params, inside(radius)).items():
            ck(f"{name} {k} stays inside its panel", not up and not down, f"{up} {down}")
        centres.append(float(K.exact([0] + params)))
        ck(f"{name} K at the centre is 3 r_s^2(P^2 + 4)/R^6", abs(centres[-1] / float(
            3 * (centre_P(radius, BL_Q) ** 2 + 4) / sp.Rational(radius) ** 6) - 1) < 1e-12)
        frames.append({"label": f"$r_s/R = {float(value):.3f}$", "value": float(value), "curves": curves, "marks": []})
    ck("bowers_liang K at the centre grows with r_s/R", all(b > a for a, b in zip(centres, centres[1:])),
       f"{centres[0]:.4f} to {centres[-1]:.4f}")
    limit = 1 - sp.Integer(3) ** (-1 / sp.Rational(BL_Q))
    ck("bowers_liang the limit of a finite central pressure is r_s/R = 80/81 at Q = 1/4", limit == sp.Rational(80, 81))
    ck("bowers_liang N vanishes at the centre at r_s/R = 80/81",
       sp.simplify(3 * (1 - limit) ** sp.Rational(BL_Q) - 1) == 0)
    near = [float(K.exact([0, 1 / (float(limit) - d), 1, q0])) * d ** 2 for d in (1e-4, 1e-6, 1e-8)]
    ck("bowers_liang K at the centre grows as (80/81 - r_s/R)^-2 toward the limit",
       abs(near[2] / near[1] - 1) < 1e-3 and abs(near[1] / near[0] - 1) < 1e-1, f"{near}")
    # Where R changes sign: at the surface at r_s/R = 4/5, at the centre at 65/81, where P = 2, and
    # last of all, at c_last, where its greatest value inside the star, in a shell, falls to zero.
    surface_zero, centre_zero = Fraction(4, 5), 1 - Fraction(2, 3) ** 4
    ck("bowers_liang R vanishes at the surface at r_s/R = 4/5",
       abs(float(Ric.exact([1 / float(surface_zero), 1 / float(surface_zero), 1, q0]))) < 1e-12)
    ck("bowers_liang R vanishes at the centre at r_s/R = 65/81, where P = 2 and 3 p_r = rho c^2",
       centre_zero == Fraction(65, 81) and centre_P(1 / centre_zero, BL_Q) == 2
       and abs(float(Ric.exact([0, 1 / float(centre_zero), 1, q0]))) < 1e-12)

    def greatest(c):
        """The greatest value of R inside the star at r_s/R = c and where it lies, by the golden
        section about the greatest of a fine grid."""
        radius = 1 / c
        grid = np.linspace(0, radius, SAMPLES)
        k = int(np.argmax(Ric(grid, radius, 1, q0)))
        lo, hi = grid[max(k - 1, 0)], grid[min(k + 1, SAMPLES - 1)]
        golden = (math.sqrt(5) - 1) / 2
        for _ in range(80):
            a, b = hi - golden * (hi - lo), lo + golden * (hi - lo)
            lo, hi = (a, hi) if float(Ric(a, radius, 1, q0)) < float(Ric(b, radius, 1, q0)) else (lo, b)
        return float(Ric(0.5 * (lo + hi), radius, 1, q0)), 0.5 * (lo + hi) / radius
    lo, hi = float(centre_zero), 0.9
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if greatest(mid)[0] > 0 else (lo, mid)
    c_last = hi
    ck("bowers_liang R turns negative at the surface, then at the centre, and last in a shell between",
       float(surface_zero) < float(centre_zero) < c_last and 0.05 < greatest(c_last)[1] < 0.95,
       f"c_last = {c_last:.6f}, the shell at r = {greatest(c_last)[1]:.3f} R")
    scan = [(c, *signs(1 / c, BL_Q)) for c in np.linspace(float(BL_COMPACT[0]), float(BL_COMPACT[-1]), 1201)]
    ck("bowers_liang R is positive throughout the star below r_s/R = 4/5, negative throughout above c_last, and "
       "of both signs between",
       all((n == 0 and centre and surface) if c < float(surface_zero) - 1e-9
           else (n == 0 and not centre and not surface) if c > c_last + 1e-9
           else n > 0 if float(surface_zero) + 1e-9 < c < c_last - 1e-9 else True
           for c, n, centre, surface in scan)
       and any(float(surface_zero) < c < c_last for c, *_ in scan), f"c_last = {c_last:.6f}")
    first, last = BL_COMPACT[0], BL_COMPACT[-1]
    compact = {
        "id": "compactness", "label": f"${first} \\le r_s/R \\le {last}$", "kind": "line", "x": x_axis, "panels": panels,
        "curves": [], "marks": [], "legend": legend,
        "caption": [
            "The Kretschmann scalar $K$ and the Ricci scalar $R^\\mu{}_\\mu$ of Bowers and Liang's star against $r$ "
            f"as a star of fixed mass is made more compact ($Q = 1/4$, ${first} \\le r_s/R \\le {last}$). Outside the star "
            "both are Schwarzschild's in every frame. At the centre $K = 3r_s^2(P^2 + 4)/R^6$, which grows from "
            f"${centres[0]:.4f}/r_s^4$ at $r_s/R = {first}$ to ${centres[-1]:.1f}/r_s^4$ at $r_s/R = {last}$, and without "
            "bound as $r_s/R$ nears $1 - 3^{-1/Q} = 80/81$, where the central pressure is infinite.",
            "The Ricci scalar is $R^\\mu{}_\\mu = 8\\pi G(\\rho c^2 - p_r - 2p_\\perp)/c^4$. It is positive throughout "
            f"the star for $r_s/R < 4/5$, negative throughout for $r_s/R > {c_last:.3f}$, and of both signs between. "
            "It turns negative first at the surface, at $r_s/R = 4/5$, then at the centre, at $r_s/R = 65/81$, where "
            "$3p_r = \\rho c^2$, and last in a shell between the two.",
        ],
        "settings": f"$r_s = 1$, the unit of every length, $Q = 1/4$, and $r_s/R$ in steps of ${BL_STEP}$.",
        "clip": "The scales are those of $Q = 1/4$, $R = 9r_s/8$, and every value drawn lies within them.",
        "movie": {"variable": "$r_s/R$", "seconds": 5, "frames": frames},
    }

    # The star of fixed mass and radius as its anisotropy changes, to Schwarzschild's star.
    frames, centres, leaving = [], [], {}
    for value in BL_QS:
        name = f"bowers_liang at Q = {float(value):.4f}"
        q = float(value)
        params = [R0, 1, q]
        curves = own(BL_RADIUS, value, name)
        for k, (up, down) in cut((K, Ric), params, inside(BL_RADIUS)).items():
            if up or down:
                leaving[(float(value), k)] = (up, down)
        if value < Fraction(1, 2):
            centres.append(float(K.exact([0] + params)))
        frames.append({"label": f"$Q = {q:.3f}$", "value": q, "curves": curves, "marks": []})
    ck("bowers_liang K at the centre grows with Q below Q = 1/2", all(b > a for a, b in zip(centres, centres[1:])),
       f"{centres[0]:.4f} to {centres[-1]:.4f}")
    half = sp.Rational(1, 2)
    ck("bowers_liang N vanishes at the centre at Q = 1/2 and R = 9 r_s/8",
       sp.simplify(3 * (1 - 1 / sp.Rational(BL_RADIUS)) ** half - 1) == 0)
    rs = np.linspace(0.05, R0, 41)
    ck("bowers_liang at Q = 1/2 has equal pressures",
       np.allclose(G["r"](rs, R0, 1, 0.5), G["\\theta"](rs, R0, 1, 0.5), rtol=1e-12, atol=0))
    K_is = nr_published(src, "interior_schwarzschild", "spherical", "kretschmann", ["r"], ["r_s", "R"])
    R_is = nr_published(src, "interior_schwarzschild", "spherical", "ricci_scalar", ["r"], ["r_s", "R"])
    ck("bowers_liang at Q = 1/2 is Schwarzschild's star of uniform density",
       all(abs(float(K.exact([r, R0, 1, 0.5])) / float(K_is.exact([r, 1, R0])) - 1) < 1e-12
           and abs(float(Ric.exact([r, R0, 1, 0.5])) - float(R_is.exact([r, 1, R0]))) < 1e-12 * abs(float(R_is.exact([r, 1, R0])))
           for r in (0.05, 0.4, 0.8, R0)))
    q_centre = sp.log(sp.Rational(3, 2)) / sp.log(9)
    q_surface = Fraction(3, 8)
    ck("bowers_liang R vanishes at the centre at Q = ln(3/2)/ln 9, where P = 2",
       abs(float((centre_P(BL_RADIUS, q_centre) - 2).evalf(60))) < 1e-50
       and abs(float(Ric.exact([0, R0, 1, float(q_centre)]))) < 1e-12, f"Q = {float(q_centre):.6f}")
    ck("bowers_liang R vanishes at the surface at Q = 3/8, where p_perp = rho c^2/2",
       abs(float(Ric.exact([R0, R0, 1, float(q_surface)]))) < 1e-12
       and abs(float(G["\\theta"].exact([R0, R0, 1, float(q_surface)])) / rho - 0.5) < 1e-12)
    ck("bowers_liang R is positive about the centre alone below Q = ln(3/2)/ln 9, negative throughout up to "
       "Q = 3/8, and positive about the surface alone above it",
       all((n == 1 and centre and not surface) if q < float(q_centre) - 1e-9
           else (n == 0 and not centre and not surface) if float(q_centre) + 1e-9 < q < float(q_surface) - 1e-9
           else (n == 1 and surface and not centre) if q > float(q_surface) + 1e-9 else True
           for q, n, centre, surface in [(q, *signs(BL_RADIUS, q)) for q in np.linspace(1e-3, 0.5, 1000)]))
    ck("bowers_liang leaves its panels at Q = 1/2 alone, at the centre",
       set(leaving) == {(0.5, "K"), (0.5, "R")} and leaving[(0.5, "K")][0] and not leaving[(0.5, "K")][1]
       and leaving[(0.5, "R")][1] and not leaving[(0.5, "R")][0]
       and leaving[(0.5, "K")][0][0][0] == 0 and leaving[(0.5, "R")][1][0][0] == 0, f"{leaving}")
    K_top, R_bottom = leaving[(0.5, "K")][0][0][1], leaving[(0.5, "R")][1][0][1]
    anisotropy = {
        "id": "anisotropy", "label": f"${BL_QS[0]} \\le Q \\le {BL_QS[-1]}$", "kind": "line", "x": x_axis, "panels": panels,
        "curves": vacuum(BL_RADIUS, "bowers_liang at R = 9 r_s/8"), "marks": [], "legend": legend,
        "caption": [
            "The Kretschmann scalar $K$ and the Ricci scalar $R^\\mu{}_\\mu$ of Bowers and Liang's star against $r$ "
            f"as its pressures change at fixed mass and radius (${BL_QS[0]} \\le Q \\le {BL_QS[-1]}$, $R = 9r_s/8$). Outside the "
            "star both are Schwarzschild's in every frame. At $Q = 1/2$ the two pressures are equal and the star "
            "is Schwarzschild's star of uniform density, whose central pressure is infinite at this radius, and so "
            "is $K$ at its centre. As $Q$ falls, $K$ at the centre falls, to "
            f"${centres[0]:.3f}/r_s^4$ at $Q = {BL_QS[0]}$.",
            "The Ricci scalar is $R^\\mu{}_\\mu = 8\\pi G(\\rho c^2 - p_r - 2p_\\perp)/c^4$. It is positive about the "
            f"centre for $Q < \\ln(3/2)/\\ln 9 = {float(q_centre):.3f}$, where $3p_r < \\rho c^2$ at the centre, "
            "negative throughout the star from there to $Q = 3/8$, and positive about the surface for $Q > 3/8$, "
            "where $2p_\\perp < \\rho c^2$ at the surface.",
        ],
        "settings": f"$r_s = 1$, the unit of every length, $R = 9r_s/8$, and $Q$ in steps of ${BL_STEP}$.",
        "clip": ("The scales are those of $Q = 1/4$, $R = 9r_s/8$, and we cut each curve where it leaves its panel, "
                 f"as $K$ passes the top inside $r = {K_top:.3f}\\,r_s$ and $R^\\mu{{}}_\\mu$ the bottom inside "
                 f"$r = {R_bottom:.3f}\\,r_s$ at $Q = 1/2$."),
        "movie": {"variable": "$Q$", "seconds": 5, "frames": frames},
    }
    return [still, compact, anisotropy]


# ---------------------------------------------------------------- the Curzon-Chazy particle

CC_CELLS = (57, 114)              # the cells of the map across and up, each a nineteenth of m on a side
CC_RANGE = (1e-3, 1e3)            # the colours of K m^4 run from the first to the second, as the Bach-Weyl ring's
CC_LINE = (1e-6, 1e6)             # the axis of K m^4 on both panels of the line view


def crossing(fn, lo, hi):
    """The point between lo and hi where fn changes sign, by bisection to the last digit."""
    f_lo = fn(lo)
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if mid in (lo, hi):
            break
        lo, hi, f_lo = (mid, hi, fn(mid)) if np.sign(fn(mid)) == np.sign(f_lo) else (lo, mid, f_lo)
    return 0.5 * (lo + hi)


def curzon_chazy(ck, src):
    """The Curzon-Chazy particle at m = 1, its Kretschmann scalar over the half plane of Weyl's
    rho and z and against r on the axis and in the plane theta = pi/2.

    Weyl's chart publishes K in rho, z and R = sqrt(rho^2 + z^2), which the reader writes out, so
    K is a function of rho and z alone, the same at every t and phi; the particle is a vacuum and
    publishes R = 0 in both charts, so the figure is K alone. The first view is a map over the box
    0 <= rho <= 3m, |z| <= 3m, the window of the research list, coloured on the Bach-Weyl ring's
    logarithmic scale from 10^-3 to 10^3 in units of 1/m^4, with the singularity marked at the
    origin and the two points of the axis where K vanishes, z = +-m. Its cells are a nineteenth of
    m on a side, coarse enough that K at the centre of the cell beside the singularity, near
    10^282/m^4, is a double, as it is not at a twentieth, where it passes e^709. The second view
    is the spherical chart's K against r, from 0 to 3m, held at theta = 0 on one panel and at
    theta = pi/2 on the other, each on a logarithmic axis from 10^-6 to 10^6; g^rr is positive at
    every point, so nothing is marked.

    What the captions state is checked here, in sixty digits and in eighty where K passes the
    doubles: R = 0 in both charts; K is 16m^2(3R^4(R - m)^2 + m^2 rho^2(3R^2 - 3mR + m^2))
    e^(2m^2 rho^2/R^4 - 4m/R)/R^12, whose last factor has no real root, so K vanishes at
    rho = 0, z = +-m alone, and every published component of the Riemann tensor vanishes there
    with it; the spherical chart's K is Weyl's at the same point; on the axis K is
    48m^2(r - m)^2 e^(-4m/r)/r^8, with its two maxima at r = (1 +- 1/sqrt 3)m and the values
    stated, and r ln K tends to -4m toward r = 0; in the plane K falls from the singularity
    outward and r^2 ln K tends to 2m^2; along theta = 0.1, pi/4 and 3pi/4, r^2 ln K tends to
    2m^2 sin^2 theta; far out K r^6 tends to 48m^2 in every direction; the circumference
    2 pi sqrt(g_phiphi) of the circle about the axis grows without bound toward the origin along
    every direction drawn but the axis; and the clips stated are where K leaves the colours, on a
    grid four times as fine as the map's, and where each curve leaves its panel."""
    m = 1.0
    params = {"m": m}
    K = nr_published(src, "curzon_chazy", "weyl", "kretschmann", ["rho", "z"], ["m"])
    Ks = nr_published(src, "curzon_chazy", "spherical", "kretschmann", ["r", "theta"], ["m"])
    for system, coords in (("weyl", ["rho", "z"]), ("spherical", ["r", "theta"])):
        ricci = nr_published(src, "curzon_chazy", system, "ricci_scalar", coords, ["m"])
        ck(f"curzon_chazy {system} publishes R = 0", ricci.expr == 0, f"{ricci.expr}")

    # The scalar in closed form, and its zeros.
    def closed(rho, z):
        rho, z, mm = (sp.Float(float(v), 70) for v in (rho, z, m))
        R = sp.sqrt(rho ** 2 + z ** 2)
        bracket = 3 * R ** 4 * (R - mm) ** 2 + mm ** 2 * rho ** 2 * (3 * R ** 2 - 3 * mm * R + mm ** 2)
        return (16 * mm ** 2 * bracket * sp.exp(2 * mm ** 2 * rho ** 2 / R ** 4 - 4 * mm / R) / R ** 12).evalf(60)
    rng = np.random.default_rng(1924)
    points = [(float(r * math.sin(t)), float(r * math.cos(t)))
              for r, t in zip(rng.uniform(0.05, 4, 20), rng.uniform(0, math.pi, 20))]
    worst = max(abs(float(K.exact([rho, z, m]) / closed(rho, z)) - 1) for rho, z in points)
    ck("curzon_chazy K is 16m^2(3R^4(R - m)^2 + m^2 rho^2(3R^2 - 3mR + m^2))e^(2m^2 rho^2/R^4 - 4m/R)/R^12",
       worst < 1e-50, f"worst {worst:.1e} at {len(points)} points")
    R_, m_ = sp.symbols("R m", positive=True)
    ck("curzon_chazy 3R^2 - 3mR + m^2 has no real root, so K vanishes on the axis alone",
       sp.discriminant(3 * R_ ** 2 - 3 * m_ * R_ + m_ ** 2, R_) == -3 * m_ ** 2)
    worst = max(abs(float(Ks.exact([math.hypot(rho, z), math.atan2(rho, z), m]) / K.exact([rho, z, m])) - 1)
                for rho, z in points)
    ck("curzon_chazy the spherical chart's K is Weyl's at the same point", worst < 1e-12, f"worst {worst:.1e}")
    axis_rho = 1e-12
    entry = json.loads((build.METRICS_DIR / "curzon_chazy.json").read_text(encoding="utf-8"))
    weyl = next(s for s in entry["coordinates"] if s["id"] == "weyl")
    worst = 0.0
    for comp in weyl["riemann"]["variants"]["ulll"]["nonzero"]:
        riemann = nr_published(src, "curzon_chazy", "weyl", "riemann", ["rho", "z"], ["m"], comp["indices"])
        worst = max(worst, abs(float(riemann.exact([axis_rho, m, m]))), abs(float(riemann.exact([axis_rho, -m, m]))))
    ck("curzon_chazy every published Riemann component vanishes at z = +-m on the axis", worst < 1e-9, f"{worst:.1e}")
    ck("curzon_chazy K vanishes at z = +-m on the axis",
       all(abs(float(K.exact([axis_rho, z, m]))) < 1e-20 for z in (m, -m)))

    # The singularity, approached along the axis, in the plane and along three other directions.
    def logs(theta, lengths):
        return [float(sp.log(Ks.exact([r, theta, m], 80))) for r in lengths]
    tiny = (1e-6, 1e-7, 1e-8)
    along = [r * v for r, v in zip(tiny, logs(0.0, tiny))]
    ck("curzon_chazy K falls to zero as e^(-4m/r) toward the singularity along the axis",
       all(abs(v / (-4 * m) - 1) < 1e-3 for v in along), f"r ln K {along}")
    for theta in (0.1, math.pi / 4, math.pi / 2, 3 * math.pi / 4):
        rates = [r * r * v / (2 * m ** 2 * math.sin(theta) ** 2) for r, v in zip(tiny, logs(theta, tiny))]
        ck(f"curzon_chazy K grows as e^(2m^2 sin^2 theta/r^2) toward the singularity at theta = {theta:.4f}",
           all(abs(q - 1) < 1e-3 for q in rates), f"{rates}")
    for theta in (0.0, math.pi / 4, math.pi / 2):
        ratios = [float(Ks.exact([r, theta, m])) * r ** 6 / (48 * m ** 2) for r in (1e6, 1e7)]
        ck(f"curzon_chazy K falls off as 48m^2/r^6 far out at theta = {theta:.4f}",
           all(abs(q - 1) < 1e-5 for q in ratios), f"{ratios}")
    gpp = nr_published(src, "curzon_chazy", "weyl", "metric_components", ["rho", "z"], ["m"], ("\\phi", "\\phi"))
    growth = [[math.sqrt(float(gpp.exact([R * math.sin(theta), R * math.cos(theta), m]))) for R in (1e-1, 1e-2, 1e-3)]
              for theta in (0.1, math.pi / 4, math.pi / 2, 3 * math.pi / 4)]
    ck("curzon_chazy the circles about the axis grow without bound toward the singularity off the axis",
       all(c[0] < c[1] < c[2] and c[2] > 10 * m for c in growth), f"{growth}")
    ck("curzon_chazy the circumference of the circle through (rho, z) is 2 pi rho e^(m/R)",
       all(abs(math.sqrt(float(gpp.exact([rho, z, m]))) / (rho * math.exp(m / math.hypot(rho, z))) - 1) < 1e-12
           for rho, z in points))

    # The map.
    x_axis, y_axis = axis("$\\rho/m$", 0, 3), axis("$z/m$", -3, 3)
    aspect = (x_axis["to"] - x_axis["from"]) / (y_axis["to"] - y_axis["from"])
    shape_checks(ck, "curzon_chazy map", [aspect])
    scale = {"kind": "log", "from": CC_RANGE[0], "to": CC_RANGE[1]}
    nx, ny = CC_CELLS
    xs, ys, values = cells(K, x_axis, y_axis, nx, ny, [m])
    ck("curzon_chazy K is finite and positive at every cell", np.all(np.isfinite(values)) and np.all(values > 0),
       f"from {np.min(values):.3e} to {np.max(values):.3e}")
    _, _, finer = cells(K, x_axis, y_axis, 60, 120, [m])
    ck("curzon_chazy K at the cell beside the singularity passes the doubles at cells a twentieth of m on a side",
       not np.all(np.isfinite(finer)))
    rng = np.random.default_rng(1925)
    picks = {(int(i), int(j)) for i, j in zip(rng.integers(0, nx, 40), rng.integers(0, ny, 40))}
    picks |= {(i, j) for i in (0, 1) for j in (ny // 2 - 2, ny // 2 - 1, ny // 2, ny // 2 + 1)}
    picks = sorted(picks)
    agrees(ck, "curzon_chazy map", K, [[xs[i], ys[j], m] for i, j in picks], [values[j, i] for i, j in picks])
    fx, fy, fine = cells(K, x_axis, y_axis, 4 * nx, 4 * ny, [m])
    FX, FY = np.meshgrid(fx, fy)
    high, low = ~(fine <= CC_RANGE[1]), fine < CC_RANGE[0]
    from_origin = np.hypot(FX, FY)
    from_zero = np.minimum(np.hypot(FX, FY - m), np.hypot(FX, FY + m))
    near_zero = low & (from_zero < from_origin)
    near_origin = low & ~near_zero
    d_high, d_zero, d_origin = (float(from_origin[high].max()), float(from_zero[near_zero].max()),
                                float(from_origin[near_origin].max()))
    ck("curzon_chazy K passes the top of its colours only beside the singularity", d_high < 0.7, f"within {d_high:.4f} m")
    ck("curzon_chazy K passes the bottom of its colours only beside its zeros and the singularity",
       d_zero < 0.2 and d_origin < 0.4, f"within {d_zero:.4f} m and {d_origin:.4f} m")
    ck("curzon_chazy K falls below its colours beside the singularity only about the axis",
       float(np.max(np.arctan2(FX[near_origin], np.abs(FY[near_origin])))) < math.pi / 4)
    rows = [[rounded(v) for v in row] for row in values]
    meridian = {
        "id": "meridian", "label": "the half plane of $\\rho$ and $z$", "kind": "map",
        "x": x_axis, "y": y_axis,
        "colour": {"name": "$K\\,m^4$", "scale": scale, "ticks": value_ticks(scale)},
        "map": {"metric": "curzon_chazy", "system": "weyl", "scalar": "kretschmann", "x": "\\rho", "y": "z",
                "params": params, "nx": nx, "ny": ny, "cells": rows},
        "marks": [{"class": "singular", "at": [0.0, 0.0]}, {"class": "zero", "at": [0.0, m]},
                  {"class": "zero", "at": [0.0, -m]}],
        "legend": [["point", "singular", "the singularity, $\\rho = z = 0$"],
                   ["point", "zero", "where $K = 0$"]],
        "caption": [
            "The Kretschmann scalar $K$ of the Curzon-Chazy particle over the half plane of Weyl's $\\rho$ and $z$ "
            "($m = 1$), the same at every $t$ and $\\phi$. The particle is a vacuum, so its Ricci scalar vanishes, "
            "and $K$ is finite everywhere off the singularity at $\\rho = z = 0$.",
            "At the singularity $K$ depends on the direction it is approached from, with $R = \\sqrt{\\rho^2 + z^2}$ "
            "and $\\theta$ the angle from the axis. Along the axis it falls to zero as $e^{-4m/R}$, and along every "
            "other direction it grows as $e^{2m^2\\sin^2\\theta/R^2}$, faster than any power of $R$.",
            "Everywhere $K$ is $16m^2e^{2m^2\\rho^2/R^4 - 4m/R}/R^{12}$ times the sum of $3R^4(R - m)^2$ and "
            "$m^2\\rho^2(3R^2 - 3mR + m^2)$, and $3R^2 - 3mR + m^2$ is positive at every $R$, so $K$ vanishes only on "
            "the axis at $z = \\pm m$, and there every component of the Riemann tensor vanishes with it.",
            "Weyl's $\\rho$ and $z$ are the cylindrical coordinates of the flat space in which $c^2\\psi$ is Newton's "
            "potential of a point mass, and the metric stretches them: the circle about the axis through "
            "$(\\rho, z)$ has circumference $2\\pi\\rho\\,e^{m/R}$, which grows without bound as $R \\to 0$ at every "
            "$\\theta$ off the axis, so the singularity, drawn at one point, is a ring of infinite circumference.",
        ],
        "settings": "$m = GM/c^2 = 1$, the unit of every length.",
        "clip": ("We colour $K\\,m^4$ on a logarithmic scale from $10^{-3}$ to $10^{3}$ and give a value beyond it "
                 f"the colour of the end it passes: $K$ passes $10^{{3}}/m^4$ within ${d_high:.2f}\\,m$ of the "
                 f"singularity and falls below $10^{{-3}}/m^4$ within ${d_zero:.2f}\\,m$ of the two points where it "
                 f"vanishes and, about the axis, within ${d_origin:.2f}\\,m$ of the singularity."),
    }

    # K against r on the axis and in the plane, from the spherical chart.
    x_r = axis("$r/m$", 0, 3)
    line_scale = {"kind": "log", "from": CC_LINE[0], "to": CC_LINE[1]}
    panels = [{"id": "axis", "name": "$K(r, 0)\\,m^4$", "aspect": 2, "scale": line_scale,
               "ticks": value_ticks(line_scale)},
              {"id": "plane", "name": "$K(r, \\pi/2)\\,m^4$", "aspect": 2, "scale": line_scale,
               "ticks": value_ticks(line_scale)}]
    shape_checks(ck, "curzon_chazy line", [p["aspect"] for p in panels])
    height = PLOT_W / 2
    on_axis, in_plane = Held(Ks, {"\\theta": 0.0}), Held(Ks, {"\\theta": math.pi / 2})
    curves = (line_curves(ck, "curzon_chazy axis", "axis", "K", on_axis, "r", x_r, line_scale, height, params)
              + line_curves(ck, "curzon_chazy plane", "plane", "K", in_plane, "r", x_r, line_scale, height, params))
    grr = nr_published(src, "curzon_chazy", "spherical", "inverse_metric_components", ["r", "theta"], ["m"], ("r", "r"))
    # In sixty digits, since e^(-2m/r) on the axis is below the doubles inside r = 0.003 m.
    ck("curzon_chazy g^rr is positive at every r drawn on the axis and in the plane, so no horizon is marked",
       all(grr.exact([r, theta, m]) > 0 for r in np.linspace(0, 3, 301)[1:] for theta in (0.0, math.pi / 2)))
    worst = max(abs(float(on_axis.exact([r, m])) / (48 * m ** 2 * (r - m) ** 2 * math.exp(-4 * m / r) / r ** 8) - 1)
                for r in (0.2, 0.5, 0.9, 1.1, 2.0, 3.0))
    ck("curzon_chazy K on the axis is 48m^2(r - m)^2 e^(-4m/r)/r^8", worst < 1e-12, f"worst {worst:.1e}")
    peaks = [m * (1 - 1 / math.sqrt(3)), m * (1 + 1 / math.sqrt(3))]
    fine_r = np.linspace(1e-3, 3, 600001)
    K_axis = on_axis(fine_r, m)
    inner, outer = fine_r < m, fine_r > m
    found = [float(fine_r[inner][np.argmax(K_axis[inner])]), float(fine_r[outer][np.argmax(K_axis[outer])])]
    peak_K = [float(on_axis.exact([r, m])) for r in peaks]
    # It rises and falls between them, checked where it is a double well above rounding: from
    # r = 0.05 m, where it is 10^-26, and outside 10^-3 m of its zero.
    def rising(a, b):
        return np.all(np.diff(K_axis[(fine_r > a) & (fine_r < b)]) > 0)

    def falling(a, b):
        return np.all(np.diff(K_axis[(fine_r > a) & (fine_r < b)]) < 0)
    ck("curzon_chazy K on the axis is largest on either side of r = m at r = (1 +- 1/sqrt 3)m",
       all(abs(a - b) < 1e-4 for a, b in zip(found, peaks)) and rising(0.05, peaks[0]) and falling(peaks[0], m - 1e-3)
       and rising(m + 1e-3, peaks[1]) and falling(peaks[1], 3), f"r = {found}, K = {peak_K}")
    K_plane = in_plane(fine_r, m)
    finite = np.isfinite(K_plane)
    ck("curzon_chazy K in the plane falls from the singularity outward", np.all(np.diff(K_plane[finite]) < 0))
    # Where each curve leaves its panel, at the crossings of 10^-6 and 10^6 found by bisection.
    def level(fn, value):
        return lambda r: math.log10(float(fn(r, m))) - math.log10(value)
    up_a, down_a = beyond(on_axis, x_r, line_scale, [m])
    up_p, down_p = beyond(in_plane, x_r, line_scale, [m])
    # The formula is not finite at r = 0 itself, and K is below the doubles inside r = 0.006 m.
    ck("curzon_chazy K on the axis leaves its panel by the bottom alone, about r = 0 and about r = m",
       all(a == b == 0 for a, b in up_a) and len(down_a) == 2 and down_a[0][0] < 0.01
       and down_a[1][0] < m < down_a[1][1], f"{up_a} {down_a}")
    ck("curzon_chazy K in the plane leaves its panel by the top alone, about r = 0",
       not down_p and len(up_p) == 1 and up_p[0][0] == 0, f"{up_p}")
    cut_axis = crossing(level(on_axis, CC_LINE[0]), 0.05, 0.3)
    dip = max(m - crossing(level(on_axis, CC_LINE[0]), 0.9, m - 1e-6),
              crossing(level(on_axis, CC_LINE[0]), m + 1e-6, 1.1) - m)
    cut_plane = crossing(level(in_plane, CC_LINE[1]), 0.2, 0.6)
    ck("curzon_chazy the cuts stated bound the samples beyond each panel",
       down_a[0][1] < cut_axis < down_a[0][1] + 3 / (SAMPLES - 1) and up_p[0][1] < cut_plane < up_p[0][1] + 3 / (SAMPLES - 1)
       and m - dip <= down_a[1][0] and down_a[1][1] <= m + dip,
       f"axis {cut_axis:.4f}, about r = m by {dip:.2e}, plane {cut_plane:.4f}")
    line = {
        "id": "directions", "label": "the axis and the plane $\\theta = \\pi/2$", "kind": "line", "x": x_r,
        "panels": panels, "curves": curves, "marks": [],
        "legend": [["line", "K", "$K$, the Kretschmann scalar"]],
        "caption": [
            "The Kretschmann scalar $K$ of the Curzon-Chazy particle against $r$ on the axis, $\\theta = 0$, and in "
            "the plane $\\theta = \\pi/2$ ($m = 1$), the same at every $t$ and $\\phi$. Far out both fall off as "
            "$48m^2/r^6$, as Schwarzschild's does at the same mass.",
            "On the axis $K = 48m^2(r - m)^2e^{-4m/r}/r^8$. It falls to zero as $e^{-4m/r}$ toward the singularity at "
            "$r = 0$, vanishes at $r = m$, and has a maximum on either side of that zero, "
            f"${peak_K[0]:.3f}/m^4$ at $r = (1 - 1/\\sqrt{{3}})\\,m = {peaks[0]:.3f}\\,m$ and "
            f"${peak_K[1]:.4f}/m^4$ at $r = (1 + 1/\\sqrt{{3}})\\,m = {peaks[1]:.3f}\\,m$.",
            "In the plane $K$ falls from the singularity outward, and toward $r = 0$ it grows without bound, as "
            "$e^{2m^2/r^2}$.",
        ],
        "settings": "$m = GM/c^2 = 1$, the unit of every length.",
        "clip": ("We draw $K\\,m^4$ on a logarithmic scale from $10^{-6}$ to $10^{6}$ on both panels and cut each curve "
                 f"where it leaves its panel: on the axis inside $r = {cut_axis:.3f}\\,m$ and within "
                 f"${tex_number(dip, 2)}\\,m$ of $r = m$, and in the plane inside $r = {cut_plane:.3f}\\,m$."),
    }
    return [meridian, line]


# ---------------------------------------------------------------- Damour and Solodukhin's wormhole

DS_LAMBDA = Fraction(1, 5)        # lambda, as every drawing of the wormhole takes it
DS_STEPS = 8                      # the frames of the movie from lambda = 1/5 to 1, and as many from 1/25 to 1/5
DS_REACH = 2.2                    # the plot runs over |u| <= 2.2 sqrt(r_s), the box of the Einstein-Rosen diagram


def damour_solodukhin(ck, src):
    """Damour and Solodukhin's wormhole at r_s = 1, its Kretschmann and Ricci scalars against
    Einstein and Rosen's u, through the throat.

    Every chart publishes the same K and R, functions of the areal r alone, and the spherical
    chart covers one side of the throat, so the figure reads them in the Einstein-Rosen chart,
    r = r_s + u^2, whose u runs through the throat at u = 0 from one side to the other. Two panels
    share the axis of u, from -2.2 to 2.2 sqrt(r_s), the box of the Einstein-Rosen spacetime
    diagram: K on a logarithmic axis from 10^-4 to 10^6, and R, negative at every u, on a signed
    one, linear within 10^-6 of zero, from -10^3 to 0, both ticked at every second decade, which
    hold every curve of both views whole, the highest K a decade below the top. The comparison is
    the Einstein-Rosen bridge's K, 12 r_s^2/(u^2 + r_s)^6, Schwarzschild's black hole of the same
    r_s in the same u, read from einstein_rosen_bridge.json; its R vanishes and has no place on
    the signed axis's decades. No
    horizon is marked, since the published g^uu = 1/4(u^2 + r_s) is positive at every u. The first
    view is the wormhole at lambda = 1/5, the value every other drawing of it takes; the second
    steps lambda from 1/25 to 1 by the factor 5^(1/8), 17 frames, through 1/5.

    What the captions state is checked here: the Einstein-Rosen chart's K and R are the spherical
    chart's at r = r_s + u^2, and the bridge's K is Schwarzschild's there and the wormhole's at
    lambda = 0; K and R are even in u; at the throat K = 6/r_s^4 + 1/4 lambda^4 r_s^4 and
    R = -1/2 lambda^2 r_s^2, in every frame; K falls and R rises toward zero from the throat
    outward, and R is negative at every u; |R| is half its throat value at the u stated; K crosses
    the bridge's once on each side, at the u and r stated, above it inside and below it outside,
    and at that r alone of all r > r_s, the one root there of the numerator of their difference;
    far out K tends to (6 + 6/(1 + lambda^2)^2) r_s^2/r^6 while R falls off as
    -lambda^2 r_s^2/2(1 + lambda^2)^2 r^4; in the movie K at the throat stands above the bridge's
    12/r_s^4 for lambda < 24^(-1/4) alone, K approaches the bridge's outside |u| = 1/2 as lambda
    falls, and at lambda = 1/25 is within a percent of it outside the u stated; and no curve of
    either view leaves its panel."""
    x_axis = axis("$u/\\sqrt{r_s}$", -DS_REACH, DS_REACH)
    K_scale = {"kind": "log", "from": 1e-4, "to": 1e6}
    R_scale = {"kind": "symlog", "from": -1e3, "to": 0, "linear": 1e-6}
    # Both axes span about ten decades, so each is ticked at every second one, as the Brans-Dicke
    # sphere's R is, which the page fits at every size.
    K_ticks = [t for t in value_ticks(K_scale) if round(math.log10(t["at"])) % 2 == 0]
    R_ticks = [t for t in value_ticks(R_scale) if t["at"] == 0 or round(math.log10(-t["at"])) % 2 == 0]
    panels = [{"id": "K", "name": "$K\\,r_s^4$", "aspect": 2, "scale": K_scale, "ticks": K_ticks},
              {"id": "R", "name": "$R\\,r_s^2$", "aspect": 2, "scale": R_scale, "ticks": R_ticks}]
    shape_checks(ck, "damour_solodukhin", [p["aspect"] for p in panels])
    height = {p["id"]: PLOT_W / p["aspect"] for p in panels}
    names = ["r_s", "\\lambda"]
    K = nr_published(src, "damour_solodukhin", "einstein_rosen", "kretschmann", ["u"], names)
    R = nr_published(src, "damour_solodukhin", "einstein_rosen", "ricci_scalar", ["u"], names)
    guu = nr_published(src, "damour_solodukhin", "einstein_rosen", "inverse_metric_components", ["u"], names, ("u", "u"))
    K_r = nr_published(src, "damour_solodukhin", "spherical", "kretschmann", ["r"], names)
    R_r = nr_published(src, "damour_solodukhin", "spherical", "ricci_scalar", ["r"], names)
    K_b = nr_published(src, "einstein_rosen_bridge", "bridge", "kretschmann", ["u"], ["r_s"])
    K_s = nr_published(src, "schwarzschild", "spherical", "kretschmann", ["r"], ["r_s"])
    lam = float(DS_LAMBDA)
    us = np.linspace(-DS_REACH, DS_REACH, SAMPLES)
    out = us[us > 0]

    # The charts and the comparison agree.
    # Sixty-fourths of sqrt(r_s), whose r = r_s + u^2 is a double exactly.
    rng = np.random.default_rng(2007)
    points = rng.integers(-140, 141, 20) / 64
    worst = max(max(abs(float(K.exact([u, 1, lam]) / K_r.exact([1 + u * u, 1, lam])) - 1),
                    abs(float(R.exact([u, 1, lam]) / R_r.exact([1 + u * u, 1, lam])) - 1)) for u in points)
    ck("damour_solodukhin the Einstein-Rosen chart's K and R are the spherical chart's at r = r_s + u^2",
       worst < 1e-50, f"worst {worst:.1e}")
    worst = max(max(abs(float(K_b.exact([u, 1]) / K_s.exact([1 + u * u, 1])) - 1),
                    abs(float(K.exact([u, 1, 0]) / K_b.exact([u, 1])) - 1)) for u in points)
    ck("damour_solodukhin the bridge's K is Schwarzschild's at r = r_s + u^2 and the wormhole's at lambda = 0",
       worst < 1e-50, f"worst {worst:.1e}")
    ck("damour_solodukhin K and R are even in u",
       np.array_equal(K(us, 1, lam), K(-us, 1, lam)) and np.array_equal(R(us, 1, lam), R(-us, 1, lam)))
    ck("the bridge's K at u = 0 is 12/r_s^4", float(K_b.exact([0, 1])) == 12)

    def reference():
        return line_curves(ck, "damour_solodukhin bridge", "K", "reference", K_b, "u", x_axis, K_scale, height["K"],
                           {"r_s": 1})

    def throat(value):
        """K and R at the throat, against 6/r_s^4 + 1/4 lambda^4 r_s^4 and -1/2 lambda^2 r_s^2."""
        k, r = float(K.exact([0, 1, value])), float(R.exact([0, 1, value]))
        ck(f"damour_solodukhin at lambda = {value:.4f} K at the throat is 6/r_s^4 + 1/4 lambda^4 r_s^4",
           abs(k / (6 + 1 / (4 * value ** 4)) - 1) < 1e-12, f"{k}")
        ck(f"damour_solodukhin at lambda = {value:.4f} R at the throat is -1/2 lambda^2 r_s^2",
           abs(r / (-1 / (2 * value ** 2)) - 1) < 1e-12, f"{r}")
        return k, r

    def own(value, name):
        params = {"r_s": 1, "lambda": float(value)}
        curves = (line_curves(ck, name, "K", "K", K, "u", x_axis, K_scale, height["K"], params)
                  + line_curves(ck, name, "R", "R", R, "u", x_axis, R_scale, height["R"], params))
        ck(f"{name} has no horizon on the range", not zeros(ck, f"{name} g^uu", guu, -DS_REACH, DS_REACH, [1, value])
           and np.all(guu(us, 1, value) > 0))
        for label, fn, scale in (("K", K, K_scale), ("R", R, R_scale)):
            up, down = beyond(fn, x_axis, scale, [1, value])
            ck(f"{name} {label} stays inside its panel", not up and not down, f"{up} {down}")
        return curves

    def apart(value):
        """K over the bridge's, less one, as a function of u."""
        return lambda u, *_: K(u, 1, value) / K_b(u, 1) - 1

    curves = reference() + own(DS_LAMBDA, "damour_solodukhin at lambda = 1/5")
    ck("the bridge's K stays inside its panel", not any(beyond(K_b, x_axis, K_scale, [1])))
    # The figures of the captions.
    K0, R0 = throat(lam)
    ck("damour_solodukhin the throat values the caption writes are 649/4 and -25/2",
       6 + 1 / (4 * DS_LAMBDA ** 4) == Fraction(649, 4) and -1 / (2 * DS_LAMBDA ** 2) == Fraction(-25, 2))
    ck("damour_solodukhin K falls and R rises from the throat outward, and R is negative at every u",
       np.all(np.diff(K(out, 1, lam)) < 0) and np.all(np.diff(R(out, 1, lam)) > 0) and np.all(R(us, 1, lam) < 0))
    half = zeros(ck, "damour_solodukhin |R| = |R_0|/2", lambda u, *_: R(u, 1, lam) / R0 - 0.5, 0, DS_REACH, [])
    ck("damour_solodukhin |R| falls to half its throat value once on each side", len(half) == 1, f"{half}")
    u_half = half[0]
    crosses = zeros(ck, "damour_solodukhin K = K_bridge", apart(lam), 0, DS_REACH, [])
    ck("damour_solodukhin K crosses the bridge's once on each side", len(crosses) == 1, f"{crosses}")
    u_cross = crosses[0]
    ck("damour_solodukhin K stands above the bridge's inside the crossing and below it outside",
       np.all(apart(lam)(us[np.abs(us) < u_cross]) > 0) and np.all(apart(lam)(us[np.abs(us) > u_cross]) < 0))
    # Beyond the range as well: in r = r_s + u^2 the difference of the two is a polynomial in r over
    # a positive denominator, with one root above r_s alone.
    r = sp.Symbol("r", positive=True)
    u_sym, rs_sym, lam_sym = K.coords + K.params
    gap = (K.expr.subs({u_sym: sp.sqrt(r - 1), rs_sym: 1, lam_sym: DS_LAMBDA})
           - K_b.expr.subs({K_b.coords[0]: sp.sqrt(r - 1), K_b.params[0]: 1}))
    numerator, _ = sp.fraction(sp.factor(sp.together(gap)))
    roots = sp.Poly(numerator, r).count_roots(1, None)
    ck("damour_solodukhin K crosses the bridge's at one r > r_s alone", roots == 1
       and abs(float(sp.Poly(numerator, r).intervals(inf=1, eps=1e-15)[0][0][0]) - (1 + u_cross ** 2)) < 1e-9,
       f"{roots} roots")
    far = math.sqrt(1e4 - 1)
    tail_K, tail_R = 6 + 6 / (1 + lam ** 2) ** 2, -lam ** 2 / (2 * (1 + lam ** 2) ** 2)
    ck("damour_solodukhin K tends to (6 + 6/(1 + lambda^2)^2) r_s^2/r^6 far out",
       abs(float(K.exact([far, 1, lam])) * 1e4 ** 6 / tail_K - 1) < 1e-3, f"{float(K.exact([far, 1, lam])) * 1e24:.6f}")
    ck("damour_solodukhin R falls off as -lambda^2 r_s^2/2(1 + lambda^2)^2 r^4",
       abs(float(R.exact([far, 1, lam])) * 1e4 ** 4 / tail_R - 1) < 1e-3, f"{float(R.exact([far, 1, lam])) * 1e16:.6f}")
    legend = [["line", "K", "$K$, the Kretschmann scalar"], ["line", "R", "$R$, the Ricci scalar"],
              ["line", "reference", "$K$ of Schwarzschild's black hole of the same $r_s$"]]
    clip = ("We draw $K\\,r_s^4$ on a logarithmic scale from $10^{-4}$ to $10^{6}$ and $R\\,r_s^2$ on a signed "
            "logarithmic scale, linear within $10^{-6}$ of zero, from $-10^{3}$ to $0$; no curve leaves its panel.")
    still = {
        "id": "throat", "label": "$\\lambda = 1/5$", "kind": "line", "x": x_axis, "panels": panels,
        "curves": curves, "marks": [], "legend": legend,
        "caption": [
            "The Kretschmann scalar $K$ and the Ricci scalar $R$ of the Damour-Solodukhin wormhole against Einstein "
            "and Rosen's $u$ ($\\lambda = 1/5$), the same at every $t$, $\\theta$, and $\\phi$, beside Schwarzschild's "
            "black hole of the same $r_s$. The areal radius is $r = r_s + u^2$, so $u$ runs from one side of the "
            "throat $u = 0$ to the other. At the throat $K = 6/r_s^4 + 1/4\\lambda^4r_s^4 = "
            f"{K0:g}/r_s^4$, where Schwarzschild's horizon has $K = 12/r_s^4$, and $K$ stands above Schwarzschild's "
            f"only within $|u| < {u_cross:.3f}\\sqrt{{r_s}}$, inside the sphere $r = {1 + u_cross ** 2:.3f}\\,r_s$.",
            f"At the throat $R = -1/2\\lambda^2r_s^2 = {R0:g}/r_s^2$, and $|R|$ falls to half of that at "
            f"$|u| = {u_half:.3f}\\sqrt{{r_s}}$. $R$ is negative at every $u$. Far out $K$ approaches "
            "$(6 + 6/(1 + \\lambda^2)^2)r_s^2/r^6$, "
            f"${tail_K / 12:.3f}$ of Schwarzschild's, and $R$ falls off as $-\\lambda^2r_s^2/2(1 + \\lambda^2)^2r^4$.",
        ],
        "settings": "$r_s = 1$, the unit of every length, and $\\lambda = 1/5$.",
        "clip": clip,
    }

    # The movie in lambda, at fixed r_s, from a wormhole near Schwarzschild's black hole to lambda = 1.
    frames, values = [], []
    for k in range(-DS_STEPS, DS_STEPS + 1):
        value = 5.0 ** (k / DS_STEPS - 1)
        values.append(value)
        frames.append({"label": f"$\\lambda = {value:.3f}$", "value": value,
                       "curves": own(value, f"damour_solodukhin at lambda = {value:.4f}"), "marks": []})
        throat(value)
    first, last = values[0], values[-1]
    ck("damour_solodukhin the movie runs from 1/25 to 1 through 1/5",
       first == 1 / 25 and last == 1.0 and values[DS_STEPS] == lam, f"{first!r} {last!r}")
    meets = 24 ** -0.25
    ck("damour_solodukhin K at the throat stands above the bridge's 12/r_s^4 for lambda < 24^(-1/4) alone",
       all((float(K.exact([0, 1, v])) > 12) == (v < meets) for v in values)
       and abs(float(K.exact([0, 1, meets])) / 12 - 1) < 1e-12, f"{meets:.6f}")
    outer = us[np.abs(us) >= 0.5]
    gaps = [float(np.max(np.abs(apart(v)(outer)))) for v in values]
    ck("damour_solodukhin K approaches the bridge's outside |u| = 1/2 as lambda falls",
       all(a < b for a, b in zip(gaps, gaps[1:])), f"{gaps[0]:.4f} to {gaps[-1]:.4f}")
    fine = np.linspace(0, DS_REACH, 400001)
    far_off = np.abs(apart(first)(fine)) > 0.01
    u_percent = float(fine[np.flatnonzero(far_off)[-1] + 1])
    wide = np.sqrt(np.logspace(math.log10(DS_REACH ** 2 + 1), 8, 2001) - 1)
    ck(f"damour_solodukhin at lambda = 1/25 K is within a percent of the bridge's outside |u| = {u_percent:.4f}",
       np.all(~far_off[fine >= u_percent]) and np.all(np.abs(apart(first)(wide)) <= 0.01))
    K_first, R_first = float(K.exact([0, 1, first])), float(R.exact([0, 1, first]))
    K_last, R_last = float(K.exact([0, 1, last])), float(R.exact([0, 1, last]))
    movie_view = {
        "id": "closing", "label": "$1/25 \\le \\lambda \\le 1$", "kind": "line", "x": x_axis, "panels": panels,
        "curves": reference(), "marks": [], "legend": legend,
        "caption": [
            "The Kretschmann scalar $K$ and the Ricci scalar $R$ of the Damour-Solodukhin wormhole against Einstein "
            "and Rosen's $u$ as $\\lambda$ grows at fixed $r_s$ ($1/25 \\le \\lambda \\le 1$), beside "
            "Schwarzschild's black hole of the same $r_s$. At the throat $K = 6/r_s^4 + 1/4\\lambda^4r_s^4$ falls "
            f"from ${tex_number(K_first, 4)}/r_s^4$ to ${K_last:g}/r_s^4$ and $R = -1/2\\lambda^2r_s^2$ from "
            f"${R_first:g}/r_s^2$ to ${R_last:g}/r_s^2$. $K$ there stands above the $12/r_s^4$ of Schwarzschild's "
            f"horizon only for $\\lambda < 24^{{-1/4}} = {meets:.3f}$.",
            "As $\\lambda$ falls, $K$ approaches Schwarzschild's at every $u$ but the throat, where it grows without "
            "bound, though at $\\lambda = 0$ the throat is Schwarzschild's horizon with $K = 12/r_s^4$. At "
            f"$\\lambda = 1/25$, $K$ is within a percent of Schwarzschild's outside $|u| = {u_percent:.3f}\\sqrt{{r_s}}$, "
            f"the sphere $r = {1 + u_percent ** 2:.3f}\\,r_s$.",
        ],
        "settings": f"$r_s = 1$, the unit of every length, and $\\lambda$ in steps of the factor $5^{{1/{DS_STEPS}}}$.",
        "clip": clip,
        "movie": {"variable": "$\\lambda$", "seconds": 5, "frames": frames},
    }
    return [still, movie_view]


# ---------------------------------------------------------------- Einstein's cluster

EC_SPEED = Fraction(1, 2)         # V, the speed of the particles of Einstein's own cluster, as every drawing takes it
EC_RADIUS = 3                     # R in units of r_s, the radius of both clusters with a surface, as every drawing takes it
EC_REACH = 9                      # the plot runs from the centre to r = 3R
EC_K = (1e-5, 1e2)                # the axis of K r_s^4 on every view
EC_R = (1e-4, 1e1)                # the axis of R^mu_mu r_s^2 on every view
EC_CENTRE = 1e-6                  # the areal chart's scalars are 0/0 at r = 0, so its curves start here, under a pixel out


def einstein_cluster(ck, src):
    """Einstein's cluster at r_s = 1, its Kretschmann and Ricci scalars against the areal radius r,
    for the three clusters its other drawings take: Einstein's own of constant speed at V = 1/2 and
    R = 3 r_s, Florides's of uniform density of the same mass and radius, and the cluster with no
    surface of the areal chart, at the core radius b = 2 r_s declared as null_rays.EC_CORE.

    Every chart is static and spherical, so both scalars depend on the radius alone and each
    view is a line view; a map over r and an angle would be the same colour along every circle.
    The isotropic and hyperspherical charts are the first two clusters in another radius, so the
    clusters are drawn once, against the areal radius, in the constant_speed, uniform, and areal
    charts, and each view belongs to every chart, as the embedding diagram's three views do.
    Two panels share the axis of r from 0 to 3R: K on a logarithmic axis from 10^-5 to 10^2 and
    R, written R^mu_mu since R is the cluster's radius, on a logarithmic one from 10^-4 to 10,
    since it is positive inside every cluster and vanishes outside one. The clusters with a
    surface are drawn from the centre to the surface, where their charts end, and Schwarzschild's
    scalars of the same mass, read from schwarzschild.json, from the surface out, in the same
    classes, since that vacuum is their exterior; its R vanishes, so the R panel has no curve
    there. Each of the two is dashed beside the other inside the surface, and the cluster with no
    surface is drawn across the whole axis with Schwarzschild's K of the same mass dashed.

    What the captions state is checked here: g^rr is positive on every stretch drawn, so no
    horizon is marked; r_s = 2V^2R/(1 + 2V^2) is 1; the density and the pressure along the
    circles are read from the published Einstein tensor, G^t_t = -8 pi G rho/c^2 and G^theta_theta
    8 pi G/c^4 times p_perp, the pressure along the radius vanishing in every chart, and the
    published R is 8 pi G(rho c^2 - 2 p_perp)/c^4, minus the trace; Einstein's cluster has
    K = 35/48r^4, R = 1/4r^2, a density growing as 1/r^2 and p_perp = rho c^2 V^2/2; Florides's
    has a uniform density, K rising from 15 r_s^2/R^6 at the centre to the value stated at the
    surface and R falling from 3 r_s/R^3 to 3 r_s(2R - 3r_s)/(2R^3(R - r_s)); Schwarzschild's
    K = 12 r_s^2/R^6 at the surface and the jump of each scalar there; the cluster with no
    surface has K and R falling from the centre outward, below Schwarzschild's K at every r
    drawn, approaching it far out, and R falling off as 3 r_s b^2/r^5; at each regular centre
    p_perp vanishes and K = 5 R^2/3; and the clip stated is where each curve leaves its panel."""
    x_axis = axis("$r/r_s$", 0, EC_REACH)
    K_scale = {"kind": "log", "from": EC_K[0], "to": EC_K[1]}
    R_scale = {"kind": "log", "from": EC_R[0], "to": EC_R[1]}
    panels = [{"id": "K", "name": "$K\\,r_s^4$", "aspect": 2, "scale": K_scale, "ticks": value_ticks(K_scale)},
              {"id": "R", "name": "$R^\\mu{}_\\mu\\,r_s^2$", "aspect": 2, "scale": R_scale,
               "ticks": value_ticks(R_scale)}]
    shape_checks(ck, "einstein_cluster", [p["aspect"] for p in panels])
    height = {p["id"]: PLOT_W / p["aspect"] for p in panels}
    V, R0, b = float(EC_SPEED), float(EC_RADIUS), 2.0
    speed = {"V": V, "R": R0}
    uniform = {"r_s": 1, "R": R0}
    core = {"Phi": nr.EC_CORE}

    def chart(system, params, functions=None):
        """The published K, R, g^rr, G^t_t and G^theta_theta of one chart, each a function of r."""
        names = list(params)
        out = {field: nr_published(src, "einstein_cluster", system, field, ["r"], names, functions=functions)
               for field in ("kretschmann", "ricci_scalar")}
        out["grr"] = nr_published(src, "einstein_cluster", system, "inverse_metric_components", ["r"], names,
                                  ("r", "r"), functions)
        for i in ("t", "\\theta"):
            out[i] = nr_published(src, "einstein_cluster", system, "einstein_tensor", ["r"], names, (i, i), functions)
        return out

    S, U, C = chart("constant_speed", speed), chart("uniform", uniform), chart("areal", {}, core)
    K_s = nr_published(src, "schwarzschild", "spherical", "kretschmann", ["r"], ["r_s"])
    grr_s = nr_published(src, "schwarzschild", "spherical", "inverse_metric_components", ["r"], ["r_s"], ("r", "r"))
    at = {"constant_speed": [V, R0], "uniform": [1, R0], "areal": []}

    for system in ("constant_speed", "uniform", "areal"):
        _, entry, _ = nr.load("einstein_cluster", system)
        listed = [c["indices"] for c in entry["einstein_tensor"]["variants"]["ul"]["nonzero"]]
        ck(f"einstein_cluster {system} has no pressure along the radius", ["r", "r"] not in listed, f"{listed}")
    ck("einstein_cluster the cluster of constant speed has r_s = 2V^2R/(1 + 2V^2) = 1",
       2 * EC_SPEED ** 2 * EC_RADIUS / (1 + 2 * EC_SPEED ** 2) == 1)

    inside, outside, whole = (0.0, R0), (R0, float(EC_REACH)), (EC_CENTRE, float(EC_REACH))
    rs_in, rs_out = np.linspace(*inside, SAMPLES), np.linspace(*outside, SAMPLES)
    ck("einstein_cluster has no horizon inside either cluster with a surface",
       np.all(S["grr"](rs_in, *at["constant_speed"]) > 0) and np.all(U["grr"](rs_in, *at["uniform"]) > 0))
    ck("einstein_cluster has no horizon outside the surface", np.all(grr_s(rs_out, 1) > 0))
    ck("einstein_cluster has no horizon in the cluster with no surface",
       np.all(C["grr"](np.linspace(*whole, SAMPLES)) > 0))

    def scalars(name, fns, system, params, cls, within, functions=None):
        """The K and R a chart draws over `within`, in the classes K and R or both as `reference`."""
        published = {"K": fns["kretschmann"], "R": fns["ricci_scalar"]}
        return [curve for panel, scale in (("K", K_scale), ("R", R_scale))
                for curve in line_curves(ck, f"{name} {panel}", panel, cls if cls == "reference" else panel,
                                         published[panel], "r", x_axis, scale, height[panel], params, within)]

    def vacuum(name):
        """Schwarzschild's K of the same mass from the surface to the edge of the plot; its R vanishes."""
        return line_curves(ck, f"{name} outside", "K", "K", K_s, "r", x_axis, K_scale, height["K"], {"r_s": 1},
                           outside)

    def cut(name, fn, params, within, scale):
        """The radius inside which a curve passes the top of its panel, found by bisection to the last
        digit, or None where it stays inside; a curve that passes the bottom, or the top anywhere but
        about the centre, fails the check."""
        up, down = beyond(fn, x_axis, scale, params, within)
        ck(f"einstein_cluster {name} leaves its panel about the centre alone",
           not down and (not up or (len(up) == 1 and up[0][0] == within[0])), f"{up} {down}")
        if not up:
            return None
        top = scale["to"]
        return crossing(lambda r: float(fn(r, *params)) - top, up[0][1], up[0][1] + (within[1] - within[0]) / (SAMPLES - 1))

    def exact(fn, system, r):
        return float(fn.exact([r] + at[system]))

    # Einstein's own cluster, of constant speed.
    rho_S = [-exact(S["t"], "constant_speed", r) for r in (0.1, 1, R0)]
    ck("einstein_cluster the density of the cluster of constant speed grows as 1/r^2",
       all(abs(rho * r ** 2 * 3 - 1) < 1e-12 for rho, r in zip(rho_S, (0.1, 1, R0))), f"{rho_S}")
    ck("einstein_cluster the pressure along the circles is rho c^2 V^2/2 in the cluster of constant speed",
       all(abs(exact(S["\\theta"], "constant_speed", r) / -exact(S["t"], "constant_speed", r) - V ** 2 / 2) < 1e-12
           for r in (0.1, 1, R0)))
    ck("einstein_cluster K of the cluster of constant speed is 35/48r^4 and R is 1/4r^2",
       all(abs(exact(S["kretschmann"], "constant_speed", r) * r ** 4 * 48 / 35 - 1) < 1e-12
           and abs(exact(S["ricci_scalar"], "constant_speed", r) * r ** 2 * 4 - 1) < 1e-12 for r in (0.01, 0.5, 2, R0)))
    surface_S = exact(S["kretschmann"], "constant_speed", R0)
    vacuum_K = float(K_s.exact([R0, 1]))
    ck("einstein_cluster Schwarzschild's K at the surface is 12 r_s^2/R^6", abs(vacuum_K * R0 ** 6 / 12 - 1) < 1e-12)
    ck("einstein_cluster K of the cluster of constant speed jumps up at the surface", surface_S < vacuum_K,
       f"{surface_S:.6f} to {vacuum_K:.6f}")

    # Florides's cluster, of uniform density.
    rho_U = -exact(U["t"], "uniform", 0)
    ck("einstein_cluster the density of the cluster of uniform density is 3 r_s/R^3 throughout",
       abs(rho_U * R0 ** 3 / 3 - 1) < 1e-12
       and np.allclose(-U["t"](rs_in, *at["uniform"]), rho_U, rtol=1e-13, atol=0), f"{rho_U}")
    centre_KU, centre_RU = exact(U["kretschmann"], "uniform", 0), exact(U["ricci_scalar"], "uniform", 0)
    surface_KU, surface_RU = exact(U["kretschmann"], "uniform", R0), exact(U["ricci_scalar"], "uniform", R0)
    ck("einstein_cluster K at the centre of the cluster of uniform density is 15 r_s^2/R^6",
       abs(centre_KU * R0 ** 6 / 15 - 1) < 1e-12, f"{centre_KU}")
    ck("einstein_cluster R at the centre of the cluster of uniform density is 3 r_s/R^3",
       abs(centre_RU * R0 ** 3 / 3 - 1) < 1e-12, f"{centre_RU}")
    ck("einstein_cluster R at the surface of the cluster of uniform density is 3 r_s(2R - 3r_s)/(2R^3(R - r_s))",
       abs(surface_RU / (3 * (2 * R0 - 3) / (2 * R0 ** 3 * (R0 - 1))) - 1) < 1e-12, f"{surface_RU}")
    ck("einstein_cluster K rises and R falls from the centre of the cluster of uniform density to its surface",
       np.all(np.diff(U["kretschmann"](rs_in, *at["uniform"])) > 0)
       and np.all(np.diff(U["ricci_scalar"](rs_in, *at["uniform"])) < 0))
    ck("einstein_cluster K of the cluster of uniform density jumps down at the surface", surface_KU > vacuum_K,
       f"{surface_KU:.6f} to {vacuum_K:.6f}")

    # The cluster with no surface.
    centre_KC, centre_RC = float(C["kretschmann"].exact([1e-20])), float(C["ricci_scalar"].exact([1e-20]))
    ck("einstein_cluster K at the centre of the cluster with no surface is 15/16r_s^4",
       abs(centre_KC * 16 / 15 - 1) < 1e-12, f"{centre_KC}")
    ck("einstein_cluster R at the centre of the cluster with no surface is 3/4r_s^2",
       abs(centre_RC * 4 / 3 - 1) < 1e-12, f"{centre_RC}")
    rs_whole = np.linspace(*whole, SAMPLES)
    ck("einstein_cluster K and R of the cluster with no surface fall from the centre outward",
       np.all(np.diff(C["kretschmann"](rs_whole)) < 0) and np.all(np.diff(C["ricci_scalar"](rs_whole)) < 0))
    ck("einstein_cluster K of the cluster with no surface lies below Schwarzschild's at every r drawn",
       np.all(C["kretschmann"](rs_whole) < K_s(rs_whole, 1)))
    ratio_edge = float(C["kretschmann"].exact([float(EC_REACH)])) / float(K_s.exact([float(EC_REACH), 1]))
    far = [float(C["kretschmann"].exact([r])) / float(K_s.exact([r, 1])) for r in (1e2, 1e3, 1e4)]
    ck("einstein_cluster K of the cluster with no surface approaches Schwarzschild's far out",
       abs(far[2] - 1) < 1e-6 and abs(far[2] - 1) < abs(far[1] - 1) < abs(far[0] - 1), f"{far}")
    tail = [float(C["ricci_scalar"].exact([r])) * r ** 5 / (3 * b ** 2) for r in (1e2, 1e3, 1e4, 1e5, 1e6)]
    ck("einstein_cluster R of the cluster with no surface falls off as 3 r_s b^2/r^5",
       abs(tail[-1] - 1) < 1e-5 and all(abs(c - 1) < abs(a - 1) for a, c in zip(tail, tail[1:])), f"{tail}")

    # At each regular centre the pressure vanishes and the Weyl tensor with it, so K = 5 R^2/3.
    ck("einstein_cluster the pressure along the circles vanishes at each regular centre",
       abs(exact(U["\\theta"], "uniform", 0)) < 1e-30 and abs(float(C["\\theta"].exact([1e-20]))) < 1e-30)
    ck("einstein_cluster K = 5 R^2/3 at each regular centre",
       abs(centre_KU / (5 * centre_RU ** 2 / 3) - 1) < 1e-12 and abs(centre_KC / (5 * centre_RC ** 2 / 3) - 1) < 1e-12)

    for fns, system, params, within in ((S, "constant_speed", speed, inside), (U, "uniform", uniform, inside),
                                        (C, "areal", {}, whole)):
        for name, fn in (("kretschmann", fns["kretschmann"]), ("ricci_scalar", fns["ricci_scalar"])):
            ck(f"einstein_cluster {system} {name} is positive on the stretch drawn",
               np.all(fn(np.linspace(within[0] + 1e-9, within[1], SAMPLES), *at[system]) > 0))
    top_KS = cut("constant speed K", S["kretschmann"], at["constant_speed"], inside, K_scale)
    top_RS = cut("constant speed R", S["ricci_scalar"], at["constant_speed"], inside, R_scale)
    top_Ks = cut("Schwarzschild K", K_s, [1], (0.0, float(EC_REACH)), K_scale)
    ck("einstein_cluster the curves of constant speed and Schwarzschild's K pass the top at "
       "(35/4800)^(1/4), 1/sqrt 40, and (12/100)^(1/6)",
       top_KS is not None and top_RS is not None and top_Ks is not None
       and abs(top_KS / (35 / 4800) ** 0.25 - 1) < 1e-12 and abs(top_RS * math.sqrt(40) - 1) < 1e-12
       and abs(top_Ks / 0.12 ** (1 / 6) - 1) < 1e-12, f"{top_KS} {top_RS} {top_Ks}")
    for fns, system, within in ((U, "uniform", inside), (C, "areal", whole)):
        for name, scale in (("kretschmann", K_scale), ("ricci_scalar", R_scale)):
            ck(f"einstein_cluster {system} {name} stays inside its panel",
               cut(f"{system} {name}", fns[name], at[system], within, scale) is None)
    ck("einstein_cluster Schwarzschild's K outside the surface stays inside its panel",
       cut("Schwarzschild K outside", K_s, [1], outside, K_scale) is None)

    K_clip = (f"We draw $K\\,r_s^4$ on a logarithmic scale from $10^{{{round(math.log10(EC_K[0]))}}}$ to "
              f"$10^{{{round(math.log10(EC_K[1]))}}}$ and $R^\\mu{{}}_\\mu\\,r_s^2$ on one from "
              f"$10^{{{round(math.log10(EC_R[0]))}}}$ to ${EC_R[1]:g}$")
    own = [["line", "K", "$K$, the Kretschmann scalar"], ["line", "R", "$R^\\mu{}_\\mu$, the Ricci scalar"]]
    settings = "$r_s = 1$, the unit of every length, $V = 1/2$, and $R = 3\\,r_s$."
    speed_view = {
        "id": "speed", "label": "Constant speed", "kind": "line", "x": x_axis, "panels": panels,
        "curves": (scalars("einstein_cluster uniform beside", U, "uniform", uniform, "reference", inside)
                   + scalars("einstein_cluster constant speed", S, "constant_speed", speed, "own", inside)
                   + vacuum("einstein_cluster constant speed")),
        "marks": [], "legend": own + [["line", "reference", "Florides's cluster of uniform density, of the same mass "
                                                            "and radius"]],
        "caption": [
            "The Kretschmann scalar $K$ and the Ricci scalar $R^\\mu{}_\\mu$ of Einstein's cluster of constant speed "
            "against $r$ ($V = 1/2$, $R = 3\\,r_s$), the same at every $t$, $\\theta$, and $\\phi$, inside the cluster "
            "and in Schwarzschild's vacuum outside it, beside Florides's cluster of uniform density of the same mass "
            "and radius. Inside, $K = 35/(48r^4)$, which grows without bound toward the centre, where the density "
            "grows as $1/r^2$, and is "
            f"${surface_S:#.3g}/r_s^4$ at the surface. Outside it is Schwarzschild's $12r_s^2/r^6$, which is "
            f"${vacuum_K:#.3g}/r_s^4$ at $r = R$.",
            "The Ricci scalar is $R^\\mu{}_\\mu = 8\\pi G(\\rho c^2 - 2p_\\perp)/c^4$, with no pressure along the "
            "radius and the pressure $p_\\perp = \\rho c^2V^2/2$ along the circles of the particles, so "
            f"$R^\\mu{{}}_\\mu = 1/(4r^2)$ inside, which is ${1 / (4 * R0 ** 2):#.3g}/r_s^2$ at the surface. Outside it "
            "vanishes. Both scalars jump at the surface, where the density falls to zero.",
        ],
        "settings": settings,
        "clip": (K_clip + ", and cut each curve where it leaves its panel, as $K$ does inside "
                 f"$r = {top_KS:.3f}\\,r_s$ and $R^\\mu{{}}_\\mu$ inside $r = {top_RS:.3f}\\,r_s$."),
    }
    uniform_view = {
        "id": "uniform", "label": "Uniform density", "kind": "line", "x": x_axis, "panels": panels,
        "curves": (scalars("einstein_cluster constant speed beside", S, "constant_speed", speed, "reference", inside)
                   + scalars("einstein_cluster uniform", U, "uniform", uniform, "own", inside)
                   + vacuum("einstein_cluster uniform")),
        "marks": [], "legend": own + [["line", "reference", "Einstein's cluster of constant speed, of the same mass "
                                                            "and radius"]],
        "caption": [
            "The Kretschmann scalar $K$ and the Ricci scalar $R^\\mu{}_\\mu$ of Florides's cluster of uniform density "
            "against $r$ ($R = 3\\,r_s$), the same at every $t$, $\\theta$, and $\\phi$, inside the cluster and in "
            "Schwarzschild's vacuum outside it, beside Einstein's cluster of constant speed of the same mass and "
            "radius. Inside, $K$ rises from "
            f"$15r_s^2/R^6 = {centre_KU:#.3g}/r_s^4$ at the centre to ${surface_KU:#.3g}/r_s^4$ at the surface. Outside "
            f"it is Schwarzschild's $12r_s^2/r^6$, which is ${vacuum_K:#.3g}/r_s^4$ at $r = R$.",
            "The Ricci scalar is $R^\\mu{}_\\mu = 8\\pi G(\\rho c^2 - 2p_\\perp)/c^4$, with no pressure along the "
            "radius and the pressure $p_\\perp$ along the circles of the particles. It falls from "
            f"$3r_s/R^3 = {centre_RU:#.3g}/r_s^2$ at the centre, where $p_\\perp = 0$ and "
            f"$K = 5(R^\\mu{{}}_\\mu)^2/3$, to ${surface_RU:#.3g}/r_s^2$ at the surface. Outside it vanishes. Both "
            "scalars jump at the surface, where the density falls to zero.",
        ],
        "settings": "$r_s = 1$, the unit of every length, and $R = 3\\,r_s$.",
        "clip": (K_clip + ", and cut each curve where it leaves its panel, as the $K$ of Einstein's cluster does "
                 f"inside $r = {top_KS:.3f}\\,r_s$ and its $R^\\mu{{}}_\\mu$ inside $r = {top_RS:.3f}\\,r_s$."),
    }
    core_view = {
        "id": "core", "label": "No surface", "kind": "line", "x": x_axis, "panels": panels,
        "curves": (line_curves(ck, "einstein_cluster Schwarzschild beside", "K", "reference", K_s, "r", x_axis,
                               K_scale, height["K"], {"r_s": 1})
                   + scalars("einstein_cluster no surface", C, "areal", {}, "own", whole)),
        "marks": [], "legend": own + [["line", "reference", "Schwarzschild's $K$ for the same mass"]],
        "caption": [
            "The Kretschmann scalar $K$ and the Ricci scalar $R^\\mu{}_\\mu$ of a cluster with no surface against $r$ "
            "($b = 2\\,r_s$), the same at every $t$, $\\theta$, and $\\phi$, beside Schwarzschild's $K$ for the same "
            "mass. Both fall from the centre outward with no jump, and $K$ lies below Schwarzschild's $12r_s^2/r^6$ "
            f"at every $r$, ${ratio_edge:.2f}$ of it at $r = 9\\,r_s$, and approaches it far out.",
            "At the centre the pressure $p_\\perp$ along the circles of the particles vanishes, as it does at the "
            "centre of Florides's cluster, and $R^\\mu{}_\\mu = 3/(4r_s^2)$ and $K = 5(R^\\mu{}_\\mu)^2/3 = 15/(16r_s^4)$. "
            "Far out $R^\\mu{}_\\mu$ falls off as $3r_sb^2/r^5$.",
        ],
        "settings": "$r_s = 1$, the unit of every length, and $b = 2\\,r_s$.",
        "input": nr.EC_CORE_INPUT,
        "clip": (K_clip + ", and cut each curve where it leaves its panel, as Schwarzschild's $K$ does inside "
                 f"$r = {top_Ks:.3f}\\,r_s$."),
    }
    return [speed_view, uniform_view, core_view]


# ---------------------------------------------------------------- the files

DRAWN = {
    "bach_weyl_ring": bach_weyl_ring,
    "bardeen": bardeen,
    "bonnor_vaidya": bonnor_vaidya,
    "born_infeld_charge": born_infeld_charge,
    "boulware_deser": boulware_deser,
    "bowers_liang": bowers_liang,
    "brans_dicke_sphere": brans_dicke_sphere,
    "brill_charged_taub_nut": brill_charged_taub_nut,
    "c_metric": c_metric,
    "curzon_chazy": curzon_chazy,
    "damour_solodukhin": damour_solodukhin,
    "einstein_cluster": einstein_cluster,
}


def draw(metric_id, ck):
    src = Sources()
    views = DRAWN[metric_id](ck, src)
    ids = [v["id"] for v in views]
    ck(f"{metric_id} names each view once", len(set(ids)) == len(ids), f"{ids}")
    return {"metric": metric_id, "source": src.stamps(), "views": views}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--metric", action="append", default=[], help="redraw only this spacetime, repeatable")
    parser.add_argument("--verify", action="store_true", help="run and print every check instead of writing")
    args = parser.parse_args(argv)
    unknown = set(args.metric) - set(DRAWN)
    if unknown:
        parser.error(f"no curvature figure is drawn for {sorted(unknown)}")
    wanted = [m for m in DRAWN if not args.metric or m in args.metric]
    start = time.time()
    ck = Checks()
    files = {metric_id: draw(metric_id, ck) for metric_id in wanted}
    if args.verify:
        ck.report()
    failed = ck.failures()
    print(f"{len(ck.made)} checks in {time.time() - start:.1f} s, {len(failed)} failed",
          file=sys.stderr if failed else sys.stdout)
    for name in failed:
        print(f"FAILED {name}", file=sys.stderr)
    if failed or args.verify:
        return 1 if failed else 0
    CURVATURE_DIR.mkdir(parents=True, exist_ok=True)
    if not args.metric:
        # A spacetime that is not drawn has no file, so a full redraw removes one left behind.
        for path in sorted(CURVATURE_DIR.glob("*.json")):
            if path.stem not in DRAWN and not build.CONFLICT_COPY.search(path.stem):
                path.unlink()
                print(f"removed {path.relative_to(build.ROOT)}")
    for metric_id, data in files.items():
        path = CURVATURE_DIR / f"{metric_id}.json"
        path.write_text(json.dumps(data, ensure_ascii=False, separators=(",", ":")) + "\n", encoding="utf-8")
        print(f"wrote {path.relative_to(build.ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
