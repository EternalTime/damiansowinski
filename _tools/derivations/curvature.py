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
       each cell of an even grid over the box, with points and curves marked on it.
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
    """A published value of one chart, a curvature scalar, a component of the inverse metric or
    one of the Riemann tensor with its first index raised, as a numpy function of the
    coordinates drawn and then the parameters, in the order named.

    The value is read through the checker's Reader, every name the reader holds as a function
    written out, and c set to 1. A coordinate or a parameter the value depends on and the figure
    does not name stops the script, since the figure would then hold only at one value of it."""

    def __init__(self, sources, metric_id, system_id, field, coords, params, indices=None):
        _, entry, reader = nr.load(metric_id, system_id)
        self.metric, self.system, self.field = metric_id, system_id, field
        if field in ("inverse_metric_components", "riemann"):
            listed = entry[field] if field == "inverse_metric_components" else entry[field]["variants"]["ulll"]["nonzero"]
            value = next(c["value"] for c in listed if c["indices"] == list(indices))
            expr = reader(value)
        else:
            expr = reader(nr.strip_lhs(entry[field]))
        sources.note(metric_id, system_id, ["coords", "parameters", field])
        if reader.held:
            expr = expr.subs(reader.held).doit()
        expr = expr.subs(reader.c, 1)
        by_plain = {reader._plain(name): symbol for name, symbol in reader.symbol.items()}
        self.coords = [by_plain[name] for name in coords]
        self.params = [reader.parameters[name] for name in params]
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


def sampled(ck, name, published, x_axis, scale, height, params):
    """A published scalar along the drawn coordinate, as the polylines a panel draws.

    The scalar is sampled at SAMPLES points across the axis, and kept where it is finite and
    within a quarter of the panel's height of its box, each stretch carried one sample further
    each way where that sample is finite, so that the line drawn runs to the edge it leaves by.
    Each stretch is thinned in the panel's own units, PLOT_W across and `height` up, so that the
    line drawn lies within TOLERANCE of the samples."""
    lo, hi = x_axis["from"], x_axis["to"]
    xs = np.linspace(lo, hi, SAMPLES)
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
        out.append([[float(x), float(v)] for x, v in zip(xs[keep], vs[keep])])
    ck(f"{name} has a curve to draw", out, f"{len(out)} stretches")
    if scale["kind"] == "log":
        ck(f"{name} is positive where a logarithmic axis draws it",
           all(v > 0 for line in out for _, v in line))
    return out


def beyond(published, x_axis, scale, params):
    """Where a scalar lies above its axis and where below it, as the intervals of the drawn
    coordinate, from the same samples sampled() takes."""
    xs = np.linspace(x_axis["from"], x_axis["to"], SAMPLES)
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


def cells(published, x_axis, y_axis, nx, ny, params):
    """A published scalar at the centres of an even grid of nx by ny cells over a map's box, a
    row for each cell up, from the bottom, each with a value for each cell across."""
    xs = x_axis["from"] + (np.arange(nx) + 0.5) * (x_axis["to"] - x_axis["from"]) / nx
    ys = y_axis["from"] + (np.arange(ny) + 0.5) * (y_axis["to"] - y_axis["from"]) / ny
    X, Y = np.meshgrid(xs, ys)
    return xs, ys, published(X, Y, *params)


def rounded(v, digits=6):
    """A value written to `digits` significant figures, which the page and the checks read as
    the value at its point."""
    return float(f"{v:.{digits}g}")


def shape_checks(ck, where, aspects):
    for k, a in enumerate(aspects):
        ck(f"{where} panel {k} is between 1:2 and 2:1", TALLEST <= a <= WIDEST, f"{a:.3f}")
    whole = 1 / sum(1 / a for a in aspects)
    ck(f"{where} is between 1:2 and 2:1 as a whole", TALLEST <= whole <= WIDEST, f"{whole:.3f}")


def curve_entry(panel, cls, published, coordinate, params, points):
    """A curve as the file holds it: what it is drawn from, so that every value can be evaluated
    again, and its points [x, value]."""
    return {"panel": panel, "class": cls, "metric": published.metric, "system": published.system,
            "scalar": published.field, "coordinate": coordinate,
            "params": {name: float(value) for name, value in params.items()}, "points": points}


def spread(points, n=24):
    """About n of a curve's points, evenly through it, its two ends among them."""
    flat = [p for line in points for p in line]
    picks = sorted({round(k * (len(flat) - 1) / (n - 1)) for k in range(n)})
    return [flat[k] for k in picks]


def line_curves(ck, name, panel, cls, published, coordinate, x_axis, scale, height, params):
    """The curves one published scalar draws on a panel, each checked in sixty digits."""
    values = [params[p] for p in params]
    lines = sampled(ck, f"{name} {cls}", published, x_axis, scale, height, values)
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
            "($g = r_s/3$), the same at every $t$, $\\theta$ and $\\phi$, beside Schwarzschild's of the same mass. "
            "Toward the centre Schwarzschild's $K = 12r_s^2/r^6$ grows without bound, and Bardeen's levels off "
            f"inside $r_-$ and reaches $24r_s^2/g^6 = {centre_K:.0f}/r_s^4$ at $r = 0$, the curvature of de Sitter "
            f"space of radius $\\sqrt{{g^3/r_s}} = {math.sqrt(g ** 3):.3f}\\,r_s$.",
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
                  f"$g < {R_top:.3f}\\,r_s$, and $R$ passes the bottom for $g < {R_bottom:.3f}\\,r_s$.")
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


def nr_published(src, metric_id, system_id, field, coords, params, indices=None):
    """A published value read as Published reads it, its parameters named as the chart spells them."""
    return Published(src, metric_id, system_id, field, coords, [p.replace("\\", "") for p in params], indices)


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
            "vanishes with it. Weyl's $\\rho$ and $z$ are lengths in a flat space that the metric stretches, and the "
            "ring, the circle $\\rho = a$, has an infinite circumference.",
        ],
        "settings": "$m = GM/c^2 = a/2$, with the radius $a$ of the ring in Weyl's coordinates the unit of every "
                    "length.",
        "clip": ("We colour $K\\,a^4$ on a logarithmic scale from $10^{-3}$ to $10^{3}$ and give a value beyond it "
                 f"the colour of the end it passes: $K$ passes $10^{{3}}/a^4$ within ${d_high:.2f}\\,a$ of the ring "
                 f"and falls below $10^{{-3}}/a^4$ within ${d_low:.2f}\\,a$ of the two points where it vanishes."),
    }
    return [view]


# ---------------------------------------------------------------- the files

DRAWN = {
    "bach_weyl_ring": bach_weyl_ring,
    "bardeen": bardeen,
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
