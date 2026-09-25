#!/usr/bin/env python3
"""Light cones and light rays in three coordinates at once, drawn from the published metrics
and projected onto the page from one fixed point of view.

Where the causal structure a reader comes for turns in a direction no plane of two
coordinates holds, as the light cones of van Stockum's cylinder tip over about its axis, a
plane is not enough. A figure here is a slice of three coordinates of one coordinate system,
one time and two of space, every other coordinate held fixed, placed in Euclidean
coordinates (X, Y, T) of the drawing with T up:

    polar      (t, r, phi) as X = r cos phi, Y = r sin phi, T = ct
    cartesian  (t, x, y)   as X = x, Y = y, T = ct

and projected orthographically from a camera at a fixed azimuth and elevation. What reaches
the file is the projection: polylines, polygons and points in the plane of the page, in the
order they are painted, and TeX labels at points, the form a conformal diagram takes. So the
page draws a figure as it draws every other one, it prints, and the application draws the
same data with its own renderer; nothing is left for a reader to turn.

null_rays.py owns the diagram files and writes each figure into its spacetime's file under
`projections`, keyed by coordinate system beside the flat views under `systems`, which an
application that reads only the flat views passes over.


The light cones
---------------

At a point of the slice the published metric on its three coordinates, g, is carried to the
drawing's coordinates, G = J^-T g J^-1 with J = d(X, Y, T)/d(slice), and the future null
directions there are n(beta) = e_0 + cos(beta) e_1 + sin(beta) e_2, for a frame that G makes
orthonormal with e_0 future timelike. Every generator is drawn to one Euclidean length in the
drawing, so a cone's size says nothing and its shape and tilt are the metric's. The future is
chosen as the null rays choose it: 'vector' along the chart's time direction, where d_t is
timelike, and 'tau' along minus the gradient of t, where dt is. Every generator drawn is
checked null against the published metric, and the published inverse metric on the slice
is checked to be the inverse of the published metric there, which holds exactly when the
slice is orthogonal to the coordinates held fixed.

A cone is painted as the convex hull of its apex and its rim, filled faintly, with its rim
and the two generators that bound it in the projection drawn over the fill, and RIBS of its
generators drawn faintly from the apex to the rim. A wide cone seen from inside its opening
projects to an oval with its apex inside, and the ribs are what show where the apex is and
which way the cone opens. Cones are painted farthest first.


Output
------

A figure is written as

  id, label       its place and the name of its button, TeX in $...$;
  box             [Xmin, Xmax, Ymin, Ymax] of the page's plane, drawn at one scale;
  layers          fills, lines and points, each with a class, painted in the order given;
  labels          TeX at a point, with an anchor and an offset (dx, dy) in units of a
                  figure 628 wide, as a conformal diagram's are;
  legend          [kind, class, text] for every class it names, kind being fill, line,
                  point or cone;
  camera          the azimuth and elevation it was projected from, in degrees;
  caption, settings, input and source, as a flat view carries them.
"""

from dataclasses import dataclass, field

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp

import build_mfs_data as build
import null_rays as nr

LENGTH = 0.32               # a cone's generators, as a fraction of the figure's slice radius
RIM = 96                    # generators per cone
RIBS = 8                    # of them drawn from the apex to the rim
NULL = 1e-12                # how far a drawn generator may miss null, against |g| |k|^2
INVERSE = 1e-12             # how far g^-1 g may miss the identity on the slice


class Camera:
    """An orthographic view of (X, Y, T) from azimuth `azimuth`, measured from X toward Y,
    and elevation `elevation` above the plane T = 0, both in degrees."""

    def __init__(self, azimuth, elevation):
        a, e = np.radians(azimuth), np.radians(elevation)
        self.azimuth, self.elevation = azimuth, elevation
        self.toward = np.array([np.cos(e) * np.cos(a), np.cos(e) * np.sin(a), np.sin(e)])
        self.right = np.array([-np.sin(a), np.cos(a), 0.0])
        self.up = np.cross(self.toward, self.right)

    def screen(self, P):
        P = np.asarray(P, dtype=float)
        return np.stack([P @ self.right, P @ self.up], -1)

    def depth(self, P):
        """How near the viewer a point is: larger is nearer."""
        return np.asarray(P, dtype=float) @ self.toward


@dataclass
class Projection:
    """One figure of one coordinate system: its place, the parameter values and the
    coordinates held fixed it is drawn at, and the function that draws it."""

    metric: str
    system: str
    view: str
    label: str
    build: object                   # spec -> (Figure.done() dict, the Slice it read)
    params: dict = field(default_factory=dict)
    fixed: dict = field(default_factory=dict)
    input: str = None
    fields: tuple = ()              # published fields read beyond FIELDS
    functions: dict = field(default_factory=dict)   # a declared function -> its expression


FIELDS = ["coords", "parameters", "metric_components", "inverse_metric_components"]


class Slice:
    """The published metric of one coordinate system on a slice of three of its coordinates,
    every other held fixed, as numpy functions of those three.

    coords    the three coordinates as the entry spells them, time first;
    layout    'polar' for (t, r, phi) or 'cartesian' for (t, x, y);
    params    parameter -> value; fixed, coordinate held fixed (plain name) -> value;
    functions a declared function -> an expression for it, in plain names.
    """

    def __init__(self, metric_id, system_id, coords, layout, params=None, fixed=None, functions=None):
        _, entry, reader = nr.load(metric_id, system_id)
        self.entry, self.reader, self.layout = entry, reader, layout
        self.metric_id, self.system_id = metric_id, system_id
        names = entry["coords"]
        self.index = [names.index(c) for c in coords]
        self.symbols = [reader.symbol[c] for c in coords]
        by_plain = {reader._plain(n): s for n, s in reader.symbol.items()}
        subs = {reader.c: 1}
        subs.update({reader.parameters[k]: nr.number(v) for k, v in (params or {}).items()})
        subs.update({by_plain[k]: nr.number(v) for k, v in (fixed or {}).items()})
        held = sorted(set(range(len(names))) - set(self.index))
        if sorted(reader._plain(names[i]) for i in held) != sorted(fixed or {}):
            raise SystemExit(f"{metric_id}/{system_id}: every coordinate off the slice must be held fixed")

        def prep(expr):
            expr = sp.sympify(expr)
            for name, rep in (functions or {}).items():
                expr = expr.replace(reader.parameters[name].func, nr._as_lambda(reader, name, rep)).doit()
            return expr.subs(subs)

        g = nr.published_matrix(reader, entry, "metric_components")
        gi = nr.published_matrix(reader, entry, "inverse_metric_components")
        block = [[prep(g[i, j]) for j in self.index] for i in self.index]
        inverse = [[prep(gi[i, j]) for j in self.index] for i in self.index]
        self._g = sp.lambdify(self.symbols, block, "numpy")
        self._gi = sp.lambdify(self.symbols, inverse, "numpy")
        self.prep = prep
        self.rho = None

    def proper_radius(self, r_max, n=200001):
        """Draw a polar slice with the proper distance from the axis as its radius,
        rho = integral of sqrt(g_rr) dr at fixed t and phi, from the published g_rr, so that
        light moving straight out runs at 45 degrees as it does at the axis. g_rr must not
        depend on t or phi, and the slice's r must be orthogonal to them there."""
        r = np.linspace(0.0, r_max, n)
        g = np.array([self.metric((0.0, x, 0.0)) for x in r[:: n // 400]])
        if np.abs(g[:, 1, [0, 2]]).max() > 0 or np.abs(np.array(
                [self.metric((1.3, x, 0.7))[1, 1] for x in r[:: n // 400]]) - g[:, 1, 1]).max() > 0:
            raise SystemExit(f"{self.metric_id}/{self.system_id}: r is not a proper radius here")
        speed = np.sqrt(np.array([self.metric((0.0, x, 0.0))[1, 1] for x in r]))
        rho = np.concatenate([[0.0], np.cumsum(0.5 * (speed[1:] + speed[:-1]) * np.diff(r))])
        self.rho = (r, rho, speed)

    def metric(self, x):
        """The published metric on the slice at one point x of its three coordinates."""
        return np.array(self._g(*x), dtype=float)

    def inverse(self, x):
        return np.array(self._gi(*x), dtype=float)

    def to_drawing(self, x):
        """A point of the slice in the drawing's (X, Y, T)."""
        t, a, b = (np.asarray(v, dtype=float) for v in x)
        if self.layout == "polar":
            a = self.radius(a)
            return np.stack([a * np.cos(b), a * np.sin(b), t], -1)
        return np.stack([a, b, t], -1)

    def radius(self, r):
        """The drawn radius of a polar slice: r itself, or the proper distance from the axis."""
        if self.rho is None:
            return r
        return np.interp(r, self.rho[0], self.rho[1])

    def jacobian(self, x):
        """d(X, Y, T)/d(slice) at one point."""
        t, a, b = x
        if self.layout == "polar":
            d = 1.0 if self.rho is None else float(np.sqrt(self.metric((t, a, b))[1, 1]))
            a = self.radius(a)
            return np.array([[0.0, d * np.cos(b), -a * np.sin(b)],
                             [0.0, d * np.sin(b), a * np.cos(b)],
                             [1.0, 0.0, 0.0]])
        return np.array([[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]])


def future_cone(sl, x, length, orient="vector", n=RIM):
    """The rim of the future light cone at the slice point x: n points of the drawing, each
    at Euclidean distance `length` from the apex along one future null generator."""
    gen, _ = generators(sl, x, orient, n)
    apex = sl.to_drawing(x)
    return apex, apex[None, :] + length * gen / np.linalg.norm(gen, axis=1, keepdims=True)


def generators(sl, x, orient="vector", n=RIM):
    """n future null generators at the slice point x, in the drawing's (X, Y, T) and in the
    slice's own coordinates.

    The generators are checked null against the published metric and the published inverse
    against the published metric's inverse, and a failure stops the figure."""
    g, gi = sl.metric(x), sl.inverse(x)
    miss = np.abs(gi @ g - np.eye(3)).max()
    if not miss <= INVERSE:
        raise SystemExit(f"{sl.metric_id}/{sl.system_id}: the published inverse metric is not the "
                         f"inverse of the published metric on the slice at {x}, by {miss:.1e}")
    J = sl.jacobian(x)
    Ji = np.linalg.inv(J)
    G = Ji.T @ g @ Ji
    if orient == "vector":
        e0 = J @ np.array([1.0, 0.0, 0.0])
    elif orient == "tau":
        e0 = -np.linalg.inv(G) @ (Ji.T @ np.array([1.0, 0.0, 0.0]))
    else:
        raise ValueError(orient)
    norm = e0 @ G @ e0
    if not norm < 0:
        raise SystemExit(f"{sl.metric_id}/{sl.system_id}: the rule '{orient}' gives no future at {x}")
    e0 = e0 / np.sqrt(-norm)
    frame = [e0]
    for trial in np.eye(3):
        # Gram-Schmidt in G, the timelike member dividing by its own negative norm.
        v = trial.copy()
        for f in frame:
            v = v - (f @ G @ trial) / (f @ G @ f) * f
        size = v @ G @ v
        if size > 1e-9 * np.abs(G).max() * (trial @ trial):
            frame.append(v / np.sqrt(size))
        if len(frame) == 3:
            break
    e0, e1, e2 = frame
    beta = np.linspace(0, 2 * np.pi, n, endpoint=False)
    gen = e0[None, :] + np.cos(beta)[:, None] * e1[None, :] + np.sin(beta)[:, None] * e2[None, :]
    k = gen @ Ji.T
    miss = np.abs(np.einsum("na,ab,nb->n", k, g, k)) / (np.abs(g).max() * np.einsum("na,na->n", k, k))
    if not miss.max() <= NULL:
        raise SystemExit(f"{sl.metric_id}/{sl.system_id}: a cone generator misses null by {miss.max():.1e}")
    return gen, k


def hull(points):
    """The convex hull of points of the plane, counterclockwise, by Andrew's monotone chain."""
    P = sorted(map(tuple, np.round(np.asarray(points, dtype=float), 12)))
    if len(P) < 3:
        return np.array(P)

    def cross(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])
    lower, upper = [], []
    for p in P:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)
    for p in reversed(P):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)
    return np.array(lower[:-1] + upper[:-1])


def rounded(points):
    return np.round(np.asarray(points, dtype=float), 4).tolist()


class Figure:
    """A projected figure: layers painted in the order they are added, labels, a legend."""

    def __init__(self, vid, label, camera):
        self.id, self.name, self.camera = vid, label, camera
        self.layers, self.labels, self.legend_items = [], [], []

    def line(self, cls, P, closed=False):
        """A polyline of the drawing's (X, Y, T), projected."""
        S = self.camera.screen(P)
        if closed:
            S = np.vstack([S, S[:1]])
        self.layers.append({"kind": "line", "class": cls, "points": rounded(nr.thin(S, 0.0005))})

    def fill(self, cls, S):
        """A polygon already in the plane of the page."""
        self.layers.append({"kind": "fill", "class": cls, "points": rounded(S)})

    def point(self, cls, P):
        self.layers.append({"kind": "point", "class": cls, "at": rounded(self.camera.screen(P))})

    def cone(self, apex, rim, cls="cone"):
        """A cone of its apex and rim, both in the drawing: the hull filled, the rim, RIBS
        generators and the two generators that bound it in the projection drawn over it."""
        a, R = self.camera.screen(apex), self.camera.screen(rim)
        H = hull(np.vstack([a[None, :], R]))
        self.fill(cls, H)
        self.layers.append({"kind": "line", "class": cls + "-rim",
                            "points": rounded(np.vstack([R, R[:1]]))})
        for i in range(0, len(R), len(R) // RIBS):
            self.layers.append({"kind": "line", "class": cls + "-rib", "points": rounded([a, R[i]])})
        at = np.flatnonzero(np.all(np.isclose(H, np.round(a, 12)), axis=1))
        if at.size:
            i = int(at[0])
            sides = [H[i - 1], a, H[(i + 1) % len(H)]]
            self.layers.append({"kind": "line", "class": cls, "points": rounded(sides)})
        self.layers.append({"kind": "point", "class": cls + "-apex", "at": rounded(a)})

    def label(self, P, text, anchor="c", cls="lab", dx=0, dy=0):
        self.labels.append({"at": rounded(self.camera.screen(P)), "text": text, "anchor": anchor,
                            "class": cls, "dx": dx, "dy": dy})

    def legend(self, kind, cls, text):
        self.legend_items.append([kind, cls, text])

    def done(self, pad=0.06):
        """The figure as it is written, boxed with a margin of `pad` of its larger side."""
        pts = [p for layer in self.layers for p in (layer["points"] if "points" in layer else [layer["at"]])]
        pts += [label["at"] for label in self.labels]
        P = np.array(pts)
        lo, hi = P.min(0), P.max(0)
        m = pad * float(np.max(hi - lo))
        drawn = {layer["class"] for layer in self.layers}
        for kind, cls, _ in self.legend_items:
            if cls not in drawn:
                raise AssertionError(f"figure {self.id}: the legend names {cls}, which is not drawn")
        return {"id": self.id, "label": self.name,
                "box": [round(float(v), 4) for v in (lo[0] - m, hi[0] + m, lo[1] - m, hi[1] + m)],
                "camera": {"azimuth": self.camera.azimuth, "elevation": self.camera.elevation},
                "layers": self.layers, "labels": self.labels, "legend": self.legend_items}


def circle(sl, t, r, n=240, phi=(0.0, 2 * np.pi)):
    """The circle of constant t and r of a polar slice, in the drawing."""
    ph = np.linspace(phi[0], phi[1], n)
    return sl.to_drawing((np.full_like(ph, t), np.full_like(ph, r), ph))


def root(f, lo, hi, tol=1e-13):
    """A sign change of f on [lo, hi], by bisection; f must change sign there."""
    flo = f(lo)
    if not flo * f(hi) < 0:
        raise SystemExit(f"no sign change on [{lo}, {hi}]")
    while hi - lo > tol * max(1.0, abs(hi)):
        mid = 0.5 * (lo + hi)
        if flo * f(mid) <= 0:
            hi = mid
        else:
            lo, flo = mid, f(mid)
    return 0.5 * (lo + hi)


# ---------------------------------------------------------------- the figures

def arrow(sl, t, r, phi, size, sense=1):
    """A chevron on the circle of constant t and r at phi, pointing along +phi, or -phi."""
    ph = np.array([phi - sense * size / r, phi, phi - sense * size / r])
    tip = sl.to_drawing((np.full(3, t), np.full(3, r), ph))
    out = sl.to_drawing((t, r, phi))
    side = out / np.linalg.norm(out)
    return np.array([tip[0] + 0.5 * size * side, tip[1], tip[2] - 0.5 * size * side])


def about_axis(spec, sl, bracket, names, camera=Camera(-90, 30)):
    """Future light cones about the axis of a polar slice (t, r, phi) whose circles of
    constant t and r turn from spacelike to timelike at a critical radius r_c.

    r_c is where the published g_phiphi vanishes, found by bisection on `bracket`, and the
    circle drawn at 3 r_c/2 is checked timelike there by the published g_phiphi < 0. Cones
    stand on the axis and at four places around each of the circles r_c/2, r_c and 3 r_c/2,
    those on r_c turned by 45 degrees from the others so that no two meet, and `names` are
    the TeX names of r_c and 3 r_c/2. The floor reaches out to 2 r_c, and the cones, the axis
    and the arrows are sized against the drawn radius of r_c."""
    g_phiphi = lambda r: sl.metric((0.0, r, 0.0))[2, 2]
    critical = root(g_phiphi, *bracket)
    beyond = 1.5 * critical
    if not g_phiphi(beyond) < 0:
        raise SystemExit(f"{key(spec)}: the circle beyond the critical radius is not timelike")
    unit = float(sl.radius(critical))
    fig = Figure(spec.view, spec.label, camera)
    for k in range(12):
        ph = k * np.pi / 6
        fig.line("floor", sl.to_drawing((np.zeros(2), np.array([0.0, 2 * critical]), np.full(2, ph))))
    for r, cls in ((0.5, "floor"), (2.0, "floor"), (1.0, "critical"), (1.5, "ctc")):
        fig.line(cls, circle(sl, 0.0, r * critical), closed=True)
    for k in range(4):
        fig.line("ctc", arrow(sl, 0.0, beyond, (k + 0.75) * np.pi / 2, 0.058 * unit))
    fig.line("axis", np.array([[0, 0, -0.4], [0, 0, 0.875]]) * unit)
    cones = [(0.0, 1e-6, 0.0)]
    for r, shift in ((0.5 * critical, 0.5), (critical, 0.0), (beyond, 0.5)):
        cones += [(0.0, r, (k + shift) * np.pi / 2) for k in range(4)]
    drawn = [future_cone(sl, x, 0.187 * unit) for x in cones]
    for apex, rim in sorted(drawn, key=lambda c: camera.depth(c[0])):
        fig.cone(apex, rim)
    fig.label(np.array([0, 0, 0.875 * unit]), "$t$", "b", dy=-4)
    for r, name in ((critical, names[0]), (beyond, names[1])):
        fig.label(sl.to_drawing((0.0, r, -7 * np.pi / 18)), f"${name}$", "tl", dx=6, dy=4)
    fig.legend("cone", "cone", "future light cone")
    fig.legend("line", "critical", f"${names[0]}$, where the circle of fixed $t$ and $r$ is null")
    fig.legend("line", "ctc", f"${names[1]}$, where it is timelike")
    return fig.done(), sl


def stockum(spec):
    """Van Stockum's light cones about the axis of the dust, on the slice z = 0 at R = 1,
    drawn polar with the proper distance from the axis as its radius."""
    sl = Slice(spec.metric, spec.system, ("t", "r", "\\phi"), "polar", spec.params, spec.fixed)
    sl.proper_radius(2.5)
    return about_axis(spec, sl, (0.5, 1.5), ("r = R", "r = 3R/2"))


def godel(spec):
    """Godel's light cones about one world line of the dust, on the slice z = 0 at omega = 1,
    drawn polar with r itself as its radius: g_rr = -g_tt, so dt = +-dr runs at 45 degrees."""
    sl = Slice(spec.metric, spec.system, ("t", "r", "\\phi"), "polar", spec.params, spec.fixed)
    return about_axis(spec, sl, (0.5, 1.5), ("r = r_c", "r = 3r_c/2"))


def ergoregion(spec, camera=Camera(-90, 30)):
    """The light cones of a rotating hole on its equator, the slice theta = pi/2 of t, r and
    phi, drawn polar with r itself as its radius, down to the horizon.

    The horizon r_+ is the outer zero of the published g^rr, the ergosurface r_E the zero of
    the published g_tt outside it, both found by bisection. Cones stand at four places around
    each of the circles halfway through the ergoregion, on the ergosurface and at 3 r_E/2, those
    on r_E turned by 45 degrees from the others, oriented by the time function t, which the
    published g^tt < 0 makes one outside r_+. Before anything is drawn, g_tt is checked positive
    on the inner circle, so that its cones hold no curve of fixed r and phi, and negative on the
    outer. The floor reaches out to 7 r_E/4 and stops at r_+, where the chart does, and the
    cones and the axis are sized against r_E. The cones are narrower than van Stockum's, since
    t runs fast against proper time near the hole, so they are drawn larger. With a cone every
    45 degrees round the ergosurface no label fits beside its circles, so the legend names them."""
    sl = Slice(spec.metric, spec.system, ("t", "r", "\\phi"), "polar", spec.params, spec.fixed)
    g_tt = lambda r: sl.metric((0.0, r, 0.0))[0, 0]
    g_rr_up = lambda r: sl.inverse((0.0, r, 0.0))[1, 1]
    horizon = root(g_rr_up, 1.0, 1.9)
    ergo = root(g_tt, horizon * (1 + 1e-9), 3.0)
    inner, outer = 0.5 * (horizon + ergo), 1.5 * ergo
    if not (g_tt(inner) > 0 and g_tt(outer) < 0 and sl.inverse((0.0, inner, 0.0))[0, 0] < 0):
        raise SystemExit(f"{key(spec)}: the ergoregion is not where the published g_tt puts it")
    # Every future generator turns toward +phi inside, none turns back on r_E, and some do
    # outside, measured by dphi against the size of the generator.
    turn = {}
    for r in (inner, ergo, outer):
        _, k = generators(sl, (0.0, r, 0.0), "tau", 3600)
        turn[r] = float((k[:, 2] * r / np.linalg.norm(k * [1, 1, r], axis=1)).min())
    if not (turn[inner] > 0 and abs(turn[ergo]) < 1e-5 and turn[outer] < 0):
        raise SystemExit(f"{key(spec)}: the cones do not turn as the ergoregion says: {turn}")
    unit = ergo
    fig = Figure(spec.view, spec.label, camera)
    for k in range(12):
        ph = k * np.pi / 6
        fig.line("floor", sl.to_drawing((np.zeros(2), np.array([horizon, 1.75 * ergo]), np.full(2, ph))))
    for r, cls in ((outer, "floor"), (1.75 * ergo, "floor"), (horizon, "horizon"), (ergo, "ergo")):
        fig.line(cls, circle(sl, 0.0, r), closed=True)
    fig.line("axis", np.array([[0, 0, -0.4], [0, 0, 1.1]]) * unit)
    cones = []
    for r, shift in ((inner, 0.5), (ergo, 0.0), (outer, 0.5)):
        cones += [(0.0, r, (k + shift) * np.pi / 2) for k in range(4)]
    drawn = [future_cone(sl, x, 0.32 * unit, orient="tau") for x in cones]
    for apex, rim in sorted(drawn, key=lambda c: camera.depth(c[0])):
        fig.cone(apex, rim)
    fig.label(np.array([0, 0, 1.1 * unit]), "$t$", "b", dy=-4)
    fig.legend("cone", "cone", "future light cone")
    fig.legend("line", "horizon", "$r_+$, the horizon")
    fig.legend("line", "ergo", "$r_E$, the ergosurface, where $g_{tt} = 0$")
    return fig.done(), sl


def bubble(spec, camera=Camera(-90, 30), later=0.75):
    """A warp bubble's light cones, on the slice z = 0 of t, x and y, drawn cartesian with t up,
    at t = 0, when the declared profile centres the bubble on x = 0, with the centre's world
    line up to ct = `later`.

    The circle v_s f = 1, where the published g_tt vanishes, is found by bisection along y
    through the centre at t = 0, and the same circle about x = v_s t at ct = `later` is checked
    to lie where the published g_tt vanishes there too, which is the check that the world line
    drawn is the bubble's centre. The centre's world line x = v_s t is checked against the
    published metric to have g(u, u) = -1 for u = (1, v_s, 0), so that it is timelike and its
    proper time is t. Cones stand at the centre and at four places around each of the circles
    r_s = R/2, the circle v_s f = 1 and r_s = 2R, those on the middle circle on the axes and the
    others turned by 45 degrees from them, oriented by t, which the published g^tt = -1 makes a
    time function."""
    sl = Slice(spec.metric, spec.system, ("t", "x", "y"), "cartesian", spec.params, spec.fixed,
               spec.functions)
    g_tt = lambda t, x, y: sl.metric((t, x, y))[0, 0]
    edge = root(lambda y: g_tt(0.0, 0.0, y), 0.5, 1.5)
    v = float(sl.prep(sl.reader.parameters["v_s"].func(sl.reader.symbol["t"])))
    ph = np.linspace(0, 2 * np.pi, 13)[:-1]
    miss = max(abs(g_tt(later, v * later + edge * np.cos(a), edge * np.sin(a))) for a in ph)
    if not miss < 1e-12:
        raise SystemExit(f"{key(spec)}: the bubble is not where x = v_s t puts it, by {miss:.1e}")
    for t in np.linspace(-0.4, later, 9):
        u = np.array([1.0, v, 0.0])
        if not abs(u @ sl.metric((t, v * t, 0.0)) @ u + 1) < 1e-12:
            raise SystemExit(f"{key(spec)}: the centre's clock does not keep t at ct = {t}")

    def ring(t, r, n=240):
        a = np.linspace(0, 2 * np.pi, n)
        return np.column_stack([v * t + r * np.cos(a), r * np.sin(a), np.full(n, t)])

    fig = Figure(spec.view, spec.label, camera)
    for k in range(12):
        a = k * np.pi / 6
        fig.line("floor", np.array([[0, 0, 0], [2.5 * np.cos(a), 2.5 * np.sin(a), 0]]))
    for r in (0.5, 2.0, 2.5):
        fig.line("floor", ring(0.0, r), closed=True)
    fig.line("ergo", ring(0.0, edge), closed=True)
    fig.line("axis", np.array([[0, 0, -0.4], [0, 0, 1.2]]))
    fig.line("world", np.array([[-0.4 * v, 0, -0.4], [v * later, 0, later]]))
    cones = [(0.0, 0.0, 0.0)]
    for r, shift in ((0.5, 0.5), (edge, 0.0), (2.0, 0.5)):
        cones += [(0.0, r * np.cos((k + shift) * np.pi / 2), r * np.sin((k + shift) * np.pi / 2)) for k in range(4)]
    drawn = [future_cone(sl, x, 0.3, orient="tau") for x in cones]
    for apex, rim in sorted(drawn, key=lambda c: camera.depth(c[0])):
        fig.cone(apex, rim)
    fig.label(np.array([0, 0, 1.2]), "$t$", "b", dy=-4)
    fig.label(np.array([v * later, 0, later]), f"$x = {v:g}ct$", "bl", dx=4, dy=-2)
    fig.legend("cone", "cone", "future light cone")
    fig.legend("line", "ergo", "$v_sf = 1$, where $g_{tt} = 0$")
    fig.legend("line", "world", f"the centre of the bubble, $x = {v:g}ct$")
    return fig.done(), sl


CAPTIONS = {
    ("cosmic_string", "conical", "beam"): [
        "This is the plane $z = 0$ around the string seen from above, $t$ left out, drawn with $r$ as "
        "the radius and the angle $(1 - 4G\\mu/c^2)\\phi$, in which the plane is flat and every light "
        "ray straight. That angle runs short of a full turn by the deficit $\\delta = 8\\pi G\\mu/c^2$, "
        "so a wedge of $\\delta$ is missing from the plane. With $\\phi$ measured from the direction the "
        "beam is heading, the wedge lies behind the string, and its two edges, $\\phi = 0$ and "
        "$\\phi = 2\\pi$, are one line. A beam of parallel light arrives from the left, each ray a null "
        "geodesic of the plane, whose Christoffel symbols are $\\Gamma^r{}_{\\phi\\phi}$ and "
        "$\\Gamma^\\phi{}_{r\\phi}$.",
        "Every ray stays straight: one that reaches an edge of the wedge carries on from the same point "
        "of the other edge, in the same direction relative to it, so the rays that passed above the "
        "string cross those that passed below. In the two sectors beside the wedge, $\\delta$ across "
        "together, light from both sides arrives, and an observer there sees a source far to the left "
        "twice, at equal brightness and $\\delta$ apart.",
    ],
    ("stockum_dust", "cylindrical", "tipping"): [
        "This is the slice $z = 0$ of $t$, $r$ and $\\phi$, with $t$ up and the proper distance from "
        "the axis, $\\int e^{-r^2/2R^2}dr$, as the radius, which puts the null directions straight out "
        "from the axis at 45°. The cones stand at $t = 0$ on the axis and around the circles "
        "$r = R/2$, $R$ and $3R/2$. On the axis they are upright, and farther out the cross term "
        "$g_{t\\phi} = -r^2/R$ tips them over toward $+\\phi$, counterclockwise seen from above.",
        "At $r = R$, where $g_{\\phi\\phi}$ vanishes, one edge of every cone lies along the circle of "
        "constant $t$ and $r$, which is a closed null curve. Beyond it the cones have tipped past the "
        "horizontal, and the circle $r = 3R/2$, run counterclockwise as its arrows point, lies inside "
        "every one of them: a closed timelike curve through each of its events.",
    ],
    ("godel", "cylindrical", "tipping"): [
        "This is the slice $z = 0$ of $t$, $r$ and $\\phi$ about the axis $r = 0$, the world line of one "
        "particle of the dust, with $t$ up and $r$ as the radius, which puts the null directions "
        "$dt = \\pm dr$ at 45°. The cones stand at $t = 0$ on the axis and around the circles "
        "$r = r_c/2$, $r_c$ and $3r_c/2$, with $\\sinh r_c = 1$. On the axis they are upright, and "
        "farther out the cross term $g_{t\\phi} = -2\\sqrt{2}\\sinh^2 r/\\omega^2$ tips them over toward "
        "$+\\phi$, counterclockwise seen from above.",
        "At $r = r_c$, where $g_{\\phi\\phi}$ vanishes, one edge of every cone lies along the circle of "
        "constant $t$ and $r$, which is a closed null curve. Beyond it the cones have tipped past the "
        "horizontal, and the circle $r = 3r_c/2$, run counterclockwise as its arrows point, lies inside "
        "every one of them: a closed timelike curve through each of its events. Every world line of the "
        "dust is equivalent to every other, so the cones tip over in the same way about each one.",
    ],
    ("alcubierre", "cartesian", "bubble"): [
        "This is the slice $z = 0$ of $t$, $x$ and $y$ through a bubble moving at twice the speed of "
        "light along $x$, with $t$ up, $ct$ and $x$ drawn at one scale, at the moment $t = 0$ when the "
        "bubble is centred on $x = 0$. Far from the bubble $f = 0$ and the cones stand upright, as "
        "Minkowski's do. Inside it $f$ is close to 1, and the shift $v_sf$ tips every cone forward along "
        "$x$ so far that the vertical lies outside it: nothing inside can stay at fixed $x$. The tilt is "
        "along $x$ everywhere, since the shift points along $x$, and it depends only on the distance "
        "from the centre, as $f$ does, so the cones tip over in a ball about the centre.",
        "On the dotted circle $v_sf = 1$, where $g_{tt} = -(1 - v_s^2f^2)$ vanishes, one edge of every "
        "cone stands vertical. The centre runs along $x = 2ct$, at 63° to the vertical, outside the "
        "upright cones far from the bubble and through the middle of the cone at the centre, where "
        "$f = 1$ and the metric gives $ds^2 = -c^2dt^2$: its world line is timelike, and its clock "
        "keeps $t$.",
    ],
    ("kerr", "boyer_lindquist", "dragging"): [
        "This is the equatorial plane $\\theta = \\pi/2$ with $t$ up and $r$ and $\\phi$ as polar "
        "coordinates about the axis, for $a = 0.9\\,GM/c^2$, down to the horizon $r_+ = "
        "1.436\\,GM/c^2$, where the chart ends. Light moving in this plane stays in it, since the "
        "reflection $\\theta \\to \\pi - \\theta$ leaves it fixed. The cones stand at $t = 0$ at four "
        "places around each of three circles: $r = 3GM/c^2$, the ergosurface $r_E = 2GM/c^2$, and "
        "halfway between $r_E$ and $r_+$. On the outer circle they stand nearly upright, and closer "
        "in the cross term $g_{t\\phi} = -2GMa/c^2r$ tips them toward $+\\phi$, counterclockwise seen "
        "from above, the way the hole turns.",
        "On the ergosurface $g_{tt} = -(1 - 2GM/c^2r)$ vanishes, so $\\partial_t$ is null and one edge "
        "of every cone stands vertical, along a curve of fixed $r$ and $\\phi$. Inside it $g_{tt}$ is "
        "positive and the cones have tipped past the vertical: every future direction, timelike or "
        "null, moves toward $+\\phi$, and nothing can stay at fixed $\\phi$. Toward $r_+$ the cones "
        "also close in $r$, as $g_{rr} = r^2/\\Delta$ grows without bound where $\\Delta = r^2 - "
        "2GMr/c^2 + a^2$ falls to zero.",
    ],
    ("kerr_newman", "boyer_lindquist", "dragging"): [
        "This is the equatorial plane $\\theta = \\pi/2$ with $t$ up and $r$ and $\\phi$ as polar "
        "coordinates about the axis, for $a = 0.6\\,GM/c^2$ and $r_Q = 0.5\\,GM/c^2$, down to the "
        "horizon $r_+ = 1.625\\,GM/c^2$, where the chart ends. As in Kerr, light moving in this plane "
        "stays in it. The cones stand at $t = 0$ at four places around each of three circles: $r = "
        "3r_E/2$, the outer ergosurface $r_E = 1.866\\,GM/c^2$, and halfway between $r_E$ and $r_+$. "
        "On the outer circle they stand nearly upright, and closer in the cross term $g_{t\\phi} = "
        "-a(2GMr/c^2 - r_Q^2)/r^2$ tips them toward $+\\phi$, counterclockwise seen from above, the "
        "way the hole turns.",
        "On the ergosurface $g_{tt} = -\\left(1 - (2GMr/c^2 - r_Q^2)/r^2\\right)$ vanishes, so one "
        "edge of every cone stands vertical, along a curve of fixed $r$ and $\\phi$. Inside it the "
        "cones have tipped past the vertical, and every future direction, timelike or null, moves "
        "toward $+\\phi$. Toward $r_+$ the cones also close in $r$, as $g_{rr} = r^2/\\Delta$ grows "
        "without bound where $\\Delta = r^2 - 2GMr/c^2 + a^2 + r_Q^2$ falls to zero.",
    ],
}


FIGURES = [
    Projection("stockum_dust", "cylindrical", "tipping", "light cones about the axis", stockum,
               {"R": 1}, {"z": "0"}),
    Projection("godel", "cylindrical", "tipping", "light cones about the axis", godel,
               {"omega": 1}, {"z": "0"}),
    # The deficit the conformal diagram draws with, 4 G mu/c^2 = 0.1, so delta = 36 degrees.
    Projection("cosmic_string", "conical", "beam", "light passing the string", lambda spec: string_rays(spec),
               {"mu": "1/40", "G": 1, "delta": "pi/5"}, {"z": "0"}, fields=("christoffel",)),
    # Alcubierre's bubble at the speed and with the profile its flat view declares.
    Projection("alcubierre", "cartesian", "bubble", "the bubble in three dimensions", bubble, {}, {"z": "0"},
               input="$v_s = 2$, and Alcubierre's own profile, "
                     "$f = [\\tanh\\sigma(r_s + R) - \\tanh\\sigma(r_s - R)]/(2\\tanh\\sigma R)$ "
                     "with $R = 1$ and $\\sigma = 4$.",
               functions={"v_s": "2", "f": nr._alcubierre_profile()}),
    # The equator, where the ergoregion is widest, at the spins and charge the flat views use.
    Projection("kerr", "boyer_lindquist", "dragging", "light cones on the equator", ergoregion,
               {"G": 1, "M": 1, "a": "9/10"}, {"theta": "pi/2"}),
    Projection("kerr_newman", "boyer_lindquist", "dragging", "light cones on the equator", ergoregion,
               {"G": 1, "M": 1, "a": "3/5", "r_Q": "1/2"}, {"theta": "pi/2"}),
]


def key(spec):
    return f"{spec.metric}/{spec.system}/{spec.view}"


def draw(spec):
    """A figure as its spacetime's diagram file carries it."""
    view, sl = spec.build(spec)
    view["caption"] = CAPTIONS[(spec.metric, spec.system, spec.view)]
    view["settings"] = nr.settings(spec, sl.entry)
    view["input"] = spec.input
    fields = FIELDS + list(spec.fields)
    view["source"] = {"fields": fields, "version": build.diagram_source_version(sl.entry, fields)}
    return view


# ---------------------------------------------------------------- light passing a cosmic string

def string_rays(spec, n=12, half_width=1.8, left=-3.0, right=3.0, height=2.1):
    """A parallel beam of light passing a cosmic string, on the plane z = 0 seen from above.

    The plane's metric, dr^2 + g_phiphi dphi^2 with g_phiphi = k^2 r^2, is flat, and the angle
    psi = pi + k (phi - pi) unrolls it isometrically onto the page, with the wedge
    |psi| < pi (1 - k) around the direction phi = 0, away from the beam, missing: its two edges
    are phi = 0 and phi = 2 pi, one line. k is read from the published g_phiphi. Each ray is a
    null geodesic launched in the plane toward +X from X = `left`, integrated in the chart with
    the published Christoffel symbols, and carried across the edges when phi leaves [0, 2 pi);
    every stretch of it is checked straight on the page, which is the check that the plane is
    flat and the unrolling right. The rays are the spatial paths of light rays, t left out,
    since g_tt = -1 makes a static light ray's path a geodesic of the plane."""
    sl = Slice(spec.metric, spec.system, ("t", "r", "\\phi"), "polar", spec.params, spec.fixed)
    k = float(np.sqrt(sl.metric((0.0, 1.0, 0.0))[2, 2]))
    if not (0 < k < 1 and abs(np.sqrt(sl.metric((0.0, 2.5, 1.0))[2, 2]) / 2.5 - k) < 1e-14):
        raise SystemExit(f"{key(spec)}: g_phiphi is not k^2 r^2 with 0 < k < 1")
    _, entry, reader = nr.load(spec.metric, spec.system)
    gamma = {tuple(c["indices"]): c["value"] for c in entry["christoffel"]["variants"]["ull"]["nonzero"]}
    r_ = reader.symbol["r"]
    G_r = sp.lambdify(r_, sl.prep(reader(gamma[("r", "\\phi", "\\phi")])), "numpy")
    G_phi = sp.lambdify(r_, sl.prep(reader(gamma[("\\phi", "r", "\\phi")])), "numpy")
    wedge = np.pi * (1 - k)

    def page(r, phi):
        psi = np.pi + k * (phi - np.pi)
        return np.stack([r * np.cos(psi), r * np.sin(psi)], -1)

    def trace(b):
        X0 = np.array([left, b])
        r0, psi0 = np.hypot(*X0), np.arctan2(X0[1], X0[0]) % (2 * np.pi)
        phi0 = np.pi + (psi0 - np.pi) / k
        # Unit speed toward +X on the page: dr = cos psi, r dpsi = -sin psi, dphi = dpsi / k.
        y0 = [r0, phi0, np.cos(psi0), -np.sin(psi0) / (k * r0)]

        def rhs(_, y):
            r, phi, rd, phid = y
            return [rd, phid, -G_r(r) * phid ** 2, -2 * G_phi(r) * rd * phid]

        def leave(_, y):
            X = page(y[0], y[1])
            return min(X[0] - left + 1e-9, right - X[0], height - abs(X[1]))
        leave.terminal = True
        sol = solve_ivp(rhs, (0, 50), y0, events=leave, rtol=1e-12, atol=1e-12, method="DOP853", max_step=0.01)
        r, phi = sol.y[0], sol.y[1]
        runs, turns = [], np.floor(phi / (2 * np.pi))
        for turn in np.unique(turns):
            keep = turns == turn
            P = page(r[keep], phi[keep] - 2 * np.pi * turn)
            A, B = P[0], P[-1]
            d = (B - A) / np.linalg.norm(B - A)
            off = np.abs((P - A) @ np.array([-d[1], d[0]])).max()
            if not off < 1e-8:
                raise SystemExit(f"{key(spec)}: a ray bends on the flat page, by {off:.1e}")
            runs.append(P)
        return runs

    fig = Figure(spec.view, spec.label, Camera(-90, 90))
    reach = np.hypot(max(abs(left), right), height) * 1.05
    edge = [np.array([[0, 0], [reach * np.cos(s * wedge), reach * np.sin(s * wedge)]]) for s in (1, -1)]
    two = np.array([[0, 0], [reach * np.cos(2 * wedge), reach * np.sin(2 * wedge)],
                    [reach * np.cos(wedge), reach * np.sin(wedge)]])
    flat = lambda P: np.column_stack([P, np.zeros(len(P))])
    fig.fill("wedge", clip_box(np.array([[0, 0], edge[0][1], [reach, 0], edge[1][1]]), left, right, height))
    for sign in (1, -1):
        fig.fill("double", clip_box(two * np.array([1, sign]), left, right, height))
    for E in edge:
        fig.line("edge", flat(clip_segment(E, left, right, height)))
    for b in np.linspace(-half_width, half_width, n):
        for P in trace(b):
            fig.line("above" if b > 0 else "below", flat(P))
    fig.point("string", np.zeros(3))
    fig.label(np.array([2.4, 0.0, 0.0]), "$\\delta$", "c")
    fig.label(np.array([0.0, -0.12, 0.0]), "the string", "t", cls="small", dy=4)
    fig.legend("line", "above", "light passing above the string")
    fig.legend("line", "below", "light passing below it")
    fig.legend("fill", "wedge", f"the missing wedge, $\\delta = {round(np.degrees(2 * wedge))}°$; its two edges are one line")
    fig.legend("fill", "double", "where light from both sides arrives")
    return fig.done(pad=0.0), sl


def clip_box(P, left, right, height):
    """A convex polygon cut to the box [left, right] x [-height, height], by Sutherland-Hodgman."""
    out = np.asarray(P, dtype=float)
    for axis, bound, keep_below in ((0, right, True), (0, left, False), (1, height, True), (1, -height, False)):
        inside = (lambda p: p[axis] <= bound) if keep_below else (lambda p: p[axis] >= bound)
        new = []
        for i in range(len(out)):
            a, b = out[i - 1], out[i]
            if inside(b):
                if not inside(a):
                    new.append(a + (b - a) * (bound - a[axis]) / (b[axis] - a[axis]))
                new.append(b)
            elif inside(a):
                new.append(a + (b - a) * (bound - a[axis]) / (b[axis] - a[axis]))
        out = np.array(new)
    return out


def clip_segment(E, left, right, height):
    return clip_box(np.vstack([E, E[::-1]]), left, right, height)[:2]
