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
same data with its own renderer. A figure of light cones, every piece of which stands in
(X, Y, T), also carries those pieces under `turn`, from which the page, with
MFS/assets/turn.js, draws it from the side a reader turns it to, by the rules its projection
was made by.

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
  turn            for a figure of light cones, its lines, cones, labels and slices in
                  (X, Y, T) at SOLID decimals, the number of ribs and the centre it turns
                  about, as _tools/README.md, "Turning a figure of light cones", defines;
  slices          the moment the spacetime's embedding diagram is cut from, where it lies
                  on the figure's floor, as a flat view carries it, its regions and lines in
                  the page's plane: the floor tinted where the embedding reaches, and a rim
                  where the embedding stops short of the floor's edge;
  caption, settings, input and source, as a flat view carries them.
"""

from dataclasses import dataclass, field
import math

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp

import build_mfs_data as build
import null_rays as nr
import slices

LENGTH = 0.32               # a cone's generators, as a fraction of the figure's slice radius
RIM = 96                    # generators per cone
RIBS = 8                    # of them drawn from the apex to the rim
NULL = 1e-12                # how far a drawn generator may miss null, against |g| |k|^2
INVERSE = 1e-12             # how far g^-1 g may miss the identity on the slice
SOLID = 6                   # decimals of the (X, Y, T) a figure that turns publishes


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
    elif not isinstance(orient, str):
        # A vector of the slice's own coordinates, for a spacetime with no arrow of time of its
        # own, as Tippett and Tsang's has none: the figure says which way it takes the future.
        e0 = J @ np.asarray(orient, dtype=float)
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


def solid(points):
    """Points of the drawing's (X, Y, T) as a figure that turns publishes them."""
    return np.round(np.asarray(points, dtype=float), SOLID).tolist()


class Figure:
    """A projected figure: layers painted in the order they are added, labels, a legend, and,
    for a figure every piece of which stands in the drawing's (X, Y, T), those pieces
    themselves, which is what drawing it from another camera takes."""

    def __init__(self, vid, label, camera):
        self.id, self.name, self.camera = vid, label, camera
        self.layers, self.labels, self.legend_items, self.slices = [], [], [], []
        self.turn = {"lines": [], "cones": [], "labels": [], "slices": []}

    def flat(self):
        """Something is drawn straight onto the page, so the figure has no other side."""
        self.turn = None

    def slice(self, moment, fills=(), lines=()):
        """A moment of the embedding diagram on the figure, painted under everything else:
        regions of the drawing's (X, Y, T), each a list of rings filled by the even odd rule,
        and the lines where the embedding stops short of the figure's edge, projected."""
        self.slices.append({**moment.json(),
                            "lines": [rounded(nr.thin(self.camera.screen(L), 0.0005)) for L in lines],
                            "points": [],
                            "fills": [[rounded(nr.thin(self.camera.screen(R), 0.0005)) for R in rings]
                                      for rings in fills]})
        if self.turn is not None:
            self.turn["slices"].append({"lines": [solid(L) for L in lines],
                                        "fills": [[solid(R) for R in rings] for rings in fills]})

    def line(self, cls, P, closed=False):
        """A polyline of the drawing's (X, Y, T), projected. Every line is painted before the
        cones, as a figure that turns paints them."""
        S = self.camera.screen(P)
        if closed:
            S = np.vstack([S, S[:1]])
        self.layers.append({"kind": "line", "class": cls, "points": rounded(nr.thin(S, 0.0005))})
        if self.turn is not None:
            if self.turn["cones"]:
                raise AssertionError(f"figure {self.id}: the line {cls} is drawn after a cone")
            P = np.asarray(P, dtype=float)
            self.turn["lines"].append({"class": cls, "points": solid(np.vstack([P, P[:1]]) if closed else P)})

    def fill(self, cls, S):
        """A polygon already in the plane of the page."""
        self.flat()
        self.layers.append({"kind": "fill", "class": cls, "points": rounded(S)})

    def point(self, cls, P):
        self.flat()
        self.layers.append({"kind": "point", "class": cls, "at": rounded(self.camera.screen(P))})

    def cone(self, apex, rim, cls="cone"):
        """A cone of its apex and rim, both in the drawing: the hull filled, the rim, RIBS
        generators and the two generators that bound it in the projection drawn over it.
        Cones are painted in the order given, which is farthest first."""
        a, R = self.camera.screen(apex), self.camera.screen(rim)
        H = hull(np.vstack([a[None, :], R]))
        self.layers.append({"kind": "fill", "class": cls, "points": rounded(H)})
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
        if self.turn is not None:
            if len(rim) % RIBS:
                raise AssertionError(f"figure {self.id}: a rim of {len(rim)} generators has no {RIBS} ribs")
            self.turn["cones"].append({"class": cls, "apex": solid(apex), "rim": solid(rim)})

    def label(self, P, text, anchor="c", cls="lab", dx=0, dy=0):
        """A label at the point P of the drawing, which it stays at from every camera."""
        self._label(P, text, anchor, cls, dx, dy, {"at": solid(P)})

    def circle_label(self, rho, t, phi, text, anchor="c", cls="lab", dx=0, dy=0):
        """A label naming the circle of radius rho about the axis at height t, at its point at
        the angle phi: from another camera it stands at the point of the circle as far round
        from the camera's azimuth, so it stays as near the reader."""
        P = np.array([rho * np.cos(phi), rho * np.sin(phi), t])
        angle = round(float(np.degrees(phi)) - self.camera.azimuth, SOLID)
        self._label(P, text, anchor, cls, dx, dy, {"circle": solid([rho, t]), "angle": angle})

    def _label(self, P, text, anchor, cls, dx, dy, place):
        self.labels.append({"at": rounded(self.camera.screen(P)), "text": text, "anchor": anchor,
                            "class": cls, "dx": dx, "dy": dy})
        if self.turn is not None:
            self.turn["labels"].append(place)

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
        box = [round(float(v), 4) for v in (lo[0] - m, hi[0] + m, lo[1] - m, hi[1] + m)]
        for mark in self.slices:
            for P in mark["lines"] + [ring for rings in mark["fills"] for ring in rings]:
                P = np.asarray(P)
                if not (P[:, 0].min() >= box[0] and P[:, 0].max() <= box[1]
                        and P[:, 1].min() >= box[2] and P[:, 1].max() <= box[3]):
                    raise AssertionError(f"figure {self.id}: the slice {mark['label']} leaves the box")
        shape = build.aspect_problems("figure", {"projections": {"": [{"id": self.id, "box": box}]}})
        if shape:
            raise AssertionError(shape[0])
        out = {"id": self.id, "label": self.name, "box": box,
               "camera": {"azimuth": self.camera.azimuth, "elevation": self.camera.elevation},
               "layers": self.layers, "labels": self.labels, "legend": self.legend_items}
        if self.slices:
            out["slices"] = self.slices
        if self.turn is not None:
            out["turn"] = {"centre": self.centre(), "ribs": RIBS, **self.turn}
        return out

    def centre(self):
        """The point the figure turns about: on the axis, halfway between the lowest and the
        highest point of everything it draws."""
        T = [p[2] for line in self.turn["lines"] for p in line["points"]]
        T += [p[2] for cone in self.turn["cones"] for p in [cone["apex"]] + cone["rim"]]
        T += [place["at"][2] if "at" in place else place["circle"][1] for place in self.turn["labels"]]
        T += [p[2] for mark in self.turn["slices"]
              for P in mark["lines"] + [ring for rings in mark["fills"] for ring in rings] for p in P]
        return [0.0, 0.0, round((min(T) + max(T)) / 2, SOLID)]


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
    # The moment t = 0 is the floor out to where the embedding stops, inside r_c, and the rim
    # stands there; beyond it no surface in flat space carries the slice.
    m = slices.moments(spec.metric)[0]
    reach = m.reach(spec.system, "r")[1]
    if not reach < critical:
        raise SystemExit(f"{key(spec)}: the embedding reaches past r_c")
    fig.slice(m, fills=[[circle(sl, 0.0, reach)]], lines=[circle(sl, 0.0, reach)])
    fig.label(np.array([0, 0, 0.875 * unit]), "$t$", "b", dy=-4)
    # Each circle is named at its own angle, so the two names, at the size the page sets them
    # at, stand clear of each other.
    for r, name, phi in ((critical, names[0], -7 * np.pi / 18), (beyond, names[1], -4 * np.pi / 18)):
        fig.circle_label(float(sl.radius(r)), 0.0, phi, f"${name}$", "tl", dx=6, dy=4)
    fig.legend("cone", "cone", "future light cone")
    fig.legend("line", "critical", f"${names[0]}$, where the circle of fixed $t$ and $r$ is null")
    fig.legend("line", "ctc", f"${names[1]}$, where it is timelike")
    return fig.done(), sl


def about_string(spec, sl, bracket, names, camera=Camera(-90, 30), sense=1, reach=None, axis="the string",
                 radius="r"):
    """Future light cones about a spinning string on a polar slice (t, r, phi), whose circles
    of constant t and r are timelike inside a critical radius r_c and spacelike outside it.

    r_c is where the published g_phiphi vanishes, found by bisection on `bracket`, and the
    circle drawn at r_c/2 is checked timelike there by the published g_phiphi < 0. Cones stand
    at four places around each of the circles r_c/2, r_c and 3 r_c/2, those on r_c turned by 45
    degrees from the others so that no two meet, and none on the string, where the metric is
    singular. `names` are the TeX names of r_c and r_c/2. The floor reaches out to 2 r_c, and
    the moment t = 0 the embedding diagram draws is the floor from r_c out. `sense` is the way
    round the timelike circle is future directed, +phi or -phi, which its arrows point; `reach`
    hands back how far the embedding runs along the slice's own radius, where that is not the
    spinning string's circumference radius; `axis` names the line r = 0 in the legend, and
    `radius` is the TeX name of the slice's radius."""
    g_phiphi = lambda r: sl.metric((0.0, r, 0.0))[2, 2]
    critical = root(g_phiphi, *bracket)
    within = 0.5 * critical
    if not g_phiphi(within) < 0:
        raise SystemExit(f"{key(spec)}: the circle inside the critical radius is not timelike")
    if not g_phiphi(1.5 * critical) > 0:
        raise SystemExit(f"{key(spec)}: the circle outside the critical radius is not spacelike")
    unit = float(sl.radius(critical))
    fig = Figure(spec.view, spec.label, camera)
    for k in range(12):
        ph = k * np.pi / 6
        fig.line("floor", sl.to_drawing((np.zeros(2), np.array([0.0, 2 * critical]), np.full(2, ph))))
    for r, cls in ((1.5, "floor"), (2.0, "floor"), (1.0, "critical"), (0.5, "ctc")):
        fig.line(cls, circle(sl, 0.0, r * critical), closed=True)
    for k in range(4):
        fig.line("ctc", arrow(sl, 0.0, within, (k + 0.75) * np.pi / 2, 0.058 * unit, sense))
    fig.line("axis", np.array([[0, 0, -0.4], [0, 0, 0.875]]) * unit)
    cones = []
    for r, shift in ((within, 0.5), (critical, 0.0), (1.5 * critical, 0.5)):
        cones += [(0.0, r, (k + shift) * np.pi / 2) for k in range(4)]
    drawn = [future_cone(sl, x, 0.187 * unit) for x in cones]
    for apex, rim in sorted(drawn, key=lambda c: camera.depth(c[0])):
        fig.cone(apex, rim)
    # The moment t = 0 is the floor from the null circle out, as far as the embedding reaches
    # or the floor does; inside r_c no surface of constant t is a moment of space.
    m = slices.moments(spec.metric)[0]
    if reach:
        lo, hi = reach(m)
    else:
        a, b = (float(nr.number(spec.params[k])) for k in ("a", "b"))
        lo, hi = (math.sqrt(R * R + a * a) / b for R in m.reach("circumference_radius", "R"))
    if not critical <= lo < critical * (1 + 1e-3) or not hi >= 2 * critical:
        raise SystemExit(f"{key(spec)}: the embedding does not run from r_c to the floor's edge")
    fig.slice(m, fills=[[circle(sl, 0.0, 2 * critical), circle(sl, 0.0, lo)]], lines=[circle(sl, 0.0, lo)])
    fig.label(np.array([0, 0, 0.875 * unit]), "$t$", "b", dy=-4)
    for r, name, phi in ((critical, names[0], -7 * np.pi / 18), (within, names[1], -4 * np.pi / 18)):
        fig.circle_label(float(sl.radius(r)), 0.0, phi, f"${name}$", "tl", dx=6, dy=4)
    fig.legend("cone", "cone", "future light cone")
    fig.legend("line", "axis", axis)
    fig.legend("line", "critical", f"${names[0]}$, where the circle of fixed $t$ and ${radius}$ is null")
    fig.legend("line", "ctc", f"${names[1]}$, where it is timelike")
    return fig.done(), sl


def spinning_string(spec):
    """The spinning string's light cones about the string, on the slice z = 0 of the proper
    radius chart at b = 0.9 and a = 0.9, drawn polar with r itself as its radius, since
    g_rr = 1 and g_tt = -1 put the null directions straight out from the string at 45 degrees."""
    sl = Slice(spec.metric, spec.system, ("t", "r", "\\phi"), "polar", spec.params, spec.fixed)
    return about_string(spec, sl, (0.5, 1.5), ("r = r_c", "r = r_c/2"))


def bonnor_rotating_dust(spec):
    """The light cones of Bonnor's dust cloud about its axis, on the slice z = 0 of the cylindrical
    chart at a = 1, drawn polar with rho itself as its radius. The circles of constant t and rho
    are timelike inside rho = a and future directed toward -phi there, since the twist enters as
    -(c dt - (a^2/rho) dphi)^2, so the arrows point clockwise; no cone stands on the axis, whose
    point rho = 0 of this slice is the singularity."""
    sl = Slice(spec.metric, spec.system, ("t", "\\rho", "\\phi"), "polar", spec.params, spec.fixed)
    return about_string(spec, sl, (0.5, 1.5), ("\\rho = a", "\\rho = a/2"), sense=-1,
                        reach=lambda m: m.reach("cylindrical", "\\rho"), axis="the axis of rotation",
                        radius="\\rho")


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


def som_raychaudhuri(spec):
    """Som and Raychaudhuri's light cones about one world line of the dust, on the slice z = 0 at
    Omega = 1, drawn polar with r itself as its radius: g_rr = -g_tt = 1, so c dt = +-dr runs at
    45 degrees."""
    sl = Slice(spec.metric, spec.system, ("t", "r", "\\phi"), "polar", spec.params, spec.fixed)
    return about_axis(spec, sl, (0.5, 1.5), ("r = r_c", "r = 3r_c/2"))


def ergoregion(spec, camera=Camera(-90, 30), horizon_between=(1.0, 1.9), ergo_below=3.0):
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
    45 degrees round the ergosurface no label fits beside its circles, so the legend names them.
    `horizon_between` brackets r_+ and `ergo_below` lies beyond r_E, in the units of the row, Kerr's
    GM/c^2 by default; Kerr-de Sitter, drawn in units of r_s with a second ergosurface next to its
    cosmological horizon, names its own."""
    sl = Slice(spec.metric, spec.system, ("t", "r", "\\phi"), "polar", spec.params, spec.fixed)
    g_tt = lambda r: sl.metric((0.0, r, 0.0))[0, 0]
    g_rr_up = lambda r: sl.inverse((0.0, r, 0.0))[1, 1]
    horizon = root(g_rr_up, *horizon_between)
    ergo = root(g_tt, horizon * (1 + 1e-9), ergo_below)
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
    m = slices.moments(spec.metric)[0]
    lo, hi = m.reach(spec.system, "r")
    if not (abs(lo - horizon) < 1e-9 and hi > 1.75 * ergo):
        raise SystemExit(f"{key(spec)}: the embedding does not run from r_+ past the floor's edge")
    fig.slice(m, fills=[[circle(sl, 0.0, 1.75 * ergo), circle(sl, 0.0, horizon)]])
    fig.label(np.array([0, 0, 1.1 * unit]), "$t$", "b", dy=-4)
    fig.legend("cone", "cone", "future light cone")
    fig.legend("line", "horizon", "$r_+$, the horizon")
    fig.legend("line", "ergo", "$r_E$, the ergosurface, where $g_{tt} = 0$")
    return fig.done(), sl


class ThroatSlice(Slice):
    """The equator of the extreme Kerr throat in Bardeen and Horowitz's global chart, (tau, y,
    phi), drawn polar with each circle of constant y at the radius 2 + arsinh(y)/sqrt 2 in units
    of r_0: the published g_yy = r_0^2/(2(1 + y^2)) on the equator, which is checked, makes
    arsinh(y)/sqrt 2 the proper distance along the throat from y = 0, and the 2 keeps the circles
    of negative y off the axis."""

    def radius(self, y):
        return 2.0 + np.arcsinh(y) / np.sqrt(2.0)

    def jacobian(self, x):
        t, y, b = x
        if not abs(self.metric(x)[1, 1] - 0.5 / (1 + y * y)) < 1e-12:
            raise SystemExit(f"{self.metric_id}/{self.system_id}: g_yy is not r_0^2/(2(1 + y^2)) on the equator")
        d, a = 1.0 / np.sqrt(2.0 * (1 + y * y)), self.radius(y)
        return np.array([[0.0, d * np.cos(b), -a * np.sin(b)],
                         [0.0, d * np.sin(b), a * np.cos(b)],
                         [1.0, 0.0, 0.0]])


def throat(spec, camera=Camera(-90, 30)):
    """The light cones of the extreme Kerr throat on its equator, the slice theta = pi/2 of tau,
    y and phi in the global chart, drawn polar with the proper distance along the throat as the
    radius, ThroatSlice.

    The circles y = +-1/sqrt 3 are the zeros of the published g_tautau, found by bisection and
    checked against 1/sqrt 3. Cones stand at four places around each of the circles y = 0,
    y = +-1/sqrt 3 and y = +-3/2, those on the dotted circles turned by 45 degrees from the
    others, oriented by the time function tau, which the published g^tautau < 0 makes one at
    every y. Before anything is drawn, g_tautau is checked negative at y = 0 and positive at
    y = +-3/2, and the cones are checked to turn as that says: some generators toward each sense
    of phi at y = 0, none toward +phi on y = 1/sqrt 3 and every one toward -phi at y = 3/2, and
    the mirror image where y is negative. The floor runs from y = -2 to 2, inside the reach of
    the embedding's cylinder, whose moment tau = 0 it is, and the circle y = 0 on it is the
    equator of the horizon the embedding's second view draws."""
    sl = ThroatSlice(spec.metric, spec.system, ("\\tau", "y", "\\phi"), "polar", spec.params, spec.fixed)
    g_tt = lambda y: sl.metric((0.0, y, 0.0))[0, 0]
    edge = root(g_tt, 0.1, 1.0)
    outer = 1.5
    if not (abs(edge - 1 / np.sqrt(3)) < 1e-10 and abs(g_tt(-edge)) < 1e-10 and g_tt(0.0) < 0
            and g_tt(outer) > 0 and g_tt(-outer) > 0
            and all(sl.inverse((0.0, y, 0.0))[0, 0] < 0 for y in np.linspace(-2, 2, 41))):
        raise SystemExit(f"{key(spec)}: d/dtau is not null where y^2 = 1/3, or tau is no time function")
    turn = {}
    for y in (-outer, -edge, 0.0, edge, outer):
        _, k = generators(sl, (0.0, y, 0.0), "tau", 3600)
        size = sl.radius(y)
        swing = k[:, 2] * size / np.linalg.norm(k * [1, 1, size], axis=1)
        turn[y] = (float(swing.min()), float(swing.max()))
    if not (turn[0.0][0] < 0 < turn[0.0][1] and abs(turn[edge][1]) < 1e-5 and turn[outer][1] < 0
            and abs(turn[-edge][0]) < 1e-5 and turn[-outer][0] > 0):
        raise SystemExit(f"{key(spec)}: the cones do not turn as g_tautau says: {turn}")
    fig = Figure(spec.view, spec.label, camera)
    reach = 2.0
    for k in range(12):
        ph = k * np.pi / 6
        fig.line("floor", sl.to_drawing((np.zeros(2), np.array([-reach, reach]), np.full(2, ph))))
    for y, cls in ((-reach, "floor"), (-outer, "floor"), (0.0, "floor"), (outer, "floor"), (reach, "floor"),
                   (-edge, "ergo"), (edge, "ergo")):
        fig.line(cls, circle(sl, 0.0, y), closed=True)
    fig.line("axis", np.array([[0, 0, -0.4], [0, 0, 1.3]]))
    cones = []
    for y, shift in ((-outer, 0.0), (-edge, 0.5), (0.0, 0.0), (edge, 0.5), (outer, 0.0)):
        cones += [(0.0, y, (k + shift) * np.pi / 2) for k in range(4)]
    drawn = [future_cone(sl, x, 0.3, orient="tau") for x in cones]
    for apex, rim in sorted(drawn, key=lambda c: camera.depth(c[0])):
        fig.cone(apex, rim)
    cylinder = slices.moments(spec.metric, "throat", label="$\\tau = 0$")[0]
    sphere = slices.moments(spec.metric, "horizon", label="$\\tau = 0$, $y = 0$")[0]
    lo, hi = cylinder.reach(spec.system, "y")
    if not (lo < -reach and hi > reach):
        raise SystemExit(f"{key(spec)}: the embedding's cylinder does not pass the floor's edges")
    fig.slice(cylinder, fills=[[circle(sl, 0.0, reach), circle(sl, 0.0, -reach)]])
    fig.slice(sphere, lines=[circle(sl, 0.0, 0.0)])
    fig.label(np.array([0, 0, 1.3]), "$\\tau$", "b", dy=-4)
    fig.legend("cone", "cone", "future light cone")
    fig.legend("line", "ergo", "$y = \\pm 1/\\sqrt{3}$, where $g_{\\tau\\tau} = 0$")
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
    m = slices.moments(spec.metric)[0]
    if not m.grid()["u"][-1] > 2.5:
        raise SystemExit(f"{key(spec)}: the embedding's plane does not pass the floor's edge")
    fig.slice(m, fills=[[ring(0.0, 2.5)]])
    fig.label(np.array([0, 0, 1.2]), "$t$", "b", dy=-4)
    fig.label(np.array([v * later, 0, later]), f"$x = {v:g}ct$", "bl", dx=4, dy=-2)
    fig.legend("cone", "cone", "future light cone")
    fig.legend("line", "ergo", "$v_sf = 1$, where $g_{tt} = 0$")
    fig.legend("line", "world", f"the centre of the bubble, $x = {v:g}ct$")
    return fig.done(), sl


def time_machine_ring(spec, camera=Camera(-62, 20), size=0.26):
    """Tippett and Tsang's bubble on the slice z = 0 of t, x and y, drawn cartesian with t up: the
    wall h = 1/2, a closed timelike curve inside it and future light cones.

    On this slice the declared bubble is the tube y^4 + (xi^2 - A^2)^2 < R^4 round the circle
    xi = A of the plane of ct and x, with xi^2 = x^2 + c^2t^2, a ring standing on its edge. Its wall
    is drawn as the four circles where it meets y = 0 and xi = A and as its sections at twelve
    angles lambda of that plane, each checked to lie where the declared h is 1/2. The circle
    xi = A, y = 0 is checked against the published metric to be timelike all the way round, with
    g(u, u) = -xi^2 (1 - 2(1 - h) sin^2(lambda)) for u = d/d lambda, within a hundredth of -xi^2
    since the declared h is 0.9993 at the centre of the box. Cones stand at eight
    places round that circle, oriented by d/d lambda, counterclockwise, at four places on the axis
    of the ring and at four outside it, oriented by t: the choice of Tippett and Tsang's figure 3,
    and each is a time direction where it is used, which generators() checks."""
    sl = Slice(spec.metric, spec.system, ("t", "x", "y"), "cartesian", spec.params, spec.fixed,
               spec.functions)
    h = sp.lambdify(sl.symbols, sl.prep(sl.reader.parameters["h"]), "numpy")
    # The sizes are read off the declared h: it is 1/2 on y = 0 at xi^2 = A^2 -+ R^2.
    inner = root(lambda x: float(h(0.0, x, 0.0)) - 0.5, 0.3, 1.0)
    outer = root(lambda x: float(h(0.0, x, 0.0)) - 0.5, 1.0, 1.5)
    A = float(np.sqrt((inner ** 2 + outer ** 2) / 2))
    R = float(np.sqrt((outer ** 2 - inner ** 2) / 2))
    if not abs(float(h(0.0, A, R)) - 0.5) < 1e-12:
        raise SystemExit(f"{key(spec)}: the declared bubble is not the box y^4 + (xi^2 - A^2)^2 = R^4")

    def ring(xi, y, n=241):
        a = np.linspace(0, 2 * np.pi, n)
        return np.column_stack([xi * np.cos(a), np.full(n, y), xi * np.sin(a)])

    def section(lam, n=161):
        b = np.linspace(0, 2 * np.pi, n)
        u, v = np.cos(b), np.sign(np.sin(b)) * np.sqrt(np.abs(np.sin(b)))      # u^2 + v^4 = 1
        xi = np.sqrt(A * A + R * R * u)
        return np.column_stack([xi * np.cos(lam), R * v, xi * np.sin(lam)])

    wall = [ring(inner, 0.0), ring(outer, 0.0), ring(A, R), ring(A, -R)] + [section(k * np.pi / 6) for k in range(12)]
    miss = max(float(np.abs(h(P[:, 2], P[:, 0], P[:, 1]) - 0.5).max()) for P in wall)
    if not miss < 1e-9:
        raise SystemExit(f"{key(spec)}: a line of the wall misses h = 1/2 by {miss:.1e}")
    for lam in np.linspace(0, 2 * np.pi, 97):
        at = (A * np.sin(lam), A * np.cos(lam), 0.0)
        u = np.array([at[1], -at[0], 0.0])              # d/d lambda = x d_ct - ct d_x
        centre = float(h(*at))
        wanted = -A * A * (1 - 2 * (1 - centre) * np.sin(lam) ** 2)
        if not (abs(u @ sl.metric(at) @ u - wanted) < 1e-12 * A * A and wanted < -0.99 * A * A):
            raise SystemExit(f"{key(spec)}: the circle xi = A is not timelike at lambda = {lam}")

    fig = Figure(spec.view, spec.label, camera)
    reach = 1.9 * A
    for y in (-reach, 0.0, reach):
        fig.line("floor", np.array([[-reach, y, 0], [reach, y, 0]]))
    for x in (-reach, 0.0, reach):
        fig.line("floor", np.array([[x, -reach, 0], [x, reach, 0]]))
    fig.line("axis", np.array([[0, 0, -1.7 * A], [0, 0, 1.7 * A]]))
    for P in wall:
        fig.line("edge", P)
    fig.line("ctc", ring(A, 0.0))
    cones = []
    for k in range(8):
        lam = (k + 0.5) * np.pi / 4
        at = (A * np.sin(lam), A * np.cos(lam), 0.0)
        cones.append(future_cone(sl, at, size, (at[1], -at[0], 0.0)))
    for t in (-1.5 * A, -0.5 * A, 0.5 * A, 1.5 * A):
        cones.append(future_cone(sl, (t, 0.0, 0.0), size, "tau"))
    for x in (-1.7 * A, 1.7 * A):
        for t in (-0.9 * A, 0.9 * A):
            cones.append(future_cone(sl, (t, x, 0.0), size, "tau"))
    for apex, rim in sorted(cones, key=lambda c: camera.depth(c[0])):
        fig.cone(apex, rim)
    # The moment the embedding diagram draws, ct = A/2 over the plane of its height.
    m = slices.moments(spec.metric, label="$ct = A/2$")[0]
    u, v = m.grid()["u"], m.grid()["v"]
    if not (u[-1] <= reach and v[-1] <= reach):
        raise SystemExit(f"{key(spec)}: the embedding's plane passes the floor's edge")
    fig.slice(m, fills=[[np.array([[u[0], v[0], slices.TIPPETT_TSANG_T], [u[-1], v[0], slices.TIPPETT_TSANG_T],
                                   [u[-1], v[-1], slices.TIPPETT_TSANG_T], [u[0], v[-1], slices.TIPPETT_TSANG_T],
                                   [u[0], v[0], slices.TIPPETT_TSANG_T]])]])
    fig.label(np.array([0, 0, 1.7 * A]), "$t$", "b", dy=-4)
    fig.label(np.array([reach, 0, 0]), "$x$", "l", dx=4)
    fig.label(np.array([0, reach, 0]), "$y$", "l", dx=4)
    fig.legend("line", "ctc", "a closed timelike curve, $x^2 + c^2t^2 = A^2$")
    fig.legend("line", "edge", "the wall of the bubble, $h = 1/2$")
    fig.legend("cone", "cone", "future light cone")
    return fig.done(), sl


CAPTIONS = {
    ("gott_time_machine", "centre_of_momentum", "loop"): [
        "The slice $z = 0$ of $t$, $x$, and $y$ in the centre of momentum frame, $t$ up, for two strings "
        "with $4G\\mu/c^2 = 1/3$, each removing a wedge of $120°$, moving at $v = 4c/5$ on the lines "
        "$y = \\pm\\ell/2$, so that $\\gamma\\sin\\alpha = 1.44$. Each string's wedge opens away from the other "
        "string, and its two faces are identified at equal times of that string's rest frame: the dashed lines "
        "are the two faces at one such time, one rising in $t$ and the other falling, and events at equal "
        "distances along them are one event.",
        "A rocket leaves the event $A$ at $t = 0$ and $x = 3\\ell/2$, crosses the upper string's wedge from the "
        "event $C$ to the same event on the other face, at an earlier $t$, and reaches the event $B$ at $x = "
        "-3\\ell/2$ at $t = 0$. It returns below the lower string in the same way, through the event $E$, and "
        "arrives at $A$ as it leaves. Each of the four stretches lies inside the future light cone at its start, "
        "so the whole is a closed timelike curve.",
    ],
    ("tippett_tsang", "cartesian", "ring"): [
        "The slice $z = 0$ of $t$, $x$, and $y$, $t$ up, with Tippett and Tsang's bubble. Its wall is a tube bent "
        "round the circle $x^2 + c^2t^2 = A^2$, a ring standing on its edge in the plane of $t$ and $x$. A cut "
        "of the ring at one $t$ is what someone outside sees at that moment: no bubble before $ct = "
        "-\\sqrt{A^2 + R^2}$, then one that splits in two, two at rest at $t = 0$, and one again before it "
        "disappears.",
        "Inside the tube the light cones turn with the angle $\\lambda$ of the plane of $t$ and $x$, so a rider "
        "who stays at the centre of the box follows the circle drawn, a closed timelike curve. We take the future "
        "counterclockwise round that circle and upward outside the bubble, as Tippett and Tsang did. On the side "
        "$x < 0$ the cones inside the tube point down in $t$ beside cones outside that point up.",
    ],
    ("lifshitz_spacetime", "poincare", "rays"): [
        "The plane $y = 0$ seen from the side, $t$ left out, with $x$ across and the depth $u$ down from the "
        "boundary along the top edge, at $z = 2$ in units of $L$. Nine light rays leave one event at $u = 2L$, "
        "$15°$ apart, each a null geodesic integrated with the Christoffel symbols of this chart.",
        "Only the ray sent straight up reaches the boundary. A ray that leaves at the angle $\\alpha$ from "
        "that one turns back at the depth $u_0 = 2L\\sin\\alpha$ and follows the catenary "
        "$u = u_0\\cosh((x - x_0)/u_0)$ through the point $x_0$ where it turns, so the boundary sees no light "
        "that left the event at an angle. At $z = 1$ the same rays are the straight lines of anti-de Sitter space, "
        "and every one of them arrives.",
    ],
    ("light_beam", "cartesian", "lens"): [
        "The plane $y = 0$ through the axis of a uniform beam of light seen from the side, $t$ left out, with $z$ "
        "across and $x$ up, in units of the beam's radius $R$ ($\\pi G\\epsilon R^2/c^4 = 1/32$). Light sent with "
        "the beam, from the left, runs straight: its $ct - z$ is constant, and every Christoffel symbol multiplies "
        "the rate of $ct - z$.",
        "Light sent against the beam, from the right, falls toward the axis. Inside the beam the pull grows in "
        "proportion to $x$, so every ray that starts there parallel to the axis reaches it when its $ct - z$ has "
        "grown by $\\sqrt{2}\\,\\pi R$, and the beam is a lens for light going the other way. "
        "A ray that starts outside the beam falls through the logarithm of the vacuum field, crosses the beam, "
        "and climbs as far out on the other side.",
    ],
    ("schrodinger_spacetime", "global", "trap"): [
        "The surface $X = 0$ of the global chart seen from the side, $V$ left out, with $cT$ across and $R$ up, "
        "in units of $\\beta$ ($L = \\beta$, $\\omega = c/\\beta$). Each curve is a null geodesic launched at "
        "$T = 0$ from the depth $R_0$ with $\\dot{R} = 0$, along the curve moving left of the plane of $T$ and "
        "$V$ there. The trap turns it back: its depth swings between $R_0$ and $\\beta^2/R_0$ as "
        "$R^2 = R_0^2\\cos^2\\omega T + (\\beta^4/R_0^2)\\sin^2\\omega T$, once in every $\\pi/\\omega$ of $T$, "
        "and the ray launched at $R_0 = \\beta$ keeps its depth.",
        "No ray reaches the boundary $R = 0$. The dashed lines are $\\omega T = \\pm\\pi/2$, the edges of the "
        "Poincaré chart, whose $t = \\tan(\\omega T)/\\omega$ and $r = R/\\cos\\omega T$ both run to infinity there: "
        "a ray leaves that chart at a finite value of its affine parameter, and in the global chart it swings on.",
    ],
    ("wormhole_time_machine", "lorentz", "trip"): [
        "The slice $Y = 0$ of $T$, $Z$, and $X$ in the Lorentz frame of the left mouth, $T$ up, with the mouths "
        "drawn as their world lines. The left mouth rests at $Z = 0$, and the right one leaves $Z = 10\\,r_0$ "
        "at $\\tau = 0$, reaches $0.905\\,c$ and $Z = 31.5\\,r_0$, and is home after $P = 54\\,r_0/c$ of its own "
        "time and $75.7\\,r_0/c$ of the left mouth's. Each dotted line joins the two mouths at one proper time "
        "$\\tau$, every $9\\,r_0/c$, and its two ends are one event.",
        "At $\\tau = 41.8\\,r_0/c$, with the right mouth on its way home, a light ray from the left mouth reaches "
        "the right one at that same $\\tau$: the closed null geodesic. For mouths small beside the distance "
        "between them, the future light cone of the event it leaves from is the Cauchy horizon, and a closed "
        "timelike curve passes through every event inside it. "
        "The one drawn leaves the left mouth at $\\tau = 60\\,r_0/c$, crosses to the right mouth at $0.46\\,c$, and "
        "steps through the throat to the event it left.",
    ],
    ("cosmic_string", "conical", "beam"): [
        "The plane $z = 0$ around the string seen from above, $t$ left out, drawn with $r$ as "
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
    ("point_particle_2plus1", "conical", "beam"): [
        "Space around the particle seen from above, $t$ left out, drawn with $r$ as the radius and the angle "
        "$\\alpha\\phi$, in which space is flat and every light ray straight. That angle runs short of a full turn "
        "by the deficit $\\delta = 2\\pi(1 - \\alpha)$, so a wedge of $\\delta$ is missing. With $\\phi$ measured "
        "from the direction the beam is heading, the wedge lies behind the particle, and its two edges, "
        "$\\phi = 0$ and $\\phi = 2\\pi$, are one line. A beam of parallel light arrives from the left, each ray a "
        "null geodesic, turned by $\\Gamma^r{}_{\\phi\\phi}$ and $\\Gamma^\\phi{}_{r\\phi}$ alone.",
        "Every ray stays straight: one that reaches an edge of the wedge carries on from the same point of the "
        "other edge, in the same direction relative to it, so the rays that passed above the particle cross "
        "those that passed below, here at a right angle. In the two sectors beside the wedge light from both "
        "sides arrives, and an observer there sees a source far to the left twice, $\\delta$ apart.",
    ],
    ("kundt_waves", "kerr_schild", "fronts"): [
        "The slice $Y = 0$ of $T$, $X$, and $Z$, with $T$ up. The envelope $x = 0$ of the wave fronts is the "
        "cone $X^2 + Z^2 = c^2T^2$, a circle growing at the speed of light, drawn at $cT = \\ell$ and $2\\ell$. "
        "A wave front is a half plane touching the cone along one null line, and on this slice a half line "
        "touching the circle of its time. Six fronts are drawn, each named by the angle $\\alpha$ between its "
        "direction of travel and the $Z$ axis, $\\alpha = \\pm 30°$, $\\pm 90°$, and $\\pm 150°$, with "
        "$u = \\tan(\\alpha/2)$.",
        "The rays are null geodesics of the wave, for any $\\ell$. On one front they are parallel, each "
        "pointing the way the front travels, and from one front to the next that direction turns with "
        "$\\alpha$. No front enters the cone, and every front ends on it.",
    ],
    ("point_particle_2plus1", "moving", "wedge"): [
        "Space and time around a particle moving along $+x$ at $v = 3c/5$ ($\\alpha = 3/4$, $\\gamma = 5/4$), "
        "$t$ up. The line element is Minkowski's, and the particle shows in the wedge it trails, drawn at three "
        "times: the two faces of the wedge are one line, and an event on one face is the event straight across "
        "from it on the other, at the same $t$.",
        "At rest the wedge is $\\delta = 2\\pi(1 - \\alpha) = 90°$ wide. In this frame lengths along $x$ are contracted and lengths "
        "along $y$ are untouched, so its half angle opens from $\\delta/2$ to "
        "$\\arctan(\\gamma\\tan(\\delta/2))$, and the wedge is $103°$ wide.",
    ],
    ("spinning_string", "proper_radius", "tipping"): [
        "The slice $z = 0$ of $t$, $r$, and $\\phi$, with $t$ up and $r$ as the radius, which puts the null "
        "directions straight out from the string at 45°. The cones stand at $t = 0$ around the circles "
        "$r = r_c/2$, $r_c$, and $3r_c/2$, with $r_c = a/b$. The cross term $g_{t\\phi} = -a$ is the same at "
        "every radius, and it tips the cones over toward $+\\phi$, counterclockwise seen from above, the more "
        "the smaller the circle they stand on.",
        "At $r = r_c$, where $g_{\\phi\\phi} = b^2r^2 - a^2$ vanishes, one edge of every cone lies along the "
        "circle of constant $t$ and $r$, which is a closed null curve. Inside it the cones have tipped past the "
        "horizontal, and the circle $r = r_c/2$, run counterclockwise as its arrows point, lies inside every "
        "one of them: a closed timelike curve through each of its events.",
    ],
    ("bonnor_rotating_dust", "cylindrical", "tipping"): [
        "The slice $z = 0$ of $t$, $\\rho$, and $\\phi$, with $ct$ up and $\\rho$ as the radius. The cones stand at "
        "$t = 0$ around the circles $\\rho = a/2$, $a$, and $3a/2$. The cross term $g_{t\\phi} = a^2/\\rho$ grows "
        "toward the centre, and it tips the cones over toward $-\\phi$, clockwise seen from above, the more the "
        "smaller the circle they stand on.",
        "At $\\rho = a$, where $g_{\\phi\\phi} = \\rho^2 - a^4/\\rho^2$ vanishes, one edge of every cone lies along "
        "the circle of constant $t$ and $\\rho$, which is a closed null curve. Inside it the cones have tipped past "
        "the horizontal, and the circle $\\rho = a/2$, run clockwise as its arrows point, lies inside every one of "
        "them: a closed timelike curve through each of its events. The centre of the slice is the singularity.",
    ],
    ("stockum_dust", "cylindrical", "tipping"): [
        "The slice $z = 0$ of $t$, $r$, and $\\phi$, with $t$ up and the proper distance from "
        "the axis, $\\int e^{-r^2/2R^2}dr$, as the radius, which puts the null directions straight out "
        "from the axis at 45°. The cones stand at $t = 0$ on the axis and around the circles $r = R/2$, "
        "$R$, and $3R/2$. On the axis they are upright, and farther out the cross term $g_{t\\phi} = "
        "-r^2/R$ tips them over toward $+\\phi$, counterclockwise seen from above.",
        "At $r = R$, where $g_{\\phi\\phi}$ vanishes, one edge of every cone lies along the circle of "
        "constant $t$ and $r$, which is a closed null curve. Beyond it the cones have tipped past the "
        "horizontal, and the circle $r = 3R/2$, run counterclockwise as its arrows point, lies inside "
        "every one of them: a closed timelike curve through each of its events.",
    ],
    ("godel", "cylindrical", "tipping"): [
        "The slice $z = 0$ of $t$, $r$, and $\\phi$ about the axis $r = 0$, the world line of one "
        "particle of the dust, with $t$ up and $r$ as the radius, which puts the null directions $dt = "
        "\\pm dr$ at 45°. The cones stand at $t = 0$ on the axis and around the circles $r = r_c/2$, "
        "$r_c$, and $3r_c/2$, with $\\sinh r_c = 1$. On the axis they are upright, and farther out the "
        "cross term $g_{t\\phi} = -2\\sqrt{2}\\sinh^2 r/\\omega^2$ tips them over toward $+\\phi$, "
        "counterclockwise seen from above.",
        "At $r = r_c$, where $g_{\\phi\\phi}$ vanishes, one edge of every cone lies along the circle of "
        "constant $t$ and $r$, which is a closed null curve. Beyond it the cones have tipped past the "
        "horizontal, and the circle $r = 3r_c/2$, run counterclockwise as its arrows point, lies inside "
        "every one of them: a closed timelike curve through each of its events. Every world line of the "
        "dust is equivalent to every other, so the cones tip over in the same way about each one.",
    ],
    ("som_raychaudhuri", "cylindrical", "tipping"): [
        "The slice $z = 0$ of $t$, $r$, and $\\phi$ about the axis $r = 0$, the world line of one "
        "particle of the dust, with $ct$ up and $r$ as the radius, which puts the null directions "
        "$c\\,dt = \\pm dr$ at 45°. The cones stand at $t = 0$ on the axis and around the circles "
        "$r = r_c/2$, $r_c$, and $3r_c/2$, with $r_c = c/\\Omega$. On the axis they are upright, and farther "
        "out the cross term $g_{t\\phi} = -\\Omega r^2/c$ tips them over toward $+\\phi$, counterclockwise "
        "seen from above.",
        "At $r = r_c$, where $g_{\\phi\\phi}$ vanishes, one edge of every cone lies along the circle of "
        "constant $t$ and $r$, which is a closed null curve. Beyond it the cones have tipped past the "
        "horizontal, and the circle $r = 3r_c/2$, run counterclockwise as its arrows point, lies inside "
        "every one of them: a closed timelike curve through each of its events. Every world line of the "
        "dust is equivalent to every other, so the cones tip over in the same way about each one.",
    ],
    ("alcubierre", "cartesian", "bubble"): [
        "The slice $z = 0$ of $t$, $x$, and $y$ through a bubble moving at twice the speed of "
        "light along $x$, with $t$ up, $ct$ and $x$ drawn at one scale, at the moment $t = 0$ when the "
        "bubble is centred on $x = 0$. Far from the bubble $f = 0$ and the cones stand upright, as "
        "Minkowski's do. Inside it $f$ is close to 1, and the shift $v_sf$ tips every cone forward along "
        "$x$ so far that the vertical lies outside it: nothing inside can stay at fixed $x$. The tilt is "
        "along $x$ everywhere, since the shift points along $x$, and it depends only on the distance from "
        "the centre, as $f$ does, so the cones tip over in a ball about the centre.",
        "On the dotted circle $v_sf = 1$, where $g_{tt} = -(1 - v_s^2f^2)$ vanishes, one edge of every "
        "cone stands vertical. The centre runs along $x = 2ct$, at 63° to the vertical, outside the "
        "upright cones far from the bubble and through the middle of the cone at the centre, where "
        "$f = 1$, $ds^2 = -c^2dt^2$: its world line is timelike, and its clock "
        "keeps $t$.",
    ],
    ("near_horizon_extreme_kerr", "global", "dragging"): [
        "The equatorial plane ($\\theta = \\pi/2$) of the global chart with $\\tau$ up and $\\phi$ the angle "
        "about the axis, each circle of constant $y$ drawn at the radius "
        "$2r_0 + r_0\\,\\mathrm{arsinh}(y)/\\sqrt{2}$, so that the distance between two circles is the proper "
        "distance along the throat. Light moving in this plane stays in it, since the reflection "
        "$\\theta \\to \\pi - \\theta$ leaves it fixed. The cones stand at $\\tau = 0$ at four places around "
        "each of five circles: $y = 0$, $y = \\pm 1/\\sqrt{3}$, and $y = \\pm 3/2$. On the middle circle they "
        "stand upright, and the cross term $g_{\\tau\\phi} = 2r_0^2y$ tips them toward $-\\phi$ where $y$ is "
        "positive and toward $+\\phi$ where it is negative.",
        "On the dotted circles $g_{\\tau\\tau} = \\tfrac{1}{2}r_0^2(3y^2 - 1)$ vanishes, so $\\partial_\\tau$ "
        "is null and one edge of every cone stands vertical. Beyond them the cones have tipped past the "
        "vertical: every future direction, timelike or null, turns about the axis, clockwise seen from above "
        "on the outer side and counterclockwise on the inner one, and nothing can stay at fixed $\\phi$. "
        "Bardeen and Horowitz likened this region to Kerr's ergosphere, and their vector "
        "$\\partial_\\tau - y\\,\\partial_\\phi$, which turns with the cones, is timelike at every point.",
    ],
    ("kerr", "boyer_lindquist", "dragging"): [
        "The equatorial plane ($\\theta = \\pi/2$) with $t$ up and $r$ and $\\phi$ as polar "
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
    ("kerr_de_sitter", "boyer_lindquist", "dragging"): [
        "The equatorial plane ($\\theta = \\pi/2$) with $t$ up and $r$ and $\\phi$ as polar coordinates about the "
        "axis, for $a = 0.45\\,r_s$ and $\\Lambda = 0.2/r_s^2$, down to the event horizon $r_+ = 0.785\\,r_s$, where "
        "the chart ends. Light moving in this plane stays in it, since the reflection $\\theta \\to \\pi - \\theta$ "
        "leaves it fixed. The cones stand at $t = 0$ at four places around each of three circles: $r = 3r_E/2$, the "
        "black hole's ergosurface $r_E = 1.105\\,r_s$, and halfway between $r_E$ and $r_+$. On the outer circle they "
        "stand nearly upright, and closer in the cross term $g_{t\\phi}$ tips them toward $+\\phi$, counterclockwise "
        "seen from above, the way the hole turns.",
        "On the ergosurface $g_{tt}$ vanishes, so $\\partial_t$ is null and one edge of every cone stands vertical, "
        "along a curve of fixed $r$ and $\\phi$. Inside it every future direction, timelike or null, moves toward "
        "$+\\phi$, and nothing can stay at fixed $\\phi$. A second ergosurface lies farther out, at $r = 3.173\\,r_s$, "
        "past the edge of the floor, and from it to the cosmological horizon $r_c = 3.232\\,r_s$ nothing can stay "
        "at fixed $\\phi$ either.",
    ],
    ("kerr_taub_nut", "boyer_lindquist", "dragging"): [
        "The equatorial plane ($\\theta = \\pi/2$) with $t$ up and $r$ and $\\phi$ as polar coordinates about the "
        "axis, for $a = m$ and $l = 5m/4$, down to the horizon $r_+ = 9m/4$, where the chart ends. The cones stand at "
        "$t = 0$ at four places around each of three circles: $r = 3r_E/2$, the ergosurface $r_E = 2.601\\,m$, and "
        "halfway between $r_E$ and $r_+$. On the outer circle they stand nearly upright, and closer in the cross term "
        "$g_{t\\phi} = a(\\Delta - r^2 - a^2 - l^2)/\\Sigma$ tips them toward $+\\phi$, counterclockwise seen from "
        "above, the way the hole turns.",
        "On the ergosurface $g_{tt} = -(\\Delta - a^2)/\\Sigma$ vanishes, so $\\partial_t$ is null and one edge of "
        "every cone stands vertical, along a curve of fixed $r$ and $\\phi$. Inside it every future direction, "
        "timelike or null, moves toward $+\\phi$, and nothing can stay at fixed $\\phi$. The twist enters the "
        "equator through $\\Sigma = r^2 + l^2$ and $\\Delta$ alone. With $l \\neq 0$ the reflection "
        "$\\theta \\to \\pi - \\theta$ is no symmetry, so of the light rays that start in this plane only some "
        "stay in it, the principal null rays among them.",
    ],
    ("kerr_newman", "boyer_lindquist", "dragging"): [
        "The equatorial plane ($\\theta = \\pi/2$) with $t$ up and $r$ and $\\phi$ as polar "
        "coordinates about the axis, for $a = 0.6\\,GM/c^2$ and $r_Q = 0.5\\,GM/c^2$, down to the "
        "horizon $r_+ = 1.624\\,GM/c^2$, where the chart ends. As in Kerr, light moving in this plane "
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
    Projection("som_raychaudhuri", "cylindrical", "tipping", "light cones about the axis", som_raychaudhuri,
               {"Omega": 1}, {"z": "0"}),
    Projection("bonnor_rotating_dust", "cylindrical", "tipping", "light cones about the axis", bonnor_rotating_dust,
               {"a": 1}, {"z": "0"}),
    # The spinning string at the values its cylinders are drawn at, r_c = a/b = 1.
    Projection("spinning_string", "proper_radius", "tipping", "light cones about the string", spinning_string,
               nr.SPINNING, {"z": "0"}),
    # The deficit the conformal diagram draws with, 4 G mu/c^2 = 0.1, so delta = 36 degrees.
    Projection("cosmic_string", "conical", "beam", "light passing the string", lambda spec: string_rays(spec),
               {"mu": "1/40", "G": 1, "delta": "pi/5"}, {"z": "0"}, fields=("christoffel",)),
    # The particle its cone is embedded at, alpha = 3/4, so delta = 90 degrees.
    Projection("point_particle_2plus1", "conical", "beam", "light passing the particle",
               lambda spec: string_rays(spec, body="particle"), {"alpha": "3/4"}, {}, fields=("christoffel",)),
    # The same particle at three fifths of the speed of light, its wedge trailing it.
    Projection("point_particle_2plus1", "moving", "wedge", "the wedge behind a moving particle",
               lambda spec: moving_wedge(spec), {"alpha": "3/4", "v": "3/5", "gamma": "5/4"}, {}),
    # Bonnor's uniform beam at the profile its flat views declare.
    Projection("light_beam", "cartesian", "lens", "light sent with the beam and against it",
               lambda spec: beam_rays(spec), {}, {"y": "0"}, input=nr.LB_ONE_INPUT, fields=("christoffel",),
               functions={"A": nr.LB_ONE}),
    # Lifshitz spacetime at z = 2, as its flat views are drawn: light from one event at u = 2L.
    Projection("lifshitz_spacetime", "poincare", "rays", "light sent toward the boundary",
               lambda spec: lifshitz_rays(spec), nr.LIFSHITZ, {"y": "0"}, fields=("christoffel",)),
    # Schrodinger spacetime's global chart at the values its flat views are drawn at, omega = c/beta.
    Projection("schrodinger_spacetime", "global", "trap", "light in the trap", lambda spec: trap_rays(spec),
               {**nr.SCHRODINGER, "omega": 1}, {"X": "0"}, fields=("christoffel",)),
    # Kundt's simplest wave on the flat space it crosses: the fronts rolled round the null cone.
    Projection("kundt_waves", "kerr_schild", "fronts", "the wave fronts round their envelope",
               lambda spec: kundt_fronts(spec), nr.KUNDT, {"Y": "0"}),
    # Gott's closed timelike curve round both strings, at the values the flat views of Grant's
    # charts are drawn at: half deficit angle pi/3, v = 4c/5 and d = l/2.
    Projection("gott_time_machine", "centre_of_momentum", "loop", "a closed timelike curve round both strings",
               lambda spec: gott_loop(spec), nr.GOTT_STRINGS, {"z": "0"}),
    # Morris, Thorne and Yurtsever's round trip, as the flat views of their chart declare it, with
    # the mouths small beside the distance between them.
    Projection("wormhole_time_machine", "lorentz", "trip", "the round trip of the right mouth",
               lambda spec: wormhole_trip(spec), {"b": 1}, {"Y": "0"},
               input="Mouths of radius $b = r_0$ that start $D = 10\\,r_0$ apart, the right one on a round trip "
                     "with rapidity $\\eta = \\tfrac{3}{2}\\sin^3(2\\pi\\tau/P)$ and $P = 54\\,r_0/c$."),
    # Tippett and Tsang's bubble as its flat views declare it, in units of A.
    Projection("tippett_tsang", "cartesian", "ring", "the bubble in three dimensions", time_machine_ring, {},
               {"z": "0"}, input=nr.TT_INPUT, functions={"h": nr.TT_H}),
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
    # Kerr-de Sitter at the values of its flat views, in units of r_s: r_+ = 0.785 and the black
    # hole's ergosurface r_E = 1.105, with the cosmological horizon's ergosurface, at 3.173, beyond
    # the floor.
    Projection("kerr_de_sitter", "boyer_lindquist", "dragging", "$\\Lambda > 0$, light cones on the equator",
               lambda spec: ergoregion(spec, horizon_between=(0.6, 1.0), ergo_below=2.0), nr.KDS, {"theta": "pi/2"}),
    # Kerr-Taub-NUT at the values of its flat views, in units of m: r_+ = 9/4 and r_E = 2.601.
    Projection("kerr_taub_nut", "boyer_lindquist", "dragging", "light cones on the equator",
               lambda spec: ergoregion(spec, horizon_between=(2.0, 2.39), ergo_below=4.0), nr.KTN, {"theta": "pi/2"}),
    # The throat of extreme Kerr on its equator, at r_0 = 1, as its flat views are drawn.
    Projection("near_horizon_extreme_kerr", "global", "dragging", "light cones on the equator", throat,
               {"r_0": 1}, {"theta": "pi/2"}),
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

def string_rays(spec, n=12, half_width=1.8, left=-3.0, right=3.0, height=2.1, body="string"):
    """A parallel beam of light passing a cosmic string, on the plane z = 0 seen from above, or
    passing a point particle in three dimensions, whose space is that plane, as `body` names it.

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
    # The ideal string's moment is the plane itself, the embedding's reference cone and the
    # sheet outside Gott's core together, out to where the embedding stops: on the page the
    # disc of that radius less the wedge the deficit removes, cut to the figure.
    # The cone unrolled is the same moment, the ideal string's cone alone.
    for m in slices.moments(spec.metric):
        reach = m.reach("conical", "r", reference=True)
        if not reach[0] == 0:
            raise SystemExit(f"{key(spec)}: the embedding's cone does not reach the {body}")
        a = np.linspace(wedge, 2 * np.pi - wedge, 721)
        arc = np.column_stack([reach[1] * np.cos(a), reach[1] * np.sin(a)])
        disc = clip_box(np.vstack([[0.0, 0.0], arc]), left, right, height)
        rims = [flat(run) for run in slices.clip_runs(arc, (left, -height), (right, height))]
        fig.slice(m, fills=[[flat(disc)]], lines=rims)
    fig.label(np.array([2.4, 0.0, 0.0]), "$\\delta$", "c")
    fig.label(np.array([0.0, -0.12, 0.0]), f"the {body}", "t", cls="small", dy=4)
    fig.legend("line", "above", f"light passing above the {body}")
    fig.legend("line", "below", "light passing below it")
    fig.legend("fill", "wedge", f"the missing wedge, $\\delta = {round(np.degrees(2 * wedge))}°$; its two edges are one line")
    fig.legend("fill", "double", "where light from both sides arrives")
    return fig.done(pad=0.0), sl


def gott_loop(spec, x0=1.5, camera=Camera(-65, 24)):
    """Gott's closed timelike curve round both strings, on the slice z = 0 of the centre of
    momentum chart, drawn cartesian with t up.

    The upper string moves along +x on y = d and the lower along -x on y = -d. In the rest frame
    of the upper string, x1 = gamma (x - beta ct), ct1 = gamma (ct - beta x), its wedge is
    |x1| < (y - d) tan(alpha) and the faces are identified at equal ct1, a rotation by 2 alpha about
    the string. A = (0, x0, 0) and B = (0, -x0, 0) are in that frame (-+gamma beta x0, +-gamma x0, 0),
    and the straight path from one to the other round the string meets each face at the foot of the
    perpendicular from the event, at the distance s = gamma x0 sin(alpha) - d cos(alpha) along the
    face and at ct1 = 0: the events C and D, which are one event. The return from B to A round the
    lower string is the same path turned by pi about the t axis. Every stretch is checked timelike
    and future directed against the published metric, C and D are checked to lie on the faces at
    one rest time and to be carried onto one another by the rotation, and the curve is checked to
    close."""
    sl = Slice(spec.metric, spec.system, ("t", "x", "y"), "cartesian", spec.params, spec.fixed)
    g = sl.metric((0.0, 0.0, 0.0))
    if not (np.array_equal(g, np.diag([-1.0, 1.0, 1.0])) and np.array_equal(sl.metric((0.7, -1.3, 2.1)), g)):
        raise SystemExit(f"{key(spec)}: the published metric on the slice is not Minkowski's")
    value = {k: float(nr.number(v)) for k, v in spec.params.items()}
    alpha, beta, d, gamma = value["alpha"], value["v"], value["d"], value["gamma"]
    if not (abs(gamma - 1 / np.sqrt(1 - beta ** 2)) < 1e-14 and abs(alpha - 4 * np.pi * value["G"] * value["mu"]) < 1e-14):
        raise SystemExit(f"{key(spec)}: gamma or alpha is not the one v and mu fix")
    if not gamma * np.sin(alpha) > 1:
        raise SystemExit(f"{key(spec)}: gamma sin(alpha) <= 1, so the strings make no closed timelike curve")

    def lab(ct1, x1, y):
        """An event of the upper string's rest frame in the chart's (ct, x, y)."""
        return np.array([gamma * (ct1 + beta * x1), gamma * (x1 + beta * ct1), y])

    def rest(e):
        return np.array([gamma * (e[0] - beta * e[1]), gamma * (e[1] - beta * e[0]), e[2]])

    w = gamma * x0
    s_ = w * np.sin(alpha) - d * np.cos(alpha)
    if not s_ > 0:
        raise SystemExit(f"{key(spec)}: the path round the string misses its wedge")
    A, B = np.array([0.0, x0, 0.0]), np.array([0.0, -x0, 0.0])
    C = lab(0.0, s_ * np.sin(alpha), d + s_ * np.cos(alpha))
    D = lab(0.0, -s_ * np.sin(alpha), d + s_ * np.cos(alpha))
    turn = np.diag([1.0, -1.0, -1.0])
    E, F = turn @ C, turn @ D
    # C and D on the two faces at one rest time, and one the rotation of the other about the string.
    c1, d1 = rest(C), rest(D)
    rot = np.array([[np.cos(2 * alpha), -np.sin(2 * alpha)], [np.sin(2 * alpha), np.cos(2 * alpha)]])
    faces = [abs(abs(e[1]) - (e[2] - d) * np.tan(alpha)) for e in (c1, d1)]
    carried = rot @ (c1[1:] - [0, d]) + [0, d]
    if not (max(faces) < 1e-12 and abs(c1[0] - d1[0]) < 1e-12 and np.abs(carried - d1[1:]).max() < 1e-12):
        raise SystemExit(f"{key(spec)}: C and D are not identified events of the wedge's faces")
    stretches = [(A, C), (D, B), (B, E), (F, A)]
    speeds = []
    for a_, b_ in stretches:
        k = b_ - a_
        if not (k @ g @ k < 0 and k[0] > 0):
            raise SystemExit(f"{key(spec)}: a stretch of the curve is not timelike and future directed")
        speeds.append(float(np.hypot(k[1], k[2]) / k[0]))
    if not (np.array_equal(stretches[-1][1], stretches[0][0]) and np.array_equal(stretches[1][1], stretches[2][0])):
        raise SystemExit(f"{key(spec)}: the curve does not close")

    draw = lambda e: sl.to_drawing((e[0], e[1], e[2]))
    fig = Figure(spec.view, spec.label, camera)
    top = float(C[0]) * 1.12
    wide, deep = float(C[1]) * 1.08, float(C[2]) * 1.15
    fig.line("floor", np.array([[-wide, -deep, 0], [wide, -deep, 0], [wide, deep, 0], [-wide, deep, 0]]), closed=True)
    fig.line("floor", np.array([[-wide, 0, 0], [wide, 0, 0]]))
    fig.line("floor", np.array([[0, -deep, 0], [0, deep, 0]]))
    for sign in (1, -1):
        fig.line("world", np.array([[-sign * beta * top, sign * d, -top], [sign * beta * top, sign * d, top]]))
    reach = 1.25 * s_
    for sign in (1, -1):
        arm = draw(lab(0.0, sign * reach * np.sin(alpha), d + reach * np.cos(alpha)))
        for flip in (np.eye(3), np.diag([-1.0, -1.0, 1.0])):
            fig.line("edge", np.array([flip @ np.array([0.0, d, 0.0]), flip @ arm]))
    for a_, b_ in ((C, D), (E, F)):
        fig.line("axis", np.array([draw(a_), draw(b_)]))
    for a_, b_ in stretches:
        fig.line("ctc", np.array([draw(a_), draw(b_)]))
    cones = [future_cone(sl, tuple(e), 0.5, "tau") for e in (A, D, B, F)]
    for apex, rim in sorted(cones, key=lambda c: camera.depth(c[0])):
        fig.cone(apex, rim)
    for e, text, anchor, dx, dy in ((A, "$A$", "tl", 6, 4), (B, "$B$", "tr", -6, 4), (C, "$C$", "b", 0, -6),
                                    (D, "$C$", "t", 0, 6), (E, "$E$", "b", 0, -6), (F, "$E$", "t", 0, 6)):
        fig.label(draw(e), text, anchor, dx=dx, dy=dy)
    fig.label(np.array([0.0, 0.0, top]), "$t$", "b", dy=-4)
    fig.legend("line", "world", "the two strings")
    fig.legend("line", "edge", "the two faces of each wedge at one time of its string's rest frame")
    fig.legend("line", "axis", "from an event on one face to the same event on the other")
    fig.legend("line", "ctc", f"a closed timelike curve, at ${max(speeds):.2f}\\,c$")
    fig.legend("cone", "cone", "future light cone")
    return fig.done(), sl


def moving_wedge(spec, top=2.0, reach=2.4, camera=Camera(-65, 24)):
    """A point particle in three dimensions moving along +x, in the frame where it moves, drawn
    cartesian with t up: its world line, and the wedge it trails at three times of that frame.

    In the particle's rest frame, x1 = gamma (x - beta ct), ct1 = gamma (ct - beta x), the wedge is
    |y| < -x1 tan(delta/2) with delta = 2 pi (1 - alpha), about the negative x1 axis, and its faces are
    identified at equal ct1 by the rotation through delta about the particle. Two identified events
    share ct1 and x1, so they share t and x as well and differ by the sign of y: in the frame where
    the particle moves the faces are identified at equal t, across the wedge, and the wedge's half
    angle there is arctan(gamma tan(delta/2)), wider than at rest. Every one of those statements is
    checked on the events drawn, against the published metric, which is checked to be Minkowski's."""
    sl = Slice(spec.metric, spec.system, ("t", "x", "y"), "cartesian", spec.params, spec.fixed)
    g = sl.metric((0.0, 0.0, 0.0))
    if not (np.array_equal(g, np.diag([-1.0, 1.0, 1.0])) and np.array_equal(sl.metric((0.7, -1.3, 2.1)), g)):
        raise SystemExit(f"{key(spec)}: the published metric on the slice is not Minkowski's")
    value = {k: float(nr.number(v)) for k, v in spec.params.items()}
    alpha, beta, gamma = value["alpha"], value["v"], value["gamma"]
    if not abs(gamma - 1 / np.sqrt(1 - beta ** 2)) < 1e-14:
        raise SystemExit(f"{key(spec)}: gamma is not the one v fixes")
    half = np.pi * (1 - alpha)
    if not 0 < half < np.pi / 2:
        raise SystemExit(f"{key(spec)}: the wedge is no narrower than a half plane")
    wide = float(np.arctan(gamma * np.tan(half)))

    def lab(ct1, x1, y):
        """An event of the particle's rest frame in the chart's (ct, x, y)."""
        return np.array([gamma * (ct1 + beta * x1), gamma * (x1 + beta * ct1), y])

    def rest(e):
        return np.array([gamma * (e[0] - beta * e[1]), gamma * (e[1] - beta * e[0]), e[2]])

    def face(t, s, sign):
        """The event of the face y = sign |y| at the chart's time t and the rest distance s from the particle."""
        x1 = -s * np.cos(half)
        return lab(t / gamma - beta * x1, x1, sign * s * np.sin(half))

    rot = np.array([[np.cos(2 * half), -np.sin(2 * half)], [np.sin(2 * half), np.cos(2 * half)]])
    times = (-top / 2, 0.0, top / 2)
    for t in times:
        for s_ in (0.5, reach):
            up, down = face(t, s_, 1), face(t, s_, -1)
            u1, d1 = rest(up), rest(down)
            carried = rot @ u1[1:]
            angle = np.arctan2(up[2], beta * t - up[1])
            if not (abs(up[0] - t) < 1e-12 and abs(down[0] - t) < 1e-12 and abs(up[1] - down[1]) < 1e-12
                    and abs(u1[0] - d1[0]) < 1e-12 and np.abs(carried - d1[1:]).max() < 1e-12
                    and abs(angle - wide) < 1e-12):
                raise SystemExit(f"{key(spec)}: the faces drawn are not the identified faces of the wedge")
    k = np.array([1.0, beta, 0.0])
    if not k @ g @ k < 0:
        raise SystemExit(f"{key(spec)}: the particle's world line is not timelike")

    draw = lambda e: sl.to_drawing((e[0], e[1], e[2]))
    fig = Figure(spec.view, spec.label, camera)
    far = max(abs(draw(face(t, reach, 1))[0]) for t in times) * 1.05
    deep = abs(draw(face(0.0, reach, 1))[1]) * 1.1
    fig.line("floor", np.array([[-far, -deep, 0], [far, -deep, 0], [far, deep, 0], [-far, deep, 0]]), closed=True)
    fig.line("floor", np.array([[-far, 0, 0], [far, 0, 0]]))
    fig.line("floor", np.array([[0, -deep, 0], [0, deep, 0]]))
    fig.line("world", np.array([draw(np.array([-top, -beta * top, 0.0])), draw(np.array([top, beta * top, 0.0]))]))
    for t in times:
        here = np.array([t, beta * t, 0.0])
        for sign in (1, -1):
            fig.line("edge", np.array([draw(here), draw(face(t, reach, sign))]))
        for s_ in (reach / 2, reach):
            fig.line("axis", np.array([draw(face(t, s_, 1)), draw(face(t, s_, -1))]))
    ahead = [np.array([t, beta * t + 1.0, 0.0]) for t in times[:2]]
    cones = [future_cone(sl, tuple(e), 0.5, "tau") for e in ahead]
    for apex, rim in sorted(cones, key=lambda c: camera.depth(c[0])):
        fig.cone(apex, rim)
    fig.label(draw(np.array([top, beta * top, 0.0])), "the particle", "b", cls="small", dy=-4)
    fig.label(np.array([0.0, 0.0, top]), "$t$", "b", dy=-4)
    fig.legend("line", "world", f"the particle, moving along $+x$ at ${beta:.1f}\\,c$")
    fig.legend("line", "edge", f"the two faces of its wedge at three times, ${2 * np.degrees(wide):.0f}°$ apart")
    fig.legend("line", "axis", "from an event on one face to the same event on the other, at one $t$")
    fig.legend("cone", "cone", "future light cone")
    return fig.done(), sl


def kundt_fronts(spec, top=2.0, reach=1.6, camera=Camera(-62, 26)):
    """The wave fronts of the simplest Kundt wave on the slice Y = 0 of the flat space it crosses,
    drawn cartesian with T up: the envelope x = 0, the null cone X^2 + Z^2 = c^2T^2, as its circles
    at three times, and six fronts, each the half plane cT = X sin(alpha) + Z cos(alpha), x >= 0,
    which is u = tan(alpha/2).

    A point of a front is cT (sin, cos) + x (cos, -sin) in (X, Z), with x the name the chart
    defines, and its rays are the lines of constant x along (cT, X, Z) = (1, sin, cos). Every
    point drawn is checked against the chart's own x and u, read from its definitions, every ray
    to be null by the published metric, and every front to touch the cone on the generator
    where its x vanishes."""
    # The reader holds the chart's three names as functions, and its `held` is each one's
    # definition; u's and v's hold x, so each is written out twice over, and handed to the slice
    # as the published values' own functions.
    _, _, reader = nr.load(spec.metric, spec.system)
    written = {n: sp.sympify(reader.parameters[n]).subs(reader.held).subs(reader.held).subs(reader.c, 1)
               for n in ("x", "u", "v")}
    sl = Slice(spec.metric, spec.system, ("T", "X", "Z"), "cartesian", spec.params, spec.fixed,
               functions={n: str(e) for n, e in written.items()})
    names = [sp.lambdify(sl.symbols, written[n], "numpy") for n in ("x", "u")]
    angles = np.radians([-150.0, -90.0, -30.0, 30.0, 90.0, 150.0])
    times = (top / 2, top)
    depths = (reach / 2, reach)

    def event(alpha, t, s):
        """(cT, X, Z) on the front alpha at the time t and the distance s from the envelope."""
        return np.array([t, t * np.sin(alpha) + s * np.cos(alpha), t * np.cos(alpha) - s * np.sin(alpha)])

    for alpha in angles:
        k = np.array([1.0, np.sin(alpha), np.cos(alpha)])
        for t in (0.3, *times):
            edge = event(alpha, t, 0.0)
            if abs(edge[1] ** 2 + edge[2] ** 2 - edge[0] ** 2) > 1e-12:
                raise SystemExit(f"{key(spec)}: a front does not touch the cone")
            for s_ in depths:
                e = event(alpha, t, s_)
                x_here, u_here = (float(f(*e)) for f in names)
                if not (abs(x_here - s_) < 1e-10 and abs(u_here - np.tan(alpha / 2)) < 1e-9):
                    raise SystemExit(f"{key(spec)}: the event drawn is not on the front u = tan(alpha/2) at x = {s_}")
                if abs(k @ sl.metric(tuple(e)) @ k) > 1e-9:
                    raise SystemExit(f"{key(spec)}: a ray drawn is not null by the published metric")

    draw = lambda e: sl.to_drawing((e[0], e[1], e[2]))
    fig = Figure(spec.view, spec.label, camera)
    far = top + reach
    fig.line("floor", np.array([[-far, -far, 0], [far, -far, 0], [far, far, 0], [-far, far, 0]]), closed=True)
    fig.line("floor", np.array([[-far, 0, 0], [far, 0, 0]]))
    fig.line("floor", np.array([[0, -far, 0], [0, far, 0]]))
    fig.line("floor", np.array([[0, 0, 0], [0, 0, top]]))
    turn = np.linspace(0.0, 2 * np.pi, 181)
    for t in (top / 2, top):
        fig.line("critical", np.array([draw(np.array([t, t * np.sin(a), t * np.cos(a)])) for a in turn]))
    for alpha in angles:
        fig.line("critical", np.array([draw(event(alpha, 0.0, 0.0)), draw(event(alpha, top, 0.0))]))
        for s_ in depths:
            fig.line("above", np.array([draw(event(alpha, 0.0, s_)), draw(event(alpha, top, s_))]))
    for alpha in angles:
        for t in times:
            fig.line("world", np.array([draw(event(alpha, t, 0.0)), draw(event(alpha, t, reach))]))
    fig.label(np.array([0.0, 0.0, top]), "$T$", "b", dy=-4)
    fig.label(np.array([far, 0.0, 0.0]), "$X$", "l", dx=4)
    fig.label(np.array([0.0, far, 0.0]), "$Z$", "l", dx=4)
    fig.legend("line", "critical", "the envelope $x = 0$, a circle growing at the speed of light, and the null "
                                   "line where each front touches it")
    fig.legend("line", "world", "six wave fronts at $cT = \\ell$ and $2\\ell$, each a half line touching the circle")
    fig.legend("line", "above", "rays, two on each front")
    return fig.done(), sl


def wormhole_trip(spec, loop=60.0, camera=Camera(-72, 20)):
    """Morris, Thorne and Yurtsever's round trip on the slice Y = 0 of the Lorentz chart, drawn
    cartesian with T up and Z across: the two mouths' world lines, the lines joining them at equal
    proper times, the closed null geodesic, the Cauchy horizon as the future light cone of the
    event that geodesic leaves from, and a closed timelike curve.

    The mouths are drawn as their world lines, their radius small beside the distance between
    them. nr.WormholeTrip is the right mouth's world line. Every step of both world lines is
    checked timelike against the published metric and of the proper time it is marked with, the
    closed null geodesic is checked null, the pairs of identified events are checked spacelike
    separated before it and timelike after, and the closed timelike curve, from the left mouth
    at tau = `loop` to the right mouth at the same tau, is checked timelike and future directed."""
    sl = Slice(spec.metric, spec.system, ("T", "Z", "X"), "cartesian", spec.params, spec.fixed)
    g = sl.metric((0.0, 0.0, 0.0))
    if not (np.array_equal(g, np.diag([-1.0, 1.0, 1.0])) and np.array_equal(sl.metric((3.1, -1.3, 2.1)), g)):
        raise SystemExit(f"{key(spec)}: the published metric on the slice is not Minkowski's")
    trip = nr.WormholeTrip()
    P, D = trip.period, trip.distance

    def event(mouth, tau):
        T, Z = mouth(tau)
        return np.array([float(T), float(Z), 0.0])

    # Both world lines, each step timelike and of the proper time it spans.
    tau = np.linspace(-4.0, P + 10.0, 3401)
    for mouth in (trip.left, trip.right):
        T, Z = mouth(tau)
        k = np.stack([np.diff(T), np.diff(Z), np.zeros(len(tau) - 1)], 1)
        elapsed = np.sqrt(-np.einsum("na,ab,nb->n", k, g, k))
        if not (np.all(k[:, 0] > 0) and np.abs(elapsed / np.diff(tau) - 1).max() < 1e-4):
            raise SystemExit(f"{key(spec)}: a mouth's world line is not marked with its proper time")
    L_c, R_c = event(trip.left, trip.tau_c), event(trip.right, trip.tau_c)
    k = R_c - L_c
    if not (abs(k @ g @ k) < 1e-9 * k[0] ** 2 and k[0] > 0):
        raise SystemExit(f"{key(spec)}: the closed geodesic is not null")
    marks = np.arange(0.0, P + 10.0, 9.0)
    for m in marks:
        k = event(trip.right, m) - event(trip.left, m)
        if (k @ g @ k < 0) != (m > trip.tau_c):
            raise SystemExit(f"{key(spec)}: the mouths at tau = {m} are on the wrong side of the horizon")
    A, B = event(trip.left, loop), event(trip.right, loop)
    k = B - A
    if not (k @ g @ k < 0 and k[0] > 0):
        raise SystemExit(f"{key(spec)}: the closed curve is not timelike and future directed")
    speed = float(abs(k[1]) / k[0])

    draw = lambda e: sl.to_drawing((e[0], e[1], e[2]))
    fig = Figure(spec.view, spec.label, camera)
    far = float(trip.right(np.linspace(0, P, 2001))[1].max())
    low, top = float(tau[0]), float(B[0]) + 3.0
    fig.line("floor", np.array([[-8.0, -8.0, 0], [far + 4.0, -8.0, 0], [far + 4.0, 8.0, 0], [-8.0, 8.0, 0]]), closed=True)
    fig.line("floor", np.array([[-8.0, 0, 0], [far + 4.0, 0, 0]]))
    fig.line("world", np.array([[0.0, 0.0, low], [0.0, 0.0, top]]))
    T, Z = trip.right(np.linspace(low, loop + 3.0, 1201))
    fig.line("world", np.stack([Z, np.zeros_like(Z), T], 1))
    for m in marks:
        fig.line("axis", np.array([draw(event(trip.left, m)), draw(event(trip.right, m))]))
    # The horizon: the future light cone of the event the closed null geodesic leaves from, its
    # generators checked null by generators(), out to the height the figure is drawn to.
    height = top - L_c[0]
    apex, rim = future_cone(sl, tuple(L_c), height * np.sqrt(2.0), "tau")
    fig.line("horizon", rim, closed=True)
    for i in range(0, len(rim), len(rim) // 12):
        fig.line("horizon", np.array([apex, rim[i]]))
    fig.line("ergo", np.array([draw(L_c), draw(R_c)]))
    fig.line("ctc", np.array([draw(A), draw(B)]))
    cones = [future_cone(sl, tuple(e), 3.0, "tau") for e in (A, L_c)]
    for apex_, rim_ in sorted(cones, key=lambda c: camera.depth(c[0])):
        fig.cone(apex_, rim_)
    for m in marks[::2]:
        fig.label(draw(event(trip.left, m)), f"${m:g}$", "r", "small", dx=-6)
    fig.label(draw(event(trip.left, marks[0])) + np.array([0, 0, -2.5]), "$c\\tau/r_0$", "r", "small", dx=-6)
    fig.label((draw(L_c) + draw(R_c)) / 2, "$\\mathcal{C}$", "tl", dx=4, dy=6)
    fig.label(np.array([0.0, 0.0, top]), "$T$", "b", dy=-4)
    fig.label(np.array([far + 4.0, 0.0, 0.0]), "$Z$", "l", dx=4)
    fig.legend("line", "world", "the two mouths, the left at rest and the right on its round trip")
    fig.legend("line", "axis", "from one mouth to the other at one proper time $\\tau$, the two ends one event")
    fig.legend("line", "ergo", "the closed null geodesic $\\mathcal{C}$")
    fig.legend("line", "horizon", "the Cauchy horizon, the future light cone of the event $\\mathcal{C}$ leaves from")
    fig.legend("line", "ctc", f"a closed timelike curve, at ${speed:.2f}\\,c$ from left mouth to right")
    fig.legend("cone", "cone", "future light cone")
    return fig.done(), sl


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


# ---------------------------------------------------------------- light sent against a beam of light

def beam_rays(spec, left=-5.0, right=5.0, height=3.0):
    """Light sent along Bonnor's uniform beam and against it, on the plane y = 0 through the beam's
    axis, seen from the side with t left out: z across the page and x up it.

    The beam's reflection in y keeps every geodesic launched in the plane in it. Each ray is a null
    geodesic of the Cartesian chart, integrated in t, x and z with the published Christoffel symbols
    and the declared profile. The rays sent with the beam leave x = b at z = `left` with
    dz = c dt, and each is checked to keep its x. The rays sent against it leave x = b at
    z = `right` with dx = 0, along the other null direction of the plane of t and z,
    c dt = -(1 - A) dz/(1 + A), and each is checked null against the published metric all the way,
    to keep ct - z's rate, which is Bonnor's (8.14), and to keep x'^2 + A (ct' - z')^2, his (8.16);
    one that starts inside the beam and has not left it is checked against b cos(u/2R), with
    sqrt(2) u = ct - z measured from its start, the pendulum his interior solution is."""
    sl = Slice(spec.metric, spec.system, ("t", "x", "z"), "cartesian", spec.params, spec.fixed, spec.functions)
    _, entry, reader = nr.load(spec.metric, spec.system)
    names = entry["coords"]
    keep = [names.index(c) for c in ("t", "x", "z")]
    x_ = reader.symbol["x"]
    gamma = {}
    for c in entry["christoffel"]["variants"]["ull"]["nonzero"]:
        ix = tuple(names.index(n) for n in c["indices"])
        value = sl.prep(reader(c["value"]))
        if value == 0:
            continue
        if not all(i in keep for i in ix):
            raise SystemExit(f"{key(spec)}: Gamma^{c['indices'][0]}_{c['indices'][1]}{c['indices'][2]} "
                             "does not vanish on the plane, so a ray launched in it leaves it")
        if value.free_symbols - {x_}:
            raise SystemExit(f"{key(spec)}: a Christoffel symbol depends on more than x on the plane")
        gamma[tuple(keep.index(i) for i in ix)] = sp.lambdify(x_, value, "numpy")
    profile = sp.lambdify(x_, sl.prep(reader.parameters["A"]), "numpy")

    def rhs(_, w):
        k = w[3:]
        acc = np.zeros(3)
        for (a, b, c), f in gamma.items():
            acc[a] -= float(f(w[1])) * k[b] * k[c]
        return np.concatenate([k, acc])

    def trace(b, against):
        A = float(profile(b))
        if against:
            start = [0.0, b, right, (1 - A) / (1 + A), 0.0, -1.0]
        else:
            start = [0.0, b, left, 1.0, 0.0, 1.0]

        def leave(_, w):
            return min(w[2] - left + 1e-9, right - w[2] + 1e-9, height - abs(w[1]) + 1e-9)
        leave.terminal = True
        sol = solve_ivp(rhs, (0, 60), start, events=leave, rtol=1e-12, atol=1e-12, method="DOP853", max_step=0.01)
        t, x, z, kt, kx, kz = sol.y
        null = max(abs(float(k @ sl.metric((0.0, xi, 0.0)) @ k)) for xi, k in zip(x[::20], sol.y[3:].T[::20]))
        if not null < 1e-9:
            raise SystemExit(f"{key(spec)}: a ray misses null by {null:.1e}")
        if not against:
            if not np.abs(x - b).max() < 1e-12:
                raise SystemExit(f"{key(spec)}: a ray sent with the beam is deflected")
            return np.column_stack([z, x])
        rate = kt - kz
        energy = kx ** 2 + profile(x) * rate ** 2
        if not (np.ptp(rate) < 1e-9 and np.ptp(energy) < 1e-9):
            raise SystemExit(f"{key(spec)}: a ray sent against the beam loses Bonnor's first integrals, "
                             f"by {np.ptp(rate):.1e} and {np.ptp(energy):.1e}")
        if np.abs(x).max() <= 1:
            u = ((t - z) - (t[0] - z[0])) / np.sqrt(2)
            miss = np.abs(x - b * np.cos(u / 2)).max()
            if not miss < 1e-8:
                raise SystemExit(f"{key(spec)}: a ray inside the beam misses b cos(u/2R) by {miss:.1e}")
        return np.column_stack([z, x])

    fig = Figure(spec.view, spec.label, Camera(-90, 90))
    flat = lambda P: np.column_stack([P, np.zeros(len(P))])
    fig.fill("double", np.array([[left, -1.0], [right, -1.0], [right, 1.0], [left, 1.0]]))
    for edge in (-1.0, 1.0):
        fig.line("edge", flat(np.array([[left, edge], [right, edge]])))
    fig.line("axis", flat(np.array([[left, 0.0], [right, 0.0]])))
    for b in (-2.5, -1.5, -0.5, 0.5, 1.5, 2.5):
        fig.line("below", flat(trace(b, False)))
    for b in np.linspace(-2.75, 2.75, 12):
        fig.line("above", flat(trace(float(b), True)))
    fig.label(np.array([left + 0.1, 1.0, 0.0]), "$x = R$", "bl", cls="small", dy=-4)
    fig.label(np.array([left + 0.1, -1.0, 0.0]), "$x = -R$", "tl", cls="small", dy=4)
    fig.legend("line", "above", "light sent against the beam, from the right, each ray starting parallel to the axis")
    fig.legend("line", "below", "light sent with the beam, from the left")
    fig.legend("fill", "double", "the beam, shining toward $+z$, to the right")
    fig.legend("line", "edge", "the edge of the beam")
    fig.legend("line", "axis", "the axis of the beam")
    return fig.done(pad=0.0), sl


# ---------------------------------------------------------------- light sent toward Lifshitz spacetime's boundary

LIFSHITZ_SOURCE = 2.0           # the depth u of the event the rays leave, in L
LIFSHITZ_ANGLES = (15, 30, 45, 60)      # degrees from the straight way to the boundary, on either side


def lifshitz_rays(spec, half_width=4.0, depth=4.5):
    """Light rays from one event of Lifshitz spacetime at z = 2, on the plane y = 0 of the inverse
    radius chart seen from the side with t left out: x across the page and the depth u down it,
    the boundary u = 0 along the top.

    The reflection of y keeps every geodesic launched in the plane in it, which is checked on the
    published Christoffel symbols. Each ray is a null geodesic of the chart, integrated in t, x
    and u with those symbols from the event x = 0, u = 2L, leaving at an angle alpha from the
    straight way up, measured in the orthonormal frame of an observer at rest there. The affine
    parameter lambda is exchanged for sigma with d(lambda) = d(sigma)/u^2, since the boundary lies
    at an infinite affine distance. Each ray is checked null against the published metric all the way,
    to keep its energy -g_tt k^t and its momentum g_xx k^x, and against its closed form, which the
    drawing never uses: the catenary u = u_0 cosh((x - x_0)/u_0) with u_0 = 2L sin(alpha), Keeler,
    Knodel and Liu's turning point, where E^2 (u/L)^{2z - 2} = p^2. The same rays at z = 1, from the
    same published metric, are the straight lines x = (2L - u) tan(alpha) of anti-de Sitter space,
    which is conformally flat, checked too, and every one of them arrives."""
    sl = Slice(spec.metric, spec.system, ("t", "x", "u"), "cartesian", spec.params, spec.fixed)
    ads = Slice(spec.metric, spec.system, ("t", "x", "u"), "cartesian", {**spec.params, "z": 1}, spec.fixed)
    _, entry, reader = nr.load(spec.metric, spec.system)
    names = entry["coords"]
    keep = [names.index(c) for c in ("t", "x", "u")]
    u_ = reader.symbol["u"]

    def symbols_of(slice_):
        gamma = {}
        for c in entry["christoffel"]["variants"]["ull"]["nonzero"]:
            ix = tuple(names.index(n) for n in c["indices"])
            value = slice_.prep(reader(c["value"]))
            if value == 0 or not all(i in keep for i in ix[1:]):
                continue
            if ix[0] not in keep:
                raise SystemExit(f"{key(spec)}: Gamma^{c['indices'][0]}_{c['indices'][1]}{c['indices'][2]} "
                                 "does not vanish on the plane, so a ray launched in it leaves it")
            if value.free_symbols - {u_}:
                raise SystemExit(f"{key(spec)}: a Christoffel symbol depends on more than u on the plane")
            gamma[tuple(keep.index(i) for i in ix)] = sp.lambdify(u_, value, "numpy")
        return gamma

    def trace(slice_, gamma, alpha):
        g = slice_.metric((0.0, 0.0, LIFSHITZ_SOURCE))
        # A null vector at the source: unit energy in the static frame there, the spatial part at
        # alpha from the direction of decreasing u.
        start = [0.0, 0.0, LIFSHITZ_SOURCE, 1 / np.sqrt(-g[0, 0]), np.sin(alpha) / np.sqrt(g[1, 1]),
                 -np.cos(alpha) / np.sqrt(g[2, 2])]

        def rhs(_, w):
            k = w[3:]
            acc = np.zeros(3)
            for (a, b, c), f in gamma.items():
                acc[a] -= float(f(w[2])) * k[b] * k[c]
            return np.concatenate([k, acc]) / w[2] ** 2

        def leave(_, w):
            return min(w[2] - 1e-3, depth - w[2] + 1e-9, half_width - abs(w[1]) + 1e-9)
        leave.terminal = True
        sol = solve_ivp(rhs, (0, 400), start, events=leave, rtol=1e-12, atol=1e-13, method="DOP853", max_step=0.01)
        if sol.status != 1:
            raise SystemExit(f"{key(spec)}: a ray neither arrives nor leaves the drawing")
        t, x, u, kt, kx, ku = sol.y
        metrics = np.array([slice_.metric((0.0, 0.0, ui)) for ui in u])
        k = sol.y[3:].T
        null = np.abs(np.einsum("ni,nij,nj->n", k, metrics, k)) / np.abs(metrics[:, 0, 0] * kt ** 2)
        energy, momentum = -metrics[:, 0, 0] * kt, metrics[:, 1, 1] * kx
        # Beside the boundary the two terms of the norm are each without bound, so it is judged from u = L/20 on.
        null = null[u >= 0.05]
        if not (np.abs(null).max() < 1e-8 and np.ptp(energy) < 1e-8 * energy[0] and np.ptp(momentum) < 1e-8):
            raise SystemExit(f"{key(spec)}: a ray loses its null norm, its energy or its momentum, by "
                             f"{np.abs(null).max():.1e}, {np.ptp(energy):.1e} and {np.ptp(momentum):.1e}")
        return x, u

    gamma, gamma_ads = symbols_of(sl), symbols_of(ads)
    fig = Figure(spec.view, spec.label, Camera(-90, 90))
    fig.flat()
    page = lambda x, u: np.column_stack([x, -np.asarray(u), np.zeros(len(x))])
    fig.line("edge", page(np.array([-half_width, half_width]), np.zeros(2)))
    arrived = 0
    for degrees in LIFSHITZ_ANGLES:
        for side in (-1, 1):
            alpha = side * np.radians(degrees)
            x, u = trace(ads, gamma_ads, alpha)
            miss = np.abs(x - (LIFSHITZ_SOURCE - u) * np.tan(alpha)).max()
            if not (miss < 1e-7 and u[-1] < 2e-3):
                raise SystemExit(f"{key(spec)}: a ray at z = 1 misses its straight line by {miss:.1e} or the boundary")
            arrived += 1
            fig.line("below", page(x, u))
    for degrees in (0,) + LIFSHITZ_ANGLES:
        for side in ((1,) if degrees == 0 else (-1, 1)):
            alpha = side * np.radians(degrees)
            x, u = trace(sl, gamma, alpha)
            if degrees == 0:
                if not (np.abs(x).max() < 1e-12 and u[-1] < 2e-3):
                    raise SystemExit(f"{key(spec)}: the ray sent straight up does not reach the boundary")
            else:
                u0 = LIFSHITZ_SOURCE * abs(np.sin(alpha))
                x0 = side * u0 * np.arccosh(LIFSHITZ_SOURCE / u0)
                miss = np.abs(u - u0 * np.cosh((x - x0) / u0)).max()
                if not (miss < 1e-7 and abs(u.min() - u0) < 1e-6 and u[-1] > LIFSHITZ_SOURCE):
                    raise SystemExit(f"{key(spec)}: a ray misses its catenary by {miss:.1e} or does not turn back at "
                                     f"u = {u0:.3f}")
            fig.line("above", page(x, u))
    if arrived != 2 * len(LIFSHITZ_ANGLES):
        raise SystemExit(f"{key(spec)}: not every ray at z = 1 arrives")
    fig.point("string", np.array([0.0, -LIFSHITZ_SOURCE, 0.0]))
    fig.label(np.array([-half_width + 0.1, 0.0, 0.0]), "$u = 0$", "tl", cls="small", dy=4)
    fig.label(np.array([0.0, -LIFSHITZ_SOURCE, 0.0]), "$u = 2L$", "t", cls="small", dy=8)
    fig.legend("line", "above", "light at $z = 2$, leaving one event every $15°$")
    fig.legend("line", "below", "the same rays at $z = 1$, anti-de Sitter space")
    fig.legend("line", "edge", "the boundary, $u = 0$")
    fig.legend("point", "string", "the event the rays leave, at $x = 0$ and $u = 2L$")
    return fig.done(pad=0.03), sl


# ---------------------------------------------------------------- light in Schrodinger spacetime's trap

def trap_rays(spec, depths=(1.0, 0.8, 0.6, 0.45, 0.35, 0.3), span=math.pi):
    """Null geodesics of Schrodinger spacetime's global chart on the surface X = 0, seen from the
    side with V left out: cT across the page and R up it, at L = beta = 1 and omega = c/beta.

    The reflection of X keeps every geodesic launched in the surface in it. Each ray is integrated
    in T, V and R with the published Christoffel symbols, from T = 0 at the depth R_0 with no
    velocity along R, along the null direction of the plane of T and V that is not d/dV,
    dV = -(1/R^2 + R^2) dT/2, both ways until |T| = `span`. It is checked null against the
    published metric all the way, to keep its momentum along V, g_TV dT/dlambda, and to follow
    R^2 = R_0^2 cos^2 T + sin^2 T/R_0^2, which is the Poincare chart's r^2 = R_0^2 + t^2/R_0^2
    carried along t = tan T and r = R/cos T; the ray at R_0 = 1 is checked to keep its depth."""
    sl = Slice(spec.metric, spec.system, ("T", "V", "R"), "cartesian", spec.params, spec.fixed, spec.functions)
    _, entry, reader = nr.load(spec.metric, spec.system)
    names = entry["coords"]
    keep = [names.index(c) for c in ("T", "V", "R")]
    R_ = reader.symbol["R"]
    gamma = {}
    for c in entry["christoffel"]["variants"]["ull"]["nonzero"]:
        ix = tuple(names.index(n) for n in c["indices"])
        value = sl.prep(reader(c["value"]))
        if value == 0:
            continue
        if not all(i in keep for i in ix[1:]):
            # A symbol with X below multiplies the velocity along X, which a ray in the surface has none of.
            continue
        if ix[0] not in keep:
            raise SystemExit(f"{key(spec)}: Gamma^{c['indices'][0]}_{c['indices'][1]}{c['indices'][2]} "
                             "does not vanish on the surface, so a ray launched in it leaves it")
        if value.free_symbols - {R_}:
            raise SystemExit(f"{key(spec)}: a Christoffel symbol depends on more than R on the surface")
        gamma[tuple(keep.index(i) for i in ix)] = sp.lambdify(R_, value, "numpy")

    def rhs(_, w):
        k = w[3:]
        acc = np.zeros(3)
        for (a, b, c), f in gamma.items():
            acc[a] -= float(f(w[2])) * k[b] * k[c]
        return np.concatenate([k, acc])

    def trace(R0, sense):
        start = [0.0, 0.0, R0, sense, -sense * (1 / R0 ** 2 + R0 ** 2) / 2, 0.0]

        def leave(_, w):
            return span - abs(w[0])
        leave.terminal = True
        sol = solve_ivp(rhs, (0, 400), start, events=leave, rtol=1e-12, atol=1e-12, method="DOP853", max_step=0.01)
        T, V, R, kT, kV, kR = sol.y
        if not abs(abs(T[-1]) - span) < 1e-9:
            raise SystemExit(f"{key(spec)}: the ray from R = {R0} stops at T = {T[-1]:.3f}")
        null = max(abs(float(k @ sl.metric((0.0, 0.0, r)) @ k)) for r, k in zip(R[::20], sol.y[3:].T[::20]))
        if not null < 1e-9:
            raise SystemExit(f"{key(spec)}: a ray misses null by {null:.1e}")
        momentum = kT / R ** 2
        if not np.ptp(momentum) < 1e-9:
            raise SystemExit(f"{key(spec)}: a ray loses its momentum along V, by {np.ptp(momentum):.1e}")
        miss = np.abs(R ** 2 - (R0 ** 2 * np.cos(T) ** 2 + np.sin(T) ** 2 / R0 ** 2)).max()
        if not miss < 1e-8:
            raise SystemExit(f"{key(spec)}: the ray from R = {R0} misses its closed form by {miss:.1e}")
        return np.column_stack([T, R])

    fig = Figure(spec.view, spec.label, Camera(-90, 90))
    fig.flat()
    flat = lambda P: np.column_stack([P, np.zeros(len(P))])
    top = 1 / min(depths) + 0.15
    fig.line("axis", flat(np.array([[-span, 0.0], [span, 0.0]])))
    for edge in (-math.pi / 2, math.pi / 2):
        fig.line("edge", flat(np.array([[edge, 0.0], [edge, top]])))
    for R0 in depths:
        back, on = trace(R0, -1), trace(R0, 1)
        fig.line("below" if R0 == 1.0 else "above", flat(np.vstack([back[::-1], on[1:]])))
    fig.label(np.array([-math.pi / 2, top, 0.0]), "$\\omega T = -\\pi/2$", "bc", cls="small", dy=-4)
    fig.label(np.array([math.pi / 2, top, 0.0]), "$\\omega T = \\pi/2$", "bc", cls="small", dy=-4)
    fig.label(np.array([span, 0.0, 0.0]), "$R = 0$", "br", cls="small", dy=-4)
    fig.legend("line", "above", "light launched at $T = 0$ from $R_0 = 0.8$, $0.6$, $0.45$, $0.35$ and $0.3\\,\\beta$")
    fig.legend("line", "below", "light launched from $R_0 = \\beta$, the bottom of the trap, which keeps its depth")
    fig.legend("line", "edge", "the edges of the Poincaré chart, $\\omega T = \\pm\\pi/2$")
    fig.legend("line", "axis", "the boundary $R = 0$")
    return fig.done(pad=0.02), sl
