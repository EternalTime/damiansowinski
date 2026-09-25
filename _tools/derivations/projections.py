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
and the two generators that bound it in the projection drawn over the fill. Cones are
painted farthest first.


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
    at Euclidean distance `length` from the apex along one future null generator.

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
    apex = sl.to_drawing(x)
    return apex, apex[None, :] + length * gen / np.linalg.norm(gen, axis=1, keepdims=True)


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
        """A cone of its apex and rim, both in the drawing: the hull filled, the rim and the
        two generators that bound it in the projection drawn over it."""
        a, R = self.camera.screen(apex), self.camera.screen(rim)
        H = hull(np.vstack([a[None, :], R]))
        self.fill(cls, H)
        self.layers.append({"kind": "line", "class": cls + "-rim",
                            "points": rounded(np.vstack([R, R[:1]]))})
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


def stockum(spec, camera=Camera(-90, 30), length=0.2, rings=(4, 6, 6)):
    """Van Stockum's light cones about the axis of the dust, on the slice z = 0 at R = 1.

    The slice is (t, r, phi), drawn polar with the proper distance from the axis as its
    radius. The critical radius is where the published g_phiphi vanishes, found by
    bisection, and the circle drawn beyond it is checked timelike there by the published
    g_phiphi < 0. Cones stand on the axis and around the circles r = R/2, R and 3R/2."""
    sl = Slice(spec.metric, spec.system, ("t", "r", "\\phi"), "polar", spec.params, spec.fixed)
    sl.proper_radius(2.5)
    g_phiphi = lambda r: sl.metric((0.0, r, 0.0))[2, 2]
    critical = root(g_phiphi, 0.5, 1.5)
    beyond = 1.5 * critical
    if not g_phiphi(beyond) < 0:
        raise SystemExit("stockum: the circle beyond the critical radius is not timelike")
    fig = Figure(spec.view, spec.label, camera)
    for k in range(12):
        ph = k * np.pi / 6
        fig.line("floor", sl.to_drawing((np.zeros(2), np.array([0.0, 2.0]), np.full(2, ph))))
    for r, cls in ((0.5, "floor"), (2.0, "floor"), (1.0, "critical"), (1.5, "ctc")):
        fig.line(cls, circle(sl, 0.0, r * critical), closed=True)
    for k in range(4):
        fig.line("ctc", arrow(sl, 0.0, beyond, (k + 0.25) * np.pi / 2, 0.05))
    fig.line("axis", np.array([[0, 0, -0.35], [0, 0, 0.75]]))
    cones = [(0.0, 1e-6, 0.0)]
    for n, r, shift in zip(rings, (0.5 * critical, critical, beyond), (0.5, 0.0, 0.5)):
        cones += [(0.0, r, (k + shift) * 2 * np.pi / n) for k in range(n)]
    drawn = [future_cone(sl, x, length) for x in cones]
    for apex, rim in sorted(drawn, key=lambda c: camera.depth(c[0])):
        fig.cone(apex, rim)
    fig.label(np.array([0, 0, 0.75]), "$t$", "b", dy=-4)
    fig.label(sl.to_drawing((0.0, critical, -np.pi / 5)), "$r = R$", "tl", dx=4, dy=4)
    fig.label(sl.to_drawing((0.0, beyond, -np.pi / 5)), "$r = 3R/2$", "tl", dx=6, dy=6)
    fig.legend("cone", "cone", "future light cone")
    fig.legend("line", "critical", "$r = R$, where the circle of fixed $t$ and $r$ is null")
    fig.legend("line", "ctc", "$r = 3R/2$, where it is timelike")
    return fig.done(), sl


CAPTIONS = {
    ("cosmic_string", "conical", "beam"): [
        "This is the plane $z = 0$ around the string seen from above, $t$ left out, drawn with $r$ as "
        "the radius and the angle $(1 - 4G\\mu/c^2)\\phi$, in which the plane is flat and every light "
        "ray straight. That angle runs short of a full turn by the deficit $\\delta = 8\\pi G\\mu/c^2$, "
        "so a wedge of $\\delta$ is missing from the plane. With $\\phi$ measured from the direction the "
        "beam is heading, the wedge lies behind the string, and its two edges, $\\phi = 0$ and "
        "$\\phi = 2\\pi$, are one line. A beam of parallel light arrives from the left, each ray run "
        "with the Christoffel symbols $\\Gamma^r{}_{\\phi\\phi}$ and $\\Gamma^\\phi{}_{r\\phi}$.",
        "Every ray stays straight: one that reaches an edge of the wedge carries on from the same point "
        "of the other edge, in the same direction relative to it, so the rays that passed above the "
        "string cross those that passed below. In the two sectors beside the wedge, $\\delta$ across "
        "together, light from both sides arrives, and an observer there sees a source far to the left "
        "twice, at equal brightness and $\\delta$ apart.",
    ],
    ("stockum_dust", "cylindrical", "tipping"): [
        "This is the slice $z = 0$ of $t$, $r$ and $\\phi$, with $t$ up and the proper distance from "
        "the axis, $\\int e^{-r^2/2R^2}dr$, as the radius, so that light moving straight out runs at "
        "45° as it does on the axis. The cones stand at $t = 0$ on the axis and around the circles "
        "$r = R/2$, $R$ and $3R/2$. On the axis they are upright, and farther out the cross term "
        "$g_{t\\phi} = -r^2/R$ tips them over toward $+\\phi$, counterclockwise seen from above.",
        "At $r = R$, where $g_{\\phi\\phi}$ vanishes, one edge of every cone lies along the circle of "
        "constant $t$ and $r$, which is a closed null curve. Beyond it the cones have tipped past the "
        "horizontal, and the circle $r = 3R/2$, run counterclockwise as its arrows point, lies inside "
        "every one of them: a closed timelike curve through each of its events.",
    ],
}


FIGURES = [
    Projection("stockum_dust", "cylindrical", "tipping", "light cones about the axis", stockum,
               {"R": 1}, {"z": "0"}),
    # The deficit the conformal diagram draws with, 4 G mu/c^2 = 0.1, so delta = 36 degrees.
    Projection("cosmic_string", "conical", "beam", "light passing the string", lambda spec: string_rays(spec),
               {"mu": "1/40", "G": 1, "delta": "pi/5"}, {"z": "0"}, fields=("christoffel",)),
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
