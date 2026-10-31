"""Euclidean meridian sphere and explicitly identified constant-field continuations."""
import math

import numpy as np

import hawking_moss as hm

METRIC = hm.METRIC
NOT_DRAWN = {system: "Positive definite Euclidean sphere: no nonzero null tangent or causal boundary."
             for system in hm.EUCLIDEAN}


def register_null(nr):
    """Every Lorentzian section has its own independently checked null coordinates."""
    rows = [
        nr.Diagram(METRIC, "static_continuation", "radial", "Static constant field continuation",
                   ("t", "r"), (0, 1, -1, 1), "$r/\\ell$", "$t/\\ell$", {"ell": 1}, nr.EQUATOR,
                   mirror=True, areal=True, inside=True),
        nr.Diagram(METRIC, "wall_continuation", "radial", "Wall constant field continuation",
                   (r"\psi", r"\xi"), (0, math.pi, -math.pi/2, math.pi/2), r"$\xi/\ell$", r"$\psi$",
                   {"ell": 1}, nr.EQUATOR, tau="psi", families=nr.SIDEWAYS, inside=True),
        nr.Diagram(METRIC, "closed_continuation", "radial", "Closed constant field continuation",
                   ("t", r"\chi"), (0, math.pi, -math.pi/2, math.pi/2), r"$\chi$", r"$t/\ell$",
                   {"ell": 1}, nr.EQUATOR, families=nr.SIDEWAYS),
    ]
    captions = {
        "static_continuation": [r"The plane of $t$ and $r$ ($\theta = \pi/2$, $\phi = 0$), through the centre of the static constant field continuation ($\ell = 1$)."],
        "wall_continuation": [r"The plane of $\psi$ and $\xi$ ($\theta = \pi/2$, $\phi = 0$), of the wall constant field continuation ($\ell = 1$)."],
        "closed_continuation": [r"The plane of $t$ and $\chi$ ($\theta = \pi/2$, $\phi = 0$), of the closed constant field continuation ($\ell = 1$)."],
    }
    forms = {
        "static_continuation": (lambda t, r: t+np.arctanh(r), lambda t, r: t-np.arctanh(r), None),
        "wall_continuation": (lambda t, x: t+np.log(np.tan(x/2)), lambda t, x: t-np.log(np.tan(x/2)), None),
        "closed_continuation": (lambda t, x: np.arctan(np.sinh(t))+x,
                                lambda t, x: np.arctan(np.sinh(t))-x, None),
    }
    for row in rows:
        key = (METRIC, row.system, row.view)
        nr.DIAGRAMS.append(row)
        nr.CAPTIONS[key] = captions[row.system]
        nr.CLOSED_FORMS[key] = forms[row.system]


def embedding(ck, src):
    """Whole Euclidean meridian S2, rho=sin(xi), z=cos(xi), at ell=1."""
    import embedding as eb
    sl = eb.Slice(src, METRIC, "four_sphere", r"\xi", r"\phi", {"chi": "pi/2", "theta": "pi/2"}, {"ell": 1})
    size = 2.0
    sphere = eb.Piece("sphere", "sheet", sl, 0, math.pi, 1.0, -1,
                      (("axis", r"the pole $\xi = 0$"), ("axis", r"the pole $\xi = \pi\ell$")),
                      [(math.pi/3, "r", None), (math.pi/2, "surface", r"$\xi = \pi\ell/2$"),
                       (2*math.pi/3, "r", None)], size)
    ck.isometry("Hawking-Moss Euclidean meridian sphere", sphere)
    ck.radius("Hawking-Moss rho=sin(xi)", sphere, np.sin, size)
    ck.form("Hawking-Moss z=cos(xi)", sphere, np.cos, size)
    surface = eb.Surface([sphere])
    fig = eb.figure_of([surface], {"sheet": "cover"}, size)
    fig.legend("fill", "cover", "Euclidean meridian sphere")
    fig.legend("line", "r", r"$\xi/\ell = \pi/3,2\pi/3$")
    fig.legend("line", "surface", r"the equator $\xi = \pi\ell/2$")
    fig.legend("line", "meridian", r"$\phi$ constant")
    return [eb.view("sphere", "Euclidean meridian sphere", r"$\ell$", [surface], fig.done(), settings=r"$\ell = 1$.")]


def register_embedding(eb):
    eb.DRAWN[METRIC] = embedding
    eb.CAPTIONS[(METRIC, "sphere")] = [r"The Euclidean meridian sphere of $\xi$ and $\phi$ ($\chi = \theta = \pi/2$), including both regular poles."]


def conformal(ck, src):
    """Only the identified static constant-field continuation, not a decay geometry."""
    import conformal as cd
    pl = cd.Plane(src, METRIC, "static_continuation", ("t", "r"), cd.EQUATOR, {"ell": 1})
    pq = cd.ds_static_pq
    ck.chart("Hawking-Moss static continuation", pl, pq, ck.uniform(-8, 8), ck.uniform(0.001, 0.999),
             lambda t, r: (1, 0))
    ck.finite("Hawking-Moss static continuation regular centre", pl.kretschmann(0, 1e-6))
    ck.limit("Hawking-Moss static continuation future horizon", cd.point(*pq(30, 0.5)), [0, cd.HALF], 1e-6)
    ck.limit("Hawking-Moss static continuation past horizon", cd.point(*pq(-30, 0.5)), [0, -cd.HALF], 1e-6)
    v = cd.View("static_continuation", "Static constant field continuation", [-0.35, 1.95, -1.9, 1.9], "static_continuation")
    triangle = [[0, -cd.HALF], [cd.HALF, 0], [0, cd.HALF]]
    v.fill("region", triangle)
    v.fill("cover", triangle)
    v.line("centre", [[[0, -cd.HALF], [0, cd.HALF]]])
    v.line("horizon", [[[0, -cd.HALF], [cd.HALF, 0], [0, cd.HALF]]])
    cd.grid(v, "r", lambda r, t: pq(t, r), (0.3, 0.6, 0.9), cd.S_ALL)
    cd.grid(v, "t", pq, (-3, -1, 0, 1, 3), cd.spread(0, 1, 400, 14))
    v.label_xt([0, 0], "$r=0$", "r", dx=-5)
    v.legend("cover", "static constant field continuation")
    v.legend("centre", "regular centre")
    v.legend("horizon", r"cosmological horizons $r=\ell$")
    v.legend("r", "$r$ constant")
    v.legend("t", "$t$ constant")
    v.set(settings=r"$\ell=1$.")
    return [v]


def register_conformal(cd):
    cd.DRAWN[METRIC] = conformal
    cd.CAPTIONS[(METRIC, "static_continuation")] = [r"The static constant field continuation ($\theta=\pi/2$, $\phi=0$), bounded by its cosmological horizons."]


def curvature(ck, src):
    """Read the stored Euclidean scalars without another curvature calculation."""
    import curvature as cv
    x = cv.axis(r"$\xi/\ell$", 0, math.pi)
    params = {"ell": 1}
    panels, curves = [], []
    for field, pid, power, top in (("kretschmann", "K", 4, 30), ("ricci_scalar", "R", 2, 15)):
        stored = cv.nr_published(src, METRIC, "four_sphere", field, ["xi"], [r"\ell"])
        scale = {"kind": "linear", "from": 0, "to": top}
        panels.append({"id": pid, "name": rf"${pid}\ell^{power}$", "aspect": 2.0,
                       "scale": scale, "ticks": cv.value_ticks(scale)})
        curves.extend(cv.line_curves(ck, METRIC, pid, pid, stored, "xi", x, scale, cv.PLOT_W/2, params))
    cv.shape_checks(ck, METRIC, [2.0, 2.0])
    return [{"id": "sphere", "label": "Barrier maximum sphere", "system": "four_sphere", "kind": "line", "x": x,
             "panels": panels, "curves": curves, "marks": [],
             "legend": [["line", "K", "Kretschmann scalar"], ["line", "R", "Ricci scalar"]],
             "caption": [r"The Kretschmann scalar $K$ and the Ricci scalar $R$ of the Euclidean sphere.",
                         r"The constant Kretschmann scalar $K\ell^4=24$.",
                         r"The constant Ricci scalar $R\ell^2=12$."],
             "ricci_sentences": [[2, 0]],
             "per_scalar": cv.scalar_clips("Linear scale; no clipping.", "Linear scale; no clipping."),
             "settings": r"$\ell=1$.", "clip": "Linear scales; both curves within their panels."}]


def register_curvature(cv):
    cv.DRAWN[METRIC] = curvature
