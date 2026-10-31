"""Drawings of the uninterrupted flat false vacuum background at p=2.

All geometric input is read from the published Jordan metric.
The conformal null coordinates compactify eta<0, r>=0 to a triangle.
"""
import math

import numpy as np


METRIC = "extended_inflation"
PARAMS = {"p": 2, "t_0": 1}
SETTINGS = "$p=2$, $t_0=1$, $c=1$."
EMBEDDING_CAPTION = r"The flat equator at $t=t_0$ ($\theta=\pi/2$), with proper radius $r$."
CONFORMAL_CAPTION = r"The uninterrupted false vacuum background, each point a 2-sphere; $T\pm X=2\arctan((c\eta\pm r)/(ct_0))$."


def null_rows(nr):
    rows, captions, forms = [], {}, {}
    for system, time, box, ylabel in (
            ("cosmic", "t", (0, 2.5, 0.2, 2.7), "$t/t_0$"),
            ("conformal", r"\eta", (0, 2.5, -2.5, -0.05), r"$\eta/t_0$")):
        rows.append(nr.Diagram(METRIC, system, "radial", rf"${time}$ and $r$", (time, "r"), box,
                               "$r/(ct_0)$", ylabel, PARAMS, nr.EQUATOR, tau=nr.vm.Reader._plain(time),
                               areal=True))
        captions[(METRIC, system, "radial")] = [
            rf"The plane of ${time}$ and $r$ ($\theta=\pi/2$, $\phi=0$), with radial null geodesics."]
        forms[(METRIC, system, "radial")] = ((lambda t, r: -1/t+r, lambda t, r: -1/t-r, None)
                                                if system == "cosmic" else
                                                (lambda e, r: e+r, lambda e, r: e-r, None))
    return rows, captions, forms


def embedding(ck, src, em):
    sl = em.Slice(src, METRIC, "cosmic", "r", r"\phi", {"t": "1", "theta": "pi/2"}, PARAMS)
    size, top = 4.0, 2.0
    plane = em.Piece("plane", "sheet", sl, 0.0, top, 0.0, 1,
                     (("axis", "the centre $r=0$"), ("edge", r"the plane extends to $r\to\infty$")),
                     [(1.0, "r", None), (top, "r", None)], size)
    ck.isometry("extended inflation: the flat equator", plane)
    ck.form("extended inflation: zero embedding height", plane, lambda r: 0*r, size)
    ck.radius("extended inflation: proper radius at the reference surface", plane, lambda r: r, size)
    surface = em.Surface([plane], label="$t=t_0$", time=1.0)
    fig = em.figure_of([surface], {"sheet": "cover"}, size, meridians=12)
    fig.legend("fill", "cover", "the spatial equator")
    fig.legend("line", "r", "$r=ct_0$ and $2ct_0$")
    fig.legend("line", "meridian", r"$\phi$ constant")
    return [em.view("plane", "The spatial equator", "$ct_0$", [surface], fig.done(), settings=SETTINGS)]


def conformal(ck, src, cf):
    conf = cf.Plane(src, METRIC, "conformal", (r"\eta", "r"), cf.EQUATOR, PARAMS)
    cosmic = cf.Plane(src, METRIC, "cosmic", ("t", "r"), cf.EQUATOR, PARAMS)
    moving = lambda t, r: cf.mink_pq(-1/np.asarray(t), r)
    ck.chart("extended inflation conformal", conf, cf.mink_pq,
             ck.uniform(-20, -0.05), ck.uniform(0.001, 20), lambda e, r: (1, 0))
    ck.chart("extended inflation cosmic", cosmic, moving,
             ck.uniform(0.1, 8), ck.uniform(0.001, 20), lambda t, r: (1, 0))
    ck.diverges("extended inflation: past curvature singularity", cosmic.kretschmann(1e-3, 0.5),
                cosmic.kretschmann(1e-4, 0.5))
    ck.finite("extended inflation: regular centre at positive time",
              cosmic.kretschmann(np.ones(50), np.geomspace(1e-9, 0.1, 50)))
    pi = math.pi
    views = []
    rr = cf.spread(0, np.inf, 500, 12)
    for system in ("cosmic", "conformal"):
        v = cf.View(system, "The false vacuum background", [-0.45, pi+0.45, -pi-0.45, 0.45], system)
        tri = [[0, 0], [pi, 0], [0, -pi]]
        v.fill("region", tri)
        v.fill("cover", tri)
        for t in (0.5, 1, 2):
            v.curve("t", *cf.mink_pq(np.full_like(rr, -1/t), rr))
        for r in (0.5, 1, 2):
            v.curve("r", *cf.mink_pq(-rr, np.full_like(rr, r)))
        v.line("chartedge", [[[0, 0], [pi, 0]]])
        v.line("singular", [[[0, -pi], [pi, 0]]], zig=True)
        v.line("centre", [[[0, -pi], [0, 0]]])
        v.line("horizon", [[[0, 0], [pi/2, -pi/2]]])
        v.label_xt([pi/2, 0], "future infinity", "b", dy=-8)
        v.label_xt([pi/2, -pi/2], "$t=0$", "t", dy=12)
        v.legend("t", "$t/t_0=1/2,1,2$")
        v.legend("r", "$r/(ct_0)=1/2,1,2$")
        v.legend("centre", "$r=0$")
        v.legend("singular", "$t=0$")
        v.legend("chartedge", r"$t\to\infty$")
        v.legend("horizon", "the central observer's event horizon")
        v.legend("cover", "the false vacuum region")
        for moment in cf.slices.moments(METRIC):
            reach = moment.reach("cosmic", "r")
            r = np.linspace(*reach, 101)
            v.slice(moment, lines=[cf.mink_pq(np.full_like(r, -1.0), r)])
        v.set(settings=SETTINGS)
        views.append(v)
    return views


def curvature(ck, src, cr):
    """Sample stored scalars only, using the common scalar-curve renderer."""
    x_axis = cr.axis("$t/t_0$", 0.2, 3.0)
    panels, curves = [], []
    for field, letter, power, limits in (("kretschmann", "K", 4, (0.1, 1e6)),
                                         ("ricci_scalar", "R", 2, (1, 1e4))):
        scale = {"kind": "log", "from": limits[0], "to": limits[1]}
        panels.append({"id": letter, "name": f"${letter}(ct_0)^{power}$", "aspect": 2,
                       "scale": scale, "ticks": cr.value_ticks(scale)})
        value = cr.Published(src, METRIC, "cosmic", field, ["t"], list(PARAMS))
        curves += cr.line_curves(ck, METRIC, letter, letter, value, "t", x_axis,
                                 scale, cr.PLOT_W/2, PARAMS)
    cr.shape_checks(ck, METRIC, [p["aspect"] for p in panels])
    return [{"id": "cosmic", "label": "Cosmic time", "kind": "line", "x": x_axis,
             "panels": panels, "curves": curves, "marks": [],
             "legend": [["line", "K", "$K$"], ["line", "R", "$R$"]],
             "caption": ["The Kretschmann scalar $K$ and the Ricci scalar $R$ of the accelerating false vacuum background, "
                         "divergent at $t=0$ and decreasing with cosmic time."],
             "ricci_sentences": [], "per_scalar": cr.scalar_clips("Log scale.", "Log scale."),
             "settings": SETTINGS, "clip": "Log scales."}]
