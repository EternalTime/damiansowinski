"""Causal planes, the conformal rectangle, and a maximum-radius sphere.

All lengths use ell=1, with the magnetic-source relation retained by the physics checks.
"""
import math

import numpy as np

METRIC = "hosotani_universe"
COSMIC = "cosmic_hyperspherical"
CONFORMAL = "conformal_hyperspherical"
PI = math.pi


def register_null(nr):
    for system, time, hi, ylabel in ((COSMIC, "t", 2, r"$t/\ell$"),
                                    (CONFORMAL, r"\eta", PI, r"$\eta$")):
        row = nr.Diagram(METRIC, system, "radial", "Radial light rays", (time, r"\chi"),
                         (-0.15, PI + 0.15, -0.15, hi + 0.15), r"$\chi$", ylabel,
                         {"ell": 1}, nr.EQUATOR, tau="t" if time == "t" else "eta",
                         families=nr.SIDEWAYS, areal=True, no_throat=True,
                         # dt/dchi vanishes at the singular endpoints in cosmic time.
                         step=0.0001 if system == COSMIC else 0.0025)
        key = (METRIC, system, "radial")
        nr.DIAGRAMS.append(row)
        nr.CAPTIONS[key] = [
            (r"The plane of $t$ and $\chi$ ($\theta=\pi/2$, $\phi=0$), with $t$ in units of $\ell$. "
             r"Radial null geodesics across expansion and contraction." if system == COSMIC else
             r"The plane of $\eta$ and $\chi$ ($\theta=\pi/2$, $\phi=0$), with radial null geodesics $d\chi/d\eta=\pm1$.")]
        eta_of = (lambda t: np.arccos(np.clip(1 - np.asarray(t), -1, 1))) if system == COSMIC else np.asarray
        nr.CLOSED_FORMS[key] = (lambda t, chi, f=eta_of: f(t) + chi,
                                lambda t, chi, f=eta_of: f(t) - chi, None)


def embedding(ck, src, eb):
    sl = eb.Slice(src, METRIC, COSMIC, r"\chi", r"\phi", {"t": 1, **eb.EQUATOR}, {"ell": 1})
    size = 2.0
    sphere = eb.Piece("sphere", "sheet", sl, 0.0, PI, -1.0, 1,
                      (("axis", r"$\chi=0$"), ("axis", r"$\chi=\pi$")),
                      [(PI / 2, "r", None)], size)
    ck.isometry("Hosotani: maximum-radius sphere", sphere)
    ck.radius("Hosotani: rho=ell sin chi", sphere, np.sin, size)
    ck.form("Hosotani: z=-ell cos chi", sphere, lambda chi: -np.cos(chi), size)
    surface = eb.Surface([sphere], label=r"$t=\ell$", time=1.0)
    fig = eb.figure_of([surface], {"sheet": "cover"}, size)
    fig.legend("fill", "cover", r"sphere of radius $\ell$")
    fig.legend("line", "r", r"$\chi=\pi/2$")
    fig.legend("line", "meridian", r"$\phi$ constant")
    return [eb.view("maximum", "Maximum expansion", r"$\ell$", [surface], fig.done(),
                    settings=r"$\ell=1$, $t=\ell$.")]


def register_embedding(eb):
    eb.DRAWN[METRIC] = lambda ck, src: embedding(ck, src, eb)
    eb.CAPTIONS[(METRIC, "maximum")] = [
        r"The equatorial 2-sphere at maximum expansion ($t=\ell$, $\theta=\pi/2$), with radius $\ell$."]


def register_slices(sl):
    for system, time in ((COSMIC, 1.0), (CONFORMAL, PI / 2)):
        sl.FLAT[(METRIC, system, "radial")] = lambda t=time: sl.one(
            METRIC, lambda m: sl.along(t, *m.reach(COSMIC, r"\chi")), label=r"$t=\ell$")
    sl.FLAT_METRICS.add(METRIC)


def conformal(ck, src, cd):
    def pq(eta, chi):
        return (np.asarray(eta) - chi) / 2, (np.asarray(eta) + chi) / 2

    eta_plane = cd.Plane(src, METRIC, CONFORMAL, (r"\eta", r"\chi"), cd.EQUATOR, {"ell": 1})
    cosmic_plane = cd.Plane(src, METRIC, COSMIC, ("t", r"\chi"), cd.EQUATOR, {"ell": 1})
    ck.chart("Hosotani: conformal chart", eta_plane, pq,
             ck.uniform(0.02, PI - 0.02), ck.uniform(0.02, PI - 0.02), lambda eta, chi: (1, 0))
    ck.chart("Hosotani: cosmic chart", cosmic_plane,
             lambda t, chi: pq(np.arccos(1 - np.asarray(t)), chi),
             ck.uniform(0.02, 1.98), ck.uniform(0.02, PI - 0.02), lambda t, chi: (1, 0))
    for label, t1, t2 in (("bang", 0.01, 0.001), ("crunch", 1.99, 1.999)):
        ck.diverges(f"Hosotani: {label}", cosmic_plane.kretschmann(t1, 0.7), cosmic_plane.kretschmann(t2, 0.7))
    v = cd.View("closed", "Closed magnetic universe", [-0.3, PI + 0.3, -0.3, PI + 0.3], CONFORMAL)
    rectangle = [[0, 0], [PI, 0], [PI, PI], [0, PI]]
    v.fill("region", rectangle)
    v.fill("cover", rectangle)
    for eta in (PI / 4, PI / 2, 3 * PI / 4):
        v.line("t", [[[0, eta], [PI, eta]]])
    for chi in (PI / 4, PI / 2, 3 * PI / 4):
        v.line("r", [[[chi, 0], [chi, PI]]])
    v.line("singular", [[[0, 0], [PI, 0]], [[0, PI], [PI, PI]]], zig=True)
    v.line("centre", [[[0, 0], [0, PI]], [[PI, 0], [PI, PI]]])
    v.line("null", [[[0, 0], [PI, PI]]])
    v.label_xt([PI / 2, 0], "big bang", "t", dy=8)
    v.label_xt([PI / 2, PI], "big crunch", "b", dy=-8)
    v.legend("cover", r"$0<\eta<\pi$, $0\leq\chi\leq\pi$")
    v.legend("t", r"$\eta$ constant")
    v.legend("r", r"$\chi$ constant")
    v.legend("singular", "curvature singularities")
    v.legend("centre", "regular centres")
    v.legend("null", "radial light ray")
    for moment in cd.slices.moments(METRIC):
        v.slice(moment, xt=[[[0, PI / 2], [PI, PI / 2]]])
    v.set(settings=r"$\ell=1$.")
    return [v]


def register_conformal(cd):
    cd.DRAWN[METRIC] = lambda ck, src: conformal(ck, src, cd)
    cd.CAPTIONS[(METRIC, "closed")] = [r"The closed magnetic universe, each point a 2-sphere of radius $a\sin\chi$."]


def curvature(ck, src, cv):
    scale = {"kind": "log", "from": 10, "to": 1e8}
    xaxis = cv.axis(r"$t/\ell$", 0.02, 1.98)
    K = cv.Published(src, METRIC, COSMIC, "kretschmann", ["t"], ["ell"])
    curves = cv.line_curves(ck, METRIC, "K", "K", K, "t", xaxis, scale, cv.PLOT_W, {"ell": 1})
    return [{"id": "lifetime", "label": "Expansion and contraction", "system": COSMIC,
             "kind": "line", "x": xaxis,
             "panels": [{"id": "K", "name": r"$K\ell^4$", "aspect": 1,
                         "scale": scale, "ticks": cv.value_ticks(scale)}],
             "curves": curves, "marks": [], "legend": [["line", "K", "$K$"]],
             "caption": [r"The Kretschmann scalar $K$ through expansion and contraction."],
             "settings": r"$\ell=1$.", "clip": "Logarithmic, none clipped."}]


def register_curvature(cv):
    cv.DRAWN[METRIC] = lambda ck, src: curvature(ck, src, cv)
