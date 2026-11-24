"""Drawings of the uninterrupted flat false vacuum background at p=2.

All geometric input is read from the published Jordan metric.
The conformal null coordinates compactify eta<0, r>=0 to a triangle.

The embedding diagram draws the equator theta = pi/2 of the cosmic chart from t = t_0/2 to
3t_0/2.  Every moment is flat, a^2(dr^2 + r^2 dphi^2), so a moment drawn alone is Euclid's plane,
and what the metric adds is the scale factor's growth: a ring of free particles at rest on
r = c t_0, which is the circle of proper radius (t/t_0)^2 c t_0, and the event horizon of the
observer at the centre, the circle of comoving radius c t_0^2/t, the line r = -c eta of the
conformal diagram, whose proper radius is c t.  The ring is inside the horizon before t = t_0,
where a = 1 and both have the proper radius c t_0, and outside it after.  It is drawn in the form
of the Big Rip's: the ring's world tube, and beside it the flat moments as a movie.
extended_inflation.md beside this file records each choice.
"""
import math

import numpy as np


METRIC = "extended_inflation"
PARAMS = {"p": 2, "t_0": 1}
SETTINGS = "$p=2$, $t_0=1$, $c=1$."

# The embedding: the key moments t/t_0 = 1/2, 1 and 3/2, symmetric about t_0, where a = 1, and
# within the spacetime diagram's run of t; between them a frame every t_0/40.
NAMED = (0.5, 1.0, 1.5)
STEP = 0.025
RING = 1.0          # the ring's comoving radius, c t_0, the conformal diagram's line r = c t_0
TOP = 2.5           # the disc's proper radius, past the ring's widest, (3/2)^2 c t_0 = 2.25 c t_0
SIZE = 2 * TOP
LIFT = 3.0          # the height of t_0 up the stack, in c t_0
EMBEDDING_CAPTIONS = {
    (METRIC, "tube"): [
        r"The world tube of a ring of free particles on $r = ct_0$ ($t$ from $t_0/2$ to $3t_0/2$), each "
        r"circle at the height of its $t$, of proper radius $(t/t_0)^2ct_0$."],
    (METRIC, "ring"): [
        r"A ring of free particles on $r = ct_0$ ($t$ from $t_0/2$ to $3t_0/2$), of proper radius $(t/t_0)^2ct_0$, inside the central observer's event horizon, of proper radius $ct$, "
        r"until $t = t_0$ and outside it after."],
}
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


def published_radial(em, src, system):
    """The coordinate speed of a radial light ray, sqrt(-g_00/g_rr), of the published chart as a
    numpy function of its time, at the drawn p and t_0 with c = 1, and sqrt(g_rr), the scale."""
    import sympy as sp

    _, entry, R = em.nr.load(METRIC, system)
    src.note(METRIC, system, em.FIELDS)
    g = em.nr.published_matrix(R, entry, "metric_components")
    at = {R.c: 1, **{R.parameters[k]: sp.sympify(v) for k, v in PARAMS.items()}}
    i, j = entry["coords"].index(entry["coords"][0]), entry["coords"].index("r")
    clock = R.symbol[entry["coords"][0]]
    g00, grr = (sp.simplify(g[k, k].subs(R.held).doit().subs(at)) for k in (i, j))
    if (g00.free_symbols | grr.free_symbols) - {clock}:
        raise SystemExit(f"{METRIC}/{system}: g_00 = {g00} and g_rr = {grr} depend on more than its time")
    return (sp.lambdify(clock, sp.sqrt(-g00 / grr), "numpy"), sp.lambdify(clock, sp.sqrt(grr), "numpy"))


def embedding(ck, src, em):
    """The equator theta = pi/2 of the cosmic chart from t = t_0/2 to 3t_0/2, flat at every
    moment, a^2(dr^2 + r^2 dphi^2) with a = (t/t_0)^2, with a ring of free particles at rest on
    r = c t_0, which the published Christoffel symbols keep at rest: the circle of proper radius
    (t/t_0)^2 c t_0. Beside it is the event horizon of the observer at the centre, the circle of
    comoving radius c times the integral of the radial light ray's dr/dt = c/a from t to infinity,
    read from the published metric and checked to be c t_0^2/t, the line r = -c eta of the
    conformal chart, whose published scale factor at that eta is checked to be the cosmic chart's
    at t. Its proper radius c t grows at a steady rate while the ring's grows faster and faster,
    and the ring leaves the horizon at t = t_0, where a = 1. Drawn first as the ring's world tube,
    each circle at the height of its time, and beside it as a movie of the flat moments, a frame
    every t_0/40."""
    from scipy.integrate import quad

    em.at_rest(ck, src, METRIC, "cosmic")
    speed, scale = published_radial(em, src, "cosmic")
    conformal_speed, conformal_scale = published_radial(em, src, "conformal")
    alpha = np.linspace(0, 2 * math.pi, 361)
    moments, at = em.ring_moments(list(NAMED), STEP,
                                  lambda k: (f"$t/t_0 = {NAMED[k]:g}$", NAMED[k], {"t": repr(NAMED[k]), "theta": "pi/2"}),
                                  lambda t: (f"$t/t_0 = {t:g}$", t, {"t": repr(t), "theta": "pi/2"}))
    frames, rings, horizons = [], [], []
    for label, time, fixed_at in moments:
        where = f"extended inflation, {label}"
        sl = em.Slice(src, METRIC, "cosmic", "r", r"\phi", fixed_at, PARAMS)
        ck.plane(f"{where}: the equator", sl, np.linspace(1e-3, 20, 200))
        a = math.sqrt(float(sl.gxx))
        ck.add(f"{where}: the scale factor is (t/t_0)^2", abs(a - time ** 2), 1e-12)
        ck.add(f"{where}: the scale factor read along t is the slice's", abs(float(scale(time)) - a), 1e-12)
        plane = em.disc(sl, "plane", TOP, "the centre $r=0$", r"the plane extends to $r\to\infty$", [], SIZE)
        ck.isometry(where, plane)
        ck.form(f"{where}: level", plane, lambda r: 0 * r, SIZE)
        ck.radius(f"{where}: rho = a r", plane, lambda r: time ** 2 * r, SIZE)
        rho = float(sl.rho_at(RING))
        ring = np.column_stack([rho * np.cos(alpha), rho * np.sin(alpha), np.zeros_like(alpha)])
        ck.on_piece(f"{where}, the ring", plane, ring)
        ck.add(f"{where}: the ring's proper radius is (t/t_0)^2 c t_0", abs(rho - time ** 2 * RING), 1e-12)
        # The comoving distance light covers from t on, the horizon of the observer at r = 0.
        comoving, _ = quad(lambda s: float(speed(s)), time, np.inf, epsabs=1e-13, epsrel=1e-12, limit=200)
        ck.add(f"{where}: the horizon's comoving radius is c t_0^2/t", abs(comoving - 1 / time), 1e-10)
        ck.add(f"{where}: the horizon is the conformal chart's r = -c eta, at the eta whose scale factor is a",
               abs(float(conformal_scale(-comoving)) - a), 1e-9)
        ck.add(f"{where}: radial light keeps r + c eta in the conformal chart",
               abs(float(conformal_speed(-comoving)) - 1), 1e-12)
        horizon = sl.rho_at(comoving) * np.column_stack([np.cos(alpha), np.sin(alpha)])
        horizon = np.column_stack([horizon, np.zeros_like(alpha)])
        ck.on_piece(f"{where}, the event horizon", plane, horizon)
        ck.add(f"{where}: the horizon's proper radius is c t", abs(float(sl.rho_at(comoving)) - time), 1e-10)
        rings.append(rho)
        horizons.append(float(sl.rho_at(comoving)))
        frames.append(em.Surface([plane], label=label, time=time,
                                 curves=[em.Curve(plane, "particles", ring, closed=True),
                                         em.Curve(plane, "horizon", horizon, closed=True)],
                                 dots=[(plane, "particles", Q) for Q in ring[:-1:30]]))
    surfaces = [frames[i] for i in at]
    ck.exact("extended inflation: the ring lies inside the horizon at t = t_0/2, on it at t_0 and outside it at 3t_0/2",
             rings[at[0]] < horizons[at[0]] and rings[at[2]] > horizons[at[2]])
    ck.add("extended inflation: the ring meets the horizon at t = t_0", abs(rings[at[1]] - horizons[at[1]]), 1e-10)
    # The frames stand a steady step of t apart: the horizon's steps are equal and the ring's lengthen.
    ck.add("extended inflation: the horizon grows by c times each step of t",
           float(np.max(np.abs(np.diff(horizons) - np.diff([f.time for f in frames])))), 1e-9)
    ck.exact("extended inflation: the ring grows faster and faster", bool(np.all(np.diff(rings, 2) > 0)))
    ck.exact("extended inflation: the disc holds the ring at its widest", max(rings) < TOP)

    def plane_at(u, a, b):
        sl = em.Slice(src, METRIC, "cosmic", "r", r"\phi", {"t": repr(float(u)), "theta": "pi/2"}, PARAMS)
        return float(np.max(np.abs(float(sl.rho_at(RING)) - np.array([a, b]))))
    tube = em.stack(ck, "extended inflation", surfaces, lambda u: (RING * scale(u),) * 2, LIFT, plane_at, SIZE,
                    "the world tube continues before $t = t_0/2$ and after $t = 3t_0/2$")
    tube_fig = em.stack_figure(tube, SIZE, "$t$", [
        ("fill", "cover", "the ring at each $t$"),
        ("line", "particles", "the ring at three moments, twelve particles marked"),
        ("line", "worldline", "their world lines"),
        ("line", "axis", "the axis of time")])
    fig, played = em.ring_movie(frames, SIZE, "$t/t_0$")
    fig.legend("fill", "cover", r"the equatorial plane $\theta = \pi/2$")
    fig.legend("line", "particles", "the ring $r = ct_0$, twelve particles marked")
    fig.legend("line", "horizon", "the central observer's event horizon")
    fig.legend("line", "meridian", "lines from the centre every $30°$")
    return [em.view("tube", "The ring's world tube", "$ct_0$", [tube], tube_fig, settings=SETTINGS,
                    height=f"$t$, a height of ${LIFT:g}\\,ct_0$ for each $t_0$"),
            em.view("ring", "A ring of particles", "$ct_0$", surfaces, fig.done(), movie=played, settings=SETTINGS)]


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
        # Every moment of t is the whole line eta = -t_0^2/t at p = 2.
        for moment in cf.slices.moments(METRIC):
            v.slice(moment, lines=[cf.mink_pq(np.full_like(rr, -1 / moment.time), rr)])
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
