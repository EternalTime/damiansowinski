"""Draw the linearly expanding Einstein frame at q = 1.

Every moment of cosmic time is flat, q^2 tau^2 (dx^2 + dy^2 + dz^2), and the dilaton
Phi = -q eta = -ln(q tau) takes one value over the whole of it, so a moment drawn alone is
Euclid's plane and no height over it varies.  What the metric adds is the scale factor q tau:
a ring of free particles at rest in the chart on x^2 + y^2 = 1/q^2 is the circle of proper
radius tau.  The embedding draws it in the form of the Big Rip's ring of galaxies, the ring's
world tube, each circle at the height of its tau, part of the cone with its apex at tau = 0,
and beside it the flat moments as a movie, q tau from 1/2 to 2.  In the conformal chart the
same moments are eta = ln(q tau)/q.  linear_dilaton.md beside this file records each choice.
"""
import math

METRIC = "linear_dilaton"
SYSTEMS = ("cosmic", "conformal")
PARAMS = {"q": 1}
SETTINGS = r"$q = 1$, with $1/q$ the unit of $\tau$ and of every length."

# The embedding: the key moments q tau = 1/2, 1, 3/2 and 2, the ring of radius 1/q at q tau = 1,
# where the scale factor is one, and between them a frame every 1/40 of q tau.
NAMED = ("1/2", "1", "3/2", "2")
STEP = 0.025
RING = 1.0          # the ring's comoving radius, 1/q, so that its proper radius is tau
TOP = 2.5           # the disc's proper radius, past the ring's widest, 2/q at q tau = 2
SIZE = 2 * TOP
LIFT = 2.0          # the height of 1/q of tau up the stack, in 1/q


def chart_time(system, tau):
    """The time of the chart at the moment tau of cosmic time, at q = 1: tau itself, or the
    conformal time eta = ln(q tau)/q."""
    return tau if system == "cosmic" else math.log(tau)


def at_rest(ck, src, em, system):
    """The published Christoffel symbols of the chart have no Gamma^i_00 for the chart's time,
    so a free particle at rest in the chart stays at rest, with the chart's time as its clock."""
    gamma, _ = em.published_christoffel(src, METRIC, system)
    time = r"\tau" if system == "cosmic" else r"\eta"
    ck.exact(f"{METRIC}/{system}: the chart's time is {time}, and its Gamma^{time}_xx is published",
             (time, "x", "x") in gamma)
    ck.exact(f"{METRIC}/{system}: no published Gamma^i_{time}{time}, so particles at rest in the chart stay there",
             not any(ix[1:] == (time, time) and ix[0] != time for ix in gamma))


def embedding_views(ck, src, em):
    """The plane z = 0 of the cosmic chart from q tau = 1/2 to 2, flat at every moment,
    q^2 tau^2 (dx^2 + dy^2), with a ring of free particles at rest on x^2 + y^2 = 1/q^2, which
    the published Christoffel symbols of both charts keep at rest: the circle of proper radius
    tau, growing at a steady rate. The conformal chart's plane at eta = ln(q tau)/q is checked
    to be the same plane, its scale factor e^(q eta) the cosmic chart's q tau. Drawn first as
    the ring's world tube, each circle at the height of its tau, which puts every row on the cone
    with its apex at tau = 0, where the published Ricci scalar R = 6/tau^2 diverges, which is
    checked, and beside it as a movie of the flat moments, a frame every 1/(40q) of tau."""
    import numpy as np
    import sympy as sp

    for system in SYSTEMS:
        at_rest(ck, src, em, system)
    # The tube's caption names the published Ricci scalar of the cosmic chart at the cone's apex.
    src.note(METRIC, "cosmic", ["ricci_scalar"])
    _, entry, reader = em.nr.load(METRIC, "cosmic")
    scalar = reader(entry["ricci_scalar"].split("=", 1)[1]).subs(reader.c, 1)
    tau = reader.symbol[r"\tau"]
    ck.exact("linear dilaton: the published Ricci scalar of the cosmic chart is 6/tau^2, as the tube's caption says",
             sp.simplify(scalar - 6 / tau ** 2) == 0)
    ck.exact("linear dilaton: the published Ricci scalar diverges as tau falls to zero, the apex of the tube's cone",
             sp.limit(scalar, tau, 0, "+") == sp.oo)
    keys = [float(sp.Rational(s)) for s in NAMED]
    moments, at = em.ring_moments(keys, STEP,
                                  lambda k: (f"$q\\tau = {keys[k]:g}$", keys[k], NAMED[k]),
                                  lambda t: (f"$q\\tau = {t:g}$", t, repr(t)))
    alpha = np.linspace(0, 2 * math.pi, 361)
    frames, radii = [], []
    for label, time, fixed in moments:
        where = f"linear dilaton, {label}"
        sl = em.FlatPlane(src, METRIC, "cosmic", "x", "y", {"tau": fixed, "z": 0}, PARAMS)
        ck.add(f"{where}: the plane's scale factor is q tau", float(np.max(np.abs(sl.scale - time))), 1e-12)
        conformal = em.FlatPlane(src, METRIC, "conformal", "x", "y", {"eta": repr(chart_time("conformal", time)), "z": 0}, PARAMS)
        ck.add(f"{where}: the conformal chart's plane at eta = ln(q tau)/q has the same scale factor",
               float(np.max(np.abs(conformal.scale - sl.scale))), 1e-12)
        plane = em.disc(sl, "plane", TOP, "the centre of the ring", "the plane continues, flat, to infinity", [], SIZE)
        ck.isometry(where, plane)
        ring, dots = em.particles(ck, where, sl, plane, RING * np.cos(alpha), RING * np.sin(alpha))
        radius = np.hypot(ring.points[:, 0], ring.points[:, 1])
        ck.add(f"{where}: the ring's proper radius is tau", float(np.max(np.abs(radius - RING * time))), 1e-12)
        radii.append(float(radius[0]))
        frames.append(em.Surface([plane], label=label, time=time, curves=[ring], dots=dots))
    surfaces = [frames[i] for i in at]
    # The frames stand a steady step of tau apart, and the ring's radius grows by the same step each time.
    ck.add("linear dilaton: the ring grows at a steady rate, by each step of tau",
           float(np.max(np.abs(np.diff(radii) - RING * np.diff([f.time for f in frames])))), 1e-12)
    ck.exact("linear dilaton: the disc holds the ring at its widest", max(radii) < TOP)

    def plane_at(u, a, b):
        sl = em.FlatPlane(src, METRIC, "cosmic", "x", "y", {"tau": repr(float(u)), "z": 0}, PARAMS)
        return float(np.max(np.abs(RING * sl.scale - [a, b])))
    tube = em.stack(ck, "linear dilaton", surfaces, lambda u: (RING * np.asarray(u, dtype=float),) * 2, LIFT,
                    plane_at, SIZE, r"the world tube continues before $q\tau = 1/2$ and after $q\tau = 2$")
    piece = tube.pieces[0]
    # Every row stands at the height LIFT tau with the radius tau: on the cone rho = z/LIFT, whose apex
    # on the axis is the moment tau = 0.
    ck.add("linear dilaton: every row of the tube lies on the cone with its apex at tau = 0",
           float(max(np.max(np.abs(piece.A - piece.Zr / LIFT)), np.max(np.abs(piece.B - piece.Zr / LIFT)))),
           10.0 ** -piece.decimals)
    tube_fig = em.stack_figure(tube, SIZE, r"$\tau$", [
        ("fill", "cover", r"the ring at each $\tau$"),
        ("line", "particles", "the ring at four moments, twelve particles marked"),
        ("line", "worldline", "their world lines"),
        ("line", "axis", "the axis of time")])
    fig, played = em.ring_movie(frames, SIZE, r"$q\tau$")
    fig.legend("fill", "cover", "the plane $z = 0$")
    fig.legend("line", "particles", "the ring $x^2 + y^2 = 1/q^2$, twelve particles marked")
    fig.legend("line", "meridian", "lines from the centre every $30°$")
    return [em.view("tube", "The ring's world tube", "$1/q$", [tube], tube_fig, settings=SETTINGS,
                    height=f"$\\tau$, a height of ${LIFT:g}/q$ for each $1/q$"),
            em.view("ring", "A ring of particles", "$1/q$", surfaces, fig.done(), movie=played, settings=SETTINGS)]


def embedding_captions():
    span = r"($q\tau$ from $1/2$ to $2$)"
    return {
        (METRIC, "tube"): [
            rf"The world tube of a ring of free particles {span}, each circle at the "
            r"height of its $\tau$ and of proper radius $\tau$, on a cone with its apex at $\tau = 0$, where the Ricci "
            r"scalar $R = 6/\tau^2$ diverges."],
        (METRIC, "ring"): [
            rf"A ring of free particles at rest in the chart on $x^2 + y^2 = 1/q^2$, in the flat plane $z = 0$ "
            rf"{span}, of proper radius $\tau$."],
    }


def slice_marks(slices, system):
    """Each moment tau of the embedding is the whole plane z = 0, so on the plane of the chart's
    time and x it is the line of constant time across the drawing: tau in the cosmic chart and
    eta = ln(q tau)/q in the conformal one, at q = 1."""
    return slices.one(METRIC, lambda m: slices.across(chart_time(system, m.time), 0.0, slices.BIG))


def diagram_specs(Diagram, sideways):
    """Null rays through the Cartesian plane in each Einstein-frame chart."""
    rows = []
    for system in SYSTEMS:
        clock = r"\tau" if system == "cosmic" else r"\eta"
        box = (-2.2, 2.2, 0.2, 4.6) if system == "cosmic" else (-2.2, 2.2, -2.2, 2.2)
        rows.append(Diagram(METRIC, system, "plane", "Cartesian plane", (clock, "x"),
                            box, "$qx$", "$q"+clock+"$", {"q": 1}, {"y": "0", "z": "0"},
                            tau="tau" if system == "cosmic" else "eta", families=sideways,
                            cones=(5, 5), step=0.003))
    return rows


def diagram_captions():
    return {(METRIC, system, "plane"): [
        (r"Null geodesics in the Einstein frame on the $\tau,x$ plane ($y=z=0$, $q=1$)."
         if system == "cosmic" else
         r"Null geodesics in the Einstein frame on the $\eta,x$ plane ($y=z=0$, $q=1$).")]
            for system in SYSTEMS}


def null_forms():
    import numpy as np
    return {(METRIC, "cosmic", "plane"): (lambda t, x: np.log(t)+x, lambda t, x: np.log(t)-x, None),
            (METRIC, "conformal", "plane"): (lambda t, x: t+x, lambda t, x: t-x, None)}


def conformal_views(ck, src, cf):
    """The Cartesian plane compactified by its two null coordinates."""
    import numpy as np

    views = []
    for system in SYSTEMS:
        cosmic = system == "cosmic"
        clock = r"\tau" if cosmic else r"\eta"
        plane = cf.Plane(src, METRIC, system, (clock, "x"), {"y": "0", "z": "0"}, {"q": 1})
        def pq(t, x):
            return cf.mink_pq(np.log(t) if cosmic else t, x)
        times = ck.uniform(0.1, 5) if cosmic else ck.uniform(-2, 2)
        ck.chart(METRIC+"/"+system, plane, pq, times, ck.uniform(-4, 4), lambda t, x: (1, 0))
        near, nearer = (0.01, 0.001) if cosmic else (np.log(0.01), np.log(0.001))
        ck.diverges(system+": past curvature", plane.kretschmann(near, 0), plane.kretschmann(nearer, 0))
        future = 1e6 if cosmic else np.log(1e6)
        ck.limit(system+": future curvature", [plane.kretschmann(future, 0)], [0], 1e-12)
        pi = np.pi
        v = cf.View(system, "Cartesian plane", [-pi-0.25, pi+0.25, -pi-0.25, pi+0.25], system)
        v.fill("region", cf.DIAMOND)
        v.fill("cover", cf.DIAMOND)
        v.line("singular", [[[-pi, 0], [0, -pi]], [[0, -pi], [pi, 0]]], zig=True)
        v.line("scri", [[[-pi, 0], [0, pi]], [[0, pi], [pi, 0]]])
        for x in (-2, -1, 0, 1, 2):
            eta = np.linspace(-8, 8, 401)
            v.curve("r", *cf.mink_pq(eta, np.full_like(eta, x)))
        for eta in (-2, -1, 0, 1, 2):
            x = np.linspace(-20, 20, 401)
            v.curve("t", *cf.mink_pq(np.full_like(x, eta), x))
        v.legend("r", "$x$ constant")
        v.legend("t", r"$\eta$ constant")
        v.legend("singular", "Past curvature singularity")
        v.legend("scri", "Future null infinity")
        # Each moment tau of the embedding is the whole plane z = 0, the line eta = ln(q tau)/q at every x.
        across = cf.spread(-np.inf, np.inf, 500, 12)
        for moment in cf.slices.moments(METRIC):
            v.slice(moment, [cf.mink_pq(np.full_like(across, chart_time("conformal", moment.time)), across)])
        v.set(settings=r"$q=1$, $p=\arctan(\eta-x)$, $v=\arctan(\eta+x)$.",
              restriction=r"The plane $y=z=0$, with each point a single event.")
        views.append(v)
    return views


def conformal_captions():
    return {(METRIC, system): [
        (r"The Einstein frame on the $\tau,x$ plane ($y=z=0$, $q=1$)." if system == "cosmic" else
         r"The Einstein frame on the $\eta,x$ plane ($y=z=0$, $q=1$).")] for system in SYSTEMS}


def curvature_views(ck, src, cv):
    """The stored Einstein-frame scalars through cosmic time."""
    xaxis = cv.axis(r"$q\tau$", 0.2, 4.4)
    scales = {"K": {"kind": "log", "from": 0.01, "to": 10000},
              "R": {"kind": "log", "from": 0.1, "to": 1000}}
    panels = [{"id": field, "name": "$"+field+("/q^4$" if field == "K" else "/q^2$"),
               "aspect": 2, "scale": scales[field], "ticks": cv.value_ticks(scales[field])}
              for field in ("K", "R")]
    cv.shape_checks(ck, METRIC, [panel["aspect"] for panel in panels])
    curves = []
    for field, scalar in (("K", "kretschmann"), ("R", "ricci_scalar")):
        stored = cv.nr_published(src, METRIC, "cosmic", scalar, ["tau"], ["q"])
        curves += cv.line_curves(ck, METRIC, field, field, stored, "tau", xaxis,
                                 scales[field], cv.PLOT_W/2, {"q": 1})
    return [{"id": "cosmic", "label": "Cosmic time", "kind": "line", "x": xaxis,
             "panels": panels, "curves": curves, "marks": [],
             "legend": [["line", "K", "$K$"], ["line", "R", "$R$"]],
             "caption": [r"The Kretschmann scalar $K$ and the Ricci scalar $R$ of the Einstein frame through cosmic time ($q=1$)."],
             "ricci_sentences": [], "settings": "$q=1$.", "clip": "Logarithmic vertical scales.",
             "per_scalar": cv.scalar_clips("Logarithmic scale.", "Logarithmic scale.")}]
