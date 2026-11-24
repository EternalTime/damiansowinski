"""Draw the pre-big bang dilaton vacuum.

The spacetime diagrams, the conformal diagram and the curvature figure draw the expanding
isotropic member, every exponent -1/sqrt3.  In its conformal time, k = -1/(1 + sqrt3), the moment
t = -1 is eta = -(1 + k).  Along x a light ray keeps x -+ (-t)^s/s with s = 1 - beta = 1 + 1/sqrt3,
which is x +- eta.

The embedding diagram draws the member (-2/7, -3/7, -6/7), every exponent negative and no two
equal, since a moment of the isotropic member is Euclid's plane with one scale and a ring on it
stays a circle.  The plane y = 0 holds the slowest and the fastest stretch, and a ring of free
particles at rest on x^2 + z^2 = l^2 at t = -1 is the ellipse of semi-axes (-t)^(-2/7) l and
(-t)^(-6/7) l from t = -1 to -1/4, drawn as its world tube and as a movie.  Its moments lie on
none of the other drawings, which draw another member.  pre_big_bang.md beside this file records
each.
"""
import math

METRIC = "pre_big_bang"
SYSTEMS = ("synchronous", "conformal")
BETA = "-1/sqrt(3)"
K = "-1/(1 + sqrt(3))"
PARAMS = {"synchronous": {"beta_1": BETA, "beta_2": BETA, "beta_3": BETA}, "conformal": {"k": K}}
CLOCKS = {"synchronous": ("t", "t"), "conformal": ("\\eta", "eta")}
VIEWS = {"synchronous": "tx", "conformal": "etax"}
S = 1 + 1 / math.sqrt(3)

# The member embedded: every exponent negative, so space expands along every axis, and the two
# exponents of the plane y = 0 the farthest apart, the slowest and the fastest stretch.
RING = {"beta_1": "-2/7", "beta_2": "-3/7", "beta_3": "-6/7"}
RING_BETAS = (-2 / 7, -3 / 7, -6 / 7)
# The ring's moments, from t = -1, where every scale factor is one and the ring is the circle of
# radius l, to t = -1/4, a frame every 1/40 of t; the stack's axis of time runs a quarter of its
# length past the last, to t = -1/16, short of the singularity at t = 0.
NAMED = ("-1", "-3/4", "-1/2", "-1/4")
STEP = 0.025
LIFT = 4.0          # the height of a unit of t up the stack
TOP = 3.5           # the disc's proper radius, past the ring's widest reach, (1/4)^(-6/7) = 3.28 l
SIZE = 2 * TOP

# The drawings on which no moment of the embedding lies, and why.
HIDDEN_REASON = ("the expanding isotropic member, every exponent -1/sqrt3, another member of the family than "
                 "(-2/7, -3/7, -6/7), whose ring is embedded")
HIDDEN = [(METRIC, system, VIEWS[system]) for system in SYSTEMS] + [(METRIC, "conformal")]


def eta_of(t):
    """The conformal time of the synchronous time t < 0."""
    return -(-t) ** S / S


def diagram_specs(Diagram, sideways):
    """Null rays through the plane y = z = 0 in each chart."""
    rows = []
    for system in SYSTEMS:
        clock, plain = CLOCKS[system]
        rows.append(Diagram(METRIC, system, VIEWS[system], f"${clock}$ and $x$", (clock, "x"),
                            (-2.2, 2.2, -4.4, 0.0), "$x$", f"${clock}$", PARAMS[system], {"y": "0", "z": "0"},
                            tau=plain, families=sideways, cones=(5, 5), step=0.003))
    return rows


def diagram_captions():
    return {
        (METRIC, "synchronous", "tx"): [
            r"Null geodesics on the $t,x$ plane of the expanding isotropic member ($y = z = 0$)."],
        (METRIC, "conformal", "etax"): [
            r"Null geodesics on the $\eta,x$ plane of the expanding isotropic member ($y = z = 0$)."],
    }


def null_forms():
    import numpy as np
    return {
        (METRIC, "synchronous", "tx"): (lambda t, x: x - (-t) ** S / S, lambda t, x: x + (-t) ** S / S,
                                        lambda t, x: t < -1e-3),
        (METRIC, "conformal", "etax"): (lambda t, x: x + t, lambda t, x: x - t, lambda t, x: np.asarray(t) < -1e-6),
    }


def ring_axes(u):
    """The semi-axes (-t)^beta_1 and (-t)^beta_3 of the ring at the time u, in units of l."""
    import numpy as np

    minus_t = -np.asarray(u, dtype=float)
    return minus_t ** RING_BETAS[0], minus_t ** RING_BETAS[2]


def embedding_views(ck, src, em):
    """The plane y = 0 of the member (-2/7, -3/7, -6/7) from t = -1 to t = -1/4, flat at every
    moment, (-t)^(-4/7) dx^2 + (-t)^(-12/7) dz^2, with a ring of free particles at rest on
    x^2 + z^2 = l^2, which the published Christoffel symbols keep at rest: the ellipse of semi-axes
    (-t)^beta_1 l along x and (-t)^beta_3 l along z, the circle of radius l at t = -1, where every
    scale factor is one. Drawn first as the ring's world tube, each ellipse at the height of its
    time, and beside it as a movie of the flat moments, a frame every 1/40 of t."""
    import numpy as np
    import sympy as sp

    betas = [sp.Rational(RING[f"beta_{i}"]) for i in (1, 2, 3)]
    ck.exact("pre-big bang: the exponents embedded lie on the sphere sum_i beta_i^2 = 1, every one negative",
             sum(b * b for b in betas) == 1 and all(b < 0 for b in betas))
    ck.add("pre-big bang: the floats of the exponents are the published fractions",
           max(abs(float(b) - f) for b, f in zip(betas, RING_BETAS)), 1e-15)
    em.at_rest(ck, src, METRIC, "synchronous")
    keys = [float(sp.Rational(s)) for s in NAMED]
    moments, at = em.ring_moments(keys, STEP,
                                  lambda k: (f"$t = {NAMED[k]}$", keys[k], {"t": NAMED[k], "y": 0}, None),
                                  lambda t: (f"$t = {t:g}$", t, {"t": repr(t), "y": 0}, None))
    frames = em.ring_sequence(ck, src, "pre-big bang", METRIC, "synchronous", moments, TOP, RING)
    surfaces = [frames[i] for i in at]
    alpha = np.linspace(0, 2 * math.pi, 361)
    for s, (_, time, _, _) in zip(frames, moments):
        a, b = ring_axes(time)
        want = np.column_stack([a * np.cos(alpha), b * np.sin(alpha)])
        ck.add(f"pre-big bang, t = {time}: the ellipse of semi-axes (-t)^beta_1 and (-t)^beta_3",
               float(np.max(np.abs(s.curves[0].points[:, :2] - want))), 1e-12)
    # Every semi-axis grows faster and faster: (-t)^beta has the second derivative
    # beta (beta - 1) (-t)^(beta - 2), positive for a negative beta.
    t = sp.Symbol("t", negative=True)
    ck.exact("pre-big bang: both semi-axes of the ring grow, and grow faster, as t rises toward zero",
             all(sp.simplify(sp.diff((-t) ** b, t) * (-t) ** (1 - b)) == -b and
                 sp.simplify(sp.diff((-t) ** b, t, 2) * (-t) ** (2 - b)) == b * (b - 1) for b in (betas[0], betas[2])))

    def plane(u, a, b):
        sl = em.FlatPlane(src, METRIC, "synchronous", "x", "z", {"t": repr(float(u)), "y": 0}, RING)
        return float(np.max(np.abs(sl.scale - [a, b])))
    tube = em.stack(ck, "pre-big bang", surfaces, ring_axes, LIFT, plane, SIZE,
                    "the world tube runs on before $t = -1$ and widens without bound toward $t = 0$")
    ck.exact("pre-big bang: the stack's axis of time ends short of the singularity at t = 0",
             tube.data()["axis"]["to"] < 0)
    tube_fig = em.stack_figure(tube, SIZE, "$t$", [
        ("fill", "cover", "the ring at each $t$"),
        ("line", "particles", "the ring at four moments, twelve particles marked"),
        ("line", "worldline", "their world lines"),
        ("line", "axis", "the axis of time")])
    fig, played = em.ring_movie(frames, SIZE, "$t$")
    fig.legend("fill", "cover", "the plane $y = 0$")
    fig.legend("line", "particles", "the ring $x^2 + z^2 = \\ell^2$, twelve particles marked")
    fig.legend("line", "meridian", "lines from the centre every $30°$")
    settings = "$(\\beta_1, \\beta_2, \\beta_3) = (-2/7, -3/7, -6/7)$, $\\ell$ the unit of $t$ and of every length."
    return [em.view("tube", "The ring's world tube", "$\\ell$", [tube], tube_fig, settings=settings,
                    height="$t$, a height of $4\\,\\ell$ for each unit of $t$"),
            em.view("ring", "A ring of particles", "$\\ell$", surfaces, fig.done(), movie=played, settings=settings)]


def embedding_captions():
    return {
        (METRIC, "tube"): [
            r"The world tube of a ring of free particles on the plane $y = 0$ ($t$ from $-1$ to $-1/4$), each "
            r"ellipse at the height of its $t$, widening faster and faster toward $t = 0$, fastest along $z$."],
        (METRIC, "ring"): [
            r"A ring of free particles ($t$ from $-1$ to $-1/4$), an ellipse of semiaxes "
            r"$(-t)^{\beta_1}\ell$ and $(-t)^{\beta_3}\ell$ along $x$ and $z$, growing faster and faster toward "
            r"$t = 0$."],
    }


def conformal_views(ck, src, cf):
    """The plane y = z = 0 of the expanding member in conformal time, the half of Minkowski's
    diamond below eta = 0."""
    import numpy as np

    pi = np.pi
    plane = cf.Plane(src, METRIC, "conformal", ("\\eta", "x"), {"y": "0", "z": "0"}, PARAMS["conformal"])
    ck.chart(METRIC + "/conformal", plane, cf.mink_pq, ck.uniform(-6, -0.01), ck.uniform(-4, 4),
             lambda t, x: (1, 0))
    ck.diverges("conformal: curvature at eta = 0", plane.kretschmann(-1e-3, 0), plane.kretschmann(-1e-4, 0))
    ck.limit("conformal: past curvature", [plane.kretschmann(-1e6, 0)], [0], 1e-12)
    v = cf.View("conformal", "$\\eta$ and $x$", [-pi - 0.25, pi + 0.25, -pi - 0.25, 0.25], "conformal")
    half = [[pi, 0], [0, -pi], [-pi, 0]]
    v.fill("region", half)
    v.fill("cover", half)
    v.line("singular", [[[-pi, 0], [pi, 0]]], zig=True)
    v.line("scri", [[[-pi, 0], [0, -pi]], [[0, -pi], [pi, 0]]])
    for x in (-2, -1, 0, 1, 2):
        eta = -np.geomspace(1e-6, 1e4, 401)
        v.curve("r", *cf.mink_pq(eta, np.full_like(eta, x)))
    for eta in (-3, -2, -1, -0.5):
        x = np.linspace(-40, 40, 801)
        v.curve("t", *cf.mink_pq(np.full_like(x, eta), x))
    v.legend("r", "$x$ constant")
    v.legend("t", r"$\eta$ constant")
    v.legend("singular", "Curvature singularity")
    v.legend("scri", "Past null infinity")
    v.set(settings=r"$k = -1/(1 + \sqrt3)$, $p = \arctan(\eta - x)$, $q = \arctan(\eta + x)$.",
          restriction=r"The plane $y = z = 0$, with each point a single event.")
    return [v]


def conformal_captions():
    return {(METRIC, "conformal"): [
        r"The expanding isotropic member on the $\eta,x$ plane ($y = z = 0$), its curvature singularity along the top."]}


def curvature_views(ck, src, cv):
    """The stored scalars of the expanding isotropic member as t rises to zero."""
    xaxis = cv.axis("$t$", -4.4, -0.2)
    scales = {"K": {"kind": "log", "from": 0.01, "to": 10000},
              "R": {"kind": "log", "from": 0.1, "to": 1000}}
    panels = [{"id": field, "name": f"${field}$", "aspect": 2, "scale": scales[field],
               "ticks": cv.value_ticks(scales[field])} for field in ("K", "R")]
    cv.shape_checks(ck, METRIC, [panel["aspect"] for panel in panels])
    curves = []
    for field, scalar in (("K", "kretschmann"), ("R", "ricci_scalar")):
        stored = cv.nr_published(src, METRIC, "synchronous", scalar, ["t"], ["beta_1", "beta_2", "beta_3"])
        curves += cv.line_curves(ck, METRIC, field, field, stored, "t", xaxis,
                                 scales[field], cv.PLOT_W / 2, {name: -1 / math.sqrt(3) for name in PARAMS["synchronous"]})
    return [{"id": "synchronous", "label": "Synchronous time", "kind": "line", "x": xaxis,
             "panels": panels, "curves": curves, "marks": [],
             "legend": [["line", "K", "$K$"], ["line", "R", "$R$"]],
             "caption": [r"The Kretschmann scalar $K$ and the Ricci scalar $R$ of the expanding isotropic member as $t$ rises to zero."],
             "ricci_sentences": [], "settings": r"$\beta_i = -1/\sqrt3$.", "clip": "Logarithmic vertical scales.",
             "per_scalar": cv.scalar_clips("Logarithmic scale.", "Logarithmic scale.")}]
