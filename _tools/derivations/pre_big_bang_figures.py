"""Draw the pre-big bang dilaton vacuum at its expanding isotropic member, every exponent -1/sqrt3.

The moment t = -1, where every scale factor is one, is a flat space.  In the conformal time of
the same member, k = -1/(1 + sqrt3), it is eta = -(1 + k).  Along x a light ray keeps
x -+ (-t)^s/s with s = 1 - beta = 1 + 1/sqrt3, which is x +- eta.  pre_big_bang.md beside this
file records each.
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
# The moment t = -1 in each chart's time: eta = -(1 + k) = -1/s.
MOMENT = {"synchronous": -1.0, "conformal": -1 / S}
REACH = 2.0


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


def slice_marks(slices, system):
    """The moment t = -1 on the plane of each chart."""
    def lines(moment):
        hi = moment.reach("synchronous", "x")[1]
        return slices.across(MOMENT[system], 0, hi)
    return slices.one(METRIC, lines, label="$t = -1$", view_id="plane")


def embedding_views(ck, src, em):
    """The plane z = 0 at t = -1, where every scale factor is one, isometric to Euclid's plane."""
    import numpy as np

    sl = em.Slice(src, METRIC, "synchronous", "x", None, {"t": -1, "z": 0}, PARAMS["synchronous"], turn="y")
    ck.plane("pre-big bang: the spatial plane", sl, np.linspace(0.001, REACH, 80))
    plane = em.Piece("plane", "sheet", sl, 0, REACH, 0, 1,
                     (("axis", "$x = y = 0$"), ("edge", "The plane extends to infinity.")),
                     [(v, "r", None) for v in (0.5, 1, 1.5, 2)], 5.0)
    ck.isometry("pre-big bang: spatial plane", plane)
    ck.radius("pre-big bang: radius", plane, lambda x: x, 5.0)
    ck.form("pre-big bang: height", plane, lambda x: 0 * x, 5.0)
    surface = em.Surface([plane])
    fig = em.figure_of([surface], {"sheet": "cover"}, 5.0, em.FLAT_CAMERA)
    fig.legend("fill", "cover", "$t = -1$, $z = 0$")
    fig.legend("line", "r", "Circles of constant radius")
    fig.legend("line", "meridian", "Radial lines")
    data = fig.done()
    data["box"] = [-2.5, 2.5, -2.5, 2.5]
    return [em.view("plane", "A spatial plane", "$\\ell$", [surface], data,
                    settings="$\\beta_i = -1/\\sqrt3$, $t = -1$, with $\\ell$ the unit of $t$ and of every length.")]


def embedding_captions():
    return {(METRIC, "plane"): [r"A flat spatial plane at $t = -1$, where every scale factor is one ($z = 0$)."]}


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
    moment = cf.slices.moments(METRIC, "plane", label="$t = -1$")[0]
    reach = moment.reach("synchronous", "x")[1]
    x = np.linspace(-reach, reach, 81)
    v.slice(moment, [cf.mink_pq(np.full_like(x, MOMENT["conformal"]), x)])
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
