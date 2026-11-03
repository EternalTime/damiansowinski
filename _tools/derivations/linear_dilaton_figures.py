"""Draw the linearly expanding Einstein frame at q = 1.

The moment tau = 1/q is a flat plane.  The conformal time at that moment
is eta = 0.  Both follow from the metric and time map in linear_dilaton.md.
"""

METRIC = "linear_dilaton"
SYSTEMS = ("cosmic", "conformal")


def embedding_views(ck, src, em):
    """An isometric plane at the moment tau = 1/q."""
    import numpy as np

    sl = em.Slice(src, METRIC, "cosmic", "x", None, {"tau": 1, "z": 0}, {"q": 1}, turn="y")
    ck.plane("linear dilaton: the spatial plane", sl, np.linspace(0.001, 2, 80))
    plane = em.Piece("plane", "sheet", sl, 0, 2, 0, 1,
                     (("axis", "$x=y=0$"), ("edge", "The plane extends to infinity.")),
                     [(v, "r", None) for v in (0.5, 1, 1.5, 2)], 5.0)
    ck.isometry("linear dilaton: spatial plane", plane)
    ck.radius("linear dilaton: radius", plane, lambda x: x, 5.0)
    ck.form("linear dilaton: height", plane, lambda x: 0*x, 5.0)
    surface = em.Surface([plane])
    fig = em.figure_of([surface], {"sheet": "cover"}, 5.0, em.FLAT_CAMERA)
    fig.legend("fill", "cover", r"$\tau=1/q$, $z=0$")
    fig.legend("line", "r", "Circles of constant radius")
    fig.legend("line", "meridian", "Radial lines")
    data = fig.done()
    data["box"] = [-2.5, 2.5, -2.5, 2.5]
    return [em.view("plane", "A spatial plane", "$1/q$", [surface], data,
                    settings=r"$q=1$, $\tau=1/q$.")]


def embedding_captions():
    return {(METRIC, "plane"): [r"A flat spatial plane ($\tau=1/q$, $z=0$)."]}


def slice_marks(slices, system):
    """The same moment on both null planes."""
    clock = 1 if system == "cosmic" else 0
    def lines(moment):
        hi = moment.reach("cosmic", "x")[1]
        return slices.across(clock, 0, hi)
    return slices.one(METRIC, lines, label=r"$\tau=1/q$", view_id="plane")


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
        moment = cf.slices.moments(METRIC, "plane", label=r"$\tau=1/q$")[0]
        reach = moment.reach("cosmic", "x")[1]
        x = np.linspace(-reach, reach, 81)
        v.slice(moment, [cf.mink_pq(np.zeros_like(x), x)])
        v.set(settings=r"$q=1$, $p=\arctan(\eta-x)$, $v=\arctan(\eta+x)$.",
              restriction=r"The plane $y=z=0$, with each point a single event.")
        views.append(v)
    return views


def conformal_captions():
    return {(METRIC, system): [r"The Einstein frame on the plane $y=z=0$ ($q=1$)."] for system in SYSTEMS}


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
