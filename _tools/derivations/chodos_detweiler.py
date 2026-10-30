"""The compact five-dimensional vacuum of Chodos and Detweiler.

Chyba (1985), page 868, (6.5)-(6.7), fixes the metric, the four periods and
the spatial rescaling. The derivation and source provenance are in
chodos_detweiler.md. Shared generators register these page-local definitions.
"""

METRIC = "chodos_detweiler"
CHARTS = ("synchronous", "rescaled")
DIMENSIONS = {
    (METRIC, "synchronous"): {"t": "T", "x": "L", "y": "L", "z": "L", "w": "L",
                             "t_0": "T", "L": "L"},
    (METRIC, "rescaled"): {"t": "T", "X": "L", "Y": "L", "Z": "L", "W": "L",
                          "t_0": "T", "L": "L", "\\tau": "T"},
}


def chart_spec(system):
    """Both five-dimensional charts, with tensor time x^0 = ct."""
    import verify_metrics as vm

    rescaled = system == "rescaled"
    space = ["X", "Y", "Z", "W"] if rescaled else ["x", "y", "z", "w"]
    coords = ["t", *space]
    parameters = ["t_0", "L"] + (["\\tau"] if rescaled else [])
    scale = "\\tau" if rescaled else "t_0"
    probe = vm.Reader(coords, parameters, ())
    time = probe.symbol["t"]

    def line(clock, lapse):
        terms = " + ".join(f"d{x}^2" for x in space[:3])
        return (f"ds^2 = -{lapse}dt^2 + \\dfrac{{{clock}}}{{{scale}}}"
                f"\\left({terms}\\right) + \\dfrac{{{scale}}}{{{clock}}}d{space[3]}^2")

    if rescaled:
        growing = "L\\sqrt{\\tau/t_0}"
        shrinking = "L\\sqrt{t_0/\\tau}"
    else:
        growing = shrinking = "L"
    domains = ["t \\in (0, \\infty)"]
    domains += [f"{x} \\in [0, {growing})" for x in space[:3]]
    domains += [f"{space[3]} \\in [0, {shrinking})", "t \\to 0 \\;\\text{(the singularity)}"]
    return {
        "metric_id": METRIC,
        "system": {"id": system, "name": "Rescaled" if rescaled else "Synchronous",
                   "coords": coords, "domains": domains, "parameters": parameters,
                   "line_element": line("t", "c^2")},
        "chart_line_element": line("\\left(t/c\\right)", ""),
        "time": "t",
        "printer": {"lead": [time], "factors": [probe.c, time, *probe.parameters.values()]},
        "check": lambda chart: chart_check(chart, system),
    }


def chart_check(chart, system):
    """The source coefficients and physical periods at distinct rational points.

    Only metric values are sampled here; tensor identities belong to the
    generated mathematics and the scoped metric checker.
    """
    import sympy as sp

    R = chart.reader
    clock = R.symbol["t"]
    for c, t, t0, tau, length in ((2, 3, 2, 5, 7), (3, 5, 7, 2, 11), (5, 2, 3, 7, 13)):
        c, t, t0, tau, length = map(sp.Rational, (c, t, t0, tau, length))
        scale = tau if system == "rescaled" else t0
        at = {clock: c*t, R.c: c, R.parameters["t_0"]: t0, R.parameters["L"]: length}
        if system == "rescaled":
            at[R.parameters["tau"]] = tau
        want = sp.diag(-1, t/scale, t/scale, t/scale, scale/t)
        got = chart.geo.g.subs(at)
        if got != want:
            raise AssertionError(f"{METRIC}/{system}: source coefficients disagree at {at}")
        periods = ([length*sp.sqrt(tau/t0)]*3 + [length*sp.sqrt(t0/tau)]
                   if system == "rescaled" else [length]*4)
        squares = [got[i, i]*periods[i-1]**2 for i in range(1, 5)]
        if squares != [length**2*t/t0]*3 + [length**2*t0/t]:
            raise AssertionError(f"{METRIC}/{system}: physical circumferences disagree at {at}")


def diagram_specs(Diagram, sideways):
    """One complete spatial period on each of the two distinct null planes."""
    import math

    rows = []
    for system in CHARTS:
        space = ["X", "Y", "Z", "W"] if system == "rescaled" else ["x", "y", "z", "w"]
        params = {"t_0": 1, "L": "2*pi"}
        if system == "rescaled":
            params["tau"] = 4
        for coordinate, period in ((space[0], 4*math.pi if system == "rescaled" else 2*math.pi),
                                   (space[3], math.pi if system == "rescaled" else 2*math.pi)):
            rows.append(Diagram(METRIC, system, "t"+coordinate, f"$t$ and ${coordinate}$",
                                ("t", coordinate), (0, period, 0, 2*math.pi),
                                f"${coordinate}/(ct_0)$", "$t/t_0$", params,
                                {x: "0" for x in space if x != coordinate},
                                families=sideways, periodic=(coordinate,), cones=(5, 6)))
    return rows


def diagram_captions():
    """Compact plane labels, with each coordinate identification explicit."""
    return {
        (METRIC, "synchronous", "tx"): [
            "The plane of $t$ and $x$ ($y=z=w=0$), with $x=0$ identified with $x=L$. "
            "Null geodesics along an expanding circle."],
        (METRIC, "synchronous", "tw"): [
            "The plane of $t$ and $w$ ($x=y=z=0$), with $w=0$ identified with $w=L$. "
            "Null geodesics along the contracting circle."],
        (METRIC, "rescaled", "tX"): [
            "The plane of $t$ and $X$ ($Y=Z=W=0$), with $X=0$ identified with $X=2L$. "
            "Null geodesics along an expanding circle."],
        (METRIC, "rescaled", "tW"): [
            "The plane of $t$ and $W$ ($X=Y=Z=0$), with $W=0$ identified with $W=L/2$. "
            "Null geodesics along the contracting circle."],
    }


def null_forms():
    """Analytic integrals of the null slopes, used only by the ray checker."""
    import numpy as np

    forms = {}
    for system, scale, axes in (("synchronous", 1, ("x", "w")), ("rescaled", 4, ("X", "W"))):
        for coordinate, growing in ((axes[0], True), (axes[1], False)):
            def distance(t, growing=growing, scale=scale):
                return (2*np.sqrt(scale*t) if growing else 2*t**1.5/(3*np.sqrt(scale)))
            forms[(METRIC, system, "t"+coordinate)] = (
                lambda t, x, distance=distance: x-distance(t),
                lambda t, x, distance=distance: x+distance(t),
                lambda t, x: t > 0.02,
            )
    return forms


def embedding_views(ck, src, e):
    """Forty-one uniformly spaced moments of a local cut of the flat two-torus.

    At fixed t the cut 0 <= x <= pi, y=z=0 has metric t dx^2 + dw^2/t,
    w periodic with period 2 pi. It embeds as a cylinder of radius 1/sqrt(t)
    and axial height sqrt(t) x. Its ends remain separate because x is cut.
    """
    import math
    import numpy as np

    frames = []
    times = np.linspace(0.25, 4, 41)
    size = 7.0
    for time in times:
        sl = e.Slice(src, METRIC, "synchronous", "x", "w",
                     {"t": repr(float(time)), "y": 0, "z": 0}, {"t_0": 1, "L": "2*pi"})
        tube = e.Piece("cylinder", "sheet", sl, 0, math.pi, -math.pi*math.sqrt(time)/2, 1,
                       (("edge", "$x=0$"), ("edge", "$x=L/2$")),
                       [(0, "r", None), (math.pi, "r", None)], size)
        ck.isometry(f"{METRIC} t={time:g}", tube)
        ck.radius(f"{METRIC} t={time:g}: fifth-circle radius", tube,
                  lambda x, time=time: np.ones_like(x)/math.sqrt(time), size)
        ck.form(f"{METRIC} t={time:g}: expanding axial length", tube,
                lambda x, time=time: math.sqrt(time)*(x-math.pi/2), size)
        frames.append(e.Surface([tube], label=f"$t/t_0 = {time:.5f}$", time=float(time)))
    fig = e.movie_figure(frames, {"sheet": "cover"}, size)
    fig.legend("fill", "cover", "the cut $y=z=0$, $0\\le x\\le L/2$")
    fig.legend("line", "r", "the two cut ends")
    fig.legend("line", "meridian", "$w$ constant")
    return [e.view("cylinder", "Expanding cut and shrinking circle", "$ct_0$", [frames[0], frames[-1]],
                   fig.done(), movie=e.movie(frames, "$t/t_0$", times),
                   settings="$ct_0=1$, $L=2\\pi ct_0$.")]


def embedding_captions():
    return {(METRIC, "cylinder"): [
        "A cylinder cut from the $x,w$ torus ($y=z=0$), with the full fifth circle. "
        "Axial length $\\pi ct_0\\sqrt{t/t_0}$ and radius $ct_0\\sqrt{t_0/t}$."]}


def slice_marks(slices, system, view):
    """The two displayed moments of the cylinder on each spacetime plane."""
    import math

    scale = 2 if system == "rescaled" and view == "tX" else 1
    reach = math.pi*scale if view in ("tx", "tX") else (math.pi if system == "rescaled" else 2*math.pi)
    return slices.one(METRIC, lambda moment: slices.along(moment.time, 0, reach))


def curvature_views(ck, src, cv):
    """A line figure read only from the stored Kretschmann scalar."""
    import math

    x_axis = cv.axis("$t/t_0$", 0.25, 4)
    scale = {"kind": "log", "from": 0.01, "to": 2000}
    panels = [{"id": "K", "name": "$Kc^4t_0^4$", "aspect": 1,
               "scale": scale, "ticks": cv.value_ticks(scale)}]
    stored = cv.nr_published(src, METRIC, "synchronous", "kretschmann", ["t"], ["t_0", "L"])
    curves = cv.line_curves(ck, METRIC, "K", "K", stored, "t", x_axis, scale, cv.PLOT_W,
                            {"t_0": 1, "L": 2*math.pi})
    return [{"id": "cosmic_time", "label": "Cosmic time", "kind": "line", "x": x_axis,
             "panels": panels, "curves": curves, "marks": [],
             "legend": [["line", "K", "five-dimensional vacuum"]],
             "caption": ["The Kretschmann scalar $K$, growing as $t^{-4}$ toward the initial singularity. "
                         "The Ricci scalar vanishes."],
             "settings": "$ct_0=1$, $L=2\\pi ct_0$.",
             "clip": "Logarithmic vertical scale; the full curve lies inside the panel."}]
