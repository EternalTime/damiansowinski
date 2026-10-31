"""The flat torus cosmologies and their two commuting holonomies.

Carlip (1992), (1.7)-(1.11), and Louko and Marolf (1994), (3.3),
(3.6) and (4.1)-(4.3), fix the families and coordinate maps.
"""

METRIC = "torus_universe_2plus1"
CHARTS = ("proper_time", "mean_curvature", "minkowski", "null_boundary", "static")
HOLONOMIES = ["a", "b", r"\lambda", r"\mu"]
DIMENSIONS = {
    (METRIC, system): {**{p: ("L" if p in ("a", "b") else "1") for p in HOLONOMIES},
                       **coords}
    for system, coords in (
        ("proper_time", {"t": "T", "x": "1", "y": "1"}),
        ("mean_curvature", {r"\tau": "1/L", "x": "1", "y": "1"}),
        ("minkowski", {"T": "T", "X": "L", "Y": "L"}),
        ("null_boundary", {"u": "L**2", "x": "1", "y": "1"}),
    )
}
DIMENSIONS[(METRIC, "static")] = {"t": "T", "x": "1", "y": "1", "r_1": "L", "r_2": "L", r"\phi": "1"}


def chart_spec(system):
    """The evolving quotient, its covering chart, and the static translation sector."""
    import verify_metrics as vm

    cell = [r"x \in [0,1)", r"y \in [0,1)"]
    orientation = r"a\mu-b\lambda>0"
    translation = r"\left(a\,dx+b\,dy\right)^2"
    rapidity = r"\left(\lambda\,dx+\mu\,dy\right)"
    parameters = list(HOLONOMIES)
    if system == "proper_time":
        coords = ["t", "x", "y"]
        domains = [r"t \in (0,\infty)", *cell, orientation]
        line = r"ds^2=-c^2dt^2+" + translation + "+c^2t^2" + rapidity + "^2"
        chart_line = r"ds^2=-dt^2+" + translation + "+t^2" + rapidity + "^2"
        name = "Proper Time"
    elif system == "mean_curvature":
        coords = [r"\tau", "x", "y"]
        domains = [r"\tau \in (0,\infty)", *cell, orientation]
        line = r"ds^2=-\frac{d\tau^2}{\tau^4}+" + translation + r"+\frac{1}{\tau^2}" + rapidity + "^2"
        chart_line, name = line, "Mean Curvature Time"
    elif system == "minkowski":
        coords = ["T", "X", "Y"]
        domains = [r"T \in (0,\infty)", r"X \in (-cT,cT)", r"Y \in (-\infty,\infty)", orientation]
        line, chart_line, name = r"ds^2=-c^2dT^2+dX^2+dY^2", r"ds^2=-dT^2+dX^2+dY^2", "Minkowski Cover"
    elif system == "null_boundary":
        coords = ["u", "x", "y"]
        domains = [r"u \in (0,\infty)", *cell, orientation]
        line = r"ds^2=-du" + rapidity + "+u" + rapidity + "^2+" + translation
        chart_line, name = line, "Null Boundary Coordinates"
    elif system == "static":
        coords, parameters = ["t", "x", "y"], ["r_1", "r_2", r"\phi"]
        domains = [r"t \in (-\infty,\infty)", *cell, "r_1>0", "r_2>0", r"0<\phi<\pi"]
        space = r"\left(r_1\,dx+r_2\cos\phi\,dy\right)^2+r_2^2\sin^2\phi\,dy^2"
        line, chart_line, name = r"ds^2=-c^2dt^2+" + space, r"ds^2=-dt^2+" + space, "Static Torus"
    else:
        raise ValueError(system)
    probe = vm.Reader(coords, parameters, ())
    return {
        "metric_id": METRIC,
        "system": {"id": system, "name": name, "coords": coords, "domains": domains,
                   "parameters": parameters, "line_element": line},
        "chart_line_element": chart_line,
        "time": coords[0] if system in ("proper_time", "minkowski", "static") else None,
        "printer": {"factors": list(probe.parameters.values())},
        "ricci_scalar": "0", "kretschmann": "0",
        "check": lambda chart: chart_check(chart, system),
    }


def chart_check(chart, system):
    """Numerical pullbacks of the source's developing maps at independent points."""
    import numpy as np
    import sympy as sp

    R = chart.reader
    clock, x, y = chart.symbols
    if system == "static":
        r1, r2, phi = (R.parameters[R._plain(p)] for p in ("r_1", "r_2", r"\phi"))
        image = [clock, r1*x+r2*sp.cos(phi)*y, r2*sp.sin(phi)*y]
    elif system == "minkowski":
        image = [clock, x, y]
    else:
        a, b, lam, mu = (R.parameters[R._plain(p)] for p in HOLONOMIES)
        q = lam*x+mu*y
        if system == "null_boundary":
            # Unit reference length in the map; it cancels from the pullback.
            image = [(sp.exp(q)+clock*sp.exp(-q))/2,
                     (sp.exp(q)-clock*sp.exp(-q))/2, a*x+b*y]
        else:
            t = 1/clock if system == "mean_curvature" else clock
            image = [t*sp.cosh(q), t*sp.sinh(q), a*x+b*y]
    J = sp.Matrix(image).jacobian(chart.symbols)
    parameters = list(R.parameters.values())
    metric_at = sp.lambdify([*chart.symbols, *parameters, R.c], chart.geo.g, "numpy")
    jacobian_at = sp.lambdify([*chart.symbols, *parameters, R.c], J, "numpy")
    for k, (t, u, v) in enumerate(((0.7, 0.13, 0.31), (1.4, 0.67, 0.21), (2.3, 0.29, 0.74))):
        values = {"a": 1.2+k/7, "b": 0.4-k/17, "lambda": 0.8+k/9, "mu": 1.7+k/11,
                  "r_1": 1.1+k/7, "r_2": 1.3+k/9, "phi": 0.7+k/8}
        arguments = [t, u, v, *[values[str(p)] for p in parameters], 1.0]
        jac = np.asarray(jacobian_at(*arguments), dtype=float)
        expected = jac.T @ np.diag([-1.0, 1.0, 1.0]) @ jac
        actual = np.asarray(metric_at(*arguments), dtype=float)
        if not np.allclose(actual, expected, rtol=2e-12, atol=2e-12):
            raise AssertionError(f"{METRIC}/{system}: the developing map disagrees at sample {k}")


EVOLVING = {"a": 1, "b": "1/3", "lambda": 1, "mu": "sqrt(2)"}
STATIC = {"r_1": 1, "r_2": "11/10", "phi": "pi/3"}


def projection_specs(Projection):
    """Both cuts of one complete torus cell, in every published chart."""
    heights = {"proper_time": "$z=ct$", "mean_curvature": r"$z=-\ell^2\tau$",
               "null_boundary": r"$z=u/(2\ell)$", "static": "$z=ct$"}
    return [Projection(METRIC, system, "cell", "The torus identifications", quotient_figure,
                       STATIC if system == "static" else EVOLVING,
                       input=(r"$\ell=1$, the length unit; Minkowski covering coordinates."
                              if system == "minkowski" else
                              r"$\ell=1$, the length unit; flat spatial developments at height " + heights[system] + "."))
            for system in CHARTS]


def projection_captions():
    return {(METRIC, system, "cell"): [caption] for system, caption in (
        ("proper_time", "The torus cell in proper time, with two pairs of identified sides."),
        ("mean_curvature", "The torus cell in mean curvature time, with the expanding future upward."),
        ("minkowski", "The Minkowski covering wedge, with paired boost and translation cuts."),
        ("null_boundary", "The torus cell in null boundary coordinates, with spacelike torus slices at $u>0$."),
        ("static", "The static torus cell, with paired spatial translations."),
    )}


def quotient_figure(spec):
    """Light cones from the full published 2+1 metric, with both deck generators.

    The Euclidean development at one time is intrinsic to that spatial slice.
    Its time dependent placement is only a drawing map of spacetime.
    The Minkowski view uses the actual developing map and boost-paired cuts.
    """
    import numpy as np
    import projections as pj

    system = spec.system
    coords = chart_spec(system)["system"]["coords"]

    class Development(pj.Slice):
        def to_drawing(self, point):
            t, x, y = np.broadcast_arrays(*[np.asarray(v, dtype=float) for v in point])
            if system == "minkowski":
                return np.stack([x, y, t], -1)
            if system == "static":
                X, Y, z = x+0.55*y, 1.1*np.sin(np.pi/3)*y, t
            else:
                proper = 1/t if system == "mean_curvature" else np.sqrt(t) if system == "null_boundary" else t
                X, Y = x+y/3, proper*(x+np.sqrt(2)*y)
                z = -t if system == "mean_curvature" else t/2 if system == "null_boundary" else t
            return np.stack([X, Y, z], -1)

        def jacobian(self, point):
            t, x, y = point
            if system == "minkowski":
                return super().jacobian(point)
            if system == "static":
                return np.array([[0, 1, 0.55], [0, 0, 1.1*np.sin(np.pi/3)], [1, 0, 0]])
            q = x+np.sqrt(2)*y
            proper = 1/t if system == "mean_curvature" else np.sqrt(t) if system == "null_boundary" else t
            slope = -1/t**2 if system == "mean_curvature" else 1/(2*np.sqrt(t)) if system == "null_boundary" else 1
            height = -1 if system == "mean_curvature" else 0.5 if system == "null_boundary" else 1
            return np.array([[0, 1, 1/3], [slope*q, proper, proper*np.sqrt(2)], [height, 0, 0]])

    sl = Development(METRIC, system, coords, "cartesian", spec.params)
    camera = pj.Camera(-45, 25)
    fig = pj.Figure(spec.view, spec.label, camera)

    def native(t, x, y):
        t, x, y = np.broadcast_arrays(t, x, y)
        if system == "minkowski":
            q = (x-0.5)+np.sqrt(2)*(y-0.5)
            return t*np.cosh(q), t*np.sinh(q), (x-0.5)+(y-0.5)/3
        clock = 1/t if system == "mean_curvature" else t*t if system == "null_boundary" else t
        return clock, x, y

    edge = np.linspace(0, 1, 161)
    # The first pair is x=0 and x=1, the second y=0 and y=1.
    for t in (0.5, 2.0):
        for side in (0, 1):
            fig.line("above", sl.to_drawing(native(t, side, edge)))
            fig.line("below", sl.to_drawing(native(t, edge, side)))
    times = np.linspace(0.5, 2, 161)
    for x, y in ((0, 0), (1, 0), (0, 1), (1, 1)):
        fig.line("floor", sl.to_drawing(native(times, x, y)))
    labels = ((r"$\ell\tau=2$", r"$\ell\tau=1/2$") if system == "mean_curvature" else
              (r"$u/\ell^2=1/4$", r"$u/\ell^2=4$") if system == "null_boundary" else
              (r"$ct/\ell=1/2$", r"$ct/\ell=2$"))
    for t, label in zip((0.5, 2.0), labels):
        fig.label(sl.to_drawing(native(t, 0, 0)), label, "l", dx=8)
    orientation = [-1, 0, 0] if system == "mean_curvature" else "tau" if system == "null_boundary" else "vector"
    cones = [pj.future_cone(sl, native(t, x, y), 0.19, orient=orientation)
             for t in (0.8, 1.5) for x, y in ((0.25, 0.25), (0.75, 0.75))]
    for apex, rim in sorted(cones, key=lambda c: camera.depth(c[0])):
        fig.cone(apex, rim)
    fig.legend("line", "above", "$x=0\\sim1$")
    fig.legend("line", "below", "$y=0\\sim1$")
    fig.legend("cone", "cone", "future light cones")
    points = np.array([p for layer in fig.layers for p in layer.get("points", [layer.get("at")])])
    width, height = np.ptp(points, axis=0)
    # Enlarge the window wherever a narrow view needs more room.
    pad = max(0.10, (max(width, height)-1.7*min(width, height))/(1.4*max(width, height)))
    return fig.done(pad=pad), sl


def curvature_views(ck, src, c):
    """Read the two identically vanishing stored scalars, with no curvature derivation."""
    views = []
    x_axis = c.axis(r"$ct/\ell$", 0.25, 2.0)
    scale = {"kind": "linear", "from": -1, "to": 1}
    for field, symbol, power, name in (("kretschmann", "K", 4, "Kretschmann scalar"),
                                      ("ricci_scalar", "R", 2, "Ricci scalar")):
        scalar = c.nr_published(src, METRIC, "proper_time", field, ["t"], [])
        panel = {"id": symbol, "name": f"${symbol}\\ell^{power}$", "aspect": 1,
                 "scale": scale, "ticks": c.value_ticks(scale)}
        curves = c.line_curves(ck, f"torus {symbol}", symbol, symbol, scalar, "t", x_axis,
                               scale, c.PLOT_W, {})
        c.shape_checks(ck, f"torus {symbol}", [1])
        views.append({"id": symbol, "label": f"${symbol}=0$", "kind": "line", "x": x_axis,
                      "panels": [panel], "curves": curves, "marks": [],
                      "legend": [["line", symbol, f"${symbol}=0$"]],
                      "caption": [f"The {name} ${symbol}$ along proper time, zero throughout the torus universe."],
                      "settings": r"$\ell=1$, the length unit.", "clip": "Linear scale."})
    return views
