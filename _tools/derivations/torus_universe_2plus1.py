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
DIMENSIONS[(METRIC, "static")] = {"t": "T", "x": "1", "y": "1", "r_1": "L", "r_2": "L", r"\varphi": "1"}


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
        coords, parameters = ["t", "x", "y"], ["r_1", "r_2", r"\varphi"]
        domains = [r"t \in (-\infty,\infty)", *cell, "r_1>0", "r_2>0", r"0<\varphi<\pi"]
        space = r"\left(r_1\,dx+r_2\cos\varphi\,dy\right)^2+r_2^2\sin^2\varphi\,dy^2"
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
        r1, r2, phi = (R.parameters[R._plain(p)] for p in ("r_1", "r_2", r"\varphi"))
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
                  "r_1": 1.1+k/7, "r_2": 1.3+k/9, "varphi": 0.7+k/8}
        arguments = [t, u, v, *[values[str(p)] for p in parameters], 1.0]
        jac = np.asarray(jacobian_at(*arguments), dtype=float)
        expected = jac.T @ np.diag([-1.0, 1.0, 1.0]) @ jac
        actual = np.asarray(metric_at(*arguments), dtype=float)
        if not np.allclose(actual, expected, rtol=2e-12, atol=2e-12):
            raise AssertionError(f"{METRIC}/{system}: the developing map disagrees at sample {k}")
