"""Bowen-York maximal initial data, with intrinsic three-dimensional geometry only.

The free conformal factor is retained in every printed tensor.
Source equations and the two boundary prescriptions are recorded in bowen_york.md.
"""
import sympy as sp

CHARTS = ("cartesian", "spherical", "cylindrical")
DIMENSIONS = {
    ("bowen_york", "cartesian"): {"x": "L", "y": "L", "z": "L", r"\psi": "1"},
    ("bowen_york", "spherical"): {"r": "L", r"\theta": "1", r"\phi": "1", r"\psi": "1"},
    ("bowen_york", "cylindrical"): {r"\rho": "L", "z": "L", r"\phi": "1", r"\psi": "1"},
}


def chart_spec(system):
    """Cartesian and spherical charts, and the axisymmetric cylindrical pullback."""
    import chart_printer as cp
    import verify_metrics as vm

    if system == "cartesian":
        coords, function = ["x", "y", "z"], r"\psi = \psi(x,y,z)"
        domains = [c + r" \in (-\infty, \infty)" for c in coords]
        flat = "dx^2 + dy^2 + dz^2"
        lap = r"\partial_x^2\psi + \partial_y^2\psi + \partial_z^2\psi"
    elif system == "spherical":
        coords, function = ["r", r"\theta", r"\phi"], r"\psi = \psi(r,\theta,\phi)"
        domains = [r"r \in (0, \infty)", r"\theta \in (0, \pi)", r"\phi \in [0, 2\pi)"]
        flat = r"dr^2 + r^2\left(d\theta^2 + \sin^2\theta\,d\phi^2\right)"
        lap = (r"\partial_r^2\psi + \frac{2\partial_r\psi}{r}"
               r" + \frac{\partial_\theta^2\psi + \cot\theta\,\partial_\theta\psi}{r^2}"
               r" + \frac{\partial_\phi^2\psi}{r^2\sin^2\theta}")
    elif system == "cylindrical":
        coords, function = [r"\rho", "z", r"\phi"], r"\psi = \psi(\rho,z)"
        domains = [r"\rho \in (0, \infty)", r"z \in (-\infty, \infty)", r"\phi \in [0, 2\pi)"]
        flat = r"d\rho^2 + dz^2 + \rho^2d\phi^2"
        lap = r"\partial_\rho^2\psi + \frac{\partial_\rho\psi}{\rho} + \partial_z^2\psi"
    else:
        raise ValueError(system)
    line = r"ds^2 = \psi^4\left(" + flat + r"\right)"
    reader = vm.Reader(coords, [function], ())
    psi = reader.parameters["psi"]
    return {
        "metric_id": "bowen_york",
        "system": {"id": system, "name": system.title(), "coords": coords,
                   "domains": domains + [r"\text{Puncture centres or excised throat interiors are excluded}"],
                   "parameters": [function], "line_element": line},
        "chart_line_element": line,
        "printer": {"lead": [psi], "flip": False,
                    "collect": lambda poly, printer: cp.collect_by(poly, [psi], printer)},
        "ricci_scalar": r"-\frac{8\left(" + lap + r"\right)}{\psi^5}",
        "check": lambda chart: chart_check(chart, system),
    }


def chart_check(chart, system):
    """Compare the conformal scalar and Cartesian pullbacks at finite sample points.

    A nonharmonic polynomial supplies an arbitrary positive conformal factor.
    No vacuum equation or Hamiltonian constraint is imposed on the generic tensors.
    """
    import numpy as np

    X = chart.symbols
    if system == "cartesian":
        image = list(X)
    elif system == "spherical":
        r, t, p = X
        image = [r * sp.sin(t) * sp.cos(p), r * sp.sin(t) * sp.sin(p), r * sp.cos(t)]
    else:
        rho, z, p = X
        image = [rho * sp.cos(p), rho * sp.sin(p), z]
    x, y, z = image
    # Axisymmetry is enough for all three charts, with angular dependence retained.
    f = 2 + (x*x + y*y) / 20 + z / 30 + z*z / 40
    psi = chart.reader.parameters["psi"]
    scalar = chart.geo.ricci_scalar().subs(psi, f).doit()
    expected = -2 / f**5  # -8 times the Cartesian Laplacian, 1/4.
    J = sp.Matrix(image).jacobian(X)
    # Each finite Jacobian is used directly, with no symbolic determinant or simplification.
    metric = chart.geo.g.subs(psi, f)
    at_scalar = sp.lambdify(X, scalar - expected, "numpy")
    at_metric = sp.lambdify(X, metric, "numpy")
    at_J = sp.lambdify(X, J, "numpy")
    at_f = sp.lambdify(X, f, "numpy")
    for point in ((1.3, 0.7, 0.4), (2.2, 1.1, 1.8), (0.8, 2.1, 2.4)):
        if abs(float(at_scalar(*point))) > 1e-12:
            raise AssertionError(f"bowen_york/{system}: the conformal scalar misses at {point}")
        jac = np.asarray(at_J(*point), dtype=float)
        pulled = float(at_f(*point))**4 * (jac.T @ jac)
        if np.max(np.abs(np.asarray(at_metric(*point), dtype=float) - pulled)) > 1e-10:
            raise AssertionError(f"bowen_york/{system}: the Cartesian pullback misses at {point}")


def embedding(ck, src):
    """Two copies of the equatorial exterior joined at the inversion throat."""
    import numpy as np
    import embedding as em
    from bowen_york_data import solution

    w = solution()
    sl = em.Slice(src, "bowen_york", "spherical", "r", r"\phi", {"theta": "pi/2"},
                  numeric={"psi": (lambda r: w.equator(r)[0], lambda r: w.equator(r)[1])})
    reach, size = 8.0, 24.0
    radicand = sl.defect_at(np.linspace(1, reach, 1001))
    ck.exact("Bowen-York: the equatorial embedding radicand is nonnegative", bool(np.min(radicand) >= -1e-12))
    marks = [(1.0, "throat", None), (2.0, "r", None), (4.0, "r", None), (reach, "r", None)]
    pieces = []
    for name, sense in (("outer_end", 1), ("other_end", -1)):
        piece = em.Piece(name, "sheet", sl, 1.0, reach, 0.0, sense,
                         (("throat", "the inversion throat $r=a$"),
                          ("edge", "an asymptotically flat end")), marks, size)
        ck.isometry(f"Bowen-York, {name}", piece)
        pieces.append(piece)
    surface = em.Surface(pieces)
    fig = em.figure_of([surface], {"sheet": "cover"}, size)
    fig.legend("fill", "cover", "equatorial sheets")
    fig.legend("line", "throat", "inversion throat")
    fig.legend("line", "r", "$r/a=2,4,8$")
    fig.legend("line", "meridian", "$\\phi$ constant")
    # A square enclosing the entire rotating surface avoids angle-dependent narrow framing.
    result = fig.done()
    x0, x1, y0, y1 = result["box"]
    width = max(x1-x0, y1-y0, 2*reach*w.field(reach, 0)**2) * 1.15
    cx, cy = (x0+x1)/2, (y0+y1)/2
    result["box"] = [cx-width/2, cx+width/2, cy-width/2, cy+width/2]
    return [em.view("spin", "Spinning throat", "$a$", [surface], result,
                    caption=["Two isometric equatorial sheets of a spinning throat."],
                    settings="$a=1$, $S=a^2$, $P=0$.",
                    input="Numerical solution of the axisymmetric Hamiltonian constraint.")]


def curvature(ck, src):
    """Evaluate the stored intrinsic Ricci scalar on the numerical spatial slice."""
    import math
    import numpy as np
    import curvature as cv
    from bowen_york_data import solution

    functions = {"psi": solution().formula()}
    scalar = cv.Published(src, "bowen_york", "spherical", "ricci_scalar", ["r", "theta"], [], functions=functions)
    xaxis = cv.axis("$r\\sin\\theta/a$", 0, 4)
    yaxis = cv.axis("$r\\cos\\theta/a$", -4, 4)
    nx, ny = 80, 160
    xs = (np.arange(nx)+0.5)*4/nx
    ys = -4+(np.arange(ny)+0.5)*8/ny
    X, Y = np.meshgrid(xs, ys)
    r, theta = np.hypot(X,Y), np.arctan2(X,Y)
    inside = r >= 1
    values = scalar(np.where(inside, r, 1), theta)
    cells = [[cv.rounded(values[j,i]) if inside[j,i] else None for i in range(nx)] for j in range(ny)]
    scale = {"kind": "linear", "from": 0, "to": 0.005}
    throat = [[float(math.sin(t)), float(math.cos(t))] for t in np.linspace(0, math.pi, 181)]
    return [{
        "id": "spin", "label": "Spinning throat", "kind": "map", "x": xaxis, "y": yaxis,
        "colour": {"name": "$R\\,a^2$", "scale": scale, "ticks": cv.value_ticks(scale)},
        "map": {"metric": "bowen_york", "system": "spherical", "scalar": "ricci_scalar",
                "x": "r", "y": r"\theta", "frame": "meridian", "params": {}, "functions": functions,
                "domain": "r >= 1", "nx": nx, "ny": ny, "cells": cells},
        "marks": [{"class": "horizon", "points": throat}],
        "legend": [["line", "horizon", "inversion throat"]],
        "caption": ["Intrinsic scalar curvature outside a spinning throat."],
        "settings": "$a=1$, $S=a^2$, $P=0$.",
        "input": "Numerical solution of the axisymmetric Hamiltonian constraint.",
        "clip": "Linear, $0\\le R a^2\\le0.005$.",
    }]
