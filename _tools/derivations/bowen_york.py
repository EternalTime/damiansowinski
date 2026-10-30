"""Bowen-York maximal initial data, with intrinsic three-dimensional geometry only.

The free conformal factor is retained in every printed tensor.
Source equations and the two boundary prescriptions are recorded in bowen_york.md.
"""
import sympy as sp

CHARTS = ("cartesian", "spherical", "cylindrical")


def chart_spec(system):
    """Cartesian and spherical charts, and the axisymmetric cylindrical pullback."""
    import chart_printer as cp
    import verify_metrics as vm

    if system == "cartesian":
        coords, function = ["x", "y", "z"], r"\psi = \psi(x,y,z)"
        flat = "dx^2 + dy^2 + dz^2"
        lap = r"\partial_x^2\psi + \partial_y^2\psi + \partial_z^2\psi"
    elif system == "spherical":
        coords, function = ["r", r"\theta", r"\phi"], r"\psi = \psi(r,\theta,\phi)"
        flat = r"dr^2 + r^2\left(d\theta^2 + \sin^2\theta\,d\phi^2\right)"
        lap = (r"\partial_r^2\psi + \frac{2\partial_r\psi}{r}"
               r" + \frac{\partial_\theta^2\psi + \cot\theta\,\partial_\theta\psi}{r^2}"
               r" + \frac{\partial_\phi^2\psi}{r^2\sin^2\theta}")
    elif system == "cylindrical":
        coords, function = [r"\rho", "z", r"\phi"], r"\psi = \psi(\rho,z)"
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
    expected = -sp.Rational(6, 5) / f**5  # -8 times the Cartesian Laplacian, 3/20.
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
