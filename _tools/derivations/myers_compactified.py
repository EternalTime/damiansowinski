"""Charts of the extremal black hole with one compact spatial direction."""

import sympy as sp

METRIC = "myers_compactified"
CHARTS = ("cylindrical", "spherical")
DIMENSIONS = {
    (METRIC, "cylindrical"): {"t": "L", "r": "L", r"\theta": "1", r"\phi": "1", "z": "L",
                               "U": "1", r"\mu": "L**2", "L": "L"},
    (METRIC, "spherical"): {"t": "L", r"\rho": "L", r"\chi": "1", r"\theta": "1", r"\phi": "1",
                             "U": "1", r"\mu": "L**2", "L": "L"},
}


def chart_spec(system):
    if system == "cylindrical":
        coords = ["t", "r", r"\theta", r"\phi", "z"]
        params = ["U = U(r,z)", r"\mu", "L"]
        space = r"dr^2+r^2\left(d\theta^2+\sin^2\theta\,d\phi^2\right)+dz^2"
        domains = [r"t \in (-\infty,\infty)", r"r \in (0,\infty)", r"\theta \in [0,\pi]",
                   r"\phi \in [0,2\pi)", r"z \in [-L/2,L/2)",
                   r"r=z=0 \;\text{(the horizon)}"]
        name = "Cylindrical"
    else:
        coords = ["t", r"\rho", r"\chi", r"\theta", r"\phi"]
        params = [r"U = U(\rho,\chi)", r"\mu", "L"]
        space = r"d\rho^2+\rho^2\left(d\chi^2+\sin^2\chi\left(d\theta^2+\sin^2\theta\,d\phi^2\right)\right)"
        domains = [r"t \in (-\infty,\infty)", r"\rho \in (0,\infty)", r"\chi \in [0,\pi]",
                   r"\theta \in [0,\pi]", r"\phi \in [0,2\pi)",
                   r"-L/2 < \rho\cos\chi < L/2", r"\rho=0 \;\text{(the horizon)}"]
        name = "Spherical"
    line = r"ds^2=-\dfrac{dt^2}{U^2}+U\left(" + space + r"\right)"
    return {"metric_id": METRIC,
            "system": {"id": system, "name": name, "coords": coords, "parameters": params,
                       "domains": domains, "line_element": line},
            "chart_line_element": line, "printer": {},
            "check": lambda chart: check_chart(chart, system)}


def check_chart(chart, system):
    """Check the sourced metric, harmonic equation, and Maxwell stress at points."""
    t, a, b, c, d = chart.symbols
    U = chart.reader.parameters["U"]
    mu, length = chart.reader.parameters["mu"], chart.reader.parameters["L"]
    r, z = (a, d) if system == "cylindrical" else (a*sp.sin(b), a*sp.cos(b))
    H = 1 + sp.pi*mu/(length*r)*sp.sinh(2*sp.pi*r/length)/(sp.cosh(2*sp.pi*r/length)-sp.cos(2*sp.pi*z/length))
    g, gi = chart.geo.g, chart.geo.ginv
    field = sp.Matrix(5, 5, lambda i, j: sp.diff(1/U, chart.symbols[i]) * int(j == 0)
                      - sp.diff(1/U, chart.symbols[j]) * int(i == 0))
    raised = gi*field*gi
    field2 = sum(field[i, j]*raised[i, j] for i in range(5) for j in range(5))
    stress = sp.Matrix(5, 5, lambda i, j: sp.Rational(3, 2) *
                       (sum(field[i, k]*field[j, l]*gi[k, l] for k in range(5) for l in range(5))
                        - g[i, j]*field2/4))
    einstein = sp.Matrix(chart.geo.einstein_ll())
    if system == "cylindrical":
        laplace = sp.diff(H, a, 2) + 2*sp.diff(H, a)/a + sp.diff(H, d, 2)
    else:
        laplace = sp.diff(H, a, 2) + 3*sp.diff(H, a)/a + (sp.diff(H, b, 2)
                    + 2*sp.cos(b)*sp.diff(H, b)/sp.sin(b))/a**2
    expressions = list(einstein-stress)
    derivatives = set().union(*(expr.atoms(sp.Derivative) for expr in expressions))
    for scale, radius, angle in ((sp.Rational(3, 4), sp.Rational(7, 10), sp.Rational(2, 5)),
                                  (sp.Rational(5, 4), sp.Rational(11, 10), sp.Rational(4, 5))):
        at = {t: 0, a: radius, b: angle, c: sp.Rational(7, 10), d: sp.Rational(1, 3),
              mu: scale, length: 3}
        vals = {U: H.evalf(35, subs=at)}
        vals.update({term: sp.diff(H, *term.variables).evalf(35, subs=at) for term in derivatives})
        if abs(complex(laplace.evalf(35, subs=at))) > 1e-25:
            raise AssertionError(f"{METRIC}/{system}: the sourced potential is not harmonic")
        for expr in expressions:
            if abs(complex(expr.xreplace(vals).evalf(30, subs=at))) > 1e-20:
                raise AssertionError(f"{METRIC}/{system}: Einstein-Maxwell equation")
