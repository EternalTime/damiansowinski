"""The closed cosmology of Nappi and Witten, hep-th/9206078 (16), (21).

The metric is the large-level string-frame geometry, with its common
length scale restored.  The 1993 plane wave is a separate spacetime.
"""

METRIC = "nappi_witten"
CHARTS = ("cosmological", "bang_regular", "invariant", "whisker")
CLOCKS = {"cosmological": r"\psi", "bang_regular": "x",
          "invariant": "r", "whisker": r"\chi"}
DIMENSIONS = {(METRIC, system): {CLOCKS[system]: "1", "s": "1", r"\rho": "1",
                               ("w" if system == "bang_regular" else r"\lambda"): "1",
                               "L": "L", "b": "1", "D": "1"}
              for system in CHARTS}


def chart_spec(system):
    """Four coordinate forms of the same classical coset geometry."""
    angular = "w" if system == "bang_regular" else r"\lambda"
    coords = [CLOCKS[system], "s", r"\rho", angular]
    common = [r"s \in (0,\pi/2)", r"\rho \in [0,2\pi)",
              angular+r" \in [0,2\pi)", "L>0", "b>0"]
    if system == "cosmological":
        denominator = r"\sin^2s\cos^2\psi+b^2\cos^2s\sin^2\psi"
        line = (r"ds^2=L^2\left(-d\psi^2+ds^2+\dfrac{\cos^2s\cos^2\psi}{D}d\rho^2"
                r"+\dfrac{b^2\sin^2s\sin^2\psi}{D}d\lambda^2\right)")
        domains, name = [r"\psi \in (0,\pi/2)"]+common, "Cosmological"
    elif system == "bang_regular":
        denominator = r"\sin^2s+b^2x\cos^2s"
        line = (r"ds^2=L^2\left(ds^2+\dfrac{(2+x)\sin^2s-b^2\cos^2s}{4(1+x)^2D}dx^2"
                r"+\dfrac{\cos^2s}{D}d\rho^2+\dfrac{b^2x\sin^2s}{D}dw^2"
                r"+\dfrac{b\sin^2s}{D}dx\,dw\right)")
        domains, name = [r"x \in (-1,\infty)", r"D\neq0"]+common, "Across the Big Bang"
    elif system == "invariant":
        denominator = r"(1-r)\sin^2s+b^2r\cos^2s"
        line = (r"ds^2=L^2\left(-\dfrac{dr^2}{4r(1-r)}+ds^2"
                r"+\dfrac{(1-r)\cos^2s}{D}d\rho^2+\dfrac{b^2r\sin^2s}{D}d\lambda^2\right)")
        domains, name = [r"r \in (-\infty,\infty)", r"r\neq0,1", r"D\neq0"]+common, "Invariant Coordinate"
    elif system == "whisker":
        denominator = r"\sin^2s\cosh^2\chi-b^2\cos^2s\sinh^2\chi"
        line = (r"ds^2=L^2\left(d\chi^2+ds^2+\dfrac{\cos^2s\cosh^2\chi}{D}d\rho^2"
                r"-\dfrac{b^2\sin^2s\sinh^2\chi}{D}d\lambda^2\right)")
        domains, name = [r"\chi \in (0,\infty)", r"D\neq0"]+common, "Static Extension"
    else:
        raise ValueError(system)
    parameters = ["L", "b", "D = "+denominator]
    return {"metric_id": METRIC, "system": {"id": system, "name": name,
            "coords": coords, "domains": domains, "parameters": parameters,
            "line_element": line}, "chart_line_element": line, "printer": {},
            "check": lambda chart: check_chart(chart, system)}


def background_fields(chart, system):
    """Dilaton gradient and a gauge potential for the sourced three-form.

    The potential differs from 1992 (14) by a closed two-form.
    Its derivative is unchanged, including in the regular gauge.
    """
    import sympy as sp
    u, s, rho, angle = chart.symbols
    L, b, denominator = (chart.reader.parameters[key] for key in ("L", "b", "D"))
    effective = denominator/(1+u) if system == "bang_regular" else denominator
    if system == "cosmological":
        numerator = sp.sin(s)**2*sp.cos(u)**2
    elif system == "bang_regular":
        numerator = sp.sin(s)**2
    elif system == "invariant":
        numerator = (1-u)*sp.sin(s)**2
    else:
        numerator = sp.sin(s)**2*sp.cosh(u)**2
    potential = sp.zeros(4)
    potential[2, 3] = L**2*numerator/denominator
    potential[3, 2] = -potential[2, 3]
    if system == "bang_regular":
        # Subtract L^2 d(rho) wedge d(log x)/(2b) to cross x=0 smoothly.
        potential[2, 0] = -L**2*b*sp.cos(s)**2/(2*denominator)
        potential[0, 2] = -potential[2, 0]
    gradient = [-sp.diff(effective, x)/(2*effective) for x in chart.symbols]
    three_form = {(i, j, k): sp.diff(potential[j, k], chart.symbols[i])
                  + sp.diff(potential[k, i], chart.symbols[j])
                  + sp.diff(potential[i, j], chart.symbols[k])
                  for i in range(4) for j in range(4) for k in range(4)}
    return gradient, three_form


def check_chart(chart, system):
    """Check (16), the coordinate pullbacks, and the metric beta equation at points.

    No symbolic determinant or simplification enters these comparisons.
    The geometry has a dilaton and a three-form: R + 2 Hess(Phi) - H^2/4 = 0.
    """
    import sympy as sp
    import verify_metrics as vm
    u, s, rho, angle = chart.symbols
    L, b = (chart.reader.parameters[key] for key in ("L", "b"))
    q = {"cosmological": sp.tan(u)**2, "bang_regular": u,
         "invariant": u/(1-u), "whisker": -sp.tanh(u)**2}[system]
    alpha_sine = (1-b*b)/(1+b*b)
    twice_time_cosine = (1-q)/(1+q)
    source_denominator = (1-twice_time_cosine*sp.cos(2*s)
                          + alpha_sine*(twice_time_cosine-sp.cos(2*s)))
    source = L**2*sp.diag(-1, 1,
        2*sp.cos(s)**2*(1+alpha_sine)/((1+q)*source_denominator),
        2*sp.sin(s)**2*q*(1-alpha_sine)/((1+q)*source_denominator))
    if system == "bang_regular":
        jacobian = sp.eye(4)
        jacobian[0, 0] = 1/(2*sp.sqrt(u)*(1+u))
        jacobian[3, 0] = 1/(2*b*u)
        source = jacobian.T*source*jacobian
    elif system == "invariant":
        jacobian = sp.diag(1/(2*sp.sqrt(u)*sp.sqrt(1-u)), 1, 1, 1)
        source = jacobian.T*source*jacobian
    elif system == "whisker":
        source[0, 0] = L**2
    gradient, three_form = background_fields(chart, system)
    gamma, ricci = chart.geo.christoffel_ull(), chart.geo.ricci_ll()
    hessian = sp.Matrix(4, 4, lambda i, j: sp.diff(gradient[j], chart.symbols[i])
                         - sum(vm._at(gamma, (k, i, j))*gradient[k] for k in range(4)))
    for time, latitude, deformation in ((sp.Rational(2, 5), sp.Rational(7, 10), sp.Rational(3, 5)),
                                         (sp.Rational(3, 5), sp.Rational(4, 5), sp.Rational(7, 5))):
        at = {u: time, s: latitude, b: deformation, L: sp.Rational(6, 5),
              rho: sp.Rational(1, 7), angle: sp.Rational(2, 9)}
        for value in chart.geo.g-source:
            if abs(sp.N(value.subs(at), 40)) > 1e-27:
                raise AssertionError(f"{METRIC}/{system}: sourced coordinate pullback")
        inverse = chart.geo.ginv.subs(at).evalf(40)
        forms = [sp.Matrix(4, 4, lambda i, j: three_form[k, i, j].subs(at).evalf(40))
                 for k in range(4)]
        raised = [inverse*h*inverse.T for h in forms]
        for i in range(4):
            for j in range(i, 4):
                stress = sum(forms[i][a, c]*raised[j][a, c] for a in range(4) for c in range(4))/4
                residual = sp.N((vm._at(ricci, (i, j))+2*hessian[i, j]).subs(at), 40)-stress
                if abs(residual) > 1e-24:
                    raise AssertionError(f"{METRIC}/{system}: metric beta equation ({i},{j}): {residual}")
