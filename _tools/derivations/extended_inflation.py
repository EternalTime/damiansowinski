"""The flat Jordan false vacuum branch of Kolitch and Eardley (5.20)-(5.21).

The coordinate time in the tensor calculation is c times the published time.
The held scale factor has its exact power-law derivative in both charts.
"""
import sympy as sp

import verify_metrics as vm


CHARTS = ("cosmic", "conformal")
SCALE = {
    "cosmic": r"a = \left(\dfrac{t}{t_0}\right)^p",
    "conformal": r"a = \left(-\dfrac{(p-1)\eta}{t_0}\right)^{-p/(p-1)}",
}
DIMENSIONS = {
    ("extended_inflation", system): {
        ("t" if system == "cosmic" else r"\eta"): "T",
        "r": "L", r"\theta": "1", r"\phi": "1", "p": "1", "t_0": "T", "a": "1",
    } for system in CHARTS
}
RATES = {
    ("extended_inflation", "cosmic"): {"a": {"t": "pa/t"}},
    ("extended_inflation", "conformal"): {"a": {r"\eta": r"-pa/((p-1)\eta)"}},
}


def chart_spec(system):
    time = "t" if system == "cosmic" else r"\eta"
    coords = [time, "r", r"\theta", r"\phi"]
    parameters = ["p", "t_0", SCALE[system]]
    internal = ["p", "t_0", SCALE[system].replace("{t_0}", "{ct_0}")]
    probe = vm.Reader(coords, internal, (), held=("a",), rates=RATES[("extended_inflation", system)])
    space = r"dr^2 + r^2d\theta^2 + r^2\sin^2\theta\,d\phi^2"
    line = (rf"ds^2 = -c^2dt^2 + a^2\left({space}\right)" if system == "cosmic"
            else rf"ds^2 = a^2\left(-c^2d\eta^2 + {space}\right)")
    return {
        "metric_id": "extended_inflation",
        "system": {"id": system, "name": "Cosmic Spherical" if system == "cosmic" else "Conformal Spherical",
                   "coords": coords, "parameters": parameters,
                   "domains": [rf"{time} \in " + (r"(0, \infty)" if system == "cosmic" else r"(-\infty, 0)"),
                               r"r \in [0, \infty)", r"\theta \in [0, \pi]", r"\phi \in [0, 2\pi)",
                               r"p > 1,\quad t_0 > 0"], "line_element": line},
        "chart_parameters": internal,
        "chart_line_element": line.replace("c^2dt", "dt").replace(r"c^2d\eta", r"d\eta"),
        "printer": {}, "time": time, "reduce": probe.surface,
        "check": lambda chart: check_chart(chart, system),
    }


def check_chart(chart, system):
    """Pointwise scalar-tensor equations and the cosmic-to-conformal pullback.

    The positive Jordan vacuum energy is constant in the action.
    No symbolic determinant or simplify call is used here.
    """
    X, P, geo = chart.symbols, chart.reader.parameters, chart.geo
    x, r, theta, _ = X
    p, t0, c = P["p"], P["t_0"], chart.reader.c
    omega = p - sp.Rational(1, 2)
    lam = (6*omega+5)*(2*omega+3)/(4*c**2*t0**2)
    tau = x if system == "cosmic" else c*t0*(-(p-1)*x/(c*t0))**(1/(1-p))
    Phi = (tau/(c*t0))**2
    gradient = sp.Matrix([sp.diff(Phi, v) for v in X])
    gamma = geo.christoffel_ull()
    hessian = sp.Matrix(4, 4, lambda i, j: sp.diff(Phi, X[i], X[j])
                        - sum(vm._at(gamma, (k, i, j))*gradient[k] for k in range(4)))
    box = sum(geo.ginv[i, j]*hessian[i, j] for i in range(4) for j in range(4))
    norm = (gradient.T*geo.ginv*gradient)[0]
    source = (-lam*geo.g/Phi + omega*(gradient*gradient.T-geo.g*norm/2)/Phi**2
              + (hessian-geo.g*box)/Phi)
    einstein = sp.Matrix(geo.einstein_ll())
    held = chart.reader.held
    for power, reference, cosmic, radius in ((2, 1.2, 0.8, 0.7), (3, 0.9, 1.8, 1.4), (sp.Rational(5, 4), 1.7, 3.1, 2.2)):
        light = 1.3
        conformal = -reference*(cosmic/reference)**(1-float(power))/(float(power)-1)
        at = {p: power, t0: reference, c: light, x: light*(cosmic if system == "cosmic" else conformal),
              r: radius, theta: 1.1, X[3]: 0.4}

        def number(value):
            return complex(sp.N(sp.sympify(value).subs(held).doit().subs(at), 35))

        def near(values, label):
            if max((abs(number(v)) for v in sp.flatten(values)), default=0) > 2e-8:
                raise AssertionError(f"extended_inflation/{system}: {label}")

        near(list(einstein-source), "Jordan false vacuum field equation")
        near([box+4*lam/(2*omega+3)], "scalar wave equation and normalization")
        # Cosmic x^0 as a function of conformal x^0 has derivative a.
        scale = (cosmic/reference)**float(power)
        pulled = sp.diag(-1, scale**2, scale**2*radius**2, scale**2*radius**2*sp.sin(1.1)**2)
        if system == "conformal":
            J = sp.diag(scale, 1, 1, 1)
            pulled = J.T*pulled*J
        near(list(geo.g-pulled), "cosmic metric pulled back")

