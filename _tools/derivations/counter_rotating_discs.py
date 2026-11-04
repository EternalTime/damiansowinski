"""Reflected Weyl fields and their material junction, BLK (1993), BL (1993), LB (2019)."""
import math

import numpy as np
import sympy as sp

METRIC = "counter_rotating_discs"
SYSTEMS = ("static", "stationary")


def chart_spec(pc, system):
    coords = ["t", r"\rho", r"\phi", "z"]
    if system == "static":
        parameters = [r"U = U(\rho,z)", "k = k(\\rho,z)", "b", "m", r"\delta", "a"]
        line = (r"ds^2 = -e^{2U}dt^2 + e^{2k-2U}(d\rho^2 + dz^2)"
                r" + \rho^2e^{-2U}d\phi^2")
    else:
        parameters = [r"f = f(\rho,z)", r"g = g(\rho,z)", r"A = A(\rho,z)", "b", "m", "a"]
        line = (r"ds^2 = -f(dt + A\,d\phi)^2 + \dfrac{g}{f}(d\rho^2 + dz^2)"
                r" + \dfrac{\rho^2}{f}d\phi^2")
    probe = pc.vm.Reader(coords, parameters, ())
    rho, z = probe.symbol[r"\rho"], probe.symbol["z"]
    funcs = [probe.parameters[n] for n in (("U", "k") if system == "static" else ("f", "g", "A"))]
    lead = funcs + [sp.diff(f, q) for f in funcs for q in (rho, z)]
    lead += [sp.diff(f, q, r) for f in funcs for q, r in ((rho, rho), (rho, z), (z, z))]
    return {
        "metric_id": METRIC,
        "system": {"id": system, "name": "Reflected Static Weyl" if system == "static" else "Reflected Stationary Weyl",
                   "coords": coords, "parameters": parameters,
                   "domains": [r"t \in (-\infty, \infty)", r"\rho \in [0, \infty)",
                               r"\phi \in [0, 2\pi)", r"z \in (-\infty, \infty)", "b > 0",
                               r"z = 0\;\text{(the matter junction)}"], "line_element": line},
        "chart_line_element": line,
        "printer": {"lead": lead + [rho], "factors": lead + [rho]},
        "pretty": lambda v: sp.powsimp(sp.factor(sp.sympify(v)), combine="exp"),
        "check": lambda chart: check_chart(chart, system),
    }


def check_chart(chart, system):
    """Compare the printed chart with the source line element at nonsingular points."""
    R = chart.reader
    rho, z = R.symbol[r"\rho"], R.symbol["z"]
    for rv, zv in ((0.4, 0.7), (1.8, 1.2), (4.2, -0.6)):
        u, k = -1 / math.hypot(rv, 3 + abs(zv)), -rv * rv / (2 * (rv * rv + (3 + abs(zv)) ** 2) ** 2)
        if system == "static":
            values = {R.parameters["U"]: u, R.parameters["k"]: k}
            target = np.diag([-math.exp(2*u), math.exp(2*k-2*u), rv*rv*math.exp(-2*u), math.exp(2*k-2*u)])
        else:
            f, g, A = math.exp(2*u), math.exp(2*k), 0.13
            values = {R.parameters["f"]: f, R.parameters["g"]: g, R.parameters["A"]: A}
            target = np.diag([-f, g/f, rv*rv/f-f*A*A, g/f])
            target[0, 2] = target[2, 0] = -f*A
        values.update({rho: rv, z: zv})
        got = np.array(chart.geo.g.subs(values), dtype=float)
        if not np.allclose(got, target, atol=1e-12, rtol=1e-12):
            raise AssertionError(f"{METRIC}/{system}: source line element at {(rv, zv)}")


def curzon_functions(side=1):
    """Smooth extension of either retained half; never differentiate abs(z) at the seam."""
    Z = f"(b+({side})*z)"
    return {"U": f"-m/sqrt(rho**2+{Z}**2)",
            "k": f"-m**2*rho**2/(2*(rho**2+{Z}**2)**2)"}


def static_matter(rho, b=3.0, m=1.0):
    radius = math.hypot(rho, b)
    u, k = -m/radius, -m*m*rho*rho/(2*radius**4)
    ur, uz = m*rho/radius**3, m*b/radius**3
    weight = math.exp(u-k)*uz/(2*math.pi)
    energy, pressure = weight*(1-rho*ur), weight*rho*ur
    return energy, pressure, math.sqrt(pressure/energy), (energy-pressure)/2
