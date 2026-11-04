"""Reflected Weyl fields and their material junction, BLK (1993), BL (1993), LB (2019)."""
import math

import numpy as np
import sympy as sp

METRIC = "counter_rotating_discs"
SYSTEMS = ("static", "stationary")
DIMENSIONS = {
    (METRIC, "static"): {"t": "L", r"\rho": "L", r"\phi": "1", "z": "L", "U": "1", "k": "1",
                        "b": "L", "m": "L", r"\delta": "1", "a": "L"},
    (METRIC, "stationary"): {"t": "L", r"\rho": "L", r"\phi": "1", "z": "L", "f": "1", "g": "1",
                            "A": "L", "b": "L", "m": "L", "a": "L"},
}


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


def seed_matrix(kind, rho, z, m=1.0, a=0.8, b=3.0, delta=1.0):
    """Analytic metric on the upper retained half, in (t,rho,phi,z); complex-step compatible."""
    Z = b + z
    matrix = np.zeros((4, 4), dtype=np.result_type(rho, z, float))
    if kind in ("curzon", "rod"):
        if kind == "curzon":
            R = np.sqrt(rho*rho + Z*Z)
            U, k = -m/R, -m*m*rho*rho/(2*R**4)
        else:
            plus, minus = np.sqrt(rho*rho+(Z+a)**2), np.sqrt(rho*rho+(Z-a)**2)
            total = plus + minus
            U = delta/2*np.log((total-2*a)/(total+2*a))
            k = delta*delta/2*np.log((total*total-4*a*a)/(4*plus*minus))
        matrix[0, 0] = -np.exp(2*U)
        matrix[1, 1] = matrix[3, 3] = np.exp(2*k-2*U)
        matrix[2, 2] = rho*rho*np.exp(-2*U)
    else:
        sigma = m*m-a*a
        s = rho*rho + Z*Z + sigma
        u = (s + np.sqrt(s*s - 4*sigma*Z*Z))/2
        r, y = m + np.sqrt(u), Z/np.sqrt(u)
        Sigma = r*r+a*a*y*y
        d = Sigma-2*m*r
        f = d/Sigma
        A = 2*m*a*r*(1-y*y)/d
        matrix[0, 0], matrix[0, 2], matrix[2, 0] = -f, -f*A, -f*A
        matrix[1, 1] = matrix[3, 3] = Sigma/(u-sigma*y*y)
        matrix[2, 2] = rho*rho/f-f*A*A
    return matrix


def seam_stress(kind, rho, **parameters):
    """Israel stress with one normal pointing from the lower half to the upper half."""
    h = seed_matrix(kind, rho, 0.0, **parameters)
    derivative = seed_matrix(kind, rho, 1j*1e-25, **parameters).imag/1e-25
    induced = h[:3, :3]
    jump = derivative[:3, :3]/np.sqrt(h[3, 3])
    trace = np.trace(np.linalg.solve(induced, jump))
    stress = -(jump-induced*trace)/(8*np.pi)
    return h, stress
