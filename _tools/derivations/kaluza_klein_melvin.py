"""The magnetic universe of Kaluza and Klein: five dimensional Minkowski space identified by a
step of 2 pi R around the fifth dimension together with a turn of 2 pi R B about an axis.

The adapted chart is Dowker, Gauntlett, Giddings and Horowitz's (5.6), hep-th/9312172, which is
Dowker, Gauntlett, Gibbons and Horowitz's (2.3), hep-th/9507143, with their x^5 written y and
their tilde phi written varphi; its two angles are identified independently, phi with period
2 pi and y with period 2 pi R.  The Minkowski chart is the 1995 (2.1) with the identifications
(2.2), phi = varphi + B y.  Both charts are flat.  kaluza_klein_melvin.md records the derivation.
"""

METRIC = "kaluza_klein_melvin"
CHARTS = ("adapted", "minkowski")
ANGLES = {"adapted": r"\varphi", "minkowski": r"\phi"}
DIMENSIONS = {(METRIC, system): {"t": "T", "z": "L", r"\rho": "L", ANGLES[system]: "1", "y": "L",
                                 "R": "L", "B": "L**-1"}
              for system in CHARTS}

REALS = r"(-\infty, \infty)"
AXIS = r"\rho = 0 \;\text{(the axis, a regular line)}"
RANGE = r"B \in (-1/2R, 1/2R] \;\text{(distinct spacetimes)}"


def line_element(system, time=r"c^2dt^2"):
    """The line element of either chart, with the time term given."""
    if system == "adapted":
        return (r"ds^2 = -" + time + r" + dz^2 + d\rho^2 + \rho^2\left(d\varphi + B\,dy\right)^2 + dy^2")
    if system == "minkowski":
        return r"ds^2 = -" + time + r" + dz^2 + d\rho^2 + \rho^2d\phi^2 + dy^2"
    raise ValueError(system)


def chart_spec(system):
    """One chart of the twisted space, its line element as its source prints it with c restored."""
    angle = ANGLES[system]
    domains = ["t \\in " + REALS, "z \\in " + REALS, r"\rho \in [0, \infty)",
               angle + r" \in [0, 2\pi)", r"y \in [0, 2\pi R)"]
    if system == "minkowski":
        domains.append(r"(\phi, y) \sim (\phi + 2\pi RB, y + 2\pi R)")
    domains += [RANGE, AXIS]
    name = {"adapted": "Adapted Angle", "minkowski": "Minkowski Angle"}[system]
    return {"metric_id": METRIC,
            "system": {"id": system, "name": name, "coords": ["t", "z", r"\rho", angle, "y"],
                       "domains": domains, "parameters": ["R", "B"], "line_element": line_element(system)},
            "chart_line_element": line_element(system, "dt^2"),
            "printer": {},
            "check": lambda chart: check_chart(chart, system)}


POINTS = ({"t": (1, 3), "z": (-2, 5), "rho": (7, 4), "angle": (2, 9), "y": (5, 3), "B": (1, 4)},
          {"t": (-3, 2), "z": (1, 7), "rho": (2, 5), "angle": (11, 7), "y": (-4, 9), "B": (-2, 5)})


def check_chart(chart, system):
    """At rational points: the adapted chart is the Minkowski chart pulled back by
    phi = varphi + B y, and every Riemann component of either chart vanishes."""
    import itertools

    import sympy as sp
    import verify_metrics as vm
    t, z, rho, angle, y = chart.symbols
    B = chart.reader.parameters["B"]
    minkowski = sp.diag(-1, 1, 1, rho**2, 1)
    if system == "adapted":
        jacobian = sp.eye(5)
        jacobian[3, 4] = B
        expected = jacobian.T * minkowski * jacobian
    else:
        expected = minkowski
    riemann = chart.geo.riemann_ulll()
    for point in POINTS:
        at = {t: sp.Rational(*point["t"]), z: sp.Rational(*point["z"]), rho: sp.Rational(*point["rho"]),
              angle: sp.Rational(*point["angle"]), y: sp.Rational(*point["y"]), B: sp.Rational(*point["B"])}
        for value in chart.geo.g - expected:
            if value.subs(at) != 0:
                raise AssertionError(f"{METRIC}/{system}: the metric is not (2.1) pulled back by phi = varphi + By")
        for index in itertools.product(range(5), repeat=4):
            if vm._at(riemann, index).subs(at) != 0:
                raise AssertionError(f"{METRIC}/{system}: Riemann component {index} does not vanish")
