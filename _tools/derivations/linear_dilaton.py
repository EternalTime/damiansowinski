"""The Einstein-frame charts of the timelike linear dilaton vacuum.

De Vega and Egusquiza (1994), equation (3.2), gives the cosmic and
conformal forms.  Angius, Delgado and Uranga (2022), section 3.2, gives
the frame transformation and the relation between their time coordinates.
"""

METRIC = "linear_dilaton"
CHARTS = ("cosmic", "conformal")
DIMENSIONS = {(METRIC, system): {clock: "L", "x": "L", "y": "L", "z": "L", "q": "1/L"}
              for system, clock in (("cosmic", r"\tau"), ("conformal", r"\eta"))}


def chart_spec(system):
    """Two Cartesian charts of the same Einstein-frame metric."""
    clock = r"\tau" if system == "cosmic" else r"\eta"
    coords = [clock, "x", "y", "z"]
    spatial = r"\left(dx^2+dy^2+dz^2\right)"
    line = (r"ds^2 = -d\tau^2+q^2\tau^2"+spatial if system == "cosmic" else
            r"ds^2 = e^{2q\eta}\left(-d\eta^2+dx^2+dy^2+dz^2\right)")
    domains = [clock + (r" \in (0,\infty)" if system == "cosmic" else r" \in (-\infty,\infty)")]
    domains += [name+r" \in (-\infty,\infty)" for name in ("x", "y", "z")]
    return {"metric_id": METRIC, "system": {"id": system,
            "name": "Cosmic Time" if system == "cosmic" else "Conformal Time",
            "coords": coords, "domains": domains, "parameters": ["q"], "line_element": line},
            "chart_line_element": line, "printer": {},
            "check": lambda chart: check_chart(chart, system)}


def check_chart(chart, system):
    """Check the sourced metric and scalar curvature at independent points."""
    import sympy as sp
    t, x, y, z = chart.symbols
    q = chart.reader.parameters["q"]
    factor = q*q*t*t if system == "cosmic" else sp.exp(2*q*t)
    expected = sp.diag(-1 if system == "cosmic" else -factor, factor, factor, factor)
    scalar = 6/t**2 if system == "cosmic" else 6*q*q/sp.exp(2*q*t)
    for at in ({t: sp.Rational(1, 2), q: sp.Rational(2, 3)},
               {t: sp.Rational(3, 2), q: sp.Rational(4, 5)}):
        at.update({x: sp.Rational(1, 5), y: sp.Rational(2, 7), z: sp.Rational(3, 8)})
        for value in list(chart.geo.g-expected)+[chart.geo.ricci_scalar()-scalar]:
            if abs(complex(sp.N(value.subs(at), 30))) > 1e-20:
                raise AssertionError(f"{METRIC}/{system}: frame metric or Ricci scalar")
