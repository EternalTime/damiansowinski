"""Null charts for a radiating spherical source with positive cosmological constant."""
import sympy as sp

METRIC = "mallett_radiating_star"
CHARTS = ("eddington_finkelstein_outgoing", "eddington_finkelstein_ingoing")
DIMENSIONS = {
    (METRIC, chart): {time: "L", "r": "L", r"\theta": "1", r"\phi": "1",
                      "m": "L", r"\Lambda": "L**-2"}
    for chart, time in zip(CHARTS, ("u", "v"))
}


def chart_spec(system):
    outgoing = system == CHARTS[0]
    time, sign = ("u", "-") if outgoing else ("v", "+")
    f = r"1 - \dfrac{2m}{r} - \dfrac{\Lambda r^2}{3}"
    line = (rf"ds^2=-\left({f}\right)d{time}^2 {sign} 2\,d{time}\,dr"
            r" + r^2\left(d\theta^2+\sin^2\theta\,d\phi^2\right)")
    coords = [time, "r", r"\theta", r"\phi"]
    cross = "-1" if outgoing else "1"
    return {"metric_id": METRIC, "system": {
        "id": system, "name": ("Outgoing" if outgoing else "Ingoing") + " Eddington-Finkelstein",
        "coords": coords, "parameters": [f"m = m({time})", r"\Lambda"],
        "domains": [time + r" \in (-\infty, \infty)", r"r \in (0, \infty)",
                    r"\theta \in [0, \pi]", r"\phi \in [0, 2\pi)",
                    r"\Lambda > 0,\quad m \ge 0", r"0 < 9\Lambda m^2 < 1 \;\text{(two positive marginal radii)}",
                    r"9\Lambda m^2 = 1 \;\text{(one double positive radius)}",
                    r"9\Lambda m^2 > 1 \;\text{(no positive marginal radius)}"],
        "line_element": line}, "chart_line_element": line, "printer": {},
        "components": {"metric_components": {(time, time): r"-\left(" + f + r"\right)",
                       (time, "r"): cross, ("r", time): cross},
                       "inverse_metric_components": {(time, "r"): cross, ("r", time): cross,
                                                     ("r", "r"): f}},
        "check": lambda chart: check_chart(chart, outgoing)}


def check_chart(chart, outgoing):
    """Compare source components and inverse at independent points."""
    w, r, theta, phi = chart.symbols
    mass = chart.reader.parameters["m"]
    lam = chart.reader.parameters["Lambda"]
    sign = -1 if outgoing else 1
    for radius, m0, cosm in ((2, sp.Rational(1, 5), sp.Rational(1, 20)),
                              (3, sp.Rational(2, 5), sp.Rational(1, 50))):
        at = {w: 0, r: radius, theta: sp.Rational(6, 5), phi: 1,
              mass: m0, lam: cosm}
        f = 1 - 2*m0/radius - cosm*radius**2/3
        expected = sp.Matrix([[-f, sign, 0, 0], [sign, 0, 0, 0],
                              [0, 0, radius**2, 0],
                              [0, 0, 0, radius**2 * sp.sin(at[theta])**2]])
        g = chart.geo.g.subs(at).evalf(30)
        inv = chart.geo.ginv.subs(at).evalf(30)
        for label, diff in (("source components", g - expected.evalf(30)),
                            ("inverse", g * inv - sp.eye(4))):
            if max(abs(complex(x)) for x in diff) > 1e-20:
                raise AssertionError(f"{METRIC}: {label} at r={radius}")
