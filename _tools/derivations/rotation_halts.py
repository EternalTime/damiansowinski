"""The cylindrical dust junction, from Gleiser and Ramirez (1106.3122), equations (1)-(15)."""

import copy
import json
from pathlib import Path

import sympy as sp

METRIC = "rotation_halts"
SYSTEMS = ("interior", "exterior", "shell")
ROOT = Path(__file__).resolve().parents[2]


def register_charts(pc):
    pc.CHARTS[METRIC] = [lambda side=side: chart_spec(pc, side) for side in SYSTEMS]


def chart_spec(pc, side):
    data = json.loads((ROOT / "MFS/assets/data/metrics/rotation_halts.json").read_text())
    prose = next(c for c in data["coordinates"] if c["id"] == side)
    params = [p["symbol"] for p in prose["parameters"]]
    if side != "shell":
        spec = copy.deepcopy(pc.einstein_rosen("cylindrical"))
        spec["metric_id"] = METRIC
        spec["system"].update(id=side, name=prose["name"], parameters=params,
                              line_element=spec["chart_line_element"])
        radial = "\\rho \\in (0, R]" if side == "interior" else "\\rho \\in [R, \\infty)"
        spec["system"]["domains"] = ["t \\in (-\\infty, \\infty)", radial,
            "\\phi \\in [0, 2\\pi)", "z \\in (-\\infty, \\infty)"]
        spec["check"] = lambda chart: bulk_check(chart, side)
        return spec
    line = "ds^2 = -d\\tau^2 + R^2e^{-2\\psi}d\\phi^2 + e^{2\\psi}dz^2"
    return {"metric_id": METRIC, "system": {
        "id": side, "name": prose["name"], "coords": prose["coords"],
        "domains": ["\\tau \\in (-\\infty, \\infty)", "\\phi \\in [0, 2\\pi)",
                    "z \\in (-\\infty, \\infty)", "R > 0"],
        "parameters": params, "line_element": line},
        "chart_line_element": line, "printer": {}, "check": shell_check}


def residual_at(expression, replacements, point):
    """Evaluate an independent expression at a point, without symbolic cancellation."""
    value = expression.subs(replacements).doit().subs(point).evalf(35)
    return abs(complex(value))


def bulk_check(chart, side):
    t, rho = chart.symbols[:2]
    psi, gamma = chart.reader.parameters["psi"], chart.reader.parameters["gamma"]
    profiles = [{psi: -sp.Rational(1, 3) * sp.log(rho / 2),
                 gamma: sp.Rational(2, 5) + sp.log(rho / 2) / 9},
                {psi: t / 7, gamma: rho**2 / 98}]
    ricci = chart.geo.ricci_ll()
    worst = 0.0
    for functions in profiles:
        for T, r in [(sp.Rational(1, 5), sp.Rational(7, 5)),
                     (sp.Rational(-2, 3), sp.Rational(11, 4))]:
            point = {t: T, rho: r}
            expected = sp.diag(-sp.exp(2 * (gamma - psi)), sp.exp(2 * (gamma - psi)),
                               rho**2 * sp.exp(-2 * psi), sp.exp(2 * psi))
            for a in range(4):
                for b in range(4):
                    worst = max(worst, residual_at(chart.geo.g[a, b] - expected[a, b], functions, point),
                                residual_at(ricci[a][b], functions, point))
    if worst > 1e-25:
        raise AssertionError(f"{side}: metric or vacuum residual {worst:.3g}")
    if side == "interior":
        flat = {psi: sp.S.Zero, gamma: sp.S.Zero}
        for entry in sp.flatten(chart.geo.riemann_ulll()):
            worst = max(worst, residual_at(entry, flat, {t: 1, rho: 2}))
        if worst > 1e-25:
            raise AssertionError(f"restricted flat interior: Riemann residual {worst:.3g}")


def shell_check(chart):
    tau = chart.symbols[0]
    R, psi = chart.reader.parameters["R"], chart.reader.parameters["psi"]
    expected = sp.diag(-1, R**2 * sp.exp(-2 * psi), sp.exp(2 * psi))
    functions = {R: 2 + tau**2 / 10, psi: tau / 8}
    worst = max(residual_at(chart.geo.g[a, b] - expected[a, b], functions, {tau: sp.Rational(1, 3)})
                for a in range(3) for b in range(3))
    if worst > 1e-25:
        raise AssertionError(f"induced shell metric residual {worst:.3g}")


def acceleration(R, V, alpha, lam):
    """Equation (15), independently obtained from (6) and (9), with J renamed alpha."""
    return ((1 + V*V) * alpha*alpha / (R * (R*R + alpha*alpha))
            - lam * (1 + V*V)**0.5 * (R*R + 2*alpha*alpha)**2
            / (R*R * (R*R + alpha*alpha)**1.5))
