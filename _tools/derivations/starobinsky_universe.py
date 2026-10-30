"""The homogeneous semiclassical geometry of Starobinsky (1980), Eqs. (4), (7), (12).

All times in the numerical solution are H_0 tau and all lengths are H_0^-1.
The numerical initial displacement is a linear approximation at the de Sitter saddle.
The evolution retains the full homogeneous equation, including the anomaly term.
"""
from functools import lru_cache

import numpy as np
import sympy as sp

import chart_printer as cp
import verify_metrics as vm

ID = "starobinsky_universe"
CHART_IDS = ("cosmic", "conformal", "closed_de_sitter", "flat_de_sitter", "open_de_sitter")
ANGLES = r"\left(d\theta^2 + \sin^2\theta\,d\phi^2\right)"


def chart_spec(system):
    """Two evolving FLRW charts and the three exact special solutions of Eq. (7)."""
    evolving = system in ("cosmic", "conformal")
    first = r"\eta" if system == "conformal" else r"\tau"
    radial = r"\chi" if system in ("closed_de_sitter", "open_de_sitter") else "r"
    coords = [first, radial, r"\theta", r"\phi"]
    if evolving:
        parameters = [f"a = a({first})", "k", "H_0", "M"]
        space = r"\dfrac{dr^2}{1-kr^2} + r^2" + ANGLES
        line = (r"ds^2 = a^2\left(-d\eta^2 + " + space + r"\right)" if system == "conformal"
                else r"ds^2 = -d\tau^2 + a^2\left(" + space + r"\right)")
        domains = [first + r"\in I", r"a>0", r"r\in[0,\infty)\;\text{for}\;k\le0",
                   r"r\in[0,1/\sqrt{k})\;\text{for}\;k>0"]
    else:
        parameters = ["H_0"]
        if system == "flat_de_sitter":
            space = r"dr^2 + r^2" + ANGLES
            factor = r"e^{2H_0\tau}"
        else:
            trig = r"\sin" if system == "closed_de_sitter" else r"\sinh"
            scale = r"\cosh" if system == "closed_de_sitter" else r"\sinh"
            space = r"d\chi^2 + " + trig + r"^2\chi" + ANGLES
            factor = scale + r"^2(H_0\tau)/H_0^2"
        line = r"ds^2 = -d\tau^2 + " + factor + r"\left(" + space + r"\right)"
        domains = [r"\tau\in(0,\infty)" if system == "open_de_sitter" else r"\tau\in(-\infty,\infty)",
                   (r"\chi\in[0,\pi]" if system == "closed_de_sitter" else radial + r"\in[0,\infty)")]
    domains += [r"\theta\in[0,\pi]", r"\phi\in[0,2\pi)"]
    probe = vm.Reader(coords, parameters, ())
    x = probe.symbol[first]
    pretty = cp.hyperbolic(probe.parameters["H_0"] * x) if system in ("closed_de_sitter", "open_de_sitter") else sp.factor
    return {"metric_id": ID,
            "system": {"id": system, "name": {"cosmic": "Cosmic Time", "conformal": "Conformal Time",
                       "closed_de_sitter": "Closed de Sitter", "flat_de_sitter": "Flat de Sitter",
                       "open_de_sitter": "Open de Sitter"}[system], "coords": coords, "domains": domains,
                       "parameters": parameters, "line_element": line},
            "chart_line_element": line, "printer": {"primed": ["a"]}, "pretty": pretty,
            "check": lambda chart: chart_check(chart, system)}


def chart_check(chart, system):
    """Check the independent FLRW Ricci formula at regular points, with no symbolic determinant."""
    x, r, theta, phi = chart.symbols
    p = chart.reader.parameters
    g = chart.geo.g
    ric = chart.geo.ricci_ll()
    evolving = system in ("cosmic", "conformal")
    for time, radius, curvature in ((sp.Rational(2, 5), sp.Rational(1, 3), -1),
                                    (sp.Rational(7, 5), sp.Rational(1, 4), 1)):
        point = {x: time, r: radius, theta: sp.Rational(4, 5), phi: sp.Rational(1, 5),
                 p["H_0"]: sp.Rational(3, 2), chart.reader.c: 1}
        if evolving:
            a = p["a"]
            values = {a: 2, sp.diff(a, x): sp.Rational(3, 7), sp.diff(a, x, 2): sp.Rational(-1, 9)}
            point.update({p["k"]: curvature, p["M"]: sp.Rational(1, 2)})
            A, B, C = (values[a], values[sp.diff(a, x)], values[sp.diff(a, x, 2)])
            if system == "cosmic":
                temporal, spatial = -3*C/A, C/A+2*(B*B+curvature)/A**2
            else:
                temporal, spatial = -3*(C/A-B*B/A**2), (C/A+B*B/A**2+2*curvature)/A**2
        else:
            values = {}
            temporal, spatial = -3*p["H_0"]**2, 3*p["H_0"]**2
        for i in range(4):
            for j in range(4):
                expected = temporal if i == j == 0 else spatial*g[i,j] if i == j else 0
                delta = (sp.sympify(vm._at(ric, (i,j))) - expected).subs(values).subs(point)
                if abs(complex(sp.N(delta, 30))) > 1e-24:
                    raise AssertionError(f"{ID}/{system}: FLRW Ricci component {i}{j} at {time}")
