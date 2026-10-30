"""Eleven dimensional wave: the two sourced charts and their static drawings.

Sources and normalization are recorded in kowalski_glikman.md.
"""

ID = "kowalski_glikman"
CHARTS = ("brinkmann", "rosen")

import itertools
import math

import numpy as np
import sympy as sp


def chart_spec(system):
    """Figueroa-O'Farrill and Papadopoulos (12), (14); Blau et al. (14)-(17)."""
    import chart_printer as cp
    import verify_metrics as vm

    rosen = system == "rosen"
    null, letter = ("w", "y") if rosen else ("v", "x")
    coords = ["u", null] + [f"{letter}_{i}" for i in range(1, 10)]
    blocks = [" + ".join(f"{letter}_{i}^2" for i in block) for block in (range(1, 4), range(4, 10))]
    differentials = [" + ".join(f"d{letter}_{i}^2" for i in block) for block in (range(1, 4), range(4, 10))]
    if rosen:
        line = (r"ds^2 = -2{c}\,du\,dw + \sin^2\left(\mu {c}u/3\right)\left(" + differentials[0]
                + r"\right) + \sin^2\left(\mu {c}u/6\right)\left(" + differentials[1] + r"\right)")
    else:
        line = (r"ds^2 = -\mu^2\left(\frac{" + blocks[0] + r"}{9} + \frac{" + blocks[1]
                + r"}{36}\right){c2}du^2 - 2{c}\,du\,dv + " + " + ".join(differentials))
    domains = [c + r" \in (-\infty,\infty)" for c in coords]
    if rosen:
        domains[0] = r"u \in (0,3\pi/(\mu c))"
    reader = vm.Reader(coords, [r"\mu"], ())
    return {"metric_id": ID, "system": {"id": system, "name": system.title(), "coords": coords,
            "domains": domains, "parameters": [r"\mu"],
            "line_element": line.replace("{c2}", "c^2").replace("{c}", "c")},
            "chart_line_element": line.replace("{c2}", "").replace("{c}", ""),
            "printer": {"lead": [reader.parameters["mu"], *reader.symbol.values()], "flip": False},
            "check": lambda chart: check_chart(chart, system), **({"time": "u"} if rosen else {})}


def check_chart(chart, system):
    """Source curvature and the Rosen pullback at three regular points, without symbolic reduction."""
    import verify_metrics as vm

    mu = chart.reader.parameters["mu"]
    x = chart.symbols
    points = [{**{s: sp.Rational((i + 2) * n, 13) for i, s in enumerate(x)},
               x[0]: sp.Rational(n + 1, 3), mu: sp.Rational(n + 2, 4)} for n in (1, 2, 3)]

    def near(value, label):
        value = sp.sympify(value)
        if value == 0:
            return
        for p in points:
            if abs(complex(value.subs(p).evalf(45))) > 1e-30:
                raise AssertionError(f"{ID}/{system}: {label}")

    if system == "brinkmann":
        # Source (5): all independent R_uiuj, including their vanishing off-diagonal entries.
        r = chart.geo.riemann_llll()
        C = chart.geo.weyl_llll()
        for a in range(11):
            for b in range(a + 1, 11):
                for c in range(11):
                    for d in range(c + 1, 11):
                        tidal = mu**2 / (9 if b < 5 else 36) if a == c == 0 and b == d and b >= 2 else 0
                        weyl = mu**2 / (18 if b < 5 else -36) if a == c == 0 and b == d and b >= 2 else 0
                        near(vm._at(r, (a, b, c, d)) - tidal, f"R_{a}{b}{c}{d}")
                        near(vm._at(C, (a, b, c, d)) - weyl, f"C_{a}{b}{c}{d}")
        ricci = chart.geo.ricci_ll()
        for a in range(11):
            for b in range(11):
                near(vm._at(ricci, (a, b)) - (mu**2/2 if a == b == 0 else 0), f"Ricci {a},{b}")
        return

    # The map is evaluated before the matrices are multiplied: no symbolic determinant or simplify.
    u, w, *ys = x
    omega = [mu / (3 if i < 3 else 6) for i in range(9)]
    scale = [sp.sin(o*u) for o in omega]
    image = sp.Matrix([u, w + sum(a*sp.diff(a, u)*y*y/2 for a, y in zip(scale, ys))]
                      + [a*y for a, y in zip(scale, ys)])
    J = image.jacobian(x)
    H = sum(o*o*X*X for o, X in zip(omega, image[2:]))
    gB = sp.diag(0, 0, *([1]*9))
    gB[0, 0], gB[0, 1], gB[1, 0] = -H, -1, -1
    for point in points:
        j = J.subs(point).evalf(45)
        pulled = j.T * gB.subs(point).evalf(45) * j
        target = chart.geo.g.subs(point).evalf(45)
        if max(abs(v) for v in pulled - target) > 1e-30:
            raise AssertionError(f"{ID}/{system}: Brinkmann metric pulled back")
