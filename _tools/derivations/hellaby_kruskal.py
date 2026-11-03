"""The comoving chart of the symmetric dust wormhole."""

import sympy as sp


METRIC = "hellaby_kruskal"
SYSTEM = "comoving"
DIMENSIONS = {(METRIC, SYSTEM): {
    "t": "L", "r": "L", "\\theta": "1", "\\phi": "1",
    "Y": "L", "f": "1", "M": "L", "a": "L", "M_0": "L",
    "M_1": "L**-2", "r_s": "L**2", "\\eta": "1",
}}


def chart_spec():
    """The sourced spherical chart, with radius and energy held for tensor printing."""
    coords = ["t", "r", "\\theta", "\\phi"]
    parameters = ["Y = Y(t,r)", "f = f(r)", "M = M(r)", "a = a(r)",
                  "M_0", "M_1", "r_s", "\\eta"]
    line = (r"ds^2=-dt^2+\dfrac{(\partial_r Y)^2}{1+f}dr^2"
            r"+Y^2(d\theta^2+\sin^2\theta\,d\phi^2)")
    return {"metric_id": METRIC, "system": {
        "id": SYSTEM, "name": "Comoving Synchronous", "coords": coords,
        "domains": [r"r \in (-\infty, \infty)", r"a(r) < t < -a(r)",
                    r"\theta \in [0, \pi]", r"\phi \in [0, 2\pi)"],
        "parameters": parameters, "line_element": line},
        "chart_line_element": line, "printer": {"primed": ["f"]},
        "check": dust_check}


def dust_check(chart):
    """Check the dust Einstein tensor at points of the sourced example."""
    t, r, theta, _ = chart.symbols
    Y, f = (chart.reader.parameters[k] for k in ("Y", "f"))
    x, eta = sp.symbols("_x _eta", real=True)
    mass = (1 + x**2 / 10)**3
    u = sp.exp(-x**2)
    radius = mass * (1 - sp.cos(eta)) / u
    time = mass * (eta - sp.sin(eta) - sp.pi) / u**sp.Rational(3, 2)
    radial_radius = sp.diff(radius, x) - sp.diff(radius, eta) * sp.diff(time, x) / sp.diff(time, eta)
    tensor = chart.geo.einstein_ll()
    for xx, ee in ((sp.Rational(1, 2), sp.pi / 2), (sp.Integer(1), 3 * sp.pi / 2)):
        point = {x: xx, eta: ee}
        m = mass.subs(point)
        mr = sp.diff(mass, x).subs(point)
        ff = -u.subs(point)
        fp = -sp.diff(u, x).subs(point)
        rr = radius.subs(point)
        rp = radial_radius.subs(point)
        velocity = sp.sin(ee) * sp.sqrt(u.subs(point)) / (1 - sp.cos(ee))
        rtr = (2 * mr / rr - 2 * m * rp / rr**2 + fp) / (2 * velocity)
        values = {Y: rr, f: ff, sp.diff(Y, r): rp, sp.diff(Y, t): velocity,
                  sp.diff(Y, t, 2): -m / rr**2, sp.diff(Y, r, t): rtr,
                  sp.diff(Y, r, t, 2): -mr / rr**2 + 2 * m * rp / rr**3,
                  sp.diff(f, r): fp, theta: sp.Rational(7, 10)}
        density = 2 * mr / (rr**2 * rp)
        for i in range(4):
            for j in range(4):
                residual = sp.sympify(tensor[i][j]).xreplace(values) - (density if i == j == 0 else 0)
                if abs(complex(sp.N(residual, 30))) > 1e-20:
                    raise AssertionError(f"{METRIC}: dust Einstein tensor at {i},{j}")


def declared_functions(lw_r):
    """Regular forms of the chosen cycloid and (1 - exp(-s))/s, with s=r^2."""
    from functools import lru_cache
    import numpy as np

    t, s = sp.symbols("_t _s", real=True)
    mass = (1 + s / 10)**3
    radius = 2 * mass * sp.exp(s) * lw_r(t * sp.exp(-3 * s / 2) / (2 * mass), 1)
    made = {}

    @lru_cache(maxsize=None)
    def radius_jet(i, j):
        return sp.lambdify((t, s), sp.diff(radius, t, i, s, j), "numpy")

    def radius_type(i, j):
        name = "hk_radius" if i == j == 0 else f"hk_radius_{i}{j}"
        if name not in made:
            def value(tval, sval):
                return radius_jet(i, j)(tval, sval)

            def fdiff(self, argindex=1):
                return radius_type(i + (argindex == 1), j + (argindex == 2))(*self.args)

            made[name] = type(name, (sp.Function,), {
                "nargs": 2, "is_real": True, "_imp_": staticmethod(value), "fdiff": fdiff})
        return made[name]

    def h_value(sval, order):
        sval = np.asarray(sval, dtype=float)
        with np.errstate(all="ignore"):
            if order == 0:
                regular = -np.expm1(-sval) / sval
                series = 1 - sval / 2 + sval**2 / 6 - sval**3 / 24 + sval**4 / 120
            else:
                regular = (np.exp(-sval) * (sval + 1) - 1) / sval**2
                series = -0.5 + sval / 3 - sval**2 / 8 + sval**3 / 30
            return np.where(np.abs(sval) < 1e-4, series, regular)

    class hk_h1(sp.Function):
        nargs = 1
        is_real = True
        _imp_ = staticmethod(lambda value: h_value(value, 1))

    class hk_h(sp.Function):
        nargs = 1
        is_real = True
        _imp_ = staticmethod(lambda value: h_value(value, 0))

        def fdiff(self, argindex=1):
            return hk_h1(self.args[0])

    return {"hk_radius": radius_type(0, 0), "hk_h": hk_h, "hk_h1": hk_h1}


def curvature_views(ck, src, cv):
    """Map the two stored curvature scalars across the regular neck."""
    import math
    import numpy as np

    xa = cv.axis(r"$r/M_0$", -1.5, 1.5)
    ya = cv.axis(r"$t/M_0$", -1.5, 1.5)
    functions = {"Y": "hk_radius(t,r**2)", "f": "-1+r**2*hk_h(r**2)",
                 "M": "(1+r**2/10)**3", "a": "-pi*(1+r**2/10)**3*exp(3*r**2/2)"}
    views = []
    for field, letter, power, name in (("kretschmann", "K", 4, "Kretschmann"),
                                        ("ricci_scalar", "R", 2, "Ricci")):
        stored = cv.Published(src, METRIC, SYSTEM, field, ["r", "t"], [], functions=functions)
        _, _, values = cv.cells(stored, xa, ya, 80, 80, [])
        ck(f"{METRIC} {name} values are finite and positive", np.all(np.isfinite(values)) and np.all(values > 0))
        low = math.floor(math.log10(float(np.min(values))))
        high = math.ceil(math.log10(float(np.max(values))))
        scale = {"kind": "log", "from": 10.0**low, "to": 10.0**high}
        cv.shape_checks(ck, f"{METRIC} {name}", [1])
        views.append({"id": letter, "label": f"{name} scalar", "kind": "map", "system": SYSTEM,
            "x": xa, "y": ya,
            "colour": {"name": f"${letter}M_0^{power}$", "scale": scale, "ticks": cv.value_ticks(scale)},
            "map": {"metric": METRIC, "system": SYSTEM, "scalar": field, "x": "r", "y": "t",
                    "params": {}, "functions": functions, "nx": 80, "ny": 80,
                    "cells": [[cv.rounded(v) for v in row] for row in values]},
            "marks": [], "legend": [],
            "caption": [f"The {name} scalar ${letter}$ across the neck in the plane of $t$ and $r$."],
            "settings": r"$M_1M_0^2=1/10$ and $r_s/M_0^2=1$; lengths in units of $M_0$.",
            "clip": f"Logarithmic colours from $10^{{{low}}}$ to $10^{{{high}}}$."})
    return views
