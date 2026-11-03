"""The comoving chart of the symmetric dust wormhole."""

import sympy as sp


METRIC = "hellaby_kruskal"
SYSTEM = "comoving"
DIMENSIONS = {(METRIC, SYSTEM): {
    "t": "L", "r": "L", "\\theta": "1", "\\phi": "1",
    "R": "L", "f": "1", "M": "L", "a": "L", "M_0": "L",
    "M_1": "L**-2", "r_s": "L**2", "\\eta": "1",
}}


def chart_spec():
    """The sourced spherical chart, with radius and energy held for tensor printing."""
    coords = ["t", "r", "\\theta", "\\phi"]
    parameters = ["R = R(t,r)", "f = f(r)", "M = M(r)", "a = a(r)",
                  "M_0", "M_1", "r_s", "\\eta"]
    line = (r"ds^2=-dt^2+\dfrac{(\partial_r R)^2}{1+f}dr^2"
            r"+R^2(d\theta^2+\sin^2\theta\,d\phi^2)")
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
    R, f = (chart.reader.parameters[k] for k in ("R", "f"))
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
        values = {R: rr, f: ff, sp.diff(R, r): rp, sp.diff(R, t): velocity,
                  sp.diff(R, t, 2): -m / rr**2, sp.diff(R, r, t): rtr,
                  sp.diff(R, r, t, 2): -mr / rr**2 + 2 * m * rp / rr**3,
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
