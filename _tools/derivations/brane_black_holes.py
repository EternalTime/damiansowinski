"""The charts of the black holes on a wall in anti-de Sitter space of four dimensions, Emparan,
Horowitz and Myers, hep-th/9911043, JHEP 2000(01) 007, and its sequel hep-th/9912135,
JHEP 2000(01) 021. brane_black_holes.md beside this file records each chart's source and every
drawing decision.

    flat_wall   the first paper's section 3.1, its equations \\adsc and \\gx as its TeX labels them:
                ds^2 = l^2/(x - y)^2 [-(y^2 + 2 mu y^3) dt^2 + dy^2/(y^2 + 2 mu y^3) + dx^2/G + G dphi^2],
                G = 1 - x^2 - 2 mu x^3, on the side 0 <= x <= x_2 that is glued to a copy of itself
                along the wall x = 0, where G'(0) = 0 and G(0) = 1 give K = g/l.
    two_walls   the sequel's section 2, its equations \\cads, \\hy, \\gx and \\defell:
                ds^2 = 1/(A^2 (x - y)^2) [H dt^2 - dy^2/H + dx^2/G + G dphi^2],
                H = -lambda + k y^2 - 2 m A y^3, G = 1 + k x^2 - 2 m A x^3, l_4 = 1/(A sqrt(lambda + 1)),
                on the quadrant x >= 0, y <= 0 between the walls x = 0 and y = 0, for lambda > 0.
                The first paper's chart is this one at lambda = 0, k = -1, mu = m A, l = 1/A, as the
                sequel's section 5.1 says.
    rotating    the sequel's section 5, its equation \\adsd with its \\hy and \\gx:
                ds^2 = 1/(A^2 (x - y)^2) [H/W (dt + a x^2 dphi)^2 - W/H dy^2 + W/G dx^2 + G/W (dphi - a y^2 dt)^2],
                W = 1 + a^2 x^2 y^2, H = -lambda + k y^2 - 2 m A y^3 - a^2 y^4, G = 1 + k x^2 - 2 m A x^3 + a^2 lambda x^4.

Every coordinate is a pure number, as the papers write them (the first paper's footnote beside \\tohors), so the
line elements carry no c and the lengths stand in front of the bracket. Each chart is held, at
rational points and exactly, to R_ab = -(3/l_4^2) g_ab, to a vanishing Weyl tensor at m = 0, to its
walls having extrinsic curvature proportional to their metric, and to the chart it reduces to.
"""
import sympy as sp

METRIC = "brane_black_holes"
CHARTS = ("flat_wall", "two_walls", "rotating")
NAMES = {"flat_wall": "Flat Wall", "two_walls": "Two Walls", "rotating": "Rotating"}
COORDS = ["t", "y", "x", "\\varphi"]
PARAMETERS = {"flat_wall": ["\\ell", "\\mu"], "two_walls": ["A", "m", "\\lambda", "k"],
              "rotating": ["A", "m", "\\lambda", "k", "a"]}
FLAT_G = "1 - x^2 - 2\\mu x^3"
FLAT_H = "y^2 + 2\\mu y^3"
WALL_H = "-\\lambda + ky^2 - 2mAy^3"
WALL_G = "1 + kx^2 - 2mAx^3"
ROT_H = "-\\lambda + ky^2 - 2mAy^3 - a^2y^4"
ROT_G = "1 + kx^2 - 2mAx^3 + a^2\\lambda x^4"
ROT_W = "1 + a^2x^2y^2"
TIME = "t \\in (-\\infty, \\infty)"
DOMAINS = {
    "flat_wall": [TIME, "y \\in (-\\infty, x)", "x \\in [0, x_2]", "\\varphi \\in [0, \\Delta\\varphi)",
                  "x = 0 \\;\\text{(the wall, where this side meets its mirror image)}",
                  "y = -1/(2\\mu) \\;\\text{(the horizon of the black hole)}",
                  "y = 0 \\;\\text{(the horizon of anti-de Sitter space)}",
                  "y \\to -\\infty \\;\\text{(the curvature singularity)}"],
    "two_walls": [TIME, "y \\in (-\\infty, 0]", "x \\in [0, x_2]", "\\varphi \\in [0, \\Delta\\varphi)", "\\lambda > 0",
                  "x = 0 \\;\\text{(the wall that holds the black hole)}",
                  "y = 0 \\;\\text{(the second wall)}",
                  "y = y_0 \\;\\text{(the horizon of the black hole)}",
                  "y \\to -\\infty \\;\\text{(the curvature singularity)}"],
    "rotating": [TIME, "y \\in (-\\infty, 0]", "x \\in [0, x_2]", "\\varphi \\in [0, \\Delta\\varphi)", "\\lambda \\ge 0",
                 "x = 0 \\;\\text{(the wall that holds the black hole)}",
                 "y = 0 \\;\\text{(the second wall for } \\lambda > 0\\text{)}",
                 "1/y^2 + a^2x^2 = 0 \\;\\text{(the ring singularity)}"],
}
DIMENSIONS = {
    (METRIC, "flat_wall"): {"t": "1", "y": "1", "x": "1", "\\varphi": "1", "\\ell": "L", "\\mu": "1"},
    (METRIC, "two_walls"): {"t": "1", "y": "1", "x": "1", "\\varphi": "1", "A": "1/L", "m": "L", "\\lambda": "1",
                            "k": "1"},
    (METRIC, "rotating"): {"t": "1", "y": "1", "x": "1", "\\varphi": "1", "A": "1/L", "m": "L", "\\lambda": "1",
                           "k": "1", "a": "1"},
}
HELD = {}
RATES = {}
KRETSCHMANN = {
    "flat_wall": "\\dfrac{24}{\\ell^4} + \\dfrac{48\\mu^2\\left(x - y\\right)^6}{\\ell^4}",
    "two_walls": "24A^4\\left(\\lambda + 1\\right)^2 + 48m^2A^6\\left(x - y\\right)^6",
    "rotating": ("24A^4\\left(\\lambda + 1\\right)^2 + \\dfrac{48m^2A^6\\left(x - y\\right)^6"
                 "\\left(1 - 15a^2x^2y^2 + 15a^4x^4y^4 - a^6x^6y^6\\right)}{\\left(1 + a^2x^2y^2\\right)^6}"),
}
RICCI_SCALAR = {"flat_wall": "-\\dfrac{12}{\\ell^2}", "two_walls": "-12A^2\\left(\\lambda + 1\\right)",
                "rotating": "-12A^2\\left(\\lambda + 1\\right)"}

_PRINT = {}


def register(module):
    _PRINT["module"] = module
    module.CHARTS[METRIC] = [lambda system=system: chart_spec(system) for system in CHARTS]


def line(system):
    """The line element of a chart; every coordinate is a pure number, so it is its own chart."""
    if system == "flat_wall":
        return (f"ds^2 = \\dfrac{{\\ell^2}}{{\\left(x - y\\right)^2}}\\left(-\\left({FLAT_H}\\right)dt^2"
                f" + \\dfrac{{dy^2}}{{{FLAT_H}}} + \\dfrac{{dx^2}}{{{FLAT_G}}} + \\left({FLAT_G}\\right)d\\varphi^2\\right)")
    if system == "two_walls":
        return (f"ds^2 = \\dfrac{{1}}{{A^2\\left(x - y\\right)^2}}\\left(\\left({WALL_H}\\right)dt^2"
                f" - \\dfrac{{dy^2}}{{{WALL_H}}} + \\dfrac{{dx^2}}{{{WALL_G}}} + \\left({WALL_G}\\right)d\\varphi^2\\right)")
    return (f"ds^2 = \\dfrac{{1}}{{A^2\\left(x - y\\right)^2}}\\left(\\dfrac{{{ROT_H}}}{{{ROT_W}}}\\left(dt + ax^2d\\varphi\\right)^2"
            f" - \\dfrac{{{ROT_W}}}{{{ROT_H}}}dy^2 + \\dfrac{{{ROT_W}}}{{{ROT_G}}}dx^2"
            f" + \\dfrac{{{ROT_G}}}{{{ROT_W}}}\\left(d\\varphi - ay^2dt\\right)^2\\right)")


def reader(system):
    import verify_metrics as vm
    return vm.Reader(COORDS, PARAMETERS[system], ())


def polynomials(system, probe):
    """The sums each chart's values are printed around, as (placeholder, polynomial, text)."""
    P, X = probe.parameters, probe.symbol
    x, y = X["x"], X["y"]
    named = [(sp.Symbol("BXY"), x - y, "x - y")]
    if system == "flat_wall":
        mu = P["mu"]
        named += [(sp.Symbol("BFG"), 1 - x ** 2 - 2 * mu * x ** 3, FLAT_G),
                  (sp.Symbol("BFH"), 1 + 2 * mu * y, "1 + 2\\mu y")]
        return named
    A, m, lam, k = P["A"], P["m"], P["lambda"], P["k"]
    if system == "two_walls":
        named += [(sp.Symbol("BWG"), 1 + k * x ** 2 - 2 * m * A * x ** 3, WALL_G),
                  (sp.Symbol("BWH"), lam - k * y ** 2 + 2 * m * A * y ** 3, "\\lambda - ky^2 + 2mAy^3")]
        return named
    a = P["a"]
    named += [(sp.Symbol("BRG"), 1 + k * x ** 2 - 2 * m * A * x ** 3 + a ** 2 * lam * x ** 4, ROT_G),
              (sp.Symbol("BRH"), lam - k * y ** 2 + 2 * m * A * y ** 3 + a ** 2 * y ** 4,
               "\\lambda - ky^2 + 2mAy^3 + a^2y^4"),
              (sp.Symbol("BRW"), 1 + a ** 2 * x ** 2 * y ** 2, ROT_W)]
    return named


def shifted(named, x, y):
    """The named factors of each value, and every other sum in x and y written in x - y where that
    leaves it fewer terms, as the Riemann tensor's 1 - 2 mu (x - y)^3."""
    factored = _PRINT["module"].named_factors(named)
    difference = named[0][0]

    def pretty(value):
        out = []
        for f in sp.Mul.make_args(factored(value)):
            base, k = (f.base, f.exp) if f.is_Pow else (f, sp.Integer(1))
            if base.is_Add and base.has(x) and base.has(y):
                moved = sp.expand(base.subs(x, difference + y))
                if len(sp.Add.make_args(moved)) < len(sp.Add.make_args(base)):
                    base = moved
            out.append(base ** k)
        return sp.Mul(*out)
    return pretty


def chart_spec(system):
    """One chart: its line element, its printer and the check at points."""
    probe = reader(system)
    P, X = probe.parameters, probe.symbol
    named = polynomials(system, probe)
    # m before A, so that a term prints as 2mA does in the papers.
    order = sorted(PARAMETERS[system], key=lambda n: n != "m")
    lead = [P[n.lstrip("\\")] for n in order] + [X["y"], X["x"]]
    return {
        "metric_id": METRIC,
        "system": {"id": system, "name": NAMES[system], "coords": COORDS, "domains": DOMAINS[system],
                   "parameters": PARAMETERS[system], "line_element": line(system)},
        "chart_line_element": line(system),
        "printer": {"lead": lead, "named": {p: text for p, _, text in named}},
        "pretty": shifted(named, X["x"], X["y"]),
        "check": lambda chart: check_chart(chart, system),
        "ricci_scalar": RICCI_SCALAR[system],
        "kretschmann": KRETSCHMANN[system],
        **({"components": {"metric_components": ROTATING_METRIC}} if system == "rotating" else {}),
    }


# The rotating chart's metric written as its line element groups it, in H and G.
_ROT_OVER = "}{A^2\\left(x - y\\right)^2\\left(" + ROT_W + "\\right)}"
ROTATING_METRIC = {
    ("t", "t"): "\\dfrac{\\left(" + ROT_H + "\\right) + a^2y^4\\left(" + ROT_G + "\\right)" + _ROT_OVER,
    ("t", "\\varphi"): ("\\dfrac{a\\left(x^2\\left(" + ROT_H + "\\right) - y^2\\left(" + ROT_G + "\\right)\\right)"
                        + _ROT_OVER),
    ("\\varphi", "t"): ("\\dfrac{a\\left(x^2\\left(" + ROT_H + "\\right) - y^2\\left(" + ROT_G + "\\right)\\right)"
                        + _ROT_OVER),
    ("\\varphi", "\\varphi"): ("\\dfrac{a^2x^4\\left(" + ROT_H + "\\right) + \\left(" + ROT_G + "\\right)"
                                 + _ROT_OVER),
}


# Rational points (x, y) inside each chart, outside the horizon and inside it.
POINTS = [(sp.Rational(1, 5), -sp.Rational(3, 4)), (sp.Rational(1, 3), -sp.Integer(5)),
          (sp.Rational(1, 7), -sp.Rational(1, 9))]
VALUES = {"flat_wall": {"ell": sp.Rational(5, 2), "mu": sp.Rational(1, 8)},
          "two_walls": {"A": sp.Rational(3, 2), "m": sp.Rational(1, 5), "lambda": sp.Rational(1, 3), "k": 1},
          "rotating": {"A": sp.Rational(3, 2), "m": sp.Rational(1, 5), "lambda": sp.Rational(1, 3), "k": 1,
                       "a": sp.Rational(2, 7)}}
SIGNS = (-1, 0, 1)


def metric_of(system, probe=None):
    """A chart's metric as a matrix of sympy expressions, and its reader."""
    import verify_metrics as vm
    probe = probe or reader(system)
    return vm.metric_from_line_element(probe, line(system), COORDS), probe


def check_chart(chart, system):
    """The chart at rational points against the field equations of its source, its walls against
    the first paper's \\extcur and the sequel's section 2, and against the chart it reduces to."""
    import verify_metrics as vm
    geo, g = chart.geo, chart.geo.g
    t, y, x, ph = chart.symbols
    P = chart.reader.parameters
    name = f"{METRIC}/{system}"
    if system == "flat_wall":
        inverse_square = 1 / P["ell"] ** 2
        cases = [VALUES[system]]
        massless = {P["mu"]: 0}
    else:
        inverse_square = P["A"] ** 2 * (P["lambda"] + 1)
        cases = [dict(VALUES[system], k=sign) for sign in SIGNS]
        massless = {P["m"]: 0}
    ricci = geo.ricci_ll()

    def at(expression, case, point):
        values = {P[key]: value for key, value in case.items()}
        values.update({x: point[0], y: point[1]})
        return sp.nsimplify(sp.sympify(expression).subs(values))

    for case in cases:
        for point in POINTS:
            for i in range(4):
                for j in range(i, 4):
                    if at(vm._at(ricci, (i, j)) + 3 * inverse_square * g[i, j], case, point) != 0:
                        raise AssertionError(f"{name}: R_ab = -(3/l_4^2) g_ab fails in slot "
                                             f"{chart.coords_tex[i]}{chart.coords_tex[j]} at {case}, {point}")
    weyl = geo.weyl_llll()
    for case in cases:
        for point in POINTS:
            for index in ((0, 1, 0, 1), (0, 2, 0, 2), (1, 2, 1, 2), (2, 3, 2, 3), (0, 3, 0, 3), (0, 1, 2, 3)):
                if at(vm._at(weyl, index).subs(massless), case, point) != 0:
                    raise AssertionError(f"{name}: the Weyl tensor does not vanish at m = 0 in {index} at {point}")

    # Each wall: the derivative of the tangential metric across it is the conformal factor's alone,
    # so its extrinsic curvature is a multiple of its metric, the first paper's \\extcur.
    tangent = {"x": (0, 1, 3), "y": (0, 2, 3)}
    walls = [(x, "x")] + ([] if system == "flat_wall" else [(y, "y")])
    for coordinate, label in walls:
        # d(x - y)^-2/dx = -2(x - y)^-3 and d(x - y)^-2/dy = +2(x - y)^-3.
        sign = 1 if label == "x" else -1
        for i in tangent[label]:
            for j in tangent[label]:
                slope = sp.diff(g[i, j], coordinate) + sign * 2 * g[i, j] / (x - y)
                if vm.norm(slope.subs(coordinate, 0)) != 0:
                    raise AssertionError(f"{name}: the wall {label} = 0 bends its slot "
                                         f"{chart.coords_tex[i]}{chart.coords_tex[j]} beyond the conformal factor")

    if system == "flat_wall":
        # The metric on the wall at x = 0, with t-hat = l t and r = -l/y, is the first paper's \\bhbrane.
        ell, mu = P["ell"], P["mu"]
        r = sp.Symbol("r", positive=True)
        on = {x: 0, y: -ell / r}
        f = 1 - 2 * mu * ell / r
        induced = [(g[0, 0] / ell ** 2, -f), (g[1, 1] * (ell / r ** 2) ** 2, 1 / f), (g[3, 3], r ** 2)]
        for value, want in induced:
            if vm.norm(value.subs(on) - want) != 0:
                raise AssertionError(f"{name}: the metric on the wall is not the first paper's \\bhbrane")
        # The first paper's chart is the sequel's at lambda = 0, k = -1, m = mu l and A = 1/l.
        other, probe = metric_of("two_walls")
        Q = probe.parameters
        same = {Q["A"]: 1 / ell, Q["m"]: mu * ell, Q["lambda"]: 0, Q["k"]: -1}
        same.update(dict(zip([probe.symbol[c] for c in COORDS], chart.symbols)))
        if any(vm.norm(v) != 0 for v in (other.subs(same) - g)):
            raise AssertionError(f"{name}: the chart is not the sequel's at lambda = 0, k = -1")
        for horizon in (-1 / (2 * mu), 0):
            if vm.norm(g[0, 0].subs(y, horizon)) != 0:
                raise AssertionError(f"{name}: g_tt does not vanish on the horizon y = {horizon}")
    elif system == "two_walls":
        # The metric on the wall at x = 0, with rho = -1/y, is the sequel's \\wallgeo.
        A, m, lam, k = P["A"], P["m"], P["lambda"], P["k"]
        rho = sp.Symbol("rho", positive=True)
        f = lam * rho ** 2 - k - 2 * m * A / rho
        on = {x: 0, y: -1 / rho}
        induced = [(g[0, 0], -f / A ** 2), (g[1, 1] / rho ** 4, 1 / (A ** 2 * f)), (g[3, 3], rho ** 2 / A ** 2)]
        for value, want in induced:
            if vm.norm(value.subs(on) - want) != 0:
                raise AssertionError(f"{name}: the metric on the wall is not the sequel's \\wallgeo")
    else:
        # At a = 0 the rotating chart is the two-wall chart, and on the wall it is the sequel's \\branegeo.
        other, probe = metric_of("two_walls")
        Q = probe.parameters
        same = {Q[n]: P[n] for n in ("A", "m", "lambda", "k")}
        same.update(dict(zip([probe.symbol[c] for c in COORDS], chart.symbols)))
        if any(vm.norm(v) != 0 for v in (other.subs(same) - g.subs(P["a"], 0))):
            raise AssertionError(f"{name}: the rotating chart at a = 0 is not the two-wall chart")
        A, m, lam, k, a = P["A"], P["m"], P["lambda"], P["k"], P["a"]
        rho = sp.Symbol("rho", positive=True)
        f = lam * rho ** 2 - k - 2 * m * A / rho + a ** 2 / rho ** 2
        on = {x: 0, y: -1 / rho}
        induced = [(g[0, 0], (-f + a ** 2 / rho ** 2) / A ** 2), (g[0, 3], -a / A ** 2),
                   (g[1, 1] / rho ** 4, 1 / (A ** 2 * f)), (g[3, 3], rho ** 2 / A ** 2)]
        for value, want in induced:
            if vm.norm(value.subs(on) - want) != 0:
                raise AssertionError(f"{name}: the metric on the wall is not the sequel's \\branegeo")
