"""The charts of the Koutras-McIntosh spacetimes, the conformally flat pure radiation of Wils, of
Koutras and McIntosh and of Edgar and Ludwig. koutras_mcintosh.md beside this file records each
chart's source.

Every chart is in the signature (-,+,+,+) of Edgar and Ludwig, CQG 14 L65 (1997), gr-qc/9612059:

    koutras_mcintosh  their eq. (3), Koutras and McIntosh's family with constants a and b,
                      g_uu = -2f(ax + b)(x^2 + y^2) + a^2 w^2, g_uw = -(ax + b), g_ux = aw,
                      which is Griffiths and Podolsky's sec. 1 and Kruglikov and Steneker's (1.1),
                      both in the signature (+,-,-,-), with their v for w and, in the second, f for 2f;
    edgar_ludwig      their eq. (8), the whole class that is not a plane wave, in Wils's coordinates
                      (u, v, x, y): g_uu = -2fx(x^2 + y^2 + gy + h) + v^2/x^2, g_uv = -1, g_ux = 2v/x,
                      Wils's metric, their eq. (1), at g = h = 0;
    kerr_schild       Griffiths and Podolsky's eq. (12), the inertial coordinates of the flat space
                      the radiation crosses, for Wils's metric with f = 1/l constant: the flat metric
                      plus (G/x) k k with G = -(x^2 + Y^2)/(2l), the profile of the Kundt form, and
                      k = (1 + u^2) c dT - 2u dX - (1 - u^2) dZ, which is 2x du.

u is a pure number labelling the wave fronts and x and y lengths on them; w is a length, v an
area, a a pure number, b a length, f an inverse length, g a length and h an area.
"""
import sympy as sp

METRIC = "koutras_mcintosh"
CHARTS = ("koutras_mcintosh", "edgar_ludwig", "kerr_schild")
NAMES = {"koutras_mcintosh": "Koutras-McIntosh", "edgar_ludwig": "Edgar-Ludwig", "kerr_schild": "Kerr-Schild"}

DIMENSIONS = {
    (METRIC, "koutras_mcintosh"): {"u": "1", "w": "L", "x": "L", "y": "L", "a": "1", "b": "L", "f": "1/L"},
    (METRIC, "edgar_ludwig"): {"u": "1", "v": "L^2", "x": "L", "y": "L", "f": "1/L", "g": "L", "h": "L^2"},
    (METRIC, "kerr_schild"): {"T": "T", "X": "L", "Y": "L", "Z": "L", "\\ell": "L", "x": "L", "u": "1", "v": "1"},
}

# The Kerr-Schild chart holds the x, u and v of the Kundt form as functions of cT, X and Z, as the
# Kerr-Schild chart of Kundt's waves does, with the same rates: 2x du is the null covector of the
# Kerr-Schild form, and every derivative of the three is rational in them.
HELD = {(METRIC, "kerr_schild"): ("x", "u", "v")}
RATES = {(METRIC, "kerr_schild"): {
    "x": {"T": "-c\\left(v + u\\left(1 + uv\\right)\\right)", "X": "1 + 2uv",
          "Z": "v - u\\left(1 + uv\\right)"},
    "u": {"T": "\\dfrac{c\\left(1 + u^2\\right)}{2x}", "X": "-\\dfrac{u}{x}",
          "Z": "-\\dfrac{1 - u^2}{2x}"},
    "v": {"T": "\\dfrac{c\\left(1 + 2v\\left(v + u\\left(1 + uv\\right)\\right)\\right)}{2x}",
          "X": "-\\dfrac{v\\left(1 + 2uv\\right)}{x}",
          "Z": "\\dfrac{1 - 2v\\left(v - u\\left(1 + uv\\right)\\right)}{2x}"},
}}

REALS = " \\in (-\\infty, \\infty)"
ENVELOPE_X = "x = 0 \\;\\text{(the envelope of the wave fronts)}"
KS_X = "x = \\sqrt{X^2 + Z^2 - c^2T^2}"
KS_U = "u = \\dfrac{X - x}{cT + Z}"
KS_V = "v = \\dfrac{cT + Z}{2x}"
KS_NULL = "\\left(1 + u^2\\right)c\\,dT - 2u\\,dX - \\left(1 - u^2\\right)dZ"
KS_WAVE = "\\dfrac{x^2 + Y^2}{2\\ell x}"

SYSTEMS = {
    "koutras_mcintosh": {
        "coords": ["u", "w", "x", "y"], "parameters": ["a", "b", "f = f(u)"],
        "domains": ["u" + REALS, "w" + REALS, "x" + REALS, "y" + REALS,
                    "ax + b > 0 \\;\\text{(the side of the envelope the radiation fills)}",
                    "ax + b = 0 \\;\\text{(the envelope of the wave fronts, for } a \\neq 0\\text{)}"],
        "line_element": ("ds^2 = dx^2 + dy^2 - 2\\left(ax + b\\right)du\\,dw + 2aw\\,du\\,dx"
                         " + \\left(a^2w^2 - 2f\\left(ax + b\\right)\\left(x^2 + y^2\\right)\\right)du^2")},
    "edgar_ludwig": {
        "coords": ["u", "v", "x", "y"], "parameters": ["f = f(u)", "g = g(u)", "h = h(u)"],
        "domains": ["u" + REALS, "v" + REALS, "x \\in (0, \\infty)", "y" + REALS, ENVELOPE_X],
        "line_element": ("ds^2 = dx^2 + dy^2 - 2\\,du\\,dv + \\dfrac{4v}{x}du\\,dx"
                         " + \\left(\\dfrac{v^2}{x^2} - 2fx\\left(x^2 + y^2 + gy + h\\right)\\right)du^2")},
    "kerr_schild": {
        "coords": ["T", "X", "Y", "Z"], "parameters": ["\\ell", KS_X, KS_U, KS_V],
        "domains": ["T" + REALS, "X" + REALS, "Y" + REALS, "Z" + REALS,
                    "X^2 + Z^2 > c^2T^2 \\;\\text{(outside the envelope of the wave fronts)}",
                    "cT + Z = 0,\\; X < 0 \\;\\text{(the wave front } u \\to \\pm\\infty\\text{)}"],
        "line_element": "ds^2 = -c^2dT^2 + dX^2 + dY^2 + dZ^2 - " + KS_WAVE + "\\left(" + KS_NULL + "\\right)^2",
        "chart": ("ds^2 = -dT^2 + dX^2 + dY^2 + dZ^2 - " + KS_WAVE
                  + "\\left(" + KS_NULL.replace("c\\,dT", "dT") + "\\right)^2")},
}

PRINT = None  # print_charts, handed over by register, for the Kerr-Schild forms of Kundt's waves


def register(module):
    global PRINT
    PRINT = module
    module.CHARTS[METRIC] = [lambda system=system: chart_spec(system) for system in CHARTS]


def chart_spec(system):
    """One chart: its line element, the same in the chart x^0 = cT where it has a time, the printer's
    order of names, and the check of the chart at points."""
    import chart_printer as cp
    import verify_metrics as vm
    chart = SYSTEMS[system]
    probe = vm.Reader(chart["coords"], chart["parameters"], (), held=HELD.get((METRIC, system), ()))
    spec = {
        "metric_id": METRIC,
        "system": {"id": system, "name": NAMES[system], "coords": chart["coords"], "domains": chart["domains"],
                   "parameters": chart["parameters"], "line_element": chart["line_element"]},
        "chart_line_element": chart.get("chart", chart["line_element"]),
        # Constants and functions before the coordinates, so that a value reads as the line element
        # writes it, ax + b and 2fx, with the bracket ax + b kept whole.
        "printer": {"lead": ([*probe.symbol.values(), *probe.parameters.values()] if system == "edgar_ludwig"
                             else [*probe.parameters.values(), *probe.symbol.values()]),
                    "factors": [*probe.parameters.values(), *probe.symbol.values()], "flip": False,
                    **({"collect": lambda poly, pr: cp.collect_by(poly, amplitudes(probe), pr)} if system != "kerr_schild" else {})},
        "check": lambda built: check_chart(built, system),
    }
    if system == "kerr_schild":
        forms = PRINT.KundtFlatForms(probe)
        spec["pretty"], spec["reduce"] = forms.pretty, forms.reduce
        spec["printer"] = {"lead": [*PRINT.KUNDT_NAMES, probe.symbol["Y"], probe.parameters["ell"]], "flip": False}
    return spec


def amplitudes(probe):
    """The amplitude f and its first two rates, by which a sum is collected, so that the flat part
    stands apart from the radiation's and g and h stay inside the bracket f multiplies."""
    u = probe.symbol["u"]
    out = []
    for name in ("f",):
        if name in probe.parameters:
            fn = probe.parameters[name]
            out += [sp.Derivative(fn, (u, 2)), sp.Derivative(fn, u), fn]
    return out


# Points at which every check is made, and test functions of u that are no special case.
POINTS = [(sp.Rational(3, 10), sp.Rational(17, 10), sp.Rational(2, 5), sp.Rational(-7, 10)),
          (sp.Rational(-11, 10), sp.Rational(1, 3), sp.Rational(9, 10), sp.Rational(6, 5)),
          (sp.Rational(7, 4), sp.Rational(-13, 20), sp.Rational(6, 5), sp.Rational(1, 4))]
# Events outside the cylinder X^2 + Z^2 = T^2, off the front cT + Z = 0, in the chart x^0 = cT.
EVENTS = [(sp.Rational(3, 10), sp.Rational(17, 10), sp.Rational(2, 5), sp.Rational(9, 10)),
          (sp.Rational(-1, 2), sp.Rational(21, 10), sp.Rational(-4, 5), sp.Rational(7, 10)),
          (sp.Rational(11, 10), sp.Rational(3, 5), sp.Rational(13, 10), sp.Rational(19, 10))]


def test_functions(u):
    return {"f": sp.Rational(1, 2) + u / 5 + sp.sin(u) / 7,
            "g": sp.Rational(2, 3) - u ** 2 / 11,
            "h": sp.Rational(1, 4) + sp.cos(u) / 3}


def eq7(u, w, x, y, f, g, h):
    """Edgar and Ludwig's eq. (7), the class in Koutras and McIntosh's w = v/x, as a matrix in
    (u, w, x, y): g_uu = -2fx(x^2 + y^2 + gy + h) + w^2, g_uw = -x, g_ux = w, g_xx = g_yy = 1."""
    return sp.Matrix([[-2 * f * x * (x ** 2 + y ** 2 + g * y + h) + w ** 2, -x, w, 0],
                      [-x, 0, 0, 0], [w, 0, 1, 0], [0, 0, 0, 1]])


def check_chart(chart, system):
    """Each chart at points against what its sources state.

    Koutras and McIntosh's chart is conformally flat, every component of its Weyl tensor zero, with
    the Ricci tensor R_uu = 4(ax + b)f and nothing else, Griffiths and Podolsky's Phi_22 = 2f/(ax + b)
    in their frame, at four pairs (a, b), among them the plane wave (0, 1) and Wils's metric (1, 0).
    Edgar and Ludwig's chart is the same with R_uu = 4xf, and it is their eq. (7) pulled back along
    w = v/x, which is their eq. (8); at a = 1, b = 0 eq. (7) is Koutras and McIntosh's chart with
    g = h = 0. The Kerr-Schild chart is eq. (7) at f = 1/l and g = h = 0 pulled back along Griffiths and
    Podolsky's eq. (7), w = cT + Z, x = sqrt(X^2 + Z^2 - c^2T^2) and u = (X - x)/(cT + Z)."""
    import verify_metrics as vm
    name = f"{METRIC}/{system}"
    s = chart.symbols
    geo, g = chart.geo, chart.geo.g
    P = chart.reader.parameters

    def vanishes(value, what, at):
        number = complex(sp.N(sp.sympify(value).subs(at), 30))
        if abs(number) > 1e-20:
            raise AssertionError(f"{name}: {what} at {at}, off by {number}")

    if system == "kerr_schild":
        T, X, Y, Z = s
        xe = sp.sqrt(X ** 2 + Z ** 2 - T ** 2)
        ue = (X - xe) / (T + Z)
        put = {P["x"]: xe, P["u"]: ue, P["v"]: (T + Z) / (2 * xe)}
        ell = P["ell"]
        image = [ue, T + Z, xe, Y]
        J = sp.Matrix(4, 4, lambda i, j: sp.diff(image[i], s[j]))
        pulled = J.T * eq7(*image, 1 / ell, 0, 0) * J
        here = g.subs(put).doit()
        for e in EVENTS:
            at = {**dict(zip(s, e)), ell: sp.Rational(7, 3)}
            for a in range(4):
                for b in range(a, 4):
                    vanishes(here[a, b] - pulled[a, b],
                             f"the chart is not Wils's metric pulled back in slot {chart.coords_tex[a]}{chart.coords_tex[b]}", at)
        return

    u = s[0]
    tests = test_functions(u)
    weyl, ricci = geo.weyl_llll(), geo.ricci_ll()
    if system == "koutras_mcintosh":
        cases = [(1, 0), (0, 1), (sp.Rational(3, 5), sp.Rational(2, 7)), (sp.Rational(-2, 3), sp.Rational(5, 4))]
        for a0, b0 in cases:
            put = lambda e: sp.sympify(e).subs(P["f"], tests["f"]).doit().subs({P["a"]: a0, P["b"]: b0})  # noqa: E731
            want = 4 * (a0 * s[2] + b0) * tests["f"]
            for p in POINTS:
                at = dict(zip(s, p))
                for a in range(4):
                    for b in range(4):
                        vanishes(put(vm._at(ricci, (a, b))) - (want if (a, b) == (0, 0) else 0),
                                 f"R_ab is not R_uu = 4(ax + b)f alone at (a, b) = ({a0}, {b0}), slot {a}{b}", at)
                        for c in range(4):
                            for d in range(c + 1, 4):
                                vanishes(put(vm._at(weyl, (a, b, c, d))), f"the Weyl tensor is not zero at (a, b) = ({a0}, {b0})", at)
        # At a = 1, b = 0 it is Edgar and Ludwig's eq. (7) with g = h = 0.
        here = g.subs({P["a"]: 1, P["b"]: 0})
        there = eq7(*s, P["f"], 0, 0)
        for p in POINTS:
            at = dict(zip(s, p))
            for a in range(4):
                for b in range(a, 4):
                    vanishes((here[a, b] - there[a, b]).subs(P["f"], tests["f"]).doit(),
                             "at a = 1, b = 0 the chart is not Wils's metric in w", at)
        return

    put = lambda e: sp.sympify(e).subs({P[k]: v for k, v in tests.items()}).doit()  # noqa: E731
    want = 4 * s[2] * tests["f"]
    for p in POINTS:
        at = dict(zip(s, p))
        for a in range(4):
            for b in range(4):
                vanishes(put(vm._at(ricci, (a, b))) - (want if (a, b) == (0, 0) else 0),
                         f"R_ab is not R_uu = 4xf alone, slot {a}{b}", at)
                for c in range(4):
                    for d in range(c + 1, 4):
                        vanishes(put(vm._at(weyl, (a, b, c, d))), "the Weyl tensor is not zero", at)
    # Edgar and Ludwig's eq. (7) pulled back along w = v/x is the chart, their eq. (8).
    uu, v, x, y = s
    image = [uu, v / x, x, y]
    J = sp.Matrix(4, 4, lambda i, j: sp.diff(image[i], s[j]))
    pulled = J.T * eq7(*image, P["f"], P["g"], P["h"]) * J
    for p in POINTS:
        at = dict(zip(s, p))
        for a in range(4):
            for b in range(a, 4):
                vanishes(put(g[a, b] - pulled[a, b]), "the chart is not eq. (7) pulled back along w = v/x", at)
