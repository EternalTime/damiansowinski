"""The charts of the generalised Vaidya spacetime, Vaidya's metric with a mass function of the null
time and the areal radius. generalised_vaidya.md beside this file records each chart's source.

Every chart is written in a null time, a time, and printed in the chart x^0 = cv or x^0 = cu, so a
derivative along the null time is the derivative along x^0. The ingoing and outgoing charts are Wang
and Wu's line element, gr-qc/9803038 eq. (1) with psi = 0 and epsilon = +1 and -1, which is Husain's
eq. (1) with psi = 0 in the ingoing case and Glass and Krisch's eq. (1), gr-qc/9803040, in the outgoing
one, there in the signature (+,-,-,-):

    ds^2 = -(1 - 2m/r) c^2 dw^2 + 2 epsilon c dw dr + r^2 dOmega^2,   w = v or u,

with m(w, r) free, a length. Its Einstein tensor is Wang and Wu's eq. (2), G^w_w = G^r_r = -2 d_r m/r^2,
G^r_w = 2 d_w m/r^2 and G^theta_theta = G^phi_phi = -d_r^2 m/r, whatever epsilon is. The third chart
is Husain's null fluid with P = k rho, his eq. (15), the ingoing chart at m = f - g/((2k - 1) r^(2k - 1)),
his eq. (12), where his energy density and pressure, eqs. (5), (7) and (13), are rho = g/(4 pi r^(2k + 2))
and P = k rho in units G = 1. The chart writes his g as q^(2k) with q a length, so that the term in the
lapse is the ratio (q/r)^(2k) and g >= 0, the sign his energy conditions ask for, is q real.
"""
import json

import sympy as sp

METRIC = "generalised_vaidya"
CHARTS = ("ingoing", "outgoing", "husain")

DIMENSIONS = {
    (METRIC, "ingoing"): {"v": "T", "r": "L", "\\theta": "1", "\\phi": "1", "m": "L"},
    (METRIC, "outgoing"): {"u": "T", "r": "L", "\\theta": "1", "\\phi": "1", "m": "L"},
    (METRIC, "husain"): {"v": "T", "r": "L", "\\theta": "1", "\\phi": "1", "k": "1", "f": "L", "q": "L", "h": "1", "m": "L"},
}

SPHERE = " + r^2\\left(d\\theta^2 + \\sin^2\\theta\\,d\\phi^2\\right)"
ANGLES = ["\\theta \\in [0, \\pi]", "\\phi \\in [0, 2\\pi)"]
REALS = " \\in (-\\infty, \\infty)"
TRAPPED = "r = 2m \\;\\text{(the marginally trapped spheres)}"
HUSAIN_M = "m = f - \\dfrac{q}{2k - 1}\\left(\\dfrac{q}{r}\\right)^{2k - 1}"
HUSAIN_LAPSE = "1 - \\dfrac{2f}{r} + \\dfrac{2}{2k - 1}\\left(\\dfrac{q}{r}\\right)^{2k}"
HUSAIN_H = "h = \\left(\\dfrac{q}{r}\\right)^{2k}"
HUSAIN_BARE = "1 - \\dfrac{2f}{r} + \\dfrac{2h}{2k - 1}"
# Husain's term in the lapse is held as a function of v and r with its two rates of change, as the
# repulson's Delta is, so that every value is printed in it. The checker reads v as a time, so the rate
# along it is c times the derivative along x^0, as the rates along T of the Kerr-Schild chart of Kundt's
# waves are.
HELD = {(METRIC, "husain"): ("h",)}
RATES = {(METRIC, "husain"): {"h": {"r": "-\\dfrac{2k\\,h}{r}", "v": "\\dfrac{2k\\,c\\,h\\,q'}{q}"}}}
# In the chart x^0 = cv that print_charts works in, v is x^0 itself and q' its derivative along it.
PRINT_RATES = {(METRIC, "husain"): {"h": {"r": "-\\dfrac{2k\\,h}{r}", "v": "\\dfrac{2k\\,h\\,q'}{q}"}}}
NAMES = {"ingoing": "Ingoing", "outgoing": "Outgoing", "husain": "Husain's Null Fluid"}


def chart_spec(system):
    """One chart: the line element with c, the same in the chart x^0 = cw, the metric and its inverse
    as the line element writes them, and the check of the chart at points."""
    import chart_printer as cp
    import verify_metrics as vm
    null = "u" if system == "outgoing" else "v"
    coords = [null, "r", "\\theta", "\\phi"]
    sign = "-" if system == "outgoing" else "+"
    one = "-1" if system == "outgoing" else "1"
    if system == "husain":
        parameters = ["k", "f = f(v)", "q = q(v)", HUSAIN_H, HUSAIN_M]
        lapse = HUSAIN_BARE
    else:
        parameters = [f"m = m({null},r)"]
        lapse = "1 - \\dfrac{2m}{r}"
    wrapped = "\\left(" + lapse + "\\right)"
    bare = HUSAIN_BARE if system == "husain" else lapse
    line = "ds^2 = -" + wrapped + "@C2@d" + null + "^2 " + sign + " 2@C@\\,d" + null + "\\,dr" + SPHERE
    domains = [null + REALS, "r \\in (0, \\infty)"] + ANGLES + [TRAPPED]
    metric = {(null, null): "-\\left(" + bare + "\\right)", (null, "r"): one, ("r", null): one}
    inverse = {(null, "r"): one, ("r", null): one, ("r", "r"): bare}
    key = (METRIC, system)
    probe = vm.Reader(coords, parameters, (), held=HELD.get(key, ()), rates=PRINT_RATES.get(key))
    r = probe.symbol["r"]
    if system == "husain":
        f, q, k = probe.parameters["f"], probe.parameters["q"], probe.parameters["k"]
        rates = [sp.Derivative(f, probe.symbol["v"]), sp.Derivative(q, probe.symbol["v"])]
        # A product is written as the line element writes it, f and q before the radius, and the
        # rates of change of f and q stand apart from what the static hole has already.
        H = probe.parameters["h"]
        printer = {"primed": ["f", "q"], "lead": rates + [H, r, f, q, k], "factors": [k, f, q, r, H] + rates,
                   "rising": [H], "flip": False, "collect": lambda poly, pr: cp.collect_by(poly, rates + [H], pr)}
    else:
        m = probe.parameters["m"]
        w = probe.symbol[null]
        rates = [sp.Derivative(m, w), sp.Derivative(m, r), sp.Derivative(m, (r, 2))]
        printer = {"lead": rates + [r, m], "factors": [m, r] + rates, "flip": False,
                   "collect": lambda poly, pr: cp.collect_by(poly, rates, pr)}
    return {
        "metric_id": METRIC,
        "system": {"id": system, "name": NAMES[system], "coords": coords, "domains": domains,
                   "parameters": parameters, "line_element": line.replace("@C2@", "c^2").replace("@C@", "c")},
        "chart_line_element": line.replace("@C2@", "").replace("@C@", ""),
        "printer": printer,
        "components": {"metric_components": metric, "inverse_metric_components": inverse},
        "check": lambda chart: check_chart(chart, system),
        **({"reduce": lambda value: vm.norm(probe.by_rates(value)),
            # Schwarzschild's 48 f^2/r^6, the cross term, and the fluid's own.
            "kretschmann": ("\\dfrac{48f^2}{r^6} - \\dfrac{32\\left(k + 1\\right)\\left(2k + 1\\right)f\\,h}{\\left(2k - 1\\right)r^5}"
                            " + \\dfrac{16\\left(4k^4 + 4k^3 + 5k^2 + 1\\right)h^2}{\\left(2k - 1\\right)^2r^4}")}
           if system == "husain" else {}),
    }


def published_metric(vm, metric_id, system_id):
    """A live page's chart as (reader, matrix of its metric components, its coordinates)."""
    data = json.loads((vm.METRICS_DIR / f"{metric_id}.json").read_text(encoding="utf-8"))
    chart = next(c for c in data["coordinates"] if c["id"] == system_id)
    reader = vm.Reader(chart["coords"], [p["symbol"] for p in chart["parameters"]], ())
    values = {tuple(e["indices"]): e["value"] for e in chart["metric_components"]}
    n = len(chart["coords"])
    matrix = sp.Matrix(n, n, lambda i, j: reader(values.get((chart["coords"][i], chart["coords"][j]), "0")))
    return reader, matrix, chart["coords"]


# Points of the plane and a test mass function of both coordinates, at which every check is made:
# the mass function is no special case, so a chart that passes there passes for the family.
POINTS = [(sp.Rational(3, 10), sp.Rational(17, 10), sp.Rational(2, 5)),
          (sp.Rational(-11, 10), sp.Rational(29, 10), sp.Rational(9, 10)),
          (sp.Rational(7, 4), sp.Rational(13, 20), sp.Rational(6, 5))]


def check_chart(chart, system):
    """Each chart at points against what its sources state.

    The ingoing and outgoing charts have Wang and Wu's Einstein tensor, eq. (2), for a mass function
    of both coordinates; with m a function of the null time alone each is the published chart of the
    same name of Vaidya's metric at G = 1, Wang and Wu's statement that rho = P = 0 is Vaidya's
    solution; and m = Lambda r^3/6 is de Sitter space, G_ab = -Lambda g_ab, their case ii). Husain's
    chart has his rho = g/(4 pi r^(2k + 2)), P = k rho and his null flux, d_v m/(4 pi r^2) with m his
    eq. (12) and his g = q^(2k), at k = 3/4 and k = 1/3; with q = 0 it is Vaidya's ingoing chart, and at
    k = 1 it is the ingoing chart of Bonnor and Vaidya with their charge squared 2g = 2q^2."""
    import verify_metrics as vm
    name = f"{METRIC}/{system}"
    w, r, theta, phi = chart.symbols
    geo, g = chart.geo, chart.geo.g

    def vanishes(value, what, at):
        number = complex(sp.N(sp.sympify(value).subs(at), 30))
        if abs(number) > 1e-20:
            raise AssertionError(f"{name}: {what} at {at}, off by {number}")

    def mixed_einstein():
        return geo.raise_indices(geo.einstein_ll(), 2, (0,))

    if system == "husain":
        f, qq, k = (chart.reader.parameters[x] for x in ("f", "q", "k"))
        held = chart.reader.held
        mixed = mixed_einstein()
        for kk in (sp.Rational(3, 4), sp.Rational(1, 3)):
            F = sp.Rational(1, 2) + w / 7 + w ** 2 / 40
            Q = sp.Rational(1, 5) + sp.sin(w) ** 2 / 9
            Gf = Q ** (2 * kk)
            mass = F - Gf / ((2 * kk - 1) * r ** (2 * kk - 1))
            put = lambda e: sp.sympify(e).xreplace(held).doit().subs({f: F, qq: Q}).doit().subs(k, kk)  # noqa: E731
            rho = Gf / (4 * sp.pi * r ** (2 * kk + 2))
            for v0, r0, th0 in POINTS:
                at = {w: v0, r: r0, theta: th0}
                # G = 8 pi T: G^v_v = G^r_r = -8 pi rho, G^theta_theta = G^phi_phi = 8 pi k rho, Husain's
                # eqs. (5), (7) and (13), and G^r_v = 8 pi T^r_v = 2 d_v m/r^2, his eq. (3).
                want = {(0, 0): -8 * sp.pi * rho, (1, 1): -8 * sp.pi * rho, (2, 2): 8 * sp.pi * kk * rho,
                        (3, 3): 8 * sp.pi * kk * rho, (1, 0): 2 * sp.diff(mass, w) / r ** 2}
                for a in range(4):
                    for b in range(4):
                        vanishes(put(mixed[a][b]) - want.get((a, b), 0),
                                 f"Husain's null fluid at k = {kk} fails in slot {chart.coords_tex[a]}{chart.coords_tex[b]}", at)
                vanishes(put(g[0, 0]) + 1 - 2 * mass / r, "the lapse is not 1 - 2m/r with Husain's eq. (12)", at)
        # With q = 0, Vaidya's ingoing chart at G = 1 with m = f; at k = 1, Bonnor and Vaidya's with
        # m = f and their charge sqrt(2) q. Each published function is written in its own reader's v first.
        for metric_id, at_k, ours, theirs in (
                ("vaidya", sp.Rational(3, 4), {qq: 0}, lambda R, v: {R.parameters["G"]: 1, R.c: 1, R.parameters["m"]: 1 + v / 3}),
                ("bonnor_vaidya", 1, {qq: sp.Rational(1, 7) + w ** 2 / 9},
                 lambda R, v: {R.parameters["m"]: 1 + v / 3, R.parameters["q"]: sp.sqrt(2) * (sp.Rational(1, 7) + v ** 2 / 9)})):
            reader, published, coords = published_metric(vm, metric_id, "eddington_finkelstein_ingoing")
            there = published.subs(theirs(reader, reader.symbol["v"])).doit()
            there = there.subs(dict(zip((reader.symbol[x] for x in coords), chart.symbols)), simultaneous=True)
            here = g.xreplace(held).doit().subs(ours).subs({f: 1 + w / 3, k: at_k}).doit()
            for v0, r0, th0 in POINTS:
                at = {w: v0, r: r0, theta: th0}
                for a in range(4):
                    for b in range(4):
                        vanishes(here[a, b] - there[a, b], f"the chart misses the published {metric_id}", at)
        return

    m = chart.reader.parameters["m"]
    mixed = mixed_einstein()
    test = sp.Rational(1, 2) + w / 5 + r ** 2 / 9 + sp.sin(w) * r / 7 + w * r ** 3 / 31
    put = lambda e: sp.sympify(e).subs(m, test).doit()  # noqa: E731
    for v0, r0, th0 in POINTS:
        at = {w: v0, r: r0, theta: th0}
        want = {(0, 0): -2 * sp.diff(test, r) / r ** 2, (1, 1): -2 * sp.diff(test, r) / r ** 2,
                (2, 2): -sp.diff(test, r, 2) / r, (3, 3): -sp.diff(test, r, 2) / r,
                (1, 0): 2 * sp.diff(test, w) / r ** 2}
        for a in range(4):
            for b in range(4):
                vanishes(put(mixed[a][b]) - want.get((a, b), 0),
                         f"Wang and Wu's eq. (2) fails in slot {chart.coords_tex[a]}{chart.coords_tex[b]}", at)
    # With m a function of the null time alone, the published chart of Vaidya's metric at G = 1.
    reader, star, coords = published_metric(vm, "vaidya", "eddington_finkelstein_" + NAMES[system].lower())
    rv = reader.symbol[coords[0]]
    there = star.subs({reader.parameters["G"]: 1, reader.c: 1, reader.parameters["m"]: 1 + rv / 3 + sp.cos(rv) / 5}).doit()
    there = there.subs(dict(zip((reader.symbol[x] for x in coords), chart.symbols)), simultaneous=True)
    here = g.subs(m, 1 + w / 3 + sp.cos(w) / 5).doit()
    for v0, r0, th0 in POINTS:
        at = {w: v0, r: r0, theta: th0}
        for a in range(4):
            for b in range(4):
                vanishes(here[a, b] - there[a, b],
                         f"with m = m({chart.coords_tex[0]}) the chart misses the published vaidya", at)
    # Wang and Wu's case ii): m = Lambda r^3/6 is de Sitter space, G_ab = -Lambda g_ab.
    Lam = sp.Rational(3, 7)
    einstein = geo.einstein_ll()
    for v0, r0, th0 in POINTS:
        at = {w: v0, r: r0, theta: th0}
        for a in range(4):
            for b in range(4):
                vanishes((vm._at(einstein, (a, b)) + Lam * g[a, b]).subs(m, Lam * r ** 3 / 6).doit(),
                         "m = Lambda r^3/6 is not de Sitter space", at)
