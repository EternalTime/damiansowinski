"""The naked black holes of Horowitz and Ross, Phys. Rev. D 56, 2180 (1997), hep-th/9704058, and Phys.
Rev. D 57, 1098 (1998), hep-th/9709050, in units G = 1 with their Planck-unit forms read as lengths.

dilaton      HR97 eqs. (gensol), (uuf) and (R), ds^2 = -F dt^2 + dr^2/F + R^2 dOmega^2 with
             F = (r - r_+)(r - r_-)/R^2 and R = r(1 - r_-/r)^{a^2/(1+a^2)}, the dilaton black holes of
             Gibbons and Maeda and of Garfinkle, Horowitz and Strominger of the action (dbhac),
             R - 2(grad phi)^2 - e^{-2 a phi} F^2. R is held as a function of r with its slope
             R((1 + a^2) r - r_-)/((1 + a^2) r (r - r_-)), so every value is a rational function of
             r, R, a, r_+ and r_-: sympy does not finish factoring powers in 1/(1 + a^2), as the
             chan_horne_mann chart of dilaton_unusual_asymptotics found first.
two_charges  HR97 eqs. (uf), (ur) and (mcff), the same F with R^2 = r^2 - Sigma^2, the black hole of
             the action (uact), R - 2(grad phi)^2 - e^{-2 phi}(F^2 + G^2), with electric charge Q in F
             and magnetic charge P in G, Sigma = (P^2 - Q^2)/2M.
string_ns    HR97 eqs. (nsmet), (dd) and (nsdil), after Sen, the string metric of the black hole with
             winding and momentum charges, ds^2 = -Delta^{-1}(1 - r_0/r) dt^2 + dr^2/(1 - r_0/r)
             + r^2 dOmega^2 with Delta = (1 + r_0 sinh^2 gamma_1/r)(1 + r_0 sinh^2 gamma_p/r) and
             e^{2 phi} = Delta^{-1/2}, in string units alpha' = 1. It is written in the lengths
             r_1 = r_0 sinh^2 gamma_1 and r_p = r_0 sinh^2 gamma_p, so Delta = (r + r_1)(r + r_p)/r^2 and
             every value is rational in r: printed in gamma_1 and gamma_p, sympy wrote each sinh as
             exponentials and the chart ran to 390,000 characters.

Each chart is checked at points, in forty digits, independently of the printed tensors: its field
equations where its source gives the fields, the static tidal components (stcur) and the infalling one
(bhtide) of HR97 against the chart's own Riemann tensor in the boosted frame (inframea), HR97's (dbhtide),
(utide), (nsgeod) and (dd) at the horizon, the mass (dbhmass) and (mcff) from the fall of g_tt, and
each special case its source names against the live chart it is: a = 0 rn_metric, a = 1 the static
chart of dilaton_black_hole, a = sqrt(3) the Einstein chart of kaluza_klein_black_hole, r_- = 0
schwarzschild, Q = P rn_metric, one charge zero the a = 1 hole, and the Einstein metric e^{-2 phi} g
of the string chart the a = 1 hole for r_1 = r_p and the a = sqrt(3) hole for r_p = 0.
naked_black_holes.md beside this file records each chart's source and every drawing decision.
"""
import json
from pathlib import Path

import sympy as sp

METRIC = "naked_black_holes"
CHARTS = ("dilaton", "two_charges", "string_ns")
NAMES = {"dilaton": "Dilaton Coupling $a$", "two_charges": "Electric and Magnetic Charges",
         "string_ns": "String Metric, Winding and Momentum"}
COORDS = ["t", "r", "\\theta", "\\phi"]
SPHERE = "\\left(d\\theta^2 + \\sin^2\\theta\\,d\\phi^2\\right)"
LAPSE = "\\left(r - r_+\\right)\\left(r - r_-\\right)"
AREA = "\\left(r^2 - \\Sigma^2\\right)"
R_DEF = "R = r\\left(1 - \\dfrac{r_-}{r}\\right)^{a^2/(1 + a^2)}"
R_SLOPE = "\\dfrac{\\left(\\left(1 + a^2\\right)r - r_-\\right)R}{\\left(1 + a^2\\right)r\\left(r - r_-\\right)}"
PARAMETERS = {"dilaton": ["r_+", "r_-", "a", R_DEF], "two_charges": ["r_+", "r_-", "\\Sigma"],
              "string_ns": ["r_0", "r_1", "r_p"]}
HELD_NAMES = {"dilaton": ("R",)}
RATE_TEXTS = {"dilaton": {"R": {"r": R_SLOPE}}}


def line(system, c2):
    """The line element of a chart, with c2 standing before dt^2."""
    if system == "dilaton":
        return ("ds^2 = -\\dfrac{" + LAPSE + "}{R^2}" + c2 + "dt^2 + \\dfrac{R^2\\,dr^2}{" + LAPSE + "} + R^2" + SPHERE)
    if system == "two_charges":
        return ("ds^2 = -\\dfrac{" + LAPSE + "}{r^2 - \\Sigma^2}" + c2 + "dt^2 + \\dfrac{" + AREA + "dr^2}{" + LAPSE
                + "} + " + AREA + SPHERE)
    return ("ds^2 = -\\dfrac{r\\left(r - r_0\\right)}{\\left(r + r_1\\right)\\left(r + r_p\\right)}" + c2
            + "dt^2 + \\dfrac{r\\,dr^2}{r - r_0} + r^2" + SPHERE)


REALS = "t \\in (-\\infty, \\infty)"
ANGLES = ["\\theta \\in [0, \\pi]", "\\phi \\in [0, 2\\pi)"]
DOMAINS = {
    "dilaton": [REALS, "r \\in (r_-, \\infty)"] + ANGLES + [
        "0 \\le r_- < r_+,\\quad a > 0", "r = r_+ \\;\\text{(the event horizon)}", "r = r_- \\;\\text{(the singularity)}"],
    "two_charges": [REALS, "r \\in (|\\Sigma|, \\infty)"] + ANGLES + [
        "|\\Sigma| \\le r_- < r_+", "r = r_+ \\;\\text{(the event horizon)}", "r = r_- \\;\\text{(the inner horizon)}",
        "r = |\\Sigma| \\;\\text{(the singularity)}"],
    "string_ns": [REALS, "r \\in (0, \\infty)"] + ANGLES + [
        "r_0 > 0,\\quad r_1 \\ge 0,\\quad r_p \\ge 0", "r = r_0 \\;\\text{(the event horizon)}",
        "r = 0 \\;\\text{(the singularity)}"],
}
_DIMS = {"t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "r_+": "L", "r_-": "L", "a": "1", "R": "L",
         "\\Sigma": "L", "r_0": "L", "r_1": "L", "r_p": "L"}
DIMENSIONS = {(METRIC, s): {name: _DIMS[name] for name in COORDS + [p.split(" = ")[0] for p in PARAMETERS[s]]}
              for s in CHARTS}
HELD = {(METRIC, s): names for s, names in HELD_NAMES.items()}
RATES = {(METRIC, s): rates for s, rates in RATE_TEXTS.items()}
DIGITS = 40
TOLERANCE = sp.Float("1e-28")


def register(module):
    module.CHARTS[METRIC] = [lambda system=system: chart_spec(system) for system in CHARTS]


def reader(system):
    import verify_metrics as vm
    return vm.Reader(COORDS, PARAMETERS[system], (), held=HELD_NAMES.get(system, ()),
                     rates=RATE_TEXTS.get(system))


def chart_spec(system):
    """One of the three charts, its values printed in r and the chart's parameters, and the dilaton chart's in R."""
    import verify_metrics as vm

    probe = reader(system)
    r, P = probe.symbol["r"], probe.parameters
    spec = {"metric_id": METRIC,
            "system": {"id": system, "name": NAMES[system], "coords": COORDS, "domains": DOMAINS[system],
                       "parameters": PARAMETERS[system], "line_element": line(system, "c^2")},
            "chart_line_element": line(system, ""),
            "check": lambda chart: chart_check(chart, system)}
    if system == "dilaton":
        spec["printer"] = {"lead": [r, P["r_+"], P["r_-"], P["a"], P["R"]], "rising": [P["a"]], "flip": False}
        spec["reduce"] = lambda value: vm.norm(probe.by_rates(value))
        # R is the radius of the spheres here, so the Ricci scalar stands under its heading alone.
        spec["bare_scalar"] = True
    elif system == "two_charges":
        spec["printer"] = {"lead": [r, P["r_+"], P["r_-"], P["Sigma"]], "flip": False}
    else:
        spec["printer"] = {"lead": [r, P["r_0"], P["r_1"], P["r_p"]], "flip": False}
    return spec


# -- the checks ------------------------------------------------------------------------------------

def _numeric(chart, expr, at):
    """The value of a chart expression at a point: the derivatives of a held name written by its
    rates, the name written out, and the point put in, in forty digits."""
    expr = sp.sympify(expr)
    if chart.reader.rates:
        expr = chart.reader.by_rates(expr)
    if chart.reader.held:
        expr = expr.xreplace(chart.reader.held)
    return sp.N(expr.subs(at), DIGITS)


def _zero(value, what):
    if abs(complex(value)) > TOLERANCE:
        raise AssertionError(f"{METRIC}: {what} fails, by {value}")


def _frame_components(chart, at, E):
    """The static frame components R_0202, R_1212 and R_2323 of HR97's (stcur) and the infalling
    R_0'20'2 of their (bhtide) for a radial fall of energy E, all at the point, from the chart's own
    Riemann tensor. With g = diag(-A, B, C, C sin^2 theta), the static frame is e_0 = A^{-1/2} d_t,
    e_1 = B^{-1/2} d_r, e_2 = C^{-1/2} d_theta, and the fall has u^t = E/A, u^r = -(E^2 - A)^{1/2}/(AB)^{1/2}."""
    riemann = chart.geo.riemann_llll()
    g = chart.geo.g
    A, B, C = (-_numeric(chart, g[0, 0], at), _numeric(chart, g[1, 1], at), _numeric(chart, g[2, 2], at))
    comp = {key: _numeric(chart, riemann[key[0]][key[1]][key[2]][key[3]], at)
            for key in [(0, 2, 0, 2), (1, 2, 1, 2), (2, 3, 2, 3), (0, 2, 1, 2)]}
    sin2 = sp.sin(at[chart.reader.symbol["\\theta"]]) ** 2
    ut, ur = E / A, -sp.sqrt(E ** 2 - A) / sp.sqrt(A * B)
    static = {"0202": comp[(0, 2, 0, 2)] / (A * C), "1212": comp[(1, 2, 1, 2)] / (B * C),
              "2323": comp[(2, 3, 2, 3)] / (C ** 2 * sin2)}
    falling = (ut ** 2 * comp[(0, 2, 0, 2)] + 2 * ut * ur * comp[(0, 2, 1, 2)] + ur ** 2 * comp[(1, 2, 1, 2)]) / C
    return static, sp.N(falling, DIGITS)


def _live_metric(metric_id, chart_id):
    """The g of a live chart, read from its line element with its own reader."""
    import verify_metrics as vm
    there = next(c for c in json.loads((Path(vm.METRICS_DIR) / f"{metric_id}.json").read_text(encoding="utf-8"))
                 ["coordinates"] if c["id"] == chart_id)
    other = vm.Reader(there["coords"], [p["symbol"] for p in there["parameters"]], ())
    g = vm.metric_from_line_element(other, there["line_element"].replace("c^2", ""), there["coords"])
    return g, other


def _against_live(chart, at, metric_id, chart_id, radius, parameters, what):
    """The chart at the point equals a live chart whose r is the given function of this chart's r and
    whose parameters are the given values."""
    g_live, other = _live_metric(metric_id, chart_id)
    t, r, th, ph = chart.symbols
    names = dict(zip([other.symbol[c] for c in other.coords], [t, radius, th, ph]))
    names.update({other.parameters[k]: v for k, v in parameters.items()})
    pulled = g_live.subs(names, simultaneous=True)
    # dr' = (d radius / dr) dr, so g_rr carries its square.
    slope = sp.diff(radius, r)
    for i, scale in ((0, 1), (1, slope ** 2), (2, 1), (3, 1)):
        _zero(_numeric(chart, chart.geo.g[i, i] - scale * pulled[i, i], at), f"{what}, slot {i}{i}")


def _field_equations(chart, at, dilaton, fields):
    """G_ab = 2 d_a phi d_b phi - g_ab (d phi)^2 + 2 sum e^{-2 k phi}(F_ac F_b^c - g_ab F^2/4), the
    Maxwell equations d_a(sqrt(-g) e^{-2 k phi} F^ab) = 0 and the dilaton's 4 box phi
    + 2 sum k e^{-2 k phi} F^2 = 0, for the action R - 2(grad phi)^2 - sum e^{-2 k phi} F^2, at the point.
    fields is a list of (F_ab lowered, k)."""
    x = chart.symbols
    held = chart.reader.held
    g = sp.Matrix(chart.geo.g).xreplace(held)
    ginv = g.inv()
    dphi = [sp.diff(dilaton, c) for c in x]
    dphi2 = sum(ginv[a, b] * dphi[a] * dphi[b] for a in range(4) for b in range(4))
    root = sp.sqrt(-g.det())
    einstein = chart.geo.einstein_ll()
    G = sp.Matrix(4, 4, lambda a, b: vm_at(einstein, (a, b)))
    T = sp.Matrix(4, 4, lambda a, b: 2 * dphi[a] * dphi[b] - g[a, b] * dphi2)
    box = sum(sp.diff(root * ginv[a, b] * dphi[b], x[a]) for a in range(4) for b in range(4)) / root
    scalar = 4 * box
    for F, k in fields:
        Fup = ginv * F * ginv
        F2 = sum(F[a, b] * Fup[a, b] for a in range(4) for b in range(4))
        T += 2 * sp.exp(-2 * k * dilaton) * (F * ginv * F.T - g * F2 / 4)
        scalar += 2 * k * sp.exp(-2 * k * dilaton) * F2
        for b in range(4):
            maxwell = sum(sp.diff(root * sp.exp(-2 * k * dilaton) * Fup[a, b], x[a]) for a in range(4))
            _zero(sp.N(maxwell.subs(at), DIGITS), f"Maxwell's equation along {chart.coords_tex[b]}")
    for a in range(4):
        for b in range(a, 4):
            _zero(_numeric(chart, G[a, b], at) - sp.N(T[a, b].subs(at), DIGITS),
                  f"Einstein's equation in slot {chart.coords_tex[a]}{chart.coords_tex[b]}")
    _zero(sp.N(scalar.subs(at), DIGITS), "the dilaton's equation")


def vm_at(tensor, index):
    import verify_metrics as vm
    return vm._at(tensor, index)


def _mass(chart, at, mass, what):
    """-g_tt = 1 - 2M/r + O(1/r^2): r(1 + g_tt) at r = 10^30 is 2M to thirty digits."""
    r = chart.symbols[1]
    far = {**at, r: sp.Integer(10) ** 30}
    value = _numeric(chart, r * (1 + chart.geo.g[0, 0]), far)
    if abs(value - 2 * sp.N(mass.subs(at), DIGITS)) > sp.Float("1e-20"):
        raise AssertionError(f"{METRIC}: {what}: r(1 + g_tt) far away is {value}, not 2M")


def chart_check(chart, system):
    {"dilaton": _check_dilaton, "two_charges": _check_two_charges, "string_ns": _check_string}[system](chart)


def _check_dilaton(chart):
    t, r, th, ph = chart.symbols
    P = chart.reader.parameters
    rp, rm, a = P["r_+"], P["r_-"], P["a"]
    n = a ** 2 / (1 + a ** 2)
    R = r * (1 - rm / r) ** n
    F = (r - rp) * (r - rm) / R ** 2
    points = [{a: sp.Rational(7, 10), rp: 3, rm: 1, r: sp.Rational(9, 2), th: sp.Rational(6, 5), "E": sp.Rational(3, 2)},
              {a: sp.sqrt(3), rp: 5, rm: 2, r: 7, th: sp.Rational(1, 3), "E": 1},
              {a: 2, rp: 2, rm: sp.Rational(1, 2), r: sp.Rational(5, 2), th: 1, "E": sp.Rational(6, 5)}]
    clean = [{k: v for k, v in at.items() if k != "E"} for at in points]
    for at in clean:
        for i, value in enumerate((-F, 1 / F, R ** 2, R ** 2 * sp.sin(th) ** 2)):
            _zero(_numeric(chart, chart.geo.g[i, i] - value, at), f"HR97 eqs. (gensol), (uuf) and (R), slot {i}{i}")
        # The magnetic hole, F = Q sin(theta) d theta ^ d phi with Q^2 = r_+ r_-/(1 + a^2), HR97 (dbhmass),
        # and e^{-2 phi} = (1 - r_-/r)^{2a/(1 + a^2)}, HR98 (dbhdil).
        Q = sp.sqrt(rp * rm / (1 + a ** 2))
        Fm = sp.zeros(4, 4)
        Fm[2, 3], Fm[3, 2] = Q * sp.sin(th), -Q * sp.sin(th)
        _field_equations(chart, at, -(a / (1 + a ** 2)) * sp.log(1 - rm / r), [(Fm, a)])
        _mass(chart, at, rp / 2 + (1 - a ** 2) * rm / (2 * (1 + a ** 2)), "HR97 (dbhmass)")
    # The static frame (stcur) and the falling frame (bhtide), with R' and R'' from the chart's own slope.
    Rp, Rpp, Fp = sp.diff(R, r), sp.diff(R, r, 2), sp.diff(F, r)
    for at in points:
        place = {k: v for k, v in at.items() if k != "E"}
        static, falling = _frame_components(chart, place, at["E"])
        _zero(static["0202"] - sp.N((Fp * Rp / (2 * R)).subs(place), DIGITS), "HR97 (stcur), R_0202 = F'R'/(2R)")
        _zero(static["2323"] - sp.N(((1 - F * Rp ** 2) / R ** 2).subs(place), DIGITS),
              "R_2323 = (1 - F R'^2)/R^2, which is HR97 (stcur) at the horizon")
        E = at["E"]
        tide = -(Rpp * (E ** 2 - F) - Fp * Rp / 2) / R
        _zero(falling - sp.N(tide.subs(place), DIGITS), "HR97 (bhtide)")
        # R''/R at the horizon is HR97 (dbhtide), with epsilon = 1 - r_-/r_+.
        eps = 1 - rm / rp
        _zero(sp.N((Rpp / R).subs({**place, r: rp}), DIGITS)
              + sp.N((a ** 2 / (1 + a ** 2) ** 2 * (1 - eps) ** 2 / (rp ** 2 * eps ** 2)).subs(place), DIGITS),
              "HR97 (dbhtide)")
    # The special cases the Related entries name, each against its live chart.
    at = {rp: 3, rm: sp.Rational(5, 4), r: sp.Rational(17, 4), th: sp.Rational(4, 5)}
    _against_live(chart, {**at, a: 1}, "dilaton_black_hole", "static", r, {"r_s": rp, "r_d": rm},
                  "a = 1 is the static chart of dilaton_black_hole")
    _against_live(chart, {**at, a: sp.sqrt(3)}, "kaluza_klein_black_hole", "einstein", r - rm,
                  {"r_s": rp - rm, "q": rp}, "a = sqrt(3) is the Einstein chart of kaluza_klein_black_hole")
    _against_live(chart, {**at, a: 0}, "rn_metric", "spherical", r, _rn_parameters(rp, rm),
                  "a = 0 is Reissner and Nordstrom's hole")
    _against_live(chart, {**at, rm: 0, a: sp.Rational(3, 5)}, "schwarzschild", "spherical", r,
                  {"r_s": rp}, "r_- = 0 is Schwarzschild's hole")


def _rn_parameters(rp, rm):
    """Reissner and Nordstrom's 1 - r_s/r + r_q^2/r^2 is (r - r_+)(r - r_-)/r^2 at r_s = r_+ + r_-, r_q^2 = r_+ r_-."""
    return {"r_s": rp + rm, "r_q": sp.sqrt(rp * rm)}


def _check_two_charges(chart):
    t, r, th, ph = chart.symbols
    P = chart.reader.parameters
    rp, rm, S = P["r_+"], P["r_-"], P["Sigma"]
    R2 = r ** 2 - S ** 2
    F = (r - rp) * (r - rm) / R2
    points = [{rp: 4, rm: 3, S: 1, r: 6, th: sp.Rational(7, 5), "E": sp.Rational(3, 2)},
              {rp: sp.Rational(9, 2), rm: 2, S: -sp.Rational(3, 2), r: 5, th: sp.Rational(1, 2), "E": 1},
              {rp: 3, rm: sp.Rational(5, 2), S: sp.Rational(2, 5), r: sp.Rational(13, 4), th: 1, "E": sp.Rational(5, 4)}]
    M = (rp + rm) / 2
    for at in points:
        place = {k: v for k, v in at.items() if k != "E"}
        for i, value in enumerate((-F, 1 / F, R2, R2 * sp.sin(th) ** 2)):
            _zero(_numeric(chart, chart.geo.g[i, i] - value, place), f"HR97 eqs. (uf) and (ur), slot {i}{i}")
        # HR97 (mcff): M = (r_+ + r_-)/2, Q^2 + P^2 = r_+ r_- + Sigma^2 and P^2 - Q^2 = 2 M Sigma; the electric
        # field F_tr = -Q e^{2 phi}/R^2 and the magnetic G = P sin(theta) d theta ^ d phi, e^{2 phi} = (r + Sigma)/(r - Sigma).
        P2 = (rp * rm + S ** 2 + 2 * M * S) / 2
        Q2 = (rp * rm + S ** 2 - 2 * M * S) / 2
        dilaton = sp.log((r + S) / (r - S)) / 2
        Fe = sp.zeros(4, 4)
        Fe[0, 1] = -sp.sqrt(Q2) * sp.exp(2 * dilaton) / R2
        Fe[1, 0] = -Fe[0, 1]
        Gm = sp.zeros(4, 4)
        Gm[2, 3], Gm[3, 2] = sp.sqrt(P2) * sp.sin(th), -sp.sqrt(P2) * sp.sin(th)
        _field_equations(chart, place, dilaton, [(Fe, 1), (Gm, 1)])
        _mass(chart, place, M, "HR97 (mcff)")
        # r_+ and r_- are M +- (M^2 + Sigma^2 - Q^2 - P^2)^{1/2}, HR97 (mcff).
        root = sp.sqrt(M ** 2 + S ** 2 - Q2 - P2)
        _zero(sp.N((M + root - rp).subs(place), DIGITS), "HR97 (mcff), r_+")
        _zero(sp.N((M - root - rm).subs(place), DIGITS), "HR97 (mcff), r_-")
        static, falling = _frame_components(chart, place, at["E"])
        R = sp.sqrt(R2)
        Rp, Rpp, Fp = sp.diff(R, r), sp.diff(R, r, 2), sp.diff(F, r)
        _zero(static["0202"] - sp.N((Fp * Rp / (2 * R)).subs(place), DIGITS), "HR97 (stcur), R_0202 = F'R'/(2R)")
        E = at["E"]
        _zero(falling - sp.N((-(Rpp * (E ** 2 - F) - Fp * Rp / 2) / R).subs(place), DIGITS), "HR97 (bhtide)")
        # HR97 (utide): |R''/R| = (1 - delta)/(r_+^2 delta^2) at the horizon, delta = 1 - Sigma^2/r_+^2.
        delta = 1 - S ** 2 / rp ** 2
        _zero(sp.N((Rpp / R).subs({**place, r: rp}), DIGITS)
              + sp.N(((1 - delta) / (rp ** 2 * delta ** 2)).subs(place), DIGITS), "HR97 (utide)")
    at = {rp: 4, rm: 3, r: 6, th: sp.Rational(7, 5)}
    # Q = P: Sigma = 0 and Reissner and Nordstrom's hole.
    _against_live(chart, {**at, S: 0}, "rn_metric", "spherical", r, _rn_parameters(rp, rm),
                  "Q = P is Reissner and Nordstrom's hole")
    # One charge zero: |Sigma| = r_-, and the a = 1 hole with r -> r + Sigma, r_s = r_+ + r_-, r_d = 2 r_-.
    _against_live(chart, {**at, S: rm}, "dilaton_black_hole", "static", r + rm, {"r_s": rp + rm, "r_d": 2 * rm},
                  "one charge zero is the a = 1 hole")


def _check_string(chart):
    t, r, th, ph = chart.symbols
    P = chart.reader.parameters
    r0, r1, rp = P["r_0"], P["r_1"], P["r_p"]
    Delta = (1 + r1 / r) * (1 + rp / r)
    points = [{r0: 2, r1: sp.Rational(1, 2), rp: sp.Rational(13, 5), r: sp.Rational(7, 2), th: sp.Rational(3, 5), "E": 2},
              {r0: sp.Rational(3, 2), r1: sp.Rational(9, 4), rp: sp.Rational(1, 5), r: sp.Rational(5, 2), th: sp.Rational(5, 4),
               "E": 1}]
    for at in points:
        place = {k: v for k, v in at.items() if k != "E"}
        # HR97 (nsmet) and (dd), with sinh^2 gamma_1 = r_1/r_0 and sinh^2 gamma_p = r_p/r_0.
        s1, sp_ = sp.asinh(sp.sqrt(r1 / r0)), sp.asinh(sp.sqrt(rp / r0))
        source = (1 + r0 * sp.sinh(s1) ** 2 / r) * (1 + r0 * sp.sinh(sp_) ** 2 / r)
        for i, value in enumerate((-(1 - r0 / r) / source, 1 / (1 - r0 / r), r ** 2, r ** 2 * sp.sin(th) ** 2)):
            _zero(_numeric(chart, chart.geo.g[i, i] - value, place), f"HR97 eqs. (nsmet) and (dd), slot {i}{i}")
        # The falling frame of HR97 (gdevb) with F = 1 - r_0/r, G = Delta and R = r, so R' = 1 and R'' = 0,
        # and its value at r_0, (nsgeod), with Delta'(r_0) of (dd).
        static, falling = _frame_components(chart, place, at["E"])
        E = at["E"]
        tide = -(E ** 2 * sp.diff(Delta, r) - sp.diff(1 - r0 / r, r)) / (2 * r)
        _zero(falling - sp.N(tide.subs(place), DIGITS), "HR97 (gdevb) with G = Delta and R = r")
        slope = -(sp.sinh(s1) ** 2 * sp.cosh(sp_) ** 2 + sp.cosh(s1) ** 2 * sp.sinh(sp_) ** 2) / r0
        _zero(sp.N((sp.diff(Delta, r).subs(r, r0) - slope).subs(place), DIGITS), "HR97 (dd), Delta'(r_0)")
        _zero(sp.N((tide.subs(r, r0) - (-E ** 2 * slope / (2 * r0) + 1 / (2 * r0 ** 2))).subs(place), DIGITS),
              "HR97 (nsgeod)")
    # The Einstein metric e^{-2 phi} g = Delta^{1/2} g: for r_1 = r_p the a = 1 hole with r -> r + r_1,
    # r_s = r_0 + r_1 = r_0 cosh^2 gamma and r_d = r_1; for r_p = 0 the a = sqrt(3) hole with r -> r,
    # r_s = r_0 and q = r_0 + r_1, the Einstein chart of kaluza_klein_black_hole.
    at = {r0: 2, r: sp.Rational(7, 2), th: sp.Rational(3, 5)}
    k = sp.Rational(4, 5)
    _einstein_against_live(chart, {**at, r1: k, rp: k}, Delta, "dilaton_black_hole", "static", r + k,
                           {"r_s": 2 + k, "r_d": k}, "r_1 = r_p is the a = 1 hole in the Einstein metric")
    _einstein_against_live(chart, {**at, r1: k, rp: 0}, Delta, "kaluza_klein_black_hole", "einstein", r,
                           {"r_s": 2, "q": 2 + k}, "r_p = 0 is the a = sqrt(3) hole in the Einstein metric")


def _einstein_against_live(chart, at, Delta, metric_id, chart_id, radius, parameters, what):
    """The Einstein metric Delta^{1/2} g of the string chart at the point equals a live chart whose r is
    the given function of this chart's r."""
    g_live, other = _live_metric(metric_id, chart_id)
    t, r, th, ph = chart.symbols
    names = dict(zip([other.symbol[c] for c in other.coords], [t, radius, th, ph]))
    names.update({other.parameters[k]: v for k, v in parameters.items()})
    pulled = g_live.subs(names, simultaneous=True)
    slope = sp.diff(radius, r)
    for i, scale in ((0, 1), (1, slope ** 2), (2, 1), (3, 1)):
        _zero(_numeric(chart, sp.sqrt(Delta) * chart.geo.g[i, i] - scale * pulled[i, i], at), f"{what}, slot {i}{i}")
