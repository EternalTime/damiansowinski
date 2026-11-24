"""The static black holes of anti-de Sitter space in five dimensions with three charges, of the STU model
of N = 2 gauged supergravity: Behrndt, Cvetic and Sabra, Nucl. Phys. B 553, 317 (1999), hep-th/9810227.

Every chart is in the signature (-,+,+,+,+), with the charge parameters q_i and the nonextremality
parameter mu lengths squared and the radius L of anti-de Sitter space a length, L = 1/g:

    sphere      Behrndt, Cvetic and Sabra's solution, their eq. (390) with the ansatz (110) and (170), and
                Cvetic and nine coauthors' eq. (d5adsbh), hep-th/9903214, with k = 1:
                ds^2 = -(H_1H_2H_3)^{-2/3} f dt^2 + (H_1H_2H_3)^{1/3}(dr^2/f + r^2 dOmega_3^2),
                f = k - mu/r^2 + r^2 H_1H_2H_3/L^2 and H_i = 1 + q_i/r^2;
    flat        the same with k = 0 and the flat horizon r^2(dx^2 + dy^2 + dz^2)/L^2, Cvetic and nine
                coauthors' eq. (d5adsbhk0) and Son and Starinets's section 2, hep-th/0601157;
    hyperbolic  the same with k = -1 and the hyperbolic 3-space d chi^2 + sinh^2 chi dOmega_2^2;
    continued   the sphere chart in Myers and Tafjord's radius rho^2 = r^2 + q_3, their eq. (rhoeq),
                hep-th/0109127, for q_1 >= q_2 >= q_3 > 0, which reaches past r = 0 to the singularity
                at rho = 0: with P_1 = r^2 H_1, P_2 = r^2 H_2 and D = r^4 f written in rho, the metric is
                -D/(rho^2 P_1 P_2)^{2/3} dt^2 + rho^2 (rho^2 P_1 P_2)^{1/3} d rho^2/D
                + (rho^2 P_1 P_2)^{1/3} dOmega_3^2;
    brane       the flat chart in Son and Starinets's u = r_+^2/r^2, their eqs. (metric_u_3) and (identif),
                and Gladden, Ivo, Kovtun and Starinets's eq. (STU-metric), 2412.12353, with
                kappa_i = q_i/r_+^2, H_i = 1 + kappa_i u, f = H_1H_2H_3 - u^2 (1 + kappa_1)(1 + kappa_2)
                (1 + kappa_3), the horizon at u = 1 and the boundary at u = 0; their T_0 = r_+/(pi L^2)
                is written out, so (pi T_0 L)^2 = r_+^2/L^2.

Every chart holds its harmonic functions and its metric function as names with declared slopes, so
each printed value is a rational function of them, the radius and the parameters, as the charts of
string_bh_three_four_charges are. stu_black_holes.md beside this file records each chart's source and
every drawing decision.
"""
import json
from pathlib import Path

import sympy as sp

METRIC = "stu_black_holes"
CHARTS = ("sphere", "flat", "hyperbolic", "continued", "brane")
NAMES = {"sphere": "Spherical Horizon, $k = 1$",
         "flat": "Flat Horizon, $k = 0$",
         "hyperbolic": "Hyperbolic Horizon, $k = -1$",
         "continued": "Past $r = 0$, Radius $\\rho$",
         "brane": "Black Brane, Coordinate $u$"}
CURVATURE = {"sphere": 1, "flat": 0, "hyperbolic": -1}

SPHERE_3 = "\\left(d\\chi^2 + \\sin^2\\chi\\left(d\\theta^2 + \\sin^2\\theta\\,d\\phi^2\\right)\\right)"
HYPERBOLIC_3 = "\\left(d\\chi^2 + \\sinh^2\\chi\\left(d\\theta^2 + \\sin^2\\theta\\,d\\phi^2\\right)\\right)"
PLANE_3 = "\\left(dx^2 + dy^2 + dz^2\\right)"
PRODUCT = "H_1\\,H_2\\,H_3"
PAIRS = "H_1\\,H_2 + H_1\\,H_3 + H_2\\,H_3"
CHARGES = ["\\mu", "q_1", "q_2", "q_3", "L"]


def _harmonic(i, variable="r"):
    if variable == "u":
        return f"H_{i} = 1 + \\kappa_{i}u"
    return f"H_{i} = 1 + \\dfrac{{q_{i}}}{{r^2}}"


def _metric_function(k):
    lead = {1: "1", 0: "", -1: "-1"}[k]
    mu = " - \\dfrac{\\mu}{r^2}" if lead else "-\\dfrac{\\mu}{r^2}"
    return "f = " + lead + mu + " + \\dfrac{r^2}{L^2}" + PRODUCT


def _slope_f(k):
    """f' = 2(k - f)/r + (2r/L^2)(H_1H_2 + H_1H_3 + H_2H_3 - H_1H_2H_3), with H_i' = -2(H_i - 1)/r."""
    lead = {1: "\\dfrac{2\\left(1 - f\\right)}{r}", 0: "-\\dfrac{2f}{r}", -1: "-\\dfrac{2\\left(1 + f\\right)}{r}"}[k]
    return lead + " + \\dfrac{2r}{L^2}\\left(" + PAIRS + " - " + PRODUCT + "\\right)"


P_DEFS = ["P_1 = \\rho^2 + q_1 - q_3", "P_2 = \\rho^2 + q_2 - q_3"]
D_DEF = ("D = \\left(\\rho^2 - q_3\\right)^2 - \\mu\\left(\\rho^2 - q_3\\right)"
         " + \\dfrac{\\rho^2\\left(\\rho^2 + q_1 - q_3\\right)\\left(\\rho^2 + q_2 - q_3\\right)}{L^2}")
D_SLOPE = ("2\\rho\\left(2\\rho^2 - 2q_3 - \\mu\\right) + \\dfrac{2\\rho}{L^2}\\left(P_1\\,P_2 + \\rho^2P_1"
           " + \\rho^2P_2\\right)")
BRANE_F = "f = H_1\\,H_2\\,H_3 - u^2\\left(1 + \\kappa_1\\right)\\left(1 + \\kappa_2\\right)\\left(1 + \\kappa_3\\right)"
# With H_i = 1 + kappa_i u, H_i' = (H_i - 1)/u, and u^2 (1 + kappa_1)(1 + kappa_2)(1 + kappa_3) = H_1H_2H_3 - f, so
# f' = (H_1H_2H_3 - H_1H_2 - H_1H_3 - H_2H_3 + 2f)/u: every value is rational in u, f and the H_i alone.
BRANE_SLOPE = "\\dfrac{H_1\\,H_2\\,H_3 - H_1\\,H_2 - H_1\\,H_3 - H_2\\,H_3 + 2f}{u}"

COORDS = {"sphere": ["t", "r", "\\chi", "\\theta", "\\phi"],
          "flat": ["t", "r", "x", "y", "z"],
          "hyperbolic": ["t", "r", "\\chi", "\\theta", "\\phi"],
          "continued": ["t", "\\rho", "\\chi", "\\theta", "\\phi"],
          "brane": ["t", "x", "y", "z", "u"]}
PARAMETERS = {
    **{s: CHARGES + [_harmonic(i) for i in (1, 2, 3)] + [_metric_function(k)] for s, k in CURVATURE.items()},
    "continued": CHARGES + P_DEFS + [D_DEF],
    "brane": ["r_+", "L", "\\kappa_1", "\\kappa_2", "\\kappa_3"] + [_harmonic(i, "u") for i in (1, 2, 3)] + [BRANE_F],
}
HELD_NAMES = {**{s: ("f", "H_1", "H_2", "H_3") for s in CURVATURE}, "continued": ("P_1", "P_2", "D"),
              "brane": ("f", "H_1", "H_2", "H_3")}
RATE_TEXTS = {
    **{s: {"f": {"r": _slope_f(k)},
           **{f"H_{i}": {"r": f"-\\dfrac{{2\\left(H_{i} - 1\\right)}}{{r}}"} for i in (1, 2, 3)}}
       for s, k in CURVATURE.items()},
    "continued": {"P_1": {"\\rho": "2\\rho"}, "P_2": {"\\rho": "2\\rho"}, "D": {"\\rho": D_SLOPE}},
    "brane": {"f": {"u": BRANE_SLOPE}, **{f"H_{i}": {"u": f"\\dfrac{{H_{i} - 1}}{{u}}"} for i in (1, 2, 3)}},
}


def line(system, c2):
    """The line element of a chart, with c2 standing before dt^2."""
    if system in CURVATURE:
        space = {"sphere": "r^2" + SPHERE_3, "flat": "\\dfrac{r^2}{L^2}" + PLANE_3,
                 "hyperbolic": "r^2" + HYPERBOLIC_3}[system]
        return ("ds^2 = -\\left(" + PRODUCT + "\\right)^{-2/3}f\\," + c2 + "dt^2 + \\left(" + PRODUCT
                + "\\right)^{1/3}\\left(\\dfrac{dr^2}{f} + " + space + "\\right)")
    if system == "continued":
        area = "\\left(P_1\\,P_2\\right)"
        return ("ds^2 = -\\dfrac{D}{\\rho^{4/3}" + area + "^{2/3}}" + c2 + "dt^2 + \\dfrac{\\rho^{8/3}" + area
                + "^{1/3}}{D}d\\rho^2 + \\rho^{2/3}" + area + "^{1/3}" + SPHERE_3)
    return ("ds^2 = -\\dfrac{r_+^2f}{L^2u\\left(" + PRODUCT + "\\right)^{2/3}}" + c2 + "dt^2 + \\dfrac{r_+^2\\left("
            + PRODUCT + "\\right)^{1/3}}{L^2u}" + PLANE_3 + " + \\dfrac{L^2\\left(" + PRODUCT
            + "\\right)^{1/3}}{4u^2f}du^2")


REALS = " \\in (-\\infty, \\infty)"
ANGLES = ["\\theta \\in [0, \\pi]", "\\phi \\in [0, 2\\pi)"]
CHARGE_RANGE = "\\mu \\ge 0,\\quad q_1, q_2, q_3 \\ge 0"
DOMAINS = {
    "sphere": ["t" + REALS, "r \\in (0, \\infty)", "\\chi \\in [0, \\pi]"] + ANGLES
    + [CHARGE_RANGE, "f = 0 \\;\\text{(the horizons)}"],
    "flat": ["t" + REALS, "r \\in (0, \\infty)", "x" + REALS, "y" + REALS, "z" + REALS]
    + [CHARGE_RANGE, "f = 0 \\;\\text{(the horizons)}"],
    "hyperbolic": ["t" + REALS, "r \\in (0, \\infty)", "\\chi \\in [0, \\infty)"] + ANGLES
    + ["\\mu \\ge q_1, q_2, q_3 \\ge 0", "f = 0 \\;\\text{(the horizons)}"],
    "continued": ["t" + REALS, "\\rho \\in (0, \\infty)", "\\chi \\in [0, \\pi]"] + ANGLES
    + ["\\mu \\ge 0,\\quad q_1 \\ge q_2 \\ge q_3 > 0", "D = 0 \\;\\text{(the horizons)}",
       "\\rho = \\sqrt{q_3} \\;\\text{(the radius } r = 0\\text{)}", "\\rho = 0 \\;\\text{(the singularity)}"],
    "brane": ["t" + REALS, "x" + REALS, "y" + REALS, "z" + REALS, "u \\in (0, \\infty)",
              "r_+ > 0,\\quad \\kappa_1, \\kappa_2, \\kappa_3 \\ge 0", "u = 1 \\;\\text{(the event horizon)}",
              "u = 0 \\;\\text{(the boundary of anti-de Sitter space)}"],
}
_DIMS = {"t": "T", "r": "L", "\\rho": "L", "u": "1", "x": "L", "y": "L", "z": "L", "\\chi": "1", "\\theta": "1",
         "\\phi": "1", "\\mu": "L^2", "q_1": "L^2", "q_2": "L^2", "q_3": "L^2", "L": "L", "r_+": "L",
         "\\kappa_1": "1", "\\kappa_2": "1", "\\kappa_3": "1", "H_1": "1", "H_2": "1", "H_3": "1", "f": "1",
         "P_1": "L^2", "P_2": "L^2", "D": "L^4"}
DIMENSIONS = {(METRIC, s): {name: _DIMS[name] for name in COORDS[s] + [p.split(" = ")[0] for p in PARAMETERS[s]]}
              for s in CHARTS}
HELD = {(METRIC, s): names for s, names in HELD_NAMES.items()}
RATES = {(METRIC, s): rates for s, rates in RATE_TEXTS.items()}
DIGITS = 40
TOLERANCE = sp.Float("1e-26")


def register(module):
    module.CHARTS[METRIC] = [lambda system=system: chart_spec(system) for system in CHARTS]


def reader(system):
    import verify_metrics as vm
    return vm.Reader(COORDS[system], PARAMETERS[system], (), held=HELD_NAMES[system], rates=RATE_TEXTS[system])


def chart_spec(system):
    """One of the five charts, its values printed in its held names, its radius and its parameters."""
    import verify_metrics as vm

    probe = reader(system)
    P = probe.parameters
    radius = probe.symbol[COORDS[system][1] if system != "brane" else "u"]
    if system == "continued":
        held = [P["P_1"], P["P_2"], P["D"]]
        lead = held + [radius, P["q_3"], P["mu"], P["L"]]
        factors = [radius, P["D"], P["P_1"], P["P_2"]]
    else:
        H = [P[f"H_{i}"] for i in (1, 2, 3)]
        lead = H + [P["f"], radius, P["L"]]
        factors = [radius, P["f"]] + H
    return {"metric_id": METRIC,
            "system": {"id": system, "name": NAMES[system], "coords": COORDS[system], "domains": DOMAINS[system],
                       "parameters": PARAMETERS[system], "line_element": line(system, "c^2")},
            "chart_line_element": line(system, ""),
            "printer": {"lead": lead, "factors": factors, "flip": False},
            "reduce": lambda value: vm.norm(probe.by_rates(value)),
            "check": lambda chart: chart_check(chart, system)}


# -- the checks ------------------------------------------------------------------------------------
#
# Each chart is checked at events in forty digits, independently of the printed tensors.
#
# The field equations are those of Cvetic and nine coauthors' Lagrangian (d5gauged),
# e^{-1} L = R - (1/2)(d phi_1)^2 - (1/2)(d phi_2)^2 + 4 g^2 sum X_i^{-1} - (1/4) sum X_i^{-2} (F^i)^2 + CS,
# with g = 1/L, X_i = exp(-a_i . phi/2), a_i . a_j = 4 delta_ij - 4/3, so sum (d ln X_i)^2 = (d phi)^2, and the
# Chern-Simons term silent for electric fields. Slot by slot:
#   G_ab = (1/2) sum (d_a ln X_i d_b ln X_i - g_ab (d ln X_i)^2/2) + (1/2) sum X_i^{-2}(F_ac F_b^c - g_ab F^2/4)
#          + (2/L^2) g_ab sum X_i^{-1},
#   d_a(sqrt(-g) X_i^{-2} F_i^ab) = 0,
#   box ln X_j = -X_j^{-2} F_j^2/2 + sum X_i^{-2} F_i^2/6 + (4/L^2)(X_j^{-1} - sum X_i^{-1}/3),
# with X_i = (H_1H_2H_3)^{1/3}/H_i and A^i = Q_i c dt/(r^2 + q_i), Q_i^2 = q_i(mu + k q_i): Behrndt, Cvetic and
# Sabra's F_r0 = -(sqrt k/2) H_I^{-2} d_r Htilde_I with their (392), qtilde_I = (mu/k) sinh beta_I cosh beta_I
# and q_I = (mu/k) sinh^2 beta_I, and Cvetic and nine coauthors' A^i = sqrt k (1 - H_i^{-1}) coth beta_i dt.

CASES = [{"\\mu": sp.Rational(13, 10), "q_1": sp.Rational(7, 10), "q_2": sp.Rational(2, 5), "q_3": sp.Rational(1, 5),
          "L": sp.Rational(3, 2)},
         {"\\mu": sp.Rational(21, 10), "q_1": sp.Rational(9, 10), "q_2": sp.Rational(9, 10), "q_3": sp.Rational(3, 10),
          "L": sp.Rational(4, 5)}]
RADII = [sp.Rational(13, 10), sp.Rational(21, 10)]
INNER = [sp.Rational(3, 10), sp.Rational(17, 10)]       # rho: inside rho^2 = q_3 for both cases, and outside
ANGLE_POINTS = [(sp.Rational(3, 5), sp.Rational(11, 10), sp.Rational(2, 5)),
                (sp.Rational(7, 5), sp.Rational(1, 3), sp.Rational(5, 2))]
BRANE_CASES = [{"r_+": sp.Rational(6, 5), "L": sp.Rational(3, 2), "kappa_1": sp.Rational(1, 2),
                "kappa_2": sp.Rational(1, 5), "kappa_3": sp.Rational(1, 10)},
               {"r_+": sp.Rational(9, 10), "L": sp.Rational(7, 10), "kappa_1": sp.Rational(3, 2),
                "kappa_2": sp.Rational(3, 2), "kappa_3": sp.Rational(0)}]
BRANE_U = [sp.Rational(2, 5), sp.Rational(17, 10)]


def _plain(name):
    return name.replace("\\", "")


def _zero(value, what):
    if abs(complex(value)) > TOLERANCE:
        raise AssertionError(f"{METRIC}: {what} fails, by {complex(value)}")


def _explicit(chart):
    """The metric with every held name written out."""
    return sp.Matrix(chart.geo.g).xreplace(chart.reader.held).doit()


_READERS = {}


def _numeric(system, chart, expr, at):
    """A chart expression at an event: the derivatives of the held names by their declared rates, the
    names written out, the event put in, in forty digits."""
    if system not in _READERS:
        _READERS[system] = reader(system)
    expr = _READERS[system].by_rates(sp.sympify(expr)).xreplace(chart.reader.held).doit()
    return sp.N(expr.subs(at), DIGITS)


def _events(chart, system):
    X = chart.symbols
    P = chart.reader.parameters
    if system == "brane":
        for case in BRANE_CASES:
            for u in BRANE_U:
                at = {P[_plain(n)] if _plain(n) in P else P[n]: v for n, v in case.items()}
                at.update({X[0]: sp.Rational(1, 3), X[1]: sp.Rational(1, 7), X[2]: sp.Rational(2, 9),
                           X[3]: sp.Rational(5, 11), X[4]: u})
                yield at
        return
    radii = INNER if system == "continued" else RADII
    for case in CASES:
        for radius in radii:
            for angles in ANGLE_POINTS:
                at = {P[_plain(n)]: v for n, v in case.items()}
                at.update({X[0]: sp.Rational(1, 3), X[1]: radius})
                if system == "flat":
                    at.update(dict(zip(X[2:], [sp.Rational(1, 7), sp.Rational(-2, 9), sp.Rational(5, 11)])))
                else:
                    at.update(dict(zip(X[2:], angles)))
                yield at


def _fields(chart, system):
    """(ln X_i, A^i lowered) for the three charges in the chart's own coordinates, x^0 = ct."""
    X = chart.symbols
    P = chart.reader.parameters
    if system == "brane":
        u = X[4]
        rp, L = P["r_+"], P["L"]
        kappa = [P[f"kappa_{i}"] for i in (1, 2, 3)]
        H = [1 + k * u for k in kappa]
        mu = rp ** 4 * sp.Mul(*[1 + k for k in kappa]) / L ** 2
        lnX = [sp.log(sp.Mul(*H)) / 3 - sp.log(h) for h in H]
        A = [sp.sqrt(k * rp ** 2 * mu) * u / (rp ** 2 * h) for k, h in zip(kappa, H)]
        return lnX, A
    mu = P["mu"]
    q = [P[f"q_{i}"] for i in (1, 2, 3)]
    if system == "continued":
        rho = X[1]
        Pi = [rho ** 2 + q[0] - q[2], rho ** 2 + q[1] - q[2], rho ** 2]      # r^2 H_i written in rho
        area = sp.Mul(*Pi)
        lnX = [sp.log(area) / 3 - sp.log(p) for p in Pi]
        A = [sp.sqrt(qi * (mu + qi)) / p for qi, p in zip(q, Pi)]
        return lnX, A
    k = CURVATURE[system]
    r = X[1]
    H = [1 + qi / r ** 2 for qi in q]
    lnX = [sp.log(sp.Mul(*H)) / 3 - sp.log(h) for h in H]
    A = [sp.sqrt(qi * (mu + k * qi)) / (r ** 2 + qi) for qi in q]
    return lnX, A


def _field_equations(chart, system):
    """The Einstein, Maxwell and scalar equations of the gauged STU model at every event."""
    import verify_metrics as vm
    X = chart.symbols
    L = chart.reader.parameters["L"]
    g = _explicit(chart)
    ginv = sp.diag(*[1 / g[a, a] for a in range(5)])
    root = sp.sqrt(-sp.Mul(*[g[a, a] for a in range(5)]))
    lnX, A = _fields(chart, system)
    Xinv = [sp.exp(-x) for x in lnX]
    weights = [w ** 2 for w in Xinv]                                         # X_i^{-2}
    dln = [[sp.diff(x, c) for c in X] for x in lnX]
    Fs = []
    for Ai in A:
        F = sp.zeros(5, 5)
        F[1 if system != "brane" else 4, 0] = sp.diff(Ai, X[1 if system != "brane" else 4])
        F[0, 1 if system != "brane" else 4] = -F[1 if system != "brane" else 4, 0]
        Fs.append(F)
    ups = [ginv * F * ginv for F in Fs]
    squares = [sum(F[a, b] * U[a, b] for a in range(5) for b in range(5)) for F, U in zip(Fs, ups)]
    einstein = chart.geo.einstein_ll()
    potential = 2 * sum(Xinv) / L ** 2
    for at in _events(chart, system):
        for a in range(5):
            for b in range(a, 5):
                stress = potential * g[a, b]
                for i in range(3):
                    grad2 = sum(ginv[c, c] * dln[i][c] ** 2 for c in range(5))
                    FF = sum(Fs[i][a, c] * Fs[i][b, c] * ginv[c, c] for c in range(5))
                    stress += (dln[i][a] * dln[i][b] - g[a, b] * grad2 / 2) / 2
                    stress += weights[i] * (FF - g[a, b] * squares[i] / 4) / 2
                _zero(_numeric(system, chart, vm._at(einstein, (a, b)), at) - sp.N(stress.subs(at), DIGITS),
                      f"{system}: Einstein's equation in slot {chart.coords_tex[a]}{chart.coords_tex[b]}")
        for i in range(3):
            for b in range(5):
                maxwell = sum(sp.diff(root * weights[i] * ups[i][a, b], X[a]) for a in range(5))
                _zero(sp.N(maxwell.subs(at), DIGITS), f"{system}: Maxwell's equation {b} of the field {i + 1}")
            box = sum(sp.diff(root * ginv[a, a] * dln[i][a], X[a]) for a in range(5)) / root
            scalar = (box + weights[i] * squares[i] / 2 - sum(w * s for w, s in zip(weights, squares)) / 6
                      - 4 * (Xinv[i] - sum(Xinv) / 3) / L ** 2)
            _zero(sp.N(scalar.subs(at), DIGITS), f"{system}: the equation of the scalar X_{i + 1}")


def _source(system):
    """Another chart of this spacetime as print_charts builds it."""
    import chart_printer as cp
    spec = chart_spec(system)
    return cp.Chart(spec["system"]["coords"], spec["system"]["parameters"], spec["chart_line_element"],
                    held=HELD_NAMES[system])


def _pulled(chart, system, source, image, names, events, what):
    """This chart's metric against `source`'s metric pulled back along the map `image` of this chart's
    coordinates, with `names` taking the source's parameters to expressions in this chart's."""
    X = chart.symbols
    J = sp.Matrix(5, 5, lambda a, b: sp.diff(image[a], X[b]))
    theirs = sp.Matrix(source.geo.g).xreplace(source.reader.held).doit()
    at = dict(zip(source.symbols, image))
    at.update(names)
    pulled = J.T * theirs.subs(at, simultaneous=True) * J
    g = _explicit(chart)
    for e in events:
        for a in range(5):
            for b in range(a, 5):
                _zero(sp.N((pulled[a, b] - g[a, b]).subs(e), DIGITS),
                      f"{system}: {what}, slot {chart.coords_tex[a]}{chart.coords_tex[b]}")


def _live(metric_id, chart_id):
    import verify_metrics as vm
    path = Path(vm.METRICS_DIR) / f"{metric_id}.json"
    there = next(c for c in json.loads(path.read_text(encoding="utf-8"))["coordinates"] if c["id"] == chart_id)
    return there


def _against_live(chart, system, metric_id, chart_id, values, limit, what):
    """This chart in the limit `limit` (a substitution) against a live chart of the same coordinates,
    whose parameters take the expressions `values` of this chart's."""
    import verify_metrics as vm
    there = _live(metric_id, chart_id)
    symbols = [p["symbol"] for p in there["parameters"]]
    key = (metric_id, chart_id)
    other = vm.Reader(there["coords"], symbols, (), held=vm.HELD.get(key, ()), rates=vm.RATES.get(key))
    theirs = vm.metric_from_line_element(other, there["line_element"].replace("c^2", ""), there["coords"])
    theirs = sp.Matrix(theirs).xreplace(other.held).doit()
    at = dict(zip([other.symbol[c] for c in other.coords], chart.symbols))
    at.update({other.parameters[_plain(k)]: v for k, v in values.items()})
    theirs = theirs.subs(at, simultaneous=True)
    mine = _explicit(chart).subs(limit)
    for e in _events(chart, system):
        for a in range(5):
            _zero(sp.N((theirs[a, a] - mine[a, a]).subs(e), DIGITS), f"{system}: {what}, slot {a}{a}")


def chart_check(chart, system):
    """Every chart against the field equations, and each against the charts it is pulled back from
    or reduces to.

    sphere      at L -> infinity the live five_charges chart of string_bh_three_four_charges, Horowitz,
                Maldacena and Strominger's hole, with r_0^2 = mu and r_i^2 = q_i, which Behrndt, Cvetic and
                Sabra name as their ungauged limit; with no charge and L -> infinity the live spherical
                chart of tangherlini, r_h^2 = mu.
    flat, hyperbolic  their field equations.
    continued   the sphere chart pulled back along r^2 = rho^2 - q_3, where rho^2 > q_3, and its field
                equations on both sides of rho^2 = q_3; with the three charges equal it is the
                Reissner-Nordstrom-anti-de Sitter black hole of five dimensions in the areal radius,
                -V dt^2 + d rho^2/V + rho^2 dOmega_3^2 with V = 1 - (2q + mu)/rho^2 + q(q + mu)/rho^4
                + rho^2/L^2, Myers and Tafjord's statement that the scalars are then constant.
    brane       the flat chart pulled back along r = r_+/sqrt u, with q_i = kappa_i r_+^2 and
                mu = r_+^4 (1 + kappa_1)(1 + kappa_2)(1 + kappa_3)/L^2, f(1) = 0, and Son and Starinets's
                Hawking temperature T_H = T_0 (2 + kappa_1 + kappa_2 + kappa_3 - kappa_1 kappa_2 kappa_3)
                /(2 sqrt((1 + kappa_1)(1 + kappa_2)(1 + kappa_3))) with T_0 = r_+/(pi L^2), from the
                surface gravity of the chart in units hbar = k_B = c = 1."""
    _field_equations(chart, system)
    X = chart.symbols
    P = chart.reader.parameters
    if system == "sphere":
        L = P["L"]
        _against_live(chart, system, "string_bh_three_four_charges", "five_charges",
                      {"r_0": sp.sqrt(P["mu"]), "r_1": sp.sqrt(P["q_1"]), "r_2": sp.sqrt(P["q_2"]),
                       "r_3": sp.sqrt(P["q_3"])}, {L: sp.oo},
                      "the hole of flat space with three charges at L -> infinity")
        _against_live(chart, system, "tangherlini", "spherical", {"r_h": sp.sqrt(P["mu"])},
                      {L: sp.oo, P["q_1"]: 0, P["q_2"]: 0, P["q_3"]: 0}, "Tangherlini's hole with no charge")
    elif system == "continued":
        source = _source("sphere")
        Q = source.reader.parameters
        rho = X[1]
        image = [X[0], sp.sqrt(rho ** 2 - P["q_3"])] + list(X[2:])
        names = {Q[n]: P[n] for n in ("mu", "q_1", "q_2", "q_3", "L")}
        outside = [e for e in _events(chart, system) if e[rho] ** 2 > e[P["q_3"]]]
        _pulled(chart, system, source, image, names, outside, "the sphere chart pulled back along r^2 = rho^2 - q_3")
        q, mu, L = sp.Symbol("q", positive=True), P["mu"], P["L"]
        equal = _explicit(chart).subs({P["q_1"]: q, P["q_2"]: q, P["q_3"]: q})
        V = 1 - (2 * q + mu) / rho ** 2 + q * (q + mu) / rho ** 4 + rho ** 2 / L ** 2
        want = [-V, 1 / V, rho ** 2, rho ** 2 * sp.sin(X[2]) ** 2,
                rho ** 2 * sp.sin(X[2]) ** 2 * sp.sin(X[3]) ** 2]
        for e in _events(chart, system):
            e = {**e, q: e[P["q_2"]]}
            for a in range(5):
                _zero(sp.N((equal[a, a] - want[a]).subs(e), DIGITS),
                      f"continued: the Reissner-Nordstrom-anti-de Sitter hole with the charges equal, slot {a}{a}")
    elif system == "brane":
        source = _source("flat")
        Q = source.reader.parameters
        u, rp, L = X[4], P["r_+"], P["L"]
        kappa = [P[f"kappa_{i}"] for i in (1, 2, 3)]
        image = [X[0], rp / sp.sqrt(u), X[1], X[2], X[3]]
        names = {Q["L"]: L, Q["mu"]: rp ** 4 * sp.Mul(*[1 + k for k in kappa]) / L ** 2,
                 **{Q[f"q_{i}"]: kappa[i - 1] * rp ** 2 for i in (1, 2, 3)}}
        _pulled(chart, system, source, image, names, list(_events(chart, system)),
                "the flat chart pulled back along r = r_+/sqrt(u)")
        g = _explicit(chart)
        lapse, stretch = -g[0, 0], g[4, 4]
        surface = sp.diff(lapse, u) / (2 * sp.sqrt(lapse * stretch))
        T0 = rp / (sp.pi * L ** 2)
        s1, s3 = sum(kappa), sp.Mul(*kappa)
        want = T0 * (2 + s1 - s3) / (2 * sp.sqrt(sp.Mul(*[1 + k for k in kappa])))
        for e in _events(chart, system):
            horizon = {**e, u: 1}
            f = _numeric(system, chart, P["f"], horizon)
            _zero(f, "brane: f at u = 1")
            # 2 pi T = |d sqrt(-g_tt)/d(proper distance)| on the horizon, in the chart x^0 = ct with c = 1.
            near = {**e, u: 1 - sp.Rational(1, 10 ** 12)}
            temperature = abs(sp.N(surface.subs(near), DIGITS)) / (2 * sp.pi)
            if abs(sp.N(temperature - want.subs(e), DIGITS)) > sp.Float("1e-9"):
                raise AssertionError(f"{METRIC}: brane: the Hawking temperature is {temperature}, not {want.subs(e)}")
