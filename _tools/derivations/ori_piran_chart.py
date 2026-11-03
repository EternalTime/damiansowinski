"""The two charts of self similar collapse of a perfect fluid with P = k rho.

Harada and Maeda, arXiv:gr-qc/0101064, section II: the comoving line element (1) with every
dimensionless function of z = r/(-t) alone, the closed forms of e^sigma and e^omega after their
equations (domegadt) and (definitionD), the ordinary differential equations (dmdz), (dsdz) and
(constraint), and the zooming coordinates tau = -ln(-t), z = r/(-t), here with a length ell.
The checks compare at numerical points: the comoving Einstein tensor with the fluid's stresses on
local solutions of the ordinary differential equations, and the zooming metric with the comoving
one pulled back.
"""
import sympy as sp

import verify_metrics as vm

METRIC = "ori_piran"
SYSTEMS = ("comoving", "zooming")
SPHERE = "\\left(d\\theta^2 + \\sin^2\\theta\\,d\\phi^2\\right)"
ANGLES = ["\\theta \\in [0, \\pi]", "\\phi \\in [0, 2\\pi)"]
COMOVING = ["\\sigma = \\sigma(t, r)", "\\omega = \\omega(t, r)", "S = S(t, r)", "k", "D"]
ZOOMING = ["\\sigma = \\sigma(z)", "\\omega = \\omega(z)", "S = S(z)", "\\ell"]
COMOVING_LINE = "ds^2 = -e^{\\sigma}dt^2 + e^{\\omega}dr^2 + r^2S^2" + SPHERE
ZOOMING_LINE = ("ds^2 = \\ell^2e^{-2\\tau}\\left(\\left(z^2\\,e^{\\omega} - e^{\\sigma}\\right)d\\tau^2"
                " - 2z\\,e^{\\omega}d\\tau\\,dz + e^{\\omega}dz^2 + z^2S^2" + SPHERE + "\\right)")


def chart_spec(system):
    if system == "comoving":
        coords = ["t", "r", "\\theta", "\\phi"]
        probe = vm.Reader(coords, COMOVING, ())
        t, r = probe.symbol["t"], probe.symbol["r"]
        functions = [probe.parameters[n] for n in ("sigma", "omega", "S")]
        lead = ([sp.Derivative(f, *v) for f in functions
                 for v in ((t, 2), (r, 2), (t, r))]
                + [sp.Derivative(f, v) for f in functions for v in (t, r)] + functions + [r])
        return {"metric_id": METRIC,
                "system": {"id": system, "name": "Comoving", "coords": coords,
                           "domains": ["t \\in (-\\infty, 0)", "r \\in (0, \\infty)"] + ANGLES,
                           "parameters": COMOVING, "line_element": COMOVING_LINE},
                "chart_line_element": COMOVING_LINE,
                "printer": {"lead": lead, "factors": functions + [r] + lead[:-len(functions) - 1]},
                "check": check_comoving}
    coords = ["\\tau", "z", "\\theta", "\\phi"]
    probe = vm.Reader(coords, ZOOMING, ())
    z = probe.symbol["z"]
    functions = [probe.parameters[n] for n in ("sigma", "omega", "S")]
    lead = ([sp.Derivative(f, (z, 2)) for f in functions] + [sp.Derivative(f, z) for f in functions]
            + functions + [z, probe.parameters["ell"]])
    return {"metric_id": METRIC,
            "system": {"id": system, "name": "Zooming", "coords": coords,
                       "domains": ["\\tau \\in (-\\infty, \\infty)", "z \\in (0, \\infty)"] + ANGLES,
                       "parameters": ZOOMING, "line_element": ZOOMING_LINE},
            "chart_line_element": ZOOMING_LINE,
            "printer": {"lead": lead, "factors": [probe.parameters["ell"], z] + functions + lead[:6]},
            "check": check_zooming}


def equations(k):
    """The right sides of Harada and Maeda's (dmdz), (dsdz) and the dependent equation for eta,
    in x = ln z, for ln M, ln S and ln eta, in the gauge a_sigma = (2D)^(2k/(1+k)), as sympy
    expressions of those logarithms, x and D."""
    lM, lS, le, x, D = sp.symbols("lM lS le x D", real=True)
    M, S, eta, z = sp.exp(lM), sp.exp(lS), sp.exp(le), sp.exp(x)
    a = (2 * D) ** (2 * k / (1 + k))
    y = M / (eta * S ** 3)
    V2 = z ** (2 * (1 - k) / (1 + k)) * eta ** (-2 * (1 - k) / (1 + k)) * S ** -4 / a
    num = 4 * y * V2 - (1 + k) ** 2 * S ** -4 * eta ** (-(1 - k) / (1 + k))
    rates = [k / (1 + k) * (1 / y - 1), -(1 - y) / (1 + k), 2 - num / (2 * (V2 - k))]
    constraint = (1 - y) ** 2 * V2 - (k + y) ** 2 + (1 + k) ** 2 * eta ** (-2 / (1 + k)) * S ** -6 * (1 - M / S)
    return (lM, lS, le, x, D), rates, constraint


def local_states(k):
    """Two states on the constraint surface, away from the sonic point, each with ln S and ln eta
    to second order in x = ln z along the flow of the equations."""
    (lM, lS, le, x, D), rates, constraint = equations(k)
    states = []
    # Near the relativistic Larson-Penston solution at k = 1/100: z = 1/200, below the sonic
    # point, and z = 1/10, above it, with the mass function solved from the constraint.
    for z0, S0, eta0, M0 in ((sp.Rational(1, 200), 36, sp.Rational(5, 100000), sp.Rational(89, 100)),
                             (sp.Rational(1, 10), 8, sp.Rational(28, 10000), sp.Rational(92, 100))):
        xv, sv, ev, dv = sp.log(z0), sp.log(S0), sp.log(eta0), sp.Rational(1438, 1000)
        at = {x: xv, lS: sv, le: ev, D: dv}
        lMv = sp.nsolve(constraint.subs(at), lM, sp.N(sp.log(M0), 40), prec=40)
        at[lM] = lMv
        first = [sp.N(f.subs(at), 40) for f in rates]

        def along(f):
            return (sp.diff(f, x) + sum(sp.diff(f, v) * rate for v, rate in zip((lM, lS, le), rates)))
        second = [sp.N(along(f).subs(at), 40) for f in rates[1:]]
        states.append({"x": xv, "D": dv, "lS": (sv, first[1], second[0]), "le": (ev, first[2], second[1])})
    return states


def check_comoving(chart):
    """The Einstein tensor equals 8 pi times the fluid's stresses, rho = eta/(8 pi r^2) and
    P = k rho at rest in the chart, on two local solutions of Harada and Maeda's equations, with
    e^sigma and e^omega their closed forms; and the mixed G^t_r vanishes there."""
    P = chart.reader.parameters
    t, r = chart.symbols[0], chart.symbols[1]
    sigma, omega, S = P["sigma"], P["omega"], P["S"]
    k = sp.Rational(1, 100)
    for state in local_states(k):
        X = sp.log(r / (-t)) - state["x"]
        lS = state["lS"][0] + state["lS"][1] * X + state["lS"][2] * X ** 2 / 2
        le = state["le"][0] + state["le"][1] * X + state["le"][2] * X ** 2 / 2
        eta = sp.exp(le)
        model = {sigma: 2 * k / (1 + k) * (sp.log(2 * state["D"]) + 2 * sp.log(r / (-t)) - le),
                 omega: -2 / (1 + k) * le - 4 * lS, S: sp.exp(lS)}
        point = {t: -1, r: sp.exp(state["x"]), chart.symbols[2]: sp.Rational(9, 10)}
        g = chart.geo.g.subs(model).doit()
        G = chart.geo.einstein_ll()
        rho = eta / (8 * sp.pi * r ** 2)
        u = sp.Matrix([-sp.sqrt(-g[0, 0]), 0, 0, 0])
        T = (1 + k) * rho * u * u.T + k * rho * g
        for i, j in ((0, 0), (1, 1), (2, 2), (3, 3), (0, 1)):
            got = sp.N(sp.sympify(G[i][j]).subs(model).doit().subs(point), 30)
            want = sp.N((8 * sp.pi * T[i, j]).subs(point), 30)
            if abs(got - want) > 1e-9 * (1 + abs(want)):
                raise AssertionError(f"{METRIC}/comoving: Einstein equation {i}{j}, {got} against {want}")


def check_zooming(chart):
    """The zooming metric is the comoving one pulled back along t = -ell e^(-tau), r = ell z e^(-tau),
    with each function the same function of z, at numerical points."""
    P = chart.reader.parameters
    tau, z, th = chart.symbols[0], chart.symbols[1], chart.symbols[2]
    ell = P["ell"]
    s, w, a = sp.symbols("ZPsig ZPom ZPS", real=True)
    names = {P["sigma"]: s, P["omega"]: w, P["S"]: a}
    t = -ell * sp.exp(-tau)
    r = ell * z * sp.exp(-tau)
    comoving = sp.diag(-sp.exp(s), sp.exp(w), (r * a) ** 2, (r * a * sp.sin(th)) ** 2)
    J = sp.Matrix([[sp.diff(t, tau), sp.diff(t, z), 0, 0], [sp.diff(r, tau), sp.diff(r, z), 0, 0],
                   [0, 0, 1, 0], [0, 0, 0, 1]])
    pulled = J.T * comoving * J
    mine = chart.geo.g.subs(names)
    for point in ({tau: sp.Rational(1, 3), z: sp.Rational(7, 4), th: sp.Rational(4, 5), ell: 2,
                   s: sp.Rational(1, 9), w: sp.Rational(5, 2), a: sp.Rational(3, 2)},
                  {tau: -1, z: sp.Rational(1, 20), th: sp.Rational(3, 2), ell: sp.Rational(1, 3),
                   s: sp.Rational(-1, 7), w: 3, a: 11}):
        for i in range(4):
            for j in range(4):
                if abs(sp.N((mine[i, j] - pulled[i, j]).subs(point), 30)) > 1e-25:
                    raise AssertionError(f"{METRIC}/zooming: slot {i}{j} is not the comoving chart pulled back")
