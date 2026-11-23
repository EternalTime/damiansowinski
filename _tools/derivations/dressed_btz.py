"""The charts of Martinez and Zanelli's conformally dressed black hole of three dimensions,
gr-qc/9604021, Phys. Rev. D 54, 3830 (1996). dressed_btz.md beside this file records each chart's
source.

Their eq. (ds) with N = 1 and their (Fp) and (psie) give

    ds^2 = -F c^2 dt^2 + dr^2/F + r^2 dtheta^2,   F = (r + B)^2 (r - 2B)/(r l^2),
    Psi = sqrt(8B/(kappa (r + B))),

a solution of G - g/l^2 = kappa T with the stress tensor of a massless scalar coupled through
-R Psi^2/16, their (Eeq), (Feq) and (Tuv). Their footnote [footnote] gives the advanced time chart,
ds^2 = -F c^2 dv^2 + 2c dv dr + r^2 dtheta^2, with c v = c t + int dr/F. The cosmological constant
is -1/l^2 and the horizon is r = 2B.
"""
import json

import sympy as sp

METRIC = "dressed_btz"
CHARTS = ("static", "eddington_finkelstein_ingoing")

DIMENSIONS = {
    (METRIC, "static"): {"t": "T", "r": "L", "\\theta": "1", "l": "L", "B": "L"},
    (METRIC, "eddington_finkelstein_ingoing"): {"v": "T", "r": "L", "\\theta": "1", "l": "L", "B": "L"},
}

F = "\\dfrac{(r + B)^2(r - 2B)}{rl^2}"
NAMES = {"static": "Static", "eddington_finkelstein_ingoing": "Ingoing Eddington-Finkelstein"}
COORDS = {"static": ["t", "r", "\\theta"], "eddington_finkelstein_ingoing": ["v", "r", "\\theta"]}
REST = ["r \\in (0, \\infty)", "\\theta \\in [0, 2\\pi)", "r = 2B \\;\\text{(the horizon)}"]
DOMAINS = {"static": ["t \\in (-\\infty, \\infty)"] + REST,
           "eddington_finkelstein_ingoing": ["v \\in (-\\infty, \\infty)"] + REST}
KRETSCHMANN = "\\dfrac{12\\left(r^6 + 2B^6\\right)}{l^4r^6}"


def line(system, c2, c1):
    """The line element of a chart, c2 standing before the square of the time and c1 before its
    product with dr."""
    if system == "static":
        return ("ds^2 = -" + F + c2 + "dt^2 + \\dfrac{rl^2\\,dr^2}{(r + B)^2(r - 2B)} + r^2d\\theta^2")
    return "ds^2 = -" + F + c2 + "dv^2 + 2" + c1 + "\\,dv\\,dr + r^2d\\theta^2"


def chart_spec(system):
    """One chart: its line element with c, the same in x^0 = ct or cv, its metric and inverse
    written as the line element has them, the Kretschmann scalar of Martinez and Zanelli's (KR),
    and the check against their field equations and the live BTZ charts at B = 0."""
    import verify_metrics as vm
    coords = COORDS[system]
    t = coords[0]
    if system == "static":
        metric = {(t, t): "-" + F, ("r", "r"): "\\dfrac{rl^2}{(r + B)^2(r - 2B)}", ("\\theta", "\\theta"): "r^2"}
        inverse = {(t, t): "-\\dfrac{rl^2}{(r + B)^2(r - 2B)}", ("r", "r"): F,
                   ("\\theta", "\\theta"): "\\dfrac{1}{r^2}"}
    else:
        metric = {(t, t): "-" + F, (t, "r"): "1", ("r", t): "1", ("\\theta", "\\theta"): "r^2"}
        inverse = {(t, "r"): "1", ("r", t): "1", ("r", "r"): F, ("\\theta", "\\theta"): "\\dfrac{1}{r^2}"}
    probe = vm.Reader(coords, ["l", "B"], ())
    r, l, B = probe.symbol["r"], probe.parameters["l"], probe.parameters["B"]
    return {
        "metric_id": METRIC,
        "system": {"id": system, "name": NAMES[system], "coords": coords, "domains": DOMAINS[system],
                   "parameters": ["l", "B"], "line_element": line(system, "c^2", "c")},
        "chart_line_element": line(system, "", ""),
        "printer": {"lead": [r, B, l], "factors": [B, l, r], "flip": False},
        "components": {"metric_components": metric, "inverse_metric_components": inverse},
        "kretschmann": KRETSCHMANN,
        "check": lambda chart: check_chart(chart, system),
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


# Points of (l, B, r) at which the field equations are checked in forty digits: outside the horizon,
# inside it, and far from it.
POINTS = [(sp.Rational(7, 5), sp.Rational(3, 10), sp.Rational(17, 10)),
          (sp.Integer(1), sp.Rational(1, 2), sp.Rational(3, 10)),
          (sp.Rational(3, 2), sp.Rational(2, 7), sp.Integer(5))]


def check_chart(chart, system):
    """Each chart at points against Martinez and Zanelli's field equations and statements.

    With chi^2 = kappa Psi^2 = 8B/(r + B), their (psie), kappa T_ab is their (Tuv),
    d_a chi d_b chi - g_ab (d chi)^2/2 + (g_ab Box - nabla_a nabla_b + G_ab) chi^2/8, and
    G_ab - g_ab/l^2 - kappa T_ab vanishes, their (Eeq); chi obeys Box chi = R chi/8, their (Feq);
    R = -6/l^2, their (R); the Kretschmann scalar is their (KR); F is their (Fp) in its two forms
    and Martinez's gr-qc/9801091 (ds) with r_+ = 2B. At B = 0 each chart is the BTZ chart of the
    same name at M = J = 0 with phi for theta and ell for l, and the ingoing chart is the static
    one along c v = c t + int dr/F."""
    import verify_metrics as vm
    geo, g, ginv = chart.geo, chart.geo.g, chart.geo.ginv
    X = list(chart.symbols)
    n = len(X)
    P = chart.reader.parameters
    l, B = P["l"], P["B"]
    r = X[1]
    name = f"{METRIC}/{system}"
    Fx = (r + B) ** 2 * (r - 2 * B) / (r * l ** 2)
    chi2 = 8 * B / (r + B)
    chi = sp.sqrt(chi2)
    gamma = geo.christoffel_ull()

    def hessian(f):
        return sp.Matrix(n, n, lambda a, b: sp.diff(f, X[a], X[b])
                         - sum(gamma[c][a][b] * sp.diff(f, X[c]) for c in range(n)))

    def box(f):
        h = hessian(f)
        return sum(ginv[a, b] * h[a, b] for a in range(n) for b in range(n))

    def vanishes(e, what):
        for ll, bb, rr in POINTS:
            value = sp.N(sp.sympify(e).subs({l: ll, B: bb, r: rr}), 40)
            if abs(complex(value)) > 1e-30:
                raise AssertionError(f"{name}: {what} fails at l = {ll}, B = {bb}, r = {rr}")

    einstein = sp.Matrix(n, n, lambda a, b: vm._at(geo.einstein_ll(), (a, b)))
    d = [sp.diff(chi, x) for x in X]
    grad2 = sum(ginv[a, b] * d[a] * d[b] for a in range(n) for b in range(n))
    h2, b2 = hessian(chi2), box(chi2)
    for a in range(n):
        for b in range(n):
            kT = d[a] * d[b] - g[a, b] * grad2 / 2 + (g[a, b] * b2 - h2[a, b] + einstein[a, b] * chi2) / 8
            vanishes(einstein[a, b] - g[a, b] / l ** 2 - kT,
                     f"Martinez and Zanelli's (Eeq) in slot {chart.coords_tex[a]}{chart.coords_tex[b]}")
    ricci = geo.ricci_scalar()
    vanishes(box(chi) - ricci * chi / 8, "the field equation (Feq)")
    vanishes(ricci + 6 / l ** 2, "R = -6/l^2, their (R)")
    vanishes(geo.kretschmann() - 12 * (r ** 6 + 2 * B ** 6) / (l ** 4 * r ** 6), "the Kretschmann scalar (KR)")
    vanishes(Fx - (r ** 2 - 3 * B ** 2 - 2 * B ** 3 / r) / l ** 2, "the two forms of (Fp)")
    vanishes(Fx - (2 * r + 2 * B) ** 2 * (r - 2 * B) / (4 * r * l ** 2), "Martinez's (ds) with r_+ = 2B")
    vanishes(g[0, 0] + Fx, "g_00 = -F")

    # At B = 0 the chart is the live BTZ chart without spin and of zero mass.
    btz = "stationary" if system == "static" else "eddington_finkelstein_ingoing"
    reader, hole, coords = published_metric(vm, "btz", btz)
    swap = dict(zip((reader.symbol[x] for x in coords), X))
    swap.update({reader.parameters["ell"]: l, reader.parameters["M"]: 0, reader.parameters["J"]: 0})
    empty = hole.subs(swap)
    for a in range(n):
        for b in range(n):
            if vm.norm(empty[a, b] - g[a, b].subs(B, 0)) != 0:
                raise AssertionError(f"{name}: at B = 0 the chart misses the published btz/{btz} "
                                     f"in slot {chart.coords_tex[a]}{chart.coords_tex[b]}")

    if system == "eddington_finkelstein_ingoing":
        # The static chart along c t = c v - int dr/F.
        J = sp.eye(n)
        J[0, 1] = -1 / Fx
        static = sp.diag(-Fx, 1 / Fx, r ** 2)
        pulled = J.T * static * J
        for a in range(n):
            for b in range(n):
                if vm.norm(pulled[a, b] - g[a, b]) != 0:
                    raise AssertionError(f"{name}: the static chart along v misses slot "
                                         f"{chart.coords_tex[a]}{chart.coords_tex[b]}")
