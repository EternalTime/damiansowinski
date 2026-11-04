"""The charged black string of three dimensions, Horne and Horowitz's exact solution, in six charts.

The static chart is Horne and Horowitz's eq. (1) (Nucl. Phys. B 368, 444, 1992),

    ds^2 = -(1 - M/r) dt^2 + (1 - Q^2/Mr) dx^2 + (1 - M/r)^-1 (1 - Q^2/Mr)^-1 k dr^2/8r^2,

with H_rtx = Q/r^2 and the dilaton Phi = ln r + (1/2) ln(k/2), their eq. (17), in units alpha' = 1/2.
The group chart is their eqs. (9) and (10), the coordinates u and v of SL(2,R) in which they gauge the
model, written in M and Q through lambda = Q^2/(M^2 - Q^2) and r = M - (M^2 - Q^2)uv/M. The ingoing
chart is Ho, Kim and Lee's eq. (10) (JHEP 08 (2022) 018). Their interior coordinate, eta = x - rbar of
their eqs. (11) and (14), is not printed: its metric carries sqrt(M/r - 1), which the printer reads in
sympy's sign as sqrt(-(r - M)) and cannot write. The extremal charts are Horne and Horowitz's eqs. (34) and (35) at |Q| = M, with
rho^2 = r - M outside the horizon and rho^2 = M - r inside it, and the superextremal chart their
eq. (38) at |Q| > M, with rho^2 = r - Q^2/M. three_dim_black_string.md beside this file records the
sources, the maps and why each chart is drawn or not.

Every check is made at points: each chart's metric against its source form, each chart against the
static one pulled back, and each chart against the three field equations of the page's convention,
from the chart's own Ricci tensor and Christoffel symbols:

    R_ab - nabla_a nabla_b Phi - (1/4) H_acd H_b^cd = 0,
    R - 2 nabla^2 Phi - (nabla Phi)^2 - H^2/12 + 8/k = 0,
    nabla_a (e^Phi H^abc) = 0,

with H = q eps, eps the volume form and q = Q sqrt(8/k)/r, which is H_rtx = Q/r^2 in the static chart.
"""
import sympy as sp

METRIC = "three_dim_black_string"
SYSTEMS = ("static", "wzw", "eddington_finkelstein_ingoing", "extremal", "extremal_interior", "superextremal")
NAMES = {"static": "Static", "wzw": "Wess-Zumino-Witten",
         "eddington_finkelstein_ingoing": "Ingoing Eddington-Finkelstein",
         "extremal": "Extremal", "extremal_interior": "Extremal Interior", "superextremal": "Superextremal"}
COORDS = {"static": ["t", "x", "r"], "wzw": ["u", "v", "x"], "eddington_finkelstein_ingoing": ["v", "x", "r"],
          "extremal": ["t", "x", "\\rho"],
          "extremal_interior": ["t", "x", "\\rho"], "superextremal": ["t", "x", "\\rho"]}
# The extremal charts are the string at |Q| = M and carry no Q of their own.
PARAMETERS = {system: (["M", "k"] if system.startswith("extremal") else ["M", "Q", "k"]) for system in SYSTEMS}

F = "\\left(1 - \\dfrac{M}{r}\\right)"
H = "\\left(1 - \\dfrac{Q^2}{Mr}\\right)"
LINES = {
    "static": f"ds^2 = -{F}dt^2 + {H}dx^2 + \\dfrac{{k\\,dr^2}}{{8r^2{F}{H}}}",
    "wzw": ("ds^2 = -\\dfrac{k\\left(Q^2v^2\\,du^2 + Q^2u^2\\,dv^2 + 2\\left(2M^2 - (2M^2 - Q^2)uv\\right)du\\,dv\\right)}"
            "{8(1 - uv)\\left(M^2 - (M^2 - Q^2)uv\\right)} + \\dfrac{(M^2 - Q^2)(1 - uv)}{M^2 - (M^2 - Q^2)uv}\\,dx^2"),
    "eddington_finkelstein_ingoing": (f"ds^2 = -{F}dv^2 + \\sqrt{{\\dfrac{{k}}{{2}}}}\\,\\dfrac{{dv\\,dr}}"
                                      f"{{r\\sqrt{{1 - \\dfrac{{Q^2}}{{Mr}}}}}} + {H}dx^2"),
    "extremal": "ds^2 = \\dfrac{\\rho^2}{\\rho^2 + M}\\left(-dt^2 + dx^2\\right) + \\dfrac{k\\,d\\rho^2}{2\\rho^2}",
    "extremal_interior": "ds^2 = \\dfrac{\\rho^2}{M - \\rho^2}\\left(dt^2 - dx^2\\right) + \\dfrac{k\\,d\\rho^2}{2\\rho^2}",
    "superextremal": ("ds^2 = -\\dfrac{Q^2 - M^2 + M\\rho^2}{Q^2 + M\\rho^2}dt^2 + \\dfrac{M\\rho^2}{Q^2 + M\\rho^2}dx^2"
                      " + \\dfrac{Mk\\,d\\rho^2}{2\\left(Q^2 - M^2 + M\\rho^2\\right)}"),
}
LINE = "(-\\infty, \\infty)"
DOMAINS = {
    "static": ["t \\in " + LINE, "x \\in " + LINE, "r \\in (0, \\infty)", "0 < |Q| < M,\\quad k > 0",
               "r = M \\;\\text{(the event horizon)}", "r = Q^2/M \\;\\text{(the inner horizon)}"],
    "wzw": ["u \\in " + LINE, "v \\in " + LINE, "x \\in " + LINE, "uv < 1", "0 < |Q| < M,\\quad k > 0",
            "uv = 0 \\;\\text{(the event horizon)}"],
    "eddington_finkelstein_ingoing": ["v \\in " + LINE, "x \\in " + LINE, "r \\in (Q^2/M, \\infty)",
                                      "0 < |Q| < M,\\quad k > 0", "r = M \\;\\text{(the event horizon)}"],
    "extremal": ["t \\in " + LINE, "x \\in " + LINE, "\\rho \\in (0, \\infty)", "M > 0,\\quad k > 0",
                 "\\rho = 0 \\;\\text{(the event horizon)}"],
    "extremal_interior": ["t \\in " + LINE, "x \\in " + LINE, "\\rho \\in (0, \\sqrt{M})", "M > 0,\\quad k > 0",
                          "\\rho = 0 \\;\\text{(the event horizon)}"],
    "superextremal": ["t \\in " + LINE, "x \\in [0, \\pi|Q|\\sqrt{2k/(Q^2 - M^2)})", "\\rho \\in [0, \\infty)",
                      "0 < M < |Q|,\\quad k > 0"],
}
TIME = {"static": "t", "eddington_finkelstein_ingoing": "v", "extremal": "t",
        "extremal_interior": "t", "superextremal": "t"}


def register_dimensions(vm):
    """With alpha' = 1/2 a number and c = 1, t, x, v and r are lengths, u and v of the group are
    numbers, rho is the square root of a length, and k, which multiplies the dimensionless dr^2/r^2,
    stands with that factor for a squared length, so it is declared L^2."""
    for system in SYSTEMS:
        dims = {"t": "L", "x": "L", "r": "L", "\\rho": "L**(1/2)", "M": "L", "Q": "L",
                "k": "L**2"}
        if system == "wzw":
            dims.update({"u": "1", "v": "1"})
        if system == "eddington_finkelstein_ingoing":
            dims["v"] = "L"
        vm.DIMENSIONS[(METRIC, system)] = {c: dims[c] for c in COORDS[system] + PARAMETERS[system]}


def spec(system):
    import verify_metrics as vm
    probe = vm.Reader(COORDS[system], PARAMETERS[system], ())
    p = probe.parameters
    lead = [p[name] for name in ("M", "Q", "k") if name in p]
    return {
        "metric_id": METRIC,
        "system": {"id": system, "name": NAMES[system], "coords": COORDS[system], "domains": DOMAINS[system],
                   "parameters": PARAMETERS[system], "line_element": LINES[system]},
        "chart_line_element": LINES[system],
        "printer": {"lead": lead, "factors": lead},
        "check": lambda chart: check(chart, system),
    }


def register(module):
    module.CHARTS[METRIC] = [lambda system=system: spec(system) for system in SYSTEMS]


# -- Checks at points, independent of the printed tensors ---------------------------------------

def radius(system, coords, M, Q):
    """The radius r of the static chart as a function of the chart's coordinates."""
    a, b, c = coords
    return {"static": c, "eddington_finkelstein_ingoing": c,
            "wzw": M - (M ** 2 - Q ** 2) * a * b / M, "extremal": c ** 2 + M, "extremal_interior": M - c ** 2,
            "superextremal": c ** 2 + Q ** 2 / M}[system]


def static_metric(M, Q, k, r):
    """Horne and Horowitz's eq. (1) at a point, in (t, x, r)."""
    f, h = 1 - M / r, 1 - Q ** 2 / (M * r)
    return sp.diag(-f, h, k / (8 * r ** 2 * f * h))


def source_metric(system, M, Q, k, a, b, c):
    """Each chart's source form at a point, written out by hand in its own coordinates (a, b, c)."""
    if system == "static":
        return static_metric(M, Q, k, c)
    if system == "wzw":
        u, v = a, b
        d = 8 * (1 - u * v) * (M ** 2 - (M ** 2 - Q ** 2) * u * v)
        cross = -k * (2 * M ** 2 - (2 * M ** 2 - Q ** 2) * u * v) / d
        return sp.Matrix([[-k * Q ** 2 * v ** 2 / d, cross, 0], [cross, -k * Q ** 2 * u ** 2 / d, 0],
                          [0, 0, (M ** 2 - Q ** 2) * (1 - u * v) / (M ** 2 - (M ** 2 - Q ** 2) * u * v)]])
    if system == "eddington_finkelstein_ingoing":
        r = c
        s = sp.sqrt(k / 8) / (r * sp.sqrt(1 - Q ** 2 / (M * r)))
        return sp.Matrix([[-(1 - M / r), 0, s], [0, 1 - Q ** 2 / (M * r), 0], [s, 0, 0]])
    rho = c
    if system == "extremal":
        w = rho ** 2 / (rho ** 2 + M)
        return sp.diag(-w, w, k / (2 * rho ** 2))
    if system == "extremal_interior":
        w = rho ** 2 / (M - rho ** 2)
        return sp.diag(w, -w, k / (2 * rho ** 2))
    n = Q ** 2 + M * rho ** 2
    return sp.diag(-(Q ** 2 - M ** 2 + M * rho ** 2) / n, M * rho ** 2 / n,
                   M * k / (2 * (Q ** 2 - M ** 2 + M * rho ** 2)))


def jacobian(system, M, Q, k, a, b, c):
    """d(t, x, r)/d(chart coordinates) at a point, from each chart's map onto the static one."""
    if system == "static":
        return sp.eye(3)
    if system == "wzw":
        u, v = a, b
        kappa = sp.sqrt(2 * (M ** 2 - Q ** 2) / k) / M
        # t = ln|u/v|/2 kappa, r = M - (M^2 - Q^2)uv/M, x = x.
        return sp.Matrix([[1 / (2 * kappa * u), -1 / (2 * kappa * v), 0], [0, 0, 1],
                          [-(M ** 2 - Q ** 2) * v / M, -(M ** 2 - Q ** 2) * u / M, 0]])
    if system == "eddington_finkelstein_ingoing":
        r = c
        f, h = 1 - M / r, 1 - Q ** 2 / (M * r)
        # dt = dv - sqrt(k/8) dr/(r f sqrt(h)).
        return sp.Matrix([[1, 0, -sp.sqrt(k / 8) / (r * f * sp.sqrt(h))], [0, 1, 0], [0, 0, 1]])
    sign = -1 if system == "extremal_interior" else 1
    return sp.Matrix([[1, 0, 0], [0, 1, 0], [0, 0, 2 * sign * c]])


# Points (M, Q, k) and chart coordinates (a, b, c), in every region each chart covers.
POINTS = {
    "static": [((2, 1, sp.Rational(52, 23)), (sp.Rational(3, 10), sp.Rational(1, 7), r))
               for r in (sp.Rational(1, 4), 1, 3)] + [((3, 2, 5), (1, 2, sp.Rational(7, 5)))],
    "wzw": [((2, 1, sp.Rational(52, 23)), (u, v, sp.Rational(1, 7)))
            for u, v in ((sp.Rational(1, 2), -1), (sp.Rational(1, 3), 1), (-2, sp.Rational(1, 3)),
                         (-sp.Rational(1, 2), -1))] + [((3, 2, 5), (sp.Rational(1, 4), sp.Rational(1, 2), 2))],
    "eddington_finkelstein_ingoing": [((2, 1, sp.Rational(52, 23)), (sp.Rational(3, 10), sp.Rational(1, 7), r))
                                      for r in (sp.Rational(3, 5), 1, 3)] + [((3, 2, 5), (1, 2, 2))],
    "extremal": [((2, 2, sp.Rational(52, 23)), (sp.Rational(3, 10), sp.Rational(1, 7), rho))
                 for rho in (sp.Rational(1, 2), 2)] + [((3, 3, 5), (1, 2, sp.Rational(1, 3)))],
    "extremal_interior": [((2, 2, sp.Rational(52, 23)), (sp.Rational(3, 10), sp.Rational(1, 7), rho))
                          for rho in (sp.Rational(1, 2), 1)] + [((3, 3, 5), (1, 2, sp.Rational(3, 2)))],
    "superextremal": [((1, 2, sp.Rational(52, 23)), (sp.Rational(3, 10), sp.Rational(1, 7), rho))
                      for rho in (sp.Rational(1, 10), 2)] + [((2, 3, 5), (1, 2, 1))],
}
TOLERANCE = 1e-25


def small(values, where, what):
    worst = max(abs(complex(sp.N(x, 40))) for x in values)
    if worst > TOLERANCE:
        raise AssertionError(f"{METRIC}/{where}: {what} ({worst:.3g})")


def levi_civita(a, b, c):
    return sp.Integer((a - b) * (b - c) * (c - a)) / 2


def check(chart, system):
    s = chart.symbols
    p = chart.reader.parameters
    M = p["M"]
    Q = p.get("Q", M)
    k = p["k"]
    r = radius(system, s, M, Q)
    phi = sp.log(r) + sp.log(k / 2) / 2
    dphi = [sp.diff(phi, x) for x in s]
    ddphi = [[sp.diff(phi, x, y) for y in s] for x in s]
    gamma = chart.geo.christoffel_ull()
    ricci = chart.geo.ricci_ll()
    scalar = chart.geo.ricci_scalar()
    for values, point in POINTS[system]:
        (Mv, Qv, kv), (av, bv, cv) = map(sp.sympify, values), map(sp.sympify, point)
        at = {s[0]: av, s[1]: bv, s[2]: cv, M: Mv, k: kv}
        if "Q" in p:
            at[Q] = Qv
        where = f"{system} at {(av, bv, cv)}, (M, Q, k) = {(Mv, Qv, kv)}"
        g = chart.geo.g.subs(at)
        small(g - source_metric(system, Mv, Qv, kv, av, bv, cv), where, "the metric differs from its source")
        rv = radius(system, (av, bv, cv), Mv, Qv)
        J = jacobian(system, Mv, Qv, kv, av, bv, cv)
        small(g - J.T * static_metric(Mv, Qv, kv, rv) * J, where, "not the static chart pulled back")
        gi = g.inv()
        root = sp.sqrt(-g.det())
        q = Qv * sp.sqrt(8 / sp.nsimplify(kv)) / rv
        # H_abc = q sqrt(-g) [abc]; in the static chart that is H_rtx = Q/r^2, the source's field.
        Hl = [[[q * root * levi_civita(a, b, c) for c in range(3)] for b in range(3)] for a in range(3)]
        if system == "static":
            small([Hl[2][0][1] - Qv / rv ** 2], where, "H_rtx is not Q/r^2")
        hh = sp.Matrix(3, 3, lambda m, n: sum(Hl[m][a][b] * Hl[n][c][d] * gi[a, c] * gi[b, d]
                                                for a in range(3) for b in range(3) for c in range(3)
                                                for d in range(3)))
        H2 = sum(Hl[a][b][c] * Hl[d][e][f] * gi[a, d] * gi[b, e] * gi[c, f] for a in range(3) for b in range(3)
                 for c in range(3) for d in range(3) for e in range(3) for f in range(3))
        hess = sp.Matrix(3, 3, lambda m, n: (ddphi[m][n] - sum(gamma[c][m][n] * dphi[c] for c in range(3))).subs(at))
        grad = sp.Matrix([d.subs(at) for d in dphi])
        R = sp.Matrix(3, 3, lambda m, n: ricci[m][n].subs(at))
        small(R - hess - hh / 4, where, "R_ab - nabla_a nabla_b Phi - HH/4 does not vanish")
        # The same equation for phi = -Phi/2, R_ab + 2 nabla_a nabla_b phi - HH/4.
        small(R + 2 * (-hess / 2) - hh / 4, where, "R_ab + 2 nabla_a nabla_b phi - HH/4 does not vanish")
        box = sum(gi[m, n] * hess[m, n] for m in range(3) for n in range(3))
        square = (grad.T * gi * grad)[0, 0]
        small([scalar.subs(at) - 2 * box - square - H2 / 12 + 8 / sp.nsimplify(kv)], where,
              "the dilaton's equation does not hold with 8/k")
        # nabla_a(e^Phi H^abc) = (1/sqrt(-g)) d_a(sqrt(-g) e^Phi H^abc), and sqrt(-g) e^Phi H^abc is
        # e^Phi q times a constant array, so the equation holds where e^Phi q is constant.
        carried = sp.exp(phi) * Q * sp.sqrt(8 / k) / r
        small([sp.diff(carried, x).subs(at) for x in s], where, "e^Phi H is not divergence free")
