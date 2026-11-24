"""The flows of the renormalisation group of holography in closed form: walls of one scalar in the
gauged supergravity of five dimensions, ds^2 = e^{2A}(-c^2 dt^2 + dx_1^2 + dx_2^2 + dx_3^2) + dr^2,
with L the radius of the anti-de Sitter space each approaches at its boundary.

proper_distance  Girardello, Petrini, Porrati, and Zaffaroni, Nucl. Phys. B 569, 451 (2000),
                 hep-th/9909047, eqs. (m37) and (m38): phi(y) = (y + log[2 sinh(y - C_1)])/2 and
                 m(y) = (sqrt3/2) log[(1 + e^{-(y - C_1)})/(1 - e^{-(y - C_1)})] for the first order
                 equations of W = -(3/4)(cosh(2m/sqrt3) + 1), with e^{2 phi} the warp factor in their
                 ansatz ds^2 = dy^2 + e^{2 phi(y)} dx^mu dx_mu. With C_1 = 0 and the radius restored,
                 y = r/L, e^{2A} = e^{2r/L} - 1, and r is the proper distance from the singularity at r = 0.
bfs_u            Bianchi, Freedman, and Skenderis, JHEP 08 (2001) 041, hep-th/0105276,
                 Sec. 3.2, the same flow in u: e^{2A} = u/(1 - u), du/dr = 2(1 - u)/L, so u = 1 - e^{-2r/L},
                 and the scalar (sqrt3/2) log[(1 + sqrt(1 - u))/(1 - sqrt(1 - u))].
coulomb_disc     Bianchi, Freedman, and Skenderis, Sec. 3.1, the flow of the Coulomb branch whose lift is
                 a uniform disc of D3-branes: v = e^{sqrt6 Phi}, e^{2A} = (ell^2/L^2) v^{2/3}/(1 - v),
                 dv/dr = (2/L) v^{2/3}(1 - v), for W = -e^{-2 Phi/sqrt6} - e^{4 Phi/sqrt6}/2; the same
                 flow is Freedman, Gubser, Pilch, and Warner's n = 2, JHEP 07 (2000) 038, hep-th/9906194,
                 eq. (SUSYCases), with v = e^{-sqrt6 mu} and ell the radius of the disc.

Each chart is checked at points, in forty digits, independently of the printed tensors: its warp
factor and its g_rr against its source's equations, its first order flow equations, Einstein's
equations G_ab = 2 T_ab of BFS's action R/4 - (grad Phi)^2/2 - V with V = (dW/dPhi)^2/2 - 4W^2/3,
their eq. (potential), divided by L^2, the scalar's wave equation, Freedman, Gubser, Pilch, and
Warner's A'' = -(2/3) Phi'^2 of hep-th/9904017 in these units, and the Ricci scalar -20/L^2 of
anti-de Sitter space at the boundary; bfs_u is the proper distance chart pulled back as well.
holographic_rg_flow.md beside this file records each chart's source and every drawing decision.
"""
import sympy as sp

METRIC = "holographic_rg_flow"
CHARTS = ("proper_distance", "bfs_u", "coulomb_disc")
NAMES = {"proper_distance": "Confining Flow, Proper Distance", "bfs_u": "Confining Flow, Coordinate $u$",
         "coulomb_disc": "Coulomb Branch Disc"}
ALONG = ["x_1", "x_2", "x_3"]
FIFTH = {"proper_distance": "r", "bfs_u": "u", "coulomb_disc": "v"}
COORDS = {system: ["t"] + ALONG + [FIFTH[system]] for system in CHARTS}
PARAMETERS = {"proper_distance": ["L"], "bfs_u": ["L"], "coulomb_disc": ["L", "\\ell"]}
FLAT = "\\left(-%sdt^2 + dx_1^2 + dx_2^2 + dx_3^2\\right)"
LINES = {
    "proper_distance": "ds^2 = \\left(e^{2r/L} - 1\\right)" + FLAT + " + dr^2",
    "bfs_u": "ds^2 = \\dfrac{u}{1 - u}" + FLAT + " + \\dfrac{L^2\\,du^2}{4\\left(1 - u\\right)^2}",
    "coulomb_disc": ("ds^2 = \\dfrac{\\ell^2v^{2/3}}{L^2\\left(1 - v\\right)}" + FLAT
                     + " + \\dfrac{L^2\\,dv^2}{4v^{4/3}\\left(1 - v\\right)^2}"),
}
REALS = " \\in (-\\infty, \\infty)"
DOMAINS = {
    "proper_distance": ["r \\in (0, \\infty)", "r = 0 \\;\\text{(the singularity)}",
                        "r \\to \\infty \\;\\text{(the boundary)}"],
    "bfs_u": ["u \\in (0, 1)", "u = 0 \\;\\text{(the singularity)}", "u \\to 1 \\;\\text{(the boundary)}"],
    "coulomb_disc": ["v \\in (0, 1)", "v = 0 \\;\\text{(the singularity)}", "v \\to 1 \\;\\text{(the boundary)}"],
}
_DIMS = {"t": "T", "x_1": "L", "x_2": "L", "x_3": "L", "r": "L", "u": "1", "v": "1", "L": "L", "\\ell": "L"}
DIMENSIONS = {(METRIC, s): {name: _DIMS[name] for name in COORDS[s] + PARAMETERS[s]} for s in CHARTS}
HELD = {}
RATES = {}
DIGITS = 40
TOLERANCE = sp.Float("1e-28")


def register(module):
    module.CHARTS[METRIC] = [lambda system=system: chart_spec(system) for system in CHARTS]


def line(system, c2):
    return LINES[system] % c2


def chart_spec(system):
    """One of the three charts, its values printed in its fifth coordinate and its parameters."""
    import verify_metrics as vm

    probe = vm.Reader(COORDS[system], PARAMETERS[system], ())
    fifth, P = probe.symbol[FIFTH[system]], probe.parameters
    lead = [fifth] + [P[p.lstrip("\\")] for p in PARAMETERS[system]]
    spec = {"metric_id": METRIC,
            "system": {"id": system, "name": NAMES[system], "coords": COORDS[system],
                       "domains": ["t" + REALS] + [x + REALS for x in ALONG] + DOMAINS[system],
                       "parameters": PARAMETERS[system], "line_element": line(system, "c^2")},
            "chart_line_element": line(system, ""),
            "printer": {"lead": lead, "flip": False},
            "check": lambda chart: chart_check(chart, system)}
    if system == "proper_distance":
        root = sp.exp(fifth / P["L"])
        spec["pretty"] = warp_whole(root)
        spec["printer"]["rank"] = falling_exponentials(root)
    else:
        # 1 - u and 1 - v are printed as the line elements write them, and the length of the disc first.
        name = FIFTH[system]
        placeholder = sp.Symbol("one_minus_" + name)
        spec["pretty"] = named_sum(placeholder, 1 - fifth)
        spec["printer"]["named"] = {placeholder: "1 - " + name}
        spec["printer"]["factors"] = ([P["ell"]] if system == "coulomb_disc" else []) + [P["L"], fifth]
    return spec


def named_sum(placeholder, poly):
    """A `pretty` that factors a value and writes each factor that is the sum `poly`, or its
    negation, as the placeholder the printer names."""
    def pretty(value):
        out = sp.Integer(1)
        for f in sp.Mul.make_args(sp.factor(value)):
            base, k = (f.base, f.exp) if f.is_Pow else (f, sp.Integer(1))
            if base.is_Add and sp.expand(base - poly) == 0:
                base = placeholder
            elif base.is_Add and sp.expand(base + poly) == 0:
                base, out = placeholder, out * (-1) ** k
            out *= base ** k
        return out
    return pretty


def _in_root(value, root, X):
    """A value with every exponential written as a power of the symbol X standing for root."""
    unit = root.args[0]
    return sp.sympify(value).replace(lambda e: isinstance(e, sp.exp) and sp.expand(e.args[0] / unit).is_Integer,
                                     lambda e: X ** sp.expand(e.args[0] / unit))


def warp_whole(root):
    """A `pretty` that keeps the warp factor e^{2r/L} - 1 whole, as the line element writes it: the
    value is factored in X = e^{r/L}, and (X - 1)(X + 1) to one power is written back as X^2 - 1."""
    X = sp.Dummy("X", positive=True)

    def pretty(value):
        factored = sp.factor(_in_root(value, root, X))
        powers = factored.as_powers_dict()
        k = powers.get(X - 1, 0)
        if k and powers.get(X + 1, 0) == k:
            del powers[X - 1], powers[X + 1]
            powers[X ** 2 - 1] = k
        back = {X: root}
        return sp.Mul(*[sp.expand(b.subs(back)) ** e if b != X else root ** e for b, e in powers.items()],
                      evaluate=False)
    return pretty


def falling_exponentials(root):
    """A `rank` that orders the terms of a sum by their power of e^{r/L}, highest first."""
    X = sp.Dummy("X", positive=True)

    def rank(term):
        degree = sp.Poly(_in_root(term, root, X), X).degree() if _in_root(term, root, X).has(X) else 0
        return -degree
    return rank


# -- the checks ------------------------------------------------------------------------------------

def _numeric(expr, at):
    return sp.N(sp.sympify(expr).subs(at), DIGITS)


def _zero(value, what):
    if abs(complex(value)) > TOLERANCE:
        raise AssertionError(f"{METRIC}: {what} fails, by {value}")


def _gppz_superpotential(phi):
    """GPPZ99 (36), W_2, and BFS's (GPPZW): W = -(3/4)(cosh(2 phi/sqrt3) + 1)."""
    return -sp.Rational(3, 4) * (sp.cosh(2 * phi / sp.sqrt(3)) + 1)


def _coulomb_superpotential(phi):
    """BFS Sec. 3.1: W = -e^{-2 Phi/sqrt6} - e^{4 Phi/sqrt6}/2."""
    return -sp.exp(-2 * phi / sp.sqrt(6)) - sp.exp(4 * phi / sp.sqrt(6)) / 2


def _field_equations(chart, at, scalar, superpotential, L, what):
    """Einstein's equations G_ab = 2 T_ab, T_ab = d_a Phi d_b Phi - g_ab((grad Phi)^2/2 + V/L^2), and
    the wave equation box Phi = V'(Phi)/L^2, with V = (dW/dPhi)^2/2 - 4W^2/3, BFS (potential), for the
    Lorentzian form of their action (scaction), R/4 - (grad Phi)^2/2 - V/L^2, at the point."""
    import verify_metrics as vm

    x = chart.symbols
    n = len(x)
    g = sp.Matrix(chart.geo.g)
    ginv = sp.Matrix(chart.geo.ginv)
    f = sp.Symbol("Phi_field")
    W = superpotential(f)
    V = (sp.diff(W, f) ** 2 / 2 - sp.Rational(4, 3) * W ** 2) / L ** 2
    dV = sp.diff(V, f)
    dphi = [sp.diff(scalar, c) for c in x]
    dphi2 = sum(ginv[a, b] * dphi[a] * dphi[b] for a in range(n) for b in range(n))
    einstein = chart.geo.einstein_ll()
    for a in range(n):
        for b in range(a, n):
            T = dphi[a] * dphi[b] - g[a, b] * (dphi2 / 2 + V.subs(f, scalar))
            _zero(_numeric(vm._at(einstein, (a, b)) - 2 * T, at),
                  f"{what}: Einstein's equation in slot {chart.coords_tex[a]}{chart.coords_tex[b]}")
    root = sp.sqrt(-g.det())
    box = sum(sp.diff(root * ginv[a, b] * dphi[b], x[a]) for a in range(n) for b in range(n)) / root
    _zero(_numeric(box - dV.subs(f, scalar), at), f"{what}: the wave equation of the scalar")


def _flow(chart, at, scalar, superpotential, L, proper, what):
    """The first order equations dA/dr = -(2/3L) W and dPhi/dr = (1/L) dW/dPhi, BFS (1storder), with
    e^{2A} = g_{x_1 x_1} and d/dr = (g_55)^{-1/2} d/dx^5, and A'' = -(2/3) Phi'^2, FGPW99's inequality
    A'' <= 0 with its matter, at the point."""
    fifth = chart.symbols[-1]
    f = sp.Symbol("Phi_field")
    W = superpotential(f)
    A = sp.log(chart.geo.g[1, 1]) / 2
    d = (lambda e: sp.diff(e, fifth) / sp.sqrt(chart.geo.g[4, 4])) if not proper else (lambda e: sp.diff(e, fifth))
    _zero(_numeric(d(A) + 2 * W.subs(f, scalar) / (3 * L), at), f"{what}: dA/dr = -2W/3L")
    _zero(_numeric(d(scalar) - sp.diff(W, f).subs(f, scalar) / L, at), f"{what}: dPhi/dr = W'/L")
    _zero(_numeric(d(d(A)) + 2 * d(scalar) ** 2 / 3, at), f"{what}: A'' = -(2/3) Phi'^2")


def _boundary(chart, at, far, what):
    """The Ricci scalar near the boundary is anti-de Sitter's of radius L, -20/L^2."""
    L = chart.reader.parameters["L"]
    value = _numeric(chart.geo.ricci_scalar() + 20 / L ** 2, {**at, **far})
    if abs(value) > sp.Float("1e-20"):
        raise AssertionError(f"{METRIC}: {what}: R + 20/L^2 near the boundary is {value}")


def chart_check(chart, system):
    {"proper_distance": _check_proper_distance, "bfs_u": _check_bfs_u, "coulomb_disc": _check_coulomb}[system](chart)


POINTS = [{"L": sp.Rational(3, 2), "x": sp.Rational(7, 10)}, {"L": sp.Integer(2), "x": sp.Integer(3)},
          {"L": sp.Integer(1), "x": sp.Rational(1, 5)}]


def _check_proper_distance(chart):
    r = chart.symbols[-1]
    L = chart.reader.parameters["L"]
    y = r / L
    # GPPZ99 (m38) at C_1 = 0: the warp factor e^{2 phi} and the scalar m.
    warp = sp.exp(2 * (y + sp.log(2 * sp.sinh(y))) / 2)
    m = sp.sqrt(3) / 2 * sp.log((1 + sp.exp(-y)) / (1 - sp.exp(-y)))
    for p in POINTS:
        at = {L: p["L"], r: p["x"]}
        for i in range(4):
            _zero(_numeric(chart.geo.g[i, i] - (-1 if i == 0 else 1) * warp, at), f"GPPZ99 (m38), slot {i}{i}")
        _zero(_numeric(chart.geo.g[4, 4] - 1, at), "r is the proper distance")
        _flow(chart, at, m, _gppz_superpotential, L, True, "GPPZ99 (m37)")
        _field_equations(chart, at, m, _gppz_superpotential, L, "GPPZ99")
    _boundary(chart, {L: sp.Rational(3, 2)}, {r: 45}, "the proper distance chart")


def _check_bfs_u(chart):
    u = chart.symbols[-1]
    L = chart.reader.parameters["L"]
    root = sp.sqrt(1 - u)
    scalar = sp.sqrt(3) / 2 * sp.log((1 + root) / (1 - root))
    # The proper distance r = -(L/2) log(1 - u), from du/dr = 2(1 - u)/L.
    radius = -L * sp.log(1 - u) / 2
    for p in POINTS:
        at = {L: p["L"], u: 1 - sp.exp(-2 * p["x"] / p["L"])}
        for i in range(4):
            _zero(_numeric(chart.geo.g[i, i] - (-1 if i == 0 else 1) * u / (1 - u), at), f"BFS Sec. 3.2, slot {i}{i}")
        _zero(_numeric(1 / sp.sqrt(chart.geo.g[4, 4]) - 2 * (1 - u) / L, at), "BFS Sec. 3.2, du/dr = 2(1 - u)/L")
        # The proper distance chart pulled back through u = 1 - e^{-2r/L}: e^{2r/L} - 1 = u/(1 - u), g_uu = (dr/du)^2.
        _zero(_numeric(sp.exp(2 * radius / L) - 1 - u / (1 - u), at), "u = 1 - e^{-2r/L}, the warp factor")
        _zero(_numeric(chart.geo.g[4, 4] - sp.diff(radius, u) ** 2, at), "u = 1 - e^{-2r/L}, g_uu")
        # BFS's scalar is GPPZ99's m at y = r/L.
        _zero(_numeric(scalar - sp.sqrt(3) * sp.atanh(sp.exp(-radius / L)), at), "BFS's scalar is GPPZ99's m")
        _flow(chart, at, scalar, _gppz_superpotential, L, False, "BFS (1storder)")
        _field_equations(chart, at, scalar, _gppz_superpotential, L, "BFS Sec. 3.2")
    _boundary(chart, {L: sp.Rational(3, 2)}, {u: 1 - sp.exp(-60)}, "the chart of u")


def _check_coulomb(chart):
    v = chart.symbols[-1]
    L, ell = chart.reader.parameters["L"], chart.reader.parameters["ell"]
    scalar = sp.log(v) / sp.sqrt(6)
    mu = sp.Symbol("mu", real=True)
    W_fgpw = -sp.exp(2 * mu / sp.sqrt(6)) - sp.exp(-4 * mu / sp.sqrt(6)) / 2
    A_fgpw = sp.log(sp.Abs(sp.exp(2 * mu / sp.sqrt(6)) / (1 - sp.exp(sp.sqrt(6) * mu)))) / 2 + sp.log(ell / L)
    points = [{L: sp.Rational(3, 2), ell: sp.Integer(2), v: sp.Rational(1, 3)},
              {L: sp.Integer(1), ell: sp.Rational(1, 2), v: sp.Rational(4, 5)},
              {L: sp.Integer(2), ell: sp.Integer(3), v: sp.Rational(1, 20)}]
    for at in points:
        warp = ell ** 2 * v ** sp.Rational(2, 3) / (L ** 2 * (1 - v))
        for i in range(4):
            _zero(_numeric(chart.geo.g[i, i] - (-1 if i == 0 else 1) * warp, at), f"BFS Sec. 3.1, slot {i}{i}")
        _zero(_numeric(1 / sp.sqrt(chart.geo.g[4, 4]) - 2 * v ** sp.Rational(2, 3) * (1 - v) / L, at),
              "BFS Sec. 3.1, dv/dr = (2/L) v^{2/3}(1 - v)")
        # FGPW00 (SUSYCases), n = 2, at mu = -log(v)/sqrt6: their A(mu) and their first order equations
        # dmu/drho = (g/2) dW/dmu and dA/drho = -(g/3) W, g = 2/L, of (SUSYCasesstart).
        here = {mu: -sp.log(at[v]) / sp.sqrt(6), L: at[L], ell: at[ell]}
        _zero(_numeric(sp.log(chart.geo.g[1, 1]) / 2, at) - _numeric(A_fgpw, here), "FGPW00 (SUSYCases), A(mu)")
        drho = 2 * v ** sp.Rational(2, 3) * (1 - v) / L
        _zero(_numeric(-drho / (sp.sqrt(6) * v), at) - _numeric(sp.diff(W_fgpw, mu) / L, here),
              "FGPW00 (SUSYCasesstart), dmu/drho")
        A = sp.log(chart.geo.g[1, 1]) / 2
        _zero(_numeric(sp.diff(A, v) * drho, at) + _numeric(2 * W_fgpw / (3 * L), here),
              "FGPW00 (SUSYCasesstart), dA/drho")
        _flow(chart, at, scalar, _coulomb_superpotential, L, False, "BFS (1storder)")
        _field_equations(chart, at, scalar, _coulomb_superpotential, L, "BFS Sec. 3.1")
    _boundary(chart, {L: sp.Rational(3, 2), ell: 2}, {v: 1 - sp.exp(-60)}, "the disc chart")
