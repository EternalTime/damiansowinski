"""Gravitational waves in Friedmann universes filled with a stiff fluid, in units G = c = 1.

Bicak and Griffiths (1996), arXiv gr-qc/9610041, (1.1), (2.4), (2.5), (3.15), (3.16), (4.7),
(5.2) and (5.4)-(5.8), give the metric form and the five stiff backgrounds; their (7.1)-(7.4)
give the plane wave in the spatially flat universe. Alekseev and Griffiths (1995), Phys. Rev.
D 52, 4497, (2.12)-(2.13) and (4.1)-(4.4), give the profile H_2 and the closed form of M for
one term. Every line element here is theirs multiplied by -1, for the signature (-,+,+,+).
See waves_in_friedmann.md beside this file for the derivation and the choices.
"""

METRIC = "waves_in_friedmann"
CHARTS = ("plane_wave", "flat_planar", "flat_cylindrical", "open_planar", "open_cylindrical",
          "closed_toroidal")
NAMES = {"plane_wave": "Plane Wave, Flat Universe", "flat_planar": "Flat Universe, Planar",
         "flat_cylindrical": "Flat Universe, Cylindrical", "open_planar": "Open Universe, Planar",
         "open_cylindrical": "Open Universe, Cylindrical", "closed_toroidal": "Closed Universe, Toroidal"}
COORDS = {"plane_wave": ["\\eta", "z", "x", "y"], "flat_planar": ["\\eta", "z", "x", "y"],
          "flat_cylindrical": ["\\eta", "r", "\\phi", "z"], "open_planar": ["\\eta", "\\mu", "x", "y"],
          "open_cylindrical": ["\\eta", "r", "\\phi", "y"],
          "closed_toroidal": ["\\eta", "\\zeta", "\\delta", "\\sigma"]}

# The profile H = H_2(z/eta) of Alekseev and Griffiths's (4.3), normalised to H_2(-1) = 1, its
# slope K = dH_2/d(z/eta), the wave V and the function S of e^{-M} = 2 eta e^{-S}.
H_DEF = r"H = \dfrac{\left(\eta^2+2z^2\right)\arccos(z/\eta)-3z\sqrt{\eta^2-z^2}}{3\pi\eta^2}"
K_DEF = r"K = \dfrac{4\left(z\arccos(z/\eta)-\sqrt{\eta^2-z^2}\right)}{3\pi\eta}"
V_DEF = r"V = A\eta^2H"
S_DEF = r"S = -\dfrac{A^2\eta^2}{8}\left(4\eta^2H^2+\left(\eta^2-z^2\right)K^2\right)"
WAVE_PARAMETERS = ["\\gamma", "A", H_DEF, K_DEF, V_DEF, S_DEF]

DIMENSIONS = {(METRIC, s): {**{x: "1" for x in COORDS[s]}, "\\gamma": "L**2",
                            **({"A": "1", "H": "1", "K": "1", "V": "1", "S": "1"} if s == "plane_wave" else {})}
              for s in CHARTS}
HELD = {(METRIC, "plane_wave"): ("H", "K", "V", "S")}
# First derivatives in the names: H' = K and, from the Euler-Poisson-Darboux equation of
# Alekseev and Griffiths's (4.2) with n = 2, (1-w^2)K' + 3wK - 4H = 0 for w = z/eta; then
# V_eta = A(2 eta H - zK), V_z = A eta K and Bicak and Griffiths's (6.7) for S.
RATES = {(METRIC, "plane_wave"): {
    "H": {"\\eta": r"-zK/\eta^2", "z": r"K/\eta"},
    "K": {"\\eta": r"-\dfrac{z\left(4\eta H-3zK\right)}{\eta\left(\eta^2-z^2\right)}",
          "z": r"\dfrac{4\eta H-3zK}{\eta^2-z^2}"},
    "V": {"\\eta": r"A\left(2\eta H-zK\right)", "z": r"A\eta K"},
    "S": {"\\eta": r"-\dfrac{A^2\eta}{2}\left(\left(2\eta H-zK\right)^2+\eta^2K^2\right)",
          "z": r"-A^2\eta^2K\left(2\eta H-zK\right)"},
}}

LINES = {
    "plane_wave": r"ds^2 = 2\gamma\eta\left(e^{-S}\left(-d\eta^2+dz^2\right)+e^{V}dx^2+e^{-V}dy^2\right)",
    "flat_planar": r"ds^2 = 2\gamma\eta\left(-d\eta^2+dz^2+dx^2+dy^2\right)",
    "flat_cylindrical": r"ds^2 = 2\gamma\eta\left(-d\eta^2+dr^2+r^2d\phi^2+dz^2\right)",
    "open_planar": r"ds^2 = \gamma\sinh(2\eta)\left(-d\eta^2+d\mu^2+e^{-2\mu}\left(dx^2+dy^2\right)\right)",
    "open_cylindrical": r"ds^2 = \gamma\sinh(2\eta)\left(-d\eta^2+dr^2+\sinh^2r\,d\phi^2+\cosh^2r\,dy^2\right)",
    "closed_toroidal": r"ds^2 = \gamma\sin(2\eta)\left(-d\eta^2+d\zeta^2+\sin^2\zeta\,d\delta^2+\cos^2\zeta\,d\sigma^2\right)",
}
REALS = r" \in (-\infty,\infty)"
DOMAINS = {
    "plane_wave": [r"\eta \in (0,\infty)", r"z \in (-\eta,\eta)", "x"+REALS, "y"+REALS,
                   r"z=\eta\;\text{(wave front)}", r"\eta\to0\;\text{(big bang)}"],
    "flat_planar": [r"\eta \in (0,\infty)", "z"+REALS, "x"+REALS, "y"+REALS,
                    r"\eta\to0\;\text{(big bang)}"],
    "flat_cylindrical": [r"\eta \in (0,\infty)", r"r \in [0,\infty)", r"\phi \in [0,2\pi)", "z"+REALS,
                         r"\eta\to0\;\text{(big bang)}"],
    "open_planar": [r"\eta \in (0,\infty)", r"\mu"+REALS, "x"+REALS, "y"+REALS,
                    r"\eta\to0\;\text{(big bang)}"],
    "open_cylindrical": [r"\eta \in (0,\infty)", r"r \in [0,\infty)", r"\phi \in [0,2\pi)", "y"+REALS,
                         r"\eta\to0\;\text{(big bang)}"],
    "closed_toroidal": [r"\eta \in (0,\pi/2)", r"\zeta \in [0,\pi/2]", r"\delta \in [0,2\pi)",
                        r"\sigma \in [0,2\pi)", r"\eta\to0\;\text{(big bang)}", r"\eta\to\pi/2\;\text{(big crunch)}"],
}
# The scale factor squared, a^2 = gamma f(eta), of Bicak and Griffiths's (5.2), and the fluid
# potential of their (5.11)-(5.13), with Ricci = (1/2) d(sigma) d(sigma): their (5.9) reads
# 16 pi rho for what the Friedmann normalisation of their (5.1) and BG94 (3.1)-(3.4) makes 32 pi rho.
CURVATURE = {"plane_wave": "flat", "flat_planar": "flat", "flat_cylindrical": "flat",
             "open_planar": "open", "open_cylindrical": "open", "closed_toroidal": "closed"}


def chart_spec(system):
    """One chart of Bicak and Griffiths's stiff universes, or the plane wave behind its front."""
    import verify_metrics as vm

    key = (METRIC, system)
    params = WAVE_PARAMETERS if system == "plane_wave" else ["\\gamma"]
    probe = vm.Reader(COORDS[system], params, (), held=HELD.get(key, ()), rates=RATES.get(key))
    spec = {"metric_id": METRIC,
            "system": {"id": system, "name": NAMES[system], "coords": COORDS[system],
                       "domains": DOMAINS[system], "parameters": params, "line_element": LINES[system]},
            "chart_line_element": LINES[system], "printer": {},
            "check": lambda chart: chart_check(chart, system)}
    if system == "plane_wave":
        spec["reduce"] = lambda value: vm.norm(probe.by_rates(value))
    elif CURVATURE[system] == "open":
        import chart_printer as cp
        spec["pretty"] = cp.hyperbolic(probe.symbol["\\eta"])
        spec["ricci_scalar"] = r"-\dfrac{6}{\gamma\sinh^3(2\eta)}"
        spec["kretschmann"] = r"\dfrac{60}{\gamma^2\sinh^6(2\eta)}"
    elif CURVATURE[system] == "closed":
        spec["pretty"] = trigonometric(probe.symbol["\\eta"])
        spec["ricci_scalar"] = r"-\dfrac{6}{\gamma\sin^3(2\eta)}"
        spec["kretschmann"] = r"\dfrac{60}{\gamma^2\sin^6(2\eta)}"
    return spec


def trigonometric(x):
    """The closed universe's values, which the checker hands back in sin x and cos x, rewritten
    as chart_printer.hyperbolic rewrites sinh and cosh: reduced by cos^2 = 1 - sin^2, the
    denominator cleared of cos by its conjugate, factored, and 1 - sin^2 written cos^2."""
    import sympy as sp

    S, C = sp.symbols("_S _C", positive=True)
    relation = [C ** 2 - 1 + S ** 2]

    def reduce(p):
        return sp.expand(sp.reduced(sp.expand(p), relation, C, S)[1])

    def pretty(value):
        v = sp.sympify(value)
        v = v.subs({sp.tan(x): sp.sin(x)/sp.cos(x), sp.cot(x): sp.cos(x)/sp.sin(x),
                    sp.sec(x): 1/sp.cos(x), sp.csc(x): 1/sp.sin(x)}).subs({sp.sin(x): S, sp.cos(x): C})
        num, den = (reduce(q) for q in sp.fraction(sp.together(v)))
        d0, d1 = sp.Poly(den, C).coeff_monomial(1), sp.Poly(den, C).coeff_monomial(C)
        if d1 != 0:
            num, den = reduce(num * (d0 - d1 * C)), sp.expand(d0 ** 2 - d1 ** 2 * (1 - S ** 2))
        out = sp.factor(sp.cancel(sp.factor(num) / sp.factor(den)))
        # (S - 1)^k (S + 1)^k, which factor makes of 1 - S^2, is read as (-C^2)^k and cancelled again.
        powers = sp.Mul.make_args(out)
        pair = {}
        for f in powers:
            base, k = (f.base, f.exp) if f.is_Pow else (f, 1)
            if sp.expand(base - (S - 1)) == 0 or sp.expand(base - (S + 1)) == 0:
                pair[sp.expand(base)] = k
        if len(pair) == 2 and len(set(pair.values())) == 1:
            k = next(iter(pair.values()))
            rest = [f for f in powers if sp.expand(f.base if f.is_Pow else f) not in pair]
            out = sp.factor(sp.Mul(*rest) * (-C ** 2) ** k)
        return out.subs({S: sp.sin(x), C: sp.cos(x)})
    return pretty


def chart_check(chart, system):
    """Every Ricci slot against the stiff fluid, and the held slopes, at points; no simplify."""
    import sympy as sp

    R = chart.reader
    X = chart.symbols
    eta = X[0]
    kind = CURVATURE[system]
    sigma = sp.sqrt(3)*(sp.log(eta) if kind == "flat" else sp.log(sp.tanh(eta)) if kind == "open"
                        else sp.log(sp.tan(eta)))
    gradient = [sp.diff(sigma, c) for c in X]
    ricci = chart.geo.ricci_ll()
    gamma = R.parameters["gamma"]
    points = ((sp.Rational(7, 10), sp.Rational(1, 5)), (sp.Rational(11, 10), sp.Rational(-1, 2)),
              (sp.Rational(1, 2), sp.Rational(1, 3)))
    for k, (ev, sv) in enumerate(points):
        at = {eta: ev, X[1]: sv if system != "closed_toroidal" else sp.Rational(1, 2)+sv/3, X[2]: sp.Rational(1, 3),
              X[3]: sp.Rational(1, 4), gamma: (1, sp.Rational(5, 2), sp.Rational(7, 10))[k]}
        if system == "plane_wave":
            at[R.parameters["A"]] = (sp.Rational(13, 10), sp.Rational(3, 5), 2)[k]

        def number(value):
            return sp.N(sp.sympify(value).subs(R.held).doit().subs(at), 40)
        for i in range(4):
            for j in range(4):
                gap = number(ricci[i][j]) - number(gradient[i]*gradient[j]/2)
                if abs(gap) > 1e-25:
                    raise AssertionError(f"{METRIC}/{system}: Ricci slot {i},{j} misses the stiff fluid by {gap}")
        for name, rates in RATES.get((METRIC, system), {}).items():
            for coord, rate in rates.items():
                gap = number(R(rate)) - number(sp.diff(R.held[R.parameters[name]], R.symbol[coord]))
                if abs(gap) > 1e-25:
                    raise AssertionError(f"{METRIC}/{system}: held slope of {name} along {coord} is off by {gap}")
