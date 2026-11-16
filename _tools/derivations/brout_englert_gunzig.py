"""The causal universe of Brout, Englert and Gunzig: empty flat space outside the forward
light cone of one event, and inside it two stages in three coordinate systems.

The fireball is the creation stage with the correction delta to the density of created
particles set to zero, and the free expansion is the radiation era of the open Friedmann
universe that follows it. Each stage is written in the inertial coordinates of the flat
background, in the proper time of the comoving world lines and in the conformal time of
Spindel and Vandepeutte's (18). brout_englert_gunzig.md records the source of each chart and
derives the joining.

Every chart computes in x^0 = c times its time, so a parameter that is a time enters the
chart's line element multiplied by c.
"""
import sympy as sp

METRIC = "brout_englert_gunzig"
CHARTS = ("fireball_flat", "expansion_flat", "outside_flat", "fireball_proper", "expansion_proper",
          "fireball_conformal", "expansion_conformal")
NAMES = {
    "outside_flat": "Empty Space, Flat Background",
    "fireball_flat": "Fireball, Flat Background",
    "expansion_flat": "Free Expansion, Flat Background",
    "fireball_proper": "Fireball, Proper Time",
    "expansion_proper": "Free Expansion, Proper Time",
    "fireball_conformal": "Fireball, Conformal Time",
    "expansion_conformal": "Free Expansion, Conformal Time",
}
ANGLES = [r"\theta \in [0, \pi]", r"\varphi \in [0, 2\pi)"]
FLAT = r"-c^2dt^2 + dr^2 + r^2d\theta^2 + r^2\sin^2\theta\,d\varphi^2"
OPEN = r"d\chi^2 + \sinh^2\chi\,d\theta^2 + \sinh^2\chi\sin^2\theta\,d\varphi^2"
TAU = r"\tau = \sqrt{t^2 - \dfrac{r^2}{c^2}}"

# Each chart's published line element and, where its time is not the chart's x^0, the line
# element the tensors are computed from. The flat charts compute from t and r written out and
# print in the flat proper time tau, which the reader reads by its published definition.
SYSTEMS = {
    "fireball_flat": {
        "coords": ["t", "r", r"\theta", r"\varphi"],
        "parameters": [r"\tau_\infty", r"\tau_1", TAU],
        "domains": [r"t \in (0, \infty)", r"r \in [0, ct)", *ANGLES, r"\tau \in (0, \tau_1]",
                    r"\tau_1 < \tau_\infty",
                    r"\tau = 0 \;\text{(the light cone of the origin)}",
                    r"\tau = \tau_1 \;\text{(the end of creation)}"],
        "line": r"ds^2 = \dfrac{" + FLAT + r"}{\left(1 - \dfrac{\tau^2}{\tau_\infty^2}\right)^2}",
        "chart": r"ds^2 = \dfrac{" + FLAT.replace("c^2dt", "dt")
                 + r"}{\left(1 - \dfrac{t^2 - r^2}{c^2\tau_\infty^2}\right)^2}",
        "time": "t",
    },
    "expansion_flat": {
        "coords": ["t", "r", r"\theta", r"\varphi"],
        "parameters": [r"\Omega_f", r"\tau_b", r"\tau_\infty", r"\tau_1", TAU],
        "domains": [r"t \in (0, \infty)", r"r \in [0, ct)", *ANGLES, r"\tau \in [\tau_1, \infty)",
                    r"\tau_b < \tau_1 < \tau_\infty",
                    r"\tau = \tau_1 \;\text{(the end of creation)}"],
        "line": r"ds^2 = \Omega_f^2\left(1 - \dfrac{\tau_b^2}{\tau^2}\right)^2\left(" + FLAT + r"\right)",
        "chart": r"ds^2 = \Omega_f^2\left(1 - \dfrac{c^2\tau_b^2}{t^2 - r^2}\right)^2\left("
                 + FLAT.replace("c^2dt", "dt") + r"\right)",
        "time": "t",
    },
    "outside_flat": {
        "coords": ["t", "r", r"\theta", r"\varphi"],
        "parameters": [],
        "domains": [r"t \in (-\infty, \infty)", r"r \in [0, \infty)", *ANGLES,
                    r"r \ge ct",
                    r"r = ct \;\text{with}\; t > 0 \;\text{(the light cone of the origin)}"],
        "line": "ds^2 = " + FLAT,
        "chart": "ds^2 = " + FLAT.replace("c^2dt", "dt"),
    },
    "fireball_proper": {
        "coords": ["T", r"\chi", r"\theta", r"\varphi"],
        "parameters": ["H", "T_1"],
        "domains": [r"T \in (0, T_1]", r"\chi \in [0, \infty)", *ANGLES,
                    r"T = 0 \;\text{(the light cone of the origin)}",
                    r"T = T_1 \;\text{(the end of creation)}"],
        "line": r"ds^2 = -c^2dT^2 + \dfrac{c^2\sinh^2\left(HT\right)}{H^2}\left(" + OPEN + r"\right)",
        "chart": r"ds^2 = -dT^2 + \dfrac{c^2\sinh^2\left(HT/c\right)}{H^2}\left(" + OPEN + r"\right)",
        "time": "T",
    },
    "expansion_proper": {
        "coords": ["T", r"\chi", r"\theta", r"\varphi"],
        "parameters": ["T_b", "b", "H", "T_1"],
        "domains": [r"T \in [T_1, \infty)", r"\chi \in [0, \infty)", *ANGLES,
                    r"T_b + b < T_1",
                    r"T = T_1 \;\text{(the end of creation)}"],
        "line": r"ds^2 = -c^2dT^2 + c^2\left(\left(T - T_b\right)^2 - b^2\right)\left(" + OPEN + r"\right)",
        "chart": r"ds^2 = -dT^2 + \left(\left(T - cT_b\right)^2 - c^2b^2\right)\left(" + OPEN + r"\right)",
        "time": "T",
    },
    "fireball_conformal": {
        "coords": [r"\eta", r"\chi", r"\theta", r"\varphi"],
        "parameters": ["H", r"\eta_1"],
        "domains": [r"\eta \in (-\infty, \eta_1]", r"\chi \in [0, \infty)", *ANGLES, r"\eta_1 < 0",
                    r"\eta \to -\infty \;\text{(the light cone of the origin)}",
                    r"\eta = \eta_1 \;\text{(the end of creation)}"],
        "line": r"ds^2 = \dfrac{c^2}{H^2\sinh^2\eta}\left(-d\eta^2 + " + OPEN + r"\right)",
    },
    "expansion_conformal": {
        "coords": [r"\eta", r"\chi", r"\theta", r"\varphi"],
        "parameters": ["b", r"\eta_b", "H", r"\eta_1"],
        "domains": [r"\eta \in [\eta_1, \infty)", r"\chi \in [0, \infty)", *ANGLES,
                    r"\eta_b < \eta_1 < 0",
                    r"\eta = \eta_1 \;\text{(the end of creation)}"],
        "line": r"ds^2 = b^2c^2\sinh^2\left(\eta - \eta_b\right)\left(-d\eta^2 + " + OPEN + r"\right)",
    },
}

DIMENSIONS = {
    (METRIC, "fireball_flat"): {"t": "T", "r": "L", r"\theta": "1", r"\varphi": "1",
                                r"\tau_\infty": "T", r"\tau_1": "T", r"\tau": "T"},
    (METRIC, "expansion_flat"): {"t": "T", "r": "L", r"\theta": "1", r"\varphi": "1", r"\Omega_f": "1",
                                 r"\tau_b": "T", r"\tau_\infty": "T", r"\tau_1": "T", r"\tau": "T"},
    (METRIC, "outside_flat"): {"t": "T", "r": "L", r"\theta": "1", r"\varphi": "1"},
    (METRIC, "fireball_proper"): {"T": "T", r"\chi": "1", r"\theta": "1", r"\varphi": "1",
                                  "H": "1/T", "T_1": "T"},
    (METRIC, "expansion_proper"): {"T": "T", r"\chi": "1", r"\theta": "1", r"\varphi": "1",
                                   "T_b": "T", "b": "T", "H": "1/T", "T_1": "T"},
    (METRIC, "fireball_conformal"): {r"\eta": "1", r"\chi": "1", r"\theta": "1", r"\varphi": "1",
                                     "H": "1/T", r"\eta_1": "1"},
    (METRIC, "expansion_conformal"): {r"\eta": "1", r"\chi": "1", r"\theta": "1", r"\varphi": "1",
                                      "b": "T", r"\eta_b": "1", "H": "1/T", r"\eta_1": "1"},
}


# -- the joining ----------------------------------------------------------------------------

def joined(tau_infinity, tau_1):
    """Every constant of the free expansion that joins the fireball at tau_1 with the conformal
    factor and its first derivative continuous, from the two free times; exact for exact input.

    With x = tau_1/tau_infinity: Omega_1 = 1/(1 - x^2), Omega_f = Omega_1^2 and
    tau_b = x tau_1, so eta_b = 2 eta_1, b = sinh^2(H T_1)/H and T_1 - T_b = sinh(2 H T_1)/(2H).
    """
    x = sp.sympify(tau_1)/tau_infinity
    H = 2/sp.sympify(tau_infinity)
    T_1 = tau_infinity*sp.atanh(x)
    b = sp.sinh(H*T_1)**2/H
    return {"H": H, "Omega_1": 1/(1 - x**2), "Omega_f": 1/(1 - x**2)**2, "tau_b": x*tau_1,
            "T_1": T_1, "b": b, "T_b": T_1 - sp.sinh(2*H*T_1)/(2*H),
            "eta_1": sp.log(x), "eta_b": 2*sp.log(x)}


def factor(stage, tau, tau_infinity, tau_1):
    """The square root of the conformal factor on the flat background, at the flat proper time tau."""
    if stage == "fireball":
        return 1/(1 - (tau/tau_infinity)**2)
    k = joined(tau_infinity, tau_1)
    return k["Omega_f"]*(1 - (k["tau_b"]/tau)**2)


# -- printing -------------------------------------------------------------------------------

def conjugates(value, written=None):
    """A product with each pair of factors (A + B)^k (A - B)^k written as one factor to the
    power k: A^2 - B^2 multiplied out, or `written(A, B)` where the chart names the difference
    its own way."""
    items = [(f.base, f.exp) if f.is_Pow else (f, sp.Integer(1)) for f in sp.Mul.make_args(sp.sympify(value))]
    out, used = [], set()
    for i, (first, k) in enumerate(items):
        if i in used:
            continue
        match = None
        if first.is_Add:
            terms = set(sp.Add.make_args(first))
            for j in range(i + 1, len(items)):
                second, other = items[j]
                if j in used or other != k or not second.is_Add:
                    continue
                against = set(sp.Add.make_args(second))
                same, flipped = terms & against, {t for t in terms if -t in against}
                if same and flipped and len(same) + len(flipped) == len(terms) == len(against):
                    match = (j, sp.Add(*same), sp.Add(*flipped))
                    break
        if match:
            j, A, B = match
            used.add(j)
            out.append((written(A, B) if written else sp.expand(A**2 - B**2))**k)
        else:
            out.append(first**k)
    return sp.Mul(*out)


def hyperbolic(*arguments, shift=None, then=None):
    """A `pretty` for a chart with sinh and cosh of more than one argument, which the checker's
    Geometry hands back as exponentials. Each value is rewritten in S_i = sinh u_i and
    C_i = cosh u_i for the arguments u_i given, reduced by C_i^2 = 1 + S_i^2, its denominator
    cleared of each C_i by its conjugate, and factored, with every factor 1 + S_i^2 written
    C_i^2. chart_printer.hyperbolic does the same for one argument.

    `shift` is (x, x_0) for an argument x - x_0: the exponentials of x and of x_0 are then
    written as those of the difference, and the powers of e^{x_0} left over cancel.
    `then` is applied to the factored value before the functions are written back."""
    pairs = [(u, *sp.symbols(f"_S{i} _C{i}", positive=True)) for i, u in enumerate(arguments)]
    relations = [C**2 - 1 - S**2 for _, S, C in pairs]
    generators = [C for _, _, C in pairs] + [S for _, S, _ in pairs]

    def reduce(p):
        return sp.expand(sp.reduced(sp.expand(p), relations, *generators)[1])

    def power(e):
        """exp of a sum of whole multiples of the arguments as a product of (C + S) and (C - S)."""
        for u, S, C in pairs:
            k = sp.cancel(e.args[0]/u)
            if k.is_Integer:
                return (C + S)**k if k > 0 else (C - S)**(-k)
        out = sp.Integer(1)
        for term in sp.Add.make_args(sp.expand(e.args[0])):
            for u, S, C in pairs:
                k = sp.nsimplify(term/u) if (term/u).is_number else None
                if k is not None and k.is_Integer:
                    out *= (C + S)**k if k > 0 else (C - S)**(-k)
                    break
            else:
                return e
        return out

    def pretty(value):
        x = sp.sympify(value)
        back = {}
        if shift:
            moved, base, gap = *shift, sp.Dummy("u")
            scale = sp.Dummy("E", positive=True)
            x = x.subs(moved, gap + base).replace(
                lambda e: isinstance(e, sp.exp),
                lambda e: sp.Mul(*[scale**(t/base) if (t/base).is_Integer else sp.exp(t)
                                   for t in sp.Add.make_args(sp.expand(e.args[0]))]))
            x = sp.cancel(sp.together(x)).subs(gap, moved - base)
            back = {scale: sp.exp(base)}
        x = x.replace(lambda e: isinstance(e, sp.exp), power)
        num, den = (reduce(p) for p in sp.fraction(sp.together(x)))
        for _, S, C in pairs:
            poly = sp.Poly(den, C)
            d0, d1 = poly.coeff_monomial(1), poly.coeff_monomial(C)
            if d1 != 0:
                num, den = reduce(num*(d0 - d1*C)), reduce(d0**2 - d1**2*(1 + S**2))
        out = sp.factor(sp.cancel(sp.factor(num)/sp.factor(den)))
        for _, S, C in pairs:
            out = out.replace(lambda e, S=S: isinstance(e, sp.Expr) and sp.expand(e - (S**2 + 1)) == 0,
                              lambda e, C=C: C**2)
        if then:
            out = then(out)
        return out.subs({s: f(u) for u, S, C in pairs for s, f in ((S, sp.sinh), (C, sp.cosh))}).subs(back)
    return pretty


TAU_SYMBOL = sp.Symbol("_tau", positive=True)
SQUARE_SYMBOL = sp.Symbol("_square", positive=True)


def in_tau(t, r, c):
    """A `pretty` for the flat charts, whose values depend on t and r through
    tau^2 = t^2 - r^2/c^2 and its derivatives. Each value is written twice, with every even
    power of r written in t and tau and with every even power of t written in r and tau, and
    the shorter of the two factored forms is kept. TAU_SYMBOL stands for tau and prints as it."""
    s = TAU_SYMBOL
    swaps = ((r, c**2*(t**2 - s**2)), (t, s**2 + r**2/c**2))

    def written(p, letter, square):
        poly = sp.Poly(sp.expand(p), letter)
        return sp.expand(sum(coefficient*letter**(k % 2)*square**(k//2) for (k,), coefficient in poly.terms()))

    def pretty(value):
        num, den = sp.fraction(sp.together(sp.sympify(value)))
        forms = [conjugates(sp.factor(sp.cancel(written(num, letter, square)/written(den, letter, square))))
                 for letter, square in swaps]
        return min(forms, key=lambda form: (sp.count_ops(form), str(form)))
    return pretty


def chart_spec(system):
    import verify_metrics as vm

    chart = SYSTEMS[system]
    internal = chart.get("chart_parameters", chart["parameters"])
    probe = vm.Reader(chart["coords"], internal, ())
    first, chi = probe.symbol[chart["coords"][0]], probe.symbol.get(r"\chi")
    P = probe.parameters
    spec = {
        "metric_id": METRIC,
        "system": {"id": system, "name": NAMES[system], "coords": chart["coords"],
                   "domains": chart["domains"], "parameters": chart["parameters"],
                   "line_element": chart["line"]},
        "chart_parameters": internal,
        "chart_line_element": chart.get("chart", chart["line"]),
        "check": lambda c: check_chart(c, system),
    }
    if "time" in chart:
        spec["time"] = chart["time"]
    names = [P[name] for name in ("Omega_f", "tau_infty", "tau_b", "H", "b", "T_b", "eta_b") if name in P]
    overrides = {P["tau_infty"]: r"\tau_\infty"} if "tau_infty" in P else {}
    if system == "outside_flat":
        spec["printer"] = {"lead": [probe.symbol["r"]]}
        return spec
    if system.endswith("_flat"):
        r = probe.symbol["r"]
        spec["pretty"] = in_tau(first, r, probe.c)
        overrides[TAU_SYMBOL] = r"\tau"
        # tau_infinity^2 - tau^2 in the fireball and tau^2 - tau_b^2 in the free expansion, each positive.
        order = [P["tau_infty"], TAU_SYMBOL] if system == "fireball_flat" else [TAU_SYMBOL, P["tau_b"]]
        spec["printer"] = {"lead": [*order, first, r], "factors": [probe.c, *names, TAU_SYMBOL, first, r],
                           "overrides": overrides}
        return spec
    spec["printer"] = {"lead": [first, *names], "factors": [probe.c, *names, first], "overrides": overrides}
    # Every value is printed with x^0 already written as c times the time, so H T is the argument.
    if system == "fireball_proper":
        spec["pretty"] = hyperbolic(P["H"]*first, chi)
    elif system == "fireball_conformal":
        spec["pretty"] = hyperbolic(first, chi)
    elif system == "expansion_conformal":
        spec["pretty"] = hyperbolic(first - P["eta_b"], chi, shift=(first, P["eta_b"]))
    elif system == "expansion_proper":
        # (T - T_b)^2 - b^2, the square of the scale factor over c^2, as the line element writes it.
        spec["printer"]["named"] = {SQUARE_SYMBOL: r"\left(T - T_b\right)^2 - b^2"}
        square = (first - P["T_b"])**2 - P["b"]**2
        spec["pretty"] = hyperbolic(chi, then=lambda value: conjugates(
            value, lambda A, B: SQUARE_SYMBOL if sp.expand(A**2 - B**2 - square) == 0 else -SQUARE_SYMBOL))
    return spec


# -- checks, all at points ------------------------------------------------------------------

SAMPLES = (
    # tau_infinity, tau_1, a fraction of the stage's range of tau, chi, theta, c
    (sp.Rational(5, 2), sp.Rational(2), sp.Rational(1, 3), sp.Rational(3, 5), sp.Rational(7, 10), sp.Rational(13, 10)),
    (sp.Rational(7, 4), sp.Rational(3, 4), sp.Rational(4, 5), sp.Rational(11, 10), sp.Rational(9, 5), sp.Rational(1)),
)


def sample(system, tau_infinity, tau_1, fraction, chi, theta, c):
    """A point of the chart, the values of its parameters at the joining, and the flat proper
    time of the point."""
    k = joined(tau_infinity, tau_1)
    tau = fraction*tau_1 if system.startswith("fireball") else tau_1*(1 + 2*fraction)
    values = {"tau_infty": tau_infinity, "tau_1": tau_1, "Omega_f": k["Omega_f"], "tau_b": k["tau_b"],
              "H": k["H"], "T_1": k["T_1"], "T_b": k["T_b"], "b": k["b"],
              "eta_1": k["eta_1"], "eta_b": k["eta_b"]}
    kind = system.rpartition("_")[2]
    if kind == "flat":
        first, second = c*tau*sp.cosh(chi), c*tau*sp.sinh(chi)
    elif kind == "proper":
        proper = (tau_infinity*sp.atanh(tau/tau_infinity) if system.startswith("fireball") else
                  k["T_1"] + k["Omega_f"]*((tau - tau_1) + k["tau_b"]**2*(1/tau - 1/tau_1)))
        first, second = c*proper, chi
    else:
        first, second = sp.log(tau/tau_infinity), chi
    return (first, second, theta, sp.Rational(2, 5)), values, tau


def check_chart(chart, system):
    """Each chart against its source's field equations and against the flat background chart.

    The fireball is an Einstein space with Lambda = 3H^2/c^2 and H = 2/tau_infinity.
    The free expansion has no Ricci scalar and the Einstein tensor of a fluid at rest on the
    hyperboloids with pressure a third of its energy density, 3 b^2 c^2/a^4 in the chart.
    Every chart is the flat background's conformal factor pulled back by
    ct = c tau cosh chi, r = c tau sinh chi, and the two stages meet at tau_1 with the scale
    factor and its first derivative continuous. No symbolic determinant or simplify is used.
    """
    X, P, c, geo = chart.symbols, chart.reader.parameters, chart.reader.c, chart.geo
    if system == "outside_flat":
        # Empty flat space: the inertial metric itself, with no curvature.
        t, r, theta, _ = X
        if geo.g != sp.diag(-1, 1, r**2, r**2*sp.sin(theta)**2) or geo.kretschmann() != 0:
            raise AssertionError(f"{METRIC}/{system}: the region outside the cone is not flat")
        return
    stage, kind = system.split("_")
    ricci, einstein = geo.ricci_ll(), geo.einstein_ll()
    scalar = geo.ricci_scalar()
    for tau_infinity, tau_1, fraction, chi, theta, light in SAMPLES:
        point, values, tau = sample(system, tau_infinity, tau_1, fraction, chi, theta, light)
        at = {**dict(zip(X, point)), c: light, **{P[n]: v for n, v in values.items() if n in P}}

        def number(value):
            return complex(sp.N(sp.sympify(value).subs(at), 40))

        def near(found, expected, what):
            if abs(number(found) - complex(sp.N(expected, 40))) > 1e-24*(1 + abs(complex(sp.N(expected, 40)))):
                raise AssertionError(f"{METRIC}/{system}: {what}")

        Omega = factor(stage, tau, tau_infinity, tau_1)
        a = light*tau*Omega
        # The flat background's metric pulled back onto this chart.
        if kind == "flat":
            pulled = Omega**2*sp.diag(-1, 1, point[1]**2, point[1]**2*sp.sin(theta)**2)
        else:
            space = a**2*sp.diag(1, sp.sinh(chi)**2, sp.sinh(chi)**2*sp.sin(theta)**2)
            pulled = sp.diag(-1 if kind == "proper" else -a**2, *space.diagonal())
        for i in range(4):
            for j in range(4):
                near(geo.g[i, j], pulled[i, j], f"the flat background pulled back, component {i}{j}")
        if stage == "fireball":
            lam = 3*(2/tau_infinity)**2/light**2
            for i in range(4):
                for j in range(4):
                    near(ricci[i][j] if isinstance(ricci, list) else ricci[i, j], lam*pulled[i, j],
                         f"Einstein space with H = 2/tau_infinity, component {i}{j}")
        else:
            near(scalar, 0, "the Ricci scalar of the radiation era")
            b = joined(tau_infinity, tau_1)["b"]
            density = 3*b**2*light**2/a**4
            # The fluid's velocity is along the flat proper time: u_mu = -d(c tau Omega...) normalised.
            if kind == "flat":
                up = sp.Matrix([point[0], point[1], 0, 0])/(light*tau*Omega)
            else:
                up = sp.Matrix([1 if kind == "proper" else 1/a, 0, 0, 0])
            down = pulled*up
            for i in range(4):
                for j in range(4):
                    fluid = density*(sp.Rational(4, 3)*down[i]*down[j] + pulled[i, j]/3)
                    near(einstein[i][j] if isinstance(einstein, list) else einstein[i, j], fluid,
                         f"radiation at rest on the hyperboloids, component {i}{j}")
    joining()


def joining():
    """The two stages at tau_1: the scale factor and its rate along the proper time agree, in
    each of the three coordinate systems, for the constants `joined` gives."""
    for tau_infinity, tau_1 in ((sp.Rational(5, 2), sp.Rational(2)), (sp.Rational(7, 4), sp.Rational(3, 4))):
        k = joined(tau_infinity, tau_1)
        tau, T, eta = sp.symbols("tau T eta", positive=True)
        H = k["H"]
        pairs = (
            (tau/(1 - tau**2/tau_infinity**2), k["Omega_f"]*(tau - k["tau_b"]**2/tau), tau, tau_1),
            (sp.sinh(H*T)/H, sp.sqrt((T - k["T_b"])**2 - k["b"]**2), T, k["T_1"]),
            (-1/(H*sp.sinh(eta)), k["b"]*sp.sinh(eta - k["eta_b"]), eta, k["eta_1"]),
        )
        for before, after, variable, where in pairs:
            for order in (0, 1):
                gap = sp.N((sp.diff(before, variable, order) - sp.diff(after, variable, order)).subs(variable, where), 40)
                if abs(gap) > 1e-30:
                    raise AssertionError(f"{METRIC}: the stages do not join at order {order} in {variable}")


# -- drawings -------------------------------------------------------------------------------
# Every drawing takes tau_1 = 0.8 tau_infinity and is drawn in units of 1/H and c/H, in
# which tau_infinity = 2.

RATIO = 0.8
END = 2
FIXED = {"theta": "pi/2", "varphi": "0"}


def drawn(system):
    """The parameter values a chart is drawn with, by the reader's plain names."""
    import math
    x = RATIO
    T_1 = 2*math.atanh(x)
    b = math.sinh(T_1)**2
    return {
        "outside_flat": {},
        "fireball_flat": {"tau_infty": END, "tau_1": x*END},
        "expansion_flat": {"Omega_f": 1/(1 - x*x)**2, "tau_b": x*x*END, "tau_infty": END, "tau_1": x*END},
        "fireball_proper": {"H": 1, "T_1": T_1},
        "expansion_proper": {"T_b": T_1 - math.sinh(2*T_1)/2, "b": b, "H": 1, "T_1": T_1},
        "fireball_conformal": {"H": 1, "eta_1": math.log(x)},
        "expansion_conformal": {"b": b, "eta_b": 2*math.log(x), "H": 1, "eta_1": math.log(x)},
    }[system]


def conformal_time(system, first):
    """eta = ln(tau/tau_infinity) along the chart's time, with the radius at the origin."""
    import numpy as np
    first = np.asarray(first, dtype=float)
    kind = system.rpartition("_")[2]
    if kind == "conformal":
        return first
    if kind == "flat":
        return np.log(first/END)
    values = drawn("expansion_proper")
    if system == "fireball_proper":
        return np.log(np.tanh(first/2))
    return 2*np.log(RATIO) + np.arccosh((first - values["T_b"])/values["b"])


PLANES = {
    # box (radius from, to, time from, to), the axes, and where the chart's region lies on the box
    "fireball_flat": ((0, 3, 0, 3), "$Hr/c$", "$Ht$", "(t**2 - r**2)*(tau_1**2 - t**2 + r**2)"),
    "expansion_flat": ((0, 6, 1, 7), "$Hr/c$", "$Ht$", "t**2 - r**2 - tau_1**2"),
    "outside_flat": ((0, 6, -2, 4), "$Hr/c$", "$Ht$", "r - t"),
    "fireball_proper": ((0, 2.4, 0, 2.4), r"$\chi$", "$HT$", None),
    "expansion_proper": ((0, 4, 2, 9), r"$\chi$", "$HT$", None),
    "fireball_conformal": ((0, 3, -3.1, -0.1), r"$\chi$", r"$\eta$", None),
    "expansion_conformal": ((0, 3, -0.3, 2.7), r"$\chi$", r"$\eta$", None),
}
PLANE = r" ($\theta = \pi/2$, $\varphi = 0$)"
STRAIGHT = "Every null ray is a straight line at 45°"
RAY_CAPTIONS = {
    "fireball_flat": [
        "The plane of $t$ and $r$" + PLANE + r", with $\tau_1 = 0.8\,\tau_\infty$. The fireball fills the "
        r"region between the light cone of the origin, $r = ct$, and the hyperbola $\tau = \tau_1$ on which "
        "creation ends. " + STRAIGHT + "."],
    "expansion_flat": [
        "The plane of $t$ and $r$" + PLANE + r", with $\tau_1 = 0.8\,\tau_\infty$. The free expansion fills "
        r"the region above the hyperbola $\tau = \tau_1$, where it joins the fireball. " + STRAIGHT + "."],
    "outside_flat": [
        "The plane of $t$ and $r$" + PLANE + ". Flat space is empty outside the light cone of the origin, "
        "$r = ct$, and the universe grows inside it."],
    "fireball_proper": [
        r"The plane of $T$ and $\chi$" + PLANE + r", up to the end of creation ($HT_1 = 2\ln 3$). Null rays "
        r"conserve $\eta \pm \chi$ with $\eta = \ln\tanh(HT/2)$, which diverges at $T = 0$: a ray from any "
        r"$\chi$ reaches the world line $\chi = 0$."],
    "expansion_proper": [
        r"The plane of $T$ and $\chi$" + PLANE + r", from the end of creation ($HT_1 = 2\ln 3$). Null rays "
        r"conserve $\eta \pm \chi$ with $\cosh(\eta - \eta_b) = (T - T_b)/b$."],
    "fireball_conformal": [
        r"The plane of $\eta$ and $\chi$" + PLANE + r", up to the end of creation ($\eta_1 = \ln 0.8$). "
        + STRAIGHT + r", and the light cone of the origin lies at $\eta \to -\infty$."],
    "expansion_conformal": [
        r"The plane of $\eta$ and $\chi$" + PLANE + r", from the end of creation ($\eta_1 = \ln 0.8$). "
        + STRAIGHT + "."],
}


def register_null_rays(nr):
    """One radial plane for each chart, with the null integrals the rays are checked against."""
    for system in CHARTS:
        coords = SYSTEMS[system]["coords"]
        box, xlabel, ylabel, where = PLANES[system]
        extra = {"where": where} if where else {}
        nr.DIAGRAMS.append(nr.Diagram(
            METRIC, system, "radial", f"${coords[0]}$ and ${coords[1]}$", (coords[0], coords[1]), box,
            xlabel, ylabel, drawn(system), FIXED, tau=coords[0].lstrip("\\"), **extra))
        nr.CAPTIONS[(METRIC, system, "radial")] = RAY_CAPTIONS[system]
        if system.endswith("_flat"):
            nr.CLOSED_FORMS[(METRIC, system, "radial")] = (lambda t, r: t + r, lambda t, r: t - r, None)
        else:
            nr.CLOSED_FORMS[(METRIC, system, "radial")] = (
                lambda t, r, s=system: conformal_time(s, t) + r,
                lambda t, r, s=system: conformal_time(s, t) - r, None)


# The embedded moment: the equatorial plane at one time of the flat background, later than the
# end of creation on the central world line and no later than tau_infinity, the latest such
# time whose plane a surface of revolution in flat space carries; the notes derive that bound.
MOMENT = 0.9*END
MOMENT_LABEL = r"$t = 0.9\,\tau_\infty$"
EMBEDDING_CAPTIONS = {
    "background": [
        r"The equatorial plane at the background time $t = 0.9\,\tau_\infty$ "
        r"($\tau_1 = 0.8\,\tau_\infty$): the free expansion at the centre, the fireball around it, and "
        r"flat space beyond the light cone. At $\delta = 0$ the surface has a crease on the light cone, "
        "where the slope of the conformal factor jumps."],
}


def embedding_diagram(ck, src):
    """The plane theta = pi/2 at the background time t = 0.9 tau_infinity, in units of c/H.

    On it the induced metric is Omega^2 (dr^2 + r^2 dvarphi^2) with Omega the conformal factor's
    root at tau^2 = t^2 - r^2/c^2: the free expansion from the centre to the circle tau = tau_1,
    the fireball from there to the light cone r = ct, and flat space beyond. Each circle has the
    radius rho = Omega r, and the surface climbs at dz/dr = sqrt(g_rr - (drho/dr)^2), all read
    from the published metric of the three flat charts. The rim is set at z = 0.
    """
    import math

    import numpy as np

    import embedding as em

    t, reach, size = MOMENT, 3.2, 6.4
    join = math.sqrt(t*t - (RATIO*END)**2)
    fixed = {"t": repr(t), "theta": "pi/2"}
    slices = {system: em.Slice(src, METRIC, system, "r", r"\varphi", fixed, drawn(system))
              for system in ("expansion_flat", "fireball_flat", "outside_flat")}
    core = em.Piece("expansion", "star", slices["expansion_flat"], 0.0, join, 0.0, 1,
                    (("axis", "the central world line"), ("join", r"the end of creation, $\tau = \tau_1$")),
                    [(0.4, "r", None), (join, "surface", None)], size)
    fireball = em.Piece("fireball", "sheet", slices["fireball_flat"], join, t, core.z[-1], 1,
                        (("join", r"the end of creation, $\tau = \tau_1$"),
                         ("crease", "the light cone of the origin, $r = ct$")),
                        [(1.2, "r", None), (1.5, "r", None), (t, "chartedge", None)], size)
    outside = em.Piece("outside", "sheet2", slices["outside_flat"], t, reach, fireball.z[-1], 1,
                       (("crease", "the light cone of the origin, $r = ct$"), ("edge", "flat space runs on")),
                       [(2.4, "r2", None), (reach, "r2", None)], size)
    shift = -outside.z[0]
    for piece in (core, fireball, outside):
        piece.z = piece.z + shift
    where = f"{METRIC}, t = {t:g}"
    for name, piece in (("the free expansion", core), ("the fireball", fireball), ("flat space", outside)):
        ck.isometry(f"{where}, {name}", piece)
    ck.join(f"{where}, the free expansion meets the fireball", core, join, fireball, join)
    # On the cone the factor is one on both sides and its slope is not: d(rho)/dr is
    # 1 - 2t^2/tau_infinity^2 just inside and 1 outside, so the two meet in one circle at a crease.
    inside, flat = np.array(fireball.at(t)), np.array(outside.at(t))
    ck.add(f"{where}, the fireball meets flat space: one circle", float(np.max(np.abs(inside - flat))), em.JOIN)
    cosine = 2*(t/END)**2 - 1
    found, level = fireball.inward(t), outside.inward(t)
    ck.add(f"{where}, the crease on the light cone: the fireball leaves it at cos = 2t^2 - 1, and the plane level",
           float(max(abs(found[0] - cosine), abs(found[1] + math.sqrt(1 - cosine**2)),
                     abs(level[0] - 1), abs(level[1]))), em.JOIN)
    k = joined(END, RATIO*END)
    ck.radius(f"{where}, rho = Omega r in the free expansion", core,
              lambda r: float(k["Omega_f"])*(1 - float(k["tau_b"])**2/(t*t - r*r))*r, size)
    ck.radius(f"{where}, rho = Omega r in the fireball", fireball, lambda r: r/(1 - (t*t - r*r)/END**2), size)
    ck.plane(f"{where}, flat space outside the cone", slices["outside_flat"], np.linspace(t, reach, 40))
    surface = em.Surface([core, fireball, outside])
    fig = em.figure_of([surface], {"star": "star", "sheet": "cover"}, size, meridians=12)
    fig.legend("fill", "star", "free expansion")
    fig.legend("fill", "cover", "fireball")
    fig.legend("line", "surface", "end of creation")
    fig.legend("line", "chartedge", "light cone")
    fig.legend("line", "r", "$r$ constant")
    fig.legend("line", "r2", "$r$ constant, flat space")
    fig.legend("line", "meridian", r"$\varphi$ constant")
    return [em.view("background", "At one background time", "$c/H$", [surface], fig.done(),
                    settings="$c/H = 1$, the unit of every length.")]


def moment_line(system, m, window):
    """The embedded moment on a chart's radial plane, as rows of (time, radius): the line
    t = 0.9 tau_infinity on the flat charts, and on the others the curve tau cosh chi = t, as far
    as the chart's region and the window reach."""
    import math

    import numpy as np

    if system.endswith("_flat"):
        return [[(MOMENT, lo), (MOMENT, hi)] for lo, hi in [m.reach(system, "r")]]
    turn = math.acosh(MOMENT/(RATIO*END))
    chi = np.linspace(0.0, turn, 81) if system.startswith("expansion") else np.linspace(turn, window, 161)
    eta = np.log(MOMENT/END/np.cosh(chi))
    if system.endswith("_conformal"):
        return [np.column_stack([eta, chi])]
    values = drawn("expansion_proper")
    proper = (2*np.arctanh(np.exp(eta)) if system == "fireball_proper" else
              values["T_b"] + values["b"]*np.cosh(eta - 2*math.log(RATIO)))
    return [np.column_stack([proper, chi])]


def register_slices(sl):
    """The one embedded moment on every radial plane it crosses."""
    for system in CHARTS:
        window = PLANES[system][0][1]
        sl.FLAT[(METRIC, system, "radial")] = (
            lambda system=system, window=window: sl.one(
                METRIC, lambda m: moment_line(system, m, window), label=MOMENT_LABEL, view_id="background"))
    sl.FLAT_METRICS.add(METRIC)


CONFORMAL_CAPTION = (
    "The whole spacetime on the compactified flat background, each point in the diagram a 2-sphere; "
    r"$p, q = \arctan(H(t \mp r/c))$. The hyperboloid on which creation ends "
    r"($\tau_1 = 0.8\,\tau_\infty$) lies wholly in the future of the origin [gott1983].")
CONFORMAL_CAPTIONS = {system: [CONFORMAL_CAPTION] for system in CHARTS}


def background(system, first, radius):
    """A chart's point as (t, r) of the flat background, in units of 1/H and c/H."""
    import numpy as np
    if system.endswith("_flat"):
        return first, radius
    tau = END*np.exp(conformal_time(system, first))
    return tau*np.cosh(radius), tau*np.sinh(radius)


def conformal_diagram(ck, src):
    """The spacetime is conformal to flat space, so its diagram is flat space's triangle with
    p, q = arctan(t -+ r): flat space below the light cone p = 0 of the origin, the fireball
    between that cone and the hyperboloid tau = tau_1, and the free expansion above it, whose
    conformal factor tends to a constant, so that its future infinity is flat space's. One view
    for each chart tints the region the chart covers and draws its lines of constant time; every
    chart's map onto the triangle is checked to send its null rays to p and q constant."""
    import numpy as np

    import conformal as cd
    import slices

    pi = np.pi

    def pq(system, first, radius):
        t, r = background(system, first, radius)
        return np.arctan(t - r), np.arctan(t + r)

    def xt(p, q):
        return np.column_stack([q - p, p + q]).tolist()

    chi = np.concatenate([np.linspace(0, 6, 240), np.geomspace(6, 40, 60)[1:]])
    hyperbola = xt(np.arctan(RATIO*END*np.exp(-chi)), np.arctan(RATIO*END*np.exp(chi)))
    origin, corner, top, bottom, right = [0, 0], [pi/2, pi/2], [0, pi], [0, -pi], [pi, 0]
    regions = {
        "outside": [bottom, right, corner, origin],
        "fireball": [origin, corner] + hyperbola[::-1],
        "expansion": hyperbola + [corner, top],
    }
    # Sample points of each chart, inside its region, for the check of its map.
    tau_in, tau_out, angle = ck.uniform(0.1, 1.58), ck.uniform(1.62, 8.0), ck.uniform(0.02, 2.5)
    span = ck.uniform(-6, 6)
    samples = {
        "fireball_flat": (tau_in*np.cosh(angle), tau_in*np.sinh(angle)),
        "expansion_flat": (tau_out*np.cosh(angle), tau_out*np.sinh(angle)),
        "outside_flat": (span, np.abs(span) + ck.uniform(0.05, 6)),
        "fireball_proper": (ck.uniform(0.1, 2.19), angle),
        "expansion_proper": (ck.uniform(2.2, 8.0), angle),
        "fireball_conformal": (ck.uniform(-3, -0.23), angle),
        "expansion_conformal": (ck.uniform(-0.22, 2.5), angle),
    }
    # The lines of constant time each view draws, as values of its chart's time.
    times = {
        "fireball_flat": (0.6, 1.2, 1.8, 3.0, 6.0), "expansion_flat": (2.0, 3.0, 6.0, 12.0),
        "outside_flat": (-6.0, -2.0, 0.0, 2.0, 6.0),
        "fireball_proper": (0.5, 1.0, 1.5, 2.0), "expansion_proper": (3.0, 5.0, 10.0, 30.0),
        "fireball_conformal": (-2.0, -1.0, -0.5), "expansion_conformal": (0.0, 0.5, 1.0, 2.0),
    }
    moment = slices.moments(METRIC, "background", label=MOMENT_LABEL)[0]
    reach = np.linspace(0.0, moment.reach("outside_flat", "r")[1], 121)
    views = []
    for system in CHARTS:
        coords = SYSTEMS[system]["coords"]
        plane = cd.Plane(src, METRIC, system, (coords[0], coords[1]), FIXED, drawn(system))
        ck.chart(f"{METRIC} {system}", plane, lambda a, b, s=system: pq(s, a, b), *samples[system],
                 lambda a, b: (1, 0))
        v = cd.View(system, NAMES[system], [-0.3, pi + 0.3, -pi - 0.3, pi + 0.3], system)
        v.fill("region", [bottom, right, top])
        v.fill("cover", regions[system.partition("_")[0]])
        v.line("scri", [[bottom, right], [right, top]])
        v.line("centre", [[bottom, top]])
        v.line("null", [[origin, corner]])
        v.line("surface", [hyperbola])
        for value in times[system]:
            if system == "outside_flat":
                r = np.abs(value) + cd.spread(0, np.inf, 300, 9) if value > 0 else cd.spread(0, np.inf, 300, 9)
                v.curve("t", *pq(system, np.full_like(r, value), r))
            elif system.endswith("_flat"):
                # t constant inside the cone, across the part of it the chart covers
                r = np.linspace(0.0, value, 400)[:-1]
                tau = np.sqrt(value*value - r*r)
                keep = (tau <= RATIO*END) if system == "fireball_flat" else (tau >= RATIO*END)
                if keep.any():
                    v.curve("t", *pq(system, np.full(keep.sum(), value), r[keep]))
            else:
                v.curve("t", *pq(system, np.full_like(chi, value), chi))
        v.slice(moment, [pq("outside_flat", np.full_like(reach, MOMENT), reach)])
        v.legend("cover", "chart domain")
        v.legend("t", f"${coords[0]}$ constant")
        v.legend("null", "light cone of the origin")
        v.legend("surface", "end of creation")
        v.legend("centre", "central world line")
        v.legend("scri", r"flat space's $\mathscr{I}^\pm$")
        v.set(settings="$H = 1$, $c = 1$.")
        views.append(v)
    return views


def curvature_diagram(ck, src):
    """The two stored scalars along the proper time of the comoving world lines, through the end
    of creation: the fireball's from its chart up to T_1 and the free expansion's from there on.
    Both are homogeneous on the hyperboloids, so one line each says everything; the Kretschmann
    scalar is continuous at T_1, where it is 24 H^4/c^4 on both sides, and the Ricci scalar
    drops there from 12 H^2/c^2 to zero. Nothing is computed again."""
    import curvature as cv

    T_1 = drawn("fireball_proper")["T_1"]
    x_axis = cv.axis("$HT$", 0, 6)
    stages = (("fireball_proper", (0.0, T_1)), ("expansion_proper", (T_1, 6.0)))
    panels, curves = [], []
    for field, kind, power, lo, hi in (("kretschmann", "K", 4, 0, 30), ("ricci_scalar", "R", 2, -3, 15)):
        scale = {"kind": "linear", "from": lo, "to": hi}
        panels.append({"id": kind, "name": f"${kind}\\,c^{power}/H^{power}$", "aspect": 2, "scale": scale,
                       "ticks": cv.value_ticks(scale)})
        for system, within in stages:
            params = drawn(system)
            stored = cv.Published(src, METRIC, system, field, ["T"], list(params))
            curves += cv.line_curves(ck, f"{METRIC} {system} {kind}", kind, kind, stored, "T", x_axis, scale,
                                     cv.PLOT_W/2, params, within=within)
    cv.shape_checks(ck, METRIC, [2, 2])
    return [{"id": "proper_time", "label": "Along the proper time", "kind": "line", "x": x_axis,
             "panels": panels, "curves": curves, "marks": [],
             "legend": [["line", "K", "$K$"], ["line", "R", "$R$"]],
             "caption": [r"The Kretschmann scalar $K$ and the Ricci scalar $R$ of the fireball and the free "
                         r"expansion, along the proper time of the comoving world lines ($HT_1 = 2\ln 3$). "
                         "At the end of creation $R$ drops from $12H^2/c^2$ to zero. "
                         "$K$ is continuous there at $24H^4/c^4$ and then falls as the inverse eighth power "
                         "of the scale factor."],
             "ricci_sentences": [[0, 1]], "per_scalar": cv.scalar_clips("Linear scale.", "Linear scale."),
             "settings": "$H = 1$, $c = 1$.", "clip": "Linear scales."}]
