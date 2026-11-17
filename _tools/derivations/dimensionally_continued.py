"""The dimensionally continued black holes of Banados, Teitelboim and Zanelli, in five charts.

The static charts are their Eqs. (22) and (23) at n = 3 (Phys. Rev. D 49, 975),

    ds^2 = -f c^2dt^2 + dr^2/f + r^2 dOmega^2,
    f = 1 - sqrt(M + 1) + r^2/l^2   in five dimensions, the odd family,
    f = 1 - sqrt(2M/r) + r^2/l^2    in six dimensions, the even family,

with the positive square root, the branch Muniain and Piriz kept (Phys. Rev. D 53, 816), on which
a positive mass has a horizon. The advanced time charts are Nozawa and Maeda's Eq. (3.7) with a
constant mass (Class. Quantum Grav. 23, 1779), -f dv^2 + 2 dv dr + r^2 dOmega^2 with
dv = c dt + dr/f, and the conformal chart of five dimensions is Shiraishi and Maki's Eqs. (7) and
(8) (Class. Quantum Grav. 11, 1687), r = r_+ sec(rho) outside the horizon
r_+ = l sqrt(sqrt(M + 1) - 1),

    ds^2 = l^2 sec^2(rho) (-(r_+^2/l^4) sin^2(rho) c^2dt^2 + d rho^2 + (r_+^2/l^2) dOmega^2).

dimensionally_continued.md beside this file records each chart's source.

Every value of the five dimensional charts is a rational function of r, l and sqrt(M + 1), and of
the six dimensional ones of r, l and w = sqrt(2M/r), a name the charts of six dimensions declare;
pretty writes each root as a placeholder before factoring, and the printer writes it back.

Each chart is checked at points against the field equations of the theory, a product of two
factors of the anti-de Sitter curvature, which in tensors is

    H_ab + 2 (D-3)(D-4) G_ab / l^2 - (D-1)(D-2)(D-3)(D-4) g_ab / (2 l^4) = 0,

with G the Einstein tensor and H the Gauss-Bonnet tensor; each static chart against the source's
metric function, horizon and surface gravity; each other chart against the static one pulled back.
"""
import random
import re

import sympy as sp

METRIC = "dimensionally_continued"
SYSTEMS = ("spherical", "eddington_finkelstein_ingoing", "conformal",
           "spherical_six", "eddington_finkelstein_ingoing_six")
SIX = ("spherical_six", "eddington_finkelstein_ingoing_six")
NAMES = {
    "spherical": "Hyperspherical, Five Dimensions",
    "eddington_finkelstein_ingoing": "Ingoing Eddington-Finkelstein, Five Dimensions",
    "conformal": "Conformal, Five Dimensions",
    "spherical_six": "Hyperspherical, Six Dimensions",
    "eddington_finkelstein_ingoing_six": "Ingoing Eddington-Finkelstein, Six Dimensions",
}
ROOT_FIVE = "\\sqrt{M + 1}"
ROOT_SIX = "\\sqrt{\\dfrac{2M}{r}}"
W_SIX = "w = \\sqrt{2M/r}"
F_FIVE = "1 - " + ROOT_FIVE + " + \\dfrac{r^2}{\\ell^2}"
F_SIX = "1 - " + ROOT_SIX + " + \\dfrac{r^2}{\\ell^2}"
F_SIX_NAMED = "1 - w + \\dfrac{r^2}{\\ell^2}"
# r_+^2/l^2 in five dimensions, Shiraishi and Maki's M-tilde.
GAP_FIVE = ROOT_FIVE + " - 1"
SPHERE_FIVE = "d\\psi^2 + \\sin^2\\psi\\,d\\theta^2 + \\sin^2\\psi\\sin^2\\theta\\,d\\phi^2"
SPHERE_SIX = ("d\\chi^2 + \\sin^2\\chi\\,d\\psi^2 + \\sin^2\\chi\\sin^2\\psi\\,d\\theta^2"
              " + \\sin^2\\chi\\sin^2\\psi\\sin^2\\theta\\,d\\phi^2")
ANGLES = {False: ["\\psi", "\\theta", "\\phi"], True: ["\\chi", "\\psi", "\\theta", "\\phi"]}
HORIZON_FIVE = "r = r_+ = \\ell\\sqrt{" + GAP_FIVE + "} \\;\\text{(the horizon)}"
HORIZON_SIX = "r = r_+ \\;\\text{(the horizon)}"


def coords(system):
    six = system in SIX
    first = "v" if system.startswith("eddington") else "t"
    second = "\\rho" if system == "conformal" else "r"
    return [first, second] + ANGLES[six]


def parameters(system):
    return ["M", "\\ell"] + ([W_SIX] if system in SIX else [])


def line_element(system, c="c"):
    """The line element as the page prints it, with c; c = "" gives the chart's own, in x^0 = ct."""
    six = system in SIX
    f = F_SIX if six else F_FIVE
    sphere = SPHERE_SIX if six else SPHERE_FIVE
    c2 = c + "^2" if c else ""
    if system.startswith("spherical"):
        return (f"ds^2 = -\\left({f}\\right){c2}dt^2 + \\dfrac{{dr^2}}{{{f}}}"
                f" + r^2\\left({sphere}\\right)")
    if system.startswith("eddington"):
        return f"ds^2 = -\\left({f}\\right)dv^2 + 2\\,dv\\,dr + r^2\\left({sphere}\\right)"
    return (f"ds^2 = -\\left({GAP_FIVE}\\right)\\tan^2\\rho\\,{c2}dt^2 + \\ell^2\\sec^2\\rho\\,d\\rho^2"
            f" + \\ell^2\\left({GAP_FIVE}\\right)\\sec^2\\rho\\left({sphere}\\right)")


def domains(system):
    six = system in SIX
    angles = [a + " \\in [0, \\pi]" for a in ANGLES[six][:-1]] + ["\\phi \\in [0, 2\\pi)"]
    if system != "conformal":
        return ([coords(system)[0] + " \\in (-\\infty, \\infty)", "r \\in (0, \\infty)"] + angles
                + [HORIZON_SIX if six else HORIZON_FIVE])
    return ["t \\in (-\\infty, \\infty)", "\\rho \\in (0, \\pi/2)"] + angles


def register_dimensions(vm):
    """t is a time and v a length; in five dimensions M is a number, in six a length, and the
    root w = sqrt(2M/r) a number; rho is an angle."""
    for system in SYSTEMS:
        six = system in SIX
        names = coords(system)
        table = {names[0]: "L" if names[0] == "v" else "T", names[1]: "1" if names[1] == "\\rho" else "L"}
        table.update({a: "1" for a in ANGLES[six]})
        table.update({"M": "L" if six else "1", "\\ell": "L"})
        if six:
            table["w"] = "1"
        vm.DIMENSIONS[(METRIC, system)] = table


def spec(system):
    import verify_metrics as vm
    six = system in SIX
    names, given = coords(system), parameters(system)
    probe = vm.Reader(names, given, ())
    x = probe.symbol[names[1]]
    M, ell = probe.parameters["M"], probe.parameters["ell"]
    root = sp.Symbol("DCroot", positive=True)
    radius, mass = sp.Symbol("DCr", positive=True), sp.Symbol("DCm", positive=True)

    lapse = sp.Symbol("DClapse")
    # l^2 times the metric function of the static chart, which vanishes on the horizon.
    whole = ell ** 2 * (1 - root) + x ** 2

    def pretty(value):
        value = sp.sympify(value)
        if value == 0:
            return value
        if six:
            # M = w^2 r/2, taken with positive names so that every root of M and of r closes.
            value = sp.together(sp.powsimp(sp.powdenest(
                value.subs({M: root ** 2 * radius / 2, x: radius}), force=True), force=True)).subs(radius, x)
        else:
            value = sp.together(sp.powdenest(value.subs(M, mass ** 2 - 1), force=True)).subs(mass, root)
        if value.has(M) or any(not p.exp.is_Integer for p in value.atoms(sp.Pow)):
            raise AssertionError(f"{METRIC}/{system}: a root is left in {value}")
        if system == "conformal":
            return conformal(value)
        # The factor that vanishes on the horizon is written whole, in one order.
        out = sp.Integer(1)
        for factor in sp.Mul.make_args(sp.factor(value)):
            base, k = (factor.base, factor.exp) if factor.is_Pow else (factor, sp.Integer(1))
            if sp.expand(base - whole) == 0:
                out *= lapse ** k
            elif sp.expand(base + whole) == 0:
                out *= (-1) ** k * lapse ** k
            else:
                out *= factor
        return out

    gap = sp.Symbol("DCgap")
    sine, cosine = sp.Symbol("DCs"), sp.Symbol("DCc")

    def conformal(value):
        """A value of the conformal chart, a rational function of sqrt(M + 1), sin(rho) and
        cos(rho): the pair (sin - 1)(sin + 1) is written -cos^2, so that every power of the
        conformal factor stands as a power of the cosine, and sqrt(M + 1) - 1, which is
        r_+^2/l^2, is written whole and in that order."""
        value = sp.factor(value.subs({sp.sin(x): sine, sp.cos(x): cosine}))
        out, powers = sp.Integer(1), {}
        for factor in sp.Mul.make_args(value):
            base, k = (factor.base, factor.exp) if factor.is_Pow else (factor, sp.Integer(1))
            if base.is_Number:
                out *= factor
            elif sp.expand(base + sine - 1) == 0 or sp.expand(base + root - 1) == 0:
                out *= (-1) ** k
                powers[-base] = powers.get(-base, 0) + k
            else:
                powers[base] = powers.get(base, 0) + k
        below, above = powers.pop(sine - 1, 0), powers.pop(sine + 1, 0)
        k = min(below, above) if below > 0 and above > 0 else max(below, above) if below < 0 and above < 0 else 0
        out *= (-1) ** k
        for base, e in ((sine - 1, below - k), (sine + 1, above - k), (cosine, 2 * k)):
            if e:
                powers[base] = powers.get(base, 0) + e
        if root - 1 in powers:
            powers[gap] = powers.pop(root - 1)
        for base, e in powers.items():
            out *= base ** e
        return out.subs({sine: sp.sin(x), cosine: sp.cos(x)})

    f = F_SIX_NAMED if six else F_FIVE
    bracket = "\\left(" + f + "\\right)"
    if system.startswith("spherical"):
        by_hand = {"metric_components": {("t", "t"): "-" + bracket, ("r", "r"): bracket + "^{-1}"},
                   "inverse_metric_components": {("t", "t"): "-" + bracket + "^{-1}", ("r", "r"): f}}
    elif system.startswith("eddington"):
        by_hand = {"metric_components": {("v", "v"): "-" + bracket},
                   "inverse_metric_components": {("r", "r"): f}}
    else:
        factor, sphere = "\\left(" + GAP_FIVE + "\\right)", "\\ell^2\\left(" + GAP_FIVE + "\\right)\\sec^2\\rho"
        by_hand = {"metric_components": {
            ("t", "t"): "-" + factor + "\\tan^2\\rho", ("\\rho", "\\rho"): "\\ell^2\\sec^2\\rho",
            ("\\psi", "\\psi"): sphere, ("\\theta", "\\theta"): sphere + "\\sin^2\\psi",
            ("\\phi", "\\phi"): sphere + "\\sin^2\\psi\\sin^2\\theta"}}
    scalars = {
        "ricci_scalar": ("-\\dfrac{20}{\\ell^2} + \\dfrac{6" + ROOT_FIVE + "\\cos^2\\rho}"
                         "{\\ell^2\\left(" + GAP_FIVE + "\\right)}"),
        "kretschmann": ("\\dfrac{40}{\\ell^4} - \\dfrac{24" + ROOT_FIVE + "\\cos^2\\rho}"
                        "{\\ell^4\\left(" + GAP_FIVE + "\\right)} + \\dfrac{12\\left(M + 1\\right)\\cos^4\\rho}"
                        "{\\ell^4\\left(" + GAP_FIVE + "\\right)^2}"),
    } if system == "conformal" else {
        "ricci_scalar": ("-\\dfrac{30}{\\ell^2} + \\dfrac{35w}{4r^2}" if six
                         else "-\\dfrac{20}{\\ell^2} + \\dfrac{6" + ROOT_FIVE + "}{r^2}"),
        "kretschmann": ("\\dfrac{60}{\\ell^4} - \\dfrac{35w}{\\ell^2 r^2} + \\dfrac{425w^2}{16r^4}" if six
                        else "\\dfrac{40}{\\ell^4} - \\dfrac{24" + ROOT_FIVE + "}{\\ell^2 r^2}"
                             " + \\dfrac{12\\left(M + 1\\right)}{r^4}"),
    }
    return {
        "metric_id": METRIC,
        "system": {"id": system, "name": NAMES[system], "coords": names, "domains": domains(system),
                   "parameters": given, "line_element": line_element(system)},
        "chart_line_element": line_element(system, c=""),
        "printer": {"rising": [root, sp.sin(x)] if system == "conformal" else [x, root],
                    "lead": [x, root, ell], "factors": [ell, x, root], "flip": False,
                    "overrides": {root: "w" if six else ROOT_FIVE},
                    "named": {gap: GAP_FIVE} if system == "conformal" else {lapse: ("\\ell^2\\left(1 - w\\right) + r^2" if six
                                      else "r^2 - \\ell^2\\left(" + GAP_FIVE + "\\right)")}},
        "pretty": pretty,
        "components": by_hand,
        "check": lambda chart: check(chart, system),
        "rewrite": [("\\," + ROOT_FIVE, ROOT_FIVE)],
        **({"after": rho_first} if system == "conformal" else {}),
        **scalars,
    }


SPHERE_THEN_RHO = re.compile(r"((?:\\(?:sin|cos|cot|csc)(?:\^\d+)?\\psi)+)"
                             r"((?:\\(?:sin|cos|cot|csc)(?:\^\d+)?\\rho)+)")


def rho_first(math, chart):
    """The conformal chart's values with the functions of rho written before those of the sphere's
    angles, as its line element writes them; each changed text is read back against the one it
    replaces."""
    import verify_metrics as vm
    if isinstance(math, dict):
        return {k: rho_first(v, chart) for k, v in math.items()}
    if isinstance(math, list):
        return [rho_first(v, chart) for v in math]
    if not isinstance(math, str):
        return math
    text = SPHERE_THEN_RHO.sub(r"\2\1", math)
    if text != math and "\\dot" not in math:
        for side, original in zip(text.split("="), math.split("="), strict=True):
            if side != original and vm.norm(chart.reader(side) - chart.reader(original)) != 0:
                raise AssertionError(f"{METRIC}/conformal: {side!r} does not read back as {original!r}")
    return text


def register(module):
    module.CHARTS[METRIC] = [lambda system=system: spec(system) for system in SYSTEMS]


# -- Checks at points, independent of the printed tensors ---------------------------------------

def gauss_bonnet_residual(geo, n, ell, at, digits=40):
    """The largest entry of H_ab + 2 (n-3)(n-4) G_ab/l^2 - (n-1)(n-2)(n-3)(n-4) g_ab/(2 l^4) at a
    point, with H_ab = 2 (R R_ab - 2 R_ac R^c_b - 2 R^cd R_acbd + R_a^cde R_bcde) - g_ab L/2 and
    L = R^2 - 4 R_cd R^cd + R_cdef R^cdef the Gauss-Bonnet density, relative to 1/l^4."""
    import verify_metrics as vm

    def number(e):
        return sp.sympify(e).xreplace(at).evalf(digits)
    span = range(n)
    riemann, ricci = geo.riemann_llll(), geo.ricci_ll()
    low = [[number(geo.g[a, b]) for b in span] for a in span]
    up = [[number(geo.ginv[a, b]) for b in span] for a in span]
    R4 = {tuple(i): number(vm._at(riemann, i)) for i in vm._indices(n, 4)}
    R2 = [[number(ricci[a][b]) for b in span] for a in span]
    R = sum(up[a][b] * R2[a][b] for a in span for b in span)
    R2uu = [[sum(up[a][c] * up[b][d] * R2[c][d] for c in span for d in span) for b in span] for a in span]
    R2mixed = [[sum(up[c][d] * R2[a][d] for d in span) for c in span] for a in span]
    # R_a^{cde}, one index raised at a time, and R^{acde} for the density.
    step = {(a, c, d, e): sum(up[c][p] * R4[a, p, d, e] for p in span)
            for a in span for c in span for d in span for e in span}
    step = {(a, c, d, e): sum(up[d][p] * step[a, c, p, e] for p in span)
            for a in span for c in span for d in span for e in span}
    raised = {(a, c, d, e): sum(up[e][p] * step[a, c, d, p] for p in span)
              for a in span for c in span for d in span for e in span}
    full = {(a, c, d, e): sum(up[a][p] * raised[p, c, d, e] for p in span)
            for a in span for c in span for d in span for e in span}
    density = (R ** 2 - 4 * sum(R2uu[a][b] * R2[a][b] for a in span for b in span)
               + sum(full[i] * R4[i] for i in R4))
    lam = 1 / number(ell) ** 2
    worst = sp.Float(0)
    for a in span:
        for b in span:
            H = 2 * (R * R2[a][b] - 2 * sum(R2mixed[a][c] * R2[c][b] for c in span)
                     - 2 * sum(R2uu[c][d] * R4[a, c, b, d] for c in span for d in span)
                     + sum(raised[a, c, d, e] * R4[b, c, d, e] for c in span for d in span for e in span)
                     ) - low[a][b] * density / 2
            einstein = R2[a][b] - low[a][b] * R / 2
            residual = (H + 2 * (n - 3) * (n - 4) * lam * einstein
                        - (n - 1) * (n - 2) * (n - 3) * (n - 4) * lam ** 2 * low[a][b] / 2)
            # A component of the tensor scales with the coordinates' own lengths, so it is held
            # against the size of the terms that cancel in it.
            scale = abs(H) + abs(lam * einstein) + abs(lam ** 2 * low[a][b]) + lam ** 2
            worst = max(worst, abs(residual) / scale)
    return worst


def source_function(system, r, M, ell):
    """The metric function of the source, with the positive square root."""
    return 1 - sp.sqrt(2 * M / r if system in SIX else M + 1) + r ** 2 / ell ** 2


def check(chart, system):
    import verify_metrics as vm
    six = system in SIX
    n = 6 if six else 5
    geo, g = chart.geo, chart.geo.g
    x = chart.symbols[1]
    M, ell = chart.reader.parameters["M"], chart.reader.parameters["ell"]
    rng = random.Random(1994)
    points = []
    for _ in range(3):
        at = {s: sp.Rational(rng.randint(400, 1300), 1000) for s in chart.symbols}
        at.update({M: sp.Rational(rng.randint(500, 3000), 1000), ell: sp.Rational(rng.randint(700, 2500), 1000)})
        if system != "conformal":
            # Outside the horizon on one draw, inside it on the next where the chart reaches there.
            at[x] = sp.Rational(rng.randint(2600, 6000), 1000)
        points.append(at)
    if system.startswith("eddington"):
        points[1][x] = sp.Rational(rng.randint(150, 400), 1000)
    for at in points:
        worst = gauss_bonnet_residual(geo, n, ell, at)
        if worst > sp.Float(10) ** -30:
            raise AssertionError(f"{METRIC}/{system}: the field equations fail at {at}, by {worst}")
    f = source_function(system, x, M, ell)
    if system.startswith("spherical"):
        for at in points:
            if abs((g[0, 0] + f).xreplace(at).evalf(40)) > sp.Float(10) ** -30 or \
                    abs((g[1, 1] - 1 / f).xreplace(at).evalf(40)) > sp.Float(10) ** -30:
                raise AssertionError(f"{METRIC}/{system}: the chart is not -f c^2dt^2 + dr^2/f")
        horizon(system, f, x, M, ell)
        return
    # Every other chart is the static one of its dimension pulled back.
    static = "spherical_six" if six else "spherical"
    import chart_printer as cp
    source = cp.Chart(coords(static), parameters(static), line_element(static, c=""))
    r = source.symbols[1]
    same = dict(zip(source.symbols[2:], chart.symbols[2:]))
    same.update({source.reader.parameters[k]: chart.reader.parameters[k] for k in ("M", "ell")})
    same[source.symbols[0]] = chart.symbols[0]
    J = sp.eye(n)
    if system.startswith("eddington"):
        # c dt = dv - dr/f.
        same[r] = x
        J[0, 1] = -1 / f
    else:
        # r = r_+ sec(rho), with r_+ = l sqrt(sqrt(M + 1) - 1).
        radius = ell * sp.sqrt(sp.sqrt(M + 1) - 1) / sp.cos(x)
        same[r] = radius
        J[1, 1] = sp.diff(radius, x)
    old = source.geo.g.subs(same, simultaneous=True)
    pulled = J.T * old * J
    for at in points:
        for i in range(n):
            for j in range(i, n):
                if abs((pulled[i, j] - g[i, j]).xreplace(at).evalf(40)) > sp.Float(10) ** -30:
                    raise AssertionError(f"{METRIC}/{system}: the static chart pulled back misses slot "
                                         f"{chart.coords_tex[i]}{chart.coords_tex[j]} at {at}")


def horizon(system, f, r, M, ell):
    """f vanishes on the horizon the source gives, with the surface gravity f'(r_+)/2 of the
    source's temperature T = hbar c f'(r_+)/(4 pi k_B): r_+/l^2 in five dimensions, where
    r_+^2 = l^2 (sqrt(M + 1) - 1), and (1 + 5 r_+^2/l^2)/(4 r_+) in six, where
    r_+ (1 + r_+^2/l^2)^2 = 2M. Anti-de Sitter space is M = -1 in five dimensions and M = 0 in six."""
    rp, l = sp.symbols("r_p l", positive=True)
    if system in SIX:
        mass = rp * (1 + rp ** 2 / l ** 2) ** 2 / 2
        gravity = (1 + 5 * rp ** 2 / l ** 2) / (4 * rp)
        empty = 0
    else:
        mass = (1 + rp ** 2 / l ** 2) ** 2 - 1
        gravity = rp / l ** 2
        empty = -1
    on = {M: mass, ell: l}
    if sp.simplify(f.subs(on).subs(r, rp)) != 0:
        raise AssertionError(f"{METRIC}/{system}: f does not vanish on the source's horizon")
    if sp.simplify(sp.diff(f, r).subs(on).subs(r, rp) / 2 - gravity) != 0:
        raise AssertionError(f"{METRIC}/{system}: the surface gravity is not the source's")
    if sp.simplify(f.subs(M, empty) - (1 + r ** 2 / ell ** 2)) != 0:
        raise AssertionError(f"{METRIC}/{system}: the chart is not anti-de Sitter space at M = {empty}")
