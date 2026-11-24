"""The wall of Dvali, Gabadadze and Porrati, Phys. Lett. B 485, 208 (2000), hep-th/0005016: a
wall with three dimensions of space in flat space of five dimensions, with a curvature term of
its own in its action, and its cosmology, Deffayet, Phys. Lett. B 502, 199 (2001), hep-th/0010186.

Every chart is in the signature (-,+,+,+,+). The bulk solves R_mu nu = 0 off the wall, and the
wall, with the mirror symmetry across it, the junction condition of the action of Dvali,
Gabadadze and Porrati,
    K_ab - h_ab K = r_c G_ab - 4 pi G_5 T_ab / c^4,     r_c = G_5 / 2G,
with K_ab the extrinsic curvature of the wall on either side, its normal pointing away from it,
G_ab the Einstein tensor of the metric h_ab on the wall and T_ab the stress tensor of its matter:
Deffayet's (junc) with [K] = 2K, 1/kappa^2 = c^4/(8 pi G_5) and 1/mu^2 = c^4/(8 pi G).

    gaussian_normal  Deffayet's (cosmback) and (bulkmet) with b = 1 and n_0 = 1, the normal
                     coordinate y written for his epsilon|y|, so that the side y >= 0 is the side
                     kept on the branch epsilon = +1 and the side y <= 0 the side kept on the branch
                     epsilon = -1: ds^2 = -N^2 c^2dt^2 + A^2 (unit 3-space of curvature k) + dy^2,
                     A = a + y sqrt(a'^2 + k), N = 1 + y a''/sqrt(a'^2 + k), a(t) free, the dot a
                     derivative along ct;
    global           the empty wall of the branch epsilon = +1, a = r_c cosh(ct/r_c) at k = 1,
                     (1 + |y|/r_c)^2 (-c^2dt^2 + r_c^2 cosh^2(ct/r_c) dOmega_3^2) + dy^2, both sides;
    planar           the same wall in flat slices, a = e^{ct/r_c} at k = 0,
                     (1 + |y|/r_c)^2 (-c^2dt^2 + e^{2ct/r_c} dx^2) + dy^2, both sides;
    inertial         Deffayet's Minkowski chart, his (canon), in which the wall is his (hyperbol),
                     the hyperboloid R^2 - c^2T^2 = r_c^2, and the side kept is R^2 - c^2T^2 >= r_c^2.

The curvature of the global and planar charts is a delta on the wall, which the checker reads as
a distribution, and their Kretschmann scalar, the square of that delta, is stated where the wall
is not. dgp_wall.md beside this file derives every chart and records every drawing decision.
"""
import sympy as sp

METRIC = "dgp_wall"
CHARTS = ("gaussian_normal", "global", "planar", "inertial")
NAMES = {"gaussian_normal": "Deffayet, Gaussian Normal", "global": "Global", "planar": "Planar",
         "inertial": "Inertial"}

SPHERE = "\\left(d\\theta^2 + \\sin^2\\theta\\,d\\phi^2\\right)"
THREE_SPHERE = "\\left(d\\psi^2 + \\sin^2\\psi" + SPHERE + "\\right)"
COMOVING = "\\left(\\dfrac{d\\chi^2}{1 - k\\chi^2} + \\chi^2" + SPHERE + "\\right)"
FLAT = "dx_1^2 + dx_2^2 + dx_3^2"
ROOT = "\\sqrt{\\dot{a}^2 + k}"
GN_PARAMETERS = ["k", "a = a(t)", "A = a + y" + ROOT, "N = 1 + \\dfrac{y\\ddot{a}}{" + ROOT + "}"]
COORDS = {"gaussian_normal": ["t", "y", "\\chi", "\\theta", "\\phi"],
          "global": ["t", "y", "\\psi", "\\theta", "\\phi"],
          "planar": ["t", "x_1", "x_2", "x_3", "y"],
          "inertial": ["T", "R", "\\psi", "\\theta", "\\phi"]}
PARAMETERS = {"gaussian_normal": GN_PARAMETERS, "global": ["r_c"], "planar": ["r_c"], "inertial": ["r_c"]}
WARP = "\\left(1 + \\dfrac{|y|}{r_c}\\right)^2"


def line(system, c, c2):
    """The line element of a chart, with c before t and c2 before dt^2."""
    if system == "gaussian_normal":
        return "ds^2 = -N^2" + c2 + "dt^2 + A^2" + COMOVING + " + dy^2"
    if system == "global":
        return ("ds^2 = " + WARP + "\\left(-" + c2 + "dt^2 + r_c^2\\cosh^2\\left(\\dfrac{" + c + "t}{r_c}\\right)"
                + THREE_SPHERE + "\\right) + dy^2")
    if system == "planar":
        return "ds^2 = " + WARP + "\\left(-" + c2 + "dt^2 + e^{2" + c + "t/r_c}\\left(" + FLAT + "\\right)\\right) + dy^2"
    return "ds^2 = -" + c2 + "dT^2 + dR^2 + R^2" + THREE_SPHERE


REALS = " \\in (-\\infty, \\infty)"
ANGLES = ["\\psi \\in [0, \\pi]", "\\theta \\in [0, \\pi]", "\\phi \\in [0, 2\\pi)"]
DOMAINS = {
    "gaussian_normal": ["t" + REALS, "y \\in [0, \\infty) \\;\\text{on the branch}\\; \\epsilon = +1",
                        "y \\in (-\\infty, 0] \\;\\text{on the branch}\\; \\epsilon = -1",
                        "\\chi \\in [0, 1/\\sqrt{k}) \\;\\text{for}\\; k > 0", "\\chi \\in [0, \\infty) \\;\\text{for}\\; k \\le 0",
                        "\\theta \\in [0, \\pi]", "\\phi \\in [0, 2\\pi)",
                        "y = 0 \\;\\text{(the wall, where the side kept meets its mirror image)}",
                        "N = 0 \\;\\text{and}\\; A = 0 \\;\\text{(coordinate singularities)}"],
    "global": ["t" + REALS, "y" + REALS] + ANGLES + ["y = 0 \\;\\text{(the wall)}"],
    "planar": ["t" + REALS, "x_1" + REALS, "x_2" + REALS, "x_3" + REALS, "y" + REALS, "y = 0 \\;\\text{(the wall)}"],
    "inertial": ["T" + REALS, "R \\in [\\sqrt{r_c^2 + c^2T^2}, \\infty)"] + ANGLES
    + ["R = \\sqrt{r_c^2 + c^2T^2} \\;\\text{(the wall, where this side meets its mirror image)}"],
}
_DIMS = {"t": "T", "T": "T", "y": "L", "R": "L", "x_1": "L", "x_2": "L", "x_3": "L", "\\chi": "1", "\\psi": "1",
         "\\theta": "1", "\\phi": "1", "k": "1", "a": "L", "A": "L", "N": "1", "r_c": "L"}
DIMENSIONS = {(METRIC, s): {name: _DIMS[name] for name in COORDS[s] + [p.split(" = ")[0] for p in PARAMETERS[s]]}
              for s in CHARTS}
HELD = {}
RATES = {}
DIGITS = 30
TOLERANCE = sp.Float("1e-22")


def register(module):
    module.CHARTS[METRIC] = [lambda system=system: chart_spec(system) for system in CHARTS]


def chart_spec(system):
    import chart_printer as cp
    import verify_metrics as vm

    coords = COORDS[system]
    probe = vm.Reader(coords, PARAMETERS[system], ())
    spec = {"metric_id": METRIC,
            "system": {"id": system, "name": NAMES[system], "coords": coords, "domains": DOMAINS[system],
                       "parameters": PARAMETERS[system], "line_element": line(system, "c", "c^2")},
            "chart_line_element": line(system, "", ""),
            "time": coords[0],
            "check": lambda chart: chart_check(chart, system)}
    if system == "gaussian_normal":
        D = sp.Derivative
        t, y = probe.symbol["t"], probe.symbol["y"]
        A, N, k = probe.parameters["A"], probe.parameters["N"], probe.parameters["k"]
        lead = [D(N, (t, 2)), D(N, t, y), D(N, (y, 2)), D(A, (t, 2)), D(A, t, y), D(A, (y, 2)),
                D(N, t), D(N, y), D(A, t), D(A, y), N, A, k]
        spec["printer"] = {"lead": lead, "rising": [k], "flip": False}
        spec["components"] = {"metric_components": {("\\chi", "\\chi"): "\\dfrac{A^2}{1 - k\\chi^2}"},
                              "inverse_metric_components": {("\\chi", "\\chi"): "\\dfrac{1 - k\\chi^2}{A^2}"}}
        return spec
    r_c, c, time_ = probe.parameters["r_c"], probe.c, probe.symbol[coords[0]]
    if system == "inertial":
        spec["printer"] = {"lead": [r_c, c, time_]}
        return spec
    y = probe.symbol["y"]
    then = cp.hyperbolic(c * time_ / r_c) if system == "global" else sp.factor
    pretty, overrides = cp.kink(y, then)
    absolute = next(p for p in overrides if p.name.startswith("_abs"))
    spec["pretty"] = pretty
    spec["printer"] = {"lead": [r_c, c, time_], "rising": [absolute], "overrides": overrides}
    spec["kretschmann_where"] = "y"
    # ct/r_c is set as the line element sets it, and sinh over cosh as tanh.
    arg = "\\left(\\dfrac{ct}{r_c}\\right)"
    tangent = [(f"\\dfrac{{{n}\\sinh{arg}}}{{r_c\\cosh{arg}}}", f"\\dfrac{{{n or 1}}}{{r_c}}\\tanh{arg}") for n in ("", "2")]
    spec["rewrite"] = [("c\\,t/r_c", "ct/r_c"), ("\\dfrac{c\\,t}{r_c}", "\\dfrac{ct}{r_c}")] + tangent
    return spec


# -- the checks ------------------------------------------------------------------------------------
#
# Each chart is checked independently of its printed tensors: the Gaussian normal chart to a
# vanishing Riemann tensor for every a(t), with A and N written out, and to the junction condition,
# which with the density rho and the pressure p of the matter on the wall is Deffayet's Friedmann
# equation on either branch; the global and planar charts to flat space of five dimensions pulled
# back along the hyperboloid of radius r_c on both sides of the wall and to the delta the junction
# condition asks of an empty wall; the inertial chart to flat space and to the wall at
# R^2 - c^2T^2 = r_c^2; and the Gaussian normal chart at a = r_c cosh(ct/r_c), k = 1 to the global
# chart on the side y >= 0.

LIGHT = sp.Rational(3, 2)                    # c, so that a factor of c that is lost is caught
RADII = [sp.Rational(7, 5), sp.Rational(5, 8)]


def _zero(value, what):
    if abs(complex(sp.N(value, DIGITS))) > TOLERANCE:
        raise AssertionError(f"{METRIC}: {what} fails, by {complex(sp.N(value, DIGITS))}")


def _unit_three_sphere(psi, theta, phi):
    return [sp.cos(psi), sp.sin(psi) * sp.cos(theta), sp.sin(psi) * sp.sin(theta) * sp.cos(phi),
            sp.sin(psi) * sp.sin(theta) * sp.sin(phi)]


def _hyperboloid(system, X, r_c):
    """The point of the unit hyperboloid -Y_0^2 + Y_1^2 + ... + Y_4^2 = 1 a chart's wall
    coordinates stand at, x^0 = ct in units of r_c."""
    tau = X[0] / r_c
    if system == "global":
        n = _unit_three_sphere(*X[2:5])
        return [sp.sinh(tau)] + [sp.cosh(tau) * e for e in n]
    u = [x / r_c for x in X[1:4]]
    u2 = sum(e ** 2 for e in u)
    return ([sp.sinh(tau) + u2 * sp.exp(tau) / 2] + [e * sp.exp(tau) for e in u]
            + [sp.cosh(tau) - u2 * sp.exp(tau) / 2])


def _minkowski_pullback(chart, system, image, at, what):
    """The chart's metric against flat space of five dimensions pulled back along `image`,
    the Cartesian X^0 = cT, X^1, ..., X^4 of each event, at the event `at`."""
    X = chart.symbols
    J = sp.Matrix([[sp.diff(p, x) for x in X] for p in image])
    pulled = J.T * sp.diag(-1, 1, 1, 1, 1) * J
    g = sp.Matrix(chart.geo.g)
    for a in range(5):
        for b in range(a, 5):
            _zero((pulled[a, b] - g[a, b]).subs(at), f"{system}: {what}, slot {chart.coords_tex[a]}{chart.coords_tex[b]}")


def _induced_einstein(h, coords):
    """The mixed Einstein tensor G^a_b of the metric h of the wall, in the coordinates given."""
    import verify_metrics as vm
    geo = vm.Geometry(h, coords, 10 ** 6)
    lower = geo.einstein_ll()
    return sp.Matrix(4, 4, lambda a, b: sp.simplify(sum(geo.ginv[a, c] * lower[c][b] for c in range(4))))


def _gaussian_normal(chart):
    import verify_metrics as vm
    reader, geo = chart.reader, chart.geo
    t, y, chi, theta, phi = chart.symbols
    # The Riemann tensor vanishes for every a(t), k and either side, with A and N written out.
    riemann = geo.riemann_llll()
    for index in vm._indices(5, 4):
        if list(index[:2]) > list(index[2:]) or index[0] >= index[1] or index[2] >= index[3]:
            continue
        if vm.norm(reader.surface(vm._at(riemann, index))) != 0:
            raise AssertionError(f"{METRIC}: the Gaussian normal chart's Riemann tensor is not zero in "
                                 f"slot {[chart.coords_tex[i] for i in index]}")
    a = reader.functions["a"][0]
    k = reader.parameters["k"]
    g = sp.Matrix(geo.g).applyfunc(reader.surface)
    # On the wall y = 0 the chart is Friedmann's universe, with t the proper time.
    h = g[[0, 2, 3, 4], [0, 2, 3, 4]].subs(y, 0)
    want = sp.diag(-1, a ** 2 / (1 - k * chi ** 2), a ** 2 * chi ** 2, a ** 2 * chi ** 2 * sp.sin(theta) ** 2)
    if vm.norm(h - want) != sp.zeros(4):
        raise AssertionError(f"{METRIC}: the Gaussian normal chart's metric on the wall is not Friedmann's")
    hinv = h.inv()
    einstein = _induced_einstein(h, [t, chi, theta, phi])
    r_c, G = sp.symbols("r_c G", positive=True)
    G5 = 2 * r_c * G
    c = reader.c
    ad, add = sp.diff(a, t), sp.diff(a, t, 2)
    S = sp.sqrt(ad ** 2 + k)

    def residual(eps, branch):
        # K_ab on the side kept, the normal pointing away from the wall: eps d_y g_ab / 2 at y = 0.
        K = (eps * sp.diff(g, y) / 2)[[0, 2, 3, 4], [0, 2, 3, 4]].subs(y, 0)
        mixed = hinv * K
        trace = mixed.trace()
        # The density from the Friedmann equation of the branch, the pressure from its conservation,
        # d rho / d(ct) = -3 (a'/a)(rho + p/c^2), with the dot along ct.
        rho = 3 * c ** 2 / (8 * sp.pi * G) * (S ** 2 / a ** 2 - branch * S / (r_c * a))
        p = -c ** 2 * (rho + sp.diff(rho, t) * a / (3 * ad))
        T = sp.diag(-rho * c ** 2, p, p, p)
        lhs = mixed - sp.eye(4) * trace
        rhs = r_c * einstein - 4 * sp.pi * G5 * T / c ** 4
        return [sp.simplify(e) for e in (lhs - rhs)]
    for eps in (1, -1):
        if any(e != 0 for e in residual(eps, eps)):
            raise AssertionError(f"{METRIC}: the junction condition is not the Friedmann equation of the branch {eps}")
        # The check can fail: the Friedmann equation of the other branch does not satisfy it.
        if all(e == 0 for e in residual(eps, -eps)):
            raise AssertionError(f"{METRIC}: the junction condition holds on the wrong branch")
    # At a = r_c cosh(ct/r_c) and k = 1 the side y >= 0 is the global chart, with chi = sin psi.
    import chart_printer as cp
    other = chart_spec("global")
    target = cp.Chart(other["system"]["coords"], other["system"]["parameters"], other["chart_line_element"])
    Y = target.symbols
    rc = target.reader.parameters["r_c"]
    J = sp.diag(1, 1, sp.cos(Y[2]), 1, 1)
    values = {a: rc * sp.cosh(t / rc)}
    explicit = g.subs(k, 1)
    explicit = explicit.subs(values).doit()
    moved = explicit.subs({t: Y[0], y: Y[1], chi: sp.sin(Y[2]), theta: Y[3], phi: Y[4]}, simultaneous=True)
    pulled = J.T * moved * J
    theirs = sp.Matrix(target.geo.g)
    for radius in RADII:
        for side in (sp.Rational(1, 3), sp.Rational(5, 2)):
            at = {rc: radius, Y[0]: sp.Rational(-2, 7), Y[1]: side, Y[2]: sp.Rational(3, 5), Y[3]: sp.Rational(4, 7),
                  Y[4]: sp.Rational(1, 9)}
            for a_ in range(5):
                for b_ in range(a_, 5):
                    _zero((pulled[a_, b_] - theirs[a_, b_]).subs(at),
                          f"gaussian_normal: the global chart at a = r_c cosh(ct/r_c), slot {a_}{b_}")
    # The Riemann tensor, now known to vanish with A and N written out, stands as zero, so that every
    # curvature tensor is printed as the bulk has it while the Christoffel symbols keep the names.
    onshell = vm.Geometry(geo.g, chart.symbols, geo.seconds)
    onshell._cache.update({"riemann_llll": onshell._zeros(4), "riemann_ulll": onshell._zeros(4)})
    chart.geo = onshell


def _empty_wall(chart, system):
    """Both sides against flat space along the hyperboloid, and the wall's delta in G_mu nu against
    the junction condition of an empty wall, -2(K_ab - h_ab K) = -2 r_c G_ab, which is (6/r_c) h_ab."""
    X = chart.symbols
    r_c, c = chart.reader.parameters["r_c"], chart.reader.c
    y = chart.reader.symbol["y"]
    across = 1 if system == "global" else 4
    wall = [i for i in range(5) if i != across]
    for radius in RADII:
        base = {r_c: radius, chart.reader.c: LIGHT}
        for side in (sp.Rational(2, 5), sp.Rational(-3, 4)):
            if system == "global":
                at = {**base, X[0]: sp.Rational(1, 3), X[1]: side, X[2]: sp.Rational(3, 5), X[3]: sp.Rational(4, 7),
                      X[4]: sp.Rational(1, 9)}
            else:
                at = {**base, X[0]: sp.Rational(1, 3), X[1]: sp.Rational(1, 7), X[2]: sp.Rational(-2, 9),
                      X[3]: sp.Rational(5, 11), X[4]: side}
            # On either side the point at distance |y| from the wall along its normal.
            image = [(r_c + sp.Abs(y)) * Y for Y in _hyperboloid(system, X, r_c)]
            _minkowski_pullback(chart, system, image, at, f"flat space along the hyperboloid at y = {side}")
    g = sp.Matrix(chart.geo.g)
    einstein = chart.geo.einstein_ll()
    Dd = sp.Symbol("_D")
    near = sp.Rational(1, 10 ** 40)
    h = g[wall, wall].subs(y, 0)
    induced = _induced_einstein(h, [X[i] for i in wall])
    # K_ab on the side y > 0, the normal d_y pointing away from the wall; the side y < 0 is its mirror image.
    outside = g.subs(sp.sign(y), 1).subs(sp.Abs(y), y)
    K = (sp.diff(outside, y) / 2)[wall, wall].subs(y, 0)
    jump_form = h.inv() * K - sp.eye(4) * (h.inv() * K).trace()
    for i in range(4):
        for j in range(4):
            if sp.simplify(jump_form[i, j] - r_c * induced[i, j]) != 0:
                raise AssertionError(f"{METRIC}: {system}: K_ab - h_ab K is not r_c G_ab on the wall")
    for i, a_ in enumerate(wall):
        for j, b_ in enumerate(wall):
            value = sp.sympify(einstein[a_][b_]).xreplace({sp.DiracDelta(y): Dd})
            on = sp.diff(value, Dd).subs(Dd, 0)
            want = -2 * r_c * (h * induced)[i, j]
            at = {r_c: RADII[0], chart.reader.c: LIGHT, y: near}
            at.update({X[m]: v for m, v in zip(wall, (sp.Rational(1, 3), sp.Rational(3, 5), sp.Rational(4, 7),
                                                       sp.Rational(1, 9)))})
            _zero((on - want).subs(at), f"{system}: the wall's delta in G_{chart.coords_tex[a_]}{chart.coords_tex[b_]}")
            _zero((on - 6 * h[i, j] / r_c).subs(at), f"{system}: the wall's delta is not (6/r_c) h_ab")
    # The check can fail: a wall of the branch epsilon = -1, which keeps the inside, misses the delta.
    if sp.simplify(-jump_form[0, 0] - r_c * induced[0, 0]) == 0:
        raise AssertionError(f"{METRIC}: {system}: the inside satisfies the junction condition of an empty wall")


def _inertial(chart):
    X = chart.symbols
    r_c = chart.reader.parameters["r_c"]
    n = _unit_three_sphere(*X[2:5])
    image = [X[0]] + [X[1] * e for e in n]
    for radius in RADII:
        at = {r_c: radius, chart.reader.c: LIGHT, X[0]: sp.Rational(1, 3), X[1]: sp.Rational(9, 4),
              X[2]: sp.Rational(3, 5), X[3]: sp.Rational(4, 7), X[4]: sp.Rational(1, 9)}
        _minkowski_pullback(chart, "inertial", image, at, "flat space in spherical coordinates")
    # The wall of the global chart, y = 0, is R^2 - c^2T^2 = r_c^2 here, and its side y > 0 lies outside.
    tau, side = sp.symbols("_tau _side", positive=True)
    cT, R = (r_c + side) * sp.sinh(tau), (r_c + side) * sp.cosh(tau)
    gap = sp.factor(sp.simplify(R ** 2 - cT ** 2 - r_c ** 2))
    if gap.subs(side, 0) != 0 or sp.simplify(gap - side * (2 * r_c + side)) != 0:
        raise AssertionError(f"{METRIC}: inertial: the wall is not R^2 - c^2T^2 = r_c^2 with the side kept outside")


def chart_check(chart, system):
    if system == "gaussian_normal":
        _gaussian_normal(chart)
    elif system == "inertial":
        _inertial(chart)
    else:
        _empty_wall(chart, system)
