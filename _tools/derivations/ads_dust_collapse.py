"""Dust collapsing under a negative cosmological constant, in eight charts.

Three dimensions, Lambda = -1/ell^2. R. B. Mann and S. F. Ross, Phys. Rev. D 47, 3319 (1993):

    disc_interior           their (3), ds^2 = -dt^2 + a^2 (dr^2/(1 - k r^2) + r^2 dtheta^2), with their (4),
                            a'' = -a/ell^2, their (6) from rest, a = cos(t/ell), and their (9) for k;
    disc_exterior           their (2), the hole of Banados, Teitelboim and Zanelli without rotation,
                            R^2/ell^2 - M, with their (21), M = 8 pi G_3 rho_0 r_0^2 - 1;
    disc_exterior_ingoing   the same exterior on the advanced time dv = dT + dR/(R^2/ell^2 - M), the
                            form S. Carlip writes in Class. Quantum Grav. 12, 2853 (1995), section 2.

Four dimensions, Lambda = -3/L^2. W. L. Smith and R. B. Mann, Phys. Rev. D 56, 4942 (1997), and
F. C. Mena, J. Natario and P. Tod, Adv. Theor. Math. Phys. 12, 1163 (2008):

    hyperbolic_interior     Smith and Mann's (4) with their k written 1/b^2, the open Friedmann universe of
                            dust sliced by hyperbolic planes, with their Friedmann equation (6);
    hyperbolic_slab         the same cloud in r = b cosh(chi/b), Mena, Natario and Tod's Table 1 (c),
                            through the smallest surface r = b to the second face;
    hyperbolic_exterior     Smith and Mann's (2), -1 - mu/R + R^2/L^2, with their (28) for the mass;
    flat_interior           Mena, Natario and Tod's section 3.2.1, case (iii), the same universe sliced
                            by horospheres, a^2 (dchi^2 + e^(2 chi/b) (dx^2 + dy^2));
    flat_exterior           their section 3.1 at b = 0 in their letters, -mu/R + R^2/L^2 over a plane.

J. P. S. Lemos's interior of Phys. Rev. D 57, 4600 (1998), his (5), is no dust solution and is not
printed: its space is a cone over a plane, with G^l_l - G^theta_theta = 1/(l^2 a^2).

The scale factor is a pure number with a = 1 at release from rest. In x^0 = ct, with the density at
release rho_0,

    three dimensions:  a'' = -a/ell^2,  a'^2 = (1 - a^2)/ell^2,  k = 8 pi G_3 rho_0/c^2 - 1/ell^2,
    four dimensions:   a'^2 = 1/b^2 + (1/L^2 - 1/b^2)/a - a^2/L^2,  8 pi G rho_0/c^2 = 3 (1/L^2 - 1/b^2).

Every check is made at points. Each interior is held to Einstein's equations with dust and its
cosmological constant, G^t_t + Lambda = -kappa rho and G^i_j + Lambda delta^i_j = 0, the slab to being
Smith and Mann's chart pulled back, each exterior to G_ab + Lambda g_ab = 0, to the published chart of
btz or of topological_black_hole it is, and to the junction: across the surface the metric of the
surface is continuous, and the surface's areal radius R_s = r_0 a obeys

    (dR_s/d(c tau))^2 + f(R_s) = E^2,

the radial geodesic of the exterior whose energy E is the interior's sqrt(1 - k r_0^2), sqrt(r_0^2/b^2 - 1)
or R_0/b, which is the continuity of the extrinsic curvature. ads_dust_collapse.md beside this file
records each chart's source.
"""
import json
from pathlib import Path

import sympy as sp

METRIC = "ads_dust_collapse"
DISC, DISC_OUT, DISC_IN = "disc_interior", "disc_exterior", "disc_exterior_ingoing"
HYP, SLAB, HYP_OUT = "hyperbolic_interior", "hyperbolic_slab", "hyperbolic_exterior"
FLAT, FLAT_OUT = "flat_interior", "flat_exterior"
SYSTEMS = (DISC, DISC_OUT, DISC_IN, HYP, SLAB, HYP_OUT, FLAT, FLAT_OUT)
INTERIORS = (DISC, HYP, SLAB, FLAT)
NAMES = {DISC: "Comoving Disc of Dust", DISC_OUT: "Static Exterior of the Disc",
         DISC_IN: "Ingoing Exterior of the Disc", HYP: "Comoving Hyperbolic Cloud",
         SLAB: "Hyperbolic Cloud between Its Faces", HYP_OUT: "Static Exterior of the Hyperbolic Cloud",
         FLAT: "Comoving Flat Cloud", FLAT_OUT: "Static Exterior of the Flat Cloud"}
COORDS = {DISC: ["t", "r", "\\theta"], DISC_OUT: ["T", "R", "\\theta"], DISC_IN: ["v", "R", "\\theta"],
          HYP: ["t", "r", "\\theta", "\\phi"], SLAB: ["t", "\\chi", "\\theta", "\\phi"],
          HYP_OUT: ["T", "R", "\\theta", "\\phi"], FLAT: ["t", "\\chi", "x", "y"],
          FLAT_OUT: ["T", "R", "x", "y"]}
PARAMETERS = {DISC: ["a = a(t)", "k", "r_0", "\\ell"], DISC_OUT: ["\\ell", "M"], DISC_IN: ["\\ell", "M"],
              HYP: ["a = a(t)", "b", "r_0", "L"], SLAB: ["a = a(t)", "b", "\\chi_0", "L"],
              HYP_OUT: ["\\mu", "L"], FLAT: ["a = a(t)", "b", "\\chi_0", "L"], FLAT_OUT: ["\\mu", "L", "b"]}

BTZ = "\\dfrac{R^2}{\\ell^2} - M"
HOLE = "-1 - \\dfrac{\\mu}{R} + \\dfrac{R^2}{L^2}"
BRANE = "\\dfrac{R^2}{L^2} - \\dfrac{\\mu}{R}"
PLANES = "r^2\\left(d\\theta^2 + \\sinh^2\\theta\\,d\\phi^2\\right)"


def _static(f, rest, c2="c^2"):
    return f"ds^2 = -\\left({f}\\right){c2}dT^2 + \\dfrac{{dR^2}}{{{f}}} + {rest}"


# Each line element as the page prints it; the chart x^0 = ct drops the c beside the time.
LINES = {
    DISC: "ds^2 = -c^2dt^2 + a^2\\left(\\dfrac{dr^2}{1 - kr^2} + r^2d\\theta^2\\right)",
    DISC_OUT: _static(BTZ, "R^2d\\theta^2"),
    DISC_IN: f"ds^2 = -\\left({BTZ}\\right)dv^2 + 2\\,dv\\,dR + R^2d\\theta^2",
    HYP: f"ds^2 = -c^2dt^2 + a^2\\left(\\dfrac{{dr^2}}{{\\dfrac{{r^2}}{{b^2}} - 1}} + {PLANES}\\right)",
    SLAB: ("ds^2 = -c^2dt^2 + a^2\\left(d\\chi^2 + b^2\\cosh^2(\\chi/b)"
           "\\left(d\\theta^2 + \\sinh^2\\theta\\,d\\phi^2\\right)\\right)"),
    HYP_OUT: _static(HOLE, "R^2\\left(d\\theta^2 + \\sinh^2\\theta\\,d\\phi^2\\right)"),
    FLAT: "ds^2 = -c^2dt^2 + a^2\\left(d\\chi^2 + e^{2\\chi/b}\\left(dx^2 + dy^2\\right)\\right)",
    FLAT_OUT: _static(BRANE, "\\dfrac{R^2}{b^2}\\left(dx^2 + dy^2\\right)"),
}
CHART_LINES = {system: line.replace("c^2dt^2", "dt^2").replace("c^2dT^2", "dT^2") for system, line in LINES.items()}

LINE = "(-\\infty, \\infty)"
ANGLE = "\\theta \\in [0, 2\\pi)"
PLANE = ["\\theta \\in [0, \\infty)", "\\phi \\in [0, 2\\pi)"]
DOMAINS = {
    DISC: ["t \\in [0, \\pi\\ell/(2c))", "r \\in [0, r_0]", ANGLE],
    DISC_OUT: ["T \\in [0, \\infty)", "R \\in [R_s(T), \\infty)", ANGLE,
               "R = \\ell\\sqrt{M} \\;\\text{(the horizon, for}\\; M > 0\\text{)}"],
    DISC_IN: ["v \\in " + LINE, "R \\in [R_s(v), \\infty)", ANGLE, "T \\ge 0 \\;\\text{(after the moment of rest)}",
              "R = \\ell\\sqrt{M} \\;\\text{(the horizon, for}\\; M > 0\\text{)}"],
    HYP: ["t \\in [0, t_c)", "r \\in (b, r_0]"] + PLANE,
    SLAB: ["t \\in [0, t_c)", "\\chi \\in [-\\chi_0, \\chi_0]"] + PLANE,
    HYP_OUT: ["T \\in [0, \\infty)", "R \\in [R_s(T), \\infty)"] + PLANE + ["R = R_h \\;\\text{(the horizon)}"],
    FLAT: ["t \\in [0, t_c)", "\\chi \\in (-\\infty, \\chi_0]", "x \\in " + LINE, "y \\in " + LINE],
    FLAT_OUT: ["T \\in [0, \\infty)", "R \\in [R_s(T), \\infty)", "x \\in " + LINE, "y \\in " + LINE,
               "R = (\\mu L^2)^{1/3} \\;\\text{(the horizon)}"],
}
CLOUD = {"ricci_scalar": "\\dfrac{6\\left(b^2\\left(a\\ddot{a} + \\dot{a}^2\\right) - 1\\right)}{a^2b^2}",
         "kretschmann": "\\dfrac{12\\left(a^2b^4\\ddot{a}^2 + \\left(b^2\\dot{a}^2 - 1\\right)^2\\right)}{a^4b^4}"}
SCALARS = {
    DISC: {"ricci_scalar": "\\dfrac{2\\left(2a\\ddot{a} + \\dot{a}^2 + k\\right)}{a^2}",
           "kretschmann": "\\dfrac{4\\left(2a^2\\ddot{a}^2 + \\left(\\dot{a}^2 + k\\right)^2\\right)}{a^4}"},
    HYP: CLOUD, SLAB: CLOUD, FLAT: CLOUD,
    DISC_OUT: {"ricci_scalar": "-\\dfrac{6}{\\ell^2}", "kretschmann": "\\dfrac{12}{\\ell^4}"},
    DISC_IN: {"ricci_scalar": "-\\dfrac{6}{\\ell^2}", "kretschmann": "\\dfrac{12}{\\ell^4}"},
    HYP_OUT: {"ricci_scalar": "-\\dfrac{12}{L^2}", "kretschmann": "\\dfrac{12\\mu^2}{R^6} + \\dfrac{24}{L^4}"},
    FLAT_OUT: {"ricci_scalar": "-\\dfrac{12}{L^2}", "kretschmann": "\\dfrac{12\\mu^2}{R^6} + \\dfrac{24}{L^4}"},
}

DIMENSIONS = {
    (METRIC, DISC): {"t": "T", "r": "L", "\\theta": "1", "a": "1", "k": "L**-2", "r_0": "L", "\\ell": "L"},
    (METRIC, DISC_OUT): {"T": "T", "R": "L", "\\theta": "1", "\\ell": "L", "M": "1"},
    (METRIC, DISC_IN): {"v": "L", "R": "L", "\\theta": "1", "\\ell": "L", "M": "1"},
    (METRIC, HYP): {"t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "a": "1", "b": "L", "r_0": "L", "L": "L"},
    (METRIC, SLAB): {"t": "T", "\\chi": "L", "\\theta": "1", "\\phi": "1", "a": "1", "b": "L", "\\chi_0": "L",
                     "L": "L"},
    (METRIC, HYP_OUT): {"T": "T", "R": "L", "\\theta": "1", "\\phi": "1", "\\mu": "L", "L": "L"},
    (METRIC, FLAT): {"t": "T", "\\chi": "L", "x": "L", "y": "L", "a": "1", "b": "L", "\\chi_0": "L", "L": "L"},
    (METRIC, FLAT_OUT): {"T": "T", "R": "L", "x": "L", "y": "L", "\\mu": "L", "L": "L", "b": "L"},
}


def hyperbolics(*arguments):
    """A `pretty` for a chart whose values are products of powers of the hyperbolic sine and cosine of each
    argument given, which the checker's Geometry hands back as exponentials. With E = e^argument a value is
    a rational function of each E; it is factored, E^2 - 1 is written 2E sinh, E^2 + 1 is written 2E cosh,
    and the powers of E that are left have to cancel, or the value is no such product and this stops."""
    names = [sp.Symbol(f"ADC_E{i}", positive=True) for i in range(len(arguments))]

    def pretty(value):
        x = sp.sympify(value)
        for argument, E in zip(arguments, names):
            x = x.replace(lambda e, w=argument: isinstance(e, sp.exp) and sp.expand(e.args[0] / w).is_Integer,
                          lambda e, w=argument, E=E: E ** sp.expand(e.args[0] / w))
        x = sp.factor(sp.cancel(sp.together(x)))
        if not any(x.has(E) for E in names):
            return x
        out, left = sp.Integer(1), {E: sp.Integer(0) for E in names}
        halves = {}
        for factor in sp.Mul.make_args(x):
            base, power = (factor.base, factor.exp) if factor.is_Pow else (factor, sp.Integer(1))
            found = False
            for argument, E in zip(arguments, names):
                if base == E:
                    left[E] += power
                elif sp.expand(base - (E ** 2 + 1)) == 0:
                    out *= (2 * sp.cosh(argument)) ** power
                    left[E] += power
                elif sp.expand(base - (E ** 2 - 1)) == 0:
                    out *= (2 * sp.sinh(argument)) ** power
                    left[E] += power
                elif sp.expand(base - (E - 1)) == 0 or sp.expand(base - (E + 1)) == 0:
                    halves.setdefault(E, {})[sp.expand(base)] = power
                else:
                    continue
                found = True
                break
            if not found:
                if any(base.has(E) for E in names):
                    raise AssertionError(f"{METRIC}: {value} is no product of hyperbolic functions")
                out *= factor
        for argument, E in zip(arguments, names):
            pair = halves.get(E, {})
            if pair:
                if len(pair) != 2 or len(set(pair.values())) != 1:
                    raise AssertionError(f"{METRIC}: {value} holds a hyperbolic function of half an argument")
                power = next(iter(pair.values()))
                out *= (2 * sp.sinh(argument)) ** power
                left[E] += power
            if left[E] != 0:
                raise AssertionError(f"{METRIC}: {value} keeps a bare exponential of {argument}")
        return out
    return pretty


def squares(e):
    """(u + v)(u - v) written back as u^2 - v^2, for every pair of two-term factors of one power whose
    product has two terms."""
    e = sp.sympify(e)
    if not e.is_Mul:
        return e
    items = [[f.base, f.exp] if f.is_Pow else [f, sp.Integer(1)] for f in sp.Mul.make_args(e)]
    changed = True
    while changed:
        changed = False
        for i in range(len(items)):
            for j in range(i + 1, len(items)):
                (b1, k1), (b2, k2) = items[i], items[j]
                if k1 != k2 or not (b1.is_Add and b2.is_Add) or len(b1.args) != 2 or len(b2.args) != 2:
                    continue
                product = sp.expand(b1 * b2)
                if len(sp.Add.make_args(product)) == 2:
                    items[i] = [product, k1]
                    del items[j]
                    changed = True
                    break
            if changed:
                break
    return sp.Mul(*[b ** k for b, k in items])


def spec(system):
    import verify_metrics as vm
    probe = vm.Reader(COORDS[system], PARAMETERS[system], ())
    p, s = probe.parameters, probe.symbol
    out = {
        "metric_id": METRIC,
        "system": {"id": system, "name": NAMES[system], "coords": COORDS[system], "domains": DOMAINS[system],
                   "parameters": PARAMETERS[system], "line_element": LINES[system]},
        "chart_line_element": CHART_LINES[system],
        "check": lambda chart: check(chart, system),
        **SCALARS.get(system, {}),
    }
    if system in INTERIORS:
        radius = s["r"] if system in (DISC, HYP) else s["\\chi"]
        others = [p[name] for name in ("k", "b") if name in p]
        # A product is written a, then the curvature, then the radius, and the disc's sums as 1 - k r^2.
        out["printer"] = {"dotted": ["a"], "lead": [p["a"], radius] + others, "factors": [p["a"]] + others + [radius]}
        if system == DISC:
            out["printer"]["rising"] = [radius]
        inner = sp.factor
        if system == HYP:
            inner = hyperbolics(s["\\theta"])
        elif system == SLAB:
            inner = hyperbolics(s["\\chi"] / p["b"], s["\\theta"])
        # (u + v)(u - v) is written u^2 - v^2, as r^2 - b^2 and b^2 a'^2 - 1 are.
        out["pretty"] = lambda value: squares(inner(value))
        if system == SLAB:
            out["printer"]["arguments"] = {s["\\chi"] / p["b"]: "\\chi/b"}
            cosh = "\\cosh^2(\\chi/b)"
            out["rewrite"] = [(f"\\left(b^2\\,\\dot{{a}}^2\\,{cosh}{tail} - {cosh}{tail}\\right)",
                               f"\\left(b^2\\,\\dot{{a}}^2 - 1\\right){cosh}{tail}")
                              for tail in ("\\sinh^2\\theta", "")]
        if system == DISC:
            out["rewrite"] = [("a\\,\\ddot{a} + k + \\dot{a}^2", "a\\,\\ddot{a} + \\dot{a}^2 + k"),
                              ("k + \\dot{a}^2", "\\dot{a}^2 + k"),
                              ("-\\dfrac{\\left(k\\,r^2 - 1\\right)\\ddot{a}}{a^3}",
                               "\\dfrac{\\left(1 - k\\,r^2\\right)\\ddot{a}}{a^3}"),
                              ("\\dfrac{\\left(k\\,r^2 - 1\\right)\\ddot{a}}{a^3}",
                               "-\\dfrac{\\left(1 - k\\,r^2\\right)\\ddot{a}}{a^3}")]
        # A sum that has to stand in a bracket keeps its hyperbolic functions.
        out["bracketed"] = lambda value: sp.expand(out["pretty"](value))
    else:
        if system in (DISC_OUT, DISC_IN):
            lead = [s["R"], p["M"], p["ell"]]
            out["printer"] = {"lead": lead, "factors": lead[1:] + lead[:1], "flip": False}
        else:
            # Every value is printed around the cubic whose root is the horizon, as the published chart of
            # the topological black hole is.
            rest = [p["b"]] if "b" in p else []
            out["printer"] = {"rising": [p["L"], p["mu"]], "lead": [p["L"], s["R"], p["mu"]] + rest, "flip": False}
        if system == HYP_OUT:
            out["pretty"] = hyperbolics(s["\\theta"])
            out["bracketed"] = lambda value: sp.expand(out["pretty"](value))
        f = {DISC_OUT: BTZ, DISC_IN: BTZ, HYP_OUT: HOLE, FLAT_OUT: BRANE}[system]
        wrapped = f"\\left({f}\\right)"
        if system == DISC_IN:
            out["components"] = {"metric_components": {("v", "v"): "-" + wrapped},
                                 "inverse_metric_components": {("R", "R"): f}}
        else:
            out["components"] = {"metric_components": {("T", "T"): "-" + wrapped, ("R", "R"): wrapped + "^{-1}"},
                                 "inverse_metric_components": {("T", "T"): "-" + wrapped + "^{-1}", ("R", "R"): f}}
    return out


def register(module):
    module.CHARTS[METRIC] = [lambda system=system: spec(system) for system in SYSTEMS]


# -- Checks at points, independent of the printed tensors ---------------------------------------

TOLERANCE = 1e-25
METRICS = Path(__file__).resolve().parents[2] / "MFS" / "assets" / "data" / "metrics"
Q = sp.Rational
# (a, the other values) at which each chart is held; a is the scale factor of the moment.
SCALES = (Q(9, 10), Q(3, 5), Q(1, 4))


def small(values, where, what):
    worst = max(abs(complex(sp.N(x, 40))) for x in sp.flatten([values]))
    if worst > TOLERANCE:
        raise AssertionError(f"{METRIC}/{where}: {what} ({worst:.3g})")


def friedmann(a, b, L):
    """The square of da/d(ct) for the clouds of four dimensions released from rest at a = 1."""
    return 1 / b ** 2 + (1 / L ** 2 - 1 / b ** 2) / a - a ** 2 / L ** 2


def mixed_einstein(chart, at):
    geo = chart.geo
    G = geo.einstein_ll()
    g = chart.geo.g.subs(at)
    n = len(chart.symbols)
    lower = sp.Matrix(n, n, lambda i, j: sp.sympify(G[i][j]).subs(at))
    return g.inv() * lower, g


def published(vm, metric_id, system_id):
    """A published chart's metric, its reader and its coordinates' symbols."""
    entry = next(c for c in json.loads((METRICS / f"{metric_id}.json").read_text(encoding="utf-8"))["coordinates"]
                 if c["id"] == system_id)
    reader = vm.Reader(entry["coords"], [q["symbol"] for q in entry["parameters"]], ())
    n = len(entry["coords"])
    g = sp.zeros(n, n)
    for e in entry["metric_components"]:
        i, j = (entry["coords"].index(x) for x in e["indices"])
        g[i, j] = reader(e["value"])
    return reader, g, [reader.symbol[c] for c in entry["coords"]]


def check(chart, system):
    import verify_metrics as vm
    s, p = chart.symbols, chart.reader.parameters
    n = len(s)
    if system in INTERIORS:
        interior(chart, system, s, p, n)
    else:
        exterior(vm, chart, system, s, p, n)


def interior(chart, system, s, p, n):
    a = p["a"]
    t = s[0]
    first, second = sp.Derivative(a, t), sp.Derivative(a, (t, 2))
    for a0 in SCALES:
        if system == DISC:
            ell, r0 = Q(7, 5), Q(9, 7)
            for k in (Q(-3, 11), Q(2, 9)):
                velocity, pull = -sp.sqrt((1 - a0 ** 2) / ell ** 2), -a0 / ell ** 2
                at = [(second, pull), (first, velocity), (a, a0), (p["k"], k), (p["r_0"], r0), (p["ell"], ell),
                      (s[0], Q(1, 3)), (s[1], Q(4, 5)), (s[2], Q(2, 7))]
                G, g = mixed_einstein(chart, at)
                Lambda, density = -1 / ell ** 2, (k + 1 / ell ** 2) / a0 ** 2
                where = f"{system} at a = {a0}, k = {k}"
                want = sp.diag(-density - Lambda, -Lambda, -Lambda)
                small(G - want, where, "the disc is not dust with its cosmological constant")
            # Mann and Ross's (6) from rest: a = cos(ct/ell) solves both equations.
            u = sp.Symbol("u")
            cosine = sp.cos(u / ell)
            small([(sp.diff(cosine, u, 2) + cosine / ell ** 2).subs(u, Q(2, 3)),
                   (sp.diff(cosine, u) ** 2 - (1 - cosine ** 2) / ell ** 2).subs(u, Q(2, 3))], system,
                  "a = cos(ct/ell) does not solve the disc's equations")
            continue
        L, b = Q(7, 5), Q(9, 4)
        velocity = -sp.sqrt(friedmann(a0, b, L))
        # The derivative of the Friedmann equation along the time.
        pull = -(1 / L ** 2 - 1 / b ** 2) / (2 * a0 ** 2) - a0 / L ** 2
        at = [(second, pull), (first, velocity), (a, a0), (p["b"], b), (p["L"], L), (s[0], Q(1, 3))]
        if system == HYP:
            at += [(p["r_0"], Q(13, 4)), (s[1], Q(11, 4)), (s[2], Q(4, 5)), (s[3], Q(2, 7))]
        elif system == SLAB:
            at += [(p["chi_0"], Q(3, 2)), (s[1], Q(-4, 5)), (s[2], Q(4, 5)), (s[3], Q(2, 7))]
        else:
            at += [(p["chi_0"], Q(3, 2)), (s[1], Q(-4, 5)), (s[2], Q(1, 5)), (s[3], Q(2, 7))]
        G, g = mixed_einstein(chart, at)
        Lambda, density = -3 / L ** 2, 3 * (1 / L ** 2 - 1 / b ** 2) / a0 ** 3
        where = f"{system} at a = {a0}"
        want = sp.diag(-density - Lambda, -Lambda, -Lambda, -Lambda)
        small(G - want, where, "the cloud is not dust with its cosmological constant")
        if system == SLAB:
            # Smith and Mann's chart pulled back along r = b cosh(chi/b).
            chi, theta = Q(-4, 5), Q(4, 5)
            r = b * sp.cosh(chi / b)
            there = sp.diag(-1, a0 ** 2 / (r ** 2 / b ** 2 - 1), a0 ** 2 * r ** 2, a0 ** 2 * r ** 2 * sp.sinh(theta) ** 2)
            J = sp.diag(1, sp.sinh(chi / b), 1, 1)
            small(g - J.T * there * J, where, "the slab is not Smith and Mann's chart pulled back")
    if system != DISC:
        # Released from rest, and falling: the Friedmann function vanishes at a = 1 with a negative slope.
        x, L, b = sp.Symbol("x", positive=True), Q(7, 5), Q(9, 4)
        small([friedmann(1, b, L)], system, "the cloud is not at rest at a = 1")
        if not sp.diff(friedmann(x, b, L), x).subs(x, 1) < 0:
            raise AssertionError(f"{METRIC}/{system}: the cloud released at a = 1 does not fall")


def exterior(vm, chart, system, s, p, n):
    three = system in (DISC_OUT, DISC_IN)
    for a0 in SCALES:
        if three:
            ell, k, r0 = Q(7, 5), Q(2, 9), Q(9, 7)
            M = (k + 1 / ell ** 2) * r0 ** 2 - 1
            R = r0 * a0
            f = R ** 2 / ell ** 2 - M
            values = [(p["ell"], ell), (p["M"], M)]
            Lambda, energy = -1 / ell ** 2, 1 - k * r0 ** 2
            speed = r0 ** 2 * (1 - a0 ** 2) / ell ** 2
            surface = [(2, R ** 2)]
        else:
            L, b = Q(7, 5), Q(9, 4)
            excess = 1 / L ** 2 - 1 / b ** 2
            Lambda = -3 / L ** 2
            if system == HYP_OUT:
                r0, theta = Q(13, 4), Q(4, 5)
                R = r0 * a0
                mu = excess * r0 ** 3
                f, energy = -1 - mu / R + R ** 2 / L ** 2, r0 ** 2 / b ** 2 - 1
                values = [(p["mu"], mu), (p["L"], L)]
                surface = [(2, a0 ** 2 * r0 ** 2), (3, a0 ** 2 * r0 ** 2 * sp.sinh(theta) ** 2)]
            else:
                chi0 = Q(3, 2)
                r0 = b * sp.exp(chi0 / b)
                R = r0 * a0
                mu = excess * r0 ** 3
                f, energy = -mu / R + R ** 2 / L ** 2, r0 ** 2 / b ** 2
                values = [(p["mu"], mu), (p["L"], L), (p["b"], b)]
                surface = [(2, a0 ** 2 * sp.exp(2 * chi0 / b)), (3, a0 ** 2 * sp.exp(2 * chi0 / b))]
            speed = r0 ** 2 * friedmann(a0, b, L)
        where = f"{system} at a = {a0}"
        at = values + [(s[0], Q(1, 3)), (s[1], R), (s[2], Q(4, 5))] + ([] if three else [(s[3], Q(2, 7))])
        G, g = mixed_einstein(chart, at)
        small(G + Lambda * sp.eye(n), where, "the exterior is not empty with its cosmological constant")
        # The metric of the surface, the interior's at its edge.
        small([g[i, i] - value for i, value in surface], where, "the surface's metric is not the interior's")
        # The surface is the exterior's radial geodesic of the interior's energy.
        small([speed + f - energy], where, "the surface is not the exterior's geodesic of the interior's energy")
        if system == DISC_IN:
            small([g[0, 0] + f, g[0, 1] - 1, g[1, 1]], where, "the ingoing chart is not -f dv^2 + 2 dv dR")
            # The static chart pulled back along dT = dv - dR/f.
            J = sp.Matrix([[1, -1 / f, 0], [0, 1, 0], [0, 0, 1]])
            small(g - J.T * sp.diag(-f, 1 / f, R ** 2) * J, where, "the ingoing chart is not the static one pulled back")
        else:
            small([g[0, 0] + f, g[1, 1] - 1 / f], where, "the static chart is not -f dT^2 + dR^2/f")
        # The published chart it is.
        source = {DISC_OUT: ("btz", "stationary"), DISC_IN: ("btz", "eddington_finkelstein_ingoing"),
                  HYP_OUT: ("topological_black_hole", "hyperbolic"), FLAT_OUT: ("topological_black_hole", "static")}
        reader, there, coords = published(vm, *source[system])
        q = reader.parameters
        if three:
            same = [(q["ell"], ell), (q["M"], M), (q["J"], 0)] + list(zip(coords, (Q(1, 3), R, Q(4, 5))))
            small(g - there.subs(same), where, f"not the published chart of {source[system][0]} without rotation")
        elif system == HYP_OUT:
            same = [(q["mu"], mu), (q["L"], L)] + list(zip(coords, (Q(1, 3), R, Q(4, 5), Q(2, 7))))
            small(g - there.subs(same), where, "not the published hyperbolic chart of topological_black_hole")
        else:
            # x = b rho cos(phi), y = b rho sin(phi) on the published chart's plane at k = 0.
            rho, phi = Q(3, 7), Q(2, 7)
            same = [(q["mu"], mu), (q["L"], L), (q["k"], 0)] + list(zip(coords, (Q(1, 3), R, rho, phi)))
            J = sp.eye(4)
            J[2, 2], J[2, 3] = b * sp.cos(phi), -b * rho * sp.sin(phi)
            J[3, 2], J[3, 3] = b * sp.sin(phi), b * rho * sp.cos(phi)
            small(J.T * g * J - there.subs(same), where, "not the published static chart of topological_black_hole at k = 0")
