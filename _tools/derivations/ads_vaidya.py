"""The charts of a shell of light falling in anti-de Sitter space, Vaidya's metric with a negative
cosmological constant, in three, four and five dimensions. ads_vaidya.md beside this file records
each chart's source.

Every chart is written in the advanced time v, a time, and printed in the chart x^0 = cv, with the
mass function m(v) free and its derivative along x^0 printed as a prime. In D = d + 1 dimensions the
ingoing charts are

    ds^2 = -f c^2 dv^2 + 2 c dv dr + r^2 dOmega^2,   f = k + r^2/l^2 - mu(v)/r^(d - 2),

with k = 1 for a round front and 0 for a flat one, and the Poincare charts are the flat ones in the
depth z = l^2/r,

    ds^2 = (l^2/z^2)(-(1 - m z^d) c^2 dv^2 - 2 c dv dz + dx^2),   m = mu/l^(2d - 2).

The one stress component is (G + Lambda g)_vv = ((d - 1)/2) mu'/r^(d - 1), with Lambda =
-d(d - 1)/2l^2: the energy density of light falling inward along the rays of constant v.
"""
import json

import sympy as sp

METRIC = "ads_vaidya"
CHARTS = ("ingoing_three", "poincare_three", "ingoing_four", "flat_front", "poincare_four", "ingoing_five")

DIMENSIONS = {
    (METRIC, "ingoing_three"): {"v": "T", "r": "L", "\\theta": "1", "\\ell": "L", "m": "1"},
    (METRIC, "poincare_three"): {"v": "T", "z": "L", "x": "L", "\\ell": "L", "m": "L**-2"},
    (METRIC, "ingoing_four"): {"v": "T", "r": "L", "\\theta": "1", "\\phi": "1",
                               "G": "L**3/(M*T**2)", "\\ell": "L", "m": "M"},
    (METRIC, "flat_front"): {"v": "T", "r": "L", "x": "L", "y": "L", "\\ell": "L", "m": "L"},
    (METRIC, "poincare_four"): {"v": "T", "z": "L", "x": "L", "y": "L", "\\ell": "L", "m": "L**-3"},
    (METRIC, "ingoing_five"): {"v": "T", "r": "L", "\\chi": "1", "\\theta": "1", "\\phi": "1",
                               "\\ell": "L", "m": "L**2"},
}

REALS = " \\in (-\\infty, \\infty)"
SPHERE = "\\left(d\\theta^2 + \\sin^2\\theta\\,d\\phi^2\\right)"
# The lapse of each ingoing chart, as its line element writes it: (coordinates, parameters, f, what
# f multiplies beside the radius, the tangential part, the domains after v and r, the scalar K).
INGOING = {
    "ingoing_three": (["v", "r", "\\theta"], ["\\ell", "m = m(v)"], "\\dfrac{r^2}{\\ell^2} - m",
                      " + r^2d\\theta^2", ["\\theta \\in [0, 2\\pi)"],
                      ["r = \\ell\\sqrt{m} \\;\\text{(the apparent horizon, where } m > 0 \\text{)}"],
                      "\\dfrac{12}{\\ell^4}"),
    "ingoing_four": (["v", "r", "\\theta", "\\phi"], ["G", "\\ell", "m = m(v)"],
                     "1 + \\dfrac{r^2}{\\ell^2} - \\dfrac{2Gm}{c^2r}", " + r^2" + SPHERE,
                     ["\\theta \\in [0, \\pi]", "\\phi \\in [0, 2\\pi)"],
                     ["1 + \\dfrac{r^2}{\\ell^2} = \\dfrac{2Gm}{c^2r} \\;\\text{(the apparent horizon)}"],
                     "\\dfrac{24}{\\ell^4} + \\dfrac{48G^2m^2}{c^4r^6}"),
    "flat_front": (["v", "r", "x", "y"], ["\\ell", "m = m(v)"], "\\dfrac{r^2}{\\ell^2} - \\dfrac{m}{r}",
                   " + \\dfrac{r^2}{\\ell^2}\\left(dx^2 + dy^2\\right)", ["x" + REALS, "y" + REALS],
                   ["r = \\left(m\\ell^2\\right)^{1/3} \\;\\text{(the apparent horizon)}"],
                   "\\dfrac{24}{\\ell^4} + \\dfrac{12m^2}{r^6}"),
    "ingoing_five": (["v", "r", "\\chi", "\\theta", "\\phi"], ["\\ell", "m = m(v)"],
                     "1 + \\dfrac{r^2}{\\ell^2} - \\dfrac{m}{r^2}",
                     " + r^2\\left(d\\chi^2 + \\sin^2\\chi" + SPHERE + "\\right)",
                     ["\\chi \\in [0, \\pi]", "\\theta \\in [0, \\pi]", "\\phi \\in [0, 2\\pi)"],
                     ["1 + \\dfrac{r^2}{\\ell^2} = \\dfrac{m}{r^2} \\;\\text{(the apparent horizon)}"],
                     "\\dfrac{40}{\\ell^4} + \\dfrac{72m^2}{r^8}"),
}
# The Poincare charts: (coordinates, d, the flat part, the domains after v and z, the scalar K).
POINCARE = {
    "poincare_three": (["v", "z", "x"], 2, "dx^2", ["x" + REALS],
                       ["z = m^{-1/2} \\;\\text{(the apparent horizon, where } m > 0 \\text{)}"],
                       "\\dfrac{12}{\\ell^4}"),
    "poincare_four": (["v", "z", "x", "y"], 3, "dx^2 + dy^2", ["x" + REALS, "y" + REALS],
                      ["z = m^{-1/3} \\;\\text{(the apparent horizon)}"],
                      "\\dfrac{24}{\\ell^4} + \\dfrac{12m^2z^6}{\\ell^4}"),
}
NAMES = {"ingoing_three": "Ingoing, Three Dimensions", "poincare_three": "Poincaré, Three Dimensions",
         "ingoing_four": "Ingoing, Four Dimensions", "flat_front": "Flat Front, Four Dimensions",
         "poincare_four": "Poincaré, Four Dimensions", "ingoing_five": "Ingoing, Five Dimensions"}


def chart_spec(system):
    """One chart: the line element with c, the same in x^0 = cv, the components written by hand as
    the line element has them, the Kretschmann scalar as anti-de Sitter's part plus the hole's, and
    the check of the chart against its source and the live pages it reduces to."""
    import verify_metrics as vm
    if system in INGOING:
        coords, parameters, f, rest, angles, marks, scalar = INGOING[system]
        wrapped = "\\left(" + f + "\\right)"
        line = "ds^2 = -" + wrapped + "@C2@dv^2 + 2@C@\\,dv\\,dr" + rest
        domains = ["v" + REALS, "r \\in (0, \\infty)"] + angles + marks
        metric = {("v", "v"): "-" + wrapped, ("v", "r"): "1", ("r", "v"): "1"}
        inverse = {("v", "r"): "1", ("r", "v"): "1", ("r", "r"): f}
    else:
        coords, d, flat, lengths, marks, scalar = POINCARE[system]
        parameters = ["\\ell", "m = m(v)"]
        f = f"1 - m\\,z^{d}"
        line = ("ds^2 = \\dfrac{\\ell^2}{z^2}\\left(-\\left(" + f + "\\right)@C2@dv^2 - 2@C@\\,dv\\,dz + "
                + flat + "\\right)")
        domains = ["v" + REALS, "z \\in (0, \\infty)"] + lengths + marks
        metric = {("v", "v"): "-\\dfrac{\\ell^2\\left(" + f + "\\right)}{z^2}",
                  ("v", "z"): "-\\dfrac{\\ell^2}{z^2}", ("z", "v"): "-\\dfrac{\\ell^2}{z^2}"}
        inverse = {("v", "z"): "-\\dfrac{z^2}{\\ell^2}", ("z", "v"): "-\\dfrac{z^2}{\\ell^2}",
                   ("z", "z"): "\\dfrac{z^2\\left(" + f + "\\right)}{\\ell^2}"}
    import chart_printer as cp
    probe = vm.Reader(coords, parameters, ())
    radius, ell, m = probe.symbol[coords[1]], probe.parameters["ell"], probe.parameters["m"]
    rate = sp.Derivative(m, probe.symbol["v"])
    newton = [probe.parameters["G"]] if "G" in probe.parameters else []
    # A product is written as the line element writes it, the mass before the lengths, and the rate
    # of change of the mass stands apart from what the static hole has already.
    printer = {"primed": ["m"], "lead": [rate, radius, ell, m], "factors": newton + [m, probe.c, ell, radius, rate],
               "flip": False, "collect": lambda poly, pr: cp.collect_by(poly, [rate], pr)}
    return {
        "metric_id": METRIC,
        "system": {"id": system, "name": NAMES[system], "coords": coords, "domains": domains,
                   "parameters": parameters, "line_element": line.replace("@C2@", "c^2").replace("@C@", "c")},
        "chart_line_element": line.replace("@C2@", "").replace("@C@", ""),
        "printer": printer,
        "components": {"metric_components": metric, "inverse_metric_components": inverse},
        "kretschmann": scalar,
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


def check_chart(chart, system):
    """Each chart against its source's one stress component and against the live pages it becomes.

    The Einstein tensor with the cosmological term, (G + Lambda g)_ab with Lambda = -d(d - 1)/2l^2,
    has the one component vv, ((d - 1)/2) mu'/r^(d - 1) with mu what multiplies r^(2 - d) in f:
    Husain's rho/2r in three dimensions (his equations (2), (4) and (6)), Lemos's 8 pi T_vv =
    2m'/r^2 and q m'/r^2 in four ((3) and (7)), and Hubeny, Liu and Rangamani's q/r^3 with m' =
    2q/3 in five (appendix F). With the mass constant the ingoing charts are the ingoing charts of
    btz without spin, of schwarzschild_ads and of topological_black_hole at k = 0, and the static
    chart of brane_wall_universe at k = 1 carried to advanced time; without the cosmological term
    the chart of four dimensions is Vaidya's and that of five Tangherlini's; and each Poincare
    chart is the flat chart of its dimension along z = l^2/r."""
    import verify_metrics as vm
    geo, g, n = chart.geo, chart.geo.g, len(chart.symbols)
    d = n - 1
    v, radius = chart.symbols[0], chart.symbols[1]
    ell, m = chart.reader.parameters["ell"], chart.reader.parameters["m"]
    rate = sp.diff(m, v)
    c = chart.reader.c
    name = f"{METRIC}/{system}"

    def same(left, right, what):
        for i in range(left.rows):
            for j in range(left.cols):
                if vm.norm(left[i, j] - right[i, j]) != 0:
                    raise AssertionError(f"{name}: {what} in slot {chart.coords_tex[i]}{chart.coords_tex[j]}")

    # The flux of light, the one component of G + Lambda g.
    if system in INGOING:
        mu = {"ingoing_three": m, "ingoing_four": 2 * chart.reader.parameters.get("G", 0) * m / c ** 2,
              "flat_front": m, "ingoing_five": m}[system]
        flux = sp.Rational(d - 1, 2) * sp.diff(mu, v) / radius ** (d - 1)
    else:
        flux = sp.Rational(d - 1, 2) * rate * radius ** (d - 1)
    cosmological = -sp.Rational(d * (d - 1), 2) / ell ** 2
    einstein = geo.einstein_ll()
    for a in range(n):
        for b in range(n):
            want = flux if a == b == 0 else 0
            if vm.norm(vm._at(einstein, (a, b)) + cosmological * g[a, b] - want) != 0:
                raise AssertionError(f"{name}: G + Lambda g is not the flux of light in slot "
                                     f"{chart.coords_tex[a]}{chart.coords_tex[b]}")
    mass = sp.Symbol("mass", positive=True)
    still = g.subs(m, mass)

    if system in POINCARE:
        # The flat chart of the same dimension along r = l^2/z, x = x, with mu = m l^(2d - 2).
        flat = chart_spec("ingoing_three" if d == 2 else "flat_front")
        import chart_printer as cp
        other = cp.Chart(flat["system"]["coords"], flat["system"]["parameters"], flat["chart_line_element"])
        mu = sp.Function("mu")(v)
        at = {other.symbols[1]: ell ** 2 / radius, other.reader.parameters["m"]: mu,
              other.reader.parameters["ell"]: ell, other.symbols[0]: v}
        J = sp.eye(n)
        J[1, 1] = -ell ** 2 / radius ** 2
        if d == 2:
            # theta = x/l: the circle of the ingoing chart unrolled.
            J[2, 2] = 1 / ell
        else:
            at.update(dict(zip(other.symbols[2:], chart.symbols[2:])))
        pulled = (J.T * other.geo.g.subs(at) * J).subs(mu, m * ell ** (2 * d - 2))
        same(pulled, g, "the flat chart along z = l^2/r misses the Poincare chart")
        return

    if system == "ingoing_three":
        reader, hole, coords = published_metric(vm, "btz", "eddington_finkelstein_ingoing")
        swap = dict(zip((reader.symbol[x] for x in coords), chart.symbols))
        swap.update({reader.parameters["ell"]: ell, reader.parameters["M"]: mass, reader.parameters["J"]: 0})
        same(hole.subs(swap), still, "with a constant mass the chart misses the published btz without spin")
        # Anti-de Sitter space is m = -1: a regular centre, where the circle of radius r has
        # circumference 2 pi r to first order, g_vv = -1 at r = 0, and constant curvature.
        ricci = geo.ricci_ll()
        for a in range(n):
            for b in range(n):
                if vm.norm((vm._at(ricci, (a, b)) + 2 * g[a, b] / ell ** 2).subs(rate, 0).subs(m, -1)) != 0:
                    raise AssertionError(f"{name}: at m = -1 the chart is not anti-de Sitter space")
        if vm.norm(g[0, 0].subs(m, -1).subs(radius, 0) + 1) != 0:
            raise AssertionError(f"{name}: at m = -1 the centre is not regular")
    elif system == "ingoing_four":
        G = chart.reader.parameters["G"]
        reader, hole, coords = published_metric(vm, "schwarzschild_ads", "eddington_finkelstein_ingoing")
        swap = dict(zip((reader.symbol[x] for x in coords), chart.symbols))
        swap.update({reader.parameters["L"]: ell, reader.parameters["r_s"]: 2 * G * mass / c ** 2})
        same(hole.subs(swap), still, "with a constant mass the chart misses the published schwarzschild_ads")
        reader, star, coords = published_metric(vm, "vaidya", "eddington_finkelstein_ingoing")
        swap = dict(zip((reader.symbol[x] for x in coords), chart.symbols))
        swap.update({reader.parameters["G"]: G})
        there = star.subs(swap).subs(reader.parameters["m"], mass)
        bare = still.copy()
        bare[0, 0] += radius ** 2 / ell ** 2
        same(there, bare, "without the cosmological term the chart misses the published vaidya")
    elif system == "flat_front":
        # The published chart at k = 0 has the plane in polar coordinates, x = l rho cos(phi).
        reader, hole, coords = published_metric(vm, "topological_black_hole", "eddington_finkelstein_ingoing")
        rho, phi = reader.symbol["\\rho"], reader.symbol["\\phi"]
        swap = {reader.symbol["v"]: v, reader.symbol["r"]: radius, reader.parameters["L"]: ell,
                reader.parameters["mu"]: mass, reader.parameters["k"]: 0}
        J = sp.eye(n)
        J[2, 2], J[2, 3] = ell * sp.cos(phi), -ell * rho * sp.sin(phi)
        J[3, 2], J[3, 3] = ell * sp.sin(phi), ell * rho * sp.cos(phi)
        same(J.T * still * J, hole.subs(swap),
             "with a constant mass the chart misses the published topological_black_hole at k = 0")
    elif system == "ingoing_five":
        # The published static chart at k = 1 along t = v - r_*, with its comoving radius sin(chi).
        reader, hole, coords = published_metric(vm, "brane_wall_universe", "static")
        chi = chart.symbols[2]
        swap = {reader.symbol["r"]: radius, reader.symbol["\\chi"]: sp.sin(chi),
                reader.symbol["\\theta"]: chart.symbols[3], reader.symbol["\\phi"]: chart.symbols[4],
                reader.parameters["ell"]: ell, reader.parameters["mu"]: mass, reader.parameters["k"]: 1}
        old = hole.subs(swap)
        J = sp.eye(n)
        J[0, 1] = 1 / old[0, 0]
        J[2, 2] = sp.cos(chi)
        same(J.T * old * J, still, "with a constant mass the chart misses the published brane_wall_universe")
        reader, hole, coords = published_metric(vm, "tangherlini", "eddington_finkelstein_ingoing")
        swap = dict(zip((reader.symbol[x] for x in coords), chart.symbols))
        there = hole.subs(swap).subs(reader.parameters["r_h"], sp.sqrt(mass))
        bare = still.copy()
        bare[0, 0] += radius ** 2 / ell ** 2
        same(there, bare, "without the cosmological term the chart misses the published tangherlini")
