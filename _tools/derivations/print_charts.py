#!/usr/bin/env python3
"""Compute and write the coordinate systems whose mathematics is printed by machine: the
charts of tov, malament_hogarth, mixmaster, lentz, einstein_static, btz, c_metric,
schwarzschild_de_sitter, milne and einstein_rosen_waves, and Godel's cylindrical chart.

    /tmp/mfs-venv/bin/python _tools/derivations/print_charts.py [--metric <id>]...
    /tmp/mfs-venv/bin/python _tools/derivations/verify_metrics.py --system <id>/<system>

Each chart below names its coordinates, its parameters, the line element it publishes, the
same line element in the chart x^0 = ct that the components are printed in, and how its values
are to be grouped for reading, and any component written by hand. The script computes every tensor from that line element with
the checker's own Geometry, prints each component through chart_printer.py, reads every
printed value back and compares it with what it was printed from, and writes the result into
the metric file. It writes only the mathematics: the entry's prose, its conventions, the
chart's own among them, and the description of each parameter stay in the metric file and are
carried over untouched, and a parameter the file does not yet describe, or a chart it does not
yet give a convention, stops the script rather than being written without one. The Ricci scalar and the Kretschmann scalar are written by hand where a structured form
reads better than an expanded one, and each of those is checked against sympy here too.

The derivations these charts rest on, and the reason each was chosen, are in tov.md,
malament_hogarth.md, mixmaster.md, lentz.md, godel.md, btz.md and schwarzschild_de_sitter.md
beside this file.
"""
import argparse
import json
import sys
import time
from pathlib import Path

import sympy as sp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import chart_printer as cp  # noqa: E402
import verify_metrics as vm  # noqa: E402

METRICS = vm.METRICS_DIR
SYSTEM_ORDER = ["id", "name", "coords", "domains", "parameters", "convention", "line_element", "metric_components",
                "inverse_metric_components", "christoffel", "riemann", "ricci_tensor", "ricci_scalar",
                "kretschmann", "einstein_tensor", "weyl_tensor", "geodesics"]


# -- Tolman-Oppenheimer-Volkoff --------------------------------------------------------

def tov():
    coords = ["t", "r", "\\theta", "\\phi"]
    parameters = ["\\Phi = \\Phi(r)", "m = m(r)"]
    probe = vm.Reader(coords, parameters, ())
    r, Phi, m = probe.symbol["r"], probe.parameters["Phi"], probe.parameters["m"]
    dPhi, d2Phi = sp.Derivative(Phi, r), sp.Derivative(Phi, (r, 2))
    Q = sp.Symbol("Q")

    def collect(poly, printer):
        # Every Riemann component carries the combination (d_r Phi)^2 + d_r^2 Phi = e^{-Phi}(e^Phi)''.
        return cp.collect_by(poly, [Q, dPhi], printer, {d2Phi: Q - dPhi ** 2})

    return {
        "metric_id": "tov",
        "system": {"id": "spherical", "name": "Spherical", "coords": coords,
                   "domains": ["t \\in (-\\infty, \\infty)", "r \\in [0, \\infty)", "\\theta \\in [0, \\pi]",
                               "\\phi \\in [0, 2\\pi)", "2m < r"],
                   "parameters": parameters,
                   "line_element": "ds^2 = -e^{2\\Phi}c^2dt^2 + \\dfrac{dr^2}{1 - \\dfrac{2m}{r}} + r^2\\left(d\\theta^2 + \\sin^2\\theta\\,d\\phi^2\\right)"},
        "chart_line_element": "ds^2 = -e^{2\\Phi}dt^2 + \\dfrac{dr^2}{1 - \\dfrac{2m}{r}} + r^2\\left(d\\theta^2 + \\sin^2\\theta\\,d\\phi^2\\right)",
        "printer": {"lead": [Q, dPhi, r, m, sp.Derivative(m, r)], "collect": collect,
                    "overrides": {Q: "\\left(\\left(\\partial_r\\Phi\\right)^2 + \\partial_r^2\\Phi\\right)"}},
        # The orthonormal frame of a static observer: 4 A^2 + 8 B^2 + 8 C^2 + 4 D^2 over the
        # four independent frame components t r t r, t theta t theta, r theta r theta, theta phi theta phi.
        "kretschmann": (
            "4\\left(\\dfrac{\\left(r - 2m\\right)\\left(\\left(\\partial_r\\Phi\\right)^2 + \\partial_r^2\\Phi\\right)}{r}"
            " - \\dfrac{\\left(r\\,\\partial_r m - m\\right)\\partial_r\\Phi}{r^2}\\right)^2"
            " + \\dfrac{8\\left(r - 2m\\right)^2\\left(\\partial_r\\Phi\\right)^2}{r^4}"
            " + \\dfrac{8\\left(r\\,\\partial_r m - m\\right)^2}{r^6} + \\dfrac{16m^2}{r^6}"),
    }


# -- Malament-Hogarth ------------------------------------------------------------------

def malament_hogarth():
    coords = ["t", "x", "y", "z"]
    parameters = ["\\Omega = \\Omega(t,x,y,z)"]
    probe = vm.Reader(coords, parameters, ())
    Om = probe.parameters["Omega"]

    def collect(poly, printer):
        return cp.collect_by(poly, [Om], printer)

    box = "\\partial_x^2\\Omega + \\partial_y^2\\Omega + \\partial_z^2\\Omega - \\partial_t^2\\Omega"
    d = {c: f"\\partial_{c}\\Omega" for c in "txyz"}

    def A(a):
        # 2 Omega^2 times the diagonal Schouten component: 4(d_a W)^2 - 2 W d_a^2 W - eta_aa (dW)^2.
        text = ""
        for c in "txyz":
            k = (4 if c == a else 0) - (-1 if a == "t" else 1) * (-1 if c == "t" else 1)
            body = ("" if abs(k) == 1 else str(abs(k))) + f"\\left({d[c]}\\right)^2"
            text += (("-" if k < 0 else "") + body) if not text else ((" - " if k < 0 else " + ") + body)
        return "\\left(" + text + f" - 2\\Omega\\,\\partial_{a}^2\\Omega" + "\\right)^2"

    def B(a, b):
        return f"\\left(2{d[a]}\\,{d[b]} - \\Omega\\,\\partial_{a}\\partial_{b}\\Omega\\right)^2"

    return {
        "metric_id": "malament_hogarth",
        "system": {"id": "cartesian", "name": "Cartesian", "coords": coords,
                   "domains": ["t \\in (-\\infty, \\infty)", "x \\in (-\\infty, \\infty)",
                               "y \\in (-\\infty, \\infty)", "z \\in (-\\infty, \\infty)",
                               "(t, x, y, z) \\neq (0, 0, 0, 0)"],
                   "parameters": parameters,
                   "line_element": "ds^2 = \\Omega^2\\left(-c^2dt^2 + dx^2 + dy^2 + dz^2\\right)"},
        "chart_line_element": "ds^2 = \\Omega^2\\left(-dt^2 + dx^2 + dy^2 + dz^2\\right)",
        "printer": {"lead": [Om], "collect": collect},
        "ricci_scalar": "-\\dfrac{6\\left(" + box + "\\right)}{\\Omega^3}",
        # K = 8 P_ab P^ab + 4 (P^a_a)^2 for a conformally flat metric, P the Schouten tensor.
        "kretschmann": (
            "\\dfrac{2\\left(" + " + ".join(A(a) for a in "txyz") + "\\right)"
            " + 16\\left(" + " + ".join(B(a, b) for a, b in (("x", "y"), ("x", "z"), ("y", "z")))
            + " - " + " - ".join(B("t", b) for b in "xyz") + "\\right)"
            " + 4\\Omega^2\\left(" + box + "\\right)^2}{\\Omega^8}"),
    }


# -- Mixmaster -------------------------------------------------------------------------

def mixmaster():
    coords = ["t", "\\psi", "\\theta", "\\phi"]
    parameters = ["a_1 = a_1(t)", "a_2 = a_2(t)", "a_3 = a_3(t)"]
    probe = vm.Reader(coords, parameters, ())
    psi, th = probe.symbol["\\psi"], probe.symbol["\\theta"]
    a = [probe.parameters[f"a_{i}"] for i in (1, 2, 3)]
    outer, inner = (sp.sin(th), sp.cos(th)), (sp.sin(psi), sp.cos(psi))

    def collect(poly, printer):
        return cp.collect_nested(poly, outer, inner, printer)

    forms = ("\\left(\\cos\\psi\\,d\\theta + \\sin\\psi\\sin\\theta\\,d\\phi\\right)^2",
             "\\left(\\sin\\psi\\,d\\theta - \\cos\\psi\\sin\\theta\\,d\\phi\\right)^2",
             "\\left(d\\psi + \\cos\\theta\\,d\\phi\\right)^2")
    spatial = f"a_1^2{forms[0]} + a_2^2{forms[1]} + a_3^2{forms[2]}"
    return {
        "metric_id": "mixmaster",
        "system": {"id": "euler_angles", "name": "Euler Angles", "coords": coords,
                   "domains": ["t \\in (0, \\infty)", "\\psi \\in [0, 4\\pi)", "\\theta \\in [0, \\pi]",
                               "\\phi \\in [0, 2\\pi)"],
                   "parameters": parameters,
                   "line_element": "ds^2 = -c^2dt^2 + " + spatial},
        "chart_line_element": "ds^2 = -dt^2 + " + spatial,
        "printer": {"primed": ["a_1", "a_2", "a_3"], "lead": a + [inner[0], inner[1], outer[0], outer[1]],
                    "collect": collect},
        "pretty": lambda v: cp.merge_squares(sp.factor(v)),
        # Bianchi I's scalar and the curvature of the closed slice, which is 6/(2a)^2 for a
        # round three sphere of radius 2a.
        "ricci_scalar": "2\\left(\\dfrac{a_1''}{a_1} + \\dfrac{a_2''}{a_2} + \\dfrac{a_3''}{a_3} + \\dfrac{a_1'\\,a_2'}{a_1\\,a_2} + \\dfrac{a_1'\\,a_3'}{a_1\\,a_3} + \\dfrac{a_2'\\,a_3'}{a_2\\,a_3}\\right) - \\dfrac{a_1^4 + a_2^4 + a_3^4 - 2a_1^2\\,a_2^2 - 2a_1^2\\,a_3^2 - 2a_2^2\\,a_3^2}{2a_1^2\\,a_2^2\\,a_3^2}",
        # The orthonormal frame c dt, a_i omega^i: 4 times the squares of the six diagonal
        # frame components R_0i0i and R_ijij, less 8 times the squares of the three that
        # pair (0i) with (jk), which enter with one time index on each side.
        "kretschmann": "4\\left(\\left(\\dfrac{a_1''}{a_1}\\right)^2 + \\left(\\dfrac{a_2''}{a_2}\\right)^2 + \\left(\\dfrac{a_3''}{a_3}\\right)^2 + \\left(\\dfrac{4a_1\\,a_2\\,a_3^2\\,a_1'\\,a_2' + a_1^4 - 2a_1^2\\,a_2^2 + 2a_1^2\\,a_3^2 + a_2^4 + 2a_2^2\\,a_3^2 - 3a_3^4}{4a_1^2\\,a_2^2\\,a_3^2}\\right)^2 + \\left(\\dfrac{4a_1\\,a_2^2\\,a_3\\,a_1'\\,a_3' + a_1^4 + 2a_1^2\\,a_2^2 - 2a_1^2\\,a_3^2 - 3a_2^4 + 2a_2^2\\,a_3^2 + a_3^4}{4a_1^2\\,a_2^2\\,a_3^2}\\right)^2 + \\left(\\dfrac{4a_1^2\\,a_2\\,a_3\\,a_2'\\,a_3' - 3a_1^4 + 2a_1^2\\,a_2^2 + 2a_1^2\\,a_3^2 + a_2^4 - 2a_2^2\\,a_3^2 + a_3^4}{4a_1^2\\,a_2^2\\,a_3^2}\\right)^2\\right) - 8\\left(\\left(\\dfrac{2a_1\\,a_2\\,a_3\\,a_1' - a_3\\left(a_1^2 + a_2^2 - a_3^2\\right)a_2' - a_2\\left(a_1^2 - a_2^2 + a_3^2\\right)a_3'}{2a_1\\,a_2^2\\,a_3^2}\\right)^2 + \\left(\\dfrac{a_1\\left(a_1^2 - a_2^2 - a_3^2\\right)a_3' + 2a_1\\,a_2\\,a_3\\,a_2' - a_3\\left(a_1^2 + a_2^2 - a_3^2\\right)a_1'}{2a_1^2\\,a_2\\,a_3^2}\\right)^2 + \\left(\\dfrac{2a_1\\,a_2\\,a_3\\,a_3' + a_1\\left(a_1^2 - a_2^2 - a_3^2\\right)a_2' - a_2\\left(a_1^2 - a_2^2 + a_3^2\\right)a_1'}{2a_1^2\\,a_2^2\\,a_3}\\right)^2\\right)",
    }


# -- Lentz -----------------------------------------------------------------------------

def lentz():
    coords = ["t", "x", "y", "z"]
    parameters = ["\\phi = \\phi(t,x,y,z)"]
    probe = vm.Reader(coords, parameters, ())
    X = [probe.symbol[c] for c in coords]
    phi = probe.parameters["phi"]
    N = [sp.Derivative(phi, v) for v in X[1:]]
    printer = cp.Printer(X, lead=N)

    def D(*variables):
        return printer.atom(sp.Derivative(phi, *variables))

    def E(i, j):
        # The tidal block of the Eulerian observers, minus n^a n^c R_{a i c j}: the rate of the
        # Hessian along their flow, d_t + N^k d_k, plus the square of the Hessian.
        rate = [D(X[0], X[i], X[j])] + [f"{D(X[k])}\\,{D(X[k], X[i], X[j])}" for k in (1, 2, 3)]
        square = [printer.term(sp.Mul(sp.Derivative(phi, X[i], X[k]), sp.Derivative(phi, X[k], X[j])))
                  for k in (1, 2, 3)]
        return "\\left(" + " + ".join(rate + square) + "\\right)^2"

    def C(i, j):
        # A cofactor of the Hessian, the spatial Riemann component through the Gauss equation.
        p, q = [k for k in (1, 2, 3) if k != i]
        u, v = [k for k in (1, 2, 3) if k != j]
        first = printer.term(sp.Mul(sp.Derivative(phi, X[p], X[u]), sp.Derivative(phi, X[q], X[v])))
        second = printer.term(sp.Mul(sp.Derivative(phi, X[p], X[v]), sp.Derivative(phi, X[q], X[u])))
        return "\\left(" + first + " - " + second + "\\right)^2"

    diagonal, off = [(1, 1), (2, 2), (3, 3)], [(1, 2), (1, 3), (2, 3)]
    shift = "\\left(d{0} - \\partial_{0}\\phi\\,{1}dt\\right)^2"

    def line(c):
        return " + ".join(shift.format(v, c) for v in "xyz")

    return {
        "metric_id": "lentz",
        "system": {"id": "cartesian", "name": "Cartesian", "coords": coords,
                   "domains": ["t \\in (-\\infty, \\infty)", "x \\in (-\\infty, \\infty)",
                               "y \\in (-\\infty, \\infty)", "z \\in (-\\infty, \\infty)"],
                   "parameters": parameters,
                   "line_element": "ds^2 = -c^2dt^2 + " + line("c\\,")},
        "chart_line_element": "ds^2 = -dt^2 + " + line(""),
        "printer": {"lead": N, "collect": lambda poly, printer: cp.collect_by(poly, N, printer)},
        # The flow is a gradient, so the Codazzi block vanishes and K is a sum of squares:
        # four times the spatial block (the Hessian's cofactors) and four times the tidal block.
        "kretschmann": (
            "4\\left(" + " + ".join(C(i, j) for i, j in diagonal) + " + 2" + " + 2".join(C(i, j) for i, j in off) + "\\right)"
            " + 4\\left(" + " + ".join(E(i, j) for i, j in diagonal) + " + 2" + " + 2".join(E(i, j) for i, j in off) + "\\right)"),
    }


# -- Godel, cylindrical --------------------------------------------------------------

def godel():
    """Godel's own cylindrical coordinates of 1949, about one world line of the dust, which the
    Cartesian chart of this entry is carried to by e^x = cosh 2r + cos(phi) sinh 2r,
    y e^x = sqrt(2) sin(phi) sinh 2r, t_x = 2t + sqrt(2)(2 arctan(e^(-2r) tan(phi/2)) - phi) and
    z_x = 2z. Before the chart is written, the published Cartesian metric is pulled back through
    that map and checked equal to this line element's metric, slot by slot; godel.md derives it."""
    coords = ["t", "r", "\\phi", "z"]
    parameters = ["\\omega"]
    line = ("ds^2 = \\dfrac{2}{\\omega^2}\\left(-dt^2 + dr^2 + \\left(\\sinh^2 r - \\sinh^4 r\\right)d\\phi^2"
            " - 2\\sqrt{2}\\sinh^2 r\\,dt\\,d\\phi + dz^2\\right)")
    return {
        "metric_id": "godel",
        "system": {"id": "cylindrical", "name": "Cylindrical", "coords": coords,
                   "domains": ["t \\in (-\\infty, \\infty)", "r \\in [0, \\infty)", "\\phi \\in [0, 2\\pi)",
                               "z \\in (-\\infty, \\infty)",
                               "\\sinh r = 1 \\;\\text{(the circles of constant } t, r, z \\text{ are null)}"],
                   "parameters": parameters, "line_element": line},
        "chart_line_element": line,
        "printer": {"lead": []},
        "pretty": cp.hyperbolic(vm.Reader(coords, parameters, ()).symbol["r"]),
        "check": godel_pullback,
    }


def godel_pullback(chart):
    """J^T g J, with g the metric the Cartesian chart of godel.json publishes and J the Jacobian
    of the map in godel()'s docstring, minus the metric of `chart`, simplified to zero in every
    slot. About seven seconds."""
    metric = json.loads((METRICS / "godel.json").read_text(encoding="utf-8"))
    cartesian = next(c for c in metric["coordinates"] if c["id"] == "cartesian")
    reader = vm.Reader(cartesian["coords"], [p["symbol"] for p in cartesian["parameters"]], ())
    g = sp.zeros(4, 4)
    for comp in cartesian["metric_components"]:
        i, j = (cartesian["coords"].index(x) for x in comp["indices"])
        g[i, j] = g[j, i] = reader(comp["value"])
    t, r, phi, z = chart.symbols
    omega = chart.reader.parameters["omega"]
    spread = sp.cosh(2 * r) + sp.cos(phi) * sp.sinh(2 * r)
    image = [2 * t + sp.sqrt(2) * (2 * sp.atan(sp.exp(-2 * r) * sp.tan(phi / 2)) - phi), sp.log(spread),
             sp.sqrt(2) * sp.sin(phi) * sp.sinh(2 * r) / spread, 2 * z]
    J = sp.Matrix(4, 4, lambda a, b: sp.diff(image[a], chart.symbols[b]))
    at = dict(zip([reader.symbol[c] for c in cartesian["coords"]], image))
    at[reader.parameters["omega"]] = omega
    pulled = J.T * g.subs(at) * J
    for a in range(4):
        for b in range(a, 4):
            if sp.simplify((pulled[a, b] - chart.geo.g[a, b]).rewrite(sp.exp)) != 0:
                raise AssertionError(f"godel: the pullback of the Cartesian metric misses the cylindrical "
                                     f"one in slot {chart.coords_tex[a]}{chart.coords_tex[b]}")


# -- Einstein static universe ----------------------------------------------------------

def einstein_static(system):
    """The three charts of Einstein's static universe, whose one parameter is the radius R of
    its three sphere. Their Ricci scalar is written bare, as the stellar interior's is, since
    "R = 6/R^2" would read as an equation for the radius; einstein_static.md derives each chart."""
    parameters = ["R"]
    charts = {
        "hyperspherical": {
            "name": "Hyperspherical", "coords": ["t", "\\chi", "\\theta", "\\phi"],
            "domains": ["t \\in (-\\infty, \\infty)", "\\chi \\in [0, \\pi]", "\\theta \\in [0, \\pi]",
                        "\\phi \\in [0, 2\\pi)"],
            "space": "R^2\\left(d\\chi^2 + \\sin^2\\chi\\,d\\theta^2 + \\sin^2\\chi\\sin^2\\theta\\,d\\phi^2\\right)"},
        "static_areal": {
            "name": "Areal", "coords": ["t", "r", "\\theta", "\\phi"],
            "domains": ["t \\in (-\\infty, \\infty)", "r \\in [0, R)", "\\theta \\in [0, \\pi]",
                        "\\phi \\in [0, 2\\pi)",
                        "r = R \\;\\text{(the equator of } S^3\\text{; the coordinates cover one hemisphere)}"],
            "space": "\\dfrac{R^2\\,dr^2}{R^2 - r^2} + r^2\\left(d\\theta^2 + \\sin^2\\theta\\,d\\phi^2\\right)"},
        "einstein_cartesian": {
            "name": "Einstein's Cartesian", "coords": ["t", "x", "y", "z"],
            "domains": ["t \\in (-\\infty, \\infty)", "x \\in (-R, R)", "y \\in (-R, R)", "z \\in (-R, R)",
                        "x^2 + y^2 + z^2 < R^2 \\;\\text{(one hemisphere of } S^3\\text{)}"],
            "space": ("dx^2 + dy^2 + dz^2 + \\dfrac{\\left(x\\,dx + y\\,dy + z\\,dz\\right)^2}"
                      "{R^2 - x^2 - y^2 - z^2}")},
    }
    chart = charts[system]
    probe = vm.Reader(chart["coords"], parameters, ())
    return {
        "metric_id": "einstein_static",
        "system": {"id": system, "name": chart["name"], "coords": chart["coords"], "domains": chart["domains"],
                   "parameters": parameters, "line_element": "ds^2 = -c^2dt^2 + " + chart["space"]},
        "chart_line_element": "ds^2 = -dt^2 + " + chart["space"],
        "printer": {"lead": [probe.parameters["R"]]},
        "bare_scalar": True,
        # The printer factors R^2 - r^2 into its two linear factors; the areal chart keeps it whole.
        "rewrite": [("{\\left(R + r\\right)\\left(R - r\\right)}", "{R^2 - r^2}"),
                    ("\\left(R + r\\right)\\left(R - r\\right)", "\\left(R^2 - r^2\\right)")],
        **({"geodesics": GEODESICS[system]} if system in GEODESICS else {}),
    }


# The geodesic equations of the two charts whose printed Christoffel symbols collect into a
# shorter form: the areal chart's angular terms share one factor, and in Einstein's projection
# Gamma^i_{jk} = x^i gamma_{jk}/R^2, so every equation is x^i/R^2 times the spatial speed squared.
_SPEED = ("\\dot{x}^2 + \\dot{y}^2 + \\dot{z}^2 + \\dfrac{\\left(x\\dot{x} + y\\dot{y} + z\\dot{z}\\right)^2}"
          "{R^2 - x^2 - y^2 - z^2}")
GEODESICS = {
    "static_areal": [
        "\\ddot{t} = 0",
        "\\ddot{r} + \\dfrac{r}{R^2 - r^2}\\dot{r}^2 - \\dfrac{r\\left(R^2 - r^2\\right)}{R^2}"
        "\\left(\\dot{\\theta}^2 + \\sin^2\\theta\\,\\dot{\\phi}^2\\right) = 0",
        "\\ddot{\\theta} + \\dfrac{2}{r}\\dot{r}\\dot{\\theta} - \\sin\\theta\\cos\\theta\\,\\dot{\\phi}^2 = 0",
        "\\ddot{\\phi} + \\dfrac{2}{r}\\dot{r}\\dot{\\phi} + 2\\cot\\theta\\,\\dot{\\theta}\\dot{\\phi} = 0",
    ],
    "einstein_cartesian": ["\\ddot{t} = 0"] + [
        f"\\ddot{{{c}}} + \\dfrac{{{c}}}{{R^2}}\\left({_SPEED}\\right) = 0" for c in "xyz"],
}


# -- Milne universe --------------------------------------------------------------------

def milne(system):
    """The four charts of the Milne universe, the inside of the future light cone of one event
    of Minkowski space. The comoving charts depend on the cosmic time, so they name it as the
    chart's `time`, and every value carries it as ct; milne.md derives each chart."""
    hyperbolic = "\\left(d\\chi^2 + \\sinh^2\\chi\\,d\\theta^2 + \\sinh^2\\chi\\sin^2\\theta\\,d\\phi^2\\right)"
    sphere = "\\left(d\\theta^2 + \\sin^2\\theta\\,d\\phi^2\\right)"
    charts = {
        "comoving_hyperbolic": {
            "name": "Comoving Hyperbolic", "coords": ["t", "\\chi", "\\theta", "\\phi"], "parameters": [],
            "domains": ["t \\in (0, \\infty)", "\\chi \\in [0, \\infty)", "\\theta \\in [0, \\pi]",
                        "\\phi \\in [0, 2\\pi)"],
            "line_element": "ds^2 = -c^2dt^2 + c^2t^2" + hyperbolic,
            "chart": "ds^2 = -dt^2 + t^2" + hyperbolic, "time": "t"},
        "comoving_spherical": {
            "name": "Comoving Spherical", "coords": ["t", "r", "\\theta", "\\phi"], "parameters": [],
            "domains": ["t \\in (0, \\infty)", "r \\in [0, \\infty)", "\\theta \\in [0, \\pi]",
                        "\\phi \\in [0, 2\\pi)"],
            "line_element": "ds^2 = -c^2dt^2 + c^2t^2\\left(\\dfrac{dr^2}{1 + r^2} + r^2" + sphere + "\\right)",
            "chart": "ds^2 = -dt^2 + t^2\\left(\\dfrac{dr^2}{1 + r^2} + r^2" + sphere + "\\right)", "time": "t"},
        "logarithmic_time": {
            "name": "Logarithmic Time", "coords": ["\\tau", "\\chi", "\\theta", "\\phi"], "parameters": ["t_0"],
            "domains": ["\\tau \\in (-\\infty, \\infty)", "\\chi \\in [0, \\infty)", "\\theta \\in [0, \\pi]",
                        "\\phi \\in [0, 2\\pi)"],
            "line_element": "ds^2 = e^{2\\tau/t_0}\\left(-c^2d\\tau^2 + c^2t_0^2" + hyperbolic + "\\right)",
            "chart": "ds^2 = e^{2\\tau/(c t_0)}\\left(-d\\tau^2 + c^2t_0^2" + hyperbolic + "\\right)", "time": "\\tau"},
        "inertial": {
            "name": "Inertial", "coords": ["T", "R", "\\theta", "\\phi"], "parameters": [],
            "domains": ["T \\in (0, \\infty)", "R \\in [0, cT)", "\\theta \\in [0, \\pi]", "\\phi \\in [0, 2\\pi)"],
            "line_element": "ds^2 = -c^2dT^2 + dR^2 + R^2" + sphere,
            "chart": "ds^2 = -dT^2 + dR^2 + R^2" + sphere},
    }
    chart = charts[system]
    probe = vm.Reader(chart["coords"], chart["parameters"], ())
    return {
        "metric_id": "milne",
        "system": {"id": system, "name": chart["name"], "coords": chart["coords"], "domains": chart["domains"],
                   "parameters": chart["parameters"], "line_element": chart["line_element"]},
        "chart_line_element": chart["chart"],
        "printer": {"lead": [probe.c, probe.symbol[chart["coords"][0]], *probe.parameters.values()]},
        "bare_scalar": True,
        **({"time": chart["time"]} if "time" in chart else {}),
        **({"pretty": cp.hyperbolic(probe.symbol["\\chi"])} if "\\chi" in chart["coords"] else {}),
    }


def rewritten(value, substitutions, chart):
    """Every string in `value` with each (old, new) of `substitutions` applied in turn, each
    changed string read back and compared with the one it replaces, so a rewrite that changed a
    value stops the script."""
    if isinstance(value, dict):
        return {k: rewritten(v, substitutions, chart) for k, v in value.items()}
    if isinstance(value, list):
        return [rewritten(v, substitutions, chart) for v in value]
    if not isinstance(value, str):
        return value
    text = value
    for old, new in substitutions:
        text = text.replace(old, new)
    if text != value:
        # A side the rewrite left alone, as the R of "R = ...", is not read again.
        for side, original in zip(text.split("="), value.split("="), strict=True):
            if side != original:
                chart.check(side, chart.reader(original))
    return text


# -- Banados-Teitelboim-Zanelli ------------------------------------------------------

BTZ_PARAMETERS = ["\\ell", "M", "J"]
BTZ_LAPSE = "\\dfrac{r^2}{\\ell^2} - M + \\dfrac{J^2}{4r^2}"


def btz_stationary():
    """The chart of Banados, Teitelboim and Zanelli, -N^2 c^2dt^2 + dr^2/N^2 + r^2(dphi + N^phi c dt)^2
    with N^phi = -J/(2r^2), written out so that every term is a parameter or a coordinate."""
    def line(c, c2):
        return (f"ds^2 = -\\left({BTZ_LAPSE}\\right){c2}dt^2 + \\dfrac{{dr^2}}{{{BTZ_LAPSE}}}"
                f" + r^2\\left(d\\phi - \\dfrac{{J}}{{2r^2}}{c}dt\\right)^2")
    return {
        "metric_id": "btz",
        "system": {"id": "stationary", "name": "Stationary", "coords": ["t", "r", "\\phi"],
                   "domains": ["t \\in (-\\infty, \\infty)", "r \\in (0, \\infty)", "\\phi \\in [0, 2\\pi)"],
                   "parameters": BTZ_PARAMETERS, "line_element": line("c\\,", "c^2")},
        "chart_line_element": line("", ""),
        "printer": {"lead": BTZ_LEAD, "factors": BTZ_LEAD[1:] + BTZ_LEAD[:1]},
    }


def btz_null(sign):
    """An Eddington-Finkelstein chart, ingoing for sign 1 and outgoing for sign -1:
    -N^2 dw^2 + 2 sign dw dr + r^2(dphi~ - J dw/(2r^2))^2 with w = v or u a length. It is checked,
    slot by slot, to be the stationary chart pulled back through d(ct) = dw - sign dr/N^2 and
    dphi = dphi~ - sign J dr/(2r^2N^2)."""
    w, name, system = ("v", "Ingoing", "eddington_finkelstein_ingoing") if sign == 1 else \
        ("u", "Outgoing", "eddington_finkelstein_outgoing")
    line = (f"ds^2 = -\\left({BTZ_LAPSE}\\right)d{w}^2 {'+' if sign == 1 else '-'} 2\\,d{w}\\,dr"
            f" + r^2\\left(d\\tilde\\phi - \\dfrac{{J}}{{2r^2}}d{w}\\right)^2")
    return {
        "metric_id": "btz",
        "system": {"id": system, "name": f"{name} Eddington-Finkelstein", "coords": [w, "r", "\\tilde\\phi"],
                   "domains": [f"{w} \\in (-\\infty, \\infty)", "r \\in (0, \\infty)", "\\tilde\\phi \\in [0, 2\\pi)"],
                   "parameters": BTZ_PARAMETERS, "line_element": line},
        "chart_line_element": line,
        "printer": {"lead": BTZ_LEAD, "factors": BTZ_LEAD[1:] + BTZ_LEAD[:1]},
        "check": lambda chart: btz_pullback(chart, sign),
    }


def btz_pullback(chart, sign):
    spec = btz_stationary()
    source = cp.Chart(spec["system"]["coords"], BTZ_PARAMETERS, spec["chart_line_element"])
    w, r, psi = chart.symbols
    ell, M, J = (chart.reader.parameters[n] for n in ("ell", "M", "J"))
    N2 = r ** 2 / ell ** 2 - M + J ** 2 / (4 * r ** 2)
    jacobian = sp.Matrix([[1, -sign / N2, 0], [0, 1, 0], [0, -sign * J / (2 * r ** 2 * N2), 1]])
    at = dict(zip(source.symbols, (w, r, psi)))
    at.update({source.reader.parameters[n]: chart.reader.parameters[n] for n in ("ell", "M", "J")})
    pulled = jacobian.T * source.geo.g.subs(at) * jacobian
    for a in range(3):
        for b in range(a, 3):
            if vm.norm(pulled[a, b] - chart.geo.g[a, b]) != 0:
                raise AssertionError(f"btz: the pullback of the stationary metric misses the "
                                     f"{chart.coords_tex[0]} chart in slot {chart.coords_tex[a]}{chart.coords_tex[b]}")


BTZ_LEAD = [sp.Symbol("r", real=True), sp.Symbol("M", real=True), sp.Symbol("J", real=True),
            sp.Symbol("ell", real=True)]


# -- C-metric ----------------------------------------------------------------------------

def named_factors(named, merges=()):
    """A pretty printer that factors a value and prints the factors the chart names as the
    sums it writes them as. `named` is [(placeholder, polynomial, text)]; `merges` is
    [(a, b, product)], two polynomials whose product is written as one expression, as
    (1 + alpha r)(1 - alpha r) as 1 - alpha^2 r^2. Each factor is matched up to its sign."""

    def matched(base, poly):
        if sp.expand(base - poly) == 0:
            return 1
        if sp.expand(base + poly) == 0:
            return -1
        return 0

    def pretty(value):
        powers = {}
        sign = sp.Integer(1)
        for f in sp.Mul.make_args(sp.factor(value)):
            base, k = (f.base, f.exp) if f.is_Pow else (f, sp.Integer(1))
            powers[base] = powers.get(base, 0) + k
        for a, b, product in merges:
            found = {}
            for base in list(powers):
                for key, poly in (("a", a), ("b", b)):
                    s = matched(base, poly)
                    if s and powers[base].is_Integer and key not in found:
                        found[key] = (base, s)
            if len(found) < 2:
                continue
            (ba, sa), (bb, sb) = found["a"], found["b"]
            ka, kb = powers[ba], powers[bb]
            k = min(ka, kb) if ka > 0 and kb > 0 else max(ka, kb) if ka < 0 and kb < 0 else 0
            if k == 0:
                continue
            powers[ba] -= k
            powers[bb] -= k
            sign *= (sa * sb) ** k
            powers[product] = powers.get(product, 0) + k
        out = sign
        for base, k in powers.items():
            if k == 0:
                continue
            if base.is_Add:
                for placeholder, poly, _ in named:
                    s = matched(base, poly)
                    if s:
                        base, out = placeholder, out * s ** k
                        break
            out *= base ** k
        return out

    return pretty


def c_metric():
    """The C-metric in the two charts of J. B. Griffiths, P. Krtous and J. Podolsky, Class.
    Quantum Grav. 23, 6745 (2006): their spherical chart, eq. (6)-(7), and the form of K. Hong
    and E. Teo, their eq. (3)-(4), in which conformal infinity is x + y = 0. Both write the
    period of phi as 2 pi C; here phi runs over 2 pi and C stands in g_phiphi, as in their
    eq. (20). The spherical chart's cos^2 theta is kept against sin^2 theta so that
    1 + alpha r cos(theta) survives as a factor, and before the Hong-Teo chart is written its
    metric is pulled back through tau = alpha t, y = 1/(alpha r), x = cos(theta) and checked
    equal to the spherical chart's, slot by slot."""
    coords = ["t", "r", "\\theta", "\\phi"]
    parameters = ["m", "\\alpha", "C"]
    probe = vm.Reader(coords, parameters, ())
    r, th = probe.symbol["r"], probe.symbol["\\theta"]
    m, al = probe.parameters["m"], probe.parameters["alpha"]
    s, c = sp.sin(th), sp.cos(th)
    named = [(sp.Symbol("CW"), 1 + al * r * c, "1 + \\alpha r\\cos\\theta"),
             (sp.Symbol("CP"), 1 + 2 * al * m * c, "1 + 2\\alpha m\\cos\\theta"),
             (sp.Symbol("CQ"), 1 - al ** 2 * r ** 2, "1 - \\alpha^2r^2"),
             (sp.Symbol("CA"), 1 - al * r, "1 - \\alpha r"),
             (sp.Symbol("CB"), 1 + al * r, "1 + \\alpha r")]
    factors = named_factors(named, [(1 + al * r, 1 - al * r, named[2][0]), (1 + c, 1 - c, s ** 2)])

    def pretty(value):
        # sin^2 as 1 - cos^2, so that cos(theta) is the one angle every factor is written in.
        def cos_only(e):
            n = int(e.exp)
            return s ** (n % 2) * (1 - c ** 2) ** (n // 2)
        return factors(value.replace(lambda e: e.is_Pow and e.base == s and e.exp.is_Integer, cos_only))

    spherical_body = ("\\left(-\\left(1 - \\alpha^2r^2\\right)\\left(1 - \\dfrac{2m}{r}\\right){}dt^2"
                      " + \\dfrac{dr^2}{\\left(1 - \\alpha^2r^2\\right)\\left(1 - \\dfrac{2m}{r}\\right)}"
                      " + \\dfrac{r^2d\\theta^2}{1 + 2\\alpha m\\cos\\theta}"
                      " + C^2\\left(1 + 2\\alpha m\\cos\\theta\\right)r^2\\sin^2\\theta\\,d\\phi^2\\right)")
    spherical_line = "ds^2 = \\dfrac{1}{\\left(1 + \\alpha r\\cos\\theta\\right)^2}" + spherical_body
    spherical = {
        "metric_id": "c_metric",
        "system": {"id": "spherical", "name": "Spherical", "coords": coords,
                   "domains": ["t \\in (-\\infty, \\infty)", "r \\in (0, \\infty)", "\\theta \\in [0, \\pi]",
                               "\\phi \\in [0, 2\\pi)", "1 + \\alpha r\\cos\\theta > 0"],
                   "parameters": parameters,
                   "line_element": spherical_line.replace("{}dt", "c^2dt")},
        "chart_line_element": spherical_line.replace("{}dt", "dt"),
        "printer": {"lead": [al, r, m, c, s], "named": {p: text for p, _, text in named}},
        "pretty": pretty,
    }

    xy = ["\\tau", "y", "x", "\\phi"]
    xprobe = vm.Reader(xy, parameters, ())
    x, y = xprobe.symbol["x"], xprobe.symbol["y"]
    xm, xal = xprobe.parameters["m"], xprobe.parameters["alpha"]
    xnamed = [(sp.Symbol("CS"), x + y, "x + y"),
              (sp.Symbol("CG"), 1 + 2 * xal * xm * x, "1 + 2\\alpha m x"),
              (sp.Symbol("CF"), 1 - 2 * xal * xm * y, "1 - 2\\alpha m y"),
              (sp.Symbol("CX"), 1 - x ** 2, "1 - x^2"),
              (sp.Symbol("CY"), 1 - y ** 2, "1 - y^2")]
    xy_line = ("ds^2 = \\dfrac{1}{\\alpha^2\\left(x + y\\right)^2}\\left(\\left(1 - y^2\\right)\\left(1 - 2\\alpha m y\\right)d\\tau^2"
               " - \\dfrac{dy^2}{\\left(1 - y^2\\right)\\left(1 - 2\\alpha m y\\right)}"
               " + \\dfrac{dx^2}{\\left(1 - x^2\\right)\\left(1 + 2\\alpha m x\\right)}"
               " + C^2\\left(1 - x^2\\right)\\left(1 + 2\\alpha m x\\right)d\\phi^2\\right)")
    hong_teo = {
        "metric_id": "c_metric",
        "system": {"id": "hong_teo", "name": "Hong-Teo", "coords": xy,
                   "domains": ["\\tau \\in (-\\infty, \\infty)", "y \\in (-\\infty, \\infty)", "x \\in [-1, 1]",
                               "\\phi \\in [0, 2\\pi)", "x + y > 0"],
                   "parameters": parameters, "line_element": xy_line},
        "chart_line_element": xy_line,
        "printer": {"lead": [xal, xm, y, x], "named": {p: text for p, _, text in xnamed}},
        "pretty": named_factors(xnamed, [(1 + x, 1 - x, xnamed[3][0]), (1 + y, 1 - y, xnamed[4][0])]),
        "check": lambda chart: c_metric_pullback(chart, spherical["chart_line_element"]),
    }
    return [spherical, hong_teo]


def c_metric_pullback(chart, spherical_line):
    """J^T g J, with g the Hong-Teo metric of `chart` and J the Jacobian of tau = alpha t,
    y = 1/(alpha r), x = cos(theta), phi = phi, minus the spherical chart's metric, in every slot."""
    coords = ["t", "r", "\\theta", "\\phi"]
    reader = vm.Reader(coords, ["m", "\\alpha", "C"], ())
    g = vm.metric_from_line_element(reader, spherical_line, coords)
    t, r, th, ph = (reader.symbol[n] for n in coords)
    at = {chart.reader.parameters[k]: reader.parameters[k] for k in ("m", "alpha", "C")}
    al = reader.parameters["alpha"]
    image = [al * t, 1 / (al * r), sp.cos(th), ph]
    J = sp.Matrix(4, 4, lambda a, b: sp.diff(image[a], [t, r, th, ph][b]))
    at.update(dict(zip(chart.symbols, image)))
    pulled = J.T * chart.geo.g.subs(at) * J
    for a in range(4):
        for b in range(a, 4):
            if vm.norm(pulled[a, b] - g[a, b]) != 0:
                raise AssertionError(f"c_metric: the Hong-Teo metric pulled back misses the spherical one "
                                     f"in slot {coords[a]}{coords[b]}")


# -- Schwarzschild-de Sitter -------------------------------------------------------------

def schwarzschild_de_sitter(system_id):
    """Kottler's static chart and the two Eddington-Finkelstein charts built on its tortoise
    coordinate, dr_*/dr = 1/f with f = 1 - r_s/r - Lambda r^2/3. Every value is printed around
    3rf = 3r - 3r_s - Lambda r^3, in the order of f itself, so that each chart reduces to
    Schwarzschild's at Lambda = 0 and to de Sitter's static chart at r_s = 0 term by term. The
    metric and its inverse are written as the line element writes f, and the Kretschmann
    scalar as Schwarzschild's 12r_s^2/r^6 plus de Sitter's 8Lambda^2/3, which it is."""
    f = "\\left(1 - \\dfrac{r_s}{r} - \\dfrac{\\Lambda r^2}{3}\\right)"
    bare = "1 - \\dfrac{r_s}{r} - \\dfrac{\\Lambda r^2}{3}"
    sphere = " + r^2\\left(d\\theta^2 + \\sin^2\\theta\\,d\\phi^2\\right)"
    domains = ["r \\in (0, \\infty)", "\\theta \\in [0, \\pi]", "\\phi \\in [0, 2\\pi)",
               "r = r_h \\;\\text{(black hole horizon)}", "r = r_c \\;\\text{(cosmological horizon)}"]
    if system_id == "static":
        coords = ["t", "r", "\\theta", "\\phi"]
        name = "Static Spherical"
        line = "ds^2 = -" + f + "c^2dt^2 + \\dfrac{dr^2}{" + bare + "}" + sphere
        chart_line = "ds^2 = -" + f + "dt^2 + \\dfrac{dr^2}{" + bare + "}" + sphere
        metric = {("t", "t"): "-" + f, ("r", "r"): f + "^{-1}"}
        inverse = {("t", "t"): "-" + f + "^{-1}", ("r", "r"): bare}
    else:
        null, sign = ("u", "-") if system_id == "eddington_finkelstein_outgoing" else ("v", "+")
        coords = [null, "r", "\\theta", "\\phi"]
        name = ("Outgoing" if null == "u" else "Ingoing") + " Eddington-Finkelstein"
        line = "ds^2 = -" + f + "d" + null + "^2 " + sign + " 2\\,d" + null + "\\,dr" + sphere
        chart_line = line
        one = "-1" if null == "u" else "1"
        metric = {(null, null): "-" + f, (null, "r"): one, ("r", null): one}
        inverse = {(null, "r"): one, ("r", null): one, ("r", "r"): bare}
    parameters = ["r_s", "\\Lambda"]
    probe = vm.Reader(coords, parameters, ())
    r, rs, L = probe.symbol["r"], probe.parameters["r_s"], probe.parameters["Lambda"]
    return {
        "metric_id": "schwarzschild_de_sitter",
        "system": {"id": system_id, "name": name, "coords": coords,
                   "domains": [coords[0] + " \\in (-\\infty, \\infty)"] + domains,
                   "parameters": parameters, "line_element": line},
        "chart_line_element": chart_line,
        "printer": {"rising": [L, rs], "lead": [L, r, rs], "flip": False},
        "components": {"metric_components": metric, "inverse_metric_components": inverse},
        "kretschmann": "\\dfrac{12r_s^2}{r^6} + \\dfrac{8\\Lambda^2}{3}",
    }


SDS_CHARTS = ["static", "eddington_finkelstein_outgoing", "eddington_finkelstein_ingoing"]


CHARTS = {"tov": tov, "malament_hogarth": malament_hogarth, "mixmaster": mixmaster, "lentz": lentz, "godel": godel,
          "einstein_static": [lambda s=s: einstein_static(s) for s in ("hyperspherical", "static_areal", "einstein_cartesian")],
          "btz": [lambda: btz_stationary(), lambda: btz_null(1), lambda: btz_null(-1)],
          "c_metric": c_metric,
          "schwarzschild_de_sitter": [lambda s=s: schwarzschild_de_sitter(s) for s in SDS_CHARTS],
          "milne": [lambda s=s: milne(s) for s in ("comoving_hyperbolic", "comoving_spherical", "logarithmic_time", "inertial")],
          "einstein_rosen_waves": [lambda s=s: einstein_rosen(s) for s in ("cylindrical", "null")]}


# -- Einstein-Rosen ------------------------------------------------------------------

ER_PARAMETERS = {"cylindrical": ["\\psi = \\psi(t,\\rho)", "\\gamma = \\gamma(t,\\rho)"],
                 "null": ["\\psi = \\psi(u,v)", "\\gamma = \\gamma(u,v)"]}


def einstein_rosen(system):
    """The Einstein-Rosen waves in the cylindrical chart of Einstein and Rosen, with psi and
    gamma free functions of t and rho, and in the null chart u = ct - rho, v = ct + rho, both
    lengths. The null chart is checked, slot by slot, to be the cylindrical metric pulled back
    through ct = (v + u)/2 and rho = (v - u)/2, with psi and gamma carried along as functions of
    the same event. No component assumes a field equation."""
    parameters = ER_PARAMETERS[system]
    if system == "cylindrical":
        coords = ["t", "\\rho", "\\phi", "z"]

        def line(c2, rho2):
            return (f"ds^2 = e^{{2(\\gamma - \\psi)}}\\left(-{c2}dt^2 + d\\rho^2\\right)"
                    f" + {rho2}e^{{-2\\psi}}d\\phi^2 + e^{{2\\psi}}dz^2")
        spec_line, chart_line = line("c^2", "\\rho^2"), line("", "\\rho^2")
        domains = ["t \\in (-\\infty, \\infty)", "\\rho \\in [0, \\infty)", "\\phi \\in [0, 2\\pi)",
                   "z \\in (-\\infty, \\infty)",
                   "\\gamma = 0 \\;\\text{at}\\; \\rho = 0 \\;\\text{(a regular axis)}"]
        name = "Cylindrical"
    else:
        coords = ["u", "v", "\\phi", "z"]

        def line(c2, rho2):
            return (f"ds^2 = -e^{{2(\\gamma - \\psi)}}{c2}du\\,dv"
                    f" + {rho2}e^{{-2\\psi}}d\\phi^2 + e^{{2\\psi}}dz^2")
        spec_line = chart_line = line("", "\\dfrac{(v - u)^2}{4}")
        domains = ["u \\in (-\\infty, \\infty)", "v \\in [u, \\infty)", "\\phi \\in [0, 2\\pi)",
                   "z \\in (-\\infty, \\infty)",
                   "\\gamma = 0 \\;\\text{at}\\; v = u \\;\\text{(a regular axis)}"]
        name = "Null"
    probe = vm.Reader(coords, parameters, ())
    x, y = probe.symbol[coords[0]], probe.symbol[coords[1]]
    psi, gam = probe.parameters["psi"], probe.parameters["gamma"]
    D = sp.Derivative
    lead = [gam, psi, D(gam, y), D(gam, x), D(psi, y), D(psi, x), D(psi, (y, 2)), D(psi, x, y), D(psi, (x, 2)),
            D(gam, (y, 2)), D(gam, x, y), D(gam, (x, 2))]
    if system == "cylindrical":
        lead.append(probe.symbol["\\rho"])
        overrides, width = {}, {}
    else:
        # Twice the radius, v - u, is printed whole wherever it stands, and never as v and u apart.
        W = sp.Symbol("W", positive=True)
        lead.append(W)
        overrides, width = {W: "\\left(v - u\\right)"}, {y: x + W}
    spec = {
        "metric_id": "einstein_rosen_waves",
        "system": {"id": system, "name": name, "coords": coords, "domains": domains,
                   "parameters": parameters, "line_element": spec_line},
        "chart_line_element": chart_line,
        "printer": {"lead": lead, "overrides": overrides, "factors": list(overrides) + lead},
        "pretty": lambda value: sp.powsimp(sp.factor(outside_functions(sp.sympify(value), width)), combine="exp"),
        "rewrite": ER_REWRITES,
        "kretschmann_text": einstein_rosen_kretschmann,
        "geodesics": ER_GEODESICS[system],
    }
    if system == "null":
        spec["check"] = einstein_rosen_pullback
    return spec


def einstein_rosen_kretschmann(chart):
    """K in the orthonormal frame e^{psi - gamma} d_t, e^{psi - gamma} d_rho, (e^psi/rho) d_phi,
    e^{-psi} d_z, with d_t = d_u + d_v and d_rho = d_v - d_u in the null chart: every nonzero
    frame component of Riemann is e^{2(psi - gamma)} times a bracket, six pair a bivector with
    itself and two, R_{0 phi 1 phi} and R_{0 z 1 z}, pair a time bivector with a space one, so
    K = 4e^{4(psi - gamma)} times the six brackets squared less twice the two squared."""
    riemann = chart.geo.riemann_llll()
    psi, gam = chart.reader.parameters["psi"], chart.reader.parameters["gamma"]
    a, b = chart.symbols[:2]
    if chart.coords_tex[0] == "t":
        rho, plane = b, [[1, 0], [0, 1]]
    else:
        rho, plane = (b - a) / 2, [[1, 1], [-1, 1]]
    lapse = sp.exp(psi - gam)
    frame = [[lapse * plane[0][0], lapse * plane[0][1], 0, 0], [lapse * plane[1][0], lapse * plane[1][1], 0, 0],
             [0, 0, sp.exp(psi) / rho, 0], [0, 0, 0, sp.exp(-psi)]]

    def component(A, B, C, E):
        total = 0
        for i in range(4):
            for j in range(4):
                for k in range(4):
                    for m in range(4):
                        w = frame[A][i] * frame[B][j] * frame[C][k] * frame[E][m]
                        if w != 0:
                            total += w * vm._at(riemann, (i, j, k, m))
        return vm.norm(sp.expand(total * sp.exp(2 * gam - 2 * psi)))

    diagonal = [(0, 1, 0, 1), (0, 2, 0, 2), (0, 3, 0, 3), (1, 2, 1, 2), (1, 3, 1, 3), (2, 3, 2, 3)]
    mixed = [(0, 2, 1, 2), (0, 3, 1, 3)]

    def square(index):
        value = component(*index)
        text = chart.text(value)
        if text.startswith("-"):
            text = chart.text(-value)
        return "\\left(" + text + "\\right)^2"

    body = " + ".join(square(i) for i in diagonal) + " - 2" + " - 2".join(square(i) for i in mixed)
    return "4e^{4\\psi - 4\\gamma}\\left(" + body + "\\right)"


def outside_functions(value, substitutions):
    """`value` with `substitutions` made in the bare coordinates only, every function and
    derivative of one held apart, so psi(u, v) keeps its arguments."""
    held = {f: sp.Dummy() for f in value.atoms(sp.Derivative) | value.atoms(sp.core.function.AppliedUndef)}
    back = {d: f for f, d in held.items()}
    return value.xreplace(held).subs(substitutions).xreplace(back)


# e^{2 psi - 2 gamma} as the files write a difference, leading with its positive term.
ER_REWRITES = [("e^{-2\\gamma + 2\\psi}", "e^{2\\psi - 2\\gamma}"), ("e^{-4\\gamma + 4\\psi}", "e^{4\\psi - 4\\gamma}"),
               ("e^{-2\\gamma + 4\\psi}", "e^{4\\psi - 2\\gamma}"), ("e^{-4\\gamma + 2\\psi}", "e^{2\\psi - 4\\gamma}"),
               ("}{\\left(v - u\\right)}", "}{v - u}")]

# Each geodesic equation with the terms the printed Christoffel symbols share gathered.
ER_GEODESICS = {
    "cylindrical": [
        "\\ddot{t} + \\left(\\partial_t\\gamma - \\partial_t\\psi\\right)\\left(\\dot{t}^2 + \\dot{\\rho}^2\\right)"
        " + 2\\left(\\partial_\\rho\\gamma - \\partial_\\rho\\psi\\right)\\dot{t}\\dot{\\rho}"
        " - \\rho^2e^{-2\\gamma}\\partial_t\\psi\\,\\dot{\\phi}^2 + e^{4\\psi - 2\\gamma}\\partial_t\\psi\\,\\dot{z}^2 = 0",
        "\\ddot{\\rho} + \\left(\\partial_\\rho\\gamma - \\partial_\\rho\\psi\\right)\\left(\\dot{t}^2 + \\dot{\\rho}^2\\right)"
        " + 2\\left(\\partial_t\\gamma - \\partial_t\\psi\\right)\\dot{t}\\dot{\\rho}"
        " + \\rho\\left(\\rho\\,\\partial_\\rho\\psi - 1\\right)e^{-2\\gamma}\\dot{\\phi}^2 - e^{4\\psi - 2\\gamma}\\partial_\\rho\\psi\\,\\dot{z}^2 = 0",
        "\\ddot{\\phi} - 2\\partial_t\\psi\\,\\dot{t}\\dot{\\phi}"
        " + \\dfrac{2\\left(1 - \\rho\\,\\partial_\\rho\\psi\\right)}{\\rho}\\dot{\\rho}\\dot{\\phi} = 0",
        "\\ddot{z} + 2\\partial_t\\psi\\,\\dot{t}\\dot{z} + 2\\partial_\\rho\\psi\\,\\dot{\\rho}\\dot{z} = 0",
    ],
    "null": [
        "\\ddot{u} + 2\\left(\\partial_u\\gamma - \\partial_u\\psi\\right)\\dot{u}^2"
        " - \\dfrac{\\left(v - u\\right)\\left(\\left(v - u\\right)\\partial_v\\psi - 1\\right)e^{-2\\gamma}}{2}\\dot{\\phi}^2"
        " + 2e^{4\\psi - 2\\gamma}\\partial_v\\psi\\,\\dot{z}^2 = 0",
        "\\ddot{v} + 2\\left(\\partial_v\\gamma - \\partial_v\\psi\\right)\\dot{v}^2"
        " - \\dfrac{\\left(v - u\\right)\\left(\\left(v - u\\right)\\partial_u\\psi + 1\\right)e^{-2\\gamma}}{2}\\dot{\\phi}^2"
        " + 2e^{4\\psi - 2\\gamma}\\partial_u\\psi\\,\\dot{z}^2 = 0",
        "\\ddot{\\phi} - \\dfrac{2\\left(1 + \\left(v - u\\right)\\partial_u\\psi\\right)}{v - u}\\dot{u}\\dot{\\phi}"
        " + \\dfrac{2\\left(1 - \\left(v - u\\right)\\partial_v\\psi\\right)}{v - u}\\dot{v}\\dot{\\phi} = 0",
        "\\ddot{z} + 2\\partial_u\\psi\\,\\dot{u}\\dot{z} + 2\\partial_v\\psi\\,\\dot{v}\\dot{z} = 0",
    ],
}


def einstein_rosen_pullback(chart):
    spec = einstein_rosen("cylindrical")
    source = cp.Chart(spec["system"]["coords"], spec["system"]["parameters"], spec["chart_line_element"])
    u, v = chart.symbols[:2]
    P, Q = sp.symbols("P Q")
    # The metric components carry psi and gamma undifferentiated, so each is one value at an event.
    values = {source.reader.parameters["psi"]: P, source.reader.parameters["gamma"]: Q}
    target = {chart.reader.parameters["psi"]: P, chart.reader.parameters["gamma"]: Q}
    t, rho = source.symbols[:2]
    jacobian = sp.Matrix([[sp.Rational(1, 2), sp.Rational(1, 2), 0, 0], [-sp.Rational(1, 2), sp.Rational(1, 2), 0, 0],
                          [0, 0, 1, 0], [0, 0, 0, 1]])
    at = {t: (u + v) / 2, rho: (v - u) / 2}
    at.update(dict(zip(source.symbols[2:], chart.symbols[2:])))
    pulled = jacobian.T * source.geo.g.subs(values).subs(at) * jacobian
    for a in range(4):
        for b in range(a, 4):
            if vm.norm(pulled[a, b] - chart.geo.g[a, b].subs(target)) != 0:
                raise AssertionError(f"einstein_rosen_waves: the pullback of the cylindrical metric misses the "
                                     f"null chart in slot {chart.coords_tex[a]}{chart.coords_tex[b]}")



def write(spec):
    start = time.time()
    chart = cp.Chart(spec["system"]["coords"], spec["system"]["parameters"], spec["chart_line_element"],
                     spec["printer"], spec.get("pretty"), spec.get("time"))
    if "check" in spec:
        spec["check"](chart)
    math = chart.mathematics()
    # A component written by hand replaces the printed one, checked against the same value.
    tensors = {"metric_components": chart.geo.g, "inverse_metric_components": chart.geo.ginv}
    for field, by_hand in spec.get("components", {}).items():
        for entry in math[field]:
            text = by_hand.get(tuple(entry["indices"]))
            if text is not None:
                i, j = (spec["system"]["coords"].index(x) for x in entry["indices"])
                entry["value"] = chart.check(text, tensors[field][i, j])
    for field in ("ricci_scalar", "kretschmann"):
        computed = chart.geo.ricci_scalar() if field == "ricci_scalar" else chart.geo.kretschmann()
        if field in spec:
            math[field] = ("R = " if field == "ricci_scalar" else "K = ") + chart.check(spec[field], computed)
        elif field == "kretschmann" and "kretschmann_text" in spec:
            math[field] = "K = " + chart.check(spec["kretschmann_text"](chart), computed)
        elif field == "kretschmann":
            math[field] = "K = " + chart.text(computed)
    if spec.get("bare_scalar"):
        math["ricci_scalar"] = math["ricci_scalar"].removeprefix("R = ")
    if "rewrite" in spec:
        math = rewritten(math, spec["rewrite"], chart)
    if "geodesics" in spec:
        # Each equation written by hand is read back against the one printed from the Christoffel symbols.
        for text, printed in zip(spec["geodesics"], math["geodesics"], strict=True):
            chart.check(text.partition("=")[0], chart.reader(printed.partition("=")[0]))
        math["geodesics"] = list(spec["geodesics"])

    path = METRICS / f"{spec['metric_id']}.json"
    metric = json.loads(path.read_text(encoding="utf-8"))
    charts = metric.get("coordinates", [])
    # A parameter keeps the description it has in this chart, or in another chart of the same
    # spacetime where it means the same thing, as Godel's omega does in both of its charts.
    described = {p["symbol"]: p["description"] for chart in charts if chart["id"] != spec["system"]["id"]
                 for p in chart.get("parameters", [])}
    described.update({p["symbol"]: p["description"] for chart in charts if chart["id"] == spec["system"]["id"]
                      for p in chart.get("parameters", [])})
    missing = [s for s in spec["system"]["parameters"] if s not in described]
    if missing:
        raise SystemExit(f"{path.name}: describe the parameters {missing} in the file before writing it")
    # The chart keeps the convention the file gives it, which says what its own coordinates are.
    convention = next((chart.get("convention") for chart in charts if chart["id"] == spec["system"]["id"]), None)
    if not convention:
        raise SystemExit(f"{path.name}: write the convention of the chart {spec['system']['id']!r} "
                         "in the file before writing it")
    system = dict(spec["system"])
    system["parameters"] = [{"symbol": s, "description": described[s]} for s in spec["system"]["parameters"]]
    system["convention"] = convention
    system.update(math)
    # The chart replaces the one of its id, or joins the spacetime's others after them.
    written = {k: system[k] for k in SYSTEM_ORDER}
    ids = [chart["id"] for chart in charts]
    if spec["system"]["id"] in ids:
        charts[ids.index(spec["system"]["id"])] = written
    else:
        charts.append(written)
    metric["coordinates"] = charts
    path.write_text(json.dumps(metric, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"{path.name}: {time.time() - start:.0f}s, {path.stat().st_size} bytes", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--metric", action="append", choices=sorted(CHARTS), default=[])
    for metric_id in parser.parse_args().metric or sorted(CHARTS):
        print(metric_id, flush=True)
        # A spacetime with several charts printed by machine lists one builder per chart, or has
        # one builder that returns every chart, when the charts are checked against each other.
        for build in CHARTS[metric_id] if isinstance(CHARTS[metric_id], list) else [CHARTS[metric_id]]:
            specs = build()
            for spec in specs if isinstance(specs, list) else [specs]:
                write(spec)


if __name__ == "__main__":
    main()
