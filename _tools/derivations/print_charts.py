#!/usr/bin/env python3
"""Compute and write the coordinate systems whose mathematics is printed by machine: the
charts of tov, malament_hogarth, mixmaster, lentz, einstein_static, btz, c_metric,
schwarzschild_de_sitter, milne, einstein_rosen_waves, nariai, aichelburg_sexl,
khan_penrose and global_monopole, and Godel's cylindrical chart.

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
import itertools
import json
import sys
import time
from pathlib import Path

import sympy as sp
from sympy.core.mul import _keep_coeff

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


# -- Global monopole -------------------------------------------------------------------

def global_monopole(system_id):
    """The field of a global monopole, which is Letelier's cloud of strings about a mass, in
    the static chart f = 1 - Delta - r_s/r and the two Eddington-Finkelstein charts built on its
    tortoise coordinate dr_*/dr = 1/f, and Barriola and Vilenkin's conical chart, the static
    chart at r_s = 0 with t and r rescaled by sqrt(1 - Delta). Every value is printed around
    rf = (1 - Delta)r - r_s, so that each chart reduces to Schwarzschild's at Delta = 0 term by
    term. The Kretschmann scalar is written as Schwarzschild's 12r_s^2/r^6 plus the terms in
    Delta, and the conical chart is checked, slot by slot, to be the static chart at r_s = 0
    pulled back through t = T/sqrt(1 - Delta), r = sqrt(1 - Delta) R."""
    f = "\\left(1 - \\Delta - \\dfrac{r_s}{r}\\right)"
    bare = "1 - \\Delta - \\dfrac{r_s}{r}"
    sphere = " + r^2\\left(d\\theta^2 + \\sin^2\\theta\\,d\\phi^2\\right)"
    angles = ["\\theta \\in [0, \\pi]", "\\phi \\in [0, 2\\pi)"]
    horizon = ["r = r_s/(1 - \\Delta) \;\\text{(the horizon, for}\; r_s > 0\\text{)}"]
    parameters = ["\\Delta", "r_s"]
    kretschmann = "\\dfrac{12r_s^2 + 8\\Delta\\,r_s\\,r + 4\\Delta^2r^2}{r^6}"
    if system_id == "conical":
        coords = ["t", "r", "\\theta", "\\phi"]
        parameters = ["\\Delta"]
        cone = " + \\left(1 - \\Delta\\right)r^2\\left(d\\theta^2 + \\sin^2\\theta\\,d\\phi^2\\right)"
        line, chart_line = "ds^2 = -c^2dt^2 + dr^2" + cone, "ds^2 = -dt^2 + dr^2" + cone
        domains = ["t \\in (-\\infty, \\infty)", "r \\in (0, \\infty)"] + angles
        extra = {"check": global_monopole_cone}
        name = "Barriola-Vilenkin"
    elif system_id == "static":
        coords = ["t", "r", "\\theta", "\\phi"]
        name = "Static Spherical"
        line = "ds^2 = -" + f + "c^2dt^2 + \\dfrac{dr^2}{" + bare + "}" + sphere
        chart_line = "ds^2 = -" + f + "dt^2 + \\dfrac{dr^2}{" + bare + "}" + sphere
        domains = ["t \\in (-\\infty, \\infty)", "r \\in (0, \\infty)"] + angles + horizon
        extra = {"components": {"metric_components": {("t", "t"): "-" + f, ("r", "r"): f + "^{-1}"},
                                "inverse_metric_components": {("t", "t"): "-" + f + "^{-1}", ("r", "r"): bare}},
                 "kretschmann": kretschmann}
    else:
        null, sign = ("u", "-") if system_id == "eddington_finkelstein_outgoing" else ("v", "+")
        coords = [null, "r", "\\theta", "\\phi"]
        name = ("Outgoing" if null == "u" else "Ingoing") + " Eddington-Finkelstein"
        line = chart_line = "ds^2 = -" + f + "d" + null + "^2 " + sign + " 2\\,d" + null + "\\,dr" + sphere
        one = "-1" if null == "u" else "1"
        domains = [null + " \\in (-\\infty, \\infty)", "r \\in (0, \\infty)"] + angles + horizon
        extra = {"components": {"metric_components": {(null, null): "-" + f, (null, "r"): one, ("r", null): one},
                                "inverse_metric_components": {(null, "r"): one, ("r", null): one, ("r", "r"): bare}},
                 "kretschmann": kretschmann}
    probe = vm.Reader(coords, parameters, ())
    r, D = probe.symbol["r"], probe.parameters["Delta"]
    lead = [D, r] + ([probe.parameters["r_s"]] if "r_s" in parameters else [])
    return {
        "metric_id": "global_monopole",
        "system": {"id": system_id, "name": name, "coords": coords, "domains": domains,
                   "parameters": parameters, "line_element": line},
        "chart_line_element": chart_line,
        "printer": {"rising": [D] + lead[2:], "lead": lead, "flip": False},
        **extra,
    }


def global_monopole_cone(chart):
    """J^T g J, with g the static chart at r_s = 0 and J the Jacobian of ct = cT/sqrt(1 - Delta),
    r = sqrt(1 - Delta) R, against the conical chart's metric, in every slot."""
    spec = global_monopole("static")
    source = cp.Chart(spec["system"]["coords"], spec["system"]["parameters"], spec["chart_line_element"])
    D = chart.reader.parameters["Delta"]
    k = sp.sqrt(1 - D)
    T, R = chart.symbols[:2]
    at = {source.reader.parameters["Delta"]: D, source.reader.parameters["r_s"]: 0,
          source.symbols[0]: T / k, source.symbols[1]: k * R}
    at.update(dict(zip(source.symbols[2:], chart.symbols[2:])))
    J = sp.diag(1 / k, k, 1, 1)
    pulled = J.T * source.geo.g.subs(at) * J
    for a in range(4):
        for b in range(a, 4):
            if sp.simplify(pulled[a, b] - chart.geo.g[a, b]) != 0:
                raise AssertionError(f"global_monopole: the static chart pulled back misses the conical "
                                     f"chart in slot {chart.coords_tex[a]}{chart.coords_tex[b]}")


GM_CHARTS = ["static", "conical", "eddington_finkelstein_outgoing", "eddington_finkelstein_ingoing"]


# -- Nariai ----------------------------------------------------------------------------

NARIAI_SPHERE = "\\dfrac{1}{\\Lambda}\\left(d\\theta^2 + \\sin^2\\theta\\,d\\phi^2\\right)"
NARIAI_DOMAINS = ["\\theta \\in [0, \\pi]", "\\phi \\in [0, 2\\pi)"]


def nariai(system):
    """The three charts of the Nariai universe, dS2 x S2 with both radii a = 1/sqrt(Lambda):
    the static patch between the horizons r = +-a, the global chart in which the circle of chi
    has the radius a cosh(ct/a), and its conformal chart, tan eta = sinh(ct/a). Each is checked,
    slot by slot, to be the metric of the surfaces -Z0^2 + Z1^2 + Z2^2 = a^2 and
    Z3^2 + Z4^2 + Z5^2 = a^2 in flat space of six dimensions, pulled back to the chart."""
    charts = {
        "static": {
            "name": "Static", "coords": ["t", "r", "\\theta", "\\phi"],
            "domains": ["t \\in (-\\infty, \\infty)",
                        "r \\in \\left(-1/\\sqrt{\\Lambda},\\, 1/\\sqrt{\\Lambda}\\right)"] + NARIAI_DOMAINS
                       + ["r = \\pm 1/\\sqrt{\\Lambda} \\;\\text{(the two horizons)}"],
            "line": lambda c2: (f"ds^2 = -\\left(1 - \\Lambda r^2\\right){c2}dt^2 + \\dfrac{{dr^2}}{{1 - \\Lambda r^2}}"
                                f" + {NARIAI_SPHERE}")},
        "global": {
            "name": "Global", "coords": ["t", "\\chi", "\\theta", "\\phi"],
            "domains": ["t \\in (-\\infty, \\infty)", "\\chi \\in [0, 2\\pi)"] + NARIAI_DOMAINS,
            "line": lambda c2: (f"ds^2 = -{c2}dt^2 + \\dfrac{{1}}{{\\Lambda}}\\left(\\cosh^2\\left(\\sqrt{{\\Lambda}}\\,"
                                f"{'c' if c2 else ''}t\\right)d\\chi^2 + d\\theta^2 + \\sin^2\\theta\\,d\\phi^2\\right)")},
        "conformal": {
            "name": "Conformal", "coords": ["\\eta", "\\chi", "\\theta", "\\phi"],
            "domains": ["\\eta \\in (-\\pi/2, \\pi/2)", "\\chi \\in [0, 2\\pi)"] + NARIAI_DOMAINS,
            "line": lambda c2: ("ds^2 = \\dfrac{1}{\\Lambda}\\left(\\dfrac{-d\\eta^2 + d\\chi^2}{\\cos^2\\eta}"
                                " + d\\theta^2 + \\sin^2\\theta\\,d\\phi^2\\right)")},
    }
    chart = charts[system]
    probe = vm.Reader(chart["coords"], ["\\Lambda"], ())
    L = probe.parameters["Lambda"]
    printing = {
        "static": lambda: {"printer": {"lead": [L], "rising": [probe.symbol["r"]], "flip": False}},
        "global": lambda: {"printer": {"lead": [L]}, "time": "t",
                           "pretty": nariai_hyperbolic(sp.sqrt(L) * probe.c * probe.symbol["t"])},
        "conformal": lambda: {"printer": {"lead": [L]}, "pretty": nariai_conformal(probe.symbol["\\eta"])},
    }[system]()
    return {
        "metric_id": "nariai",
        "system": {"id": system, "name": chart["name"], "coords": chart["coords"], "domains": chart["domains"],
                   "parameters": ["\\Lambda"], "line_element": chart["line"]("c^2")},
        "chart_line_element": chart["line"](""),
        "check": lambda c: nariai_embedding(c, system),
        **printing,
    }


def nariai_hyperbolic(r):
    """cp.hyperbolic, with -(sinh^2 r + 1), which sympy leaves a sum, written as -cosh^2 r."""
    pretty = cp.hyperbolic(r)
    return lambda value: pretty(value).replace(lambda e: e.is_Add and sp.expand(e + sp.sinh(r) ** 2 + 1) == 0,
                                               lambda e: -sp.cosh(r) ** 2)


def nariai_conformal(eta):
    """Each value of the de Sitter factor in eta simplified by trigonometry, so
    (sin eta + 1)(sin eta - 1) reads as -cos^2 eta; the sphere's values stay as sympy factors them."""
    return lambda value: sp.factor(sp.trigsimp(sp.factor(value))) if sp.sympify(value).has(eta) else sp.factor(value)


def nariai_embedding(chart, system):
    """J^T eta J for the embedding (Z0, Z1, Z2) of the de Sitter factor, plus a^2 times the
    round sphere, against the chart's metric, with the chart's time already ct."""
    x0, x1, theta, _ = chart.symbols
    L = chart.reader.parameters["Lambda"]
    a = 1 / sp.sqrt(L)
    if system == "static":
        s = sp.sqrt(a ** 2 - x1 ** 2)
        Z = [s * sp.sinh(x0 / a), s * sp.cosh(x0 / a), x1]
    elif system == "global":
        Z = [a * sp.sinh(x0 / a), a * sp.cosh(x0 / a) * sp.cos(x1), a * sp.cosh(x0 / a) * sp.sin(x1)]
    else:
        Z = [a * sp.tan(x0), a * sp.cos(x1) / sp.cos(x0), a * sp.sin(x1) / sp.cos(x0)]
    J = sp.Matrix([[sp.diff(z, v) for v in (x0, x1)] for z in Z])
    pulled = sp.zeros(4, 4)
    pulled[:2, :2] = J.T * sp.diag(-1, 1, 1) * J
    pulled[2, 2], pulled[3, 3] = a ** 2, a ** 2 * sp.sin(theta) ** 2
    positive = sp.Symbol("Lambda_", positive=True)
    for i in range(4):
        for j in range(i, 4):
            diff = (pulled[i, j] - chart.geo.g[i, j]).subs(L, positive)
            if sp.simplify(sp.expand_trig(diff.rewrite(sp.exp))) != 0:
                raise AssertionError(f"nariai: the embedding pulled back misses the {system} chart "
                                     f"in slot {chart.coords_tex[i]}{chart.coords_tex[j]}")


# -- Aichelburg-Sexl ------------------------------------------------------------------

def aichelburg_sexl(system):
    """The three charts of the Aichelburg-Sexl ultraboost, a pp-wave whose profile is
    F(rho) delta(u) with F = -(8GE/c^4) ln(rho/rho_0), harmonic off the axis rho = 0. The
    Cartesian chart is Aichelburg and Sexl's own, with u = ct - z, and the null charts write
    u and v = ct + z as coordinates, with the transverse plane in Cartesian or polar form;
    aichelburg_sexl.md derives each. The profile is written with ln((x^2 + y^2)/rho_0^2) in
    the Cartesian planes, so no radical appears in any component."""
    log_xy = "\\ln\\left(\\dfrac{x^2 + y^2}{\\rho_0^2}\\right)"
    charts = {
        "cartesian": {
            "name": "Cartesian", "coords": ["t", "x", "y", "z"],
            "domains": ["t \\in (-\\infty, \\infty)", "x \\in (-\\infty, \\infty)", "y \\in (-\\infty, \\infty)",
                        "z \\in (-\\infty, \\infty)", "(x, y) \\neq (0, 0)"],
            "line_element": "ds^2 = -c^2dt^2 + dx^2 + dy^2 + dz^2 - \\dfrac{4GE}{c^4}" + log_xy
                            + "\\delta(ct - z)\\left(c\\,dt - dz\\right)^2",
            "chart": "ds^2 = -dt^2 + dx^2 + dy^2 + dz^2 - \\dfrac{4GE}{c^4}" + log_xy
                     + "\\delta(t - z)\\left(dt - dz\\right)^2",
            "time": "t"},
        "null_cartesian": {
            "name": "Null Cartesian", "coords": ["u", "v", "x", "y"],
            "domains": ["u \\in (-\\infty, \\infty)", "v \\in (-\\infty, \\infty)", "x \\in (-\\infty, \\infty)",
                        "y \\in (-\\infty, \\infty)", "(x, y) \\neq (0, 0)"],
            "line_element": "ds^2 = -du\\,dv + dx^2 + dy^2 - \\dfrac{4GE}{c^4}" + log_xy + "\\delta(u)\\,du^2"},
        "null_cylindrical": {
            "name": "Null Cylindrical", "coords": ["u", "v", "\\rho", "\\phi"],
            "domains": ["u \\in (-\\infty, \\infty)", "v \\in (-\\infty, \\infty)", "\\rho \\in (0, \\infty)",
                        "\\phi \\in [0, 2\\pi)"],
            "line_element": "ds^2 = -du\\,dv + d\\rho^2 + \\rho^2d\\phi^2 - \\dfrac{8GE}{c^4}"
                            "\\ln\\left(\\dfrac{\\rho}{\\rho_0}\\right)\\delta(u)\\,du^2"},
    }
    # The Cartesian chart's metric and inverse are written as Minkowski's plus the shock, each
    # term of the shock carrying (c dt - dz)^2 or its dual, which is how they read.
    shock = "\\dfrac{4GE}{c^4}" + log_xy + "\\delta(ct - z)"
    components = {
        "cartesian": {
            "metric_components": {("t", "t"): "-1 - " + shock, ("t", "z"): shock, ("z", "t"): shock,
                                  ("z", "z"): "1 - " + shock},
            "inverse_metric_components": {("t", "t"): "-1 + " + shock, ("t", "z"): shock, ("z", "t"): shock,
                                          ("z", "z"): "1 + " + shock}},
    }
    # The Cartesian chart's geodesics grouped around ct - z, whose rate is the one the shock sees.
    kick = "\\dfrac{8GE\\,\\delta(ct - z)}{c^4\\left(x^2 + y^2\\right)}"
    along = ("\\dfrac{2GE}{c^4}" + log_xy + "\\delta'(ct - z)\\left(\\dot{t} - \\dot{z}\\right)^2 + " + kick
             + "\\left(x\\dot{x} + y\\dot{y}\\right)\\left(\\dot{t} - \\dot{z}\\right) = 0")
    across = "\\dfrac{4GE\\,{0}\\,\\delta(ct - z)}{c^4\\left(x^2 + y^2\\right)}\\left(\\dot{t} - \\dot{z}\\right)^2 = 0"
    null_kick = "\\dfrac{16GE\\,\\delta(u)}{c^4\\left(x^2 + y^2\\right)}"
    geodesics = {
        "cartesian": ["\\ddot{t} + " + along, "\\ddot{x} + " + across.replace("{0}", "x"),
                      "\\ddot{y} + " + across.replace("{0}", "y"), "\\ddot{z} + " + along],
        "null_cartesian": ["\\ddot{u} = 0",
                           "\\ddot{v} + \\dfrac{4GE}{c^4}" + log_xy + "\\delta'(u)\\,\\dot{u}^2 + " + null_kick
                           + "\\left(x\\dot{x} + y\\dot{y}\\right)\\dot{u} = 0",
                           "\\ddot{x} + \\dfrac{4GE\\,x\\,\\delta(u)}{c^4\\left(x^2 + y^2\\right)}\\dot{u}^2 = 0",
                           "\\ddot{y} + \\dfrac{4GE\\,y\\,\\delta(u)}{c^4\\left(x^2 + y^2\\right)}\\dot{u}^2 = 0"],
    }
    chart = charts[system]
    parameters = ["G", "E", "\\rho_0"]
    probe = vm.Reader(chart["coords"], parameters, ())
    G, E = probe.parameters["G"], probe.parameters["E"]
    # G and E enter every component as the product GE, printed together as the line element has it.
    GE = sp.Symbol("GE")
    return {
        **({"components": components[system]} if system in components else {}),
        **({"geodesics": geodesics[system]} if system in geodesics else {}),
        "metric_id": "aichelburg_sexl",
        "system": {"id": system, "name": chart["name"], "coords": chart["coords"], "domains": chart["domains"],
                   "parameters": parameters, "line_element": chart["line_element"]},
        "chart_line_element": chart.get("chart", chart["line_element"]),
        "printer": {"lead": [probe.c, GE, *probe.symbol.values()], "flip": False},
        "pretty": lambda value: sp.factor(value).subs(E, GE / G),
        **({"time": chart["time"]} if "time" in chart else {}),
    }


AS_CHARTS = ["cartesian", "null_cartesian", "null_cylindrical"]


# -- Khan-Penrose ------------------------------------------------------------------------

KP_NULL_LINE = (
    "ds^2 = -\\dfrac{2L^2\\left(1 - u^2 - v^2\\right)^{3/2}}{\\sqrt{1 - u^2}\\sqrt{1 - v^2}"
    "\\left(uv + \\sqrt{1 - u^2}\\sqrt{1 - v^2}\\right)^2}du\\,dv"
    " + \\left(1 - u^2 - v^2\\right)\\dfrac{1 + u\\sqrt{1 - v^2} + v\\sqrt{1 - u^2}}"
    "{1 - u\\sqrt{1 - v^2} - v\\sqrt{1 - u^2}}dx^2"
    " + \\left(1 - u^2 - v^2\\right)\\dfrac{1 - u\\sqrt{1 - v^2} - v\\sqrt{1 - u^2}}"
    "{1 + u\\sqrt{1 - v^2} + v\\sqrt{1 - u^2}}dy^2")
KP_COSMOLOGICAL_LINE = (
    "ds^2 = \\dfrac{L^2\\left(\\cos\\tau\\right)^{3/2}}{2\\sqrt{\\cos\\sigma}}\\left(-d\\tau^2 + d\\sigma^2\\right)"
    " + \\dfrac{\\cos\\sigma}{\\cos\\tau}\\left(\\left(1 + \\sin\\tau\\right)^2dx^2"
    " + \\left(1 - \\sin\\tau\\right)^2dy^2\\right)")


def _factor_powers(expr, factor=True):
    """(numeric coefficient, {base: exponent}) of a product, factored first unless told not to."""
    coefficient, powers = sp.Integer(1), {}
    for f in sp.Mul.make_args(sp.factor(expr) if factor else expr):
        base, k = (f.base, f.exp) if f.is_Pow else (f, sp.Integer(1))
        if base.is_Number:
            coefficient *= f
        else:
            powers[base] = powers.get(base, 0) + k
    return coefficient, powers


def _merge_pair(powers, one, other, product):
    """one^p other^q with one * other = product: the common power written as product^k, and
    the rest left as it is. Returns k, so the caller can carry a sign."""
    p, q = powers.pop(one, 0), powers.pop(other, 0)
    k = min(p, q) if p > 0 and q > 0 else max(p, q) if p < 0 and q < 0 else 0
    for base, e in ((one, p - k), (other, q - k), (product, k)):
        if e != 0:
            powers[base] = powers.get(base, 0) + e
    return k


def trig_pairs(expr, pairs):
    """expr factored, with each (s, c, plus, minus) of pairs, where c^2 = 1 - s^2: every pair
    (1 + s)(1 - s) written as c^2, every 1 + s or 1 - s left in a denominator cleared by its
    conjugate, the survivors written as the placeholders plus and minus, and every factor that
    is a sum written with s^2 as 1 - c^2 wherever that has fewer terms. So the cosmological
    chart's g_xx prints as its line element writes it, (1 + sin tau)^2 cos sigma/cos tau."""
    coefficient, powers = _factor_powers(expr)
    for s, c, plus, minus in pairs:
        p = q = 0
        for base in list(powers):
            if sp.expand(base - (1 + s)) == 0:
                p += powers.pop(base)
            elif sp.expand(base - (1 - s)) == 0:
                q += powers.pop(base)
            elif sp.expand(base + (1 - s)) == 0:
                k = powers.pop(base)
                q += k
                coefficient *= (-1) ** k
        k = min(p, q) if p > 0 and q > 0 else max(p, q) if p < 0 and q < 0 else 0
        p, q, ce = p - k, q - k, 2 * k
        if p < 0:
            q, ce, p = q - p, ce + 2 * p, 0
        if q < 0:
            p, ce, q = p - q, ce + 2 * q, 0
        for base, e in ((plus, p), (minus, q), (c, ce)):
            if e != 0:
                powers[base] = powers.get(base, 0) + e
    sines = [s for s, _, _, _ in pairs]
    for base in list(powers):
        if not base.is_Add or not base.has(*sines):
            continue
        best = base
        for chosen in itertools.product((False, True), repeat=len(pairs)):
            alt = base
            for use, (s, c, _, _) in zip(chosen, pairs):
                if use and alt.has(s):
                    alt = sp.expand(sum(k * s ** (m[0] % 2) * (1 - c ** 2) ** (m[0] // 2)
                                        for m, k in sp.Poly(alt, s).terms()))
            if len(sp.Add.make_args(alt)) < len(sp.Add.make_args(best)):
                best = alt
        if best is not base:
            e = powers.pop(base)
            b_coefficient, b_powers = _factor_powers(best)
            coefficient *= b_coefficient ** e
            for b, k in b_powers.items():
                powers[b] = powers.get(b, 0) + k * e
    return _keep_coeff(coefficient, sp.Mul(*[b ** k for b, k in powers.items()]))


class KhanPenroseForms:
    """The double null chart's values built from the cosmological chart's by the law of
    transformation, and a pretty printer that prints a computed value as the one it equals.

    With u = sin(alpha) and v = sin(beta), tau = alpha + beta and sigma = alpha - beta, so
    cos tau = sqrt(1 - u^2) sqrt(1 - v^2) - uv, cos sigma = sqrt(1 - u^2) sqrt(1 - v^2) + uv,
    sin tau = u sqrt(1 - v^2) + v sqrt(1 - u^2) and sin sigma = u sqrt(1 - v^2) - v sqrt(1 - u^2),
    and the cosmological values, compact in these, carry over through the Jacobian
    d(tau, sigma)/d(u, v), whose entries are 1/sqrt(1 - u^2) and +-1/sqrt(1 - v^2), and for the
    Christoffel symbols the second derivatives u/(1 - u^2)^{3/2} and +-v/(1 - v^2)^{3/2}. Each of
    those four is a named placeholder, and so are 1 +- sin tau and 1 +- sin sigma. A half power
    of cos tau comes with one of cos sigma and their product is 1 - u^2 - v^2, and a sum of at
    most six terms in u, v and the two roots is written out and factored over them.

    The printed value is matched to its form at two points to 1e-25 and read back by the chart
    exactly, so a form that differs is refused twice. The checker reads sqrt(1 - u^2) as
    i sqrt(u - 1) sqrt(u + 1), which is -sqrt(1 - u^2) in principal roots, so the value is read
    in the real roots before it is matched."""

    POINTS = ((sp.Rational(3, 10), sp.Rational(11, 20), sp.Rational(7, 5)),
              (sp.Rational(-1, 7), sp.Rational(2, 5), sp.Rational(9, 4)))

    def __init__(self, null_symbols, L):
        u, v = null_symbols[:2]
        self.u, self.v, self.L = u, v, L
        Cm, Cp, S, Sp, Tp, Tm, Up, Um = sp.symbols("KPCm KPCp KPS KPSp KPTp KPTm KPUp KPUm", positive=True)
        self.PA, self.PB, self.PC = sp.symbols("KPA KPB KPC", positive=True)
        self.atoms = (Cm, Cp, S, Sp, Tp, Tm, Up, Um)
        a, b = sp.sqrt(1 - u ** 2), sp.sqrt(1 - v ** 2)
        self.value_of = {Cm: a * b - u * v, Cp: a * b + u * v, S: u * b + v * a, Sp: u * b - v * a,
                         Tp: 1 + u * b + v * a, Tm: 1 - u * b - v * a, Up: 1 + u * b - v * a, Um: 1 - u * b + v * a,
                         self.PA: 1 - u ** 2, self.PB: 1 - v ** 2, self.PC: 1 - u ** 2 - v ** 2}
        roots = "\\sqrt{1 - u^2}\\sqrt{1 - v^2}"
        self.named = {Cm: roots + " - u\\,v", Cp: roots + " + u\\,v", S: "u\\sqrt{1 - v^2} + v\\sqrt{1 - u^2}",
                      Sp: "u\\sqrt{1 - v^2} - v\\sqrt{1 - u^2}", Tp: "1 + u\\sqrt{1 - v^2} + v\\sqrt{1 - u^2}",
                      Tm: "1 - u\\sqrt{1 - v^2} - v\\sqrt{1 - u^2}", Up: "1 + u\\sqrt{1 - v^2} - v\\sqrt{1 - u^2}",
                      Um: "1 - u\\sqrt{1 - v^2} + v\\sqrt{1 - u^2}", self.PA: "1 - u^2", self.PB: "1 - v^2",
                      self.PC: "1 - u^2 - v^2"}
        self.pairs = [(S, Cm, Tp, Tm), (Sp, Cp, Up, Um)]
        self.real_roots = {sp.sqrt(u - 1): -sp.I * a / sp.sqrt(u + 1), sp.sqrt(v - 1): -sp.I * b / sp.sqrt(v + 1),
                           sp.sqrt(u ** 2 + v ** 2 - 1): -sp.I * sp.sqrt(1 - u ** 2 - v ** 2)}
        cosmological = cp.Chart(["\\tau", "\\sigma", "x", "y"], ["L"], KP_COSMOLOGICAL_LINE)
        tau, sigma = cosmological.symbols[:2]
        to_atoms = {sp.sin(tau): S, sp.cos(tau): Cm, sp.sin(sigma): Sp, sp.cos(sigma): Cp,
                    cosmological.reader.parameters["L"]: L}
        geo = cosmological.geo

        def at(e):
            return sp.sympify(e).subs(to_atoms)
        P, Q, n, zero = sp.sqrt(self.PA), sp.sqrt(self.PB), 4, sp.Integer(0)
        # J[a][b] = d(tau, sigma, x, y)^a/d(u, v, x, y)^b, K its inverse, H the second derivatives.
        J = [[zero] * n for _ in range(n)]
        J[0][0], J[0][1], J[1][0], J[1][1], J[2][2], J[3][3] = 1 / P, 1 / Q, 1 / P, -1 / Q, 1, 1
        K = [[zero] * n for _ in range(n)]
        K[0][0], K[0][1], K[1][0], K[1][1], K[2][2], K[3][3] = P / 2, P / 2, Q / 2, -Q / 2, 1, 1
        H = [[[zero] * n for _ in range(n)] for _ in range(n)]
        H[0][0][0], H[0][1][1], H[1][0][0], H[1][1][1] = u / P ** 3, v / Q ** 3, u / P ** 3, -v / Q ** 3
        N = range(n)
        g = [[at(geo.g[i, j]) for j in N] for i in N]
        gi = [[at(geo.ginv[i, j]) for j in N] for i in N]
        G = [[[at(c) for c in row] for row in plane] for plane in geo.christoffel_ull()]
        lower = geo.riemann_llll()
        upper = geo.raise_indices(lower, 4, (0,))
        gn = [[sum(J[a][b] * J[d][c] * g[a][d] for a in N for d in N) for c in N] for b in N]
        gin = [[sum(K[b][a] * K[c][d] * gi[a][d] for a in N for d in N) for c in N] for b in N]
        Gn = [[[sum(K[a][d] * (sum(G[d][e][f] * J[e][b] * J[f][c] for e in N for f in N) + H[d][b][c]) for d in N)
                for c in N] for b in N] for a in N]
        Gl = [[[sum(gn[a][d] * Gn[d][b][c] for d in N) for c in N] for b in N] for a in N]
        forms = [T[i][j] for T in (gn, gin) for i in N for j in N]
        # A geodesic equation carries twice each Christoffel symbol with two different indices.
        forms += [T[i][j][k] * m for T in (Gn, Gl) for i in N for j in N for k in N for m in (1, 2)]
        down = [[i for i in N if J[i][b] != 0] for b in N]
        for a, b, c, d in itertools.product(N, repeat=4):
            forms.append(sum(J[p][a] * J[q][b] * J[r][c] * J[s][d] * at(vm._at(lower, (p, q, r, s)))
                             for p in down[a] for q in down[b] for r in down[c] for s in down[d]))
            forms.append(sum(K[a][p] * J[q][b] * J[r][c] * J[s][d] * at(vm._at(upper, (p, q, r, s)))
                             for p in N if K[a][p] != 0 for q in down[b] for r in down[c] for s in down[d]))
        forms.append(at(geo.kretschmann()))
        self.forms = [(self.numbers(f, self.value_of), f) for f in {sp.together(f) for f in forms if f != 0}]

    def numbers(self, e, values):
        e = e.subs(values)
        return [complex(sp.N(e.subs({self.u: pu, self.v: pv, self.L: pL}), 40)) for pu, pv, pL in self.POINTS]

    def pretty(self, value):
        if value == 0:
            return value
        want = self.numbers(value, self.real_roots)
        for got, form in self.forms:
            for sign in (1, -1):
                if all(abs(w - sign * m) <= 1e-25 * (1 + abs(w)) for w, m in zip(want, got)):
                    return self.shape(sign * form)
        raise AssertionError(f"khan_penrose: no transformed form equals {value}")

    def shape(self, form):
        Cm, Cp = self.atoms[:2]
        coefficient, powers = _factor_powers(trig_pairs(form, self.pairs), factor=False)
        km, kp = powers.get(Cm, 0), powers.get(Cp, 0)
        if (2 * km) % 2 or (2 * kp) % 2:
            if not ((2 * km) % 2 and (2 * kp) % 2):
                raise AssertionError(f"khan_penrose: an unpaired half power of a cosine in {form}")
            powers[Cm], powers[Cp] = km - sp.Rational(1, 2), kp - sp.Rational(1, 2)
            powers[self.PC] = powers.get(self.PC, 0) + sp.Rational(1, 2)
        _merge_pair(powers, Cm, Cp, self.PC)
        rest = []
        for base, k in powers.items():
            written = self.written_out(base) if base.is_Add else None
            if written is None:
                rest.append(base ** k)
                continue
            w_coefficient, w_powers = written
            coefficient *= w_coefficient ** k
            rest += [b ** (e * k) for b, e in w_powers]
        return _keep_coeff(coefficient, sp.Mul(*rest))

    def written_out(self, base):
        """A sum in u, v, sqrt(1 - u^2) and sqrt(1 - v^2), factored over them, as (coefficient,
        [(base, exponent)]), or None where it has more than six terms."""
        u, v = self.u, self.v
        a, b = sp.symbols("KPa KPb", positive=True)
        values = {s: self.value_of[s].subs({sp.sqrt(1 - u ** 2): a, sp.sqrt(1 - v ** 2): b}) for s in self.atoms}
        values.update({self.PA: a ** 2, self.PB: b ** 2, self.PC: 1 - u ** 2 - v ** 2})
        e = sp.expand(base.subs(values))
        e = sp.expand(sum(k * (1 - u ** 2) ** (i // 2) * a ** (i % 2) * (1 - v ** 2) ** (j // 2) * b ** (j % 2)
                          for (i, j), k in sp.Poly(e, a, b).terms()))
        if len(sp.Add.make_args(e)) > 6:
            return None
        coefficient, powers = _factor_powers(e)
        for one, other, name in ((u + 1, u - 1, self.PA), (v + 1, v - 1, self.PB)):
            coefficient *= (-1) ** _merge_pair(powers, one, other, name)
        for base in list(powers):
            if sp.expand(base - (u ** 2 + v ** 2 - 1)) == 0:
                k = powers.pop(base)
                powers[self.PC] = powers.get(self.PC, 0) + k
                coefficient *= (-1) ** k
        back = {a: sp.sqrt(self.PA), b: sp.sqrt(self.PB)}
        return coefficient, [(base.subs(back), k) for base, k in powers.items()]


def khan_penrose():
    """The Khan-Penrose spacetime where both waves have passed, in two charts. The first is
    Khan and Penrose's null chart, in the form J. Frauendiener, C. Stevens and B. Whale give as
    their eq. (11), Phys. Rev. D 89, 104026 (2014), with the signature flipped and the focal
    length L restored on the plane of u and v; its g_xx factor (R + Q)(W + P)/((R - Q)(W - P))
    is written as (1 + S)/(1 - S) with S = u sqrt(1 - v^2) + v sqrt(1 - u^2), which it equals.
    The second is the chart of tau = arcsin u + arcsin v and sigma = arcsin u - arcsin v, whose
    sines are the t and z of J. B. Griffiths and M. Santano-Roco, eq. (42), Class. Quantum
    Grav. 19, 4273 (2002); before it is written, its metric pulled back through that map is
    checked equal to the null chart's at random points.

    The cosmological chart's values are printed by trig_pairs, and the null chart's as the
    transformed cosmological ones by KhanPenroseForms, whose docstring says how."""
    null_coords = ["u", "v", "x", "y"]
    parameters = ["L"]
    probe = vm.Reader(null_coords, parameters, ())
    u, v, L = probe.symbol["u"], probe.symbol["v"], probe.parameters["L"]
    forms = KhanPenroseForms([u, v], L)
    domains = ["x \\in (-\\infty, \\infty)", "y \\in (-\\infty, \\infty)"]
    null = {
        "metric_id": "khan_penrose",
        "system": {"id": "double_null", "name": "Double Null", "coords": null_coords,
                   "domains": ["u \\in [0, 1)", "v \\in [0, \\sqrt{1 - u^2})"] + domains,
                   "parameters": parameters, "line_element": KP_NULL_LINE},
        "chart_line_element": KP_NULL_LINE,
        "printer": {"lead": [L, u, v], "named": forms.named},
        "pretty": forms.pretty,
    }
    cosmological_coords = ["\\tau", "\\sigma", "x", "y"]
    cr = vm.Reader(cosmological_coords, parameters, ())
    tau, sigma = cr.symbol["\\tau"], cr.symbol["\\sigma"]
    Tp, Tm, Up, Um = sp.symbols("KPTp KPTm KPUp KPUm", positive=True)
    pairs = [(sp.sin(tau), sp.cos(tau), Tp, Tm), (sp.sin(sigma), sp.cos(sigma), Up, Um)]
    cosmological = {
        "metric_id": "khan_penrose",
        "system": {"id": "cosmological", "name": "Cosmological", "coords": cosmological_coords,
                   "domains": ["\\tau \\in [0, \\pi/2)", "\\sigma \\in [-\\tau, \\tau]"] + domains,
                   "parameters": parameters, "line_element": KP_COSMOLOGICAL_LINE},
        "chart_line_element": KP_COSMOLOGICAL_LINE,
        "printer": {"lead": [cr.parameters["L"]], "named": {Tp: "1 + \\sin\\tau", Tm: "1 - \\sin\\tau",
                                                             Up: "1 + \\sin\\sigma", Um: "1 - \\sin\\sigma"}},
        "pretty": lambda value: trig_pairs(value, pairs),
        "check": khan_penrose_pullback,
    }
    return [null, cosmological]


def khan_penrose_pullback(chart, points=24):
    """J^T g J, with g the cosmological metric and J the Jacobian of tau = arcsin u + arcsin v,
    sigma = arcsin u - arcsin v, against the null chart's metric, at random points where both
    waves have passed, to forty digits. Both line elements are read as printed, with sympy's
    principal roots, which are the real roots there."""
    null_coords, cosmological_coords = ["u", "v", "x", "y"], ["\\tau", "\\sigma", "x", "y"]
    rn, rc = vm.Reader(null_coords, ["L"], ()), vm.Reader(cosmological_coords, ["L"], ())

    def matrix(reader, line, coords):
        form = sp.expand(reader(line.partition("=")[2]))
        d = [reader.differential[name] for name in coords]
        return sp.Matrix(4, 4, lambda i, j: form.coeff(d[i], 2) if i == j
                         else form.coeff(d[i], 1).coeff(d[j], 1) / 2)

    gn, gc = matrix(rn, KP_NULL_LINE, null_coords), matrix(rc, KP_COSMOLOGICAL_LINE, cosmological_coords)
    u, v, x, y = (rn.symbol[n] for n in null_coords)
    image = [sp.asin(u) + sp.asin(v), sp.asin(u) - sp.asin(v), x, y]
    J = sp.Matrix(4, 4, lambda i, j: sp.diff(image[i], [u, v, x, y][j]))
    at = dict(zip((rc.symbol[n] for n in cosmological_coords), image))
    at[rc.parameters["L"]] = rn.parameters["L"]
    difference = J.T * gc.subs(at) * J - gn
    rng = __import__("random").Random(9)
    checked = 0
    while checked < points:
        p = {u: sp.Rational(rng.randint(1, 999), 1000), v: sp.Rational(rng.randint(1, 999), 1000),
             rn.parameters["L"]: sp.Rational(rng.randint(1, 50), 10)}
        if p[u] ** 2 + p[v] ** 2 >= 1:
            continue
        worst = max(abs(sp.N(difference[i, j].subs(p), 40)) for i in range(4) for j in range(4))
        if worst > sp.Float("1e-30"):
            raise AssertionError(f"khan_penrose: the cosmological metric pulled back misses the null chart's "
                                 f"by {worst} at {p}")
        checked += 1


CHARTS = {"tov": tov, "malament_hogarth": malament_hogarth, "mixmaster": mixmaster, "lentz": lentz, "godel": godel,
          "einstein_static": [lambda s=s: einstein_static(s) for s in ("hyperspherical", "static_areal", "einstein_cartesian")],
          "btz": [lambda: btz_stationary(), lambda: btz_null(1), lambda: btz_null(-1)],
          "c_metric": c_metric,
          "schwarzschild_de_sitter": [lambda s=s: schwarzschild_de_sitter(s) for s in SDS_CHARTS],
          "milne": [lambda s=s: milne(s) for s in ("comoving_hyperbolic", "comoving_spherical", "logarithmic_time", "inertial")],
          "einstein_rosen_waves": [lambda s=s: einstein_rosen(s) for s in ("cylindrical", "null")],
          "nariai": [lambda s=s: nariai(s) for s in ("static", "global", "conformal")],
          "aichelburg_sexl": [lambda s=s: aichelburg_sexl(s) for s in AS_CHARTS],
          "khan_penrose": khan_penrose,
          "global_monopole": [lambda s=s: global_monopole(s) for s in GM_CHARTS]}


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
        # Both are printed texts, so both already carry a named time as c times it.
        for text, printed in zip(spec["geodesics"], math["geodesics"], strict=True):
            if vm.norm(chart.reader(text.partition("=")[0]) - chart.reader(printed.partition("=")[0])) != 0:
                raise AssertionError(f"the geodesic {text!r} does not read back as {printed!r}")
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
