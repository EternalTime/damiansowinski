#!/usr/bin/env python3
"""Compute and write the coordinate systems whose mathematics is printed by machine: the
charts of tov, malament_hogarth, mixmaster, lentz, einstein_static, btz, c_metric,
schwarzschild_de_sitter, schwarzschild_ads, milne, einstein_rosen_waves, nariai, aichelburg_sexl,
khan_penrose, global_monopole, domain_wall, majumdar_papapetrou, melvin, thin_shell_wormhole, levi_civita, curzon_chazy,
robinson_trautman, string_black_hole, mcvittie, tangherlini, gott_time_machine, zipoy_voorhees, szekeres,
kaluza_klein_monopole, bell_szekeres, spinning_string, photon_rocket, light_beam, wormhole_time_machine,
damour_solodukhin, ori_time_machine, reissner_nordstrom_de_sitter, string_wave, simpson_visser, bardeen, hayward and fisher_jnw, and
Godel's cylindrical chart.

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
malament_hogarth.md, mixmaster.md, lentz.md, godel.md, btz.md, schwarzschild_de_sitter.md,
majumdar_papapetrou.md, robinson_trautman.md, tangherlini.md, szekeres.md, photon_rocket.md and fisher_jnw.md beside this file.
"""
import argparse
import itertools
import json
import random
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
        # Q leads each sum and is left out of `factors`, so that in a product the bracket it stands
        # for closes the term, r(r - 2m)((d_r Phi)^2 + d_r^2 Phi), as a bracketed sum does.
        "printer": {"lead": [Q, dPhi, r, m, sp.Derivative(m, r)], "factors": [dPhi, r, m, sp.Derivative(m, r)],
                    "collect": collect,
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
            # Both are printed texts, so both already carry a named time as c times it.
            if side != original and vm.norm(chart.reader(side) - chart.reader(original)) != 0:
                raise AssertionError(f"the rewrite {side!r} does not read back as {original!r}")
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


# -- Schwarzschild-anti-de Sitter --------------------------------------------------------

def schwarzschild_ads(system_id):
    """Hawking and Page's black hole in anti-de Sitter space: Kottler's metric with
    Lambda = -3/L^2, f = 1 - r_s/r + r^2/L^2, in the static chart and in the two
    Eddington-Finkelstein charts built on its tortoise coordinate, dr_*/dr = 1/f. The parameters
    are Schwarzschild's r_s and anti-de Sitter's L, so each chart reduces to Schwarzschild's as
    L grows without bound and to anti-de Sitter's static chart at r_s = 0. Every value is
    printed around L^2 r f = r^3 + L^2 r - L^2 r_s, the cubic whose positive root is the horizon. The metric and
    its inverse are written as the line element writes f, and the Kretschmann scalar as
    Schwarzschild's 12r_s^2/r^6 plus anti-de Sitter's 24/L^4, which it is."""
    f = "\\left(1 - \\dfrac{r_s}{r} + \\dfrac{r^2}{L^2}\\right)"
    bare = "1 - \\dfrac{r_s}{r} + \\dfrac{r^2}{L^2}"
    sphere = " + r^2\\left(d\\theta^2 + \\sin^2\\theta\\,d\\phi^2\\right)"
    domains = ["r \\in (0, \\infty)", "\\theta \\in [0, \\pi]", "\\phi \\in [0, 2\\pi)",
               "r = r_h \\;\\text{(the horizon)}"]
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
    parameters = ["r_s", "L"]
    probe = vm.Reader(coords, parameters, ())
    r, rs, L = probe.symbol["r"], probe.parameters["r_s"], probe.parameters["L"]
    return {
        "metric_id": "schwarzschild_ads",
        "system": {"id": system_id, "name": name, "coords": coords,
                   "domains": [coords[0] + " \\in (-\\infty, \\infty)"] + domains,
                   "parameters": parameters, "line_element": line},
        "chart_line_element": chart_line,
        "printer": {"rising": [L, rs], "lead": [L, r, rs], "flip": False},
        "components": {"metric_components": metric, "inverse_metric_components": inverse},
        "ricci_scalar": "-\\dfrac{12}{L^2}",
        "kretschmann": "\\dfrac{12r_s^2}{r^6} + \\dfrac{24}{L^4}",
    }


SADS_CHARTS = ["static", "eddington_finkelstein_outgoing", "eddington_finkelstein_ingoing"]


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
    horizon = ["r = r_s/(1 - \\Delta) \\;\\text{(the horizon, for}\\; r_s > 0\\text{)}"]
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


# -- McVittie ----------------------------------------------------------------------------

MCV_SPHERE = "\\left(d\\theta^2 + \\sin^2\\theta\\,d\\phi^2\\right)"
MCV_ANGLES = ["\\theta \\in [0, \\pi]", "\\phi \\in [0, 2\\pi)"]


def mcvittie(system_id):
    """McVittie's point mass in a spatially flat expanding universe, in his own isotropic
    comoving chart and in the areal chart R = ar(1 + r_s/4ar)^2 of Nolan and of Kaloper,
    Kleban and Martin, whose metric depends on the scale factor only through the Hubble
    rate H = (da/dt)/a, a frequency as de Sitter's flat slicing has it. Every curvature value is collected by the derivatives of a, or by
    H and its derivative, so that each reads as its Schwarzschild value plus the terms the
    expansion adds, and the Ricci scalar is written as 12H^2/c^2 plus the term in Hdot that
    diverges on the singular sphere R = r_s. The areal chart is checked, slot by slot, to be
    the isotropic chart pulled back; mcvittie.md derives both."""
    if system_id == "isotropic":
        coords = ["t", "r", "\\theta", "\\phi"]
        parameters = ["r_s", "a = a(t)"]
        mu = "\\dfrac{r_s}{4ar}"
        lapse = "\\left(\\dfrac{1 - " + mu + "}{1 + " + mu + "}\\right)^2"
        space = "a^2\\left(1 + " + mu + "\\right)^4\\left(dr^2 + r^2" + MCV_SPHERE + "\\right)"
        line = "ds^2 = -" + lapse + "c^2dt^2 + " + space
        chart_line = "ds^2 = -" + lapse + "dt^2 + " + space
        probe = vm.Reader(coords, parameters, ())
        r, rs, a = probe.symbol["r"], probe.parameters["r_s"], probe.parameters["a"]
        t = probe.symbol["t"]
        rates = [sp.Derivative(a, (t, 2)), sp.Derivative(a, t)]
        conformal = "a^2\\left(1 + " + mu + "\\right)^4"
        inverse = "\\dfrac{1}{a^2}\\left(1 + " + mu + "\\right)^{-4}"
        extra = {
            "components": {
                "metric_components": {("t", "t"): "-\\left(\\dfrac{4ar - r_s}{4ar + r_s}\\right)^2",
                                      ("r", "r"): conformal, ("\\theta", "\\theta"): conformal[:3] + "r^2" + conformal[3:],
                                      ("\\phi", "\\phi"): conformal[:3] + "r^2\\sin^2\\theta" + conformal[3:]},
                "inverse_metric_components": {("t", "t"): "-\\left(\\dfrac{4ar + r_s}{4ar - r_s}\\right)^2",
                                              ("r", "r"): inverse,
                                              ("\\theta", "\\theta"): "\\dfrac{1}{a^2r^2}\\left(1 + " + mu + "\\right)^{-4}",
                                              ("\\phi", "\\phi"): "\\dfrac{1}{a^2r^2\\sin^2\\theta}\\left(1 + " + mu + "\\right)^{-4}"}},
            "ricci_scalar": ("\\dfrac{12\\dot{a}^2}{a^2} + 6\\left(\\dfrac{\\ddot{a}}{a} - \\dfrac{\\dot{a}^2}{a^2}\\right)"
                             "\\dfrac{4ar + r_s}{4ar - r_s}"),
            # The Weyl tensor's 12r_s^2/R^6 at the areal radius R = ar(1 + r_s/4ar)^2, and the
            # Ricci tensor's share, 12H^4 + 12(H^2 + Hdot(1 + mu)/(1 - mu))^2.
            "kretschmann": ("\\dfrac{12r_s^2}{a^6r^6}\\left(1 + " + mu + "\\right)^{-12} + \\dfrac{12\\dot{a}^4}{a^4}"
                            " + 12\\left(\\dfrac{\\dot{a}^2}{a^2} + \\left(\\dfrac{\\ddot{a}}{a} - \\dfrac{\\dot{a}^2}{a^2}\\right)"
                            "\\dfrac{4ar + r_s}{4ar - r_s}\\right)^2"),
        }
        domains = (["t \\in (0, \\infty)", "r \\in (r_s/(4a), \\infty)"] + MCV_ANGLES
                   + ["r = r_s/(4a) \\;\\text{(curvature singularity unless}\\; \\dot{a}/a \\;\\text{is constant)}"])
        name = "Isotropic Comoving"
        printer = {"lead": [a, r, rs], "dotted": ["a"],
                   "collect": lambda poly, printer: cp.collect_by(poly, rates, printer)}
    else:
        coords = ["t", "R", "\\theta", "\\phi"]
        parameters = ["r_s", "H = H(t)"]
        f = "\\left(1 - \\dfrac{r_s}{R} - \\dfrac{H^2R^2}{c^2}\\right)"
        root = "\\sqrt{1 - \\dfrac{r_s}{R}}"
        rest = " + \\dfrac{dR^2}{1 - \\dfrac{r_s}{R}} + R^2" + MCV_SPHERE
        line = "ds^2 = -" + f + "c^2dt^2 - \\dfrac{2HR}{" + root + "}dt\\,dR" + rest
        chart_line = "ds^2 = -" + f + "dt^2 - \\dfrac{2HR}{c" + root + "}dt\\,dR" + rest
        probe = vm.Reader(coords, parameters, ())
        R, rs, H = probe.symbol["R"], probe.parameters["r_s"], probe.parameters["H"]
        rates = [sp.Derivative(H, probe.symbol["t"]), H]
        extra = {
            "check": mcvittie_areal,
            "components": {
                "metric_components": {("t", "t"): "-" + f, ("R", "R"): "\\left(1 - \\dfrac{r_s}{R}\\right)^{-1}",
                                      ("t", "R"): "-\\dfrac{HR}{c" + root + "}", ("R", "t"): "-\\dfrac{HR}{c" + root + "}"},
                "inverse_metric_components": {("t", "t"): "-\\left(1 - \\dfrac{r_s}{R}\\right)^{-1}",
                                              ("t", "R"): "-\\dfrac{HR}{c" + root + "}",
                                              ("R", "t"): "-\\dfrac{HR}{c" + root + "}",
                                              ("R", "R"): f[6:-7]}},
            "ricci_scalar": "\\dfrac{12H^2}{c^2} + \\dfrac{6\\dot{H}}{c" + root + "}",
            # Schwarzschild's 12r_s^2/R^6 from the Weyl tensor, and the Ricci tensor's share,
            # which is de Sitter's 24H^4/c^4 where Hdot = 0.
            "kretschmann": ("\\dfrac{12r_s^2}{R^6} + \\dfrac{12H^4}{c^4}"
                            " + 12\\left(\\dfrac{H^2}{c^2} + \\dfrac{\\dot{H}}{c" + root + "}\\right)^2"),
        }
        domains = (["t \\in (0, \\infty)", "R \\in (r_s, \\infty)"] + MCV_ANGLES
                   + ["R = r_s \\;\\text{(curvature singularity where}\\; \\dot{H} \\neq 0\\text{)}",
                      "1 - r_s/R - H^2R^2/c^2 = 0 \\;\\text{(apparent horizons)}"])
        name = "Areal Radius"
        printer = {"lead": [H, R, rs, probe.c], "factors": [probe.c, H, R, rs], "dotted": ["H"],
                   "collect": lambda poly, printer: cp.collect_by(poly, rates, printer)}
    return {
        "metric_id": "mcvittie",
        "system": {"id": system_id, "name": name, "coords": coords, "domains": domains,
                   "parameters": parameters, "line_element": line},
        "chart_line_element": chart_line,
        "printer": printer,
        **extra,
    }


def mcvittie_areal(chart):
    """J^T g J, with g the areal chart's metric at R = ar(1 + r_s/4ar)^2 and H = c adot/a, and J
    the Jacobian of (t, R) with respect to (t, r), against the isotropic chart's metric, in
    every slot. The areal metric is first written with a symbol S for sqrt(1 - r_s/R) and
    checked against the chart's own; S is then (4ar - r_s)/(4ar + r_s), which is positive on
    the isotropic chart's domain r > r_s/4a."""
    spec = mcvittie("isotropic")
    target = cp.Chart(spec["system"]["coords"], spec["system"]["parameters"], spec["chart_line_element"])
    t, r = target.symbols[:2]
    a, rs = target.reader.parameters["a"], target.reader.parameters["r_s"]
    T, R = chart.symbols[:2]
    H, Rs = chart.reader.parameters["H"], chart.reader.parameters["r_s"]
    S = sp.Symbol("S", positive=True)
    c = chart.reader.c
    f = 1 - Rs / R - H ** 2 * R ** 2 / c ** 2
    areal = sp.Matrix([[-f, -H * R / (c * S)], [-H * R / (c * S), 1 / (1 - Rs / R)]])
    for i in range(2):
        for j in range(2):
            if vm.norm(areal[i, j].subs(S, sp.sqrt(1 - Rs / R)) - chart.geo.g[i, j]) != 0:
                raise AssertionError("mcvittie: the areal metric written with S misses the chart's own")
    radius = a * r * (1 + rs / (4 * a * r)) ** 2
    at = {S: (4 * a * r - rs) / (4 * a * r + rs), Rs: rs}
    plane = areal.subs(at).subs(H, c * sp.Derivative(a, t) / a).subs(R, radius).subs(T, t)
    J = sp.Matrix([[1, 0], [sp.diff(radius, t), sp.diff(radius, r)]])
    pulled = J.T * plane * J
    for i in range(2):
        for j in range(i, 2):
            if vm.norm(pulled[i, j] - target.geo.g[i, j]) != 0:
                raise AssertionError(f"mcvittie: the areal chart pulled back misses the isotropic chart "
                                     f"in slot {target.coords_tex[i]}{target.coords_tex[j]}")
    for k in (2, 3):
        if vm.norm(chart.geo.g[k, k].subs(R, radius).subs(chart.symbols[2], target.symbols[2])
                   - target.geo.g[k, k]) != 0:
            raise AssertionError(f"mcvittie: the areal sphere misses the isotropic one in slot {k}")


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


# -- Domain wall -----------------------------------------------------------------------

def domain_wall(system):
    """The four charts of the domain wall of Vilenkin and of Ipser and Sikivie, flat on either
    side of the wall and with a kink across it: the planar chart with the de Sitter world
    volume of the wall in flat slices, the global chart with it in closed slices, the global
    chart in the conformal form of Cvetic and Soleng, and the inertial chart of one side, in
    which the wall is the hyperboloid R^2 - c^2T^2 = 1/k^2. The curvature of the first three
    is a delta on the wall, which the checker reads as a distribution, and their Kretschmann
    scalar, the square of that delta, is stated where the wall is not; domain_wall.md derives
    each chart."""
    sphere = "\\left(d\\theta^2 + \\sin^2\\theta\\,d\\phi^2\\right)"
    charts = {
        "planar": {"coords": ["t", "x", "y", "z"], "kink": "z",
                   "chart": "ds^2 = \\left(1 - k|z|\\right)^2\\left(-dt^2 + e^{2kt}\\left(dx^2 + dy^2\\right)\\right) + dz^2"},
        "global": {"coords": ["t", "z", "\\theta", "\\phi"], "kink": "z",
                   "chart": "ds^2 = \\left(1 - k|z|\\right)^2\\left(-dt^2 + \\dfrac{\\cosh^2(kt)}{k^2}" + sphere + "\\right) + dz^2"},
        "conformal": {"coords": ["t", "w", "\\theta", "\\phi"], "kink": "w",
                      "chart": "ds^2 = e^{-2k|w|}\\left(-dt^2 + dw^2 + \\dfrac{\\cosh^2(kt)}{k^2}" + sphere + "\\right)"},
        "inertial": {"coords": ["T", "R", "\\theta", "\\phi"],
                     "chart": "ds^2 = -dT^2 + dR^2 + R^2" + sphere},
    }
    chart = charts[system]
    path = METRICS / "domain_wall.json"
    published = next(c for c in json.loads(path.read_text(encoding="utf-8"))["coordinates"] if c["id"] == system)
    probe = vm.Reader(chart["coords"], ["k"], ())
    k, c, time_ = probe.parameters["k"], probe.c, probe.symbol[chart["coords"][0]]
    spec = {
        "metric_id": "domain_wall",
        "system": {"id": system, "name": published["name"], "coords": chart["coords"], "domains": published["domains"],
                   "parameters": ["k"], "line_element": published["line_element"]},
        "chart_line_element": chart["chart"],
        "time": chart["coords"][0],
        "printer": {"lead": [k, c, time_]},
    }
    if "kink" in chart:
        x = probe.symbol[chart["kink"]]
        then = cp.hyperbolic(k * c * time_) if system != "planar" else sp.factor
        pretty, overrides = cp.kink(x, then)
        absolute = next(p for p in overrides if p.name.startswith("_abs"))
        spec["pretty"] = pretty
        spec["printer"] = {"lead": [k, c, time_], "rising": [absolute], "overrides": overrides}
        spec["kretschmann_where"] = chart["kink"]
    # kct and k|z| are set as the line element sets them, and sinh over cosh as tanh.
    hyperbolic_tangent = [(f"\\dfrac{{{n}k\\sinh\\left(kct\\right)}}{{\\cosh\\left(kct\\right)}}",
                           f"{n}k\\tanh\\left(kct\\right)") for n in ("", "2")]
    spec["rewrite"] = [("k\\,c\\,t", "kct"), ("k\\,c\\,T", "kcT"), ("k\\,|", "k|")] + hyperbolic_tangent
    return spec


DW_CHARTS = ["planar", "global", "conformal", "inertial"]


# -- Majumdar-Papapetrou -----------------------------------------------------------------

def majumdar_papapetrou(system_id):
    """Majumdar's and Papapetrou's static metric, -U^{-2}c^2dt^2 + U^2 times flat space, in
    the three charts its literature uses: the Cartesian chart of Majumdar, Papapetrou and Hartle
    and Hawking with U = U(x,y,z) left free, the cylindrical chart of holes strung along one
    axis with U = U(rho,z), and the isotropic chart of one hole, U = 1 + m/r. The free U is left
    free in every tensor, so no component assumes Laplace's equation, which the Einstein-Maxwell
    equations impose away from the holes, as majumdar_papapetrou_check confirms before anything is written;
    majumdar_papapetrou.md records the choices."""
    flat = {"cartesian": "dx^2 + dy^2 + dz^2", "cylindrical": "d\\rho^2 + \\rho^2d\\phi^2 + dz^2",
            "isotropic": "dr^2 + r^2\\left(d\\theta^2 + \\sin^2\\theta\\,d\\phi^2\\right)"}[system_id]
    if system_id == "isotropic":
        U = "\\left(1 + \\dfrac{m}{r}\\right)"
        line = "ds^2 = -" + U + "^{-2}{}dt^2 + " + U + "^2\\left(" + flat + "\\right)"
    else:
        line = "ds^2 = -\\dfrac{{}dt^2}{U^2} + U^2\\left(" + flat + "\\right)"
    if system_id == "cartesian":
        coords, name, parameters = ["t", "x", "y", "z"], "Cartesian", ["U = U(x,y,z)"]
        domains = ["x \\in (-\\infty, \\infty)", "y \\in (-\\infty, \\infty)", "z \\in (-\\infty, \\infty)",
                   "(x, y, z) \\neq \\mathbf{x}_i \\;\\text{(the horizons)}"]
    elif system_id == "cylindrical":
        coords, name, parameters = ["t", "\\rho", "\\phi", "z"], "Cylindrical", ["U = U(\\rho,z)"]
        domains = ["\\rho \\in [0, \\infty)", "\\phi \\in [0, 2\\pi)", "z \\in (-\\infty, \\infty)",
                   "(\\rho, z) \\neq (0, z_i) \\;\\text{(the horizons)}"]
    else:
        coords, name, parameters = ["t", "r", "\\theta", "\\phi"], "Isotropic", ["m"]
        domains = ["r \\in (0, \\infty)", "\\theta \\in [0, \\pi]", "\\phi \\in [0, 2\\pi)",
                   "r = 0 \\;\\text{(the horizon)}"]
    probe = vm.Reader(coords, parameters, ())
    spec = {
        "metric_id": "majumdar_papapetrou",
        "system": {"id": system_id, "name": name, "coords": coords,
                   "domains": ["t \\in (-\\infty, \\infty)"] + domains,
                   "parameters": parameters, "line_element": line.replace("{}dt", "c^2dt")},
        "chart_line_element": line.replace("{}dt", "dt"),
        "check": majumdar_papapetrou_check,
    }
    if system_id == "isotropic":
        r, m = probe.symbol["r"], probe.parameters["m"]
        spec["printer"] = {"lead": [r, m], "factors": [m, r]}
        spec["components"] = {"metric_components": {("t", "t"): "-" + U + "^{-2}", ("r", "r"): U + "^2"},
                              "inverse_metric_components": {("t", "t"): "-" + U + "^2", ("r", "r"): U + "^{-2}"}}
    else:
        Uf = probe.parameters["U"]
        spec["printer"] = {"lead": [Uf], "flip": False,
                           "collect": lambda poly, printer: cp.collect_by(poly, [Uf], printer)}
    return spec


def majumdar_papapetrou_check(chart):
    """The Cartesian chart against the Einstein-Maxwell equations: with A = U^{-1} dt in units
    where G = c = 4 pi epsilon_0 = 1, the Einstein tensor equals 2(F_ma F_n^a - g_mn F^2/4) once
    d_z^2 U is replaced by -d_x^2 U - d_y^2 U, and Maxwell's equations div F = 0 reduce to
    Laplace's for U. Every other chart against the Cartesian one."""
    if chart.coords_tex[1] == "x":
        X, U = chart.symbols, chart.reader.parameters["U"]
        g, gi = chart.geo.g, chart.geo.ginv
        A = [1 / U, 0, 0, 0]
        F = sp.Matrix(4, 4, lambda a, b: sp.diff(A[b], X[a]) - sp.diff(A[a], X[b]))
        Fu = gi * F * gi
        F2 = sum(F[a, b] * Fu[a, b] for a in range(4) for b in range(4))
        stress = sp.Matrix(4, 4, lambda a, b: 2 * (sum(F[a, c] * F[b, d] * gi[c, d] for c in range(4) for d in range(4))
                                                  - g[a, b] * F2 / 4))
        einstein = sp.Matrix(chart.geo.einstein_ll())
        laplace = sp.Derivative(U, (X[1], 2)) + sp.Derivative(U, (X[2], 2)) + sp.Derivative(U, (X[3], 2))
        harmonic = {sp.Derivative(U, (X[3], 2)): sp.Derivative(U, (X[3], 2)) - laplace}
        for a in range(4):
            for b in range(a, 4):
                if vm.norm((einstein[a, b] - stress[a, b]).subs(harmonic)) != 0:
                    raise AssertionError(f"majumdar_papapetrou: G_{a}{b} is not the Maxwell stress for harmonic U")
        root = sp.sqrt(-g.det())
        for b in range(4):
            divergence = sum(sp.diff(root * Fu[a, b], X[a]) for a in range(4)) / root
            if vm.norm(divergence - (laplace / U ** 2 if b == 0 else 0)) != 0:
                raise AssertionError(f"majumdar_papapetrou: Maxwell's equation {b} is not Laplace's for U")
        return
    majumdar_papapetrou_pullback(chart)


def majumdar_papapetrou_pullback(chart):
    """Each chart against the Cartesian one: the Cartesian metric with U = U(rho cos phi,
    rho sin phi, z) pulled back to cylindrical coordinates, and with U = 1 + m/sqrt(x^2 + y^2 + z^2)
    pulled back to spherical ones, minus the chart's own metric, vanish in every slot."""
    cart = ["t", "x", "y", "z"]
    reader = vm.Reader(cart, ["U = U(x,y,z)"], ())
    g = vm.metric_from_line_element(
        reader, "ds^2 = -\\dfrac{dt^2}{U^2} + U^2\\left(dx^2 + dy^2 + dz^2\\right)", cart)
    X = [reader.symbol[c] for c in cart]
    Uc = reader.parameters["U"]
    t, a, b, c = chart.symbols
    if chart.coords_tex[1] == "\\rho":
        image = [t, a * sp.cos(b), a * sp.sin(b), c]
        U_on_chart = chart.reader.parameters["U"]
        g = g.subs(Uc, sp.Function("V")(*X))
        at = {sp.Function("V")(*image): U_on_chart}
    else:
        image = [t, a * sp.sin(b) * sp.cos(c), a * sp.sin(b) * sp.sin(c), a * sp.cos(b)]
        at = {}
        g = g.subs(Uc, 1 + chart.reader.parameters["m"] / sp.sqrt(X[1] ** 2 + X[2] ** 2 + X[3] ** 2))
    J = sp.Matrix(4, 4, lambda i, j: sp.diff(image[i], chart.symbols[j]))
    pulled = (J.T * g.subs(dict(zip(X, image))) * J).subs(at)
    # The radius and rho are positive on the chart, which is what takes sqrt(r^2) to r.
    positive = sp.Symbol("positive_radius", positive=True)
    for i in range(4):
        for j in range(i, 4):
            miss = sp.simplify((pulled[i, j] - chart.geo.g[i, j]).subs(a, positive)).subs(positive, a)
            if vm.norm(miss) != 0:
                raise AssertionError(f"majumdar_papapetrou: the Cartesian metric pulled back misses the "
                                     f"{chart.coords_tex[1]} chart in slot {chart.coords_tex[i]}{chart.coords_tex[j]}")


MP_CHARTS = ["cartesian", "cylindrical", "isotropic"]


# -- The thin shell wormhole -------------------------------------------------------------

TSW_SPHERE = " + r^2\\left(d\\theta^2 + \\sin^2\\theta\\,d\\phi^2\\right)"
TSW_KINK = "\\left(a + |\\ell|\\right)"
TSW_KINK_BARE = "a + |\\ell|"
TSW_ANGLES = ["\\theta \\in [0, \\pi]", "\\phi \\in [0, 2\\pi)"]


def thin_shell_wormhole(system):
    """Visser's thin shell wormhole, two copies of Schwarzschild's exterior r >= a > r_s
    joined at r = a, in two charts. The first is Schwarzschild's own on either side, where
    the spacetime is vacuum and every value is Schwarzschild's. The second runs through the
    throat with r = a + |ell|, ell = r - a on one side and a - r on the other, so the metric is
    continuous at ell = 0 with a jump in its first derivative: the Christoffel symbols carry
    sgn ell, and the curvature a delta at the throat, whose Einstein tensor is the surface
    energy and pressure of the Israel junction conditions. The Kretschmann scalar squares
    that delta, so it is stated off the throat, where it is Schwarzschild's.
    thin_shell_wormhole.md derives each."""
    parameters = ["r_s", "a"]
    if system == "spherical":
        coords = ["t", "r", "\\theta", "\\phi"]
        f = "\\left(1 - \\dfrac{r_s}{r}\\right)"
        bare = "1 - \\dfrac{r_s}{r}"
        line = "ds^2 = -" + f + "c^2dt^2 + \\dfrac{dr^2}{" + bare + "}" + TSW_SPHERE
        probe = vm.Reader(coords, parameters, ())
        r, rs = probe.symbol["r"], probe.parameters["r_s"]
        return {
            "metric_id": "thin_shell_wormhole",
            "system": {"id": "spherical", "name": "Schwarzschild", "coords": coords,
                       "domains": ["t \\in (-\\infty, \\infty)", "r \\in [a, \\infty)"] + TSW_ANGLES
                       + ["r = a \\;\\text{(throat)}"],
                       "parameters": parameters, "line_element": line},
            "chart_line_element": "ds^2 = -" + f + "dt^2 + \\dfrac{dr^2}{" + bare + "}" + TSW_SPHERE,
            "printer": {"rising": [rs], "lead": [r, rs], "flip": False},
            "components": {"metric_components": {("t", "t"): "-" + f, ("r", "r"): f + "^{-1}"},
                           "inverse_metric_components": {("t", "t"): "-" + f + "^{-1}", ("r", "r"): bare}},
            "kretschmann": "\\dfrac{12r_s^2}{r^6}",
        }
    coords = ["t", "\\ell", "\\theta", "\\phi"]
    f = "\\left(1 - \\dfrac{r_s}{" + TSW_KINK_BARE + "}\\right)"
    bare = "1 - \\dfrac{r_s}{" + TSW_KINK_BARE + "}"
    sphere = " + " + TSW_KINK + "^2\\left(d\\theta^2 + \\sin^2\\theta\\,d\\phi^2\\right)"
    line = "ds^2 = -" + f + "c^2dt^2 + \\dfrac{d\\ell^2}{" + bare + "}" + sphere
    probe = vm.Reader(coords, parameters, ())
    ell, rs, a = probe.symbol["\\ell"], probe.parameters["r_s"], probe.parameters["a"]
    # chart_printer.kink prints each value as A(|ell|) + sgn(ell) B(|ell|); the smooth part and the
    # coefficient of the delta at the throat are each factored on their own, the delta's term last.
    delta = sp.Symbol("_delta" + ell.name)

    def then(e):
        return sp.factor(e.subs(delta, 0)) + sp.factor(sp.diff(e, delta)) * delta

    pretty, overrides = cp.kink(ell, then)
    absolute = next(p for p in overrides if p.name.startswith("_abs"))
    return {
        "metric_id": "thin_shell_wormhole",
        "system": {"id": "throat", "name": "Through the Throat", "coords": coords,
                   "domains": ["t \\in (-\\infty, \\infty)", "\\ell \\in (-\\infty, \\infty)"] + TSW_ANGLES
                   + ["\\ell = 0 \\;\\text{(throat)}"],
                   "parameters": parameters, "line_element": line},
        "chart_line_element": "ds^2 = -" + f + "dt^2 + \\dfrac{d\\ell^2}{" + bare + "}" + sphere,
        "printer": {"rising": [rs], "lead": [a, absolute, rs], "flip": False, "last": [delta], "overrides": overrides},
        "pretty": pretty,
        "bracketed": pretty,
        "components": {"metric_components": {("t", "t"): "-" + f, ("\\ell", "\\ell"): f + "^{-1}"},
                       "inverse_metric_components": {("t", "t"): "-" + f + "^{-1}", ("\\ell", "\\ell"): bare}},
        "kretschmann_where": "\\ell",
        "check": thin_shell_pullback,
    }


def thin_shell_pullback(chart):
    """The chart through the throat with ell > 0 is Schwarzschild's with r = a + ell."""
    ell = chart.reader.symbol["\\ell"]
    spherical = vm.Reader(["t", "r", "\\theta", "\\phi"], ["r_s", "a"], ())
    r = spherical.symbol["r"]
    g = vm.metric_from_line_element(
        spherical, "ds^2 = -\\left(1 - \\dfrac{r_s}{r}\\right)dt^2 + \\dfrac{dr^2}{1 - \\dfrac{r_s}{r}}" + TSW_SPHERE,
        ["t", "r", "\\theta", "\\phi"])
    a = chart.reader.parameters["a"]
    names = {spherical.parameters[n]: chart.reader.parameters[n] for n in ("r_s", "a")}
    names[spherical.symbol["t"]] = chart.reader.symbol["t"]
    names[spherical.symbol["\\theta"]] = chart.reader.symbol["\\theta"]
    for side in (1, -1):
        pulled = g.subs(names).subs(r, a + side * ell)
        mine = chart.geo.g.subs(sp.sign(ell), side)
        if vm.norm(pulled - mine) != sp.zeros(4, 4):
            raise AssertionError(f"thin_shell_wormhole: the chart through the throat misses Schwarzschild's "
                                 f"on the side sgn ell = {side}")
# -- Kantowski-Sachs -------------------------------------------------------------------

def kantowski_sachs(system):
    """The three charts of the Kantowski-Sachs cosmologies, whose moments are cylinders of
    spheres: the comoving chart with both scale factors free, the vacuum member, which is the
    inside of Schwarzschild's horizon with the areal radius T as its time, and the dust
    solution in the parametric time eta, c dt = 2 b_0 cos^2 eta d eta. The last two are checked
    to be the first with a and b as their conventions state them."""
    sphere = "\\left(d\\theta^2 + \\sin^2\\theta\\,d\\phi^2\\right)"
    angles = ["\\theta \\in [0, \\pi]", "\\phi \\in [0, 2\\pi)"]
    if system == "comoving":
        coords, parameters = ["t", "r", "\\theta", "\\phi"], ["a = a(t)", "b = b(t)"]
        probe = vm.Reader(coords, parameters, ())
        a, b = probe.parameters["a"], probe.parameters["b"]
        return {
            "metric_id": "kantowski_sachs",
            "system": {"id": system, "name": "Comoving", "coords": coords,
                       "domains": ["t \\in (-\\infty, \\infty)", "r \\in (-\\infty, \\infty)"] + angles,
                       "parameters": parameters,
                       "line_element": "ds^2 = -c^2dt^2 + a^2dr^2 + b^2" + sphere},
            "chart_line_element": "ds^2 = -dt^2 + a^2dr^2 + b^2" + sphere,
            "printer": {"primed": ["a", "b"], "lead": [a, b]},
            # The sphere's own curvature, which the printer leaves as a bracketed difference.
            "rewrite": [("-\\left(-\\left(b'\\right)^2\\,\\sin^2\\theta - \\sin^2\\theta\\right)",
                         "\\left(\\left(b'\\right)^2 + 1\\right)\\sin^2\\theta"),
                        ("\\left(-\\left(b'\\right)^2\\,\\sin^2\\theta - \\sin^2\\theta\\right)",
                         "-\\left(\\left(b'\\right)^2 + 1\\right)\\sin^2\\theta")],
            # The scalars term by term over the four frame curvatures a''/a, b''/b, a'b'/(ab)
            # and (b'^2 + 1)/b^2, the last the sphere's own curvature plus its expansion.
            "ricci_scalar": ("2\\left(\\dfrac{a''}{a} + \\dfrac{2b''}{b} + \\dfrac{2a'\\,b'}{a\\,b}"
                             " + \\dfrac{\\left(b'\\right)^2 + 1}{b^2}\\right)"),
            "kretschmann": ("4\\left(\\dfrac{a''}{a}\\right)^2 + 8\\left(\\dfrac{b''}{b}\\right)^2"
                            " + 8\\left(\\dfrac{a'\\,b'}{a\\,b}\\right)^2"
                            " + 4\\left(\\dfrac{\\left(b'\\right)^2 + 1}{b^2}\\right)^2"),
        }
    if system == "schwarzschild_interior":
        coords, parameters = ["T", "r", "\\theta", "\\phi"], ["r_s"]
        f = "\\left(\\dfrac{r_s}{T} - 1\\right)"
        line = "ds^2 = -\\dfrac{dT^2}{\\dfrac{r_s}{T} - 1} + " + f + "dr^2 + T^2" + sphere
        probe = vm.Reader(coords, parameters, ())
        return {
            "metric_id": "kantowski_sachs",
            "system": {"id": system, "name": "Vacuum (Schwarzschild Interior)", "coords": coords,
                       "domains": ["T \\in (0, r_s)", "r \\in (-\\infty, \\infty)"] + angles,
                       "parameters": parameters, "line_element": line},
            "chart_line_element": line,
            "printer": {"lead": [probe.parameters["r_s"], probe.symbol["T"]]},
            "components": {"metric_components": {("T", "T"): "-" + f + "^{-1}", ("r", "r"): "\\dfrac{r_s}{T} - 1"},
                           "inverse_metric_components": {("T", "T"): "-" + f, ("r", "r"): f + "^{-1}"}},
            "check": lambda c: kantowski_sachs_member(c, lambda T, p: (
                sp.sqrt(p["r_s"] / T - 1), T, 1 / sp.sqrt(p["r_s"] / T - 1))),
        }
    coords, parameters = ["\\eta", "r", "\\theta", "\\phi"], ["b_0", "\\kappa"]
    a = "\\left(1 + \\left(\\eta + \\kappa\\right)\\tan\\eta\\right)"
    line = "ds^2 = b_0^2\\cos^4\\eta\\left(-4\\,d\\eta^2 + d\\theta^2 + \\sin^2\\theta\\,d\\phi^2\\right) + " + a + "^2dr^2"
    probe = vm.Reader(coords, parameters, ())
    eta = probe.symbol["\\eta"]
    return {
        "metric_id": "kantowski_sachs",
        "system": {"id": system, "name": "Dust, Parametric Time", "coords": coords,
                   "domains": ["\\eta \\in (-\\pi/2, \\pi/2)", "r \\in (-\\infty, \\infty)"] + angles,
                   "parameters": parameters, "line_element": line},
        "chart_line_element": line,
        "printer": {"lead": [probe.parameters["b_0"], KS_Q, sp.tan(eta)],
                    "named": {KS_A: "1 + \\left(\\eta + \\kappa\\right)\\tan\\eta", KS_Q: "\\eta + \\kappa"}},
        "pretty": kantowski_sachs_dust(eta, probe.parameters["kappa"]),
        "check": lambda c: kantowski_sachs_member(c, lambda eta, p: (
            1 + (eta + p["kappa"]) * sp.tan(eta), p["b_0"] * sp.cos(eta) ** 2, 2 * p["b_0"] * sp.cos(eta) ** 2)),
    }


KS_A, KS_Q = sp.symbols("KSA KSQ")


def kantowski_sachs_dust(eta, kappa):
    """A `pretty` for the dust chart, every value of which is a rational function of
    t = tan eta and Q = eta + kappa, since cos^4 eta = 1/(1 + t^2)^2, dt/d eta = 1 + t^2 and
    a = 1 + Q t. Each value is factored there, which is unique, and printed with 1 + t^2 as
    1/cos^2 eta, the factor 1 + Q t as the scale factor the line element writes, and Q as
    eta + kappa."""
    t, C = sp.symbols("_t _C")

    def pretty(value):
        x = sp.sympify(value).subs({sp.tan(eta): t, sp.sin(eta): t * C, sp.cos(eta): C})
        num, den = (sp.expand(p) for p in sp.fraction(sp.together(x.subs(kappa, KS_Q - eta))))
        if num.has(eta) or den.has(eta):
            raise AssertionError(f"kantowski_sachs: eta stands outside eta + kappa in {value}")

        def even(p):
            # C stands only in even powers, and C^2 = 1/(1 + t^2).
            poly = sp.Poly(sp.expand(p), C)
            if any(k % 2 for (k,) in poly.monoms()):
                raise AssertionError(f"kantowski_sachs: an odd power of cos eta in {value}")
            return sum(c * (1 + t ** 2) ** (-k // 2) for (k,), c in poly.terms())

        out = sp.factor(sp.cancel(even(num) / even(den)))
        result = sp.Integer(1)
        for f in sp.Mul.make_args(out):
            base, k = (f.base, f.exp) if f.is_Pow else (f, sp.Integer(1))
            # sympy leaves -(t^2 + 1) alone as a sum, so each factor is matched up to its sign.
            named = [(s, to) for poly, to in ((t ** 2 + 1, sp.cos(eta) ** -2), (KS_Q * t + 1, KS_A))
                     for s in (1, -1) if sp.expand(base - s * poly) == 0]
            if named:
                (s, to), = named
                result *= s ** k * to ** k
            else:
                result *= base.subs(t, sp.tan(eta)) ** k
        return result

    return pretty


def kantowski_sachs_member(chart, member):
    """The chart's metric is -(lapse d time)^2 + a^2 dr^2 + b^2 dOmega^2 with (a, b, lapse) as
    `member` states them, so it is the comoving chart with c dt = lapse d time."""
    time, r, theta = chart.symbols[0], chart.symbols[1], chart.symbols[2]
    a, b, lapse = member(time, chart.reader.parameters)
    expected = sp.diag(-lapse ** 2, a ** 2, b ** 2, b ** 2 * sp.sin(theta) ** 2)
    for i, j in itertools.product(range(4), repeat=2):
        if vm.norm(chart.geo.g[i, j] - expected[i, j]) != 0:
            raise AssertionError(f"kantowski_sachs: slot {i}{j} is {chart.geo.g[i, j]}, not {expected[i, j]}")


KS_CHARTS = ["comoving", "schwarzschild_interior", "dust"]


# -- Robinson-Trautman -------------------------------------------------------------------

def robinson_trautman(system_id):
    """Robinson and Trautman's metric, -2H c^2du^2 - 2c du dr + r^2 times a metric on the wave
    fronts, in the two charts its literature uses: the chart of their 1960 letter, with the
    fronts' metric (dx^2 + dy^2)/P^2 and zeta = (x + iy)/sqrt 2 the complex coordinate of
    later work, and the axisymmetric chart of the numerical work, with the fronts' metric
    (dtheta^2 + sin^2 theta dphi^2)/f^2. H and P, or f, are left free in every tensor, so
    no component assumes a field equation; robinson_trautman_check confirms before anything
    is written that the vacuum H leaves the Robinson-Trautman equation as the one field
    equation, and robinson_trautman.md records the choices."""
    if system_id == "stereographic":
        coords, name = ["u", "r", "x", "y"], "Robinson-Trautman"
        parameters = ["P = P(u,x,y)", "H = H(u,r,x,y)"]
        fronts = "\\dfrac{r^2}{P^2}\\left(dx^2 + dy^2\\right)"
        domains = ["x \\in (-\\infty, \\infty)", "y \\in (-\\infty, \\infty)"]
    else:
        coords, name = ["u", "r", "\\theta", "\\phi"], "Axisymmetric"
        parameters = ["f = f(u,\\theta)", "H = H(u,r,\\theta)"]
        fronts = "\\dfrac{r^2}{f^2}\\left(d\\theta^2 + \\sin^2\\theta\\,d\\phi^2\\right)"
        domains = ["\\theta \\in [0, \\pi]", "\\phi \\in [0, 2\\pi)"]
    line = "ds^2 = -2H{2}du^2 - 2{1}du\\,dr + " + fronts
    probe = vm.Reader(coords, parameters, ())
    shape = probe.parameters["P" if system_id == "stereographic" else "f"]
    return {
        "metric_id": "robinson_trautman",
        "system": {"id": system_id, "name": name, "coords": coords,
                   "domains": ["u \\in (-\\infty, \\infty)", "r \\in (0, \\infty)"] + domains,
                   "parameters": parameters,
                   "line_element": line.replace("{2}", "\\,c^2").replace("{1}", "c\\,")},
        "chart_line_element": line.replace("{2}", "\\,").replace("{1}", ""),
        "printer": {"lead": [shape, probe.parameters["H"]], "flip": False,
                    "collect": lambda poly, printer: cp.collect_by(poly, [shape], printer)},
        "check": robinson_trautman_check,
    }


def robinson_trautman_vacuum(chart):
    """The chart's front function, the Gaussian curvature K of the fronts at r = 1, the
    Laplacian of the fronts, and the vacuum 2H = K - 2r d_u ln P - 2m/r, with m a length."""
    u, r, a, b = chart.symbols
    m = sp.Symbol("m", positive=True)
    if chart.coords_tex[2] == "x":
        P = chart.reader.parameters["P"]

        def laplacian(F):
            return P ** 2 * (sp.diff(F, a, 2) + sp.diff(F, b, 2))
        K = laplacian(sp.log(P))
    else:
        P = chart.reader.parameters["f"]

        def laplacian(F):
            return P ** 2 * sp.diff(sp.sin(a) * sp.diff(F, a), a) / sp.sin(a)
        K = P ** 2 + laplacian(sp.log(P))
    H = (K - 2 * r * sp.diff(P, u) / P - 2 * m / r) / 2
    return P, K, laplacian, H, m


def robinson_trautman_check(chart):
    """With the vacuum H every component of the Ricci tensor vanishes but R_uu, which is
    (Laplacian K + 12m d_u ln P)/(2r^2), the Robinson-Trautman equation; with P = 1 + (x^2 + y^2)/4
    or f = 1 the fronts have K = 1 and the Kretschmann scalar is Schwarzschild's 48m^2/r^6; and the
    axisymmetric chart is the chart of x and y pulled back through x + iy = 2 tan(theta/2) e^{i phi},
    P = f (1 + (x^2 + y^2)/4)."""
    u, r, a, b = chart.symbols
    P, K, laplacian, H, m = robinson_trautman_vacuum(chart)
    free = chart.reader.parameters["H"]
    ricci = chart.geo.ricci_ll()
    equation = (laplacian(K) + 12 * m * sp.diff(P, u) / P) / (2 * r ** 2)
    for i in range(4):
        for j in range(i, 4):
            value = sp.sympify(ricci[i][j]).subs(free, H).doit()
            if vm.norm(value - (equation if i == j == 0 else 0)) != 0:
                raise AssertionError(f"robinson_trautman: R_{chart.coords_tex[i]}{chart.coords_tex[j]} "
                                     "with the vacuum H is not the Robinson-Trautman equation")
    round_front = 1 + (a ** 2 + b ** 2) / 4 if chart.coords_tex[2] == "x" else sp.Integer(1)
    schwarzschild = {free: (1 - 2 * m / r) / 2, P: round_front}
    if vm.norm(K.subs(P, round_front).doit() - 1) != 0:
        raise AssertionError("robinson_trautman: the round front has not K = 1")
    if vm.norm(sp.sympify(chart.geo.kretschmann()).subs(schwarzschild).doit() - 48 * m ** 2 / r ** 6) != 0:
        raise AssertionError("robinson_trautman: the round front is not Schwarzschild")
    if chart.coords_tex[2] == "x":
        return
    # The pullback, with t = tan(theta/2) so that every entry is rational: x = 2t cos phi and
    # y = 2t sin phi take (dx^2 + dy^2)/(1 + (x^2 + y^2)/4)^2 to 4dt^2/(1 + t^2)^2 + 4t^2dphi^2/(1 + t^2)^2,
    # which is dtheta^2 + sin^2 theta dphi^2, since dtheta = 2dt/(1 + t^2) and sin theta = 2t/(1 + t^2).
    t = sp.Symbol("t", positive=True)
    image = [2 * t * sp.cos(b), 2 * t * sp.sin(b)]
    J = sp.Matrix(2, 2, lambda i, j: sp.diff(image[i], (t, b)[j]))
    flat = (J.T * J) / (1 + (image[0] ** 2 + image[1] ** 2) / 4) ** 2
    sphere = sp.Matrix([[4 / (1 + t ** 2) ** 2, 0], [0, 4 * t ** 2 / (1 + t ** 2) ** 2]])
    for i in range(2):
        for j in range(2):
            if vm.norm(flat[i, j] - sphere[i, j]) != 0:
                raise AssertionError("robinson_trautman: the stereographic fronts pulled back are not the sphere's")


RT_CHARTS = ["stereographic", "axisymmetric"]


# -- Dilaton black hole ----------------------------------------------------------------

DILATON_ANGLES = ["\\theta \\in [0, \\pi]", "\\phi \\in [0, 2\\pi)"]
DILATON_MARKS = ["r = r_s \\;\\text{(the horizon)}", "r = r_d \\;\\text{(the singularity)}"]
DILATON_F = "\\left(1 - \\dfrac{r_s}{r}\\right)"
DILATON_H = "\\left(1 - \\dfrac{r_d}{r}\\right)"
DILATON_SPHERE = "\\left(d\\theta^2 + \\sin^2\\theta\\,d\\phi^2\\right)"


def dilaton_black_hole(system_id):
    """The charged black hole of Gibbons and Maeda and of Garfinkle, Horowitz and Strominger.
    Its Einstein metric is Schwarzschild's on the plane of t and r, f = 1 - r_s/r, with spheres
    of area 4 pi r(r - r_d), r_d = Q^2/M in units G = c = 1; the static chart and the two
    Eddington-Finkelstein charts on Schwarzschild's tortoise coordinate print it. The two
    string charts print the string metric e^{2 phi} g of the magnetically charged hole,
    e^{-2 phi} = 1 - r_d/r, and of the electrically charged one, e^{2 phi} = 1 - r_d/r, each in
    the same r. dilaton_field_equations checks the Einstein charts against Einstein's equations
    with the Maxwell field and the dilaton, Maxwell's equations and the dilaton's equation, and
    dilaton_string_frame checks each string chart to be that multiple of the static chart."""
    f, h, sphere = DILATON_F, DILATON_H, DILATON_SPHERE
    bare = "1 - \\dfrac{r_s}{r}"
    parameters = ["r_s", "r_d"]
    radial = ["r \\in (r_d, \\infty)"] + DILATON_ANGLES + DILATON_MARKS
    coords = ["t", "r", "\\theta", "\\phi"]
    time = ["t \\in (-\\infty, \\infty)"]
    extra = {}
    if system_id == "static":
        name = "Static Spherical"

        def line(c2):
            return ("ds^2 = -" + f + c2 + "dt^2 + \\dfrac{dr^2}{" + bare + "} + r\\left(r - r_d\\right)" + sphere)
        spec_line, chart_line = line("c^2"), line("")
        extra = {"components": {"metric_components": {("t", "t"): "-" + f, ("r", "r"): f + "^{-1}"},
                                "inverse_metric_components": {("t", "t"): "-" + f + "^{-1}", ("r", "r"): bare}},
                 "check": dilaton_field_equations}
    elif system_id.startswith("eddington_finkelstein"):
        null, sign = ("u", "-") if system_id.endswith("outgoing") else ("v", "+")
        coords = [null, "r", "\\theta", "\\phi"]
        time = [null + " \\in (-\\infty, \\infty)"]
        name = ("Outgoing" if null == "u" else "Ingoing") + " Eddington-Finkelstein"
        spec_line = chart_line = ("ds^2 = -" + f + "d" + null + "^2 " + sign + " 2\\,d" + null + "\\,dr"
                                  + " + r\\left(r - r_d\\right)" + sphere)
        one = "-1" if null == "u" else "1"
        extra = {"components": {"metric_components": {(null, null): "-" + f, (null, "r"): one, ("r", null): one},
                                "inverse_metric_components": {(null, "r"): one, ("r", null): one, ("r", "r"): bare}},
                 "check": dilaton_field_equations}
    elif system_id == "string_magnetic":
        name = "String Metric, Magnetic Charge"

        def line(c2):
            return ("ds^2 = -\\dfrac{1 - \\dfrac{r_s}{r}}{1 - \\dfrac{r_d}{r}}" + c2 + "dt^2 + \\dfrac{dr^2}{"
                    + f + h + "} + r^2" + sphere)
        spec_line, chart_line = line("c^2"), line("")
        extra = {"check": lambda chart: dilaton_string_frame(chart, -1)}
    else:
        name = "String Metric, Electric Charge"

        def line(c2):
            return ("ds^2 = " + h + "\\left(-" + f + c2 + "dt^2 + \\dfrac{dr^2}{" + bare + "}\\right)"
                    " + \\left(r - r_d\\right)^2" + sphere)
        spec_line, chart_line = line("c^2"), line("")
        extra = {"check": lambda chart: dilaton_string_frame(chart, 1)}
    probe = vm.Reader(coords, parameters, ())
    r, rs, rd = probe.symbol["r"], probe.parameters["r_s"], probe.parameters["r_d"]
    return {
        "metric_id": "dilaton_black_hole",
        "system": {"id": system_id, "name": name, "coords": coords, "domains": time + radial,
                   "parameters": parameters, "line_element": spec_line},
        "chart_line_element": chart_line,
        "printer": {"rising": [rs, rd], "lead": [r, rs, rd], "flip": False},
        **extra,
    }


def dilaton_field_equations(chart):
    """The equations of the action R - 2(d phi)^2 - e^{-2 phi} F^2 for the magnetically charged
    hole, e^{-2 phi} = 1 - r_d/r and F = Q sin(theta) d theta ^ d phi with Q^2 = r_s r_d/2, slot
    by slot: G_mu nu = 2 d_mu phi d_nu phi - g_mu nu (d phi)^2 + 2 e^{-2 phi}(F_mu a F_nu^a
    - g_mu nu F^2/4), d_mu(sqrt(-g) e^{-2 phi} F^mu nu) = 0, and box phi = -e^{-2 phi} F^2/2."""
    x = chart.symbols
    r, th = chart.reader.symbol["r"], chart.reader.symbol["\\theta"]
    rs, rd = chart.reader.parameters["r_s"], chart.reader.parameters["r_d"]
    g, ginv = chart.geo.g, chart.geo.ginv
    e = 1 - rd / r                                   # e^{-2 phi}
    dphi = [sp.diff(-sp.log(e) / 2, c) for c in x]
    F = sp.zeros(4, 4)
    F[2, 3], F[3, 2] = sp.sin(th), -sp.sin(th)       # F / Q
    Q2 = rs * rd / 2
    Fup = ginv * F * ginv
    F2 = Q2 * sum(F[a, b] * Fup[a, b] for a in range(4) for b in range(4))
    dphi2 = sum(ginv[a, b] * dphi[a] * dphi[b] for a in range(4) for b in range(4))
    G = chart.geo.einstein_ll()
    for a in range(4):
        for b in range(a, 4):
            FF = Q2 * sum(F[a, c] * F[b, d] * ginv[c, d] for c in range(4) for d in range(4))
            T = 2 * dphi[a] * dphi[b] - g[a, b] * dphi2 + 2 * e * (FF - g[a, b] * F2 / 4)
            if vm.norm(vm._at(G, (a, b)) - T) != 0:
                raise AssertionError(f"dilaton_black_hole: the Einstein tensor misses the stress of the Maxwell "
                                     f"field and the dilaton in slot {chart.coords_tex[a]}{chart.coords_tex[b]}")
    # sqrt(-g) is r(r - r_d) sin(theta) in all three charts.
    root = r * (r - rd) * sp.sin(th)
    if sp.simplify(root ** 2 + g.det()) != 0:
        raise AssertionError("dilaton_black_hole: sqrt(-g) is not r(r - r_d) sin(theta)")
    for b in range(4):
        if sp.simplify(sum(sp.diff(root * e * Fup[a, b], x[a]) for a in range(4))) != 0:
            raise AssertionError(f"dilaton_black_hole: Maxwell's equations fail along {chart.coords_tex[b]}")
    box = sum(sp.diff(root * ginv[a, b] * dphi[b], x[a]) for a in range(4) for b in range(4)) / root
    if sp.simplify(box + e * F2 / 2) != 0:
        raise AssertionError("dilaton_black_hole: the dilaton's equation fails")


def dilaton_string_frame(chart, sign):
    """The string chart is e^{2 phi} times the static chart, with e^{2 phi} = (1 - r_d/r)^sign:
    sign = 1 for the electrically charged hole and -1 for the magnetically charged one."""
    spec = dilaton_black_hole("static")
    source = cp.Chart(spec["system"]["coords"], spec["system"]["parameters"], spec["chart_line_element"])
    names = dict(zip(source.symbols, chart.symbols))
    names.update({source.reader.parameters[n]: chart.reader.parameters[n] for n in ("r_s", "r_d")})
    r, rd = chart.reader.symbol["r"], chart.reader.parameters["r_d"]
    scaled = (1 - rd / r) ** sign * source.geo.g.subs(names)
    if vm.norm(scaled - chart.geo.g) != sp.zeros(4, 4):
        raise AssertionError(f"dilaton_black_hole: the string chart of sign {sign} is not e^(2 phi) times the "
                             "Einstein metric")


DILATON_CHARTS = ["static", "eddington_finkelstein_outgoing", "eddington_finkelstein_ingoing",
                  "string_magnetic", "string_electric"]


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
          "global_monopole": [lambda s=s: global_monopole(s) for s in GM_CHARTS],
          "domain_wall": [lambda s=s: domain_wall(s) for s in DW_CHARTS],
          "majumdar_papapetrou": [lambda s=s: majumdar_papapetrou(s) for s in MP_CHARTS],
          "thin_shell_wormhole": [lambda s=s: thin_shell_wormhole(s) for s in ("spherical", "throat")],
          "kantowski_sachs": [lambda s=s: kantowski_sachs(s) for s in KS_CHARTS],
          "robinson_trautman": [lambda s=s: robinson_trautman(s) for s in RT_CHARTS],
          "mcvittie": [lambda s=s: mcvittie(s) for s in ("isotropic", "areal")],
          "dilaton_black_hole": [lambda s=s: dilaton_black_hole(s) for s in DILATON_CHARTS]}


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



# -- Melvin ------------------------------------------------------------------------------

def melvin(system):
    """Melvin's magnetic universe in its cylindrical chart, with Lambda = 1 + B^2 rho^2/4, and
    Ernst's Schwarzschild black hole inside it in Schwarzschild's coordinates, with
    Lambda = 1 + B^2 r^2 sin^2(theta)/4. Every value is printed around 4 + B^2 rho^2, or
    4 + B^2 r^2 sin^2(theta), which is 4 Lambda, and 4 - B^2 rho^2, which vanishes on the widest circle,
    and the metric and its inverse as the line element writes them. melvin_maxwell checks both
    against the Einstein-Maxwell equations, and melvin.md beside this file is the derivation."""
    if system == "cylindrical":
        coords, parameters = ["t", "\\rho", "\\phi", "z"], ["B"]
        L = "\\left(1 + \\dfrac{B^2\\rho^2}{4}\\right)"

        def line(c2):
            return f"ds^2 = {L}^2\\left(-{c2}dt^2 + d\\rho^2 + dz^2\\right) + \\dfrac{{\\rho^2}}{{{L}^2}}d\\phi^2"
        probe = vm.Reader(coords, parameters, ())
        rho, B = probe.symbol["\\rho"], probe.parameters["B"]
        named = [(sp.Symbol("MP"), 4 + B ** 2 * rho ** 2, "4 + B^2\\rho^2"),
                 (sp.Symbol("MM"), 4 - B ** 2 * rho ** 2, "4 - B^2\\rho^2")]
        merges = [(2 + B * rho, 2 - B * rho, named[1][0])]
        metric = {("t", "t"): "-" + L + "^2", ("\\rho", "\\rho"): L + "^2", ("z", "z"): L + "^2",
                  ("\\phi", "\\phi"): "\\dfrac{\\rho^2}{" + L + "^2}"}
        inverse = {("t", "t"): "-" + L + "^{-2}", ("\\rho", "\\rho"): L + "^{-2}", ("z", "z"): L + "^{-2}",
                   ("\\phi", "\\phi"): "\\dfrac{" + L + "^2}{\\rho^2}"}
        domains = ["t \\in (-\\infty, \\infty)", "\\rho \\in [0, \\infty)", "\\phi \\in [0, 2\\pi)",
                   "z \\in (-\\infty, \\infty)"]
        lead = [B, rho]
        name = "Cylindrical"
    else:
        coords, parameters = ["t", "r", "\\theta", "\\phi"], ["r_s", "B"]
        L = "\\left(1 + \\dfrac{B^2r^2\\sin^2\\theta}{4}\\right)"
        f = "\\left(1 - \\dfrac{r_s}{r}\\right)"

        def line(c2):
            return (f"ds^2 = {L}^2\\left(-{f}{c2}dt^2 + \\dfrac{{dr^2}}{{1 - \\dfrac{{r_s}}{{r}}}} + r^2d\\theta^2\\right)"
                    f" + \\dfrac{{r^2\\sin^2\\theta}}{{{L}^2}}d\\phi^2")
        probe = vm.Reader(coords, parameters, ())
        r, th = probe.symbol["r"], probe.symbol["\\theta"]
        B, rs = probe.parameters["B"], probe.parameters["r_s"]
        s = sp.sin(th)
        named = [(sp.Symbol("MP"), 4 + B ** 2 * r ** 2 * s ** 2, "4 + B^2r^2\\sin^2\\theta"),
                 (sp.Symbol("MM"), 4 - B ** 2 * r ** 2 * s ** 2, "4 - B^2r^2\\sin^2\\theta")]
        merges = [(2 + B * r * s, 2 - B * r * s, named[1][0])]
        metric = {("t", "t"): "-" + L + "^2" + f, ("r", "r"): "\\dfrac{" + L + "^2}{1 - \\dfrac{r_s}{r}}",
                  ("\\theta", "\\theta"): L + "^2r^2", ("\\phi", "\\phi"): "\\dfrac{r^2\\sin^2\\theta}{" + L + "^2}"}
        inverse = {("t", "t"): "-\\dfrac{1}{" + L + "^2" + f + "}", ("r", "r"): "\\dfrac{1 - \\dfrac{r_s}{r}}{" + L + "^2}",
                   ("\\theta", "\\theta"): "\\dfrac{1}{" + L + "^2r^2}",
                   ("\\phi", "\\phi"): "\\dfrac{" + L + "^2}{r^2\\sin^2\\theta}"}
        domains = ["t \\in (-\\infty, \\infty)", "r \\in (0, \\infty)", "\\theta \\in [0, \\pi]",
                   "\\phi \\in [0, 2\\pi)", "r = r_s \\;\\text{(horizon)}"]
        lead = [B, r, rs, sp.cos(th), s]
        name = "Ernst"
    return {
        "metric_id": "melvin",
        "system": {"id": system, "name": name, "coords": coords, "domains": domains, "parameters": parameters,
                   "line_element": line("c^2")},
        "chart_line_element": line(""),
        "printer": {"lead": lead, "named": {p: text for p, _, text in named}},
        "pretty": named_factors(named, merges),
        "components": {"metric_components": metric, "inverse_metric_components": inverse},
        **({"geodesics": MELVIN_GEODESICS} if system == "cylindrical" else {}),
        "check": melvin_maxwell,
    }


def melvin_maxwell(chart):
    """G_mu nu = 2(F_mu a F_nu^a - g_mu nu F^2/4) with F = dA and A_phi = 2B s^2/(4 + B^2 s^2), s the
    distance from the axis, rho or r sin(theta), in every slot, and d(sqrt(-g) F^mu nu) = 0: Einstein's
    equations with Maxwell's field in units where B is an inverse length, and Maxwell's equations."""
    x = chart.symbols
    B = chart.reader.parameters["B"]
    s = x[1] if chart.coords_tex[1] == "\\rho" else x[1] * sp.sin(x[2])
    A = [0, 0, 0, 0]
    A[chart.coords_tex.index("\\phi")] = 2 * B * s ** 2 / (4 + B ** 2 * s ** 2)
    g, ginv = chart.geo.g, chart.geo.ginv
    F = sp.Matrix(4, 4, lambda a, b: sp.diff(A[b], x[a]) - sp.diff(A[a], x[b]))
    Fup = ginv * F * ginv
    F2 = sum(F[a, b] * Fup[a, b] for a in range(4) for b in range(4))
    G = chart.geo.einstein_ll()
    for a in range(4):
        for b in range(a, 4):
            T = 2 * (sum(F[a, c] * F[b, d] * ginv[c, d] for c in range(4) for d in range(4)) - g[a, b] * F2 / 4)
            if vm.norm(vm._at(G, (a, b)) - T) != 0:
                raise AssertionError(f"melvin: the Einstein tensor misses Maxwell's stress in slot "
                                     f"{chart.coords_tex[a]}{chart.coords_tex[b]}")
    # sqrt(-g) is rho Lambda^2 or r^2 sin(theta) Lambda^2, each factor positive on the chart.
    root = sp.powdenest(sp.sqrt(sp.factor(-g.det())), force=True)
    for b in range(4):
        if sp.simplify(sum(sp.diff(root * Fup[a, b], x[a]) for a in range(4))) != 0:
            raise AssertionError(f"melvin: Maxwell's equations fail along {chart.coords_tex[b]}")


# The cylindrical geodesics with the three terms that share 2B^2 rho/(4 + B^2 rho^2) gathered.
MELVIN_GEODESICS = [
    "\\ddot{t} + \\dfrac{4B^2\\,\\rho}{4 + B^2\\rho^2}\\dot{t}\\dot{\\rho} = 0",
    "\\ddot{\\rho} + \\dfrac{2B^2\\,\\rho}{4 + B^2\\rho^2}\\left(\\dot{t}^2 + \\dot{\\rho}^2 - \\dot{z}^2\\right)"
    " - \\dfrac{256\\rho\\left(4 - B^2\\rho^2\\right)}{\\left(4 + B^2\\rho^2\\right)^5}\\dot{\\phi}^2 = 0",
    "\\ddot{\\phi} + \\dfrac{2\\left(4 - B^2\\rho^2\\right)}{\\rho\\left(4 + B^2\\rho^2\\right)}\\dot{\\rho}\\dot{\\phi} = 0",
    "\\ddot{z} + \\dfrac{4B^2\\,\\rho}{4 + B^2\\rho^2}\\dot{\\rho}\\dot{z} = 0",
]


CHARTS["melvin"] = [lambda s=s: melvin(s) for s in ("cylindrical", "ernst")]
CHARTS["schwarzschild_ads"] = [lambda s=s: schwarzschild_ads(s) for s in SADS_CHARTS]


# -- Levi-Civita -------------------------------------------------------------------------

def radial_powers(radius, radius_tex, state):
    """A pretty printer for a chart whose every value is one monomial in a radius raised to an
    exponent that holds a parameter, as Levi-Civita's rho^{4 sigma}. The printer writes rational
    exponents only, so the power is handed to it as a placeholder whose text is the radius with
    its exponent: in the numerator where the whole power of the radius is positive, as
    rho^{2 - 4 sigma}, or with the whole power below the line and the rest above it, as
    2 sigma rho^{8 sigma - 8 sigma^2}/rho, and all of it below the line where every term of the
    exponent is negative, as 1/rho^{16 sigma^2 - 8 sigma + 4}. `state["printer"]` is the chart's
    printer, which the chart's check puts there before anything is printed."""
    table = {}

    def pretty(value):
        # The radius is positive on the chart, which is what lets its powers be gathered.
        positive = sp.Dummy(positive=True)
        value = sp.powsimp(sp.powdenest(sp.together(value.subs(radius, positive)), force=True), force=True)
        value = value.subs(positive, radius)
        exponent, rest = sp.Integer(0), sp.Integer(1)
        for f in sp.Mul.make_args(value):
            base, k = f.as_base_exp()
            if base == radius:
                exponent += k
            else:
                rest *= f
        if rest.has(radius):
            raise AssertionError(f"{value} is not one monomial in {radius}")
        whole, symbolic = sp.expand(exponent).as_coeff_Add()
        out = sp.factor(rest)
        if symbolic == 0:
            return out * radius ** whole
        printer = state["printer"]
        negative = all(c.is_negative for c in sp.Add.make_args(symbolic) for c in [c.as_coeff_Mul()[0]])
        if negative and whole <= 0:
            shown, side = -(symbolic + whole), -1
        elif whole > 0:
            shown, side = symbolic + whole, 1
        else:
            shown, side = symbolic, 1
            out *= radius ** whole
        text = radius_tex + "^{" + printer.positive_first(shown) + "}"
        placeholder = table.setdefault(text, sp.Symbol(f"RADIALPOWER{len(table)}", positive=True))
        printer.overrides[placeholder] = text
        return out * placeholder ** side

    return pretty


LC_DOMAINS = ["t \\in (-\\infty, \\infty)", "{r} \\in (0, \\infty)", "\\phi \\in [0, 2\\pi)", "z \\in (-\\infty, \\infty)",
              "{r} = 0 \\;\\text{(the axis, a curvature singularity unless } {flat}\\text{)}"]
LC_FLAT = {"weyl": "\\sigma = 0 \\text{ or } \\sigma = 1/2", "kasner": "p_0\\,p_2\\,p_3 = 0"}


def levi_civita(system):
    """Levi-Civita's static cylinder in Weyl's coordinates, with the mass parameter sigma and the
    conicity C, and in its Kasner form, whose radial coordinate is the proper distance from the
    axis and whose exponents lie on Kasner's circle, p_0 + p_2 + p_3 = p_0^2 + p_2^2 + p_3^2 = 1.
    Every value is one monomial in the radius, printed by radial_powers, and the metric and its
    inverse as the line element writes them. levi_civita.md beside this file derives the Kasner
    form from Weyl's."""
    state = {}

    def check(chart):
        state["printer"] = chart.printer

    if system == "weyl":
        coords, parameters = ["t", "\\rho", "\\phi", "z"], ["\\sigma", "C"]

        def line(c2):
            return (f"ds^2 = -\\rho^{{4\\sigma}}{c2}dt^2 + \\rho^{{4\\sigma(2\\sigma - 1)}}\\left(d\\rho^2 + dz^2\\right)"
                    " + \\dfrac{\\rho^{2(1 - 2\\sigma)}}{C^2}d\\phi^2")
        probe = vm.Reader(coords, parameters, ())
        radius, sigma, C = probe.symbol["\\rho"], probe.parameters["sigma"], probe.parameters["C"]
        metric = {("t", "t"): "-\\rho^{4\\sigma}", ("\\rho", "\\rho"): "\\rho^{4\\sigma(2\\sigma - 1)}",
                  ("z", "z"): "\\rho^{4\\sigma(2\\sigma - 1)}", ("\\phi", "\\phi"): "\\dfrac{\\rho^{2(1 - 2\\sigma)}}{C^2}"}
        inverse = {("t", "t"): "-\\rho^{-4\\sigma}", ("\\rho", "\\rho"): "\\rho^{-4\\sigma(2\\sigma - 1)}",
                   ("z", "z"): "\\rho^{-4\\sigma(2\\sigma - 1)}", ("\\phi", "\\phi"): "C^2\\rho^{-2(1 - 2\\sigma)}"}
        spec = {
            "system": {"id": "weyl", "name": "Weyl", "coords": coords,
                       "domains": [d.replace("{r}", "\\rho").replace("{flat}", LC_FLAT["weyl"]) for d in LC_DOMAINS], "parameters": parameters,
                       "line_element": line("c^2")},
            "chart_line_element": line(""),
            "printer": {"lead": [sigma, C, radius]},
            "pretty": radial_powers(radius, "\\rho", state),
            "components": {"metric_components": metric, "inverse_metric_components": inverse},
            "kretschmann": "\\dfrac{64\\sigma^2\\left(2\\sigma - 1\\right)^2\\left(4\\sigma^2 - 2\\sigma + 1\\right)}"
                           "{\\rho^{4(4\\sigma^2 - 2\\sigma + 1)}}",
            "check": check,
        }
    else:
        coords, parameters = ["t", "r", "\\phi", "z"], ["p_0", "p_2", "p_3", "\\ell"]

        def line(c2):
            return f"ds^2 = -r^{{2p_0}}{c2}dt^2 + dr^2 + \\ell^2r^{{2p_2}}d\\phi^2 + r^{{2p_3}}dz^2"
        probe = vm.Reader(coords, parameters, ())
        radius = probe.symbol["r"]
        p0, p2, p3, ell = (probe.parameters[name] for name in ("p_0", "p_2", "p_3", "ell"))
        metric = {("t", "t"): "-r^{2p_0}", ("\\phi", "\\phi"): "\\ell^2r^{2p_2}", ("z", "z"): "r^{2p_3}"}
        inverse = {("t", "t"): "-r^{-2p_0}", ("\\phi", "\\phi"): "\\dfrac{r^{-2p_2}}{\\ell^2}", ("z", "z"): "r^{-2p_3}"}
        spec = {
            "system": {"id": "kasner", "name": "Kasner", "coords": coords,
                       "domains": [d.replace("{r}", "r").replace("{flat}", LC_FLAT["kasner"]) for d in LC_DOMAINS], "parameters": parameters,
                       "line_element": line("c^2")},
            "chart_line_element": line(""),
            "printer": {"lead": [p0, p2, p3, ell, radius]},
            "pretty": radial_powers(radius, "r", state),
            "components": {"metric_components": metric, "inverse_metric_components": inverse},
            "check": check,
            "after": levi_civita_on_kasner_circle,
        }
    return {"metric_id": "levi_civita", **spec}


LC_KASNER_KRETSCHMANN = "-\\dfrac{16p_0\\,p_2\\,p_3}{r^4}"


def levi_civita_on_kasner_circle(math, chart):
    """The Kasner form is a vacuum only on Kasner's circle, and the file claims its values only
    there, as Kasner's cosmology does. The Christoffel symbols and the Riemann tensor are printed
    for free exponents, which holds on the circle too; here the Ricci and Einstein tensors are
    checked to vanish on the checker's own parametrisation of the circle and written as
    vanishing, the Weyl tensor is the Riemann tensor, and the Kretschmann scalar is checked there."""
    relations = vm.PARAMETER_RELATIONS[("levi_civita", "kasner")]
    surface = {chart.reader.parameters[name]: sp.sympify(value) for name, value in relations.items()}

    def on_circle(value):
        return vm.norm(sp.sympify(value).subs(surface))

    ricci = chart.geo.ricci_ll()
    n = len(chart.symbols)
    for index in vm._indices(n, 2):
        if on_circle(vm._at(ricci, index)) != 0:
            raise AssertionError(f"levi_civita: the Kasner form's Ricci tensor does not vanish in slot {index}")
    if on_circle(chart.reader(LC_KASNER_KRETSCHMANN) - chart.geo.kretschmann()) != 0:
        raise AssertionError("levi_civita: the Kasner form's Kretschmann scalar is misprinted")
    for field in ("ricci_tensor", "einstein_tensor"):
        for variant in math[field]["variants"].values():
            variant["nonzero"] = []
    math["ricci_scalar"] = "R = 0"
    math["kretschmann"] = "K = " + LC_KASNER_KRETSCHMANN
    notes = {"ulll": "C^\\mu{}_{\\nu\\rho\\sigma}", "llll": "C_{\\mu\\nu\\rho\\sigma}"}
    math["weyl_tensor"] = {"default": math["riemann"]["default"], "variants": {
        v: {"note": notes[v], "nonzero": [dict(entry) for entry in block["nonzero"]]}
        for v, block in math["riemann"]["variants"].items()}}
    return math


CHARTS["levi_civita"] = [lambda s=s: levi_civita(s) for s in ("weyl", "kasner")]


# -- Curzon-Chazy ------------------------------------------------------------------------

def curzon_chazy(system):
    """The Curzon-Chazy particle, the member of Weyl's static axisymmetric class whose potential
    is Newton's for a point mass, psi = -m/R and gamma = -m^2 rho^2/(2R^4), in Weyl's canonical
    chart, where R = sqrt(rho^2 + z^2) is a name the chart defines and every value is printed
    around it, and in the spherical chart rho = r sin(theta), z = r cos(theta), where R is r.
    curzon_chazy_pullback checks the second against the first, and curzon_chazy.md beside this
    file is the derivation."""
    if system == "weyl":
        coords, parameters = ["t", "\\rho", "\\phi", "z"], ["m", "R = \\sqrt{\\rho^2 + z^2}"]

        def line(c2):
            return (f"ds^2 = -e^{{-2m/R}}{c2}dt^2 + e^{{2m/R}}\\left(e^{{-m^2\\rho^2/R^4}}\\left(d\\rho^2 + dz^2\\right)"
                    " + \\rho^2d\\phi^2\\right)")
        probe = vm.Reader(coords, parameters, ())
        rho, z, m = probe.symbol["\\rho"], probe.symbol["z"], probe.parameters["m"]
        R = sp.Symbol("R", positive=True)
        domains = ["t \\in (-\\infty, \\infty)", "\\rho \\in [0, \\infty)", "\\phi \\in [0, 2\\pi)",
                   "z \\in (-\\infty, \\infty)", "(\\rho, z) \\neq (0, 0) \\;\\text{(the singularity)}"]
        lead, pretty, name = [m, rho, R, z], weyl_distance(rho, z, R), "Weyl"
        metric = {("\\rho", "\\rho"): "e^{2m/R}e^{-m^2\\rho^2/R^4}", ("z", "z"): "e^{2m/R}e^{-m^2\\rho^2/R^4}"}
        inverse = {("\\rho", "\\rho"): "e^{-2m/R}e^{m^2\\rho^2/R^4}", ("z", "z"): "e^{-2m/R}e^{m^2\\rho^2/R^4}"}
        rewrite = [("\\dfrac{2R^3 - 2m\\,\\rho^2}{\\rho\\,R^3}", "\\dfrac{2\\left(R^3 - m\\,\\rho^2\\right)}{\\rho\\,R^3}")]
    else:
        coords, parameters = ["t", "r", "\\theta", "\\phi"], ["m"]

        def line(c2):
            return (f"ds^2 = -e^{{-2m/r}}{c2}dt^2 + e^{{2m/r}}\\left(e^{{-m^2\\sin^2\\theta/r^2}}"
                    "\\left(dr^2 + r^2d\\theta^2\\right) + r^2\\sin^2\\theta\\,d\\phi^2\\right)")
        probe = vm.Reader(coords, parameters, ())
        r, th, m = probe.symbol["r"], probe.symbol["\\theta"], probe.parameters["m"]
        domains = ["t \\in (-\\infty, \\infty)", "r \\in (0, \\infty)", "\\theta \\in [0, \\pi]",
                   "\\phi \\in [0, 2\\pi)", "r = 0 \\;\\text{(the singularity)}"]
        lead, pretty, name = [m, r, sp.cos(th), sp.sin(th)], weyl_distance(None, None, None), "Spherical"
        metric = {}
        inverse = {("\\phi", "\\phi"): "\\dfrac{e^{-2m/r}}{r^2\\sin^2\\theta}"}
        rewrite = [("\\left(r\\,e^{2m/r}\\sin^2\\theta - m\\,e^{2m/r}\\sin^2\\theta\\right)",
                    "\\left(r - m\\right)e^{2m/r}\\sin^2\\theta"),
                   ("\\dfrac{2r - 2m}{r^2}", "\\dfrac{2\\left(r - m\\right)}{r^2}"),
                   ("\\dfrac{2m^2\\sin^2\\theta - 2m\\,r + 2r^2}{r^3}",
                    "\\dfrac{2\\left(m^2\\sin^2\\theta - m\\,r + r^2\\right)}{r^3}")]
    return {
        "metric_id": "curzon_chazy",
        "system": {"id": system, "name": name, "coords": coords, "domains": domains, "parameters": parameters,
                   "line_element": line("c^2")},
        "chart_line_element": line(""),
        "printer": {"lead": lead},
        "pretty": pretty,
        "components": {"metric_components": metric, "inverse_metric_components": inverse},
        "rewrite": rewrite,
        "check": curzon_chazy_check,
    }


def weyl_distance(rho, z, R):
    """A pretty printer for a Weyl chart that names R = sqrt(rho^2 + z^2): every power of
    rho^2 + z^2 is written as a power of R, every even power of z left over as R^2 - rho^2, so a
    value is a polynomial in rho and R with at most one z in front, and the exponentials are
    gathered and written one for each term of their exponent, as the line element writes them.
    Without a chart's rho, z and R it only gathers the exponentials."""

    def rational(e):
        if R is None:
            return sp.factor(e)
        square = rho ** 2 + z ** 2
        e = e.replace(lambda p: p.is_Pow and sp.expand(p.base - square) == 0, lambda p: R ** (2 * p.exp))
        sides = []
        for side in sp.fraction(sp.together(e)):
            poly = sp.Poly(sp.expand(side), z)
            low = min(k for (k,), _ in poly.terms())
            rest = sum(c * z ** ((k - low) % 2) * (R ** 2 - rho ** 2) ** ((k - low) // 2) for (k,), c in poly.terms())
            sides.append(z ** low * sp.factor(sp.expand(rest)))
        return sides[0] / sides[1]

    def pretty(value):
        value = sp.sympify(value)
        if value == 0:
            return value
        exponent, rest = sp.Integer(0), sp.Integer(1)
        for f in sp.Mul.make_args(sp.powsimp(sp.factor(value))):
            base, k = (f.base, f.exp) if f.is_Pow else (f, sp.Integer(1))
            if isinstance(base, sp.exp):
                exponent += base.args[0] * k
            else:
                rest *= f
        if rest.atoms(sp.exp):
            raise AssertionError(f"an exponential stands inside a sum of {value}")
        return sp.Mul(rational(rest), *[sp.exp(rational(term)) for term in sp.Add.make_args(sp.expand(exponent))
                                        if term != 0])

    return pretty


def curzon_chazy_check(chart):
    """Both charts solve the vacuum equations, and Weyl's two functions solve Weyl's: psi = -m/R
    is harmonic in the flat space of rho, phi and z, and gamma is its quadrature,
    d_rho gamma = rho((d_rho psi)^2 - (d_z psi)^2) and d_z gamma = 2 rho d_rho psi d_z psi."""
    ricci = chart.geo.ricci_ll()
    if any(vm.norm(ricci[a][b]) != 0 for a in range(4) for b in range(4)):
        raise AssertionError("curzon_chazy: the Ricci tensor does not vanish")
    if chart.coords_tex[1] != "\\rho":
        return
    rho, z = chart.symbols[1], chart.symbols[3]
    m, R = chart.reader.parameters["m"], chart.reader.parameters["R"]
    psi, gamma = -m / R, -m ** 2 * rho ** 2 / (2 * R ** 4)
    checks = {"Laplace": sp.diff(rho * sp.diff(psi, rho), rho) / rho + sp.diff(psi, z, 2),
              "d_rho gamma": sp.diff(gamma, rho) - rho * (sp.diff(psi, rho) ** 2 - sp.diff(psi, z) ** 2),
              "d_z gamma": sp.diff(gamma, z) - 2 * rho * sp.diff(psi, rho) * sp.diff(psi, z),
              "g_tt": chart.geo.g[0, 0] + sp.exp(2 * psi),
              "g_rhorho": chart.geo.g[1, 1] - sp.exp(2 * gamma - 2 * psi)}
    for label, value in checks.items():
        if vm.norm(value) != 0:
            raise AssertionError(f"curzon_chazy: Weyl's {label} fails")


def curzon_chazy_pullback():
    """The spherical chart is Weyl's with rho = r sin(theta), z = r cos(theta): the pullback of
    Weyl's metric is the spherical one in every slot."""
    weyl, sph = (cp.Chart(spec["system"]["coords"], spec["system"]["parameters"], spec["chart_line_element"])
                 for spec in (curzon_chazy("weyl"), curzon_chazy("spherical")))
    t, r, th, ph = sph.symbols
    T, rho, phi, z = weyl.symbols
    at = {T: t, rho: r * sp.sin(th), phi: ph, z: r * sp.cos(th)}
    images = [at[x] for x in weyl.symbols]
    jacobian = sp.Matrix(4, 4, lambda a, b: sp.diff(images[a], sph.symbols[b]))
    source = weyl.geo.g.subs(at, simultaneous=True).applyfunc(lambda e: sp.powdenest(sp.simplify(e), force=True))
    pulled = jacobian.T * source * jacobian
    for a in range(4):
        for b in range(a, 4):
            if sp.simplify(pulled[a, b] - sph.geo.g[a, b]) != 0:
                raise AssertionError(f"curzon_chazy: the pullback of Weyl's chart misses the spherical chart in slot "
                                     f"{sph.coords_tex[a]}{sph.coords_tex[b]}")


def curzon_chazy_charts():
    curzon_chazy_pullback()
    return [curzon_chazy("weyl"), curzon_chazy("spherical")]


CHARTS["curzon_chazy"] = curzon_chazy_charts


# -- Zipoy-Voorhees ----------------------------------------------------------------------

def named_powers(names, state, merges=()):
    """A pretty printer for a chart whose values are each a rational function times powers, with
    exponents that hold a parameter, of a few positive factors, where the chart names the ratios
    those factors come in: Zipoy and Voorhees's metric carries (1 - 2m/r)^q and
    (1 + m^2 sin^2(theta)/(r^2 - 2mr))^{q(2 + q)}, and names the two bases f and h. `names` maps
    the text of each name to its powers of the irreducible factors, as f = (r - 2m)/r is
    {r - 2m: 1, r: -1}. The part of every exponent that holds a parameter is written on the names,
    as f^{2q}h^{q(2 + q)}, below the line where every term of the exponent is negative, and the
    whole part stays in the rational function, factored, with each pair of `merges` written as
    its product, as named_factors writes them. `state["printer"]` is the chart's printer, which
    the chart's check puts there before anything is printed, and `pretty.tidy` is the factoring,
    for a chart that collects its numerators."""
    bases = []
    for powers in names.values():
        bases += [b for b in powers if b not in bases]
    unknowns = {name: sp.Dummy() for name in names}
    table = {}
    tidy = named_factors([], merges)

    def base_of(expression):
        for b in bases:
            if sp.expand(expression - b) == 0:
                return b
        return None

    def exponent_text(printer, e):
        factored = sp.factor(e)
        sums = [f for f in sp.Mul.make_args(factored) if f.is_Add]
        if len(sums) == 1 and len(sp.Mul.make_args(factored)) > 1 and not factored.could_extract_minus_sign():
            return printer.term(factored).replace("\\left(", "(").replace("\\right)", ")")
        return printer.positive_first(sp.expand(e))

    def pretty(value):
        value = sp.factor(sp.sympify(value))
        symbolic, rest = {b: sp.Integer(0) for b in bases}, sp.Integer(1)
        for f in sp.Mul.make_args(value):
            base, k = f.as_base_exp()
            if k.is_Number:
                rest *= f
                continue
            b = base_of(base)
            if b is None:
                raise AssertionError(f"{f} is a power of none of {bases}")
            whole, part = sp.expand(k).as_coeff_Add()
            symbolic[b] += part
            rest *= b ** whole
        if any(a.func is sp.Pow and not a.exp.is_Number for a in sp.preorder_traversal(rest)):
            raise AssertionError(f"{value} is not one product of powers of {bases}")
        equations = [sum(unknowns[name] * powers.get(b, 0) for name, powers in names.items()) - symbolic[b] for b in bases]
        solution = sp.solve(equations, list(unknowns.values()), dict=True)
        if len(solution) != 1 or set(solution[0]) != set(unknowns.values()):
            raise AssertionError(f"the powers of {value} are not powers of {list(names)}")
        printer = state["printer"]
        out = tidy(rest)
        # The names on one side of the line are written together, in the order the chart gives them.
        sides = {1: "", -1: ""}
        for name in names:
            e = sp.expand(solution[0][unknowns[name]])
            if e == 0:
                continue
            below = all(c.as_coeff_Mul()[0].is_negative for c in sp.Add.make_args(e))
            sides[-1 if below else 1] += name + "^{" + exponent_text(printer, -e if below else e) + "}"
        for side, text in sides.items():
            if text:
                placeholder = table.setdefault(text, sp.Symbol(f"NAMEDPOWER{len(table)}", positive=True))
                printer.overrides[placeholder] = text
                out *= placeholder ** side
        return out

    pretty.tidy = tidy
    return pretty


def zipoy_voorhees(system):
    """Zipoy and Voorhees's deformed mass, the member of Weyl's class whose potential is that of
    a uniform rod of length 2m and mass per unit length delta/2: in Quevedo's coordinates, which
    are Schwarzschild's at q = delta - 1 = 0, and in the prolate spheroidal coordinates x = r/m - 1,
    y = cos(theta) of Voorhees. Both charts name f = 1 - 2m/r and h = e^{-2 gamma/delta^2} and are
    printed around their powers by named_powers; zipoy_voorhees_pullback checks the first against
    the second, and zipoy_voorhees.md beside this file is the derivation."""
    state = {}
    if system == "spherical":
        coords = ["t", "r", "\\theta", "\\phi"]
        parameters = ["m", "q", "f = 1 - \\dfrac{2m}{r}", "h = 1 + \\dfrac{m^2\\sin^2\\theta}{r^2 - 2mr}"]

        def line(c2):
            return (f"ds^2 = -f^{{1+q}}{c2}dt^2 + f^{{-q}}\\left(h^{{-q(2+q)}}\\left(\\dfrac{{dr^2}}{{f}} + r^2d\\theta^2\\right)"
                    " + r^2\\sin^2\\theta\\,d\\phi^2\\right)")
        probe = vm.Reader(coords, parameters, ())
        r, th, m, q = probe.symbol["r"], probe.symbol["\\theta"], probe.parameters["m"], probe.parameters["q"]
        sigma = r ** 2 - 2 * m * r + m ** 2 * sp.sin(th) ** 2
        names = {"f": {r - 2 * m: 1, r: -1}, "h": {sigma: 1, r: -1, r - 2 * m: -1}}
        domains = ["t \\in (-\\infty, \\infty)", "r \\in (2m, \\infty)", "\\theta \\in [0, \\pi]", "\\phi \\in [0, 2\\pi)",
                   "r = 2m \\;\\text{(the singularity, for } q \\neq 0\\text{)}"]
        name, parameter = "Spherical", q
        printer = {"lead": [r, m, sp.cos(th), sp.sin(th)], "factors": [m, q, r], "rising": [q],
                   "collect": lambda poly, pr: cp.collect_by(poly, [sp.sin(th)], pr)}
        metric = {("t", "t"): "-f^{1+q}", ("r", "r"): "\\dfrac{1}{f^{1+q}h^{q(2+q)}}",
                  ("\\theta", "\\theta"): "\\dfrac{r^2}{f^{q}h^{q(2+q)}}", ("\\phi", "\\phi"): "\\dfrac{r^2\\sin^2\\theta}{f^{q}}"}
        inverse = {("t", "t"): "-\\dfrac{1}{f^{1+q}}", ("r", "r"): "f^{1+q}h^{q(2+q)}",
                   ("\\theta", "\\theta"): "\\dfrac{f^{q}h^{q(2+q)}}{r^2}", ("\\phi", "\\phi"): "\\dfrac{f^{q}}{r^2\\sin^2\\theta}"}
        merges = []
        # Quevedo's grouping of the Kretschmann scalar.
        kretschmann = ("\\dfrac{16m^2\\left(1 + q\\right)^2\\left(3\\left(r - 2m - m\\,q\\right)^2"
                       "\\left(r^2 - 2m\\,r + m^2\\sin^2\\theta\\right) + m^2\\,q\\left(2 + q\\right)"
                       "\\left(m^2\\,q\\left(2 + q\\right) + 3\\left(r - m\\right)\\left(r - 2m - m\\,q\\right)\\right)"
                       "\\sin^2\\theta\\right)f^{2q}h^{2q(2 + q)}}"
                       "{r^6\\left(r - 2m\\right)^2\\left(r^2 - 2m\\,r + m^2\\sin^2\\theta\\right)}")
    else:
        coords = ["t", "x", "y", "\\phi"]
        parameters = ["m", "\\delta", "f = \\dfrac{x - 1}{x + 1}", "h = \\dfrac{x^2 - y^2}{x^2 - 1}"]

        def line(c2):
            return (f"ds^2 = -f^{{\\delta}}{c2}dt^2 + \\dfrac{{m^2}}{{f^{{\\delta}}}}\\left(\\dfrac{{x^2 - y^2}}{{h^{{\\delta^2}}}}"
                    "\\left(\\dfrac{dx^2}{x^2 - 1} + \\dfrac{dy^2}{1 - y^2}\\right) + \\left(x^2 - 1\\right)\\left(1 - y^2\\right)d\\phi^2\\right)")
        probe = vm.Reader(coords, parameters, ())
        x, y, m, delta = probe.symbol["x"], probe.symbol["y"], probe.parameters["m"], probe.parameters["delta"]
        names = {"f": {x - 1: 1, x + 1: -1}, "h": {x - y: 1, x + y: 1, x - 1: -1, x + 1: -1}}
        domains = ["t \\in (-\\infty, \\infty)", "x \\in (1, \\infty)", "y \\in [-1, 1]", "\\phi \\in [0, 2\\pi)",
                   "x = 1 \\;\\text{(the singularity, for } \\delta \\neq 1\\text{)}"]
        name, parameter = "Prolate spheroidal", delta
        merges = [(x - 1, x + 1, x ** 2 - 1), (1 - y, 1 + y, 1 - y ** 2), (x - y, x + y, x ** 2 - y ** 2),
                  (delta - 1, delta + 1, delta ** 2 - 1)]

        def collect(poly, pr):
            # Each power of delta with its coefficient factored as the denominators are.
            return cp.Sum([(c, pretty.tidy(rest)) for c, rest in cp.collect_by(poly, [delta], pr).terms])
        printer = {"lead": [delta, x, y], "factors": [m, delta, x, y], "rising": [y], "collect": collect}
        metric, inverse, kretschmann = {}, {}, None

    def check(chart):
        state["printer"] = chart.printer
        ricci = chart.geo.ricci_ll()
        if any(vm.norm(ricci[a][b]) != 0 for a in range(4) for b in range(4)):
            raise AssertionError("zipoy_voorhees: the Ricci tensor does not vanish")

    pretty = named_powers(names, state, merges)
    return {
        **({"kretschmann": kretschmann} if kretschmann else {}),
        "metric_id": "zipoy_voorhees",
        "system": {"id": "spherical" if system == "spherical" else "prolate_spheroidal", "name": name, "coords": coords,
                   "domains": domains, "parameters": parameters, "line_element": line("c^2")},
        "chart_line_element": line(""),
        "printer": printer,
        "pretty": pretty,
        "components": {"metric_components": metric, "inverse_metric_components": inverse},
        "check": check,
    }


def zipoy_voorhees_pullback():
    """The spherical chart is the prolate spheroidal one with x = r/m - 1, y = cos(theta), and
    delta = 1 + q: the pullback of the second metric is the first in every slot."""
    sph, pro = (cp.Chart(spec["system"]["coords"], spec["system"]["parameters"], spec["chart_line_element"])
                for spec in (zipoy_voorhees("spherical"), zipoy_voorhees("prolate_spheroidal")))
    t, r, th, ph = sph.symbols
    T, x, y, phi = pro.symbols
    m, q = sph.reader.parameters["m"], sph.reader.parameters["q"]
    at = {T: t, x: r / m - 1, y: sp.cos(th), phi: ph, pro.reader.parameters["delta"]: 1 + q}
    images = [at[u] for u in pro.symbols]
    jacobian = sp.Matrix(4, 4, lambda a, b: sp.diff(images[a], sph.symbols[b]))
    # x - y and x + y each come to a factor linear in cos(theta), and the spherical chart writes
    # their product, x^2 - y^2 = (r^2 - 2mr + m^2 sin^2(theta))/m^2, as one factor: every power of
    # x - y is written through the product W before the coordinates are replaced.
    W = sp.Dummy("W", positive=True)
    source = pro.geo.g.applyfunc(lambda e: vm.norm(e.replace(
        lambda u: u.func is sp.Pow and u.base == x - y and not u.exp.is_Number,
        lambda u: (W / (x + y)) ** u.exp)))
    at[W] = (r ** 2 - 2 * m * r + m ** 2 * sp.sin(th) ** 2) / m ** 2
    source = source.subs(at, simultaneous=True)
    pulled = jacobian.T * source * jacobian
    for a in range(4):
        for b in range(a, 4):
            if vm.norm(pulled[a, b] - sph.geo.g[a, b]) != 0:
                raise AssertionError(f"zipoy_voorhees: the pullback of the prolate spheroidal chart misses the "
                                     f"spherical chart in slot {sph.coords_tex[a]}{sph.coords_tex[b]}")


def zipoy_voorhees_charts():
    zipoy_voorhees_pullback()
    return [zipoy_voorhees("spherical"), zipoy_voorhees("prolate_spheroidal")]


CHARTS["zipoy_voorhees"] = zipoy_voorhees_charts


# -- A black hole threaded by a cosmic string ------------------------------------------

def string_black_hole(system_id):
    """Aryal, Ford and Vilenkin's black hole on a cosmic string: Schwarzschild's metric with
    g_phiphi multiplied by b^2, b = 1 - 4G mu/c^2, in the static chart, in the chart whose angle
    b phi runs over 2 pi b, where the line element is Schwarzschild's own, and in the two
    Eddington-Finkelstein charts built on Schwarzschild's tortoise coordinate. Every value is
    printed around r - r_s, as Schwarzschild's are. The wedge chart is checked, slot by slot,
    to be the static chart pulled back through phi = phi~/b."""
    f = "\\left(1 - \\dfrac{r_s}{r}\\right)"
    bare = "1 - \\dfrac{r_s}{r}"
    parameters = ["r_s", "b"]
    kretschmann = "\\dfrac{12r_s^2}{r^6}"
    string = ["\\theta = 0, \\pi \\;\\text{(the string, a conical singularity)}"]
    horizon = ["r = r_s \\;\\text{(the horizon)}"]
    extra = {}
    if system_id == "wedge":
        angle = "\\tilde\\phi"
        sphere = " + r^2\\left(d\\theta^2 + \\sin^2\\theta\\,d\\tilde\\phi^2\\right)"
        angles = ["\\theta \\in [0, \\pi]", "\\tilde\\phi \\in [0, 2\\pi b)"]
        extra["check"] = string_black_hole_wedge
    else:
        angle = "\\phi"
        sphere = " + r^2\\left(d\\theta^2 + b^2\\sin^2\\theta\\,d\\phi^2\\right)"
        angles = ["\\theta \\in [0, \\pi]", "\\phi \\in [0, 2\\pi)"]
    if system_id in ("static", "wedge"):
        coords = ["t", "r", "\\theta", angle]
        name = "Static" if system_id == "static" else "Static with the Wedge Removed"
        line = "ds^2 = -" + f + "c^2dt^2 + \\dfrac{dr^2}{" + bare + "}" + sphere
        chart_line = "ds^2 = -" + f + "dt^2 + \\dfrac{dr^2}{" + bare + "}" + sphere
        domains = ["t \\in (-\\infty, \\infty)", "r \\in (r_s, \\infty)"] + angles + string
        components = {"metric_components": {("t", "t"): "-" + f, ("r", "r"): f + "^{-1}"},
                      "inverse_metric_components": {("t", "t"): "-" + f + "^{-1}", ("r", "r"): bare}}
    else:
        null, sign = ("u", "-") if system_id == "eddington_finkelstein_outgoing" else ("v", "+")
        coords = [null, "r", "\\theta", angle]
        name = ("Outgoing" if null == "u" else "Ingoing") + " Eddington-Finkelstein"
        line = chart_line = "ds^2 = -" + f + "d" + null + "^2 " + sign + " 2\\,d" + null + "\\,dr" + sphere
        one = "-1" if null == "u" else "1"
        domains = [null + " \\in (-\\infty, \\infty)", "r \\in (0, \\infty)"] + angles + string + horizon
        components = {"metric_components": {(null, null): "-" + f, (null, "r"): one, ("r", null): one},
                      "inverse_metric_components": {(null, "r"): one, ("r", null): one, ("r", "r"): bare}}
    probe = vm.Reader(coords, parameters, ())
    r, rs, b = probe.symbol["r"], probe.parameters["r_s"], probe.parameters["b"]
    return {
        "metric_id": "string_black_hole",
        "system": {"id": system_id, "name": name, "coords": coords, "domains": domains,
                   "parameters": parameters, "line_element": line},
        "chart_line_element": chart_line,
        "printer": {"rising": [rs], "lead": [b, r, rs], "flip": False},
        "components": components,
        "kretschmann": kretschmann,
        **extra,
    }


def string_black_hole_wedge(chart):
    """J^T g J, with g the static chart and J the Jacobian of phi = phi~/b, against the wedge
    chart's metric, which is Schwarzschild's, in every slot."""
    spec = string_black_hole("static")
    source = cp.Chart(spec["system"]["coords"], spec["system"]["parameters"], spec["chart_line_element"])
    b = chart.reader.parameters["b"]
    at = {source.reader.parameters[n]: chart.reader.parameters[n] for n in ("r_s", "b")}
    at.update(dict(zip(source.symbols[:3], chart.symbols[:3])))
    at[source.symbols[3]] = chart.symbols[3] / b
    J = sp.diag(1, 1, 1, 1 / b)
    pulled = J.T * source.geo.g.subs(at) * J
    for i in range(4):
        for j in range(i, 4):
            if sp.simplify(pulled[i, j] - chart.geo.g[i, j]) != 0:
                raise AssertionError(f"string_black_hole: the static chart pulled back misses the wedge "
                                     f"chart in slot {chart.coords_tex[i]}{chart.coords_tex[j]}")
    if chart.geo.g.has(b):
        raise AssertionError("string_black_hole: the wedge chart's metric is not Schwarzschild's")


SBH_CHARTS = ["static", "wedge", "eddington_finkelstein_outgoing", "eddington_finkelstein_ingoing"]
CHARTS["string_black_hole"] = [lambda s=s: string_black_hole(s) for s in SBH_CHARTS]

# -- Schwarzschild-Tangherlini ---------------------------------------------------------

TANGHERLINI_CHARTS = ["spherical", "eddington_finkelstein_outgoing", "eddington_finkelstein_ingoing", "spherical_six"]


def tangherlini(system_id):
    """Tangherlini's black hole, f = 1 - (r_h/r)^(D-3), in D = 5 in the static chart and the two
    Eddington-Finkelstein charts built on its tortoise coordinate
    r_* = r + (r_h/2) ln|(r - r_h)/(r + r_h)|, and in D = 6 in the static chart. The angles are
    the hyperspherical ones, each sphere's line element the next angle's plus its sine squared
    times the sphere below. The metric and its inverse are written as the line element writes f,
    and the Kretschmann scalar is (D-1)(D-2)^2(D-3) r_h^(2D-6)/r^(2D-2), Schwarzschild's
    12r_s^2/r^6 at D = 4. tangherlini.md derives each chart."""
    six = system_id == "spherical_six"
    n = 3 if six else 2
    f = "\\left(1 - \\dfrac{r_h^%d}{r^%d}\\right)" % (n, n)
    bare = "1 - \\dfrac{r_h^%d}{r^%d}" % (n, n)
    if six:
        angles = ["\\chi", "\\psi", "\\theta", "\\phi"]
        sphere = (" + r^2\\left(d\\chi^2 + \\sin^2\\chi\\,d\\psi^2 + \\sin^2\\chi\\sin^2\\psi\\,d\\theta^2"
                  " + \\sin^2\\chi\\sin^2\\psi\\sin^2\\theta\\,d\\phi^2\\right)")
    else:
        angles = ["\\psi", "\\theta", "\\phi"]
        sphere = (" + r^2\\left(d\\psi^2 + \\sin^2\\psi\\,d\\theta^2"
                  " + \\sin^2\\psi\\sin^2\\theta\\,d\\phi^2\\right)")
    domains = [a + " \\in [0, \\pi]" for a in angles[:-1]] + ["\\phi \\in [0, 2\\pi)"]
    if system_id.startswith("spherical"):
        coords = ["t", "r"] + angles
        name = "Hyperspherical, " + ("Six" if six else "Five") + " Dimensions"
        radial = "r \\in (r_h, \\infty)"
        line = "ds^2 = -" + f + "c^2dt^2 + \\dfrac{dr^2}{" + bare + "}" + sphere
        chart_line = "ds^2 = -" + f + "dt^2 + \\dfrac{dr^2}{" + bare + "}" + sphere
        metric = {("t", "t"): "-" + f, ("r", "r"): f + "^{-1}"}
        inverse = {("t", "t"): "-" + f + "^{-1}", ("r", "r"): bare}
    else:
        null, sign = ("u", "-") if system_id == "eddington_finkelstein_outgoing" else ("v", "+")
        coords = [null, "r"] + angles
        name = ("Outgoing" if null == "u" else "Ingoing") + " Eddington-Finkelstein, Five Dimensions"
        radial = "r \\in (0, \\infty)"
        line = "ds^2 = -" + f + "d" + null + "^2 " + sign + " 2\\,d" + null + "\\,dr" + sphere
        chart_line = line
        one = "-1" if null == "u" else "1"
        metric = {(null, null): "-" + f, (null, "r"): one, ("r", null): one}
        inverse = {(null, "r"): one, ("r", null): one, ("r", "r"): bare}
    parameters = ["r_h"]
    probe = vm.Reader(coords, parameters, ())
    r, rh = probe.symbol["r"], probe.parameters["r_h"]
    return {
        "metric_id": "tangherlini",
        "system": {"id": system_id, "name": name, "coords": coords,
                   "domains": [coords[0] + " \\in (-\\infty, \\infty)", radial] + domains,
                   "parameters": parameters, "line_element": line},
        "chart_line_element": chart_line,
        "printer": {"rising": [rh], "lead": [r, rh], "flip": False},
        "components": {"metric_components": metric, "inverse_metric_components": inverse},
        "kretschmann": "\\dfrac{240r_h^6}{r^{10}}" if six else "\\dfrac{72r_h^4}{r^8}",
        # The printer factors r^n - r_h^n; every value keeps it whole, as f itself is written.
        "rewrite": [("\\left(r - r_h\\right)\\left(r^2 + r\\,r_h + r_h^2\\right)", "\\left(r^3 - r_h^3\\right)")] if six
        else [("\\left(r + r_h\\right)\\left(r - r_h\\right)", "\\left(r^2 - r_h^2\\right)")],
    }


CHARTS["tangherlini"] = [lambda s=s: tangherlini(s) for s in TANGHERLINI_CHARTS]


# -- Gott's time machine ---------------------------------------------------------------

GOTT_CHARTS = ("centre_of_momentum", "string_rest", "grant_rindler", "grant_milne")


def gott_time_machine(system):
    """The four charts of Gott's two strings, every one flat: the inertial chart of the centre of
    momentum frame, whose wedges and identifications are its domains; the conical chart of one
    string's rest frame, the cosmic string's own; and the Rindler and Milne charts of the
    Minkowski space James Grant showed the spacetime to be away from the strings, identified
    under a boost of rapidity a and a shift b along its axis. gott_time_machine.md derives a
    and b from the two strings' rotations."""
    reals = "(-\\infty, \\infty)"
    flat = "\\text{flat: every curvature tensor vanishes}"
    deficit = "\\left(1 - \\dfrac{4G\\mu}{c^2}\\right)"
    charts = {
        "centre_of_momentum": {
            "name": "Centre of Momentum", "coords": ["t", "x", "y", "z"],
            "parameters": ["\\mu", "G", "v", "d", "\\alpha", "\\gamma"],
            "domains": ["t \\in " + reals, "x \\in " + reals, "y \\in " + reals, "z \\in " + reals,
                        "x = \\pm vt,\\; y = \\pm d \\;\\text{(the strings)}",
                        "\\gamma|x \\mp vt| < \\pm(y \\mp d)\\tan\\alpha \\;\\text{(the wedges removed)}",
                        "\\gamma\\sin\\alpha > 1 \\;\\text{(closed timelike curves)}", flat],
            "line_element": "ds^2 = -c^2dt^2 + dx^2 + dy^2 + dz^2",
            "chart": "ds^2 = -dt^2 + dx^2 + dy^2 + dz^2"},
        "string_rest": {
            "name": "String Rest Frame", "coords": ["t", "r", "\\phi", "z"], "parameters": ["\\mu", "G", "d"],
            "domains": ["t \\in " + reals, "r \\in (0, \\infty)", "\\phi \\in [0, 2\\pi)", "z \\in " + reals,
                        "r = 0 \\;\\text{(the string, a conical singularity)}",
                        "r\\cos\\left[" + deficit + "(\\phi - \\pi)\\right] \\le d "
                        "\\;\\text{(the half space } y \\ge 0\\text{)}", flat],
            "line_element": "ds^2 = -c^2dt^2 + dr^2 + " + deficit + "^2 r^2 d\\phi^2 + dz^2",
            "chart": "ds^2 = -dt^2 + dr^2 + " + deficit + "^2 r^2 d\\phi^2 + dz^2"},
        "grant_rindler": {
            "name": "Grant Rindler", "coords": ["\\eta", "\\xi", "Y", "z"], "parameters": ["a", "b"],
            "domains": ["\\eta \\in [0, a)", "\\xi \\in (0, \\infty)", "Y \\in " + reals, "z \\in " + reals,
                        "(\\eta, Y) \\sim (\\eta + a, Y + b)", "\\xi = 0 \\;\\text{(chronology horizon)}",
                        "\\xi = \\dfrac{nb}{2\\sinh(na/2)} \\;\\text{(the } n\\text{th polarised hypersurface)}",
                        flat],
            "line_element": "ds^2 = -\\xi^2d\\eta^2 + d\\xi^2 + dY^2 + dz^2",
            "chart": "ds^2 = -\\xi^2d\\eta^2 + d\\xi^2 + dY^2 + dz^2"},
        "grant_milne": {
            "name": "Grant Milne", "coords": ["\\tau", "\\chi", "Y", "z"], "parameters": ["a", "b"],
            "domains": ["\\tau \\in " + reals, "\\chi \\in [0, a)", "Y \\in " + reals, "z \\in " + reals,
                        "(\\chi, Y) \\sim (\\chi + a, Y + b)", "\\tau = 0 \\;\\text{(chronology horizon)}", flat],
            "line_element": "ds^2 = -c^2d\\tau^2 + c^2\\tau^2d\\chi^2 + dY^2 + dz^2",
            "chart": "ds^2 = -d\\tau^2 + \\tau^2d\\chi^2 + dY^2 + dz^2", "time": "\\tau"},
    }
    chart = charts[system]
    probe = vm.Reader(chart["coords"], chart["parameters"], ())
    return {
        "metric_id": "gott_time_machine",
        "system": {"id": system, "name": chart["name"], "coords": chart["coords"], "domains": chart["domains"],
                   "parameters": chart["parameters"], "line_element": chart["line_element"]},
        "chart_line_element": chart["chart"],
        "printer": {"lead": [probe.c, probe.symbol[chart["coords"][0]], *probe.parameters.values()]},
        **({"time": chart["time"]} if "time" in chart else {}),
        # The printer collects the deficit over c^4; the cone keeps it as the cosmic string writes it.
        **({"rewrite": [("\\dfrac{c^4}{r^2\\left(c^2 - 4\\mu\\,G\\right)^2}", "\\dfrac{1}{" + deficit + "^2 r^2}"),
                        ("\\dfrac{r^2\\left(c^2 - 4\\mu\\,G\\right)^2}{c^4}", deficit + "^2 r^2"),
                        ("\\dfrac{r\\left(c^2 - 4\\mu\\,G\\right)^2}{c^4}", deficit + "^2 r")]}
           if system == "string_rest" else {}),
    }


CHARTS["gott_time_machine"] = [lambda s=s: gott_time_machine(s) for s in GOTT_CHARTS]

# -- Gowdy -----------------------------------------------------------------------------

GOWDY_CHARTS = ["areal", "logarithmic", "sphere"]
GOWDY_PARAMETERS = {
    "areal": ["L", "P = P(t,\\theta)", "Q = Q(t,\\theta)", "\\lambda = \\lambda(t,\\theta)"],
    "logarithmic": ["L", "P = P(\\tau,\\theta)", "Q = Q(\\tau,\\theta)", "\\lambda = \\lambda(\\tau,\\theta)"],
    "sphere": ["L", "P = P(t,\\theta)", "Q = Q(t,\\theta)", "a = a(t,\\theta)"],
}
GOWDY_TORUS = "\\left(e^{P}\\left(d\\sigma + Q\\,d\\delta\\right)^2 + e^{-P}d\\delta^2\\right)"
GOWDY_LINES = {
    "areal": "ds^2 = L^2\\left(\\dfrac{e^{\\lambda/2}}{\\sqrt{t}}\\left(-dt^2 + d\\theta^2\\right) + t" + GOWDY_TORUS + "\\right)",
    "logarithmic": ("ds^2 = L^2\\left(e^{(\\lambda + \\tau)/2}\\left(-e^{-2\\tau}d\\tau^2 + d\\theta^2\\right)"
                    " + e^{-\\tau}" + GOWDY_TORUS + "\\right)"),
    "sphere": ("ds^2 = L^2\\left(e^{2a}\\left(-dt^2 + d\\theta^2\\right) + \\sin t\\,\\sin\\theta"
               + GOWDY_TORUS + "\\right)"),
}


def gowdy(system):
    """Gowdy's cosmologies in his own notation, ds^2 = L^2(e^{2a}(-dt^2 + dtheta^2) +
    R(e^P(dsigma + Q ddelta)^2 + e^{-P} ddelta^2)), the orbits of the two Killing vectors of
    area 4 pi^2 L^2 R: the torus, R = t and e^{2a} = e^{lambda/2}/sqrt(t), in the areal time t
    and in tau = -ln t, as Rendall and Weaver's equation (1) writes it, and the sphere and the
    handle S^2 x S^1, R = sin t sin theta. The length L multiplies each line element whole and
    every coordinate is a pure number. P, Q and lambda, or a, are free functions of the time and
    theta, and no component assumes a field equation; gowdy_check holds the field equations each
    chart's parameters state to making every Ricci component vanish, the logarithmic chart to
    being the areal one pulled back through t = e^{-tau}, and the sphere chart to holding the
    inside of Schwarzschild's horizon. gowdy.md beside this file is the derivation."""
    parameters = GOWDY_PARAMETERS[system]
    period = " \\in [0, 2\\pi)"
    if system == "logarithmic":
        coords = ["\\tau", "\\theta", "\\sigma", "\\delta"]
        domains = ["\\tau \\in (-\\infty, \\infty)"] + [x + period for x in coords[1:]]
        name = "Torus, Logarithmic Time"
    else:
        coords = ["t", "\\theta", "\\sigma", "\\delta"]
        if system == "areal":
            domains = ["t \\in (0, \\infty)"] + [x + period for x in coords[1:]]
            name = "Torus, Areal Time"
        else:
            domains = ["t \\in (0, \\pi)", "\\theta \\in [0, \\pi]", "\\sigma" + period, "\\delta" + period]
            name = "Sphere and Handle"
    probe = vm.Reader(coords, parameters, ())
    x, y = probe.symbol[coords[0]], probe.symbol[coords[1]]
    P, Q, third = (probe.parameters[n] for n in ("P", "Q", "a" if system == "sphere" else "lambda"))
    D = sp.Derivative
    lead = [probe.parameters["L"], third, P, Q]
    for f in (third, P, Q):
        lead += [D(f, x), D(f, y)]
    for f in (P, Q, third):
        lead += [D(f, (x, 2)), D(f, x, y), D(f, (y, 2))]
    lead.append(x)
    return {
        "metric_id": "gowdy",
        "system": {"id": system, "name": name, "coords": coords, "domains": domains,
                   "parameters": parameters, "line_element": GOWDY_LINES[system]},
        "chart_line_element": GOWDY_LINES[system],
        "printer": {"lead": lead, "factors": lead},
        "pretty": lambda value: sp.powsimp(sp.factor(sp.sympify(value)), combine="exp"),
        "check": gowdy_check,
    }


def gowdy_vacuum(chart):
    """The field equations each chart's parameters state, as replacements for the second time
    derivatives of the two wave amplitudes and for both first derivatives of the third function."""
    x, y = chart.symbols[:2]
    D = sp.Derivative
    P, Q = chart.reader.parameters["P"], chart.reader.parameters["Q"]
    e = sp.exp(2 * P)
    if "a" in chart.reader.parameters:
        cx, cy = sp.cos(x) / sp.sin(x), sp.cos(y) / sp.sin(y)
        waves = {D(P, (x, 2)): D(P, (y, 2)) - cx * D(P, x) + cy * D(P, y) + e * (D(Q, x) ** 2 - D(Q, y) ** 2),
                 D(Q, (x, 2)): D(Q, (y, 2)) - cx * D(Q, x) + cy * D(Q, y) - 2 * (D(P, x) * D(Q, x) - D(P, y) * D(Q, y))}
        m = (cx * cy + D(P, x) * D(P, y) + e * D(Q, x) * D(Q, y)) / 2
        n = (D(P, x) ** 2 + D(P, y) ** 2 + e * (D(Q, x) ** 2 + D(Q, y) ** 2) - cx ** 2 - cy ** 2) / 4 - 1
        det = cy ** 2 - cx ** 2
        return waves, chart.reader.parameters["a"], {x: (m * cy - n * cx) / det, y: (n * cy - m * cx) / det}
    lam = chart.reader.parameters["lambda"]
    if chart.coords_tex[0] == "t":
        waves = {D(P, (x, 2)): D(P, (y, 2)) - D(P, x) / x + e * (D(Q, x) ** 2 - D(Q, y) ** 2),
                 D(Q, (x, 2)): D(Q, (y, 2)) - D(Q, x) / x - 2 * (D(P, x) * D(Q, x) - D(P, y) * D(Q, y))}
        first = {x: x * (D(P, x) ** 2 + D(P, y) ** 2 + e * (D(Q, x) ** 2 + D(Q, y) ** 2)),
                 y: 2 * x * (D(P, x) * D(P, y) + e * D(Q, x) * D(Q, y))}
    else:
        w = sp.exp(-2 * x)
        waves = {D(P, (x, 2)): w * D(P, (y, 2)) + e * (D(Q, x) ** 2 - w * D(Q, y) ** 2),
                 D(Q, (x, 2)): w * D(Q, (y, 2)) - 2 * (D(P, x) * D(Q, x) - w * D(P, y) * D(Q, y))}
        first = {x: -(D(P, x) ** 2 + w * D(P, y) ** 2 + e * (D(Q, x) ** 2 + w * D(Q, y) ** 2)),
                 y: -2 * (D(P, x) * D(P, y) + e * D(Q, x) * D(Q, y))}
    return waves, lam, first


def gowdy_on_shell(chart, value):
    """`value` where the chart's field equations hold: the third function's derivatives replaced
    by the constraints and their derivatives, then every time derivative of the two amplitudes
    beyond the first by the wave equations, until none is left."""
    x, y = chart.symbols[:2]
    waves, third, first = gowdy_vacuum(chart)
    D = sp.Derivative
    value = sp.sympify(value).subs({D(third, (x, 2)): sp.diff(first[x], x), D(third, (y, 2)): sp.diff(first[y], y),
                                    D(third, x, y): sp.diff(first[x], y)}).doit()
    value = value.subs({D(third, x): first[x], D(third, y): first[y]}).doit()
    for _ in range(4):
        higher = {}
        for d in value.atoms(D):
            counts = dict(d.variable_count)
            for target, replacement in waves.items():
                if d.expr == target.expr and counts.get(x, 0) >= 2:
                    rest = [(v, n - (2 if v == x else 0)) for v, n in d.variable_count]
                    rest = [(v, n) for v, n in rest if n]
                    higher[d] = sp.diff(replacement, *rest) if rest else replacement
        if not higher:
            break
        value = value.subs(higher).doit()
    return value


def gowdy_vanishes(value):
    """Whether a value is zero: by the checker's canonical form, and by sympy's own
    simplification where that form leaves something standing."""
    try:
        if vm.norm(sp.together(value)) == 0:
            return True
    except Exception:  # noqa: BLE001  the canonical form refuses some shapes, which simplify takes
        pass
    return sp.simplify(value) == 0


def gowdy_check(chart):
    ricci = chart.geo.ricci_ll()
    for i in range(4):
        for j in range(i, 4):
            if not gowdy_vanishes(gowdy_on_shell(chart, ricci[i][j])):
                raise AssertionError(f"gowdy: the stated field equations leave R_{chart.coords_tex[i]}"
                                     f"{chart.coords_tex[j]} standing")
    x, y = chart.symbols[:2]
    if chart.coords_tex[0] == "\\tau":
        # The areal chart pulled back through t = e^{-tau}, the three functions one value at an event.
        spec = gowdy("areal")
        source = cp.Chart(spec["system"]["coords"], spec["system"]["parameters"], spec["chart_line_element"])
        values = sp.symbols("P0 Q0 l0")
        at = {source.reader.parameters[n]: v for n, v in zip(("P", "Q", "lambda"), values)}
        here = {chart.reader.parameters[n]: v for n, v in zip(("P", "Q", "lambda"), values)}
        here[chart.reader.parameters["L"]] = source.reader.parameters["L"]
        J = sp.diag(-sp.exp(-x), 1, 1, 1)
        move = dict(zip(source.symbols[1:], chart.symbols[1:]))
        move[source.symbols[0]] = sp.exp(-x)
        pulled = J.T * source.geo.g.subs(at).subs(move) * J
        if sp.simplify(pulled - chart.geo.g.subs(here)) != sp.zeros(4, 4):
            raise AssertionError("gowdy: the areal chart pulled back through t = e^(-tau) misses the logarithmic chart")
    if "a" in chart.reader.parameters:
        # The inside of Schwarzschild's horizon, r = L(1 - cos t) and T = L delta with r_s = 2L:
        # e^{2a} = (1 - cos t)^2 and e^P = (1 - cos t)^2 sin(theta)/sin(t), Q = 0.
        P, Q, a = (chart.reader.parameters[n] for n in ("P", "Q", "a"))
        hole = {P: sp.log((1 - sp.cos(x)) ** 2 * sp.sin(y) / sp.sin(x)), Q: sp.Integer(0), a: sp.log(1 - sp.cos(x))}
        for i in range(4):
            for j in range(i, 4):
                if sp.simplify(ricci[i][j].subs(hole).doit()) != 0:
                    raise AssertionError("gowdy: the inside of Schwarzschild's horizon is not a vacuum in the sphere chart")
        L, r = chart.reader.parameters["L"], 1 - sp.cos(x)
        inside = sp.diag(-L ** 2 * r ** 2, L ** 2 * r ** 2, L ** 2 * r ** 2 * sp.sin(y) ** 2,
                         L ** 2 * (2 / r - 1))
        if sp.simplify(chart.geo.g.subs(hole).doit() - inside) != sp.zeros(4, 4):
            raise AssertionError("gowdy: the sphere chart's member is not Schwarzschild's interior")


CHARTS["gowdy"] = [lambda s=s: gowdy(s) for s in GOWDY_CHARTS]


# -- Szekeres ----------------------------------------------------------------------------

SZ_E = ("\\dfrac{S}{2}\\left(\\dfrac{\\left(p - P\\right)^2 + \\left(q - Q\\right)^2}{S^2}"
        " + \\epsilon\\right)")
SZ_ORDERS = 4


class SzekeresForms:
    """Szekeres's function E = (S/2)(((p - P)^2 + (q - Q)^2)/S^2 + epsilon), held as a function
    of r, p and q while the tensors are built, and the two ways a value in it is written.

    E is quadratic in p and q, so d_p^2 E = d_q^2 E = 1/S, d_p d_q E = 0 and every third
    derivative along p and q vanishes, and 2E = S((d_p E)^2 + (d_q E)^2 + epsilon). `reduce`
    writes those in and removes E and its derivatives along r alone by the last relation, which
    leaves d_p E, d_q E and their derivatives along r, with no relation among them, so a value
    that vanishes for Szekeres's E is exactly zero. `pretty` factors a reduced value and writes
    each factor back in E and its derivatives along r, by the same relation, in whichever of
    four orders of elimination is shortest."""

    def __init__(self, reader):
        self.r, self.p, self.q = (reader.symbol[name] for name in ("r", "p", "q"))
        self.E, self.S, self.epsilon = (reader.parameters[name] for name in ("E", "S", "epsilon"))
        a, b = sp.Function("_a")(self.r), sp.Function("_b")(self.r)
        self.along_r = {}      # the k-th derivative of E along r, in d_p E, d_q E and their derivatives
        self.symbol = {}       # a generator as the plain symbol Poly takes
        relations = []
        for k in range(SZ_ORDERS):
            for letter, x in (("a", self.p), ("b", self.q)):
                self.symbol[self.derivative(x, k)] = sp.Symbol(f"_{letter}{k}")
            self.symbol[self.derivative(None, k)] = sp.Symbol(f"_e{k}")
            value = sp.diff(self.S / 2 * (a ** 2 + b ** 2 + self.epsilon), self.r, k)
            orders = {d: d.variable_count[0][1] for d in value.atoms(sp.Derivative) if d.expr in (a, b)}
            value = value.xreplace({d: self.derivative(self.p if d.expr == a else self.q, n)
                                    for d, n in orders.items()})
            value = value.xreplace({a: self.derivative(self.p, 0), b: self.derivative(self.q, 0)})
            self.along_r[k] = value
            relations.append(sp.expand(value.xreplace(self.symbol) - sp.Symbol(f"_e{k}")))
        self.back = {symbol: generator for generator, symbol in self.symbol.items()}
        self.eliminations = []
        for letter in "ba":
            eliminated = [sp.Symbol(f"_{letter}{k}") for k in reversed(range(SZ_ORDERS))]
            self.eliminations += [(eliminated, relations), (eliminated, relations[::-1])]

    def derivative(self, x, k):
        """d_r^k of E, or of d_x E where x is p or q."""
        variables = ([x] if x is not None else []) + ([(self.r, k)] if k else [])
        return sp.Derivative(self.E, *variables) if variables else self.E

    def reduce(self, value):
        value = sp.sympify(value)
        if not value.has(self.E):
            return vm.norm(value)
        written = {}
        for d in value.atoms(sp.Derivative):
            if d.expr != self.E:
                continue
            count = {x: 0 for x in (self.r, self.p, self.q)}
            for x, n in d.variable_count:
                count[x] += n
            i, j, k = count[self.p], count[self.q], count[self.r]
            if k >= SZ_ORDERS:
                raise AssertionError(f"szekeres: {d} is of a higher order along r than SZ_ORDERS allows")
            if i + j >= 3 or (i, j) == (1, 1):
                written[d] = sp.Integer(0)
            elif i + j == 2:
                written[d] = sp.diff(1 / self.S, self.r, k)
            elif i + j == 0:
                written[d] = self.along_r[k]
        # E itself is replaced only where it stands outside a derivative.
        held = {d: sp.Dummy() for d in value.atoms(sp.Derivative) if d not in written}
        for image in written.values():
            for d in image.atoms(sp.Derivative):
                held.setdefault(d, sp.Dummy())
        for d in self.along_r[0].atoms(sp.Derivative):
            held.setdefault(d, sp.Dummy())
        value = value.xreplace({d: image.xreplace(held) for d, image in written.items()}).xreplace(held)
        value = value.xreplace({self.E: self.along_r[0].xreplace(held)})
        return vm.norm(value.xreplace({dummy: d for d, dummy in held.items()}))

    def factor(self, polynomial):
        polynomial = sp.expand(polynomial.xreplace(self.symbol))
        best = polynomial
        for eliminated, relations in self.eliminations:
            if not any(polynomial.has(x) for x in eliminated):
                continue
            remainder = sp.reduced(polynomial, relations, *eliminated, order="lex")[1]
            if sp.count_ops(remainder) < sp.count_ops(best):
                best = remainder
        return sp.factor(best.xreplace(self.back))

    def pretty(self, value):
        out = sp.Integer(1)
        for f in sp.Mul.make_args(sp.factor(sp.sympify(value))):
            base, k = (f.base, f.exp) if f.is_Pow else (f, sp.Integer(1))
            out *= self.factor(base) ** k if base.is_Add else f
        return out


def szekeres_polar(theta):
    """A `pretty` for the axisymmetric chart. Its metric has a term in dr dtheta, so a value comes
    back with cos^2(theta) and sin^2(theta) mixed and neither side of its fraction factors. Here
    every even power of the sine is written in the cosine, the fraction cancelled, and each side
    factored as it stands, all in the cosine and all in the sine, the shortest of the three kept."""
    c, s = sp.Symbol("_c"), sp.Symbol("_s")
    held = {sp.cos(theta): c, sp.sin(theta): s}
    back = {c: sp.cos(theta), s: sp.sin(theta)}
    circle = s ** 2 + c ** 2 - 1

    def in_cosine(value):
        return sp.expand(value).xreplace(held).replace(
            lambda e: e.is_Pow and e.base == s and e.exp.is_Integer and e.exp > 1,
            lambda e: (1 - c ** 2) ** (e.exp // 2) * s ** (e.exp % 2))

    def side(polynomial):
        polynomial = sp.expand(polynomial)
        forms = [polynomial, sp.reduced(polynomial, [circle], s, c)[1], sp.reduced(polynomial, [circle], c, s)[1]]
        return min((sp.factor(form) for form in forms), key=sp.count_ops).xreplace(back)

    def pretty(value):
        numerator, denominator = sp.fraction(sp.together(sp.sympify(value)))
        numerator, denominator = sp.fraction(sp.cancel(in_cosine(numerator) / in_cosine(denominator)))
        return side(numerator) / side(denominator)
    return pretty


def szekeres(system_id):
    """Szekeres's dust cosmologies of the class that holds the Lemaitre-Tolman models, in the
    two charts its literature uses: Hellaby and Krasinski's, whose surfaces of constant t and r
    carry the stereographic coordinates p and q and the sign epsilon of their curvature, and the
    chart of the quasispherical case with an axis of symmetry, P and Q constant, in the polar
    angles p - P = S cot(theta/2) cos(phi), q - Q = S cot(theta/2) sin(phi). R, f, S, P and Q
    are left free in every tensor, so no component assumes a field equation; szekeres_check
    confirms before anything is written that with (d_t R)^2 = 2M/R + f the matter is dust at
    rest in the chart, and szekeres.md is the derivation."""
    if system_id == "stereographic":
        coords, name = ["t", "r", "p", "q"], "Stereographic"
        parameters = ["R = R(t,r)", "f = f(r)", "S = S(r)", "P = P(r)", "Q = Q(r)", "\\epsilon", "E = " + SZ_E]
        line = ("ds^2 = -{c2}dt^2 + \\dfrac{\\left(\\partial_r R - \\dfrac{R\\,\\partial_r E}{E}\\right)^2}"
                "{\\epsilon + f}dr^2 + \\dfrac{R^2}{E^2}\\left(dp^2 + dq^2\\right)")
        domains = ["p \\in (-\\infty, \\infty)", "q \\in (-\\infty, \\infty)"]
        probe = vm.Reader(coords, parameters, ())
        forms = SzekeresForms(probe)
        lead = [probe.parameters["epsilon"], probe.parameters["f"], probe.parameters["E"], probe.parameters["R"]]
        slope = sp.Derivative(probe.parameters["f"], probe.symbol["r"])
        extra = {"reduce": forms.reduce, "pretty": forms.pretty}
        collect = {"collect": lambda poly, printer: cp.collect_by(poly, [slope], printer)}
    else:
        coords, name = ["t", "r", "\\theta", "\\phi"], "Axisymmetric"
        parameters = ["R = R(t,r)", "f = f(r)", "S = S(r)"]
        line = ("ds^2 = -{c2}dt^2 + \\dfrac{\\left(\\partial_r R + \\dfrac{R\\,S'\\cos\\theta}{S}\\right)^2}"
                "{1 + f}dr^2 + R^2\\left(d\\theta - \\dfrac{S'\\sin\\theta}{S}dr\\right)^2 + R^2\\sin^2\\theta\\,d\\phi^2")
        domains = ["\\theta \\in [0, \\pi]", "\\phi \\in [0, 2\\pi)"]
        probe = vm.Reader(coords, parameters, ())
        lead = [probe.parameters["f"], probe.parameters["S"], probe.parameters["R"]]
        extra, collect = {"pretty": szekeres_polar(probe.symbol["\\theta"])}, {}
    return {
        "metric_id": "szekeres",
        "system": {"id": system_id, "name": name, "coords": coords,
                   "domains": ["t \\in (-\\infty, \\infty)", "r \\in [0, \\infty)"] + domains,
                   "parameters": parameters, "line_element": line.replace("{c2}", "c^2")},
        "chart_line_element": line.replace("{c2}", ""),
        "printer": {"lead": lead, "primed": ["f", "S", "P", "Q"], **collect},
        "check": szekeres_check,
        **extra,
    }


def szekeres_check(chart):
    """With (d_t R)^2 = 2M/R + f, M a function of r and a length, the Einstein tensor is that
    of dust at rest in the chart: its one component is G^t_t = -8 pi G rho/c^2, with the density
    Szekeres found, 2(M' - 3M E'/E)/(R^2(R' - R E'/E)) in the stereographic chart and the same
    with E'/E = -S' cos(theta)/S in the axisymmetric one."""
    t, r = chart.symbols[:2]
    R, f = chart.reader.parameters["R"], chart.reader.parameters["f"]
    M = sp.Function("M", real=True)(r)
    w = sp.Symbol("_w", positive=True)
    square = 2 * M / R + f
    rates = {(1, 0): w, (1, 1): sp.diff(square, r) / (2 * w), (2, 0): -M / R ** 2, (2, 1): sp.diff(-M / R ** 2, r)}

    def dust(value):
        """The value with every derivative of R along the time written by the dust's equation."""
        written = {}
        for d in value.atoms(sp.Derivative):
            count = dict(d.variable_count)
            if d.expr == R and count.get(t):
                written[d] = rates[count[t], count.get(r, 0)]
        return value.xreplace(written)
    if chart.coords_tex[2] == "p":
        ratio = sp.Derivative(chart.reader.parameters["E"], r) / chart.reader.parameters["E"]
    else:
        S = chart.reader.parameters["S"]
        ratio = -sp.Derivative(S, r) * sp.cos(chart.symbols[2]) / S
    density = 2 * (sp.diff(M, r) - 3 * M * ratio) / (R ** 2 * (sp.Derivative(R, r) - R * ratio))
    mixed = chart.geo.raise_indices(chart.geo.einstein_ll(), 2, (0,))
    reduce = chart.reduce or vm.norm
    for a in range(4):
        for b in range(4):
            value = dust(sp.sympify(mixed[a][b]))
            value = sp.together(value).subs(w ** 2, square)
            value = sp.expand(sp.numer(sp.together(value))).subs(w ** 2, square) / sp.denom(sp.together(value))
            if value.has(w) or reduce(value + (density if a == b == 0 else 0)) != 0:
                raise AssertionError(f"szekeres: with dust G^{chart.coords_tex[a]}_{chart.coords_tex[b]} "
                                     "is not that of dust at rest")


def szekeres_pullback():
    """The axisymmetric chart is the stereographic one with epsilon = 1 and P and Q constant,
    pulled back through p - P = S u cos(phi), q - Q = S u sin(phi), where u = cot(theta/2):
    both metrics are compared in the coordinates t, r, u and phi, in which each is rational."""
    stereo, polar = (cp.Chart(spec["system"]["coords"], spec["system"]["parameters"], spec["chart_line_element"])
                     for spec in (szekeres("stereographic"), szekeres("axisymmetric")))
    t, r, p, q = stereo.symbols
    T, rr, theta, phi = polar.symbols
    u = sp.Symbol("u", positive=True)
    S, P0, Q0 = polar.reader.parameters["S"], sp.Symbol("P_0", real=True), sp.Symbol("Q_0", real=True)
    held = {stereo.reader.parameters["P"]: P0, stereo.reader.parameters["Q"]: Q0,
            stereo.reader.parameters["epsilon"]: 1}
    at = {t: T, r: rr, p: P0 + S * u * sp.cos(phi), q: Q0 + S * u * sp.sin(phi)}
    source = stereo.reader.surface(stereo.geo.g).subs(held).doit().subs(at, simultaneous=True)
    images = [at[x] for x in stereo.symbols]
    jacobian = sp.Matrix(4, 4, lambda a, b: sp.diff(images[a], (T, rr, u, phi)[b]))
    pulled = jacobian.T * source * jacobian
    # theta as a function of u: cos(theta) = (u^2 - 1)/(u^2 + 1), sin(theta) = 2u/(u^2 + 1), dtheta = -2du/(u^2 + 1).
    angle = {sp.cos(theta): (u ** 2 - 1) / (u ** 2 + 1), sp.sin(theta): 2 * u / (u ** 2 + 1)}
    change = sp.diag(1, 1, -2 / (u ** 2 + 1), 1)
    target = change.T * polar.geo.g.subs(angle) * change
    for a in range(4):
        for b in range(a, 4):
            if vm.norm(pulled[a, b] - target[a, b]) != 0:
                raise AssertionError("szekeres: the pullback of the stereographic chart misses the axisymmetric "
                                     f"chart in slot {polar.coords_tex[a]}{polar.coords_tex[b]}")


def szekeres_charts():
    szekeres_pullback()
    return [szekeres("stereographic"), szekeres("axisymmetric")]


CHARTS["szekeres"] = szekeres_charts


# -- Teo's rotating wormhole ---------------------------------------------------------------

def teo_wormhole(system):
    """The example of E. Teo, Phys. Rev. D 58, 024014 (1998), his eq. (26) in the canonical
    metric (19): N = K = 1 + (4a cos theta)^2/r, b = 1 and omega = 2a/r^3 in units of the throat,
    here with the throat radius b_0 restored, N = K = 1 + 16a^2 b_0 cos^2(theta)/r, b = b_0 and
    omega = 2a b_0^2 c/r^3, and a = GJ/(c^3 b_0^2) a pure number. The first chart is his
    (t, r, theta, phi) on one side of the throat; the second is his eq. (16), in the proper
    radial distance l of his (15) and (28), with r a function of l the chart leaves free, so its
    values hold for every r(l). Every value is printed around rN = r + 16a^2 b_0 cos^2 theta,
    with cos^2 theta kept against sin^2 theta so that it survives as a factor.
    teo_wormhole_pullback checks the second chart against the first."""
    proper = system == "proper_radial"
    x = "l" if proper else "r"
    coords = ["t", x, "\\theta", "\\phi"]
    parameters = ["b_0", "a"] + (["r = r(l)"] if proper else [])
    probe = vm.Reader(coords, parameters, ())
    th, b0, a = probe.symbol["\\theta"], probe.parameters["b_0"], probe.parameters["a"]
    r = probe.parameters["r"] if proper else probe.symbol["r"]
    s, c = sp.sin(th), sp.cos(th)
    W = r + 16 * a ** 2 * b0 * c ** 2
    rN, rN_text = sp.Symbol("TEO_rN", positive=True), "r + 16a^2b_0\\cos^2\\theta"

    def in_cosine(e):
        return e.replace(lambda p: p.is_Pow and p.base == s and p.exp.is_Integer and p.exp > 1,
                         lambda p: s ** (int(p.exp) % 2) * (1 - c ** 2) ** (int(p.exp) // 2))

    def in_sine(e):
        return e.replace(lambda p: p.is_Pow and p.base == c and p.exp.is_Integer and p.exp > 1,
                         lambda p: c ** (int(p.exp) % 2) * (1 - s ** 2) ** (int(p.exp) // 2))

    def pretty(value):
        # Factored with sin^2 as 1 - cos^2, the angle N is written in, so that r + 16a^2 b_0 cos^2 theta
        # comes out as a factor, and (cos + 1)(cos - 1) as -sin^2; every other sum is then written
        # in whichever of sin^2 and cos^2 leaves it fewer terms, as r^4 - 4a^2 b_0^4 sin^2 theta,
        # the cosine where they tie.
        powers = {}
        for f in sp.Mul.make_args(sp.factor(in_cosine(sp.sympify(value)))):
            base, k = (f.base, f.exp) if f.is_Pow else (f, sp.Integer(1))
            powers[base] = powers.get(base, 0) + k
        k = min(powers.get(c + 1, 0), powers.get(c - 1, 0))
        if k > 0:
            powers[c + 1] -= k
            powers[c - 1] -= k
            powers[s] = powers.get(s, 0) + 2 * k
        out = sp.Integer(-1) ** max(k, 0)
        for base, k in powers.items():
            if base.is_Add and sp.expand(base - W) == 0:
                base = rN
            elif base.is_Add and sp.expand(base + W) == 0:
                base = -rN
            elif base.is_Add:
                sine = sp.expand(in_sine(base))
                if len(sp.Add.make_args(sine)) < len(sp.Add.make_args(base)):
                    base = sine
            out *= base ** k
        return out

    lapse = "\\left(1 + \\dfrac{16a^2b_0\\cos^2\\theta}{r}\\right)^2"

    def line(cdt):
        radial = "dl^2" if proper else "\\dfrac{dr^2}{1 - \\dfrac{b_0}{r}}"
        return ("ds^2 = -" + lapse + ("c^2" if cdt else "") + "dt^2 + " + radial + " + \\left(" + rN_text
                + "\\right)^2\\left(d\\theta^2 + \\sin^2\\theta\\left(d\\phi - \\dfrac{2ab_0^2}{r^3}" + cdt
                + "dt\\right)^2\\right)")

    domains = ["t \\in (-\\infty, \\infty)", "l \\in (-\\infty, \\infty)" if proper else "r \\in [b_0, \\infty)",
               "\\theta \\in [0, \\pi]", "\\phi \\in [0, 2\\pi)",
               ("l = 0" if proper else "r = b_0") + " \\;\\text{(throat)}"]
    return {
        "metric_id": "teo_wormhole",
        "system": {"id": system, "name": "Proper Radial Distance" if proper else "Spherical", "coords": coords,
                   "domains": domains, "parameters": parameters, "line_element": line("c\\,")},
        "chart_line_element": line(""),
        "printer": {"lead": [a, r, b0, c, s], "named": {rN: rN_text}},
        "pretty": pretty,
        "bracketed": pretty,
        # The metric as the line element writes it: omega^2 r^2 N^2 sin^2 theta taken from N^2 in g_tt.
        "components": {
            "metric_components": {
                ("t", "t"): "-" + lapse + "\\left(1 - \\dfrac{4a^2b_0^4\\sin^2\\theta}{r^4}\\right)",
                **({} if proper else {("r", "r"): "\\left(1 - \\dfrac{b_0}{r}\\right)^{-1}"})},
            "inverse_metric_components": {} if proper else {("r", "r"): "1 - \\dfrac{b_0}{r}"}},
    }


def teo_wormhole_pullback(chart):
    """The proper distance chart with dr/dl = sqrt(1 - b_0/r) is the spherical one: J^T g J with
    J = diag(1, dr/dl, 1, 1) and g the spherical metric at r = r(l), in every slot."""
    spec = teo_wormhole("spherical")
    source = cp.Chart(spec["system"]["coords"], spec["system"]["parameters"], spec["chart_line_element"])
    r = chart.reader.parameters["r"]
    at = dict(zip(source.symbols, [chart.symbols[0], r, chart.symbols[2], chart.symbols[3]]))
    at.update({source.reader.parameters[k]: chart.reader.parameters[k] for k in ("b_0", "a")})
    b0 = chart.reader.parameters["b_0"]
    jacobian = sp.diag(1, sp.sqrt(1 - b0 / r), 1, 1)
    pulled = jacobian.T * source.geo.g.subs(at, simultaneous=True) * jacobian
    for i in range(4):
        for j in range(i, 4):
            if sp.simplify(pulled[i, j] - chart.geo.g[i, j]) != 0:
                raise AssertionError(f"teo_wormhole: the spherical chart pulled back misses the proper distance "
                                     f"chart in slot {chart.coords_tex[i]}{chart.coords_tex[j]}")


def teo_wormhole_charts():
    spherical, proper = teo_wormhole("spherical"), teo_wormhole("proper_radial")
    proper["check"] = teo_wormhole_pullback
    return [spherical, proper]


CHARTS["teo_wormhole"] = teo_wormhole_charts


# -- Kaluza-Klein monopole -------------------------------------------------------------

KK_CHARTS = ["gross_perry", "hopf", "taub_nut"]


def kaluza_klein_monopole(system_id):
    """Sorkin's and Gross and Perry's monopole, Euclidean Taub-NUT space with -c^2dt^2 added, in
    three charts of its five dimensions. Gross and Perry's keeps the fifth coordinate x_5, of
    period 16 pi m, with the potential 4m(1 - cos theta) d phi, regular on theta = 0. The Hopf
    chart trades x_5 for the Euler angle psi = -phi - x_5/4m of period 4 pi, which leaves
    16m^2 (d psi + cos theta d phi)^2/V on every 3-sphere of constant t and r. The Taub-NUT
    chart moves the radius to rho = r + 2m, Hawking's form of the self-dual instanton with NUT
    parameter 2m. Each chart after the first is checked, slot by slot, to be the first pulled
    back through its map. The Kretschmann scalar is 384 m^2/(r + 4m)^6, finite
    at the nut r = 0. kaluza_klein_monopole.md derives each chart."""
    V = "\\left(1 + \\dfrac{4m}{r}\\right)"
    bareV = "1 + \\dfrac{4m}{r}"
    reals = "(-\\infty, \\infty)"
    time = "t \\in " + reals
    nut = "r = 0 \\;\\text{(the nut, a regular point)}"
    kretschmann = "\\dfrac{384m^2}{\\left(r + 4m\\right)^6}"
    parameters = ["m"]
    components = {}
    if system_id == "gross_perry":
        coords, name = ["t", "r", "\\theta", "\\phi", "x_5"], "Gross-Perry"
        domains = [time, "r \\in [0, \\infty)", "\\theta \\in [0, \\pi]", "\\phi \\in [0, 2\\pi)",
                   "x_5 \\in [0, 16\\pi m)", nut,
                   "\\theta = \\pi \\;\\text{(the Dirac string, a coordinate singularity)}"]
        space = (V + "\\left(dr^2 + r^2d\\theta^2 + r^2\\sin^2\\theta\\,d\\phi^2\\right) + " + V
                 + "^{-1}\\left(dx_5 + 4m\\left(1 - \\cos\\theta\\right)d\\phi\\right)^2")
        components = {"metric_components": {("r", "r"): bareV, ("x_5", "x_5"): V + "^{-1}"}}
    elif system_id == "hopf":
        coords, name = ["t", "r", "\\theta", "\\phi", "\\psi"], "Hopf Angle"
        domains = [time, "r \\in [0, \\infty)", "\\theta \\in [0, \\pi]", "\\phi \\in [0, 2\\pi)",
                   "\\psi \\in [0, 4\\pi)", nut]
        space = (V + "\\left(dr^2 + r^2d\\theta^2 + r^2\\sin^2\\theta\\,d\\phi^2\\right) + 16m^2" + V
                 + "^{-1}\\left(d\\psi + \\cos\\theta\\,d\\phi\\right)^2")
        components = {"metric_components": {("r", "r"): bareV, ("\\psi", "\\psi"): "16m^2" + V + "^{-1}"}}
    else:
        coords, name = ["t", "\\rho", "\\theta", "\\phi", "\\psi"], "Taub-NUT Radius"
        domains = [time, "\\rho \\in [2m, \\infty)", "\\theta \\in [0, \\pi]", "\\phi \\in [0, 2\\pi)",
                   "\\psi \\in [0, 4\\pi)", "\\rho = 2m \\;\\text{(the nut, a regular point)}"]
        space = ("\\dfrac{\\rho + 2m}{\\rho - 2m}d\\rho^2 + \\left(\\rho^2 - 4m^2\\right)\\left(d\\theta^2"
                 " + \\sin^2\\theta\\,d\\phi^2\\right) + \\dfrac{16m^2\\left(\\rho - 2m\\right)}{\\rho + 2m}"
                 "\\left(d\\psi + \\cos\\theta\\,d\\phi\\right)^2")
        kretschmann = "\\dfrac{384m^2}{\\left(\\rho + 2m\\right)^6}"
    probe = vm.Reader(coords, parameters, ())
    radius = probe.symbol[coords[1]]
    spec = {
        "metric_id": "kaluza_klein_monopole",
        "system": {"id": system_id, "name": name, "coords": coords, "domains": domains,
                   "parameters": parameters, "line_element": "ds^2 = -c^2dt^2 + " + space},
        "chart_line_element": "ds^2 = -dt^2 + " + space,
        "printer": {"lead": [radius, probe.parameters["m"]], "flip": False},
        "components": components,
        "kretschmann": kretschmann,
    }
    if system_id != "gross_perry":
        spec["check"] = lambda chart, s=system_id: kaluza_klein_pullback(chart, s)
    return spec


def kaluza_klein_pullback(chart, system_id):
    """J^T g J, with g Gross and Perry's chart and J the Jacobian of the map from `chart` into
    it, against the metric of `chart`, in every slot: x_5 = -4m(psi + phi) for the Hopf angle,
    and r = rho - 2m besides for the Taub-NUT radius."""
    spec = kaluza_klein_monopole("gross_perry")
    source = cp.Chart(spec["system"]["coords"], spec["system"]["parameters"], spec["chart_line_element"])
    m = chart.reader.parameters["m"]
    t, a, b, c, w = chart.symbols
    if system_id == "hopf":
        image = [t, a, b, c, -4 * m * (w + c)]
    else:
        image = [t, a - 2 * m, b, c, -4 * m * (w + c)]
    J = sp.Matrix(5, 5, lambda i, j: sp.diff(image[i], chart.symbols[j]))
    at = dict(zip(source.symbols, image))
    at[source.reader.parameters["m"]] = m
    pulled = J.T * source.geo.g.subs(at, simultaneous=True) * J
    own = chart.geo.g
    for i in range(5):
        for j in range(i, 5):
            if sp.simplify(pulled[i, j] - own[i, j]) != 0:
                raise AssertionError(f"kaluza_klein_monopole: Gross and Perry's chart pulled back misses the "
                                     f"{system_id} chart in slot {chart.coords_tex[i]}{chart.coords_tex[j]}")


CHARTS["kaluza_klein_monopole"] = [lambda s=s: kaluza_klein_monopole(s) for s in KK_CHARTS]
# -- Bell-Szekeres -----------------------------------------------------------------------

BS_CHARTS = ("double_null", "time_space", "regular", "global", "kruskal_szekeres", "bertotti_robinson")


def bell_szekeres(system):
    """The Bell-Szekeres spacetime where both electromagnetic waves have passed, in six charts,
    every one with the signature flipped from its source's. The double null chart is Bell and
    Szekeres's, Gen. Rel. Grav. 5, 275 (1974), eq. (15.7) of J. B. Griffiths, Colliding Plane
    Waves in General Relativity (Oxford, 1991). The chart of xi = au + bv and eta = bv - au is
    eq. (9) of A. Feinstein and M. A. Perez Sebastian, Class. Quantum Grav. 12, 2723 (1995). The
    regular chart is the one Bell and Szekeres removed the coordinate singularity with, Griffiths's
    eqs. (15.11) and (15.12). The global chart is C. J. S. Clarke and S. A. Hayward's, Class.
    Quantum Grav. 6, 615 (1989), Griffiths's eqs. (15.13) to (15.15). The Kruskal-Szekeres chart
    is Feinstein and Perez Sebastian's eqs. (10) to (13), after M. Dorca and E. Verdaguer, and the
    conformally flat chart of the Bertotti-Robinson universe is Griffiths's eq. (15.9).
    bell_szekeres_check pulls each of the last five back onto the first and checks the
    Einstein-Maxwell equations in every chart, and bell_szekeres.md beside this file is the
    derivation."""
    reals = "(-\\infty, \\infty)"
    charts = {
        "double_null": {
            "name": "Double Null", "coords": ["u", "v", "x", "y"],
            "domains": ["u \\in [0, \\pi/(2a))", "v \\in [0, (\\pi/2 - au)/b)", "x \\in " + reals, "y \\in " + reals,
                        "au + bv = \\pi/2 \\;\\text{(Killing-Cauchy horizon)}"],
            "line": "ds^2 = -2\\,du\\,dv + \\cos^2\\left(au - bv\\right)dx^2 + \\cos^2\\left(au + bv\\right)dy^2"},
        "time_space": {
            "name": "Time and Space", "coords": ["\\xi", "\\eta", "x", "y"],
            "domains": ["\\xi \\in [0, \\pi/2)", "\\eta \\in [-\\xi, \\xi]", "x \\in " + reals, "y \\in " + reals,
                        "\\xi = \\pi/2 \\;\\text{(Killing-Cauchy horizon)}"],
            "line": "ds^2 = \\dfrac{-d\\xi^2 + d\\eta^2}{2ab} + \\cos^2\\eta\\,dx^2 + \\cos^2\\xi\\,dy^2"},
        "regular": {
            "name": "Regular", "coords": ["T", "Z", "X", "Y"],
            "domains": ["T \\in \\left[-\\sqrt{1 + Z^2},\\, -|Z|\\right)", "Z \\in " + reals, "X \\in (-1, 1)",
                        "Y \\in \\left(-\\sqrt{1 - X^2},\\, \\sqrt{1 - X^2}\\right)",
                        "T^2 - Z^2 \\le X^2 + Y^2 \\;\\text{(behind both wave fronts)}",
                        "T^2 = Z^2 \\;\\text{(Killing-Cauchy horizon)}"],
            "line": "ds^2 = \\dfrac{-\\left(1 + Z^2\\right)dT^2 + 2TZ\\,dT\\,dZ + \\left(1 - T^2\\right)dZ^2}"
                    "{2ab\\left(1 - T^2 + Z^2\\right)}"
                    " + \\dfrac{\\left(1 - Y^2\\right)dX^2 + 2XY\\,dX\\,dY + \\left(1 - X^2\\right)dY^2}"
                    "{2ab\\left(1 - X^2 - Y^2\\right)}"},
        "global": {
            "name": "Global", "coords": ["\\chi", "\\rho", "\\theta", "\\phi"],
            "domains": ["\\chi \\in \\left[-\\pi/2,\\, -\\arcsin|\\tanh\\rho|\\right)", "\\rho \\in " + reals,
                        "\\theta \\in \\left[\\arccos(\\cos\\chi\\cosh\\rho),\\, \\pi - \\arccos(\\cos\\chi\\cosh\\rho)\\right]",
                        "\\phi \\in " + reals,
                        "\\cos\\chi\\cosh\\rho = 1 \\;\\text{(Killing-Cauchy horizon)}"],
            "line": "ds^2 = \\dfrac{1}{2ab}\\left(-\\cosh^2\\rho\\,d\\chi^2 + d\\rho^2 + d\\theta^2 + \\sin^2\\theta\\,d\\phi^2\\right)"},
        "kruskal_szekeres": {
            "name": "Kruskal-Szekeres", "coords": ["U", "V", "\\eta", "x"],
            "domains": ["U \\in (-\\infty, 0)", "V \\in \\left[\\dfrac{1}{2abU},\\, 0\\right)",
                        "\\eta \\in (-\\pi/2, \\pi/2)", "x \\in " + reals,
                        "2abUV \\le \\dfrac{1 - \\sin|\\eta|}{1 + \\sin|\\eta|} \\;\\text{(behind both wave fronts)}",
                        "UV = 0 \\;\\text{(Killing-Cauchy horizon)}"],
            "line": "ds^2 = -\\dfrac{4\\,dU\\,dV}{\\left(1 + 2abUV\\right)^2} + \\dfrac{d\\eta^2}{2ab} + \\cos^2\\eta\\,dx^2"},
        "bertotti_robinson": {
            "name": "Bertotti-Robinson", "coords": ["t", "r", "\\theta", "\\phi"],
            "domains": ["t \\in [0, r)", "r \\in (0, \\infty)",
                        "\\theta \\in \\left[\\pi/2 - \\arcsin(t/r),\\, \\pi/2 + \\arcsin(t/r)\\right]",
                        "\\phi \\in " + reals, "t = r \\;\\text{(Killing-Cauchy horizon)}"],
            "line": "ds^2 = \\dfrac{1}{2ab\\,r^2}\\left(-dt^2 + dr^2 + r^2d\\theta^2 + r^2\\sin^2\\theta\\,d\\phi^2\\right)"},
    }
    chart = charts[system]
    probe = vm.Reader(chart["coords"], ["a", "b"], ())
    a, b = probe.parameters["a"], probe.parameters["b"]
    if system == "double_null":
        printing = {"printer": {"lead": [a, b]},
                    "pretty": bell_szekeres_sum_and_difference(a * probe.symbol["u"], b * probe.symbol["v"])}
    elif system == "regular":
        # Each sum is printed whole and in the order the line element writes it.
        T, Z, X, Y = (probe.symbol[name] for name in chart["coords"])
        named = [(sp.Symbol("BSh"), 1 - T ** 2 + Z ** 2, "1 - T^2 + Z^2"),
                 (sp.Symbol("BSs"), 1 - X ** 2 - Y ** 2, "1 - X^2 - Y^2"),
                 (sp.Symbol("BSz"), 1 + Z ** 2, "1 + Z^2"), (sp.Symbol("BSt"), 1 - T ** 2, "1 - T^2"),
                 (sp.Symbol("BSx"), 1 - X ** 2, "1 - X^2"), (sp.Symbol("BSy"), 1 - Y ** 2, "1 - Y^2")]
        merges = [(1 + w, 1 - w, placeholder) for w, (placeholder, _, _) in zip((T, X, Y), named[3:])]
        printing = {"printer": {"lead": [a, b], "named": {placeholder: text for placeholder, _, text in named}},
                    "pretty": named_factors(named, merges)}
    else:
        # (1 + sin)(1 - sin) of xi, eta or theta reads as the cosine squared the line element writes.
        angles = [name for name in ("\\xi", "\\eta", "\\theta") if name in chart["coords"]]
        pairs = [(sp.sin(probe.symbol[name]), sp.cos(probe.symbol[name]),
                  sp.Symbol(f"BSp{i}", positive=True), sp.Symbol(f"BSm{i}", positive=True))
                 for i, name in enumerate(angles)]
        named = {}
        for (_, _, plus, minus), name in zip(pairs, angles):
            named[plus], named[minus] = f"1 + \\sin{name}", f"1 - \\sin{name}"
        # The Kruskal-Szekeres chart's 1 + 2abUV keeps the order the line element writes it in.
        rising = [probe.symbol["U"]] if system == "kruskal_szekeres" else []
        hyperbolic = nariai_hyperbolic(probe.symbol["\\rho"]) if system == "global" else (lambda value: value)
        printing = {"printer": {"lead": [a, b], "named": named, "rising": rising, "flip": False},
                    "pretty": lambda value: trig_pairs(hyperbolic(value), pairs)}
    return {
        "metric_id": "bell_szekeres",
        "system": {"id": system, "name": chart["name"], "coords": chart["coords"], "domains": chart["domains"],
                   "parameters": ["a", "b"], "line_element": chart["line"]},
        "chart_line_element": chart["line"],
        "check": lambda c: bell_szekeres_check(c, system),
        **printing,
    }


def bell_szekeres_sum_and_difference(p, q):
    """A pretty printer for the double null chart, whose every value is a function of the sum
    p + q = au + bv and the difference p - q = au - bv. The checker's canonical form spreads
    each over the sines and cosines of p and of q; this writes it back with
        cos^2 p = (1 + CP CM - SP SM)/2,    sin p cos p = (SP CM + CP SM)/2,
        cos^2 q = (1 + CP CM + SP SM)/2,    sin q cos q = (SP CM - CP SM)/2,
    where SP, CP are the sine and cosine of the sum and SM, CM those of the difference, which
    reaches every monomial of even degree in p's pair and in q's. A numerator and a denominator
    of odd degree are both multiplied by cos p or cos q first. Of the four ways of writing the
    squares, sines or cosines for each angle, the one that factors shortest is printed."""
    s1, c1, s2, c2 = sp.sin(p), sp.cos(p), sp.sin(q), sp.cos(q)
    SP, CP, SM, CM = sp.symbols("BSsp BScp BSsm BScm")
    square = {0: (1 + CP * CM - SP * SM) / 2, 1: (1 + CP * CM + SP * SM) / 2}
    cross = {0: (SP * CM + CP * SM) / 2, 1: (SP * CM - CP * SM) / 2}
    back = {SP: sp.sin(p + q), CP: sp.cos(p + q), SM: sp.sin(p - q), CM: sp.cos(p - q)}

    def odd(poly, pair):
        degrees = {sum(m[k] for k in pair) % 2 for m in poly.monoms()}
        if len(degrees) > 1:
            raise AssertionError(f"bell_szekeres: {poly.as_expr()} mixes parities in one angle")
        return degrees == {1}

    def carried(expr):
        poly = sp.Poly(sp.expand(expr), s1, c1, s2, c2)
        total = sp.Integer(0)
        for (i, j, k, m), coefficient in poly.terms():
            term = coefficient
            for which, (sine, cosine) in enumerate(((i, j), (k, m))):
                term *= (1 - square[which]) ** (sine // 2) * cross[which] ** (sine % 2) \
                    * square[which] ** ((cosine - sine % 2) // 2)
            total += term
        return sp.expand(total)

    def reduced(expr, use_sines):
        for (sine, cosine), use_sine in zip(((SP, CP), (SM, CM)), use_sines):
            keep, drop = (sine, cosine) if use_sine else (cosine, sine)
            expr = sp.expand(sum(k * drop ** (d % 2) * (1 - keep ** 2) ** (d // 2) * rest
                                 for (d,), k in sp.Poly(expr, drop).terms()
                                 for rest in [1]))
        return expr

    def pretty(value):
        value = sp.sympify(value)
        if not value.has(s1, c1, s2, c2):
            return sp.factor(value)
        numerator, denominator = sp.fraction(sp.together(value))
        for pair, cosine in (((0, 1), c1), ((2, 3), c2)):
            if odd(sp.Poly(sp.expand(numerator), s1, c1, s2, c2), pair):
                numerator, denominator = numerator * cosine, denominator * cosine
        numerator, denominator = carried(numerator), carried(denominator)
        best = None
        for use_sines in itertools.product((False, True), repeat=2):
            candidate = sp.factor(reduced(numerator, use_sines)) / sp.factor(reduced(denominator, use_sines))
            candidate = sp.factor(sp.cancel(candidate))
            if best is None or sp.count_ops(candidate) < sp.count_ops(best):
                best = candidate
        return best.subs(back, simultaneous=True)

    return pretty


def bell_szekeres_field(chart, system):
    """The field where both waves have passed, F = sqrt(2ab) times the area form of the anti-de
    Sitter factor, uniform, in units where G_mu nu = 2(F_mu a F_nu^a - g_mu nu F^2/4). In the
    double null chart and the chart of xi and eta it is dA with A = sin(au + bv) dy, and in the
    others, whose first two coordinates span that factor, it is read off the metric."""
    x, g = chart.symbols, chart.geo.g
    a, b = chart.reader.parameters["a"], chart.reader.parameters["b"]
    if system in ("double_null", "time_space"):
        A = [0, 0, 0, sp.sin(a * x[0] + b * x[1]) if system == "double_null" else sp.sin(x[0])]
        return sp.Matrix(4, 4, lambda i, j: sp.diff(A[j], x[i]) - sp.diff(A[i], x[j]))
    area = sp.powdenest(sp.sqrt(sp.factor(-g[:2, :2].det())), force=True).replace(sp.Abs, lambda e: e)
    F = sp.zeros(4, 4)
    F[0, 1], F[1, 0] = sp.sqrt(2 * a * b) * area, -sp.sqrt(2 * a * b) * area
    return F


def bell_szekeres_maps(u, v, x, y, a, b):
    """Each chart's coordinates as functions of the double null chart's."""
    q = 1 / sp.sqrt(2 * a * b)
    xi, eta = a * u + b * v, b * v - a * u
    fall = q * sp.cos(xi) / (1 + sp.sin(xi))
    c = sp.sqrt(2 * a * b)
    T, Z = -sp.cos(xi) * sp.cosh(c * y), sp.cos(xi) * sp.sinh(c * y)
    return {"time_space": [xi, eta, x, y],
            "regular": [T, Z, sp.cos(eta) * sp.cos(c * x), sp.cos(eta) * sp.sin(c * x)],
            "global": [sp.asin(T / sp.sqrt(1 + Z ** 2)), sp.asinh(Z), sp.pi / 2 - eta, c * x],
            "kruskal_szekeres": [-fall * sp.exp(y / q), -fall * sp.exp(-y / q), eta, x],
            "bertotti_robinson": [sp.exp(y / q) * sp.tan(xi), sp.exp(y / q) / sp.cos(xi), sp.pi / 2 - eta, x / q]}


def bell_szekeres_check(chart, system, points=12):
    """The Einstein tensor is the stress of bell_szekeres_field in every slot and that field
    solves Maxwell's equations; the Weyl tensor vanishes; and J^T g J, with J the Jacobian of
    the chart's coordinates as functions of u, v, x and y, is the double null metric at random
    points where both waves have passed, to forty digits."""
    x = chart.symbols
    a, b = chart.reader.parameters["a"], chart.reader.parameters["b"]
    g, ginv = chart.geo.g, chart.geo.ginv
    F = bell_szekeres_field(chart, system)
    for i, j, k in itertools.combinations(range(4), 3):
        if sp.simplify(sp.diff(F[i, j], x[k]) + sp.diff(F[j, k], x[i]) + sp.diff(F[k, i], x[j])) != 0:
            raise AssertionError(f"bell_szekeres: the field is not closed in the {system} chart")
    Fup = ginv * F * ginv
    F2 = sum(F[i, j] * Fup[i, j] for i in range(4) for j in range(4))
    G = chart.geo.einstein_ll()
    for i in range(4):
        for j in range(i, 4):
            T = 2 * (sum(F[i, c] * F[j, d] * ginv[c, d] for c in range(4) for d in range(4)) - g[i, j] * F2 / 4)
            if sp.simplify(vm._at(G, (i, j)) - T) != 0:
                raise AssertionError(f"bell_szekeres: the Einstein tensor misses Maxwell's stress in slot "
                                     f"{chart.coords_tex[i]}{chart.coords_tex[j]} of the {system} chart")
    # sqrt(-g) is a product of factors that are each positive on the chart.
    root = sp.powdenest(sp.sqrt(sp.factor(-g.det())), force=True).replace(sp.Abs, lambda e: e)
    for j in range(4):
        if sp.simplify(sum(sp.diff(root * Fup[i, j], x[i]) for i in range(4)) / root) != 0:
            raise AssertionError(f"bell_szekeres: Maxwell's equations fail along {chart.coords_tex[j]} "
                                 f"in the {system} chart")
    weyl = chart.geo.weyl_llll()
    if any(vm.norm(vm._at(weyl, index)) != 0 for index in vm._indices(4, 4)):
        raise AssertionError(f"bell_szekeres: the {system} chart is not conformally flat")
    if system == "double_null":
        return
    null = cp.Chart(["u", "v", "x", "y"], ["a", "b"], bell_szekeres("double_null")["chart_line_element"])
    u, v, xx, y = null.symbols
    na, nb = null.reader.parameters["a"], null.reader.parameters["b"]
    image = bell_szekeres_maps(u, v, xx, y, na, nb)[system]
    J = sp.Matrix(4, 4, lambda i, j: sp.diff(image[i], null.symbols[j]))
    at = dict(zip(x, image))
    at.update({a: na, b: nb})
    difference = J.T * g.subs(at, simultaneous=True) * J - null.geo.g
    rng = __import__("random").Random(5)
    for _ in range(points):
        p = {na: sp.Rational(rng.randint(5, 30), 10), nb: sp.Rational(rng.randint(5, 30), 10),
             xx: sp.Rational(rng.randint(-20, 20), 10), y: sp.Rational(rng.randint(-20, 20), 10)}
        total, share = sp.Rational(rng.randint(1, 150), 100), sp.Rational(rng.randint(1, 99), 100)
        p[u], p[v] = total * share / p[na], total * (1 - share) / p[nb]
        worst = max(abs(sp.N(difference[i, j].subs(p), 40)) for i in range(4) for j in range(4))
        if worst > sp.Float("1e-30"):
            raise AssertionError(f"bell_szekeres: the {system} chart pulled back misses the double null chart's "
                                 f"metric by {worst} at {p}")


CHARTS["bell_szekeres"] = [lambda s=s: bell_szekeres(s) for s in BS_CHARTS]
# -- Kinnersley's photon rocket ----------------------------------------------------------

ROCKET_P = ("\\sqrt{1 + U_1^2 + U_2^2 + U_3^2} - U_1\\cos\\theta - U_2\\sin\\theta\\cos\\phi"
            " - U_3\\sin\\theta\\sin\\phi")
ROCKET_ORDERS = 4


class RocketForms:
    """The function p = U_0 - U_1 cos(theta) - U_2 sin(theta) cos(phi) - U_3 sin(theta) sin(phi)
    of the rocket's Robinson-Trautman chart, with U_0 = sqrt(1 + U_1^2 + U_2^2 + U_3^2), held
    as a function of u, theta and phi while the tensors are built, and the way a value in it
    is written.

    p and each of its derivatives along u is a constant on the sphere less a function linear in
    the direction of the ray, q = A - B.n, and every such q obeys d_theta^2 q = A - q,
    d_theta d_phi q = cot(theta) d_phi q and d_phi^2 q = (A - q) sin^2(theta)
    - sin(theta) cos(theta) d_theta q. The constant of p is U_0, and since the four-velocity is
    a unit vector, 2 U_0 p = p^2 + 1 + (d_theta p)^2 + (d_phi p)^2/sin^2(theta), so the constant
    of the k-th derivative of p along u is the k-th derivative of that. `reduce` writes those
    in, which leaves p, d_theta p, d_phi p and their derivatives along u, with no relation
    among them, so a value that vanishes for the rocket's p is exactly zero."""

    def __init__(self, reader):
        self.u, self.theta, self.phi = (reader.symbol[name] for name in ("u", "\\theta", "\\phi"))
        self.p = reader.parameters["p"]
        a, b, c = (sp.Function(name)(self.u) for name in ("_a", "_b", "_c"))
        which = {a: None, b: self.theta, c: self.phi}
        constant = (a ** 2 + 1 + b ** 2 + c ** 2 / sp.sin(self.theta) ** 2) / (2 * a)
        self.constant = {}     # the constant of the k-th derivative of p along u, in the generators
        for k in range(ROCKET_ORDERS):
            value = sp.diff(constant, self.u, k)
            value = value.xreplace({d: self.generator(which[d.expr], d.variable_count[0][1])
                                    for d in value.atoms(sp.Derivative)})
            self.constant[k] = value.xreplace({f: self.generator(x, 0) for f, x in which.items()})

    def generator(self, x, k):
        """d_u^k of p, or of d_x p where x is theta or phi."""
        variables = ([x] if x is not None else []) + ([(self.u, k)] if k else [])
        return sp.Derivative(self.p, *variables) if variables else self.p

    def angular(self, i, j, k):
        """d_theta^i d_phi^j d_u^k p with i + j >= 2, in the generators."""
        theta, phi = self.theta, self.phi
        q, A = sp.Function("_q")(theta, phi), sp.Symbol("_A")
        second = {(2, 0): A - q, (1, 1): sp.cos(theta) / sp.sin(theta) * sp.Derivative(q, phi),
                  (0, 2): (A - q) * sp.sin(theta) ** 2 - sp.sin(theta) * sp.cos(theta) * sp.Derivative(q, theta)}

        def lowered(value):
            written = {}
            for d in value.atoms(sp.Derivative):
                count = dict(d.variable_count)
                m, n = count.get(theta, 0), count.get(phi, 0)
                if m + n < 2:
                    continue
                first = (2, 0) if m >= 2 else (1, 1) if m == 1 else (0, 2)
                rest = [theta] * (m - first[0]) + [phi] * (n - first[1])
                written[d] = lowered(sp.diff(second[first], *rest) if rest else second[first])
            return value.xreplace(written)
        value = lowered(sp.Derivative(q, *([theta] * i + [phi] * j)))
        return value.xreplace({sp.Derivative(q, theta): self.generator(theta, k),
                               sp.Derivative(q, phi): self.generator(phi, k)}
                              ).xreplace({q: self.generator(None, k), A: self.constant[k]})

    def reduce(self, value):
        value = sp.sympify(value)
        if not value.has(self.p):
            return vm.norm(value)
        written = {}
        for d in value.atoms(sp.Derivative):
            if d.expr != self.p:
                continue
            count = {x: 0 for x in (self.u, self.theta, self.phi)}
            for x, n in d.variable_count:
                count[x] += n
            if count[self.u] >= ROCKET_ORDERS:
                raise AssertionError(f"photon_rocket: {d} is of a higher order along u than ROCKET_ORDERS allows")
            if count[self.theta] + count[self.phi] >= 2:
                written[d] = self.angular(count[self.theta], count[self.phi], count[self.u])
        # A generator is held as a symbol while the others are written in, so that none is
        # replaced inside another.
        held = {}
        for expression in (value, *written.values()):
            for d in expression.atoms(sp.Derivative):
                if d not in written:
                    held.setdefault(d, sp.Dummy())
        value = value.xreplace({d: image.xreplace(held) for d, image in written.items()}).xreplace(held)
        return vm.norm(value.xreplace({dummy: d for d, dummy in held.items()}))


def photon_rocket(system_id):
    """Kinnersley's photon rocket in the two charts its literature uses: the chart of a rocket
    that moves along a line, whose angle theta is measured in the rocket's rest frame from the
    direction opposite to its acceleration alpha (von der Gonna and Kramer's (58), Podolsky's
    (12) of 2008 without the cosmological constant), and the Robinson-Trautman chart of a
    rocket in arbitrary motion, whose angles are those of the background inertial frame
    (Podolsky's (52) and (53) of 2011 in four dimensions). m, alpha and the four-velocity are
    left free in every tensor; photon_rocket_check confirms before anything is written that
    the matter is null dust along d_r, and photon_rocket.md is the derivation."""
    coords = ["u", "r", "\\theta", "\\phi"]
    extra = {}
    if system_id == "rectilinear":
        name, parameters = "Rectilinear", ["m = m(u)", "\\alpha = \\alpha(u)"]
        line = ("ds^2 = -\\left(1 - \\dfrac{2m}{r} - 2\\alpha r\\cos\\theta\\right){c2}du^2 - 2{c}du\\,dr"
                " + r^2\\left(d\\theta + \\alpha\\sin\\theta\\,{c}du\\right)^2 + r^2\\sin^2\\theta\\,d\\phi^2")
        probe = vm.Reader(coords, parameters, ())
        lead = [probe.parameters["alpha"], probe.parameters["m"]]
    else:
        name = "Robinson-Trautman"
        parameters = ["m = m(u)", "U_1 = U_1(u)", "U_2 = U_2(u)", "U_3 = U_3(u)", "p = " + ROCKET_P]
        line = ("ds^2 = -\\left(1 - \\dfrac{2m}{r} - \\dfrac{2r\\,\\partial_u p}{p}\\right){c2}du^2 - 2{c}du\\,dr"
                " + \\dfrac{r^2}{p^2}\\left(d\\theta^2 + \\sin^2\\theta\\,d\\phi^2\\right)")
        probe = vm.Reader(coords, parameters, ())
        lead = [probe.parameters["p"], probe.parameters["m"]]
        extra = {"reduce": RocketForms(probe).reduce}
    return {
        "metric_id": "photon_rocket",
        "system": {"id": system_id, "name": name, "coords": coords,
                   "domains": ["u \\in (-\\infty, \\infty)", "r \\in (0, \\infty)", "\\theta \\in [0, \\pi]",
                               "\\phi \\in [0, 2\\pi)"],
                   "parameters": parameters,
                   "line_element": line.replace("{c2}", "c^2").replace("{c}", "c\\,")},
        "chart_line_element": line.replace("{c2}", "").replace("{c}", ""),
        "printer": {"lead": lead},
        "check": photon_rocket_check,
        **extra,
    }


def photon_rocket_check(chart):
    """The Ricci tensor has the one component R_uu = 2(3m d_u ln p - d_u m)/r^2, with
    d_u ln p = alpha cos(theta) in the rectilinear chart, which is null dust streaming along
    d_r; the Kretschmann scalar is 48m^2/r^6; and with m = 0 the Riemann tensor vanishes, so the
    chart is one of Minkowski space."""
    u, r, theta, _ = chart.symbols
    m = chart.reader.parameters["m"]
    if "p" in chart.reader.parameters:
        p = chart.reader.parameters["p"]
        rate = sp.Derivative(p, u) / p
    else:
        rate = chart.reader.parameters["alpha"] * sp.cos(theta)
    reduce = chart.reduce or vm.norm
    dust = 2 * (3 * m * rate - sp.Derivative(m, u)) / r ** 2
    ricci = chart.geo.ricci_ll()
    for i in range(4):
        for j in range(i, 4):
            if reduce(sp.sympify(ricci[i][j]) - (dust if i == j == 0 else 0)) != 0:
                raise AssertionError(f"photon_rocket: R_{chart.coords_tex[i]}{chart.coords_tex[j]} is not that of "
                                     "null dust along d_r")
    if reduce(sp.sympify(chart.geo.kretschmann()) - 48 * m ** 2 / r ** 6) != 0:
        raise AssertionError("photon_rocket: the Kretschmann scalar is not 48m^2/r^6")
    riemann = chart.geo.riemann_llll()
    for index in itertools.product(range(4), repeat=4):
        if reduce(sp.sympify(vm._at(riemann, index)).subs(m, 0).doit()) != 0:
            raise AssertionError("photon_rocket: the chart with m = 0 is not flat")


def photon_rocket_pullback():
    """The rectilinear chart is the Robinson-Trautman one for a rocket of rapidity w(u) along the
    first axis, U = (sinh w, 0, 0), with alpha = d_u w, pulled back through the aberration of
    the rays, tan(theta'/2) = e^{-w} cot(theta/2), where theta' is the angle in the background
    frame from the direction of motion and theta the angle in the rest frame from the direction
    opposite to the acceleration. Both metrics are compared in the coordinates u, r,
    tau = tan(theta/2) and phi, with E = e^w, in which each is rational."""
    general, straight = (cp.Chart(spec["system"]["coords"], spec["system"]["parameters"], spec["chart_line_element"])
                         for spec in (photon_rocket("robinson_trautman"), photon_rocket("rectilinear")))
    u, r, theta, phi = straight.symbols
    U, R, Theta, Phi = general.symbols
    tau = sp.Symbol("tau", positive=True)
    E = sp.Function("E", positive=True)(u)
    m = straight.reader.parameters["m"]
    T = 1 / (E * tau)                                   # tan(theta'/2)
    velocity = {general.reader.parameters[f"U_{i}"]: value
                for i, value in ((1, (E - 1 / E) / 2), (2, sp.Integer(0)), (3, sp.Integer(0)))}
    source = general.reader.surface(general.geo.g).subs(general.reader.parameters["m"], m)
    source = source.subs(velocity).doit()
    energy = (E + 1 / E) / 2                            # U_0 = cosh w, the root of 1 + sinh^2 w

    def rooted(power):
        if sp.simplify(power.base - energy ** 2) != 0:
            raise AssertionError("photon_rocket: a square root in the pullback is not U_0")
        return energy ** (2 * power.exp)
    source = source.replace(lambda e: e.is_Pow and abs(e.exp) == sp.Rational(1, 2), rooted)
    source = source.subs({sp.cos(Theta): (1 - T ** 2) / (1 + T ** 2), sp.sin(Theta): 2 * T / (1 + T ** 2)})
    source = source.subs({U: u, R: r, Phi: phi})
    images = [u, r, 2 * sp.atan(T), phi]
    jacobian = sp.Matrix(4, 4, lambda a, b: sp.diff(images[a], (u, r, tau, phi)[b]))
    pulled = jacobian.T * source * jacobian
    angle = {sp.cos(theta): (1 - tau ** 2) / (1 + tau ** 2), sp.sin(theta): 2 * tau / (1 + tau ** 2)}
    change = sp.diag(1, 1, 2 / (1 + tau ** 2), 1)
    target = straight.geo.g.subs(straight.reader.parameters["alpha"], sp.diff(E, u) / E).doit().subs(angle)
    target = change.T * target * change
    for a in range(4):
        for b in range(a, 4):
            if sp.simplify(pulled[a, b] - target[a, b]) != 0:
                raise AssertionError("photon_rocket: the pullback of the Robinson-Trautman chart misses the "
                                     f"rectilinear chart in slot {straight.coords_tex[a]}{straight.coords_tex[b]}")


def photon_rocket_charts():
    photon_rocket_pullback()
    return [photon_rocket("rectilinear"), photon_rocket("robinson_trautman")]


CHARTS["photon_rocket"] = photon_rocket_charts


# -- The spinning cosmic string --------------------------------------------------------

SPINNING_CHARTS = ("proper_radius", "rescaled_radius", "circumference_radius", "helical", "extended_source")


def spinning_string(system):
    """The spinning cosmic string of Deser, Jackiw and 't Hooft, with a = 4GJ/c^3 and
    b = 1 - 4G mu/c^2, in four charts of its flat exterior and one of an extended source:

    proper_radius         Deser and Jackiw's form of 1992, (c dt + a dphi)^2 with r the proper radius;
    rescaled_radius       rho = b r, the form of Cornish and Frankel's (1.2) and of Mena, Natario and Tod;
    circumference_radius  R^2 = b^2 r^2 - a^2, Cornish and Frankel's (1.6), whose g_phiphi is R^2;
    helical               tau = t + a phi/c and phi~ = b phi, Deser and Jackiw's (9), Minkowski's line
                          element with a circuit of the string identifying (tau, phi~) with
                          (tau + 2 pi a/c, phi~ + 2 pi b);
    extended_source       Soleng's frame of 1994, (c dt + M dphi)^2 and rho^2 dphi^2 with M(r) and
                          rho(r) free, which is the proper radius chart at M = a and rho = b r.

    Each exterior chart is checked, slot by slot, to be the proper radius chart pulled back, and
    the extended source's Einstein tensor is checked in Soleng's frame against his equation (4)
    without torsion. spinning_string.md is the derivation."""
    reals = "(-\\infty, \\infty)"
    flat = "\\text{flat: every curvature tensor vanishes}"
    string = "= 0 \\;\\text{(the string)}"
    twist = "-\\left(c\\,dt + a\\,d\\phi\\right)^2"
    charts = {
        "proper_radius": {
            "name": "Proper Radius", "coords": ["t", "r", "\\phi", "z"], "parameters": ["a", "b"],
            "domains": ["t \\in " + reals, "r \\in (0, \\infty)", "\\phi \\in [0, 2\\pi)", "z \\in " + reals,
                        "r " + string,
                        "r < a/b \\;\\text{(the circles of constant } t, r, z \\text{ are closed timelike curves)}",
                        flat],
            "line_element": "ds^2 = " + twist + " + dr^2 + b^2r^2d\\phi^2 + dz^2"},
        "rescaled_radius": {
            "name": "Rescaled Radius", "coords": ["t", "\\rho", "\\phi", "z"], "parameters": ["a", "b"],
            "domains": ["t \\in " + reals, "\\rho \\in (0, \\infty)", "\\phi \\in [0, 2\\pi)", "z \\in " + reals,
                        "\\rho " + string,
                        "\\rho < a \\;\\text{(the circles of constant } t, \\rho, z \\text{ are closed timelike curves)}",
                        flat],
            "line_element": "ds^2 = " + twist + " + \\dfrac{d\\rho^2}{b^2} + \\rho^2d\\phi^2 + dz^2"},
        "circumference_radius": {
            "name": "Circumference Radius", "coords": ["t", "R", "\\phi", "z"], "parameters": ["a", "b"],
            "domains": ["t \\in " + reals, "R \\in (0, \\infty)", "\\phi \\in [0, 2\\pi)", "z \\in " + reals,
                        "R = 0 \\;\\text{(the circle of constant } t, R, z \\text{ is null)}", flat],
            "line_element": ("ds^2 = -c^2dt^2 - 2ac\\,dt\\,d\\phi + \\dfrac{R^2dR^2}{b^2\\left(R^2 + a^2\\right)}"
                             " + R^2d\\phi^2 + dz^2")},
        "helical": {
            "name": "Helical Time", "coords": ["\\tau", "r", "\\tilde\\phi", "z"], "parameters": ["a", "b"],
            "domains": ["\\tau \\in " + reals, "r \\in (0, \\infty)", "\\tilde\\phi \\in [0, 2\\pi b)", "z \\in " + reals,
                        "(\\tau, \\tilde\\phi) \\sim (\\tau + 2\\pi a/c, \\tilde\\phi + 2\\pi b)", "r " + string, flat],
            "line_element": "ds^2 = -c^2d\\tau^2 + dr^2 + r^2d\\tilde\\phi^2 + dz^2"},
        "extended_source": {
            "name": "Extended Source", "coords": ["t", "r", "\\phi", "z"],
            "parameters": ["M = M(r)", "\\rho = \\rho(r)", "r_0"],
            "domains": ["t \\in " + reals, "r \\in [0, r_0]", "\\phi \\in [0, 2\\pi)", "z \\in " + reals,
                        "M(0) = 0,\\; M(r_0) = a",
                        "M^2 > \\rho^2 \\;\\text{(the circles of constant } t, r, z \\text{ are closed timelike curves)}"],
            "line_element": "ds^2 = -\\left(c\\,dt + M\\,d\\phi\\right)^2 + dr^2 + \\rho^2d\\phi^2 + dz^2"},
    }
    chart = charts[system]
    probe = vm.Reader(chart["coords"], chart["parameters"], ())
    lead = [probe.symbol[chart["coords"][1]], *probe.parameters.values()]
    return {
        "metric_id": "spinning_string",
        "system": {"id": system, "name": chart["name"], "coords": chart["coords"], "domains": chart["domains"],
                   "parameters": chart["parameters"], "line_element": chart["line_element"]},
        # The chart coordinate is x^0 = ct, so c dt is the chart's dt and 2ac dt is its 2a dt.
        "chart_line_element": chart["line_element"].replace("2ac\\,dt", "2a\\,dt").replace("c\\,dt", "dt").replace("c^2d", "d"),
        "printer": {"lead": lead},
        "check": spinning_string_source if system == "extended_source" else
                 (lambda c, system=system: spinning_string_pullback(c, system)) if system != "proper_radius" else
                 (lambda c: None),
    }


def spinning_string_pullback(chart, system):
    """J^T g J, with g the proper radius chart's metric and J the Jacobian of the map onto it,
    against this chart's metric in every slot: r = rho/b, r = sqrt(R^2 + a^2)/b, and
    t = tau - a phi~/b with phi = phi~/b, the time already the chart's x^0."""
    spec = spinning_string("proper_radius")
    source = cp.Chart(spec["system"]["coords"], spec["system"]["parameters"], spec["chart_line_element"])
    a, b = (chart.reader.parameters[n] for n in ("a", "b"))
    x0, x1, x2, x3 = chart.symbols
    image = {"rescaled_radius": [x0, x1 / b, x2, x3],
             "circumference_radius": [x0, sp.sqrt(x1 ** 2 + a ** 2) / b, x2, x3],
             "helical": [x0 - a * x2 / b, x1, x2 / b, x3]}[system]
    at = {source.reader.parameters[n]: chart.reader.parameters[n] for n in ("a", "b")}
    at.update(dict(zip(source.symbols, image)))
    J = sp.Matrix(4, 4, lambda i, j: sp.diff(image[i], chart.symbols[j]))
    pulled = J.T * source.geo.g.subs(at, simultaneous=True) * J
    for i in range(4):
        for j in range(i, 4):
            if sp.simplify(pulled[i, j] - chart.geo.g[i, j]) != 0:
                raise AssertionError(f"spinning_string: the proper radius chart pulled back misses the {system} "
                                     f"chart in slot {chart.coords_tex[i]}{chart.coords_tex[j]}")
    if system == "helical" and chart.geo.g != sp.diag(-1, 1, x1 ** 2, 1):
        raise AssertionError("spinning_string: the helical chart's metric is not Minkowski's")


def spinning_string_source(chart):
    """The extended source against Soleng's equation (4) of 1994 without torsion, in his frame
    e^0 = dt + M dphi, e^1 = dr, e^2 = rho dphi, e^3 = dz, with Omega = M'/(2 rho):
    G^0_0 = -3 Omega^2 + rho''/rho, G^1_1 = G^2_2 = Omega^2, G^3_3 = -Omega^2 + rho''/rho and a
    heat flow G^0_2 of magnitude Omega'. At M = a and rho = b r the metric is the proper radius
    chart's."""
    t, r, phi, z = chart.symbols
    M, rho = chart.reader.parameters["M"], chart.reader.parameters["rho"]
    e = sp.Matrix([[1, 0, M, 0], [0, 1, 0, 0], [0, 0, rho, 0], [0, 0, 0, 1]])
    eta = sp.diag(-1, 1, 1, 1)
    if (e.T * eta * e - chart.geo.g).applyfunc(sp.simplify) != sp.zeros(4, 4):
        raise AssertionError("spinning_string: Soleng's frame is not orthonormal in the extended source chart")
    lowered = chart.geo.einstein_ll()
    G = chart.geo.ginv * sp.Matrix(4, 4, lambda i, j: lowered[i][j])
    frame = (e * G * e.inv()).applyfunc(sp.simplify)
    Om = sp.diff(M, r) / (2 * rho)
    ddrho = sp.diff(rho, r, 2) / rho
    want = {(0, 0): -3 * Om ** 2 + ddrho, (1, 1): Om ** 2, (2, 2): Om ** 2, (3, 3): -Om ** 2 + ddrho}
    for i in range(4):
        for j in range(4):
            if (i, j) in ((0, 2), (2, 0)):
                if sp.simplify(frame[i, j] ** 2 - sp.diff(Om, r) ** 2) != 0:
                    raise AssertionError(f"spinning_string: the heat flow G^{i}_{j} is not Omega'")
            elif sp.simplify(frame[i, j] - want.get((i, j), 0)) != 0:
                raise AssertionError(f"spinning_string: G^{i}_{j} in Soleng's frame is {frame[i, j]}")
    spec = spinning_string("proper_radius")
    source = cp.Chart(spec["system"]["coords"], spec["system"]["parameters"], spec["chart_line_element"])
    a, b = source.reader.parameters["a"], source.reader.parameters["b"]
    exterior = chart.geo.g.subs({M: a, rho: b * r}).doit()
    if (exterior - source.geo.g.subs(dict(zip(source.symbols, chart.symbols)))).applyfunc(sp.simplify) != sp.zeros(4, 4):
        raise AssertionError("spinning_string: the extended source at M = a, rho = b r is not the exterior")


CHARTS["spinning_string"] = [lambda s=s: spinning_string(s) for s in SPINNING_CHARTS]
# -- Bonnor's beam of light ---------------------------------------------------------------

LIGHT_BEAM_CHARTS = ("cartesian", "null_cartesian", "null_cylindrical_interior", "null_cylindrical_exterior")


def light_beam(system):
    """The four charts of William Bonnor's steady beams of light, a pp-wave whose profile A does
    not depend on u and is not harmonic where the light is. The Cartesian chart is his equation
    (2.13) and the null Cartesian chart his (2.1), with sqrt(2) u = ct - z and sqrt(2) v = ct + z,
    both with A a free function of x and y, written here in the signature (-,+,+,+), so that his
    A >= 0 enters with a minus sign. The null cylindrical charts are his uniform beam of circular
    cross section, equations (5.6) and (5.5), in the polar coordinates of his (8.10):
    A = 4 pi G epsilon rho^2/c^4 inside the radius R and 4 pi G epsilon R^2 (1 + 2 ln(rho/R))/c^4
    outside it, his m being pi G epsilon R^2/c^4. light_beam.md derives each."""
    reals = "(-\\infty, \\infty)"
    k = "\\dfrac{8\\pi G\\epsilon}{c^4}"
    log = "\\left(1 + 2\\ln\\left(\\dfrac{\\rho}{R}\\right)\\right)"
    charts = {
        "cartesian": {
            "name": "Cartesian", "coords": ["t", "x", "y", "z"], "parameters": ["A = A(x,y)"],
            "domains": ["t \\in " + reals, "x \\in " + reals, "y \\in " + reals, "z \\in " + reals],
            "line_element": "ds^2 = -c^2dt^2 + dx^2 + dy^2 + dz^2 - A\\left(c\\,dt - dz\\right)^2",
            "chart": "ds^2 = -dt^2 + dx^2 + dy^2 + dz^2 - A\\left(dt - dz\\right)^2"},
        "null_cartesian": {
            "name": "Null Cartesian", "coords": ["u", "v", "x", "y"], "parameters": ["A = A(x,y)"],
            "domains": ["u \\in " + reals, "v \\in " + reals, "x \\in " + reals, "y \\in " + reals],
            "line_element": "ds^2 = -2\\,du\\,dv + dx^2 + dy^2 - 2A\\,du^2"},
        "null_cylindrical_interior": {
            "name": "Null Cylindrical, Inside the Beam", "coords": ["u", "v", "\\rho", "\\phi"],
            "parameters": ["G", "\\epsilon", "R"],
            "domains": ["u \\in " + reals, "v \\in " + reals, "\\rho \\in [0, R]", "\\phi \\in [0, 2\\pi)"],
            "line_element": "ds^2 = -2\\,du\\,dv + d\\rho^2 + \\rho^2d\\phi^2 - " + k + "\\rho^2\\,du^2"},
        "null_cylindrical_exterior": {
            "name": "Null Cylindrical, Outside the Beam", "coords": ["u", "v", "\\rho", "\\phi"],
            "parameters": ["G", "\\epsilon", "R"],
            "domains": ["u \\in " + reals, "v \\in " + reals, "\\rho \\in [R, \\infty)", "\\phi \\in [0, 2\\pi)"],
            "line_element": "ds^2 = -2\\,du\\,dv + d\\rho^2 + \\rho^2d\\phi^2 - " + k + "R^2" + log + "du^2"},
    }
    chart = charts[system]
    probe = vm.Reader(chart["coords"], chart["parameters"], ())
    spec = {
        "metric_id": "light_beam",
        "system": {"id": system, "name": chart["name"], "coords": chart["coords"], "domains": chart["domains"],
                   "parameters": chart["parameters"], "line_element": chart["line_element"]},
        "chart_line_element": chart.get("chart", chart["line_element"]),
        "check": lambda chart: light_beam_null_dust(chart, system),
    }
    if system == "cartesian":
        # Minkowski's metric and inverse plus the beam, each term of it carrying (c dt - dz)^2 or its dual.
        spec["components"] = {
            "metric_components": {("t", "t"): "-1 - A", ("z", "z"): "1 - A"},
            "inverse_metric_components": {("t", "t"): "-1 + A", ("z", "z"): "1 + A"}}
        spec["printer"] = {"lead": [probe.parameters["A"]]}
        spec["time"] = "t"
        # The geodesics grouped around ct - z, whose rate is the one the beam sees.
        along = ("\\left(\\partial_x A\\,\\dot{x} + \\partial_y A\\,\\dot{y}\\right)"
                 "\\left(\\dot{t} - \\dot{z}\\right) = 0")
        across = "\\dfrac{\\partial_{0} A}{2}\\left(\\dot{t} - \\dot{z}\\right)^2 = 0"
        spec["geodesics"] = ["\\ddot{t} + " + along, "\\ddot{x} + " + across.replace("{0}", "x"),
                             "\\ddot{y} + " + across.replace("{0}", "y"), "\\ddot{z} + " + along]
    elif system == "null_cartesian":
        spec["printer"] = {"lead": [probe.parameters["A"]]}
    else:
        G, eps = probe.parameters["G"], probe.parameters["epsilon"]
        # pi, G and epsilon enter every component as their product, printed together as the line element has it.
        piGe = sp.Symbol("piGe")
        spec["printer"] = {"lead": [probe.c, piGe, probe.parameters["R"], *probe.symbol.values()], "flip": False,
                           "overrides": {piGe: "\\pi G\\epsilon"}}
        spec["pretty"] = lambda value: sp.factor(value).subs(eps, piGe / (sp.pi * G))
        if system == "null_cylindrical_exterior":
            profile = k + "R^2" + log
            spec["components"] = {"metric_components": {("u", "u"): "-" + profile},
                                  "inverse_metric_components": {("v", "v"): profile}}
    return spec


def light_beam_null_dust(chart, system):
    """Einstein's equations with light for the source: G_mu nu = (8 pi G/c^4) epsilon k_mu k_nu, with k the
    null covector c dt - dz, which is sqrt(2) du, and the energy density epsilon = c^4 (A_xx + A_yy)/(16 pi G)
    for a free profile, Bonnor's (2.10), the declared epsilon inside the uniform beam and zero outside it.
    k^mu is checked null, a Killing vector, and covariantly constant, which makes the metric a pp-wave."""
    x, names = chart.symbols, chart.coords_tex
    g, ginv = chart.geo.g, chart.geo.ginv
    if "t" in names:
        k = sp.Matrix([1 if n == "t" else -1 if n == "z" else 0 for n in names])
    else:
        k = sp.Matrix([sp.sqrt(2) if n == "u" else 0 for n in names])
    if "A" in chart.reader.parameters:
        A = chart.reader.parameters["A"]
        density = sum(sp.diff(A, s, 2) for n, s in zip(names, x) if n in ("x", "y")) / 2
    elif system == "null_cylindrical_interior":
        density = 8 * sp.pi * chart.reader.parameters["G"] * chart.reader.parameters["epsilon"] / chart.reader.c ** 4
    else:
        density = 0
    G = chart.geo.einstein_ll()
    for a in range(4):
        for b in range(a, 4):
            if vm.norm(vm._at(G, (a, b)) - density * k[a] * k[b]) != 0:
                raise AssertionError(f"light_beam: the Einstein tensor misses the light's stress in slot "
                                     f"{names[a]}{names[b]}")
    up = ginv * k
    if vm.norm((k.T * up)[0]) != 0:
        raise AssertionError("light_beam: the rays of the beam are not null")
    gamma = chart.geo.christoffel_ull()
    for a in range(4):
        for b in range(4):
            nabla = sp.diff(up[a], x[b]) + sum(vm._at(gamma, (a, b, c)) * up[c] for c in range(4))
            if vm.norm(nabla) != 0:
                raise AssertionError(f"light_beam: the rays of the beam are not covariantly constant, slot "
                                     f"{names[a]}{names[b]}")


CHARTS["light_beam"] = [lambda s=s: light_beam(s) for s in LIGHT_BEAM_CHARTS]


# -- The wormhole time machine ---------------------------------------------------------

WTM_CHARTS = ("lorentz", "wormhole", "short_throat")
WTM_SPHERE = "\\left(d\\theta^2 + \\sin^2\\theta\\,d\\phi^2\\right)"


def wormhole_time_machine(system):
    """Morris, Thorne and Yurtsever's wormhole with one mouth accelerating, in three charts.

    The Lorentz chart of the flat space outside the mouths, where the physics is in the domains:
    a ball of radius b removed about each mouth's world line and the two spheres identified at
    equal proper times, the model of Friedman, Morris, Novikov, Echeverria, Klinkhammer, Thorne
    and Yurtsever (1990). Morris, Thorne and Yurtsever's own chart through the wormhole (1988),
    the static wormhole's Phi(l) and r(l) with the lapse N = 1 + g l F cos(theta)/c^2 of the right
    mouth's accelerated frame\\; N is held as a function of t, l and theta while the tensors are
    built and every value is written in N and its derivatives, so no component assumes a form of
    g, F, Phi or r. And the mouth of the 1990 model in its rest frame, r = b + |l|, whose
    curvature is the delta function of their equation (2). wormhole_time_machine.md derives each."""
    reals = "(-\\infty, \\infty)"
    angles = ["\\theta \\in [0, \\pi]", "\\phi \\in [0, 2\\pi)"]
    if system == "lorentz":
        coords, parameters = ["T", "X", "Y", "Z"], ["b"]
        probe = vm.Reader(coords, parameters, ())
        return {
            "metric_id": "wormhole_time_machine",
            "system": {"id": system, "name": "External Lorentz", "coords": coords,
                       "domains": [f"{x} \\in {reals}" for x in coords] + [
                           "|\\mathbf{X} - \\mathbf{X}_A(\\tau)| \\le b \\;\\text{in the rest frame of mouth } A = 1, 2 "
                           "\\;\\text{(the balls removed)}",
                           "\\mathcal{B}_1(\\tau) \\sim \\mathcal{B}_2(\\tau) "
                           "\\;\\text{(the two spheres at equal proper time)}",
                           "\\text{flat: every curvature tensor vanishes}"],
                       "parameters": parameters, "line_element": "ds^2 = -c^2dT^2 + dX^2 + dY^2 + dZ^2"},
            "chart_line_element": "ds^2 = -dT^2 + dX^2 + dY^2 + dZ^2",
            "printer": {"lead": [probe.c, probe.symbol["T"]]},
            "time": "T",
        }
    if system == "wormhole":
        coords = ["t", "l", "\\theta", "\\phi"]
        parameters = ["g = g(t)", "F = F(l)", "\\Phi = \\Phi(l)", "r = r(l)", "N = 1 + \\dfrac{glF\\cos\\theta}{c^2}"]
        rest = " + dl^2 + r^2" + WTM_SPHERE
        probe = vm.Reader(coords, parameters, ())
        P = probe.parameters
        return {
            "metric_id": "wormhole_time_machine",
            "system": {"id": system, "name": "Morris-Thorne-Yurtsever", "coords": coords,
                       "domains": ["t \\in " + reals, "l \\in " + reals] + angles + [
                           "l = 0 \\;\\text{(throat)}",
                           "g_{\\max}S \\ll c^2 \\;\\text{(a wormhole of length } S \\text{ from mouth to mouth)}"],
                       "parameters": parameters,
                       "line_element": "ds^2 = -\\left(1 + \\dfrac{glF\\cos\\theta}{c^2}\\right)^2e^{2\\Phi}c^2dt^2" + rest},
            "chart_line_element": "ds^2 = -N^2e^{2\\Phi}dt^2" + rest,
            "printer": {"lead": [P["N"], P["Phi"], P["r"]]},
            "check": wormhole_accelerated_frame,
        }
    coords, parameters = ["t", "l", "\\theta", "\\phi"], ["b"]
    sphere = " + \\left(b + |l|\\right)^2" + WTM_SPHERE
    probe = vm.Reader(coords, parameters, ())
    ell, b = probe.symbol["l"], probe.parameters["b"]
    # chart_printer.kink prints each value as A(|l|) + sgn(l) B(|l|), the delta's term last.
    delta = sp.Symbol("_delta" + ell.name)

    def then(e):
        return sp.factor(e.subs(delta, 0)) + sp.factor(sp.diff(e, delta)) * delta

    pretty, overrides = cp.kink(ell, then)
    absolute = next(p for p in overrides if p.name.startswith("_abs"))
    return {
        "metric_id": "wormhole_time_machine",
        "system": {"id": system, "name": "Short Throat", "coords": coords,
                   "domains": ["t \\in " + reals, "l \\in " + reals] + angles + ["l = 0 \\;\\text{(throat)}"],
                   "parameters": parameters, "line_element": "ds^2 = -c^2dt^2 + dl^2" + sphere},
        "chart_line_element": "ds^2 = -dt^2 + dl^2" + sphere,
        "printer": {"lead": [b, absolute], "flip": False, "last": [delta], "overrides": overrides},
        "pretty": pretty,
        "bracketed": pretty,
        "kretschmann_where": "l",
    }


def wormhole_accelerated_frame(chart):
    """Just outside the right mouth, where F = 1, Phi = 0 and r = l, the chart is flat space in
    the right mouth's accelerated frame. Morris, Thorne and Yurtsever's transformation to the
    Lorentz coordinates outside, with c = 1 and the mouth's rapidity eta, v = tanh(eta),

        T = T_R + v gamma l cos(theta),    Z = Z_R + gamma l cos(theta),
        X = l sin(theta) cos(phi),         Y = l sin(theta) sin(phi),

    where dT_R/dt = gamma and dZ_R/dt = v gamma along the mouth's world line, pulls
    -dT^2 + dX^2 + dY^2 + dZ^2 back to -(1 + g l cos(theta))^2 dt^2 + dl^2 + l^2 dOmega^2 with
    g = d eta/dt = gamma^2 dv/dt, in every slot. And the chart's own metric is that one with
    N^2 e^(2 Phi) for the lapse squared and r for l in the angular part."""
    t, l, theta, phi = sp.symbols("t l theta phi", real=True)
    eta, TR, ZR = (sp.Function(name)(t) for name in ("eta", "T_R", "Z_R"))
    x = [t, l, theta, phi]
    X = [TR + sp.sinh(eta) * l * sp.cos(theta), l * sp.sin(theta) * sp.cos(phi),
         l * sp.sin(theta) * sp.sin(phi), ZR + sp.cosh(eta) * l * sp.cos(theta)]
    J = sp.Matrix(4, 4, lambda i, j: sp.diff(X[i], x[j])).subs({sp.Derivative(TR, t): sp.cosh(eta),
                                                               sp.Derivative(ZR, t): sp.sinh(eta)})
    pulled = (J.T * sp.diag(-1, 1, 1, 1) * J).applyfunc(sp.simplify)
    lapse = 1 + sp.Derivative(eta, t) * l * sp.cos(theta)
    if (pulled - sp.diag(-lapse ** 2, 1, l ** 2, l ** 2 * sp.sin(theta) ** 2)).applyfunc(sp.simplify) != sp.zeros(4, 4):
        raise AssertionError("the accelerated frame outside the right mouth is not the Lorentz chart pulled back")
    R, P = chart.reader, chart.reader.parameters
    th = R.symbol["\\theta"]
    expected = sp.diag(-P["N"] ** 2 * sp.exp(2 * P["Phi"]), 1, P["r"] ** 2, P["r"] ** 2 * sp.sin(th) ** 2)
    if (chart.geo.g - expected).applyfunc(vm.norm) != sp.zeros(4, 4):
        raise AssertionError("the wormhole chart's metric is not -N^2 e^(2 Phi) dt^2 + dl^2 + r^2 dOmega^2")


CHARTS["wormhole_time_machine"] = [lambda s=s: wormhole_time_machine(s) for s in WTM_CHARTS]


# -- The Damour-Solodukhin wormhole ------------------------------------------------------

DS_SPHERE = "\\left(d\\theta^2 + \\sin^2\\theta\\,d\\phi^2\\right)"
DS_ANGLES = ["\\theta \\in [0, \\pi]", "\\phi \\in [0, 2\\pi)"]
DS_REALS = "(-\\infty, \\infty)"
DS_THROAT = "\\;\\text{(throat)}"
DS_FAR = "\\;\\text{(the far end of the other side)}"
DS_CHARTS = ["spherical", "rescaled", "throat", "isotropic", "einstein_rosen"]


def damour_solodukhin(system):
    """Damour and Solodukhin's wormhole, Schwarzschild's metric with the constant lambda^2 added
    to -g_tt, in five charts, each with r_s the areal radius of the throat. Their own keeps
    Schwarzschild's t and r on either side of the throat r = r_s, where g_rr diverges and g_tt =
    -lambda^2 does not vanish. Bueno, Cano, Goelen, Hertog and Vercnocke's time T = sqrt(1 +
    lambda^2) t is the proper time of a clock at rest far away, their (35) with their mass
    M = r_s c^2/2G(1 + lambda^2). Their auxiliary variable rho, r = r_s(2 + lambda^2(1 + cosh
    rho))/2(1 + lambda^2), runs over the whole line through the throat rho = 0, and in it the
    plane of T and rho is conformally flat with dr_* = r d rho, their (40) and (41). Nandi,
    Karimov, Izmailov and Potapov's isotropic radius, with the areal radius r(1 + r_s/4r)^2, has
    the throat at r_s/4 and the other side inside it, their (21) and (22). Einstein and Rosen's
    coordinate u, u^2 = r - r_s, runs through the throat with a polynomial. Each chart after
    the first is checked to be the first pulled back, on both sides of the throat where it
    covers both. The energy density vanishes, since g_rr is Schwarzschild's, which is checked as
    G^t_t = 0, and at the throat G^r_r = -1/r_s^2, G^theta_theta = (1 + 2 lambda^2)/(4 lambda^2
    r_s^2) and K = (1 + 24 lambda^4)/(4 lambda^4 r_s^4). The Kretschmann scalar is written as
    4A^2 + 8B^2 + 8C^2 + 4D^2 over the four frame components of a static observer, t r t r,
    t theta t theta, r theta r theta and theta phi theta phi. damour_solodukhin.md derives each."""
    parameters = ["r_s", "\\lambda"]
    one = sp.Symbol("DSone", positive=True)
    overrides = {one: "\\left(1 + \\lambda^2\\right)"}
    named, wholes, hyper, kretschmann = {}, [], None, None
    time = "t"
    if system in ("spherical", "rescaled"):
        coords = [time, "r", "\\theta", "\\phi"]
        domains = [time + " \\in " + DS_REALS, "r \\in [r_s, \\infty)"] + DS_ANGLES + ["r = r_s \\;\\text{(throat)}"]
        bare = "1 - \\dfrac{r_s}{r}"
        space = "\\dfrac{dr^2}{" + bare + "} + r^2" + DS_SPHERE
        W = "\\left(\\left(1 + \\lambda^2\\right)r - r_s\\right)"
        if system == "spherical":
            name = "Damour-Solodukhin"
            f = "\\left(1 - \\dfrac{r_s}{r} + \\lambda^2\\right)"
            weight = ""
        else:
            name = "Rescaled Time"
            f = "\\left(1 - \\dfrac{r_s}{\\left(1 + \\lambda^2\\right)r}\\right)"
            weight = ""
        components = {"metric_components": {(time, time): "-" + f, ("r", "r"): "\\left(" + bare + "\\right)^{-1}"},
                      "inverse_metric_components": {(time, time): "-" + f + "^{-1}", ("r", "r"): bare}}
        kretschmann = ("\\dfrac{r_s^2\\left(4\\left(r - r_s\\right)" + W + " - \\lambda^2\\,r\\,r_s\\right)^2}{4r^6" + W + "^4}"
                       " + \\dfrac{2r_s^2\\left(r - r_s\\right)^2}{r^6" + W + "^2} + \\dfrac{6r_s^2}{r^6}")
    elif system == "einstein_rosen":
        coords, name = ["t", "u", "\\theta", "\\phi"], "Einstein-Rosen"
        domains = ["t \\in " + DS_REALS, "u \\in " + DS_REALS] + DS_ANGLES + ["u = 0 \\;\\text{(throat)}"]
        f = "\\left(\\dfrac{u^2}{u^2 + r_s} + \\lambda^2\\right)"
        space = "4\\left(u^2 + r_s\\right)du^2 + \\left(u^2 + r_s\\right)^2" + DS_SPHERE
        components = {"metric_components": {("t", "t"): "-" + f, ("u", "u"): "4\\left(u^2 + r_s\\right)"},
                      "inverse_metric_components": {("t", "t"): "-" + f + "^{-1}"}}
        W = "\\left(\\left(1 + \\lambda^2\\right)u^2 + \\lambda^2\\,r_s\\right)"
        kretschmann = ("\\dfrac{r_s^2\\left(4u^2" + W + " - \\lambda^2\\,r_s\\left(u^2 + r_s\\right)\\right)^2}"
                       "{4\\left(u^2 + r_s\\right)^6" + W + "^4}"
                       " + \\dfrac{2u^4\\,r_s^2}{\\left(u^2 + r_s\\right)^6" + W + "^2}"
                       " + \\dfrac{6r_s^2}{\\left(u^2 + r_s\\right)^6}")
    elif system == "isotropic":
        coords, name = ["t", "r", "\\theta", "\\phi"], "Isotropic"
        domains = (["t \\in " + DS_REALS, "r \\in (0, \\infty)"] + DS_ANGLES
                   + ["r = r_s/4 " + DS_THROAT, "r \\to 0 " + DS_FAR])
        f = "\\left(\\left(\\dfrac{4r - r_s}{4r + r_s}\\right)^2 + \\lambda^2\\right)"
        space = "\\left(1 + \\dfrac{r_s}{4r}\\right)^4\\left(dr^2 + r^2" + DS_SPHERE + "\\right)"
        components = {"metric_components": {("t", "t"): "-" + f, ("r", "r"): "\\left(1 + \\dfrac{r_s}{4r}\\right)^4"},
                      "inverse_metric_components": {("t", "t"): "-" + f + "^{-1}"}}
        W = "\\left(\\left(4r - r_s\\right)^2 + \\lambda^2\\left(4r + r_s\\right)^2\\right)"
        kretschmann = ("\\dfrac{16777216r^6\\,r_s^2}{\\left(4r + r_s\\right)^{12}}\\left(\\dfrac{4\\left(\\left(4r - r_s\\right)^2"
                       + W + " - 4\\lambda^2\\,r\\,r_s\\left(4r + r_s\\right)^2\\right)^2}{" + W + "^4}"
                       " + \\dfrac{2\\left(4r - r_s\\right)^4}{" + W + "^2} + 6\\right)")
    else:
        coords, name = ["t", "\\rho", "\\theta", "\\phi"], "Through the Throat"
        domains = ["t \\in " + DS_REALS, "\\rho \\in " + DS_REALS] + DS_ANGLES + ["\\rho = 0 \\;\\text{(throat)}"]
        up = "\\left(1 + \\cosh\\rho\\right)"
        whole = "2 + \\lambda^2" + up
        f = "\\dfrac{\\lambda^2" + up + "}{" + whole + "}"
        scale = "\\dfrac{r_s^2}{4\\left(1 + \\lambda^2\\right)^2}"
        space = (scale + "\\left(\\lambda^2" + up + "\\left(" + whole + "\\right)d\\rho^2 + \\left(" + whole
                 + "\\right)^2" + DS_SPHERE + "\\right)")
        components = {}
    probe = vm.Reader(coords, parameters, ())
    x, rs, lam = probe.symbol[coords[1]], probe.parameters["r_s"], probe.parameters["lambda"]
    lead = [one, x, rs, lam]
    if system in ("einstein_rosen", "isotropic"):
        placeholder = sp.Symbol("DSw")
        named[placeholder] = W.removeprefix("\\left(").removesuffix("\\right)")
        wholes = [(placeholder, (1 + lam ** 2) * x ** 2 + lam ** 2 * rs if system == "einstein_rosen"
                   else (4 * x - rs) ** 2 + lam ** 2 * (4 * x + rs) ** 2)]
    elif system == "throat":
        # Every value is a rational function of C = cosh rho with at most one S = sinh rho, and its
        # factors are 1 + cosh rho, cosh rho - 1 and the areal radius's 2 + lambda^2(1 + cosh rho).
        S, C = sp.symbols("DSs DSc", positive=True)
        hyper = (S, C)
        placeholder, lift = sp.Symbol("DSw"), sp.Symbol("DSu")
        named[placeholder], named[lift] = whole, "1 + \\cosh\\rho"
        wholes = [(placeholder, 2 + lam ** 2 * (1 + C)), (lift, 1 + C)]
        overrides.update({S: "\\sinh\\rho", C: "\\cosh\\rho"})
        lead = [one, C, S, rs, lam]

    def in_hyperbolics(value):
        S, C = hyper
        value = sp.sympify(value).subs({sp.sinh(x): S, sp.cosh(x): C})
        value = value.replace(lambda e: isinstance(e, sp.exp) and sp.expand(e.args[0] / x).is_Integer,
                              lambda e: (S + C) ** sp.expand(e.args[0] / x))

        def reduced(poly):
            return sp.expand(sp.reduced(sp.expand(poly), [S ** 2 - C ** 2 + 1], S, C)[1])

        num, den = (reduced(part) for part in sp.fraction(sp.together(value)))
        d0, d1 = sp.Poly(den, S).coeff_monomial(1), sp.Poly(den, S).coeff_monomial(S)
        if d1 != 0:
            num, den = reduced(num * (d0 - d1 * S)), sp.expand(d0 ** 2 - d1 ** 2 * (C ** 2 - 1))
        return sp.cancel(sp.factor(num) / sp.factor(den))

    def pretty(value):
        # 1 + lambda^2 is written whole wherever that shortens a factor, as (1 + lambda^2) r - r_s.
        if hyper:
            value = in_hyperbolics(value)
        out = sp.Integer(1)
        for factor in sp.Mul.make_args(sp.factor(value)):
            base, k = (factor.base, factor.exp) if factor.is_Pow else (factor, sp.Integer(1))
            if base.is_Add:
                for placeholder, poly in wholes:
                    if sp.expand(base - poly) == 0:
                        base = placeholder
                    elif sp.expand(base + poly) == 0 and k.is_Integer:
                        base, out = placeholder, out * (-1) ** k
            if base.is_Add and base.has(lam):
                short = sp.expand(base.subs(lam ** 2, one - 1))
                if not short.has(lam) and len(sp.Add.make_args(short)) < len(sp.Add.make_args(sp.expand(base))):
                    base = short
            out *= base ** k
        return out

    spec = {
        "metric_id": "damour_solodukhin",
        "system": {"id": system, "name": name, "coords": coords, "domains": domains, "parameters": parameters,
                   "line_element": "ds^2 = -" + f + "c^2dt^2 + " + space},
        "chart_line_element": "ds^2 = -" + f + "dt^2 + " + space,
        "printer": {"lead": lead, "factors": [one, lam] + lead[1:-1], "flip": False,
                    "overrides": overrides, "named": named},
        "pretty": pretty,
        "components": components,
        "check": lambda chart: damour_solodukhin_check(chart, system),
    }
    if hyper:
        back = {hyper[0]: sp.sinh(x), hyper[1]: sp.cosh(x)}
        spec["pretty"] = lambda value: pretty(value).subs(back)
        spec["printer"]["overrides"] = {one: overrides[one]}
        spec["printer"]["lead"] = [one, rs, lam]
        spec["printer"]["factors"] = [one, lam, rs]
    if kretschmann:
        spec["kretschmann"] = kretschmann
    return spec


def damour_solodukhin_check(chart, system):
    """G^t_t = 0 in every chart, the energy density of the matter at the throat; at the throat
    G^r_r = -1/r_s^2, G^theta_theta = (1 + 2 lambda^2)/(4 lambda^2 r_s^2) and
    K = (1 + 24 lambda^4)/(4 lambda^4 r_s^4); at lambda = 0 the Ricci tensor vanishes; and each
    chart after the first is Damour and Solodukhin's pulled back: t = T/sqrt(1 + lambda^2) for
    the rescaled time, r = r_s(2 + lambda^2(1 + cosh rho))/2(1 + lambda^2) besides through the
    throat, the areal radius (4r + r_s)^2/16r for the isotropic radius, and r = r_s + u^2 for
    Einstein and Rosen's u, each on both sides of the throat."""
    rs, lam = chart.reader.parameters["r_s"], chart.reader.parameters["lambda"]
    t, x = chart.symbols[:2]
    throat = {x: rs} if system in ("spherical", "rescaled") else {x: rs / 4} if system == "isotropic" else {x: 0}
    G = chart.geo.raise_indices(chart.geo.einstein_ll(), 2, (0,))
    if vm.norm(vm._at(G, (0, 0))) != 0:
        raise AssertionError(f"damour_solodukhin: the energy density does not vanish in the {system} chart")
    expected = {(1, 1): -1 / rs ** 2, (2, 2): (1 + 2 * lam ** 2) / (4 * lam ** 2 * rs ** 2)}
    for slot, value in expected.items():
        if sp.simplify(sp.simplify(vm._at(G, slot)).subs(throat) - value) != 0:
            raise AssertionError(f"damour_solodukhin: G^a_b in slot {slot} at the throat is not {value}")
    K = sp.simplify(chart.geo.kretschmann())
    if sp.simplify(K.subs(throat) - (1 + 24 * lam ** 4) / (4 * lam ** 4 * rs ** 4)) != 0:
        raise AssertionError("damour_solodukhin: the Kretschmann scalar at the throat is not "
                             "(1 + 24 lambda^4)/(4 lambda^4 r_s^4)")
    ricci = chart.geo.ricci_ll()
    if system != "throat" and any(sp.simplify(sp.cancel(vm._at(ricci, index)).subs(lam, 0)) != 0
                                  for index in vm._indices(4, 2)):
        raise AssertionError(f"damour_solodukhin: the {system} chart is not vacuum at lambda = 0")
    if system == "spherical":
        return
    source = damour_solodukhin("spherical")
    areal = cp.Chart(source["system"]["coords"], source["system"]["parameters"], source["chart_line_element"])
    names = dict(zip(areal.symbols[2:], chart.symbols[2:]))
    names.update({areal.reader.parameters[n]: chart.reader.parameters[n] for n in ("r_s", "lambda")})
    image = {"rescaled": [t / sp.sqrt(1 + lam ** 2), x],
             "throat": [t / sp.sqrt(1 + lam ** 2), rs * (2 + lam ** 2 * (1 + sp.cosh(x))) / (2 * (1 + lam ** 2))],
             "isotropic": [t, (4 * x + rs) ** 2 / (16 * x)],
             "einstein_rosen": [t, rs + x ** 2]}[system]
    J = sp.diag(sp.diff(image[0], t), sp.diff(image[1], x), 1, 1)
    at = dict(zip(areal.symbols[:2], image))
    pulled = J.T * areal.geo.g.subs(names).subs(at, simultaneous=True) * J
    difference = (pulled - chart.geo.g).applyfunc(lambda e: sp.simplify(e.rewrite(sp.exp)))
    if difference != sp.zeros(4, 4):
        raise AssertionError(f"damour_solodukhin: the {system} chart is not the wormhole's own pulled back")


CHARTS["damour_solodukhin"] = [lambda s=s: damour_solodukhin(s) for s in DS_CHARTS]

# -- The Simpson-Visser black bounce -----------------------------------------------------

SV_SPHERE = "\\left(d\\theta^2 + \\sin^2\\theta\\,d\\phi^2\\right)"
SV_ANGLES = ["\\theta \\in [0, \\pi]", "\\phi \\in [0, 2\\pi)"]
SV_REALS = "(-\\infty, \\infty)"
SV_RADIUS = "\\rho = \\sqrt{r^2 + a^2}"
SV_CHARTS = ["spherical", "areal", "eddington_finkelstein_outgoing", "eddington_finkelstein_ingoing"]


def simpson_visser(system):
    """Simpson and Visser's black bounce, Schwarzschild's metric with sqrt(r^2 + a^2) written for
    the areal radius and dr left as it is, in four charts, each with r_s = 2Gm/c^2. Their own
    keeps t and lets r run over the whole line; the sphere r = 0 has the least area,
    4 pi a^2, and f = 1 - r_s/sqrt(r^2 + a^2) vanishes at r = +-sqrt(r_s^2 - a^2) when a < r_s.
    Tsukamoto's standard radial coordinate rho = sqrt(r^2 + a^2), the areal radius, covers one
    side of r = 0 at a time, and in it g_tt and the spheres are Schwarzschild's. Simpson,
    Martin-Moruno and Visser's two Eddington-Finkelstein charts, u = ct - r_* and
    v = ct + r_* with dr_*/dr = 1/f, are regular on the horizons. The charts with r
    name rho = sqrt(r^2 + a^2) and print every value in r, a and rho. Each chart after the first
    is checked to be the first pulled back, and every chart against Simpson and Visser's
    G^t_t = a^2(rho - 2r_s)/rho^5, G^r_r = -a^2/rho^4, G^theta_theta = a^2(2 rho - r_s)/(2 rho^5),
    R = a^2(3r_s - 2 rho)/rho^5 and Kretschmann scalar, against Schwarzschild's vacuum at a = 0
    and against Ellis and Bronnikov's K = 12a^4/rho^8 at r_s = 0. simpson_visser.md derives each."""
    areal = system == "areal"
    parameters = ["r_s", "a"] + ([] if areal else [SV_RADIUS])
    if system == "spherical":
        coords, name = ["t", "r", "\\theta", "\\phi"], "Simpson-Visser"
        f = "\\left(1 - \\dfrac{r_s}{\\sqrt{r^2 + a^2}}\\right)"
        bare = "1 - \\dfrac{r_s}{\\sqrt{r^2 + a^2}}"
        rest = "\\dfrac{dr^2}{" + bare + "} + \\left(r^2 + a^2\\right)" + SV_SPHERE
        line, chart_line = "ds^2 = -" + f + "c^2dt^2 + " + rest, "ds^2 = -" + f + "dt^2 + " + rest
        domains = (["t \\in " + SV_REALS, "r \\in " + SV_REALS] + SV_ANGLES
                   + ["r = 0 \\;\\text{(the sphere of least area)}",
                      "r = \\pm\\sqrt{r_s^2 - a^2} \\;\\text{(the horizons, where } a < r_s\\text{)}"])
        g = "\\left(1 - \\dfrac{r_s}{\\rho}\\right)"
        components = {"metric_components": {("t", "t"): "-" + g, ("r", "r"): g + "^{-1}",
                                            ("\\theta", "\\theta"): "\\rho^2"},
                      "inverse_metric_components": {("t", "t"): "-" + g + "^{-1}", ("r", "r"): "1 - \\dfrac{r_s}{\\rho}",
                                                    ("\\theta", "\\theta"): "\\dfrac{1}{\\rho^2}"}}
    elif areal:
        coords, name = ["t", "\\rho", "\\theta", "\\phi"], "Areal Radius"
        f = "\\left(1 - \\dfrac{r_s}{\\rho}\\right)"
        rest = ("\\dfrac{d\\rho^2}{\\left(1 - \\dfrac{r_s}{\\rho}\\right)\\left(1 - \\dfrac{a^2}{\\rho^2}\\right)}"
                " + \\rho^2" + SV_SPHERE)
        line, chart_line = "ds^2 = -" + f + "c^2dt^2 + " + rest, "ds^2 = -" + f + "dt^2 + " + rest
        domains = (["t \\in " + SV_REALS, "\\rho \\in [a, \\infty)"] + SV_ANGLES
                   + ["\\rho = a \\;\\text{(the sphere of least area)}",
                      "\\rho = r_s \\;\\text{(the horizon, where } a < r_s\\text{)}"])
        components = {"metric_components": {
                          ("t", "t"): "-" + f,
                          ("\\rho", "\\rho"): f + "^{-1}\\left(1 - \\dfrac{a^2}{\\rho^2}\\right)^{-1}"},
                      "inverse_metric_components": {
                          ("t", "t"): "-" + f + "^{-1}",
                          ("\\rho", "\\rho"): f + "\\left(1 - \\dfrac{a^2}{\\rho^2}\\right)"}}
    else:
        null, sign = ("u", "-") if system == "eddington_finkelstein_outgoing" else ("v", "+")
        coords = [null, "r", "\\theta", "\\phi"]
        name = ("Outgoing" if null == "u" else "Ingoing") + " Eddington-Finkelstein"
        f = "\\left(1 - \\dfrac{r_s}{\\sqrt{r^2 + a^2}}\\right)"
        line = ("ds^2 = -" + f + "d" + null + "^2 " + sign + " 2\\,d" + null + "\\,dr + \\left(r^2 + a^2\\right)"
                + SV_SPHERE)
        chart_line = line
        domains = ([null + " \\in " + SV_REALS, "r \\in " + SV_REALS] + SV_ANGLES
                   + ["r = 0 \\;\\text{(the sphere of least area)}",
                      "r = \\pm\\sqrt{r_s^2 - a^2} \\;\\text{(the horizons, where } a < r_s\\text{)}"])
        one = "-1" if null == "u" else "1"
        g = "\\left(1 - \\dfrac{r_s}{\\rho}\\right)"
        components = {"metric_components": {(null, null): "-" + g, (null, "r"): one, ("r", null): one,
                                            ("\\theta", "\\theta"): "\\rho^2"},
                      "inverse_metric_components": {(null, "r"): one, ("r", null): one,
                                                    ("r", "r"): "1 - \\dfrac{r_s}{\\rho}",
                                                    ("\\theta", "\\theta"): "\\dfrac{1}{\\rho^2}"}}
    probe = vm.Reader(coords, parameters, ())
    x, rs, a = probe.symbol[coords[1]], probe.parameters["r_s"], probe.parameters["a"]
    if areal:
        # rho^2 - a^2 is written whole, the square of Simpson and Visser's r.
        whole = sp.Symbol("SVw")
        named = {whole: "\\rho^2 - a^2"}

        def pretty(value):
            out, powers = _factor_powers(sp.sympify(value))
            flipped = powers.pop(a - x, 0)
            down, up = powers.pop(x - a, 0) + flipped, powers.pop(x + a, 0)
            if down != up:
                raise AssertionError(f"rho - a and rho + a stand to different powers in {value}")
            out *= (-1) ** flipped * whole ** down
            for base, k in powers.items():
                out *= base ** k
            return out

        printer = {"lead": [x, rs, a], "factors": [rs, a, x], "flip": False, "named": named}
        kretschmann = ("\\dfrac{3r_s^2\\left(4\\rho^4 - 12a^2\\rho^2 + 11a^4\\right)"
                       " + 16a^2\\,r_s\\,\\rho\\left(\\rho^2 - 2a^2\\right) + 12a^4\\rho^2}{\\rho^{10}}")
    else:
        rho = sp.Symbol("rho", positive=True)
        pretty = simpson_visser_radius(x, a, rho)
        printer = {"lead": [rho, x, rs, a], "factors": [rs, a, x, rho], "flip": False}
        kretschmann = ("\\dfrac{3r_s^2\\left(4r^4 - 4a^2r^2 + 3a^4\\right) + 16a^2\\,r_s\\left(r^2 - a^2\\right)\\rho"
                       " + 12a^4\\rho^2}{\\rho^{10}}")
    spec = {
        "metric_id": "simpson_visser",
        "system": {"id": system, "name": name, "coords": coords, "domains": domains, "parameters": parameters,
                   "line_element": line},
        "chart_line_element": chart_line,
        "printer": printer,
        "pretty": pretty,
        "components": components,
        "check": lambda chart: simpson_visser_check(chart, system),
    }
    if kretschmann:
        spec["kretschmann"] = kretschmann
    return spec


def simpson_visser_radius(r, a, rho):
    """A pretty printer for a chart that names rho = sqrt(r^2 + a^2): every power of r^2 + a^2 is a
    power of rho, a^2 is written rho^2 - r^2 while the value is factored, so that r and rho are
    generators with no relation between them, and each factor is then written in whichever of its
    three forms has the fewest terms: in r, a and at most one rho, as Simpson and Visser write
    2r^2 - a^2, in r and rho, or in rho and a."""
    square = r ** 2 + a ** 2

    def forms(base):
        linear = sp.Poly(base, rho)
        with_a = sp.expand(sum(c * rho ** (k % 2) * square ** (k // 2) for (k,), c in linear.terms()))
        even = sp.Poly(base, r)
        candidates = [with_a, sp.expand(base)]
        if all(k % 2 == 0 for (k,), _ in even.terms()):
            candidates.append(sp.expand(sum(c * (rho ** 2 - a ** 2) ** (k // 2) for (k,), c in even.terms())))
        return candidates

    def pretty(value):
        value = sp.sympify(value)
        if value == 0:
            return value
        value = value.replace(lambda p: p.is_Pow and sp.expand(p.base - square) == 0, lambda p: rho ** (2 * p.exp))
        value = sp.factor(sp.together(value.subs(a, sp.sqrt(rho ** 2 - r ** 2))))
        if value.has(sp.sqrt(rho ** 2 - r ** 2)):
            raise AssertionError(f"an odd power of a stands in {value}")
        out, powers = _factor_powers(value)
        flipped = powers.pop(r - rho, 0)
        down, up = powers.pop(rho - r, 0) + flipped, powers.pop(rho + r, 0)
        if down != up:
            raise AssertionError(f"rho - r and rho + r stand to different powers in {value}")
        out *= (-1) ** flipped * a ** (2 * down)
        for base, k in powers.items():
            if base.is_Add:
                base = min(forms(base), key=lambda e: len(sp.Add.make_args(e)))
                if base.could_extract_minus_sign() and k.is_Integer:
                    base, out = -base, out * (-1) ** k
            out *= base ** k
        return out

    return pretty


def simpson_visser_check(chart, system):
    """Simpson and Visser's curvature in every chart, with rho the areal radius: G^t_t =
    a^2(rho - 2r_s)/rho^5, G^r_r = -a^2/rho^4 and G^theta_theta = a^2(2 rho - r_s)/(2 rho^5) in the
    diagonal charts, theirs with m = r_s/2, R = a^2(3r_s - 2 rho)/rho^5 and their Kretschmann
    scalar. At a = 0 the Ricci tensor vanishes, Schwarzschild's vacuum, and at
    r_s = 0 the Kretschmann scalar is Ellis and Bronnikov's 12a^4/rho^8. Each chart after the
    first is Simpson and Visser's pulled back, the areal chart on both sides of r = 0."""
    rs, a = chart.reader.parameters["r_s"], chart.reader.parameters["a"]
    x = chart.symbols[1]
    rho = x if system == "areal" else sp.sqrt(x ** 2 + a ** 2)
    r2 = rho ** 2 - a ** 2
    G = chart.geo.raise_indices(chart.geo.einstein_ll(), 2, (0,))
    expected = {(2, 2): a ** 2 * (2 * rho - rs) / (2 * rho ** 5), (3, 3): a ** 2 * (2 * rho - rs) / (2 * rho ** 5)}
    if system in ("spherical", "areal"):
        expected.update({(0, 0): a ** 2 * (rho - 2 * rs) / rho ** 5, (1, 1): -a ** 2 / rho ** 4})
    for slot, value in expected.items():
        if vm.norm(vm._at(G, slot) - value) != 0:
            raise AssertionError(f"simpson_visser: G^a_b in slot {slot} of the {system} chart is not {value}")
    if vm.norm(chart.geo.ricci_scalar() - a ** 2 * (3 * rs - 2 * rho) / rho ** 5) != 0:
        raise AssertionError(f"simpson_visser: the Ricci scalar of the {system} chart is not Simpson and Visser's")
    K = chart.geo.kretschmann()
    theirs = (3 * rs ** 2 * (4 * r2 ** 2 - 4 * a ** 2 * r2 + 3 * a ** 4) + 16 * a ** 2 * rs * (r2 - a ** 2) * rho
              + 12 * a ** 4 * rho ** 2) / rho ** 10
    if vm.norm(K - theirs) != 0:
        raise AssertionError(f"simpson_visser: the Kretschmann scalar of the {system} chart is not Simpson and Visser's")
    if vm.norm(theirs.subs(rs, 0) - 12 * a ** 4 / rho ** 8) != 0:
        raise AssertionError("simpson_visser: the Kretschmann scalar at r_s = 0 is not Ellis and Bronnikov's")
    ricci = chart.geo.ricci_ll()
    if any(sp.simplify(sp.sympify(vm._at(ricci, index)).subs(a, 0)) != 0 for index in vm._indices(4, 2)):
        raise AssertionError(f"simpson_visser: the {system} chart is not vacuum at a = 0")
    if system == "spherical":
        return
    source = simpson_visser("spherical")
    own = cp.Chart(source["system"]["coords"], source["system"]["parameters"], source["chart_line_element"])
    t, r = own.symbols[:2]
    names = dict(zip(own.symbols[2:], chart.symbols[2:]))
    names.update({own.reader.parameters[n]: chart.reader.parameters[n] for n in ("r_s", "a")})
    w = chart.symbols[0]
    if system == "areal":
        # r = +-sqrt(rho^2 - a^2), one sign on each side of the sphere of least area.
        maps = [([w, side * sp.sqrt(x ** 2 - a ** 2)], sp.diag(1, side * x / sp.sqrt(x ** 2 - a ** 2), 1, 1))
                for side in (1, -1)]
    else:
        # ct = u + r_* or v - r_*, with dr_*/dr = 1/f.
        side = 1 if system == "eddington_finkelstein_outgoing" else -1
        f = 1 - rs / sp.sqrt(x ** 2 + a ** 2)
        maps = [([w, x], sp.Matrix([[1, side / f, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]))]
    for image, J in maps:
        at = {r: image[1]}
        pulled = J.T * own.geo.g.subs(names).subs(at, simultaneous=True) * J
        # The areal radius is positive, so sqrt(rho^2) is rho.
        difference = (pulled - chart.geo.g).applyfunc(lambda e: sp.simplify(sp.powdenest(sp.simplify(e), force=True)))
        if difference != sp.zeros(4, 4):
            raise AssertionError(f"simpson_visser: the {system} chart is not Simpson and Visser's pulled back")


CHARTS["simpson_visser"] = [lambda s=s: simpson_visser(s) for s in SV_CHARTS]


# -- Ori's time machine ------------------------------------------------------------------

ORI_CHARTS = ("vacuum_core", "foliation", "brinkmann")


def ori_time_machine(system):
    """The vacuum core of Amos Ori's time machine of 2005 in three charts:

    vacuum_core  his equation (1), -2 dz dT + (f - T) dz^2 with f(x, y, z) free and z periodic;
    foliation    his equation (5), the core with f = a(x^2 - y^2)/2 in the time
                 t = T - a(x^2 - y^2)/2 + e(x^2 + y^2), whose surfaces t < 0 are spacelike;
    brinkmann    the same core with u = -2 exp(-z/2) and v = T exp(z/2), the plane wave
                 2a(x^2 - y^2)/u^2 du^2 - 2 du dv, on which a circuit of z is a boost.

    ori_core_check holds the first to his equations (2) and (3), and ori_pullback checks each of
    the other two, slot by slot, to be the first pulled back. ori_time_machine.md is the
    derivation."""
    reals = "(-\\infty, \\infty)"
    plane = "\\left(x^2 + y^2\\right)"
    charts = {
        "vacuum_core": {
            "name": "Vacuum Core", "coords": ["T", "x", "y", "z"], "parameters": ["f = f(x,y,z)", "L"],
            "domains": ["T \\in " + reals, "x \\in " + reals, "y \\in " + reals, "z \\in [0, L)",
                        "T = f \\;\\text{(the circle of constant } T, x, y \\text{ is null)}",
                        "T > f \\;\\text{(the circles of constant } T, x, y \\text{ are closed timelike curves)}"],
            "line_element": "ds^2 = dx^2 + dy^2 - 2\\,dz\\,dT + \\left(f - T\\right)dz^2"},
        "foliation": {
            "name": "Spacelike Foliation", "coords": ["t", "x", "y", "z"], "parameters": ["a", "e", "L"],
            "domains": ["t \\in " + reals, "x \\in " + reals, "y \\in " + reals, "z \\in [0, L)",
                        "x = y = t = 0 \\;\\text{(the closed null geodesic } N\\text{)}",
                        "t > e" + plane + " \\;\\text{(the circles of constant } t, x, y "
                        "\\text{ are closed timelike curves)}"],
            "line_element": ("ds^2 = dx^2 + dy^2 - 2\\,dz\\,dt + \\left(e" + plane + " - t\\right)dz^2"
                             " + 2\\left(\\left(2e - a\\right)x\\,dx + \\left(2e + a\\right)y\\,dy\\right)dz")},
        "brinkmann": {
            "name": "Brinkmann", "coords": ["u", "v", "x", "y"], "parameters": ["a", "L"],
            "domains": ["u \\in (-\\infty, 0)", "v \\in " + reals, "x \\in " + reals, "y \\in " + reals,
                        "(u, v) \\sim \\left(u\\exp(-L/2),\\, v\\exp(L/2)\\right)",
                        "uv < -a\\left(x^2 - y^2\\right) \\;\\text{(the closed curves of constant } uv, x, y "
                        "\\text{ are timelike)}"],
            "line_element": "ds^2 = \\dfrac{2a\\left(x^2 - y^2\\right)}{u^2}du^2 - 2\\,du\\,dv + dx^2 + dy^2"},
    }
    chart = charts[system]
    probe = vm.Reader(chart["coords"], chart["parameters"], ())
    spec = {
        "metric_id": "ori_time_machine",
        "system": {"id": system, "name": chart["name"], "coords": chart["coords"], "domains": chart["domains"],
                   "parameters": chart["parameters"], "line_element": chart["line_element"]},
        # No coordinate is a time, so the chart's line element is the published one.
        "chart_line_element": chart["line_element"],
        "printer": {"lead": [*probe.parameters.values()]},
        "check": ori_core_check if system == "vacuum_core" else (lambda c: ori_pullback(c, system)),
    }
    # The printer writes a - 2e and a + 2e; the values keep Ori's 2e - a and 2e + a, and his g^tt.
    down, up = "\\left(2e - a\\right)x", "\\left(2e + a\\right)y"
    if system == "vacuum_core":
        spec["components"] = {"inverse_metric_components": {("T", "T"): "T - f"}}
    elif system == "foliation":
        spec["components"] = {
            "metric_components": {("x", "z"): down, ("z", "x"): down, ("y", "z"): up, ("z", "y"): up,
                                  ("z", "z"): "e" + plane + " - t"},
            "inverse_metric_components": {
                ("t", "x"): down, ("x", "t"): down, ("t", "y"): up, ("y", "t"): up,
                ("t", "t"): "t + \\left(2e - a\\right)^2x^2 + \\left(2e + a\\right)^2y^2 - e" + plane}}
        tzz = "\\dfrac{t - \\left(e + 2ae - a^2\\right)x^2 - \\left(e - 2ae - a^2\\right)y^2}{2}"
        spec["rewrite"] = [
            ("-\\dfrac{a\\,x\\left(a - 2e\\right)}{2}", "\\dfrac{a\\left(2e - a\\right)x}{2}"),
            ("\\dfrac{a\\,x\\left(a - 2e\\right)}{2}", "-\\dfrac{a\\left(2e - a\\right)x}{2}"),
            ("\\dfrac{a\\,y\\left(a + 2e\\right)}{2}", "\\dfrac{a\\left(2e + a\\right)y}{2}"),
            ("\\dfrac{a^2\\,x^2 + a^2\\,y^2 - 2a\\,e\\,x^2 + 2a\\,e\\,y^2 - e\\,x^2 - e\\,y^2 + t}{2}", tzz),
            ("-\\left(a - 2e\\right)", "2e - a"), ("a - 2e", "-\\left(2e - a\\right)"),
            ("-\\left(a + 2e\\right)", "-\\left(2e + a\\right)"), ("a + 2e", "2e + a")]
        spec["geodesics"] = [
            "\\ddot{t} + \\dot{t}\\dot{z} - \\left(2e - a\\right)\\dot{x}^2 - \\left(2e + a\\right)\\dot{y}^2"
            " - 2e\\left(x\\dot{x} + y\\dot{y}\\right)\\dot{z} + " + tzz + "\\dot{z}^2 = 0",
            "\\ddot{x} - \\dfrac{a\\,x}{2}\\dot{z}^2 = 0", "\\ddot{y} + \\dfrac{a\\,y}{2}\\dot{z}^2 = 0",
            "\\ddot{z} - \\dfrac{1}{2}\\dot{z}^2 = 0"]
    else:
        spec["rewrite"] = [("\\left(x + y\\right)\\left(x - y\\right)", "\\left(x^2 - y^2\\right)")]
    return spec


def ori_core_check(chart):
    """Ori's equations (2) and (3) and his footnote 1: det g = -1, R_izjz = -f,ij/2 with every
    other component zero up to the symmetries, and a Ricci tensor whose one component is
    R_zz = -(f,xx + f,yy)/2, so the core is a vacuum exactly where f is harmonic in x and y."""
    T, x, y, z = chart.symbols
    f = chart.reader.parameters["f"]
    if sp.simplify(chart.geo.g.det() + 1) != 0:
        raise AssertionError("ori_time_machine: det g is not -1 in the vacuum core chart")
    R = chart.geo.riemann_llll()
    plane = {1: x, 2: y}
    for i, j, k, l in itertools.product(range(4), repeat=4):
        want = 0
        if j == 3 and l == 3 and i in plane and k in plane:
            want = -sp.diff(f, plane[i], plane[k]) / 2
        elif i == 3 and k == 3 and j in plane and l in plane:
            want = -sp.diff(f, plane[j], plane[l]) / 2
        elif i == 3 and l == 3 and j in plane and k in plane:
            want = sp.diff(f, plane[j], plane[k]) / 2
        elif j == 3 and k == 3 and i in plane and l in plane:
            want = sp.diff(f, plane[i], plane[l]) / 2
        if sp.simplify(R[i][j][k][l] - want) != 0:
            raise AssertionError(f"ori_time_machine: R_{i}{j}{k}{l} misses Ori's equation (3)")
    ricci = chart.geo.ricci_ll()
    for i, j in itertools.product(range(4), repeat=2):
        want = -(sp.diff(f, x, 2) + sp.diff(f, y, 2)) / 2 if i == j == 3 else 0
        if sp.simplify(ricci[i][j] - want) != 0:
            raise AssertionError(f"ori_time_machine: R_{i}{j} is not Ori's null dust term")


def ori_pullback(chart, system):
    """J^T g J, with g the vacuum core chart's metric at f = a(x^2 - y^2)/2 and J the Jacobian
    of the map onto it, against this chart's metric in every slot: T = t + a(x^2 - y^2)/2 -
    e(x^2 + y^2) for the foliation, and T = -uv/2 with z = -2 ln(-u/2) for the Brinkmann chart.
    Then the vacuum, det g = -1, and for the foliation Ori's g^tt and the closed null geodesic
    N, the circle x = y = t = 0, whose tangent the published connection carries into itself."""
    spec = ori_time_machine("vacuum_core")
    source = cp.Chart(spec["system"]["coords"], spec["system"]["parameters"], spec["chart_line_element"])
    a = chart.reader.parameters["a"]
    if system == "foliation":
        e = chart.reader.parameters["e"]
        t, x, y, z = chart.symbols
        image = [t + a * (x ** 2 - y ** 2) / 2 - e * (x ** 2 + y ** 2), x, y, z]
    else:
        u, v, x, y = chart.symbols
        image = [-u * v / 2, x, y, -2 * sp.log(-u / 2)]
    at = dict(zip(source.symbols, image))
    profile = a * (x ** 2 - y ** 2) / 2
    g = source.geo.g.subs(source.reader.parameters["f"], profile).subs(at, simultaneous=True)
    J = sp.Matrix(4, 4, lambda i, j: sp.diff(image[i], chart.symbols[j]))
    pulled = J.T * g * J
    for i in range(4):
        for j in range(i, 4):
            if sp.simplify(pulled[i, j] - chart.geo.g[i, j]) != 0:
                raise AssertionError(f"ori_time_machine: the vacuum core chart pulled back misses the {system} "
                                     f"chart in slot {chart.coords_tex[i]}{chart.coords_tex[j]}")
    if sp.simplify(chart.geo.g.det() + 1) != 0:
        raise AssertionError(f"ori_time_machine: det g is not -1 in the {system} chart")
    ricci = chart.geo.ricci_ll()
    if any(sp.simplify(ricci[i][j]) != 0 for i in range(4) for j in range(4)):
        raise AssertionError(f"ori_time_machine: the {system} chart is not a vacuum")
    if system == "foliation":
        want = t + (2 * e - a) ** 2 * x ** 2 + (2 * e + a) ** 2 * y ** 2 - e * (x ** 2 + y ** 2)
        if sp.simplify(chart.geo.ginv[0, 0] - want) != 0:
            raise AssertionError("ori_time_machine: g^tt is not Ori's")
        on_N = {t: 0, x: 0, y: 0}
        if chart.geo.g[3, 3].subs(on_N) != 0:
            raise AssertionError("ori_time_machine: the circle x = y = t = 0 is not null")
        gamma = chart.geo.christoffel_ull()
        for i in range(3):
            if sp.simplify(gamma[i][3][3].subs(on_N)) != 0:
                raise AssertionError("ori_time_machine: the circle x = y = t = 0 is not a geodesic")


CHARTS["ori_time_machine"] = [lambda s=s: ori_time_machine(s) for s in ORI_CHARTS]

# -- Fisher, Janis, Newman and Winicour's scalar field -----------------------------------

FJNW_SPHERE = "\\left(d\\theta^2 + \\sin^2\\theta\\,d\\phi^2\\right)"
FJNW_ANGLES = ["\\theta \\in [0, \\pi]", "\\phi \\in [0, 2\\pi)"]
FJNW_CHARTS = ["spherical", "jnw", "isotropic", "harmonic"]


def fjnw_hyperbolic(u, m, k):
    """A `pretty` for the harmonic chart, whose values carry e^{2mu}, sinh(ku) and cosh(ku) and
    come back from the checker's Geometry as exponentials: each is rewritten in S = sinh(ku),
    C = cosh(ku) and E = e^{mu}, reduced by C^2 = 1 + S^2, its denominator cleared of C by its
    conjugate, and factored, as chart_printer.hyperbolic does for one argument."""
    S, C, E = sp.symbols("FJNWs FJNWc FJNWe", positive=True)

    def reduce(p):
        return sp.expand(sp.reduced(sp.expand(p), [C ** 2 - 1 - S ** 2], C, S)[1])

    def powers(e):
        a = sp.Poly(sp.expand(e.args[0]), u)
        if a.degree() != 1 or a.coeff_monomial(1) != 0:
            raise AssertionError(f"{e} is no exponential of a multiple of u")
        rate = sp.Poly(a.coeff_monomial(u), m, k)
        n, j = rate.coeff_monomial(k), rate.coeff_monomial(m)
        if not (n.is_Integer and j.is_Integer) or sp.expand(rate.as_expr() - n * k - j * m) != 0:
            raise AssertionError(f"{e} is no power of e^(ku) and e^(mu)")
        return (S + C) ** n * E ** j

    def pretty(value):
        x = sp.sympify(value).replace(sp.sinh, lambda a: (sp.exp(a) - sp.exp(-a)) / 2)
        x = x.replace(sp.cosh, lambda a: (sp.exp(a) + sp.exp(-a)) / 2).replace(lambda e: isinstance(e, sp.exp), powers)
        num, den = (reduce(part) for part in sp.fraction(sp.together(x)))
        d0, d1 = sp.Poly(den, C).coeff_monomial(1), sp.Poly(den, C).coeff_monomial(C)
        if d1 != 0:
            num, den = reduce(num * (d0 - d1 * C)), sp.expand(d0 ** 2 - d1 ** 2 * (1 + S ** 2))
        out = sp.factor(sp.cancel(sp.factor(num) / sp.factor(den)))
        return out.subs({S: sp.sinh(k * u), C: sp.cosh(k * u)}).replace(
            lambda e: e.is_Pow and e.base == E, lambda e: sp.exp(e.exp * m * u)).subs(E, sp.exp(m * u))

    return pretty


def fisher_jnw(system):
    """The static, spherically symmetric, asymptotically flat solution of Einstein's equations
    with a massless scalar field, in four charts. Wyman's, in the form Virbhadra gives it, his
    (11), has g_tt g_rr = -1 and is Schwarzschild's at gamma = 1, with b = 2 sqrt(m^2 + q^2) and
    gamma = 2m/b. Janis, Newman and Winicour's own radius is R = r - (1 + gamma) b/2, Virbhadra's
    (14), with their r_0 = gamma b and mu = 1/gamma; it is Fisher's Z/c, his (14) and (25), and in
    it Gamma^R_theta theta = -R. The isotropic radius has r = rho (1 + b/4 rho)^2, the chart
    Fisher names in his (35) and (36) and Xanthopoulos and Zannias integrate in. Bronnikov's
    harmonic coordinate u has e^{-2ku} = 1 - b/r with k = b/2 and m = gamma k, in which the
    scalar field is linear in u; it is (54) of Bronnikov, Fabris and Zhidenko. The first three
    name the ratio f or h and are printed around its powers by named_powers; the fourth is
    printed in sinh(ku), cosh(ku) and e^{mu} by fjnw_hyperbolic. fisher_jnw_check holds each to
    the field equations and each after the first to being the first pulled back, and
    fisher_jnw.md beside this file is the derivation."""
    state = {}
    named, wholes, kretschmann = {}, [], None
    if system == "harmonic":
        coords, name = ["t", "u", "\\theta", "\\phi"], "Harmonic"
        parameters = ["m", "k"]
        domains = (["t \\in (-\\infty, \\infty)", "u \\in (0, \\infty)"] + FJNW_ANGLES
                   + ["u \\to 0 \\;\\text{(spatial infinity)}", "u \\to \\infty \\;\\text{(the singularity, for } m < k\\text{)}"])

        def line(c2):
            return (f"ds^2 = -e^{{-2mu}}{c2}dt^2 + \\dfrac{{k^2e^{{2mu}}}}{{\\sinh^2(ku)}}\\left(\\dfrac{{k^2du^2}}{{\\sinh^2(ku)}}"
                    " + d\\theta^2 + \\sin^2\\theta\\,d\\phi^2\\right)")
        probe = vm.Reader(coords, parameters, ())
        x, m, k = probe.symbol["u"], probe.parameters["m"], probe.parameters["k"]
        printer = {"lead": [k, m], "factors": [k, m], "flip": False}
        pretty = fjnw_hyperbolic(x, m, k)
        components = {}
        kretschmann = ("\\dfrac{4e^{-4mu}\\sinh^6(ku)\\left(\\left(3k^4 + 14k^2m^2 + 7m^4\\right)\\sinh^2(ku)"
                       " - 8km\\left(k^2 + 2m^2\\right)\\sinh(ku)\\cosh(ku) + 12k^2m^2\\right)}{k^8}")
    else:
        if system == "spherical":
            coords, name = ["t", "r", "\\theta", "\\phi"], "Spherical"
            parameters = ["b", "\\gamma", "f = 1 - \\dfrac{b}{r}"]
            domains = (["t \\in (-\\infty, \\infty)", "r \\in (b, \\infty)"] + FJNW_ANGLES
                       + ["r = b \\;\\text{(the singularity, for } \\gamma < 1\\text{)}"])

            def line(c2):
                return f"ds^2 = -f^{{\\gamma}}{c2}dt^2 + \\dfrac{{dr^2}}{{f^{{\\gamma}}}} + f^{{1-\\gamma}}r^2" + FJNW_SPHERE
            probe = vm.Reader(coords, parameters, ())
            x, b, g = probe.symbol["r"], probe.parameters["b"], probe.parameters["gamma"]
            names = {"f": {x - b: 1, x: -1}}
            sums = [(2 * x - (1 + g) * b, "2r - b\\left(1 + \\gamma\\right)")]
            area = "r^2f^{1-\\gamma}"
            ratio = "f"
        elif system == "jnw":
            coords, name = ["t", "R", "\\theta", "\\phi"], "Janis-Newman-Winicour"
            parameters = ["b", "\\gamma", "f = \\dfrac{2R - b(1-\\gamma)}{2R + b(1+\\gamma)}"]
            domains = (["t \\in (-\\infty, \\infty)", "R \\in (b(1-\\gamma)/2, \\infty)"] + FJNW_ANGLES
                       + ["R = b(1-\\gamma)/2 \\;\\text{(the singularity, for } \\gamma < 1\\text{)}"])

            def line(c2):
                return (f"ds^2 = -f^{{\\gamma}}{c2}dt^2 + \\dfrac{{dR^2}}{{f^{{\\gamma}}}}"
                        " + f^{1-\\gamma}\\left(R + \\dfrac{b(1+\\gamma)}{2}\\right)^2" + FJNW_SPHERE)
            probe = vm.Reader(coords, parameters, ())
            x, b, g = probe.symbol["R"], probe.parameters["b"], probe.parameters["gamma"]
            names = {"f": {2 * x - (1 - g) * b: 1, 2 * x + (1 + g) * b: -1}}
            sums = [(2 * x + (1 + g) * b, "2R + b\\left(1 + \\gamma\\right)"),
                    (2 * x - (1 - g) * b, "2R - b\\left(1 - \\gamma\\right)")]
            area = "\\left(R + \\dfrac{b(1+\\gamma)}{2}\\right)^2f^{1-\\gamma}"
            ratio = "f"
        else:
            coords, name = ["t", "\\rho", "\\theta", "\\phi"], "Isotropic"
            parameters = ["b", "\\gamma", "h = \\dfrac{4\\rho - b}{4\\rho + b}"]
            domains = (["t \\in (-\\infty, \\infty)", "\\rho \\in (b/4, \\infty)"] + FJNW_ANGLES
                       + ["\\rho = b/4 \\;\\text{(the singularity, for } \\gamma < 1\\text{)}"])

            def line(c2):
                return (f"ds^2 = -h^{{2\\gamma}}{c2}dt^2 + \\left(1 + \\dfrac{{b}}{{4\\rho}}\\right)^4h^{{2(1-\\gamma)}}"
                        "\\left(d\\rho^2 + \\rho^2" + FJNW_SPHERE + "\\right)")
            probe = vm.Reader(coords, parameters, ())
            x, b, g = probe.symbol["\\rho"], probe.parameters["b"], probe.parameters["gamma"]
            # The reader holds h^(2 gamma) as a power of b - 4 rho, so that is the base the powers are counted on.
            names = {"h": {b - 4 * x: 1, 4 * x + b: -1}}
            sums = []
            area = None
            ratio = "h"
        for i, (poly, text) in enumerate(sums):
            placeholder = sp.Symbol(f"FJNWsum{i}")
            named[placeholder] = text
            wholes.append((placeholder, poly))
        printer = {"lead": [x, b, g], "factors": [b, g, x], "rising": [g], "flip": False, "named": named,
                   "collect": lambda poly, pr: cp.collect_by(poly, [x], pr)}
        powers = named_powers(names, state)

        def pretty(value):
            # A base written with the other sign under a power that holds gamma is the same positive
            # ratio, and each sum the chart names is written the way the chart writes it.
            value = sp.factor(sp.sympify(value)).replace(
                lambda e: e.is_Pow and e.base == -1 and not e.exp.is_Number, lambda e: sp.Integer(1))
            out = sp.Integer(1)
            for factor in sp.Mul.make_args(powers(value)):
                base, e = (factor.base, factor.exp) if factor.is_Pow else (factor, sp.Integer(1))
                if base.is_Add and e.is_Integer:
                    # sympy carries a rational coefficient into a lone sum, so the sum is matched by
                    # its primitive part and the coefficient put back in front.
                    content, primitive = base.as_content_primitive()
                    for placeholder, poly in wholes:
                        if sp.expand(primitive - poly) == 0:
                            base, out = placeholder, out * content ** e
                        elif sp.expand(primitive + poly) == 0:
                            base, out = placeholder, out * (-content) ** e
                out *= base ** e
            return out

        if system == "isotropic":
            conformal = "\\left(1 + \\dfrac{b}{4\\rho}\\right)^4h^{2(1-\\gamma)}"
            components = {"metric_components": {("t", "t"): "-h^{2\\gamma}", ("\\rho", "\\rho"): conformal,
                                                ("\\theta", "\\theta"): "\\rho^2" + conformal,
                                                ("\\phi", "\\phi"): "\\rho^2" + conformal + "\\sin^2\\theta"},
                          "inverse_metric_components": {("t", "t"): "-\\dfrac{1}{h^{2\\gamma}}"}}
        else:
            components = {
                "metric_components": {("t", "t"): "-f^{\\gamma}", (coords[1], coords[1]): "\\dfrac{1}{f^{\\gamma}}",
                                      ("\\theta", "\\theta"): area, ("\\phi", "\\phi"): area + "\\sin^2\\theta"},
                "inverse_metric_components": {("t", "t"): "-\\dfrac{1}{f^{\\gamma}}", (coords[1], coords[1]): "f^{\\gamma}"}}

    def check(chart):
        state["printer"] = chart.printer
        fisher_jnw_check(chart, system)

    spec = {
        "metric_id": "fisher_jnw",
        "system": {"id": system, "name": name, "coords": coords, "domains": domains, "parameters": parameters,
                   "line_element": line("c^2")},
        "chart_line_element": line(""),
        "printer": printer,
        "pretty": pretty,
        "components": components,
        "check": check,
    }
    if kretschmann:
        spec["kretschmann"] = kretschmann
    return spec


def fisher_jnw_check(chart, system):
    """In every chart the Ricci tensor is 2 d phi d phi for the field phi = sqrt(1 - gamma^2) ln(f)/2,
    which is -sqrt(k^2 - m^2) u in the harmonic chart, and that field solves the wave equation;
    without the field, gamma = 1 or m = k, the Ricci tensor vanishes; the harmonic coordinate is
    harmonic; and each chart after the first is the first pulled back: r = R + (1 + gamma) b/2,
    r = rho (1 + b/4 rho)^2, and r = 2k/(1 - e^{-2ku}) with b = 2k and gamma = m/k."""
    t, x, th, ph = chart.symbols
    P = chart.reader.parameters
    if system == "harmonic":
        m, k = P["m"], P["k"]
        b, gamma = 2 * k, m / k
        radius = 2 * k / (1 - sp.exp(-2 * k * x))
        vacuum = {m: k}
    else:
        b, gamma = P["b"], P["gamma"]
        radius = {"spherical": x, "jnw": x + (1 + gamma) * b / 2, "isotropic": x * (1 + b / (4 * x)) ** 2}[system]
        vacuum = {gamma: 1}
    field = sp.sqrt(1 - gamma ** 2) * sp.log(1 - b / radius) / 2
    g, ginv = chart.geo.g, chart.geo.ginv
    ricci = chart.geo.ricci_ll()
    for a, c in vm._indices(4, 2):
        source = 2 * sp.diff(field, chart.symbols[a]) * sp.diff(field, chart.symbols[c])
        if vm.norm(vm._at(ricci, (a, c)) - source) != 0:
            raise AssertionError(f"fisher_jnw: R_ab is not 2 d_a phi d_b phi in slot {(a, c)} of the {system} chart")
        if vm.norm(vm.norm(vm._at(ricci, (a, c))).subs(vacuum)) != 0:
            raise AssertionError(f"fisher_jnw: the {system} chart is not vacuum without the field")
    # d_x(sqrt(-g) g^xx d_x phi) = 0, written through logarithmic derivatives so that no root is taken.
    density = sp.diff(g.det(), x) / (2 * g.det()) + sp.diff(ginv[1, 1], x) / ginv[1, 1]
    if vm.norm(density + sp.diff(field, x, 2) / sp.diff(field, x)) != 0:
        raise AssertionError(f"fisher_jnw: the field does not solve the wave equation in the {system} chart")
    if system == "harmonic" and vm.norm(density) != 0:
        raise AssertionError("fisher_jnw: the coordinate u is not harmonic")
    if system == "spherical":
        return
    source = fisher_jnw("spherical")
    own = cp.Chart(source["system"]["coords"], source["system"]["parameters"], source["chart_line_element"])
    names = dict(zip(own.symbols[2:], chart.symbols[2:]))
    names.update({own.symbols[0]: t, own.reader.parameters["b"]: b, own.reader.parameters["gamma"]: gamma})
    J = sp.diag(1, sp.diff(radius, x), 1, 1)
    pulled = J.T * own.geo.g.subs(names, simultaneous=True).subs(own.symbols[1], radius) * J
    # Compared on positive symbols, the distance X from the singularity among them, so that
    # (h^2)^gamma and h^(2 gamma) are one power of a positive base.
    X, A, B = sp.symbols("FJNWx FJNWa FJNWb", positive=True)
    if system == "harmonic":
        positive = {x: X, P["m"]: A, P["k"]: B}
    else:
        edge = {"jnw": (1 - A) * B / 2, "isotropic": B / 4}[system]
        positive = {x: edge + X, P["gamma"]: A, P["b"]: B}
    for a, c in vm._indices(4, 2):
        if g[a, c] == 0:
            same = pulled[a, c] == 0
        else:
            ratio = (pulled[a, c] / g[a, c]).subs(positive, simultaneous=True)
            if system == "harmonic":
                ratio = ratio.rewrite(sp.exp)
            ratio = sp.powdenest(sp.factor(ratio), force=True)
            # The reader holds h^(2 gamma) as a power of b - 4 rho with (-1)^(2 gamma) beside it: a
            # phase of its generator, and no sign of the component.
            ratio = ratio.replace(lambda e: e.is_Pow and e.base.is_negative and not e.exp.is_Number,
                                  lambda e: (-e.base) ** e.exp)
            same = sp.simplify(sp.powsimp(sp.powdenest(ratio, force=True), force=True)) == 1
        if not same:
            raise AssertionError(f"fisher_jnw: the {system} chart is not the spherical chart pulled back in slot {(a, c)}")


CHARTS["fisher_jnw"] = [lambda s=s: fisher_jnw(s) for s in FJNW_CHARTS]



# -- Van Den Broeck's warp drive -----------------------------------------------------

VDB_CHARTS = ["cartesian", "pocket", "proper_radial"]


def van_den_broeck(system_id):
    """Van Den Broeck's warp drive of 1999, Alcubierre's metric with the flat space of each
    slice multiplied by B^2, in the two charts his paper uses. The Cartesian chart is his
    equation (1), with the shape function f and the expansion factor B left as arbitrary
    functions of the four coordinates, as Alcubierre's f is, so that B = 1 gives Alcubierre's
    chart slot for slot. The comoving spherical chart is the one he computes in: inside the
    wall of the bubble f = 1, x' = x - x_s(t) turns dx - v_s c dt into dx' whatever v_s(t) is,
    and the metric is diag(-1, B^2, B^2, B^2), static and spherically symmetric about the ship.
    That chart is checked to be the Cartesian one pulled back at f = 1, and its G_tt and its
    frame component of Riemann on a sphere against his equations for the energy density and
    for the largest curvature. van_den_broeck.md is the derivation."""
    reals = "(-\\infty, \\infty)"
    if system_id == "cartesian":
        coords = ["t", "x", "y", "z"]
        parameters = ["v_s = v_s(t)", "f = f(t,x,y,z)", "B = B(t,x,y,z)"]
        probe = vm.Reader(coords, parameters, ())
        lead = [probe.parameters["v_s"], probe.parameters["B"], probe.parameters["f"]]
        space = "B^2\\left(\\left(dx - v_s f\\,{0}dt\\right)^2 + dy^2 + dz^2\\right)"
        return {
            "metric_id": "van_den_broeck",
            "system": {"id": "cartesian", "name": "Cartesian", "coords": coords,
                       "domains": [f"{c} \\in {reals}" for c in coords], "parameters": parameters,
                       "line_element": "ds^2 = -c^2dt^2 + " + space.format("c\\,")},
            "chart_line_element": "ds^2 = -dt^2 + " + space.format(""),
            "printer": {"lead": lead, "primed": ["v_s"]},
            # Compared once the tensors are printed: reading Alcubierre's chart first leaves the
            # lowered Riemann tensor of this one printing for tens of minutes where it takes seconds.
            "after": van_den_broeck_alcubierre,
            # sympy does not finish factoring the numerator, so it is written out over 4B^8.
            "kretschmann_text": van_den_broeck_kretschmann,
        }
    if system_id == "proper_radial":
        coords = ["t", "l", "\\theta", "\\phi"]
        parameters = ["r = r(l)", "l_0"]
        probe = vm.Reader(coords, parameters, ())
        space = "dl^2 + r^2\\left(d\\theta^2 + \\sin^2\\theta\\,d\\phi^2\\right)"
        return {
            "metric_id": "van_den_broeck",
            "system": {"id": "proper_radial", "name": "Proper Radial Distance", "coords": coords,
                       "domains": ["t \\in " + reals, "l \\in [-l_0, \\infty)", "\\theta \\in [0, \\pi]",
                                   "\\phi \\in [0, 2\\pi)",
                                   "l = -l_0 \\;\\text{(the centre of the pocket, where } r = 0\\text{)}"],
                       "parameters": parameters, "line_element": "ds^2 = -c^2dt^2 + " + space},
            "chart_line_element": "ds^2 = -dt^2 + " + space,
            "printer": {"lead": [probe.parameters["r"]], "primed": ["r"]},
            "check": van_den_broeck_proper,
            # Two squares: the frame components l theta l theta and theta phi theta phi.
            "kretschmann": ("\\dfrac{8\\left(r''\\right)^2}{r^2} + "
                            "\\dfrac{4\\left(1 - \\left(r'\\right)^2\\right)^2}{r^4}"),
        }
    coords = ["t", "r", "\\theta", "\\phi"]
    parameters = ["B = B(r)", "R"]
    probe = vm.Reader(coords, parameters, ())
    space = "B^2\\left(dr^2 + r^2d\\theta^2 + r^2\\sin^2\\theta\\,d\\phi^2\\right)"
    return {
        "metric_id": "van_den_broeck",
        "system": {"id": "pocket", "name": "Comoving Spherical", "coords": coords,
                   "domains": ["t \\in " + reals, "r \\in [0, R)", "\\theta \\in [0, \\pi]",
                               "\\phi \\in [0, 2\\pi)"],
                   "parameters": parameters, "line_element": "ds^2 = -c^2dt^2 + " + space},
        "chart_line_element": "ds^2 = -dt^2 + " + space,
        "printer": {"lead": [probe.parameters["B"], probe.symbol["r"]], "primed": ["B"]},
        "check": van_den_broeck_pocket,
    }


def van_den_broeck_kretschmann(chart):
    numerator, denominator = sp.fraction(sp.together(chart.geo.kretschmann()))
    return chart.printer(sp.expand(numerator) / sp.factor(denominator))


def van_den_broeck_alcubierre(math, chart):
    """At B = 1 the Cartesian chart is the published Alcubierre metric, slot by slot."""
    B = chart.reader.parameters["B"]
    entry = json.loads((METRICS / "alcubierre.json").read_text(encoding="utf-8"))["coordinates"][0]
    reader = vm.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]], ())
    published = vm.metric_from_line_element(reader, entry["line_element"].replace("c^2dt^2", "dt^2")
                                            .replace("c\\,dt", "dt"), entry["coords"])
    same = {reader.parameters[k]: chart.reader.parameters[k] for k in ("v_s", "f")}
    same.update(zip((reader.symbol[c] for c in entry["coords"]), chart.symbols))
    own = chart.geo.g.subs(B, 1).doit()
    for i in range(4):
        for j in range(i, 4):
            if sp.simplify(published[i, j].subs(same, simultaneous=True) - own[i, j]) != 0:
                raise AssertionError(f"van_den_broeck: at B = 1 the slot {entry['coords'][i]}{entry['coords'][j]} "
                                     "is not Alcubierre's")
    return math


def van_den_broeck_pocket(chart):
    """The comoving spherical chart against the Cartesian one and against Van Den Broeck's
    equations. With f = 1, x - x_s = r cos(theta), y = r sin(theta) cos(phi) and
    z = r sin(theta) sin(phi), J^T g J is the chart's metric in every slot, for any v_s(t).
    His energy density of the Eulerian observers is G_tt/8 pi = (B'^2/B^4 - 2B''/B^3 - 4B'/(B^3 r))/8 pi,
    and his largest frame component of Riemann is R_{theta phi theta phi}/(B^4 r^4 sin^2(theta))
    = B'^2/B^4 - B''/B^3 - B'/(B^3 r) for B a function of r alone... see van_den_broeck.md."""
    spec = van_den_broeck("cartesian")
    source = cp.Chart(spec["system"]["coords"], spec["system"]["parameters"], spec["chart_line_element"])
    t, r, th, ph = chart.symbols
    B = chart.reader.parameters["B"]
    xs = sp.Function("x_s")(t)
    image = [t, xs + r * sp.cos(th), r * sp.sin(th) * sp.cos(ph), r * sp.sin(th) * sp.sin(ph)]
    J = sp.Matrix(4, 4, lambda i, j: sp.diff(image[i], chart.symbols[j]))
    g = source.geo.g
    at = {source.reader.parameters["f"]: 1, source.reader.parameters["B"]: B,
          source.reader.parameters["v_s"]: sp.diff(xs, t)}
    pulled = J.T * g.subs(at, simultaneous=True) * J
    for i in range(4):
        for j in range(i, 4):
            if sp.simplify(pulled[i, j] - chart.geo.g[i, j]) != 0:
                raise AssertionError(f"van_den_broeck: the Cartesian chart pulled back at f = 1 misses the "
                                     f"comoving chart in slot {chart.coords_tex[i]}{chart.coords_tex[j]}")
    B1, B2 = sp.diff(B, r), sp.diff(B, r, 2)
    density = B1 ** 2 / B ** 4 - 2 * B2 / B ** 3 - 4 * B1 / (B ** 3 * r)
    if vm.norm(chart.geo.einstein_ll()[0][0] - density) != 0:
        raise AssertionError("van_den_broeck: G_tt is not Van Den Broeck's energy density times 8 pi")
    riemann = chart.geo.riemann_llll()
    # His frame at a point of the line y = z = 0 has e_1 along x, the radial direction there,
    # and e_2 across it, so his R_1212 is the radial-tangential frame component.
    radial = riemann[1][2][1][2] / (B ** 4 * r ** 2)
    if vm.norm(radial - (B1 ** 2 / B ** 4 - B2 / B ** 3 - B1 / (B ** 3 * r))) != 0:
        raise AssertionError("van_den_broeck: the frame component r theta r theta of Riemann is not his R_1212")


def van_den_broeck_proper(chart):
    """Krasnikov's chart of the pocket against the comoving spherical one and against his
    Einstein tensor. With dl = B d rho and r(l) = B rho, rho the comoving chart's radius, the
    chart pulls back to diag(-1, B^2, B^2 rho^2, B^2 rho^2 sin^2(theta)), and in the static
    frame G_tt = (1 - r'^2 - 2 r r'')/r^2, G_tt + G_ll = -2 r''/r and
    G_tt + G_thetatheta/r^2 = (1 - r'^2 - r r'')/r^2."""
    spec = van_den_broeck("pocket")
    pocket = cp.Chart(spec["system"]["coords"], spec["system"]["parameters"], spec["chart_line_element"])
    rho, B = pocket.symbols[1], pocket.reader.parameters["B"]
    l, r = chart.symbols[1], chart.reader.parameters["r"]
    pulled = chart.geo.g.subs(r, B * rho)
    pulled[1, 1] = pulled[1, 1] * B ** 2
    at = dict(zip(chart.symbols[2:], pocket.symbols[2:]))
    for i in range(4):
        for j in range(i, 4):
            if sp.simplify(pulled[i, j].subs(at) - pocket.geo.g[i, j]) != 0:
                raise AssertionError(f"van_den_broeck: Krasnikov's chart pulled back misses the comoving chart "
                                     f"in slot {chart.coords_tex[i]}{chart.coords_tex[j]}")
    G = chart.geo.einstein_ll()
    r1, r2 = sp.diff(r, l), sp.diff(r, l, 2)
    for value, his in ((G[0][0], (1 - r1 ** 2 - 2 * r * r2) / r ** 2), (G[0][0] + G[1][1], -2 * r2 / r),
                       (G[0][0] + G[2][2] / r ** 2, (1 - r1 ** 2 - r * r2) / r ** 2)):
        if vm.norm(value - his) != 0:
            raise AssertionError("van_den_broeck: the Einstein tensor of the proper distance chart is not Krasnikov's")


CHARTS["van_den_broeck"] = [lambda s=s: van_den_broeck(s) for s in VDB_CHARTS]


# -- Reissner-Nordström-de Sitter ------------------------------------------------------

RNDS_CHARTS = ["static", "eddington_finkelstein_outgoing", "eddington_finkelstein_ingoing", "cosmological"]
RNDS_SPHERE = "\\left(d\\theta^2 + \\sin^2\\theta\\,d\\phi^2\\right)"
RNDS_ANGLES = ["\\theta \\in [0, \\pi]", "\\phi \\in [0, 2\\pi)"]


def reissner_nordstrom_de_sitter(system_id):
    """The charged black hole in de Sitter space, f = 1 - r_s/r + r_q^2/r^2 - Lambda r^2/3, in the
    static chart and the two Eddington-Finkelstein charts built on its tortoise coordinate,
    dr_*/dr = 1/f, and the lukewarm hole r_q = r_s/2 in the cosmological chart of Kastor and
    Traschen, ds^2 = -c^2 dtau^2/W^2 + W^2 (drho^2 + rho^2 dOmega^2) with W = H tau + r_s/(2 rho)
    and Lambda = 3H^2/c^2, whose areal radius is r = H tau rho + r_s/2. Every value of the first
    three is printed around 3r^2 f = 3r^2 - 3r_s r + 3r_q^2 - Lambda r^4, in the order of f itself,
    so that each chart reduces to Kottler's at r_q = 0 and to Reissner and Nordstrom's at
    Lambda = 0 term by term. The Kretschmann scalar is Reissner and Nordstrom's plus de Sitter's
    8 Lambda^2/3, which it is. The cosmological chart is checked, slot by slot, to be the static
    chart pulled back; reissner_nordstrom_de_sitter.md derives the map."""
    if system_id == "cosmological":
        coords = ["\\tau", "\\rho", "\\theta", "\\phi"]
        parameters = ["r_s", "H"]

        def line(c2, rate):
            W = "\\left(" + rate + " + \\dfrac{r_s}{2\\rho}\\right)"
            return ("ds^2 = -" + W + "^{-2}" + c2 + "d\\tau^2 + " + W + "^2\\left(d\\rho^2 + \\rho^2"
                    + RNDS_SPHERE + "\\right)")
        probe = vm.Reader(coords, parameters, ())
        tau, rho = probe.symbol["\\tau"], probe.symbol["\\rho"]
        rs, H = probe.parameters["r_s"], probe.parameters["H"]
        D = sp.Symbol("RD", positive=True)
        return {
            "metric_id": "reissner_nordstrom_de_sitter",
            "system": {"id": system_id, "name": "Cosmological", "coords": coords,
                       "domains": ["\\tau \\in (-r_s/(2H\\rho), \\infty)", "\\rho \\in (0, \\infty)"] + RNDS_ANGLES
                       + ["r = H\\tau\\rho + r_s/2 \\;\\text{(areal radius)}",
                          "2H\\tau\\rho + r_s = 0 \\;\\text{(the singularity)}"],
                       "parameters": parameters, "line_element": line("c^2", "H\\tau")},
            # The symbol tau is x^0 = c tau in the geometry, so the rate that multiplies it is H/c there.
            "chart_line_element": line("", "\\dfrac{H\\tau}{c}"),
            "time": "\\tau",
            "printer": {"lead": [H, D, tau, rho, rs, probe.c],
                        "named": {D: "2H\\tau\\rho + r_s"}},
            "pretty": rnds_cosmological_pretty(2 * H * tau * rho + rs, tau, D),
            "ricci_scalar": "\\dfrac{12H^2}{c^2}",
            # Reissner and Nordstrom's scalar at r_q = r_s/2 and the areal radius, plus de Sitter's.
            "kretschmann": ("\\dfrac{128r_s^2\\left(24H^2\\tau^2\\rho^2 + r_s^2\\right)}"
                            "{\\left(2H\\tau\\rho + r_s\\right)^8} + \\dfrac{24H^4}{c^4}"),
            "check": reissner_nordstrom_de_sitter_pullback,
        }
    f = "\\left(1 - \\dfrac{r_s}{r} + \\dfrac{r_q^2}{r^2} - \\dfrac{\\Lambda r^2}{3}\\right)"
    bare = f[6:-7]
    sphere = " + r^2" + RNDS_SPHERE
    domains = ["r \\in (0, \\infty)"] + RNDS_ANGLES + [
        "r = r_- \\;\\text{(inner horizon)}", "r = r_+ \\;\\text{(black hole horizon)}",
        "r = r_c \\;\\text{(cosmological horizon)}"]
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
    parameters = ["r_s", "r_q", "\\Lambda"]
    probe = vm.Reader(coords, parameters, ())
    r, rs, rq, L = probe.symbol["r"], probe.parameters["r_s"], probe.parameters["r_q"], probe.parameters["Lambda"]
    return {
        "metric_id": "reissner_nordstrom_de_sitter",
        "system": {"id": system_id, "name": name, "coords": coords,
                   "domains": [coords[0] + " \\in (-\\infty, \\infty)"] + domains,
                   "parameters": parameters, "line_element": line},
        "chart_line_element": chart_line,
        "printer": {"rising": [L, rq, rs], "lead": [L, r, rs, rq], "flip": False},
        "components": {"metric_components": metric, "inverse_metric_components": inverse},
        "ricci_scalar": "4\\Lambda",
        "kretschmann": ("\\dfrac{12r_s^2}{r^6} - \\dfrac{48r_s\\,r_q^2}{r^7} + \\dfrac{56r_q^4}{r^8}"
                        " + \\dfrac{8\\Lambda^2}{3}"),
    }


def rnds_cosmological_pretty(twice, tau, D):
    """A pretty printer for the cosmological chart, whose values are built on twice the areal
    radius, 2 H tau rho + r_s: each value is factored, that sum is written as the placeholder D
    wherever it is a factor, and any other sum is written as a polynomial in D where that
    leaves it fewer terms, as 3H^2 D^4 + 4 r_s^2 c^2 in the Einstein tensor."""
    back = sp.solve(twice - D, tau)[0]

    def pretty(value):
        out = sp.Integer(1)
        for factor in sp.Mul.make_args(sp.factor(value)):
            base, k = (factor.base, factor.exp) if factor.is_Pow else (factor, sp.Integer(1))
            if base.is_Add:
                if sp.expand(base - twice) == 0:
                    base = D
                elif sp.expand(base + twice) == 0:
                    base, out = D, out * (-1) ** k
                else:
                    inside = sp.expand(base.subs(tau, back))
                    if not sp.fraction(sp.together(inside))[1].free_symbols and len(inside.args) < len(base.args):
                        base = inside
            out *= base ** k
        return out

    return pretty


def reissner_nordstrom_de_sitter_pullback(chart):
    """J^T g J, with g the static chart's metric at r_q = r_s/2 and Lambda = 3H^2/c^2, on its plane
    of the time T and r, and J the Jacobian of r = H tau rho + r_s/2 and
    cT = (c/H) ln(H tau) + G(r), where dG/dr = H r^2/(c (r - r_s/2) f), against the cosmological
    chart's metric in every slot, and the two spheres against one another."""
    spec = reissner_nordstrom_de_sitter("static")
    static = cp.Chart(spec["system"]["coords"], spec["system"]["parameters"], spec["chart_line_element"])
    r = static.symbols[1]
    tau, rho = chart.symbols[:2]                                 # the symbol tau is x^0 = c tau here
    rs, H, c = chart.reader.parameters["r_s"], chart.reader.parameters["H"], chart.reader.c
    h = H / c
    radius = h * tau * rho + rs / 2
    at = {static.reader.parameters["r_s"]: rs, static.reader.parameters["r_q"]: rs / 2,
          static.reader.parameters["Lambda"]: 3 * h ** 2}
    f = (-static.geo.g[0, 0]).subs(at)
    slope = (h * r ** 2 / ((r - rs / 2) * f)).subs(r, radius)    # dT/dr at fixed tau, in the chart x^0 = cT
    J = sp.Matrix([[1 / (h * tau) + slope * sp.diff(radius, tau), slope * sp.diff(radius, rho)],
                   [sp.diff(radius, tau), sp.diff(radius, rho)]])
    plane = static.geo.g[:2, :2].subs(at).subs(r, radius)
    pulled = J.T * plane * J
    for i in range(2):
        for j in range(i, 2):
            if vm.norm(pulled[i, j] - chart.geo.g[i, j]) != 0:
                raise AssertionError(f"reissner_nordstrom_de_sitter: the static chart pulled back misses the "
                                     f"cosmological chart in slot {chart.coords_tex[i]}{chart.coords_tex[j]}")
    for k in (2, 3):
        if vm.norm(static.geo.g[k, k].subs(r, radius).subs(static.symbols[2], chart.symbols[2])
                   - chart.geo.g[k, k]) != 0:
            raise AssertionError(f"reissner_nordstrom_de_sitter: the static sphere misses the cosmological one in slot {k}")


CHARTS["reissner_nordstrom_de_sitter"] = [lambda s=s: reissner_nordstrom_de_sitter(s) for s in RNDS_CHARTS]


# -- The travelling wave on a cosmic string ----------------------------------------------

STRING_WAVE_CHARTS = ("null_conical", "isotropic", "moving_string")
STRING_WAVE_FACTOR = "\\left(\\dfrac{\\rho}{\\ell}\\right)^{2b - 2}"


def string_wave(system):
    """Garfinkle and Vachaspati's travelling wave on a cosmic string, with b = 1 - 4G mu/c^2 and
    the null coordinates u = ct - z and v = ct + z, both lengths, in three charts:

    null_conical   the pp-wave on the cone, -du dv + F du^2 + dr^2 + b^2 r^2 dphi^2 with F(u, r, phi)
                   free, Frolov and Garfinkle's family and the form of Garriga and Peter's (14)\\;
    isotropic      Garfinkle's own chart, as Anderson's (1.5) and (1.6) write it, the cone
                   conformally flat, (rho/ell)^(2b - 2)(dx^2 + dy^2), and the wave the one term
                   -2(x A'' + y B'') du^2\\;
    moving_string  the chart of Anderson's (1.10) and of Dabholkar, Gauntlett, Harvey and Waldram's
                   (2.14), in which the string lies along X = A(u), Y = B(u).

    The isotropic chart is checked to be the null conical chart pulled back at the profile
    F = -2 ell (b r/ell)^(1/b)(A'' cos phi + B'' sin phi), and the moving string chart to be the
    isotropic chart pulled back\\; the isotropic chart's Ricci tensor is checked to vanish.
    string_wave.md is the derivation."""
    reals = "(-\\infty, \\infty)"
    state = {}
    probe_parameters = ["b", "\\ell", "A = A(u)", "B = B(u)"]
    if system == "null_conical":
        coords, parameters = ["u", "v", "r", "\\phi"], ["b", "F = F(u,r,\\phi)"]
        probe = vm.Reader(coords, parameters, ())
        return {
            "metric_id": "string_wave",
            "system": {"id": system, "name": "Null Conical", "coords": coords, "parameters": parameters,
                       "domains": ["u \\in " + reals, "v \\in " + reals, "r \\in (0, \\infty)", "\\phi \\in [0, 2\\pi)",
                                   "r = 0 \\;\\text{(the string)}"],
                       "line_element": "ds^2 = -du\\,dv + F\\,du^2 + dr^2 + b^2r^2d\\phi^2"},
            "chart_line_element": "ds^2 = -du\\,dv + F\\,du^2 + dr^2 + b^2r^2d\\phi^2",
            "printer": {"lead": [probe.parameters["b"], probe.symbol["r"]]},
        }
    if system == "isotropic":
        coords = ["u", "v", "x", "y"]
        parameters = probe_parameters + ["\\rho = \\sqrt{x^2 + y^2}"]
        line = ("ds^2 = -du\\,dv - 2\\left(x\\,A'' + y\\,B''\\right)du^2 + " + STRING_WAVE_FACTOR
                + "\\left(dx^2 + dy^2\\right)")
        name, string = "Isotropic", "(x, y) = (0, 0) \\;\\text{(the string)}"
        components = {"metric_components": {("u", "u"): "-2\\left(x\\,A'' + y\\,B''\\right)"}}
    else:
        coords = ["u", "V", "X", "Y"]
        parameters = probe_parameters + ["\\rho = \\sqrt{(X - A)^2 + (Y - B)^2}"]
        line = ("ds^2 = -du\\,dV + dX^2 + dY^2 + \\left(" + STRING_WAVE_FACTOR + " - 1\\right)"
                "\\left(\\left(dX - A'\\,du\\right)^2 + \\left(dY - B'\\,du\\right)^2\\right)")
        name, string = "Moving String", "(X, Y) = (A, B) \\;\\text{(the string)}"
        # The metric as the line element writes it, around (rho/ell)^(2b - 2) - 1.
        lift = "\\left(" + STRING_WAVE_FACTOR + " - 1\\right)"
        components = {"metric_components": {
            ("u", "u"): lift + "\\left(\\left(A'\\right)^2 + \\left(B'\\right)^2\\right)",
            **{pair: "-" + lift + f"{f}'" for f, c in (("A", "X"), ("B", "Y")) for pair in (("u", c), (c, "u"))}}}
    probe = vm.Reader(coords, parameters, ())
    first, second = probe.symbol[coords[2]], probe.symbol[coords[3]]
    b, ell, A, B = (probe.parameters[n] for n in ("b", "ell", "A", "B"))
    moving = system == "moving_string"
    # Across the string every value is written in the two distances from it, x and y, or X - A
    # and Y - B, which the moving string chart prints as the differences they are.
    xi, eta = (sp.Symbol("_xi", real=True), sp.Symbol("_eta", real=True)) if moving else (first, second)
    rho = sp.Symbol("rho", positive=True)
    powers = {}

    def pretty(value):
        """The value with every power of (x^2 + y^2)/ell^2 whose exponent holds b written as a
        power of (rho/ell)^(2b - 2), the conformal factor of the cone, and then factored, with
        every factor x^2 + y^2 written rho^2."""
        L, S, T, W = (sp.Dummy(positive=True) for _ in range(4))
        value = sp.sympify(value)
        if moving:
            value = state["reader"].surface(value)
            value = value.subs({first: xi + A, second: eta + B}, simultaneous=True).doit()
        value = value.subs(ell, L)

        def split(power):
            # The base is (x^2 + y^2)^p ell^q; with x^2 + y^2 = T ell^2 the power is T^(p e) ell^((2p + q) e),
            # and the part of ell's exponent that holds b cancels between the factors of a term.
            base = sp.powsimp(sp.expand_power_base(sp.factor(power.base).subs(xi ** 2 + eta ** 2, S), force=True))
            found = base.as_powers_dict()
            p, q = found.get(S, 0), found.get(L, 0)
            if sp.simplify(base - S ** p * L ** q) != 0:
                raise AssertionError(f"string_wave: {power} is no power of the distance from the string")
            return T ** sp.expand(p * power.exp) * L ** sp.expand((2 * p + q) * power.exp)

        value = value.replace(lambda e: e.is_Pow and e.exp.has(b), split)
        value = sp.factor(sp.powsimp(sp.together(value), force=True))

        def conformal(power):
            # T^(k b + n) is W^k (T)^(n + k), with W = T^(b - 1) the conformal factor.
            whole, symbolic = sp.expand(power.exp).as_coeff_Add()
            k = sp.cancel(symbolic / b)
            if not k.is_Integer:
                raise AssertionError(f"string_wave: {power} is no power of the conformal factor")
            return W ** k * ((xi ** 2 + eta ** 2) / L ** 2) ** (whole + k)

        value = value.replace(lambda e: e.is_Pow and e.base == T and e.exp.has(b), conformal)
        value = value.subs(T, (xi ** 2 + eta ** 2) / L ** 2)
        if any(e.is_Pow and e.exp.has(b) for e in sp.preorder_traversal(value)):
            raise AssertionError(f"string_wave: {value} keeps a power of b")

        def named(k):
            name = powers.setdefault(k, sp.Symbol(f"STRINGWAVEPOWER{len(powers)}", positive=True))
            exponent = state["printer"].positive_first(sp.expand(2 * k * (b - 1)))
            state["printer"].overrides[name] = "\\left(\\dfrac{\\rho}{\\ell}\\right)^{" + exponent + "}"
            return name

        def finished(part):
            # Factored, with the factors that hold the conformal factor multiplied out into one
            # bracket, as (rho/ell)^(2 - 2b) - 1, and every factor x^2 + y^2 written rho^2.
            part = sp.factor(part)
            holding = [f for f in sp.Mul.make_args(part) if f.has(W)]
            rest = sp.cancel(part / sp.Mul(*holding))
            bracket = sp.expand(sp.Mul(*holding))
            bracket = bracket.replace(lambda e: e.is_Pow and e.base == W, lambda e: named(e.exp))
            bracket = bracket.subs(W, named(sp.Integer(1))) if bracket.has(W) else bracket
            return sp.factor(rest).subs(xi ** 2 + eta ** 2, rho ** 2).subs(L, ell) * bracket

        # A value that holds the string's acceleration beside terms without it is written as the
        # two parts it is made of, the cone's own and the wave's, each factored.
        accelerations = {sp.diff(A, state["u"], 2): 0, sp.diff(B, state["u"], 2): 0}
        still = value.xreplace(accelerations)
        if still == 0 or sp.expand(value - still) == 0:
            return finished(value)
        return sp.Add(finished(still), finished(sp.together(value - still)), evaluate=False)

    def check(chart):
        state["printer"], state["reader"], state["u"] = chart.printer, chart.reader, chart.symbols[0]
        if moving:
            chart.printer.named.update({xi: "X - A", eta: "Y - B"})
        string_wave_check(chart, system)

    return {
        "metric_id": "string_wave",
        "system": {"id": system, "name": name, "coords": coords, "parameters": parameters,
                   "domains": [f"{c} \\in " + reals for c in coords] + [string], "line_element": line},
        "chart_line_element": line,
        "printer": {"lead": [b, rho, xi, eta, first, second, ell], "primed": ["A", "B"]},
        "pretty": pretty,
        # The moving string chart's rho holds A and B, so the checker holds it as a function, and
        # every value is written out in X - A and Y - B before it is printed or compared.
        **({"reduce": lambda value: vm.norm(state["reader"].surface(sp.sympify(value)))} if moving else {}),
        "components": components,
        "check": check,
    }


def string_wave_chart(system):
    spec = string_wave(system)
    return cp.Chart(spec["system"]["coords"], spec["system"]["parameters"], spec["chart_line_element"])


def string_wave_check(chart, system):
    """The isotropic chart against the null conical chart at the travelling wave's own profile,
    x = ell (b r/ell)^(1/b) cos(phi) and y = ell (b r/ell)^(1/b) sin(phi), checked at random points
    since the map carries the power 1/b, and its Ricci tensor against zero\\; the moving string
    chart against the isotropic chart, x = X - A, y = Y - B and
    v = V - 2 A' (X - A) - 2 B' (Y - B) - int (A'^2 + B'^2) du, checked exactly."""
    u, v, x, y = chart.symbols
    p = chart.reader.parameters
    b, ell, A, B = (p[n] for n in ("b", "ell", "A", "B"))
    surface = chart.reader.surface
    if system == "isotropic":
        ricci = chart.geo.ricci_ll()
        if any(vm.norm(surface(ricci[i][j])) != 0 for i in range(4) for j in range(4)):
            raise AssertionError("string_wave: the isotropic chart's Ricci tensor does not vanish")
        source = string_wave_chart("null_conical")
        r, phi = source.symbols[2:]
        sb, F = source.reader.parameters["b"], source.reader.parameters["F"]
        radius = sp.sqrt(x ** 2 + y ** 2)
        image = [u, v, ell * (radius / ell) ** b / b, sp.atan2(y, x)]
        profile = -2 * (x * sp.diff(A, u, 2) + y * sp.diff(B, u, 2))
        g = source.geo.g.subs({F: profile, sb: b}).subs(dict(zip(source.symbols, image)), simultaneous=True)
        J = sp.Matrix(4, 4, lambda i, j: sp.diff(image[i], chart.symbols[j]))
        difference = J.T * g * J - chart.geo.g.applyfunc(surface)
        # The profile is harmonic on the cone: the null conical chart's only Ricci component.
        lap = (sp.diff(F, r, 2) + sp.diff(F, r) / r + sp.diff(F, phi, 2) / (sb ** 2 * r ** 2))
        if vm.norm(source.geo.ricci_ll()[0][0] + lap / 2) != 0:
            raise AssertionError("string_wave: the null conical chart's R_uu is not minus half the Laplacian of F")
        wave = -2 * ell * (sb * r / ell) ** (1 / sb) * (sp.Symbol("a2") * sp.cos(phi) + sp.Symbol("b2") * sp.sin(phi))
        on_cone = sp.diff(wave, r, 2) + sp.diff(wave, r) / r + sp.diff(wave, phi, 2) / (sb ** 2 * r ** 2)
        if sp.simplify(on_cone) != 0:
            raise AssertionError("string_wave: the travelling wave's profile is not harmonic on the cone")
    else:
        source = string_wave_chart("isotropic")
        sp_ = source.reader.parameters
        dA, dB = sp.diff(A, u), sp.diff(B, u)
        image = [u, v - 2 * dA * (x - A) - 2 * dB * (y - B) - sp.Integral(dA ** 2 + dB ** 2, u), x - A, y - B]
        g = source.geo.g.applyfunc(source.reader.surface)
        g = g.subs({sp_["A"]: A, sp_["B"]: B, sp_["b"]: b, sp_["ell"]: ell})
        g = g.subs(dict(zip(source.symbols, image)), simultaneous=True).doit()
        J = sp.Matrix(4, 4, lambda i, j: sp.diff(image[i], chart.symbols[j]))
        difference = J.T * g * J - chart.geo.g.applyfunc(surface)
    rng = random.Random(0)
    for _ in range(6):
        atoms = set()
        for entry in difference:
            atoms |= entry.atoms(sp.Derivative) | entry.atoms(sp.core.function.AppliedUndef) | entry.free_symbols
        at = {a: sp.Rational(rng.randint(10 ** 3, 2 * 10 ** 3), 10 ** 3)
              for a in sorted(atoms, key=lambda a: (-len(str(a)), str(a)))}
        for i in range(4):
            for j in range(4):
                if abs(complex(difference[i, j].xreplace(at).evalf(30))) > 1e-20:
                    raise AssertionError(f"string_wave: the {system} chart is not its source pulled back in slot "
                                         f"{chart.coords_tex[i]}{chart.coords_tex[j]}")


CHARTS["string_wave"] = [lambda s=s: string_wave(s) for s in STRING_WAVE_CHARTS]


# -- Bardeen ---------------------------------------------------------------------------

def bardeen_pretty(r, g):
    """A pretty printer for Bardeen's charts, whose every value is a rational function of r,
    g^2, r_s and one radical, Q = sqrt(r^2 + g^2). Written in r, Q and r_s, with g^2 = Q^2 - r^2,
    no relation is left among the generators, so the value factors there, and each factor is
    then written back: a factor even in Q as a polynomial in r and g, factored, and any other
    with each odd power of Q as that power of r^2 + g^2, so that the zero of f stands in a
    denominator as (r^2 + g^2)^{3/2} - r_s r^2 and no radical is ever rationalised away. Two
    factors that Q -> -Q exchanges are multiplied out first, since r^2 - Q^2 is -g^2."""
    Q = sp.Symbol("Q", positive=True)
    square = r ** 2 + g ** 2

    def back(poly):
        terms = sp.Poly(poly, Q).terms()
        if all(k % 2 == 0 for (k,), _ in terms):
            return sp.factor(sp.expand(sum(c * square ** (k // 2) for (k,), c in terms)))
        return sum(c * (sp.expand(square ** (k // 2)) if k % 2 == 0 else square ** sp.Rational(k, 2))
                   for (k,), c in terms)

    def pretty(value):
        value = sp.sympify(value)
        if value == 0:
            return value
        value = value.replace(lambda p: p.is_Pow and sp.expand(p.base - square) == 0, lambda p: Q ** (2 * p.exp))
        value = sp.together(value).subs(g, sp.sqrt(Q ** 2 - r ** 2))
        if value.has(g) or any(not p.exp.is_Integer for p in value.atoms(sp.Pow)):
            raise AssertionError(f"bardeen: an odd power of g is left in {value}")
        powers = {}
        for f in sp.Mul.make_args(sp.factor(value)):
            base, k = (f.base, f.exp) if f.is_Pow else (f, sp.Integer(1))
            powers[base] = powers.get(base, 0) + k
        # A factor and its image under Q -> -Q, as r - Q and r + Q, are one factor even in Q.
        for base in list(powers):
            mirror = sp.expand(base.subs(Q, -Q))
            other = next((b for b in powers if b != base and sp.expand(b - mirror) * sp.expand(b + mirror) == 0), None)
            if base != Q and base.has(Q) and other is not None and powers.get(base) == powers.get(other):
                k = powers.pop(base)
                del powers[other]
                powers[sp.expand(base * other)] = k
        out = sp.Integer(1)
        for base, k in powers.items():
            out *= (square ** (k / 2) if base == Q else back(base) ** k if base.has(Q) else base ** k)
        return out

    return pretty


def bardeen(system_id):
    """Bardeen's regular black hole, f = 1 - r_s r^2/(r^2 + g^2)^{3/2}, in the static chart and
    in the two Eddington-Finkelstein charts built on its tortoise coordinate, dr_*/dr = 1/f.
    The parameters are Schwarzschild's r_s and the length g, the charge of Ayon-Beato and
    Garcia's magnetic monopole, so each chart reduces to Schwarzschild's at g = 0. Every value
    is printed with the radical kept as a power of r^2 + g^2, by bardeen_pretty, the metric and
    its inverse are written as the line element writes f, and bardeen_source checks the
    Einstein tensor against their nonlinear electrodynamics; bardeen.md beside this file is
    the derivation."""
    f = "\\left(1 - \\dfrac{r_s\\,r^2}{\\left(r^2 + g^2\\right)^{3/2}}\\right)"
    bare = "1 - \\dfrac{r_s\\,r^2}{\\left(r^2 + g^2\\right)^{3/2}}"
    sphere = " + r^2\\left(d\\theta^2 + \\sin^2\\theta\\,d\\phi^2\\right)"
    domains = ["r \\in [0, \\infty)", "\\theta \\in [0, \\pi]", "\\phi \\in [0, 2\\pi)",
               "r = r_\\pm \\;\\text{(the horizons, for}\\; g < 2r_s/3\\sqrt{3}\\text{)}"]
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
    parameters = ["r_s", "g"]
    probe = vm.Reader(coords, parameters, ())
    r, rs, g = probe.symbol["r"], probe.parameters["r_s"], probe.parameters["g"]
    return {
        "metric_id": "bardeen",
        "system": {"id": system_id, "name": name, "coords": coords,
                   "domains": [coords[0] + " \\in (-\\infty, \\infty)"] + domains,
                   "parameters": parameters, "line_element": line},
        "chart_line_element": chart_line,
        "printer": {"lead": [r ** 2 + g ** 2, r, g, rs], "factors": [rs, g, r], "flip": False},
        "pretty": bardeen_pretty(r, g),
        "components": {"metric_components": metric, "inverse_metric_components": inverse},
        "ricci_scalar": "\\dfrac{3r_s\\,g^2\\left(4g^2 - r^2\\right)}{\\left(r^2 + g^2\\right)^{7/2}}",
        "kretschmann": ("\\dfrac{3r_s^2\\left(4r^8 - 12g^2\\,r^6 + 47g^4\\,r^4 - 4g^6\\,r^2 + 8g^8\\right)}"
                        "{\\left(r^2 + g^2\\right)^7}"),
        "check": bardeen_source,
    }


def bardeen_source(chart):
    """Ayon-Beato and Garcia's source (Phys. Lett. B 493, 149, 2000): the monopole field
    F = g sin(theta) dtheta ^ dphi of the nonlinear electrodynamics with Lagrangian
    L(F) = (3/(2 s g^2)) (sqrt(2 g^2 F)/(1 + sqrt(2 g^2 F)))^{5/2}, s = g/r_s, whose invariant
    F = F_{mu nu}F^{mu nu}/4 is g^2/(2r^4) on the monopole. Their Einstein equations are
    G_mu^nu = 2(L_F F_{mu lambda}F^{nu lambda} - delta_mu^nu L), which for a magnetic field read
    G^t_t = G^r_r = -2L and G^theta_theta = G^phi_phi = 2(2F L_F - L). The field is closed,
    and since L_F depends on r alone and sqrt(-g) F^{theta phi} = g/r^2 on it, the field
    equation d(L_F *F) = 0 holds."""
    r = chart.reader.symbol["r"]
    rs, g = chart.reader.parameters["r_s"], chart.reader.parameters["g"]
    F = sp.Symbol("F", positive=True)
    root = sp.sqrt(2 * g ** 2 * F)
    lagrangian = 3 * rs / (2 * g ** 3) * (root / (1 + root)) ** sp.Rational(5, 2)
    on = {F: g ** 2 / (2 * r ** 4)}
    L = sp.powdenest(sp.simplify(lagrangian.subs(on)), force=True)
    LF = sp.powdenest(sp.simplify(sp.diff(lagrangian, F).subs(on)), force=True)
    mixed = chart.geo.raise_indices(chart.geo.einstein_ll(), 2, (0,))
    wanted = [-2 * L, -2 * L, 2 * (2 * on[F] * LF - L), 2 * (2 * on[F] * LF - L)]
    for a in range(4):
        for b in range(4):
            if sp.simplify(mixed[a][b] - (wanted[a] if a == b else 0)) != 0:
                raise AssertionError(f"bardeen: the Einstein tensor misses the monopole's stress in slot "
                                     f"{chart.coords_tex[a]}{chart.coords_tex[b]}")
    # With p_r = -rho the strong energy condition is p_t >= 0, G^theta_theta >= 0.
    # The energy density, -G^t_t/(8 pi), is positive everywhere, and G^theta_theta changes sign
    # at r^2 = 2g^2/3, inside which the strong energy condition fails.
    tangential = 3 * rs * g ** 2 * (3 * r ** 2 - 2 * g ** 2) / (2 * (r ** 2 + g ** 2) ** sp.Rational(7, 2))
    if sp.simplify(wanted[2] - tangential) != 0:
        raise AssertionError("bardeen: the tangential pressure does not change sign at r^2 = 2g^2/3")


BARDEEN_CHARTS = ["static", "eddington_finkelstein_outgoing", "eddington_finkelstein_ingoing"]
CHARTS["bardeen"] = [lambda s=s: bardeen(s) for s in BARDEEN_CHARTS]


# -- Hayward ---------------------------------------------------------------------------

HAYWARD_CHARTS = ("static", "eddington_finkelstein_ingoing", "eddington_finkelstein_outgoing", "evaporating")


def hayward(system_id):
    """Hayward's regular black hole, F = 1 - 2 m r^2/(r^3 + 2 m ell^2), his (1) and (5), in his
    own two lengths: the static chart, the two Eddington-Finkelstein charts built on its tortoise
    coordinate, dr_*/dr = 1/F, of which the ingoing one is his (10), and the same ingoing chart
    with the mass a function m(v) of advanced time, the chart his black hole forms and evaporates
    in. Every value is printed around r^3 + 2 m ell^2 and r^3 - 2 m r^2 + 2 m ell^2, the
    denominator and the numerator of F, so that each chart is Schwarzschild's at ell = 0 term by
    term. hayward_check holds each chart to his Einstein tensor, to de Sitter space at the centre
    and to Schwarzschild's metric far away."""
    f = "\\left(1 - \\dfrac{2mr^2}{r^3 + 2m\\ell^2}\\right)"
    bare = "1 - \\dfrac{2mr^2}{r^3 + 2m\\ell^2}"
    sphere = " + r^2\\left(d\\theta^2 + \\sin^2\\theta\\,d\\phi^2\\right)"
    domains = ["r \\in [0, \\infty)", "\\theta \\in [0, \\pi]", "\\phi \\in [0, 2\\pi)"]
    marks = ["r = r_+ \\;\\text{(the outer horizon)}", "r = r_- \\;\\text{(the inner horizon)}"]
    parameters = ["m", "\\ell"]
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
        if system_id == "evaporating":
            name = "Forming and Evaporating"
            parameters = ["m = m(v)", "\\ell"]
            marks = ["r^3 - 2mr^2 + 2m\\ell^2 = 0 \\;\\text{(the trapping horizons)}"]
    probe = vm.Reader(coords, parameters, ())
    r, m, ell = probe.symbol["r"], probe.parameters["m"], probe.parameters["ell"]
    # A product is written as the line element writes it, the mass and the length before the
    # radius, 2 m ell^2 and 2 m r^2, and a sum by falling powers of r.
    printer = {"lead": [r, m, ell], "factors": [m, ell, r], "flip": False}
    if system_id == "evaporating":
        dm = sp.Derivative(m, probe.symbol["v"])
        printer["lead"] = [dm, r, m, ell]
        printer["factors"] = [m, ell, r, dm]
        # The rate of change of the mass stands apart from what the static chart has already.
        printer["collect"] = lambda poly, pr: cp.collect_by(poly, [dm], pr)
    return {
        "metric_id": "hayward",
        "system": {"id": system_id, "name": name, "coords": coords,
                   "domains": [coords[0] + " \\in (-\\infty, \\infty)"] + domains + marks,
                   "parameters": parameters, "line_element": line},
        "chart_line_element": chart_line,
        "printer": printer,
        "components": {"metric_components": metric, "inverse_metric_components": inverse},
        "check": hayward_check,
    }


def hayward_check(chart):
    """Hayward's Einstein tensor, his (7), (8) and (12), as the mixed components he writes:
    G^t_t = G^r_r = -12 ell^2 m^2/(r^3 + 2 ell^2 m)^2 and G^theta_theta = G^phi_phi =
    24(r^3 - ell^2 m) ell^2 m^2/(r^3 + 2 ell^2 m)^3, and with m a function of v the one further
    component G^r_v = 2 r^4 m'/(r^3 + 2 ell^2 m)^2, the ingoing flux of the radiation, which in
    retarded time is G^r_u with the opposite sign. At the centre F = 1 - r^2/ell^2 + O(r^5), his
    (3), de Sitter space of radius ell for every mass, and far away F = 1 - 2m/r + O(1/r^4), his
    (2), Schwarzschild's."""
    r = chart.reader.symbol["r"]
    m, ell = chart.reader.parameters["m"], chart.reader.parameters["ell"]
    D = r ** 3 + 2 * m * ell ** 2
    G = chart.geo.raise_indices(chart.geo.einstein_ll(), 2, (0,))
    want = sp.zeros(4, 4)
    want[0, 0] = want[1, 1] = -12 * ell ** 2 * m ** 2 / D ** 2
    want[2, 2] = want[3, 3] = 24 * ell ** 2 * m ** 2 * (r ** 3 - ell ** 2 * m) / D ** 3
    x0 = chart.symbols[0]
    if chart.coords_tex[0] != "t":
        sign = -1 if chart.coords_tex[0] == "u" else 1
        want[1, 0] = sign * 2 * r ** 4 * sp.diff(m, x0) / D ** 2
    for a in range(4):
        for b in range(4):
            if vm.norm(vm._at(G, (a, b)) - want[a, b]) != 0:
                raise AssertionError(f"hayward: the mixed Einstein tensor misses Hayward's in slot "
                                     f"{chart.coords_tex[a]}{chart.coords_tex[b]}")
    mass = sp.Symbol("mass", positive=True)
    F = (1 - 2 * mass * r ** 2 / (r ** 3 + 2 * mass * ell ** 2))
    if sp.series(F, r, 0, 5).removeO() != 1 - r ** 2 / ell ** 2:
        raise AssertionError("hayward: F is not 1 - r^2/ell^2 + O(r^5) at the centre")
    if sp.limit((F - 1 + 2 * mass / r) * r ** 3, r, sp.oo) != 0:
        raise AssertionError("hayward: F is not 1 - 2m/r + O(1/r^4) far away")


CHARTS["hayward"] = [lambda s=s: hayward(s) for s in HAYWARD_CHARTS]


# -- Kerr-de Sitter ----------------------------------------------------------------------

KDS_CHARTS = ["boyer_lindquist", "nonrotating", "kerr_ingoing", "kerr_outgoing", "kerr_schild"]
# The four names every chart defines, as Carter's metric is always written.
KDS_PARAMETERS = ["r_s", "a", "\\Lambda",
                  "\\Xi = 1 + \\dfrac{\\Lambda a^2}{3}",
                  "\\rho = \\sqrt{r^2 + a^2\\cos^2\\theta}",
                  "\\Delta_r = \\left(r^2 + a^2\\right)\\left(1 - \\dfrac{\\Lambda r^2}{3}\\right) - r_s r",
                  "\\Delta_\\theta = 1 + \\dfrac{\\Lambda a^2}{3}\\cos^2\\theta"]


class KerrDeSitterForms:
    """How a value of Kerr-de Sitter is written: in the four names the charts define.

    The checker hands every value back as a rational function of r, a, Lambda, r_s, sin(theta)
    and cos(theta). Here it is factored with the even powers of the sine written in the cosine,
    which is the angle rho and Delta_theta are written in; each factor that is a multiple of
    rho^2, Delta_r, Delta_theta or Xi is written by its name, (u + v)(u - v) is written
    u^2 - v^2, and every sum left over is written in the shortest of the forms `regrouped`
    tries. `joined` puts two such values over one denominator without multiplying them out,
    which is how the Riemann tensor is written as the Weyl tensor plus the part that carries
    Lambda."""

    def __init__(self, reader):
        self.r, self.theta = reader.symbol["r"], reader.symbol["\\theta"]
        self.rs, self.a, self.L = (reader.parameters[k] for k in ("r_s", "a", "Lambda"))
        self.s, self.c = sp.sin(self.theta), sp.cos(self.theta)
        self.C = sp.Symbol("KDS_cos")
        self.Xi, self.rho = sp.Symbol("Xi", positive=True), sp.Symbol("rho", positive=True)
        self.Dr, self.Dth = sp.Symbol("Delta_r"), sp.Symbol("Delta_theta", positive=True)
        r, a, L, C = self.r, self.a, self.L, self.C
        self.radial = (r ** 2 + a ** 2) * (1 - L * r ** 2 / 3) - self.rs * r
        self.named = [(self.Xi, 1 + L * a ** 2 / 3), (self.rho ** 2, r ** 2 + a ** 2 * C ** 2),
                      (self.Dr, self.radial), (self.Dth, 1 + L * a ** 2 * C ** 2 / 3)]

    def printer(self):
        return {"lead": [self.r, self.a, self.c, self.s], "rising": [self.L, self.rs], "flip": False,
                "overrides": {self.Dth: "\\Delta_\\theta"},
                "factors": [self.L, self.a, self.r, self.rs, self.Xi, self.rho, self.Dr, self.Dth]}

    def in_cosine(self, e):
        s = self.s
        e = e.replace(lambda p: p.is_Pow and p.base == s and p.exp.is_Integer and p.exp > 1,
                      lambda p: s ** (int(p.exp) % 2) * (1 - self.c ** 2) ** (int(p.exp) // 2))
        return e.subs(self.c, self.C)

    def in_sine(self, e):
        C = self.C
        return sp.expand(e.replace(lambda p: p.is_Pow and p.base == C and p.exp.is_Integer and p.exp > 1,
                                   lambda p: C ** (int(p.exp) % 2) * (1 - self.s ** 2) ** (int(p.exp) // 2)))

    def name_of(self, base):
        for placeholder, polynomial in self.named:
            ratio = sp.cancel(base / polynomial)
            if ratio.is_Number:
                return ratio * placeholder
        return None

    def product(self, value, depth=0):
        """A rational function of r, a, Lambda, r_s, the cosine and sin(theta), factored."""
        out, powers, C, s = sp.Integer(1), {}, self.C, self.s
        for f in sp.Mul.make_args(sp.factor(value)):
            base, k = (f.base, f.exp) if f.is_Pow else (f, sp.Integer(1))
            powers[base] = powers.get(base, 0) + k
        pairs = [b for b in powers if b.is_Add and len(b.args) == 2]
        for i, one in enumerate(pairs):
            for other in pairs[i + 1:]:
                both = sp.expand(one * other)
                k = min(powers[one], powers[other], key=abs) if powers[one] * powers[other] > 0 else 0
                if k and len(sp.Add.make_args(both)) == 2:
                    powers[one] -= k
                    powers[other] -= k
                    powers[both] = powers.get(both, 0) + k
        for base, sign in ((sp.expand(1 - C ** 2), 1), (sp.expand(C ** 2 - 1), -1)):
            k = powers.pop(base, 0)
            if k:
                powers[s] = powers.get(s, 0) + 2 * k
                out *= sp.Integer(sign) ** k
        factors = []
        for base, k in powers.items():
            if k == 0:
                continue
            if base.is_Add:
                named = self.name_of(base)
                base = named if named is not None else self.regrouped(base, depth)
            if (base ** k).is_Number:
                out *= base ** k
            else:
                factors.append(base ** k)
        # Every factor multiplied in one step, the number last and held in front of a lone sum, as
        # the 3 of 3(r^2 - a^2 cos^2(theta)): a number multiplied onto a sum alone is carried into it.
        rest = sp.Mul(*factors)
        number, rest = (out * rest.as_coeff_Mul()[0], rest.as_coeff_Mul()[1])
        return sp.Mul(number, rest, evaluate=False) if rest.is_Add and number != 1 else sp.Mul(number, rest)

    @staticmethod
    def summed(terms):
        """The terms as one sum, their common numerical factor in front of it."""
        terms = [t.as_coeff_Mul() + (monomial,) for t, monomial in terms if t != 0]
        content = sp.gcd_list([number for number, _, _ in terms])
        if all(number < 0 for number, _, _ in terms):
            content = -content
        # The number, the factored coefficient and the monomial multiplied in one step: a number
        # multiplied onto a sum alone would be carried into it.
        return sp.Mul(content, sp.Add(*[sp.Mul(number / content, rest, monomial) for number, rest, monomial in terms],
                                      evaluate=False), evaluate=False)

    def regrouped(self, base, depth):
        """A sum that is none of the named ones, in the shortest of: as it stands in the cosine
        or in the sine; grouped by the powers of r_s, of Lambda or of both, each coefficient
        factored in its turn; and, with r_s r written as (r^2 + a^2)(1 - Lambda r^2/3) - Delta_r,
        grouped by the powers of Delta_r."""
        candidates = [base, self.in_sine(base)]
        if depth == 0:
            for generators in ((self.rs,), (self.rs, self.L), (self.L,)):
                polynomial = sp.Poly(base, *generators)
                if len(polynomial.terms()) > 1:
                    candidates.append(self.summed([
                        (self.product(coefficient, 1), sp.Mul(*[g ** e for g, e in zip(generators, monomial)]))
                        for monomial, coefficient in polynomial.terms()]))
            if base.has(self.rs):
                swapped = sp.expand(base.subs(self.rs, (self.radial + self.rs * self.r - self.Dr) / self.r))
                if sp.denom(sp.together(swapped)).is_Number:
                    candidates.append(self.summed([(self.product(coefficient, 1), self.Dr ** k)
                                                   for (k,), coefficient in sp.Poly(swapped, self.Dr).terms()]))
        return min(candidates, key=lambda e: len(str(e)))

    def pretty(self, value):
        value = sp.sympify(value)
        if value == 0:
            return value
        return self.product(self.in_cosine(sp.together(value))).subs(self.C, self.c)

    @staticmethod
    def _powers(e):
        number, rest = e.as_coeff_Mul()
        powers = {}
        for f in sp.Mul.make_args(rest):
            if f == 1:
                continue
            base, k = (f.base, f.exp) if f.is_Pow else (f, sp.Integer(1))
            powers[base] = powers.get(base, 0) + k
        return number, powers

    def joined(self, pieces):
        """The sum of the pieces, each a product as `pretty` writes one, over their common
        denominator, with whatever every numerator shares taken out in front."""
        fractions = [sp.fraction(piece) for piece in pieces]
        bottoms = [self._powers(d) for _, d in fractions]
        number = sp.ilcm(*[n for n, _ in bottoms])
        low = {}
        for _, powers in bottoms:
            for base, k in powers.items():
                low[base] = max(low.get(base, 0), k)
        tops = []
        for (top, _), (n, powers) in zip(fractions, bottoms):
            c, own = self._powers(top)
            for base, k in low.items():
                if k != powers.get(base, 0):
                    own[base] = own.get(base, 0) + k - powers.get(base, 0)
            tops.append((c * number / n, own))
        content = sp.gcd_list([c for c, _ in tops])
        shared = {base: min(own.get(base, 0) for _, own in tops) for base in tops[0][1]}
        terms = [sp.Mul(c / content, *[base ** (k - shared.get(base, 0)) for base, k in own.items()])
                 for c, own in tops]
        return (content * sp.Mul(*[base ** k for base, k in shared.items()]) * sp.Add(*terms, evaluate=False)
                / (number * sp.Mul(*[base ** k for base, k in low.items()])))


def kerr_de_sitter_line(system, c):
    """The line element of a chart, with c the text that stands before the differential of its
    time: empty in the chart x^0 = ct, and c where the time is the named one."""
    tetrad = ("\\dfrac{\\rho^2}{\\Delta_\\theta}d\\theta^2 + \\dfrac{\\Delta_\\theta\\sin^2\\theta}{\\rho^2}"
              "\\left(a\\,{T} - \\dfrac{r^2 + a^2}{\\Xi}d{P}\\right)^2")
    principal = "\\left({T} - \\dfrac{a\\sin^2\\theta}{\\Xi}d{P}\\right)"
    if system == "boyer_lindquist":
        line = ("ds^2 = -\\dfrac{\\Delta_r}{\\rho^2}" + principal + "^2 + \\dfrac{\\rho^2}{\\Delta_r}dr^2 + " + tetrad)
        return line.replace("{T}", c + "dt").replace("{P}", "\\phi")
    if system == "nonrotating":
        return ("ds^2 = -\\dfrac{\\Delta_r}{\\Xi^2\\rho^2}\\left(\\Delta_\\theta\\," + c + "dt - a\\sin^2\\theta\\,d\\Phi\\right)^2"
                " + \\dfrac{\\rho^2}{\\Delta_r}dr^2 + \\dfrac{\\rho^2}{\\Delta_\\theta}d\\theta^2"
                " + \\dfrac{\\Delta_\\theta\\sin^2\\theta}{\\Xi^2\\rho^2}\\left(a\\left(1 - \\dfrac{\\Lambda r^2}{3}\\right)" + c
                + "dt - \\left(r^2 + a^2\\right)d\\Phi\\right)^2")
    if system in ("kerr_ingoing", "kerr_outgoing"):
        null, sign = ("v", "+") if system == "kerr_ingoing" else ("u", "-")
        line = ("ds^2 = -\\dfrac{\\Delta_r}{\\rho^2}" + principal + "^2 " + sign + " 2\\,dr" + principal + " + " + tetrad)
        return line.replace("{T}", "d" + null).replace("{P}", "\\tilde\\phi")
    static = "\\left(1 - \\dfrac{\\Lambda r^2}{3}\\right)"
    radial = "\\dfrac{\\rho^2}{" + static + "\\left(r^2 + a^2\\right)}"
    return ("ds^2 = -\\dfrac{\\Delta_\\theta}{\\Xi}" + static + c.replace("c", "c^2") + "d\\tau^2 + " + radial + "dr^2"
            " + \\dfrac{\\rho^2}{\\Delta_\\theta}d\\theta^2 + \\dfrac{r^2 + a^2}{\\Xi}\\sin^2\\theta\\,d\\psi^2"
            " + \\dfrac{r_s r}{\\rho^2}\\left(\\dfrac{\\Delta_\\theta}{\\Xi}" + c + "d\\tau + " + radial + "dr"
            " - \\dfrac{a\\sin^2\\theta}{\\Xi}d\\psi\\right)^2")


def kerr_de_sitter(system):
    """Carter's rotating black hole with a cosmological constant of either sign, in five charts:
    his own, which is Boyer and Lindquist's at Lambda = 0; the same chart turned at the angular
    velocity a Lambda c/3, in which the frame at infinity does not rotate; the ingoing and
    outgoing Kerr charts, along the two principal null congruences; and the Kerr-Schild chart
    of Gibbons, Lu, Page and Pope, the de Sitter metric in spheroidal coordinates plus r_s r/rho^2
    times the square of a null one form. kerr_de_sitter_check holds each chart to
    R_mu_nu = Lambda g_mu_nu, to Riemann = Weyl + (Lambda/3)(g g - g g), and each chart after
    the first to being the first pulled back. kerr_de_sitter.md beside this file is the
    derivation."""
    time, angle, name = {"boyer_lindquist": ("t", "\\phi", "Boyer-Lindquist"),
                         "nonrotating": ("t", "\\Phi", "Nonrotating at Infinity"),
                         "kerr_ingoing": ("v", "\\tilde\\phi", "Ingoing Kerr"),
                         "kerr_outgoing": ("u", "\\tilde\\phi", "Outgoing Kerr"),
                         "kerr_schild": ("\\tau", "\\psi", "Kerr-Schild")}[system]
    coords = [time, "r", "\\theta", angle]
    forms = KerrDeSitterForms(vm.Reader(coords, KDS_PARAMETERS, ()))
    timed = time in ("t", "\\tau")
    # Carter's chart and the one turned from it hold between two roots of Delta_r, and are published
    # for the region outside the event horizon; the Kerr charts run through every root; and the
    # Kerr-Schild chart's own time and angle fail where the de Sitter background's static chart does.
    if system in ("boyer_lindquist", "nonrotating"):
        radial = ["r \\in (r_+, r_c) \\;\\text{for}\\; \\Lambda > 0", "r \\in (r_+, \\infty) \\;\\text{for}\\; \\Lambda \\le 0"]
        notes = ["\\Delta_r = 0 \\;\\text{(the event horizon } r_+ \\text{ and, for } \\Lambda > 0 \\text{, the cosmological horizon } r_c\\text{)}"]
    else:
        radial = (["r \\in (-\\infty, \\infty)"] if system != "kerr_schild" else
                  ["r \\in \\left(-\\sqrt{3/\\Lambda}, \\sqrt{3/\\Lambda}\\right) \\;\\text{for}\\; \\Lambda > 0",
                   "r \\in (-\\infty, \\infty) \\;\\text{for}\\; \\Lambda \\le 0"])
        notes = ["\\Delta_r = 0 \\;\\text{(the horizons)}", "\\rho = 0 \\;\\text{(the ring singularity)}"]
    domains = [time + " \\in (-\\infty, \\infty)"] + radial + ["\\theta \\in [0, \\pi]", angle + " \\in [0, 2\\pi)"] + notes
    spec = {
        "metric_id": "kerr_de_sitter",
        "system": {"id": system, "name": name, "coords": coords, "domains": domains,
                   "parameters": KDS_PARAMETERS, "line_element": kerr_de_sitter_line(system, "c\\," if timed else "")},
        "chart_line_element": kerr_de_sitter_line(system, ""),
        "printer": forms.printer(),
        "pretty": forms.pretty,
        "bracketed": forms.pretty,
        "check": lambda chart: kerr_de_sitter_check(chart, system),
        "after": lambda math, chart: kerr_de_sitter_riemann(math, chart, forms),
        "ricci_scalar": "4\\Lambda",
        "kretschmann": ("\\dfrac{12r_s^2\\left(r^2 - a^2\\cos^2\\theta\\right)\\left(\\rho^4 - 16a^2r^2\\cos^2\\theta\\right)}"
                        "{\\rho^{12}} + \\dfrac{8\\Lambda^2}{3}"),
    }
    if system == "boyer_lindquist":
        block = "\\left(\\left(r^2 + a^2\\right)\\Delta_\\theta - \\Delta_r\\right)"
        square = "\\left(\\left(r^2 + a^2\\right)^2\\Delta_\\theta - a^2\\Delta_r\\sin^2\\theta\\right)"
        lapse = "\\left(\\Delta_r - a^2\\Delta_\\theta\\sin^2\\theta\\right)"
        cross = "-\\dfrac{a" + block + "\\sin^2\\theta}{\\Xi\\rho^2}"
        up = "-\\dfrac{a\\,\\Xi" + block + "}{\\rho^2\\Delta_r\\Delta_\\theta}"
        spec["components"] = {
            "metric_components": {("t", "t"): "-\\dfrac{\\Delta_r - a^2\\Delta_\\theta\\sin^2\\theta}{\\rho^2}",
                                  ("t", "\\phi"): cross, ("\\phi", "t"): cross,
                                  ("\\phi", "\\phi"): "\\dfrac{" + square + "\\sin^2\\theta}{\\Xi^2\\rho^2}"},
            "inverse_metric_components": {
                ("t", "t"): "-\\dfrac{\\left(r^2 + a^2\\right)^2\\Delta_\\theta - a^2\\Delta_r\\sin^2\\theta}"
                            "{\\rho^2\\Delta_r\\Delta_\\theta}",
                ("t", "\\phi"): up, ("\\phi", "t"): up,
                ("\\phi", "\\phi"): "\\dfrac{\\Xi^2" + lapse + "}{\\rho^2\\Delta_r\\Delta_\\theta\\sin^2\\theta}"}}
    return spec


def kerr_de_sitter_jacobian(chart, system):
    """d(Boyer-Lindquist coordinate)/d(the chart's), in the chart x^0 = ct of both."""
    r = chart.symbols[1]
    rs, a, L, Xi, Dr = (chart.reader.parameters[k] for k in ("r_s", "a", "Lambda", "Xi", "Delta_r"))
    if system == "nonrotating":
        # phi = Phi + (a Lambda/3) ct
        return sp.Matrix([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [a * L / 3, 0, 0, 1]])
    if system in ("kerr_ingoing", "kerr_outgoing"):
        # dv = c dt + (r^2 + a^2) dr/Delta_r and d(phi~) = d(phi) + a Xi dr/Delta_r, and their mirror for u
        sign = 1 if system == "kerr_ingoing" else -1
        return sp.Matrix([[1, -sign * (r ** 2 + a ** 2) / Dr, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0],
                          [0, -sign * a * Xi / Dr, 0, 1]])
    # c d(tau) = c dt + r_s r dr/((1 - Lambda r^2/3) Delta_r) and
    # d(psi) = d(phi) - (a Lambda/3) c dt + a r_s r dr/((r^2 + a^2) Delta_r)
    static = 1 - L * r ** 2 / 3
    time = -rs * r / (static * Dr)
    return sp.Matrix([[1, time, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0],
                      [a * L / 3, a * L * time / 3 - a * rs * r / ((r ** 2 + a ** 2) * Dr), 0, 1]])


def kerr_de_sitter_check(chart, system):
    g, geo, L = chart.geo.g, chart.geo, chart.reader.parameters["Lambda"]
    ricci, riemann, weyl = geo.ricci_ll(), geo.riemann_llll(), geo.weyl_llll()
    for i in range(4):
        for j in range(4):
            if vm.norm(ricci[i][j] - L * g[i, j]) != 0:
                raise AssertionError(f"kerr_de_sitter/{system}: R_mu_nu is not Lambda g_mu_nu in slot {i}{j}")
    for i, j, k, l in itertools.product(range(4), repeat=4):
        constant = L * (g[i, k] * g[j, l] - g[i, l] * g[j, k]) / 3
        if vm.norm(riemann[i][j][k][l] - weyl[i][j][k][l] - constant) != 0:
            raise AssertionError(f"kerr_de_sitter/{system}: Riemann is not Weyl plus the constant curvature part")
    if system == "boyer_lindquist":
        return
    spec = kerr_de_sitter("boyer_lindquist")
    source = cp.Chart(spec["system"]["coords"], KDS_PARAMETERS, spec["chart_line_element"])
    at = dict(zip(source.symbols, chart.symbols))
    jacobian = kerr_de_sitter_jacobian(chart, system)
    pulled = jacobian.T * source.geo.g.subs(at, simultaneous=True) * jacobian
    for i in range(4):
        for j in range(i, 4):
            if vm.norm(pulled[i, j] - g[i, j]) != 0:
                raise AssertionError(f"kerr_de_sitter/{system}: the Boyer-Lindquist chart pulled back misses slot "
                                     f"{chart.coords_tex[i]}{chart.coords_tex[j]}")


def kerr_de_sitter_riemann(math, chart, forms):
    """The Riemann tensor written as the Weyl tensor, which is Kerr's with the four names in
    place of Kerr's, plus (Lambda/3)(g g - g g), the two over one denominator."""
    geo = chart.geo
    riemann, weyl = geo.riemann_llll(), geo.weyl_llll()
    blocks = {"llll": (riemann, weyl),
              "ulll": (geo.raise_indices(riemann, 4, (0,)), geo.raise_indices(weyl, 4, (0,)))}
    parts = {}
    for whole, traceless in blocks.values():
        for index in vm._indices(4, 4):
            value, conformal = vm._at(whole, index), vm._at(traceless, index)
            constant = vm.norm(value - conformal)
            if value != 0 and conformal != 0 and constant != 0:
                parts[value], parts[-value] = (conformal, constant), (-conformal, -constant)

    def split(value):
        found = parts.get(sp.sympify(value))
        return forms.joined([forms.pretty(part) for part in found]) if found else forms.pretty(value)

    chart.pretty = chart.bracketed = split
    for variance, (whole, _) in blocks.items():
        math["riemann"]["variants"][variance]["nonzero"] = chart.block(whole, 4)
    chart.pretty = chart.bracketed = forms.pretty
    return math


CHARTS["kerr_de_sitter"] = [lambda s=s: kerr_de_sitter(s) for s in KDS_CHARTS]


def write(spec):
    start = time.time()
    chart = cp.Chart(spec["system"]["coords"], spec["system"]["parameters"], spec["chart_line_element"],
                     spec["printer"], spec.get("pretty"), spec.get("time"), spec.get("bracketed"), spec.get("reduce"))
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
        elif field == "kretschmann" and "kretschmann_where" in spec:
            # A square of a delta is no distribution, so the scalar is stated where the kink is not.
            where = spec["kretschmann_where"]
            off = vm.off_support(computed, chart.reader.symbol[where])
            math[field] = ("K = " + chart.text(off)
                           + " \\;\\text{for}\\; " + where + " \\neq 0")
        elif field == "kretschmann" and "kretschmann_text" in spec:
            math[field] = "K = " + chart.check(spec["kretschmann_text"](chart), computed)
        elif field == "kretschmann":
            math[field] = "K = " + chart.text(computed)
    if spec.get("bare_scalar"):
        math["ricci_scalar"] = math["ricci_scalar"].removeprefix("R = ")
    if "rewrite" in spec:
        math = rewritten(math, spec["rewrite"], chart)
    if "after" in spec:
        # A chart whose values hold only on a surface of its parameters states them there.
        math = spec["after"](math, chart)
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
