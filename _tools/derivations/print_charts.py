#!/usr/bin/env python3
"""Compute and write the coordinate systems whose mathematics is printed by machine: the
charts of tov, malament_hogarth, mixmaster and lentz, and Godel's cylindrical chart.

    /tmp/mfs-venv/bin/python _tools/derivations/print_charts.py [--metric <id>]...
    /tmp/mfs-venv/bin/python _tools/derivations/verify_metrics.py --system <id>/<system>

Each chart below names its coordinates, its parameters, the line element it publishes, the
same line element in the chart x^0 = ct that the components are printed in, and how its values
are to be grouped for reading. The script computes every tensor from that line element with
the checker's own Geometry, prints each component through chart_printer.py, reads every
printed value back and compares it with what it was printed from, and writes the result into
the metric file. It writes only the mathematics: the entry's prose, its convention, and the
description of each parameter stay in the metric file and are carried over untouched, and a
parameter the file does not yet describe stops the script rather than being written without
one. The Ricci scalar and the Kretschmann scalar are written by hand where a structured form
reads better than an expanded one, and each of those is checked against sympy here too.

The derivations these charts rest on, and the reason each was chosen, are in tov.md,
malament_hogarth.md, mixmaster.md, lentz.md and godel.md beside this file.
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
SYSTEM_ORDER = ["id", "name", "coords", "domains", "parameters", "line_element", "metric_components",
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


CHARTS = {"tov": tov, "malament_hogarth": malament_hogarth, "mixmaster": mixmaster, "lentz": lentz, "godel": godel}


def write(spec):
    start = time.time()
    chart = cp.Chart(spec["system"]["coords"], spec["system"]["parameters"], spec["chart_line_element"],
                     spec["printer"], spec.get("pretty"))
    if "check" in spec:
        spec["check"](chart)
    math = chart.mathematics()
    for field in ("ricci_scalar", "kretschmann"):
        computed = chart.geo.ricci_scalar() if field == "ricci_scalar" else chart.geo.kretschmann()
        if field in spec:
            math[field] = ("R = " if field == "ricci_scalar" else "K = ") + chart.check(spec[field], computed)
        elif field == "kretschmann":
            math[field] = "K = " + chart.text(computed)

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
    system = dict(spec["system"])
    system["parameters"] = [{"symbol": s, "description": described[s]} for s in spec["system"]["parameters"]]
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
        write(CHARTS[metric_id]())


if __name__ == "__main__":
    main()
