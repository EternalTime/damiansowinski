"""Charts and cold-disc boundary equations of Lynden-Bell and Pineault (1978).

The Weyl chart includes the finite discs of Paper I section 3.2.
The other charts describe infinite discs, on the open upper vacuum side.
Reflection across the equator supplies the lower side and its thin matter source.
"""
import math

import numpy as np
import sympy as sp

import chart_printer as cp
import verify_metrics as vm

METRIC = "lynden_bell_pineault"
SYSTEMS = ("weyl", "spherical", "photons", "rotating")


def spec(system):
    """One chart in the conventions of Paper I (3.1), (3.16), (4.2), or Paper II (2.2)."""
    import print_charts as pc

    if system == "weyl":
        coords = ["t", "R", "z", "\\phi"]
        parameters = ["\\nu = \\nu(R,z)", "\\xi = \\xi(R,z)", "n", "A"]
        line = ("ds^2 = -e^{2\\nu}{c2}dt^2 + e^{2\\xi - 2\\nu}"
                "\\left(dR^2 + dz^2\\right) + R^2e^{-2\\nu}d\\phi^2")
        domains = ["t \\in (-\\infty, \\infty)", "R \\in (0, \\infty)",
                   "z \\in (0, \\infty)", "\\phi \\in [0, 2\\pi)", "0 < n \\le 1/2", "A > 0"]
        name = "Weyl, Finite and Infinite Counter-Rotating Discs"
    else:
        coords = ["t", "r", "\\theta", "\\phi"]
        domains = ["t \\in (-\\infty, \\infty)", "r \\in (0, \\infty)",
                   "\\theta \\in (0, \\pi/2)", "\\phi \\in [0, 2\\pi)", "A > 0"]
        if system == "rotating":
            parameters = ["n", "A", "N = N(\\theta)", "P = P(\\theta)",
                          "Q = Q(\\theta)", "Z = Z(\\theta)", "f = r/A"]
            domains += ["0 < n < 1/2"]
            line = ("ds^2 = -f^{2n}e^N{c2}dt^2 + A^2f^{2 - 2n}e^{2P - N}"
                    "\\left(d\\phi - \\dfrac{f^{2n - 1}e^{N - P}Q}{A}{c}dt\\right)^2"
                    " + f^{-2n}e^{Z - N}\\left(dr^2 + r^2d\\theta^2\\right)")
            name = "Rotating Self-Similar Disc"
        elif system == "photons":
            parameters = ["A", "h = (1 + \\cos\\theta)/2"]
            line = ("ds^2 = -\\dfrac{rh}{A}{c2}dt^2 + \\dfrac{A}{r}"
                    "\\left(dr^2 + r^2d\\theta^2\\right) + \\dfrac{Ar\\sin^2\\theta}{h}d\\phi^2")
            name = "Counter-Rotating Massless Particles"
        else:
            parameters = ["n", "A", "f = r/A", "h = (1 + \\cos\\theta)/2"]
            domains += ["0 \\le n \\le 1/2"]
            line = ("ds^2 = -f^{2n}h^{2n}{c2}dt^2 + f^{-2n}h^{4n^2 - 2n}"
                    "\\left(dr^2 + r^2d\\theta^2\\right)"
                    " + f^{-2n}h^{-2n}r^2\\sin^2\\theta\\,d\\phi^2")
            name = "Infinite Counter-Rotating Disc"
    probe = vm.Reader(coords, parameters, ())
    state = {}
    if system in ("spherical", "rotating"):
        r, A = probe.symbol["r"], probe.parameters["A"]
        names = {"f": {r: 1, A: -1}}
        if system == "spherical":
            th = probe.symbol["\\theta"]
            names["h"] = {1 + sp.cos(th): 1, sp.Integer(2): -1}
        pretty = pc.named_powers(names, state)
    else:
        pretty = lambda value: sp.powsimp(sp.factor(sp.sympify(value)), combine="exp")

    def check(chart):
        state["printer"] = chart.printer
        check_chart(chart, system)

    return {"metric_id": METRIC,
            "system": {"id": system, "name": name, "coords": coords, "domains": domains,
                       "parameters": parameters, "line_element": line.replace("{c2}", "c^2").replace("{c}", "c")},
            "chart_line_element": line.replace("{c2}", "").replace("{c}", ""),
            "printer": {"lead": [*probe.parameters.values(), *probe.symbol.values()], "flip": False},
            "pretty": pretty, "check": check}


def check_chart(chart, system):
    """Pointwise source comparisons, completed beside each chart's implementation."""
    if system not in SYSTEMS:
        raise AssertionError(system)


def register(pc):
    """Register only this spacetime with the shared chart printer."""
    pc.CHARTS[METRIC] = [lambda s=s: spec(s) for s in SYSTEMS]
