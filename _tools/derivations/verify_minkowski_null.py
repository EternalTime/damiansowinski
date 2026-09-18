#!/usr/bin/env python3
"""Check the two null coordinate systems of minkowski.json against sympy.

Builds the Double Null and Spherical Null metrics from the Cartesian line element
that the same file publishes, computes the inverse metric, the Christoffel symbols
in both variants, the Riemann tensor, the Ricci tensor, the Ricci scalar, the
Kretschmann scalar, the Einstein tensor and the Weyl tensor, then reads every value
the two entries print and asserts that it is what sympy got. Exits non-zero on any
mismatch and names what disagreed.

    python3 -m venv .venv && .venv/bin/pip install sympy
    .venv/bin/python _tools/derivations/verify_minkowski_null.py

The companion derivation is _tools/derivations/minkowski_null.md, which carries the
same algebra by hand.
"""

import json
import re
import sys
from pathlib import Path

import sympy as sp
from sympy.parsing.sympy_parser import (
    implicit_multiplication,
    parse_expr,
    standard_transformations,
)

METRIC_FILE = Path(__file__).resolve().parents[2] / "MFS/assets/data/metrics/minkowski.json"

TRANSFORMS = standard_transformations + (implicit_multiplication,)

c, u, v, y, z, theta, phi = sp.symbols("c u v y z theta phi", real=True)
t, x, r = sp.symbols("t x r", real=True)

COORD_SYMBOL = {"t": t, "x": x, "y": y, "z": z, "r": r, "u": u, "v": v, "\\theta": theta, "\\phi": phi}
DIFFERENTIAL = {name: sp.Symbol("d" + str(symbol)) for name, symbol in COORD_SYMBOL.items()}
DOT = {name: sp.Symbol(str(symbol) + "_dot") for name, symbol in COORD_SYMBOL.items()}
DDOT = {name: sp.Symbol(str(symbol) + "_ddot") for name, symbol in COORD_SYMBOL.items()}

ALLOWED = (
    set(COORD_SYMBOL.values())
    | set(DIFFERENTIAL.values())
    | set(DOT.values())
    | set(DDOT.values())
    | {c}
)

LOCAL = {str(symbol): symbol for symbol in ALLOWED}


class LatexError(Exception):
    pass


def _matching_brace(text, opening):
    depth = 0
    for position in range(opening, len(text)):
        if text[position] == "{":
            depth += 1
        elif text[position] == "}":
            depth -= 1
            if depth == 0:
                return position
    raise LatexError(f"unbalanced braces in {text!r}")


def _expand_fractions(text):
    while True:
        found = re.search(r"\\d?frac\s*\{", text)
        if not found:
            return text
        first_open = text.index("{", found.start())
        first_close = _matching_brace(text, first_open)
        second_open = text.index("{", first_close)
        second_close = _matching_brace(text, second_open)
        numerator = text[first_open + 1:first_close]
        denominator = text[second_open + 1:second_close]
        text = f"{text[:found.start()]}(({numerator})/({denominator})){text[second_close + 1:]}"


def to_sympy(latex):
    """Translate the LaTeX dialect the metric files are written in into a sympy expression."""
    text = latex
    text = re.sub(r"\\ddot\s*\{\s*\\?([A-Za-z]+)\s*\}", r" \1_ddot ", text)
    text = re.sub(r"\\dot\s*\{\s*\\?([A-Za-z]+)\s*\}", r" \1_dot ", text)
    text = text.replace("\\left", " ").replace("\\right", " ")
    text = re.sub(r"\\[,;:!> ]", " ", text)
    text = _expand_fractions(text)
    text = text.replace("d\\theta", "dtheta").replace("d\\phi", "dphi")
    text = text.replace("\\theta", "theta").replace("\\phi", "phi")
    text = re.sub(r"\\(sin|cos|tan|cot)\s*\^\s*(\d+)\s*(theta|phi)", r"\1(\3)**\2", text)
    text = re.sub(r"\\(sin|cos|tan|cot)\s*(theta|phi)", r"\1(\2)", text)
    text = text.replace("^", "**")
    if "\\" in text:
        raise LatexError(f"{latex!r} still carries LaTeX this reader does not know: {text!r}")
    expression = parse_expr(text, local_dict=LOCAL, transformations=TRANSFORMS, evaluate=True)
    unknown = expression.free_symbols - ALLOWED
    if unknown:
        raise LatexError(f"{latex!r} parsed to unknown symbols {sorted(map(str, unknown))}")
    return sp.simplify(expression)


def metric_from_line_element(line_element, coords):
    """Read g_{mu nu} off a line element, in the chart its own coords name."""
    _, _, right = line_element.partition("=")
    form = sp.expand(to_sympy(right))
    differentials = [DIFFERENTIAL[name] for name in coords]
    polynomial = sp.Poly(form, *differentials)
    for powers in polynomial.monoms():
        if sum(powers) != 2:
            raise LatexError(f"{line_element!r} has a term of degree {sum(powers)} in the differentials")
    size = len(coords)
    g = sp.zeros(size, size)
    for i in range(size):
        for j in range(size):
            if i == j:
                g[i, j] = form.coeff(differentials[i], 2)
            else:
                mixed = sp.expand(form).coeff(differentials[i], 1).coeff(differentials[j], 1)
                g[i, j] = mixed / 2
    return sp.simplify(g)


def pull_back(g_old, old_coords, new_coords, substitution):
    """g_new[mu][nu] = (dx^a/dy^mu)(dx^b/dy^nu) g_old[a][b], with x^a given by substitution."""
    jacobian = sp.Matrix(
        len(old_coords),
        len(new_coords),
        lambda a, mu: sp.diff(substitution[old_coords[a]], new_coords[mu]),
    )
    replaced = g_old.subs({symbol: substitution[symbol] for symbol in old_coords}, simultaneous=True)
    return sp.simplify(jacobian.T * replaced * jacobian)


class Geometry:
    def __init__(self, g, coords):
        self.n = len(coords)
        self.coords = coords
        self.g = sp.simplify(g)
        self.ginv = sp.simplify(self.g.inv())
        self.christoffel_ull = self._christoffel_ull()
        self.christoffel_lll = self._christoffel_lll()
        self.riemann_ulll = self._riemann_ulll()
        self.riemann_llll = self._lower_first(self.riemann_ulll)
        self.ricci_ll = self._ricci_ll()
        self.ricci_scalar = sp.simplify(
            sum(self.ginv[a, b] * self.ricci_ll[a][b] for a in range(self.n) for b in range(self.n))
        )
        self.kretschmann = self._kretschmann()
        self.einstein_ll = self._einstein_ll()
        self.weyl_llll = self._weyl_llll()
        self.weyl_ulll = self._raise_first(self.weyl_llll)

    def _zeros(self, rank):
        if rank == 1:
            return [sp.Integer(0)] * self.n
        return [self._zeros(rank - 1) for _ in range(self.n)]

    def _christoffel_ull(self):
        out = self._zeros(3)
        for mu in range(self.n):
            for nu in range(self.n):
                for rho in range(self.n):
                    total = sum(
                        self.ginv[mu, alpha]
                        * (
                            sp.diff(self.g[alpha, rho], self.coords[nu])
                            + sp.diff(self.g[alpha, nu], self.coords[rho])
                            - sp.diff(self.g[nu, rho], self.coords[alpha])
                        )
                        for alpha in range(self.n)
                    )
                    out[mu][nu][rho] = sp.simplify(total / 2)
        return out

    def _christoffel_lll(self):
        out = self._zeros(3)
        for mu in range(self.n):
            for nu in range(self.n):
                for rho in range(self.n):
                    out[mu][nu][rho] = sp.simplify(
                        sum(self.g[mu, alpha] * self.christoffel_ull[alpha][nu][rho] for alpha in range(self.n))
                    )
        return out

    def _riemann_ulll(self):
        out = self._zeros(4)
        gamma = self.christoffel_ull
        for mu in range(self.n):
            for nu in range(self.n):
                for rho in range(self.n):
                    for sigma in range(self.n):
                        out[mu][nu][rho][sigma] = sp.simplify(
                            sp.diff(gamma[mu][nu][sigma], self.coords[rho])
                            - sp.diff(gamma[mu][nu][rho], self.coords[sigma])
                            + sum(
                                gamma[mu][rho][lam] * gamma[lam][nu][sigma]
                                - gamma[mu][sigma][lam] * gamma[lam][nu][rho]
                                for lam in range(self.n)
                            )
                        )
        return out

    def _lower_first(self, tensor):
        out = self._zeros(4)
        for mu in range(self.n):
            for nu in range(self.n):
                for rho in range(self.n):
                    for sigma in range(self.n):
                        out[mu][nu][rho][sigma] = sp.simplify(
                            sum(self.g[mu, alpha] * tensor[alpha][nu][rho][sigma] for alpha in range(self.n))
                        )
        return out

    def _raise_first(self, tensor):
        out = self._zeros(4)
        for mu in range(self.n):
            for nu in range(self.n):
                for rho in range(self.n):
                    for sigma in range(self.n):
                        out[mu][nu][rho][sigma] = sp.simplify(
                            sum(self.ginv[mu, alpha] * tensor[alpha][nu][rho][sigma] for alpha in range(self.n))
                        )
        return out

    def _ricci_ll(self):
        out = self._zeros(2)
        for nu in range(self.n):
            for sigma in range(self.n):
                out[nu][sigma] = sp.simplify(
                    sum(self.riemann_ulll[mu][nu][mu][sigma] for mu in range(self.n))
                )
        return out

    def _kretschmann(self):
        upper = self._zeros(4)
        for a in range(self.n):
            for b in range(self.n):
                for cc in range(self.n):
                    for d in range(self.n):
                        upper[a][b][cc][d] = sp.simplify(
                            sum(
                                self.ginv[b, q] * self.ginv[cc, s] * self.ginv[d, w] * self.riemann_ulll[a][q][s][w]
                                for q in range(self.n)
                                for s in range(self.n)
                                for w in range(self.n)
                            )
                        )
        return sp.simplify(
            sum(
                self.riemann_llll[a][b][cc][d] * upper[a][b][cc][d]
                for a in range(self.n)
                for b in range(self.n)
                for cc in range(self.n)
                for d in range(self.n)
            )
        )

    def _einstein_ll(self):
        out = self._zeros(2)
        for a in range(self.n):
            for b in range(self.n):
                out[a][b] = sp.simplify(self.ricci_ll[a][b] - self.ricci_scalar * self.g[a, b] / 2)
        return out

    def _weyl_llll(self):
        n = self.n
        out = self._zeros(4)
        for a in range(n):
            for b in range(n):
                for cc in range(n):
                    for d in range(n):
                        correction = (
                            self.g[a, cc] * self.ricci_ll[d][b]
                            - self.g[a, d] * self.ricci_ll[cc][b]
                            - self.g[b, cc] * self.ricci_ll[d][a]
                            + self.g[b, d] * self.ricci_ll[cc][a]
                        ) / (n - 2)
                        trace = (
                            self.ricci_scalar
                            * (self.g[a, cc] * self.g[d, b] - self.g[a, d] * self.g[cc, b])
                            / ((n - 1) * (n - 2))
                        )
                        out[a][b][cc][d] = sp.simplify(
                            self.riemann_llll[a][b][cc][d] - correction + trace
                        )
        return out


class Report:
    def __init__(self):
        self.failures = []

    def check(self, condition, message):
        if not condition:
            self.failures.append(message)

    def equal(self, got, expected, message):
        difference = sp.simplify(sp.expand(sp.together(got - expected)))
        if difference != 0:
            self.failures.append(f"{message}: file says {got}, sympy says {expected}")

    def note(self, message):
        print(message)


def index_of(coords, name, where, report):
    if name not in coords:
        report.check(False, f"{where}: index {name!r} is not one of {coords}")
        return None
    return coords.index(name)


def check_rank2(report, published, matrix, coords, where):
    seen = set()
    for entry in published:
        i = index_of(coords, entry["indices"][0], where, report)
        j = index_of(coords, entry["indices"][1], where, report)
        if i is None or j is None:
            continue
        seen.add((i, j))
        report.equal(to_sympy(entry["value"]), matrix[i, j], f"{where} [{entry['indices'][0]},{entry['indices'][1]}]")
    for i in range(len(coords)):
        for j in range(len(coords)):
            if (i, j) not in seen and sp.simplify(matrix[i, j]) != 0:
                report.check(False, f"{where}: [{coords[i]},{coords[j]}] is missing, sympy says {matrix[i, j]}")


def check_rank3(report, published, tensor, coords, where):
    seen = set()
    for entry in published:
        position = [index_of(coords, name, where, report) for name in entry["indices"]]
        if any(p is None for p in position):
            continue
        seen.add(tuple(position))
        i, j, k = position
        report.equal(to_sympy(entry["value"]), tensor[i][j][k], f"{where} {entry['indices']}")
    for i in range(len(coords)):
        for j in range(len(coords)):
            for k in range(len(coords)):
                if (i, j, k) not in seen and sp.simplify(tensor[i][j][k]) != 0:
                    report.check(
                        False,
                        f"{where}: [{coords[i]},{coords[j]},{coords[k]}] is missing, "
                        f"sympy says {tensor[i][j][k]}",
                    )


def check_all_vanish(report, variants, where):
    for name, variant in variants.items():
        report.check(
            variant["nonzero"] == [],
            f"{where}.{name} lists {len(variant['nonzero'])} components, sympy says every one vanishes",
        )


def check_geodesics(report, published, geometry, coords, where):
    report.check(
        len(published) == len(coords),
        f"{where}: {len(published)} equations for {len(coords)} coordinates",
    )
    for position, equation in enumerate(published):
        left, _, right = equation.partition("=")
        name = coords[position]
        report.check(
            to_sympy(left) == DDOT[name],
            f"{where}: equation {position} is for {to_sympy(left)}, expected {DDOT[name]}",
        )
        expected = -sum(
            geometry.christoffel_ull[position][nu][rho] * DOT[coords[nu]] * DOT[coords[rho]]
            for nu in range(len(coords))
            for rho in range(len(coords))
        )
        report.equal(to_sympy(right), sp.simplify(expected), f"{where}: {name} equation right hand side")


def system(metric, wanted):
    for entry in metric["coordinates"]:
        if entry["id"] == wanted:
            return entry
    raise SystemExit(f"minkowski.json carries no coordinate system {wanted!r}")


def check_system(report, entry, geometry, transformed):
    where = entry["id"]
    coords = entry["coords"]
    report.note(f"\n{where}: {entry['name']}")

    published = metric_from_line_element(entry["line_element"], coords)
    for i in range(len(coords)):
        for j in range(len(coords)):
            report.equal(
                published[i, j],
                transformed[i, j],
                f"{where}.line_element [{coords[i]},{coords[j]}] against the transformation from Cartesian",
            )
    report.note("  line element reproduced by the transformation from Cartesian")

    check_rank2(report, entry["metric_components"], geometry.g, coords, f"{where}.metric_components")
    check_rank2(report, entry["inverse_metric_components"], geometry.ginv, coords, f"{where}.inverse_metric_components")
    product = sp.simplify(geometry.g * geometry.ginv)
    report.check(product == sp.eye(len(coords)), f"{where}: g.ginv is {product}, not the identity")
    report.note("  metric, inverse metric and g^{mu alpha} g_{alpha nu} = delta^mu_nu")

    check_rank3(
        report,
        entry["christoffel"]["variants"]["ull"]["nonzero"],
        geometry.christoffel_ull,
        coords,
        f"{where}.christoffel.ull",
    )
    check_rank3(
        report,
        entry["christoffel"]["variants"]["lll"]["nonzero"],
        geometry.christoffel_lll,
        coords,
        f"{where}.christoffel.lll",
    )
    report.note("  Christoffel symbols in both variants")

    check_all_vanish(report, entry["riemann"]["variants"], f"{where}.riemann")
    check_all_vanish(report, entry["ricci_tensor"]["variants"], f"{where}.ricci_tensor")
    check_all_vanish(report, entry["einstein_tensor"]["variants"], f"{where}.einstein_tensor")
    check_all_vanish(report, entry["weyl_tensor"]["variants"], f"{where}.weyl_tensor")
    for label, tensor in (
        ("riemann", geometry.riemann_ulll),
        ("weyl", geometry.weyl_ulll),
    ):
        for a in range(len(coords)):
            for b in range(len(coords)):
                for d in range(len(coords)):
                    for e in range(len(coords)):
                        report.check(
                            sp.simplify(tensor[a][b][d][e]) == 0,
                            f"{where}: sympy makes {label} component "
                            f"[{coords[a]},{coords[b]},{coords[d]},{coords[e]}] nonzero",
                        )
    for a in range(len(coords)):
        for b in range(len(coords)):
            report.check(
                sp.simplify(geometry.ricci_ll[a][b]) == 0,
                f"{where}: sympy makes Ricci [{coords[a]},{coords[b]}] nonzero",
            )
            report.check(
                sp.simplify(geometry.einstein_ll[a][b]) == 0,
                f"{where}: sympy makes Einstein [{coords[a]},{coords[b]}] nonzero",
            )
    report.note("  Riemann, Ricci, Einstein and Weyl vanish in sympy and are published empty")

    report.equal(to_sympy(entry["ricci_scalar"].split("=")[1]), geometry.ricci_scalar, f"{where}.ricci_scalar")
    report.equal(to_sympy(entry["kretschmann"].split("=")[1]), geometry.kretschmann, f"{where}.kretschmann")
    report.note("  Ricci scalar and Kretschmann scalar")

    check_geodesics(report, entry["geodesics"], geometry, coords, f"{where}.geodesics")
    report.note("  geodesic equations")


def main():
    metric = json.loads(METRIC_FILE.read_text(encoding="utf-8"))
    report = Report()

    cartesian = system(metric, "cartesian")
    eta = metric_from_line_element(cartesian["line_element"], cartesian["coords"])
    report.equal(eta[0, 0], -c ** 2, "cartesian line element g_tt")
    report.note(f"Cartesian seed read from the file: diag{tuple(eta[i, i] for i in range(4))}")

    half = sp.Rational(1, 2)
    radius = c * (v - u) / 2

    double_null = pull_back(
        eta,
        [t, x, y, z],
        [u, v, y, z],
        {t: (u + v) * half, x: radius, y: y, z: z},
    )
    spherical_null = pull_back(
        eta,
        [t, x, y, z],
        [u, v, theta, phi],
        {
            t: (u + v) * half,
            x: radius * sp.sin(theta) * sp.cos(phi),
            y: radius * sp.sin(theta) * sp.sin(phi),
            z: radius * sp.cos(theta),
        },
    )

    spherical = system(metric, "spherical")
    g_spherical = metric_from_line_element(spherical["line_element"], spherical["coords"])
    by_way_of_spherical = pull_back(
        g_spherical,
        [t, r, theta, phi],
        [u, v, theta, phi],
        {t: (u + v) * half, r: radius, theta: theta, phi: phi},
    )
    report.check(
        sp.simplify(by_way_of_spherical - spherical_null) == sp.zeros(4, 4),
        "the route through the Spherical entry disagrees with the route straight from Cartesian",
    )
    report.note("The two routes to the Spherical Null metric agree")

    check_system(report, system(metric, "double_null"), Geometry(double_null, [u, v, y, z]), double_null)
    check_system(
        report,
        system(metric, "spherical_null"),
        Geometry(spherical_null, [u, v, theta, phi]),
        spherical_null,
    )

    if report.failures:
        print(f"\n{len(report.failures)} disagreements:", file=sys.stderr)
        for failure in report.failures:
            print(f"  {failure}", file=sys.stderr)
        return 1
    print("\nEverything the two null systems print agrees with sympy.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
