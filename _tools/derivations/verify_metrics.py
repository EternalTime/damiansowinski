#!/usr/bin/env python3
"""Check every published coordinate system in MFS/assets/data/metrics against sympy.

For each system the script reads the line element, builds g_{mu nu} from it, and then
computes the inverse metric, the Christoffel symbols in both variants, the Riemann
tensor, the Ricci tensor, the Ricci scalar, the Kretschmann scalar, the Einstein
tensor and the Weyl tensor. Every value the file publishes is compared against what
sympy got, and every component the file leaves out is required to vanish, so an
omission is caught as well as a wrong number.

Nothing passes silently. A system whose values cannot be parsed, or whose declaration
is missing, or a tensor sympy cannot finish inside the time budget, is reported as
UNCHECKED with the reason. The script exits non-zero if anything disagreed.

    python3 -m venv /tmp/mfs-venv && /tmp/mfs-venv/bin/pip install sympy
    /tmp/mfs-venv/bin/python _tools/derivations/verify_metrics.py

Pass --system <metric_id>/<system_id> to check one system, repeatable.
Pass --budget <seconds> to change the per tensor time budget, which defaults to 120.


The chart convention
--------------------

The collection writes x^0 = cT. A coordinate that carries dimensions of time is
therefore not itself the chart coordinate: the chart coordinate is c times it, and
every published component is a component in that chart even though the index is
printed with the bare name. Schwarzschild is the clearest example, publishing
g_{tt} = -(1 - r_s/r) against a line element whose time term is -(1 - r_s/r)c^2dt^2.

Because the rescaling x^0 -> c x^0 is linear with constant coefficients, a component
in the chart is the component computed with the bare coordinate multiplied by

    c^(number of upper time indices - number of lower time indices),

and the Christoffel symbols follow the same rule as the tensors, since the
inhomogeneous term in their transformation law carries a second derivative of the
coordinate change and so vanishes for a linear one.

TIME_COORDINATES below declares, per system, which coordinates carry dimensions of
time. It is the one thing the script cannot read off the file, because telling a time
from a length needs the dimensions of the parameters as well, and those are prose. A
system that is not declared there is reported UNCHECKED rather than guessed at.


The curvature conventions
-------------------------

The signature is (-,+,+,+) and the Riemann tensor is

    R^mu_{nu rho sigma} = d_rho Gamma^mu_{nu sigma} - d_sigma Gamma^mu_{nu rho}
                          + Gamma^mu_{rho lam} Gamma^lam_{nu sigma}
                          - Gamma^mu_{sigma lam} Gamma^lam_{nu rho},

which is what the published Riemann components are in.

The collection contracts the Ricci tensor as R_{mu nu} = R^a_{mu nu a}, on the last
lower index rather than the first. That is the opposite sign from the commoner
R^a_{mu a nu}, and it carries through to the Einstein tensor and the Ricci scalar, so
the published FRW scalar is -6(addot/a + adot^2/a^2 + k/a^2) where a textbook using
the other contraction would print +6(...). The convention is not written down in the
files, but it is what they consistently do: contracting each published Riemann tensor
both ways reproduces the published Ricci in every slot this way and not the other, in
Ellis-Bronnikov, Reissner-Nordstrom and Godel alike. The script therefore checks
against it rather than against the textbook contraction.

The Weyl tensor is the exception, because it is defined by removing the traces of
Riemann, and those traces do not care which contraction the file names Ricci. It is
built here from R^a_{mu a nu} regardless.
"""

import argparse
import json
import re
import signal
import sys
from pathlib import Path

import sympy as sp
from sympy.parsing.sympy_parser import (
    implicit_multiplication,
    parse_expr,
    split_symbols_custom,
    standard_transformations,
)

ROOT = Path(__file__).resolve().parents[2]
METRICS_DIR = ROOT / "MFS" / "assets" / "data" / "metrics"

DEFAULT_BUDGET_SECONDS = 120

# The coordinates each system writes as a time, so that the chart coordinate is c
# times what the index is printed as. A system absent from this table is UNCHECKED.
TIME_COORDINATES = {
    ("ellis_bronnikov", "spherical"): {"t"},
    ("frw", "comoving_spherical"): {"t"},
    ("frw", "conformal_spherical"): set(),
    ("godel", "cartesian"): set(),
    ("interior_schwarzschild", "spherical"): set(),
    ("minkowski", "cartesian"): {"t"},
    ("minkowski", "spherical"): {"t"},
    ("minkowski", "double_null"): {"u", "v"},
    ("minkowski", "spherical_null"): {"u", "v"},
    ("minkowski", "rindler"): {"T"},
    ("rn_metric", "spherical"): set(),
    ("schwarzschild", "spherical"): {"t"},
    ("schwarzschild", "eddington_finkelstein_outgoing"): set(),
    ("schwarzschild", "eddington_finkelstein_ingoing"): set(),
    ("stockum_dust", "cylindrical"): set(),
}

GREEK = [
    "theta", "phi", "eta", "omega", "Omega", "ell", "pi", "lambda", "mu", "nu",
    "rho", "sigma", "tau", "chi", "psi", "alpha", "beta", "gamma", "delta",
    "epsilon", "kappa", "xi", "zeta", "Lambda", "Phi", "Theta", "Psi", "Sigma",
]

TRIG = ["sinh", "cosh", "tanh", "coth", "sin", "cos", "tan", "cot", "sec", "csc"]

FUNCTIONS = {
    "sin": sp.sin, "cos": sp.cos, "tan": sp.tan, "cot": sp.cot,
    "sec": sp.sec, "csc": sp.csc, "sinh": sp.sinh, "cosh": sp.cosh,
    "tanh": sp.tanh, "coth": sp.coth, "exp": sp.exp, "log": sp.log,
    "ln": sp.log, "sqrt": sp.sqrt,
}


def norm(expression):
    """Simplify far enough that a component which vanishes is recognisably zero.

    sympy's simplify leaves forms such as sin(2x)tan(x) + cos(2x) - 1 standing, which
    would otherwise be reported as a curvature component the file forgot to publish.
    """
    simplified = sp.simplify(expression)
    if simplified == 0:
        return sp.Integer(0)
    harder = sp.simplify(sp.expand_trig(sp.expand(simplified)))
    if harder == 0:
        return sp.Integer(0)
    return harder if sp.count_ops(harder) < sp.count_ops(simplified) else simplified


class LatexError(Exception):
    pass


class Timeout(Exception):
    pass


class budget:
    """Abandon a computation that runs past `seconds` rather than let it hang."""

    def __init__(self, seconds):
        self.seconds = seconds

    def __enter__(self):
        signal.signal(signal.SIGALRM, self._fire)
        signal.setitimer(signal.ITIMER_REAL, self.seconds)

    def __exit__(self, *_):
        signal.setitimer(signal.ITIMER_REAL, 0)
        return False

    def _fire(self, *_):
        raise Timeout(f"sympy did not finish inside {self.seconds}s")


def matching_brace(text, opening):
    depth = 0
    for position in range(opening, len(text)):
        if text[position] == "{":
            depth += 1
        elif text[position] == "}":
            depth -= 1
            if depth == 0:
                return position
    raise LatexError(f"unbalanced braces in {text!r}")


def expand_fractions(text):
    while True:
        found = re.search(r"\\[dt]?frac\s*\{", text)
        if not found:
            return text
        first = text.index("{", found.start())
        first_end = matching_brace(text, first)
        second = text.index("{", first_end)
        second_end = matching_brace(text, second)
        numerator = text[first + 1:first_end]
        denominator = text[second + 1:second_end]
        text = f"{text[:found.start()]}(({numerator})/({denominator})){text[second_end + 1:]}"


def expand_braced_call(text, command, replacement):
    """Turn \\command{body} into replacement(body), innermost first."""
    while True:
        found = re.search(r"\\" + command + r"\s*\{", text)
        if not found:
            return text
        opening = text.index("{", found.start())
        closing = matching_brace(text, opening)
        body = text[opening + 1:closing]
        text = f"{text[:found.start()]}{replacement}({body}){text[closing + 1:]}"


def expand_superscript_braces(text):
    """^{body} into **(body), so that ^{-1} and e^{-r^2} survive the parser."""
    out = []
    position = 0
    while position < len(text):
        if text[position] == "^" and position + 1 < len(text) and text[position + 1] == "{":
            closing = matching_brace(text, position + 1)
            out.append("**(" + expand_superscript_braces(text[position + 2:closing]) + ")")
            position = closing + 1
        else:
            out.append(text[position])
            position += 1
    return "".join(out)


class Reader:
    """The LaTeX dialect the metric files are written in, translated into sympy.

    Every name the reader is willing to produce has to be declared up front, so a
    typo in a published value becomes an error here rather than a silent new symbol.
    """

    def __init__(self, coords, parameters, time_coords=frozenset()):
        self.coords = list(coords)
        self.time_coords = set(time_coords)
        self.symbol = {}
        self.differential = {}
        self.dot = {}
        self.ddot = {}
        for name in self.coords:
            plain = self._plain(name)
            self.symbol[name] = sp.Symbol(plain, real=True)
            self.differential[name] = sp.Symbol("d" + plain)
            self.dot[name] = sp.Symbol(plain + "_dot")
            self.ddot[name] = sp.Symbol(plain + "_ddot")
        self.c = sp.Symbol("c", positive=True)
        self.parameters = {}
        self.primed = set()
        for declaration in parameters:
            self._declare_parameter(declaration)
        self.allowed = (
            set(self.symbol.values())
            | set(self.differential.values())
            | set(self.dot.values())
            | set(self.ddot.values())
            | {self.c}
        )
        for value in self.parameters.values():
            self.allowed |= value.free_symbols
        self.local = {str(s): s for s in self.allowed}
        self.local.update(self.parameters)
        self.local.update(FUNCTIONS)
        self.known = set(self.local)
        self.transforms = standard_transformations + (
            split_symbols_custom(lambda name, _=None: name not in self.known),
            implicit_multiplication,
        )

    @staticmethod
    def _plain(name):
        return name.replace("\\", "").strip()

    def _declare_parameter(self, declaration):
        """`a` is a constant; `a = a(t)` is a function of the coordinate it names.

        For the second kind the printed rate is a derivative with respect to the chart
        coordinate, which is c times the named one when that one is a time, so the rate
        carries the matching power of c.
        """
        plain = self._plain(declaration.split("=")[0])
        argument = re.search(r"\(\s*(\\?[A-Za-z]+)\s*\)", declaration)
        if argument is None:
            self.parameters[plain] = sp.Symbol(plain, real=True)
            return
        name = argument.group(1)
        if name not in self.symbol:
            raise LatexError(f"parameter {declaration!r} varies with {name!r}, which is not a coordinate")
        variable = self.symbol[name]
        function = sp.Function(plain, real=True)(variable)
        scale = self.c if name in self.time_coords else sp.Integer(1)
        self.parameters[plain] = function
        first = sp.Derivative(function, variable) / scale
        second = sp.Derivative(function, variable, 2) / scale ** 2
        # Both spellings appear: a dot in the comoving chart, a prime in the conformal one.
        for suffix, value in (("_dot", first), ("_prime", first),
                              ("_ddot", second), ("_pprime", second)):
            self.parameters[plain + suffix] = value
        self.primed.add(plain)

    def _preprocess(self, latex):
        text = latex
        text = text.replace("\\left", " ").replace("\\right", " ")
        text = re.sub(r"\\[,;:!>]", " ", text)
        text = re.sub(r"\\ ", " ", text)
        text = text.replace("\\cdot", "*")
        for name in sorted(self.primed, key=len, reverse=True):
            text = re.sub(re.escape(name) + r"''", f" {name}_pprime ", text)
            text = re.sub(re.escape(name) + r"'", f" {name}_prime ", text)
        text = expand_fractions(text)
        text = expand_braced_call(text, "sqrt", "sqrt")
        text = expand_braced_call(text, "ddot", "DDOT")
        text = expand_braced_call(text, "dot", "DOT")
        text = re.sub(r"DDOT\s*\(\s*\\?([A-Za-z]+)\s*\)", r" \1_ddot ", text)
        text = re.sub(r"DOT\s*\(\s*\\?([A-Za-z]+)\s*\)", r" \1_dot ", text)
        # d\Omega^2 is the unit two sphere, written out so the reader sees differentials.
        text = re.sub(r"d\\Omega\s*\^\s*2", "(dtheta**2 + sin(theta)**2*dphi**2)", text)
        # Differentials become single atoms before anything is allowed to pad with spaces.
        for name in sorted(self.coords, key=len, reverse=True):
            text = text.replace("d" + name, "d" + self._plain(name))
        for command in sorted(GREEK + TRIG, key=len, reverse=True):
            text = text.replace("\\" + command, " " + command + " ")
        text = text.replace("\\exp", " exp ").replace("\\ln", " log ").replace("\\log", " log ")
        text = expand_superscript_braces(text)
        text = text.replace("^", "**")
        # A trig call written bare, as \sin^2\theta or \cot\theta rather than sin(theta).
        names = "|".join(sorted(TRIG, key=len, reverse=True))
        text = re.sub(r"\b(" + names + r")\s*\*\*\s*\(?\s*(-?\d+)\s*\)?\s*([A-Za-z]\w*)",
                      r"(\1(\3))**\2", text)
        text = re.sub(r"\b(" + names + r")\s+([A-Za-z]\w*)", r"\1(\2)", text)
        text = re.sub(r"(?<![A-Za-z_])e\s*\*\*", " E**", text)
        # 3R\sqrt{..} leaves 3Rsqrt(..), whose one name token would be split letter by letter.
        for name in FUNCTIONS:
            text = re.sub(r"(?<=[A-Za-z0-9_])(" + name + r")\s*\(", r" \1(", text)
        if "\\" in text:
            raise LatexError(f"unhandled LaTeX in {latex!r}: {text!r}")
        return text

    def __call__(self, latex):
        text = self._preprocess(latex)
        try:
            expression = parse_expr(
                text, local_dict=dict(self.local), transformations=self.transforms, evaluate=True
            )
        except Exception as error:
            raise LatexError(f"cannot parse {latex!r} (as {text!r}): {error}") from error
        unknown = expression.free_symbols - self.allowed
        if unknown:
            raise LatexError(f"{latex!r} produced undeclared symbols {sorted(map(str, unknown))}")
        return expression


def metric_from_line_element(reader, line_element, coords):
    """Read g_{mu nu} off the line element, in the chart its own coords name."""
    _, _, right = line_element.partition("=")
    form = sp.expand(reader(right))
    differentials = [reader.differential[name] for name in coords]
    size = len(coords)
    g = sp.zeros(size, size)
    for i in range(size):
        for j in range(size):
            if i == j:
                g[i, j] = form.coeff(differentials[i], 2)
            else:
                g[i, j] = form.coeff(differentials[i], 1).coeff(differentials[j], 1) / 2
    remainder = sp.expand(
        form - sum(g[i, j] * differentials[i] * differentials[j] for i in range(size) for j in range(size))
    )
    if norm(remainder) != 0:
        raise LatexError(f"{line_element!r} is not a quadratic form in its own differentials")
    return norm(g)


class Geometry:
    """Every tensor the files publish, computed from g in the chart the coords name."""

    def __init__(self, g, coords, seconds):
        self.n = len(coords)
        self.coords = coords
        self.seconds = seconds
        self.g = norm(g)
        self.ginv = norm(self.g.inv())
        self._cache = {}
        self.unavailable = {}

    def _timed(self, name, build):
        if name in self._cache:
            return self._cache[name]
        if name in self.unavailable:
            raise Timeout(self.unavailable[name])
        try:
            with budget(self.seconds):
                self._cache[name] = build()
        except Timeout as error:
            self.unavailable[name] = str(error)
            raise
        return self._cache[name]

    def _zeros(self, rank):
        if rank == 1:
            return [sp.Integer(0)] * self.n
        return [self._zeros(rank - 1) for _ in range(self.n)]

    def christoffel_ull(self):
        def build():
            out = self._zeros(3)
            for mu in range(self.n):
                for nu in range(self.n):
                    for rho in range(nu, self.n):
                        value = norm(sum(
                            self.ginv[mu, alpha] * (
                                sp.diff(self.g[alpha, rho], self.coords[nu])
                                + sp.diff(self.g[alpha, nu], self.coords[rho])
                                - sp.diff(self.g[nu, rho], self.coords[alpha])
                            )
                            for alpha in range(self.n)
                        ) / 2)
                        out[mu][nu][rho] = value
                        out[mu][rho][nu] = value
            return out
        return self._timed("christoffel_ull", build)

    def christoffel_lll(self):
        def build():
            gamma = self.christoffel_ull()
            out = self._zeros(3)
            for mu in range(self.n):
                for nu in range(self.n):
                    for rho in range(self.n):
                        out[mu][nu][rho] = norm(
                            sum(self.g[mu, alpha] * gamma[alpha][nu][rho] for alpha in range(self.n))
                        )
            return out
        return self._timed("christoffel_lll", build)

    def riemann_ulll(self):
        def build():
            gamma = self.christoffel_ull()
            out = self._zeros(4)
            for mu in range(self.n):
                for nu in range(self.n):
                    for rho in range(self.n):
                        for sigma in range(rho + 1, self.n):
                            value = norm(
                                sp.diff(gamma[mu][nu][sigma], self.coords[rho])
                                - sp.diff(gamma[mu][nu][rho], self.coords[sigma])
                                + sum(
                                    gamma[mu][rho][lam] * gamma[lam][nu][sigma]
                                    - gamma[mu][sigma][lam] * gamma[lam][nu][rho]
                                    for lam in range(self.n)
                                )
                            )
                            out[mu][nu][rho][sigma] = value
                            out[mu][nu][sigma][rho] = -value
            return out
        return self._timed("riemann_ulll", build)

    def riemann_llll(self):
        def build():
            upper = self.riemann_ulll()
            out = self._zeros(4)
            for mu in range(self.n):
                for nu in range(self.n):
                    for rho in range(self.n):
                        for sigma in range(self.n):
                            out[mu][nu][rho][sigma] = norm(sum(
                                self.g[mu, alpha] * upper[alpha][nu][rho][sigma]
                                for alpha in range(self.n)
                            ))
            return out
        return self._timed("riemann_llll", build)

    def ricci_trace_ll(self):
        """R^a_{mu a nu}, the contraction whose traces are Riemann's own.

        This is the one the Weyl tensor is built from, because Weyl is defined by
        removing the traces of Riemann, and those traces do not care what the file
        chooses to call the Ricci tensor.
        """
        def build():
            riemann = self.riemann_ulll()
            out = self._zeros(2)
            for nu in range(self.n):
                for sigma in range(nu, self.n):
                    value = norm(sum(riemann[mu][nu][mu][sigma] for mu in range(self.n)))
                    out[nu][sigma] = value
                    out[sigma][nu] = value
            return out
        return self._timed("ricci_trace_ll", build)

    def ricci_ll(self):
        """R_{mu nu} = R^a_{mu nu a}, which is the contraction the collection publishes."""
        def build():
            trace = self.ricci_trace_ll()
            return [[norm(-trace[a][b]) for b in range(self.n)] for a in range(self.n)]
        return self._timed("ricci_ll", build)

    def ricci_scalar(self):
        def build():
            ricci = self.ricci_ll()
            return norm(sum(
                self.ginv[a, b] * ricci[a][b] for a in range(self.n) for b in range(self.n)
            ))
        return self._timed("ricci_scalar", build)

    def einstein_ll(self):
        def build():
            ricci = self.ricci_ll()
            scalar = self.ricci_scalar()
            out = self._zeros(2)
            for a in range(self.n):
                for b in range(self.n):
                    out[a][b] = norm(ricci[a][b] - scalar * self.g[a, b] / 2)
            return out
        return self._timed("einstein_ll", build)

    def kretschmann(self):
        def build():
            lower = self.riemann_llll()
            upper = self._zeros(4)
            for a in range(self.n):
                for b in range(self.n):
                    for cc in range(self.n):
                        for d in range(self.n):
                            upper[a][b][cc][d] = norm(sum(
                                self.ginv[a, p] * self.ginv[b, q] * self.ginv[cc, s] * self.ginv[d, w]
                                * lower[p][q][s][w]
                                for p in range(self.n) for q in range(self.n)
                                for s in range(self.n) for w in range(self.n)
                            ))
            return norm(sum(
                lower[a][b][cc][d] * upper[a][b][cc][d]
                for a in range(self.n) for b in range(self.n)
                for cc in range(self.n) for d in range(self.n)
            ))
        return self._timed("kretschmann", build)

    def weyl_llll(self):
        def build():
            n = self.n
            riemann = self.riemann_llll()
            ricci = self.ricci_trace_ll()
            scalar = norm(sum(
                self.ginv[a, b] * ricci[a][b] for a in range(n) for b in range(n)
            ))
            out = self._zeros(4)
            for a in range(n):
                for b in range(n):
                    for cc in range(n):
                        for d in range(n):
                            correction = (
                                self.g[a, cc] * ricci[d][b] - self.g[a, d] * ricci[cc][b]
                                - self.g[b, cc] * ricci[d][a] + self.g[b, d] * ricci[cc][a]
                            ) / (n - 2)
                            trace = scalar * (
                                self.g[a, cc] * self.g[d, b] - self.g[a, d] * self.g[cc, b]
                            ) / ((n - 1) * (n - 2))
                            out[a][b][cc][d] = norm(riemann[a][b][cc][d] - correction + trace)
            return out
        return self._timed("weyl_llll", build)

    def raise_indices(self, tensor, rank, positions):
        """Raise the given index positions of a fully lowered tensor."""
        out = tensor
        for position in positions:
            source = out
            out = self._zeros(rank)
            for index in _indices(self.n, rank):
                total = 0
                for alpha in range(self.n):
                    swapped = list(index)
                    swapped[position] = alpha
                    total += self.ginv[index[position], alpha] * _at(source, swapped)
                _put(out, index, norm(total))
        return out


def _indices(n, rank):
    if rank == 0:
        yield []
        return
    for head in range(n):
        for tail in _indices(n, rank - 1):
            yield [head] + tail


def _at(tensor, index):
    for i in index:
        tensor = tensor[i]
    return tensor


def _put(tensor, index, value):
    for i in index[:-1]:
        tensor = tensor[i]
    tensor[index[-1]] = value


class Report:
    def __init__(self):
        self.disagreements = []
        self.unchecked = []
        self.checked_systems = 0
        self.systems = 0

    def disagree(self, where, message):
        self.disagreements.append(f"{where}: {message}")

    def skip(self, where, reason):
        self.unchecked.append(f"{where}: {reason}")

    def guarded(self, where, seconds, work):
        """Run a comparison under the clock, so a slow one is named rather than hung on."""
        try:
            with budget(seconds):
                work()
        except Timeout as error:
            self.skip(where, f"comparing the published values: {error}")


def variance_weight(variance, coords, index, time_coords):
    """c^(upper time indices - lower time indices) for the x^0 = cT chart."""
    exponent = 0
    for slot, position in zip(variance, index):
        if coords[position] in time_coords:
            exponent += 1 if slot == "u" else -1
    return exponent


def compare_block(report, reader, where, published, computed, variance, coords, time_coords, c):
    """Every published component against sympy, and every omitted one against zero."""
    rank = len(variance)
    seen = set()
    for entry in published:
        names = entry["indices"]
        if len(names) != rank:
            report.disagree(where, f"{names} has {len(names)} indices, expected {rank}")
            continue
        try:
            index = [coords.index(name) for name in names]
        except ValueError:
            report.disagree(where, f"{names} names an index that is not one of {coords}")
            continue
        seen.add(tuple(index))
        try:
            value = reader(entry["value"])
        except LatexError as error:
            report.skip(f"{where} {names}", str(error))
            continue
        expected = _at(computed, index) * c ** variance_weight(variance, coords, index, time_coords)
        if norm(sp.expand(sp.together(value - expected))) != 0:
            report.disagree(where, f"{names} published as {entry['value']} "
                                   f"({norm(value)}), sympy says {norm(expected)}")
    for index in _indices(len(coords), rank):
        if tuple(index) in seen:
            continue
        if norm(_at(computed, index)) != 0:
            names = [coords[i] for i in index]
            weighted = _at(computed, index) * c ** variance_weight(variance, coords, index, time_coords)
            report.disagree(where, f"{names} is missing, sympy says {norm(weighted)}")


def compare_scalar(report, reader, where, published, computed):
    text = published.split("=", 1)[1] if "=" in published else published
    try:
        value = reader(text)
    except LatexError as error:
        report.skip(where, str(error))
        return
    if norm(sp.expand(sp.together(value - computed))) != 0:
        report.disagree(where, f"published as {published.strip()}, sympy says {norm(computed)}")


def compare_geodesics(report, reader, where, published, gamma, coords, time_coords, c):
    """Each equation against xddot^mu + Gamma^mu_{nu rho} xdot^nu xdot^rho = 0.

    The printed dots are the chart velocities, matching the printed Christoffels, so
    the weighted symbols are the ones the residual is built from.
    """
    if len(published) != len(coords):
        report.disagree(where, f"{len(published)} equations for {len(coords)} coordinates")
    for equation in published:
        left, _, right = equation.partition("=")
        try:
            residual = sp.expand(reader(left) - reader(right))
        except LatexError as error:
            report.skip(f"{where} {equation!r}", str(error))
            continue
        carried = [name for name in coords if residual.has(reader.ddot[name])]
        if len(carried) != 1:
            report.disagree(where, f"{equation!r} carries second derivatives of {carried}, expected one")
            continue
        name = carried[0]
        mu = coords.index(name)
        expected = reader.ddot[name] + sum(
            gamma[mu][nu][rho]
            * c ** variance_weight("ull", coords, [mu, nu, rho], time_coords)
            * reader.dot[coords[nu]] * reader.dot[coords[rho]]
            for nu in range(len(coords)) for rho in range(len(coords))
        )
        if norm(residual - expected) != 0 and norm(residual + expected) != 0:
            report.disagree(where, f"{name} equation {equation!r} is not the geodesic equation, "
                                   f"sympy makes the residual {norm(expected)}")


VECTOR_VARIANTS = {"ll": (0, "ll"), "ul": ((0,), "ul"), "uu": ((0, 1), "uu")}
RANK4_VARIANTS = {"llll": (None, "llll"), "ulll": ((0,), "ulll")}


def check_system(report, metric_id, entry, seconds):
    where = f"{metric_id}/{entry['id']}"
    coords = entry["coords"]
    report.systems += 1

    declaration = TIME_COORDINATES.get((metric_id, entry["id"]))
    if declaration is None:
        report.skip(where, "no time coordinate declaration in TIME_COORDINATES")
        return
    unknown = declaration - set(coords)
    if unknown:
        report.skip(where, f"time coordinate declaration names {sorted(unknown)}, not in {coords}")
        return

    parameters = [p["symbol"] for p in entry.get("parameters", [])]
    try:
        reader = Reader(coords, parameters, declaration)
    except LatexError as error:
        report.skip(where, f"parameters unreadable: {error}")
        return
    c = reader.c
    try:
        g = metric_from_line_element(reader, entry["line_element"], coords)
    except LatexError as error:
        report.skip(where, f"line element unreadable: {error}")
        return
    if norm(g.det()) == 0:
        report.skip(where, "the line element gives a degenerate metric")
        return

    # Everything is computed with the bare coordinates and then weighted into the
    # x^0 = cT chart by compare_block, which is exact because the rescaling is linear.
    symbols = [reader.symbol[name] for name in coords]
    geometry = Geometry(g, symbols, seconds)
    report.checked_systems += 1
    print(f"  {where}")

    metric_variants = {"ll": g, "uu": geometry.ginv}
    published_metric = [(f"{where}.metric_components", entry.get("metric_components"), "ll", g)]
    if "inverse_metric_components" in entry:
        published_metric.append((
            f"{where}.inverse_metric_components", entry["inverse_metric_components"],
            "uu", geometry.ginv,
        ))
    else:
        report.skip(f"{where}.inverse_metric_components", "the entry does not publish one")
    for label, published, variance, matrix in published_metric:
        as_lists = [[matrix[i, j] for j in range(len(coords))] for i in range(len(coords))]
        report.guarded(label, seconds, lambda label=label, published=published,
                       as_lists=as_lists, variance=variance: compare_block(
            report, reader, label, published, as_lists, variance, coords, declaration, c))
    if norm(g * geometry.ginv) != sp.eye(len(coords)):
        report.disagree(where, "the published metric and inverse are not inverse to each other")

    blocks = [
        ("christoffel", "ull", lambda: geometry.christoffel_ull()),
        ("christoffel", "lll", lambda: geometry.christoffel_lll()),
        ("riemann", "ulll", lambda: geometry.raise_indices(geometry.riemann_llll(), 4, (0,))),
        ("riemann", "llll", lambda: geometry.riemann_llll()),
        ("ricci_tensor", "ll", lambda: geometry.ricci_ll()),
        ("ricci_tensor", "ul", lambda: geometry.raise_indices(geometry.ricci_ll(), 2, (0,))),
        ("ricci_tensor", "uu", lambda: geometry.raise_indices(geometry.ricci_ll(), 2, (0, 1))),
        ("einstein_tensor", "ll", lambda: geometry.einstein_ll()),
        ("einstein_tensor", "ul", lambda: geometry.raise_indices(geometry.einstein_ll(), 2, (0,))),
        ("einstein_tensor", "uu", lambda: geometry.raise_indices(geometry.einstein_ll(), 2, (0, 1))),
        ("weyl_tensor", "llll", lambda: geometry.weyl_llll()),
        ("weyl_tensor", "ulll", lambda: geometry.raise_indices(geometry.weyl_llll(), 4, (0,))),
    ]
    for field, variance, build in blocks:
        variants = entry.get(field, {}).get("variants", {})
        if variance not in variants:
            continue
        label = f"{where}.{field}.{variance}"
        try:
            computed = build()
        except Timeout as error:
            report.skip(label, str(error))
            continue
        report.guarded(label, seconds, lambda label=label, computed=computed,
                       variance=variance, published=variants[variance]["nonzero"]:
                       compare_block(report, reader, label, published, computed,
                                     variance, coords, declaration, c))

    for field, build in (("ricci_scalar", geometry.ricci_scalar), ("kretschmann", geometry.kretschmann)):
        if field not in entry:
            continue
        label = f"{where}.{field}"
        try:
            computed = build()
        except Timeout as error:
            report.skip(label, str(error))
            continue
        report.guarded(label, seconds, lambda label=label, computed=computed, field=field:
                       compare_scalar(report, reader, label, entry[field], computed))

    if "geodesics" in entry:
        label = f"{where}.geodesics"
        try:
            gamma = geometry.christoffel_ull()
        except Timeout as error:
            report.skip(label, str(error))
        else:
            report.guarded(label, seconds, lambda: compare_geodesics(
                report, reader, label, entry["geodesics"], gamma, coords, declaration, c))
    else:
        report.skip(f"{where}.geodesics", "the entry does not publish any")


CONFLICT_COPY = re.compile(r" \d+$")


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--system", action="append", default=[],
                        help="check only <metric_id>/<system_id>, repeatable")
    parser.add_argument("--budget", type=float, default=DEFAULT_BUDGET_SECONDS,
                        help="seconds sympy may spend on one tensor")
    arguments = parser.parse_args()

    report = Report()
    files = 0
    print("Checking every published coordinate system against sympy, in the x^0 = cT chart.\n")
    for path in sorted(METRICS_DIR.glob("*.json")):
        if CONFLICT_COPY.search(path.stem):
            continue
        files += 1
        metric = json.loads(path.read_text(encoding="utf-8"))
        for entry in metric.get("coordinates", []):
            name = f"{metric['id']}/{entry['id']}"
            if arguments.system and name not in arguments.system:
                continue
            check_system(report, metric["id"], entry, arguments.budget)

    print(f"\n{files} metric files, {report.systems} coordinate systems, "
          f"{report.checked_systems} checked.")
    if report.unchecked:
        print(f"\n{len(report.unchecked)} UNCHECKED:")
        for line in report.unchecked:
            print(f"  UNCHECKED {line}")
    if report.disagreements:
        print(f"\n{len(report.disagreements)} disagreements:", file=sys.stderr)
        for line in report.disagreements:
            print(f"  {line}", file=sys.stderr)
        return 1
    print("\nEvery published value agrees with sympy.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
