#!/usr/bin/env python3
r"""Hold every symbol a spacetime uses to a place that defines it.

    python3 _tools/symbols.py [spacetime ...]

A reader meets a symbol in a chart's mathematics, on a drawing or in a caption, and looks for what
it means. This check refuses a symbol that nothing the reader can see defines. It prints each one
with the place that uses it and exits non-zero if there is any; `build_mfs_data.py` holds the same
rule and refuses to write while it fails, and `_tools/README.md` says how to mend a finding.

What a chart defines:

- its `coords`;
- the `symbol` of each of its `parameters`, and what a parameter's description defines in one of
  the three ways below;
- every symbol in the mathematics of its `convention` and of the spacetime's shared `convention`;
- what one of its `domains` names in its words, as $N$ is in "(the closed null geodesic $N$)".

What is held to those definitions:

- a chart's `domains`, `line_element`, components, curvature and `geodesics`, held to that chart
  alone, since the page shows one chart at a time;
- every text of every spacetime diagram, conformal diagram and embedding diagram, its labels,
  axes, legends, settings, inputs, restrictions and captions, held to all the spacetime's charts
  together, since a caption may set one chart beside another, as in "the ingoing chart's $g_{vr}$".

A drawing may also define a symbol of its own, one that belongs to the picture and not to the
metric, such as the scale of a conformal map. Any of its texts may do so, in one of three ways:

- by an equation with the symbol alone on its left, "$r_* = r + r_s\ln(r/r_s - 1)$";
- by a noun of NAMING right before it, "the height $z$";
- by saying what it is, "with $\ell$ any length".

`defined_in` gives each of the three in full. A value alone, "$a = 1$", says how much and not what,
so it defines a symbol only in a drawing's `settings` and `input`, the two lines that say what was
drawn.

Histories, descriptions and related entries are prose about other people's notation and are not
read; each defines what it uses in its own sentences.

A symbol is a letter, Latin or Greek, with its subscript and any accent that changes its meaning:
$r$, $r_s$, $\tilde\phi$ and $r_*$ are four symbols. A dot is a derivative and a power is a power,
so neither makes a new symbol. A subscript made only of coordinates or of index letters marks a
component, as in $g_{tt}$ or $u_\mu$, and the symbol is the letter that carries it.

UNIVERSAL, DRAWN and TENSORS below are the whole list of symbols that need no definition. Keep
them short: a symbol belongs there only if it means the same thing on every page of the collection.
"""

import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "MFS" / "assets" / "data"

# Symbols that mean the same thing in every spacetime and are defined in none.
UNIVERSAL = {
    "c": "the speed of light",
    "G": "Newton's constant",
    "\\hbar": "Planck's constant",
    "k_B": "Boltzmann's constant",
    "M_\\odot": "the mass of the Sun",
    "\\pi": "the number",
    "e": "the base of the natural logarithm",
    "i": "the imaginary unit, and with a superscript the points at infinity $i^0$ and $i^\\pm$",
    "d": "the differential",
    "s": "the interval, in $ds^2$",
    "\\partial": "a partial derivative",
    "\\nabla": "the covariant derivative",
    "\\infty": "infinity",
    "g": "the metric",
    "\\Gamma": "the Christoffel symbols",
    "\\mathscr{I}": "null infinity",
    "\\mathbb{R}": "the real numbers",
    "\\mathbb{Z}": "the integers",
    "S^1": "the circle",
    "S^2": "the 2-sphere",
    "S^3": "the 3-sphere",
}
# Two more are universal by the letter before or after them: $d\Omega^2$ is the metric of the unit
# 2-sphere, and $\Delta v$ is a change in $v$. Standing alone, $\Omega$ and $\Delta$ need defining.
SOLID_ANGLE, CHANGE = "\\Omega", "\\Delta"
# The coordinates of a drawing's own chart, which every drawing of its kind shares.
DRAWN = {
    "conformal": {
        "T": "the diagram's coordinate up the page",
        "X": "the diagram's coordinate across the page",
        "p": "the diagram's null coordinate, constant along the rays that run up to the right",
        "q": "the diagram's null coordinate, constant along the rays that run up to the left",
    },
    "embedding": {
        "X": "a Cartesian coordinate of the flat space the surface stands in",
        "Y": "a Cartesian coordinate of the flat space the surface stands in",
        "Z": "the Cartesian coordinate of that space along the axis, where the space is Minkowski's",
        "z": "the height of the surface, where the space is Euclid's",
    },
}
# Symbols that are universal only while they carry indices: the curvature tensors, the Einstein
# tensor, the stress-energy tensor, the extrinsic curvature, the metric, the Christoffel symbols,
# the Kronecker delta, the flat metric and the alternating tensor. Bare, each is an ordinary letter.
TENSORS = {"R", "G", "C", "T", "K", "g", "\\Gamma", "\\delta", "\\eta", "\\epsilon", "\\varepsilon"}
# Letters that are indices when a subscript or a superscript holds nothing else.
INDICES = {"\\mu", "\\nu", "\\rho", "\\sigma", "\\alpha", "\\beta", "\\gamma", "\\delta", "\\lambda", "\\kappa",
           "a", "b", "i", "j", "k"}
# The scalars a chart's own fields name on the left of their equals sign.
SCALAR_FIELDS = {"ricci_scalar": "R", "kretschmann": "K"}

# The nouns that name a symbol standing alone right after them, as in "the height $z$" or "for
# integer $k$". A caption that names a symbol with another noun adds the noun here.
NAMING = (
    "angle", "charge", "circle", "coordinate", "coordinates", "cosmological constant", "curve", "density",
    "distance", "event", "events", "factor", "function", "functions", "geodesic", "height", "integer", "integers",
    "length", "map", "maps", "mass", "number", "parameter", "point", "points", "pressure", "radius", "rate", "scale",
    "speed", "tension", "time", "unit", "width",
)

GREEK = {
    "alpha", "beta", "gamma", "delta", "epsilon", "varepsilon", "zeta", "eta", "theta", "vartheta", "iota", "kappa",
    "lambda", "mu", "nu", "xi", "pi", "varpi", "rho", "varrho", "sigma", "varsigma", "tau", "upsilon", "phi",
    "varphi", "chi", "psi", "omega", "Gamma", "Delta", "Theta", "Lambda", "Xi", "Pi", "Sigma", "Upsilon", "Phi",
    "Psi", "Omega", "ell", "hbar", "partial", "nabla", "infty",
}
# Accents that leave the symbol what it was: derivatives along a curve.
TRANSPARENT = {"dot", "ddot"}
# Accents and alphabets that make a new symbol of the letter they hold.
MARKS = {"tilde", "widetilde", "bar", "overline", "hat", "widehat", "vec", "check", "mathcal", "mathscr", "mathbb",
         "mathbf", "boldsymbol"}
# Commands whose argument is words, a unit or the name of a function, never a symbol.
WORDS = {"text", "mathrm", "operatorname", "textit", "mathit"}
# The old switch to upright letters, which runs to the end of its group.
UPRIGHT = re.compile(r"\\rm\b[^{}]*")

MATH = re.compile(r"\$\$(.+?)\$\$|\$([^$]+)\$", re.DOTALL)
TOKEN = re.compile(r"\\[A-Za-z]+|\\.|[A-Za-z]|[_^{}']|[^\\A-Za-z_^{}'\s]")
# A list of mathematics in a sentence: "$U$", "$U$ and $V$", "$Z_0$, $Z_1$, and $Z_2$".
LIST = r"\$[^$]+\$(?:(?:, | and |, and )\$[^$]+\$)*"
NAMED = re.compile(rf"\b(?:{'|'.join(NAMING)}) ({LIST})")
SAID = re.compile(rf"({LIST}),? (?:is |are |being )?(?:the|a|an|any|its|their|his|her) ")
EQUALS = re.compile(r"(?<![<>!\\])=")
RELATION = re.compile(r"\\(?:equiv|sim|approx|le|ge|leq|geq|ne|neq|in|to)\b|[<>]")
ARGUMENTS = re.compile(r"\([^()]*\)$")
LIGHT_SPEED = re.compile(r"^c(?:\\,)?\s*(?=\\|[A-Za-z])")
DIFFERENTIAL = re.compile(r"^d(?=\\|[A-Za-z])([^/]*?)(?:/d.+)?$")
POWER = re.compile(r"\^(?:\d|\{\d+\})$")
NOT_A_SYMBOL = re.compile(r"_\d|_\{[^{}]*\}|\\[A-Za-z]+|[A-Za-z_^{}'+\-*,\s]|\\[,;!]")
WORDED = re.compile(r"\\text\{([^{}]*)\}")


def mathematics(text):
    """The mathematics of a sentence: what stands between its dollar signs."""
    return [a or b for a, b in MATH.findall(text)]


class Reader:
    """Read the symbols out of one TeX expression."""

    def __init__(self, tex="", tokens=None):
        self.tokens = TOKEN.findall(UPRIGHT.sub("", tex)) if tokens is None else tokens
        self.at = 0
        self.found = []

    def peek(self):
        return self.tokens[self.at] if self.at < len(self.tokens) else None

    def take(self):
        token = self.peek()
        self.at += 1
        return token

    def group(self):
        """The tokens of the next argument: a braced group without its braces, or one token."""
        if self.peek() != "{":
            return [self.take()] if self.peek() is not None else []
        self.take()
        depth, inside = 1, []
        while self.peek() is not None:
            token = self.take()
            depth += (token == "{") - (token == "}")
            if depth == 0:
                break
            inside.append(token)
        return inside

    def letter(self):
        """The next symbol's letter with its accents, or None once something else is passed over."""
        token = self.take()
        name = token[1:] if token.startswith("\\") else None
        if name in WORDS:
            self.group()
        elif name in TRANSPARENT:
            self.found += Reader(tokens=self.group()).read()
        elif name in MARKS:
            return "\\" + name + "{" + "".join(self.group()) + "}"
        elif name in GREEK or (name is None and token.isalpha()):
            return token
        return None

    def scripts(self):
        """The subscript and the superscripts after a symbol, each as its tokens. Primes and the
        empty group that spaces a tensor's indices, as in $R^\\mu{}_\\nu$, are passed over."""
        sub, sups = [], []
        while self.peek() in ("_", "^", "'") or self.tokens[self.at:self.at + 2] == ["{", "}"]:
            mark = self.take()
            if mark == "{":
                self.take()
            elif mark == "_":
                sub += self.group()
            elif mark == "^":
                sups.append(self.group())
        return sub, sups

    def read(self):
        """Every symbol of the expression as (letter, subscript tokens, whether it carries indices)."""
        while self.peek() is not None:
            if not (self.peek()[0] == "\\" or self.peek().isalpha()):
                self.take()
                continue
            previous = self.tokens[self.at - 1] if self.at else None
            letter = self.letter()
            if letter is None:
                continue
            following = self.peek() or ""
            if letter == SOLID_ANGLE and previous == "d":
                continue
            if letter == CHANGE and (following.isalpha() or following[1:] in GREEK):
                continue
            sub, sups = self.scripts()
            indexed = is_component(sub, INDICES)
            for sup in sups:
                if sup and all(token in INDICES for token in sup):
                    indexed = True
                else:
                    self.found += Reader(tokens=sup).read()
            if len(sups) == 1 and not sub and f"{letter}^{''.join(sups[0])}" in UNIVERSAL:
                continue
            self.found.append((letter, sub, indexed))
        return self.found


def is_component(sub, letters):
    """Whether a subscript holds only letters of `letters` beside its punctuation, and at least one."""
    return any(token in letters for token in sub) and \
        all(token in letters or not token.strip("\\").isalpha() for token in sub)


def symbols(tex, coords=()):
    """The symbols an expression uses. Each is a name such as `r_s`, or for a letter that carries
    a subscript of coordinates or indices a pair such as (`p_\\phi`, `p`): the whole and the letter,
    either of which may be what is defined. A tensor that carries indices is left out."""
    letters, found = set(coords) | INDICES, set()
    for letter, sub, indexed in Reader(tex).read():
        component = indexed or is_component(sub, letters)
        if letter in TENSORS and component:
            continue
        whole = letter + "_" + "".join(token for token in sub if token not in ("{", "}", "\\,")) if sub else letter
        found.add((whole, letter) if component else whole)
    return found


def names(found):
    r"""Every name a set of symbols answers to: a component $p_\phi$ answers to itself and to $p$,
    and $r_\pm$ to $r_+$ and $r_-$."""
    answer = set()
    for symbol in found:
        for name in (symbol if isinstance(symbol, tuple) else (symbol,)):
            answer.add(name)
            for sign in ("\\pm", "\\mp"):
                if name.endswith("_" + sign):
                    answer |= {name[:-len(sign)] + "+", name[:-len(sign)] + "-"}
    return answer


def alone(tex, coords=()):
    r"""The symbols of an expression that is nothing but symbols in a list, "p, q" or "f(r)", and
    the empty set for any other expression. A symbol's power, its differential and its derivative
    are still the symbol alone, "\rho^2", "dr_*" and "dz/dr", and so is a time written as a
    length, "c\tau"."""
    tex = ARGUMENTS.sub("", tex.strip())
    tex = LIGHT_SPEED.sub("", tex)
    tex = DIFFERENTIAL.sub(r"\1", tex)
    tex = POWER.sub("", tex)
    found = Reader(tex).read()
    if NOT_A_SYMBOL.sub("", tex) or not found or len(found) > tex.count(",") + 1:
        return set()
    return names(symbols(tex, coords))


def defined_in(text, coords=(), values=False):
    r"""The symbols a sentence defines, in one of three ways.

    By an equation with the symbol alone on its left, as $r_*$ and $p, q$ are in "$r_* = r +
    r_s\ln(r/r_s - 1)$" and "$p, q = \arctan(u, v)$", or as a function is in "$f(r) = 1 - r_s/r$".
    Setting a symbol to a number, "$a = 1$", gives its value and not its meaning, so it defines
    nothing unless `values` says the sentence is one that says what was drawn. A function holding
    the symbol, as in "$\sinh r_c = 1$", makes an equation for it, which does define it. An
    expression on the left, "$ex^2 = 2$", names nothing.

    By a noun of NAMING right before it, or before a list it stands in: "the height $z$", "for
    integer $k$", "the Kruskal coordinates $U$ and $V$".

    By saying what it is, alone or in a list: "with $\ell$ the unit of length", "$\bar\rho$ is the
    mean density", "$Y_0$ and $Y_1$ the Bessel functions"."""
    found = set()
    for tex in mathematics(text):
        sides = EQUALS.split(tex)
        for left, right in zip(sides, sides[1:]):
            left = RELATION.split(left)[-1]
            solved = any(command[1:] not in GREEK for command in re.findall(r"\\[A-Za-z]+", left))
            if values or solved or names(symbols(right, coords)) - set(UNIVERSAL):
                found |= alone(left, coords)
    for listed in NAMED.findall(text) + SAID.findall(text):
        for tex in mathematics(listed):
            found |= alone(tex, coords)
    return found


def as_sentence(domain):
    r"""A domain as a sentence. A domain is mathematics with its words inside \text{...}, the other
    way round from a sentence, which has its mathematics between dollar signs."""
    return re.sub(r"\$\s*\$", "", "$" + WORDED.sub(lambda words: f"${words.group(1)}$", domain) + "$")


def chart_definitions(metric, chart):
    """The symbols a chart defines: its coordinates, its parameters, its conventions and what its
    domains name."""
    coords = chart.get("coords") or []
    found = set()
    for coord in coords:
        found |= names(symbols(coord, coords))
    for parameter in chart.get("parameters") or []:
        found |= names(symbols(parameter.get("symbol") or "", coords))
        found |= defined_in(parameter.get("description") or "", coords)
    for text in (metric.get("convention") or "", chart.get("convention") or ""):
        for tex in mathematics(text):
            found |= names(symbols(tex, coords))
    for domain in chart.get("domains") or []:
        found |= defined_in(as_sentence(domain), coords)
    return found


def chart_mathematics(chart):
    """Yield (field, TeX) for every expression of a chart's own mathematics."""
    for domain in chart.get("domains") or []:
        yield "domains", domain
    if chart.get("line_element"):
        yield "line_element", chart["line_element"]
    for field in ("metric_components", "inverse_metric_components"):
        for component in chart.get(field) or []:
            yield field, component["value"]
    for field in ("christoffel", "riemann", "ricci_tensor", "einstein_tensor", "weyl_tensor"):
        for variant in ((chart.get(field) or {}).get("variants") or {}).values():
            for component in variant.get("nonzero") or []:
                yield field, component["value"]
    for field, scalar in SCALAR_FIELDS.items():
        if chart.get(field):
            yield field, re.sub(rf"^\s*{scalar}\s*=", "", chart[field])
    for equation in chart.get("geodesics") or []:
        yield "geodesics", equation


def drawings(diagram=None, conformal=None, embedding=None):
    """Yield (place, kind, view) for every drawing of a spacetime."""
    for part in ("systems", "projections"):
        for system, views in ((diagram or {}).get(part) or {}).items():
            for view in views:
                yield f"diagrams {system}/{view['id']}", "diagrams", view
    for kind, data in (("conformal", conformal), ("embedding", embedding)):
        for view in (data or {}).get("views") or []:
            yield f"{kind} {view['id']}", kind, view


def view_texts(view):
    """Yield (field, text) for every text of a drawing."""
    figure, movie = view.get("figure") or {}, view.get("movie") or {}
    for field in ("label", "xlabel", "ylabel", "unit", "settings", "input", "restriction", "cone", "height"):
        if isinstance(view.get(field), str):
            yield field, view[field]
    for field in ("caption", "families", "stops"):
        for text in view.get(field) or []:
            yield field, text
    for holder in (view, figure):
        for entry in holder.get("legend") or []:
            yield "legend", entry[2]
        for label in holder.get("labels") or []:
            yield "labels", label["text"]
    for marker in view.get("markers") or []:
        if marker.get("legend"):
            yield "legend", marker["legend"]
    for mark in view.get("slices") or []:
        yield "slices", mark["label"]
    for ticks in (view.get("ticks") or {}).values():
        for tick in ticks:
            yield "ticks", tick["label"]
    for shade in view.get("shades") or []:
        yield "shades", shade["legend"][2]
        yield "shades", shade["caption"]
    if movie.get("variable"):
        yield "movie", movie["variable"]
    for surface in (view.get("surfaces") or []) + (movie.get("frames") or []):
        texts = [surface.get("label")] + [curve.get("label") for curve in surface.get("curves") or []] + \
            [ring.get("label") for ring in surface.get("rings") or []] + \
            [piece.get(end, {}).get("text") for piece in surface.get("pieces") or []
             for end in ("start", "end", "edge")]
        for text in texts:
            if text:
                yield "surfaces", text


def undefined(tex, defined, coords):
    """The symbols of an expression that neither `defined` nor UNIVERSAL covers."""
    missing = set()
    for symbol in symbols(tex, coords):
        answers = symbol if isinstance(symbol, tuple) else (symbol,)
        if not any(name in defined or name in UNIVERSAL for name in answers):
            missing.add(answers[-1])
    return missing


def undefined_uses(metric, diagram=None, conformal=None, embedding=None):
    """Yield (place, who could define it, field, symbol, text) for every use of a symbol that
    nothing defines."""
    charts = {chart["id"]: chart for chart in metric.get("coordinates") or []}
    definitions = {chart_id: chart_definitions(metric, chart) for chart_id, chart in charts.items()}
    every_definition = set().union(*definitions.values())
    every_coord = [coord for chart in charts.values() for coord in chart.get("coords") or []]
    for chart_id, chart in charts.items():
        for field, tex in chart_mathematics(chart):
            for symbol in sorted(undefined(tex, definitions[chart_id], chart.get("coords") or [])):
                yield (f"the chart '{chart_id}'", "neither its coordinates, its parameters nor its convention",
                       field, symbol, tex)
    for place, kind, view in drawings(diagram, conformal, embedding):
        defined = every_definition | set(DRAWN.get(kind, ()))
        texts = list(view_texts(view))
        for field, text in texts:
            defined |= defined_in(text, every_coord, values=field in ("settings", "input"))
        for field, text in texts:
            for tex in mathematics(text):
                for symbol in sorted(undefined(tex, defined, every_coord)):
                    yield place, "neither a chart of the spacetime nor a definition of its own", field, symbol, text


def symbol_problems(metric, diagram=None, conformal=None, embedding=None):
    """Every symbol used where nothing defines it, one sentence for each symbol of each chart and
    of each drawing, naming the fields that use it."""
    found = {}
    for place, definer, field, symbol, _ in undefined_uses(metric, diagram, conformal, embedding):
        found.setdefault((place, definer, symbol), set()).add(field)
    return [f"{metric['id']}.json: {place} uses ${symbol}$ in {', '.join(sorted(fields))}, and {definer} defines it"
            for (place, definer, symbol), fields in found.items()]


def load(folder, name):
    path = DATA / folder / f"{name}.json"
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else None


def main(argv=None):
    given = sys.argv[1:] if argv is None else argv
    chosen = given or sorted(path.stem for path in (DATA / "metrics").glob("*.json"))
    problems = []
    for name in chosen:
        files = [load(folder, name) for folder in ("metrics", "diagrams", "conformal", "embedding")]
        problems += symbol_problems(*files)
    for problem in problems:
        print(problem)
    print(f"{len(problems)} symbols without a definition in {len(chosen)} spacetimes", file=sys.stderr)
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
