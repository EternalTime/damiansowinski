"""Print sympy expressions in the LaTeX dialect the metric files use, and build their blocks.

This is how the four entries added on 25 September 2026 (`tov`, `malament_hogarth`,
`mixmaster` and `lentz`) were written: `print_charts.py` beside this file defines each chart
and calls into here. Everything is computed in the chart the collection prints, whose time
coordinate is x^0 = ct already, so a derivative along it is the chart derivative and no factor
of c appears. Every printed value is read back through the checker's own Reader and compared
with the value it was printed from before it is written, so a printing slip cannot reach a
file; `verify_metrics.py` then checks the file as it checks every other entry.

The printer knows a small class of expressions: sums, products, rational powers, symbols,
undefined functions and their derivatives, exp and the trigonometric functions. A value is
printed as a sign, a rational coefficient, and numerator and denominator factors. A sum prints
its terms in an order the chart chooses, (r - 2m) rather than (-2m + r), and a chart may pass a
`collect` function that regroups the numerator of a value, which is how TOV's curvature is
written around (d_r Phi)^2 + d_r^2 Phi and Bianchi IX's around a_1^2 cos^2 psi + a_2^2 sin^2 psi.
"""
import time

import sympy as sp

import verify_metrics as vm

GREEK = {"theta", "phi", "psi", "chi", "eta", "tau", "Phi", "Omega", "omega", "lambda", "mu", "nu", "rho"}
TRIG = (sp.sin, sp.cos, sp.tan, sp.cot, sp.csc, sp.sec, sp.sinh, sp.cosh)


def tex_name(name):
    return "\\" + name if name in GREEK else name


class Sum:
    """A sum whose terms print in the order given, each held as (coefficient, rest) so that
    sympy never distributes a number over a bracket the entry wants kept."""

    def __init__(self, terms):
        self.terms = []
        for t in terms:
            if isinstance(t, tuple):
                self.terms.append(t)
            else:
                self.terms.append(sp.sympify(t).as_coeff_Mul())

    def negated(self):
        return Sum([(-c, rest) for c, rest in self.terms])

    def scaled(self, k):
        return Sum([(k * c, rest) for c, rest in self.terms])


class Printer:
    def __init__(self, coords, primed=(), lead=(), overrides=None, collect=None):
        """coords: coordinate symbols in chart order.
        primed: names of functions of one variable printed with primes.
        lead: generators, most significant first, that order the terms of a sum.
        overrides: {placeholder symbol: its printed text}.
        collect: a function turning the polynomial numerator of a value into a Sum.
        """
        self.coords = list(coords)
        self.primed = set(primed)
        self.lead = list(lead)
        self.overrides = dict(overrides or {})
        self.collect = collect

    # -- ordering ----------------------------------------------------------------------
    def degrees(self, term):
        _, rest = term.as_coeff_Mul()
        powers = rest.as_powers_dict()
        out = []
        for g in self.lead:
            d = sp.sympify(powers.get(g, 0))
            out.append(-d if d.is_Number else 0)
        return out

    def ordered(self, terms):
        return sorted(terms, key=lambda t: (self.degrees(t[1]), -sp.count_ops(t[1]), self.plain(t[1])))

    def plain(self, rest):
        try:
            return self.term(rest)
        except Exception:
            return str(rest)

    def sum_of(self, e):
        return Sum(self.ordered([t.as_coeff_Mul() for t in sp.Add.make_args(e)]))

    def positive_first(self, e):
        """A bare sum printed with its first positively signed term leading."""
        s = self.sum_of(e)
        for i, pair in enumerate(s.terms):
            if not self.leads_negative(pair):
                s = Sum([pair] + s.terms[:i] + s.terms[i + 1:])
                break
        return self.sum_text(s)

    def leads_negative(self, pair):
        """Whether a term of a Sum prints with a minus in front."""
        c, rest = pair
        if rest.is_Add:
            inner = self.sum_of(rest)
            return (c * inner.terms[0][0]).is_negative if abs(c) == 1 else c.is_negative
        body = self.term(abs(c) * rest)
        return c.is_negative != body.startswith("-")

    def collected(self, base):
        s = self.collect(base, self)
        return s if isinstance(s, Sum) else self.sum_of(s)

    # -- atoms -------------------------------------------------------------------------
    def atom(self, e):
        if e in self.overrides:
            return self.overrides[e]
        if isinstance(e, sp.Symbol):
            return tex_name(e.name)
        if isinstance(e, sp.core.function.AppliedUndef):
            return tex_name(e.func.__name__)
        if isinstance(e, sp.Derivative):
            name = tex_name(e.expr.func.__name__)
            counts = {}
            for v, k in e.variable_count:
                counts[v] = counts.get(v, 0) + k
            if e.expr.func.__name__ in self.primed:
                (n,) = counts.values()
                return name + "'" * n
            parts = []
            for c in self.coords:
                if c in counts:
                    parts.append(f"\\partial_{tex_name(c.name)}" + (f"^{counts[c]}" if counts[c] > 1 else ""))
            return "".join(parts) + ("" if name.startswith("\\") else " ") + name
        raise ValueError(f"not an atom: {e!r}")

    # -- the recursion -----------------------------------------------------------------
    def __call__(self, e):
        e = sp.sympify(e)
        if e.is_Add:
            s = self.collected(e) if self.collect else self.sum_of(e)
            if self.leads_negative(s.terms[0]) and len(s.terms) > 1:
                return "-\\left(" + self.sum_text(s.negated()) + "\\right)"
            return self.sum_text(s)
        return self.term(e, top=True)

    def expr(self, e):
        if isinstance(e, Sum):
            return self.sum_text(e)
        if e.is_Add:
            return self.sum_text(self.sum_of(e))
        return self.term(e)

    def sum_text(self, s):
        out = []
        for c, rest in s.terms:
            if rest.is_Add:
                inner = self.sum_of(rest)
                if abs(c) == 1:
                    for cc, rr in inner.scaled(c).terms:
                        body = self.term(abs(cc) * rr)
                        if body.startswith("-"):
                            cc, body = -cc, body[1:]
                        out.append((cc, body))
                    continue
                out.append((c, _num(abs(c)) + "\\left(" + self.sum_text(inner) + "\\right)"))
                continue
            body = self.term(abs(c) * rest)
            if body.startswith("-"):
                c, body = -c, body[1:]
            out.append((c, body))
        text = ""
        for i, (c, body) in enumerate(out):
            if c.is_negative:
                text += ("-" if i == 0 else " - ") + body
            else:
                text += ("" if i == 0 else " + ") + body
        return text

    def parts(self, e, top=False):
        """(sign, coefficient, numerator factors, denominator factors) of a product, with
        every sum a Sum whose leading term is positive and the trigonometry tidied."""
        coefficient, rest = e.as_coeff_Mul()
        sign = -1 if coefficient.is_negative else 1
        coefficient = abs(coefficient)
        numerator, denominator = [], []
        for f in sp.Mul.make_args(rest):
            if f == 1:
                continue
            base, exponent = (f.base, f.exp) if f.is_Pow else (f, sp.Integer(1))
            if isinstance(base, sp.exp):
                numerator.append((sp.exp(base.args[0] * exponent), sp.Integer(1)))
            elif exponent.is_negative:
                denominator.append((base, -exponent))
            else:
                numerator.append((base, exponent))
        numerator, denominator = self.tidy_trig(numerator, denominator)
        out = []
        for side in (numerator, denominator):
            done = []
            for base, exponent in side:
                if base.is_Add:
                    s = self.collected(base) if (self.collect and top and side is numerator) else self.sum_of(base)
                    if self.leads_negative(s.terms[0]):
                        s = s.negated()
                        if exponent.is_Integer and exponent % 2 == 1:
                            sign = -sign
                    base = s
                done.append((base, exponent))
            out.append(done)
        # -(A - B)/D reads better as (B - A)/D.
        if top and sign < 0 and out[1]:
            sums = [i for i, (b, k) in enumerate(out[0]) if isinstance(b, Sum) and k == 1]
            if len(sums) == 1:
                i = sums[0]
                s = out[0][i][0]
                if len(s.terms) > 1 and self.leads_negative(s.terms[-1]):
                    out[0][i] = (Sum(list(reversed(s.negated().terms))), sp.Integer(1))
                    sign = -sign
        return sign, coefficient, out[0], out[1]

    def tidy_trig(self, numerator, denominator):
        """cos^k/sin^k as cot^k, and what sine is left below the line as csc above it."""
        def split(side):
            powers, rest = {}, []
            for base, exponent in side:
                if isinstance(base, (sp.sin, sp.cos)) and exponent.is_Integer:
                    key = (base.func, base.args[0])
                    powers[key] = powers.get(key, 0) + exponent
                else:
                    rest.append((base, exponent))
            return powers, rest
        num, rest_num = split(numerator)
        den, rest_den = split(denominator)
        extra = []
        for (func, argument), k in list(den.items()):
            if func is sp.sin:
                c = num.get((sp.cos, argument), 0)
                pair = min(c, k)
                if pair:
                    extra.append((sp.cot(argument), sp.Integer(pair)))
                    num[(sp.cos, argument)] = c - pair
                    den[(func, argument)] = k - pair
                if den[(func, argument)]:
                    extra.append((sp.csc(argument), sp.Integer(den[(func, argument)])))
                    den[(func, argument)] = 0
        numerator = rest_num + [(f(a), sp.Integer(k)) for (f, a), k in num.items() if k] + extra
        denominator = rest_den + [(f(a), sp.Integer(k)) for (f, a), k in den.items() if k]
        return numerator, denominator

    def rank_of(self, g):
        return self.lead.index(g) if g in self.lead else 99

    def factor_key(self, item):
        base, _ = item
        if isinstance(base, Sum):
            return (3, len(base.terms))
        if base.is_Number:
            return (0,)
        if base in self.overrides:
            return (3, 100)
        if isinstance(base, sp.Symbol):
            return (1, self.rank_of(base), base.name)
        if isinstance(base, sp.core.function.AppliedUndef):
            return (1, self.rank_of(base), base.func.__name__)
        if isinstance(base, sp.exp):
            return (4,)
        if isinstance(base, sp.Derivative):
            return (5, self.rank_of(base), self.atom(base))
        if isinstance(base, TRIG):
            order = {sp.sin: 0, sp.cos: 1, sp.cot: 2, sp.csc: 3, sp.tan: 4, sp.sec: 5, sp.sinh: 6, sp.cosh: 7}
            return (6, str(base.args[0]), order[base.func])
        return (7, str(base))

    def product(self, factors, alone=None):
        factors = sorted(factors, key=self.factor_key)
        if alone is None:
            alone = len(factors) == 1
        out = ""
        previous = None
        for base, exponent in factors:
            piece = self.factor(base, exponent, alone)
            if out and _needs_space(out, piece, previous):
                out += "\\,"
            out += piece
            previous = base
        return out

    def factor(self, base, exponent, alone=False):
        if isinstance(base, Sum):
            text = self.sum_text(base)
            if exponent == 1:
                return text if alone else "\\left(" + text + "\\right)"
            if exponent == sp.Rational(1, 2):
                return "\\sqrt{" + text + "}"
            return "\\left(" + text + "\\right)^" + _sup(exponent)
        if exponent == sp.Rational(1, 2):
            return "\\sqrt{" + (_num(base) if base.is_Number else self.expr(base)) + "}"
        if isinstance(base, sp.exp):
            return "e^{" + self.expr(base.args[0]) + "}"
        if isinstance(base, TRIG):
            head = "\\" + base.func.__name__
            if exponent != 1:
                head += "^" + _sup(exponent)
            return head + self.trig_argument(base.args[0])
        text = _num(base) if base.is_Number else self.atom(base)
        if exponent == 1:
            return text
        if isinstance(base, sp.Derivative) or "'" in text:
            return "\\left(" + text + "\\right)^" + _sup(exponent)
        return text + "^" + _sup(exponent)

    def trig_argument(self, argument):
        if isinstance(argument, sp.Symbol):
            return tex_name(argument.name) if argument.name in GREEK else " " + argument.name
        return "\\left(" + self.expr(argument) + "\\right)"

    def term(self, e, top=False):
        if isinstance(e, Sum):
            return self.sum_text(e)
        if e.is_Add:
            return self.sum_text(self.sum_of(e))
        sign, coefficient, numerator, denominator = self.parts(e, top)
        head = "-" if sign < 0 else ""
        p, q = coefficient.p, coefficient.q
        # A coefficient against a lone sum is carried into it.
        if p != 1 and len(numerator) == 1 and isinstance(numerator[0][0], Sum) and numerator[0][1] == 1:
            numerator = [(numerator[0][0].scaled(p), sp.Integer(1))]
            p = 1
        if not denominator and q == 1:
            if not numerator:
                return head + str(p)
            if len(numerator) == 1 and isinstance(numerator[0][0], Sum) and numerator[0][1] == 1 and head:
                # A lone bracket with a minus in front: the sign goes into the sum as well.
                return self.sum_text(numerator[0][0].negated())
            body = self.product(numerator)
            return head + (str(p) if p != 1 else "") + body
        top_text = self.product(numerator) if numerator else ""
        if p != 1 or not top_text:
            top_text = str(p) + top_text
        bottom = self.product(denominator, alone=(q == 1 and len(denominator) == 1)) if denominator else ""
        if q != 1:
            bottom = str(q) + bottom
        return head + "\\dfrac{" + top_text + "}{" + bottom + "}"


def _num(n):
    n = sp.Rational(n)
    return str(n.p) if n.q == 1 else f"{n.p}/{n.q}"


def _sup(exponent):
    text = _num(exponent)
    return text if len(text) == 1 else "{" + text + "}"


TRIG_TEX = ("\\sin", "\\cos", "\\cot", "\\csc", "\\sinh", "\\cosh")


def _needs_space(left, right, previous):
    """Whether two adjacent factors want a thin space between them, as the files set them."""
    if right.startswith("e^"):
        return left[-1].isalpha() or left.endswith("'")
    if right.startswith("\\left") or right.startswith("\\dfrac"):
        return False
    if isinstance(previous, sp.exp) or left.endswith("\\right)"):
        return False
    if left.endswith("}") and not left.endswith("'}"):
        if not right.startswith("\\partial") and not right[0].isalpha():
            return False
    if left.lstrip("-").isdigit():
        return False
    if right.startswith(TRIG_TEX):
        return isinstance(previous, sp.Derivative) or left.endswith("'")
    return True


# -- hyperbolic functions --------------------------------------------------------------

def hyperbolic(r):
    """A `pretty` for a chart whose metric carries sinh r and cosh r, which the checker's
    Geometry hands back as exponentials: each value is rewritten in S = sinh r and C = cosh r,
    reduced by C^2 = 1 + S^2, its denominator cleared of C by its conjugate, factored, and every
    factor 1 + S^2 written as C^2, as Godel's cylindrical chart is printed."""
    S, C = sp.symbols("_S _C", positive=True)
    relation = [C ** 2 - 1 - S ** 2]

    def reduce(p):
        return sp.expand(sp.reduced(sp.expand(p), relation, C, S)[1])

    def pretty(value):
        x = sp.sympify(value)
        x = x.replace(lambda e: isinstance(e, sp.exp) and sp.expand(e.args[0] / r).is_Integer,
                      lambda e: (S + C) ** sp.expand(e.args[0] / r))
        num, den = (reduce(p) for p in sp.fraction(sp.together(x)))
        d0, d1 = sp.Poly(den, C).coeff_monomial(1), sp.Poly(den, C).coeff_monomial(C)
        if d1 != 0:
            num, den = reduce(num * (d0 - d1 * C)), sp.expand(d0 ** 2 - d1 ** 2 * (1 + S ** 2))
        out = sp.factor(sp.cancel(sp.factor(num) / sp.factor(den)))
        out = out.replace(lambda e: sp.expand(e - (S ** 2 + 1)) == 0, lambda e: C ** 2)
        # (S - 1)(S + 1) is read as the one factor S^2 - 1 it came from.
        powers = sp.Mul.make_args(out)
        pair = {}
        for f in powers:
            base, k = (f.base, f.exp) if f.is_Pow else (f, 1)
            if sp.expand(base - (S - 1)) == 0 or sp.expand(base - (S + 1)) == 0:
                pair[sp.expand(base)] = k
        if len(pair) == 2 and len(set(pair.values())) == 1:
            k = next(iter(pair.values()))
            out = sp.Mul(*[f for f in powers if sp.expand((f.base if f.is_Pow else f)) not in pair]) * (S ** 2 - 1) ** k
        return out.subs({S: sp.sinh(r), C: sp.cosh(r)})
    return pretty


# -- regrouping a numerator ------------------------------------------------------------

def symbolize(expr):
    """The expression with every derivative and function application replaced by a plain
    symbol, derivatives first, and the map back, so Poly can take them as generators."""
    forward, back = {}, {}
    derivatives = sorted(expr.atoms(sp.Derivative), key=lambda d: -sum(n for _, n in d.variable_count))
    for i, d in enumerate(derivatives):
        s = sp.Symbol(f"_D{i}")
        forward[d] = s
        back[s] = d
    expr = expr.xreplace(forward)
    functions = expr.atoms(sp.core.function.AppliedUndef)
    for i, f in enumerate(functions):
        s = sp.Symbol(f"_F{i}")
        forward[f] = s
        back[s] = f
    return expr.xreplace({f: forward[f] for f in functions}), forward, back


def collect_by(poly, generators, printer, substitutions=None):
    """poly as a Sum grouped by monomials in the generators, each coefficient factored."""
    poly = sp.expand(poly.subs(substitutions) if substitutions else poly)
    e, forward, back = symbolize(poly)
    gens = [forward.get(g, g) for g in generators]
    P = sp.Poly(e, *gens)
    if P.total_degree() == 0:
        return printer.sum_of(poly)
    terms = []
    for monomial, coeff in sorted(P.terms(), key=lambda mc: tuple(-k for k in mc[0])):
        c, rest = sp.factor(coeff.as_expr().xreplace(back)).as_coeff_Mul()
        terms.append((c, rest * sp.Mul(*[g ** k for g, k in zip(generators, monomial)])))
    return Sum(terms)


def merge_squares(e):
    """(u - v)(u + v) written back as u^2 - v^2, for every such pair of factors."""
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
                terms = sp.Add.make_args(product)
                if len(terms) == 2 and all(t.as_coeff_Mul()[1].is_Pow and t.as_coeff_Mul()[1].exp == 2 for t in terms):
                    items[i] = [product, k1]
                    del items[j]
                    changed = True
                    break
            if changed:
                break
    return sp.Mul(*[b ** k for b, k in items])


def homogenize(e, pairs):
    """Restore sin^2 + cos^2 = 1 where a canonical form took it away: within each parity of
    total degree in one (sin, cos) pair, raise every term to the highest degree present."""
    e = sp.expand(e)
    for s, c in pairs:
        degree = {}
        for term in sp.Add.make_args(e):
            degree[term] = sp.degree(term, s) + sp.degree(term, c)
        top = {}
        for d in degree.values():
            top[d % 2] = max(top.get(d % 2, 0), d)
        e = sp.expand(sp.Add(*[term * (s ** 2 + c ** 2) ** ((top[d % 2] - d) // 2) for term, d in degree.items()]))
    return e


_inner_count = [0]


def _trig_part(e):
    return sp.Mul(*[f for f in sp.Mul.make_args(e) if (f.base if f.is_Pow else f).func in (sp.sin, sp.cos)])


def collect_nested(poly, outer_pair, inner_pair, printer):
    """poly grouped first by monomials in the outer (sin, cos) pair, homogenized there, and
    each coefficient homogenized in the inner pair and grouped by its monomials in turn.

    For Bianchi IX this keeps a_1^2 cos^2 psi + a_2^2 sin^2 psi together as the unit it is,
    the metric of the rotated 1-2 frame, inside a sum over the powers of sin theta."""
    e, forward, back = symbolize(sp.expand(poly))
    so, co = [forward.get(g, g) for g in outer_pair]
    si, ci = [forward.get(g, g) for g in inner_pair]
    e = homogenize(e, [(so, co)])
    P = sp.Poly(e, so, co)
    groups = []
    for (ks, kc), C in sorted(P.terms(), key=lambda mc: (-mc[0][0], -mc[0][1])):
        Q = sp.Poly(homogenize(C.as_expr(), [(si, ci)]), si, ci)
        inner = []
        for (js, jc), D in Q.terms():
            c, rest = merge_squares(sp.factor(D.as_expr().xreplace(back))).as_coeff_Mul()
            inner.append((c, rest * inner_pair[0] ** js * inner_pair[1] ** jc, rest))
        inner.sort(key=lambda t: printer.term(t[2]))
        groups.append((outer_pair[0] ** ks * outer_pair[1] ** kc, inner))
    if len(groups) == 1:
        mono, inner = groups[0]
        return Sum([(c, rest * mono) for c, rest, _ in inner])
    terms = []
    for mono, inner in groups:
        if len(inner) == 1:
            c, rest, _ = inner[0]
            terms.append((c, rest * mono))
            continue
        inner_expr = sp.Add(*[c * rest for c, rest, _ in inner])
        content, primitive = sp.factor_terms(inner_expr).as_coeff_Mul()
        common, body = sp.Integer(1), primitive
        if primitive.is_Mul:
            adds = [f for f in sp.Mul.make_args(primitive) if f.is_Add]
            if len(adds) == 1:
                body = adds[0]
                common = sp.Mul(*[f for f in sp.Mul.make_args(primitive) if not f.is_Add])
        if not body.is_Add:
            terms.append((content, primitive * mono))
            continue
        s = printer.sum_of(body)
        s = Sum(sorted(s.terms, key=lambda t: printer.term(t[1] / _trig_part(t[1]))))
        if printer.leads_negative(s.terms[0]):
            positive = [i for i, pair in enumerate(s.terms) if not printer.leads_negative(pair)]
            if positive:
                i = positive[0]
                s = Sum([s.terms[i]] + s.terms[:i] + s.terms[i + 1:])
            else:
                s = s.negated()
                content = -content
        _inner_count[0] += 1
        placeholder = sp.Symbol(f"_S{_inner_count[0]}")
        printer.overrides[placeholder] = "\\left(" + printer.sum_text(s) + "\\right)"
        terms.append((content, common * placeholder * mono))
    return Sum(terms)


# -- a chart and its blocks ------------------------------------------------------------

def negate(text):
    return text[1:] if text.startswith("-") else "-" + text


def is_single_term(text):
    """True when no + or - sits at the top level after the first character."""
    depth = 0
    i = 0
    while i < len(text):
        if text.startswith("\\left", i):
            depth += 1
            i += 5
            continue
        if text.startswith("\\right", i):
            depth -= 1
            i += 6
            continue
        ch = text[i]
        if ch == "{":
            depth += 1
        elif ch == "}":
            depth -= 1
        elif ch in "+-" and depth == 0 and i > 0:
            return False
        i += 1
    return True


class Chart:
    def __init__(self, coords_tex, parameters, chart_line_element, printer_options=None, pretty=None):
        # Time is already the chart coordinate here, so no coordinate is scaled by c.
        self.coords_tex = coords_tex
        self.reader = vm.Reader(coords_tex, parameters, ())
        self.symbols = [self.reader.symbol[name] for name in coords_tex]
        g = vm.metric_from_line_element(self.reader, chart_line_element, coords_tex)
        self.geo = vm.Geometry(g, self.symbols, 10 ** 6)
        self.printer = Printer(self.symbols, **(printer_options or {}))
        self.pretty = pretty or sp.factor

    def check(self, text, value):
        if vm.norm(self.reader(text) - value) != 0:
            raise AssertionError(f"printed {text!r} does not read back as {value}")
        return text

    def text(self, value):
        """The printed form of a value, checked by reading it back."""
        return self.check(self.printer(self.pretty(value)), value)

    def single_term(self, value):
        """A printed form of the value with no top level sum, so a leading minus negates it.
        A sum that has to be bracketed is printed as minus the bracketed negation, so what the
        page shows reads as a negated quantity rather than as a stray bracket."""
        text = self.text(value)
        if is_single_term(text) and not text.startswith("\\left("):
            return text
        candidate = negate(self.text(-value))
        if is_single_term(candidate) and not candidate.startswith("\\left("):
            return candidate
        return self.check("-\\left(" + self.printer.positive_first(sp.expand(-value)) + "\\right)", value)

    def block(self, tensor, rank):
        """Nonzero components, each class of values listed together, and a value and its
        negation printed so that the second is the first with a minus in front, which is
        what lets the page merge them onto one line."""
        n = len(self.symbols)
        classes = []
        for index in vm._indices(n, rank):
            value = vm._at(tensor, index)
            if value == 0:
                continue
            for entry in classes:
                if vm.norm(value - entry[0]) == 0:
                    entry[1].append((index, 1))
                    break
                if vm.norm(value + entry[0]) == 0:
                    entry[1].append((index, -1))
                    break
            else:
                classes.append([value, [(index, 1)]])
        out = []
        for value, members in classes:
            text = self.single_term(value) if len({s for _, s in members}) > 1 else self.text(value)
            for index, sign in members:
                out.append({"indices": [self.coords_tex[i] for i in index],
                            "value": text if sign == 1 else negate(text)})
        for entry in out:
            index = [self.coords_tex.index(name) for name in entry["indices"]]
            self.check(entry["value"], vm._at(tensor, index))
        return out

    def mathematics(self, log=print):
        geo = self.geo
        n = len(self.symbols)

        def t(name, build):
            start = time.time()
            result = build()
            log(f"  {name}: {time.time() - start:.1f}s")
            return result

        g = [[geo.g[i, j] for j in range(n)] for i in range(n)]
        ginv = [[geo.ginv[i, j] for j in range(n)] for i in range(n)]
        riemann = geo.riemann_llll()
        ricci = geo.ricci_ll()
        einstein = geo.einstein_ll()
        weyl = geo.weyl_llll()

        def variants(field, specs):
            return {"default": specs[0][0], "variants": {
                v: {"note": note, "nonzero": t(f"{field} {v}", build)} for v, note, build in specs}}

        return {
            "metric_components": t("metric", lambda: self.block(g, 2)),
            "inverse_metric_components": t("inverse", lambda: self.block(ginv, 2)),
            "christoffel": variants("christoffel", [
                ("ull", "\\Gamma^\\mu{}_{\\nu\\rho}", lambda: self.block(geo.christoffel_ull(), 3)),
                ("lll", "\\Gamma_{\\mu\\nu\\rho}", lambda: self.block(geo.christoffel_lll(), 3))]),
            "riemann": variants("riemann", [
                ("ulll", "R^\\mu{}_{\\nu\\rho\\sigma}", lambda: self.block(geo.raise_indices(riemann, 4, (0,)), 4)),
                ("llll", "R_{\\mu\\nu\\rho\\sigma}", lambda: self.block(riemann, 4))]),
            "ricci_tensor": variants("ricci", [
                ("ll", "R_{\\mu\\nu}", lambda: self.block(ricci, 2)),
                ("ul", "R^\\mu{}_\\nu", lambda: self.block(geo.raise_indices(ricci, 2, (0,)), 2)),
                ("uu", "R^{\\mu\\nu}", lambda: self.block(geo.raise_indices(ricci, 2, (0, 1)), 2))]),
            "ricci_scalar": "R = " + t("ricci scalar", lambda: self.text(geo.ricci_scalar())),
            "einstein_tensor": variants("einstein", [
                ("ll", "G_{\\mu\\nu}", lambda: self.block(einstein, 2)),
                ("ul", "G^\\mu{}_\\nu", lambda: self.block(geo.raise_indices(einstein, 2, (0,)), 2)),
                ("uu", "G^{\\mu\\nu}", lambda: self.block(geo.raise_indices(einstein, 2, (0, 1)), 2))]),
            "weyl_tensor": variants("weyl", [
                ("ulll", "C^\\mu{}_{\\nu\\rho\\sigma}", lambda: self.block(geo.raise_indices(weyl, 4, (0,)), 4)),
                ("llll", "C_{\\mu\\nu\\rho\\sigma}", lambda: self.block(weyl, 4))]),
            "geodesics": t("geodesics", self.geodesics),
        }

    def geodesics(self):
        """x''^mu + Gamma^mu_{nu rho} x'^nu x'^rho = 0, one equation per coordinate."""
        gamma = self.geo.christoffel_ull()
        n = len(self.symbols)
        lines = []
        for mu in range(n):
            text = "\\ddot{" + self.coords_tex[mu] + "}"
            for nu in range(n):
                for rho in range(nu, n):
                    value = gamma[mu][nu][rho] * (1 if nu == rho else 2)
                    if value == 0:
                        continue
                    velocity = ("\\dot{" + self.coords_tex[nu] + "}^2" if nu == rho else
                                "\\dot{" + self.coords_tex[nu] + "}\\dot{" + self.coords_tex[rho] + "}")
                    coefficient = self.text(value)
                    negative = coefficient.startswith("-") and is_single_term(coefficient)
                    if negative:
                        coefficient = coefficient[1:]
                    if not is_single_term(coefficient):
                        coefficient = "\\left(" + coefficient + "\\right)"
                    if coefficient == "1":
                        coefficient = ""
                    elif coefficient[-1].isalpha() or coefficient.endswith("'") or "\\partial" in coefficient[-14:]:
                        coefficient += "\\,"
                    text += (" - " if negative else " + ") + coefficient + velocity
            lines.append(text + " = 0")
        return lines
