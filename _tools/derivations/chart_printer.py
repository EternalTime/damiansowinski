"""Print sympy expressions in the LaTeX dialect the metric files use, and build their blocks.

This is how the four entries added on 25 September 2026 (`tov`, `malament_hogarth`,
`mixmaster` and `lentz`) were written: `print_charts.py` beside this file defines each chart
and calls into here. Everything is computed in the chart the collection prints, whose time
coordinate is x^0 = ct already, so a derivative along it is the chart derivative and no factor
of c appears. Every printed value is read back through the checker's own Reader and compared
with the value it was printed from before it is written, so a printing slip cannot reach a
file; `verify_metrics.py` then checks the file as it checks every other entry.

The printer knows a small class of expressions: sums, products, rational powers, symbols,
undefined functions and their derivatives, exp, the logarithm, the trigonometric functions,
and the Dirac delta with its derivatives, as \\delta'(u). A value is printed as a sign, a
rational coefficient, and numerator and denominator factors. A sum prints
its terms in an order the chart chooses, (r - 2m) rather than (-2m + r), and a chart may pass a
`collect` function that regroups the numerator of a value, which is how TOV's curvature is
written around (d_r Phi)^2 + d_r^2 Phi and Bianchi IX's around a_1^2 cos^2 psi + a_2^2 sin^2 psi.
"""
import time

import sympy as sp

import verify_metrics as vm

GREEK = {"theta", "phi", "psi", "chi", "eta", "tau", "Phi", "Omega", "omega", "lambda", "mu", "nu", "rho", "ell", "alpha",
         "Lambda", "gamma", "sigma", "Delta", "kappa", "xi", "delta"}
# A name the reader spells from an accented command, as it reads \tilde\phi as tildephi.
ACCENTED = {"tildephi": "\\tilde\\phi"}
TRIG = (sp.sin, sp.cos, sp.tan, sp.cot, sp.csc, sp.sec, sp.sinh, sp.cosh)


def tex_name(name):
    if name in ACCENTED:
        return ACCENTED[name]
    # A Greek letter with a subscript, as rho_0, keeps its command.
    return "\\" + name if name in GREEK or name.partition("_")[0] in GREEK else name


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
    def __init__(self, coords, primed=(), lead=(), overrides=None, collect=None, factors=None, named=None,
                 rising=(), flip=True, last=(), dotted=()):
        """coords: coordinate symbols in chart order.
        primed: names of functions of one variable printed with primes.
        dotted: names of functions of the time printed with dots, as \\dot{a} and \\ddot{a},
            the way FRW writes its scale factor; each dot is a derivative along the chart's x^0.
        lead: generators, most significant first, that order the terms of a sum.
        overrides: {placeholder symbol: its printed text}.
        named: {placeholder symbol: the text of a sum it stands for}, bracketed wherever a sum
            would be, so a factor such as 1 + alpha r cos(theta) keeps the order it is written in.
        collect: a function turning the polynomial numerator of a value into a Sum.
        factors: generators in the order they are written within a product, lead by default.
        rising: generators, most significant first, that order the terms of a sum before
            lead does, lowest power first, so that 3r - 3r_s - Lambda r^3 keeps the order of
            1 - r_s/r - Lambda r^2/3.
        flip: whether -(A - B)/D is printed as (B - A)/D; a chart whose sums keep one order
            throughout passes False.
        last: expressions whose terms close a sum, as a chart with a kink puts the delta at
            the kink after the smooth part of every value.
        """
        self.coords = list(coords)
        self.primed = set(primed)
        self.dotted = set(dotted)
        self.lead = list(lead)
        self.factors = list(lead if factors is None else factors)
        self.rising = list(rising)
        self.flip = flip
        self.last = list(last)
        self.overrides = dict(overrides or {})
        self.named = dict(named or {})
        self.collect = collect

    # -- ordering ----------------------------------------------------------------------
    def degrees(self, term):
        _, rest = term.as_coeff_Mul()
        powers = rest.as_powers_dict()
        out = []
        for g in self.rising:
            d = sp.sympify(powers.get(g, 0))
            out.append(d if d.is_Number else 0)
        for g in self.lead:
            d = sp.sympify(powers.get(g, 0))
            out.append(-d if d.is_Number else 0)
        return out

    def ordered(self, terms):
        return sorted(terms, key=lambda t: (any(t[1].has(g) for g in self.last), self.degrees(t[1]),
                                            -sp.count_ops(t[1]), self.plain(t[1])))

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
        if isinstance(e, sp.DiracDelta):
            # A prime for each derivative, taken with respect to the delta's own argument.
            order = e.args[1] if len(e.args) > 1 else 0
            # Its argument is set tight, as ct - z is in the line element.
            return "\\delta" + "'" * order + "(" + self.positive_first(e.args[0]).replace("\\,", "") + ")"
        if isinstance(e, sp.Derivative):
            name = tex_name(e.expr.func.__name__)
            counts = {}
            for v, k in e.variable_count:
                counts[v] = counts.get(v, 0) + k
            if e.expr.func.__name__ in self.dotted:
                (n,) = counts.values()
                return ("\\dot{", "\\ddot{")[n - 1] + name + "}"
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

    def exponent(self, e):
        """An exponent set inline, as the files write $e^{r^2/R^2}$: a single term over a
        denominator is printed with a slash, since a built fraction in a superscript is unreadable."""
        num, den = sp.fraction(e)
        if den == 1 or sp.sympify(num).is_Add or sp.sympify(den).is_Add:
            return self.expr(e)
        head = "-" if sp.sympify(num).could_extract_minus_sign() else ""
        num = -num if head else num
        bottom = self.expr(den)
        if len(sp.Mul.make_args(den)) > 1:
            bottom = "(" + bottom + ")"
        return head + self.expr(num) + "/" + bottom

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
            if rest in self.named and c != 1:
                # A named sum under a coefficient or a minus sign keeps its brackets.
                out.append((c, ("" if abs(c) == 1 else _num(abs(c))) + "\\left(" + self.named[rest] + "\\right)"))
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
        if self.flip and top and sign < 0 and out[1]:
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
        return self.factors.index(g) if g in self.factors else 99

    def factor_key(self, item):
        base, _ = item
        if isinstance(base, Sum):
            return (3, len(base.terms))
        if base.is_Number:
            return (0,)
        if base in self.overrides:
            # An override the chart lists among its factors is written where a symbol would be.
            return (1, self.rank_of(base), "") if base in self.factors else (3, 100)
        if base in self.named:
            return (3, 50 + list(self.named).index(base))
        if isinstance(base, sp.Symbol):
            return (1, self.rank_of(base), base.name)
        if isinstance(base, sp.core.function.AppliedUndef):
            return (1, self.rank_of(base), base.func.__name__)
        if isinstance(base, sp.exp):
            return (4,)
        if isinstance(base, sp.log):
            return (6.5,)
        if isinstance(base, sp.DiracDelta):
            return (8,)
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
        if isinstance(base, Sum) or base in self.named:
            text = self.named[base] if base in self.named else self.sum_text(base)
            if exponent == 1:
                return text if alone else "\\left(" + text + "\\right)"
            if exponent == sp.Rational(1, 2):
                return "\\sqrt{" + text + "}"
            return "\\left(" + text + "\\right)^" + _sup(exponent)
        if exponent == sp.Rational(1, 2):
            return "\\sqrt{" + (_num(base) if base.is_Number else self.expr(base)) + "}"
        if isinstance(base, sp.exp):
            return "e^{" + self.exponent(base.args[0]) + "}"
        if isinstance(base, sp.log) and exponent == 1:
            # The argument is one fraction, as ln((x^2 + y^2)/rho_0^2), however sympy expanded it.
            return "\\ln\\left(" + self.expr(sp.factor(base.args[0])) + "\\right)"
        if isinstance(base, TRIG) and not exponent.is_Integer:
            # The reader takes \cos^2\tau but not \cos^{3/2}\tau, so a fractional power is bracketed.
            return "\\left(\\" + base.func.__name__ + self.trig_argument(base.args[0]) + "\\right)^" + _sup(exponent)
        if isinstance(base, TRIG):
            head = "\\" + base.func.__name__
            if exponent != 1:
                head += "^" + _sup(exponent)
            return head + self.trig_argument(base.args[0])
        text = _num(base) if base.is_Number else self.atom(base)
        if exponent == 1:
            return text
        if (isinstance(base, sp.Derivative) and base.expr.func.__name__ not in self.dotted) or "'" in text:
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
            # A named sum standing alone is bracketed once a coefficient or a sign stands before it.
            body = self.product(numerator, alone=(len(numerator) == 1 and p == 1 and not head))
            return head + (str(p) if p != 1 else "") + body
        top_text = self.product(numerator, alone=(len(numerator) == 1 and p == 1)) if numerator else ""
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


# -- a kink -----------------------------------------------------------------------------

def kink(x, then=sp.factor):
    """A `pretty` for a metric with a kink at x = 0, as the domain wall's (1 - k|x|)^2, and the
    placeholders it prints with, as (pretty, overrides).

    The checker hands every value back in x and sgn(x), its denominators cleared of the sign;
    here x is written sgn(x)|x|, the sign's powers are reduced by sgn(x)^2 = 1, inside every
    exponential too, and the value is split into A(|x|) + sgn(x) B(|x|), each part put in
    the form `then` gives it, so that 1/(1 - k|x|) is printed as it is written and not as
    (1 + k|x|)/(1 - k^2x^2). The delta of x is held aside while x is replaced."""
    a, s, d = sp.Symbol("_abs" + x.name, positive=True), sp.Symbol("_sgn" + x.name), sp.Symbol("_delta" + x.name)
    name = tex_name(x.name)
    overrides = {a: "|" + name + "|", s: "\\mathrm{sgn}(" + name + ")", d: "\\delta(" + name + ")"}

    def reduce(e):
        return e.replace(lambda f: f.is_Pow and f.base == s and f.exp.is_Integer, lambda f: s ** (int(f.exp) % 2))

    def pretty(value):
        v = sp.sympify(value).replace(lambda e: isinstance(e, sp.Abs) and e.args[0] == x, lambda e: a)
        v = v.xreplace({sp.DiracDelta(x): d, sp.sign(x): s}).subs(x, s * a)
        v = reduce(v.replace(lambda e: isinstance(e, sp.exp), lambda e: sp.exp(reduce(sp.expand(e.args[0])))))
        numerator, denominator = (reduce(sp.expand(side)) for side in sp.fraction(sp.together(v)))
        if denominator.has(s):
            raise ValueError(f"{value} keeps sgn({x}) below the line")
        even, odd = numerator.subs(s, 0), sp.expand((numerator - numerator.subs(s, 0)) / s)
        return then(even / denominator) + s * then(odd / denominator)
    return pretty, overrides


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
    def __init__(self, coords_tex, parameters, chart_line_element, printer_options=None, pretty=None, time=None,
                 bracketed=None):
        # Time is already the chart coordinate here, so no coordinate is scaled by c. A chart
        # whose components depend on the time names it as `time`: its symbol then stands for
        # x^0 = ct in the geometry, and every value is printed and read back with it written
        # as c times the time the file prints, as Milne's comoving c^2t^2 is.
        self.coords_tex = coords_tex
        self.reader = vm.Reader(coords_tex, parameters, ())
        self.symbols = [self.reader.symbol[name] for name in coords_tex]
        self.bare = {self.reader.symbol[time]: self.reader.c * self.reader.symbol[time]} if time else {}
        g = vm.metric_from_line_element(self.reader, chart_line_element, coords_tex)
        self.geo = vm.Geometry(g, self.symbols, 10 ** 6)
        self.printer = Printer(self.symbols, **(printer_options or {}))
        self.pretty = pretty or sp.factor
        # The sum a value is printed as when it has to stand in a bracket with a minus in front:
        # expanded, unless the chart writes its sums its own way, as a chart with a kink does.
        self.bracketed = bracketed or sp.expand
        self.own_brackets = bracketed is not None

    def check(self, text, value):
        if vm.norm(self.reader(text) - sp.sympify(value).subs(self.bare, simultaneous=True)) != 0:
            raise AssertionError(f"printed {text!r} does not read back as {value}")
        return text

    def text(self, value):
        """The printed form of a value, checked by reading it back."""
        return self.check(self.printer(self.pretty(sp.sympify(value).subs(self.bare, simultaneous=True))), value)

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
        if self.own_brackets:
            # A chart that writes its own sums keeps their order, and brackets whichever of the
            # value and its negation leads with a positive term, minus first.
            for sign, inner in ((-1, -value), (1, value)):
                s = self.printer.sum_of(self.bracketed(inner))
                if not self.printer.leads_negative(s.terms[0]):
                    text = "\\left(" + self.printer.sum_text(s) + "\\right)"
                    return self.check(text if sign == 1 else "-" + text, value)
        return self.check("-\\left(" + self.printer.positive_first(self.bracketed(-value)) + "\\right)", value)

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
