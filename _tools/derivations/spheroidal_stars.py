"""Complete neutral stellar families, with source locations in spheroidal_stars.md."""
import functools

import sympy as sp

METRIC = "spheroidal_stars"
SYSTEMS = ("vaidya_tikekar", "durgapal_iii", "durgapal_iv", "durgapal_v", "finch_skea")
NAMES = ("Vaidya-Tikekar and Durgapal-Bannerji", "Durgapal III (Heintzmann IIa)",
         "Durgapal IV", "Durgapal V", "Finch and Skea")
COORDS = ["t", "r", "\\theta", "\\phi"]
SURFACES = (sp.Rational(4, 11), sp.Rational(1, 8), sp.Rational(1, 11),
            sp.Rational(1, 14), sp.Rational(2, 3))
LINE = (r"ds^2 = -y^2{c2}dt^2 + \dfrac{dr^2}{Z}"
        r" + r^2\left(d\theta^2 + \sin^2\theta\,d\phi^2\right)")


def parameters(system):
    common = ["C", "r_b"]
    x = r"x = Cr^2"
    if system == "vaidya_tikekar":
        return common + ["A", "B", x, r"Z = \dfrac{2-x}{2(1+x)}",
                         r"y = A(1+x)^{3/2}+B(5+2x)\sqrt{2-x}"]
    if system == "finch_skea":
        return common + ["A", "B", x, r"s = \sqrt{1+x}", r"Z = \dfrac{1}{1+x}",
                         r"y = A(\cos s+s\sin s)+B(\sin s-s\cos s)"]
    n = {"durgapal_iii": 3, "durgapal_iv": 4, "durgapal_v": 5}[system]
    radial = {
        3: r"Z = 1-\dfrac{3x}{2(1+x)}\left(1+\dfrac{k}{\sqrt{1+4x}}\right)",
        4: r"Z = \dfrac{7-10x-x^2}{7(1+x)^2}+\dfrac{kx}{(1+x)^2(1+5x)^{2/5}}",
        5: r"Z = \dfrac{1-x(309+54x+8x^2)/112+kx/(1+6x)^{1/3}}{(1+x)^3}",
    }[n]
    return common + ["a", "k", x, radial, f"y = a(1+x)^{{{n}/2}}"]


def register_dimensions(vm):
    for system in SYSTEMS:
        vm.HELD[(METRIC, system)] = ("y", "Z")
        vm.DIMENSIONS[(METRIC, system)] = {
            "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "C": "L^-2", "r_b": "L",
            **{value.split(" = ")[0]: "1" for value in parameters(system)[2:]},
        }


@functools.lru_cache(None)
def potentials(system):
    """Explicit source potentials in x, independently of the tensor printer."""
    x = sp.Symbol("x", real=True)
    A, B, a, k = sp.symbols("A B a k", real=True)
    if system == "vaidya_tikekar":
        return x, (2-x)/(2*(1+x)), A*(1+x)**sp.Rational(3, 2)+B*(5+2*x)*sp.sqrt(2-x)
    if system == "finch_skea":
        s = sp.sqrt(1+x)
        return x, 1/(1+x), A*(sp.cos(s)+s*sp.sin(s))+B*(sp.sin(s)-s*sp.cos(s))
    n = SYSTEMS.index(system)+2
    Z = {
        3: 1-3*x/(2*(1+x))*(1+k/sp.sqrt(1+4*x)),
        4: (7-10*x-x*x)/(7*(1+x)**2)+k*x/((1+x)**2*(1+5*x)**sp.Rational(2, 5)),
        5: (1-x*(309+54*x+8*x*x)/112+k*x/(1+6*x)**sp.Rational(1, 3))/(1+x)**3,
    }[n]
    return x, Z, a*(1+x)**sp.Rational(n, 2)


@functools.lru_cache(None)
def drawing_parameters(system):
    """r_b=1 and u=1/5, with all constants fixed by p(q)=0 and y(q)^2=Z(q)."""
    x, Z, y = potentials(system)
    q = SURFACES[SYSTEMS.index(system)]
    target = sp.Rational(3, 5)
    by_name = {str(v): v for v in (Z.free_symbols | y.free_symbols)-{x}}
    if "k" in by_name:
        k, a = by_name["k"], by_name["a"]
        values = {k: (target-Z.subs(k, 0).subs(x, q))/sp.diff(Z, k).subs(x, q),
                  a: sp.sqrt(target)/(1+q)**sp.Rational(SYSTEMS.index(system)+2, 2)}
    else:
        A, B = by_name["A"], by_name["B"]
        f, g = sp.diff(y, A), sp.diff(y, B)
        f0, g0, f1, g1 = [e.subs(x, q) for e in (f, g, sp.diff(f, x), sp.diff(g, x))]
        slope = sp.sqrt(target)*(1-target)/(4*q*target)
        det = f0*g1-f1*g0
        values = {A: (sp.sqrt(target)*g1-slope*g0)/det,
                  B: (slope*f0-sp.sqrt(target)*f1)/det}
    return {"C": str(sp.N(q, 45)), "r_b": "1", **{str(k): str(sp.N(v, 45)) for k, v in values.items()}}


def spec(system):
    import verify_metrics as vm
    probe = vm.Reader(COORDS, parameters(system), (), held=("y", "Z"))
    r, y, Z = probe.symbol["r"], probe.parameters["y"], probe.parameters["Z"]
    return {
        "metric_id": METRIC,
        "system": {"id": system, "name": NAMES[SYSTEMS.index(system)], "coords": COORDS,
                   "domains": [r"t \in (-\infty,\infty)", r"r \in [0,r_b]",
                               r"\theta \in [0,\pi]", r"\phi \in [0,2\pi)",
                               r"C>0,\quad y>0,\quad Z>0"],
                   "parameters": parameters(system), "line_element": LINE.replace("{c2}", "c^2")},
        "chart_line_element": LINE.replace("{c2}", ""),
        "printer": {"lead": [r, y, Z]},
        "check": lambda chart: check_chart(chart, system),
    }


def check_chart(chart, system):
    """Point comparisons of both potentials and the full Einstein tensor with sourced matter.

    Centre, first pressure zero, junction, monotonicity, and sound speed are checked here once.
    No symbolic determinant or simplify is used.
    """
    import mpmath as mp
    import verify_metrics as vm
    x, z, y = potentials(system)
    C = chart.reader.parameters["C"]
    r = chart.reader.symbol["r"]
    explicit_z = z.subs(x, C*r*r)
    explicit_y = y.subs(x, C*r*r)
    # Potentials constructed independently from the source, compared to the public definitions.
    source = {chart.reader.parameters["Z"]: explicit_z, chart.reader.parameters["y"]: explicit_y}
    # Rebind the source's constants to the reader's symbols.
    constants = {v: chart.reader.parameters[str(v)] for e in source.values() for v in e.free_symbols
                 if str(v) in chart.reader.parameters}
    source = {k: v.subs(constants) for k, v in source.items()}
    pars = drawing_parameters(system)
    at = {chart.reader.parameters[k]: sp.Float(v, 45) for k, v in pars.items()}
    at.update({chart.reader.c: 1, chart.reader.symbol["t"]: 0,
               chart.reader.symbol["\\theta"]: sp.Rational(7, 5), chart.reader.symbol["\\phi"]: 0})
    def out(e):
        return sp.sympify(e).subs(source).doit().subs(at)
    for held, value in chart.reader.held.items():
        residual = out(held)-value.subs(at)
        for radius in (sp.Rational(1, 9), sp.Rational(3, 7), sp.Rational(8, 9)):
            if abs(sp.N(residual.subs(r, radius), 35)) > sp.Float("1e-30"):
                raise AssertionError(f"{system}: a potential differs from its source")
    G = chart.geo.raise_indices(chart.geo.einstein_ll(), 2, (0,))
    zr, yr = source[chart.reader.parameters["Z"]], source[chart.reader.parameters["y"]]
    density = (1-zr)/r**2-sp.diff(zr, r)/r
    pressure = 2*zr*sp.diff(yr, r)/(r*yr)+(zr-1)/r**2
    for radius in (sp.Rational(1, 9), sp.Rational(3, 7), sp.Rational(8, 9)):
        for i in range(4):
            expected = -density if i == 0 else pressure
            residual = out(vm._at(G, (i, i))-expected).subs(r, radius)
            if abs(sp.N(residual, 35)) > sp.Float("1e-28"):
                raise AssertionError(f"{system}: sourced density or isotropic pressure differs at {radius}")
    # The matter and junction audit uses x directly, including exact Taylor limits at the centre.
    mp.mp.dps = 40
    values = {v: sp.Float(pars[str(v)], 45) for v in z.free_symbols | y.free_symbols if str(v) in pars}
    z, y = z.subs(values), y.subs(values)
    rho = (1-z)/x-2*sp.diff(z, x)
    p = 4*z*sp.diff(y, x)/y+(z-1)/x
    # At x=0 use Taylor coefficients, never a divided expression.
    z1, z2, y0, y1, y2 = [mp.mpf(str(sp.N(e.subs(x, 0), 38))) for e in
                           (sp.diff(z, x), sp.diff(z, x, 2), y, sp.diff(y, x), sp.diff(y, x, 2))]
    assert y0 > 0 and -3*z1 > 0 and 4*y1/y0+z1 > 0
    speed0 = (4*(z1*y1/y0+y2/y0-(y1/y0)**2)+z2/2)/(-5*z2/2)
    assert 0 < speed0 < 1
    funcs = sp.lambdify(x, [z, y, rho, p, sp.diff(rho, x), sp.diff(p, x)], "mpmath")
    q = mp.mpf(str(sp.N(SURFACES[SYSTEMS.index(system)], 38)))
    for i in range(1, 101):
        Z, Y, R, P, dR, dP = funcs(q*i/100)
        assert Z > 0 and Y > 0 and R > 0 and dR < 0 and 0 < dP/dR < 1
        assert P > 0 if i < 100 else abs(P) < mp.mpf("1e-30")
    Z, Y, *_ = funcs(q)
    assert abs(Z-mp.mpf("0.6")) < mp.mpf("1e-30")
    assert abs(Y*Y-Z) < mp.mpf("1e-30")


def register(pc):
    pc.CHARTS[METRIC] = [lambda s=s: spec(s) for s in SYSTEMS]
