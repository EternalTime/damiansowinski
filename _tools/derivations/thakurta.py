"""The spherical Thakurta family, with the scale factor tied to static time.

Chart sources and the restrictions on the expanding exterior are in thakurta.md.
The checks compare numeric pullbacks and the published fluid source at points.
"""
import sympy as sp

import verify_metrics as vm


CHART_IDS = ("conformal", "cosmic", "areal", "isotropic", "tortoise", "kruskal")
ANGLES = [r"\theta \in [0, \pi]", r"\phi \in [0, 2\pi)"]
SPHERE = r"\left(d\theta^2 + \sin^2\theta\,d\phi^2\right)"
RADIUS_X = r"r = 2m\left(1 + \mathrm{W}\left(e^{x/(2m)-1}\right)\right)"
RADIUS_UV = r"r = 2m\left(1 + \mathrm{W}\left(-UV e^{-1}\right)\right)"
ETA_UV = r"\eta = 2m\ln\left(-V/U\right)"


def chart_spec(system):
    """Mello (2)-(4), Sato (2.77), and exterior tortoise and Kruskal maps."""
    parameters = ["m"]
    if system in ("conformal", "cosmic"):
        time = r"\eta" if system == "conformal" else "t"
        coords = [time, "r", r"\theta", r"\phi"]
        parameters += [f"a = a({time})"]
        f = r"\left(1-\dfrac{2m}{r}\right)"
        radial = rf"\dfrac{{dr^2}}{{{f}}} + r^2{SPHERE}"
        line = (rf"ds^2 = a^2\left(-{f}d\eta^2 + {radial}\right)" if system == "conformal"
                else rf"ds^2 = -{f}dt^2 + a^2\left({radial}\right)")
        domains = [rf"{time} \in (0, \infty)", r"r \in (2m, \infty)"]
    elif system == "areal":
        coords = ["t", "R", r"\theta", r"\phi"]
        parameters += ["a = a(t)"]
        f = r"\left(1-\dfrac{2ma}{R}\right)"
        line = (rf"ds^2 = -{f}dt^2 + \dfrac{{\left(dR-\dfrac{{\partial_t a}}{{a}}R\,dt\right)^2}}{{{f}}}"
                rf" + R^2{SPHERE}")
        domains = [r"t \in (0, \infty)", r"R \in (2ma, \infty)"]
    elif system == "isotropic":
        coords = [r"\eta", r"\sigma", r"\theta", r"\phi"]
        parameters += [r"a = a(\eta)"]
        line = (r"ds^2 = a^2\left(-\dfrac{\left(1-m/(2\sigma)\right)^2}{\left(1+m/(2\sigma)\right)^2}d\eta^2"
                rf" + \left(1+\dfrac{{m}}{{2\sigma}}\right)^4\left(d\sigma^2+\sigma^2{SPHERE}\right)\right)")
        domains = [r"\eta \in (0, \infty)", r"\sigma \in (m/2, \infty)"]
    elif system == "tortoise":
        coords = [r"\eta", "x", r"\theta", r"\phi"]
        parameters += [r"a = a(\eta)", RADIUS_X]
        f = r"\left(1-\dfrac{2m}{r}\right)"
        line = rf"ds^2 = a^2\left({f}\left(-d\eta^2+dx^2\right)+r^2{SPHERE}\right)"
        domains = [r"\eta \in (0, \infty)", r"x \in (-\infty, \infty)"]
    elif system == "kruskal":
        coords = ["U", "V", r"\theta", r"\phi"]
        parameters += [r"\eta_0", RADIUS_UV, ETA_UV]
        # UV = -(r/(2m)-1) exp(r/(2m)); this rational form avoids cancellations.
        line = (r"ds^2 = \dfrac{\eta^4}{\eta_0^4}\left(\dfrac{16m^2(r-2m)}{rUV}dU\,dV"
                rf" + r^2{SPHERE}\right)")
        domains = [r"U \in (-\infty, 0)", r"V \in (-U, \infty)"]
    else:
        raise ValueError(system)
    held = vm.HELD.get(("thakurta", system), ())
    rates = vm.RATES.get(("thakurta", system), {})
    probe = vm.Reader(coords, parameters, (), held=held, rates=rates)
    reduce = probe.surface if rates else None
    return {"metric_id": "thakurta", "system": {"id": system,
            "name": {"conformal": "Conformal Spherical", "cosmic": "Cosmic Spherical", "areal": "Cosmic Areal",
                     "isotropic": "Conformal Isotropic", "tortoise": "Tortoise", "kruskal": "Kruskal Null, Dust Expansion"}[system],
            "coords": coords, "parameters": parameters, "domains": domains + ANGLES, "line_element": line},
            "chart_line_element": line, "printer": {}, "reduce": reduce,
            "check": lambda chart: check_chart(chart, system)}


def check_chart(chart, system):
    """Numeric pullbacks from the conformal metric, and its isotropic fluid source."""
    X, P = chart.symbols, chart.reader.parameters
    m = P["m"]
    # A positive, nonconstant test scale, with nonzero first and second derivatives.
    z = sp.Symbol("z", positive=True)
    A = 1 + z / 5 + z ** 2 / 20
    for eta0, r0, mass in ((2.3, 5.2, 0.7), (4.1, 8.5, 1.2), (7.2, 12.1, 0.9)):
        eta, r, a = eta0, r0, float(A.subs(z, eta0))
        numbers = {m: mass, X[2]: 1.1, X[3]: 0.3}
        J = sp.eye(4)
        if system in ("conformal", "cosmic"):
            numbers.update({X[0]: eta, X[1]: r})
            scale = {P["a"]: A.subs(z, X[0])}
            if system == "cosmic":
                J[0, 0] = 1 / a
        elif system == "areal":
            numbers.update({X[0]: eta, X[1]: a * r})
            scale = {P["a"]: A.subs(z, X[0])}
            adot = float(sp.diff(A, z).subs(z, eta))
            J[0, 0], J[1, 0], J[1, 1] = 1 / a, -r * adot / a, 1 / a
        elif system == "isotropic":
            sigma = (r - mass + (r * (r - 2 * mass)) ** 0.5) / 2
            numbers.update({X[0]: eta, X[1]: sigma})
            scale = {P["a"]: A.subs(z, X[0])}
            J[1, 1] = 1 - mass ** 2 / (4 * sigma ** 2)
        elif system == "tortoise":
            import math
            numbers.update({X[0]: eta, X[1]: r + 2 * mass * math.log(r / (2 * mass) - 1), P["r"]: r})
            scale = {P["a"]: A.subs(z, X[0])}
            J[1, 1] = 1 - 2 * mass / r
        else:
            import math
            x = r + 2 * mass * math.log(r / (2 * mass) - 1)
            U, V = -math.exp(-(eta - x) / (4 * mass)), math.exp((eta + x) / (4 * mass))
            a = eta ** 2
            numbers.update({X[0]: U, X[1]: V, P["r"]: r, P["eta"]: eta, P["eta_0"]: 1})
            scale = {}
            f = 1 - 2 * mass / r
            J[0, 0], J[0, 1] = -2 * mass / U, 2 * mass / V
            J[1, 0], J[1, 1] = 2 * mass * f / U, 2 * mass * f / V
        f = 1 - 2 * mass / r
        old = sp.diag(-a ** 2 * f, a ** 2 / f, a ** 2 * r ** 2, a ** 2 * r ** 2 * sp.sin(1.1) ** 2)
        def at(value):
            return float(sp.sympify(value).subs(scale).doit().subs(numbers).evalf())
        got = chart.geo.g.applyfunc(at)
        wanted = J.T * old * J
        if max(abs(float(v)) for v in got - wanted) > 1e-8 * max(1, max(abs(float(v)) for v in wanted)):
            raise AssertionError(f"thakurta/{system}: conformal pullback at {eta0}, {r0}")
        if system == "conformal":
            ap = float(sp.diff(A, z).subs(z, eta))
            app = float(sp.diff(A, z, 2).subs(z, eta))
            G = sp.Matrix(chart.geo.einstein_ll()).applyfunc(at)
            wanted_G = sp.diag(3 * ap ** 2 / a ** 2, (ap ** 2 - 2 * a * app) / (a ** 2 * f ** 2),
                              r ** 2 * (ap ** 2 - 2 * a * app) / (a ** 2 * f),
                              r ** 2 * sp.sin(1.1) ** 2 * (ap ** 2 - 2 * a * app) / (a ** 2 * f))
            wanted_G[0, 1] = wanted_G[1, 0] = 2 * mass * ap / (a * r ** 2 * f)
            if max(abs(float(v)) for v in G - wanted_G) > 1e-9:
                raise AssertionError("thakurta: Mello's density, isotropic pressure and radial flux")
