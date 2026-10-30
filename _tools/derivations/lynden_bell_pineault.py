"""Charts and cold-disc boundary equations of Lynden-Bell and Pineault (1978).

The Weyl chart includes the finite discs of Paper I section 3.2.
The other charts describe infinite discs, on the open upper vacuum side.
Reflection across the equator supplies the lower side and its thin matter source.
"""
import math
from functools import lru_cache

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
            line = ("ds^2 = -f^{2n}e^{N}{c2}dt^2 + A^2f^{2 - 2n}e^{2P - N}"
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
        powers = pc.named_powers(names, state)
        if system == "rotating":
            def pretty(value):
                exponential, rest = sp.Integer(1), sp.Integer(1)
                for factor in sp.Mul.make_args(sp.factor(sp.sympify(value))):
                    if factor.func == sp.exp:
                        exponential *= factor
                    else:
                        rest *= factor
                return exponential*powers(rest)
        else:
            pretty = powers
    else:
        pretty = lambda value: sp.powsimp(sp.factor(sp.sympify(value)), combine="exp")

    def check(chart):
        state["printer"] = chart.printer
        check_chart(chart, system)

    out = {"metric_id": METRIC,
            "system": {"id": system, "name": name, "coords": coords, "domains": domains,
                       "parameters": parameters, "line_element": line.replace("{c2}", "c^2").replace("{c}", "c")},
            "chart_line_element": line.replace("{c2}", "").replace("{c}", ""),
            "printer": {"lead": [*probe.parameters.values(), *probe.symbol.values()], "flip": False},
            "pretty": pretty, "check": check}
    if system == "spherical":
        out["kretschmann"] = ("\\dfrac{16n^2(2n - 1)^2f^{4n}h^{4n - 8n^2}"
                              "\\left(2 + \\cos\\theta + n(2n - 1)(1 - \\cos\\theta)\\right)}"
                              "{r^4(1 + \\cos\\theta)^2}")
    return out


def check_chart(chart, system):
    """Compare the source line elements and exterior Ricci tensors at separate events.

    No determinant or symbolic simplification belongs in these source checks.
    Free-function charts are evaluated on an exact static member or the cold rotating solution.
    """
    x, p = chart.symbols, chart.reader.parameters
    A = p["A"]
    for radius, angle, exponent in [(1.3, 0.4, 0.2), (2.1, 0.9, 0.12), (0.8, 1.2, 0.31)]:
        at = {A: 1.0, x[0]: 0.0, x[3]: 0.2}
        if system == "weyl":
            R, z, n = x[1], x[2], p["n"]
            rr = sp.sqrt(R**2 + z**2)
            nu = n * sp.log((rr + z) / (2*A))
            xi = 2*n**2 * sp.log((1 + z/rr)/2)
            at.update({R: radius*math.sin(angle), z: radius*math.cos(angle), n: exponent})
            jets = function_jets({p["nu"]: nu, p["xi"]: xi}, [R, z], at)
        elif system == "rotating":
            solution = rotating_solution()
            n = solution["n"]
            at.update({x[1]: radius, x[2]: angle, p["n"]: n})
            jets = rotating_jets(chart.reader, angle, solution)
            N, P, Q, Z = (jets[p[name]] for name in ("N", "P", "Q", "Z"))
            f = radius
            # Paper II (2.2), k = 1-n, independently expanded in the source tetrad.
            want = np.diag([f**(2*n)*math.exp(N)*(Q*Q-1),
                            f**(-2*n)*math.exp(Z-N),
                            f**(2-2*n)*math.exp(Z-N),
                            f**(2-2*n)*math.exp(2*P-N)])
            want[0, 3] = want[3, 0] = -f*math.exp(P)*Q
            actual = np.array(chart.geo.g.xreplace(jets).subs(at), dtype=float)
            if np.max(np.abs(actual-want)) > 2e-11:
                raise AssertionError("the rotating chart differs from Paper II (2.2)")
        else:
            n = 0.5 if system == "photons" else exponent
            at.update({x[1]: radius, x[2]: angle})
            if system == "spherical":
                at[p["n"]] = n
            jets = {}
            h = math.cos(angle/2)**2
            # Paper I (3.16), cos(theta/2)^(8n^2), with (4.2) at n=1/2.
            want = np.diag([-(radius*h)**(2*n),
                            (radius*h)**(-2*n)*math.cos(angle/2)**(8*n*n),
                            radius**2*(radius*h)**(-2*n)*math.cos(angle/2)**(8*n*n),
                            radius**2*math.sin(angle)**2*(radius*h)**(-2*n)])
            if np.max(np.abs(np.array(chart.geo.g.subs(at), dtype=float)-want)) > 2e-11:
                raise AssertionError("the static chart differs from Paper I (3.16)")
        ricci = sp.Matrix(chart.geo.ricci_ll()).xreplace(jets).subs(at)
        if max(abs(float(v.evalf())) for v in ricci) > (2e-8 if system == "rotating" else 2e-10):
            raise AssertionError(f"{system}: smooth exterior Ricci tensor is nonzero")


def function_jets(functions, variables, at):
    """Values and derivatives through second order for a pointwise field-equation check."""
    out = {}
    for function, expression in functions.items():
        out[function] = expression.subs(at).evalf(35)
        for a in variables:
            out[sp.diff(function, a)] = sp.diff(expression, a).subs(at).evalf(35)
            for b in variables:
                out[sp.diff(function, a, b)] = sp.diff(expression, a, b).subs(at).evalf(35)
    return out


def rotating_rhs(theta, y, n):
    """Paper II vacuum equations with k+n=1 and e^P=sin(theta).

    Cai-Shu (4-10),(4-11) have Theta=theta_LP and their exponent n_CS=n/(1-n).
    Z' is Paper II (2.21); no WKB expression is used.
    """
    N, dN, Q, dQ, Z = y
    cot = 1/math.tan(theta)
    b = 1-2*n
    S = dQ-Q*cot+Q*dN
    ddN = Q*Q*b*b+S*S-2*n-dN*cot
    ddQ = Q*((1-Q*Q)*b*b+(dN-cot)**2-S*S)-dQ*cot
    dZ = b*Q*S+2*n*dN
    return np.array([dN, ddN, dQ, ddQ, dZ])


def integrate_rotating(n, eta, tolerance=2e-12):
    """Regular-axis expansion followed by adaptive integration to the equatorial sheet."""
    from scipy.integrate import solve_ivp
    eps = 1e-5
    y0 = [-n*eps*eps/2, -n*eps, eta*eps, eta, -n*n*eps*eps]
    result = solve_ivp(lambda theta, y: rotating_rhs(theta, y, n), (eps, math.pi/2), y0,
                       rtol=tolerance, atol=tolerance/100, dense_output=True)
    if not result.success:
        raise AssertionError(result.message)
    return result


@lru_cache(maxsize=4)
def rotating_solution(v=0.2):
    """Shoot on n and Q'(0) using Paper II (2.25) and the cold-disc jump conditions."""
    from scipy.optimize import root

    def residual(values):
        n, eta = values
        y = integrate_rotating(n, eta).y[:, -1]
        N, dN, Q, dQ, Z = y
        return [(1-n)*v*v+(1-2*n)*Q*v-n, dQ+dN*(Q+2*v/(1+v*v))]

    fitted = root(residual, [v*v, 2*v**3], tol=1e-10)
    if not fitted.success or max(abs(x) for x in residual(fitted.x)) > 1e-10:
        raise AssertionError(f"cold-disc boundary solve failed: {fitted.message}")
    n, eta = map(float, fitted.x)
    return {"n": n, "eta": eta, "v": v, "solution": integrate_rotating(n, eta)}


def rotating_jets(reader, theta, solution=None):
    """Angular values and first two derivatives, in the original 1978 radial chart."""
    solution = solution or rotating_solution()
    n = solution["n"]
    y = solution["solution"].sol(theta)
    d = rotating_rhs(theta, y, n)
    N, dN, Q, dQ, Z = y
    cot = 1/math.tan(theta)
    S = dQ-Q*cot+Q*dN
    dS = d[3]-dQ*cot+Q/(math.sin(theta)**2)+dQ*dN+Q*d[1]
    ddZ = (1-2*n)*(dQ*S+Q*dS)+2*n*d[1]
    vals = {"N": (N, dN, d[1]), "P": (math.log(math.sin(theta)), cot, -1/math.sin(theta)**2),
            "Q": (Q, dQ, d[3]), "Z": (Z, d[4], ddZ)}
    out = {}
    th = reader.symbol["\\theta"]
    for name, jet in vals.items():
        function = reader.parameters[name]
        for order, value in enumerate(jet):
            out[sp.diff(function, th, order)] = float(value)
    return out


def register(pc):
    """Register only this spacetime with the shared chart printer."""
    pc.CHARTS[METRIC] = [lambda s=s: spec(s) for s in SYSTEMS]


def register_dimensions():
    """Dimensions used by the published-value readers and the scoped metric checker."""
    for system in SYSTEMS:
        values = {"t": "T", "A": "L", "n": "1"}
        if system == "weyl":
            values.update({"R": "L", "z": "L", "\\phi": "1", "\\nu": "1", "\\xi": "1"})
        else:
            values.update({"r": "L", "\\theta": "1", "\\phi": "1", "f": "1", "h": "1"})
            if system == "rotating":
                values.update({name: "1" for name in ("N", "P", "Q", "Z")})
        vm.DIMENSIONS[(METRIC, system)] = values
