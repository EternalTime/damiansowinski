"""The first-order collapsing texture, Nötzold (10) and Durrer et al. (0.46).

The displayed tensors belong exactly to the truncated metric.
The texture stress matches its Einstein tensor only through first order.
"""
import sympy as sp

ANGULAR = r"H = 1 - 2\epsilon + \dfrac{2\epsilon t}{r}\arctan\left(\dfrac{r}{t}\right)"
LINE = r"ds^2=-dt^2+dr^2+r^2H\left(d\theta^2+\sin^2\theta\,d\phi^2\right)"


def chart_spec():
    import verify_metrics as vm
    coords = ["t", "r", r"\theta", r"\phi"]
    probe = vm.Reader(coords, [r"\epsilon", ANGULAR], (), held=("H",))
    H = probe.parameters["H"]
    r = probe.symbol["r"]
    return {
        "metric_id": "global_texture",
        "system": {"id": "collapse", "name": "Before Collapse", "coords": coords,
                   "domains": [r"t \in (-\infty, 0)", r"r \in (0, \infty)",
                               r"\theta \in [0, \pi]", r"\phi \in [0, 2\pi)"],
                   "parameters": [r"\epsilon", ANGULAR], "line_element": LINE},
        "chart_line_element": LINE,
        "printer": {"lead": [H, r], "flip": False},
        "check": chart_check,
    }


def chart_check(chart):
    """Compare the linear Einstein coefficient to the flat-space sigma-model stress.

    All comparisons are numerical at rational points, without symbolic determinants
    or simplification of a check expression.
    """
    t, r, theta, _ = chart.symbols
    epsilon = chart.reader.parameters["epsilon"]
    explicit = chart.reader.held
    g0 = sp.diag(-1, 1, r**2, r**2 * sp.sin(theta)**2)
    chi = sp.pi + 2 * sp.atan(t/r)
    dt, dr = sp.diff(chi, t), sp.diff(chi, r)
    inner = sp.zeros(4)
    inner[0, 0], inner[0, 1], inner[1, 0], inner[1, 1] = dt**2, dt*dr, dt*dr, dr**2
    inner[2, 2] = sp.sin(chi)**2
    inner[3, 3] = sp.sin(chi)**2 * sp.sin(theta)**2
    trace = -dt**2 + dr**2 + 2*sp.sin(chi)**2/r**2
    stress = inner - g0*trace/2
    G = chart.geo.einstein_ll()
    points = [(-3, sp.Rational(2, 3)), (-1, 2), (sp.Rational(-1, 4), 3)]
    for a in range(4):
        for b in range(a, 4):
            exact = G[a][b].subs(explicit).doit()
            leading = sp.diff(exact, epsilon).subs(epsilon, 0)
            for time, radius in points:
                at = {t: time, r: radius, theta: sp.Rational(4, 5)}
                if abs(sp.N((leading-stress[a, b]).subs(at), 40)) > sp.Float('1e-28'):
                    raise AssertionError(f"global_texture: first-order Einstein equation fails in {a}{b}")
                if abs(sp.N(exact.subs(epsilon, 0).subs(at), 40)) > sp.Float('1e-28'):
                    raise AssertionError("global_texture: flat limit is not vacuum")
    gamma = chart.geo.christoffel_ull()
    for a in (0, 1):
        for b in (0, 1):
            for c in (0, 1):
                if gamma[a][b][c] != 0:
                    raise AssertionError("global_texture: radial geodesics are not straight")
