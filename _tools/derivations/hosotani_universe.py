"""The fixed magnetic SU(2) cosmology, Gibbons and Steif (2.1)-(2.12).

The source normalization and the original full-text debt are in hosotani_universe.md.
The two charts describe one expanding and recollapsing solution.
"""
import sympy as sp


SYSTEMS = ("cosmic_hyperspherical", "conformal_hyperspherical")
SPATIAL = (r"d\chi^2+\sin^2\chi\left(d\theta^2+\sin^2\theta\,d\phi^2\right)")


def chart_spec(system):
    cosmic = system == SYSTEMS[0]
    time = "t" if cosmic else r"\eta"
    line = (rf"ds^2=-dt^2+t(2\ell-t)\left({SPATIAL}\right)" if cosmic else
            rf"ds^2=\ell^2\sin^2\eta\left(-d\eta^2+{SPATIAL}\right)")
    return {"metric_id": "hosotani_universe", "system": {
        "id": system, "name": ("Cosmic" if cosmic else "Conformal") + " Hyperspherical",
        "coords": [time, r"\chi", r"\theta", r"\phi"], "parameters": [r"\ell"],
        "domains": [r"t \in (0, 2\ell)" if cosmic else r"\eta \in (0, \pi)",
                    r"\chi \in [0, \pi]", r"\theta \in [0, \pi]", r"\phi \in [0, 2\pi)"],
        "line_element": line}, "chart_line_element": line, "printer": {},
        "check": lambda chart: check_chart(chart, system)}


def check_chart(chart, system):
    """Point checks of the radiation source, curvature, and cosmic/conformal pullback."""
    X, P = chart.symbols, chart.reader.parameters
    time, chi, theta, phi = X
    ell = P["ell"]
    cosmic = system == SYSTEMS[0]
    a2 = time * (2 * ell - time) if cosmic else ell**2 * sp.sin(time)**2
    mixed = chart.geo.ginv * sp.Matrix(chart.geo.einstein_ll())
    expected = sp.diag(-3, 1, 1, 1) * ell**2 / a2**2
    # The conformal chart pulled back with dt/deta=a gives the cosmic chart.
    eta = sp.Symbol("eta_probe", real=True)
    tc = ell * (1 - sp.cos(eta))
    scale2 = ell**2 * sp.sin(eta)**2
    spatial = sp.diag(1, sp.sin(chi)**2, sp.sin(chi)**2 * sp.sin(theta)**2)
    target = sp.diag(-scale2, *(scale2 * spatial[i, i] for i in range(3)))
    pulled = chart.geo.g
    if cosmic:
        J = sp.diag(ell * sp.sin(eta), 1, 1, 1)
        pulled = J.T * pulled.subs(time, tc) * J
    else:
        pulled = pulled.subs(time, eta)
    curvature = [chart.geo.ricci_scalar(), chart.geo.kretschmann() - 24 * ell**4 / a2**4]
    for ev, lv in ((0.43, 0.8), (1.57, 1.3), (2.6, 2.1)):
        at = {ell: lv, time: lv * (1 - sp.cos(ev)) if cosmic else ev,
              eta: ev, chi: 0.87, theta: 1.1, phi: 0.2}
        for label, gap in (("Einstein radiation source", list(mixed - expected)),
                           ("chart pullback", list(pulled - target)), ("curvature", curvature)):
            error = max(abs(complex(sp.N(v.subs(at), 35))) for v in gap)
            if error > 1e-8:
                raise AssertionError(f"hosotani_universe/{system}: {label}: {error}")
