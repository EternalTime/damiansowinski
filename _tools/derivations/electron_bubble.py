"""The rotating geometric shell, from gr-qc/0109085v2, (94) and (108).

The shell is outside every smooth chart's open domain.
The electromagnetic qualification and the Israel tensor are derived in
electron_bubble.md; these bulk tensors do not contain either surface layer.
"""

METRIC = "electron_bubble"
CHARTS = ("exterior", "interior", "exterior_kerr", "interior_kerr", "cartesian")
PARAMETERS = ["m", "q", "a", "r_0 = q^2/(2m)"]
SIGMA = r"\Sigma = r^2 + a^2\cos^2\theta"
B = "B = r^2 + a^2"
DELTA = r"\Delta = r^2 + a^2 - 2mr + q^2"
DIMENSIONS = {
    (METRIC, system): {
        **({"T": "T", "x": "L", "y": "L", "z": "L"} if system == "cartesian" else
           {("T" if system.endswith("kerr") else "t"): "T", "r": "L", r"\theta": "1",
            (r"\psi" if system.endswith("kerr") else r"\phi"): "1",
            r"\Sigma": "L**2", "B": "L**2",
            **({r"\Delta": "L**2"} if system == "exterior" else {})}),
        "m": "L", "q": "L", "a": "L", "r_0": "L",
    } for system in CHARTS
}


def line_element(system, clock):
    """The source's f=mr-q²/2 or f=0, with x^0=c times the printed time."""
    if system == "cartesian":
        return f"ds^2 = -{clock}dT^2 + dx^2 + dy^2 + dz^2"
    kerr = system.endswith("kerr")
    exterior = system.startswith("exterior")
    time, angle = ("T", r"\psi") if kerr else ("t", r"\phi")
    c = "c\\," if clock else ""
    metric = f"ds^2 = -{clock}d{time}^2"
    radial = "dr^2" if kerr else (r"\dfrac{\Sigma}{\Delta}dr^2" if exterior else r"\dfrac{\Sigma}{B}dr^2")
    metric += f" + {radial}" + r" + \Sigma d\theta^2 + B\sin^2\theta\,d" + angle + "^2"
    if kerr:
        metric += r" - 2a\sin^2\theta\,dr\,d\psi"
    if exterior:
        metric += r" + \dfrac{2mr - q^2}{\Sigma}\left(" + c + f"d{time}"
        metric += " + dr" if kerr else ""
        metric += r" - a\sin^2\theta\,d" + angle + r"\right)^2"
    return metric


def chart_spec(system):
    """Five source-defined charts and a numerical comparison with its matrix."""
    import verify_metrics as vm

    cartesian = system == "cartesian"
    kerr = system.endswith("kerr")
    exterior = system.startswith("exterior")
    time, angle = ("T", r"\psi") if kerr else ("t", r"\phi")
    coords = ["T", "x", "y", "z"] if cartesian else [time, "r", r"\theta", angle]
    parameters = PARAMETERS + ([] if cartesian else [SIGMA, B] + ([DELTA] if system == "exterior" else []))
    names = {"exterior": "Boyer-Lindquist Exterior", "interior": "Flat Oblate Interior",
             "exterior_kerr": "Kerr Angular Exterior", "interior_kerr": "Flat Kerr Angular Interior",
             "cartesian": "Flat Cartesian Interior"}
    if cartesian:
        domains = [f"{v} \\in (-\\infty, \\infty)" for v in coords]
        domains += [r"\dfrac{x^2 + y^2}{r_0^2 + a^2} + \dfrac{z^2}{r_0^2} < 1"]
    else:
        domains = [time + r" \in (-\infty, \infty)",
                   r"r \in (r_0, \infty)" if exterior else r"r \in (0, r_0)",
                   r"\theta \in [0, \pi]", angle + r" \in [0, 2\pi)"]
    probe = vm.Reader(coords, parameters, ())
    return {
        "metric_id": METRIC,
        "system": {"id": system, "name": names[system], "coords": coords,
                   "parameters": parameters, "domains": domains,
                   "line_element": line_element(system, "c^2")},
        "chart_line_element": line_element(system, ""),
        "printer": {"lead": list(probe.symbol.values()), "flip": False},
        "ricci_scalar": "0",
        "check": lambda chart: chart_check(chart, system),
    }


def chart_check(chart, system):
    """Source matrices at points, including both charge signs and a static control.

    No symbolic determinant or tensor recomputation is needed here.
    The angular-chart Jacobian is the Appendix C transformation.
    """
    import numpy as np

    reader = chart.reader
    exterior = system.startswith("exterior")
    for m, q, a, theta in ((1, 2, 3, .7), (2, -3, 1, 1.1), (1, 1.5, 0, .9)):
        r0 = q*q/(2*m)
        r = r0 * (1.4 if exterior else .6)
        at = {reader.parameters[k]: v for k, v in (("m", m), ("q", q), ("a", a))}
        if system == "cartesian":
            want = np.diag([-1., 1., 1., 1.])
        else:
            at.update({reader.symbol["r"]: r, reader.symbol[r"\theta"]: theta})
            sine = np.sin(theta)**2
            sigma, b = r*r + a*a*np.cos(theta)**2, r*r + a*a
            f2 = 2*m*r - q*q if exterior else 0
            want = np.diag([-1., sigma/(b-f2), sigma, b*sine])
            v = np.array([1., 0., 0., -a*sine])
            want += f2/sigma * np.outer(v, v)
            if system.endswith("kerr"):
                jac = np.eye(4)
                jac[0, 1], jac[3, 1] = -f2/(b-f2), -a/(b-f2)
                want = jac.T @ want @ jac
        got = np.array(chart.geo.g.subs(at).evalf(), dtype=float)
        if not np.allclose(got, want, atol=2e-12, rtol=2e-12):
            raise AssertionError(f"{METRIC}/{system}: source matrix or Appendix C pullback disagrees")


def register_charts(pc):
    pc.CHARTS[METRIC] = [lambda system=system: chart_spec(system) for system in CHARTS]
