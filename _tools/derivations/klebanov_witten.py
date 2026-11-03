"""The Klebanov-Witten throat, from hep-th/9807080, equations (4) and (10).

The conifold and its base are separate Riemannian geometries, used to derive
and test the ten-dimensional throat, not additional charts of spacetime.
"""

METRIC = "klebanov_witten"
ANGLES = [r"\theta_1", r"\phi_1", r"\theta_2", r"\phi_2", r"\psi"]
BASE_LINE = (r"\dfrac{1}{9}\left(d\psi+\cos\theta_1\,d\phi_1+\cos\theta_2\,d\phi_2\right)^2"
             r"+\dfrac{1}{6}\left(d\theta_1^2+\sin^2\theta_1\,d\phi_1^2"
             r"+d\theta_2^2+\sin^2\theta_2\,d\phi_2^2\right)")
DIMENSIONS = {(METRIC, "poincare"): {"t": "L", "x": "L", "y": "L", "z": "L", "r": "L", "L": "L",
                                      **{a: "1" for a in ANGLES}}}


def base_metric(angles):
    """The unit Einstein metric in the angular order theta1, phi1, theta2, phi2, psi."""
    import sympy as sp
    a, p, b, q, w = angles
    connection = sp.Matrix([0, sp.cos(a), 0, sp.cos(b), 1])
    return sp.diag(sp.Rational(1, 6), sp.sin(a)**2/6, sp.Rational(1, 6), sp.sin(b)**2/6, 0) + connection*connection.T/9


def chart_spec():
    """The near-horizon limit of (4), with the normalized angular metric (10)."""
    line = (r"ds^2 = \dfrac{r^2}{L^2}\left(-{}dt^2+dx^2+dy^2+dz^2\right)"
            r"+\dfrac{L^2}{r^2}dr^2+L^2\left(" + BASE_LINE + r"\right)")
    coords = ["t", "x", "y", "z", "r"] + ANGLES
    domains = [a+r" \in (-\infty,\infty)" for a in coords[:4]]
    domains += [r"r \in (0,\infty)", r"\theta_1 \in [0,\pi]", r"\phi_1 \in [0,2\pi)",
                r"\theta_2 \in [0,\pi]", r"\phi_2 \in [0,2\pi)", r"\psi \in [0,4\pi)"]
    return {"metric_id": METRIC, "system": {"id": "poincare", "name": "Poincaré",
            "coords": coords, "domains": domains, "parameters": ["L"], "line_element": line.replace("{}", "c^2")},
            "chart_line_element": line.replace("{}", ""), "printer": {},
            "ricci_scalar": "0", "kretschmann": r"\dfrac{176}{L^4}", "check": check_chart}


def check_chart(chart):
    """Check the sourced product and its opposite Ricci eigenvalues at three regular points."""
    import sympy as sp
    r, L = chart.symbols[4], chart.reader.parameters["L"]
    expected = sp.diag(-r*r/L**2, r*r/L**2, r*r/L**2, r*r/L**2, L*L/(r*r),
                       L*L*base_metric(chart.symbols[5:]))
    ricci = sp.Matrix(chart.geo.ricci_ll())
    targets = sp.diag(*([-4/L**2]*5 + [4/L**2]*5))*expected
    expressions = list(chart.geo.g-expected) + list(ricci-targets)
    expressions += [chart.geo.ricci_scalar(), chart.geo.kretschmann()-176/L**4]
    for k in range(3):
        at = {q: sp.Rational(j+k+2, j+3) for j, q in enumerate(chart.symbols)}
        at[L] = sp.Rational(k+3, 2)
        for value in expressions:
            if abs(complex(sp.N(value.subs(at), 35))) > 1e-24:
                raise AssertionError(f"{METRIC}/poincare: sourced product or curvature mismatch")
