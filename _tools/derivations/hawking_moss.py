"""The barrier-top sphere and its exact constant-field continuations.

Brown and Weinberg (2007), equations (3.6), (5.3), (5.6), (5.7), and footnote 7.
All coordinates with dimensions of length use c = hbar = 1.
"""
import sympy as sp

METRIC = "hawking_moss"
EUCLIDEAN = ("four_sphere", "thermal", "thermal_angle")
LORENTZIAN = ("static_continuation", "wall_continuation", "closed_continuation")
SYSTEMS = EUCLIDEAN + LORENTZIAN
ANGLES = [r"\theta \in [0, \pi]", r"\phi \in [0, 2\pi)"]
S2 = r"\left(d\theta^2 + \sin^2\theta\,d\phi^2\right)"
S3 = r"\left(d\chi^2 + \sin^2\chi\,d\theta^2 + \sin^2\chi\sin^2\theta\,d\phi^2\right)"


def spec(system):
    """One sourced coordinate system, with its numerical ambient-space pullback check."""
    sphere = system in ("four_sphere", "closed_continuation")
    if sphere:
        first = r"\xi" if system == "four_sphere" else "t"
        coords = [first, r"\chi", r"\theta", r"\phi"]
        domains = ([r"\xi \in [0, \pi\ell]"] if system == "four_sphere" else
                   [r"t \in (-\infty, \infty)"]) + [r"\chi \in [0, \pi]"] + ANGLES
        line = (r"ds^2 = d\xi^2 + \ell^2\sin^2(\xi/\ell)" + S3 if system == "four_sphere" else
                r"ds^2 = -dt^2 + \ell^2\cosh^2(t/\ell)" + S3)
    elif system == "wall_continuation":
        coords = [r"\xi", r"\psi", r"\theta", r"\phi"]
        domains = [r"\xi \in (0, \pi\ell)", r"\psi \in (-\infty, \infty)"] + ANGLES
        line = r"ds^2 = d\xi^2 + \ell^2\sin^2(\xi/\ell)\left(-d\psi^2 + \cosh^2\psi" + S2 + r"\right)"
    else:
        angle = system == "thermal_angle"
        first = "t" if system == "static_continuation" else r"\tau"
        radial = r"\eta" if angle else "r"
        coords = [first, radial, r"\theta", r"\phi"]
        domains = ([r"t \in (-\infty, \infty)"] if first == "t" else
                   [r"\tau \in [0, 2\pi\ell)"])
        domains += ([r"\eta \in [0, \pi/2]"] if angle else
                    [r"r \in [0, \ell)" if first == "t" else r"r \in [0, \ell]"]) + ANGLES
        sign = "-" if first == "t" else ""
        line = (r"ds^2 = \cos^2\eta\,d\tau^2 + \ell^2d\eta^2 + \ell^2\sin^2\eta" + S2 if angle else
                rf"ds^2 = {sign}\left(1 - r^2/\ell^2\right)d{first}^2"
                + r" + \dfrac{dr^2}{1 - r^2/\ell^2} + r^2" + S2)
    names = {"four_sphere": "Euclidean 4-Sphere", "thermal": "Euclidean Thermal Coordinates",
             "thermal_angle": "Euclidean Thermal Angle", "static_continuation": "Static Constant Field Continuation",
             "wall_continuation": "Wall Constant Field Continuation",
             "closed_continuation": "Closed Constant Field Continuation"}
    return {"metric_id": METRIC, "system": {"id": system, "name": names[system], "coords": coords,
            "domains": domains + [r"\ell > 0"], "parameters": [r"\ell"], "line_element": line},
            "chart_line_element": line, "printer": {},
            "ricci_scalar": r"\dfrac{12}{\ell^2}", "kretschmann": r"\dfrac{24}{\ell^4}",
            "check": lambda chart: check_chart(chart, system)}


def ambient(chart, system):
    """Five-dimensional sphere or de Sitter hyperboloid, with the ambient time first."""
    x, y, th, ph = chart.symbols
    ell = chart.reader.parameters["ell"]
    n2 = sp.Matrix([sp.cos(th), sp.sin(th)*sp.cos(ph), sp.sin(th)*sp.sin(ph)])
    if system in ("four_sphere", "closed_continuation"):
        n3 = sp.Matrix([sp.cos(y), *(sp.sin(y)*n2)])
        if system == "four_sphere":
            return sp.Matrix([ell*sp.cos(x/ell), *(ell*sp.sin(x/ell)*n3)])
        return sp.Matrix([ell*sp.sinh(x/ell), *(ell*sp.cosh(x/ell)*n3)])
    if system == "wall_continuation":
        return sp.Matrix([ell*sp.sin(x/ell)*sp.sinh(y), ell*sp.cos(x/ell),
                          *(ell*sp.sin(x/ell)*sp.cosh(y)*n2)])
    r = ell*sp.sin(y) if system == "thermal_angle" else y
    root = sp.sqrt(ell**2-r**2)
    if system == "static_continuation":
        return sp.Matrix([root*sp.sinh(x/ell), root*sp.cosh(x/ell), *(r*n2)])
    return sp.Matrix([root*sp.cos(x/ell), root*sp.sin(x/ell), *(r*n2)])


def check_chart(chart, system):
    """Point checks, including the radius constraint and Ricci equation, without symbolic reduction."""
    ell = chart.reader.parameters["ell"]
    image = ambient(chart, system)
    flat = sp.diag(-1 if system in LORENTZIAN else 1, 1, 1, 1, 1)
    jac = image.jacobian(chart.symbols)
    pulled = jac.T*flat*jac
    ricci = sp.Matrix(chart.geo.ricci_ll())
    for length, u, v, theta, phi in ((0.8, 0.32, 0.43, 0.7, 0.4), (1.3, 0.81, 0.62, 1.1, 2.1),
                                   (2.1, 1.1, 0.27, 2.3, 4.4)):
        # Both thermal coordinates stay inside the regular open part of the sphere.
        if system in ("thermal", "static_continuation"):
            v *= length
        if system in ("four_sphere", "wall_continuation"):
            u *= length
        at = {**dict(zip(chart.symbols, (u, v, theta, phi))), ell: length, chart.reader.c: 1}
        for label, values in (("ambient pullback", list(chart.geo.g-pulled)),
                              ("ambient radius", [(image.T*flat*image)[0]-ell**2]),
                              ("Einstein equation", list(ricci-3*chart.geo.g/ell**2))):
            error = max(abs(complex(sp.N(value.subs(at), 35))) for value in values)
            if error > 1e-9:
                raise AssertionError(f"{METRIC}/{system}: {label}: {error}")


def register(pc):
    pc.CHARTS[METRIC] = [lambda s=s: spec(s) for s in SYSTEMS]


def register_dimensions(vm):
    for system in SYSTEMS:
        coords = spec(system)["system"]["coords"]
        vm.DIMENSIONS[(METRIC, system)] = {**{c: "L" if c in ("t", "r", r"\tau", r"\xi") else "1"
                                               for c in coords}, r"\ell": "L"}
