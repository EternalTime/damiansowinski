"""Source-defined localization backgrounds from 1982 to 1987.

Flat defects, the six dimensional vacuum, and the electric solution retain
their distinct field equations. Source locators are in early_braneworlds.md.
"""
import math

import numpy as np
import sympy as sp

METRIC = "early_braneworlds"
SYSTEMS = ("vortex_cartesian", "vortex_polar", "scalar_wall", "vacuum_six",
           "visser_lapse", "visser_electric", "squires_three", "squires_family",
           "squires_axis", "magnetic_flat", "magnetic_general", "magnetic_curved")
REAL = r" \in (-\infty,\infty)"
POSITIVE = r" \in (0,\infty)"
DOTS = (r" We take components in $x^0=ct$; a dot denotes an affine derivative, "
        r"so $\dot t=d(ct)/d\lambda$ for affine parameter $\lambda$.")


def definition(system):
    """The source line element and its physical assumptions, before tensor generation."""
    coords = ["t", "x", "y", "z", "w"]
    parameters = {}
    domains = [s + REAL for s in coords]
    convention = ("We use time $t$ and Cartesian lengths $x$, $y$, and $z$ along the defect, "
                  "with transverse length $w$.")
    if system.startswith("vortex"):
        if system == "vortex_cartesian":
            coords = ["t", "x", "y", "z", "u", "v"]
            line = r"-c^2dt^2+dx^2+dy^2+dz^2+du^2+dv^2"
            domains = [s + REAL for s in coords]
            name = "Akama, Cartesian Background"
            convention = ("We use time $t$ and Cartesian lengths $x$, $y$, and $z$ along the vortex, "
                          "with transverse lengths $u$ and $v$. The six dimensional background is flat.")
        else:
            coords = ["t", "x", "y", "z", r"\rho", r"\phi"]
            line = r"-c^2dt^2+dx^2+dy^2+dz^2+d\rho^2+\rho^2d\phi^2"
            domains = [s + REAL for s in coords[:4]] + [r"\rho \in (0,\infty)", r"\phi \in [0,2\pi)"]
            name = "Akama, Transverse Polar Coordinates"
            convention = (r"We use time $t$, Cartesian lengths $x$, $y$, and $z$ along the vortex, "
                          r"transverse radius $\rho$, and angle $\phi$, with $u=\rho\cos\phi$ and "
                          r"$v=\rho\sin\phi$. The axis is regular.")
    elif system == "scalar_wall":
        name = "Scalar Wall Background"
        line = r"-c^2dt^2+dx^2+dy^2+dz^2+dw^2"
        convention += " Gravitational backreaction of the scalar wall is neglected."
    elif system == "vacuum_six":
        name = "Rubakov and Shaposhnikov, Six Dimensional Vacuum"
        coords = ["t", "x", "y", "z", r"\rho", r"\phi"]
        parameters = {
            "L": r"Radial length scale, with cosmological constant $\Lambda=8/(5L^2)$ in $G_{AB}+\Lambda g_{AB}=0$",
            r"f = \cos^{2/5}(\rho/L)": "Scale factor of the four dimensional Minkowski factor",
            r"h = L\sin(\rho/L)\cos^{-3/5}(\rho/L)": "Circumferential radius of a transverse circle, a length",
        }
        line = r"f^2\left(-c^2dt^2+dx^2+dy^2+dz^2\right)+d\rho^2+h^2d\phi^2"
        domains = [s + REAL for s in coords[:4]] + [r"\rho \in (0,\pi L/2)", r"\phi \in [0,2\pi)", "L>0"]
        convention = (r"We use time $t$, Cartesian lengths $x$, $y$, and $z$, proper radial length $\rho$, "
                      r"and angle $\phi$. The regular axis is at $\rho=0$ and the singular endpoint at "
                      r"$\rho=\pi L/2$. We reverse Lavrelashvili's overall signature.")
    elif system.startswith("squires"):
        parameters = {"L": r"Positive radial length scale; the five dimensional vacuum has $\Lambda=-3/(2L^2)$ and the three dimensional vacuum has $\Lambda=-1/L^2$"}
        if system == "squires_three":
            name = "Squires, Three Dimensional Vacuum"
            coords = ["t", "w", r"\phi"]
            line = r"-\cosh^2(w/L)c^2dt^2+dw^2+L^2\sinh^2(w/L)d\phi^2"
            convention = r"We use time $t$, proper radial length $w$, and angle $\phi$, scaled to a regular axis."
        else:
            coords = ["t", "x", "y", "w", r"\phi"]
            convention = r"We use time $t$, Cartesian lengths $x$ and $y$, proper radial length $w$, and angle $\phi$."
            if system == "squires_axis":
                name = "Squires, Regular Axis"
                line = r"\cosh(w/L)\left(-c^2dt^2+dx^2+dy^2\right)+dw^2+\dfrac{L^2\sinh^2(w/L)}{\cosh(w/L)}d\phi^2"
            else:
                name = "Squires, Constrained Vacuum Family"
                for i in range(4):
                    parameters[f"q_{i}"] = (r"Dimensionless exponent; $q_0+q_1+q_2+q_3=0$ and "
                                             r"$q_0^2+q_1^2+q_2^2+q_3^2=3/4$")
                    parameters[rf"A_{i} = \left(\sinh(w/L)\right)^{{1/4+q_{i}}}\left(\cosh(w/L)\right)^{{1/4-q_{i}}}"] = "Dimensionless tangential scale factor"
                line = r"-A_0^2c^2dt^2+A_1^2dx^2+A_2^2dy^2+dw^2+L^2A_3^2d\phi^2"
        domains = [s + (REAL if s in ("t", "x", "y") else POSITIVE if s == "w" else r" \in [0,2\pi)") for s in coords] + ["L>0"]
    elif system.startswith("magnetic"):
        if system == "magnetic_flat":
            name = "Gibbons and Wiltshire, Flat Magnetic Membrane"
            coords = ["t", "x", "y", "z", r"\rho", r"\phi"]
            parameters = {"a": r"Magnetic length, with central field $B_0$ fixed by $a^2=2/(3\pi G B_0^2)$",
                          r"H = 1+\rho^2/a^2": "Dimensionless magnetic warp factor"}
            line = r"H^{2/3}\left(-c^2dt^2+dx^2+dy^2+dz^2+d\rho^2\right)+\rho^2H^{-2}d\phi^2"
            domains = [s + REAL for s in coords[:4]] + [r"\rho \in (0,\infty)", r"\phi \in [0,2\pi)", "a>0"]
            convention = (r"We use time $t$, Cartesian lengths $x$, $y$, and $z$, transverse radius $\rho$, "
                          r"and angle $\phi$. The magnetic two form is $F=B_0\rho H^{-2}d\rho\wedge d\phi$ "
                          r"with Maxwell stress normalized by $8\pi G$.")
        else:
            curved = system == "magnetic_curved"
            name = "Gibbons and Wiltshire, " + ("Curved Membrane" if curved else "Radial Warp Coordinates")
            coords = ["t", "r", r"\theta", r"\psi", r"\beta", r"\phi"] if curved else ["t", "x", "y", "z", r"\beta", r"\phi"]
            parameters = {"C": r"Angular scale, a squared length; at a simple regular root $\beta_0$ of $D$ we set $C=2/|\partial_\beta D(\beta_0)|$",
                          r"\mu": r"Integration constant, an inverse squared length, $\mu=2Gk$ with integration constant $k$",
                          "Q": r"Magnetic parameter, an inverse length, $Q^2=2\pi G B_0^2/3$ with magnetic field parameter $B_0$",
                          r"\Lambda": "Bulk cosmological constant, an inverse squared length"}
            if curved:
                parameters["K"] = r"Membrane sectional curvature, an inverse squared length; its Ricci tensor is $\bar R_{ab}=3K\bar g_{ab}$"
                leaf = r"-\left(1-Kr^2\right)c^2dt^2+\dfrac{dr^2}{1-Kr^2}+r^2d\theta^2+r^2\sin^2\theta\,d\psi^2"
            else:
                leaf = r"-c^2dt^2+dx^2+dy^2+dz^2"
            parameters[r"D = " + ("K+" if curved else "") + r"\mu/\beta^3-Q^2/\beta^6-\Lambda\beta^2/10"] = "Radial function, an inverse squared length"
            line = r"\beta^2\left(" + leaf + r"\right)+\dfrac{d\beta^2}{D}+C^2D\,d\phi^2"
            domains = (["t"+REAL, "r"+POSITIVE, r"\theta \in (0,\pi)", r"\psi \in [0,2\pi)", "1-Kr^2>0"] if curved else [s+REAL for s in coords[:4]])
            domains += [r"\beta>0", "D>0", r"\phi \in [0,2\pi)", "C>0"]
            convention = (r"We use time $t$, " + (r"areal length $r$ and angles $\theta$ and $\psi$ on the maximally symmetric membrane, " if curved else r"Cartesian lengths $x$, $y$, and $z$ on the flat membrane, ") +
                          r"dimensionless warp coordinate $\beta$, and angle $\phi$, with $\tilde\phi=C\phi$. "
                          r"The magnetic two form is $F=CB_0\beta^{-4}d\beta\wedge d\phi$ with Maxwell stress normalized by $8\pi G$.")
    else:
        if system == "visser_lapse":
            name = "Visser, Lapse Ansatz"
            parameters = {"N = N(w)": "Positive dimensionless lapse; the supporting stress is fixed by Einstein's equations"}
            line = r"-N^2c^2dt^2+dx^2+dy^2+dz^2+dw^2"
            convention += " We leave the lapse free and impose no matter equation."
        else:
            name = "Visser, Electric Trapping"
            parameters = {"L": r"Curvature length, with normalized electric field $E=1/L$, potential $A_t=\sinh(w/L)$, and cosmological constant $\Lambda=-1/(2L^2)$ in $G_{AB}+\Lambda g_{AB}=T_{AB}$"}
            line = r"-\cosh^2(w/L)c^2dt^2+dx^2+dy^2+dz^2+dw^2"
            domains += ["L>0"]
            convention += " A constant electric field and a cosmological constant support the varying clock rate."
    return {"id": system, "name": name, "coords": coords, "domains": domains,
            "parameters": [{"symbol": k, "description": v} for k, v in parameters.items()],
            "convention": convention + DOTS, "line_element": "ds^2=" + line}


def spec(system, pc=None):
    """A chart with source checks made only at numerical points."""
    import verify_metrics as vm
    d = definition(system)
    parameters = [p["symbol"] for p in d["parameters"]]
    reader = vm.Reader(d["coords"], parameters, ())
    return {"metric_id": METRIC,
            "system": {k: v for k, v in d.items() if k != "convention" and k != "parameters"}
                      | {"parameters": parameters},
            "chart_line_element": d["line_element"].replace("c^2", ""),
            "printer": {"lead": [*reader.parameters.values(), *reader.symbol.values()], "flip": False},
            "check": lambda chart: check_chart(chart, system)}


def jets(reader, variables, point):
    """Held functions and their derivatives, evaluated from their own definitions."""
    out = {}
    for f, expression in reader.held.items():
        out[f] = expression.subs(point).evalf(40)
        for a in variables:
            out[sp.diff(f, a)] = sp.diff(expression, a).subs(point).evalf(40)
            for b in variables:
                out[sp.diff(f, a, b)] = sp.diff(expression, a, b).subs(point).evalf(40)
    return out


def check_chart(chart, system):
    """Source Ricci and stress tensors, including the flat polar pullback."""
    x, p = chart.symbols, chart.reader.parameters
    ricci = sp.Matrix(chart.geo.ricci_ll())
    for a in (sp.Rational(2, 5), sp.Rational(4, 5), sp.Rational(6, 5)):
        point = {s: a for s in x}
        if "L" in p:
            point[p["L"]] = sp.Rational(7, 5)
        for i, value in enumerate((-sp.Rational(1,4), -sp.Rational(1,4), -sp.Rational(1,4), sp.Rational(3,4))
                                  if a == sp.Rational(2,5) else
                                  (sp.sqrt(3)/4, sp.sqrt(3)/4, -sp.sqrt(3)/4, -sp.sqrt(3)/4)):
            if f"q_{i}" in p:
                point[p[f"q_{i}"]] = value
        if system.startswith("magnetic"):
            values = {"a": sp.Rational(7,5), "C": 1, "mu": 2, "Q": sp.Rational(2,5),
                      "Lambda": -sp.Rational(1,10), "K": sp.Rational(1,5)}
            point.update({p[k]: v for k,v in values.items() if k in p})
            if system != "magnetic_flat":
                point[x[4]] = 1+a
        substitutions = jets(chart.reader, x, point)
        if system == "visser_lapse":
            N, w = p["N"], x[-1]
            substitutions.update({N: sp.cosh(a), sp.diff(N, w): sp.sinh(a), sp.diff(N, w, 2): sp.cosh(a)})
        def at(value):
            return np.array(sp.Matrix(value).xreplace(substitutions).subs(point).evalf(30), dtype=float)
        g, R = at(chart.geo.g), at(ricci)
        want = np.zeros_like(g)
        if system == "vacuum_six":
            want = 4 / (5 * float(point[p["L"]])**2) * g
        elif system.startswith("squires"):
            want = -(2 if system == "squires_three" else 1) / float(point[p["L"]])**2 * g
        elif system.startswith("magnetic"):
            # The source Maxwell form, rescaled by sqrt(8*pi*G), with d=6.
            radial, angular = 4, 5
            F = np.zeros((6,6))
            if system == "magnetic_flat":
                radius, scale = float(a), float(point[p["a"]])
                F[radial, angular] = 4/math.sqrt(3)/scale * radius/(1+radius**2/scale**2)**2
                cosmological = 0
            else:
                F[radial, angular] = math.sqrt(12)*float(point[p["Q"]])*float(point[p["C"]])/float(point[x[4]])**4
                cosmological = float(point[p["Lambda"]])
            F[angular, radial] = -F[radial, angular]
            inverse = np.linalg.inv(g)
            contraction = float(np.sum(F*(inverse @ F @ inverse)))
            # Trace-reversed Einstein equation R_ab = Lambda*g_ab/2 + F_ac F_b^c - g_ab F^2/8.
            want = cosmological*g/2 + F @ inverse @ F.T - g*contraction/8
        elif system in ("visser_lapse", "visser_electric"):
            L = float(point.get(p.get("L"), 1))
            want[0, 0] = -g[0, 0] / L**2
            want[-1, -1] = -1 / L**2
        if np.max(np.abs(R-want)) > 1e-11 * max(1, np.max(np.abs(want))):
            raise AssertionError(f"{system}: source Ricci tensor at {a}: {np.max(np.abs(R-want))}")
        if system == "vortex_polar":
            rho, phi = float(a), float(a)
            J = np.eye(6)
            J[4:, 4:] = [[math.cos(phi), -rho*math.sin(phi)], [math.sin(phi), rho*math.cos(phi)]]
            if np.max(np.abs(J.T @ np.diag([-1,1,1,1,1,1]) @ J-g)) > 1e-12:
                raise AssertionError("Akama: Cartesian to polar pullback")


def register(pc):
    pc.CHARTS[METRIC] = [lambda s=s: spec(s, pc) for s in SYSTEMS]


def register_dimensions(vm):
    for system in SYSTEMS:
        d = definition(system)
        dims = {s: ("T" if s == "t" else "1" if s in (r"\phi",r"\theta",r"\psi",r"\beta") else "L") for s in d["coords"]}
        dims.update({"L": "L", "f": "1", "h": "L", "N": "1", "a": "L", "H": "1",
                     "C": "L**2", r"\mu": "1/L**2", "Q": "1/L", r"\Lambda": "1/L**2", "K": "1/L**2", "D": "1/L**2"})
        dims.update({f"q_{i}": "1" for i in range(4)})
        dims.update({f"A_{i}": "1" for i in range(4)})
        used = {vm.Reader._plain(s) for s in d["coords"]}
        used.update(vm.Reader._plain(p["symbol"].partition("=")[0].strip()) for p in d["parameters"])
        vm.DIMENSIONS[(METRIC, system)] = {k:v for k,v in dims.items() if vm.Reader._plain(k) in used}
    vm.HELD[(METRIC, "vacuum_six")] = ("f", "h")
    vm.HELD[(METRIC, "squires_family")] = tuple(f"A_{i}" for i in range(4))
    for system in ("magnetic_general", "magnetic_curved"):
        vm.HELD[(METRIC, system)] = ("D",)
