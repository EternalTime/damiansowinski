"""Static de Sitter regions of the charged membrane bubble.

Hirano, arXiv:1804.09985, equations (2.28)-(2.29) and (3.3).
The wall is excluded from the smooth bulk tensors and is matched by Israel's condition.
"""
import sympy as sp

SYSTEMS = ("static_inside", "static_outside", "sliced_inside", "sliced_outside",
           "flat_inside", "flat_outside", "conformal_inside", "conformal_outside")


def chart_spec(system):
    import verify_metrics as vm

    if system.startswith("sliced_"):
        return sliced_spec(system)
    if system.startswith(("flat_", "conformal_")):
        return expanding_spec(system)
    inside = system == "static_inside"
    lam = r"\Lambda_i" if inside else r"\Lambda_o"
    coords = ["t", "r", r"\theta", r"\phi"]
    parameters = [lam, "r_w = r_w(t)"]
    f = rf"1 - \dfrac{{{lam} r^2}}{{3}}"
    sphere = r"r^2\left(d\theta^2 + \sin^2\theta\,d\phi^2\right)"
    line = rf"ds^2 = -\left({f}\right)c^2dt^2 + \dfrac{{dr^2}}{{{f}}} + {sphere}"
    probe = vm.Reader(coords, parameters, ())
    return {
        "metric_id": "brown_teitelboim",
        "system": {
            "id": system,
            "name": "Static, Inside the Membrane" if inside else "Static, Outside the Membrane",
            "coords": coords,
            "parameters": parameters,
            "domains": [r"t \in (-\infty, \infty)",
                        r"r \in [0, r_w)" if inside else r"r \in (r_w, \sqrt{3/\Lambda_o})",
                        r"r < \sqrt{3/\Lambda_i}" if inside else r"\Lambda_o > 0",
                        r"\theta \in [0, \pi]", r"\phi \in [0, 2\pi)",
                        r"r = r_w \;\text{(the membrane)}"],
            "line_element": line,
        },
        "chart_line_element": line.replace("c^2dt^2", "dt^2"),
        "printer": {"lead": [probe.parameters[vm.Reader._plain(lam)], probe.symbol["r"]]},
        "kretschmann": rf"\dfrac{{8{lam}^2}}{{3}}",
        "check": lambda chart: check_chart(chart, system),
    }


def check_chart(chart, system):
    """Check the de Sitter Einstein equation at independent points."""
    lam = chart.reader.parameters["Lambda_i" if system == "static_inside" else "Lambda_o"]
    _, r, theta, phi = chart.symbols
    ricci = chart.geo.ricci_ll()
    for curvature, radius, angle in ((sp.Rational(3, 4), sp.Rational(1, 2), sp.Rational(2, 3)),
                                     (sp.Rational(1, 3), sp.Rational(4, 5), sp.Rational(7, 8))):
        point = {lam: curvature, r: radius, theta: angle, phi: sp.Rational(1, 4)}
        for i in range(4):
            for j in range(4):
                diff = ricci[i][j] - lam * chart.geo.g[i, j]
                if abs(complex(sp.N(diff.subs(point), 35))) > 1e-25:
                    raise AssertionError(f"brown_teitelboim/{system}: Ricci equation {i},{j}")


def sliced_spec(system):
    """Hirano's de Sitter slices, equations (2.16) and (2.19)."""
    import chart_printer as cp
    import verify_metrics as vm

    inside = system == "sliced_inside"
    suffix = "i" if inside else "o"
    ell = rf"\ell_{suffix}"
    wall = rf"\chi_{suffix}"
    coords = [r"\eta", r"\chi", r"\theta", r"\phi"]
    parameters = [ell, wall]
    factor = rf"{ell}^2\sin^2\left(\dfrac{{\chi}}{{{ell}}}\right)"
    sphere = r"d\theta^2 + \sin^2\theta\,d\phi^2"
    line = (rf"ds^2 = d\chi^2 + {factor}\left(-d\eta^2 + "
            rf"\cosh^2\eta\left({sphere}\right)\right)")
    probe = vm.Reader(coords, parameters, ())
    return {
        "metric_id": "brown_teitelboim",
        "system": {"id": system,
                   "name": "de Sitter Slices, " + ("Inside" if inside else "Outside") + " the Membrane",
                   "coords": coords, "parameters": parameters,
                   "domains": [r"\eta \in (-\infty, \infty)",
                               rf"\chi \in (0, {wall})" if inside else rf"\chi \in ({wall}, \pi {ell})",
                               r"\theta \in [0, \pi]", r"\phi \in [0, 2\pi)",
                               rf"\chi = {wall} \;\text{{(the membrane)}}"],
                   "line_element": line},
        "chart_line_element": line,
        "printer": {"lead": [probe.parameters[vm.Reader._plain(ell)]]},
        "pretty": cp.hyperbolic(probe.symbol[r"\eta"]),
        "kretschmann": rf"\dfrac{{24}}{{{ell}^4}}",
        "check": lambda chart: check_sliced(chart, system),
    }


def check_sliced(chart, system):
    ell = chart.reader.parameters["ell_i" if system == "sliced_inside" else "ell_o"]
    eta, chi, theta, phi = chart.symbols
    ricci = chart.geo.ricci_ll()
    for radius, distance, rapidity, angle in ((2, sp.Rational(1, 2), sp.Rational(1, 3), sp.Rational(2, 3)),
                                              (3, 1, sp.Rational(2, 5), sp.Rational(4, 5))):
        point = {ell: radius, chi: distance, eta: rapidity, theta: angle, phi: sp.Rational(1, 5)}
        for i in range(4):
            for j in range(4):
                diff = ricci[i][j] - 3 * chart.geo.g[i, j] / ell**2
                if abs(complex(sp.N(diff.subs(point), 35))) > 1e-25:
                    raise AssertionError(f"brown_teitelboim/{system}: Ricci equation {i},{j}")


def expanding_spec(system):
    """Hirano's expanding and conformally flat charts, equations (2.5)-(2.6)."""
    import verify_metrics as vm

    inside = system.endswith("inside")
    conformal = system.startswith("conformal_")
    suffix = "i" if inside else "o"
    ell = rf"\ell_{suffix}"
    time = r"\zeta" if conformal else "T"
    coords = [time, r"\rho", r"\theta", r"\phi"]
    parameters = [ell, rf"s_w = s_w({time})"]
    spatial = r"d\rho^2 + \rho^2\left(d\theta^2 + \sin^2\theta\,d\phi^2\right)"
    line = (rf"ds^2 = \dfrac{{{ell}^2}}{{\zeta^2}}\left(-d\zeta^2 + {spatial}\right)" if conformal
            else rf"ds^2 = -dT^2 + e^{{2T/{ell}}}\left({spatial}\right)")
    probe = vm.Reader(coords, parameters, ())
    return {
        "metric_id": "brown_teitelboim",
        "system": {"id": system,
                   "name": ("Conformal" if conformal else "Expanding") + ", "
                           + ("Inside" if inside else "Outside") + " the Membrane",
                   "coords": coords, "parameters": parameters,
                   "domains": [rf"{time} \in " + (r"(-\infty, 0)" if conformal else r"(-\infty, \infty)"),
                               r"\rho \in [0, s_w)" if inside else r"\rho \in (s_w, \infty)",
                               r"\theta \in [0, \pi]", r"\phi \in [0, 2\pi)",
                               r"\rho = s_w \;\text{(the membrane)}"],
                   "line_element": line},
        "chart_line_element": line,
        "printer": {"lead": [probe.parameters[vm.Reader._plain(ell)]]},
        "kretschmann": rf"\dfrac{{24}}{{{ell}^4}}",
        "check": lambda chart: check_expanding(chart, system),
    }


def check_expanding(chart, system):
    ell = chart.reader.parameters["ell_i" if system.endswith("inside") else "ell_o"]
    time, radius, theta, phi = chart.symbols
    ricci = chart.geo.ricci_ll()
    for length, instant, r, angle in ((2, -1, sp.Rational(1, 2), sp.Rational(2, 3)),
                                       (3, -2, 1, sp.Rational(4, 5))):
        point = {ell: length, time: instant, radius: r, theta: angle, phi: sp.Rational(1, 5)}
        for i in range(4):
            for j in range(4):
                diff = ricci[i][j] - 3 * chart.geo.g[i, j] / ell**2
                if abs(complex(sp.N(diff.subs(point), 35))) > 1e-25:
                    raise AssertionError(f"brown_teitelboim/{system}: Ricci equation {i},{j}")
