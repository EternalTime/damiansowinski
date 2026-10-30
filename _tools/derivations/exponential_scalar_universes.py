"""Exact flat exponential scalar cosmologies, in units c = 1.

Copeland, Liddle and Wands (1998), (1)-(4), Table I and (A7)-(A8), and
Heard and Wands (2002), (10)-(11) and (22)-(24), fix the two matter branches.
Abbott and Wise (1984), (2.1)-(2.3), fix the inflationary conformal map.
The continuation to deceleration and the logarithmic coasting map follow
by integrating dt/a. See exponential_scalar_universes.md for the domains.
"""

METRIC = "exponential_scalar_universes"
CHARTS = ("cosmic_cartesian", "cosmic_spherical", "conformal", "coasting")
SCALES = {
    "cosmic_cartesian": r"a = a_0\left(t/t_0\right)^p",
    "cosmic_spherical": r"a = a_0\left(t/t_0\right)^p",
    "conformal": r"a = a_0\left((1-p)a_0\eta/t_0\right)^{p/(1-p)}",
    "coasting": r"a = a_0e^{a_0\eta/t_0}",
}
COORDS = {s: (["t", "r", r"\theta", r"\phi"] if s == "cosmic_spherical" else
              ["t" if s.startswith("cosmic") else r"\eta", "x", "y", "z"]) for s in CHARTS}
DIMENSIONS = {(METRIC, s): {**{x: ("1" if x in (r"\theta", r"\phi") else "L") for x in COORDS[s]},
                            "p": "1", "a_0": "1", "t_0": "L", "a": "1"} for s in CHARTS}
HELD = {(METRIC, s): ("a",) for s in CHARTS}
RATES = {(METRIC, s): {"a": {COORDS[s][0]: ("pa/t" if s.startswith("cosmic") else
                    r"pa/((1-p)\eta)" if s == "conformal" else "a_0a/t_0")}} for s in CHARTS}


def chart_spec(system):
    """Four source-defined charts with a positive scale factor held by its slope."""
    import verify_metrics as vm

    cosmic = system.startswith("cosmic")
    params = (["p"] if system != "coasting" else []) + ["a_0", "t_0", SCALES[system]]
    if system == "cosmic_spherical":
        line = r"ds^2 = -dt^2+a^2\left(dr^2+r^2d\theta^2+r^2\sin^2\theta\,d\phi^2\right)"
        space = [r"r \in [0,\infty)", r"\theta \in [0,\pi]", r"\phi \in [0,2\pi)"]
    else:
        line = (r"ds^2 = -dt^2+a^2\left(dx^2+dy^2+dz^2\right)" if cosmic else
                r"ds^2 = a^2\left(-d\eta^2+dx^2+dy^2+dz^2\right)")
        space = [x+r" \in (-\infty,\infty)" for x in ("x", "y", "z")]
    domain = ([r"t \in (0,\infty)"] if cosmic else
              [r"\eta \in (-\infty,0)\;\text{for}\;p>1",
               r"\eta \in (0,\infty)\;\text{for}\;p<1"] if system == "conformal" else
              [r"\eta \in (-\infty,\infty)"])
    key = (METRIC, system)
    probe = vm.Reader(COORDS[system], params, (), held=HELD[key], rates=RATES[key])
    return {"metric_id": METRIC,
            "system": {"id": system, "name": {"cosmic_cartesian": "Cosmic Cartesian",
                        "cosmic_spherical": "Cosmic Spherical", "conformal": "Conformal Time",
                        "coasting": "Coasting Conformal Time"}[system],
                       "coords": COORDS[system], "parameters": params, "domains": domain+space,
                       "line_element": line},
            "chart_line_element": line, "printer": {},
            "reduce": lambda value: vm.norm(probe.by_rates(value)),
            "check": lambda chart: chart_check(chart, system)}


def chart_check(chart, system):
    """Point checks of the coordinate maps and the two exact matter branches."""
    # The implementation is completed before the charts are generated.
    raise NotImplementedError(system)
