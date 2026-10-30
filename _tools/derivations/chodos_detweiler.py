"""The compact five-dimensional vacuum of Chodos and Detweiler.

Chyba (1985), page 868, (6.5)-(6.7), fixes the metric, the four periods and
the spatial rescaling. The derivation and source provenance are in
chodos_detweiler.md. Shared generators register these page-local definitions.
"""

METRIC = "chodos_detweiler"
CHARTS = ("synchronous", "rescaled")


def chart_spec(system):
    """Both five-dimensional charts, with tensor time x^0 = ct."""
    import verify_metrics as vm

    rescaled = system == "rescaled"
    space = ["X", "Y", "Z", "W"] if rescaled else ["x", "y", "z", "w"]
    coords = ["t", *space]
    parameters = ["t_0", "L"] + (["\\tau"] if rescaled else [])
    scale = "\\tau" if rescaled else "t_0"
    probe = vm.Reader(coords, parameters, ())
    time = probe.symbol["t"]

    def line(clock, lapse):
        terms = " + ".join(f"d{x}^2" for x in space[:3])
        return (f"ds^2 = -{lapse}dt^2 + \\dfrac{{{clock}}}{{{scale}}}"
                f"\\left({terms}\\right) + \\dfrac{{{scale}}}{{{clock}}}d{space[3]}^2")

    if rescaled:
        growing = "L\\sqrt{\\tau/t_0}"
        shrinking = "L\\sqrt{t_0/\\tau}"
    else:
        growing = shrinking = "L"
    domains = ["t \\in (0, \\infty)"]
    domains += [f"{x} \\in [0, {growing})" for x in space[:3]]
    domains += [f"{space[3]} \\in [0, {shrinking})", "t \\to 0 \\;\\text{(the singularity)}"]
    return {
        "metric_id": METRIC,
        "system": {"id": system, "name": "Rescaled" if rescaled else "Synchronous",
                   "coords": coords, "domains": domains, "parameters": parameters,
                   "line_element": line("t", "c^2")},
        "chart_line_element": line("\\left(t/c\\right)", ""),
        "time": "t",
        "printer": {"lead": [time], "factors": [probe.c, time, *probe.parameters.values()]},
        "check": lambda chart: chart_check(chart, system),
    }


def chart_check(chart, system):
    """The source coefficients and physical periods at distinct rational points.

    Only metric values are sampled here; tensor identities belong to the
    generated mathematics and the scoped metric checker.
    """
    import sympy as sp

    R = chart.reader
    clock = R.symbol["t"]
    for c, t, t0, tau, length in ((2, 3, 2, 5, 7), (3, 5, 7, 2, 11), (5, 2, 3, 7, 13)):
        c, t, t0, tau, length = map(sp.Rational, (c, t, t0, tau, length))
        scale = tau if system == "rescaled" else t0
        at = {clock: c*t, R.c: c, R.parameters["t_0"]: t0, R.parameters["L"]: length}
        if system == "rescaled":
            at[R.parameters["tau"]] = tau
        want = sp.diag(-1, t/scale, t/scale, t/scale, scale/t)
        got = chart.geo.g.subs(at)
        if got != want:
            raise AssertionError(f"{METRIC}/{system}: source coefficients disagree at {at}")
        periods = ([length*sp.sqrt(tau/t0)]*3 + [length*sp.sqrt(t0/tau)]
                   if system == "rescaled" else [length]*4)
        squares = [got[i, i]*periods[i-1]**2 for i in range(1, 5)]
        if squares != [length**2*t/t0]*3 + [length**2*t0/t]:
            raise AssertionError(f"{METRIC}/{system}: physical circumferences disagree at {at}")
