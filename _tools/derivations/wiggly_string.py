"""Camps and Emparan (2010), equation (2.8), without higher-order completion."""

METRIC = "wiggly_string"
SYSTEM = "cylindrical"
DIMENSIONS = {(METRIC, SYSTEM): {"t": "L", "r": "L", r"\phi": "1", "z": "L",
                                 "U": "1", "T": "1", "r_0": "L"}}


def chart_spec():
    """Print exact tensors of the explicitly truncated, sourced line element."""
    line = (r"ds^2 = -\left(1+4(U-T)\ln(r/r_0)\right)dt^2"
            r"+\left(1-4(U+T)\ln(r/r_0)\right)\left(dr^2+r^2d\phi^2\right)"
            r"+\left(1-4(U-T)\ln(r/r_0)\right)dz^2")
    import sympy as sp
    import verify_metrics as vm
    probe = vm.Reader(["t", "r", r"\phi", "z"], ["U", "T", "r_0"], ())
    logarithm = sp.log(probe.symbol["r"]/probe.parameters["r_0"])
    return {"metric_id": METRIC, "system": {
        "id": SYSTEM, "name": "Cylindrical Transverse Gauge",
        "coords": ["t", "r", r"\phi", "z"],
        "domains": [r"t \in (-\infty,\infty)", r"r \in (0,\infty)",
                    r"\phi \in [0,2\pi)", r"z \in (-\infty,\infty)",
                    r"0<T\le U\ll1", r"\left|4(U+T)\ln(r/r_0)\right|\ll1"],
        "parameters": ["U", "T", "r_0"], "line_element": line},
        "chart_line_element": line, "printer": {"overrides": {logarithm: r"\ln(r/r_0)"}}, "check": check_chart}


def check_chart(chart):
    """Compare with the source coefficient by coefficient at independent points."""
    import sympy as sp
    t, r, phi, z = chart.symbols
    U, T, r0 = (chart.reader.parameters[key] for key in ("U", "T", "r_0"))
    f, q = 4*(U-T)*sp.log(r/r0), 4*(U+T)*sp.log(r/r0)
    expected = sp.diag(-1-f, 1-q, r*r*(1-q), 1-f)
    for radius, length, energy, tension in ((2, 1, sp.Rational(1, 100), sp.Rational(1, 400)),
                                           (3, 2, sp.Rational(1, 200), sp.Rational(1, 300))):
        at = {t: 0, r: radius, phi: sp.Rational(2, 7), z: 0,
              r0: length, U: energy, T: tension}
        for value in chart.geo.g-expected:
            if abs(complex(sp.N(value.subs(at), 30))) > 1e-22:
                raise AssertionError("wiggly_string: equation (2.8) coefficient differs")
