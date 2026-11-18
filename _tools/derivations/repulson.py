"""Kallosh and Linde's repulsons, hep-th/9507022, with the string coupling at infinity g = 1 and
G = 1, in their signature (+,-,-,-) turned to (-,+,+,+).

The Einstein metric is their eq. (36), TeX label II, ds^2 = -dt^2/sqrt(Delta) + sqrt(Delta) dx^2
with Delta = 1 + 4M/r - 4q^2/r^2, and the string metric is their eq. (43), TeX label stringyelectric,
ds^2 = -dt^2/Delta + dx^2, the Einstein metric times e^{2 phi} with e^{-2 phi} = sqrt(Delta), their
four-dimensional dilaton. Delta is held as a function of r by its slope. repulson.md beside this
file is the derivation, and records why Cvetic and Youm's form of the same solution is the
isotropic chart itself and not another chart.
"""

METRIC = "repulson"
CHARTS = ("isotropic", "string")
COORDS = ["t", "r", r"\theta", r"\phi"]
DELTA = r"\Delta = 1 + \dfrac{4M}{r} - \dfrac{4q^2}{r^2}"
DIMENSIONS = {(METRIC, s): {"t": "T", "r": "L", r"\theta": "1", r"\phi": "1", "M": "L", "q": "L",
                            r"\Delta": "1"} for s in CHARTS}
HELD = {(METRIC, s): ("Delta",) for s in CHARTS}
RATES = {(METRIC, s): {"Delta": {"r": r"-\dfrac{4\left(Mr - 2q^2\right)}{r^3}"}} for s in CHARTS}
SPHERE = r"r^2\left(d\theta^2 + \sin^2\theta\,d\phi^2\right)"
NAMES = {"isotropic": "Einstein Metric, Isotropic", "string": "String Metric"}
DOMAINS = [r"t \in (-\infty, \infty)", r"r \in \left(2\sqrt{M^2 + q^2} - 2M, \infty\right)", r"\theta \in [0, \pi]",
           r"\phi \in [0, 2\pi)", r"r = r_0 \;\text{(the singularity)}"]


def line(system, c2):
    """The line element of a chart, with c2 standing before dt^2."""
    if system == "isotropic":
        return (r"ds^2 = -\dfrac{" + c2 + r"dt^2}{\sqrt{\Delta}} + \sqrt{\Delta}\left(dr^2 + " + SPHERE + r"\right)")
    return r"ds^2 = -\dfrac{" + c2 + r"dt^2}{\Delta} + dr^2 + " + SPHERE


def chart_spec(system):
    """One of the two charts, its values printed in Delta and r."""
    import verify_metrics as vm

    key = (METRIC, system)
    parameters = ["M", "q", DELTA]
    probe = vm.Reader(COORDS, parameters, (), held=HELD[key], rates=RATES[key])
    r, Delta = probe.symbol["r"], probe.parameters["Delta"]
    return {"metric_id": METRIC,
            "system": {"id": system, "name": NAMES[system], "coords": COORDS, "domains": DOMAINS,
                       "parameters": parameters, "line_element": line(system, "c^2")},
            "chart_line_element": line(system, ""),
            "printer": {"lead": [Delta, r], "factors": [r, Delta], "flip": False},
            "reduce": lambda value: vm.norm(probe.by_rates(value)),
            "check": lambda chart: chart_check(chart, system)}


def chart_check(chart, system):
    """Each chart at points against what its sources state.

    The Einstein metric is Cvetic and Youm's extreme hole of four charges, hep-th/9507160 eq. (2),
    published as the four_extreme chart of string_bh_three_four_charges, with the two magnetic charge
    radii zero and the two electric ones 2M + 2 sqrt(M^2 + q^2) and 2M - 2 sqrt(M^2 + q^2), the second
    negative, which is their condition for a singular, repelling state. With M = 0 its Ricci scalar is
    Kallosh's, hep-th/9506113 App. B, 4q^2(2q^2 + r^2)/(r (r^2 - 4q^2)^{5/2}), and the string metric's is
    her -8q^2(8q^2 + r^2)/(r^2 (r^2 - 4q^2)^2). The string metric is the Einstein metric divided by
    sqrt(Delta), and Delta vanishes at r_0 = 2(sqrt(M^2 + q^2) - M), Kallosh and Linde's eq. (41)."""
    import json
    from pathlib import Path

    import sympy as sp
    import verify_metrics as vm

    X = chart.symbols
    P = chart.reader.parameters
    t, r, th, ph = X
    M, q = P["M"], P["q"]
    held = chart.reader.held
    g = sp.Matrix(chart.geo.g).xreplace(held)
    Delta = 1 + 4 * M / r - 4 * q ** 2 / r ** 2
    root = sp.sqrt(M ** 2 + q ** 2)
    points = [{M: sp.Rational(3, 10), q: sp.Rational(7, 10), r: sp.Rational(29, 10), th: sp.Rational(1, 3)},
              {M: sp.Rational(1, 2), q: sp.Rational(1, 5), r: sp.Rational(41, 10), th: sp.Rational(9, 10)},
              {M: 0, q: sp.Rational(3, 4), r: sp.Rational(37, 10), th: sp.Rational(2, 5)}]

    def vanishes(e, what):
        for at in points:
            if abs(complex(sp.N(sp.sympify(e).subs(at), 40))) > 1e-30:
                raise AssertionError(f"{METRIC}/{system}: {what} fails at {at}")

    einstein = sp.diag(-1 / sp.sqrt(Delta), sp.sqrt(Delta), sp.sqrt(Delta) * r ** 2,
                       sp.sqrt(Delta) * r ** 2 * sp.sin(th) ** 2)
    if system == "isotropic":
        for i in range(4):
            vanishes(g[i, i] - einstein[i, i], f"Kallosh and Linde's eq. (36), slot {i}{i}")
        there = next(c for c in json.loads((Path(vm.METRICS_DIR) / "string_bh_three_four_charges.json")
                                           .read_text(encoding="utf-8"))["coordinates"] if c["id"] == "four_extreme")
        other = vm.Reader(there["coords"], [p["symbol"] for p in there["parameters"]], ())
        theirs = vm.metric_from_line_element(other, there["line_element"].replace("c^2", ""), there["coords"])
        names = dict(zip([other.symbol[c] for c in other.coords], X))
        names.update({other.parameters["r_1"]: 2 * M + 2 * root, other.parameters["r_2"]: 2 * M - 2 * root,
                      other.parameters["r_3"]: 0, other.parameters["r_4"]: 0})
        pulled = theirs.subs(names, simultaneous=True)
        for i in range(4):
            vanishes(pulled[i, i] - g[i, i], f"Cvetic and Youm's extreme hole with opposite charges, slot {i}{i}")
        scalar = sp.sympify(chart.geo.ricci_scalar()).xreplace(held)
        kallosh = 4 * q ** 2 * (2 * q ** 2 + r ** 2) / (r * (r ** 2 - 4 * q ** 2) ** sp.Rational(5, 2))
    else:
        for i in range(4):
            vanishes(g[i, i] - einstein[i, i] / sp.sqrt(Delta), f"the Einstein metric over sqrt(Delta), slot {i}{i}")
        scalar = sp.sympify(chart.geo.ricci_scalar()).xreplace(held)
        kallosh = -8 * q ** 2 * (8 * q ** 2 + r ** 2) / (r ** 2 * (r ** 2 - 4 * q ** 2) ** 2)
    at = points[2]
    if abs(complex(sp.N((scalar - kallosh).subs(at), 40))) > 1e-30:
        raise AssertionError(f"{METRIC}/{system}: the Ricci scalar at M = 0 is not Kallosh's")
    for at in points:
        r0 = (2 * (root - M)).subs(at)
        if abs(sp.N(Delta.subs({**at, r: r0}), 40)) > 1e-30:
            raise AssertionError(f"{METRIC}/{system}: Delta does not vanish at r_0")
