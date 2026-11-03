"""The static magnetic exterior in Eiroa and Simeone (36)-(39)."""
import sympy as sp

METRIC = "superconducting_string"
SYSTEM = "cylindrical"
COORDS = ["t", "r", r"\phi", "z"]
F = r"F = \dfrac{(r/\ell)^m+k(r/\ell)^{-m}}{1+k}"
H = r"H = \dfrac{(r/\ell)^m-k(r/\ell)^{-m}}{1+k}"
P = r"P = (r/\ell)^{2m^2}"
PARAMETERS = ["m", "k", "W", r"\ell", "r_c", F, H, P]
LINE = r"ds^2 = F^2\left(P\left(-dt^2+dr^2\right)+W^2r^2d\phi^2\right)+F^{-2}dz^2"
DIMENSIONS = {(METRIC, SYSTEM): {"t":"L", "r":"L", r"\phi":"1", "z":"L",
    "m":"1", "k":"1", "W":"1", r"\ell":"L", "F":"1", "H":"1", "P":"1", "r_c":"L"}}
HELD = {(METRIC, SYSTEM): ("F", "H", "P")}
RATES = {(METRIC, SYSTEM): {"F": {"r":"mH/r"}, "H": {"r":"mF/r"}, "P": {"r":"2m^2P/r"}}}


def chart_spec():
    """The normalized cylindrical exterior, with radial rates held explicitly."""
    import verify_metrics as vm
    probe = vm.Reader(COORDS, PARAMETERS, (), held=HELD[(METRIC, SYSTEM)], rates=RATES[(METRIC, SYSTEM)])
    return {"metric_id": METRIC, "system": {
        "id": SYSTEM, "name": "Cylindrical Exterior", "coords": COORDS,
        "parameters": PARAMETERS,
        "domains": [r"t\in(-\infty,\infty)", r"r\in(r_c,\infty)",
                    r"\phi\in[0,2\pi)", r"z\in(-\infty,\infty)",
                    r"\ell>0,\quad r_c>0,\quad W>0,\quad k\geq 0,\quad -1<m<1"],
        "line_element": LINE}, "chart_line_element": LINE,
        "printer": {}, "reduce": probe.surface, "check": check_chart}


def check_chart(chart):
    """Compare the published metric with the source at sample points."""
    t, r, phi, z = chart.symbols
    p = chart.reader.parameters
    for rv, ellv, mv, kv, wv in ((1.5, 1.0, 0.2, 0.7, 0.8),
                                 (3.2, 1.3, -0.3, 1.4, 0.9),
                                 (2.0, 0.8, 0.4, 0.2, 0.7)):
        x = rv / ellv
        f = (x**mv + kv*x**(-mv))/(1+kv)
        a = x**(2*mv*mv)*f*f
        expected = (-a, a, wv*wv*rv*rv*f*f, 1/(f*f))
        at = {t: 0, r: rv, phi: 0, z: 0, p['m']: mv, p['k']: kv,
              p['W']: wv, p['ell']: ellv, p['r_c']: 0.5}
        for i in range(4):
            got = chart.geo.g[i,i].subs(chart.reader.held).doit().subs(at).evalf(30)
            if abs(float(got)-expected[i]) > 1e-10:
                raise AssertionError(f"{METRIC}/{SYSTEM}: source metric differs at {rv}, component {i}")
