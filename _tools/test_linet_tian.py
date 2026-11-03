"""Linet and Tian's static cylinder with a cosmological constant: what its History states, checked
on the published charts at points of Kasner's circle.

The published curvature is held to the line element by verify_metrics.py and the tags by the data
build; these tests read the published values themselves and hold them to the physics the prose
states: the Einstein space condition R_mu nu = Lambda g_mu nu, the exponents of the History adding
to zero, the hyperbolic form for a negative constant, Levi-Civita's cylinder as Lambda -> 0, the
split into Petrov types I and D, and a Kretschmann scalar that never takes the value 8 Lambda^2/3 of
de Sitter or anti-de Sitter space.
"""
import json
import sys
import unittest
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "_tools" / "derivations"))
import verify_metrics as vm

METRIC = json.loads((ROOT / "MFS" / "assets" / "data" / "metrics" / "linet_tian.json").read_text())
CHARTS = {chart["id"]: chart for chart in METRIC["coordinates"]}

# Parameters and radii at which each chart is read: Lambda = +-1/3 puts the far axis of the toroidal
# chart at rho = pi, and L = 3 gives Sarioglu and Tekin's chart Lambda = -1/3.
POINTS = {
    "toroidal": ({"Lambda": sp.Rational(1, 3), "B": sp.Rational(5, 4), "C": sp.Rational(3, 2)}, [0.4, 1.3, 2.9]),
    "cylindrical": ({"Lambda": -sp.Rational(1, 3), "C": sp.Rational(3, 2)}, [0.4, 1.3, 4.0]),
    "gleiser": ({"Lambda": sp.Rational(1, 3), "C": sp.Rational(3, 2)}, [0.3, 1.5, 2.8]),
    "sarioglu_tekin": ({"L": sp.Integer(3), "r_0": sp.Integer(1), "C": sp.Rational(3, 2)}, [1.1, 2.0, 6.0]),
}
# Levi-Civita's mass parameter at the five members of type D and at two of type I.
TYPE_D = [sp.Rational(-1, 2), sp.Integer(0), sp.Rational(1, 4), sp.Rational(1, 2), sp.Integer(1)]
TYPE_I = [sp.Rational(1, 5), sp.Rational(37, 100)]
TOLERANCE = sp.Float("1e-25")


def kasner(sigma):
    """Levi-Civita's Kasner exponents at mass parameter sigma."""
    S = 4 * sigma ** 2 - 2 * sigma + 1
    return {"p_0": 2 * sigma / S, "p_2": (1 - 2 * sigma) / S, "p_3": 2 * sigma * (2 * sigma - 1) / S}


class Chart:
    """One published chart, read through the checker's reader with its defined names written out."""

    def __init__(self, chart_id):
        self.entry = CHARTS[chart_id]
        self.id = chart_id
        self.reader = vm.Reader(self.entry["coords"], [p["symbol"] for p in self.entry["parameters"]], ())
        self.symbols = [self.reader.symbol[name] for name in self.entry["coords"]]
        n = len(self.symbols)
        self.g = sp.zeros(n)
        for component in self.entry["metric_components"]:
            i, j = (self.entry["coords"].index(x) for x in component["indices"])
            self.g[i, j] = self.reader(component["value"])

    def tensor(self, field, variant, rank):
        out = {}
        for component in self.entry[field]["variants"][variant]["nonzero"]:
            index = tuple(self.entry["coords"].index(x) for x in component["indices"])
            out[index] = self.reader(component["value"])
        return out

    def scalar(self, field):
        return self.reader(self.entry[field].partition("=")[2])

    def at(self, sigma, radius, **overrides):
        values, _ = POINTS[self.id]
        point = {self.reader.parameters[name]: value for name, value in {**values, **overrides}.items()}
        point.update({self.reader.parameters[name]: value for name, value in kasner(sigma).items()})
        point.update({self.symbols[0]: 0, self.symbols[1]: sp.Float(radius, 40), self.symbols[2]: 0, self.symbols[3]: 0})
        return point

    def value(self, expression, point):
        return sp.N(sp.sympify(expression).subs(self.reader.c, 1).subs(point), 40)

    def lam(self, point):
        if self.id == "sarioglu_tekin":
            return -3 / point[self.reader.parameters["L"]] ** 2
        return point[self.reader.parameters["Lambda"]]


class EinsteinSpace(unittest.TestCase):
    def test_the_published_ricci_tensor_is_lambda_times_the_metric(self):
        for chart_id in CHARTS:
            chart = Chart(chart_id)
            ricci = chart.tensor("ricci_tensor", "ll", 2)
            for sigma in (sp.Integer(0), sp.Rational(1, 4), sp.Rational(1, 3)):
                for radius in POINTS[chart_id][1]:
                    point = chart.at(sigma, radius)
                    for a in range(4):
                        gap = chart.value(ricci.get((a, a), 0) - chart.lam(point) * chart.g[a, a], point)
                        self.assertLess(abs(complex(gap)), 1e-25, (chart_id, sigma, radius, a))
                    self.assertLess(abs(complex(chart.value(chart.scalar("ricci_scalar") - 4 * chart.lam(point), point))),
                                    1e-25, chart_id)

    def test_the_exponents_of_the_history_add_to_zero(self):
        sigma = sp.Symbol("sigma")
        S = 4 * sigma ** 2 - 2 * sigma + 1
        history = {"t": -(4 * sigma ** 2 - 8 * sigma + 1) / (3 * S), "z": (8 * sigma ** 2 - 4 * sigma - 1) / (3 * S),
                   "phi": -2 * (2 * sigma ** 2 + 2 * sigma - 1) / (3 * S)}
        self.assertEqual(sp.simplify(sum(history.values())), 0)
        p = kasner(sigma)
        # Q^{2/3} P^{2p} = cos^{4/3}(sqrt(3 Lambda) rho/2) P^{2(p + 1/3)}: each is a Kasner exponent less 1/3.
        for direction, name in (("t", "p_0"), ("z", "p_3"), ("phi", "p_2")):
            self.assertEqual(sp.simplify(history[direction] - (p[name] - sp.Rational(1, 3))), 0, direction)
        self.assertEqual(sp.simplify(sum(p.values()) - 1), 0)
        self.assertEqual(sp.simplify(sum(v ** 2 for v in p.values()) - 1), 0)


class Forms(unittest.TestCase):
    def test_the_negative_constant_takes_the_hyperbolic_functions(self):
        """The cylindrical chart's P is the toroidal chart's continued to Lambda < 0, and the two charts
        write the same line element in it."""
        toroidal, cylindrical = Chart("toroidal"), Chart("cylindrical")
        P_tan = toroidal.reader("P")
        P_tanh = cylindrical.reader("P")
        lam = -sp.Rational(1, 3)
        for radius in (0.4, 1.3, 4.0):
            continued = sp.N(P_tan.subs({toroidal.reader.parameters["Lambda"]: lam, toroidal.symbols[1]: radius}), 40)
            direct = sp.N(P_tanh.subs({cylindrical.reader.parameters["Lambda"]: lam, cylindrical.symbols[1]: radius}), 40)
            self.assertLess(abs(complex(continued - direct)), 1e-30)
        for sigma in (sp.Rational(1, 5), sp.Rational(1, 3)):
            point_t = toroidal.at(sigma, 1.3, Lambda=lam, B=1)
            point_c = cylindrical.at(sigma, 1.3)
            for a, b in ((0, 0), (3, 3), (2, 2)):
                self.assertLess(abs(complex(toroidal.value(toroidal.g[a, a], point_t)
                                            - cylindrical.value(cylindrical.g[b, b], point_c))), 1e-25)

    def test_levi_civita_as_lambda_tends_to_zero(self):
        for chart_id, sign in (("toroidal", 1), ("cylindrical", -1)):
            chart = Chart(chart_id)
            sigma = sp.Rational(1, 5)
            p = kasner(sigma)
            point = chart.at(sigma, 0.7, Lambda=sign * sp.Rational(1, 10 ** 14), C=1, **({"B": 1} if sign > 0 else {}))
            rho = 0.7
            expected = {0: -rho ** (2 * p["p_0"]), 2: rho ** (2 * p["p_3"]), 3: rho ** (2 * p["p_2"])}
            for a, value in expected.items():
                ratio = complex(chart.value(chart.g[a, a], point)) / complex(sp.N(value, 40))
                self.assertLess(abs(ratio - 1), 1e-10, (chart_id, a))


class Curvature(unittest.TestCase):
    def electric(self, chart, sigma, radius):
        """The eigenvalues of the electric part of the Weyl tensor in the static frame, C_{t i t i}
        over g_tt g_ii, from the published Weyl tensor: the metric is diagonal, so the frame is."""
        weyl = chart.tensor("weyl_tensor", "llll", 4)
        point = chart.at(sigma, radius)
        return [complex(chart.value(weyl.get((0, i, 0, i), 0) / (chart.g[0, 0] * chart.g[i, i]), point)).real
                for i in (1, 2, 3)]

    def test_five_members_are_of_type_d_and_the_rest_of_type_i(self):
        """A static metric's Weyl tensor is electric; two equal eigenvalues make it type D, three
        different ones type I."""
        for chart_id in ("toroidal", "cylindrical"):
            chart = Chart(chart_id)
            for radius in POINTS[chart_id][1][:2]:
                for sigma in TYPE_D:
                    e = self.electric(chart, sigma, radius)
                    gaps = [abs(e[a] - e[b]) for a, b in ((0, 1), (0, 2), (1, 2))]
                    self.assertLess(min(gaps), 1e-20 * max(map(abs, e)), (chart_id, sigma, e))
                    self.assertGreater(max(map(abs, e)), 1e-6, (chart_id, sigma))
                for sigma in TYPE_I:
                    e = self.electric(chart, sigma, radius)
                    gaps = [abs(e[a] - e[b]) for a, b in ((0, 1), (0, 2), (1, 2))]
                    self.assertGreater(min(gaps), 1e-6 * max(map(abs, e)), (chart_id, sigma, e))

    def test_the_kretschmann_scalar_never_takes_the_maximally_symmetric_value(self):
        """K - 8 Lambda^2/3 is the square of the Weyl tensor, so a member that were de Sitter or anti-de
        Sitter space would read 8 Lambda^2/3; every member reads more, at every radius sampled."""
        for chart_id in CHARTS:
            chart = Chart(chart_id)
            K = chart.scalar("kretschmann")
            low, high = POINTS[chart_id][1][0], POINTS[chart_id][1][-1]
            radii = [low + (high - low) * k / 12 for k in range(13)]
            for sigma in TYPE_D + TYPE_I:
                for radius in radii:
                    point = chart.at(sigma, radius)
                    excess = complex(chart.value(K - 8 * chart.lam(point) ** 2 / 3, point)).real
                    self.assertGreater(excess, 0, (chart_id, sigma, radius))

    def test_the_axis_and_the_far_end_are_regular_only_where_the_history_says(self):
        """Toward rho = 0 the Kretschmann scalar grows without bound unless sigma is 0 or 1/2, and toward
        rho = pi/sqrt(3 Lambda) unless sigma is -1/2, 1/4 or 1."""
        chart = Chart("toroidal")
        K = chart.scalar("kretschmann")
        for sigma, axis, far in ((sp.Integer(0), True, False), (sp.Rational(1, 2), True, False),
                                 (sp.Rational(-1, 2), False, True), (sp.Rational(1, 4), False, True),
                                 (sp.Integer(1), False, True), (sp.Rational(1, 5), False, False)):
            near = [complex(chart.value(K, chart.at(sigma, e))).real for e in (1e-3, 1e-4)]
            end = [complex(chart.value(K, chart.at(sigma, float(sp.pi) - e))).real for e in (1e-3, 1e-4)]
            self.assertEqual(abs(near[1]) < 2 * abs(near[0]) + 1, axis, (sigma, near))
            self.assertEqual(abs(end[1]) < 2 * abs(end[0]) + 1, far, (sigma, end))


if __name__ == "__main__":
    unittest.main()
