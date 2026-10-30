"""Penney's charged scalar solution, held on the published file with the checker's own reader:
.venv.noindex/bin/python -m unittest _tools/test_penney_scalar.py

The published Ricci tensor of the spherical chart is that of the scalar field
psi = sqrt(1 - gamma^2) ln((r - a)/(r - b))/2 and the Coulomb field of the charge q = gamma sqrt(ab);
at gamma = 1 the chart is Reissner and Nordstrom's published chart with r_s = a + b and r_q^2 = ab,
and at b = 0 Fisher's published spherical chart; far out -g_tt = 1 - gamma(a + b)/r; the spheres
shrink to a point at r = a and g^rr stays positive outside it, so no horizon stands there; and the
harmonic chart is the spherical one pulled back by e^(-2ku) = (r - a)/(r - b), with a - b = 2k,
gamma = h/k and a/b = e^(2 h u_1), its H at spatial infinity the mass. Every comparison is made at
points, to thirty digits. Needs sympy, and is skipped where it is absent.
"""
import importlib.util
import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "MFS" / "assets" / "data"
HAS_SYMPY = importlib.util.find_spec("sympy") is not None
sys.path.insert(0, str(Path(__file__).resolve().parent / "derivations"))

# Points of (a, b, gamma, r), each outside r = a.
POINTS = [(1, "1/4", "1/2", "2"), (3, 1, "1/5", "7/2"), ("5/2", 2, "9/10", "11/4")]


def chart(metric_id, system):
    import verify_metrics
    entry = next(c for c in json.loads((DATA / "metrics" / f"{metric_id}.json").read_text())["coordinates"]
                 if c["id"] == system)
    reader = verify_metrics.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]], ())

    def read(text):
        value = reader(text)
        return value.subs(reader.held).subs(reader.held).doit() if reader.held else value
    return entry, reader, read


def matrix(entry, read, field, sp):
    out = sp.zeros(4, 4)
    for e in entry[field]:
        i, j = (entry["coords"].index(c) for c in e["indices"])
        out[i, j] = read(e["value"])
    return out


@unittest.skipUnless(HAS_SYMPY, "sympy is not installed")
class ThePublishedCharts(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import sympy
        cls.sp = sympy
        entry, reader, read = chart("penney_scalar", "spherical")
        cls.entry, cls.reader, cls.read = entry, reader, staticmethod(read)
        cls.X = [reader.symbol[c] for c in entry["coords"]]
        cls.g = matrix(entry, read, "metric_components", sympy)
        cls.ginv = matrix(entry, read, "inverse_metric_components", sympy)
        cls.a, cls.b, cls.gamma = (reader.parameters[n] for n in ("a", "b", "gamma"))

    def at(self, point):
        sp = self.sp
        a, b, gamma, r = (sp.Rational(v) for v in point)
        return {self.a: a, self.b: b, self.gamma: gamma, self.X[1]: r, self.X[2]: sp.Rational(7, 10), self.reader.c: 1}

    def test_the_published_ricci_tensor_is_that_of_the_scalar_and_coulomb_fields(self):
        sp = self.sp
        _, r, theta, _ = self.X
        a, b, gamma = self.a, self.b, self.gamma
        psi = sp.sqrt(1 - gamma ** 2) * sp.log((r - a) / (r - b)) / 2
        F = sp.zeros(4, 4)
        F[0, 1] = gamma * sp.sqrt(a * b) * (-self.g[0, 0]) * self.g[1, 1] / self.g[2, 2]
        F[1, 0] = -F[0, 1]
        ricci = sp.zeros(4, 4)
        for e in self.entry["ricci_tensor"]["variants"]["ll"]["nonzero"]:
            i, j = (self.entry["coords"].index(c) for c in e["indices"])
            ricci[i, j] = self.read(e["value"])
        for point in POINTS:
            at = self.at(point)
            gi = self.ginv.subs(at).evalf(40)
            g = self.g.subs(at).evalf(40)
            Fn = F.subs(at).evalf(40)
            square = sum(Fn[i, j] * Fn[k, l] * gi[i, k] * gi[j, l]
                         for i in range(4) for j in range(4) for k in range(4) for l in range(4))
            for i in range(4):
                for j in range(4):
                    source = (2 * (sp.diff(psi, self.X[i]) * sp.diff(psi, self.X[j])).subs(at)
                              + 2 * sum(Fn[i, k] * Fn[j, l] * gi[k, l] for k in range(4) for l in range(4))
                              - g[i, j] * square / 2)
                    self.assertLess(abs(sp.N(ricci[i, j].subs(at) - source, 30)), 1e-25, (point, i, j))

    def test_gamma_one_is_reissner_and_nordstrom_and_b_zero_is_fisher(self):
        sp = self.sp
        for other, system, names in (("rn_metric", "spherical", ("r_s", "r_q")), ("fisher_jnw", "spherical", ("b", "gamma"))):
            entry, reader, read = chart(other, system)
            theirs = matrix(entry, read, "metric_components", sp)
            Y = [reader.symbol[c] for c in entry["coords"]]
            for point in POINTS:
                at = self.at(point)
                if other == "rn_metric":
                    at[self.gamma] = 1
                    same = {reader.parameters["r_s"]: at[self.a] + at[self.b], reader.parameters["r_q"]: sp.sqrt(at[self.a] * at[self.b])}
                else:
                    at[self.b] = 0
                    same = {reader.parameters["b"]: at[self.a], reader.parameters["gamma"]: at[self.gamma]}
                same.update({Y[1]: at[self.X[1]], Y[2]: at[self.X[2]], reader.c: 1})
                for i in range(4):
                    self.assertLess(abs(sp.N(self.g[i, i].subs(at) - theirs[i, i].subs(same), 30)), 1e-25, (other, point, i))

    def test_the_mass_is_gamma_a_plus_b_over_two(self):
        sp = self.sp
        for point in POINTS:
            at = self.at(point)
            at[self.X[1]] = sp.Integer(10) ** 15
            mass = sp.N(((1 + self.g[0, 0]) * self.X[1]).subs(at), 30)
            self.assertLess(abs(mass - at[self.gamma] * (at[self.a] + at[self.b])), 1e-10, point)

    def test_the_spheres_shrink_to_a_point_at_r_equal_a_with_no_horizon_outside_it(self):
        sp = self.sp
        for point in POINTS:
            at = self.at(point)
            if at[self.gamma] == 1:
                continue
            areas = []
            for k in (20, 30, 40):
                at[self.X[1]] = at[self.a] + sp.Rational(1, 10 ** k)
                areas.append(sp.N(self.g[2, 2].subs(at), 30))
            # g_thetatheta goes as (r - a)^(1 - gamma), falling by 10^(10 (1 - gamma)) at each step.
            fall = sp.Integer(10) ** (-10 * (1 - at[self.gamma]))
            for near, nearer in zip(areas, areas[1:]):
                self.assertLess(abs(nearer / near / fall - 1), 1e-3, (point, areas))
            for step in (sp.Rational(1, 10 ** 6), sp.Rational(1, 100), 1, 100):
                at[self.X[1]] = at[self.a] + step
                self.assertGreater(sp.N(self.ginv[1, 1].subs(at), 30), 0, (point, step))

    def test_the_harmonic_chart_is_the_spherical_chart_pulled_back(self):
        sp = self.sp
        entry, reader, read = chart("penney_scalar", "harmonic")
        theirs = matrix(entry, read, "metric_components", sp)
        Y = [reader.symbol[c] for c in entry["coords"]]
        k, h, u1 = (reader.parameters[n] for n in ("k", "h", "u_1"))
        u = Y[1]
        for kv, hv, u1v, uv in ((1, sp.Rational(1, 2), sp.log(4), sp.Rational(1, 3)),
                                (sp.Rational(3, 2), sp.Rational(1, 5), 2, sp.Rational(5, 2))):
            bv = 2 * kv / (sp.exp(2 * hv * u1v) - 1)
            av = bv + 2 * kv
            radius = bv + 2 * kv / (1 - sp.exp(-2 * kv * u))
            here = {k: kv, h: hv, u1: u1v, Y[2]: sp.Rational(7, 10), reader.c: 1}
            at = {self.a: av, self.b: bv, self.gamma: hv / kv, self.X[2]: sp.Rational(7, 10), self.reader.c: 1}
            for i in range(4):
                mine = self.g[i, i].subs(at).subs(self.X[1], radius) * (sp.diff(radius, u) ** 2 if i == 1 else 1)
                self.assertLess(abs(sp.N((mine - theirs[i, i].subs(here)).subs(u, uv), 30)), 1e-25, (kv, hv, i))
            # H at u = 0 is the mass, gamma (a + b)/2.
            H = reader.parameters["H"]
            self.assertLess(abs(sp.N(H.subs(here).subs(u, 0) - (hv / kv) * (av + bv) / 2, 30)), 1e-25)


if __name__ == "__main__":
    unittest.main()
