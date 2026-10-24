"""The perturbed Friedmann universe in Bardeen's potentials, kept to first order in them:
.venv.noindex/bin/python -m unittest discover -s _tools

What its texts state, held on the published files: with the potentials set to zero each
longitudinal chart is the published conformal chart of FRW; the published synchronous chart is the
published flat longitudinal chart pulled back along eta -> eta + d_eta E and x^i -> x^i + d_i E,
with Bardeen's potentials the ones its convention states; the space components off the diagonal
of the Einstein tensor are d_i d_j (Psi - Phi); and the nine tenths the History gives. It also
holds the rule by which the checker counts a potential of a system whose order names no slow time. The
tests need sympy, and are skipped where it is absent.
"""
import importlib.util
import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
METRICS = ROOT / "MFS" / "assets" / "data" / "metrics"
HAS_SYMPY = importlib.util.find_spec("sympy") is not None
sys.path.insert(0, str(Path(__file__).resolve().parent / "derivations"))


def flat(tensor):
    if isinstance(tensor, list):
        for item in tensor:
            yield from flat(item)
    else:
        yield tensor


def chart(metric_id, system):
    data = json.loads((METRICS / f"{metric_id}.json").read_text(encoding="utf-8"))
    return next(c for c in data["coordinates"] if c["id"] == system)


@unittest.skipUnless(HAS_SYMPY, "needs sympy")
class Orders(unittest.TestCase):
    """The first order the checker keeps where the small quantities are potentials of all four
    coordinates and the line's third item, None, says no time is slow."""

    @classmethod
    def setUpClass(cls):
        import sympy
        import verify_metrics
        cls.sp, cls.vm = sympy, verify_metrics

    def test_every_derivative_of_a_potential_counts_as_the_potential(self):
        sp, vm = self.sp, self.vm
        reader = vm.Reader(["\\eta", "x"], ["a = a(\\eta)", "\\Phi = \\Phi(\\eta,x)"], (), kept=({"Phi": 1}, 1, None))
        eta, x = reader.symbol["\\eta"], reader.symbol["x"]
        a, Phi = reader.parameters["a"], reader.parameters["Phi"]
        D = sp.Derivative
        value = a ** 2 + Phi * a + D(Phi, eta) + D(Phi, (eta, 3)) + D(Phi, x) * D(Phi, eta) + Phi ** 2
        kept = a ** 2 + Phi * a + D(Phi, eta) + D(Phi, (eta, 3))
        self.assertEqual(sp.expand(reader.through(value, 1) - kept), 0)
        self.assertEqual(sp.expand(reader.through(value, 0) - a ** 2), 0)
        self.assertIsNone(reader.slow)
        # A line of two items still counts constants alone and refuses a function.
        with self.assertRaises(vm.LatexError):
            vm.Reader(["\\eta", "x"], ["\\Phi = \\Phi(\\eta,x)"], (), kept=({"Phi": 1}, 1))

    def test_terms_of_the_second_order_in_the_metric_move_nothing_that_is_kept(self):
        """With an arbitrary function of the four coordinates put in each of the ten components at
        second order, every tensor the checker hands over for the flat longitudinal chart is what
        it was."""
        import print_charts
        vm = self.vm
        spec = print_charts.bardeen_perturbations("longitudinal")
        coords = spec["system"]["coords"]
        weights = dict(vm.ORDERS[("bardeen_perturbations", "longitudinal")][0])

        def geometry(open_terms):
            names = [f"F_{i}{j} = F_{i}{j}(\\eta,x,y,z)" for i in range(4) for j in range(i, 4)] if open_terms else []
            extra = {f"F_{i}{j}": 2 for i in range(4) for j in range(i, 4)} if open_terms else {}
            reader = vm.Reader(coords, spec["system"]["parameters"] + names, (), kept=({**weights, **extra}, 1))
            g = vm.metric_from_line_element(reader, spec["chart_line_element"], coords)
            symbols = [reader.symbol[name] for name in coords]
            for i in range(4):
                for j in range(i, 4):
                    if open_terms:
                        g[i, j] += reader.parameters[f"F_{i}{j}"]
                        g[j, i] = g[i, j]
            return vm.geometry_of(g, symbols, 10 ** 6, reader)

        fixed, opened = geometry(False), geometry(True)
        for name in ("christoffel_ull", "riemann_llll", "ricci_ll", "einstein_ll", "weyl_llll", "ricci_scalar"):
            ours, theirs = getattr(fixed, name)(), getattr(opened, name)()
            self.assertTrue(all(vm.norm(a - b) == 0 for a, b in zip(flat(ours), flat(theirs))), name)

    def test_a_post_newtonian_system_still_counts_its_slow_time(self):
        sp, vm = self.sp, self.vm
        reader = vm.Reader(["t", "x"], ["U = U(t,x)"], (), kept=({"U": 2}, 2, "t"))
        t, U = reader.symbol["t"], reader.parameters["U"]
        self.assertEqual(sp.expand(reader.through(U + sp.Derivative(U, t), 2) - U), 0)


@unittest.skipUnless(HAS_SYMPY, "needs sympy")
class Published(unittest.TestCase):
    """The published charts read through the checker's own Reader."""

    @classmethod
    def setUpClass(cls):
        import sympy
        import verify_metrics
        cls.sp, cls.vm = sympy, verify_metrics

    def reader(self, system):
        vm = self.vm
        key = ("bardeen_perturbations", system)
        c = chart(*key)
        reader = vm.Reader(c["coords"], [p["symbol"] for p in c["parameters"]], (),
                           kept=vm.ORDERS[key], held=vm.HELD.get(key, ()), rates=vm.RATES.get(key))
        return c, reader

    def scalar(self, reader, text):
        value = reader(text.partition("=")[2])
        return value.subs(reader.held).doit() if reader.held else value

    def test_with_no_potentials_each_longitudinal_chart_is_friedmanns(self):
        """R and K at Phi = Psi = 0 against the published conformal chart of FRW, at k = 0 for the
        flat chart and at the same k for the curved one."""
        sp, vm = self.sp, self.vm
        frw = chart("frw", "conformal_spherical")
        theirs = vm.Reader(frw["coords"], [p["symbol"] for p in frw["parameters"]], ())
        for system in ("longitudinal", "longitudinal_curved"):
            c, reader = self.reader(system)
            zero = {reader.parameters["Phi"]: 0, reader.parameters["Psi"]: 0}
            k = reader.parameters["k"] if system == "longitudinal_curved" else 0
            for field in ("ricci_scalar", "kretschmann"):
                ours = self.scalar(reader, c[field]).subs(zero).doit()
                friedmann = theirs(frw[field].partition("=")[2])
                friedmann = friedmann.subs(theirs.parameters["k"], k).subs(
                    theirs.parameters["a"], reader.parameters["a"]).subs(theirs.symbol["\\eta"], reader.symbol["\\eta"]).doit()
                self.assertEqual(vm.norm(ours - friedmann), 0, f"{system} {field}")

    def test_the_synchronous_chart_is_the_longitudinal_chart_pulled_back(self):
        sp, vm = self.sp, self.vm
        c, reader = self.reader("synchronous")
        long_c, long_reader = self.reader("longitudinal")
        symbols = [reader.symbol[name] for name in c["coords"]]
        eta = symbols[0]
        a, psi, E = (reader.parameters[k] for k in ("a", "psi", "E"))
        H = sp.diff(a, eta) / a
        shift = sp.diff(E, eta)
        potentials = {long_reader.parameters["Phi"]: -sp.diff(E, eta, 2) - H * shift,
                      long_reader.parameters["Psi"]: psi + H * shift}
        theirs = sp.zeros(4, 4)
        for entry in long_c["metric_components"]:
            i, j = (long_c["coords"].index(name) for name in entry["indices"])
            value = long_reader(entry["value"]).subs(dict(zip([long_reader.symbol[n] for n in long_c["coords"]], symbols)))
            theirs[i, j] = value.subs(potentials).subs(a, a + sp.diff(a, eta) * shift).doit()
        moved = [eta + shift] + [s + sp.diff(E, s) for s in symbols[1:]]
        jacobian = sp.Matrix(4, 4, lambda i, j: sp.diff(moved[i], symbols[j]))
        pulled = jacobian.T * theirs * jacobian
        ours = sp.zeros(4, 4)
        for entry in c["metric_components"]:
            i, j = (c["coords"].index(name) for name in entry["indices"])
            ours[i, j] = reader(entry["value"])
        for i in range(4):
            for j in range(4):
                self.assertEqual(reader.through(pulled[i, j] - ours[i, j], 1), 0, f"slot {i}{j}")

    def test_the_space_components_off_the_diagonal_are_the_difference_of_the_potentials(self):
        sp, vm = self.sp, self.vm
        c, reader = self.reader("longitudinal")
        Phi, Psi = reader.parameters["Phi"], reader.parameters["Psi"]
        x, y = reader.symbol["x"], reader.symbol["y"]
        entry = next(e for e in c["einstein_tensor"]["variants"]["ll"]["nonzero"] if e["indices"] == ["x", "y"])
        self.assertEqual(vm.norm(reader(entry["value"]) - sp.diff(Psi - Phi, x, y)), 0)


@unittest.skipUnless(HAS_SYMPY, "needs sympy")
class History(unittest.TestCase):
    def test_a_constant_potential_falls_to_nine_tenths_from_radiation_to_matter(self):
        """Mukhanov, Feldman and Brandenberger's (6.58) with the potential constant in each era:
        zeta = (2/3) Phi/(1 + w) + Phi is 3 Phi/2 at w = 1/3 and 5 Phi/3 at w = 0."""
        import sympy as sp
        Phi_r, Phi_m = sp.symbols("Phi_r Phi_m")

        def zeta(Phi, w):
            return sp.Rational(2, 3) * Phi / (1 + w) + Phi

        solved = sp.solve(zeta(Phi_r, sp.Rational(1, 3)) - zeta(Phi_m, 0), Phi_m)[0]
        self.assertEqual(solved, sp.Rational(9, 10) * Phi_r)
        data = json.loads((METRICS / "bardeen_perturbations.json").read_text(encoding="utf-8"))
        self.assertIn("nine tenths", data["history"])


if __name__ == "__main__":
    unittest.main()
