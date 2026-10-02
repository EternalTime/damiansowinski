"""Tests for _tools/derivations/print_charts.py: python3 -m unittest discover -s _tools

A chart the script wrote is printed again and held to the published file byte for byte, so a
change to the printer that reorders a published value is caught here and not by a reader. The
tests need sympy, which the checker's own environment has, and are skipped where it is absent.
"""
import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path

DERIVATIONS = Path(__file__).resolve().parent / "derivations"
HAS_SYMPY = importlib.util.find_spec("sympy") is not None


@unittest.skipUnless(HAS_SYMPY, "print_charts.py needs sympy")
class Reprint(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        sys.path.insert(0, str(DERIVATIONS))
        import print_charts
        cls.pc = print_charts

    @classmethod
    def tearDownClass(cls):
        sys.path.remove(str(DERIVATIONS))

    def reprinted(self, metric_id):
        """The file print_charts.py writes for one spacetime, printed into a copy of its own."""
        pc = self.pc
        published = pc.METRICS / f"{metric_id}.json"
        with tempfile.TemporaryDirectory() as folder:
            copy = Path(folder) / published.name
            copy.write_bytes(published.read_bytes())
            kept, pc.METRICS = pc.METRICS, Path(folder)
            try:
                builders = pc.CHARTS[metric_id]
                for build in builders if isinstance(builders, list) else [builders]:
                    specs = build()
                    for spec in specs if isinstance(specs, list) else [specs]:
                        pc.write(spec)
            finally:
                pc.METRICS = kept
            return published.read_bytes(), copy.read_bytes()

    def test_tov_prints_again_as_published(self):
        # Its collected bracket closes each product, r(r - 2m)((d_r Phi)^2 + d_r^2 Phi), which a
        # change to where the printer writes a named factor moved to the front on 30 September 2026.
        published, printed = self.reprinted("tov")
        self.assertEqual(printed, published)

    def test_a_chart_prints_the_same_twice(self):
        # The random point block() tells values apart at is seeded, so two prints agree.
        first = self.reprinted("einstein_static")[1]
        self.assertEqual(self.reprinted("einstein_static")[1], first)


@unittest.skipUnless(HAS_SYMPY, "verify_metrics.py needs sympy")
class Orders(unittest.TestCase):
    """A system kept to an order, as Hartle and Thorne's exterior is to the second order of its spin:
    every tensor is a Taylor polynomial to that order, and a value that carries more is refused."""

    @classmethod
    def setUpClass(cls):
        sys.path.insert(0, str(DERIVATIONS))
        import sympy
        import verify_metrics
        cls.sp, cls.vm = sympy, verify_metrics

    @classmethod
    def tearDownClass(cls):
        sys.path.remove(str(DERIVATIONS))

    def reader(self, weights, highest):
        return self.vm.Reader(["t", "r", "\\theta", "\\phi"], ["m", "a", "q"], {"t"}, kept=(weights, highest))

    def test_a_value_is_cut_at_the_order_with_each_parameter_counted_by_its_weight(self):
        reader = self.reader({"a": 1, "q": 2}, 2)
        r, (m, a, q) = reader.symbol["r"], (reader.parameters[k] for k in "maq")
        value = 1 / (1 - 2 * m / r + a ** 2 / r ** 2 + q * a / r ** 3) + a ** 3 + q ** 2
        cut = 1 / (1 - 2 * m / r) - a ** 2 / (r ** 2 * (1 - 2 * m / r) ** 2)
        self.assertEqual(self.vm.norm(reader.truncated(value) - cut), 0)
        self.assertTrue(self.vm.beyond_order(reader, value))
        self.assertFalse(self.vm.beyond_order(reader, cut))

    def test_a_system_with_no_order_is_left_whole(self):
        reader = self.vm.Reader(["t", "r"], ["m", "a"], {"t"})
        value = reader.parameters["a"] ** 5 / reader.symbol["r"]
        self.assertEqual(reader.truncated(value), value)
        self.assertFalse(self.vm.beyond_order(reader, value))

    def test_an_order_is_refused_for_a_name_that_is_no_constant_of_the_system(self):
        with self.assertRaises(self.vm.LatexError):
            self.vm.Reader(["t", "r"], ["m"], {"t"}, kept=({"a": 1}, 1))

    def test_the_inverse_and_the_curvature_of_the_weak_field_are_those_of_linearised_gravity(self):
        sp, vm = self.sp, self.vm
        reader = self.reader({"m": 1}, 1)
        t, r, th, ph = (reader.symbol[k] for k in ("t", "r", "\\theta", "\\phi"))
        m, a = reader.parameters["m"], reader.parameters["a"]
        g = sp.diag(-(1 - 2 * m / r), 1 + 2 * m / r, r ** 2 * (1 + 2 * m / r), r ** 2 * sp.sin(th) ** 2 * (1 + 2 * m / r))
        g[0, 3] = g[3, 0] = -2 * a * m * sp.sin(th) ** 2 / r
        geometry = vm.Geometry(g, [t, r, th, ph], 120, reader)
        self.assertEqual(vm.norm(geometry.ginv[0, 0] + 1 + 2 * m / r), 0)
        self.assertEqual(vm.norm(geometry.ginv[0, 3] + 2 * a * m / r ** 3), 0)
        self.assertEqual(reader.truncated(g * geometry.ginv), sp.eye(4))
        self.assertTrue(all(value == 0 for row in geometry.ricci_ll() for value in row))
        self.assertEqual(vm.norm(geometry.riemann_llll()[0][1][0][1] + 2 * m / r ** 3), 0)
        # The Kretschmann scalar is of second order in the mass, so to first order it vanishes.
        self.assertEqual(geometry.kretschmann(), 0)

    def test_hartle_and_thornes_exterior_is_kerr_to_second_order_at_kerrs_quadrupole_moment(self):
        sys.path.insert(0, str(DERIVATIONS))
        import chart_printer
        import print_charts
        spec = print_charts.hartle_thorne("hartle_thorne")
        chart = chart_printer.Chart(spec["system"]["coords"], spec["system"]["parameters"], spec["chart_line_element"],
                                    spec["printer"], order=self.vm.ORDERS[("hartle_thorne", "hartle_thorne")])
        print_charts.hartle_thorne_kerr(chart)


@unittest.skipUnless(HAS_SYMPY and importlib.util.find_spec("scipy") is not None, "null_rays.py needs sympy and scipy")
class PointParticles(unittest.TestCase):
    """Point particles in three dimensions, held to numbers the drawings never use."""

    @classmethod
    def setUpClass(cls):
        sys.path.insert(0, str(DERIVATIONS))
        import null_rays
        import print_charts
        cls.nr, cls.pc = null_rays, print_charts

    @classmethod
    def tearDownClass(cls):
        sys.path.remove(str(DERIVATIONS))

    def test_two_particles_of_alpha_three_quarters_stand_a_beta_function_apart(self):
        """The proper distance between two particles at x = +-d is d times the integral of
        (1 - s^2)^(-1/4) from -1 to 1, which is B(1/2, 3/4), and beyond a particle the distance
        grows as the square root far away, 2 sqrt(x), as on a cone of alpha = 1/2."""
        import math
        beta = math.gamma(0.5) * math.gamma(0.75) / math.gamma(1.25)
        distance = self.nr._two_particles_distance
        self.assertAlmostEqual(2 * float(distance(1.0)), beta, places=12)
        self.assertAlmostEqual(float(distance(-1.0)), -float(distance(1.0)), places=14)
        far = float(distance(4.0e6)) - float(distance(1.0e6))
        self.assertAlmostEqual(far / (2 * math.sqrt(4.0e6) - 2 * math.sqrt(1.0e6)), 1.0, places=6)

    def test_every_chart_passes_its_own_check(self):
        """Each flat chart is the conical chart pulled back, the two bodies are one particle when
        the other has no mass, and the planet is dust at rest that meets the cone."""
        import chart_printer as cp
        for system in self.pc.PARTICLE_CHARTS:
            spec = self.pc.point_particle_2plus1(system)
            chart = cp.Chart(spec["system"]["coords"], spec["system"]["parameters"], spec["chart_line_element"])
            self.pc.point_particle_check(chart, system)

    def test_a_misstated_deficit_is_refused(self):
        """The wedge chart against a cone of another alpha: the pullback misses."""
        import chart_printer as cp
        import sympy as sp
        spec = self.pc.point_particle_2plus1("wedge")
        chart = cp.Chart(spec["system"]["coords"], spec["system"]["parameters"],
                         spec["chart_line_element"].replace("r^2d\\theta^2", "4r^2d\\theta^2"))
        with self.assertRaises(AssertionError):
            self.pc.point_particle_check(chart, "wedge")


if __name__ == "__main__":
    unittest.main()
