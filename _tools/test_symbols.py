"""Every symbol a spacetime uses is defined where a reader can find it.

The captain found Ori's time machine drawn in units of a length $\\ell$ that nothing on the page
defined, beside the $L$ its metric uses for the period of $z$ (1 October 2026). `symbols.py`
holds the rule that came of it, and these tests hold the rule to the files on disk and to what
it is meant to catch and to pass.
"""
import copy
import unittest

import build_mfs_data as build
import symbols


def metric(**chart):
    """A spacetime of one chart, Schwarzschild's, with any of the chart's fields replaced."""
    system = {
        "id": "spherical",
        "coords": ["t", "r", "\\theta", "\\phi"],
        "parameters": [{"symbol": "r_s", "description": "Schwarzschild radius, $r_s = 2GM/c^2$ for the mass $M$"}],
        "convention": "We use coordinates $(t, r, \\theta, \\phi)$, with $t$ carrying dimensions of time.",
        "domains": ["r \\in (r_s, \\infty)"],
        "line_element": "ds^2 = -\\left(1 - \\dfrac{r_s}{r}\\right)c^2 dt^2 + \\dfrac{dr^2}{1 - \\dfrac{r_s}{r}} "
                        "+ r^2 d\\Omega^2",
        "metric_components": [{"indices": ["t", "t"], "value": "-\\left(1 - \\dfrac{r_s}{r}\\right)"}],
        "christoffel": {"variants": {"ull": {"nonzero": [{"indices": ["r", "t", "t"], "value": "\\dfrac{r_s}{2r^2}"}]}}},
        "kretschmann": "K = \\dfrac{12 r_s^2}{r^6}",
        "geodesics": ["\\ddot{t} + \\dfrac{r_s}{r(r - r_s)}\\dot{t}\\dot{r} = 0"],
    }
    system.update(chart)
    return {"id": "x", "convention": "We keep factors of $c$ explicit.", "coordinates": [system]}


def diagram(**view):
    """A spacetime diagram file of one view of that chart, with any of the view's fields replaced."""
    drawn = {"id": "radial", "label": "$t$ and $r$", "xlabel": "$r/r_s$", "ylabel": "$ct/r_s$", "families": [],
             "markers": [], "settings": "$\\theta = \\frac{\\pi}{2}$, $\\phi = 0$",
             "caption": ["The plane of $t$ and $r$, where the cones close at $r = r_s$."]}
    drawn.update(view)
    return {"metric": "x", "systems": {"spherical": [drawn]}}


class Reading(unittest.TestCase):
    """What the check takes a symbol to be."""

    def test_a_symbol_is_a_letter_with_its_subscript_and_its_accent(self):
        self.assertEqual(symbols.symbols("r + r_s + \\tilde\\phi + r_* + \\ell"),
                         {"r", "r_s", "\\tilde{\\phi}", "r_*", "\\ell"})

    def test_a_dot_and_a_power_make_no_new_symbol(self):
        self.assertEqual(symbols.symbols("\\ddot{r} + \\dot{t}^2 + e^{-z/2}"), {"r", "t", "e", "z"})

    def test_words_units_and_functions_are_not_symbols(self):
        self.assertEqual(symbols.symbols("\\sin\\theta\\,\\mathrm{sgn}(x) \\;\\text{(the horizon)} + r_{\\rm EB}"),
                         {"\\theta", "x", "r"})

    def test_a_tensor_that_carries_indices_is_left_out(self):
        self.assertEqual(symbols.symbols("R^\\mu{}_{\\nu\\rho\\sigma} + g_{tt} + \\Gamma^\\mu{}_{tt} + K_{ij}K^{ij}",
                                         ["t", "r"]), set())

    def test_a_bare_tensor_letter_is_an_ordinary_symbol(self):
        self.assertEqual(symbols.symbols("R = 2T", ["t", "r"]), {"R", "T"})

    def test_a_component_answers_to_itself_and_to_its_letter(self):
        self.assertEqual(symbols.symbols("p_\\phi", ["t", "\\phi"]), {("p_\\phi", "p")})
        self.assertEqual(symbols.undefined("p_\\phi", {"p"}, ["t", "\\phi"]), set())
        self.assertEqual(symbols.undefined("p_\\phi", set(), ["t", "\\phi"]), {"p"})

    def test_a_subscript_that_is_no_coordinate_makes_a_symbol_of_its_own(self):
        self.assertEqual(symbols.undefined("r_h", {"r"}, ["t", "r"]), {"r_h"})

    def test_both_signs_are_defined_with_the_pair(self):
        self.assertEqual(symbols.names({"r_\\pm"}), {"r_\\pm", "r_+", "r_-"})

    def test_the_solid_angle_and_a_change_are_universal_only_as_written(self):
        self.assertEqual(symbols.symbols("r^2 d\\Omega^2 + \\Delta v"), {"r", "d", "v"})
        self.assertEqual(symbols.symbols("\\Omega^2 + \\Delta/\\Sigma"), {"\\Omega", "\\Delta", "\\Sigma"})

    def test_a_sphere_is_universal_and_a_length_named_like_it_is_not(self):
        self.assertEqual(symbols.undefined("S^1 \\times S^2", set(), []), set())
        self.assertEqual(symbols.undefined("g S \\ll c^2", set(), []), {"S"})


class Defining(unittest.TestCase):
    """The three ways a sentence defines a symbol, and what does not."""

    def defined(self, text, **options):
        return symbols.defined_in(text, ["t", "r"], **options)

    def test_an_equation_defines_the_symbol_alone_on_its_left(self):
        self.assertEqual(self.defined("with $r_* = r + r_s\\ln(r/r_s - 1)$"), {"r_*"})
        self.assertEqual(self.defined("so $p, q = \\arctan((ct \\mp r)/\\ell)$ bring it in"), {"p", "q"})
        self.assertEqual(self.defined("with $f(r) = 1 - r_s/r$"), {"f"})
        self.assertEqual(self.defined("the horizons $r_\\pm = m \\pm \\sqrt{m^2 - a^2}$"), {"r_\\pm", "r_+", "r_-"})

    def test_its_power_its_differential_and_its_derivative_are_the_symbol_still(self):
        self.assertEqual(self.defined("with $\\rho^2 = x^2 + y^2$"), {"\\rho"})
        self.assertEqual(self.defined("with $dr_* = dr/Q$"), {"r_*"})
        self.assertEqual(self.defined("the surface climbs at $dz/dr = \\sqrt{r_s/(r - r_s)}$"), {"z"})
        self.assertEqual(self.defined("and $c\\,d\\tau = dT/\\sqrt{U}$"), {"\\tau"})

    def test_an_expression_on_the_left_defines_nothing(self):
        self.assertEqual(self.defined("with $ex^2 = 2\\,\\ell^2$"), set())
        self.assertEqual(self.defined("where $1 + 2E = 1 - r_s/r$"), set())

    def test_a_value_says_how_much_and_not_what(self):
        self.assertEqual(self.defined("drawn at $\\ell = 1$"), set())
        self.assertEqual(self.defined("$L = 2\\pi$, $x = 4$"), set())

    def test_a_value_defines_where_the_sentence_says_what_was_drawn(self):
        self.assertEqual(self.defined("$w = 0.05$", values=True), {"w"})

    def test_an_equation_to_solve_defines_the_symbol_it_holds(self):
        self.assertEqual(self.defined("with $\\sinh r_c = 1$"), {"r_c"})

    def test_a_naming_noun_defines_the_symbol_after_it(self):
        self.assertEqual(self.defined("the height $z$ of the surface"), {"z"})
        self.assertEqual(self.defined("for integer $k$"), {"k"})
        self.assertEqual(self.defined("the Kruskal coordinates $U$ and $V$"), {"U", "V"})
        self.assertEqual(self.defined("coordinates $Z_0$, $Z_1$, and $Z_2$"), {"Z_0", "Z_1", "Z_2"})

    def test_a_noun_before_an_expression_names_nothing(self):
        self.assertEqual(self.defined("the radius $2GM/c^2$"), set())
        self.assertEqual(self.defined("lines of constant $\\chi$"), set())

    def test_saying_what_it_is_defines_it(self):
        self.assertEqual(self.defined("with $\\ell$ any length, bring it in"), {"\\ell"})
        self.assertEqual(self.defined("$\\bar\\rho$ is the mean density"), {"\\bar{\\rho}"})
        self.assertEqual(self.defined("with $Y_0$ and $Y_1$ the Bessel functions"), {"Y_0", "Y_1"})

    def test_a_domain_names_a_symbol_in_its_words(self):
        domain = "x = y = t = 0 \\;\\text{(the closed null geodesic } N\\text{)}"
        self.assertIn("N", symbols.defined_in(symbols.as_sentence(domain), ["t", "x", "y"]))


class Charts(unittest.TestCase):
    """A chart's own mathematics is held to that chart's coordinates, parameters and conventions."""

    def test_a_chart_that_defines_all_it_uses_passes(self):
        self.assertEqual(symbols.symbol_problems(metric()), [])

    def test_a_symbol_in_the_line_element_that_nothing_defines_is_named_with_its_chart_and_field(self):
        problems = symbols.symbol_problems(metric(line_element="ds^2 = -c^2 dt^2 + \\dfrac{r^2}{b^2}dr^2"))
        self.assertEqual(problems, ["x.json: the chart 'spherical' uses $b$ in line_element, and neither its "
                                    "coordinates, its parameters nor its convention defines it"])

    def test_every_kind_of_mathematics_is_read(self):
        for field, value in (
            ("domains", ["r \\in (r_h, \\infty)"]),
            ("metric_components", [{"indices": ["t", "t"], "value": "-N^2"}]),
            ("christoffel", {"variants": {"ull": {"nonzero": [{"indices": ["r", "t", "t"], "value": "\\kappa"}]}}}),
            ("kretschmann", "K = 12 m^2/r^6"),
            ("geodesics", ["\\ddot{t} + \\omega\\dot{t}\\dot{r} = 0"]),
        ):
            with self.subTest(field):
                self.assertEqual(len(symbols.symbol_problems(metric(**{field: value}))), 1)

    def test_a_parameter_a_convention_or_a_coordinate_defines(self):
        used = dict(line_element="ds^2 = -N^2c^2 dt^2 + dr^2")
        radius = {"symbol": "r_s", "description": "Schwarzschild radius"}
        self.assertEqual(len(symbols.symbol_problems(metric(**used))), 1)
        for definition in (
            dict(parameters=[radius, {"symbol": "N = N(r)", "description": "Lapse"}]),
            dict(parameters=[dict(radius, description="Schwarzschild radius, with $N^2 = 1 - r_s/r$")]),
            dict(convention="We use coordinates $(t, r)$ and write $N$ for the lapse."),
        ):
            with self.subTest(definition):
                self.assertEqual(symbols.symbol_problems(metric(**used, **definition)), [])

    def test_the_shared_convention_defines_for_every_chart(self):
        spacetime = metric(line_element="ds^2 = -N^2c^2 dt^2 + dr^2")
        spacetime["convention"] = "We keep factors of $c$ explicit and write $N$ for the lapse."
        self.assertEqual(symbols.symbol_problems(spacetime), [])

    def test_one_chart_does_not_define_for_another(self):
        spacetime = metric()
        other = copy.deepcopy(spacetime["coordinates"][0])
        other.update(id="tortoise", coords=["t", "r_*", "\\theta", "\\phi"], domains=[], christoffel={}, geodesics=[],
                     kretschmann="K = 0", metric_components=[], line_element="ds^2 = -c^2 dt^2 + dr_*^2",
                     convention="We use the tortoise coordinate $r_*$.")
        spacetime["coordinates"][0]["line_element"] = "ds^2 = -c^2 dt^2 + dr_*^2"
        spacetime["coordinates"].append(other)
        self.assertEqual([problem.split(" uses ")[0] for problem in symbols.symbol_problems(spacetime)],
                         ["x.json: the chart 'spherical'"])

    def test_the_scalar_a_field_names_on_its_left_is_its_own(self):
        self.assertEqual(symbols.symbol_problems(metric(kretschmann="K = 0", ricci_scalar="R = 0")), [])


class Drawings(unittest.TestCase):
    """A drawing is held to the spacetime's charts and to what its own texts define."""

    def problems(self, **view):
        return symbols.symbol_problems(metric(), diagram(**view))

    def test_a_drawing_in_the_charts_symbols_passes(self):
        self.assertEqual(self.problems(), [])

    def test_a_unit_on_an_axis_that_nothing_defines_is_named_with_its_view_and_field(self):
        self.assertEqual(self.problems(xlabel="$r/\\ell$"),
                         ["x.json: diagrams spherical/radial uses $\\ell$ in xlabel, and neither a chart of the "
                          "spacetime nor a definition of its own defines it"])

    def test_every_text_of_a_drawing_is_read(self):
        for field, value in (
            ("label", "$r = 4\\,\\ell$"), ("ylabel", "$ct/\\ell$"), ("caption", ["At $r = \\ell$ the cones close."]),
            ("markers", [{"legend": "the circle $r = \\ell$"}]), ("slices", [{"label": "$ct = 2\\,\\ell$"}]),
            ("legend", [["line", "r", "$r$ constant, in units of $\\ell$"]]), ("labels", [{"text": "$r = \\ell$"}]),
            ("restriction", "The plane $r > \\ell$ only."), ("input", "A shell at $r = 3\\ell$."),
        ):
            with self.subTest(field):
                self.assertEqual(len(self.problems(**{field: value})), 1)

    def test_one_symbol_in_many_fields_of_a_view_is_one_finding(self):
        self.assertEqual(self.problems(xlabel="$r/\\ell$", ylabel="$ct/\\ell$", label="$r = 4\\,\\ell$"),
                         ["x.json: diagrams spherical/radial uses $\\ell$ in label, xlabel, ylabel, and neither a "
                          "chart of the spacetime nor a definition of its own defines it"])

    def test_a_definition_anywhere_in_the_view_serves_the_whole_view(self):
        self.assertEqual(self.problems(xlabel="$r/\\ell$", caption=["The plane of $t$ and $r$, with $\\ell$ any length."]),
                         [])

    def test_a_value_defines_in_the_settings_and_not_in_a_caption(self):
        self.assertEqual(self.problems(xlabel="$r/\\ell$", settings="$\\ell = 1$, $\\phi = 0$"), [])
        self.assertEqual(len(self.problems(xlabel="$r/\\ell$", caption=["Drawn at $\\ell = 1$."])), 1)

    def test_the_convention_defines_for_every_drawing(self):
        spacetime = metric()
        spacetime["convention"] = "We keep factors of $c$ explicit. The drawings take any length $\\ell$ as their unit."
        self.assertEqual(symbols.symbol_problems(spacetime, diagram(xlabel="$r/\\ell$")), [])

    def test_any_chart_of_the_spacetime_defines_for_a_drawing(self):
        spacetime = metric()
        other = copy.deepcopy(spacetime["coordinates"][0])
        other.update(id="ingoing", coords=["v", "r", "\\theta", "\\phi"], geodesics=[],
                     line_element="ds^2 = 2\\,dv\\,dr", convention="We use the advanced time $v$.")
        spacetime["coordinates"].append(other)
        caption = ["Its $g_{tr}$ vanishes where the ingoing chart's $g_{vr}$ is $+1$."]
        self.assertEqual(symbols.symbol_problems(spacetime, diagram(caption=caption)), [])

    def test_a_conformal_diagram_has_its_own_coordinates_and_no_other_drawing_has_them(self):
        view = {"id": "spherical", "label": "Spherical", "caption": ["The singularity lies on the lines $T = \\pm\\pi/2$."],
                "legend": [], "labels": [{"text": "$q = 0$"}]}
        self.assertEqual(symbols.symbol_problems(metric(), None, {"views": [view]}), [])
        self.assertEqual(len(symbols.symbol_problems(metric(), diagram(caption=view["caption"]))), 1)

    def test_an_embedding_diagram_has_its_height_and_the_space_it_stands_in(self):
        view = {"id": "flamm", "label": "Flamm", "unit": "$r_s$", "caption": ["The paraboloid $z^2 = 4r_s(r - r_s)$."],
                "settings": "Every length is measured with $dX^2 + dY^2 - dZ^2$.",
                "figure": {"legend": [], "labels": []}, "surfaces": []}
        self.assertEqual(symbols.symbol_problems(metric(), None, None, {"views": [view]}), [])


class Collection(unittest.TestCase):
    """The rule over every spacetime on disk, and its place in the build."""

    @classmethod
    def setUpClass(cls):
        cls.metrics = build.load_metrics()
        cls.diagrams = build.load_diagrams(cls.metrics)
        cls.conformal = build.load_conformal(cls.metrics)
        cls.embedding = build.load_embedding(cls.metrics)

    def test_every_spacetime_on_disk_defines_every_symbol_it_uses(self):
        self.assertGreater(len(self.metrics), 60, "the spacetimes were not all read")
        for spacetime in self.metrics:
            with self.subTest(spacetime["id"]):
                name = spacetime["id"]
                self.assertEqual(symbols.symbol_problems(spacetime, self.diagrams.get(name), self.conformal.get(name),
                                                         self.embedding.get(name)), [])

    def test_oris_time_machine_defines_its_unit_and_keeps_his_period(self):
        ori = next(spacetime for spacetime in self.metrics if spacetime["id"] == "ori_time_machine")
        defined = set.intersection(*(symbols.chart_definitions(ori, chart) for chart in ori["coordinates"]))
        self.assertLessEqual({"L", "\\ell"}, defined)
        self.assertIn("any length $\\ell$", ori["convention"])
        for chart in ori["coordinates"]:
            self.assertIn("L", [parameter["symbol"] for parameter in chart["parameters"]])

    def test_oris_drawings_fail_without_the_unit_the_captain_missed(self):
        ori = copy.deepcopy(next(spacetime for spacetime in self.metrics if spacetime["id"] == "ori_time_machine"))
        ori["convention"] = ori["convention"].split(" The metric fixes no length")[0]
        problems = symbols.symbol_problems(ori, self.diagrams["ori_time_machine"])
        self.assertTrue(problems)
        self.assertTrue(all(" uses $\\ell$ in " in problem for problem in problems))

    def test_the_build_refuses_an_undefined_symbol(self):
        with self.assertRaises(build.DataError) as raised:
            build.check_symbols([metric()], {"x": diagram(xlabel="$r/\\ell$")}, {}, {})
        self.assertIn("x.json: diagrams spherical/radial uses $\\ell$ in xlabel", str(raised.exception))
        self.assertIsNone(build.check_symbols([metric()], {"x": diagram()}, {}, {}))

    def test_every_symbol_that_needs_no_definition_says_what_it_is(self):
        for name, meaning in list(symbols.UNIVERSAL.items()) + \
                [item for drawn in symbols.DRAWN.values() for item in drawn.items()]:
            with self.subTest(name):
                self.assertTrue(meaning.strip())
        self.assertLessEqual(len(symbols.UNIVERSAL), 24, "the list of universal symbols is meant to stay short")


if __name__ == "__main__":
    unittest.main()
