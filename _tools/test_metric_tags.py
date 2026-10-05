"""Tests for the tags of My Favorite Spacetimes: .venv.noindex/bin/python -m unittest discover -s _tools

Every tag the metric decides is held to what _tools/derivations/metric_tags.py computed from
the published charts: vacuum, Einstein space, conformally flat, flat, stationary, static,
spherically symmetric and the dimension. The computation needs sympy and five minutes, so
its answers are kept in _tools/derivations/metric_tags.json, each under a stamp of the chart
it was computed from, and these tests read that file. A chart that is new, or whose line
element or domains changed, has no answer under its stamp and fails here until the script
has been run again, so a new spacetime cannot be tagged by guess.

The tests that compute need sympy, which the checker's own environment has, and are skipped
where it is absent.
"""
import importlib.util
import json
import re
import shutil
import subprocess
import sys
import unittest
from pathlib import Path

DERIVATIONS = Path(__file__).resolve().parent / "derivations"
HAS_SYMPY = importlib.util.find_spec("sympy") is not None

sys.path.insert(0, str(DERIVATIONS))
import metric_tags as mt  # noqa: E402
sys.path.remove(str(DERIVATIONS))

RUN = "run .venv.noindex/bin/python _tools/derivations/metric_tags.py"


def published_metrics():
    return [json.loads(path.read_text(encoding="utf-8"))
            for path in sorted(mt.METRICS_DIR.glob("*.json"))]


def chart(**facts):
    """The facts of a chart of four dimensions with nothing to say, and then `facts`."""
    return {"dimension": 4, "ricci_flat": False, "einstein": None, "conformally_flat": False,
            "killing": [], "sphere": None, **facts}


def killing(coord="t", timelike=True, orthogonal=True, translation=True):
    return {"coord": coord, "timelike": timelike, "hypersurface_orthogonal": orthogonal,
            "translation": translation}


class PublishedTags(unittest.TestCase):
    """The collection as it stands."""

    @classmethod
    def setUpClass(cls):
        cls.metrics = published_metrics()
        cls.facts = mt.load_facts()

    def test_every_published_chart_has_facts_computed_from_it_as_it_stands(self):
        stale = []
        for metric_id, entry in mt.published_charts():
            name = f"{metric_id}/{entry['id']}"
            if name not in self.facts:
                stale.append(f"{name} has never been computed")
            elif self.facts[name]["stamp"] != mt.stamp(entry):
                stale.append(f"{name} has changed since it was computed")
        self.assertEqual(stale, [], f"{RUN}")

    def test_the_facts_hold_no_chart_that_is_not_published(self):
        published = {f"{metric_id}/{entry['id']}" for metric_id, entry in mt.published_charts()}
        self.assertEqual(sorted(set(self.facts) - published), [], f"{RUN}")

    def test_every_spacetime_carries_exactly_the_tags_its_charts_decide(self):
        faults = {}
        for metric in self.metrics:
            if all(f"{metric['id']}/{e['id']}" in self.facts for e in metric["coordinates"]):
                found = mt.tag_faults(metric, self.facts)
                if found:
                    faults[metric["id"]] = found
        self.assertEqual(faults, {})

    def test_the_tag_this_was_written_for_is_on_every_conformally_flat_spacetime(self):
        """Three spacetimes carried it on 2 October 2026, which the captain saw could not be right."""
        tagged = {m["id"] for m in self.metrics if "conformally flat" in m["tags"]}
        for metric_id in ("minkowski", "de_sitter", "anti_de_sitter", "milne", "frw",
                          "einstein_static", "bertotti_robinson", "bell_szekeres",
                          "interior_schwarzschild"):
            self.assertIn(metric_id, tagged)
        for metric_id in ("schwarzschild", "kerr", "kasner", "oppenheimer_snyder"):
            self.assertNotIn(metric_id, tagged)

    def test_every_region_names_charts_its_spacetime_publishes_and_says_why(self):
        charts = {m["id"]: [e["id"] for e in m["coordinates"]] for m in self.metrics}
        for metric_id, declared in mt.REGIONS.items():
            self.assertIn(metric_id, charts)
            named = [c for region in declared["regions"] for c in region]
            self.assertEqual(len(named), len(set(named)), metric_id)
            self.assertLessEqual(set(named), set(charts[metric_id]), metric_id)
            self.assertTrue(all(declared["regions"]) and declared["why"], metric_id)
            self.assertNotEqual(declared["regions"], [charts[metric_id]],
                                f"{metric_id} declares what is taken anyway")

    def test_every_overruling_overrules_something(self):
        """One the charts have come to agree with is stale, and would hide a later change."""
        by_id = {m["id"]: m for m in self.metrics}
        for (metric_id, tag), (value, why) in mt.OVERRULED.items():
            self.assertIn(metric_id, by_id)
            self.assertIn(tag, mt.OWNED)
            self.assertTrue(why)
            kept, mt.OVERRULED = mt.OVERRULED, {}
            try:
                owned, _ = mt.decided_tags(metric_id, mt.facts_of(by_id[metric_id], self.facts))
            finally:
                mt.OVERRULED = kept
            self.assertNotEqual(owned[tag], value, f"{metric_id}: {tag}")

    def test_no_merged_tag_is_a_tag_something_was_merged_into(self):
        self.assertEqual(set(mt.RETIRED) & {into for into in mt.RETIRED.values() if into}, set())

    def test_every_tag_something_was_merged_into_is_carried(self):
        carried = {tag for metric in self.metrics for tag in metric["tags"]}
        self.assertEqual({into for into in mt.RETIRED.values() if into} - carried, set())

    def test_no_keyword_is_carried_by_one_spacetime_or_spelled_like_another(self):
        """610 keywords on 3 October 2026, 259 of them on one spacetime, and among them
        "point mass" and "point masses" and "alternative theory of gravity" and "alternative
        theories of gravity", until the captain had the whole collection revisited."""
        self.assertEqual(mt.keyword_faults(self.metrics), [])


class Rules(unittest.TestCase):
    """decided_tags() on charts made up for the purpose."""

    def owned(self, charts, metric_id="made_up"):
        return mt.decided_tags(metric_id, charts)[0]

    def test_a_curvature_tag_needs_every_chart(self):
        flat = chart(ricci_flat=True, conformally_flat=True)
        self.assertTrue(self.owned({"a": flat, "b": flat})["flat"])
        owned = self.owned({"a": flat, "b": chart(conformally_flat=True)})
        self.assertEqual((owned["vacuum"], owned["flat"], owned["conformally flat"]), (False, False, True))

    def test_a_symmetry_needs_one_chart_to_show_it(self):
        owned = self.owned({"a": chart(killing=[killing()]), "b": chart()})
        self.assertTrue(owned["stationary"] and owned["static"])

    def test_static_needs_the_time_translation_to_be_hypersurface_orthogonal(self):
        owned = self.owned({"a": chart(killing=[killing(orthogonal=False)])})
        self.assertEqual((owned["stationary"], owned["static"]), (True, False))

    def test_a_killing_vector_that_is_spacelike_or_no_translation_of_its_chart_is_no_time_translation(self):
        for vector in (killing(timelike=False), killing(translation=False)):
            self.assertFalse(self.owned({"a": chart(killing=[vector])})["stationary"])

    def test_a_sphere_has_to_fill_all_but_two_dimensions(self):
        self.assertTrue(self.owned({"a": chart(sphere=["\\theta", "\\phi"])})["spherically symmetric"])
        self.assertTrue(self.owned({"a": chart(sphere=["\\chi", "\\theta", "\\phi"])})["spherically symmetric"])
        string = chart(dimension=5, sphere=["\\theta", "\\phi"])
        self.assertFalse(self.owned({"a": string})["spherically symmetric"])

    def test_the_dimension_of_every_chart_is_a_tag(self):
        owned = self.owned({"a": chart(dimension=5), "b": chart(dimension=6)})
        self.assertEqual([tag for tag in mt.DIMENSION.values() if owned[tag]],
                         ["five-dimensional", "six-dimensional"])
        with self.assertRaises(ValueError):
            self.owned({"a": chart(dimension=12)})

    def test_an_einstein_space_implies_its_cosmological_constant_and_its_sign(self):
        negative = chart(einstein={"k": "-3/L**2", "sign": -1}, conformally_flat=True)
        owned, implied = mt.decided_tags("made_up", {"a": negative})
        self.assertTrue(owned["Einstein space"] and not owned["vacuum"])
        self.assertEqual(implied, {"cosmological constant", "negative cosmological constant"})
        either = chart(einstein={"k": "Lambda", "sign": None})
        self.assertEqual(mt.decided_tags("made_up", {"a": either})[1], {"cosmological constant"})

    def test_a_symmetry_of_a_spacetime_of_two_regions_has_to_show_in_both(self):
        """Oppenheimer and Snyder's exterior is static and their dust ball is not."""
        ball = chart(conformally_flat=True, sphere=["\\chi", "\\theta", "\\phi"])
        outside = chart(ricci_flat=True, killing=[killing()], sphere=["\\theta", "\\phi"])
        owned = self.owned({"interior_comoving": ball, "exterior_schwarzschild": outside},
                           "oppenheimer_snyder")
        self.assertTrue(owned["spherically symmetric"])
        for tag in ("static", "stationary", "vacuum", "conformally flat"):
            self.assertFalse(owned[tag], tag)

    def test_a_chart_left_out_of_its_spacetimes_regions_is_not_counted(self):
        cone = chart(ricci_flat=True, conformally_flat=True, killing=[killing()])
        owned = self.owned({"conical": cone, "interior_cap": chart()}, "cosmic_string")
        self.assertTrue(owned["flat"] and owned["static"])

    def test_tag_faults_names_what_is_missing_what_is_wrong_and_what_was_merged(self):
        facts = {"x/a": chart(ricci_flat=True)}
        metric = {"id": "x", "coordinates": [{"id": "a"}],
                  "tags": ["vacuum solution", "static", "four-dimensional", "four-dimensional"]}
        faults = mt.tag_faults(metric, facts)
        self.assertIn("lacks 'vacuum', which its charts decide it has", faults)
        self.assertIn("carries 'static', which its charts decide it has not", faults)
        self.assertIn("carries 'vacuum solution', which was merged into 'vacuum'", faults)
        self.assertIn("carries 'four-dimensional' twice", faults)
        self.assertEqual(mt.tag_faults(dict(metric, tags=["vacuum", "four-dimensional"]), facts), [])


class Keywords(unittest.TestCase):
    """keyword_faults() on collections made up for the purpose."""

    def faults(self, **tags):
        return mt.keyword_faults([{"id": metric_id, "tags": carried} for metric_id, carried in tags.items()])

    def test_a_keyword_on_one_spacetime_is_a_fault_and_on_two_is_not(self):
        self.assertEqual(self.faults(a=["vacuum", "geon"], b=["vacuum"]), ["'geon' is carried by a alone"])
        self.assertEqual(self.faults(a=["geon"], b=["geon"]), [])

    def test_a_tag_the_charts_decide_may_be_carried_alone(self):
        self.assertEqual(self.faults(a=["eleven-dimensional", "static"], b=[]), [])

    def test_a_spacetime_listing_a_keyword_twice_still_carries_it_alone(self):
        self.assertEqual(self.faults(a=["geon", "geon"], b=[]), ["'geon' is carried by a alone"])

    def test_two_keywords_differing_in_case_plural_or_punctuation_are_a_fault(self):
        for one, other in (("Black hole", "black hole"), ("point mass", "point masses"),
                           ("alternative theory of gravity", "alternative theories of gravity"),
                           ("pp-wave", "pp wave"), ("Painlevé-Gullstrand", "Painleve Gullstrand"),
                           ("2+1 dimensional gravity", "2+1-dimensional gravity")):
            with self.subTest(one=one, other=other):
                faults = self.faults(a=[one, other], b=[one, other])
                self.assertEqual(len(faults), 1)
                self.assertIn("spelled alike", faults[0])

    def test_keywords_of_different_words_are_not_spelled_alike(self):
        for one, other in (("compact", "compactness"), ("radiating", "radiation"), ("torus", "tori"),
                           ("Kasner", "Kerr"), ("dust", "null dust")):
            self.assertNotEqual(mt.spelling(one), mt.spelling(other), (one, other))


class Domains(unittest.TestCase):
    """What the script reads off a chart's domains, which needs no sympy."""

    def test_an_interval_is_read_with_its_ends_as_written(self):
        self.assertEqual(mt.interval("t \\in (-\\infty, \\infty)", "t"), ("-\\infty", "\\infty"))
        self.assertEqual(mt.interval("r \\in \\left[0,\\, \\sqrt{3/\\Lambda}\\right)", "r"),
                         ("0", "\\sqrt{3/\\Lambda}"))
        self.assertEqual(mt.interval("\\phi \\in [0, 2\\pi)", "\\phi"), ("0", "2\\pi"))
        self.assertIsNone(mt.interval("r \\in [0, \\infty) \\;\\text{for}\\; k \\le 0", "r"))
        self.assertIsNone(mt.interval("\\theta \\in [0, \\pi]", "t"))

    def test_a_translation_runs_over_the_whole_line_and_is_named_by_no_other_entry(self):
        whole = {"domains": ["t \\in (-\\infty, \\infty)", "r \\in (r_s, \\infty)"]}
        self.assertTrue(mt.is_translation(whole, "t"))
        self.assertFalse(mt.is_translation(whole, "r"))
        milne = {"domains": ["T \\in (0, \\infty)", "R \\in [0, cT)"]}
        self.assertFalse(mt.is_translation(milne, "T"))
        collapse = {"domains": ["t \\in (-\\infty, \\infty)", "r \\in [R(t), \\infty)"]}
        self.assertFalse(mt.is_translation(collapse, "t"))
        strings = {"domains": ["t \\in (-\\infty, \\infty)", "x \\in (-\\infty, \\infty)",
                               "x = \\pm vt,\\; y = \\pm d \\;\\text{(the strings)}"]}
        self.assertFalse(mt.is_translation(strings, "t"))

    def test_a_remark_in_parentheses_does_not_name_a_coordinate(self):
        remarked = {"domains": [
            "t \\in (-\\infty, \\infty)", "\\theta \\in [0, \\pi]",
            "r < a/b \\;\\text{(the circles of constant } t, r, z \\text{ are closed timelike curves)}",
        ]}
        self.assertTrue(mt.is_translation(remarked, "t"))

    def test_a_whole_sphere_has_its_last_angle_run_once_round(self):
        whole = {"domains": ["\\theta \\in [0, \\pi]", "\\phi \\in [0, 2\\pi)"]}
        self.assertTrue(mt.is_whole_sphere(whole, ["\\theta", "\\phi"]))
        wedge = {"domains": ["\\theta \\in [0, \\pi]", "\\tilde\\phi \\in [0, 2\\pi b)"]}
        self.assertFalse(mt.is_whole_sphere(wedge, ["\\theta", "\\tilde\\phi"]))
        cap = {"domains": ["\\chi \\in [0, \\chi_0]", "\\theta \\in [0, \\pi]", "\\phi \\in [0, 2\\pi)"]}
        self.assertFalse(mt.is_whole_sphere(cap, ["\\chi", "\\theta", "\\phi"]))
        self.assertTrue(mt.is_whole_sphere(cap, ["\\theta", "\\phi"]))

    def test_the_stamp_changes_with_the_line_element_and_with_the_domains(self):
        entry = {"coords": ["t", "r"], "domains": ["t \\in (-\\infty, \\infty)"],
                 "line_element": "ds^2 = -dt^2 + dr^2", "parameters": []}
        self.assertEqual(mt.stamp(entry), mt.stamp(dict(entry, name="any", convention="other")))
        self.assertNotEqual(mt.stamp(entry), mt.stamp(dict(entry, line_element="ds^2 = -dt^2 + 2dr^2")))
        self.assertNotEqual(mt.stamp(entry), mt.stamp(dict(entry, domains=["t \\in (0, \\infty)"])))


class PageFilter(unittest.TestCase):
    """The page's own rule for which tags a search finds, run in Node over the published tags."""

    LAYOUT = Path(__file__).resolve().parents[1] / "_layouts" / "mfs.html"

    def found(self, text):
        """The ids of the spacetimes a tag of which the page's rule finds for the text typed."""
        if shutil.which("node") is None:
            self.skipTest("node is not installed")
        source = self.LAYOUT.read_text(encoding="utf-8")
        function = re.search(r"function mfsTagMatches\(tag, q\) \{.*?\n      \}\n", source, re.S)
        self.assertIsNotNone(function, "the page no longer defines mfsTagMatches")
        self.assertIn("mfsTagMatches(t, q)", source.replace(function.group(0), ""),
                      "the page's search no longer asks mfsTagMatches")
        index = [{"id": m["id"], "tags": m["tags"]} for m in published_metrics()]
        script = (function.group(0) + "const [index, q] = JSON.parse(process.argv[1]);"
                  "console.log(JSON.stringify(index.filter(m => m.tags.some(t => mfsTagMatches(t, q)))"
                  ".map(m => m.id)));")
        run = subprocess.run(["node", "-e", script, json.dumps([index, text.strip().lower()])],
                             capture_output=True, text=True, check=True)
        return set(json.loads(run.stdout))

    def carrying(self, *tags):
        return {m["id"] for m in published_metrics() if set(tags) & set(m["tags"])}

    def test_a_tag_typed_whole_finds_the_spacetimes_that_carry_it_and_no_others(self):
        """No others but those with a tag the text begins a word of, as "vacuum decay" for "vacuum"."""
        for tag in ("conformally flat", "vacuum", "static", "Einstein space", "three-dimensional"):
            begins = {m["id"] for m in published_metrics()
                      if any(re.search(r"(^|[ /-])" + re.escape(tag.lower()), t.lower()) for t in m["tags"])}
            found = self.found(tag)
            self.assertLessEqual(self.carrying(tag), found, tag)
            self.assertEqual(found, begins, tag)
        self.assertEqual(self.found("conformally flat"), self.carrying("conformally flat"))

    def test_a_tag_with_capitals_is_found_however_it_is_typed(self):
        for typed, tag in (("AdS/CFT", "AdS/CFT"), ("ads/cft", "AdS/CFT"), ("Petrov type D", "Petrov type D"),
                           ("petrov type d", "Petrov type D"), ("EINSTEIN SPACE", "Einstein space")):
            self.assertEqual(self.found(typed), self.carrying(tag), typed)
            self.assertTrue(self.found(typed), typed)

    def test_a_word_inside_a_tag_finds_it_and_the_middle_of_a_word_does_not(self):
        self.assertLessEqual(self.carrying("pp-wave", "gravitational wave"), self.found("wave"))
        self.assertLessEqual(self.carrying("AdS/CFT"), self.found("cft"))
        self.assertLessEqual(self.carrying("conformally flat", "asymptotically flat", "flat"), self.found("flat"))
        self.assertEqual(self.found("vacuum") & (self.carrying("electrovacuum") - self.carrying("vacuum")), set())
        no_hole = {m["id"] for m in published_metrics()
                   if not any(word.startswith("hole") for tag in m["tags"] for word in re.split(r"[ /-]", tag))}
        self.assertTrue(self.carrying("wormhole") & no_hole)
        self.assertEqual(self.found("hole") & no_hole, set())


@unittest.skipUnless(HAS_SYMPY, "metric_tags.py computes with sympy")
class Computed(unittest.TestCase):
    """A few quick charts computed again and held to the file, so that the script and the
    answers it once wrote cannot drift apart. The whole collection takes five minutes:
    .venv.noindex/bin/python _tools/derivations/metric_tags.py --check
    """

    CHARTS = (
        "minkowski/spherical",              # flat, static, a sphere
        "schwarzschild/spherical",          # vacuum and not conformally flat
        "de_sitter/flat_slicing",           # an Einstein space with k > 0
        "btz/stationary",                   # three dimensions, where Cotton decides
        "kasner/cartesian",                 # vacuum only on the Kasner circle
        "milne/inertial",                   # a timelike Killing vector that is no translation
        "oppenheimer_snyder/interior_comoving",  # a cap of a three sphere, a whole two sphere
        "string_black_hole/wedge",          # a sphere with a wedge missing
        "tangherlini/spherical_six",        # a four sphere in six dimensions
        "spinning_string/proper_radius",    # locally static about a string that spins
        "ab_metrics/a2_kruskal",            # vacuum by Lambert's W, which sympy leaves uncancelled
    )

    @classmethod
    def setUpClass(cls):
        sys.path.insert(0, str(DERIVATIONS))
        cls.facts = mt.load_facts()
        cls.entries = {f"{metric_id}/{entry['id']}": (metric_id, entry)
                       for metric_id, entry in mt.published_charts()}

    @classmethod
    def tearDownClass(cls):
        sys.path.remove(str(DERIVATIONS))

    def test_the_charts_computed_again_are_what_the_file_holds(self):
        for name in self.CHARTS:
            with self.subTest(chart=name):
                metric_id, entry = self.entries[name]
                self.assertEqual(mt.compute(metric_id, entry, 120), self.facts[name])

    def test_a_number_is_surely_not_zero_only_if_it_is_rational_or_far_from_zero(self):
        import sympy as sp
        x = sp.Rational(133, 10) / sp.E
        self.assertFalse(mt._not_zero(sp.Integer(0)))
        self.assertTrue(mt._not_zero(sp.Rational(1, 10**60)))
        self.assertFalse(mt._not_zero(sp.LambertW(x) * sp.exp(sp.LambertW(x)) / x - 1))
        self.assertTrue(mt._not_zero(sp.LambertW(x) - 1))

    def test_the_relations_the_file_was_computed_under_are_the_checkers(self):
        import verify_metrics as vm
        for name, (metric_id, entry) in self.entries.items():
            declared = vm.PARAMETER_RELATIONS.get((metric_id, entry["id"]), {})
            self.assertEqual(self.facts[name]["relations"], declared, f"{name}: {RUN}")


if __name__ == "__main__":
    unittest.main()
