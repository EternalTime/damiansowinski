"""The command that checks one spacetime's page, _tools/check_spacetime.py:
.venv.noindex/bin/python -m unittest discover -s _tools

The steps of a real spacetime are held to commands that exist, to the option each needs and to
sections of _tools/README.md that exist. The running of steps is held to steps made up here, which
run in a second: a run that passes, one that stops at a failure, one whose command is missing, a
drawing program that draws nothing for the spacetime, and a spacetime with no test file of its own.
"""
import ast
import contextlib
import copy
import io
import re
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import build_mfs_data as build
import check_spacetime as cs

README = cs.ROOT / cs.README


def fake(name, code, **fields):
    """A step that runs `code` in this interpreter."""
    return cs.Step(name, [Path(sys.executable), "-c", code], [Path(sys.executable)], "Checking", **fields)


def run(plan, metric_id="kerr"):
    """Run `plan` into a folder of its own, returning the status, what it printed and the folder."""
    logs = Path(tempfile.mkdtemp()) / metric_id
    out, err = io.StringIO(), io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
        status = cs.check(metric_id, plan, logs)
    return status, out.getvalue() + err.getvalue(), logs


class TheSteps(unittest.TestCase):
    """The steps of a real spacetime, in the order and with the options the README gives."""

    def setUp(self):
        self.plan = cs.steps("kerr")

    def test_the_steps_run_in_order(self):
        self.assertEqual([step.name for step in self.plan],
                         ["environment", "derivations-source", "verify-metrics", "null-rays", "conformal",
                          "embedding", "curvature", "tex", "build", "agent-build", "drawing-tests", "own-tests"])

    def test_every_command_a_step_needs_exists(self):
        for step in self.plan:
            needs = [need for need in step.needs if not (isinstance(need, Path) and need == cs.PYTHON)]
            with self.subTest(step=step.name):
                self.assertIsNone(cs.missing(cs.Step(step.name, step.command, needs, step.section)))

    def test_every_section_named_is_a_heading_of_the_readme(self):
        headings = set(re.findall(r"^#+ (.+)$", README.read_text(encoding="utf-8"), re.MULTILINE))
        for step in self.plan:
            with self.subTest(step=step.name):
                self.assertIn(step.section, headings)

    def test_every_drawing_program_takes_verify_and_metric_and_refuses_what_it_does_not_draw(self):
        for step in self.plan:
            if step.drawn is None:
                continue
            source = step.command[1].read_text(encoding="utf-8")
            with self.subTest(step=step.name):
                self.assertEqual(step.command[2:], ["--verify", "--metric", "kerr"])
                self.assertIn('add_argument("--verify"', source)
                self.assertIn('add_argument("--metric"', source)
                said = re.search(r'parser\.error\(f"(no [a-z ]+ is drawn for )\{sorted\(unknown\)\}"\)', source)
                self.assertIsNotNone(said)
                self.assertTrue(cs.not_drawn(cs.Step("x", [], [], "", drawn=Path("/nowhere")), "kerr", 2,
                                             f"error: {said.group(1)}['kerr']"))

    def test_the_other_commands_take_the_options_they_are_given(self):
        sources = {step.name: step.command[1].read_text(encoding="utf-8") for step in self.plan
                   if step.name in ("verify-metrics", "tex", "build", "agent-build")}
        self.assertIn('add_argument("--metric"', sources["verify-metrics"])
        self.assertIn("'--metric='", sources["tex"])
        self.assertIn('add_argument("--check"', sources["build"])
        self.assertIn('add_argument("--check"', sources["agent-build"])

    def test_the_drawing_tests_are_classes_of_test_build_mfs_data(self):
        tree = ast.parse((cs.TOOLS / "test_build_mfs_data.py").read_text(encoding="utf-8"))
        classes = {node.name for node in tree.body if isinstance(node, ast.ClassDef)}
        self.assertEqual(set(cs.DRAWING_TESTS) - classes, set())
        step = next(step for step in self.plan if step.name == "drawing-tests")
        self.assertEqual(step.environment, {"MFS_SPACETIME": "kerr"})

    def test_a_spacetime_runs_its_own_test_file_where_it_has_one(self):
        own = next(step for step in cs.steps("kerr_newman") if step.name == "own-tests")
        self.assertEqual(own.command[-1], "test_kerr_newman")
        self.assertEqual(own.only_if, cs.TOOLS / "test_kerr_newman.py")


class Running(unittest.TestCase):
    """The running of steps, held to steps made up here."""

    def test_a_run_that_passes_writes_every_log_and_ends_on_one_summary_line(self):
        status, said, logs = run([fake("one", "print('first')"), fake("two", "print('second')")])
        self.assertEqual(status, 0)
        self.assertIn("first", (logs / "one.log").read_text())
        self.assertIn("second", (logs / "two.log").read_text())
        self.assertRegex(said.splitlines()[-1], r"^kerr: 2 of 2 steps passed, in \d+ s; logs in ")

    def test_a_run_stops_at_the_first_failure_naming_the_step_and_its_section(self):
        status, said, logs = run([fake("one", "print('first')"), fake("two", "import sys; sys.exit(3)"),
                                  fake("three", "print('third')")])
        self.assertEqual(status, 1)
        self.assertIn("FAILED  two:", said)
        self.assertIn('see "Checking" in _tools/README.md', said)
        self.assertRegex(said.splitlines()[-1], r'^kerr: FAILED at two after 1 of 3 steps passed')
        self.assertFalse((logs / "three.log").exists())

    def test_a_missing_command_is_named_before_any_step_runs(self):
        gone = cs.Step("gone", [Path("/nowhere/script.py")], [Path("/nowhere/script.py")], "Checking")
        status, said, logs = run([fake("one", "print('first')"), gone])
        self.assertEqual(status, 1)
        self.assertIn("FAILED  gone: its command needs /nowhere/script.py, which does not exist", said)
        self.assertFalse((logs / "one.log").exists())

    def test_a_drawing_program_that_draws_nothing_for_the_spacetime_is_not_drawn_and_logged(self):
        refuse = "import sys; print(\"error: no curvature figure is drawn for ['kerr']\"); sys.exit(2)"
        status, said, logs = run([fake("curvature", refuse, drawn=Path("/nowhere/kerr.json"))])
        self.assertEqual(status, 0)
        self.assertIn("is drawn for", (logs / "curvature.log").read_text())
        self.assertIn("kerr: 0 of 1 steps passed, 1 not drawn (curvature)", said)

    def test_that_refusal_is_a_failure_where_a_file_of_that_kind_is_published(self):
        refuse = "import sys; print(\"error: no curvature figure is drawn for ['kerr']\"); sys.exit(2)"
        status, said, _ = run([fake("curvature", refuse, drawn=Path(__file__))])
        self.assertEqual(status, 1)
        self.assertIn("FAILED  curvature:", said)

    def test_a_spacetime_with_no_test_file_of_its_own_says_so_in_the_log(self):
        status, said, logs = run([fake("own-tests", "print('ran')", only_if=Path("/nowhere/test_kerr.py"))])
        self.assertEqual(status, 0)
        self.assertIn("has no /nowhere/test_kerr.py", (logs / "own-tests.log").read_text())
        self.assertIn("1 with no file (own-tests)", said)

    def test_a_run_empties_the_logs_of_the_last_one(self):
        _, _, logs = run([fake("one", "print('first')")])
        (logs / "stale.log").write_text("old\n")
        with contextlib.redirect_stdout(io.StringIO()):
            cs.check("kerr", [fake("one", "print('first')")], logs)
        self.assertFalse((logs / "stale.log").exists())

    def test_a_unittest_step_says_how_many_tests_ran_and_how_many_were_passed_over(self):
        self.assertEqual(cs.tests_ran("...\nRan 12 tests in 0.3s\n\nOK (skipped=5)\n"), "12 tests, 5 passed over")
        self.assertIsNone(cs.tests_ran("every check passed\n"))


class Texts(unittest.TestCase):
    """--texts holds the spacetime's texts to the build's rules alone."""

    def texts(self, metrics):
        logs = Path(tempfile.mkdtemp()) / "kerr"
        out = io.StringIO()
        with mock.patch.object(build, "load_metrics", return_value=metrics), \
                contextlib.redirect_stdout(out), contextlib.redirect_stderr(out):
            status = cs.check_texts("kerr", logs)
        return status, out.getvalue(), (logs / "texts.log").read_text()

    def test_a_fault_in_the_metric_file_or_in_a_text_relating_another_to_it_is_named(self):
        metrics = copy.deepcopy(build.load_metrics())
        kerr = next(m for m in metrics if m["id"] == "kerr")
        kerr["description"] += " It is important to note this entry."
        other = next(m for m in metrics if m["id"] == kerr["related"][0]["id"])
        next(e for e in other["related"] if e["id"] == "kerr")["text"] = "A well-known limit."
        status, said, log = self.texts(metrics)
        self.assertEqual(status, 1)
        self.assertIn("kerr.json: description has a word for the collection's machinery: 'entry'", said)
        self.assertIn(f"{other['id']}.json: related kerr has a compound spelled with a hyphen: 'well-known'", said)
        self.assertIn("2 texts break a rule", log)

    def test_texts_that_keep_the_rules_pass(self):
        status, said, _ = self.texts(build.load_metrics())
        self.assertEqual(status, 0)
        self.assertIn("every one keeps the rules", said)


if __name__ == "__main__":
    unittest.main()
