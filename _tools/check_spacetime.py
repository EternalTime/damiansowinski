#!/usr/bin/env python3
"""Run the checks of one spacetime's own page, in order, stopping at the first that fails.

    .venv.noindex/bin/python _tools/check_spacetime.py kerr

runs, for the spacetime whose metric file is MFS/assets/data/metrics/kerr.json:

  1. the environment, and the source of the derivation scripts;
  2. verify_metrics.py over every chart of that spacetime;
  3. the checks of null_rays.py, conformal.py, embedding.py and curvature.py for its drawings;
  4. tex_check.cjs over its own files;
  5. build_mfs_data.py --check and build_agent_data.py --check, which read the whole collection;
  6. the drawing-data tests of test_build_mfs_data.py, with MFS_SPACETIME naming it;
  7. its own test file, _tools/test_kerr.py, where it has one.

Each step's output goes to .check.noindex/<id>/<step>.log at the top of the checkout, which the
run empties first, so every log there was written by that run. A failed step is named with the
section of _tools/README.md that explains it, and the run ends on one summary line. A drawing
program that draws nothing of its kind for the spacetime refuses its name, and the step records
that refusal in its log and is counted as not drawn; a refusal while a published file of that
kind exists is a failure.

    .venv.noindex/bin/python _tools/check_spacetime.py kerr --texts

runs the build's text rules alone, over the spacetime's metric file and every text in another
metric file on why that one is related to it, through the functions build_mfs_data.py calls.
"""
import argparse
import os
import re
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path

import build_mfs_data as build

TOOLS = Path(__file__).resolve().parent
ROOT = TOOLS.parent
DERIVATIONS = TOOLS / "derivations"
PYTHON = ROOT / ".venv.noindex" / "bin" / "python"
LOGS = ROOT / ".check.noindex"
README = "_tools/README.md"
DATA = ROOT / "MFS" / "assets" / "data"

# The classes of test_build_mfs_data.py that test the drawing data of the spacetimes
# MFS_SPACETIME names.
DRAWING_TESTS = ("Diagrams", "Figures", "ConformalDiagrams", "EmbeddingDiagrams", "Slices", "StacksAndMovies",
                 "TimeSlicedViewsAreMovies")
TAIL = 25


@dataclass
class Step:
    """One check: its name, which names its log, the command, the files that command needs, the
    section of _tools/README.md that explains it, and where a drawing program draws, the folder
    of its published files."""
    name: str
    command: list
    needs: list
    section: str
    drawn: Path = None
    environment: dict = field(default_factory=dict)
    only_if: Path = None


def drawing(name, script, folder, section, metric_id):
    path = DERIVATIONS / script
    return Step(name, [PYTHON, path, "--verify", "--metric", metric_id], [PYTHON, path], section,
                drawn=DATA / folder / f"{metric_id}.json")


def unittests(name, modules, section, metric_id, only_if=None):
    """A run of the named tests from _tools, with MFS_SPACETIME naming the one spacetime."""
    files = sorted({TOOLS / f"{module.split('.')[0]}.py" for module in modules} - {only_if})
    return Step(name, [PYTHON, "-m", "unittest", *modules], [PYTHON, *files], section,
                environment={"MFS_SPACETIME": metric_id}, only_if=only_if)


def steps(metric_id):
    """The steps for one spacetime, in the order they run."""
    own = TOOLS / f"test_{metric_id}.py"
    return [
        Step("environment", [PYTHON, TOOLS / "environment.py"], [PYTHON, TOOLS / "environment.py"],
             "The test environment"),
        unittests("derivations-source", ["test_derivations_source"], "Conformal diagrams", metric_id),
        Step("verify-metrics", [PYTHON, DERIVATIONS / "verify_metrics.py", "--metric", metric_id],
             [PYTHON, DERIVATIONS / "verify_metrics.py"], "Checking the physics"),
        drawing("null-rays", "null_rays.py", "diagrams", "Checking the rays", metric_id),
        drawing("conformal", "conformal.py", "conformal", "Conformal diagrams", metric_id),
        drawing("embedding", "embedding.py", "embedding", "Checking the surface", metric_id),
        drawing("curvature", "curvature.py", "curvature", "Checking the figures", metric_id),
        Step("tex", ["node", DERIVATIONS / "tex_check.cjs", f"--metric={metric_id}"],
             ["node", DERIVATIONS / "tex_check.cjs"], "Conformal diagrams"),
        Step("build", [PYTHON, TOOLS / "build_mfs_data.py", "--check"], [PYTHON, TOOLS / "build_mfs_data.py"],
             "Checking"),
        Step("agent-build", [PYTHON, TOOLS / "build_agent_data.py", "--check"],
             [PYTHON, TOOLS / "build_agent_data.py"], "Adding a spacetime to My Favorite Spacetimes"),
        unittests("drawing-tests", [f"test_build_mfs_data.{name}" for name in DRAWING_TESTS], "Checking",
                  metric_id),
        unittests("own-tests", [own.stem], "Adding a spacetime to My Favorite Spacetimes", metric_id, only_if=own),
    ]


def place(path):
    """A path as a person would type it at the top of the checkout."""
    return str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)


def missing(step):
    """The first file or program a step's command needs that is not there, or None."""
    for need in step.needs:
        if isinstance(need, Path) and not need.exists():
            return place(need)
        if isinstance(need, str) and shutil.which(need) is None:
            return need
    return None


def shown(command):
    """A command as a person would type it at the top of the checkout."""
    return " ".join(place(part) if isinstance(part, Path) else str(part) for part in command)


def not_drawn(step, metric_id, status, output):
    """Whether a drawing program refused the spacetime as one it draws nothing of its kind for,
    with no file of that kind published for it."""
    refusal = re.compile(rf"error: no [a-z ]+ is drawn for \['{re.escape(metric_id)}'\]")
    return step.drawn is not None and status == 2 and refusal.search(output) and not step.drawn.exists()


def tests_ran(output):
    """What a unittest run says it ran, as 'N tests, K passed over', or None for another command."""
    ran = re.search(r"^Ran (\d+) tests?", output, re.MULTILINE)
    if not ran:
        return None
    skipped = re.search(r"skipped=(\d+)", output)
    return f"{ran.group(1)} tests, {skipped.group(1) if skipped else 0} passed over"


def run_step(step, log):
    """Run one step with its output in `log`, returning its exit status and its output."""
    environment = {key: value for key, value in os.environ.items() if key != "MFS_SPACETIME"}
    environment.update(step.environment)
    # A Python step's two streams share the log, so they are written unbuffered, in the order said.
    environment["PYTHONUNBUFFERED"] = "1"
    cwd = TOOLS if step.environment else ROOT
    header = [f"$ {shown(step.command)}", f"  in {cwd.relative_to(ROOT)}/"]
    header += [f"  with {key}={value}" for key, value in step.environment.items()]
    start = time.time()
    done = subprocess.run([str(part) for part in step.command], cwd=cwd, env=environment,
                          stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    footer = f"exit status {done.returncode} after {time.time() - start:.1f} s"
    log.write_text("\n".join(header) + "\n\n" + done.stdout + ("" if done.stdout.endswith("\n") else "\n")
                   + footer + "\n", encoding="utf-8")
    return done.returncode, done.stdout


def failed(metric_id, step, why, log, passed, total, start):
    """Report the step that failed, the end of its log, and the summary line, and return 1."""
    print(f"\nFAILED  {step.name}: {why}", file=sys.stderr)
    print(f'        see "{step.section}" in {README}', file=sys.stderr)
    if log is not None and log.exists():
        lines = log.read_text(encoding="utf-8").splitlines()
        print(f"        the last lines of {place(log)}:", file=sys.stderr)
        for line in lines[-TAIL:]:
            print(f"        | {line}", file=sys.stderr)
    print(f"{metric_id}: FAILED at {step.name} after {passed} of {total} steps passed, in "
          f"{time.time() - start:.0f} s; see \"{step.section}\" in {README}", file=sys.stderr)
    return 1


def check(metric_id, plan, logs):
    """Run every step of `plan` in order, stopping at the first that fails, and return the exit
    status: 0 when every step passed, 1 when one failed or its command is missing."""
    start = time.time()
    shutil.rmtree(logs, ignore_errors=True)
    logs.mkdir(parents=True)
    # Every command is looked for before any runs, so a missing one is named at once.
    for step in plan:
        absent = missing(step)
        if absent:
            return failed(metric_id, step, f"its command needs {absent}, which does not exist", None,
                          0, len(plan), start)
    print(f"Checking {metric_id}, each step's output in {place(logs)}/\n")
    passed, undrawn, unwritten = 0, [], []
    for step in plan:
        log = logs / f"{step.name}.log"
        if step.only_if is not None and not step.only_if.exists():
            where = place(step.only_if)
            log.write_text(f"{metric_id} has no {where}, so it has no tests of its own to run\n", encoding="utf-8")
            unwritten.append(step.name)
            print(f"  no file    {step.name:<20} {'':>9}  {metric_id} has no {where}")
            continue
        began = time.time()
        status, output = run_step(step, log)
        took = f"{time.time() - began:.1f} s"
        if not_drawn(step, metric_id, status, output):
            undrawn.append(step.name)
            print(f"  not drawn  {step.name:<20} {took:>9}  the program draws nothing of its kind for {metric_id}")
            continue
        if status != 0:
            return failed(metric_id, step, f"{shown(step.command)} exited with status {status}", log,
                          passed, len(plan), start)
        passed += 1
        ran = tests_ran(output)
        print(f"  ok         {step.name:<20} {took:>9}" + (f"  {ran}" if ran else ""))
    others = [f"{len(names)} {what} ({', '.join(names)})"
              for names, what in ((undrawn, "not drawn"), (unwritten, "with no file")) if names]
    print(f"\n{metric_id}: {passed} of {len(plan)} steps passed" + "".join(f", {other}" for other in others)
          + f", in {time.time() - start:.0f} s; logs in {place(logs)}/")
    return 0


def check_texts(metric_id, logs):
    """Hold the spacetime's texts to the build's text rules alone, writing what they find to
    texts.log as well, and return the exit status."""
    try:
        metrics = build.load_metrics()
    except build.DataError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    problems = build.spacetime_text_problems(metric_id, metrics)
    related = sum(1 for metric in metrics if metric["id"] != metric_id
                  for _ in build.related_texts(metric, metric_id))
    count = len([text for metric in metrics if metric["id"] == metric_id for text in build.metric_texts(metric)])
    found = f"{len(problems)} {'text breaks' if len(problems) == 1 else 'texts break'} a rule"
    summary = (f"{metric_id}: {count} texts of its metric file and {related} on why it is related to another, "
               f"{found if problems else 'every one keeps the rules'}")
    logs.mkdir(parents=True, exist_ok=True)
    (logs / "texts.log").write_text("\n".join([*problems, summary]) + "\n", encoding="utf-8")
    for problem in problems:
        print(problem, file=sys.stderr)
    print(summary, file=sys.stderr if problems else sys.stdout)
    return 1 if problems else 0


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("metric", help="the id of the spacetime, the name of its metric file")
    parser.add_argument("--texts", action="store_true", help="run the build's text rules alone")
    args = parser.parse_args(argv)
    # Each step's line shows as it finishes, in order with the lines of a failure.
    sys.stdout.reconfigure(line_buffering=True)
    if not (build.METRICS_DIR / f"{args.metric}.json").exists():
        parser.error(f"{args.metric}: no metric file {(build.METRICS_DIR / args.metric).relative_to(ROOT)}.json")
    logs = LOGS / args.metric
    if args.texts:
        return check_texts(args.metric, logs)
    return check(args.metric, steps(args.metric), logs)


if __name__ == "__main__":
    sys.exit(main())
