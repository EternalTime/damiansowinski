"""The environment the tests and derivations run in, and the one check that it is whole.

_tools/setup-env.sh builds it inside the checkout: a Python venv in .venv.noindex with the
packages pinned in _tools/requirements.txt, and mathjax-full in .node.noindex at the version
pinned in _tools/node-requirements.txt. It used to live in /tmp/mfs-venv, until on 5 October
2026 macOS's cleaner of /tmp deleted that venv's pyvenv.cfg and package files and the suite
ran as bare Python, its sympy tests skipped and nothing failing.

test_environment.py runs this check as the suite is gathered and stops the suite on the
first problem, so a missing or broken environment fails loudly and never skips. Run directly,
it prints what is wrong, or that the environment is whole:

    .venv.noindex/bin/python _tools/environment.py

MFS_NODE names another folder than .node.noindex to look for mathjax-full in.
"""
import importlib
import importlib.metadata
import json
import os
import shutil
import sys
from pathlib import Path

TOOLS = Path(__file__).resolve().parent
ROOT = TOOLS.parent
VENV = ROOT / ".venv.noindex"
NODE = Path(os.environ.get("MFS_NODE", ROOT / ".node.noindex"))
SETUP = "_tools/setup-env.sh"
RUN = ".venv.noindex/bin/python -m unittest discover -s _tools"


def pins(path, separator):
    """The (name, version) pairs a requirements file pins, comments and blank lines skipped."""
    lines = (line.split("#")[0].strip() for line in path.read_text().splitlines())
    return [tuple(line.split(separator)) for line in lines if line]


def problems():
    """Everything wrong with the environment this interpreter runs in, empty when it is whole."""
    found = []
    if not (VENV / "pyvenv.cfg").is_file():
        found.append(f"{VENV}/pyvenv.cfg is missing: the venv is gone, or runs as bare Python")
    elif Path(sys.prefix).resolve() != VENV.resolve():
        found.append(f"this is {sys.executable}, not the checkout's {VENV}/bin/python")
    for name, version in pins(TOOLS / "requirements.txt", "=="):
        try:
            installed = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            found.append(f"{name} is not installed (pinned at {version})")
            continue
        if installed != version:
            found.append(f"{name} is {installed}, pinned at {version}")
            continue
        try:
            importlib.import_module(name)
        except Exception as error:
            found.append(f"{name} {version} is installed but does not import: {error!r}")
    if shutil.which("node") is None:
        found.append("node is not on the PATH")
    for name, version in pins(TOOLS / "node-requirements.txt", "@"):
        package = NODE / "node_modules" / name / "package.json"
        if not package.is_file():
            found.append(f"{name} is not installed in {NODE} (pinned at {version})")
        elif (installed := json.loads(package.read_text())["version"]) != version:
            found.append(f"{name} in {NODE} is {installed}, pinned at {version}")
    return found


def report(found):
    """The message that names each problem and the command that mends them all."""
    lines = ["The test environment is missing or broken:"]
    lines += [f"  - {problem}" for problem in found]
    lines += ["Build it with", f"  {SETUP}", "and run the suite with", f"  {RUN}"]
    return "\n".join(lines)


if __name__ == "__main__":
    found = problems()
    print(report(found) if found else f"The test environment in {VENV} and {NODE} is whole.",
          file=sys.stderr if found else sys.stdout)
    sys.exit(1 if found else 0)
