#!/bin/sh
# Builds the environment the tests and derivations run in, inside this checkout, from scratch:
#   .venv.noindex  Python with the packages pinned in _tools/requirements.txt
#   .node.noindex  mathjax-full, pinned in _tools/node-requirements.txt, for the TeX check
# Both folders are gitignored, so every worktree builds its own. PYTHON names another
# interpreter than python3.14 to build the venv from.
set -eu
ROOT=$(cd "$(dirname "$0")/.." && pwd)
PYTHON=${PYTHON:-python3.14}

"$PYTHON" -m venv --clear "$ROOT/.venv.noindex"
"$ROOT/.venv.noindex/bin/python" -m pip install --quiet --disable-pip-version-check \
    -r "$ROOT/_tools/requirements.txt"

npm install --silent --no-audit --no-fund --save-exact --prefix "$ROOT/.node.noindex" \
    $(cat "$ROOT/_tools/node-requirements.txt")

"$ROOT/.venv.noindex/bin/python" "$ROOT/_tools/environment.py"
