"""Every TeX string the spacetimes page sets is one MathJax can set.

_tools/derivations/tex_check.cjs does the typesetting; this runs it with the tests, so a page
whose mathematics would print as an error box cannot land. It needs node and mathjax-full, and
fails rather than skips without them, since a check that is skipped holds nothing:

    _tools/setup-env.sh

which installs it in .node.noindex; MFS_NODE names another folder to look for it in.
test_environment.py stops the suite before this runs when either is missing.
"""
import os
import shutil
import subprocess
import unittest
from pathlib import Path

TOOLS = Path(__file__).resolve().parent
PREFIX = Path(os.environ.get("MFS_NODE", TOOLS.parent / ".node.noindex"))
INSTALL = "_tools/setup-env.sh"


class TexCheck(unittest.TestCase):
    def test_mathjax_sets_every_string_the_page_sets(self):
        self.assertIsNotNone(shutil.which("node"), "the TeX check needs node")
        self.assertTrue((PREFIX / "node_modules" / "mathjax-full" / "package.json").exists(),
                        f"the TeX check needs mathjax-full: {INSTALL}")
        run = subprocess.run(["node", str(TOOLS / "derivations" / "tex_check.cjs"), str(PREFIX)],
                             capture_output=True, text=True, timeout=1800)
        self.assertEqual(run.returncode, 0, (run.stderr + run.stdout)[-4000:])


if __name__ == "__main__":
    unittest.main()
