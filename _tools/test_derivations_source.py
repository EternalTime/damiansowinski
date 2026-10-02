"""The source of the scripts in _tools/derivations, read without running them. A block spliced into
conformal.py once left forty functions defined twice, where Python keeps the later one in silence,
and three more spliced into the end of a drawing took its `return views`, so the whole redraw
failed while each new spacetime drawn alone passed. Needs no sympy."""
import ast
import unittest
from collections import Counter
from pathlib import Path

DERIVATIONS = Path(__file__).resolve().parent / "derivations"
DEFINITIONS = (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)


def modules():
    return {path.name: ast.parse(path.read_text(encoding="utf-8")) for path in sorted(DERIVATIONS.glob("*.py"))}


def defined_twice(tree):
    """The top-level functions and classes a module defines more than once, with their lines."""
    lines = {}
    for node in tree.body:
        if isinstance(node, DEFINITIONS):
            lines.setdefault(node.name, []).append(node.lineno)
    return {name: at for name, at in lines.items() if len(at) > 1}


def drawings_without_return(tree):
    """The top-level functions of (ck, src), which draw one spacetime's views, that can run off
    their end and hand back None."""
    return [node.name for node in tree.body
            if isinstance(node, ast.FunctionDef) and [a.arg for a in node.args.args] == ["ck", "src"]
            and not isinstance(node.body[-1], ast.Return)]


class DerivationsSource(unittest.TestCase):
    def test_there_are_modules_to_read(self):
        self.assertIn("conformal.py", modules())

    def test_no_module_defines_a_top_level_name_twice(self):
        twice = {f"{name}: {found}" for name, tree in modules().items() if (found := defined_twice(tree))}
        self.assertEqual(twice, set(), "a function or class defined twice at the top level of a module; "
                                       "Python keeps the later one, so remove the copy it does not use")

    def test_every_drawing_returns_its_views(self):
        tree = modules()["conformal.py"]
        self.assertGreater(Counter(type(n) for n in tree.body)[ast.FunctionDef], 80)
        self.assertEqual(drawings_without_return(tree), [])

    def test_the_checks_catch_what_they_are_for(self):
        tree = ast.parse("def a(ck, src):\n    views = []\n\n\ndef b(x):\n    return x\n\n\n"
                         "class C:\n    def b(self):\n        pass\n\n\ndef b(x):\n    return 2 * x\n")
        self.assertEqual(defined_twice(tree), {"b": [5, 14]})
        self.assertEqual(drawings_without_return(tree), ["a"])


if __name__ == "__main__":
    unittest.main()
