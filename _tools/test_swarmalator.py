"""Deterministic checks of the Swarmalator numerical core."""
import subprocess
import unittest
from pathlib import Path


class Swarmalator(unittest.TestCase):
    def test_numerical_model(self):
        root = Path(__file__).resolve().parents[1]
        result = subprocess.run(
            ["node", str(root / "_tools/test_swarmalator.cjs")],
            cwd=root, capture_output=True, text=True, timeout=180,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)


if __name__ == "__main__":
    unittest.main()
