"""The two deficits and the regular axis of the cylindrical wave geometry."""

import json
import sys
import unittest
from pathlib import Path

import mpmath as mp

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / '_tools' / 'derivations'))
import string_cylindrical_waves as waves


class StringCylindricalWaves(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.metric = json.loads((ROOT / 'MFS/assets/data/metrics/string_cylindrical_waves.json').read_text())


if __name__ == '__main__':
    unittest.main()
