"""Checks of the stringy cosmic string's sourced asymptotic angle and drawn moment."""
import json
import unittest
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'MFS/assets/data'
METRIC = 'stringy_cosmic_strings'


def radial_slope(n, radius):
    """For equally spaced unit poles, H is proportional to (r**n - 1)**(-1/6).

    The torus factor tends to a nonzero constant at infinity. For g=H(dr^2+r^2 dphi^2),
    this is d(r sqrt(H))/d(proper radius), with its sign retained at the far pole.
    """
    r = Fraction(radius)
    p = r**n
    return 1 - Fraction(n, 12) * p / (p - 1)


class StringyCosmicStrings(unittest.TestCase):
    def test_sourced_factor_and_asymptotic_deficit(self):
        metric = json.loads((DATA / 'metrics' / f'{METRIC}.json').read_text())
        self.assertIn(r'\prod_i|z-z_i|^{-1/6}', metric['history'])
        for n, expected_deficit in ((1, Fraction(1, 6)), (2, Fraction(1, 3)),
                                    (6, Fraction(1)), (12, Fraction(2)), (24, Fraction(4))):
            with self.subTest(n=n):
                slope = radial_slope(n, 10**6)
                self.assertLess(abs(slope - (1 - Fraction(n, 12))), Fraction(1, 10**6))
                self.assertLess(abs(2 * (1 - slope) - expected_deficit), Fraction(1, 10**6))
        self.assertEqual(1 - Fraction(12, 12), 0)  # limiting cylinder
        self.assertEqual(abs(1 - Fraction(24, 12)), 1)  # regular other pole of the sphere
        self.assertEqual(Fraction(24, 6), 4)  # total deficit in units of pi

    def test_embedding_moment_on_all_drawn_planes(self):
        embedded = json.loads((DATA / 'embedding' / f'{METRIC}.json').read_text())
        view = embedded['views'][0]
        self.assertEqual(view['id'], 'six_string_cone')
        self.assertTrue(view['surfaces'][0]['pieces'])
        diagrams = json.loads((DATA / 'diagrams' / f'{METRIC}.json').read_text())
        conformal = json.loads((DATA / 'conformal' / f'{METRIC}.json').read_text())
        def views(obj):
            if isinstance(obj, dict):
                if obj.get('id') and 'caption' in obj:
                    yield obj
                for value in obj.values():
                    yield from views(value)
            elif isinstance(obj, list):
                for value in obj:
                    yield from views(value)
        for drawn in [*views(diagrams), *views(conformal)]:
            with self.subTest(view=drawn['id']):
                slices = drawn.get('slices', [])
                self.assertTrue(any(s['view'] == view['id'] and (s.get('lines') or s.get('points'))
                                    for s in slices))


if __name__ == '__main__':
    unittest.main()
