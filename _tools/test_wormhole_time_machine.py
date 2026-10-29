"""The wormhole time machine: the trip of its moving mouth and the closed curves it makes.

The readers of the published files and MFS_SPACETIME are those of test_build_mfs_data.py."""
import unittest

from test_build_mfs_data import checks, diagram_file


@checks("wormhole_time_machine")
class WormholeTrip(unittest.TestCase):
    """The round trip of the wormhole time machine's right mouth, as its figure in three
    dimensions publishes it in (Z, X, T) with c = 1: the physics the figure states, held to the
    numbers on disk."""

    def setUp(self):
        figure = diagram_file("wormhole_time_machine")["projections"]["lorentz"][0]
        self.lines = {}
        for line in figure["turn"]["lines"]:
            self.lines.setdefault(line["class"], []).append(line["points"])

    @staticmethod
    def interval(a, b):
        return (b[0] - a[0]) ** 2 + (b[1] - a[1]) ** 2 - (b[2] - a[2]) ** 2

    def test_both_mouths_move_slower_than_light_and_the_right_one_comes_home(self):
        left, right = self.lines["world"]
        self.assertEqual({(p[0], p[1]) for p in left}, {(0.0, 0.0)})
        for a, b in zip(right, right[1:]):
            self.assertGreater(b[2], a[2])
            self.assertLessEqual(self.interval(a, b), 1e-9)
        self.assertAlmostEqual(right[0][0], 10.0, places=5)
        self.assertAlmostEqual(right[-1][0], 10.0, places=5)
        self.assertAlmostEqual(max(p[0] for p in right), 31.48, places=2)

    def test_the_closed_geodesic_is_null_and_the_first_pair_of_ends_light_can_join(self):
        (start, end), = self.lines["ergo"]
        self.assertEqual(start[:2], [0.0, 0.0])
        self.assertAlmostEqual(self.interval(start, end), 0.0, delta=1e-6)
        self.assertAlmostEqual(start[2], 41.807, places=3)
        # Each dotted line joins the mouths at one proper time, the left mouth's being its T:
        # spacelike before the closed null geodesic and timelike after.
        pairs = self.lines["axis"]
        self.assertEqual([left[2] for left, _ in pairs], [9.0 * k for k in range(len(pairs))])
        for left, right in pairs:
            self.assertEqual(self.interval(left, right) < 0, left[2] > start[2], left)

    def test_the_closed_curve_is_timelike_and_inside_the_horizon(self):
        (start, end), = self.lines["ctc"]
        self.assertLess(self.interval(start, end), 0)
        self.assertGreater(end[2], start[2])
        (apex, _), = self.lines["ergo"]
        self.assertLess(self.interval(apex, start), 0)
        # Its far end is the right mouth at the proper time it left the left mouth at, after the
        # trip, when the right mouth's T runs ahead of its proper time by the time shift.
        self.assertAlmostEqual(end[2] - start[2], 21.718, places=3)
        self.assertAlmostEqual(end[0], 10.0, places=5)
        # Every generator drawn of the horizon is null, from the event the geodesic leaves.
        generators = [line for line in self.lines["horizon"] if len(line) == 2]
        self.assertEqual(len(generators), 12)
        for a, b in generators:
            self.assertEqual(a, apex)
            self.assertAlmostEqual(self.interval(a, b), 0.0, delta=1e-4)


if __name__ == "__main__":
    unittest.main()
