"""The Poincare dodecahedral universe: python3 -m unittest discover -s _tools

What its texts state, held on the binary icosahedral group built from Aurich, Lustig and Steiner's
two generators, gamma_1 = j and gamma_2 = sigma/2 + i/(2 sigma) + j/2 with sigma the golden ratio,
acting by left multiplication on the unit quaternions w + x i + y j + z ij, their (18). The group
has 120 elements and no element but the identity fixes a point; every element moves every point the
same distance, its Clifford translation, the shortest being pi/5 and taken by twelve elements, one
for each face of the cell; the cell about an observer is a dodecahedron of inradius pi/10 whose twenty
vertices lie 0.388 of the radius of curvature away, Gott's figure; the space vibrates at the wave
numbers 1, 13, 21, 25, 31, ... and never at an even one, Aurich, Lustig and Steiner's (26); and a shift
of the toroidal chart's alpha and gamma together is the Clifford translation along w + z ij. The
published charts carry the identification, and the drawings hold what their captions say. Everything
here needs nothing but the files and the standard library.
"""
import itertools
import json
import math
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "MFS" / "assets" / "data"
SIGMA = (1 + math.sqrt(5)) / 2


def mul(p, q):
    """The quaternion product in the basis (1, i, j, ij), with i^2 = j^2 = -1 and ij = -ji."""
    a1, b1, c1, d1 = p
    a2, b2, c2, d2 = q
    return (a1 * a2 - b1 * b2 - c1 * c2 - d1 * d2, a1 * b2 + b1 * a2 + c1 * d2 - d1 * c2,
            a1 * c2 - b1 * d2 + c1 * a2 + d1 * b2, a1 * d2 + b1 * c2 - c1 * b2 + d1 * a2)


def key(q):
    return tuple(round(x, 9) + 0.0 for x in q)


def group():
    """The binary icosahedral group, the closure of Aurich, Lustig and Steiner's two generators."""
    generators = ((0.0, 0.0, 1.0, 0.0), (SIGMA / 2, 1 / (2 * SIGMA), 0.5, 0.0))
    found = {key((1.0, 0.0, 0.0, 0.0)): (1.0, 0.0, 0.0, 0.0)}
    frontier = list(found.values())
    while frontier:
        fresh = []
        for q in frontier:
            for g in generators:
                r = mul(g, q)
                if key(r) not in found:
                    found[key(r)] = r
                    fresh.append(r)
        frontier = fresh
    return list(found.values())


def distance(p, q):
    """The distance on the unit 3-sphere, cos d = p.q."""
    return math.acos(max(-1.0, min(1.0, sum(a * b for a, b in zip(p, q)))))


def null_vector(rows):
    """A vector of R^4 orthogonal to three rows, the generalised cross product."""
    out = []
    for k in range(4):
        cols = [c for c in range(4) if c != k]
        m = [[row[c] for c in cols] for row in rows]
        det = (m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
               + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]))
        out.append((-1) ** k * det)
    return out


G = group()
ONE = (1.0, 0.0, 0.0, 0.0)
NEIGHBOURS = [g for g in G if abs(distance(ONE, g) - math.pi / 5) < 1e-9]


def vertices():
    """The corners of the Dirichlet cell about 1: points equidistant from 1 and three of its nearest
    images that no image is nearer to."""
    found = {}
    for a, b, c in itertools.combinations(NEIGHBOURS, 3):
        rows = [[(1.0 if i == 0 else 0.0) - g[i] for i in range(4)] for g in (a, b, c)]
        v = null_vector(rows)
        size = math.sqrt(sum(x * x for x in v))
        if size < 1e-9:
            continue
        v = [x / size for x in v]
        if v[0] < 0:
            v = [-x for x in v]
        if all(sum(x * y for x, y in zip(v, g)) <= v[0] + 1e-9 for g in G):
            found[key(v)] = v
    return list(found.values())


def multiplicity(beta):
    """How many eigenfunctions of the Laplacian with eigenvalue beta^2 - 1 the quotient keeps: beta
    times the invariants of I* in the irreducible representation of SU(2) of dimension beta."""
    total = 0.0
    for g in G:
        theta = math.acos(max(-1.0, min(1.0, g[0])))
        if abs(math.sin(theta)) < 1e-12:
            total += beta * (1 if g[0] > 0 else (-1) ** (beta - 1))
        else:
            total += math.sin(beta * theta) / math.sin(theta)
    return round(beta * total / len(G))


class Group(unittest.TestCase):
    def test_order_120(self):
        self.assertEqual(len(G), 120)

    def test_no_fixed_point(self):
        # g q = q for a unit quaternion q only when g = 1.
        self.assertEqual(sum(1 for g in G if abs(g[0] - 1) < 1e-9), 1)

    def test_every_point_moved_the_same_distance(self):
        for g in G[:20]:
            for q in ((0.5, 0.5, 0.5, 0.5), (0.0, 0.6, 0.0, 0.8), (math.cos(1.0), 0.0, math.sin(1.0), 0.0)):
                self.assertAlmostEqual(distance(q, mul(g, q)), distance(ONE, g), places=12)

    def test_twelve_nearest_images_at_pi_over_5(self):
        self.assertAlmostEqual(min(distance(ONE, g) for g in G if g != ONE), math.pi / 5, places=12)
        self.assertEqual(len(NEIGHBOURS), 12)

    def test_minus_one_is_in_the_group(self):
        self.assertIn(key((-1.0, 0.0, 0.0, 0.0)), {key(g) for g in G})


class Cell(unittest.TestCase):
    def test_inradius_is_pi_over_10(self):
        # The face between 1 and a nearest image bisects the arc between them.
        self.assertAlmostEqual(distance(ONE, NEIGHBOURS[0]) / 2, math.pi / 10, places=12)

    def test_twenty_vertices_at_gotts_radius(self):
        corners = vertices()
        self.assertEqual(len(corners), 20)
        radii = {round(distance(ONE, v), 9) for v in corners}
        self.assertEqual(len(radii), 1)
        self.assertAlmostEqual(radii.pop(), 0.388, places=3)

    def test_volume_is_a_hundred_and_twentieth(self):
        # Gott and Luminet: 120 cells of one volume fill the 3-sphere of volume 2 pi^2.
        self.assertAlmostEqual(2 * math.pi ** 2 / len(G), math.pi ** 2 / 60, places=12)


class Spectrum(unittest.TestCase):
    def test_the_wave_numbers_of_aurich_lustig_and_steiner(self):
        theirs = [1, 13, 21, 25, 31, 33, 37, 41, 43, 45, 49, 51, 53, 55, 57] + list(range(61, 100, 2))
        self.assertEqual([b for b in range(1, 100) if multiplicity(b) > 0], theirs)

    def test_none_even(self):
        self.assertTrue(all(multiplicity(b) == 0 for b in range(2, 100, 2)))

    def test_the_three_lowest_hold_59_eigenfunctions(self):
        self.assertEqual(sum(multiplicity(b) for b in (13, 21, 25)), 59)


class Toroidal(unittest.TestCase):
    def test_shift_of_alpha_and_gamma_is_a_clifford_translation(self):
        for v, alpha, gamma, s in ((0.1, 0.3, 1.2, math.pi / 5), (0.4, 2.0, -0.7, 0.9), (0.0, 1.0, 0.0, 2.5)):
            point = (math.sqrt(1 - 2 * v) * math.cos(alpha), math.sqrt(2 * v) * math.cos(gamma),
                     math.sqrt(2 * v) * math.sin(gamma), math.sqrt(1 - 2 * v) * math.sin(alpha))
            moved = mul((math.cos(s), 0.0, 0.0, math.sin(s)), point)
            shifted = (math.sqrt(1 - 2 * v) * math.cos(alpha + s), math.sqrt(2 * v) * math.cos(gamma + s),
                       math.sqrt(2 * v) * math.sin(gamma + s), math.sqrt(1 - 2 * v) * math.sin(alpha + s))
            for a, b in zip(moved, shifted):
                self.assertAlmostEqual(a, b, places=12)

    def test_the_group_turned_holds_that_translation(self):
        # Every element of real part cos(pi/5) is conjugate to cos(pi/5) + sin(pi/5) ij, by the unit
        # quaternion that turns its axis onto ij, and conjugation is a rotation of the 3-sphere.
        g = NEIGHBOURS[0]
        axis = [x / math.sin(math.pi / 5) for x in g[1:]]
        target = (0.0, 0.0, 1.0)
        cross = (axis[1] * target[2] - axis[2] * target[1], axis[2] * target[0] - axis[0] * target[2],
                 axis[0] * target[1] - axis[1] * target[0])
        half = math.acos(max(-1.0, min(1.0, sum(a * b for a, b in zip(axis, target))))) / 2
        size = math.sqrt(sum(x * x for x in cross)) or 1.0
        u = (math.cos(half), *(math.sin(half) * x / size for x in cross))
        conjugate = mul(mul(u, g), (u[0], -u[1], -u[2], -u[3]))
        for a, b in zip(conjugate, (math.cos(math.pi / 5), 0.0, 0.0, math.sin(math.pi / 5))):
            self.assertAlmostEqual(a, b, places=12)


class Published(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.metric = json.loads((DATA / "metrics" / "poincare_dodecahedral.json").read_text(encoding="utf-8"))

    def test_every_chart_states_the_identification(self):
        for chart in self.metric["coordinates"]:
            self.assertTrue(any("I^*" in d for d in chart["domains"]), chart["id"])

    def test_the_toroidal_chart_states_the_shift(self):
        chart = next(c for c in self.metric["coordinates"] if c["id"] == "toroidal")
        self.assertTrue(any("\\alpha + \\pi/5, \\gamma + \\pi/5" in d for d in chart["domains"]))

    def test_the_comoving_chart_is_friedmanns_closed_universe(self):
        chart = next(c for c in self.metric["coordinates"] if c["id"] == "comoving")
        self.assertIn("a^2\\left(d\\chi^2 + \\sin^2\\chi", chart["line_element"])


class Sky(unittest.TestCase):
    """Luminet's universe, Omega_m = 0.28 and Omega_0 = 1.013, integrated here by Simpson's rule in
    the scale factor, against the numbers the conformal drawing and its caption state."""

    @staticmethod
    def conformal_time(a_from, a_to, steps=20000):
        om, o0 = 0.28, 1.013
        rate = 1 / math.sqrt(o0 - 1)

        def f(s):
            # a = s^2 takes the integrable 1/sqrt(a) at the bang away.
            a = s * s
            return 2 * s / (a * a * rate * math.sqrt(om / a ** 3 + (o0 - om) - (o0 - 1) / a ** 2))
        lo, hi = math.sqrt(a_from), math.sqrt(a_to)
        h = (hi - lo) / steps
        total = f(lo + 1e-15) + f(hi)
        total += sum((4 if k % 2 else 2) * f(lo + k * h) for k in range(1, steps))
        return total * h / 3

    def test_today_and_the_last_scattering(self):
        self.assertAlmostEqual(self.conformal_time(0, 1), 0.38882, places=4)
        self.assertAlmostEqual(self.conformal_time(0, 1 / 1101), 0.012988, places=5)

    def test_the_matched_circles_are_about_35_degrees(self):
        chi = self.conformal_time(0, 1) - self.conformal_time(0, 1 / 1101)
        self.assertAlmostEqual(chi, 0.3758, places=3)
        self.assertGreater(chi, math.pi / 10)
        alpha = math.degrees(math.acos(math.tan(math.pi / 10) / math.tan(chi)))
        self.assertAlmostEqual(alpha, 34.6, places=1)


if __name__ == "__main__":
    unittest.main()
