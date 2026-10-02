"""Petrov's homogeneous vacuum, the one vacuum field whose four symmetries reach every event once:
python3 -m unittest discover -s _tools

What the texts and drawings of `petrov_homogeneous` state, held on the published files. The
published components of Petrov's chart are held to the determinant -1, to light cones on the
plane of t and phi that have turned through psi/2 at the four radii the flat views draw, to the
turning of light off that plane, toward smaller r when launched toward +phi and toward larger r
when launched toward -phi, and to the fourth Killing vector of Gibbons and Gielen. The chart
outside the dust is held to the numbers its captions state, to the null circle on the surface
and to the band of closed timelike curves out to r = e^(pi/sqrt 3) R. The figure's closed
timelike curve is held to closing, to being timelike and future directed at every step by the
published metric, and to the speed its caption states, from the numbers written in the diagram
file. The embedding diagram's surface is held to the circles' radius ell e^(-r/ell) and to lying
level at r = 0. The tests of the published components need sympy and are skipped where it is
absent.
"""
import importlib.util
import json
import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "MFS" / "assets" / "data"
HAS_SYMPY = importlib.util.find_spec("sympy") is not None
sys.path.insert(0, str(Path(__file__).resolve().parent / "derivations"))

ROOT3 = math.sqrt(3)
JOINED = math.exp(-0.5)             # ell = R/sqrt(e), at R = 1
DRAWN = {"petrov": {"ell": 1.0}, "cylinder": {"R": 1.0, "ell": JOINED}}
# The radii of the four flat views of Petrov's chart and the angle their cones have turned through.
PLANES = {"upright": (0.0, 0.0), "diagonal": (math.pi / (2 * ROOT3), math.pi / 4),
          "sideways": (math.pi / ROOT3, math.pi / 2), "inverted": (2 * math.pi / ROOT3, math.pi)}


def published(system):
    metric = json.loads((DATA / "metrics" / "petrov_homogeneous.json").read_text(encoding="utf-8"))
    return next(c for c in metric["coordinates"] if c["id"] == system)


@unittest.skipUnless(HAS_SYMPY, "sympy is not installed")
class Components(unittest.TestCase):
    """The published metric and Christoffel symbols of each chart as numbers."""

    @classmethod
    def setUpClass(cls):
        import sympy as sp
        import verify_metrics as vm
        cls.sp = sp
        cls.charts = {}
        for system, values in DRAWN.items():
            entry = published(system)
            reader = vm.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]], ())
            at = {reader.c: 1}
            at.update({s: values[s.name] for s in reader.allowed if s.name in values})
            cls.charts[system] = (entry, reader, at)

    def matrix(self, system, field, point):
        entry, reader, at = self.charts[system]
        there = {tuple(e["indices"]): e["value"] for e in entry[field]}
        coords = entry["coords"]
        where = {**at, **{reader.symbol[c]: v for c, v in point.items()}}
        return [[float(reader(there.get((a, b), "0")).subs(where)) for b in coords] for a in coords]

    def christoffel(self, system, upper, point):
        entry, reader, at = self.charts[system]
        where = {**at, **{reader.symbol[c]: v for c, v in point.items()}}
        return {tuple(e["indices"][1:]): float(reader(e["value"]).subs(where))
                for e in entry["christoffel"]["variants"]["ull"]["nonzero"] if e["indices"][0] == upper}

    def test_the_determinant_of_petrovs_chart_is_minus_one(self):
        sp = self.sp
        for r in (-2.0, 0.0, 0.9, 3.6):
            g = sp.Matrix(self.matrix("petrov", "metric_components", {"t": 0.3, "r": r, "\\phi": -0.7, "z": 1.1}))
            self.assertAlmostEqual(float(g.det()), -1.0, places=9)

    def test_the_cones_have_turned_through_half_the_phase_on_each_plane_drawn(self):
        """On the plane of t and phi the null directions make the angles pi/4 - a and 3 pi/4 - a
        with the phi axis, and cos(a) d_t + sin(a) d_phi is timelike of norm -e^r, a = sqrt(3) r/2."""
        for view, (r, a) in PLANES.items():
            self.assertAlmostEqual(ROOT3 * r / 2, a, places=12, msg=view)
            g = self.matrix("petrov", "metric_components", {"t": 0.0, "r": r, "\\phi": 0.0, "z": 0.0})
            block = [[g[0][0], g[0][2]], [g[2][0], g[2][2]]]

            def norm(dt, dphi):
                return block[0][0] * dt * dt + 2 * block[0][1] * dt * dphi + block[1][1] * dphi * dphi
            for angle in (math.pi / 4 - a, 3 * math.pi / 4 - a):
                self.assertAlmostEqual(norm(math.sin(angle), math.cos(angle)) / math.exp(r), 0.0, places=12, msg=view)
            self.assertAlmostEqual(norm(math.cos(a), math.sin(a)) / math.exp(r), -1.0, places=12, msg=view)

    def test_light_on_a_plane_of_t_and_phi_is_turned_off_it(self):
        """-Gamma^r_ab k^a k^b for the two null directions: negative for the one that moves toward
        +phi where the cones stand upright, positive for the other, on every plane."""
        for view, (r, a) in PLANES.items():
            gamma = self.christoffel("petrov", "r", {"t": 0.0, "r": r, "\\phi": 0.0, "z": 0.0})
            found = []
            for angle in (math.pi / 4 - a, 3 * math.pi / 4 - a):
                k = {"t": math.sin(angle), "\\phi": math.cos(angle)}
                found.append(-sum(gamma.get((m, n), 0.0) * k[m] * k[n] for m in k for n in k) / math.exp(r))
            self.assertAlmostEqual(found[0], -ROOT3 / 2, places=12, msg=view)
            self.assertAlmostEqual(found[1], ROOT3 / 2, places=12, msg=view)

    def test_the_fourth_killing_vector_is_gibbons_and_gielens(self):
        """R = d_r + (z/ell) d_z + (sqrt 3 x^0 - phi)/(2 ell) d_phi - (x^0 + sqrt 3 phi)/(2 ell) d_0."""
        sp = self.sp
        entry, reader, at = self.charts["petrov"]
        X = [reader.symbol[c] for c in entry["coords"]]
        t, r, phi, z = X
        there = {tuple(e["indices"]): e["value"] for e in entry["metric_components"]}
        g = sp.Matrix(4, 4, lambda i, j: reader(there.get((entry["coords"][i], entry["coords"][j]), "0")).subs(at))
        xi = [-(t + sp.sqrt(3) * phi) / 2, 1, (sp.sqrt(3) * t - phi) / 2, z]
        lie = sp.Matrix(4, 4, lambda a, b: sum(xi[c] * sp.diff(g[a, b], X[c]) + g[c, b] * sp.diff(xi[c], X[a])
                                               + g[a, c] * sp.diff(xi[c], X[b]) for c in range(4)))
        self.assertEqual(lie.applyfunc(sp.simplify), sp.zeros(4, 4))

    def test_the_cylinders_captions_state_the_published_numbers(self):
        """F, M and L at r = R, 3R and 12R, and the radii where F and L change sign."""
        entry, reader, at = self.charts["cylinder"]

        def value(name, r):
            return float(reader.parameters[name].subs({**at, reader.symbol["r"]: r}))
        self.assertAlmostEqual(value("F", 1.0), 1.0, places=12)
        self.assertAlmostEqual(value("M", 1.0), 1.0, places=12)
        self.assertAlmostEqual(value("L", 1.0), 0.0, places=12)
        self.assertAlmostEqual(value("F", 3.0), -2.62, places=2)
        self.assertAlmostEqual(value("L", 3.0), -3.27, places=2)
        self.assertAlmostEqual(value("F", 12.0), 1.59, places=2)
        self.assertAlmostEqual(value("L", 12.0), 12.7, places=1)
        self.assertAlmostEqual(math.degrees(ROOT3 * math.log(12.0) / 2), 123, delta=0.5)
        edge = math.exp(math.pi / ROOT3)
        self.assertAlmostEqual(edge, 6.13, places=2)
        self.assertAlmostEqual(math.exp(2 * math.pi / ROOT3), 37.6, places=1)
        self.assertAlmostEqual(value("L", edge), 0.0, places=9)
        self.assertLess(value("L", edge - 0.01), 0.0)
        self.assertGreater(value("L", edge + 0.01), 0.0)
        for r in (1.001, 2.0, 4.0, 6.0):
            self.assertLess(value("L", r), 0.0)
        # FL + M^2 = r^2 on the plane of t and phi.
        for r in (1.0, 3.0, 12.0, 40.0):
            self.assertAlmostEqual(value("F", r) * value("L", r) + value("M", r) ** 2, r * r, places=9)


class Figure(unittest.TestCase):
    """The figure of turning cones and its closed timelike curve, from the numbers written."""

    @classmethod
    def setUpClass(cls):
        diagrams = json.loads((DATA / "diagrams" / "petrov_homogeneous.json").read_text(encoding="utf-8"))
        cls.figure = next(f for f in diagrams["projections"]["petrov"] if f["id"] == "turning")

    @staticmethod
    def metric(r):
        """Petrov's metric on the slice z = 0 in the drawing's (X, Y, T) = (r, phi, ct), at ell = 1."""
        c, s, e = math.cos(ROOT3 * r), math.sin(ROOT3 * r), math.exp(r)
        return {("T", "T"): -e * c, ("Y", "Y"): e * c, ("T", "Y"): -e * s, ("X", "X"): 1.0}

    def dot(self, r, u, v):
        g = self.metric(r)
        (ux, uy, ut), (vx, vy, vt) = u, v
        return (g[("X", "X")] * ux * vx + g[("Y", "Y")] * uy * vy + g[("T", "T")] * ut * vt
                + g[("T", "Y")] * (ut * vy + uy * vt))

    def test_the_curve_closes_and_is_timelike_and_future_directed(self):
        legs = [line["points"] for line in self.figure["turn"]["lines"] if line["class"] == "ctc"]
        self.assertEqual(len(legs), 4)
        for leg, following in zip(legs, legs[1:] + legs[:1]):
            for a, b in zip(leg[-1], following[0]):
                self.assertAlmostEqual(a, b, places=5)
        speeds = []
        for leg in legs:
            for p, q in zip(leg, leg[1:]):
                k = [q[i] - p[i] for i in range(3)]
                r = (p[0] + q[0]) / 2
                a = ROOT3 * r / 2
                future = (0.0, math.sin(a), math.cos(a))
                kk, kv, vv = self.dot(r, k, k), self.dot(r, k, future), self.dot(r, future, future)
                self.assertLess(kk, 0.0)
                self.assertLess(kv, 0.0)
                gamma = -kv / math.sqrt(kk * vv)
                speeds.append(math.sqrt(max(0.0, 1 - 1 / gamma ** 2)))
        self.assertAlmostEqual(max(speeds), 0.8, places=2)
        self.assertAlmostEqual(min(speeds), 0.8, places=2)

    def test_the_cones_are_null_and_turn_through_half_the_phase(self):
        cones = self.figure["turn"]["cones"]
        self.assertEqual(len(cones), 5)
        radii = sorted(cone["apex"][0] for cone in cones)
        for found, wanted in zip(radii, (-math.pi / ROOT3, -math.pi / (2 * ROOT3), 0.0, math.pi / (2 * ROOT3),
                                         math.pi / ROOT3)):
            self.assertAlmostEqual(found, wanted, places=5)
        for cone in cones:
            r = cone["apex"][0]
            a = ROOT3 * r / 2
            axis = [0.0, 0.0, 0.0]
            for point in cone["rim"]:
                k = [point[i] - cone["apex"][i] for i in range(3)]
                size = sum(x * x for x in k)
                self.assertAlmostEqual(self.dot(r, k, k) / (size * math.exp(abs(r))), 0.0, places=4)
                self.assertLess(self.dot(r, k, (0.0, math.sin(a), math.cos(a))), 0.0)
                axis = [axis[i] + k[i] for i in range(3)]
            # The rim's centre lies along cos(a) d_t + sin(a) d_phi from the apex.
            self.assertAlmostEqual(math.atan2(axis[1], axis[2]), a, places=3)


class Embedding(unittest.TestCase):
    """The hyperbolic plane of r and z, from the numbers written."""

    def test_the_circles_have_the_radius_of_the_horocycles_and_the_surface_lies_level_at_zero(self):
        surface = json.loads((DATA / "embedding" / "petrov_homogeneous.json").read_text(encoding="utf-8"))
        pieces = {p["id"]: p for p in surface["views"][0]["surfaces"][0]["pieces"]}
        self.assertEqual(set(pieces), {"minkowski", "pseudosphere"})
        for piece in pieces.values():
            for r, rho, _ in piece["points"]:
                self.assertAlmostEqual(rho, math.exp(-r), places=5)
        horn = pieces["pseudosphere"]["points"]
        self.assertAlmostEqual(horn[0][0], 0.0, places=9)
        # The tractrix: Z = arcosh(e^r) - sqrt(1 - e^(-2r)) above r = 0.
        for r, _, z in horn:
            self.assertAlmostEqual(z - horn[0][2], math.acosh(math.exp(r)) - math.sqrt(1 - math.exp(-2 * r)), places=4)
        sheet = pieces["minkowski"]["points"]
        self.assertAlmostEqual(sheet[-1][0], 0.0, places=9)
        self.assertAlmostEqual(sheet[0][0], -math.log(3.0), places=5)
        for r, _, z in sheet:
            root = math.sqrt(max(0.0, math.exp(-2 * r) - 1))
            self.assertAlmostEqual(z - sheet[-1][2], math.atan(root) - root, places=4)


if __name__ == "__main__":
    unittest.main()
