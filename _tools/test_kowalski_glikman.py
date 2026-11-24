"""The eleven dimensional supergravity flux, using the published tensors and connection, and the
embedding diagram's ring of free particles, held to the two rates of the published g_uu from the
numbers in its file."""
import itertools
import json
import math
from pathlib import Path
import sys
import unittest

import numpy as np
import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parent / "derivations"))
import verify_metrics as vm

ROOT = Path(__file__).resolve().parent.parent
ID = "kowalski_glikman"


class FourForm(unittest.TestCase):
    def test_parallel_null_flux_and_eleven_dimensional_einstein_equation(self):
        metric = json.loads((ROOT / "MFS/assets/data/metrics" / f"{ID}.json").read_text())
        for chart in metric["coordinates"]:
            names = chart["coords"]
            self.assertEqual(len(names), 11)
            reader = vm.Reader(names, [r"\mu"], ())
            mu, u = reader.parameters["mu"], reader.symbol["u"]
            rosen = chart["id"] == "rosen"
            coefficient = mu * (sp.sin(mu*u/3)**3 if rosen else 1)
            slots = []
            for order in itertools.permutations(range(4)):
                indices = tuple((0, 2, 3, 4)[i] for i in order)
                sign = (-1)**sum(order[i] > order[j] for i in range(4) for j in range(i+1, 4))
                slots.append((indices, sign))
            for number in (1, 2, 3):
                at = {reader.c: 1, mu: (number+2)/4,
                      **{s: (i+number+1)/7 for i, s in enumerate(reader.symbol.values())}}
                def evaluate(value):
                    return float(reader(value).subs(at).evalf())
                inverse = np.zeros((11, 11))
                ricci = np.zeros((11, 11))
                gamma = np.zeros((11, 11, 11))
                for c in chart["inverse_metric_components"]:
                    inverse[tuple(names.index(i) for i in c["indices"])] = evaluate(c["value"])
                for c in chart["ricci_tensor"]["variants"]["ll"]["nonzero"]:
                    ricci[tuple(names.index(i) for i in c["indices"])] = evaluate(c["value"])
                for c in chart["christoffel"]["variants"]["ull"]["nonzero"]:
                    gamma[tuple(names.index(i) for i in c["indices"])] = evaluate(c["value"])
                form = np.zeros((11,)*4)
                derivative = np.zeros((11,)*5)
                for indices, sign in slots:
                    form[indices] = sign*float(coefficient.subs(at))
                    derivative[(0,)+indices] = sign*float(sp.diff(coefficient, u).subs(at))
                # The connection acts on each covariant slot independently.
                nabla = derivative - np.einsum('fab,fcde->abcde', gamma, form)
                nabla -= np.einsum('fac,bfde->abcde', gamma, form)
                nabla -= np.einsum('fad,bcfe->abcde', gamma, form)
                nabla -= np.einsum('fae,bcdf->abcde', gamma, form)
                self.assertLess(np.max(np.abs(nabla)), 1e-11, chart["id"])
                raised3 = np.einsum('pi,qj,rk,aijk->apqr', inverse, inverse, inverse, form, optimize=True)
                source = np.einsum('apqr,bpqr->ab', form, raised3)/12
                self.assertLess(np.max(np.abs(source-ricci)), 1e-11)
                # A doubled flux would quadruple the source and fail this normalization.
                self.assertGreater(np.max(np.abs(4*source-ricci)), 0.1)
                squared = np.einsum('ab,apqr,bpqr', inverse, form, raised3)
                self.assertEqual(float(squared), 0.0)
                # Each nonzero term contains du, hence F wedge F vanishes.
                self.assertTrue(all(0 in indices for indices, _ in slots))
                # Parallel F is closed and co-closed; the Chern-Simons term also vanishes.


class RingOfFreeParticles(unittest.TestCase):
    def test_the_ring_closes_onto_each_axis_in_turn_at_the_two_rates_of_the_wave(self):
        # A particle at rest has d^2x_i/d(cu)^2 = -Gamma^x_i_uu = d(g_uu)/dx_i/2 = -omega_i^2 x_i, with
        # omega_1 = mu/3 on x_1, a direction the flux threads, and omega_4 = mu/6 on x_4, one of the other
        # six, so the ring at rest on x_1^2 + x_4^2 = 1/mu^2 at u = 0 is the ellipse of semi-axes
        # |cos(mu c u/3)| along x_1 and |cos(mu c u/6)| along x_4: a segment of the x_4 axis of half length
        # 1/sqrt2 at mu c u = 3 pi/2 and 9 pi/2, a segment of the x_1 axis of half length 1 at 3 pi, and the
        # circle again at 6 pi.
        metric = json.loads((ROOT / "MFS/assets/data/metrics" / f"{ID}.json").read_text())
        chart = next(c for c in metric["coordinates"] if c["id"] == "brinkmann")
        reader = vm.Reader(chart["coords"], [r"\mu"], ())
        mu = reader.parameters["mu"]
        guu = reader(next(c["value"] for c in chart["metric_components"] if c["indices"] == ["u", "u"]))
        rates = [sp.simplify(-sp.diff(guu, reader.symbol[x], 2) / 2) for x in ("x_1", "x_4")]
        self.assertEqual(rates, [mu ** 2 / 9, mu ** 2 / 36])
        embedding = json.loads((ROOT / "MFS/assets/data/embedding" / f"{ID}.json").read_text())
        views = {view["id"]: view for view in embedding["views"]}
        self.assertEqual(list(views), ["tube", "ring"])
        ring = views["ring"]
        for surface, k in zip(ring["surfaces"], (0, 15, 30, 45, 60)):
            self.assertAlmostEqual(surface["time"], k * math.pi / 10, delta=1e-6)
        frames = ring["movie"]["frames"]
        self.assertEqual(len(frames), 61)
        for k, frame in enumerate(frames):
            u = frame["value"]
            self.assertAlmostEqual(u, k * math.pi / 10, delta=1e-6)
            points = frame["curves"][0]["points"]
            a, b = max(abs(p[0]) for p in points), max(abs(p[1]) for p in points)
            self.assertAlmostEqual(a, abs(math.cos(u / 3)), delta=1e-6, msg=k)
            self.assertAlmostEqual(b, abs(math.cos(u / 6)), delta=1e-6, msg=k)
            self.assertTrue(all(Z == 0 for _, _, Z in points), k)
            if min(a, b) > 1e-3:
                for X, Y, _ in points:
                    self.assertAlmostEqual(math.hypot(X / a, Y / b), 1, delta=1e-5, msg=k)
        for k, half, on in ((15, 1 / math.sqrt(2), 1), (30, 1.0, 0), (45, 1 / math.sqrt(2), 1)):
            points = frames[k]["curves"][0]["points"]
            self.assertLess(max(abs(p[1 - on]) for p in points), 1e-6, k)
            self.assertAlmostEqual(max(abs(p[on]) for p in points), half, delta=1e-6, msg=k)
        for k in (0, 60):
            for X, Y, _ in frames[k]["curves"][0]["points"]:
                self.assertAlmostEqual(math.hypot(X, Y), 1, delta=1e-6, msg=k)
        # The tube stacks the ring at every u, a height of 1/(6 mu) for each 1/mu of cu.
        tube = views["tube"]["surfaces"][0]
        grid = tube["pieces"][0]["grid"]
        for u, a, b, z in zip(grid["u"], grid["a"], grid["b"], grid["z"]):
            self.assertAlmostEqual(a, abs(math.cos(u / 3)), delta=1e-6, msg=u)
            self.assertAlmostEqual(b, abs(math.cos(u / 6)), delta=1e-6, msg=u)
            self.assertAlmostEqual(z, u / 6, delta=1e-6, msg=u)
        self.assertEqual([c["label"] for c in tube["curves"]],
                         ["$\\mu cu = 0$", "$\\mu cu = 3\\pi/2$", "$\\mu cu = 3\\pi$", "$\\mu cu = 9\\pi/2$",
                          "$\\mu cu = 6\\pi$"])


if __name__ == "__main__":
    unittest.main()
