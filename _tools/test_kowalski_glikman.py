"""The eleven dimensional supergravity flux, using the published tensors and connection."""
import itertools
import json
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


if __name__ == "__main__":
    unittest.main()
