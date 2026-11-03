"""The printed magnetic string field, current and Gaussian Einstein-Maxwell equations."""
import json
import re
import sys
import unittest
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "_tools" / "derivations"))
import verify_metrics as vm


class MagneticExterior(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        page = json.loads((ROOT / "MFS/assets/data/metrics/superconducting_string.json").read_text())
        cls.page = page
        cls.chart, = page["coordinates"]
        cls.reader = vm.Reader(cls.chart["coords"], [p["symbol"] for p in cls.chart["parameters"]], (),
                               held=vm.HELD[("superconducting_string", "cylindrical")],
                               rates=vm.RATES[("superconducting_string", "cylindrical")])
        cls.r = cls.reader.symbol["r"]

    @classmethod
    def written(cls, value):
        return sp.sympify(value).subs(cls.reader.held).doit()

    def printed_expression(self, pattern):
        match = re.search(pattern, self.page["history"])
        self.assertIsNotNone(match, pattern)
        text = match.group(1).replace(r"\pm ", "").replace("x", r"(r/\ell)")
        return self.written(self.reader(text))

    def test_metric_field_current_and_maxwell_flux(self):
        """Compare the published Ricci tensor to Gaussian stress, using the printed field and current."""
        self.assertIn("Gaussian electromagnetic units", self.page["convention"])
        field = self.printed_expression(r"\$F_\{zr\}=([^$]+)\$")
        current = self.printed_expression(r"\$I=([^$]+)\$")
        entry = self.chart
        names = entry["coords"]
        g = sp.zeros(4)
        for item in entry["metric_components"]:
            i, j = (names.index(c) for c in item["indices"])
            g[i,j] = self.written(self.reader(item["value"]))
        ricci = sp.zeros(4)
        for item in entry["ricci_tensor"]["variants"]["ll"]["nonzero"]:
            i, j = (names.index(c) for c in item["indices"])
            ricci[i,j] = self.written(self.reader(item["value"]))
        inv = sp.diag(*(1/g[i,i] for i in range(4)))
        em = sp.zeros(4)
        em[3,1], em[1,3] = field, -field
        square = sum(inv[i,i]*inv[j,j]*em[i,j]**2 for i in range(4) for j in range(4))
        source = sp.Matrix(4, 4, lambda i,j: 2*(sum(em[i,a]*inv[a,a]*em[j,a] for a in range(4))
                                                    - g[i,j]*square/4))
        p = self.reader.parameters
        for rv, lv, mv, kv, wv in ((1.5, 1.0, 0.2, 0.7, 0.8),
                                   (3.2, 1.3, -0.3, 1.4, 0.9),
                                   (2.0, 0.8, 0.4, 0.2, 0.7)):
            at = {self.r:rv, p["ell"]:lv, p["m"]:mv, p["k"]:kv, p["W"]:wv}
            for i in range(4):
                for j in range(4):
                    self.assertAlmostEqual(float((ricci[i,j]-source[i,j]).subs(at).evalf()), 0, delta=1e-9)
            # sqrt(-g) F^{rz} = -2I and its radial derivative vanish in the exterior.
            flux = -p["W"]*self.r*self.written(p["F"])**2*field
            self.assertAlmostEqual(float((flux+2*current).subs(at).evalf()), 0, delta=1e-10)
            self.assertAlmostEqual(float(sp.diff(flux,self.r).subs(at).evalf()), 0, delta=1e-10)
            self.assertGreater(abs(float(source[0,0].subs(at).evalf())), 1e-5)


if __name__ == "__main__":
    unittest.main()
