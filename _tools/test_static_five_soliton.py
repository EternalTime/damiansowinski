"""The scalar reduction, bubble regularity, and sphere areas of the neutral static family."""
import math
import sys
import unittest
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(Path(__file__).resolve().parent / "derivations"))
import null_rays as nr


class StaticFiveSoliton(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.metric, cls.entry, cls.reader = nr.load("static_five_soliton", "spherical")
        cls.g = nr.published_matrix(cls.reader, cls.entry, "metric_components")

    def test_the_full_family_and_three_sourced_charts_are_retained(self):
        self.assertEqual(self.metric["signature"], "(-,+,+,+,+)")
        self.assertEqual([c["id"] for c in self.metric["coordinates"]],
                         ["spherical", "isotropic", "lake"])
        for chart in self.metric["coordinates"]:
            self.assertIn("a^2+ab+b^2=1", chart["domains"])
            self.assertTrue({"a", "b", "r_s"}.issubset(p["symbol"] for p in chart["parameters"]))

    def test_reduction_in_einstein_units_is_fishers_massless_scalar(self):
        R = self.reader
        _, entry, fisher = nr.load("fisher_jnw", "spherical")
        fg = nr.published_matrix(fisher, entry, "metric_components")
        for aa, bb in ((sp.Rational(8,7),-sp.Rational(3,7)), (0,1), (1,0)):
            aa,bb=sp.sympify(aa),sp.sympify(bb)
            gamma = aa + bb/2
            self.assertEqual(sp.expand(gamma**2+3*bb**2/4), 1)
            for r in (sp.Rational(5,4),2,5):
                own={R.parameters["a"]:aa,R.parameters["b"]:bb,R.parameters["r_s"]:1,
                     R.symbol["r"]:r,R.symbol[r"\theta"]:sp.Rational(3,4)}
                there={fisher.parameters["b"]:1,fisher.parameters["gamma"]:gamma,
                       fisher.symbol["r"]:r,fisher.symbol[r"\theta"]:sp.Rational(3,4)}
                for i in range(4):
                    reduced=(self.g[i,i]*self.g[4,4]**sp.Rational(1,2)).subs(own)
                    self.assertLess(abs(complex(sp.N(reduced-fg[i,i].subs(there),35))),1e-25)

    def test_four_pi_period_closes_the_bubble_smoothly(self):
        R=self.reader
        for scale in (sp.Rational(3,4),sp.Rational(5,3)):
            eps=sp.Rational(1,10**9)
            at={R.parameters["a"]:0,R.parameters["b"]:1,R.parameters["r_s"]:scale,
                R.symbol["r"]:scale*(1+eps)}
            circle=4*sp.pi*scale*sp.sqrt(self.g[4,4].subs(at))
            proper=2*scale*sp.sqrt(eps)
            self.assertAlmostEqual(float(circle/proper),2*math.pi,places=7)
            christoffel=self.entry["christoffel"]["variants"]["ull"]["nonzero"]
            radial=next(q for q in christoffel if q["indices"]==["r","t","t"])
            self.assertEqual(R(radial["value"]).subs(at),0)

    def test_sphere_area_can_turn_before_the_singular_endpoint(self):
        R=self.reader
        for aa,bb,turns in ((sp.Rational(8,7),-sp.Rational(3,7),False),
                            (1/sp.sqrt(3),1/sp.sqrt(3),True)):
            area=self.g[2,2].subs({R.parameters["a"]:aa,R.parameters["b"]:bb,R.parameters["r_s"]:1})
            derivative=sp.diff(area,R.symbol["r"])
            if turns:
                minimum=(1+aa+bb)/2
                self.assertLess(abs(float(derivative.subs(R.symbol["r"],minimum))),1e-12)
                self.assertLess(float(derivative.subs(R.symbol["r"],1.01)),0)
            else:
                self.assertGreater(float(derivative.subs(R.symbol["r"],1.01)),0)
            self.assertGreater(float(derivative.subs(R.symbol["r"],2)),0)

    def test_history_year_uses_the_original_kramer_scan(self):
        bib=(ROOT/"assets/data/references.bib").read_text()
        block=bib.split("@article{kramer1971projective,",1)[1].split("\n}",1)[0]
        self.assertIn("year    = {1971}",block)
        self.assertIn("pages   = {807--811}",block)
        self.assertIn("1971",self.metric["history"])


if __name__ == "__main__":
    unittest.main()
