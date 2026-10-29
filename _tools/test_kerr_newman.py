"""The horizons of Kerr's and Kerr and Newman's black holes as their texts print them.

The readers of the published files and MFS_SPACETIME are those of test_build_mfs_data.py."""
import decimal
import re
import unittest
from decimal import ROUND_HALF_EVEN, Decimal
from fractions import Fraction

from test_build_mfs_data import (checks, conformal_file, conformal_prose, diagram_file, diagram_prose, embedding_file,
                                 embedding_prose)


@checks("kerr", "kerr_newman")
class StatedHorizons(unittest.TestCase):
    """Every horizon radius a Kerr or Kerr-Newman text prints is the root of Delta = r^2 - 2GMr/c^2
    + a^2 + r_Q^2 at the spin and charge its diagrams are drawn for, rounded once from that root,
    so that no value rounded to more places first, as 1.6245 for r_+ = 1.62449..., can be rounded
    again the wrong way."""

    PARAMETER = r"\$({}) = (?:\\frac\{{(\d+)\}}\{{(\d+)\}}|(\d+))\$"

    def parameters(self, name):
        """The values every diagram view of `name` is drawn for, which must agree."""
        stated = set()
        for views in diagram_file(name)["systems"].values():
            for view in views:
                values = {m[1]: Fraction(int(m[2]), int(m[3])) if m[2] else Fraction(int(m[4]))
                          for m in re.finditer(self.PARAMETER.format("G|M|a|r_Q"), view["settings"])}
                stated.add(tuple(sorted(values.items())))
        self.assertEqual(len(stated), 1, f"{name}: its diagrams are drawn for different values: {stated}")
        values = dict(stated.pop())
        self.assertEqual((values["G"], values["M"]), (1, 1), name)
        return values["a"], values.get("r_Q", Fraction(0))

    def test_every_printed_horizon_radius_is_the_root_rounded_once(self):
        texts = {"kerr": [], "kerr_newman": []}
        for name in texts:
            texts[name] += diagram_prose(name, diagram_file(name))
            texts[name] += conformal_prose(name, conformal_file(name))
            texts[name] += embedding_prose(name, embedding_file(name))
        found = []
        for name, fields in texts.items():
            a, r_Q = self.parameters(name)
            with decimal.localcontext(prec=50):
                root = (1 - a * a - r_Q * r_Q)
                root = (Decimal(root.numerator) / Decimal(root.denominator)).sqrt()
                radius = {"+": 1 + root, "-": 1 - root}
            for field, text in fields:
                for sign, printed in re.findall(r"r_([+-]) = (\d+\.\d+)", text):
                    places = Decimal(1).scaleb(-len(printed.split(".")[1]))
                    self.assertEqual(printed, str(radius[sign].quantize(places, ROUND_HALF_EVEN)),
                                     f"{field} prints r_{sign} = {printed}, where r_{sign} = {radius[sign]:.12f}")
                    found.append((name, field, sign))
        self.assertTrue({name for name, _, _ in found} == set(texts), f"a spacetime prints no horizon: {found}")
        # Kerr-Newman's r_+ = 1.62449... GM/c^2 is printed in the caption of the principal rays from
        # above and of the light cones on the equator, the two that once printed 1.625.
        for view in ("boyer_lindquist/above", "boyer_lindquist/dragging"):
            self.assertTrue(any(name == "kerr_newman" and f" {view}.caption" in field and sign == "+"
                                for name, field, sign in found), view)


if __name__ == "__main__":
    unittest.main()
