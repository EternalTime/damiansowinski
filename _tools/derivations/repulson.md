# Repulsons and massless black holes: derivation notes

## Sources of the charts

Both charts are Kallosh and Linde's, hep-th/9507022 v3 (Phys. Rev. D 52, 7137), read in their signature (+,-,-,-) and printed here in (-,+,+,+).

- The Einstein metric is their eq. (36), TeX label `II`: $ds^2_{\rm can} = \Delta^{-1/2}dt^2 - \Delta^{1/2}d\vec x^2$ with $\Delta = 1 + 4M/r - 4g^2q^2/r^2$ and $e^{-2\phi_0} = 1/g^2$ (their eq. (33)).
  The page takes $g = 1$, the dilaton zero at infinity, as their Sec. 5 does, and their units $G_N = 1$; $M \ge 0$ is their supersymmetric positivity bound, eqs. (34) and (35).
- The string metric is their eq. (43), TeX label `stringyelectric`: $ds^2_{\rm str} = \Delta^{-1}dt^2 - d\vec x^2$.
  It is the Einstein metric times $e^{2\varphi}$ with $e^{-2\varphi} = \Delta^{1/2}$, their four-dimensional dilaton of eq. (32) at $g = 1$.
- The singularity is at $r_0 = 2(\sqrt{M^2 + q^2} - M)$, their eq. (41), where $\Delta = 0$ and the volume of the compact space vanishes, eq. (42).

Each chart's domain is $r > r_0$, where $\Delta > 0$.
Kallosh and Linde continue the metric to $r < r_0$ in their eq. (44), with $|\Delta|$ under the root, to discuss tunnelling; that continuation changes the signature of the determinant and is not drawn.

## Checks run by `repulson.py`

`chart_check` holds each printed chart, at three points, to:

1. Kallosh and Linde's eq. (36) or (43), slot by slot.
2. For the Einstein metric, the extreme hole of four charges of Cvetič and Youm, hep-th/9507160 eq. (2), as published in the `four_extreme` chart of `string_bh_three_four_charges`, with charge radii $r_1 = 2M + 2\sqrt{M^2 + q^2}$, $r_2 = 2M - 2\sqrt{M^2 + q^2}$ and $r_3 = r_4 = 0$.
   Then $H_1H_2 = \Delta$, and $r_2 < 0$ is the opposite relative sign of the two electric charges for which Cvetič and Youm find the solutions singular and repelling.
   That chart's own check holds it to the field equations of the heterotic string on a six-torus; the equations are polynomial identities in the charge radii, so they hold with $r_2 < 0$.
3. With $M = 0$, the Ricci scalars Kallosh gives in hep-th/9506113 App. B, $4q^2(2q^2 + r^2)/(r(r^2 - 4q^2)^{5/2})$ for the Einstein metric and $-8q^2(8q^2 + r^2)/(r^2(r^2 - 4q^2)^2)$ for the string metric, her $g$ being this page's $q$; both agree in sign as printed under (-,+,+,+) with the standard contraction.
4. $\Delta(r_0) = 0$.

## Cvetič and Youm's form is this chart

Cvetič and Youm write the same solutions as $ds^2 = -\lambda dt^2 + \lambda^{-1}dr^2 + R\,d\Omega^2$ with $\lambda = r^2/\sqrt{\prod_i(r - \eta C_i)}$ and $R = \sqrt{\prod_i(r - \eta C_i)}$ (hep-th/9507160 eq. (2)).
With the magnetic charges zero and the electric pair $\eta C_1 = r_0$, $\eta C_2 = -2M - 2\sqrt{M^2 + q^2}$, the product is $r^2(r - \eta C_1)(r - \eta C_2) = r^4\Delta$, so $\lambda = \Delta^{-1/2}$ and $R = r^2\Delta^{1/2}$ in the same $r$: their form is the isotropic chart, not another chart.
`check_points.py` in the task folder ran this at two points on 2026-11-18 before the chart was printed.

## Not charts

Ortín's quadruhole, hep-th/9602067 eq. (14), and the massless hole of Cano, Chimento, Ortín and Ruipérez, 1806.08377 eq. (6), $ds^2 = -(1 - Q^2/r^2)^{-1}dt^2 + (1 - Q^2/r^2)d\vec x^2$, are a different solution, with four charges and a constant dilaton, and are not printed here.
