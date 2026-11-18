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

## Drawings

Every drawing is of the massless member, $M = 0$, in units of $q$, so $\Delta = 1 - 4q^2/r^2$ and the singularity is $r = 2q$; `repulson_drawings.py` beside this file draws them.

- Spacetime diagrams, one view per chart on the plane of $t$ and $r$ at $\theta = \pi/2$, $\phi = 0$, for $2q < r < 6q$ and $|ct| < 2q$.
  The radial light rays of both charts obey $c\,dt = \pm\sqrt{\Delta}\,dr$, since the string metric is the Einstein metric times $1/\sqrt{\Delta}$, so $ct \mp r_*$ is constant along them with $r_* = \sqrt{r^2 - 4q^2} - 2q\arccos(2q/r)$, zero at $r = 2q$; `--verify` holds every ray to it.
  A ray moving in reaches the singularity after the finite time $r_*/c$, as Cvetič and Youm state for massless particles of zero angular momentum (hep-th/9507160).
- Conformal diagram, one view per chart: Minkowski's triangle by $\arctan(ct \mp r_*)$, the singularity a timelike line where Minkowski's centre stands, each point a 2-sphere; the Kretschmann scalar is checked to diverge on it.
- Embedding diagram, one view per chart at $t = 0$ on the equator.
  For the Einstein metric $\rho = r\Delta^{1/4}$ and $(d\rho/dr)^2 - g_{rr} = 4(1 - 3q^2/r^2)/(r^2\Delta^{3/2}) > 0$ for every $r > 2q$, so every circle grows faster than the distance out to it and the slice is drawn in Minkowski space, $dX^2 + dY^2 - dZ^2$, as anti-de Sitter space's static slice is.
  Near $r = 2q$ the surface turns toward a light cone of that space, its metric length about $(r - 2q)/q$ of its drawn length, so the profile is drawn from $r = 2.05q$, with knots crowded toward that end and eleven decimals, which keeps every chord within $2.7 \times 10^{-5}$ of its proper length.
  For a massive member the slice would be in Minkowski space inside $r_c = 2q^2/M$ and in flat space outside it, level at $r_c$; that member is not drawn.
  The string metric's slice is $dr^2 + r^2d\phi^2$, the flat plane outside the circle $r = 2q$.
- Curvature figure: the stored Kretschmann and Ricci scalars of the Einstein metric against $r$ from $2q$ to $6q$, on logarithmic scales cut near $r = 2q$, where both grow without bound; the Ricci scalar is checked against Kallosh's closed form, hep-th/9506113 App. B.
