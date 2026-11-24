# The flows of the renormalisation group of holography: derivation notes

The page holds three charts of two flows in closed form, walls of one scalar in the gauged supergravity of five dimensions with the metric $ds^2 = e^{2A}(-c^2dt^2 + dx_1^2 + dx_2^2 + dx_3^2) + dr^2$.
The founding walls of Girardello, Petrini, Porrati, and Zaffaroni (hep-th/9810126) and of Freedman, Gubser, Pilch, and Warner (hep-th/9904017) are known only numerically, so no chart of them is drawn; the two flows below are the closed forms their successors printed.
Each source sets the radius of the anti-de Sitter space at the boundary to 1; the page restores it as $L$, reads $t$ as a time and $x_i$ and $r$ as lengths, and keeps $c$ in the line element, with $g_{tt}$ dimensionless.
`holographic_rg_flow.py` beside this file holds the charts and their checks, and `holographic_rg_flow_drawings.py` the drawings.

## The charts and their sources

- `proper_distance`, "Confining Flow, Proper Distance": L. Girardello, M. Petrini, M. Porrati, and A. Zaffaroni, Nucl. Phys. B 569, 451 (2000), hep-th/9909047 (GPPZ99), eqs. (m37) and (m38), read from the arXiv TeX source.
  In their ansatz $ds^2 = dy^2 + e^{2\phi(y)}dx^\mu dx_\mu$ the first order equations of $W = -\tfrac{3}{4}(\cosh(2m/\sqrt{3}) + 1)$ give $\phi(y) = \tfrac{1}{2}(y + \log[2\sinh(y - C_1)])$ and $m(y) = \tfrac{\sqrt{3}}{2}\log[(1 + e^{-(y - C_1)})/(1 - e^{-(y - C_1)})]$.
  So $e^{2\phi} = e^{C_1}(e^{2(y - C_1)} - 1)$, and the shift $y \to y + C_1$ with the constant $e^{C_1}$ taken into $t$ and $x_i$ sets $C_1 = 0$, which the History says as "up to a constant factor".
  With $y = r/L$ the warp factor is $e^{2A} = e^{2r/L} - 1$, $r$ is the proper distance from the curvature singularity at $r = 0$, GPPZ99 (m39), and $e^{2A} \to e^{2r/L}$ is anti-de Sitter space of radius $L$ as $r \to \infty$.
  The scalar $m = \sqrt{3}\,\mathrm{artanh}(e^{-r/L})$ is matter and stands in no component.
- `bfs_u`, "Confining Flow, Coordinate $u$": M. Bianchi, D. Z. Freedman, and K. Skenderis, JHEP 08 (2001) 041, hep-th/0105276 (BFS01), Sec. 3.2, the same flow as $e^{2A} = u/(1 - u)$ with $du/dr = 2(1 - u)$ and $\varphi_B = \tfrac{\sqrt{3}}{2}\log[(1 + \sqrt{1 - u})/(1 - \sqrt{1 - u})]$.
  Restoring $L$, $du/dr = 2(1 - u)/L$, so $u = 1 - e^{-2r/L}$, $u/(1 - u) = e^{2r/L} - 1$, and $g_{uu} = L^2/(4(1 - u)^2)$; $u = 0$ is the singularity and $u \to 1$ the boundary.
  BFS01's scalar is GPPZ99's $m$, since $\sqrt{1 - u} = e^{-r/L}$.
- `coulomb_disc`, "Coulomb Branch Disc": BFS01 Sec. 3.1, $v = e^{\sqrt{6}\Phi}$, $e^{2A} = l^2v^{2/3}/(1 - v)$, $dv/dr = 2v^{2/3}(1 - v)$ for $W = -e^{-2\Phi/\sqrt{6}} - \tfrac{1}{2}e^{4\Phi/\sqrt{6}}$, with the boundary at $v = 1$ and a curvature singularity at $v = 0$, "the origin of a disc distribution of $D3$-branes".
  Restoring $L$, $e^{2A} = (\ell^2/L^2)v^{2/3}/(1 - v)$ and $dv/dr = (2/L)v^{2/3}(1 - v)$, so $g_{vv} = L^2/(4v^{4/3}(1 - v)^2)$.
  It is D. Z. Freedman, S. S. Gubser, K. Pilch, and N. P. Warner, JHEP 07 (2000) 038, hep-th/9906194 (FGPW00), the flow $n = 2$ of their eq. (SUSYCases), $A(\mu) = \tfrac{1}{2}\log|e^{2\mu/\sqrt{6}}/(1 - e^{\sqrt{6}\mu})| + \log(\ell/L)$ for $W(\mu) = -e^{2\mu/\sqrt{6}} - \tfrac{1}{2}e^{-4\mu/\sqrt{6}}$, with $v = e^{-\sqrt{6}\mu}$, $\Phi = -\mu$, and their equations (SUSYCasesstart) $d\mu/d\rho = (g/2)\,\partial W/\partial\mu$ and $dA/d\rho = -(g/3)W$ at $g = 2/L$.
  FGPW00 (TenForm) and the text after it make $\ell$ the radius of the ball, here the disc, over which the D3-branes spread in ten dimensions, a length, which the parameter's description says; BFS01's $l$ is the same constant at $L = 1$.
  In five dimensions $\ell$ only rescales $t$ and $x_i$, so every $\ell > 0$ is one geometry, and the drawings take $\ell = 3L/2$ for the reason given under the embedding diagram.

The two flows are different spacetimes, and `metric_tags.py` decides its eight tags over all three charts together: every chart is conformally flat, since $e^{2A}\eta + dr^2 = e^{2A}(\eta + dz^2)$ with $dz = e^{-A}dr$, static, stationary, and of five dimensions.

## The checks

`chart_check` in `holographic_rg_flow.py` runs in `print_charts.py` before anything is written, in forty digits at three points of each chart, with no symbolic simplification and independently of the printed tensors.
Each chart is held to its source's warp factor and $g_{55}$, to BFS01's first order equations (1storder) with $L$ restored, $dA/dr = -(2/3L)W$ and $d\Phi/dr = W'/L$, and to Einstein's equations $G_{ab} = 2T_{ab}$ of BFS01's action (scaction), $\tfrac{1}{4}R - \tfrac{1}{2}(\nabla\Phi)^2 - V/L^2$ in Lorentzian signature, with $V = \tfrac{1}{2}(\partial_\Phi W)^2 - \tfrac{4}{3}W^2$, their (potential), and to the wave equation $\Box\Phi = V'(\Phi)/L^2$.
These are GPPZ99's (m37) for the confining flow and BFS01's Sec. 3.1 for the disc, whose normalisation agrees with GPPZ99's.
Each is held as well to $A'' = -\tfrac{2}{3}\Phi'^2$, Freedman, Gubser, Pilch, and Warner's inequality $A'' \le 0$ of hep-th/9904017 in these units, and to the Ricci scalar $-20/L^2$ of anti-de Sitter space at the boundary.
The chart of $u$ is held besides to being the proper distance chart pulled back through $u = 1 - e^{-2r/L}$, and the disc to FGPW00's $A(\mu)$ and both of their first order equations.
The Ricci scalars printed agree with $R = -8A'' - 20A'^2$ for this metric: $-4(5e^{2r/L} - 4)e^{2r/L}/(L^2(e^{2r/L} - 1)^2)$, $-(16u + 4)/(L^2u^2)$, and $(4v^2 - 48v - 16)/(3L^2v^{2/3})$, each $-20/L^2$ at the boundary, and the Kretschmann scalars each tend to $40/L^4$ there.
The test points are sympy numbers throughout: a Python division turned one of them into a float of fifteen digits while the checks were written, and the forty digit comparisons failed on it at $10^{-17}$.

The proper distance chart is printed with $e^{2r/L} - 1$ kept whole and its sums in falling powers of $e^{r/L}$, `warp_whole` and `falling_exponentials`; the other two print $1 - u$ and $1 - v$ as the line elements write them, `named_sum`.

## The drawings

Every drawing takes $L = 1$ and holds $x_1 = x_2 = x_3 = 0$ but where it rolls $x_1$ up.

The spacetime diagrams draw the plane of $t$ and the fifth coordinate of each chart, with radial light rays $c\,dt = \pm e^{-A}dr$.
The proper distance chart is drawn on $0 \le r \le 2L$ and $|ct| \le L$, where its cones close toward the singularity as $\sqrt{e^{2r/L} - 1}$ and open toward the boundary.
The chart of $u$ is drawn on $0 \le u \le 1$ and $|ct| \le 0.8L$, where $du/d(ct) = \pm 2\sqrt{u(1 - u)}/L$ closes the cones at both ends; $g^{uu} = 4(1 - u)^2/L^2$ vanishes on the boundary $u = 1$, which the diagram marks.
The disc is drawn on $0 \le v \le 1$ and $|ct| \le 0.8L$, where $dv/d(ct) = \pm 2\ell v\sqrt{1 - v}/L^2$.
Its Kretschmann scalar grows only as $v^{-4/3}$, the fourth inverse power of the proper distance $r \propto v^{1/3}$, which is too slow in $v$ for the edge test at $10^{-5}$ of the drawing, so the row declares $v = 0$ singular in `singular_zero` and the generator checks the scalar there at $10^{-20}$ and $10^{-30}$ in sixty digits.

The conformal diagrams use $z = \int e^{-A}dr$, in which the plane of $t$ and the fifth coordinate is $e^{2A}(-c^2dt^2 + dz^2)$.
For the confining flow $z = L\arcsin(e^{-r/L})$, finite at both ends, so the plane is conformal as it stands to the strip $0 < z < \pi L/2$ between the boundary and the singularity, both timelike, and light crosses it in the time $\pi L/2c$.
It is drawn as $X = \arccos(e^{-r/L})$, the singularity on the left and the boundary on the right, and $T = ct/L$, over $|T| \le \pi/2$ with the strip faded above and below, since $t$ runs on.
The chart of $u$ covers the same plane, so it has no conformal view of its own.
For the disc $dz = (L^2/2\ell)\,dv/(v\sqrt{1 - v})$, so $z = (L^2/\ell)\,\mathrm{artanh}\sqrt{1 - v}$, which diverges logarithmically at $v = 0$: the plane is conformal to half of Minkowski's plane, as FGPW00 say after their (zMetric) for $n = 2$, and $p, q = \arctan((ct \mp w)/\lambda)$ with $w = -z$ and $\lambda = L^2/\ell$ draw it as the left half of Minkowski's diamond, the boundary on $X = 0$.
Radial light reaches $v = 0$ at a finite affine parameter, $\int e^{A}dr \propto \int v^{-1/3}dv$, but at infinite $t$, so the singularity is the two null edges of the half diamond and their corner $i^0$, which the view draws as singular.
Near $v = 0$ the drawing writes $z = \lambda(\ln(1 + \sqrt{1 - v}) - \tfrac{1}{2}\ln v)$, since $1 - v$ rounds to 1 there.

The embedding diagrams draw the surface of $r$, or of $v$, and $x_1$ at $t = 0$, $x_2 = x_3 = 0$, with the strip $0 \le x_1 < 2\pi L$ rolled up so that $x_1$ runs once round, as Randall and Sundrum's page draws its wall.
Its circles have the radius $\rho = Le^{A}$, and along the proper distance $d\rho/dr = \rho A'$.
For the confining flow $d\rho/dr = e^{2r/L}/\sqrt{e^{2r/L} - 1} \ge 2$, least at $e^{2r/L} = 2$: every circle grows faster than the distance out to it, for this or any period of $x_1$, since $d\rho/dr$ diverges at both ends, so no surface of revolution in Euclidean space has this geometry, and the slice is drawn in Minkowski space, $dX^2 + dY^2 - dZ^2$.
It nears the light cone at both ends, and its chords near $r = 0$ are so short in the metric against their length on the page that the drawing starts at $r = L/20$ and runs to $r = 3L/2$, with 719 knots and twelve digits, which keeps every chord within $1.1 \times 10^{-4}$ and every line across within $9 \times 10^{-5}$ of the metric.
For the disc $d\rho/dr = (2\ell/3L)(1 + v/2)/\sqrt{1 - v}$, which diverges at the boundary for every $\ell$ and is at least 1 everywhere exactly when $\ell \ge 3L/2$, so at $\ell = 3L/2$ a surface in Minkowski space has the geometry of the whole slice from the singularity $v = 0$, where the slope is 1 and the surface leaves the axis level, out to $v = 0.9$.

The curvature figure draws the stored Kretschmann and Ricci scalars of the proper distance chart against $r$, from $L/20$ to $3L$, and of the disc against $v$, from 0 to 1, each on a logarithmic axis for $K$ and a signed logarithmic one for $\mathcal{R}$, in plain boxes.
Its checks read the stored scalars only: $K \to 40/L^4$ and $\mathcal{R} \to -20/L^2$ at the boundary, $K > 0$ and $\mathcal{R} < 0$ throughout, and $K$ past $10^6$ near the singularity.

No drawing is left out: every chart has its spacetime diagram, each flow its conformal and embedding views and its curvature view, and each moment $t = 0$ is marked on the spacetime diagrams and the conformal view of its flow.
