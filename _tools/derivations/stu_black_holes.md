# The black holes of anti-de Sitter space with three charges: derivation notes

The page holds five charts of the static black holes of the STU model of $N = 2$ gauged supergravity in five dimensions, found by Klaus Behrndt, Mirjam Cvetič, and Wafic Sabra, "Non-extreme black holes of five-dimensional $N = 2$ AdS supergravity", Nucl. Phys. B 553, 317 (1999), hep-th/9810227 (BCS), read from the arXiv TeX source.
Every parameter and every radius is a length or a length squared, $L = 1/g$ is the radius of anti-de Sitter space, and $c$ stands in the line element alone.
`stu_black_holes.py` beside this file holds the charts and their checks, and `stu_black_holes_drawings.py` the drawings.

## The charts and their sources

- `sphere`, "Spherical Horizon, $k = 1$": BCS eq. (390) with the ansatz (110) and (170), and M. Cvetič and nine coauthors, Nucl. Phys. B 558, 96 (1999), hep-th/9903214, eq. (d5adsbh): $ds^2 = -(H_1H_2H_3)^{-2/3}f\,dt^2 + (H_1H_2H_3)^{1/3}(dr^2/f + r^2d\Omega_3^2)$, $f = k - \mu/r^2 + r^2H_1H_2H_3/L^2$, $H_i = 1 + q_i/r^2$, with $q_i = (\mu/k)\sinh^2\beta_i$, BCS (392).
- `flat`, "Flat Horizon, $k = 0$": Cvetič and nine coauthors' eq. (d5adsbhk0), after the rescaling $\sinh^2\beta_i \to k\sinh^2\beta_i$ and $k \to 0$ they describe, with the flat horizon $r^2(dx^2 + dy^2 + dz^2)/L^2$ of Dam T. Son and Andrei O. Starinets, JHEP 03 (2006) 052, hep-th/0601157, section 2, so that $x$, $y$, $z$ are lengths.
- `hyperbolic`, "Hyperbolic Horizon, $k = -1$": the same line element with $d\Omega_{3,-1}^2 = d\chi^2 + \sinh^2\chi\,d\Omega_2^2$, BCS's $d\Omega_{3,k}$ and Cvetič and nine coauthors' $H^3$.
- `continued`, "Past $r = 0$, Radius $\rho$": the spherical chart in Robert C. Myers and Øyvind Tafjord's radius $\rho^2 = r^2 + q_3$, JHEP 11 (2001) 009, hep-th/0109127, eq. (rhoeq), for $q_1 \ge q_2 \ge q_3 > 0$.
  With $P_1 = r^2H_1$, $P_2 = r^2H_2$, and $D = r^4f$ written in $\rho$, $H_1H_2H_3 = \rho^2P_1P_2/(\rho^2 - q_3)^3$ and the metric is $-D\,dt^2/(\rho^{4/3}(P_1P_2)^{2/3}) + \rho^{8/3}(P_1P_2)^{1/3}d\rho^2/D + \rho^{2/3}(P_1P_2)^{1/3}d\Omega_3^2$, written with $\rho^{2/3}$ in place of $(\rho^2)^{1/3}$, which sympy reads as $|\rho|^{2/3}$.
  Myers and Tafjord say $\rho \ge 0$ covers the whole geometry with the naked singularity at $\rho = 0$; $r = 0$ is the regular 3-sphere $\rho = \sqrt{q_3}$.
- `brane`, "Black Brane, Radius $u$": the flat chart in Son and Starinets's $u = r_+^2/r^2$, their eqs. (metric_u_3) and (identif), and Liam Gladden, Victor Ivo, Pavel Kovtun, and Starinets, Phys. Rev. D 111, 086030 (2025), 2412.12353, eq. (STU-metric), with $\kappa_i = q_i/r_+^2$, $H_i = 1 + \kappa_iu$, $f = H_1H_2H_3 - u^2(1 + \kappa_1)(1 + \kappa_2)(1 + \kappa_3)$.
  Their $T_0 = r_+/(\pi L^2)$, in units $\hbar = k_B = c = 1$, is written out, so $(\pi T_0L)^2 = r_+^2/L^2$ and the chart is in lengths; the conventions give $T_0 = \hbar cr_+/(\pi k_BL^2)$.

Every chart holds its harmonic functions and its metric function, or $P_1$, $P_2$, and $D$, as names with declared slopes, so every printed value is a rational function of them, the radius, and the parameters, as the charts of `string_bh_three_four_charges` are.
The slopes are $H_i' = -2(H_i - 1)/r$ and $f' = 2(k - f)/r + (2r/L^2)(H_1H_2 + H_1H_3 + H_2H_3 - H_1H_2H_3)$ in $r$; $P_i' = 2\rho$ and $D' = 2\rho(2\rho^2 - 2q_3 - \mu) + (2\rho/L^2)(P_1P_2 + \rho^2P_1 + \rho^2P_2)$ in $\rho$; and in $u$, $H_i' = (H_i - 1)/u$ and $f' = (H_1H_2H_3 - H_1H_2 - H_1H_3 - H_2H_3 + 2f)/u$, which leaves no $\kappa_i$ in any value of the brane chart: written with $\kappa_i$ the chart was 445 KB, and with these slopes it is 81 KB.
The checker holds each declared slope to its name's definition.

## The checks

Each chart is checked at points in forty digits, apart from the printed tensors, in `stu_black_holes.chart_check`:

- the Einstein, Maxwell, and scalar equations of Cvetič and nine coauthors' Lagrangian (d5gauged), $R - \frac{1}{2}(\partial\varphi)^2 + 4g^2\sum X_i^{-1} - \frac{1}{4}\sum X_i^{-2}(F^i)^2$, with $X_i = (H_1H_2H_3)^{1/3}/H_i$ and $A^i = Q_i\,c\,dt/(r^2 + q_i)$, $Q_i^2 = q_i(\mu + kq_i)$, which is BCS's $F_{r0} = -\frac{\sqrt k}{2}H_I^{-2}\partial_r\tilde H_I$ with (392) and Cvetič and nine coauthors' $A^i = \sqrt k(1 - H_i^{-1})\coth\beta_i\,dt$; the Chern-Simons term is silent for electric fields, and $\sum(\partial\ln X_i)^2 = (\partial\varphi)^2$ for $a_i\cdot a_j = 4\delta_{ij} - 4/3$; a field one percent too strong, or the scalars dropped, fails Einstein's equation at once;
- `sphere` at $L \to \infty$ against the live `five_charges` chart of `string_bh_three_four_charges` with $r_0^2 = \mu$ and $r_i^2 = q_i$, BCS's statement that only $f$ differs, and with no charge against the live `spherical` chart of `tangherlini`, $r_h^2 = \mu$;
- `continued` against `sphere` pulled back along $r^2 = \rho^2 - q_3$ where $\rho^2 > q_3$, its field equations on both sides of $\rho^2 = q_3$, and with the three charges equal the Reissner-Nordström-anti-de Sitter hole of five dimensions, $V = 1 - (2q + \mu)/\rho^2 + q(q + \mu)/\rho^4 + \rho^2/L^2$;
- `brane` against `flat` pulled back along $r = r_+/\sqrt u$ with $q_i = \kappa_ir_+^2$ and $\mu = r_+^4(1 + \kappa_1)(1 + \kappa_2)(1 + \kappa_3)/L^2$, $f(1) = 0$, and Son and Starinets's temperature $T_0(2 + \kappa_1 + \kappa_2 + \kappa_3 - \kappa_1\kappa_2\kappa_3)/(2\sqrt{(1 + \kappa_1)(1 + \kappa_2)(1 + \kappa_3)})$ from the surface gravity of the chart.

On the hyperbolic chart $Q_i^2 = q_i(\mu - q_i)$, so its gauge fields are real for $\mu \ge q_i$, which its domains and the description of $\mu$ say.

## The drawings

- Members drawn, in units of the event horizon $r_+ = 1$.
  The spherical and continued charts: $q_1 = 4$, $q_2 = 2$, $q_3 = 1$, $L = 3$, $\mu = 13/3$, where $r^4fL^2 = (x - 1)(x^2 + 17x - 8)$ in $x = r^2$, so $r_-^2 = (\sqrt{321} - 17)/2 = 0.458$, and in $\rho$ the horizons are $\sqrt 2$ and $1.2075$, $r = 0$ is $\rho = 1$, and the singularity $\rho = 0$ is timelike, $D(0) = q_3(q_3 + \mu) > 0$.
  Charges well apart and $L = 3r_+$ keep the inner horizon clear of $r = 0$ and the embedding long enough to read.
  The flat and brane charts: the same charges and $L$, with $\mu = 10/3$, the inner horizon at $x = 2\sqrt6 - 4$, $u = 1.112$.
  The hyperbolic chart: $L = 1$, $q = (1, 1/2, 1/4)$, $\mu = 11/4$, chosen so that $\mu \ge q_i$.
- Spacetime diagrams: one view of the plane of $t$ and the radius for each chart, with ingoing cones, and outgoing for the brane chart, whose $u$ falls outward.
- Embedding diagram: the equatorial plane of the slice $t = 0$ of the spherical chart, circles of radius $r(H_1H_2H_3)^{1/6}$, both exteriors through the throat $r_+$, out to the circle $r = 2.356\,r_+$, where $g_{rr} = (dR/dr)^2$ and the surface lies level; beyond it the circles grow faster than the distance out to them, as in every black hole of anti-de Sitter space, and the page draws the part that lies in flat space.
  `Slice.horizons()` reads no cube roots, so both horizons are found by `horizon_between` and the level circle by `level_radius` to eighty digits.
  Its moment lies on the spherical and continued spacetime diagrams and on both conformal views; the other charts are other members, named in `HIDDEN`.
- Conformal diagram: the tower of Reissner-Nordström-anti-de Sitter, drawn by `HyperbolicHoleDrawing` as `reissner_nordstrom_ads` is, for the spherical and the continued charts, each covering the exterior.
  The tortoise coordinate is $r_*(\rho) = \int_0^\rho s\sqrt{P(s)}/D(s)\,ds$, its two logarithms at the horizons written out and the rest tabulated by quadrature and cubic splines, beyond $\rho = 4$ in $w = 4/\rho$; $r_*(0) = 0$ puts the singularity on the vertical lines, and $r_* \to 2.5518$ on the boundary.
  The flat, hyperbolic, and brane charts get no conformal view: each ends at $r = 0$, a regular surface inside the inner horizon, and no source here continues them past it, as Myers and Tafjord continue the spherical chart.
- Curvature figure: the stored Kretschmann and Ricci scalars of the continued chart against $\rho$ from $0.25$ to $6$, both finite at the horizons and at $\rho = 1$, both growing without bound toward $\rho = 0$, and far out $K \to 40/L^4$ and $\mathcal{R} \to -20/L^2$, anti-de Sitter space of five dimensions.
