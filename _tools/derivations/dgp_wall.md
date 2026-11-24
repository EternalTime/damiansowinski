# The Dvali-Gabadadze-Porrati braneworld: derivation notes

The page holds four charts of flat space of five dimensions with a wall in it, the model of Gia Dvali, Gregory Gabadadze, and Massimo Porrati, "4D gravity on a brane in 5D Minkowski space", Phys. Lett. B 485, 208 (2000), hep-th/0005016 (DGP), and its cosmology by Cédric Deffayet, "Cosmology on a brane in Minkowski bulk", Phys. Lett. B 502, 199 (2001), hep-th/0010186 (D), read from the arXiv TeX source.
Every chart is in the signature $(-,+,+,+,+)$, $r_c$ is the crossover length, and $c$ stands in the line element alone.
`dgp_wall.py` beside this file holds the charts and their checks, and `dgp_wall_drawings.py` the drawings.

## The equations the metric solves

DGP's action is $\frac{1}{2\kappa^2}\int\sqrt{-g}\,R_5 + \frac{1}{2\mu^2}\int\sqrt{-h}\,R_4 + S_m$, the second term the curvature term of the wall itself, with $\kappa^2 = 8\pi G_5/c^4$ and $\mu^2 = 8\pi G/c^4$ for Newton's constants $G_5$ and $G$ in five and four dimensions.
Off the wall the bulk is empty and $R_{\mu\nu} = 0$.
On the wall D's junction condition, his (junc), is $[K_{ab}] - h_{ab}[K] = -\kappa^2(T_{ab} - G_{ab}/\mu^2)$, with $G_{ab}$ the Einstein tensor of the metric $h_{ab}$ on the wall.
With the mirror symmetry the jump is twice the extrinsic curvature on either side, its normal pointing away from the wall, and with $r_c = \kappa^2/2\mu^2 = G_5/2G$ the condition reads $G_{ab} - 8\pi GT_{ab}/c^4 = (K_{ab} - h_{ab}K)/r_c$, which is what the shared convention states: Einstein's equation on the wall with the extrinsic curvature as an added source, so the theory on the wall is not Einstein's.
DGP's $r_c = M_P^2/(2M^3)$, the History's, is the same length, with $M^3 = c^4/16\pi G_5$ and $M_P^2 = c^4/16\pi G$ in their normalisation of the action.

## The charts and their sources

- `gaussian_normal`, "Deffayet, Gaussian Normal": D's (cosmback) and (bulkmet) with $b = 1$ and $n_0 = 1$, $a(\tau, y) = a_0 + \epsilon|y|\sqrt{\dot a_0^2 + k}$ and $n = 1 + \epsilon|y|\ddot a_0/\sqrt{\dot a_0^2 + k}$, for the two branches $\epsilon = \pm1$, $\epsilon$ the sign of the jump of $\partial_ya$.
  The coordinate $y$ is written for D's $\epsilon|y|$, so that one chart with no $\epsilon$ in it holds both branches: the branch $\epsilon = +1$ keeps two copies of the side $y \ge 0$ and the branch $\epsilon = -1$ two copies of the side $y \le 0$, glued at $y = 0$.
  Written this way the Riemann tensor vanishes with no use of $\epsilon^2 = 1$, which the checker could not know; with $\epsilon$ kept, the spatial components hold $(1 - \epsilon^2)(\dot a^2 + k)$.
  $a$ is a length, $k \in \{-1, 0, 1\}$ and $\chi$ a pure number, as on the universe on a wall moving past a black hole, and the dot is a derivative along $x^0 = ct$.
  $A$ and $N$ are names the chart defines, held as functions while the tensors are built, so the Christoffel symbols are written in them; the Riemann tensor, checked to vanish with both written out, is printed as zero, and every curvature tensor with it.
  $N = 0$ where $\ddot a < 0$ on the side $y > 0$, and $A = 0$ at $y = -a/\sqrt{\dot a^2 + k}$ on the side $y < 0$, where the normals to the wall meet; both are coordinate singularities.
- `global`, "Global": the empty wall of the branch $\epsilon = +1$, D's (bulkmet) at $k = 1$ and $a_0 = r_c\cosh(ct/r_c)$, which solves the Friedmann equation with $\rho = 0$, $\sqrt{H^2 + c^2/a^2} = c/r_c$; it is $(1 + |y|/r_c)^2(-c^2dt^2 + r_c^2\cosh^2(ct/r_c)d\Omega_3^2) + dy^2$, both copies, with $d\Omega_3^2 = d\psi^2 + \sin^2\psi\,d\Omega_2^2$.
- `planar`, "Planar": the same wall in flat slices, $a_0 = e^{ct/r_c}$ at $k = 0$ with $x_i$ lengths, $(1 + |y|/r_c)^2(-c^2dt^2 + e^{2ct/r_c}dx_i^2) + dy^2$; it is the hyperboloid's flat slicing scaled by $1 + |y|/r_c$, $X^A = (r_c + |y|)Y^A$, and covers the half $Y^0 + Y^4 > 0$ of each copy, D's $Y^0 - Y^5 > 0$ in his labels.
- `inertial`, "Inertial": D's flat chart (canon), $-c^2dT^2 + dR^2 + R^2d\Omega_3^2$, in which his (hyperbol) at constant $H = c/r_c$ is the hyperboloid $R^2 - c^2T^2 = r_c^2$; "the SI solution covers the exterior", so the side kept is $R \ge \sqrt{r_c^2 + c^2T^2}$, glued along the wall to its mirror image.
  The global chart is this one along $cT = (r_c + |y|)\sinh(ct/r_c)$ and $R = (r_c + |y|)\cosh(ct/r_c)$.

The curvature of the global and planar charts is a delta on the wall, $G^a{}_b = (6/r_c)\delta(y)\delta^a{}_b$ along it, read by the checker as a distribution as the domain wall's is, and their Kretschmann scalar, the square of that delta, is stated for $y \neq 0$.
That delta is $-2(K_{ab} - h_{ab}K)$, $-2r_cG_{ab}$ of the empty wall's de Sitter space, which has the sign of a negative tension, as Gregory, Kaloper, Myers, and Padilla note.

## The checks

Each chart is checked apart from its printed tensors, in `dgp_wall.chart_check`.

- `gaussian_normal`: every component of the Riemann tensor vanishes with $A$ and $N$ written out, for any $a(t)$ and $k$; the metric at $y = 0$ is Friedmann's with $t$ the proper time; and with $K_{ab} = \epsilon\,\partial_yg_{ab}/2$ at $y = 0$, $\rho$ from the Friedmann equation of the branch, $H^2 + kc^2/a^2 - (\epsilon c/r_c)\sqrt{H^2 + kc^2/a^2} = 8\pi G\rho/3$, and $p$ from $\dot\rho = -3H(\rho + p/c^2)$, the junction condition holds exactly on both branches and fails with the other branch's Friedmann equation, so the check can fail.
  This is D's (refrite), $\frac{\mu^2}{3}\rho = H^2 + k/a_0^2 - 2\epsilon\frac{\mu^2}{\kappa^2}\sqrt{H^2 + k/a_0^2}$ with $2\mu^2/\kappa^2 = 1/r_c$.
  At $a = r_c\cosh(ct/r_c)$ and $k = 1$ the side $y \ge 0$ is the global chart with $\chi = \sin\psi$, checked at points in thirty digits.
- `global` and `planar`: on both sides of the wall the metric is flat space pulled back along $X^A = (r_c + |y|)Y^A$, $Y^A$ the point of the unit hyperboloid each chart stands at, at points in thirty digits for two values of $r_c$; $K_{ab} - h_{ab}K = r_cG_{ab}$ on the wall; the delta in $G_{ab}$ is $-2r_cG_{ab} = (6/r_c)h_{ab}$; and the inside of the wall fails the junction condition of an empty wall.
- `inertial`: flat space in spherical coordinates pulled back, and the wall $y = 0$ of the global chart at $R^2 - c^2T^2 = r_c^2$ with the side $y > 0$ outside it.

`verify_metrics.py --metric dgp_wall` checks every printed value of the four charts against sympy and the dimensions of every term; it took 11 seconds on 24 November 2026, and the four charts printed in under ten.
`_tools/test_dgp_wall.py` holds each line element to its source equation at points, the Friedmann equation of both branches to the History's form, $r_c = 1.4\,c/H_0$ to a flat universe with a density parameter of matter of 0.3, and the drawings to the member they draw.

## The drawings

Every drawing is in units of $r_c$.

- Spacetime diagrams: the Gaussian normal chart at $a = r_c\cosh(ct/r_c)$ and $k = 1$, $\chi = 1/2$, over $-1.1 \le y/r_c \le 1.4$, its rays keeping $ct/r_c \mp \ln(1 + y/r_c)$; the planar chart at $x_i = 0$, its rays keeping $ct/r_c \mp \mathrm{sgn}(y)\ln(1 + |y|/r_c)$; and the inertial chart along a line through the centre, $\psi = \theta = \pi/2$, mirrored, its rays at 45°, hatched inside the wall.
  At that member the Gaussian normal chart's side $y > 0$ is the empty wall of the branch $\epsilon = +1$, and its side $y < 0$ the inside of the same hyperboloid, kept by a wall of tension $3c^4/4\pi Gr_c^2$ on the branch $\epsilon = -1$, since there $\sqrt{H^2 + c^2/a^2} = c/r_c$ makes $8\pi G\rho/3 = 2c^2/r_c^2$; its cones close toward the centre $y = -r_c$, where $N = A = 0$, and it is hatched beyond.
  No Christoffel symbol turns a ray out of any of the three planes, so every curve drawn is a null geodesic.
  The global chart's plane of $t$ and $y$ at $\psi = \theta = \pi/2$ has the metric $-(1 + |y|/r_c)^2c^2dt^2 + dy^2$, the planar chart's plane, so it would only repeat that drawing and is not drawn, as the domain wall's global and conformal planes are not.
- Conformal diagrams: Minkowski's maps $p, q = \arctan(cT \mp R)$ of flat space's plane of $T$ and $R$, which send the wall to the line $X = \pi/2$.
  The global, planar, and inertial charts each draw the two copies of the outside, the side $y < 0$ the mirror image of the side $y > 0$ in the wall, a diamond with $i^0$ at its two side corners and the wall running from the middle of $\mathscr{I}^-$ to the middle of $\mathscr{I}^+$.
  The Gaussian normal chart at its member draws one copy of flat space, Minkowski's triangle, its side $y > 0$ outside the wall and its side $y < 0$ inside it out to the light cone of the centre, which $t$ reaches only at either end.
- Embedding diagram: the equator $\psi = \theta = \pi/2$ of the global chart's moments from $ct = -r_c$ to $r_c$, out to $|y| = 3r_c/2$ on either side, as the domain wall's moments are drawn: the circles grow at $|d\rho/dy| = \cosh(ct/r_c)$, faster than the distance out to them but at $ct = 0$, so each moment lies in three dimensional Minkowski space, each side the cone $\rho = (r_c + |y|)\cosh(ct/r_c)$, $Z = y\sinh(ct/r_c)$, the two meeting at their narrowest circle, the wall.
  At $ct = 0$ the moment is the flat plane outside the circle of radius $r_c$ taken twice, which `exceptions/dgp_wall.json` records; the movie has 41 frames, one every $0.05\,r_c$ of $ct$.
  Each moment is the line $t = t_k$ on the planar and Gaussian normal planes, the latter on its side $y > 0$ alone, $cT = |R|\tanh(ct_k/r_c)$ on the inertial plane, and the same curve under Minkowski's maps on every conformal view.
- Curvature figure: not drawn.
  Every curvature scalar vanishes off the wall in every chart, $K = 0$ and $R = 0$ there, and the one curvature is the delta on the wall, $R = -16\delta(y)/r_c$, which no sampled curve draws, so a figure would draw two zero lines; the flat bulk is stated by the scalars themselves.

Every view's text is within the soft 400 characters of the shared count in `view_text.py` but the inertial spacetime diagram's, 408 with its approved caption.
