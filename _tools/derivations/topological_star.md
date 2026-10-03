# The topological star

Ibrahima Bah and Pierre Heidmann, Phys. Rev. Lett. 126, 151101 (2021), arXiv:2011.08851, their (1.2), (1.5) and (1.6): the solution of Einstein's equations with Maxwell's field in five dimensions

$$ds^2 = -f_S\,dt^2 + f_B\,dy^2 + \frac{dr^2}{f_Sf_B} + r^2d\Omega^2, \qquad f_S = 1 - \frac{r_S}{r}, \quad f_B = 1 - \frac{r_B}{r},$$

with $F = P\sin\theta\,d\theta\wedge d\phi$ and $\kappa_5^2P^2 = 3r_Sr_B/2$, and $y$ a circle of period $2\pi R_y$.
For $r_B > r_S$ it is the star, which ends on the bubble $r = r_B$; for $r_S > r_B$ it is a black string with that bubble behind its horizon.
The line element is the class of solutions of Miyamoto and Kudoh, JHEP 2006(12), 048, their (2.4) at $d = 4$, and of Stotyn and Mann, Phys. Lett. B 705, 269 (2011), their (1) and (6), which Bah and Heidmann review.
Its five charts are written by `_tools/derivations/print_charts.py --metric topological_star`, and `verify_metrics.py --system topological_star/<chart>` checks each in seconds.

## Step 1. The parameters

$r_S$ and $r_B$ are lengths, and so is $y$, as in Bah and Heidmann's paper; $t$ is a time, so $c$ stands where they set it to one.
The mass seen from four dimensions is their (1.9), $M = 2\pi(2r_S + r_B)/\kappa_4^2$, which with $\kappa_4^2 = 8\pi G/c^4$ is $M = c^2(2r_S + r_B)/(4G)$; at $r_B = 0$ it is Schwarzschild's $c^2r_S/(2G)$.
The square of the magnetic charge is their $Q_m^2 = 3r_Br_S/(2\kappa_4^2)$, proportional to $r_Sr_B$.

## Step 2. The charts and their sources

- `bah_heidmann`: their (1.2), (1.5) and (1.6), the star, with $R_y = 2\sqrt{r_B^3/(r_B - r_S)}$, their (2.2) at $k = 1$, the period that makes the bubble smooth.
- `bubble`: their $\rho^2 = 4(r - r_B)/(r_B - r_S)$, under their (2.1), taken over the whole chart, and the angle $\psi = y/R_y$. With $r = r_B + (r_B - r_S)\rho^2/4$, $f_S = (r_B - r_S)(4 + \rho^2)/(4r)$, $f_B = (r_B - r_S)\rho^2/(4r)$ and $dr^2/(f_Sf_B) = 4r^2d\rho^2/(4 + \rho^2)$, so $$ds^2 = -\frac{(r_B - r_S)(4 + \rho^2)}{4r}c^2dt^2 + \frac{4r^2}{4 + \rho^2}d\rho^2 + r^2d\Omega^2 + \frac{r_B^3\rho^2}{r}d\psi^2,$$ whose plane of $\rho$ and $\psi$ is $r_B^2(d\rho^2 + \rho^2d\psi^2)$ at the bubble, their (2.1). Bah and Heidmann use $\rho$ near the bubble only; the chart over the whole star is our own, derived here, and `topological_star_check` holds it to theirs pulled back and to $g_{\psi\psi}/g_{\rho\rho} = \rho^2(1 + O(\rho^2))$. The chart names $r$ for the expression, as Weyl's chart of the Curzon-Chazy particle names $R$, and `bubble_radius` writes every power of $4r_B + (r_B - r_S)\rho^2$ as one of $4r$.
- `eddington_finkelstein_ingoing`: the black string, $r_S > r_B$, in $v = ct + r_*$ with $dr_*/dr = 1/(f_S\sqrt{f_B})$, so that $$ds^2 = -f_S\,dv^2 + 2\sqrt{\frac{r}{r - r_B}}\,dv\,dr + r^2d\Omega^2 + f_B\,dy^2.$$ Bah and Heidmann draw the string's Penrose diagram, their figure 3, and say the whole interior wants Kruskal's coordinates; this chart, which crosses the horizon and reaches the bubble, is our own, and `topological_star_check` holds it to theirs pulled back.
- `extremal`: $r_S = r_B = m$ in the isotropic radius $\rho = r - m$, as the companion paper writes it, JHEP 2021(09), 147, arXiv:2012.13407, section 3.4.2; the letter's (3.1) is the same metric at $r = \rho^2 + m$.
- `einstein`: the Einstein metric of four dimensions, their (1.8), $f_B^{1/2}$ times the rest of their (1.2) on a surface of constant $y$, with $e^{2\Phi} = f_B^{-1/2}$.

`topological_star_check` in `print_charts.py` holds every chart of five dimensions to $R_{ab} = \kappa_5^2(F_{ac}F_b{}^c - g_{ab}F^2/6)$ and to Maxwell's equations; Bah and Heidmann's chart to the published static black string at $r_B = 0$, to Schwarzschild's published metric with its time made the circle at $r_S = 0$, and to the symmetry $g_{tt} \leftrightarrow -g_{yy}$ under $r_S \leftrightarrow r_B$; and the Einstein metric to the reduction and to the equations of the action of their footnote, $G_{ab} = 6(\partial_a\Phi\partial_b\Phi - \tfrac12g_{ab}(\partial\Phi)^2) + \tfrac32r_Sr_Be^{-2\Phi}(F_{ac}F_b{}^c - \tfrac14g_{ab}F^2)$ with $F = \sin\theta\,d\theta\wedge d\phi$, Maxwell's equation with $e^{-2\Phi}$, and $\Box\Phi = -\tfrac18r_Sr_Be^{-2\Phi}F^2$, at six points in forty digits.

## Step 3. The tortoise coordinate

With $u = \sqrt{r - r_B}$, $\int dr/(f_S\sqrt{f_B}) = \int 2(u^2 + r_B)^{3/2}du/(u^2 + r_B - r_S)$, which is
$$r_* = \sqrt{r(r - r_B)} + (r_B + 2r_S)\,\mathrm{arsinh}\sqrt{r/r_B - 1} + \frac{2r_S^{3/2}}{\sqrt{r_B - r_S}}\arctan\sqrt{\frac{r_S(r - r_B)}{(r_B - r_S)r}}$$
for the star and, for the black string, the same with the last term $\frac{r_S^{3/2}}{\sqrt{r_S - r_B}}\ln\left|\frac{a - b}{a + b}\right|$, $a = \sqrt{r_S(r - r_B)}$, $b = \sqrt{(r_S - r_B)r}$.
Both vanish on the bubble.
The logarithm's coefficient is $1/(2\kappa)$ with the surface gravity $\kappa = \sqrt{r_S - r_B}/(2r_S^{3/2})$, the temperature the companion paper gives.
`slices.topological_star_rstar` is that closed form, checked against quadrature, and the extremal string's $\rho_* = \sqrt{\rho(\rho + m)} + 3m\,\mathrm{arsinh}\sqrt{\rho/m} - 2m\sqrt{1 + m/\rho}$ is `slices.extremal_string_rstar`.

## Step 4. What the diagrams draw

The star is drawn at $r_B = 1$ and $r_S = 3/4$, where $R_y = 4r_B$, inside the range $r_S < r_B < 2r_S$ in which Bah, Dey and Heidmann, JHEP 2022(04), 168, find it stable, and below $3r_S/2$, so it has two photon spheres, Heidmann, Speeney, Berti and Bah's second kind.
The black string is drawn at $r_S = 1$ and $r_B = 3/4$, the same two numbers exchanged, inside the range $r_S/2 < r_B < r_S$ in which Miyamoto and Kudoh find no Gregory-Laflamme mode; its surface gravity is then $1/4$.
The extremal string is drawn at $m = 1$.

On the star's plane of $t$ and $r$ light from the bubble reaches $r = 2r_B$ at $ct = r_*(2r_B) = 5.92\,r_B$, and through the bubble from $\rho = 4$ to $\rho = 4$ on the other side in $11.84\,r_B$.
A clock at rest on the bubble runs at $\sqrt{1 - r_S/r_B} = 1/2$, and on the bubble $c\,dt/d\rho = r_B^{3/2}/\sqrt{r_B - r_S} = 2r_B$.
The Kretschmann scalar is greatest on the bubble, $(48r_B^2 - 72r_Br_S + 43r_S^2)/(4r_B^6) = 4.55/r_B^4$, and on the string's bubble it is $22.5/r_S^4$; the extremal string's is $19/(4m^4)$ on its horizon, and the Einstein metric's diverges at $r_B$ as $(r - r_B)^{-3}$.

The conformal diagrams are Minkowski's triangle for the star, $p, q = \arctan((ct \mp r_*)/\ell)$ with $\ell = 6r_B$, Bah and Heidmann's figure 1, and the whole diamond for the chart about the bubble with $r_*$ taken negative on the side $\psi = \pi$; Kruskal and Szekeres's square for the black string, $U = \mp e^{-\kappa(v - 2r_*)}$, $V = e^{\kappa v}$, in which $r_*(r_B) = 0$ puts the bubble on $UV = 1$, their figure 3; the whole diamond for the extremal string's exterior, whose left edges are the degenerate horizon; and the triangle again for the Einstein metric, with a singular left edge.
The free particles on the diamond about the bubble, released from rest at $\rho_0 = 2$ and $4$, obey $\ddot\rho = -(h'/2h)\dot\rho^2 - E^2(N^2)'/(2hN^4)$ and $\dot t = E/N^2$ with $N^2 = (4 + \rho^2)/(16 + \rho^2)$, $h = (16 + \rho^2)^2/(64(4 + \rho^2))$ and $E = N(\rho_0)$; their periods in $ct$ are $35.56$ and $52.83\,r_B$, above the $16\pi/\sqrt3 = 29.02\,r_B$ of a small swing, and the quadrature $4\int_0^{\rho_0}E\sqrt{h}\,d\rho/(N^2\sqrt{E^2/N^2 - 1})$ gives the same.

## Step 5. The embeddings

The cigar is the surface of $\rho$ and $\psi$ of the chart about the bubble at one moment and one point of the sphere: $g_{\rho\rho} = 4r^2/(4 + \rho^2)$ and the circles have the radius $R_y\sqrt{f_B} = 4\rho/\sqrt{16 + \rho^2}$.
Along it $dR/ds = R_y r_B\sqrt{f_S}/(2r^2)$, which is $1$ on the bubble and falls from there wherever $r > 5r_S/4$, so the whole surface stands in flat space for $r_B > 5r_S/4$, as at the drawn $r_B = 4r_S/3$.
From the bubble to $\rho$ the proper distance is $\sqrt{(r - r_S)(r - r_B)} + (r_S + r_B)\,\mathrm{arsinh}(\rho/2)$, $7.79\,r_B$ to $\rho = 8$, where $r = 5r_B$.
The equator is the surface of $r$ and $\phi$ of Bah and Heidmann's chart at $y = 0$, carried through the bubble's equator onto $y = \pi R_y$, the line of $\rho$ through the origin of the plane of $\rho$ and $\psi$; its height is the quadrature of $\sqrt{((r_S + r_B)r - r_Sr_B)/((r - r_S)(r - r_B))}$.
