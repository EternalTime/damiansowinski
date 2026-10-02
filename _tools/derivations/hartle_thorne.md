# Hartle and Thorne's slowly rotating star

Why each chart is the one published, what an order is to the checker, and what `print_charts.py`, `null_rays.py`, `conformal.py` and `embedding.py` check before they write.

## Step 1: Hartle and Thorne's chart

The line element is equation (A1) of Hartle and Thorne (1968), the field outside a star that rotates slowly, through second order in its angular velocity.
Their constants are the mass $\mathfrak{M}$, the angular momentum $J$ and the mass quadrupole moment $Q$, in units with $G = c = 1$.
The chart writes them as three lengths: $m = GM/c^2$, $a = J/(Mc)$ and $q = Q/M$, an area, so that $J^2/(\mathfrak{M}r^3)$ is $a^2m/r^3$, $(Q - J^2/\mathfrak{M})/\mathfrak{M}^3$ is $(q - a^2)/m^2$, and the dragging $2J/r^3$ is $2amc/r^3$.
Hartle and Thorne do not print their Legendre functions; Hartle (1967) does, in his (137) and (141), with $\zeta = r/m - 1$:
$Q_2^{\,2}(\zeta) = \tfrac{3}{2}(\zeta^2 - 1)\ln\frac{\zeta + 1}{\zeta - 1} - \frac{3\zeta^3 - 5\zeta}{\zeta^2 - 1}$ and $Q_2^{\,1}(\zeta) = \sqrt{\zeta^2 - 1}\left(\frac{3\zeta^2 - 2}{\zeta^2 - 1} - \tfrac{3}{2}\zeta\ln\frac{\zeta + 1}{\zeta - 1}\right)$.
With $\zeta^2 - 1 = r(r - 2m)/m^2$ and $(\zeta + 1)/(\zeta - 1) = r/(r - 2m)$ the square root in $Q_2^{\,1}$ cancels the one in front of it in (A1), and the two functions the chart names are rational functions of $r$ and $m$ and one logarithm, $L = \ln(r/(r - 2m))$:
$A = \tfrac{5}{8}Q_2^{\,2} = \frac{15r(r - 2m)}{16m^2}L - \frac{5(r - m)(3r^2 - 6mr - 2m^2)}{8mr(r - 2m)}$ and $B = \tfrac{5}{8}\frac{2m}{\sqrt{r(r - 2m)}}Q_2^{\,1} = \frac{5(3r^2 - 6mr + m^2)}{4r(r - 2m)} - \frac{15(r - m)}{8m}L$.
Far from the star $A \to m^3/r^3$ and $B \to m^4/2r^4$, so $g_{tt} \to -(1 - 2m/r + 2qm\,P_2/r^3)$ and $q$ is the quadrupole moment per unit mass of the Newtonian potential, positive for an oblate star.

## Step 2: what an order is to the checker

The line element solves the vacuum equations through second order in the spin and no further, so its exact Ricci tensor is not zero.
`ORDERS` in `verify_metrics.py` declares that $a$ counts as first order and $q$ as second and that the chart keeps the second.
`Reader.truncated` scales each small parameter by its power of one symbol, differentiates the expression as it stands and sets the symbol to zero, so only denominators of zeroth order are ever factored; the canonical form of the whole metric would have to factor $1 + 2h_2P_2$, which did not finish in ten minutes.
`Reader.inverse` builds the inverse as the series $g_0^{-1} - g_0^{-1}\delta g_0^{-1} + \dots$ about Schwarzschild's, and `Geometry` cuts every tensor after each product, which is the same polynomial as cutting the exact tensor and takes about a minute for the whole chart.
A published connection or curvature component that carries a term beyond the order is a disagreement.
The metric and its inverse are published as the line element writes them, $g_{rr} = (1 - 2j_2P_2)/F$ and $g^{rr} = F/(1 - 2j_2P_2)$, each compared to the order, so the drawings read a metric and its exact inverse.
Every value is printed order by order: Schwarzschild's part, the part linear in $a$, the part in $a^2$ that Kerr has, and the part in $q - a^2$, which alone holds the logarithm.

## Step 3: Kerr

`hartle_thorne_kerr` pulls Boyer and Lindquist's line element with $J = Ma$ back through $r \to r - \frac{a^2}{2r^3}\left((r + 2m)(r - m) - \cos^2\theta\,(r - 2m)(r + 3m)\right)$ and $\theta \to \theta - \frac{a^2}{2r^3}(r + 2m)\sin\theta\cos\theta$ and finds the chart at $q = a^2$ in every slot, to second order.
That is Hartle and Thorne's (A5) with the sign of its $\cos^2\theta$ term as Berti, White, Maniopoulou and Bruni (2005) correct it; with the sign as printed the slots $tt$, $rr$, $r\theta$, $\theta\theta$ and $\phi\phi$ miss, which the function also checks.
Hartle and Thorne print $J = -\mathfrak{M}a$ for the form of Boyer and Lindquist they compare with; with the sign of $g_{t\phi}$ in the collection's Kerr it is $J = Ma$.

## Step 4: Lense and Thirring's chart

The line element is the field of Lense and Thirring (1918) as Baines, Berry, Simpson and Visser (2021) write it in their equation (1), in coordinates isotropic far from the body: $-(1 - 2m/r)c^2dt^2 + (1 + 2m/r)(dr^2 + r^2d\Omega^2) - (4am/r)\sin^2\theta\,c\,dt\,d\phi$.
It is linear in the source, so `ORDERS` counts $m$ as first order and keeps the first; the spin enters as $am$ and comes with it.
The inverse metric, the connection and the curvature are those of linearised gravity, the Ricci tensor vanishes, and the Kretschmann scalar, of second order in the mass, is zero to the order kept.
`hartle_thorne_check` holds the chart to being Hartle and Thorne's at first order in $a$ pulled back through $r \to r + m$, the isotropic radius to first order in the mass.
The chart holds where $m/r$ is small, and there its planes and its equator are flat to the eye, so nothing is drawn of it.

## Step 5: the Painlevé-Gullstrand chart

The line element is equation (6) of Baines, Berry, Simpson and Visser (2021), equation (2.4) of their arXiv version 2, and it is published exactly, with no order: they offer it as a spacetime in its own right, and its tensors are short.
`hartle_thorne_check` holds it to two things.
With $dt \to dt + \sqrt{2m/r}\,dr/(1 - 2m/r)$ and $d\phi \to d\phi + (2am/r^3)\sqrt{2m/r}\,dr/(1 - 2m/r)$ it is exactly $-(1 - 2m/r)c^2dt^2 + dr^2/(1 - 2m/r) + r^2(d\theta^2 + \sin^2\theta(d\phi - (2am/r^3)c\,dt)^2)$, Hartle and Thorne's line element with every term of second order dropped from $F$, $h_2$, $j_2$ and $k_2$.
Its Ricci scalar is $18a^2m^2\sin^2\theta/r^6$, their (12) with $J = am$.
Its Kretschmann scalar comes out as the Weyl scalar of their (14) plus $1836\,a^4m^4\sin^4\theta/r^{12}$, which is $\tfrac{17}{3}R^2$; their (15) prints that term both as $1728\,J^4\sin^4\theta/r^{12}$ and as $\tfrac{17}{3}R^2$, and the second is the one that holds.
Every half power of $m$ and $r$ comes from the speed $\sqrt{2m/r}$, so each value is printed as a rational function plus a rational function times that root.

## Step 6: the declared star

Every diagram is drawn for one star in units of $m$: the surface at $R = 6\,m$, which is 12.4 km for 1.4 solar masses, the spin $a = m/4$, and $q = 4a^2$, four times Kerr's quadrupole moment, inside the range of 2.0 to 12.1 that Laarakkers and Poisson (1999) found for neutron stars.
At these values $h_2$, $j_2$ and $k_2$ are $0.0019$, $0.0016$ and $-0.0018$ at the surface, so every drawing is Schwarzschild's exterior to the eye, and the dragging at the surface is a speed of $2amc/R^2 = 0.014\,c$ around the axis.
No figure in three dimensions is drawn: a cone tipped by $0.014\,c$ cannot be told from an upright one.

## Step 7: the diagrams

Spacetime diagrams: the plane of $t$ and $r$ on the axis, where $\sin\theta = 0$ removes the dragging term and the rays stay on the plane, and the equator with $\phi$ divided out, in Hartle and Thorne's chart and in the Painlevé-Gullstrand one, from the surface out to $14\,m$.
Hartle and Thorne's line element completes the square in $d\phi$, so the metric orthogonal to the circles is $-F(1 + 2h_2P_2)c^2dt^2 + (1 - 2j_2P_2)dr^2/F$ at both latitudes, and `--verify` checks the rays against $ct \pm r_*$ with `ht_star`, the quadrature of $\sqrt{(1 - 2j_2P_2)/(1 + 2h_2P_2)}/F$ from `ht_functions`, which writes (A1) out again independently of the file.
The Painlevé-Gullstrand plane is $-c^2dt^2 + (dr + \sqrt{2m/r}\,c\,dt)^2$ on the axis and on the equator for every spin, and its rays are checked against $ct + r - 2\sqrt{2mr} + 4m\ln(\sqrt{r/2m} + 1)$ and $ct - r - 2\sqrt{2mr} - 4m\ln(\sqrt{r/2m} - 1)$.
The lifted rays of the divided out views are checked null against the published metric and geodesic by the published Christoffel symbols, which Hartle and Thorne's chart holds to second order; the miss, of third order, is below the check's $10^{-4}$.

Conformal diagram: the equatorial plane with $\phi$ divided out, one view for each of the two charts, Minkowski's triangle by $p, q = \arctan((ct \mp r_*)/R)$ with the surface of the star on $X = 0$, and for the Painlevé-Gullstrand chart by the arctangents of the two quadratures above, each zero at the surface.

Embedding diagram: the equator of a moment of Hartle and Thorne's $t$, $(1 + j_2)dr^2/F + r^2(1 - k_2)d\phi^2$, from the surface out to $3R$, each circle checked against $r\sqrt{1 - k_2}$ from `ht_functions`, the same slice at $a = q = 0$ checked to be Flamm's paraboloid, and the two checked to lie within $m/50$ of each other.
A moment of the Painlevé-Gullstrand time is flat, which the view states and the script checks, and the moment embedded is marked on Hartle and Thorne's equatorial drawings alone: the two line elements agree to first order in the spin and no further.
