# The black holes on a wall of Emparan, Horowitz and Myers

Roberto Emparan, Gary T. Horowitz and Robert C. Myers, "Exact description of black holes on branes", Journal of High Energy Physics 2000(01), 007, arXiv hep-th/9911043, read in full from its TeX source.
Their sequel, "Exact description of black holes on branes II: comparison with BTZ black holes and black strings", Journal of High Energy Physics 2000(01), 021, arXiv hep-th/9912135, read in full from its TeX source.
The papers number their equations by harvmac labels, so each is named here by its label in the TeX.

## Step 1. The charts

`flat_wall` is the first paper's section 3.1, its \adsc with \gx: $ds^2 = \ell^2/(x - y)^2\left(-(y^2 + 2\mu y^3)dt^2 + dy^2/(y^2 + 2\mu y^3) + dx^2/G + G\,d\varphi^2\right)$ with $G = 1 - x^2 - 2\mu x^3$, which satisfies $R_{AB} = -(3/\ell^2)g_{AB}$.
Every coordinate is a pure number, as the footnote beside their \tohors says, so the line element has no $c$ and the length $\ell$ stands in front of the bracket; `DIMENSIONS` declares $t$, $y$, $x$ and $\varphi$ dimensionless and $\ell$ a length.
The chart is one side of the wall: they cut the spacetime along $x = 0$, where $G'(0) = 0$ and $G(0) = 1$ make the extrinsic curvature \extcur equal to $g_{\mu\nu}/\ell$, and glue two copies of $0 \le x \le x_2$ along it, $x_2$ the positive root of $G$.
The period of $\varphi$ is \delphi, $4\pi/|G'(x_2)|$, which the convention writes $\Delta\varphi = 2\pi/(x_2(1 + 3\mu x_2))$.
`two_walls` is the sequel's section 2, its \cads with \hy, \gx and \defell: $ds^2 = 1/(A^2(x - y)^2)\left(H\,dt^2 - dy^2/H + dx^2/G + G\,d\varphi^2\right)$ with $H = -\lambda + ky^2 - 2mAy^3$, $G = 1 + kx^2 - 2mAx^3$ and $\ell_4 = 1/(A\sqrt{\lambda + 1})$.
For $\lambda > 0$ the sequel cuts it along $x = 0$ and $y = 0$ and keeps the quadrant $x \ge 0$, $y \le 0$, both walls of positive tension, and gives $G_3 = AG_4/2$ and $\ell_3^2 = 1/(A^2\lambda)$.
Its section 5.1 says the first paper's chart is this one at $\lambda = 0$, $k = -1$, $\mu = mA$ and $\ell = 1/A$, which `check_chart` holds slot by slot.
`rotating` is the sequel's section 5, its \adsd with its own \hy and \gx: $H = -\lambda + ky^2 - 2mAy^3 - a^2y^4$, $G = 1 + kx^2 - 2mAx^3 + a^2\lambda x^4$ and $W = 1 + a^2x^2y^2$, with the ring singularity \ring, $1/y^2 + a^2x^2 = 0$.
Its tensors print in 40 seconds and verify, so the chart stays, with the tags `rotating black hole` and `angular momentum`; its metric components are written as the line element groups them, $H + a^2y^4G$ over $W$ for $g_{tt}$ and so on, and checked against the computed ones.
The massless limits, the first paper's \hors and the sequel's \nice, are anti-de Sitter space and no chart of the black hole; the metrics on the wall, \bhbrane, \wallgeo and \branegeo, have three dimensions, and `quantum_btz` already prints the wall's chart.

## Step 2. The tidal-charge metric

Naresh Dadhich, Roy Maartens, Philippos Papadopoulos and Vahid Rezania's metric, Physics Letters B 487, 1-6 (2000), arXiv hep-th/0003061, solves the equations induced on a wall with three dimensions of space in a bulk of five, so it is no coordinate system of this spacetime, which has four.
They state that no warped bulk holds it for a tidal charge that is not zero, so no fetched source gives a spacetime of five dimensions to chart, and as a metric of four dimensions it is Reissner and Nordström's form with $Q^2$ replaced by the tidal charge, the geometry `rn_metric` prints.
It therefore enters the History alone.

## Step 3. The checks

`check_chart` in `brane_black_holes.py` holds each chart, exactly at three rational points and for each $k \in \{-1, 0, 1\}$ where the chart has $k$, to $R_{\mu\nu} = -(3/\ell_4^2)g_{\mu\nu}$, and its Weyl tensor to zero at $m = 0$ or $\mu = 0$.
Each wall is held to having extrinsic curvature proportional to its metric: across $x = 0$, and across $y = 0$ in the sequel's charts, the tangential metric changes only through the conformal factor $(x - y)^{-2}$, since $G'(0) = 0$ and $H'(0) = 0$.
The flat wall's chart is held to the first paper's \bhbrane on the wall, $-(1 - 2\mu\ell/r)d\hat t^2 + (1 - 2\mu\ell/r)^{-1}dr^2 + r^2d\varphi^2$ with $\hat t = \ell t$ and $r = -\ell/y$, and to $g_{tt} = 0$ on the horizon $y = -1/(2\mu)$ and on $y = 0$.
The two-wall chart is held to the sequel's \wallgeo on the wall with $\rho = -1/y$, and the rotating chart to \branegeo there and to the two-wall chart at $a = 0$.
The Kretschmann scalars are written by hand as $24/\ell^4 + 48\mu^2(x - y)^6/\ell^4$, $24A^4(\lambda + 1)^2 + 48m^2A^6(x - y)^6$ and, with rotation, $24A^4(\lambda + 1)^2 + 48m^2A^6(x - y)^6(1 - 15a^2x^2y^2 + 15a^4x^4y^4 - a^6x^6y^6)/W^6$, each checked against sympy by `print_charts.py`.
In the null tetrad of the plane of $t$ and $y$ and the plane of $x$ and $\varphi$ the static charts have $\Psi_2$ alone, $-\mu(x - y)^3/\ell^2$ and $-mA^3(x - y)^3$, so they are of Petrov type D, as the sequel's section 5 says of the whole family; `_tools/test_brane_black_holes.py` holds that from the published Weyl tensor.

## Step 4. The drawings

Every drawing of the flat wall is of one hole, $\mu = 1/8$ in units of $\ell$, on the wall $x = 0$, where the plane of $t$ and $y$ is $-(1 + 2\mu y)dt^2 + dy^2/(y^4(1 + 2\mu y))$, Schwarzschild's plane in $r = -\ell/y$ and $\hat t = \ell t$ with $r_s = 2\mu\ell = \ell/4$.
At $\mu = 1/2$ the cones would stand nearly upright across the whole exterior, since $dt/dy = 1/(y^2(1 + 2\mu y)) \ge 6.75$ there; at $\mu = 1/8$ they open to $dt/dy = 0.42$ in the middle of the exterior, so that value is drawn.
The rays keep $t \mp y_*$ with $y_* = -1/y + 2\mu\ln|(1 + 2\mu y)/y|$, which is $r_*/\ell$, the closed form `null_rays.py --verify` checks.
$g^{yy} = y^4(1 + 2\mu y)/\ell^2$ vanishes at $y \to 0$ with the diverging conformal factor, the wall's infinity, so the row declares `where_is_infinity` and no horizon is marked there.
The plane is not totally geodesic in the chart of one side, since the wall bends it, so its captions call the curves null curves.
The two-wall chart is drawn at $A = 1$, $m = 1$, $\lambda = 1/2$ and $k = 1$, where $G = (1 - x)(1 + x + 2x^2)$ has $x_2 = 1$ and $H = -(y + 1/2)(2y^2 - 2y + 1)$ the one root $y_0 = -1/2$, on the axis $x = x_2$, which the rotation holds fixed, so its null curves are null geodesics; the second wall $y = 0$ is its right edge, where $H = -\lambda$ and light arrives in a finite $t$.
Its rays keep $t \mp y_*$ with $y_* = \tfrac25\ln|y + \tfrac12| - \tfrac15\ln(2y^2 - 2y + 1) + \tfrac45\arctan(2y - 1)$.
The rotating chart is drawn on a flat wall, $\lambda = 0$ and $k = -1$, at $mA = 1/8$ and $a = 3/25$, on $x = 0$ with $\varphi$ divided out, where the metric orthogonal to the circles is $(H\,dt^2 - dy^2/H)/(A^2y^2)$ with $H = -y^2(1 + 2mAy + a^2y^2)$; its roots $y = -25/4$ and $y = -100/9$ are Kerr's $r_\pm = 4/25$ and $9/100$ in units of $1/A$, and $g_{tt}$ vanishes on the ergosurface $y = -1/(2mA) = -4$, Kerr's $r = 2m$.
Those two views draw other members of the family than the hole whose moment is embedded, which `exceptions/brane_black_holes.json` and `HIDDEN` record.
The conformal diagram is the flat wall's plane in cells I and II of Schwarzschild's tower with $r_s = \ell/4$, the exterior and the black hole the chart covers, as the first paper's Fig. 1 describes the section $x = 0$; the singularity $y \to -\infty$ lies on $T = \pi/2$ and the wall's infinity on $i^0$ and $\mathscr{I}^\pm$.
The conformal diagram draws the hole whose moment is embedded, on the section the first paper's Fig. 1 describes; the two other charts draw other members of the family, so neither has a conformal view of its own, as Vaidya's outgoing chart has none.
The embedding diagram is the flat wall's moment $t = 0$ on $x = 0$ from the horizon $y = -4$ to $y = -1/3$, $r = 3\ell$, swept as $\varphi = s\phi$ with $s = \Delta\varphi/2\pi = 1/(x_2(1 + 3\mu x_2)) = 0.827$, so its circles have the radius $sr$ and their true circumference $r\Delta\varphi$.
That surface is Aryal, Ford and Vilenkin's equator at $b = s$, whose height $r_s(w\sqrt{1 + k^2w^2} + \mathrm{arsinh}(kw)/k)$ with $w = \sqrt{r/r_s - 1}$ and $k = \sqrt{1 - s^2}$ the generator checks, steep at the horizon and opening far out into a cone of deficit angle $2\pi(1 - s) = 62°$, the cone the History says the wall approaches.
The curvature figure plots the stored Kretschmann scalar on the wall against $r = -\ell/y$ from $0.1\ell$ to $4\ell$ beside Schwarzschild's published scalar at the same $r_s$: on the wall $K = 24/\ell^4 + 12r_s^2/r^6$, which the generator checks at every sample, so $K$ is Schwarzschild's plus the constant of anti-de Sitter space.
The stored Ricci scalar is $-12/\ell^2$ everywhere, so the figure has $K$ alone and its caption states $R$.
