# Schwarzschild's black hole in a tidal field

Schwarzschild's black hole in the quadrupole field of distant matter, in the prolate spheroidal chart, in the chart that is Schwarzschild's at $q = 0$, $x = r/m - 1$, $y = \cos\theta$, and in Weyl's canonical chart, $\rho = m\sqrt{(x^2 - 1)(1 - y^2)}$, $z = mxy$, is

$$ds^2 = -\frac{x - 1}{x + 1}e^{2U}c^2dt^2 + m^2\left(x + 1\right)^2e^{-2U}\left(e^{2V}\left(\frac{dx^2}{x^2 - 1} + \frac{dy^2}{1 - y^2}\right) + \left(1 - y^2\right)d\phi^2\right),$$

$$ds^2 = -\left(1 - \frac{2m}{r}\right)e^{2U}c^2dt^2 + e^{-2U}\left(e^{2V}\left(\frac{dr^2}{1 - 2m/r} + r^2d\theta^2\right) + r^2\sin^2\theta\,d\phi^2\right),$$

$$ds^2 = -e^{2\psi}c^2dt^2 + e^{-2\psi}\left(e^{2\gamma}\left(d\rho^2 + dz^2\right) + \rho^2d\phi^2\right),$$

with

$$U = \frac{q}{4}\left(3x^2 - 1\right)\left(3y^2 - 1\right), \qquad V = -3qx\left(1 - y^2\right) - \frac{9q^2}{16}\left(x^2 - 1\right)\left(1 - y^2\right)\left(9x^2y^2 - x^2 - y^2 + 1\right),$$

$$\psi = \frac{1}{2}\ln\frac{x - 1}{x + 1} + U, \qquad \gamma = \frac{1}{2}\ln\frac{x^2 - 1}{x^2 - y^2} + V.$$

The three charts are written by `_tools/derivations/print_charts.py --metric distorted_schwarzschild`, one chart to a run with `--system`, and `verify_metrics.py --system distorted_schwarzschild/<chart>` checks them.

## Step 1. The sources of the three charts

The solution is Appendix IV of Doroshkevich, Zel'dovich, and Novikov (1965), in the prolate spheroidal $\lambda$ and $\mu$ of their Appendix I, read from the English translation in the journal's own archive.
Their $\psi$ and $\gamma$ are the two functions above, and `distorted_schwarzschild_check` and `_tools/test_distorted_schwarzschild.py` hold the published $U$ and $V$ to them.
The line element with Schwarzschild's part taken out of Weyl's two functions, $U$ and $V$ for what the distant matter adds, is Shoom, Walsh, and Booth's (1) (2016), after Geroch and Hartle (1982); the chart in $r$ and $\theta$ is Fairhurst and Krishnan's (3.7) (2001); and Weyl's chart is the one Geroch and Hartle and Fairhurst and Krishnan's (2.1) to (2.3) work in.
Geroch and Hartle's paper itself was not read: its record is doi:10.1063/1.525384, and what is said of it comes from its abstract and from Fairhurst and Krishnan, Frolov and Shoom (2007), and Shoom, Walsh, and Booth.

## Step 2. Which quadrupole parameter

Two conventions are in print.
Doroshkevich, Zel'dovich, and Novikov's $q$ multiplies $P_2(x)P_2(y)$, and it is Frolov and Shoom's $a_2$, whose potential is $\sum_n a_nP_n(\cos\psi)P_n(\cos\theta)$.
Shoom, Walsh, and Booth expand in $R^nP_n(xy/R)$ with $R = \sqrt{x^2 + y^2 - 1}$, and $R^2P_2 = \tfrac{1}{2}(3x^2y^2 - x^2 - y^2 + 1) = \tfrac{2}{3}P_2(x)P_2(y) + \tfrac{1}{3}$, so their quadrupole moment is $\tfrac{3}{2}q$ and their $U$ differs by a constant.
A constant added to $U$ rescales $t$ and $m$ and is no new solution.
The published $U$ carries none beyond what $P_2(x)P_2(y)$ has, so $U = q$ at the poles of the horizon, Frolov and Shoom's $u_0$.

## Step 3. The surface gravity, the area and the mass

On the horizon $V = 2U - 2q$, which is checked, so $g_{tt}g_{rr} \to -e^{-4q}$ times Schwarzschild's and
$$\kappa^2 = \lim_{r \to 2m}\frac{\left(\partial_rg_{tt}\right)^2}{-4g_{tt}g_{rr}} = \frac{e^{4q}}{16m^2}$$
at every $\theta$: the surface gravity is $e^{2q}/4m$, Frolov and Shoom's (21) with $u_0 = q$, and the exponent is $+2q$.
The horizon's metric is $4m^2e^{-2U}\left(e^{2V}d\theta^2 + \sin^2\theta\,d\phi^2\right)$ with $U = \tfrac{q}{2}(3\cos^2\theta - 1)$ and $V = -3q\sin^2\theta$, its area element is $4m^2e^{-2q}\sin\theta$, its area $16\pi m^2e^{-2q}$, and $\kappa A/4\pi = m$, the Komar mass Fairhurst and Krishnan evaluate on the horizon.
Its Gaussian curvature is
$$\frac{e^{2q + 3q\sin^2\theta}}{4m^2}\left(1 + 3q - 15q\cos^2\theta - 18q^2\sin^2\theta\cos^2\theta\right),$$
Frolov and Shoom's (48), $(1 - 12q)e^{2q}/4m^2$ at the poles and $(1 + 3q)e^{5q}/4m^2$ on the equator, and the Kretschmann scalar on the horizon is twelve times its square.
The English translation of Doroshkevich, Zel'dovich, and Novikov prints another polynomial for this curvature, with $e^q$ in front, which is not used.
`_tools/test_distorted_schwarzschild.py` holds the published metric to each of these.

## Step 4. Two names each chart holds, and what is left standing

Each chart declares its two functions among its parameters and the checker holds them as functions of the chart's two coordinates, `HELD` in `verify_metrics.py`, as Erez and Rosen's are.
Weyl's quadrature for $V$ divides by $x^2 - y^2$, which the line element with Schwarzschild's part taken out does not hold.
So in the prolate spheroidal and spherical charts `distorted_schwarzschild_reduce` writes only the second derivative of $U$ along the radial coordinate, by Laplace's equation, and that of $V$, by the one field equation of second order,
$$\left(x^2 - 1\right)\partial_x^2V + x\,\partial_xV + \left(1 - y^2\right)\partial_y^2V - y\,\partial_yV = -\left(x^2 - 1\right)\left(\partial_xU\right)^2 - \left(1 - y^2\right)\left(\partial_yU\right)^2 - 2\partial_xU,$$
and leaves the first derivatives of $V$ standing, so that no value divides by $x^2 - y^2$.
A value that vanishes on the quadrature as well is written as zero, tried at two random points before it is asked exactly, which is how the Ricci tensor is zero.
The Weyl tensor built apart carries the Ricci tensor's terms, which vanish only on the quadrature, so `distorted_schwarzschild_vacuum` checks each of its components against the Riemann tensor's there and writes it as the Riemann tensor is written.
The Kretschmann scalar is written as $16(a^2 + ab + b^2 + c^2)$ with $a = R^{tx}{}_{tx}$, $b = R^{ty}{}_{ty}$ and $c^2 = R^{tx}{}_{ty}R^{ty}{}_{tx}$, the square of the electric part of a static vacuum's Weyl tensor, checked against the scalar built from every component.
The spherical chart's values are factored in $\cos\theta$ as Erez and Rosen's are, and a sine left under a sum is taken into it term by term, so that $(\partial_\theta U\sin\theta - \cos\theta)/\sin\theta$ is written $\partial_\theta U - \cot\theta$.
Weyl's chart writes $\gamma$'s derivatives by the quadrature, which divides by nothing there, with `morgan_morgan_reduce`, and names the prolate spheroidal $x$ and $y$ as functions of $\rho$ and $z$ to define $\psi$ and $\gamma$ by.
With those roots written out and differentiated twice the checker left its Kretschmann scalar unchecked at 120 seconds, after 445 in all, so `RATES` declares the first derivatives of $\psi$, rational in $\rho$, $z$, $x$ and $y$, and of $\gamma$, the quadrature, and the chart then checks whole in 134 seconds; the other two take 9 and 13.

## Step 5. What the check holds

`distorted_schwarzschild_check` holds each of the first two charts, with the names written out, to a Ricci tensor that is exactly zero, $U$ to Laplace's equation, $V$ to both quadratures, to the equation of second order, to vanishing on the axis and to $2U - 2q$ on the horizon, and both to vanishing at $q = 0$.
The prolate spheroidal chart is held to being Weyl's line element with $\psi$ and $\gamma$ as above and to Doroshkevich, Zel'dovich, and Novikov's printed functions, and the spherical chart to the prolate spheroidal chart's functions, to being it pulled back, and at $q = 0$ to Schwarzschild's published metric with $r_s = 2m$.
Weyl's chart is held at six random points in forty digits to Laplace's equation, to the quadrature, to $\gamma = 0$ on the axis beyond the rod and to the prolate spheroidal chart's functions, and that chart to being Weyl's pulled back.

## Step 6. The diagrams

Every drawing is at $m = 1$ for the oblate $q = 1/12$, the largest $q$ whose horizon stands in flat space, and the prolate $q = -1/12$, Frolov and Shoom's two shapes, on the two totally geodesic planes, the axis and the equatorial plane.
On either plane the metric is conformal to $-c^2dt^2 + dx_*^2$ with $dx_*/dx = m\,e^{V - 2U}(x + 1)/(x - 1)$, where $V - 2U = -q(3x^2 - 1)$ on the axis and $\tfrac{9}{16}q^2(x^2 - 1)^2 + \tfrac{q}{2}(3x^2 - 6x - 1)$ in the plane.
Both have the one pole $x = 1$, of residue $2m\,e^{-2q}$, so `_distorted_star` in `null_rays.py` takes $x_*$ as $2e^{-2q}\ln|x - 1|$ plus the quadrature of a smooth function, zero at the singularity $x = -1$, and `null_rays.py --verify` compares the rays with it.
The prolate spheroidal and spherical charts' components are analytic through the horizon, so their spacetime diagrams are drawn on to the singularity with the future taken from the ingoing rays, as Schwarzschild's own chart is drawn; Weyl's chart ends on the horizon, the end $z = m$ of the rod on its axis and the point $\rho = 0$ of its equatorial plane.
The conformal diagrams are a `Tower` of the one root with that $x_*$ and $\kappa = e^{2q}/4m$, `DistortedTower` in `conformal.py`: Kruskal and Szekeres's hexagon on three of the four planes, and on the oblate axis, where $x_*$ tends to $R = 0.3443\,m$, the same cells with the edge $x \to \infty$ the timelike curve $\tan p\tan q = -e^{2\kappa R}$, as Schwarzschild-anti-de Sitter's boundary is drawn.
Light reaches that edge after an infinite affine distance and the Kretschmann scalar grows without bound toward it; on the other three planes it falls to zero far out.
The embedding diagrams are the equatorial plane at $t = 0$ for each shape, through the throat $r = 2m$ of radius $2m\,e^{q/2}$ into the other exterior, and the horizon itself.
At $q = 1/12$ the plane lies level at $r = 2.5388\,m$ and is drawn on in Minkowski space, as Schwarzschild-anti-de Sitter's is, and the horizon is $2.085\,m$ across its equator and $1.225\,m$ from its centre to a pole, flat at the poles; at $q = -1/12$ the plane stands in flat space as far as it is drawn and the horizon is $1.918\,m$ and $2.652\,m$.
