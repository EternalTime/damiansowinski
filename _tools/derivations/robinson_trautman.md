# Robinson-Trautman, an isolated body radiating on expanding fronts

Robinson and Trautman's metric is

$$ds^2 = -2H\,c^2du^2 - 2c\,du\,dr + \frac{r^2}{P^2}\left(dx^2 + dy^2\right),$$

with $u$ a retarded time, $r$ an affine length along the rays $u, x, y$ constant, and $x$ and $y$ pure numbers on each wave front.
Its two charts are written by `_tools/derivations/print_charts.py --metric robinson_trautman` in about 21 seconds, and `verify_metrics.py --system robinson_trautman/<chart>` checks each in under three.

## Step 1. The charts

The chart of $u$, $r$, $x$ and $y$ is the one Robinson and Trautman wrote in 1960, with their $\sigma$, $\rho$, $\xi$, $\eta$ and $p$; the complex coordinate of later work is $\zeta = (x + iy)/\sqrt2$, so that $2\,d\zeta\,d\bar\zeta = dx^2 + dy^2$, and the checker's reader takes real coordinates only.
The axisymmetric chart is the one the numerical work uses, with the fronts' metric $r^2f^{-2}\left(d\theta^2 + \sin^2\theta\,d\phi^2\right)$ and $f = f(u,\theta)$, Macedo and Saa's $Q$.
The two are one metric: $x + iy = 2\tan(\theta/2)\,e^{i\phi}$ takes $(dx^2 + dy^2)/(1 + (x^2 + y^2)/4)^2$ to the round sphere, so $P = f\left(1 + (x^2 + y^2)/4\right)$, which the printer checks with $t = \tan(\theta/2)$, where every entry is rational.

## Step 2. The field equations

Both charts leave $H$ and $P$, or $f$, free, so no component assumes a field equation.
The Gaussian curvature of a front at $r = 1$ is

$$K = P^2\left(\partial_x^2 + \partial_y^2\right)\ln P = f^2\left(1 + \frac{\partial_\theta\left(\sin\theta\,\partial_\theta\ln f\right)}{\sin\theta}\right),$$

and with $m = GM/c^2$ the vacuum $H$ is $2H = K - 2r\,\partial_u\ln P - 2m/r$.
Before anything is written the printer puts that $H$ into the Ricci tensor it computed and checks that every component vanishes but

$$R_{uu} = \frac{\Delta K + 12m\,\partial_u\ln P}{2r^2},$$

with $\Delta$ the Laplacian of the front at $r = 1$, so the Robinson-Trautman equation $\Delta K + 12m\,\partial_u\ln P = 0$ is the one field equation left.
Every $\partial_u$ is taken in the chart $x^0 = cu$.
It also checks that the round front, $P = 1 + (x^2 + y^2)/4$ or $f = 1$, has $K = 1$ and the Kretschmann scalar $48m^2/r^6$, Schwarzschild's in outgoing Eddington-Finkelstein coordinates.

## Step 3. Printing

Each value is printed with the front's function first and collected by its powers, and with `flip = False`, so the Christoffel symbols read $-\partial_uP/P$ and the Ricci tensor keeps the order $P^4$, $P^2$, $P$ of its terms.

## Step 4. The fronts drawn

No vacuum member with smooth round fronts is known in closed form but Schwarzschild's, so the drawings take a declared first front and solve the equation.
`FrontSolver` in `null_rays.py` writes $f(u,\theta) = \sum_l b_l(u)P_l(\cos\theta)$ with forty Legendre polynomials, as Macedo and Saa do, in units $G = c = m = 1$, where

$$\partial_uf = -\frac{f^3\,\Delta_0K}{12},\qquad K = f^2 + f\,\Delta_0f - (1 - x^2)\left(\partial_xf\right)^2,\qquad x = \cos\theta,$$

and $\Delta_0$ is the Laplacian of the unit sphere.
The right hand side is a polynomial in $x$, projected onto the $P_l$ by Gauss-Legendre quadrature at 160 nodes, and the modes are integrated by Radau's implicit method with the Jacobian written out, since differences of a rate that carries the fourth derivative of $f$ lose to rounding the digits Newton's iteration needs.
The first front is Macedo and Saa's prolate one,

$$f(0,\theta)^2 = f_0^2\left(1 - \epsilon^2\cos^2\theta\right),\qquad f_0^2 = \frac{1}{2\epsilon}\ln\frac{1 + \epsilon}{1 - \epsilon},$$

at $\epsilon = 4/5$, which gives every front the area $4\pi r^2$, so the last front has $f = 1$ and the mass parameter $m$ is the mass of the Schwarzschild black hole left.
The equation has no solution toward the past, so nothing is drawn before $u = 0$.

`front_checks` holds the fronts to what their construction does not use, at points all over $u$, $r$ and $\theta$: every published component of the axisymmetric chart's Ricci tensor vanishes on them, with every derivative of $f$ and $H$ taken from the series itself, and the published Kretschmann scalar is $48m^2/r^6$, as it is for every vacuum metric whose Weyl tensor has only $\Psi_2 = -m/r^3$, $\Psi_3$ and $\Psi_4$.
The area of a front stays $4\pi r^2$; the Bondi mass $M_B = (m/2)\int_{-1}^{1}f^{-3}\,dx$ starts at Macedo and Saa's closed form $m/(1 - \Delta)$, with

$$\Delta = 1 - \sqrt{\frac{1 - \epsilon^2}{8\epsilon^3}\ln^3\frac{1 + \epsilon}{1 - \epsilon}},$$

which is $1.0357\,m$ at $\epsilon = 4/5$, falls at every step and ends at $m$; and the last mode to die, $l = 2$, falls as $e^{-2cu/m}$, Foster and Newman's rate.

The fronts are even about the equator, so $\partial_\theta H$ vanishes on the axis and on the equator, and $\Gamma^\theta{}_{uu} = f^2\partial_\theta H/r^2$ with it: light launched along either plane of $u$ and $r$ stays on it, and the null curves drawn there are null geodesics.
The outgoing rays are $u$ constant and the ingoing ones obey $dr/d(cu) = -H$.
A published value that divides by $\sin\theta$, as the Kretschmann scalar does, is read $10^{-3}$ beside the axis.

The embedding diagram is the front itself, the surface of $\theta$ and $\phi$ at one $u$ and one $r$, a surface of revolution with circles of radius $r\sin\theta/f$ and

$$\frac{dz}{d\theta} = \frac{r}{f}\sqrt{1 - \left(\cos\theta - \frac{\sin\theta\,\partial_\theta f}{f}\right)^2}.$$

The fronts of one $u$ differ only in size, so each is drawn with its own $r$ as the unit, and on the spacetime diagrams a moment of the embedding diagram is the whole outgoing ray $u = u_k$.
