# The Kaluza-Klein monopole

Sorkin's and Gross and Perry's monopole is Euclidean Taub-NUT space, self dual, with $-c^2dt^2$ added: a static vacuum solution in five dimensions.
It is the collection's one spacetime whose fifth dimension is a circle, so this note records the charts, the period of that circle and what each diagram draws.

## Step 1: Gross and Perry's chart

$ds^2 = -c^2dt^2 + V(dr^2 + r^2d\theta^2 + r^2\sin^2\theta\,d\phi^2) + V^{-1}(dx_5 + 4m(1 - \cos\theta)\,d\phi)^2$ with $V = 1 + 4m/r$.
Gross and Perry write $V$ for the reciprocal of this function; the line element here spells the function out, so no letter is needed.
The fifth coordinate is written $x_5$, since the checker's reader takes a subscript in a coordinate's name and not a superscript, which it reads as a power.
The one parameter $m$ is a length, and $G$ appears nowhere in the line element.

## Step 2: the period of $x_5$

Near $r = 0$, with $r = s^2/16m$, the metric of a slice of constant $t$ is $ds^2 + \tfrac{s^2}{4}\left(d\theta^2 + \sin^2\theta\,d\phi^2 + (dx_5/4m + (1 - \cos\theta)\,d\phi)^2\right)$ to leading order.
With $\psi = -\phi - x_5/4m$ the bracket is $d\theta^2 + \sin^2\theta\,d\phi^2 + (d\psi + \cos\theta\,d\phi)^2$, four times the round metric of the unit 3-sphere in Euler angles when $\psi$ has period $4\pi$.
So the slice is flat $\mathbb{R}^4$ at the nut exactly when $x_5$ has period $16\pi m$, the fifth dimension far away is a circle of radius $R = 8m$, and with any other period the nut is the apex of a cone.
The same period makes the shift $x_5 \to x_5 + 8m\phi$, which moves the Dirac string from $\theta = \pi$ to $\theta = 0$, a well defined change of chart.

## Step 3: the other two charts

The Hopf chart keeps $r$ and trades $x_5$ for $\psi = -\phi - x_5/4m$: $ds^2 = -c^2dt^2 + V(dr^2 + r^2d\Omega^2) + 16m^2V^{-1}(d\psi + \cos\theta\,d\phi)^2$.
The Taub-NUT chart moves the radius as well, $\rho = r + 2m$: $V = (\rho + 2m)/(\rho - 2m)$ and $Vr^2 = \rho^2 - 4m^2$, which is Hawking's self dual Taub-NUT instanton with NUT parameter $n = 2m$, $\frac{\rho + n}{\rho - n}d\rho^2 + (\rho^2 - n^2)d\Omega^2 + 4n^2\frac{\rho - n}{\rho + n}(d\psi + \cos\theta\,d\phi)^2$.
`print_charts.kaluza_klein_pullback` pulls Gross and Perry's metric back through each map and compares it with the chart's own, slot by slot.

Gibbons and Hawking's Cartesian form, $V\,d\mathbf{x}\cdot d\mathbf{x} + V^{-1}(dx_5 + \mathbf{A}\cdot d\mathbf{x})^2$ with $\mathbf{A}\cdot d\mathbf{x} = 4m(x\,dy - y\,dx)/(r(r + z))$, was tried and is not published.
Its metric is checked to be Gross and Perry's pulled back, but with $r = \sqrt{x^2 + y^2 + z^2}$ in every component sympy had not finished the Riemann tensor after ten minutes on 1 October 2026, where each polar chart takes five seconds.

## Step 4: the curvature

All three charts are Ricci flat, so the Weyl tensor is the Riemann tensor.
The Kretschmann scalar is $K = 384m^2/(r + 4m)^6$, checked against sympy as it is written: $3/(32m^4)$ at the nut, finite, and falling as $r^{-6}$.

## Step 5: the diagrams

On the plane of $t$ and $r$ at fixed angles the metric is $-c^2dt^2 + V\,dr^2$, and no Christoffel symbol with an upper angle and lower indices among $t$ and $r$ is nonzero, so the plane is totally geodesic and its null curves are null geodesics.
With $r_* = \int_0^r\sqrt{1 + 4m/r'}\,dr' = \sqrt{r(r + 4m)} + 4m\,\mathrm{arsinh}\sqrt{r/4m}$, the proper distance from the nut, the rays are $ct \pm r_*$ constant, which `null_rays.py --verify` checks.
$r_*$ runs from $0$ without bound, so $p, q = \arctan((ct \mp r_*)/\ell)$ give Minkowski's triangle with the nut a regular centre and no horizon; the drawing takes $\ell = 4m$.
The planes are taken on the half axis $\theta = 0$, where the embedding diagram's surface lies.

The embedding diagram is the surface of $r$ and $x_5$ at $\theta = 0$ and one $t$, where the potential $4m(1 - \cos\theta)$ vanishes: $V\,dr^2 + dx_5^2/V$.
Its angle is $x_5/8m$, so `Slice` sweeps $x_5 = 8m\phi$ with the chart's own $\phi$, which moves nothing at $\theta = 0$.
The circle at $r$ has radius $\rho = 8m\sqrt{r/(r + 4m)}$, and $d\rho/ds = 16m^2/(r + 4m)^2 \le 1$ with $s$ the proper distance, so the surface embeds in flat space everywhere: a cigar, smooth at the nut where $d\rho/ds = 1$, a cylinder of radius $8m$ far away.
Its height is the quadrature of $dz/dr = \sqrt{(r^3 + 16mr^2 + 96m^2r + 256m^3)/(r + 4m)^3}$, which is $2$ at the nut and tends to $1$.

The horizontal lift of the equator of the base, the surface $\theta = \pi/2$ at constant $\psi$, is not drawn: the Hopf connection's holonomy round the equator is half the fibre, so the circle of $\phi$ closes only after two turns and the surface covers the base's plane twice.
