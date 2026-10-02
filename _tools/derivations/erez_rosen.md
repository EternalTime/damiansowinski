# Erez and Rosen's quadrupole

Erez and Rosen's metric in their prolate spheroidal chart and in the chart that is Schwarzschild's at $q = 0$, $x = r/m - 1$, $y = \cos\theta$, is

$$ds^2 = -e^{2\psi}c^2dt^2 + m^2e^{-2\psi}\left(e^{2\gamma}\left(x^2 - y^2\right)\left(\frac{dx^2}{x^2 - 1} + \frac{dy^2}{1 - y^2}\right) + \left(x^2 - 1\right)\left(1 - y^2\right)d\phi^2\right),$$

$$ds^2 = -e^{2\psi}c^2dt^2 + e^{-2\psi}\left(e^{2\gamma}\left(r^2 - 2mr + m^2\sin^2\theta\right)\left(\frac{dr^2}{r^2 - 2mr} + d\theta^2\right) + \left(r^2 - 2mr\right)\sin^2\theta\,d\phi^2\right),$$

with, for $L = \ln\frac{x - 1}{x + 1} = \ln\left(1 - \frac{2m}{r}\right)$,

$$\psi = \frac{L}{2} + \frac{q\left(3y^2 - 1\right)}{8}\left(\left(3x^2 - 1\right)L + 6x\right),$$

$$\gamma = \frac{(1 + q)^2}{2}\ln\frac{x^2 - 1}{x^2 - y^2} - \frac{3q\left(1 - y^2\right)}{2}\left(xL + 2\right) + \frac{9q^2\left(1 - y^2\right)}{64}\left(\left(x^2 - 1\right)\left(x^2 + y^2 - 9x^2y^2 - 1\right)L^2 + 4x\left(x^2 + 7y^2 - 9x^2y^2 - \tfrac{5}{3}\right)L + 4\left(x^2 + 4y^2 - 9x^2y^2 - \tfrac{4}{3}\right)\right).$$

Both charts are written by `_tools/derivations/print_charts.py --metric erez_rosen`, one chart to a run with `--system`, and `verify_metrics.py --system erez_rosen/prolate_spheroidal` and `--system erez_rosen/spherical` check them.

## Step 1. The sources of the two charts

The prolate spheroidal coordinates are Erez and Rosen's own, and Doroshkevich, Zel'dovich, and Novikov (1965, Appendix I), Young and Coulter (1969), and Quevedo (1989) write the metric in them.
The 1959 paper itself was not read: its record is OSTI 4201189, and $\psi$ and $\gamma$ are taken from Quevedo's "Multipolar Solutions" (arXiv:1201.1608) and Boshkayev and coauthors (2019), who print them identically.
The chart in $r$ and $\theta$ is how Bini and coauthors (2013) and Memmen and Perlick (2021) pass to "Schwarzschild-like coordinates", and at $q = 0$ it is Schwarzschild's published metric with $r_s = 2m$, which the check holds it to.
$\psi = -Q_0(x) - q\,P_2(y)Q_2(x)$ in Legendre's functions, Schwarzschild's potential plus the quadrupole harmonic that dies away at infinity.

## Step 2. The coefficient that was misprinted

Doroshkevich, Zel'dovich, and Novikov's Appendix I, in the English translation, prints the first term of $\gamma$ as $\tfrac{1}{2}(1 + q + q^2)\ln\frac{\lambda^2 - 1}{\lambda^2 - \mu^2}$, and every other term as above.
Weyl's quadrature for $\gamma$,

$$\partial_x\gamma = \frac{1 - y^2}{x^2 - y^2}\left(x\left(x^2 - 1\right)\left(\partial_x\psi\right)^2 - x\left(1 - y^2\right)\left(\partial_y\psi\right)^2 - 2y\left(x^2 - 1\right)\partial_x\psi\,\partial_y\psi\right),$$

$$\partial_y\gamma = \frac{x^2 - 1}{x^2 - y^2}\left(y\left(x^2 - 1\right)\left(\partial_x\psi\right)^2 - y\left(1 - y^2\right)\left(\partial_y\psi\right)^2 + 2x\left(1 - y^2\right)\partial_x\psi\,\partial_y\psi\right),$$

holds with $\tfrac{1}{2}(1 + q)^2$ and fails with $\tfrac{1}{2}(1 + q + q^2)$: the two differ by $\tfrac{1}{2}q\ln\frac{x^2 - 1}{x^2 - y^2}$, which is not constant.
`erez_rosen_check` holds the published $\gamma$ to both quadratures and `_tools/test_erez_rosen.py` holds the misprinted coefficient to failing them.

## Step 3. Two names each chart holds

Written out, a curvature component is a polynomial of the sixth degree in $L$, so each chart declares $L$, $\psi$ and $\gamma$ among its parameters and holds $\psi$ and $\gamma$ as functions of its two coordinates while the tensors are built, `HELD` in `verify_metrics.py`, as the first Morgan-Morgan disc's are.
`erez_rosen_reduce` writes every derivative of $\gamma$ by the quadrature and the second derivative of $\psi$ along the radial coordinate by Laplace's equation, which leaves the other derivatives of $\psi$ with no relation among them, so the Ricci tensor is exactly zero.
Each value is printed as a polynomial in those derivatives, with every coefficient in the sums of the line element.
The spherical chart's values are factored in $\cos\theta$: the metric is even in $\theta$, so a value is a rational function of $r$, $m$ and $\cos\theta$, times $\sin\theta$ or not, once each derivative of $\psi$ taken an odd number of times along $\theta$ is counted as odd, and in $\cos\theta$ the sum $r^2 - 2mr + m^2\sin^2\theta = (r - m - m\cos\theta)(r - m + m\cos\theta)$ comes out as the factor it is.
The checker took 350 seconds on the prolate spheroidal chart with the definitions written out on both sides of each comparison and left its Kretschmann scalar unchecked at 120 seconds, so `RATES` declares the first derivatives of $\psi$, polynomials in $L$, and of $\gamma$, the quadrature; the reader checks each against the definition, and the chart then checks whole in 91 seconds.

## Step 4. What the check holds

`erez_rosen_check` holds each chart to a vanishing Ricci tensor, $\psi$ to Laplace's equation, $\gamma$ to both quadratures and to vanishing on the axis, and both to Schwarzschild's $\tfrac{1}{2}\ln f$ and $\tfrac{1}{2}\ln\frac{x^2 - 1}{x^2 - y^2}$ at $q = 0$.
On the axis Weyl's $z = mx$, and Geroch's potential $-\tanh\psi$ there is $m/z + \tfrac{2}{15}qm^3/z^3 + \dots$, so the mass is $m$ and the quadrupole moment $\tfrac{2}{15}qm^3$, which the check holds the prolate spheroidal chart to.
The spherical chart's $\psi$ and $\gamma$ are held to the prolate spheroidal chart's at six random points in forty digits, the chart to being the prolate spheroidal one pulled back, and at $q = 0$ to Schwarzschild's published metric.

## Step 5. The singular surface

Toward $x = 1$, $2\psi \to (1 + qP_2(y))\ln(x - 1)$ and $2\gamma \to (1 + qP_2(y))^2\ln(x - 1)$, so near each point of the surface the field is Zipoy and Voorhees's with $\delta = 1 + qP_2(y)$, which is Doroshkevich, Zel'dovich, and Novikov's Appendix II.
On the axis $\delta = 1 + q$ and on the equator $\delta = 1 - q/2$, so $g_{tt}$ vanishes on the whole surface for $-1 < q < 2$.
The Kretschmann scalar goes as $(x - 1)^{2\delta - 4}$ on the axis and as $(x - 1)^{-2(\delta^2 - \delta + 1)}$ off it: at $q = 1$ it is finite on the axis, $3e^6/(4m^4)$ at $x = 1$, and diverges as $(x - 1)^{-3/2}$ on the equator, and at $q = -1/2$ it diverges as $(x - 1)^{-3}$ and as $(x - 1)^{-21/8}$.
A ray on the axis has $c\,dt/dx = m\,e^{-2\psi} \sim (x - 1)^{-1-q}$, which is integrable only for $q < 0$, and a ray in the equatorial plane has $c\,dt/dx = m\,x\,e^{\gamma - 2\psi}/\sqrt{x^2 - 1} \sim (x - 1)^{q^2/8 - 1}$, integrable for every $q \neq 0$: those are their times $(\lambda_0 - 1)^{-q}$ and $(\lambda_0 - 1)^{q^2/8}$.
The circle about the axis in the equatorial plane has radius $m(x + 1)f^{q(3x^2 - 1)/8}e^{3qx/4}$, which goes as $(x - 1)^{q/4}$: it shrinks to zero for the prolate $q > 0$, their cucumber, and grows without bound for the oblate $q < 0$.

## Step 6. The diagrams

Every drawing is at $m = 1$ for the prolate $q = 1$ and the oblate $q = -1/2$, on the two totally geodesic planes, the axis and the equatorial plane, in each chart.
The tortoise coordinate has no closed form, so `_erez_rosen_star` in `null_rays.py` takes it by quadrature from the closed forms in floats, in the variable $u = (x - 1)^{1/n}$ that makes the integrand finite at $x = 1$, with $n = 2$ on the axis and $n = 8/q^2$ in the plane, and `null_rays.py --verify` compares the rays with it.
At $q = -1/2$ the equatorial integrand goes as $(x - 1)^{-31/32}$: light from $x = 1.1$ takes $69.9\,m/c$ to reach the singularity, so that plane's conformal maps take $\ell = 40\,m$ where the others take $4\,m$.
The conformal diagrams are the whole diamond for the prolate axis, its left edges the edge of the chart with the Kretschmann scalar checked to be $3e^6/(4m^4)$ there, and Minkowski's triangle with a timelike singularity for the other three planes.
The embedding diagram is the equatorial plane at $t = 0$: at $q = 1$ it runs in to $r = 2.0279\,m$, a circle of radius $1.371\,m$, and at $q = -1/2$ it has a neck of radius $2.2867\,m$ at $r = 2.1297\,m$ and runs in to $r = 2.0288\,m$.
