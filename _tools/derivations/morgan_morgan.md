# The Morgan-Morgan discs

Morgan and Morgan's static discs of counterrotating dust are members of Weyl's class,

$$ds^2 = -e^{2\psi}c^2dt^2 + e^{-2\psi}\left(e^{2\gamma}\left(d\rho^2 + dz^2\right) + \rho^2d\phi^2\right),$$

and in the oblate spheroidal coordinates $\rho = a\sqrt{(1 + \xi^2)(1 - \eta^2)}$, $z = a\xi\eta$, in which the disc of radius $a$ is the surface $\xi = 0$,

$$ds^2 = -e^{2\psi}c^2dt^2 + a^2e^{-2\psi}\left(e^{2\gamma}\left(\xi^2 + \eta^2\right)\left(\dfrac{d\xi^2}{1 + \xi^2} + \dfrac{d\eta^2}{1 - \eta^2}\right) + \left(1 + \xi^2\right)\left(1 - \eta^2\right)d\phi^2\right).$$

Both charts are written by `_tools/derivations/print_charts.py --metric morgan_morgan`, in about a minute and a half, and `verify_metrics.py --system morgan_morgan/weyl` and `--system morgan_morgan/oblate_spheroidal` check them in two seconds and in forty.

## Step 1. The sources of the two charts

Weyl's chart is how Lynden-Bell and Pineault (1978), Bičák, Lynden-Bell, and Katz (1993), and Semerák (2002) write the discs, and it holds for every member, so it leaves $\psi$ and $\gamma$ free, as the double Kerr chart leaves its functions, and no component assumes a field equation; `OVERRULED` in `metric_tags.py` sets its tag `vacuum`.
The oblate spheroidal coordinates are Morgan and Morgan's own, after Hunter (1963): their erratum of 1970 writes the coefficient $\tfrac{1}{2}\left(\xi/(1 + \xi^2) - \cot^{-1}\xi\right)$.
Their potential is $\phi$; the collection writes Weyl's two functions $\psi$ and $\gamma$, as on the Curzon-Chazy page, since $\phi$ is the azimuth.
The paper itself is behind the publisher's wall and was not read; its abstract, its two errata and its reference list were, and every formula here was derived and then compared with the open papers named below.

## Step 2. The first disc

The disc of order $n$ is defined by its Newtonian density $(2n + 1)M(1 - \rho^2/a^2)^{n - 1/2}/2\pi a^2$ (Semerák 2002, eq. 80; González and Reina 2006, eq. 20).
Its potential is $-(m/a)\sum_{k \le n} C_{2k}\,q_{2k}(\xi)P_{2k}(\eta)$ with $q_{2k}(\xi) = i^{2k+1}Q_{2k}(i\xi)$, and matching $\partial_z\psi = (1/a\eta)\partial_\xi\psi$ on the disc to the density gives $C = (1, 1)$ for $n = 1$:

$$\psi = -\dfrac{3m}{4a}\left(\left(1 + \eta^2 - \xi^2 + 3\xi^2\eta^2\right)\mathrm{arccot}\,\xi - \xi\left(3\eta^2 - 1\right)\right),$$

which is González and Reina's (24a) and (25a) and Gutiérrez-Piñeres and González's (47a).
The quadrature for $\gamma$ in these coordinates is the pair of equations the chart's parameters state, Kofroň, Kotlařík, and Semerák's (57) and (58), and integrating the second from the axis $\eta = 1$ gives

$$\gamma = -\dfrac{9m^2\left(1 - \eta^2\right)}{16a^2}\left(\left(1 + \xi^2\right)\left(9\xi^2\eta^2 + \eta^2 - \xi^2 - 1\right)\alpha^2 - 2\xi\left(9\xi^2\eta^2 + 7\eta^2 - \xi^2 + 1\right)\alpha + 9\xi^2\eta^2 + 4\eta^2 - \xi^2 + 4\right),\qquad \alpha = \mathrm{arccot}\,\xi,$$

their (A1) to (A5).
`morgan_morgan_check` holds $\psi$ to Laplace's equation, to the density and to $-m/r$ far away, $\gamma$ to both quadratures, to vanishing on the axis and to the Curzon-Chazy particle's $-m^2\sin^2\theta/2r^2$ far away, the Ricci tensor to vanishing, and the chart to being Weyl's pulled back.

## Step 3. Two names the chart holds

Written out, a curvature component is a polynomial of the sixth degree in $\mathrm{arccot}\,\xi$, so the chart declares $\alpha$, $\psi$ and $\gamma$ among its parameters and holds $\psi$ and $\gamma$ as functions of $\xi$ and $\eta$ while the tensors are built, `HELD` in `verify_metrics.py`.
`morgan_morgan_reduce` writes every derivative of $\gamma$ by the quadrature and $\partial_\xi^2\psi$ by Laplace's equation, which leaves $\partial_\xi\psi$, $\partial_\eta\psi$, $\partial_\xi\partial_\eta\psi$ and $\partial_\eta^2\psi$ with no relation among them, so a value that vanishes for a harmonic $\psi$ is exactly zero and the Ricci tensor is printed as vanishing.
Each value is printed as a polynomial in those four, every coefficient in the three sums of the line element, $1 + \xi^2$, $1 - \eta^2$ and $\xi^2 + \eta^2$.
The checker writes the definitions back in on both sides of every comparison, so what it checks is the first disc itself.
`null_rays.load` and `metric_tags.py` pass `HELD` to their readers as the checker does, so that a published value which writes $\partial_\xi\psi$ is read by every generator.

## Step 4. The disc

Israel's junction conditions across $z = 0$ give a layer with no stress along $\rho$, by the quadrature $\partial_z\gamma = 2\rho\,\partial_\rho\psi\,\partial_z\psi$, and with the ratio of the stress round the axis to the energy density $V^2/c^2 = \rho\,\partial_\rho\psi/(1 - \rho\,\partial_\rho\psi)$, the speed of two streams of dust, which is Lynden-Bell and Pineault's (3.6).
On the first disc $\psi = -(3\pi m/8a)(2 - \rho^2/a^2)$, so $\rho\,\partial_\rho\psi = (3\pi m/4a^3)\rho^2$, and the streams stay slower than light out to the rim for $m/a < 2/(3\pi) = 0.2122$, where light from the centre has the redshift $e^{1/2} - 1 = 0.6487$.
The same limit taken for each order rises with the order, to $1.5803$ for a Gaussian disc, the bound of Morgan and Morgan's abstract, which `_tools/test_morgan_morgan.py` reproduces.
The Kretschmann scalar is finite on the disc and on the axis and grows as $2.58\,a^{-4}/(\xi^2 + \eta^2)$ toward the rim at $m = a/5$, the inverse of the distance from it, since $\rho - a = a\xi^2/2$ in the plane; that is the curvature singularity of Semerák (2001), and the higher members do not have it.

## Step 5. The diagrams

Every drawing is the first disc at $m = a/5$ in units of $a$: the dust moves at $0.94\,c$ at the rim and the redshift from the centre is $0.60$.
The rotations fix the axis and the reflection through the disc fixes the plane $z = 0$, so both are totally geodesic.
The spacetime diagrams draw the axis and the plane in Weyl's chart, with $\psi$ and $\gamma$ declared as the published functions written in $\rho$ and $z$, and the axis above the disc, the plane outside the rim and the disc itself in the oblate spheroidal chart; `null_rays.py --verify` compares their rays with the tortoise coordinate of each, taken by quadrature from the closed forms in floats.
The conformal diagrams use $p, q = \arctan((ct \mp x_*)/\ell)$: the whole axis is Minkowski's diamond with the disc down its middle, the plane is the triangle with the rim a timelike singular line, and the oblate chart's three views are the triangle from the centre of the disc, the triangle from the rim, and the strip between centre and rim.
The embedding diagram is the plane $z = 0$ at $t = 0$, a cap for the disc and a sheet outside it that meet at the rim with one tangent, since $e^{2\gamma} \ge (1 - \rho\,\partial_\rho\psi)^2$ everywhere on the plane.
