# Bach and Weyl's ring

Bach and Weyl's ring is the member of Weyl's class whose potential is Newton's potential of a uniform circular ring of mass $M$ and radius $a$ in Weyl's coordinates,

$$ds^2 = -e^{2\psi}c^2dt^2 + e^{-2\psi}\left(e^{2\gamma}\left(d\rho^2 + dz^2\right) + \rho^2d\phi^2\right),$$

$$\psi = -\frac{2m\,\mathrm{K}(\kappa)}{\pi\,l_2}, \qquad \gamma = -\frac{m^2}{4\pi^2a^2\rho}\left(\left(\rho + a\right)\left(\mathrm{E}(\kappa) - \mathrm{K}(\kappa)\right)^2 + \frac{\left(\rho - a\right)\left(\mathrm{E}(\kappa) - \left(1 - \kappa\right)\mathrm{K}(\kappa)\right)^2}{1 - \kappa}\right),$$

with $m = GM/c^2$, $l_{1,2} = \sqrt{(\rho \mp a)^2 + z^2}$ the least and greatest distances from the point to the ring, $\kappa = 4a\rho/l_2^2 = 1 - l_1^2/l_2^2$, and $\mathrm{K}$ and $\mathrm{E}$ the complete elliptic integrals of the first and second kinds with parameter $\kappa$.

The three charts are written by `_tools/derivations/print_charts.py --metric bach_weyl_ring`, one chart to a run with `--system`, and `verify_metrics.py --system bach_weyl_ring/weyl`, `/toroidal` and `/oblate_spheroidal` check them.

## Step 1. The sources

$\psi$ and $\gamma$ are Semerák's equations (10) and (11) (Physical Review D 94, 104021, 2016, arXiv:1611.03299), with his $M$, $\nu$, $\lambda$ and modulus $k$ written $m$, $\psi$, $\gamma$ and $\kappa = k^2$, the names the collection's other Weyl fields use.
His Section II gives the three charts: Weyl's $(\rho, z)$; the toroidal $(\zeta, \psi)$ of his equation (2), whose angle is written $\sigma$ here because $\psi$ is Weyl's function; and the oblate spheroidal $(R, \vartheta)$ of his equation (4), written $\xi = R/a$ and $\eta = \cos\vartheta$ here, the coordinates the first Morgan-Morgan disc's chart uses, in which the line element holds no trigonometric function.
In the oblate spheroidal chart $l_{1,2} = a\left(\sqrt{1 + \xi^2} \mp \sqrt{1 - \eta^2}\right)$, and Landen's transformation, $\mathrm{K}(1 - l_1^2/l_2^2) = \frac{l_1 + l_2}{l_2}\mathrm{K}(\kappa)$ with $\kappa = \left(\frac{l_2 - l_1}{l_2 + l_1}\right)^2 = \frac{1 - \eta^2}{1 + \xi^2}$, makes the parameter rational in the coordinates:

$$\psi = -\frac{2m\,\mathrm{K}(\kappa)}{\pi a\sqrt{1 + \xi^2}}, \qquad \gamma = -\frac{2m^2}{\pi^2a^2(\xi^2 + \eta^2)^2}\left((1 - \eta^2)^2\mathrm{K}^2 + (1 - \eta^2)\left((1 + \xi^2)(\mathrm{E}^2 + 2\mathrm{E}\mathrm{K} - 2\mathrm{K}^2) - 2\mathrm{E}\mathrm{K}\right) + (1 + \xi^2)(\mathrm{E} - \mathrm{K})\left((1 + \xi^2)(\mathrm{E} - \mathrm{K}) - 2\mathrm{E}\right)\right),$$

with $\mathrm{K}$ and $\mathrm{E}$ at that $\kappa$, from $\mathrm{E}(1 - l_1^2/l_2^2) = \frac{2l_2}{l_1 + l_2}\mathrm{E}(\kappa) - \frac{2l_1}{l_1 + l_2}\mathrm{K}(\kappa)$.
Written with Semerák's parameter, the chart held three radicals, and the checker left its Riemann and Weyl tensors and its Kretschmann scalar unchecked at 120 seconds each; in Landen's form the one radical is the $\sqrt{1 + \xi^2}$ under $\psi$, $\gamma$ is regular on the axis as written, and `bach_weyl_ring_check` holds both functions to Semerák's at six random points.
Bach and Weyl's paper of 1922 (Mathematische Zeitschrift 13, 134) was not read: its record is Crossref's, and its English translation of 2012 (General Relativity and Gravitation 44, 817) is behind the publisher's wall.
D'Afonseca, Letelier, and Oliveira (2005, Appendix A) and Semerák's footnote 8 record that several ring solutions in the literature carry misprints, Semerák's own earlier $\lambda$ among them, which is why Step 2 checks the formula and does not take it on trust.

## Step 2. The formula is checked, not copied

The reader reads `\mathrm{K}\left(\kappa\right)` and `\mathrm{E}\left(\kappa\right)` as `EllipticK` and `EllipticE` in `verify_metrics.py`, sympy's `elliptic_k` and `elliptic_e` under names of their own, as `EllipticF` is: their argument is the parameter, they differentiate by $d\mathrm{K}/d\kappa = (\mathrm{E} - (1 - \kappa)\mathrm{K})/(2\kappa(1 - \kappa))$ and $d\mathrm{E}/d\kappa = (\mathrm{E} - \mathrm{K})/(2\kappa)$, evaluate through mpmath to any number of digits, and carry scipy's and mpmath's numbers for `lambdify`.
The derivative of each is written in the two again, so every derivative of $\psi$ and $\gamma$ is a rational function of the coordinates, the radical $l_2$, $\mathrm{K}$ and $\mathrm{E}$, among which there is no relation.
`RATES` declares the first derivatives of $\psi$ in each chart, in Weyl's chart

$$\partial_\rho\psi = \frac{m}{\pi\rho\,l_2}\left(\mathrm{K} - \frac{\left(a^2 - \rho^2 + z^2\right)\mathrm{E}}{l_1^2}\right), \qquad \partial_z\psi = \frac{2m\,z\,\mathrm{E}}{\pi\,l_1^2\,l_2},$$

and those of $\gamma$ as Weyl's quadrature, $\partial_\rho\gamma = \rho\left((\partial_\rho\psi)^2 - (\partial_z\psi)^2\right)$ and $\partial_z\gamma = 2\rho\,\partial_\rho\psi\,\partial_z\psi$.
The reader differentiates each definition and compares it with the declared rate exactly, so Semerák's $\gamma$ is held to both quadratures in each chart, as an identity in $\rho$, $z$, $\mathrm{K}$ and $\mathrm{E}$.
`bach_weyl_ring_check` holds the declared rates of $\psi$ to Laplace's equation in the same way.

## Step 3. Two names each chart holds

Each chart declares $\psi$ and $\gamma$ among its parameters and holds them as functions of its two coordinates while the tensors are built, `HELD` in `verify_metrics.py`, as Erez and Rosen's quadrupole does.
`bach_weyl_ring_reduce` writes every derivative of $\gamma$ by the quadrature and the second derivative of $\psi$ along the first coordinate by Laplace's equation, which leaves the other derivatives of $\psi$ with no relation among them, so the Ricci tensor is exactly zero.
In Weyl's chart and in the oblate spheroidal chart the quadrature and Laplace's equation are the Morgan-Morgan discs', `morgan_morgan_quadrature` and `morgan_morgan_laplace`.
In the toroidal chart, with $D = \cosh\zeta - \cos\sigma$ and $\rho = a\sinh\zeta/D$,

$$\partial_\zeta\gamma = \frac{\sinh\zeta}{D}\left(\left(1 - \cosh\zeta\cos\sigma\right)\left((\partial_\zeta\psi)^2 - (\partial_\sigma\psi)^2\right) - 2\sinh\zeta\sin\sigma\,\partial_\zeta\psi\,\partial_\sigma\psi\right),$$

$$\partial_\sigma\gamma = \frac{\sinh\zeta}{D}\left(\sinh\zeta\sin\sigma\left((\partial_\zeta\psi)^2 - (\partial_\sigma\psi)^2\right) + 2\left(1 - \cosh\zeta\cos\sigma\right)\partial_\zeta\psi\,\partial_\sigma\psi\right),$$

which is the quadrature in any coordinates conformal to $\rho$ and $z$, $\partial_u\gamma = \rho\left(\rho_u(\psi_u^2 - \psi_v^2) + 2\rho_v\psi_u\psi_v\right)/(\rho_u^2 + \rho_v^2)$ and its companion, with $\rho_\zeta = a(1 - \cosh\zeta\cos\sigma)/D^2$, $\rho_\sigma = -a\sinh\zeta\sin\sigma/D^2$ and $\rho_\zeta^2 + \rho_\sigma^2 = a^2/D^2$.
There $\kappa = 1 - e^{-2\zeta}$ depends on $\zeta$ alone, $l_2 = a\sqrt{2e^{\zeta}/D}$, and $\rho \pm a = \pm a(e^{\pm\zeta} - \cos\sigma)/D$, which gives the chart's $\psi$ and $\gamma$.
The toroidal chart's values are factored in $\cosh\zeta$ and $\cos\sigma$: the metric holds $\zeta$ in $\cosh\zeta$ and $\sinh^2\zeta$ and $\sigma$ in $\cos\sigma$, so a value is a rational function of the two cosines, times $\sinh\zeta$ or not and times $\sin\sigma$ or not, once each derivative of $\psi$ taken an odd number of times along a coordinate is counted as odd in it, `bach_weyl_ring_toroidal_factors`.
`chart_printer.py` now knows $\zeta$ as a Greek letter, which it needs to write `\partial_\zeta`.

## Step 4. What the check holds

`bach_weyl_ring_check` holds each chart to a vanishing Ricci tensor, the declared rates to Laplace's equation and Weyl's quadrature, and $\psi$ and $\gamma$ to Semerák's (10) and (11) at Weyl's $\rho$ and $z$ of six random points in fifty digits.
The toroidal and oblate spheroidal charts are held to being Weyl's chart pulled back.
Weyl's chart is held to $\psi = -m/R$ at $a = 0$ and $\gamma \to -m^2\rho^2/(2R^4)$ as $a \to 0$, Curzon and Chazy's functions; to $\psi = -m/\sqrt{z^2 + a^2}$ and $\gamma = 0$ on the axis; and to $\psi R \to -m$ far away, so the mass is $m$.

## Step 5. The ring

Toward the ring $\kappa \to 1$, $\mathrm{K} \to \ln(4/\sqrt{1 - \kappa}) = \ln(4 l_2/l_1)$ and $\mathrm{E} \to 1$, so $\psi \to (m/\pi a)\ln(l_1/8a)$ and $e^{\psi}$ vanishes as $(l_1/8a)^{m/\pi a}$, Semerák's Section III B.
In the toroidal chart $l_1/2a \to e^{-\zeta}$, and in the plane of the ring $\gamma \to \mp(m^2/2\pi^2a^2)\,e^{\zeta}$, the upper sign outside the ring, $\sigma = 0$, and the lower inside it, $\sigma = \pi$, which is $\gamma \to m^2/(\pi^2a(a - \rho))$ inside, the exponent of Semerák's Section V A; `_tools/test_bach_weyl_ring.py` holds the closed form to it at $\zeta = 8$ and $10$.
So $\gamma \to -\infty$ on the outer side, where the proper distance to the ring and the time light takes to reach it are finite, and $\gamma \to +\infty$ on the inner side, where both are infinite.
The circle of Weyl's radius $\rho$ in the plane of the ring has radius $\rho\,e^{-\psi}$, which diverges at the ring from both sides.
An affine parameter along a ray in the plane has $d\lambda \propto e^{\gamma}d\rho$, so the ring is at infinite affine distance from inside as well.

## Step 6. The diagrams

Every drawing is at $m = a/2$ in units of $a$, on the two totally geodesic surfaces, the axis and the plane of the ring.
As Weyl's and the toroidal charts write it, $\gamma$ is $0/0$ on the axis, a bracket that vanishes as $\rho^2$ over $\rho$, so their spacetime diagrams of the axis are taken at $\rho = 10^{-4}a$ and $\zeta = 10^{-4}$, where $\gamma$ is below $10^{-9}$, and `_bach_weyl_numbers` in `null_rays.py` returns zero there; the oblate spheroidal chart's axis is $\eta = 1$ itself.
The tortoise coordinate has no closed form, so `_bach_weyl_star` takes it by quadrature from the closed forms in floats, and `null_rays.py --verify` compares the rays with it.
In the plane light from the axis reaches $\rho = 0.99\,a$ at $ct = 3.67\,a$ and $0.999\,a$ at $3 \times 10^7\,a$, and light from $\rho = 2a$ outside reaches the ring at $ct = 2.00\,a$.
The circle about the axis is narrowest outside the ring at $\rho = 1.220\,a$, of radius $2.075\,a$, and the circular photon orbit, $2\rho\,\partial_\rho\psi = 1$, is at $\rho = 1.535\,a$.
The slice $t = 0$, $z = 0$ has $g_{\rho\rho} \ge (\partial_\rho\sqrt{g_{\phi\phi}})^2$ outside $\rho = 1.116\,a$, where it is drawn as a surface in flat space, and the reverse on the disc inside $\rho = 0.993\,a$, where it is drawn in Minkowski space.
