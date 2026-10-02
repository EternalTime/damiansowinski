# Initial data for two black holes

Misner's and Brill and Lindquist's data are one moment of a spacetime: a Riemannian space of three dimensions, signature $(+,+,+)$, with no time in it.
No spacetime around the two-hole slice is known in closed form, so the page states the slice and nothing of its evolution.
This note records what that changes for the checker, the five charts with their sources, the minimal surfaces, and what the embedding diagram draws.

## Step 1: one moment

At a moment of time symmetry the extrinsic curvature vanishes, the momentum constraint holds identically, and the Hamiltonian constraint is $R = 16\pi G\rho/c^2$, which in a vacuum is $R = 0$ (Cook, Living Rev. Relativ. 3, 5, section 2.3; Giulini, arXiv:1505.01403).
With the slice conformally flat, $ds^2 = \psi^4\,\delta$, the Ricci scalar is $R = -8\psi^{-5}\nabla^2\psi$, so the constraint is Laplace's equation.
No coordinate is declared a time in `DIMENSIONS`, as for the Eguchi-Hanson space, so the checker scales nothing by $c$ and a dot in a geodesic equation is a derivative with respect to arc length.
In three dimensions the Weyl tensor vanishes identically, and the charts print none.
The owned tags follow from the Cartesian chart alone, `three-dimensional` and `conformally flat`, the second by the Cotton tensor; `REGIONS` in `metric_tags.py` leaves out the isotropic chart, the single hole, which alone is spherically symmetric.
There is no light cone, so `null_rays.py` draws nothing, and the slice stands in `NOT_DRAWN` in `conformal.py`.

## Step 2: Brill and Lindquist's chart

$ds^2 = \psi^4(dx^2 + dy^2 + dz^2)$ with $\psi = 1 + \sum_i \alpha_i/|\mathbf{x} - \mathbf{x}_i|$ is Brill and Lindquist's solution (Phys. Rev. 131, 471, 1963), which Misner and Wheeler had noted in 1957 (Ann. Phys. 2, 525, their (253), as Brandt and Brügmann cite it, Phys. Rev. Lett. 78, 3606).
The papers of 1960 and 1963 were not to hand; the formulas are taken from Cook's review, where $\alpha_i = \mu_i/2$, and from Giulini's, where $\alpha_i = a_i$, and each is checked here.
$\psi$ is left a free function of $x$, $y$ and $z$, as Majumdar and Papapetrou's $U$ is, so no component assumes Laplace's equation, and the Ricci scalar is published as $-8\nabla^2\psi/\psi^5$.
`misner_brill_lindquist_check` holds the Ricci scalar to that form and the two-hole $\psi$ to being harmonic.
The cylindrical chart is the same with both holes on the axis, at $z = \pm a$, and is checked to be the Cartesian chart pulled back.

Each pole is the infinity of another sheet.
Near the pole $i$, $\psi \approx \alpha_i/r_i + A_i$ with $A_i = 1 + \sum_{j \ne i}\alpha_j/r_{ij}$, and in $r' = \alpha_i^2/r_i$ the metric is $(1 + A_i\alpha_i/r')^4$ times flat space, asymptotically flat with the mass $M_i = 2\alpha_i A_i$ as a length.
The shared sheet has $M = 2\sum_i\alpha_i$, so for two holes $M - M_1 - M_2 = -4\alpha_1\alpha_2/r_{12} = -M_1M_2/r_{12} + \dots$, Giulini's (binding energy) and Brill and Lindquist's interaction energy.
`_tools/test_misner_brill_lindquist.py` reads the three masses off $\psi$ numerically.

## Step 3: Misner's chart

Bispherical coordinates: $\rho = a\sin\eta/(\cosh\mu - \cos\eta)$, $z = a\sinh\mu/(\cosh\mu - \cos\eta)$, in which flat space is $a^2(d\mu^2 + d\eta^2 + \sin^2\eta\,d\phi^2)/(\cosh\mu - \cos\eta)^2$ and the poles $\mu \to \pm\infty$ are the points $z = \pm a$ of the axis.
So $\psi^4\,\delta = a^2\Psi^4(d\mu^2 + d\eta^2 + \sin^2\eta\,d\phi^2)$ with $\Psi = \psi/\sqrt{\cosh\mu - \cos\eta}$, the form Misner's data take (Phys. Rev. 118, 1110, 1960, as Price and Pullin write it, Phys. Rev. Lett. 72, 3297, their (1) and (2)).
The metric in brackets is the cylinder of $\mu$ times the unit 2-sphere, of scalar curvature 2, so $R = -(8/a^2)\Psi^{-5}(\partial_\mu^2\Psi + \partial_\eta^2\Psi + \cot\eta\,\partial_\eta\Psi - \Psi/4)$, and the constraint is that bracket set to zero.
$\Psi$ is left free; the check holds the Ricci scalar to that form, each term $(\cosh(\mu + 2n\mu_0) - \cos\eta)^{-1/2}$ of Misner's sum to solving the constraint, and Brill and Lindquist's $\Psi = (\cosh\mu - \cos\eta)^{-1/2} + (\alpha_1e^{\mu/2} + \alpha_2e^{-\mu/2})/\sqrt{2}\,a$ to solving it and to being their $\psi$ over $\sqrt{\cosh\mu - \cos\eta}$.
The distances from the two poles are $r_{1,2} = \sqrt{2}\,a\,e^{\mp\mu/2}/\sqrt{\cosh\mu - \cos\eta}$, which is where the exponentials come from.

Misner's sum is unchanged by $\mu \to \mu + 2\mu_0$ and by $\mu \to 2\mu_0 - \mu$.
The first lets the sphere $\mu = \mu_0$ be identified with $\mu = -\mu_0$, the wormhole of 1960, $S^1 \times S^2$ less the point at infinity; the second makes each of those spheres the fixed set of an isometry, so two copies of the region $|\mu| \le \mu_0$ join along both, the two sheets of 1963 (Ann. Phys. 24, 102; Cook, section 3.1.2; Giulini, arXiv:0910.2574, section 4.2).
In the flat coordinates the sum is $1 + \sum_{n \ge 1}(a/\sinh n\mu_0)(1/r_n^+ + 1/r_n^-)$ with $r_n^\pm$ the distance from $z = \pm a\coth n\mu_0$, the method of images, which the test holds numerically.

## Step 4: one hole and the charged chart

One hole is $\psi = 1 + r_s/4r$, the slice $t = 0$ of Schwarzschild's spacetime in isotropic coordinates, checked against the space part of the published isotropic chart of `einstein_rosen_bridge`.
The inversion $r \to r_s^2/16r$ is an isometry with the throat $r = r_s/4$ fixed, so $r \to 0$ is the infinity of the second sheet.

Brill and Lindquist's charged data are $ds^2 = (\chi\psi)^2\,\delta$ with $\chi$ and $\psi$ both harmonic and the electric field $E_i = \partial_i\ln(\chi/\psi)$.
With $\Omega = \sqrt{\chi\psi}$, $R = -8\Omega^{-5}\nabla^2\Omega = 2(\chi\psi)^{-2}|\nabla\ln(\chi/\psi)|^2$ for harmonic $\chi$ and $\psi$, which is $2E_iE^i$, the constraint with the field's energy density in units $G = c = 1$; the check holds the published Ricci scalar to it.
At $\chi = \psi$ the chart is Step 2's, and at $\psi = 1$ it is the space part of the published Cartesian chart of `majumdar_papapetrou` with $U = \chi$: every charge equals its mass and the interaction energy vanishes.
Beyond the pole $i$ the factor is $(1 + A_i\beta_i/r')(1 + B_i\alpha_i/r')$, the isotropic form of the Reissner-Nordström slice, with the mass $\alpha_i + \beta_i + \sum_{j \ne i}(\alpha_i\beta_j + \alpha_j\beta_i)/r_{ij}$.

## Step 5: the minimal surfaces

`two_holes.py` finds the minimal surfaces of revolution of Step 2's slice: a curve $(\rho(s), z(s))$ of the half plane with tangent angle $\theta$ is minimal when $d\theta/ds = n\cdot\nabla\ln(\psi^4\rho)$, the geodesic equation of $(\psi^4\rho)^2(d\rho^2 + dz^2)$.
A curve shot from the axis closes around one hole, a throat, or meets the midplane of two equal holes at right angles, a surface around both.
At a moment of time symmetry a minimal surface is marginally trapped, and the outermost is the apparent horizon (Gibbons, Commun. Math. Phys. 27, 87).
For two equal holes, $\alpha_1 = \alpha_2 = \alpha$ at $z = \pm a$, a surface around both exists for $a < 1.5324\,\alpha$, a coordinate separation of $1.532$ in units of $2\alpha$.
Brill and Lindquist found $1.56$, Bishop $1.53$ (Gen. Relativ. Gravit. 14, 717), and Alcubierre and others $1.532$ (Class. Quantum Grav. 17, 2159), who quote the first two.
The tests hold a single hole's throat to the sphere of radius $\alpha$, each throat to being stationary in area and below $16\pi M_i^2$, and the surface around both to the throat of the mass $4\alpha$ as $a \to 0$.

## Step 6: the embedding diagram

The slice is a moment already, so the three views are surfaces of it, each a surface of revolution.

The sphere through both holes, $\eta = \pi/2$, which is $\rho^2 + z^2 = a^2$: its metric is $a^2\Psi^4(d\mu^2 + d\phi^2)$ with $\Psi = (\cosh\mu)^{-1/2} + \sqrt{2}\,(\alpha/a)\cosh(\mu/2)$.
The circle at $\mu$ has the radius $a\Psi^2$ and $dz/d\mu = a\Psi\sqrt{\Psi^2 - 4\Psi'^2}$, real since $2\Psi' = -\tanh\mu\,(\cosh\mu)^{-1/2} + \sqrt{2}\,(\alpha/a)\sinh(\mu/2)$ is smaller than $\Psi$ in size.
Toward each pole the radius grows as fast as the distance, so the surface flattens into the plane of the sheet beyond that hole: it runs from one far sheet through one throat, across the shared sheet, and through the other throat to the other far sheet.
$\Psi''(0) = -1/2 + \alpha/(2\sqrt{2}\,a)$, so the circle $\mu = 0$ bulges while $a > \alpha/\sqrt{2}$ and is the narrowest on the surface below.
It is no fixed set of an isometry and so not totally geodesic; it is drawn because it is the one surface of revolution that passes through both throats.
The pieces join where each hole's throat cuts the sphere, and the apparent horizon of the pair cuts it for $1.4126 < a/\alpha < 1.5324$, below which the whole sphere lies inside it.
The movie runs $a$ from $\alpha/2$ to $2\alpha$.

The plane $z = 0$ between the holes, the fixed set of $z \to -z$: $\psi = 1 + 2\alpha/\sqrt{\rho^2 + a^2}$, the radius $\rho\psi^2$, $dz/d\rho = \psi\sqrt{\psi^2 - (\psi + 2\rho\psi')^2}$, flat on the axis and Flamm's paraboloid of the mass $4\alpha$ far out.
The apparent horizon cuts it in one circle while $a < 1.5324\,\alpha$.

One hole's equator in the isotropic chart, drawn as one piece through the throat: Flamm's paraboloid $z^2 = 4r_s(R - r_s)$ over $R = r(1 + r_s/4r)^2$ on both sheets, from $r = r_s^2/16r_{\max}$ to $r_{\max}$.

The charts took 10 seconds to print and the embedding diagram 130 seconds to draw on 2 October 2026.
