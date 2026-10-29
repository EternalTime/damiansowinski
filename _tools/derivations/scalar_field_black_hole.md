# A black hole in a universe filled by a scalar field

Schwarzschild's black hole in a massless scalar field whose gradient points along one axis far away.
This note records the three charts, where each comes from, the field equations each is held to, and what each diagram draws.
Units here have $G = c = 1$ where a source does; the metric file keeps $c$ and writes the mass as the length $m = GM/c^2$.

## Step 1: the sources and who wrote the metric first

Somasri Sen and Narayan Banerjee give the metric as the third set of their minimally coupled solutions, Pramana 56, 487 (2001), eqs. (2.33) and (2.34), arXiv:gr-qc/9809064v4.
Their (2.33) is Kerr-Newman in Boyer-Lindquist form with $d\theta^2 + dr^2/\Delta$ multiplied by $e^{-\frac{\sigma^2}{2}(r^2 - 2mr + a^2 + e^2)\sin^2\theta}$, and their (2.34) is the field $\phi = \sigma(r - m)\cos\theta$, with $R_{\mu\nu} = -\phi_{,\mu}\phi_{,\nu}$ in the signature $(+,-,-,-)$.
At $a = e = 0$ and $\sigma^2 = 2b^2$ it is the metric of the spherical chart here.
The third set is absent from v1, v2 and v3 of the preprint (September 1998 and March 1999) and first appears in v4 (31 May 2001), posted after the journal issue of April 2001 under its title and DOI; the Pramana text itself was not reached.
Vitor Cardoso and José Natário found it again from Buchdahl's reciprocal ansatz and Ernst's lead, General Relativity and Gravitation 57, 138 (2025), arXiv:2410.02851, and do not cite Sen and Banerjee.
Carlos Herdeiro placed it in Weyl's construction for a massless scalar, Physics Letters B 860, 139160 (2025), arXiv:2410.12950, as the growing dipole of the scalar sector on Schwarzschild's rod.
Cristian Stelea, Marina-Aura Dariescu and Vitalie Lungu point to Sen and Banerjee's form of the rotating version, arXiv:2501.11807.
The literature review behind these statements is kept with the task's records, `literature.md`.

## Step 2: the spherical chart

Cardoso and Natário's line element, with $f = 1 - 2m/r$:
$ds^2 = -f\,dt^2 + F(dr^2/f + r^2d\theta^2) + r^2\sin^2\theta\,d\phi^2$, $F = e^{-b^2r^2f\sin^2\theta} = e^{-b^2(r^2 - 2mr)\sin^2\theta}$, $\Phi = K(r - M)\cos\theta/\sqrt{4\pi}$ with their $K = b$ and $M = m$.
Their $R_{\mu\nu} = 8\pi\Phi_{,\mu}\Phi_{,\nu}$ is $R_{\mu\nu} = 2\,\partial_\mu\varphi\,\partial_\nu\varphi$ for $\varphi = b(r - m)\cos\theta$, and in the units of the metric file the field is $\sqrt{c^4/4\pi G}\,b(r - m)\cos\theta$ with $R_{\mu\nu} = (8\pi G/c^4)\,\partial_\mu\varphi\,\partial_\nu\varphi$.
Herdeiro's section 4.2 gives the same metric with his $b_1 = b$; Lim, Eur. Phys. J. C 86, 366 (2026), arXiv:2604.08951, writes it with his $b$ and $\varphi = \sqrt{2}\,b(r - m)\cos\theta$ for his normalisation of the action.
The letter $b$ is Lim's, since $K$ is the Kretschmann scalar on the page.
`print_charts.scalar_field_black_hole_check` holds the chart to $R_{\mu\nu} = 2\,\partial_\mu\varphi\,\partial_\nu\varphi$ in every slot, to the wave equation with $\sqrt{-g} = r^2F\sin\theta$, and to Schwarzschild's metric at $b = 0$.
The Ricci scalar is $R = 2b^2(r^2 - 2mr + m^2\sin^2\theta)/(r^2F)$, Herdeiro's section 4.2 and Lim's section on its causal structure; it is finite on the horizon and grows without bound as $r \to \infty$ off the axis.
The sum $r^2 - 2mr + m^2\sin^2\theta = (r - m)^2 - m^2\cos^2\theta$ is $m^2(x^2 - y^2)$, the product of the distances from the two ends of Schwarzschild's rod in Weyl's chart, and the printer writes it as one factor.
The Kretschmann scalar is printed by powers of $b$, $K = \left(48m^2 + 16b^2mr(\dots) + 4b^4r^2(r^2 - 2mr + m^2\sin^2\theta)(\dots)\right)/(r^6F^2)$, checked against the scalar built from every component; at $b = 0$ it is Schwarzschild's $48m^2/r^6$, and on the horizon $r = 2m$ it is $\frac{3}{4m^4}(1 - \frac{4}{3}\beta + \beta^2)$ with $\beta = b^2m^2\sin^2\theta$, finite.
Values are factored as the black hole in a tidal field's spherical chart is, by `distorted_spherical_factors`, with $F$ standing as a factor and every factor linear in $r$ written with $r$ first.

## Step 3: Weyl's chart

Herdeiro's construction: with $ds^2 = -e^{2\psi}dt^2 + e^{-2\psi}(e^{2\gamma}(d\rho^2 + dz^2) + \rho^2d\phi^2)$ and a field $\varphi(\rho, z)$, the equations are Laplace's for $\psi$ and for $\varphi$ in the flat space of $\rho$, $\phi$ and $z$, and the quadrature $\partial_\rho\gamma = \rho((\partial_\rho\psi)^2 - (\partial_z\psi)^2 + (\partial_\rho\varphi)^2 - (\partial_z\varphi)^2)$, $\partial_z\gamma = 2\rho(\partial_\rho\psi\,\partial_z\psi + \partial_\rho\varphi\,\partial_z\varphi)$, his section 2 with his field scaled by $\sqrt{4\pi}$ from Cardoso and Natário's.
The two sources of $\gamma$ add, so the field $\varphi = bz$ on Schwarzschild's rod adds $-b^2\rho^2/2$ to the rod's $\gamma$.
The rod's functions are written in the prolate spheroidal $x$ and $y$, as in the black hole in a tidal field: $\psi = \frac{1}{2}\ln\frac{x - 1}{x + 1}$, $\gamma = \frac{1}{2}\ln\frac{x^2 - 1}{x^2 - y^2} - \frac{b^2\rho^2}{2}$.
Herdeiro's change of coordinates in his section 3, $\rho = \sqrt{r^2 - 2mr}\sin\theta$, $z = (r - m)\cos\theta$, carries the chart onto the spherical one, and $\varphi = bz = b(r - m)\cos\theta$.
$\psi$ and $\gamma$ are held as functions (`HELD` in `verify_metrics.py`) with their slopes declared (`RATES`), $\psi$'s Schwarzschild's and $\gamma$'s the quadrature with the field, each checked there against its definition.
`print_charts.scalar_field_black_hole_weyl_reduce` writes $\gamma$'s derivatives by that quadrature and $\psi$'s second derivative along $\rho$ by Laplace's equation, so every value is printed in $\psi$, $\gamma$, $b$, $\rho$ and the derivatives of $\psi$.
On those equations the Ricci tensor is $R_{zz} = 2b^2$ alone, and the components that vanish only there are left out with the rest of the vanishing ones.
`scalar_field_black_hole_weyl_check` holds the chart to that, to the quadrature being integrable on Laplace's equation, and at six random points in forty digits to Laplace's equation, to the quadrature, to $\gamma = 0$ on the axis beyond the rod, and to being the spherical chart carried by that change of coordinates.

## Step 4: the universe with no black hole

Cardoso and Natário's scalar counterpart to Melvin's universe, $m = 0$: $ds^2 = -dt^2 + e^{-K^2\rho^2}(d\rho^2 + dz^2) + \rho^2d\phi^2$, $\Phi = Kz/\sqrt{4\pi}$, from Buchdahl's ansatz with $\beta = 0$ and $2\gamma = -K^2\rho^2$.
It is Weyl's chart at $m = 0$, where $\psi = 0$.
The chart is held to $R_{\mu\nu} = 2\,\partial_\mu\varphi\,\partial_\nu\varphi$ with $\varphi = bz$, to the wave equation with $\sqrt{-g} = \rho e^{-b^2\rho^2}$, and to flat space at $b = 0$.
Its curvature is $R = 2b^2e^{b^2\rho^2}$ and $K = 12b^4e^{2b^2\rho^2}$, least on the axis and growing without bound away from it.
The size of the field's gradient is $|\nabla\varphi| = be^{b^2\rho^2/2}$, Cardoso and Natário's $Ke^{K^2\rho^2/2}/\sqrt{4\pi}$, and a radial distance $\int_0^\infty e^{-b^2\rho^2/2}d\rho = \sqrt{\pi/2}/b$ reaches $\rho = \infty$, Herdeiro's affine length for light along the same line.

## Step 5: the diagrams

(written with the diagrams)
