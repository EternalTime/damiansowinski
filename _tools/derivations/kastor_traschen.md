# Kastor-Traschen, charged black holes that collide

Kastor and Traschen's metric (Phys. Rev. D 47, 5370, 1993) is

$$ds^2 = -\frac{c^2dt^2}{\Omega^2} + a^2\Omega^2\left(dx^2 + dy^2 + dz^2\right),\qquad \Omega = 1 + \sum_i \frac{m_i}{a|\mathbf{x} - \mathbf{x}_i|},\qquad a = e^{Ht},$$

with $m_i = GM_i/c^2$, every charge fixed by $GQ_i^2/(4\pi\epsilon_0c^4) = m_i^2$ and of one sign, and $H = \pm c\sqrt{\Lambda/3}$.
Its four charts are written by `_tools/derivations/print_charts.py --metric kastor_traschen` in about 35 seconds, and `verify_metrics.py --system kastor_traschen/<chart>` checks each in a few seconds.

## Step 1. The charts and their sources

The Cartesian chart is Brill, Horowitz, Kastor and Traschen's (Phys. Rev. D 49, 840, 1994, their 3.1), in the time $\tau$ with $H\tau = e^{Ht}$:

$$ds^2 = -\frac{c^2d\tau^2}{U^2} + U^2\left(dx^2 + dy^2 + dz^2\right),\qquad U = H\tau + V,\qquad V = \sum_i \frac{m_i}{|\mathbf{x} - \mathbf{x}_i|}.$$

It runs on past $\tau = 0$ to the singularity $U = 0$, which the comoving chart never reaches, and it is the page's first chart for that reason.
The cylindrical chart, $V = V(\rho,z)$, holds any number of holes strung along one axis, the setting of their section 5, whose axis and midplane the diagrams draw.
The isotropic chart is one hole alone, $V = m/r$, their (2.3), which is the lukewarm Reissner-Nordström-de Sitter hole with areal radius $R = rU = H\tau r + m$.
The comoving chart is Kastor and Traschen's own, with $a\Omega$ for $U$.
The static chart of one hole is on the page of the Reissner-Nordström-de Sitter black hole and is not repeated.

The free $V$ is left free in every tensor, so no component assumes Laplace's equation.
Each chart defines a name for the one sum its metric is built on, $U = H\tau + V$ or $\Omega = 1 + V/a$ with $a = e^{Ht}$, and the time enters every value through that sum alone.
`KastorTraschenForms` therefore writes the time out of each value, $c\tau = c(U - V)/H$ and $V = a(\Omega - 1)$, and factors what is left, so no $\tau$ and no $t$ stands in any printed value.
The checker holds $U$ and $\Omega$ as functions while it builds the tensors and writes their definitions out where two values are compared.

## Step 2. The field equations

With $A = c\,d\tau/U$ in units where $G = c = 4\pi\epsilon_0 = 1$ and $\Lambda = 3H^2$, `kastor_traschen_check` confirms before anything is written that

$$G_{\mu\nu} + \Lambda g_{\mu\nu} = 2\left(F_{\mu\alpha}F_\nu{}^\alpha - \tfrac14 g_{\mu\nu}F_{\alpha\beta}F^{\alpha\beta}\right)$$

once $\partial_z^2V$ is replaced by $-\partial_x^2V - \partial_y^2V$, that Maxwell's equations reduce to Laplace's for $V$, and that at $H = 0$ the chart is Majumdar and Papapetrou's Cartesian chart.
The Ricci scalar is $12H^2/c^2 - 2\nabla^2V/U^3$, which is $4\Lambda$ where $V$ is harmonic.
`kastor_traschen_pullback` checks every other chart against the Cartesian one, slot by slot, and the isotropic chart against the cosmological chart of the Reissner-Nordström-de Sitter black hole at $r_s = 2m$.

## Step 3. What is drawn

Every diagram draws holes falling together, $H < 0$, in units of the mass parameter $m$ of one hole.
One hole alone is drawn at $H = -3c/(16m)$, where $4m|H|/c = 3/4$ and the horizons are at the areal radii $0.861\,m$, $4m/3$ and $4m$, twice the lukewarm hole's at $r_s = 1$.
Two holes of mass parameter $m$ each stand at $z = \pm 2m$ with $H = -3c/(32m)$, so that they merge into the same lukewarm hole of mass parameter $2m$.

On the axis through two holes and on the plane midway between them, symmetry keeps every ray of the plane a null geodesic, with $dx/d(c\tau) = \pm 1/U^2$.
With $R = H\tau x$ and $y = \ln x$ an outgoing ray obeys $dR/dy = R + H(R + xV)^2$, their (4.6) with $xV$ for the mass, and far from the pair, where $xV \to 2m$, the right hand side vanishes at $R = 2m/3$ and $R = 6m$.
A ray above $2m/3$ runs on to the cosmological horizon, a ray below it falls to $U = 0$, and the one ray that tends to $2m/3$ divides them.
That ray is the event horizon on the plane, by their definition of the horizon as the boundary of the events that can reach the cosmological horizon: the horizon's null normal at a point of either plane lies in the plane, by the reflection and the rotation, so its generator there is a ray of the plane that neither escapes nor falls.
The slope of the right hand side at $2m/3$ is $+1/2$, so `slices.kt_last_ray` finds the ray by integrating inward from $x = 10^7 m$, where every ray closes on it.
On the midplane it leaves the axis at $c\tau = -6.1995\,m$, the event at which the horizons of the two holes join, and far out it keeps $c\tau x = -64m^2/9$.

The conformal diagram is one hole's, the lukewarm tower read from the static chart of the Reissner-Nordström-de Sitter black hole: with $\tau' = -\tau/2$ and $\rho = r/2$ the isotropic chart with $H < 0$ is that spacetime's cosmological chart run backward, so an event is drawn at the point of $(\tau', \rho)$ turned upside down.
Two holes have no conformal diagram: their axis is cut by the holes into three strips whose horizons are only finitely differentiable, and no pair of null coordinates in closed form covers them.
The embedding diagram draws one hole's equator at $H\tau = 1/2$, Majumdar and Papapetrou's throat over the areal radius, and the midplane between two holes as a movie from $c\tau = -8m$ to $-m$, a surface of revolution with circumference radius $\rho U$ that folds up into the floor of one throat of radius $2m$ as $\tau \to 0$.
