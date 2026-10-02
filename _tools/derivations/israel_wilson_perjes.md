# Israel-Wilson-Perjés, Majumdar and Papapetrou's charges set spinning

Israel, Wilson and Perjés's metric, as Hartle and Hawking write it (Commun. Math. Phys. 26, 87, 1972, their 4.1), is

$$ds^2 = -|U|^{-2}\left(c\,dt + \boldsymbol\omega\cdot d\mathbf{x}\right)^2 + |U|^2\,d\mathbf{x}\cdot d\mathbf{x},\qquad \nabla^2U = 0,\qquad \nabla\times\boldsymbol\omega = i\left(U\nabla\bar U - \bar U\nabla U\right),$$

with $U$ complex and both operators those of flat space.
The sources are Perjés (Phys. Rev. Lett. 27, 1668, 1971) and Israel and Wilson (J. Math. Phys. 13, 865, 1972).
Its three charts are written by `_tools/derivations/print_charts.py --metric israel_wilson_perjes` in about 55 seconds, and `verify_metrics.py --system israel_wilson_perjes/<chart>` checks each in a few seconds.

## Step 1. The charts and their sources

The cylindrical chart is Hartle and Hawking's (4.1) for an axisymmetric $U$, with $\boldsymbol\omega$ given its one component along $\phi$ as in their (4.22):

$$ds^2 = -\frac{\left(c\,dt + \omega\,d\phi\right)^2}{W^2} + W^2\left(d\rho^2 + \rho^2d\phi^2 + dz^2\right),\qquad W = |U|.$$

$W$ and $\omega$ are left free in every tensor, so no component assumes a field equation.
The general Cartesian chart, with $W$ and three components of $\boldsymbol\omega$ free in three coordinates, is not published: with four free functions of three coordinates its curvature had not printed after eight minutes on 2 October 2026, where the cylindrical chart takes twenty seconds.

The oblate spheroidal chart is one source at an imaginary place, the case Israel and Wilson worked out, which Hartle and Hawking report as "the charged Kerr metric with equal charge and mass".
With $\rho = \sqrt{r^2 + a^2}\sin\theta$ and $z = r\cos\theta$, flat space is $\Sigma_0\left(dr^2/(r^2 + a^2) + d\theta^2\right) + (r^2 + a^2)\sin^2\theta\,d\phi^2$ with $\Sigma_0 = r^2 + a^2\cos^2\theta$, and

$$U = 1 + \frac{m}{r + ia\cos\theta},\qquad |U|^2 = \frac{\Sigma}{\Sigma_0},\qquad \omega = \frac{am(2r + m)\sin^2\theta}{\Sigma_0},\qquad \Sigma = (r + m)^2 + a^2\cos^2\theta.$$

It is the published Boyer-Lindquist chart of `kerr_newman` at $GM/c^2 = r_Q = m$ in the radius $r + m$, where $\Delta = r^2 + a^2$ has no zero for $a \ne 0$.
The curvature singularity is $\Sigma = 0$, the ring $r = -m$ on the equator; $\Sigma_0 = 0$, the ring $r = 0$ on the equator where $U$ is infinite, is a regular place of the spacetime where $\partial_t$ is null.

The spherical chart is one source of complex mass, Hartle and Hawking's (4.9) to (4.12), $U = 1 + (m + il)/(r - m)$ in their radius $R = r_{\rm flat} + M$:

$$ds^2 = -\frac{(r - m)^2}{r^2 + l^2}\left(c\,dt + 2l\cos\theta\,d\phi\right)^2 + \frac{r^2 + l^2}{(r - m)^2}dr^2 + (r^2 + l^2)\left(d\theta^2 + \sin^2\theta\,d\phi^2\right).$$

It is Brill's charged NUT space (Phys. Rev. 133, B845, 1964) with the square of the charge radius equal to $m^2 + l^2$, in the notation of the published `taub_nut`, and at $l = 0$ it is one hole of `majumdar_papapetrou`.

## Step 2. The field equations

`israel_wilson_perjes_check` holds the cylindrical chart to the Einstein-Maxwell equations before anything is written, in units where $G = c = 4\pi\epsilon_0 = 1$.
With $U = P + iQ$ the chart's equation for $\omega$ is

$$\partial_z\omega = -2\rho\left(P\,\partial_\rho Q - Q\,\partial_\rho P\right),\qquad \partial_\rho\omega = 2\rho\left(P\,\partial_zQ - Q\,\partial_zP\right),$$

whose integrability condition is $P\nabla^2Q - Q\nabla^2P = 0$.
The Maxwell field is Hartle and Hawking's (4.4) to (4.6): $\Phi + i\chi = 1/U$, $F_{ti} = \partial_i\Phi$ and $F^{ij} = \epsilon^{ijk}\partial_k\chi/\sqrt{-g}$, with $\sqrt{-g} = W^2\rho$.
On a solution, which means $W = \sqrt{P^2 + Q^2}$, the two equations for $\omega$, and Laplace's equation for $P$ and $Q$, the check finds $G_{\mu\nu} = 2\left(F_{\mu\alpha}F_\nu{}^\alpha - \tfrac14 g_{\mu\nu}F_{\alpha\beta}F^{\alpha\beta}\right)$ in every slot, $\nabla_\mu F^{\mu\nu} = 0$ and $dF = 0$.
The sign of $\omega$ is a convention, undone by $\phi \to -\phi$ or $U \to \bar U$; this one makes the spherical chart's cross term the $+2l\cos\theta$ of `taub_nut` and the spheroidal chart's the $g_{t\phi}$ of `kerr_newman`.

Each explicit chart is checked to be the cylindrical one pulled back at its own $U$, with $U$ harmonic and $\omega$ solving its equation in the flat metric of its own coordinates, and then against the published metric it reduces to.

## Step 3. Printing

The cylindrical chart is printed as Majumdar and Papapetrou's is, each numerator collected by powers of $W$.
The spheroidal chart defines $\Sigma_0$ and $\Sigma$, and `IwpForms` writes every value factored, the even powers of the sine in the cosine, each factor that is one of the two by its name, and any other sum in the shortest of three forms, one of them in powers of $r + m$.
Its Kretschmann scalar is Kerr and Newman's at $GM/c^2 = r_Q = m$ and $r \to r + m$, and $g_{\phi\phi}$ and $g^{tt}$ are written around $((r + m)^2 + a^2)^2 - a^2(r^2 + a^2)\sin^2\theta$, as Boyer and Lindquist's are.

## Step 4. What is drawn

The spinning source is drawn at $a = m$: its axis, where $dr/d(ct) = \pm(r^2 + a^2)/((r + m)^2 + a^2)$ and

$$r_* = r + m\ln\frac{r^2 + a^2}{a^2} + \frac{m^2}{a}\arctan\frac{r}{a},$$

and its principal null rays on the equator, whose projection on $t$ and $r$ keeps the same $ct \mp r_*$.
On the equator $g_{\phi\phi}$ vanishes where $R^4 + a^2R^2 + 2a^2mR = a^2m^2$ with $R = r + m$, at $r = -0.595\,m$, and inside it the circles of $\phi$ are timelike, so the equatorial view begins at $r = -m/2$.

The source of complex mass is drawn at $l = m/2$ on its plane of $t$ and $r$, with $r_* = r + 2m\ln(r - m) - (m^2 + l^2)/(r - m)$.

The cylindrical chart is drawn for Hartle and Hawking's two sources, their (4.21), at $z = \pm d$ with $d = 2m$ and $l = m/2$:

$$U = 1 + \frac{m - il}{r_1} + \frac{m + il}{r_2},\qquad \omega = \frac{2l(z + d)}{r_2} - \frac{2l(z - d)}{r_1} - \frac{2ml}{d}\left(\frac{\rho^2 + z^2 - d^2}{r_1r_2} - 1\right),$$

with $r_1$ and $r_2$ the distances from $z = d$ and $z = -d$.
The constant makes $\omega$ vanish on the axis beyond the sources; between them it is $4l(1 + m/d)$, their (4.24) with $a = 2d$, and far away $\omega \to 4ld\sin^2\theta/r$, the field of an angular momentum $2ld$, their $Na$.
On the midplane $W = 1 + 2m/s$ and $\omega = 4ld/s + 4mld/s^2$ with $s = \sqrt{\rho^2 + d^2}$, the reflection $z \to -z$ keeps the null geodesics of no angular momentum in the plane, with $d\rho/d(ct) = \pm\rho/\sqrt{W^4\rho^2 - \omega^2}$, and the circles of $\phi$ are timelike inside $\rho = 0.734\,m$.
The axis itself is not drawn in this chart, whose inverse metric is not finite on it.

The embedding diagram draws three moments $t = 0$, one for each chart: the infinite throat of the source of complex mass, the equator of the spinning source from $r = -0.21\,m$ out, and the midplane of the two sources from $\rho = 1.56\,m$ out, each of the last two beginning where $g_{xx} = (d\sqrt{g_{\phi\phi}}/dx)^2$.
The conformal diagram is the spinning source's axis, the full diamond by $p, q = \arctan((ct \mp r_*)/2m)$; the source of complex mass has none, as Taub-NUT has none, and neither have the two sources.
