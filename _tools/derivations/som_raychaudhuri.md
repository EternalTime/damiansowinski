# Som and Raychaudhuri's universe: the charts, the source and what the diagrams draw

Som and Raychaudhuri's homogeneous universe has one parameter, the angular velocity $\Omega$ of its dust, and one length, $r_c = c/\Omega$.
This note records where each chart of `som_raychaudhuri.json` comes from, what its source is, and what the diagrams take from it.
Units here have $c = 1$ unless a line keeps it; the file keeps $c$.

## Step 1. The cylindrical chart

Som and Raychaudhuri (1968) solved the Einstein-Maxwell equations for charged dust in rigid rotation with cylindrical symmetry.
The homogeneous member of their family is, in the signature $(-,+,+,+)$,

$$ds^2 = -\left(dt + \Omega r^2\,d\phi\right)^2 + dr^2 + r^2d\phi^2 + dz^2,$$

with $\phi$ periodic in $2\pi$.
That is the form every later source writes: Reboucas and Tiomno's linear class of 1983, $H = \Omega r^2$ and $D = r$, which is their family at $m = 0$; Paiva and Teixeira's equation (1), which keeps $c$ as $[c\,dt - (\Omega r^2/c)\,d\phi]^2$ with the opposite sense of $\phi$; Boyda, Ganguli, Horava and Varadarajan's (2.12) in the limit $m \to 0$; Drukker, Fiol and Simon's (2.8); and Russo and Tseytlin's (2.26) and Horowitz and Tseytlin's (4.11), with their $f/2$ and $H/2$ for $\Omega$.
The file writes it with $c$, $ds^2 = -(c\,dt + \Omega r^2d\phi/c)^2 + \dots$, so that $\Omega$ is an angular velocity: with $u = \partial_t$ the vorticity of the dust is $\tfrac{1}{2}|du|$, which is $\Omega/c$ per unit length.
Its $g_{\phi\phi} = r^2(1 - \Omega^2r^2/c^2)$ changes sign at $r_c = c/\Omega$: the circle of constant $t$, $r$ and $z$ is null there and timelike beyond, Drukker, Fiol and Simon's (2.5) and (2.6) at $l = 0$.
`print_charts.som_raychaudhuri("cylindrical")` builds it.

## Step 2. The source

In the frame of the dust, $e^0 = dt + \Omega r^2d\phi$, $e^1 = dr$, $e^2 = r\,d\phi$, $e^3 = dz$, the Einstein tensor is

$$G^a{}_b = \Omega^2\,\mathrm{diag}(-3, 1, 1, -1),$$

which is Soleng's frame for a thick spinning string at $M = \Omega r^2$ and $\rho = r$, where his $\Omega = M'/2\rho$ is constant and his heat flow $\Omega'$ vanishes.
Russo and Tseytlin's (2.28) gives the same curvature in that frame, $R_{00} = R_{11} = R_{22} = 2\Omega^2$.
A magnetic field along $z$, $F = B\,e^1\wedge e^2 = Br\,dr\wedge d\phi$, has $8\pi T^a{}_b = B^2\,\mathrm{diag}(-1, 1, 1, -1)$ in Gaussian units with $G = c = 1$.
At $B = \Omega$ it accounts for every stress and leaves $\mathrm{diag}(-2\Omega^2, 0, 0, 0)$: dust at rest in the chart, of density $\rho = \Omega^2/4\pi$, which is $\Omega^2/4\pi G$ with the constants written out, and $B = \Omega c/\sqrt{G}$.
The field is closed, and $\nabla_\nu F^{\mu\nu}$ is $2\Omega B$ along $\partial_t$ and nothing else, a charge at rest with the dust of density $\Omega B/2\pi$, and $F^\mu{}_\nu J^\nu = 0$: no Lorentz force acts on the dust, which is Som and Raychaudhuri's statement.
`som_raychaudhuri_source` checks the frame orthonormal, the Einstein tensor in it, the field's stress, the remainder, the current and the force, before the chart is written.
The Ricci scalar is $2\Omega^2/c^2$ and the Kretschmann scalar $44\,\Omega^4/c^4$, from the frame components $R_{0101} = R_{0202} = \Omega^2$ and $R_{1212} = 3\Omega^2$.

## Step 3. The Cartesian chart

With $x = r\cos\phi$ and $y = r\sin\phi$, $r^2d\phi = x\,dy - y\,dx$, and the line element is

$$ds^2 = -\left(dt + \Omega\left(x\,dy - y\,dx\right)\right)^2 + dx^2 + dy^2 + dz^2,$$

Drukker, Fiol and Simon's (2.13) and Das and Gegenberg's (32), each with the other sense of rotation.
`som_raychaudhuri_pullback` pulls the cylindrical chart back through $r = \sqrt{x^2 + y^2}$, $\phi = \mathrm{atan2}(y, x)$ and checks every slot.
In this chart the geodesic equations are $\ddot x = -2\Omega\,\dot t\,\dot y + \dots$ and $\ddot y = 2\Omega\,\dot t\,\dot x + \dots$, the Lorentz force of a uniform magnetic field on a charge, which is why every geodesic projects onto the plane of $x$ and $y$ as a circle, of radius $c/2\Omega$ for light moving across $z$.

## Step 4. The spacetime diagrams

On a cylinder of $t$ and $\phi$ at radius $r$ the metric is $-(dt + \Omega r^2d\phi)^2 + r^2d\phi^2$, van Stockum's block with $R = r_c$, so its null curves are straight: $dt = r(1 - r/r_c)\,d\phi$ moving to $+\phi$ and $dt = -r(1 + r/r_c)\,d\phi$ moving to $-\phi$.
The radial acceleration of light launched along one is $-\Gamma^r{}_{ab}k^ak^b$ with $\Gamma^r{}_{t\phi} = \Omega r$ and $\Gamma^r{}_{\phi\phi} = r(2\Omega^2r^2 - 1)$: for the curve moving to $+\phi$ it is $-r(2r/r_c - 1)$ per $d\phi^2$, zero at $r = r_c/2$, the Larmor circle of light about the axis, and toward the axis beyond it; for the curve moving to $-\phi$ it is $r(2r/r_c + 1)$, away from the axis.
`null_rays.py` draws the cylinders at $r_c/2$ and $3r_c/2$, `CYLINDERS` holds the slopes and `TURNING` what the captions say of each family.
On the Cartesian plane $y = z = 0$ the metric is $-dt^2 + dx^2$, since $g_{tx} = \Omega y$ vanishes there and $g_{xx} = 1 - \Omega^2y^2$ is 1, so the null curves run at 45°, and $\Gamma^y{}_{tx} = -\Omega$ turns light launched along one into $y$.
The cones are oriented by $\partial_t$, which is timelike everywhere, since $g^{tt} = \Omega^2r^2 - 1$ is positive beyond $r_c$ and $t$ is a time function only inside it.
The figure in three dimensions is van Stockum's and Godel's, `about_axis` in `projections.py`, with $r$ itself as the radius, since $g_{rr} = -g_{tt} = 1$.

## Step 5. The embedding diagram

The plane $z = 0$ at one moment of $t$ has $g_{rr} = 1$ and $g_{\phi\phi} = r^2(1 - r^2)$ in units of $r_c$, so its circles have radius $\rho = r\sqrt{1 - r^2}$, widest at $r = 1/\sqrt{2}$, where $\rho = 1/2$.
With $d\rho/dr = (1 - 2r^2)/\sqrt{1 - r^2}$,

$$\left(\frac{dz}{dr}\right)^2 = 1 - \frac{(1 - 2r^2)^2}{1 - r^2} = \frac{r^2(3 - 4r^2)}{1 - r^2},$$

which vanishes at $r = \sqrt{3}/2$ and is negative from there to $r_c$: the surface stops at $\sqrt{3}\,r_c/2$, where it runs level and its circles shrink at one unit of radius per unit of distance.
With $s = 1 - r^2$ the height is $\tfrac{1}{2}\int_s^1\sqrt{(4s' - 1)/s'}\,ds'$, and $\int\sqrt{(4s - 1)/s}\,ds = \sqrt{s(4s - 1)} - \tfrac{1}{4}\mathrm{arcosh}(8s - 1)$, so

$$z = \frac{\sqrt{3} - \sqrt{(1 - r^2)(3 - 4r^2)}}{2} - \frac{\mathrm{arcosh}\,7 - \mathrm{arcosh}(7 - 8r^2)}{8},$$

which the script writes with $\mathrm{arcosh}(8s - 1) = 2\,\mathrm{artanh}\sqrt{(4s - 1)/4s}$.
`embedding.som_raychaudhuri` draws the surface by quadrature and checks it against that closed form, the radius against $r\sqrt{1 - r^2}$, and the stop exactly.
Boyda and his collaborators' holographic screen for this universe, their (2.13) at $m \to 0$ with their $\Omega$ equal to $\Omega/\sqrt{2}$ here, is at $r = r_c/\sqrt{2}$, the widest circle of this surface.

## Step 6. No conformal diagram

Every event's future is the whole spacetime, as for Godel and van Stockum, so `conformal.py` names the spacetime in `NOT_DRAWN`.
