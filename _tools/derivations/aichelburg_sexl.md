# The Aichelburg-Sexl ultraboost

The three charts of `aichelburg_sexl.json` are written by `print_charts.py --metric aichelburg_sexl`, which computes every tensor from each line element with the checker's own `Geometry` and reads every printed value back before it writes.
This note derives what the charts are, why their tensors come out as printed, and what the diagrams draw.

## Step 1. The Cartesian chart

Aichelburg and Sexl boost the Schwarzschild field to speed $V$ along $z$ and let $V \to c$ while the rest mass $M \to 0$ with $E = Mc^2/\sqrt{1 - V^2/c^2}$ held fixed.
Everything but a single null plane becomes flat, and the limit is

$$ds^2 = -c^2dt^2 + dx^2 + dy^2 + dz^2 - \frac{4GE}{c^4}\ln\left(\frac{x^2 + y^2}{\rho_0^2}\right)\delta(ct - z)\left(c\,dt - dz\right)^2,$$

with $\rho_0$ any length, since changing it adds a constant times $\delta(ct - z)(c\,dt - dz)^2$, which the coordinate change $v \to v + \text{const}\cdot\theta(ct - z)$ removes.
The profile is written with $\ln((x^2 + y^2)/\rho_0^2) = 2\ln(\rho/\rho_0)$ so that no component carries a radical.

Every component is taken in the chart $x^0 = ct$.
The delta's argument $ct - z$ is a function of two chart coordinates, so the checker reads it as sympy's `DiracDelta(c*t - z)`, and a prime on it as a derivative with respect to its argument: $\partial_{x^0}\delta(ct - z) = \delta'(ct - z)$ and $\partial_z\delta(ct - z) = -\delta'(ct - z)$ follow from the chain rule, so the combinations that cancel in the curvature cancel exactly.
The metric and its inverse are written by hand as Minkowski's plus the shock, and both are checked against sympy.

## Step 2. The null charts

With $u = ct - z$ and $v = ct + z$, both lengths, $-c^2dt^2 + dz^2 = -du\,dv$ and $c\,dt - dz = du$, so

$$ds^2 = -du\,dv + dx^2 + dy^2 + H\,du^2, \qquad H = -\frac{4GE}{c^4}\ln\left(\frac{x^2 + y^2}{\rho_0^2}\right)\delta(u),$$

a pp-wave in Brinkmann's form with the normalisation $g_{uv} = -1/2$ of Eardley and Giddings's equation (5).
With $x = \rho\cos\phi$ and $y = \rho\sin\phi$ the transverse plane is $d\rho^2 + \rho^2d\phi^2$ and $H = -(8GE/c^4)\ln(\rho/\rho_0)\,\delta(u)$.
No coordinate is a time, so no component carries a factor of $c$ from the chart.

## Step 3. The metric and its inverse

$g_{uu} = H$, $g_{uv} = -1/2$ and the transverse block is flat, so $\det g = -1/4$ in the Cartesian transverse plane.
The inverse has $g^{uv} = -2$, $g^{vv} = -4H$ and $g^{uu} = 0$.
$H$ enters the inverse in one slot, and linearly, so no product of deltas appears anywhere below.

## Step 4. The Christoffel symbols

The only nonconstant component is $g_{uu} = H(u, x, y)$, so the lowered symbols are $\Gamma_{uuu} = \tfrac{1}{2}\partial_uH$, $\Gamma_{uui} = \Gamma_{uiu} = \tfrac{1}{2}\partial_iH$ and $\Gamma_{iuu} = -\tfrac{1}{2}\partial_iH$ for $i = x, y$.
Raising with $g^{vu} = -2$ gives

$$\Gamma^v{}_{uu} = -\partial_uH, \qquad \Gamma^v{}_{ui} = -\partial_iH, \qquad \Gamma^i{}_{uu} = -\tfrac{1}{2}\partial_iH,$$

and with $\partial_x\ln(x^2 + y^2) = 2x/(x^2 + y^2)$ these are the printed $\Gamma^v{}_{uu} = (4GE/c^4)\ln((x^2 + y^2)/\rho_0^2)\,\delta'(u)$, $\Gamma^v{}_{ux} = 8GEx\,\delta(u)/c^4(x^2 + y^2)$ and $\Gamma^x{}_{uu} = 4GEx\,\delta(u)/c^4(x^2 + y^2)$.
No $\Gamma^u$ is nonzero, so $\ddot u = 0$ and $u$ is an affine parameter on every geodesic that crosses the shock.

## Step 5. The curvature

The Riemann tensor of a pp-wave is $R_{uiuj} = -\tfrac{1}{2}\partial_i\partial_jH$.
With $\partial_x^2\ln(x^2 + y^2) = 2(y^2 - x^2)/(x^2 + y^2)^2$ and $\partial_x\partial_y\ln(x^2 + y^2) = -4xy/(x^2 + y^2)^2$, every component is a multiple of $\delta(u)$ with the quadrupole $(x^2 - y^2)$ or $xy$ over $(x^2 + y^2)^2$, as printed.
The Ricci tensor is $R_{uu} = -\tfrac{1}{2}(\partial_x^2 + \partial_y^2)H$, and $\ln(x^2 + y^2)$ is harmonic off the axis, so every Ricci and Einstein component vanishes on the domain $(x, y) \neq (0, 0)$.
On the axis $(\partial_x^2 + \partial_y^2)\ln(x^2 + y^2) = 4\pi\,\delta(x)\delta(y)$, so $R_{uu} = (8\pi GE/c^4)\,\delta(u)\delta(x)\delta(y)$, the energy $E$ of the source carried along $u = 0$; that is why the domain leaves out the axis.
The Weyl tensor equals the Riemann tensor, since the Ricci tensor vanishes, and the Kretschmann scalar vanishes because every Riemann component carries two lower $u$ indices and $g^{uu} = 0$.

## Step 6. The geodesics across the shock

A null curve in the plane of $u$ and $v$ at fixed $x$ and $y$ has $-du\,dv + H\,du^2 = 0$, so either $u$ is constant, a ray moving right that never meets the shock, or $dv/du = H$, a ray moving left, whose $v$ jumps by

$$\Delta v = \int H\,du = -\frac{8GE}{c^4}\ln\left(\frac{\rho}{\rho_0}\right)$$

across $u = 0$.
The geodesic equations give the same jump: at $\dot u = 1$ and fixed $x$ and $y$, $\ddot v = -\Gamma^v{}_{uu} = -(4GE/c^4)\ln(\rho^2/\rho_0^2)\,\delta'(u)$, which integrated twice across $u = 0$ is $\Delta v = -(8GE/c^4)\ln(\rho/\rho_0)$, and $\ddot x = -\Gamma^x{}_{uu}$ gives each geodesic the velocity $\Delta\dot x^i = -(4GE/c^4)\,x^i/\rho^2$ toward the axis.
These are Eardley and Giddings's equations (10) and (11) with $\Phi = -8G\mu\ln\rho$.
The terms of $\ddot v$ in $\delta(u)\dot x$ multiply a delta by a velocity that jumps at the same place, a product Kunzinger and Steinbauer make sense of in Colombeau's algebra; the null curves of a fixed $x$ and $y$ never meet it, since their $\dot x$ is zero.

## Step 7. What the charts leave out

The coordinates in which every geodesic crosses the shock continuously, Eardley and Giddings's equations (9) to (11), put $\theta(u)$ and $u\,\theta(u)$ into the metric, and its curvature then holds products such as $\theta(u)\delta(u)$ and $\delta(u)^2$ that no reading of the delta as a distribution fixes.
The checker cannot compare such a product with anything, so that chart is not printed.

## Step 8. The spacetime diagrams

The plane of $u$ and $v$ at $x = \rho_0/n$, $y = 0$ for $n = 2$, $8$ and $32$, drawn with $z = (v - u)/2$ and $ct = (u + v)/2$, in units of $8GE/c^4 = \rho_0 = 1$.
`null_rays.py` integrates the rays numerically, so a view declares a smooth pulse in place of the delta, `smoothed` with $e^{-u^2/w^2}/(w\sqrt{\pi})$ at $w = 0.05$.
A pp-wave's profile may be any function of $u$, so the pulse is itself an exact solution, the field of a pulse of light of that length along the axis with the same energy, and every ray moving left crosses it with $\Delta v$ exactly $\ln n$, since the pulse integrates to one.
`--verify` holds each ray moving left to $v - \ln n\,(1 + \operatorname{erf}(u/w))/2$ and each ray moving right to $u$.
The jumps $\ln 2$, $\ln 8$ and $\ln 32$ step by $\ln 4$ between neighbours, which is the logarithm the drawings show.
The box runs from $ct = -2.5$ to $3.5$ so that the cone lattice, which is symmetric about the middle of the box, keeps every cone at least $0.36$ from $u = 0$, seven widths of the pulse.

## Step 9. The conformal diagram

On the plane $x = \rho_0/8$, $y = 0$ the metric is $-du\,dv$ off the shock, so $p = \arctan u$ and $q = \arctan(v - \Delta v\,\theta(u))$, with $\Delta v = \ln 8$, give each side of the shock half of Minkowski's diamond.
Across $u = 0$ a ray moving left keeps its $q$, so it is one straight line, and a line of constant $v$ breaks, its part behind the shock moved along the shock by $\Delta v$.
`conformal.py` checks both maps null against the published metric, the jump against $-\Gamma^v{}_{uu}$ integrated twice, and that $q$ is continuous along a ray moving left.
The plane is totally geodesic off the shock, where every Christoffel symbol vanishes, and not on it, where $\Gamma^x{}_{uu}$ turns the rays toward the axis.

## Step 10. The embedding diagram

A surface of constant $u$ and $v$ has the metric $dx^2 + dy^2$, flat, whatever $u$ is, so each moment is a flat disc, and the shock shows in what it does to free particles.
360 particles at rest on the circle of radius $8GE/c^4$ about the axis are run with the published $\Gamma^x{}_{uu}$ and $\Gamma^y{}_{uu}$ through a pulse of width $0.005$, a tenth of the spacetime diagram's, since the particles move while they cross the pulse and shift the focus by the order of its width.
By Step 6 the ring stays a circle of radius $\rho - (4GE/c^4)u/\rho$, $1 - u/2$ here, and closes on the axis at $u = c^4\rho^2/4GE = 2$; the script checks the circle to $10^{-9}$ and the radius and focus to twice the pulse's width.
A ring of radius $\rho$ focuses at a $u$ growing as $\rho^2$, so the shock focuses like a lens whose focal length grows with the square of the distance from its axis.
The moments are $u = -1$, $0.5$, $1$ and $1.5$, each the null hypersurface $u = u_k$, every $v$ on it carrying the same flat front, and they are marked as null lines on the spacetime and conformal diagrams.
