# The Einstein-Rosen bridge

Einstein and Rosen's paper of 1935, "The Particle Problem in the General Theory of Relativity", has two bridges: a neutral one made from Schwarzschild's metric and a charged one made from Reissner and Nordström's with the sign of the Maxwell stress tensor reversed.
This note records the five charts, where each comes from, what the bridge is in the extension through its horizon, and what each diagram draws.
Their $2m$ is $r_s = 2GM/c^2$ here, and their $\epsilon^2/2$ is $r_q^2$, with $r_q$ the charge radius of the Reissner-Nordström page, $r_q^2 = GQ^2/(4\pi\epsilon_0c^4)$.

## Step 1: the neutral bridge

Their equation (5) is Schwarzschild's metric for $r > 2m$, and with $u^2 = r - 2m$ it becomes their (5a), $ds^2 = -4(u^2 + 2m)du^2 - (u^2 + 2m)^2d\Omega^2 + \frac{u^2}{u^2 + 2m}dt^2$ in their signature.
In the collection's signature that is the chart `bridge`: $ds^2 = -\frac{u^2}{u^2 + r_s}c^2dt^2 + 4(u^2 + r_s)du^2 + (u^2 + r_s)^2d\Omega^2$ with $u$ over the whole line.
Every component is finite at $u = 0$, where $g_{tt}$ vanishes, and the determinant is $-4u^2(u^2 + r_s)^4\sin^2\theta$, which vanishes there and has one sign on both sheets, their footnote 2.
$u$ carries the dimension of the square root of a length, which `DIMENSIONS` declares as `L**(1/2)`.

The chart `spherical` is their (5), Schwarzschild's $t$ and $r$ on one sheet, and `isotropic` is the isotropic radius, whose areal radius is $r(1 + r_s/4r)^2$: the bridge is $r = r_s/4$, the inversion $r \to r_s^2/16r$ exchanges the sheets, and $r = 0$ is the far end of the other one.
`print_charts.einstein_rosen_bridge_check` holds each of the three to being a vacuum for $u \neq 0$ and to being the published metric of Schwarzschild's spherical chart carried along its map.

## Step 2: the charged bridge

Their equation (8) is $ds^2 = -\frac{dr^2}{1 - 2m/r - \epsilon^2/2r^2} - r^2d\Omega^2 + \left(1 - \frac{2m}{r} - \frac{\epsilon^2}{2r^2}\right)dt^2$, with the potential $\varphi_4 = \epsilon/r$.
They say why the sign is what it is: "we are forced to put the negative of the above into the gravitational equations if it is to be possible to obtain static spherically symmetric solutions of the equations, free from singularities", and their footnote 4 adds that the usual sign gives $+\epsilon^2$ and no bridge.
The chart `charged_spherical` is that metric, $f = 1 - r_s/r - r_q^2/r^2 = (r - r_+)(r - r_-)/r^2$ with $r_\pm = (r_s \pm \sqrt{r_s^2 + 4r_q^2})/2$, so $r_- < 0 < r_+$ for every mass of either sign.
The script checks that it is the published Reissner-Nordström metric with $r_q^2 \to -r_q^2$, that $R = 0$, and that $G^t{}_t = G^r{}_r = r_q^2/r^4 = -G^\theta{}_\theta$: a negative energy density $-c^4r_q^2/8\pi Gr^4$.

With $m = 0$ and $u^2 = r^2 - \epsilon^2/2$ they reach (8a), $ds^2 = -du^2 - (u^2 + \epsilon^2/2)d\Omega^2 + \frac{2u^2}{2u^2 + \epsilon^2}dt^2$.
The chart `charged_bridge` is that metric: $ds^2 = -\frac{u^2}{u^2 + r_q^2}c^2dt^2 + du^2 + (u^2 + r_q^2)d\Omega^2$, with $u$ a length, the proper radial distance from the bridge.
Its space is Ellis and Bronnikov's, $dr^2 + (r^2 + \ell^2)d\Omega^2$ with $\ell = r_q$, and its $g_{tt}$ vanishes at $u = 0$.
The script checks that it is `charged_spherical` at $r_s = 0$ pulled back along $r = \sqrt{u^2 + r_q^2}$.
Einstein and Rosen say that a bridge with both a mass and a charge can be formed as well and write no coordinate for it, so no chart through that bridge is published.

## Step 3: what the bridge is

For the neutral bridge $r_* = r + r_s\ln(r/r_s - 1) \to -\infty$ on $u = 0$, so a light ray reaches the bridge only as $t \to \pm\infty$: the bridge is the horizon of Schwarzschild's spacetime.
In Kruskal and Szekeres's extension the sheet $u > 0$ is the exterior I and the sheet $u < 0$ the exterior I', and the map from $(t, u)$ sends every point of $u = 0$ at a finite $t$ to the bifurcation sphere, the one sphere the two exteriors share.
Guendelman, Kaganovich, Nissimov and Pacheva stress the difference: in the chart of $t$ and $u$ the two sheets meet along the whole hypersurface $u = 0$, and there $\sqrt{-g_{tt}} \sim |u|$ makes $R^0{}_0 \sim \delta(u)/|u|$, so the metric written in $u$ is a vacuum only away from the bridge.
The components published are those of $u \neq 0$, which is what the checker compares.
Kruskal's time runs against $t$ on I', so the conformal diagram's check takes $-\partial_t$ for the future there.

The charged bridges have one horizon, at $r_+$.
$1/f = 1 + \frac{A_+}{r - r_+} + \frac{A_-}{r - r_-}$ with $A_\pm = r_\pm^2/(r_\pm - r_\mp)$, so $r_* = r + A_+\ln|r/r_+ - 1| + A_-\ln|r/r_- - 1|$ vanishes at $r = 0$, and the extension through $r_+$ has the diagram of Schwarzschild's, with $r = 0$ a spacelike singularity on $T = \pm\pi/2$, where $K = (12r_s^2r^2 + 48r_sr_q^2r + 56r_q^4)/r^8$ diverges.
With no mass, $r_\pm = \pm r_q$, $r_* = r + \tfrac{1}{2}r_q\ln((r - r_q)/(r + r_q))$ and the surface gravity is $c^2/r_q$.
With a mass the diagrams take $r_q = \sqrt{3}\,r_s/2$, where $r_+ = 3r_s/2$, $r_- = -r_s/2$, $A_+ = 9r_s/8$ and $A_- = -r_s/8$.

## Step 4: the diagrams

The spacetime diagrams are drawn at $r_s = 1$, at $r_s = 1$ with $r_q = \sqrt{3}/2$, and at $r_q = 1$.
`null_rays.py --verify` checks the rays of the five planes against $ct \mp r_*\,\mathrm{sgn}(u)$, $ct \mp r_*$, $ct \mp r_*\,\mathrm{sgn}(4r - r_s)$, $ct \mp r_*$ and $ct \mp r_*\,\mathrm{sgn}(u)$ with the three tortoise coordinates of Step 3.
The two charts of an areal radius are drawn from below the bridge, with `where` leaving that part hatched and empty, since no sphere of the spacetime is smaller than the bridge.

The conformal diagrams set each chart in the extension through its horizon, the tower of one positive root, and tint the exteriors the chart covers.
The embedding diagrams are the equator of a moment of $t$: Flamm's paraboloid $z^2 = 4r_s(r - r_s)$ on both sheets for the neutral bridge, on which Einstein and Rosen's coordinate is the height, $z = 2\sqrt{r_s}\,u$; the catenoid $z = r_q\,\mathrm{arsinh}(u/r_q)$ for the charged bridge with no mass; and the surface $dz/dr = \sqrt{(r_sr + r_q^2)/(r^2 - r_sr - r_q^2)}$ for the charged bridge with a mass, by quadrature.
Each is static, so each is one surface and no movie.
