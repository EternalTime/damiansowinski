# Schwarzschild's black hole in a swirling universe

Astorino, Martelli and Viganò's black hole and the swirling universe it sits in are, in Schwarzschild's coordinates and in the cylindrical chart,

$$ds^2 = F\left(-\left(1 - \frac{2m}{r}\right)c^2dt^2 + \frac{dr^2}{1 - 2m/r} + r^2\,d\theta^2\right) + \frac{r^2\sin^2\theta}{F}\left(d\phi + 4j\left(r - 2m\right)\cos\theta\,c\,dt\right)^2,\qquad F = 1 + j^2r^4\sin^4\theta,$$

$$ds^2 = F\left(-c^2dt^2 + d\rho^2 + dz^2\right) + \frac{\rho^2}{F}\left(d\phi + 4jz\,c\,dt\right)^2,\qquad F = 1 + j^2\rho^4.$$

Both charts are written by `_tools/derivations/print_charts.py --metric swirling_universe_black_hole`, and `verify_metrics.py --system swirling_universe_black_hole/spherical` and `--system swirling_universe_black_hole/cylindrical` check each in seconds.
The source is M. Astorino, R. Martelli and A. Viganò, Phys. Rev. D 106, 064014 (2022), arXiv:2205.13548, read in its TeX source (v4).

## Step 1. The metric and its constant

The source's metric carries $\left[4\jmath(r - 2m)\cos\theta + \omega_0\right]dt$ in the bracket with $d\varphi$, $\omega_0$ an integration constant "related to the frame of reference".
It is the shift $\varphi \to \varphi - \omega_0 t$, so we set $\omega_0 = 0$: the dragging rate $-g_{t\phi}/g_{\phi\phi} = -4j(r - 2m)\cos\theta$ then vanishes on the horizon at every $\theta$, and the azimuth is at rest there.
The source's $\jmath$ is written $j$, since the reader has no `\jmath`.
With no mass the metric is the source's background, written in $\rho = r\sin\theta$, $z = r\cos\theta$; it is the metric of Gibbons, Mujtaba and Pope's appendix C (arXiv:1301.3927) with their $c = -j$.

## Step 2. Units

$F$ is a pure number, so $j$ is an inverse area, $[j] = L^{-2}$, and $4j(r - 2m)\cos\theta\,c\,dt$ and $4jz\,c\,dt$ are pure numbers beside $d\phi$.
The mass enters as the length $m = GM/c^2$.
`DIMENSIONS` declares $[m] = L$ and $[j] = L^{-2}$, and the dimensional pass balances every published term.

## Step 3. The field equations

`swirling_check` in `print_charts.py` holds both charts Ricci flat in every slot before anything is written, holds the spherical chart at $j = 0$ to the published spherical chart of `schwarzschild` with $r_s = 2m$, and holds the cylindrical chart, pulled back along $\rho = r\sin\theta$, $z = r\cos\theta$, to the spherical one at $m = 0$.
It holds the spherical chart's Kretschmann scalar to the closed form of the source's appendix A, with $x = \cos\theta$, and the cylindrical chart's to Moreira, Herdeiro and Crispino's $192j^2(j^6\rho^{12} - 15j^4\rho^8 + 15j^2\rho^4 - 1)/(1 + j^2\rho^4)^6$ (arXiv:2401.05658), and the cylindrical $g_{tt}$ to Gibbons, Mujtaba and Pope's $-[(1 + c^2\rho^4)^2 - 16c^2z^2\rho^2]/(1 + c^2\rho^4)$.
$g_{tt}$ and $g^{\phi\phi}$ are printed in that grouping, $F^2$ less the square of the dragging, in both charts.

## Step 4. The hole

$g^{rr} = (r - 2m)/(rF)$, so the horizon is $r = 2m$ at every $\theta$, and $g_{t\phi} \propto r - 2m$ vanishes there.
Its area is $\int 2m\sqrt{F}\cdot 2m\sin\theta/\sqrt{F}\,d\theta\,d\phi = 16\pi m^2$ for every $j$, and its surface gravity, $\tfrac{1}{2}\partial_r(-g_{tt})/\sqrt{-g_{tt}g_{rr}}$ on the horizon, is Schwarzschild's $c^2/4m$, as the source's $T = 1/8\pi m$ says.
The ratio $g_{\phi\phi}/(g_{\theta\theta}\sin^2\theta)$ tends to $1$ at both poles, so the axis is regular with $\phi$ of period $2\pi$, and $g_{\phi\phi} > 0$ everywhere, so no circle of $\phi$ is timelike.
On the horizon's equator the circumference is $4\pi m/\sqrt{1 + 16j^2m^4}$, shorter than Schwarzschild's, while the meridian $\int_0^\pi 2m\sqrt{F}\,d\theta$ is longer, so the horizon is prolate, as Moreira, Herdeiro and Crispino find; the source's "more oblate" is not borne out.
Its Gaussian curvature on the equator is $(1 - 16j^2m^4)/(4m^2(1 + 16j^2m^4)^2)$, negative for $jm^2 > 1/4$.
`_tools/test_swirling_universe_black_hole.py` holds the published metric to all of this.

## Step 5. The planes

Reflection $\theta \to \pi - \theta$ with $\phi \to -\phi$ is an isometry, since it flips the sign of $\cos\theta$ and of $d\phi$ together, so its fixed plane $\theta = \pi/2$, $\phi = 0$ is totally geodesic.
On it $\cos\theta = 0$, the dragging term drops out, and the metric of $t$ and $r$ is $(1 + j^2r^4)$ times Schwarzschild's, so its null curves are Schwarzschild's, $ct = \pm r_* + $ const with $r_* = r + 2m\ln|r/2m - 1|$.
On the axis $F = 1$ and the plane of $t$ and $r$ is Schwarzschild's exactly.
In the universe, $z \to -z$ with $\phi \to -\phi$ fixes the plane $z = 0$, $\phi = 0$, on which the metric is $(1 + j^2\rho^4)(-c^2dt^2 + d\rho^2)$, conformally flat.
The circles of the equator, of radius $r/\sqrt{F}$, are widest where $1 - j^2r^4 = 0$, at $r = 1/\sqrt{j}$ with radius $1/\sqrt{2j}$, and the same holds in the universe at $\rho = 1/\sqrt{j}$; the distance out grows as $jr^3/3$, so the moment narrows into an unbounded spike as Melvin's does.

## Step 6. The ergoregions

$g_{tt} = -(r - 2m)[F^2 - 16j^2r^3(r - 2m)\sin^2\theta\cos^2\theta]/(rF)$ is negative on the axis, where it is $-(1 - 2m/r)$, on the equator, where it is $-F(1 - 2m/r)$, and just outside the horizon, where the bracket tends to $F^2$.
It is positive where $16j^2r^3(r - 2m)\sin^2\theta\cos^2\theta > F^2$, one region above the equator and its mirror image below, which at $jm^2 = 1/16$ come no nearer the centre than $r = 3.67\,m$ and widen to infinity beside the axis; at a fixed distance from the axis $g_{tt}$ grows as the square of the height.
In the universe the region is $|4jz\rho| > 1 + j^2\rho^4$, Gibbons, Mujtaba and Pope's, lowest at $\rho = (3j^2)^{-1/4}$, $z = 3^{-3/4}/\sqrt{j}$.
`swirl_ergo` in `projections.py` checks each of these before it draws.

## Step 7. The diagrams

The spacetime diagrams draw the hole's equator at $m = 1$, $j = 1/16$, where the widest circle is at $r = 4m$, and the universe's plane $z = 0$ at $j = 1$, each with the widest circle marked; `null_rays.py --verify` compares their rays with $ct \pm r_*$ and $ct \pm \rho$.
Each chart also has a figure of its ergoregions, the meridional plane seen from the side.
The conformal diagrams are Minkowski's half diamond for the universe, $p, q = \arctan((ct \mp \rho)\sqrt{j})$, and Kruskal and Szekeres's hexagon for the hole's equator, built by `Tower` on $f = 1 - 2m/r$ and checked against the published metric, whose conformal factor changes no null direction.
The embedding diagrams are the universe's plane $z = 0$ out to $\rho = 2/\sqrt{j}$ and the hole's equator through the bifurcation sphere out to $r = 6m$ on both sheets, each circle checked against $\rho/\sqrt{F}$ and $r/\sqrt{F}$.
The curvature figure draws the universe's Kretschmann scalar against $\rho$, $-192j^2$ on the axis with zeros at $\sqrt{j}\,\rho = \sqrt{2 - \sqrt{3}}$, $1$ and $\sqrt{2 + \sqrt{3}}$, and the hole's over its meridian plane, $48m^2/r^6 - 192j^2$ on the axis, on Kerr's signed logarithmic scale.
