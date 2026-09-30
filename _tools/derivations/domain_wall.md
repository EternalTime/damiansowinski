# The domain wall of Vilenkin and of Ipser and Sikivie

The four charts of `domain_wall.json` are written by `print_charts.py --metric domain_wall`, which computes every tensor from each line element with the checker's own `Geometry` and reads every printed value back before it writes.
This note derives what the charts are, why their tensors come out as printed, and what the diagrams draw.

## Step 1. The planar chart

A thin wall with surface energy density $\sigma$ equal to its tension has, in Gaussian normal coordinates about it,

$$ds^2 = \left(1 - k|z|\right)^2\left(-c^2dt^2 + e^{2kct}\left(dx^2 + dy^2\right)\right) + dz^2, \qquad k = \frac{2\pi G\sigma}{c^4},$$

with $z$ the proper distance from the wall and $-1/k < z < 1/k$.
Each surface of constant $z$ is $(1 - k|z|)^2$ times de Sitter space in $2 + 1$ dimensions of radius $1/k$, in its flat slicing.
Every component is taken in the chart $x^0 = ct$.

## Step 2. The two other charts of the wall's world volume

The flat slicing covers half of de Sitter space.
Its closed slicing, $-c^2dt^2 + \cosh^2(kct)\,d\Omega^2/k^2$, covers all of it, and puts the wall at each moment on a sphere; that is the global chart,

$$ds^2 = \left(1 - k|z|\right)^2\left(-c^2dt^2 + \frac{\cosh^2(kct)}{k^2}d\Omega^2\right) + dz^2 .$$

With $1 - k|z| = e^{-k|w|}$, so that $dz = e^{-k|w|}dw$, it becomes the conformal form of Cvetič and Soleng's equation (4.1), $a(w) = -k|w|$,

$$ds^2 = e^{-2k|w|}\left(-c^2dt^2 + dw^2 + \frac{\cosh^2(kct)}{k^2}d\Omega^2\right),$$

and $w$ runs over the whole line while $z$ runs from $-1/k$ to $1/k$.

## Step 3. The inertial chart

With $\zeta = 1/k - |z|$ the global chart's side $z < 0$ is $d\zeta^2 + \zeta^2\left(-k^2c^2dt^2 + \cosh^2(kct)\,d\Omega^2\right)$, which is Minkowski space outside the light cone of one event, in its de Sitter slicing: $cT = \zeta\sinh(kct)$ and $R = \zeta\cosh(kct)$ take it to $-c^2dT^2 + dR^2 + R^2d\Omega^2$.
The wall, $\zeta = 1/k$, is the hyperboloid $R^2 - c^2T^2 = 1/k^2$, and the horizons $z = \pm 1/k$, $\zeta = 0$, are the light cone $R = c|T|$ of the centre at $T = 0$.
The inertial chart covers the whole side, the inside of the hyperboloid, the region inside that light cone included; the other side is the same region of a second copy.
Every tensor of the inertial chart vanishes.

## Step 4. The kink

The metric of the first three charts is continuous with a jump in its first derivative at the wall.
The checker's `norm` writes $|z|$ as $z\,\mathrm{sgn}(z)$, so that $\partial_z|z| = \mathrm{sgn}(z)$ and $\partial_z\mathrm{sgn}(z) = 2\delta(z)$, and reduces $\mathrm{sgn}(z)^2 = 1$.
A curvature is two derivatives of the metric, so it is a bounded function plus $\delta(z)$ times a function continuous at $z = 0$, which is its value there times $\delta(z)$, and `_on_a_kink` makes that reduction term by term.
The Christoffel symbols carry $\mathrm{sgn}(z)$, as $\Gamma^t{}_{tz} = -k\,\mathrm{sgn}(z)/(1 - k|z|)$.

## Step 5. The curvature

Off the wall every chart is flat: in the planar chart the $(1 - k|z|)^2$ de Sitter factor has curvature $k^2/(1 - k|z|)^2$, cancelled exactly by $(\partial_z(1 - k|z|))^2/(1 - k|z|)^2 = k^2/(1 - k|z|)^2$.
On it $\partial_z^2(1 - k|z|) = -2k\,\delta(z)$ gives $R_{tztz} = -2k\,\delta(z)$ and $R^x{}_{zxz} = R^y{}_{zyz} = 2k\,\delta(z)$, and the Ricci scalar $12k\,\delta(z)$.
The Einstein tensor is $G^t{}_t = G^x{}_x = G^y{}_y = -4k\,\delta(z)$ and $G^z{}_z = 0$, which with $G^\mu{}_\nu = (8\pi G/c^4)T^\mu{}_\nu$ is $T^t{}_t = T^x{}_x = T^y{}_y = -\sigma\,\delta(z)$: surface energy density $\sigma$ and tension $\sigma$, and no stress across the wall.
That is the Israel junction condition, the jump $2k$ in the extrinsic curvature $-k\,h_{ab}$ of the wall on either side.
The Weyl tensor vanishes, on the wall too: each chart is conformal to de Sitter space of $2 + 1$ dimensions times a line, which is conformally flat.
The Kretschmann scalar is $48k^2\delta(z)^2$ in sympy's hands, and a square of a delta is no distribution, so it is published as $K = 0$ for $z \neq 0$, which `scalar_parts` reads and `off_support` compares where $z$ is not zero.

## Step 6. The spacetime diagrams

On the plane $x = y = 0$ of the planar chart the metric is $-(1 - k|z|)^2c^2dt^2 + dz^2$, Rindler's on each side, and a ray runs along $kct \mp \mathrm{sgn}(z)\ln(1 - k|z|) = $ const, which `--verify` holds each ray to.
No Christoffel symbol turns a ray of that plane out of it, since $\Gamma^x{}_{ab}$ and $\Gamma^y{}_{ab}$ vanish for $a, b \in \{t, z\}$ at $x = y = 0$.
The inertial chart's plane through the centre is Minkowski's, cut at the hyperbola, whose outside the domain hatches.
The global chart's plane of $t$ and $z$ is the planar chart's, and the conformal chart's is conformal to it, so neither is drawn again.

## Step 7. The conformal diagram

Minkowski's maps $p = \arctan(k(cT - R))$ and $q = \arctan(k(cT + R))$ send the hyperbola $\tan p\tan q = -1$ to the vertical line $X = q - p = \pi/2$, from $(\pi/2, -\pi/2)$ on $\mathscr{I}^-$ to $(\pi/2, \pi/2)$ on $\mathscr{I}^+$.
The side $z < 0$ is the part $0 \le X \le \pi/2$ of Minkowski's triangle, and the side $z > 0$ its mirror image $X \to \pi - X$, which is $p = \arctan(k(cT + R)) - \pi/2$ and $q = \arctan(k(cT - R)) + \pi/2$ and meets the first along the wall at the same $T$.
A ray moving right keeps its $p$ up to the wall, and on the far side, where it moves toward the centre, the mirrored map gives it the same $p$, so it is one straight line.
The whole spacetime is the hexagon with the centres on $X = 0$ and $X = \pi$ and no $i^0$, each point a 2-sphere, since the $SO(3)$ about each centre is one symmetry of both sides.
The planar, global and conformal charts cover the diamond between the wall and the horizons, the light cones of the two centres at $T = 0$, and `conformal.py` checks each map null and future directed against the published metric.

## Step 8. The embedding diagram

A moment of the global chart has, on its equator, $g_{zz} = 1$ and $\rho = \sqrt{g_{\phi\phi}} = (1 - k|z|)\cosh(kct)/k$, so $|d\rho/dz| = \cosh(kct) \geq 1$, and away from $ct = 0$ no surface of revolution in flat space carries it.
In Minkowski space it climbs at $dZ/dz = \sqrt{\cosh^2(kct) - 1} = |\sinh(kct)|$: $Z = z\sinh(kct)$, two cones joined rim to rim at the wall, each closing at the centre of its side at $T = 0$, where the horizon meets it.
The side $z < 0$ is then the cone $cT = R\tanh(kct)$ of its own inertial chart with $Z = cT$, moved so that the wall stands at $Z = 0$, and the side $z > 0$ its mirror image in the wall, an isometry of Minkowski space.
At $ct = 0$ both cones are the flat disc of radius $1/k$, the moment the wall stops, and the equator is that disc taken twice, joined at its rim.
Beyond $|kct| = 1.1$ the cones lie so near the light cone that a chord turned $0.02$ round the axis beside the apex misses the surface by more than the check across allows, so the moments run from $kct = -1$ to $1$.
