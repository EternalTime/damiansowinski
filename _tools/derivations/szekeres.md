# The Szekeres cosmologies

The two charts of `szekeres.json` are written by `print_charts.py --metric szekeres`, which computes every tensor from each line element with the checker's own `Geometry` and reads every printed value back before it writes.
This note derives what the charts are, why their tensors come out as printed, and what the diagrams draw.

## Step 1. The stereographic chart

Hellaby and Krasiński's equation (1), with $c$ restored, is

$$ds^2 = -c^2dt^2 + \frac{\left(\partial_r R - R\,\partial_r E/E\right)^2}{\epsilon + f}\,dr^2 + \frac{R^2}{E^2}\left(dp^2 + dq^2\right),$$

with $R(t,r)$ a length, $f(r)$ a pure number, and their equation (13),

$$E = \frac{S}{2}\left(\frac{(p - P)^2 + (q - Q)^2}{S^2} + \epsilon\right),$$

in three functions $S$, $P$ and $Q$ of $r$ alone and the sign $\epsilon$.
That form of $E$ solves Szekeres's constraint $4(AC - B_1^2 - B_2^2) = \epsilon$ on $E = A(p^2 + q^2) + 2B_1p + 2B_2q + C$ identically, so the chart has no relation among its parameters.
$R$, $f$, $S$, $P$ and $Q$ are left free in every tensor, as Tolman-Bondi leaves its $R$ and its energy function, so no component assumes a field equation.

## Step 2. E held as a function

$E$ is a name the chart defines, and the checker's `Reader` reads a defined name as the expression it names.
Built with $E$ written out in $p$, $q$, $S$, $P$ and $Q$, the Christoffel symbols and the curvature swell at every differentiation, and the geometry did not finish in ten minutes.
So a defined name whose definition holds a function of the coordinates is held by the `Reader` as a function of the coordinates it varies with, here $E(r,p,q)$, while the tensors are built, which takes two seconds, and `Reader.surface` writes the definition out on both sides of each comparison.
A published value therefore writes $\partial_r E$, $\partial_p E$ and $\partial_r\partial_p E$, each read as the derivative of the definition.

$E$ is quadratic in $p$ and $q$, so

$$\partial_p^2E = \partial_q^2E = \frac{1}{S}, \qquad \partial_p\partial_qE = 0, \qquad 2E = S\left((\partial_pE)^2 + (\partial_qE)^2 + \epsilon\right),$$

and every third derivative along $p$ and $q$ vanishes.
`SzekeresForms.reduce` in `print_charts.py` writes those into each value and removes $E$ and its derivatives along $r$ alone by the last relation and its derivatives along $r$.
What is left is written in $\partial_pE = (p - P)/S$, $\partial_qE = (q - Q)/S$ and their derivatives along $r$, which are independent, since $p$, $q$, $P'$, $Q'$ and the higher derivatives of $P$ and $Q$ are.
So a value that vanishes for Szekeres's $E$ is exactly zero there, and two values equal for it are equal.
`SzekeresForms.pretty` factors a reduced value and writes each factor back in $E$ and its derivatives along $r$ by the same relations, eliminating $\partial_qE$ or $\partial_pE$ first and taking the shortest of the four results.
Every printed value is read back and reduced before it is written, and `verify_metrics.py` checks the file with $E$ written out.

With those relations the Einstein tensor is diagonal, the twenty four Riemann and twenty four Weyl components are Tolman-Bondi's in number, and $S$, $P$ and $Q$ enter every tensor through $E$ alone.

## Step 3. Dust

`szekeres_check` substitutes $(\partial_tR)^2 = 2M/R + f$, with $M(r)$ a length, and its derivatives, and confirms before anything is written that the only component of $G^\mu{}_\nu$ left is

$$G^t{}_t = -\frac{2\left(M' - 3M\,\partial_rE/E\right)}{R^2\left(\partial_rR - R\,\partial_rE/E\right)} = -\frac{8\pi G\rho}{c^2},$$

Hellaby and Krasiński's equation (20), so the matter is dust at rest in the chart.

## Step 4. The axisymmetric chart

For $\epsilon = 1$ the Riemann projection $p - P = S\cot(\theta/2)\cos\phi$, $q - Q = S\cot(\theta/2)\sin\phi$, their equation (16), gives $E = S/(1 - \cos\theta)$ and, at fixed $p$ and $q$,

$$\frac{\partial_rE}{E} = -\frac{S'\cos\theta + \left(P'\cos\phi + Q'\sin\phi\right)\sin\theta}{S},$$

Buckley and Schlegel's equation (15).
With $P$ and $Q$ constant the spacetime has the axis of symmetry $\theta = 0$, $\pi$, and their equation (13) becomes

$$ds^2 = -c^2dt^2 + \frac{\left(\partial_rR + R\,S'\cos\theta/S\right)^2}{1 + f}\,dr^2 + R^2\left(d\theta - \frac{S'\sin\theta}{S}\,dr\right)^2 + R^2\sin^2\theta\,d\phi^2.$$

`szekeres_pullback` checks slot by slot that this is the stereographic chart pulled back, in the variable $u = \cot(\theta/2)$, in which both metrics are rational.

With $f = 0$ the moment of $t$ is Euclidean space: the map $x = R\sin\theta\cos\phi$, $y = R\sin\theta\sin\phi$, $z = R\cos\theta + Z(r)$ pulls $dx^2 + dy^2 + dz^2$ back onto the spatial metric when $Z' = R\,S'/S$.
So the shell $r$ is the sphere of radius $R$ about the point of the axis at height $Z$, and $R\,S'/S$ is the rate at which the centres move along the axis toward $\theta = 0$.
`szekeres_fields` in `embedding.py` checks that pullback against the published metric.

## Step 5. The declared cloud

Every diagram takes the marginally bound cloud of Tolman-Bondi's spacetime diagram, $f = 0$, $R(0,r) = r$ and $2GM/c^2 = r^3(5 - 3r^2)/4$ inside $r_b = 1$ and $1/2$ beyond, each shell falling as $R^{3/2} = r^{3/2} - \tfrac{3}{2}\sqrt{2GM/c^2}\,ct$, with

$$\frac{S'}{S} = 2r\left(1 - r^2\right), \qquad S = e^{r^2 - r^4/2},$$

inside the cloud and $S$ constant beyond, where the spacetime is Schwarzschild's in Lemaître's coordinates.
The density is positive and no shells cross where $|S'/S| < M'/3M$ and $|S'/S| < \partial_rR/R$.
Here $M'/3M = 5(1 - r^2)/(r(5 - 3r^2))$, so $(S'/S)/(M'/3M) = 2r^2(5 - 3r^2)/5 \le 5/6$, and $r\,S'/S \to 0$ at the centre, which keeps the origin regular.
With $\mu = M'/3M$, $\partial_rR/R - S'/S$ has the sign of $\sqrt{r}\,(1 - \mu r) + (1 - S'/(S\mu))\,\mu R^{3/2}$, which is positive since the density falls outward, $\mu r \le 1$.

## Step 6. The spacetime diagrams

`null_rays.py` draws the two halves of the axis, $\theta = 0$ and $\theta = \pi$, whose null curves are null geodesics by the symmetry.
On them $c\,dt = \pm(\partial_rR \pm R\,S'/S)\,dr$, the upper sign toward $\theta = 0$, so the cones are narrower in $r$ on the half the centres move toward.
The marked ray on each is the last along that half to reach infinity, traced back from the event $R = r_s$ at $r = 1.1\,r_b$ of the exterior: it leaves the centre at $ct = -0.81\,r_b$ toward $\theta = 0$ and at $-0.12\,r_b$ toward $\theta = \pi$, and both cross the surface at $0.61\,r_b$.
The first is not called the event horizon, since an event of the centre between those two times still sends light out toward $\theta = \pi$.

No conformal diagram is drawn: the general spacetime has no symmetry to divide out, and on the axis the diagram is whatever the free functions make it, as for Tolman-Bondi.

## Step 7. The embedding diagram

The surface $\theta = \pi/2$ of a moment, through the equator of every shell, has the metric $\left((\partial_rR)^2 + R^2S'^2/S^2\right)dr^2 + R^2d\phi^2$ at $f = 0$, a surface of revolution with $\rho = R$ and $dz/dr = R\,S'/S$.
By Step 4 that is the surface as it lies in the flat space of the moment: each circle is the equator of a shell at the height of its centre.
At $ct = 0$ it is $z = 2r^3/3 - 2r^5/5$, which rises by $4/15$ from the centre to the surface of the cloud, and outside the cloud it is a plane.
The movie runs $ct$ from $-0.4$ to $0.55\,r_b$, the centre being crushed at $4/(3\sqrt{5}) = 0.596\,r_b$, with the plane at $z = 0$ in every frame and every height drawn three times over.
The surface meets the two halves of the axis only at the centre, so the spacetime diagrams do not mark its moments, which `slices.HIDDEN` records.
