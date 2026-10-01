# The spinning cosmic string: the charts and what the diagrams draw

The spinning string is flat everywhere off the string, so every exterior chart of `spinning_string.json` has vanishing curvature and the physics sits in two numbers: $b = 1 - 4G\mu/c^2$, which removes a wedge, and the length $a = 4GJ/c^3$, which shifts time on each circuit.
This note records where each chart comes from and what the diagrams take from it.
Units here have $c = 1$; the file keeps $c$.

## Step 1. The proper radius chart

Deser and Jackiw's line element of 1992 (their equation 7, with $c$ and $G$ written out) is

$$ds^2 = -\left(dt + a\,d\phi\right)^2 + dr^2 + b^2r^2d\phi^2 + dz^2,$$

with $\phi$ periodic in $2\pi$.
Its $g_{\phi\phi} = b^2r^2 - a^2$ changes sign at $r_c = a/b$, which is their $r < 4GJ/(c^3 - 4G\mu c)$ for the closed timelike curves.
`print_charts.spinning_string("proper_radius")` builds it and the checker's `Geometry` finds every curvature tensor zero.

## Step 2. The helical chart

With $\tau = t + a\phi$ and $\tilde\phi = b\phi$ (their equation 9) the line element is Minkowski's, $-d\tau^2 + dr^2 + r^2d\tilde\phi^2 + dz^2$.
The identification $\phi \sim \phi + 2\pi$ at fixed $t$ becomes $(\tau, \tilde\phi) \sim (\tau + 2\pi a, \tilde\phi + 2\pi b)$: the angle runs over $2\pi b$, and $\tau$ jumps by $2\pi a = 8\pi GJ/c^4$ on each circuit.
`spinning_string_pullback` pulls the proper radius chart back through $t = \tau - a\tilde\phi/b$, $\phi = \tilde\phi/b$ and checks every slot against the helical chart, and the helical metric against $\mathrm{diag}(-1, 1, r^2, 1)$.

## Step 3. The rescaled and circumference radii

Cornish and Frankel (1994) write the same exterior with $\rho = br$, $ds^2 = -(dt + a\,d\phi)^2 + d\rho^2/b^2 + \rho^2d\phi^2$, their equation (1.2), which is also the form Mena, Natário, and Tod match their shell to, with $C = 1/b$ and $m = a$.
Their equation (1.6) takes the radius of the circle itself, $R^2 = \rho^2 - a^2$, so that $g_{\phi\phi} = R^2$ and $g_{RR} = R^2/b^2(R^2 + a^2)$.
That chart covers $r > r_c$ only, and $R = 0$ is the null circle.
Both are checked to be the proper radius chart pulled back, through $r = \rho/b$ and $r = \sqrt{R^2 + a^2}/b$.

## Step 4. The extended source

Soleng (1994) writes a thick string in the frame $\omega^0 = dt + M\,d\phi$, $\omega^1 = dr$, $\omega^2 = \rho\,d\phi$, $\omega^3 = dz$, with $M$ and $\rho$ functions of $r$.
Without torsion his Einstein tensor (his equation 4 at $\sigma = 0$) is, with $\Omega = M'/2\rho$,

$$G^0{}_0 = -3\Omega^2 + \rho''/\rho, \qquad G^1{}_1 = G^2{}_2 = \Omega^2, \qquad G^3{}_3 = -\Omega^2 + \rho''/\rho, \qquad G^0{}_2 = -\Omega'.$$

`spinning_string_source` checks that frame orthonormal for the chart's metric, rotates the chart's $G^\mu{}_\nu$ into it, and checks each of those components, the heat flow up to the sign of the frame's orientation.
At $M = a$ and $\rho = br$ the metric is the proper radius chart's, which is checked too.
The radial pressure is $\Omega^2$, so a surface with no pressure has $M' = 0$ there, and smoothness on the axis asks $M(0) = 0$; Soleng's argument is that $\Omega$ then has $\Omega' \neq 0$ somewhere, a heat flow, which breaks the weak energy condition where the closed timelike curves are.

## Step 5. The values drawn

Every diagram is drawn at $b = 0.9$, the deficit the cosmic string is drawn at, and $a = 0.9$, so that $r_c = a/b$ is the unit of length.

On a cylinder of $t$ and $\phi$ at one radius the metric is $-(dt + a\,d\phi)^2 + b^2r^2d\phi^2$, the same at every point, and its null curves are $dt = (br - a)\,d\phi$ and $dt = -(br + a)\,d\phi$: $-0.45$ and $-1.35$ at $r_c/2$, $0.45$ and $-2.25$ at $3r_c/2$, the slopes `CYLINDERS` holds the traced rays to.
In the circumference radius chart $br = \sqrt{R^2 + a^2}$, so at $R = a$ the slopes are $a(\sqrt2 - 1)$ and $-a(\sqrt2 + 1)$.
In the helical chart the cylinder is flat, the rays run at 45° against $r\tilde\phi$, and the circle of constant $t$ is the straight line $\tau = a\tilde\phi/b$, which is steeper than the rays inside $r_c$ and shallower outside.
The spacetime is flat, so $\Gamma^r{}_{\phi\phi} = -b^2r$ turns light launched along any of these curves off the cylinder toward larger $r$, which `verify_turning` checks.

## Step 6. The conformal diagram

The half plane of fixed $\phi$ and $z$ has the metric $-dt^2 + dr^2$, and it is totally geodesic, since every $\Gamma^\phi{}_{ab}$ and $\Gamma^z{}_{ab}$ with $a$, $b$ among $t$ and $r$ vanishes.
It is not orthogonal to the circles of $\phi$, because $g_{t\phi} = -a$, so the published $g^{tt} = -(b^2r^2 - a^2)/b^2r^2$ is not the inverse of the half plane's own metric, and $t \mp r$ are null coordinates of the half plane without being null coordinates of the spacetime, whose null hypersurfaces are $\tau \mp r$.
`conformal.Plane(..., induced=True)` reads such a surface: it checks the published Christoffel symbols for total geodesy and takes the inverse of the surface's own metric.
In the helical chart the half plane of fixed $\tilde\phi$ is the same surface, with $\tau = t$ on $\tilde\phi = 0$, and needs none of that.
Each view is Minkowski's half diamond, $p, q = \arctan((t \mp r)/r_c)$, each point a single event, with the region $r < r_c$ tinted.

## Step 7. The embedding diagram

The moment $t = 0$ of the plane $z = 0$ has the metric $dr^2 + (b^2r^2 - a^2)\,d\phi^2$, positive only outside $r_c$.
In the circumference radius it is $g_{RR}\,dR^2 + R^2d\phi^2$ with $g_{RR} = R^2/b^2(R^2 + a^2)$, which is below one inside $R_1 = ab/\sqrt{1 - b^2}$ and above it outside; in the proper radius $R_1$ is $r_c/\sqrt{1 - b^2} = 2.29\,r_c$.
Inside $R_1$ the circles grow faster than the distance out to them and the surface is drawn in Minkowski space, $dZ/dR = \sqrt{1 - g_{RR}}$, which is one at $R = 0$: the surface leaves the axis along the light cone.
Outside it is drawn in flat space, $dz/dR = \sqrt{g_{RR} - 1} \to \sqrt{1/b^2 - 1}$, the cosmic string's cone.
Both slopes vanish at $R_1$, so the two pieces meet with one tangent.

Near the axis the profile is nearly null, and there a straight chord of Minkowski space is longer than the arc it spans.
With $Z \approx R - R^3/6a^2b^2$ the chord from $R_1$ to $R_2$ has length $(R_2 - R_1)\sqrt{(R_1^2 + R_1R_2 + R_2^2)/3}/ab$ and the arc $(R_2^2 - R_1^2)/2ab$, a ratio of $1 + (R_2 - R_1)^2/6(R_1 + R_2)^2$ to leading order, and $2/\sqrt3$ for the chord from the axis itself however short it is.
So the profile steps toward the axis by a fortieth of $R$ at a time, is written to twelve decimals, and starts at $r = (1 + 2\times10^{-5})\,r_c$.
