# The travelling wave on a cosmic string: the charts and what the diagrams draw

Garfinkle and Vachaspati's travelling wave is a pp-wave whose transverse plane is the cone of a cosmic string, and `string_wave.json` gives it in three charts.
This note records where each chart comes from, the maps between them, and what the diagrams take from them.
Units here have $c = 1$, with $u = t - z$ and $v = t + z$ both lengths, and $b = 1 - 4G\mu/c^2$.

## Step 1. The sources

Garfinkle's paper of 1990 (Physical Review D 41, 1112) and Garfinkle and Vachaspati's (Physical Review D 42, 1960) are not free to read, so the line elements are taken from sources that quote them and are free: Malcolm Anderson's "Near-field expansion of the metric due to a cosmic string" (1999), whose equations (1.5) to (1.10) give Garfinkle's solution and its second form, and Dabholkar, Gauntlett, Harvey, and Waldram's "Strings as solitons & black holes as strings" (1996), whose section 2.3 reviews the transformation and whose (2.13) and (2.14) are the same two forms for the fundamental string.
The transformation itself is read from Kaloper, Myers, and Roussel (1997), section 2, and from Garfinkle's own account of it in "Black string traveling waves" (1992), equations (9) to (12).
The abstracts of the two papers of 1990, read at INSPIRE, carry what the History says of them.

## Step 2. The null conical chart

Garriga and Peter's equation (14) (1994) is the pp-wave on a cone,

$$ds^2 = -du\,dv + F\,du^2 + dr^2 + b^2r^2d\phi^2,$$

with their $-H$ for $F$ and their $\gamma$ for $b$.
The checker's `Geometry` finds one Ricci component,

$$R_{uu} = -\tfrac{1}{2}\left(\partial_r^2F + \frac{\partial_rF}{r} + \frac{\partial_\phi^2F}{b^2r^2}\right),$$

minus half the Laplacian of $F$ on the cone, so the metric is a vacuum solution exactly where $F$ is harmonic there, which is Frolov and Garfinkle's family.
The chart leaves $F$ free, as the pp-wave's Brinkmann chart leaves its $H$, and no component assumes a field equation.
The harmonic functions regular on the string are $r^{m/b}\cos m\phi$ and $r^{m/b}\sin m\phi$; the string's own wave is $m = 1$.

## Step 3. Garfinkle's isotropic chart

Anderson's (1.2) writes the cone conformally flat: with $x = (br)^{1/b}\cos\phi$ and $y = (br)^{1/b}\sin\phi$ the transverse metric is $\rho^{-8\mu}(dx^2 + dy^2)$, $\rho^2 = x^2 + y^2$, and $-8\mu = 2b - 2$.
The checker reads a power with a symbolic exponent with its base in a fixed unit, so the file writes the base as a pure number, $(\rho/\ell)^{2b-2}$, with $\ell$ an arbitrary length: $r = (\ell/b)(\rho/\ell)^b$.
Anderson's (1.5) and (1.6), in the signature $(+,-,-,-)$, are $ds^2 = dt^2 - dz^2 + F(dt - dz)^2 - \rho^{-8\mu}(dx^2 + dy^2)$ with $F = 2xA'' + 2yB''$, which in the collection's signature is

$$ds^2 = -du\,dv - 2\left(xA'' + yB''\right)du^2 + \left(\frac{\rho}{\ell}\right)^{2b-2}\left(dx^2 + dy^2\right).$$

It is written in $u$ and $v$ because the checker's `Reader` takes a function of a coordinate, and $A$ is a function of $t - z$.
`string_wave_check` checks that this is the null conical chart pulled back at $F = -2(xA'' + yB'') = -2\ell(br/\ell)^{1/b}(A''\cos\phi + B''\sin\phi)$, at random points since the map carries the power $1/b$, that the profile is harmonic on the cone, and that the isotropic chart's Ricci tensor vanishes for every $A$ and $B$.
Its curvature is $R_{uxux} = -R_{uyuy} = (1 - b)(xA'' - yB'')/\rho^2$ and $R_{uxuy} = (1 - b)(xB'' + yA'')/\rho^2$: nothing at $b = 1$, where the linear profile is a change of coordinates, and curved for every $b < 1$.
The Kretschmann scalar and the Ricci scalar vanish, as Kaloper, Myers, and Roussel's theorem says every scalar invariant's difference from the cone's must.

A Christoffel symbol differentiates $g_{uu}$ once more, so the chart needs $A'''$, and the `Reader` now reads a third prime.
`norm` now also splits a power of a reciprocal, since sympy turns $(\rho/\ell)^{2b-2}$ into $(1/\ell)^{2b}$ times the rest.

## Step 4. The moving string chart

Anderson's (1.7) to (1.9) and Dabholkar, Gauntlett, Harvey, and Waldram's (2.14) change coordinates so that the string is seen to move:

$$X = x + A, \qquad Y = y + B, \qquad V = v + 2xA' + 2yB' + \int\left(A'^2 + B'^2\right)du .$$

Then $-du\,dv - 2(xA'' + yB'')du^2 = -du\,dV + 2A'\,dx\,du + 2B'\,dy\,du + (A'^2 + B'^2)du^2$, and with $dx = dX - A'du$ the line element is

$$ds^2 = -du\,dV + dX^2 + dY^2 + \left(\left(\frac{\rho}{\ell}\right)^{2b-2} - 1\right)\left(\left(dX - A'du\right)^2 + \left(dY - B'du\right)^2\right),$$

with $\rho^2 = (X - A)^2 + (Y - B)^2$, which is Anderson's (1.10) in the collection's signature.
At $b = 1$ it is Minkowski's line element as it stands, and the string lies along $X = A(u)$, $Y = B(u)$.
`string_wave_check` checks it against the isotropic chart pulled back, at random points.
The chart's $\rho$ holds $A$ and $B$, so the checker holds it as a function and `reduce` writes it out, as Szekeres's $E$ is held.
Far from the string $g_{uu} = ((\rho/\ell)^{2b-2} - 1)(A'^2 + B'^2)$ falls as $-2(1 - b)(A'^2 + B'^2)\ln(\rho/\ell)$ for a light string, a logarithmic potential whose gradient is Vachaspati's $1/r$ force.

## Step 5. What the diagrams declare

Every diagram is drawn for one pulse, $A = (\ell/2)e^{-4u^2/\ell^2}$ and $B = 0$, on a string with $b = 1/2$, whose cone lacks half a turn; at the $b = 0.9$ the cosmic string is drawn at, the conformal factor differs from 1 by a few hundredths and nothing would show.
In the null conical chart that is $F = -2(r/2)^2(32u^2 - 4)e^{-4u^2}\cos\phi$ in units of $\ell$.

The spacetime diagrams are planes of $u$ and $v$ at fixed places across the string.
On each the rays of constant $u$ are null geodesics, the generators of the null Killing vector, and the other family, $dv/du = g_{uu}$, is a null curve that light crossing the wave leaves, since $\Gamma^i{}_{uu} \neq 0$ across the string.
The time function is $3u + v$, whose gradient is timelike wherever $g_{uu} > -3$; $u + v$ stops being one below $g_{uu} = -1$, and $g_{uu}$ reaches $-2$ on the far side of the string.
`--verify` checks the curves against $v \pm \tfrac{1}{2}\ell A'$ in the null conical and isotropic charts, and against $V - \int g_{uu}\,du$ by a quadrature in the moving string chart, where the shifts across the whole pulse are $0.65\,\ell$ and $0.78\,\ell$.

No conformal diagram is drawn: inside the pulse the plane of $u$ and $v$ at a fixed place across the string is not totally geodesic, so a diagram of it would be read as a picture of the causal structure, which its second family of curves does not carry.
The pp-wave has none for the same reason.

The embedding diagram is the surface of constant $u$ and $v$, the cone of half angle $\arcsin b = 30°$ at every $u$.
What the wave does is shown on it by a ring of 360 free particles at rest on $r = 2\ell$ before the pulse: `string_wave_ring` runs them with the published Christoffel symbols, $u$ being affine since no $\Gamma^u$ is published, once in the null conical chart and once in the isotropic chart, and the two rings are checked to be one through $r = (\ell/b)(\rho/\ell)^b$.
The ring crosses the cone as the string swings out and back, reaching from $r = 1.5\,\ell$ to $2.5\,\ell$ on the crest, and is left falling toward the string, its largest radius $1.70\,\ell$ at $u = \ell$ and $1.36\,\ell$ at $u = 2.5\,\ell$.
`StringWaveRing` in the tests runs twelve of the particles again with no sympy, from the isotropic chart's geodesic equations written out, and holds the file's marked particles to them.
