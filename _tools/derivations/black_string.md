# The black string and its Gregory-Laflamme instability

The black string is Schwarzschild's black hole with one more direction of space, $z$, along which nothing changes: $ds^2 = -f\,c^2dt^2 + dr^2/f + r^2d\Omega^2 + dz^2$ with $f = 1 - r_s/r$.
Gregory and Laflamme's general metric is a Schwarzschild-Tangherlini black hole of $D$ dimensions times $10 - D$ flat directions, and a chart needs its coordinates named one by one, so the collection publishes the string of five dimensions in three charts and the string of six, Tangherlini's black hole of five times a line, in one.

## Step 1: the parameter

Reduce the string of five dimensions along a circle of length $L$ and the four dimensional Newton constant is $G_4 = G_5/L$, Gregory's (1.10), while the metric of the four dimensions is Schwarzschild's with the string's whole mass $M$.
So $r_s = 2G_4M/c^2 = 2G_5M/(c^2L)$, the mass per unit length standing where the mass stood.
The linearised field agrees: a string of energy $\lambda$ per unit length and tension $\tau$ has $\nabla^2h_{tt} = -16\pi G_5(2\lambda - \tau)/3$ on the three flat directions across it, and the uniform black string has $\tau = \lambda/2$, which gives $h_{tt} = 2G_5\lambda/(c^2r)$.
In six dimensions the same reduction gives Tangherlini's $r_h^2 = 8G_5M/(3\pi c^2)$ with $G_5 = G_6/L$.

## Step 2: the charts

The static chart is Gregory and Laflamme's (1) at $D = 4$ with one flat direction.
`black_string_product` in `print_charts.py` checks that its components along $z$ are those of $dz^2$, that nothing depends on $z$, and that the rest is the published metric of `schwarzschild.json`, slot by slot, and the same for the chart of six dimensions against `tangherlini.json`.

The ingoing Eddington-Finkelstein chart is their (3.16) of 1994, $ds^2 = -f\,dv^2 + 2\,dv\,dr + r^2d\Omega^2 + dz^2$ with $v = ct + r + r_s\ln|r/r_s - 1|$, which they call $u$ and in which they find the apparent horizon of the rippled string.

The Kerr-Schild chart is the form Choptuik, Lehner, Olabarrieta, Petryk, Pretorius, and Villegas base their coordinates on, their (8): $ds^2 = -f\,c^2dT^2 + (2r_s/r)\,c\,dT\,dr + (1 + r_s/r)\,dr^2 + r^2d\Omega^2 + dz^2$.
Its time is $cT = v - r$, so $c\,dT = c\,dt + r_s\,dr/(r - r_s)$, and $g^{TT} = -(1 + r_s/r)$ is negative at every $r$: the surfaces of constant $T$ are spacelike through the horizon, which is what an evolution with the singularity excised needs.
`black_string_pullback` checks both charts to be the static one pulled back, in every slot.

Gregory and Laflamme judge a perturbation regular in Kruskal's coordinates, $R$ and $T$ of their (2.14).
The areal radius is an implicit function of those, as it is for Schwarzschild's black hole, and the collection publishes no chart of either in them; the conformal diagram is drawn with them.

## Step 3: the curvature

No curvature component has an index along $z$, so the Riemann tensor is Schwarzschild's, the Ricci tensor vanishes, the Weyl tensor is the Riemann tensor, and the Kretschmann scalar is $12r_s^2/r^6$, or $72r_h^4/r^8$ in six dimensions.

## Step 4: the instability

`gregory_laflamme.py` holds the perturbation.
With $h_{ab} = e^{\Omega t + ikz}H_{ab}(r)$, spherically symmetric and with no component along $z$, Gregory's review (arXiv:1107.5821, (1.21) to (1.23)) reduces the linearised vacuum equations to two first order equations and one algebraic one for $H = H_{tr}$ and $H_\pm = H_{tt}/V \pm VH_{rr}$, with $h_{\theta\theta} = r^2H_-/2$ from the vanishing trace.
`linearised_ricci()` checks them without assuming a gauge: it expands the Ricci tensor of $g + \epsilon h$ to first order in sympy, puts the three equations in, and every one of the fifteen components simplifies to zero.

The mode regular on the future horizon behaves as $H \sim (\Omega r_s - \tfrac12)(r - r_s)^{\Omega r_s - 1}$, and the one that decays far away as $e^{-\sqrt{\Omega^2 + k^2}\,r}$.
`mismatch()` integrates the first outward and the second inward, in $H/\Omega$ and $H_-$, which leaves only $\Omega^2$ in the equations, and takes their Wronskian where they meet; an unstable mode is a zero of it.
In units of $r_s$ and $c$:

- the threshold, where $\Omega \to 0$, is $kr_s = 0.876$, Gregory and Laflamme's value, a wavelength of $7.17\,r_s$;
- the fastest mode has $kr_s = 0.354$, a wavelength of $17.8\,r_s$, and $\Omega = 0.0923\,c/r_s$;
- on Lehner and Pretorius's circle of length $10\,r_s$, their $L = 20M$ with $r_s = 2M$, the one mode that fits has $\Omega = 0.0633\,c/r_s$.

The results do not move in the fifth digit when the start is taken from $10^{-4}$ to $10^{-10}\,r_s$ outside the horizon or the meeting point from $1.5$ to $6\,r_s$.

## Step 5: the diagrams

The plane of $t$ and $r$ at fixed angles and fixed $z$ is Schwarzschild's, so the spacetime diagrams are Schwarzschild's in each chart and the conformal diagram is Kruskal and Szekeres's, each point a sphere times the line of $z$; in six dimensions both are Tangherlini's of five.
In the Kerr-Schild chart the ingoing rays are $cT + r = $ const and the outgoing ones $cT + r - 2r_* = $ const, $c\,dT/dr = (r + r_s)/(r - r_s)$.

The embedding diagram has three views.
Across the string, the equator of a moment of $t$ at one $z$ is Flamm's paraboloid, and in six dimensions the catenoid.
Along the string, the horizon at one advanced time is the surface of $z$ and $\phi$ at $r = r_s$ in the ingoing chart, a cylinder of radius $r_s$.
Gregory and Laflamme's (3.20) puts the apparent horizon of the perturbed string at $r_s$ plus a term in $\cos kz$, and the mode has $h_{zz} = 0$, so to first order the horizon's cross section has the metric $dz^2 + r(z)^2d\phi^2$ with $r = r_s(1 + a\cos kz)$ and $a \propto e^{\Omega v/c}$.
That is the metric of the surface $r = r(z)$ of constant $v$ in the published ingoing chart, where $g_{rr} = 0$, so `Slice` reads it there with `along`, and the surface climbs at $\sqrt{1 - (dr/dz)^2}$.
The view plays it as a movie from $a = 0.02$ at $v = 0$ to $a = 0.29$ at $v = 42\,r_s$, a frame every $r_s$.
It is linear theory, and the caption says where it stops holding; the perturbed string is another spacetime than the one published, so its moments are marked on no other drawing.
