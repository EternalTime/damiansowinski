# The RP3 geon

The three charts of `rp3_geon.json` are written by `print_charts.py --metric rp3_geon`, which computes every tensor from each line element with the checker's own `Geometry` and reads every printed value back before it writes.
This note records the source of each chart, the identification, what was read in each source, and what the diagrams draw.

## Step 1. The spacetime

Kruskal's manifold is the maximal extension of Schwarzschild's solution (`kruskal1960`).
In Kruskal's coordinates $(T, X, \theta, \phi)$, with $T^2 - X^2 < 1$, the metric is

$$ds^2 = \frac{4r_s^3}{r}e^{-r/r_s}\left(-dT^2 + dX^2\right) + r^2d\Omega^2, \qquad \left(\frac{r}{r_s} - 1\right)e^{r/r_s} = X^2 - T^2,$$

which is `louko1998` (4.5) and (4.6) with $2M = r_s$.
The map

$$J: (T, X, \theta, \phi) \to (T, -X, \pi - \theta, \phi + \pi)$$

is `louko1998` (4.12a) and `giulini2018` (38).
It is an isometry, since the metric holds $X$ only through $X^2$ and $dX^2$ and the antipodal map is an isometry of the sphere.
It has no fixed point, since no point of the sphere is its own antipode.
Its Jacobian is $\mathrm{diag}(1, -1, -1, 1)$, which leaves $\partial_T$ alone and has determinant $+1$, so the quotient is orientable in time and in space.
The quotient is the RP3 geon, named in `friedman1993`.
`RP3_GEON_INVOLUTION` in `print_charts.py` holds the map, and `rp3_geon_check` holds each of these statements.

A fundamental domain is the half $X \ge 0$, whose edge $X = 0$ is glued to itself by the antipodal map of the sphere, which the Kruskal chart states as its last domain.
Each point of the edge is a sphere with opposite points one, a projective plane, of area $2\pi r^2$.
The sphere $T = X = 0$, where the horizons cross, has $r = r_s$, so its projective plane has area $2\pi r_s^2 = 8\pi M^2$, half of the horizon's $16\pi M^2$ (`louko1998`, section IV B).

## Step 2. The charts

- `kruskal`: as above, on $X \ge 0$. The areal radius is $r = r_s(1 + \mathrm{W}((X^2 - T^2)/e))$ with $\mathrm{W}$ Lambert's function, held as a function while the tensors are built (`HELD` in `verify_metrics.py`), with $\partial_T r = -2r_s^2Te^{-r/r_s}/r$ and $\partial_X r = 2r_s^2Xe^{-r/r_s}/r$. `rp3_geon_kruskal` writes every derivative of $r$ by those two, tests a value for zero with the exponential written $(r - r_s)/(r_s(X^2 - T^2))$, and writes it back with every product $(T - X)(T + X)$ as $(r_s - r)e^{r/r_s}/r_s$, so that no value divides by $T$ or $X$.
- `schwarzschild`: the one exterior, `louko1998` (4.8) and (4.9): $T = \sqrt{r/r_s - 1}\,e^{r/2r_s}\sinh(ct/2r_s)$ and $X$ the same with $\cosh$.
- `isotropic`: `misner1957` (234), with their $m^* = r_s/2$, on the one sheet $\rho > r_s/4$; the areal radius is $\rho(1 + r_s/4\rho)^2$.

The fold reverses the Killing vector $(X\partial_T + T\partial_X)/2r_s$, which is $\partial_{ct}$ in the exterior, so the geon is static only outside its horizons.
The moment $t$ of the exterior is the line $T = X\tanh(ct/2r_s)$, and its image under the fold is $T = -X\tanh(ct/2r_s)$, so only $t = 0$ meets the edge level and extends to a smooth slice of the whole spacetime (`louko1998`, section I).

## Step 3. What Misner and Wheeler wrote

The candidate list gave the attribution to `misner1957` as from memory, so the paper was read: pages 590 to 594 of the scan at liphy-annuaire.univ-grenoble-alpes.fr, the section "The Schwarzschild and Reissner-Nordstrom solutions".
Their (234) is Schwarzschild's metric in the isotropic radius $\rho$, and (236) to (238) take the moment $T = 0$ as initial data with any lapse $V$ that does not vanish.
Their (245), at no charge, is the inversion $\rho \to m^{*2}/4\rho$, which is $\rho \to r_s^2/16\rho$, and (246) to (248) complete it by $T \to T$, $\theta \to \pi - \theta$, $\phi \to \phi + \pi$.
They write that the mapping "has no fixed points, and we may identify the point $x$ with the point $Jx$ to obtain a new manifold with only one region, $\rho \to +\infty$, which is asymptotically flat", and that "The surface $T = 0$ in this manifold is topologically equivalent to projective 3-space $P^3$ (see Sec. III B) with one point removed".
So the paper bears out the attribution for the moment of time symmetry.
It also says the procedure fails for the Reissner-Nordström solution, since "the electromagnetic field changes sign", and cites a Princeton senior thesis of 1957 by Morton R. Dubman for the first order of the evolution, in which "the area of the critical sphere begins to decrease".

Kruskal's coordinates are of 1960, so the fold as a map of the whole spacetime is later.
Taken with Schwarzschild's static time on both sheets, $(t, \rho) \to (t, r_s^2/16\rho)$ is $(T, X) \to (-T, -X)$, the elliptic interpretation of `gibbons1986elliptic`, which is how `giulini2018`, footnote 11, reads Misner and Wheeler's map; with a lapse of one sign, as their (236) takes it, $T \to T$ keeps the direction of time and the development of the data is the geon.
The two maps agree on the moment $T = 0$, which is all the paper constructs, and the History claims no more of it.
`rp3_geon_check` holds the inversion to being an isometry of the isotropic chart and to being $X \to -X$ on the moment $t = 0$.
`louko2005`, footnote 1, and `louko2010`, section 2, give the same credit for the initial data.

## Step 4. Topological censorship

`friedman1993` use the geon as their example.
A ray moving left keeps $V = T + X$, and one from the past null infinity of the exterior has $V > 0$.
It meets the edge at $T = V$, where $T^2 - X^2 = V^2 > 0$, inside the black hole, and leaves it moving right with $U = T - X = V > 0$, on which $T^2 - X^2 = V(V + 2X)$ only grows, so it ends on the singularity.
The time reverse leaves the past singularity, meets the edge inside the white hole and reaches future null infinity, the passive detection of the topology that the theorem allows.
`_tools/test_rp3_geon.py` holds both.

## Step 5. The diagrams

- Spacetime diagrams: Kruskal's plane on $X \ge 0$, with the edge, both horizons and the singularities $T^2 - X^2 = 1$; Schwarzschild's plane of the exterior; and the isotropic plane from $\rho = r_s/4$ out. `--verify` checks the rays against $T \pm X$ and $ct \mp r_*$.
- Conformal diagram: the right half of Kruskal's diagram, by $p = \arctan(T - X)$, $q = \arctan(T + X)$, in three views, each with the glued edge and the two rays of Step 4.
- Embedding diagram: the moment $t = 0$ on the equatorial plane, Flamm's paraboloid from the throat out, one sheet. The throat circle is glued to itself point to opposite point, so it is a closed geodesic of length $\pi r_s$ and the surface is a Möbius band; two opposite meridians are marked as one geodesic through the throat.

## Sources as read

`louko1998`, `friedman1993`, `louko2005`, `louko2010` and `giulini2018` were read in their arXiv versions, and every quotation in the History was matched word for word against them.
`misner1957` was read in the scan named above; its record is Crossref's.
`gannon1975` and `sorkin1986` are cited for what `friedman1993` and `louko2010` say of them, and their records are Crossref's.
No relation to `misner` is stated, since no source read ties Misner space to the geon; the page relates the geon to `misner_brill_lindquist` and `einstein_rosen_bridge`, both through `misner1957`.
