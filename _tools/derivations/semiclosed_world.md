# The semiclosed world, the bag of gold

The ball $\chi \le \chi_0$ of the closed Friedmann universe of dust with $\chi_0$ past the equator of the three sphere, $\pi/2 < \chi_0 < \pi$, and Schwarzschild's vacuum outside it.
This note records the four charts, where each comes from, how the two regions are joined, and what each diagram draws.
The id is `semiclosed_world`, spelled as Frolov, Markov and Mukhanov spell it in their title of 1990, since the collection's prose carries no hyphen in an English compound.

## Step 1: the sources

Hsu and Reeb's "Unitarity and the Hilbert space of quantum gravity" (2008), section (b), writes the construction out at the moment of time symmetry: the slice $a^2(d\chi^2 + \sin^2\chi\,d\Omega^2)$ of a closed universe at its greatest expansion, cut at $\chi = \chi_{1l}$, and the slice $U + V = 0$ of Kruskal's manifold with the Einstein-Rosen bridge in it.
Their (8) is the continuity of the area, $a\sin\chi_0 = r$, their (9) the continuity of its derivative along the proper distance, "which forces $\chi_{1l} \in [\pi/2, \pi)$", and their (10) the mass, $2M = a\sin^3\chi_0$.
They say the evolved spacetime of dust is Kruskal's and Friedmann's sewn together as in Oppenheimer and Snyder's collapse, and cite Misner, Thorne and Wheeler's section 32.4 for it.
Marolf's "Black holes, AdS, and CFTs" (2009) describes the same spacetime, a closed universe sewn onto the back of Kruskal's extension and reached through the Einstein-Rosen bridge, and says where the name bag of gold comes from.
Zelmanov's report to the International Astronomical Union (1965) credits the semiclosed world to Zel'dovich (1962) and to Novikov (1962, 1963).
Zel'dovich's own paper, ZhETF 43, 1037 (1962), could not be retrieved: the journal's archive lacks that issue, so it is not in the references and nothing is taken from it.
Novikov's paper of 1964 is behind the publisher's wall, so what is said of it is what Krasiński's editor's note says of it.

## Step 2: the four charts

`comoving` is the published interior of Oppenheimer and Snyder, $ds^2 = -c^2d\tau^2 + a^2(d\chi^2 + \sin^2\chi\,d\Omega^2)$ with $a(\tau)$ left free, on $\chi \le \chi_0$ with $\chi_0 > \pi/2$.
$\tau = 0$ is the moment of greatest expansion, and the dust runs from its bang at $-\tau_s$ to its crunch at $\tau_s = \pi a_m/2c$.

`conformal` writes the cycloid out, $a = \tfrac{1}{2}a_m(1 + \cos\eta)$ and $c\tau = \tfrac{1}{2}a_m(\eta + \sin\eta)$, the form the Oppenheimer-Snyder page uses, with $\eta$ from $-\pi$ to $\pi$.
With the scale factor explicit the chart is checked to be dust: $G^\eta{}_\eta = -3a_m/a^3$ and every other component zero, so $8\pi G\rho/c^2 = 3a_m/a^3$, and the Weyl tensor vanishes.
`semiclosed_world_cycloid` prints every value as a rational function of $\cos\eta$ and at most the first power of $\sin\eta$, factored, so that the scale factor shows as powers of $\cos\eta + 1$.

`schwarzschild` is Schwarzschild's $t$ and $r$, which cover one sheet of the exterior at a time: the far sheet, $r > r_s$ out to infinity, and the sheet behind the throat, $r_s < r \le R(t)$.
`isotropic` is the isotropic radius, areal radius $r(1 + r_s/4r)^2$, which covers both sheets together, the throat at $r_s/4$ and the surface of the dust at $r_d(t) < r_s/4$.
Both are checked to be vacuum and to be the published Schwarzschild metric carried along their maps.

No chart of Kruskal's or of Novikov's is published, though the conformal diagram is drawn in Kruskal's coordinates and the embedding in Novikov's slicing.
Kruskal's $r(U, V)$ and Novikov's $r(\tau, R^*)$ are implicit functions, and a chart whose components hold an undetermined function cannot be held to being a vacuum by the checker; Schwarzschild's own page publishes neither for the same reason.

## Step 3: the junction

`semiclosed_world_matching` checks, for any $\chi_0$ in $(0, \pi)$ with $r_s = a_m\sin^3\chi_0$, that the surface $R = a\sin\chi_0$ is a radial geodesic of the exterior with energy $E = \cos\chi_0$ per unit mass, $(dR/d\tau)^2 = E^2 - (1 - r_s/R)$, and that the areal radius changes along the unit outward normal at the rate $\cos\chi_0$ on both sides.
For $\chi_0 > \pi/2$ that rate is negative: the areal radius falls as one leaves the dust, which on the moment of time symmetry is the part of Flamm's paraboloid between the surface and the throat, on the sheet of the Kruskal manifold that an observer at infinity does not live on.
The energy $\cos\chi_0$ is negative for the same reason, the surface being a geodesic of the other exterior, whose static time runs backward.

The rest mass of the dust is $M_0 = \frac{3c^2a_m}{4G}(\chi_0 - \sin\chi_0\cos\chi_0)$, from $8\pi G\rho/c^2 = 3a_m/a^3$ and the volume $2\pi a^3(\chi_0 - \sin\chi_0\cos\chi_0)$ of the ball, and the mass seen from outside is $M = \frac{c^2a_m}{2G}\sin^3\chi_0$.
At $\chi_0 = 3\pi/4$ the ratio $M_0/M$ is $12.1$, the twelfth of the embedding diagram's caption, and as $\chi_0 \to \pi$ it grows without bound.

## Step 4: the diagrams

Every diagram is drawn at $\chi_0 = 3\pi/4$, where $a_m = 2\sqrt{2}\,r_s$ and the surface reaches $2\,r_s$: the complement on the three sphere of the star the Oppenheimer-Snyder page draws at $\chi_0 = \pi/4$, with the same mass.

The embedding diagram is a movie of the dust's proper time from $\tau = 0$ to $c\tau = 1.5\,r_s$.
Inside it is the sphere of radius $a(\tau)$ from its pole to $\chi_0$.
Outside it is Novikov's slice on both sheets, `embedding.NovikovSheets`: the shell labelled $s$ rests at $R = r_s(s^2 + 1)$ and falls on its cycloid, $s < 0$ behind the throat, $s = 0$ the throat and $s > 0$ the far sheet; it is Tolman-Bondi's chart with no dust and $E = -1/(2(s^2 + 1))$, Misner, Thorne and Wheeler's (31.12) with $s$ their $R^*$.
$g_{ss} = (\partial_sR)^2(s^2 + 1)/s^2$ is finite on the throat only as a ratio, so $\partial_sR$ is written with its factor $s$ explicit.
The throat's own shell reaches $r = 0$ at $c\tau = \pi r_s/2$, where the slice pinches off, so the movie stops at $1.5\,r_s$.
At $\tau = 0$ the surface is checked against Flamm's paraboloid in the Tolman-Bondi chart, in `isotropic` through the throat, and in `schwarzschild` on each sheet.

The conformal diagram is `conformal.Bag`.
The dust is the rectangle $0 \le \chi \le \chi_0$, $-\pi \le \eta \le \pi$.
The surface is Novikov's shell $s = \cot\chi_0$, in Kruskal's coordinates in closed form, `slices.novikov_sheets`, and the exterior is drawn in $P(U)$ and $Q(V)$ fixed by the surface and by the two singularities being the lines $T = \pm\pi$.
The bifurcation sphere is at $X = 3\chi_0 - \pi$, $i^\pm$ at $(3\chi_0, \pm\pi)$ and $i^0$ at $X = \pi + 3\chi_0$.
The event horizon $U = 0$ meets the surface at $\eta = \pi - 2\chi_0$ and the bang at $\chi = 3\chi_0 - 2\pi$, so for $\chi_0 > 2\pi/3$ no light from the centre of the dust reaches the outside.
Schwarzschild's $t$ is taken future directed on each sheet; the isotropic chart's $t$, a single static coordinate, runs to the past behind the throat, and the check takes $-\partial_t$ for the future there.

The spacetime diagrams are one for each chart.
The dust is drawn from bang to crunch, with `DustSolver`'s origin `rest`, added for it: the solution kept on both sides of a moment of time symmetry.
`schwarzschild` draws the far sheet, and `isotropic` both sheets from $t = 0$ on, with the surface of the dust behind the throat integrated from rest by `Surface` and the time function $-t$ there.
`null_rays.py --verify` checks the rays against $\eta \mp \chi$, $ct \mp r_*$ and $ct \mp r_*\,\mathrm{sgn}(4r - r_s)$.
