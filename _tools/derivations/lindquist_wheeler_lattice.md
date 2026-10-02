# The lattice universe of Lindquist and Wheeler

Richard Lindquist and John Wheeler's lattice (Rev. Mod. Phys. 29, 432, 1957) is $N$ equal masses at the centres of the cells of a regular tiling of the three sphere, each cell replaced by a ball with Schwarzschild's geometry inside it.
It is an approximate spacetime, so what is published is the geometry of one cell in the three charts its literature uses, and the comparison hypersphere the cells are fitted to.
The charts are written by `_tools/derivations/print_charts.py --metric lindquist_wheeler_lattice`, which took 3 seconds on 2 October 2026, and `verify_metrics.py --system lindquist_wheeler_lattice/<chart>` checks each in under a second.

## Step 1. The source of each chart

The 1957 paper is behind a paywall, so each line element is taken from a paper that quotes it or derives it, and Lindquist's own account is the open report of the Chapel Hill conference (doi:10.34663/9783945561294-17).

- `schwarzschild_cell`: Schwarzschild's chart about the mass of one cell, Clifton and Ferreira's (3) (Phys. Rev. D 80, 103503, arXiv:0907.4109), with $2m$ written $r_s$.
- `cosmological_time`: Clifton and Ferreira's (8), reached from Schwarzschild's chart by their (7), $d\tau = \sqrt{E}\,dt - \sqrt{E - 1 + 2m/r}\,dr/(1 - 2m/r)$, in which every shell with $(dr/d\tau)^2 = E - 1 + 2m/r$ has $\tau$ for its proper time and is orthogonal to the surfaces of constant $\tau$. For the closed lattice $E = \cos^2\psi$.
- `lindquist_wheeler`: Lindquist and Wheeler's comoving coordinates as Liu writes them (Phys. Rev. D 92, 063529, arXiv:1501.05169), his (17) to (19): each shell is labelled by the largest radius it reaches, his $\tilde r_{LW}$, here $\rho$, and $\tau$ is its proper time from that moment. The line element is Tolman and Bondi's with no dust, $-c^2d\tau^2 + (\partial_\rho r)^2d\rho^2/(1 - r_s/\rho) + r^2d\Omega^2$, and the areal radius $r(\tau, \rho)$ is left free, as Tolman-Bondi's $R$ is, since the cycloid has no closed form in $\tau$.
- `comparison_hypersphere`: the three sphere of radius $a(\tau)$ the cells are tangent to, with $a$ left free. Clifton and Ferreira's (4) is its equation, $\dot a^2 = 2m/(a\sin^3\psi) - 1$.

The angle $\psi$ is the angular radius of a cell of equal volume, $N(2\psi - \sin 2\psi) = 2\pi$, Liu's (35).
The boundary of a cell has $r = a\sin\psi$, and the tangency of the cell to the hypersphere there, Liu's (40), is $\sqrt{1 - r_s/r_m} = \cos\psi$, so the boundary turns around at $r_m = r_s/\sin^2\psi$ and the hypersphere at $a_m = r_s/\sin^3\psi$, Clifton and Ferreira's (5).
Friedmann's closed universe of the same total mass $Nm$ turns around at $4Nm/3\pi$, their (6), so the ratio is $3\pi/(2N\sin^3\psi)$: 1.428, 1.277, 1.156, 1.114, 1.036 and 1.012 for $N$ = 5, 8, 16, 24, 120 and 600.

## Step 2. What the script checks

`lindquist_wheeler_check` holds the cell in Schwarzschild's chart to the published metric of `schwarzschild` slot by slot and to being a vacuum.
It holds the cosmological time chart to being that metric pulled back along Clifton and Ferreira's (7), to a vanishing Ricci tensor and $K = 12r_s^2/r^6$, and the vector $(1, \sqrt{E - 1 + r_s/r})$ to being a unit normal of the surfaces of constant $\tau$.
The comoving chart is a vacuum only where every shell obeys $(\partial_\tau r)^2 = r_s/r - r_s/\rho$, so the check writes $\partial_\tau^2 r = -r_s/2r^2$, $\partial_\tau\partial_\rho r$ and $\partial_\tau^2\partial_\rho r$ by that equation, reduces every power of $\partial_\tau r$ by it, and then finds the Ricci tensor zero and the Kretschmann scalar Schwarzschild's.
The hypersphere is held to $G_{\tau\tau} = 3r_s/(a^3\sin^3\psi)$ and no other component where $a$ obeys its equation: dust at rest, which is what Friedmann's equation with that constant says.

`metric_tags.py` counts the two exact cell charts as the lattice's region, so the page carries `vacuum`, and overrules `static`, `stationary` and `spherically symmetric`, which those charts show and the lattice of cells does not have.

## Step 3. The drawings

Every drawing is of the lattice of eight cells, the tiling by cubes that Bentivegna and Korzyński evolved, with $\psi = 0.8832$, $r_m = 1.6746\,r_s$ and $a_m = 2.1671\,r_s$; `slices.py` holds those numbers for all three scripts.

The spacetime diagrams draw one cell in each chart and the hypersphere.
`CellBoundary` in `null_rays.py` is the boundary: in Schwarzschild's chart the radial geodesic from rest at $r_m$ and its mirror image in $t = 0$, and in the cosmological time chart, which is singular at $r_m$, the geodesic of the same energy started a part in $10^8$ inside $r_m$ and run back to the bang.
The comoving chart takes its areal radius from `lw_r`, the cycloid with its derivatives along $\tau$ and $\rho$ in closed form in the cycloid's parameter, and the row's `bang` asks for the curve where the metric stops being finite below the middle of the drawing as well as above it.
The hypersphere's radius is solved from its own published $G^\chi{}_\chi = 0$, from rest, on both sides of the moment of rest.
`--verify` holds the rays of each view to closed forms the drawing does not use: $t \pm r_*$ in a cell, through Clifton and Ferreira's (7) and through Novikov's $t$ of the comoving shells (Misner, Thorne and Wheeler's (31.10)), and $\eta \pm \chi$ on the hypersphere.

The embedding diagram is the picture Lindquist showed at Chapel Hill, funnels tangent to a hypersphere: the equatorial plane through two opposite cells, each cell the slice of constant $\tau$ of the comoving chart, where $dz/dr = 1/\sqrt{\rho/r_s - 1}$, and between them the hypersphere of radius $a(\tau)$.
The two meet with one tangent at every moment, which the script checks to $10^{-9}$, and at $\tau = 0$ the cell is Flamm's paraboloid.
A shell $\rho$ reaches $r = 0$ at $c\tau = \tfrac{\pi}{2}\rho\sqrt{\rho/r_s}$, so from $c\tau = \pi r_s/2$ on each funnel ends in a point.
The movie runs from the moment of rest to $c\tau = 3.36\,r_s$; the expansion is its time reverse.

The conformal diagram of a cell is the part of Kruskal and Szekeres's diagram between the throat's world line, $X = 0$, and the boundary, with the comoving shells placed by $V = (k\cos\tfrac{\eta}{2} + \sin\tfrac{\eta}{2})\exp\left(\tfrac{1}{2}(r + k(\eta + \tfrac{\rho}{2}(\eta + \sin\eta)))\right)$, $k = \sqrt{\rho/r_s - 1}$, and $U(\eta) = -V(-\eta)$.
The cosmological time chart is placed by $u = (\tau - \tau_m)/\sqrt{E} + G(r)$ with $dG/dr = r(W/\sqrt{E} - 1)/(r - r_s)$, regular at the horizon, and covers the cell before the slice $\tau = \tau_m$, which is the gap Liu describes between the expanding and the contracting charts.
The hypersphere is the rectangle of its cycloid parameter and $\chi$.

The moments of the embedding are Novikov's slices in Schwarzschild's chart, where the two that still reach outside $r_s$ are drawn, lines of constant $\tau$ in the comoving chart and on the hypersphere, and no part of the cosmological time chart, which covers the cell while it expands.

## Step 4. What the tests hold

`_tools/test_lindquist_wheeler_lattice.py` holds the published files to the numbers the texts state: the six angles and ratios, the on-shell vacuum of the comoving chart and the dust of the hypersphere, the tangency $\sqrt{1 - r_s/r_m} = \cos\psi$ in every frame of the movie, the boundary's turn and its proper time in the drawings, and the cycloid's solver.
