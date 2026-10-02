# Bartnik and McKinnon's soliton

Bartnik and McKinnon's soliton of 1988 is a static, spherically symmetric ball of SU(2) Yang-Mills field held together by its own gravity, regular at its centre, with no horizon, and flat far away.
No solution is known in closed form.
Its four charts are written by `_tools/derivations/print_charts.py --metric bartnik_mckinnon` with their functions left free, `verify_metrics.py --system bartnik_mckinnon/<chart>` checks each in seconds, and every drawing rests on the numerical solution of `_tools/derivations/bartnik_mckinnon.py`.

## Step 1. The charts and their sources

The Physical Review Letter itself could not be read when this was written, on 2 October 2026, so no equation number below is Bartnik and McKinnon's own.
The equations are taken from Volkov and Gal'tsov's review (Phys. Rep. 319, 1, 1999, arXiv:hep-th/9810070), whose numbers these are, and Smoller, Wasserman, Yau and McLeod's announcement (Bull. Amer. Math. Soc. 27, 239, 1992) records that Bartnik and McKinnon's own line element is $ds^2 = -T^{-2}dt^2 + A^{-1}dr^2 + r^2d\Omega^2$.

- `areal`: the review's (2.50) with the signature reversed, $ds^2 = -\sigma^2N\,c^2dt^2 + dr^2/N + r^2d\Omega^2$ with $N = 1 - 2m/r$.
  The gauge potential is its (2.49) with no electric part, $A = w\,(T_2\,d\vartheta - T_1\sin\vartheta\,d\varphi) + T_3\cos\vartheta\,d\varphi$, so the field is one function $w(r)$: $w = \pm 1$ is a pure gauge and $w = 0$ the Dirac monopole.
- `isotropic`: Kleihaus and Kunz's (60) (Phys. Rev. D 57, 834, 1998, arXiv:gr-qc/9707045), $ds^2 = -f\,dt^2 + (m/f)(dr^2 + r^2d\Omega^2)$, the spherical case $l = m$ of their (7).
  Their $m$ is written $h$ and their $r$ is written $\rho$ here, since $m$ and $r$ are the mass function and the areal radius in the other charts.
  Their (61) and (62) are $h/f = r^2/\rho^2$ and $d\rho/\rho = dr/(r\sqrt{N})$, and their (65) fixes the constant by $\rho/r \to 1$ at infinity.
- `tortoise`: the coordinate of the review's (5.6), $d\xi/dr = 1/(\sigma N)$, in which it writes the equation of the soliton's small oscillations, its (5.5).
  The review calls it $\rho$; it is $\xi$ here, since $\rho$ is the isotropic radius and $x$ the line through the centre in the drawings.
  The line element is $F(-c^2dt^2 + d\xi^2) + r^2d\Omega^2$ with $F = \sigma^2N$.
- `flow`: Breitenlohner, Forgács and Maison's (48) (Commun. Math. Phys. 163, 141, 1994; equation numbers are those of the preprint MPI-Ph/93-41), $ds^2 = A^2N^2dt^2 - r^2(d\tau^2 + d\Omega^2)$ with $N = \sqrt{\mu}$, $\mu = 1 - 2m/r$, and $\tau$ defined by $dr = rN\,d\tau$.
  In it the field equations are their (50), a flow with no explicit $\tau$, which is how they classify every solution with a regular centre.
  Their $N$ is the square root of the review's, and their $A$ is its $\sigma$.
  The chart leaves $A$, $N$ and $r$ free, so its printed values hold for any three functions, and $\partial_\tau r = rN$ is stated with its coordinates.

## Step 2. The field equations

With a prime for $d/dr$ and every length in units of $\ell = \sqrt{4\pi G}/(gc^2)$, the review's (2.29), its (3.2) to (3.4) are

$$m' = Nw'^2 + \frac{(1 - w^2)^2}{2r^2},\qquad r^2Nw'' + \left(2m - \frac{(1 - w^2)^2}{r}\right)w' + w(1 - w^2) = 0,\qquad \frac{\sigma'}{\sigma} = \frac{2w'^2}{r}.$$

The mass comes in units of $c^2\ell/G = \sqrt{4\pi/G}/g$.
With $V = (1 - w^2)^2/(2r^2)$ the stress tensor is diagonal,

$$G^t{}_t = -\frac{2\ell^2(Nw'^2 + V)}{r^2},\qquad G^r{}_r = \frac{2\ell^2(Nw'^2 - V)}{r^2},\qquad G^\theta{}_\theta = G^\phi{}_\phi = \frac{2\ell^2V}{r^2},$$

and its trace vanishes, as a Yang-Mills field's does, so the Ricci scalar is zero on every solution.
The first two equations are $G^t{}_t$ and $G^r{}_r - G^t{}_t$; $G^\theta{}_\theta$ is never used to build a solution and holds by the Bianchi identity.

`bartnik_mckinnon_check` in `print_charts.py` holds each chart to this before it is written.
It runs along the chart's own radial coordinate: the areal radius, $m$, $\sigma$, $w$ and $w'$ are five functions of it whose derivatives are the equations above times $dr/dx$, the chart's functions are written by them, and what is left of each mixed Einstein component after its stress is subtracted is a rational function in those five and $\sqrt{N}$, which has to cancel to zero.
Each chart after the first is also held to the areal radius $\sqrt{g_{\theta\theta}}$, to the Misner-Sharp mass $1 - |\nabla r|^2 = 2m/r$, and to the areal chart's lapse and radial part.

## Step 3. The numerical solution

A regular centre leaves one free number, $b$, the review's (3.5), and an asymptotically flat end two, $M$ and $a$, its (3.6):

$$w = 1 - br^2 + \left(\tfrac{3}{10}b^2 - \tfrac{4}{5}b^3\right)r^4 + \dots,\quad m = 2b^2r^3 - \tfrac{8}{5}b^3r^5 + \dots;\qquad w = \pm\left(1 - \frac{a}{r} + \frac{3a(a - 2M)}{4r^2} + \dots\right),\quad m = M - \frac{a^2}{r^3} + \dots$$

`bartnik_mckinnon.py` carries the first through $r^6$ in $w$ and $r^7$ in $m$, and the second through $1/r^5$ and $1/r^7$; `_tools/test_bartnik_mckinnon.py` holds both to the field equations symbolically.

For a general $b$ the solution shot outward from the centre leaves the strip $|w| < 1$, and the number of zeros $w$ has by then steps up by one as $b$ passes each $b_n$.
`bracket` finds $b_n$ by bisection on that count.
Shooting outward cannot reach the flat end, since the perturbation of $w$ about $\pm 1$ that grows is $r^2$, so each soliton is then solved from both ends, as the review describes: from the centre outward with $b$, from far out inward with $M$ and $a$, the two meeting in $w$, $w'$ and $m$ just beyond the last zero of $w$, which Newton's method brings about in the three numbers.
The far radius is 200 times the radius of the last zero, since the zeros of the higher solitons stand a factor of about 40 apart.

| $n$ | $b_n$ | $M_n$ | $a_n$ | $\sigma_n(0)$ | least $N$ | at $r$ | zeros of $w$ |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 0.4537162727 | 0.8286469821 | 0.89338 | 0.126431 | 0.242384 | 1.5975 | 1.5457 |
| 2 | 0.6517255256 | 0.9713454943 | 8.86391 | 0.020892 | 0.035059 | 1.1302 | 1.1001, 3.6957 |
| 3 | 0.6970400503 | 0.9953164722 | 58.9326 | 0.003392 | 0.002974 | 1.0187 | 0.9898, 1.6663, 14.164 |

The review's Table 1 gives $b_n$, $M_n$, $a_n$ and $\sigma_n$ to four figures, and Breitenlohner, Forgács and Maison's Table 1 gives $b_n$ and $M_n$ to eleven: $0.45371627277$ and $0.82864698216$, $0.65172552552$ and $0.97134549426$, $0.69704005033$ and $0.99531647219$.
The test holds the three solitons to both tables, the second within $2 \times 10^{-10}$.

Three more functions ride along the integration: $\delta = -\ln\sigma$, zero at infinity; $\ln(\rho/r)$, from $d\ln\rho/dr = 1/(r\sqrt{N})$, zero at infinity, where the isotropic radius is Schwarzschild's; and $\xi$, from $d\xi/dr = 1/(\sigma N)$, zero at the centre.
At the centre $\rho/r = 0.2357$ for the first soliton, so $r = 4.24\,\rho$ there.

A drawing asks for these functions point by point, so `Soliton._tables` puts all of them on cubic splines in $\ln r$, a part in 900 apart and twenty times finer about the neck, and the test holds the splines to the integrations.
Inside $10^{-3}\,\ell$ every function and its derivatives are the series at the centre, differentiated term by term, so that nothing is divided by a small radius.

## Step 4. What the drawings declare

Every drawing is in units of $\ell$ and, unless it says otherwise, of the soliton with one zero, `BM_INPUT`.
`null_rays.py` declares the functions of each chart from the solver's, as sympy functions that carry their own derivatives, `DECLARED_FUNCTIONS`:

- `areal`: $m$ and $\sigma = e^{-\delta}$.
- `isotropic`: $f = \sigma^2N$ and $h = fk^2$ at $r = \rho k$, with $k(\rho) = r/\rho$ and its two derivatives, $k' = -2(m/r^2)k^2/(1 + \sqrt{N})$, which is finite at the centre.
- `tortoise`: $F = \sigma^2N$ at $r(\xi)$, with $dr/d\xi = \sigma N$.
- `flow`: $A = \sigma$, $N = \sqrt{1 - 2m/r}$ and $r$ at $\rho = \ell e^\tau$, which fixes the constant in $\tau$ as $\tau = \ln(\rho/\ell)$.

The test evaluates the published Einstein tensor and Ricci scalar of each chart with these functions and holds them to the stress of the Yang-Mills field on the same spheres, which needs every declared function and both of its derivatives to be right.

Each spacetime diagram marks the sphere on which $w = 0$, $r = 1.5457\,\ell$.
The conformal diagram is one view for all four charts, whose radial coordinates are functions of one another: $p, q = \arctan((ct \mp \xi)/8\ell)$, Minkowski's triangle, checked against each chart's published metric.
The embedding diagram has a view for each of the first three solitons, $dz/dr = \sqrt{2m/(r - 2m)}$, with Flamm's paraboloid of the same mass beside it.
Its steepest slope is $\sqrt{(1 - N)/N}$ at the least $N$: $1.77$, $5.2$ and $18.3$, the neck that in the limit $n \to \infty$ becomes the infinite throat of the extreme Reissner-Nordström black hole, $N = (1 - \ell/r)^2$, the review's (3.10).
Only the first soliton's moment lies on the spacetime diagrams and the conformal diagram, which draw that soliton alone.
