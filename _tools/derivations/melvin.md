# Melvin's magnetic universe and Ernst's black hole

Melvin's universe in its cylindrical chart and Ernst's black hole in Schwarzschild's coordinates are

$$ds^2 = \Lambda^2\left(-c^2dt^2 + d\rho^2 + dz^2\right) + \frac{\rho^2}{\Lambda^2}\,d\phi^2,\qquad \Lambda = 1 + \frac{B^2\rho^2}{4},$$

$$ds^2 = \Lambda^2\left(-\left(1 - \frac{r_s}{r}\right)c^2dt^2 + \frac{dr^2}{1 - r_s/r} + r^2\,d\theta^2\right) + \frac{r^2\sin^2\theta}{\Lambda^2}\,d\phi^2,\qquad \Lambda = 1 + \frac{B^2r^2\sin^2\theta}{4}.$$

Both charts are written by `_tools/derivations/print_charts.py --metric melvin`, and `verify_metrics.py --system melvin/cylindrical` and `--system melvin/ernst` check each in seconds.

## Step 1. Units

The field enters only through $B = \sqrt{G}\,B_0/c^2$, with $B_0$ the field on the axis in Gaussian units, an inverse length, as Bini and Mashhoon write it, $1 + GB_0^2\rho^2/4c^4$.
Their length $a = 2c^2/(\sqrt{G}B_0) = 2/B$ is Kastor and Traschen's Melvin radius.
The vector potential is $A_\phi = (c^2/\sqrt{G})\,2Bs^2/(4 + B^2s^2)$, with $s = \rho$ or $s = r\sin\theta$ the distance from the axis, and its limit $2c^2/(\sqrt{G}B)$ as $s \to \infty$ makes the total flux $4\pi/B$ in units with $G = c = 1$, Kastor and Traschen's eq. (6).
`DIMENSIONS` declares $[B] = L^{-1}$ in both charts, and $[r_s] = L$ in Ernst's.

## Step 2. The field equations

In units with $G = c = 1$ and Gaussian fields, $G_{\mu\nu} = 8\pi T_{\mu\nu} = 2\left(F_{\mu\alpha}F_\nu{}^\alpha - \tfrac{1}{4}g_{\mu\nu}F_{\alpha\beta}F^{\alpha\beta}\right)$.
`melvin_maxwell` in `print_charts.py` takes $F = dA$ from the potential above and checks this in every slot of both charts, and $\partial_\mu(\sqrt{-g}\,F^{\mu\nu}) = 0$, before anything is written.
On Melvin's axis $G^t{}_t = -B^2$, the energy density $B_0^2/8\pi$ of the field, and $R = 0$ in both charts, since Maxwell's stress is trace free.

## Step 3. Melvin's geometry

The circle of constant $\rho$ has circumference $2\pi\rho/\Lambda$, whose derivative $(1 - B^2\rho^2/4)/\Lambda^2$ vanishes at $\rho = 2/B$, where the radius is $1/B$.
Beyond it the radius falls as $4/(B^2\rho)$ while the proper distance $\int\Lambda\,d\rho = \rho + B^2\rho^3/12$ grows without bound.
On the slice $t = 0$, $z = 0$, $d(\rho/\Lambda)/ds = (1 - B^2\rho^2/4)/\Lambda^3$ has magnitude below 1 at every $\rho > 0$, so the slice embeds in flat space everywhere.
The Kretschmann scalar is $16384B^4(3B^4\rho^4 - 24B^2\rho^2 + 80)/(4 + B^2\rho^2)^8$, $20B^4$ on the axis and falling as $\rho^{-8}$ far from it, so the spacetime has no curvature singularity.

## Step 4. Ernst's geometry

The block of $t$ and $r$ is $\Lambda^2$ times Schwarzschild's, so its null curves are Schwarzschild's, $ct = \pm r_* + $ const with $r_* = r + r_s\ln|r/r_s - 1|$, at every $\theta$.
Reflection in the equator and in the plane $\phi = 0$ fix the plane $\theta = \pi/2$, $\phi = 0$, and the rotations about the axis fix the axis, so both are totally geodesic and their null curves are null geodesics.
The horizon $r = r_s$ has area $\int\Lambda r_s\cdot r_s\sin\theta/\Lambda\,d\theta\,d\phi = 4\pi r_s^2$ for every $B$, and the factor $\Lambda^2$, static and regular at the horizon, leaves its surface gravity at Schwarzschild's $c^2/2r_s$, as Radu's unchanged temperature and entropy require.
On the equator the circles have radius $r/\Lambda$, greatest at $r = 2/B$ when $Br_s < 2$, and at the throat when $Br_s \ge 2$.

## Step 5. The diagrams

The spacetime diagrams draw Melvin's plane of $t$ and $\rho$ at $B = 1$ and Ernst's equator at $r_s = 1$, $B = 1/2$, with the Melvin radius and the widest circle marked; `null_rays.py --verify` compares their rays with $ct \pm \rho$ and with Schwarzschild's $ct \pm r_*$.
The conformal diagrams are Minkowski's half diamond for Melvin, $p, q = \arctan((ct \mp \rho)B)$, and Kruskal and Szekeres's hexagon for Ernst's equator, built by `Tower` on $f = 1 - r_s/r$ and checked against the published metric, whose conformal factor changes no null direction.
The embedding diagrams are Melvin's plane $z = 0$ out to $\rho = 4/B$ and Ernst's equator through the bifurcation sphere out to $r = 5\,r_s$ on both sheets, each circle checked against $\rho/\Lambda$ and $r/\Lambda$; the length of either spike grows as the cube of its coordinate, so neither is drawn farther.
