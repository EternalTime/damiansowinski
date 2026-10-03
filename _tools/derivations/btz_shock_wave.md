# The shock wave in the BTZ black hole

The three charts of `btz_shock_wave.json` are written by `print_charts.py --metric btz_shock_wave`, which computes every tensor from each line element with the checker's own `Geometry` and reads every printed value back before it writes.
This note records the source of each chart, how its step and its delta are read, and what each drawing draws.
The sources are S. H. Shenker and D. Stanford, JHEP 03 (2014) 067, arXiv:1306.0622, cited below as SS, and T. Dray and G. 't Hooft, Nucl. Phys. B 253 (1985) 173, with K. Sfetsos, Nucl. Phys. B 436 (1995) 721, for the shift across a null surface.

## Step 1. The black hole and its Kruskal chart

The eternal black hole of Bañados, Teitelboim and Zanelli without rotation is SS (6) and (7), $ds^2 = -(r^2 - R^2)dt^2/\ell^2 + \ell^2dr^2/(r^2 - R^2) + r^2d\phi^2$ with $R^2 = 8GM\ell^2$ in units with $c = 1$.
Its Kruskal chart is SS (8), $ds^2 = \left(-4\ell^2du\,dv + R^2(1 - uv)^2d\phi^2\right)/(1 + uv)^2$, with the right outside at $u < 0 < v$, the two boundaries on $uv = -1$ and the two singularities on $uv = 1$.
The map between them is SS (10): $r/R = (1 - uv)/(1 + uv)$ and $(v + u)/(v - u) = \tanh(Rt/\ell^2)$, which on the right outside is $u = -\sqrt{(r - R)/(r + R)}\,e^{-Rt/\ell^2}$ and $v = \sqrt{(r - R)/(r + R)}\,e^{Rt/\ell^2}$.
`btz_shock_check` holds the exterior chart to being the Kruskal chart pulled back along that map, slot by slot.
The exterior chart is printed with $c$ explicit, $Rct/\ell^2$ in place of $Rt/\ell^2$.

## Step 2. The shock

SS glue the black hole of mass $M$ to one of mass $M + E$ along the null surface $u_w = e^{-Rt_w/\ell^2}$, the path of quanta of energy $E$ let go from the left boundary a time $t_w$ before $t = 0$.
Continuity of the circle gives their (12), and in the limit $E/M \to 0$, $t_w \to \infty$ with $\alpha = (E/4M)e^{Rt_w/\ell^2}$ fixed, their (13), the shift $\tilde v = v + \alpha$ behind the shock, which lies on $u = 0$.
The Kruskal chart is SS (14), $ds^2 = \left(-4\ell^2du\,dv + R^2\left(1 - u(v + \alpha\Theta(u))\right)^2d\phi^2\right)/\left(1 + u(v + \alpha\Theta(u))\right)^2$, continuous across the shock.
The discontinuous chart is SS (15), in $U = u$ and $V = v + \alpha\Theta(u)$, $ds^2 = \left(-4\ell^2dU\,dV + 4\ell^2\alpha\,\delta(U)\,dU^2 + R^2(1 - UV)^2d\phi^2\right)/(1 + UV)^2$.
That is the form in which Dray and 't Hooft introduce a shock by a shift in $v$ at $u = 0$, their (B.2), and in which Sfetsos extends it to a background with a cosmological constant, his (2.5) and (2.6), with $F = -2Af\delta$; the black hole of three dimensions is his section 4.2, there with a single massless particle and here with a shell spread evenly round the horizon, for which the shift $f$ is the constant $\alpha$.
On either side of the shock the Kruskal chart is the black hole's, in $u$ and $v$ ahead of it and in $u$ and $v + \alpha$ behind it, which `btz_shock_check` holds.

## Step 3. The step and the delta

The Kruskal chart defines the step as a name, $\Theta = \tfrac{1}{2}(1 + \mathrm{sgn}(u))$, as Penrose's impulsive wave does, so the checker reads its metric with $|u|$ and $\mathrm{sgn}(u)$ and reduces a delta of $u$ times a function continuous at $u = 0$ to that function's value there, `_on_a_kink`.
Every delta of the curvature multiplies a function of $u(v + \alpha\Theta)$, which is continuous, so no reduction meets a function that jumps.
`btz_shock_shifted` prints each value once, as its value ahead of the shock with $v + \alpha\Theta$ written for $v$, after checking that the value behind the shock is that function at $v + \alpha$, and closes it with its delta.
The discontinuous chart's delta multiplies $1/(1 + UV)^2$, which varies across the shock, so the chart is listed in `IMPULSES` and its delta read as a pulse even about the shock, as Hotta and Tanaka's Kruskal chart's is.

## Step 4. What is checked

Both Kruskal charts have $R_{\mu\nu} = -(2/\ell^2)g_{\mu\nu}$ off the shock and $R_{uu} + (2/\ell^2)g_{uu} = 2\alpha\,\delta(u)$ on it, which with $\Lambda = -1/\ell^2$ and $c = 1$ is SS (16), $T_{uu} = \alpha\,\delta(u)/4\pi G$, "a shell of null particles symmetrically distributed on the horizon".
`btz_shock_check` holds every slot to it.
The Ricci scalar is $-6/\ell^2$ and the Kretschmann scalar $12/\ell^4$ everywhere, the shock included, since its Ricci tensor is null, and the Weyl tensor vanishes, as it does in every spacetime of three dimensions.
The published $\Gamma^V{}_{UU} = -\alpha\left(\delta'(U) + 2V\delta(U)\right)$ integrated twice across $U = 0$ makes a ray of constant $u$ jump by $\alpha$ in $V$, which is $V = v + \alpha\Theta(u)$, and `conformal.py` checks the jump from the published value.

## Step 5. The spacetime diagrams

Every drawing is at $\ell = R = 1$ and $\alpha = 1$.
On the plane of $u$ and $v$ only $g_{uv}$ is nonzero, so the rays are $u$ and $v$ constant, and no Christoffel symbol turns a ray out of the plane in any chart.
The Kruskal view hatches $|u(v + \alpha\Theta)| \ge 1$ and marks the shock, the past horizon $v = 0$ of the right outside and the future horizon $v = -\alpha$ of the left one, which miss each other by $\alpha$.
The discontinuous view draws the delta as the pulse $p(U) = e^{-U^2/w^2}/(w\sqrt\pi)$, $w = 0.05$; a ray moving left then follows $dV/dU = \alpha p(U)(1 + UV)^2$, which `_btz_shock_before` integrates back across the pulse, so `--verify` holds each such ray to the $V$ it had before it.
The exterior view is the black hole's, and its rays are $ct \mp r_*$ with $r_* = (\ell^2/2R)\ln|(r - R)/(r + R)|$.

## Step 6. The conformal diagram

In the Kruskal chart $u$ and $v$ are null on both sides of the shock, so $p = \arctan u$ and $q = \arctan(v + \alpha/2)$ bring the plane into a finite drawing with light at 45° and each ray moving left drawn as one line across the shock, the diagonal $X = T$.
The shift $\alpha/2$ makes the drawing symmetric under $(X, T) \to (-X, -T)$ and puts the embedded moment on $T = 0$.
The map is checked on all three charts, the exterior one on both outsides, the left outside with $u$ and $v + \alpha$ for $-v$ and $-u$.

## Step 7. The moment the embedding diagram draws

On the moment $u + v = -\alpha/2$ of the Kruskal chart $u(v + \alpha\Theta(u)) = \alpha|u|/2 - u^2$, so its metric is the same function of $w = |u|$ on both sides of the shock:
$$ds^2 = \frac{4\ell^2dw^2}{(1 + \alpha w/2 - w^2)^2} + \rho^2d\phi^2, \qquad \rho = R\,\frac{1 - \alpha w/2 + w^2}{1 + \alpha w/2 - w^2}.$$
It is the black hole's Kruskal time $-\alpha/2$ ahead of the shock and its own Kruskal time $\alpha/2$ behind it, and it crosses the shock at $v = -\alpha/2$, halfway between the horizons, where the geodesic of SS (17) to (19) from $t = 0$ on the left boundary to $t = 0$ on the right one crosses it, since $d_1$ and $d_2$ there are the same functions of $v + \alpha$ and $-v$.
With $ds = 2\ell\,dw/(1 + \alpha w/2 - w^2)$, $d\rho/ds = -R(\alpha - 2w)/\ell(1 + \alpha w/2 - w^2)$: the circles shrink at $\alpha R/2\ell$ leaving the shock, so the two halves meet there at an angle, and the slice reaches a neck at $w = \alpha/4$ of radius $R(16 - \alpha^2)/(16 + \alpha^2)$.
It crosses the horizons at $w = 1/2$, where $\rho = R$, and lies level where $(d\rho/ds)^2 = 1$, at $w = (\sqrt{33} - 3)/4$ for $\alpha = 1$; beyond that it is drawn in three dimensional Minkowski space, as the black hole's $t = 0$ is beyond $r = \sqrt2\,\ell$, out to $w = 1$, where $\rho = 3R$.
For $\alpha \ge 2\ell/R$ the circles leave the shock faster than the distance out, and for $\alpha \ge 4$ the moment meets $r = 0$, so $\alpha = 1$ is drawn.
In the exterior chart the moment's part outside either horizon is the curve $r = R(1 - uv)/(1 + uv)$, $ct = (\ell^2/2R)\ln(v/(-u))$ with $u < -\alpha/2$ and $v = -\alpha/2 - u$, the same for both outsides, and `slices.py` draws it once.
