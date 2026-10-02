# The Simpson-Visser black bounce

Simpson and Visser's metric is Schwarzschild's with $\sqrt{r^2 + a^2}$ written for the areal radius and $dr$ left as it is.
This note records the four charts, where each comes from, the curvature, the tortoise coordinate in closed form, and what each diagram draws.
Every chart writes $r_s = 2Gm/c^2$ for the Schwarzschild radius of Simpson and Visser's mass parameter $m$, and $\rho = \sqrt{r^2 + a^2}$ for the areal radius.

## Step 1: Simpson and Visser's chart

Their line element is $ds^2 = -f\,c^2dt^2 + dr^2/f + (r^2 + a^2)\,d\Omega^2$ with $f = 1 - r_s/\sqrt{r^2 + a^2}$, and $t$ and $r$ both run over the whole line.
The sphere $r = 0$ has the least area, $4\pi a^2$.
For $a < r_s$, $f$ vanishes at $r = \pm h$ with $h = \sqrt{r_s^2 - a^2}$, two simple roots, and between them $r$ is the time: the black bounce.
For $a = r_s$, $f = r^2/2r_s^2 + \dots$ has a double root at $r = 0$: the one way wormhole, an extremal horizon on the sphere of least area.
For $a > r_s$, $f \ge 1 - r_s/a > 0$: the traversable wormhole.
At $a = 0$ the metric is Schwarzschild's, and at $r_s = 0$ it is Ellis and Bronnikov's wormhole with $\ell = a$.

The chart names $\rho = \sqrt{r^2 + a^2}$ as a parameter defined by the coordinates, so every value is printed in $r$, $a$ and $\rho$.
`print_charts.simpson_visser_radius` writes each power of $r^2 + a^2$ as a power of $\rho$, factors the value with $a^2$ written as $\rho^2 - r^2$, so that $r$ and $\rho$ are generators with no relation between them, and then writes each factor in whichever of its three forms has the fewest terms.
That is how $2r^2 - a^2$ and $\rho - r_s$ come out as Simpson and Visser write them.

## Step 2: the areal radius

Tsukamoto's standard radial coordinate is $\rho = \sqrt{r^2 + a^2}$, in the subsection of his paper that carries that name: $ds^2 = -(1 - r_s/\rho)\,c^2dt^2 + d\rho^2/\left((1 - r_s/\rho)(1 - a^2/\rho^2)\right) + \rho^2d\Omega^2$.
In it $g_{tt}$ and the spheres are Schwarzschild's, which is why the photon sphere, the innermost stable circular orbit and the shadow are Schwarzschild's while $a$ is smaller than their radii.
The chart covers one side of $r = 0$ at a time, $r = \pm\sqrt{\rho^2 - a^2}$, and $g_{\rho\rho}$ diverges on $\rho = a$ with the curvature finite.
It is the Morris-Thorne form with $e^{2\Phi} = 1 - r_s/\rho$ and $b = r_s + a^2/\rho - r_sa^2/\rho^2$, so $b(a) = a$.

## Step 3: the Eddington-Finkelstein charts

Simpson, Martín-Moruno and Visser write the static metric as $ds^2 = -f\,dw^2 \mp 2\,dw\,dr + (r^2 + a^2)\,d\Omega^2$, in the introduction and second section of their paper, the upper sign for the retarded time $u$ and the lower for the advanced time $v$.
With $u = ct - r_*$ and $v = ct + r_*$, $dr_*/dr = 1/f$, each is Simpson and Visser's chart pulled back, and each is regular on the horizons.
The ingoing chart runs from an exterior through the black hole, across $r = 0$ and out through the second horizon; the outgoing chart is its time reverse.

`print_charts.simpson_visser_check` pulls Simpson and Visser's metric back through each of the three maps, the areal one on both sides of $r = 0$, and compares it with the chart's own, slot by slot.

## Step 4: the curvature

Simpson and Visser's Einstein tensor, Ricci scalar and Kretschmann scalar, from their sections on the curvature tensors and invariants, with $m = r_s/2$:
$G^t{}_t = a^2(\rho - 2r_s)/\rho^5$, $G^r{}_r = -a^2/\rho^4$, $G^\theta{}_\theta = G^\phi{}_\phi = a^2(2\rho - r_s)/2\rho^5$, $R = a^2(3r_s - 2\rho)/\rho^5$, and
$K = \left(3r_s^2(4r^4 - 4a^2r^2 + 3a^4) + 16a^2r_s(r^2 - a^2)\rho + 12a^4\rho^2\right)/\rho^{10}$.
The script checks each of these in every chart, the Ricci tensor against zero at $a = 0$, and $K$ against Ellis and Bronnikov's $12a^4/\rho^8$ at $r_s = 0$.
On the sphere of least area $K = (9r_s^2 - 16ar_s + 12a^2)/a^6$, which is $256/r_s^4$, $5/r_s^4$ and $25/64r_s^4$ at $a = r_s/2$, $r_s$ and $2r_s$.
Outside a horizon $\rho_{\rm matter} + p_r = (G^r{}_r - G^t{}_t)c^4/8\pi G = -a^2(\rho - r_s)c^4/4\pi G\rho^5 < 0$, so the null energy condition fails everywhere off the horizons.

## Step 5: the tortoise coordinate

$dr_*/dr = \rho/(\rho - r_s) = 1 + r_s(\rho + r_s)/(r^2 + a^2 - r_s^2)$, and $\int \rho\,dr/(r^2 + a^2 - r_s^2)$ and $\int dr/(r^2 + a^2 - r_s^2)$ are elementary.
With $r_*(0) = 0$ where that is finite, $r_*$ is odd in $r$:

- $a < r_s$, $h = \sqrt{r_s^2 - a^2}$: $r_* = r + r_s\,\mathrm{arsinh}(r/a) + \dfrac{r_s^2}{2h}\ln\left|\dfrac{(r - h)(r_sr - h\rho)}{(r + h)(r_sr + h\rho)}\right|$.
  Both factors of the numerator vanish at $r = h$ and both of the denominator at $r = -h$, so $r_* \to (r_s^2/h)\ln|r \mp h|$ there, which is $1/2\kappa$ times the logarithm with $\kappa = h/2r_s^2$, Simpson and Visser's surface gravity, the same on both horizons.
- $a = r_s$: $r_* = r + r_s\,\mathrm{arsinh}(r/r_s) - r_s(r_s + \rho)/r$, which diverges as $-2r_s^2/r$.
- $a > r_s$, $k = \sqrt{a^2 - r_s^2}$: $r_* = r + r_s\,\mathrm{arsinh}(r/a) + \dfrac{r_s^2}{k}\left(\arctan\dfrac{r}{k} + \arctan\dfrac{r_sr}{k\rho}\right)$, finite at every $r$.

`null_rays._sv_rstar` holds the three, and `null_rays.py --verify` checks the traced rays of every view against $ct \mp r_*$, $v - 2r_*$ and $u + 2r_*$.
`conformal.py` checks $dr_*/dr$ against the published $g_{rr}$ numerically in each case.

## Step 6: the diagrams

Every diagram is drawn at $r_s = 1$ for the three geometries $a = r_s/2$, $r_s$ and $2r_s$.

The spacetime diagrams are one view for each geometry in each chart.
Where $a < r_s$ the region between the horizons takes its future from the ingoing chart, toward smaller $r$.

The conformal diagram of $a = r_s/2$ is `BounceTower`, a `Tower` whose tortoise coordinate is the closed form above: the cells are Carter's for the axis of Kerr, with the regions $r < -h$ asymptotically flat.
On $r = 0$, $r_* = 0$, so $u = v = ct$ and $G(u) + G(-v) = \arctan e^{-\kappa ct} + \arctan e^{\kappa ct} = \pi/2$: the sphere of least area is the straight line $T = \pi/2$ through the middle of the black hole.
Simpson and Visser's chart covers the cells I, II and III, the areal radius I and the half of II below $r = 0$, the ingoing chart I, II and III', and the outgoing chart III', the white hole above it and the exterior above that.
The diagram of $a = r_s$ has no surface gravity to build Kruskal's exponential from, so each region is placed by arctangents of $(ct \mp r_*)/\ell$ with $\ell = 2r_s$, as Majumdar and Papapetrou's hole is: the regions $r > 0$ on the right, the regions $r < 0$ on the left, each a whole diamond, alternating upward.
The diagram of $a = 2r_s$ is the full diamond by $\arctan((ct \mp r_*)/\ell)$ with $\ell = 4r_s$, the throat on its axis.

The embedding diagram has four views.
Outside the horizon of $a = r_s/2$ the equator of $t = 0$ is read in the areal chart, where $g_{\rho\rho}$ is rational and the horizon is the root $\rho = r_s$: $dz/d\rho = \sqrt{g_{\rho\rho} - 1}$, steeper than Flamm's paraboloid at every radius, on both sheets through the bifurcation sphere.
Between the horizons a moment of constant $r$ has the equator $(r_s/\rho - 1)c^2dt^2 + \rho^2d\phi^2$, a flat cylinder of radius $\rho$ on which the stretch $|ct| \le r_s$ is $2r_s\sqrt{r_s/\rho - 1}$ long.
The movie runs in $r$ from $0.75\,r_s$ to $-0.75\,r_s$, each frame's value the proper time since the horizon of an observer at fixed $t$, $c\tau = \int_r^h \sqrt{\rho/(r_s - \rho)}\,dr$, which is $1.831\,r_s$ at $r = 0$ and $3.663\,r_s$ from horizon to horizon.
For $a = r_s$ the equator of $t = 0$ has $g_{rr} \to 2r_s^2/r^2$, an infinitely long throat closing on the radius $r_s$, drawn from $r = r_s/50$.
For $a = 2r_s$ it is one piece through the throat, with $dz/dr = \sqrt{\rho/(\rho - r_s) - r^2/\rho^2}$, which is $\sqrt{a/(a - r_s)}$ at the throat and $a/\rho$, the catenoid's, at $r_s = 0$.
The heights have no elementary closed form, so each is checked against scipy's quadrature of its slope.
