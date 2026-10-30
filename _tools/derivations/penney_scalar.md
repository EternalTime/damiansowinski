# Penney's charged scalar solution

The static, spherically symmetric, asymptotically flat solution of Einstein's equations with a Maxwell field and a massless, minimally coupled scalar field, which R. Penney published in 1969.
This note records the two charts, where each comes from, the field equations each is held to, and what each diagram draws.
Units here have $G = c = 1$ where a source does; the metric file keeps $c$ and writes the mass and the charge as lengths.

## Step 1: the spherical chart

Banerjee and Sen's (2.14) to (2.19) and Astorino's (2.25) to (2.29) write the solution in a radius with $g_{tt}g_{rr} = -1$:
$e^{\bar\alpha} = e^{-\bar\gamma} = (r - a)^{-\Lambda}(r - b)^{-\Lambda}\left[(b(r - a)^\Lambda - a(r - b)^\Lambda)/(b - a)\right]^2$ and $e^{\bar\beta} = e^{\bar\alpha}(r - a)(r - b)$.
With $E = ((r - a)/(r - b))^\gamma$ the bracket is $(r - b)^\gamma(a - bE)/(a - b)$, so

$$ds^2 = -\frac{(a - b)^2E}{(a - bE)^2}c^2dt^2 + \frac{(a - bE)^2}{(a - b)^2E}\left(dr^2 + (r - a)(r - b)\,d\Omega^2\right).$$

The parameters are lengths $a > b > 0$ and $0 < \gamma \le 1$, Banerjee and Sen's $\Lambda$ and Astorino's $A$.
$E$ is held as a function of $r$ with the slope $dE/dr = \gamma(a - b)E/((r - a)(r - b))$, `HELD` and `RATES` in `verify_metrics.py`, so every value is a rational function of $r$, $E$ and the parameters.
A second name such as $D = (a - bE)/(a - b)$ cannot be given: the reader holds any name defined through a held one, and $D$ and $E$ would then stand in a relation the rates leave out, so the sum $a - bE$ is printed as it stands.

## Step 2: the mass and the charge

Far away $-g_{tt} = 1 - \gamma(a + b)/r + O(r^{-2})$, so $2GM/c^2 = \gamma(a + b)$.
Banerjee and Sen's (2.19b) and Astorino's text both write $a + b = 2m$ and call $m$ the mass; that holds only at $\gamma = 1$, and their own (2.15) gives $\gamma(a + b)$.
Bronnikov's (30), $m^2 = q^2 + h^2$, agrees with $\gamma(a + b)/2$ under Step 4's dictionary.
The charge enters as $q^2 = \gamma^2ab$, with $q^2 = GQ^2/(4\pi\varepsilon_0c^4)$, in the normalisation where Reissner and Nordström's $-g_{tt}$ is $1 - r_s/r + r_q^2/r^2$; Banerjee and Sen's $2\Lambda^2ab = q^2$ carries another normalisation of $q$.

## Step 3: the field equations

With $R_{\mu\nu} = 2\,\partial_\mu\psi\,\partial_\nu\psi + 2F_{\mu\alpha}F_\nu{}^\alpha - \tfrac{1}{2}g_{\mu\nu}F_{\alpha\beta}F^{\alpha\beta}$ the scalar is $\psi = \tfrac{1}{2}\sqrt{1 - \gamma^2}\ln((r - a)/(r - b))$, Banerjee and Sen's (2.17), and the Coulomb field is $F_{tr} = qE(a - b)^2/((a - bE)^2(r - a)(r - b))$, which is $q$ over the square of the areal radius since $g_{tt}g_{rr} = -1$, their (2.18).
In the units of the metric file $\varphi = \sqrt{(1 - \gamma^2)c^4/16\pi G}\,\ln((r - a)/(r - b))$, as on Fisher's page.
`print_charts.penney_scalar_check` holds every component of the Ricci tensor to that source at three points of the parameters and the radius, to forty digits, never symbolically.
At $\gamma = 1$, $E = (r - a)/(r - b)$ and the chart is the published chart of `rn_metric` with $r_s = a + b$ and $r_q^2 = ab$; at $b = 0$ it is the published spherical chart of `fisher_jnw` with its $b$ equal to $a$.
The check compares both, slot by slot, at the same points.
For $\gamma < 1$ the sphere $r = a$ has zero area, since $E \to 0$ there while $a - bE \to a$, and the Ricci scalar $R = (a - b)^4(1 - \gamma^2)E/(2(a - bE)^2(r - a)^2(r - b)^2)$ diverges as $(r - a)^{\gamma - 2}$.

## Step 4: Bronnikov's harmonic chart

Bronnikov's (36), (37) and (40) of 1973 with $F = 1$, and the same solution in Bronnikov, Bolokhov, Skvortsova, Ibadov and Shaymanova's (19) and (24) to (27), take the harmonic coordinate $u$ with $\alpha = 2\beta + \gamma$, in which $\psi = Cu$.
On the branch $k > h > 0$, $u_1 > 0$, their $s(k, u) = \sinh(ku)/k$ and $q\,s(h, u + u_1) = \sinh(h(u + u_1))/\sinh(hu_1)$ by their (29).
With $K = k\coth(ku)$, $H = h\coth(h(u + u_1))$ and $q = h/\sinh(hu_1)$, $1/s^2(k, u) = K^2 - k^2$ and the metric is

$$ds^2 = -\frac{H^2 - h^2}{q^2}c^2dt^2 + \frac{q^2(K^2 - k^2)}{H^2 - h^2}\left((K^2 - k^2)du^2 + d\Omega^2\right).$$

$K$ and $H$ are held with the slopes $k^2 - K^2$ and $h^2 - H^2$, each written in itself, so no relation between names is left out.
At $u = 0$, $H = h\coth(hu_1)$, the mass, and $q$ is the charge, so $m^2 = q^2 + h^2$, Bronnikov's (30).
The dictionary to the spherical chart is $r = b + 2k/(1 - e^{-2ku})$, $a - b = 2k$, $\gamma = h/k$ and $a/b = e^{2hu_1}$, so $e^{-2ku} = (r - a)/(r - b)$ and $\psi = -\sqrt{k^2 - h^2}\,u$.
`print_charts.penney_harmonic_check` holds each declared slope to its definition, $u$ to being harmonic, and the chart to the spherical one pulled back by that dictionary, slot by slot, at three points.
The factors $(e^{hu_1} \pm 1)$ in which the reader writes $q$ are written back as powers of $q$ by the chart's `pretty`.

## Step 5: what the drawings use

Every drawing shows one spacetime, $b = a/4$ and $\gamma = 1/2$, so $q = a/4$ and $GM/c^2 = 5a/16$: the spherical chart at $a = 1$, and the harmonic chart at $k = 1$, $h = 1/2$ and $u_1 = \ln 4$, which is $a = 8/3$ in units of $k$, with $q = 2k/3$ and $GM/c^2 = 5k/6$.
A ray moving in from $r = 2a$ reaches the singularity after $2.8\,a/c$; with $P = \sqrt{(r - a)(r - b)}$ the tortoise coordinate at $\gamma = 1/2$ is $r_* = ((a^2 + b^2)P - 2ab(r - a))/(a - b)^2 + (a + b)\ln\left((\sqrt{r - a} + \sqrt{r - b})/\sqrt{a - b}\right)$, which vanishes at $r = a$.

## Step 6: the embedding diagram

The equatorial plane at one moment is read in the spherical chart, where the circle of the sphere $r$ has the radius $\rho = (a - bE)\sqrt{(r - a)(r - b)/E}/(a - b)$, which shrinks to a point at $r = a$ as $(r - a)^{1/4}$.
At $b = a/4$ and $\gamma = 1/2$, $g_{rr} - (d\rho/dr)^2$ vanishes where $E^2 - 15E + 4 = 0$, so $E_e = (15 - \sqrt{209})/2$ and $r_e = (209 + 9\sqrt{209})a/320 = 1.0597\,a$, the larger root of $800r^2 - 1045ar + 209a^2$.
Beyond $r_e$ the surface is drawn in flat space and inside it in three dimensional Minkowski space, as on Fisher's page, and the two pieces lie level at $r_e$.
Toward the singular point the surface runs along a light cone, where a chord's proper length is the small difference of its extents in $\rho$ and $Z$, so the profile stops at $r - a = 10^{-7}a$, a circle of radius $0.019\,a$.
`embedding.penney_scalar` checks the heights against quadratures of $g_{rr} - (d\rho/dr)^2$ written by hand and the harmonic chart against the same circles and heights.

## Step 7: the conformal diagram

On the plane of $t$ and $r$ the metric is $F(-c^2dt^2 + dr_*^2)$ with $F = -g_{tt}$, and $r_*$ of Step 5 is finite at $r = a$ for every $\gamma < 1$, since $F$ vanishes there as $(r - a)^\gamma$.
So $p, q = \arctan((ct \mp r_*)/\ell)$ with $\ell = 4a$ bring the spacetime into Minkowski's triangle, with $r = a$ on its edge $X = 0$, a timelike singularity with no horizon, as for Fisher's solution.
The harmonic view draws the same spacetime, $r/a = 1/4 + (3/4)/(1 - e^{-2ku})$ with times and lengths $3/8$ of their count in $a$.
`conformal.penney_scalar` samples the harmonic chart only to $u = e^{1.5}/k$: beyond it $K^2 - k^2$ is a difference of numbers that agree to within $e^{-2ku}$, and the published $g^{uu}$, which holds its square, loses digits past the check's part in $10^8$ by $u = 6/k$.

## Step 8: the curvature figure

`curvature.penney_scalar` draws the published $K$ and $R$ of the spherical chart against $r$ from $a$ to $5a$, beside Reissner and Nordström's $K$ at $r_s = a + b$ and $r_q^2 = ab$, the case $\gamma = 1$, whose value on its horizon $r = a$ is $29/4a^4$.
Near $r = a$, $K$ diverges as $(r - a)^{2\gamma - 4}$ and $R$ as $(r - a)^{\gamma - 2}$, with corrections of order $(r - a)^{1/2}$ from $a - bE$, so the powers are measured at $r - a = 2^{-34}a$, $2^{-43}a$ and $2^{-52}a$.
Far out $K$ falls to $\gamma^2$ of the black hole's, since the mass is $\gamma(a + b)/2$, and $R$ falls off as $(1 - \gamma^2)(a - b)^2/2r^4$.
