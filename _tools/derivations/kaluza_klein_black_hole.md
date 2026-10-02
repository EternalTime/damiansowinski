# The Kaluza-Klein black holes

Gravity in five dimensions with the fifth a circle has charged black holes that are vacuum solutions: seen from four dimensions the twist and the size of the circle are a Maxwell field and a scalar field.
This note records the parameters, the source of each chart, what each is checked against, and what the diagrams draw.
Equation numbers of Horowitz and Wiseman's review are those of arXiv:1107.5563.

## Step 1: the parameters

The collection keeps Horowitz and Wiseman's three lengths, with their $2m$ written $r_s$: the horizon is at $r = r_s$, and $q \ge r_s$ and $p \ge r_s$ carry the electric and the magnetic charge.
In units with $G = c = 1$ in four dimensions, their (2.40), (2.42) and (2.43), $M = (p + q)/4$, $Q^2 = q(q^2 - r_s^2)/4(p + q)$ and $P^2 = p(p^2 - r_s^2)/4(p + q)$.
A hole with no magnetic charge has $p = r_s$, so $M = (q + r_s)/4$ and $Q = \sqrt{q(q - r_s)}/2$, and $q = r_s\cosh^2\alpha$ for a boost of rapidity $\alpha$.
Gibbons and Wiltshire, and Rasheed after them, write the family with $M$, $Q$, $P$ and the scalar charge $\Sigma$ bound by $Q^2/(\Sigma + \sqrt{3}M) + P^2/(\Sigma - \sqrt{3}M) = 2\Sigma/3$; with $\Sigma = \sqrt{3}(q - p)/4$ the three lengths solve that cubic identically, which is why they are used.

## Step 2: the general dyon

`kaluza_klein_dyon` in `print_charts.py` is Horowitz and Wiseman's (2.31) to (2.37) at $a = 0$, Larsen's form of Rasheed's solution:

$$ds^2 = \frac{H_2}{H_1}\left(dy + 2A\right)^2 - \frac{r(r - r_s)}{H_2}c^2dt^2 + H_1\left(\frac{dr^2}{r(r - r_s)} + d\Omega^2\right),$$

with $H_1 = r^2 + (p - r_s)r + p(p - r_s)(q - r_s)/2(p + q)$, $H_2$ the same with $p$ and $q$ exchanged, and $2A = -2Q\left(r + (p - r_s)/2\right)H_2^{-1}c\,dt + 2P(1 - \cos\theta)\,d\phi$.
The potential along $\phi$ is theirs moved by a constant, so that it vanishes on $\theta = 0$ as Gross and Perry's does.
`kaluza_klein_vacuum` holds it to being Ricci flat at six random points in forty digits.

It is not published as a chart.
Printed through the checker's printer it runs to 1.2 megabytes, more than all six published charts together, and its values carry six square roots of the parameters; the three charts below are its cases with closed forms a reader can use.

## Step 3: the charts

The electric chart is the dyon at $p = r_s$, Horowitz and Wiseman's (2.11): Schwarzschild's black string in a frame moving along $y$,

$$ds^2 = -\left(1 - \frac{q}{r}\right)c^2dt^2 - \frac{2\sqrt{q(q - r_s)}}{r}\,c\,dt\,dy + \left(1 + \frac{q - r_s}{r}\right)dy^2 + \frac{dr^2}{1 - r_s/r} + r^2d\Omega^2.$$

`kaluza_klein_check` holds it to the dyon at $p = r_s$ and to the published static chart of `black_string` pulled back through $t' = t\cosh\alpha - y\sinh\alpha$, $y' = y\cosh\alpha - t\sinh\alpha$.
Its Kretschmann scalar is the string's, $12r_s^2/r^6$.
$g_{tt}$ vanishes at $r = q$, outside the horizon: between $r_s$ and $q$ no observer keeps $y$ fixed, an ergoregion of the moving string.

The ingoing Eddington-Finkelstein chart is the same boost of the string's ingoing chart, $v = ct + \sqrt{q/r_s}\,r_*$ and $w = y + \sqrt{(q - r_s)/r_s}\,r_*$ with Schwarzschild's $r_* = r + r_s\ln|r/r_s - 1|$, checked to be the electric chart pulled back.

The magnetic chart is the dyon at $q = r_s$, Horowitz and Wiseman's (2.28).
At $r_s = 0$ it is the published chart of `kaluza_klein_monopole` with $p = 4m$, which is checked, and the period $4\pi\sqrt{p(p - r_s)}$ of $y$ is the one that removes the Dirac string, $8\pi P$.

The dyonic chart is the dyon at $q = p$, where $H_1 = H_2 = (r + (p - r_s)/2)^2$, written in the areal radius $\rho = r + (p - r_s)/2$.
The circle has one length, $g_{yy} = 1$, and the metric orthogonal to it is the published metric of `rn_metric` with $r_s \to p$ and $r_q^2 \to (p^2 - r_s^2)/4$, horizons at $\rho_\pm = (p \pm r_s)/2$; both are checked.
Larsen notes that this case is the standard Reissner-Nordström solution, since with equal charges the scalar field may be held constant.

The two Einstein charts print the metric of four dimensions, $ds_5^2 = e^{-4\varphi/\sqrt{3}}(dy + 2A)^2 + e^{2\varphi/\sqrt{3}}g_{\mu\nu}dx^\mu dx^\nu$, Horowitz and Wiseman's (2.15) and (2.16), and its ingoing form on $v = ct + r_*$ with $dr_*/dr = \sqrt{r(r + q - r_s)}/(r - r_s)$.
The static one is held to being the reduction of the electric chart and to the field equations of four dimensions, their (2.4) to (2.6), at random points; the magnetic hole has the same metric with $p$ for $q$.
They are no vacuum, so `REGIONS` in `metric_tags.py` leaves them out of the tags.

## Step 4: the tortoise coordinate

With the circle divided out the plane of $t$ and $r$ of the electric hole has $c\,dt/dr = \pm\sqrt{1 + (q - r_s)/r}\,(1 - r_s/r)^{-1}$, and so have the magnetic hole's plane on $\theta = 0$ and the Einstein metric.
In units of $r_s$ with $a = q - r_s$ and $s = \sqrt{r(r + a)}$,

$$r_* = s + \left(1 + \frac{a}{2}\right)\ln\frac{2s + 2r + a}{a} - \sqrt{1 + a}\,\ln\frac{(a + 2)r + a + 2\sqrt{1 + a}\,s}{a\,|r - 1|},$$

which vanishes at $r = 0$; `slices.kkbh_rstar` is this, and sympy confirms its derivative.
Its one logarithm at the horizon has the coefficient $\sqrt{1 + a} = \sqrt{q/r_s}$, so the surface gravity is $\kappa = 1/(2\sqrt{q\,r_s})$: Schwarzschild's divided by $\cosh\alpha$.
`_tools/test_kaluza_klein_black_hole.py` takes the mass, the charge, the temperature, the entropy and the potential from the published charts and holds them to Smarr's relation $M = 2TS + \Phi Q$.

## Step 5: the diagrams

The drawings are at $q = 2\,r_s$, $p = 2\,r_s$ and, for equal charges, $p = 3\,r_s$, where $\rho_+ = 2\,r_s$ and $\rho_- = r_s$.

Spacetime diagrams. The plane of $t$ and $r$ at fixed $y$ has no null direction between $r_s$ and $q$, so the electric charts and the dyonic chart are drawn with the circle divided out, `quotient`, as the rotating holes are with $\phi$: the rays are those with no momentum along the circle, which is no charge in four dimensions.
In the ingoing chart those rays do not keep $v$, and $v - r$ is no time function there, so the plane is drawn against $v - \sqrt{q/r_s}\,r$, which is one at every $r$.
The magnetic chart is drawn on the half axis $\theta = 0$, where its potential vanishes and the plane is totally geodesic.

Conformal diagrams. `KaluzaKleinTower` is a `Tower` on the $r_*$ above: $U = -e^{-\kappa u}$, $V = e^{\kappa v}$, $UV = (1 - r/r_s)e^{2\kappa D(r)}$ with $D$ the smooth part of $r_*$, and $r_*(0) = 0$ puts the singularity on $UV = 1$, the straight lines of Schwarzschild's square.
The equal charges are Reissner and Nordström's tower.

Embedding diagrams. Three moments $t = 0$, each through the horizon into the other exterior: the equator in the Einstein metric, whose throat has the radius $(q\,r_s^3)^{1/4}$; the surface of $r$ and $y$ at $\theta = 0$ of the electric hole, where the circle of $y$, drawn at $L = 2\pi r_s$, has the radius $r_s\sqrt{1 + (q - r_s)/r}$ and is widest on the horizon; and the same surface of the magnetic hole, where the circle has the radius $2\sqrt{p(p - r_s)}\sqrt{r/(r + p - r_s)}$ and is narrowest on the horizon, the monopole's cigar cut open by a horizon.
The equal charges have no view of their own: their circle has one length, and their equator is Reissner and Nordström's.
The geometry is static, so there is no movie.
