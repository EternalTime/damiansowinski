# The wormhole of Einstein-Dirac-Maxwell theory

The three charts of `einstein_dirac_maxwell_wormhole` are written by `print_charts.py`, and this note records where each comes from and what the script checks.

## Step 1: the metric

Blázquez-Salcedo, Knoll and Radu, Phys. Rev. Lett. 126, 101102 (2021), arXiv:2010.07317, solve their Einstein-Dirac-Maxwell equations exactly with the gauge coupling $q = 0$, massless spinors and frequency $w = 0$.
Their line element, read off the TeX source, is

$$ds^2 = -\left(1 - \frac{M}{r}\right)^2dt^2 + \frac{dr^2}{\left(1 - \frac{r_0}{r}\right)\left(1 - \frac{Q_e^2}{r_0r}\right)} + r^2d\Omega^2, \qquad M = \frac{2Q_e^2r_0}{Q_e^2 + r_0^2},$$

with $V = (M/Q_e)\sqrt{(1 - r_0/r)(1 - Q_e^2/r_0r)}$, $r_0$ the throat and $Q_e < r_0$.
Their units are $G = c = 1$ with the field equations $G_{\mu\nu} = 2T_{\mu\nu}$ and $T^{(M)}_{\mu\nu} = F_{\mu\alpha}F_\nu{}^\alpha - \tfrac14F^2g_{\mu\nu}$, the Gaussian geometrized units in which Reissner and Nordström's $g_{tt}$ is $-(1 - 2M/r + Q^2/r^2)$.
The collection keeps $c$ and writes $r_0$, $Q_e$ and $M$ as lengths.

## Step 2: the same metric in 2003

Bronnikov and Kim, Phys. Rev. D 67, 064027 (2003), arXiv:gr-qc/0212112, Example 3, eq. (wh3), take $e^{2\gamma} = (1 - 2m/r)^2$ and find $g_{rr} = r^2/((r - r_0)(r - r_1))$ with $r_1 = mr_0/(r_0 - m)$.
With $2m = M$, $r_1 = Q_e^2/r_0$: it is the same metric, as Bolokhov, Bronnikov, Krasnikov and Skvortsova point out, Grav. Cosmol. 27, 401 (2021), point 4.
Their substitution $r = r_0 + x^2$ gives the chart `bronnikov_kim`, written in $u$ to keep $x$ for the compact chart:

$$g_{uu} = \frac{4r_0(r_0 + u^2)^2}{r_0^2 - Q_e^2 + r_0u^2}.$$

(Their printed line has $(r_0 + x^2)\,d\Omega^2$ where the square belongs, a misprint.)
The metric is symmetric under $r_0 \leftrightarrow r_1$, so $Q_e > r_0$ gives the same spacetime with the throat at $Q_e^2/r_0$; the published domain takes $0 \le Q_e < r_0$.

## Step 3: the compact chart

Blázquez-Salcedo and Knoll, Eur. Phys. J. C 80, 174 (2020), arXiv:1910.03565, write their uncharged wormhole in $\rho = \sqrt{N}$, $N = 1 - \mu/r$, extended to $\rho \in (-1, 1)$, and Konoplya and Zhidenko, Phys. Rev. Lett. 128, 091104 (2022), arXiv:2106.05034, use $x = \kappa\sqrt{1 - r_0/r}$ for this theory, "suggested" by the analytic solution.
With $r = r_0/(1 - x^2)$, $1 - r_0/r = x^2$ and $dr = 2r_0x\,dx/(1 - x^2)^2$, so

$$g_{xx} = \frac{4r_0^4}{(1 - x^2)^4\left(r_0^2 - Q_e^2(1 - x^2)\right)}, \qquad g_{tt} = -\left(1 - \frac{M(1 - x^2)}{r_0}\right)^2.$$

Every component of the charts `bronnikov_kim` and `compact` is a function of $u^2$ or $x^2$, so the two sides of the throat join smoothly.
The chart of Kanti, Kleihaus and Kunz with $r^2 + r_0^2$ for the areal radius, which Blázquez-Salcedo, Knoll and Radu use for their numerical solutions, was never written for this metric, and the proper radial distance has no closed inverse; neither is published.

## Step 4: what `einstein_dirac_maxwell_wormhole_check` holds

In each chart: the Ricci scalar vanishes (massless Dirac and Maxwell fields are both traceless); their $V$, continued through the throat with the sign of the coordinate, solves Maxwell's equations with no source; $G^t{}_t = -Q_e^2/r^4$, the electric field's energy density; $G^a{}_b - 2T^{(M)a}{}_b$, twice the Dirac stress, has no energy density and no trace, and its radial component is $-(r_0^2 - Q_e^2)^2/\left(r_0(Q_e^2 + r_0^2)r^2(r - M)\right)$, the tension that violates the null energy condition and vanishes at $Q_e = r_0$; and the Kretschmann scalar at the throat is $2(3Q_e^4 - 2Q_e^2r_0^2 + 3r_0^4)/r_0^8$, $8/r_0^4$ at $Q_e = r_0$ as on the extreme horizon.
The areal chart at $Q_e = r_0$ is the published `rn_metric` with $r_s = 2r_0$, $r_q = r_0$, and at $Q_e = 0$ it has $g_{tt} = -1$, Blázquez-Salcedo and Knoll's uncharged wormhole.
The other two charts are the areal chart pulled back along $r = r_0 + u^2$ and $r = r_0/(1 - x^2)$.
Bronnikov and Kim's effective energy density $\rho = mr_0^2/\left(r^4(r_0 - m)\right)$ is $Q_e^2/r^4$ in these parameters.

## Step 5: the areal chart's Kretschmann scalar

For a static observer $K = 4A^2 + 8B^2 + 8C^2 + 4D^2$ with $A = R^{tr}{}_{tr}$, $B = R^{t\theta}{}_{t\theta}$, $C = R^{r\theta}{}_{r\theta}$, $D = R^{\theta\phi}{}_{\theta\phi}$.
With $W = (r_0^2 + Q_e^2)r - 2r_0Q_e^2 = (r_0^2 + Q_e^2)(r - M)$:
$A = Q_e^2\left(4r_0r^2 - 5(r_0^2 + Q_e^2)r + 6r_0Q_e^2\right)/r^4W$, $B = -2Q_e^2(r - r_0)(r_0r - Q_e^2)/r^4W$, $C = -W/2r_0r^4$ and $D = \left((r_0^2 + Q_e^2)r - r_0Q_e^2\right)/r_0r^4$, which is the published form, checked against sympy.

## Step 6: the tortoise coordinate

$dr_*/dr = r^2/\left((r - M)\sqrt{(r - r_0)(r - b)}\right)$ with $b = Q_e^2/r_0$, and $b < M < r_0$ for $0 < Q_e < r_0$.
Writing $r^2/(r - M) = r + M + M^2/(r - M)$ and integrating from the throat,

$$r_* = \sqrt{(r - r_0)(r - b)} + \frac{r_0 + b + 2M}{2}\ln\frac{2r - r_0 - b + 2\sqrt{(r - r_0)(r - b)}}{r_0 - b} + \frac{2M^2}{\sqrt{(r_0 - M)(M - b)}}\arctan\sqrt{\frac{(M - b)(r - r_0)}{(r_0 - M)(r - b)}},$$

which `_edm_rstar` in `null_rays.py` evaluates at $r_0 = 1$, $Q_e = 1/2$ and which agrees with a quadrature in $u$ to $10^{-13}$.
Far out $g_{rr} - 1 \to (r_0 + b)/r$, so the slice of constant $t$ approaches Flamm's paraboloid with $r_s = r_0 + b$, while distant orbits feel the mass $M$: $r_0 + b > 2M$ unless $Q_e = r_0$.

## Step 7: the drawings

All are at $r_0 = 1$ and $Q_e = 1/2$, so $M = 2/5$ and $b = 1/4$.
The spacetime diagram of each chart is checked against $ct \mp r_*$, with $r_*$ times the sign of $u$ or $x$ through the throat; the conformal diagram is the full diamond by $p, q = \arctan((ct \mp x)/4r_0)$, the areal chart covering the right half; the embedding diagram is the equator at one moment on both sides, its height checked against the quadrature $dz/du = 2\sqrt{(r_0^2 + (r_0 + b)u^2)/(r_0 - b + u^2)}$ in all three charts.
