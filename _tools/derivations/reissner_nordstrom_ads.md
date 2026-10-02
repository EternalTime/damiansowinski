# Reissner-Nordström-anti-de Sitter

The charged black hole in anti-de Sitter space is the static, spherically symmetric solution of the Einstein-Maxwell equations with $\Lambda = -3/L^2$,

$$ds^2 = -f\,c^2dt^2 + \frac{dr^2}{f} + r^2\,d\Omega^2,\qquad f = 1 - \frac{r_s}{r} + \frac{r_q^2}{r^2} + \frac{r^2}{L^2}.$$

Its three charts are written by `_tools/derivations/print_charts.py --metric reissner_nordstrom_ads`, and `verify_metrics.py --system reissner_nordstrom_ads/<chart>` checks each in seconds.

## Step 1. The parameters and the sources of the charts

The parameters are Reissner and Nordström's $r_s$ and $r_q$ and anti-de Sitter's radius $L$, the three the collection already uses for the parents.
Chamblin, Emparan, Johnson and Myers write $V(r) = 1 - m/r^{n-2} + q^2/r^{2n-4} + r^2/l^2$ in $n + 1$ dimensions, their (8) and (9), Phys. Rev. D 60, 064018, which at $n = 3$ is $f$ with $m = r_s$, $q = r_q$ and $l = L$.

- `static`: their (8) at $n = 3$. It is Schwarzschild-anti-de Sitter's static chart at $r_q = 0$, which `reissner_nordstrom_ads_check` holds it to, and Reissner and Nordström's as $L \to \infty$.
- `eddington_finkelstein_ingoing` and `eddington_finkelstein_outgoing`: $v, u = ct \pm r_*$ with $dr_*/dr = 1/f$ and $r_*(0) = 0$, as the Schwarzschild-anti-de Sitter and Reissner-Nordström-de Sitter pages build them; each is checked to be the static chart pulled back.

`reissner_nordstrom_ads_check` also holds every chart to the Einstein-Maxwell equations, $G^\mu{}_\nu = (3/L^2)\delta^\mu{}_\nu + (r_q^2/r^4)\,\mathrm{diag}(-1, -1, 1, 1)$, the stress of the Coulomb field.
A Kruskal chart gives $r$ only implicitly, so none is printed; the conformal diagram draws its own map through both horizons.

## Step 2. The horizons and the temperature

$L^2r^2f = r^4 + L^2r^2 - L^2r_sr + L^2r_q^2$ has at most two positive roots, $r_- \le r_+$, and a complex pair.
With $r_s = r_+(1 + r_q^2/r_+^2 + r_+^2/L^2)$, $f'(r_+) = 1/r_+ - r_q^2/r_+^3 + 3r_+/L^2$, so

$$T = \frac{\hbar c}{4\pi k_B}\left(\frac{1}{r_+} - \frac{r_q^2}{r_+^3} + \frac{3r_+}{L^2}\right),$$

Chamblin and his coauthors' (23) at $n = 3$ with $\hbar$, $c$ and $k_B$ restored.
Their critical charge, (28) at $n = 3$, is $r_q = L/6$ at $r_+ = L/\sqrt{6}$, which is Kubizňák and Mann's (3.17) with $P = 3/(8\pi L^2)$.
At $r_s = 2r_q$, $f = (1 - r_q/r)^2 + r^2/L^2 > 0$: Romans's supersymmetric solution, their (16), has no horizon.

## Step 3. What the diagrams draw

Every diagram is drawn at $L = 1$, $r_s = 27L/8$ and $r_q^2 = 11L^2/8$, where

$$r^2f = (r - 1)\left(r - \tfrac{1}{2}\right)\left(r^2 + \tfrac{3}{2}r + \tfrac{11}{4}\right):$$

with roots $a$, $b$ and a quadratic $r^2 + (a + b)r + c$, the quartic needs $c = L^2 + a^2 + ab + b^2$, $r_s = (a + b)(c - ab)/L^2$ and $r_q^2 = abc/L^2$, and $a = L$, $b = L/2$ puts the event horizon where Schwarzschild-anti-de Sitter is drawn with its own.
Then $f'(1) = 21/8$ and $f'(1/2) = -15/2$, so $\kappa_+ = 21/(16L)$ and $\kappa_- = 15/(4L)$, the complex roots are $(-3 \pm i\sqrt{35})/4$ with residues $-13/105 \mp i\sqrt{35}/35$, and

$$r_* = \mathrm{Re}\sum_i \frac{\ln(1 - r/r_i)}{f'(r_i)} = \frac{8}{21}\ln|1 - r| - \frac{2}{15}\ln|1 - 2r| + \dots$$

vanishes at $r = 0$ and tends to $R = 0.4052\,L$ as $r \to \infty$; `slices.rnads_rstar` sums it over the four roots.

The Eddington-Finkelstein rows of `null_rays.py` take the time functions $v - h(r)$ and $u + h(r)$ with $h' = r^2/((r^2 + 1)(r^2 + 2))$, $h = \sqrt{2}\arctan(r/\sqrt{2}) - \arctan r$, since $v - r$ is spacelike where $f > 2$, near the singularity and far outside, and $0 < h' < 2/f$ wherever $f > 0$.
Their planes are drawn with the box's top below the height at which the moment $t = 0$ inside $r_-$ runs off along the inner horizon, where its slope is $2/(15(1/2 - r))$ and four decimals no longer place it.

## Step 4. The conformal diagram

`CarterAdSTower` in `conformal.py`, built for Kerr-anti-de Sitter's axis, takes any $f$ with two real roots and a complex pair, and `HyperbolicHoleDrawing` draws its cells with each exterior ending on the timelike curve $(-G(t - R), G(-t - R))$ and each region inside $r_-$ on the singularity, the straight line $X = \pm\pi/2$ since $r_*(0) = 0$.
That is the tower Brecher, He and Rozali draw for the hole in five dimensions, their figure 2, JHEP 04 (2005) 004, with the singularities straight and the boundaries curved where theirs has the boundaries straight; they show that both cannot be straight.

## Step 5. The embedding

On the equator of $t = 0$ the metric is $dr^2/f + r^2d\phi^2$, and $f = 1$ at the two real roots of $8r^4 - 27r + 11$, $r = 0.4163\,L$ and $1.3275\,L$.
Outside, the surface is Schwarzschild-anti-de Sitter's: in flat space from the throat $r_+$ to $1.3275\,L$, and in three dimensional Minkowski space at $dZ/dr = \sqrt{1 - 1/f}$ beyond.
Inside $r_-$ it runs through the inner bifurcation sphere in flat space in to $0.4163\,L$, and nearer the singularity, where the charge makes $f > 1$ again, in Minkowski space: $1/f \to r^2/r_q^2$, so the slope tends to $1$ and the sheet closes on $r = 0$ along a light cone, a proper distance $\int_0^{0.4163}dr/\sqrt{f} = 0.131\,L$ from the level circle.
A chord of that sheet keeps its proper length only while it is short against $r$, so the profile is taken through radii in the ratio $1.02$ from $L/100$ and written to twelve decimals.
Stuchlík and Hledík, Acta Physica Slovaca 52, 363 (2002), draw the embedding diagrams of the equatorial plane in flat space, their section 5.

## Step 6. The names in the History

Each name is as its paper's author line prints it, and Larry Romans, whose paper prints L. J. Romans, as INSPIRE's author record 991345 gives him, Romans, Larry James.
