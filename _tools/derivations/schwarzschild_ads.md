# Schwarzschild-anti-de Sitter

Hawking and Page's black hole is Kottler's metric with $\Lambda = -3/L^2$,

$$ds^2 = -f\,c^2dt^2 + \frac{dr^2}{f} + r^2\,d\Omega^2,\qquad f = 1 - \frac{r_s}{r} + \frac{r^2}{L^2}.$$

Its three charts are written by `_tools/derivations/print_charts.py --metric schwarzschild_ads`, and `verify_metrics.py --system schwarzschild_ads/<chart>` checks each in seconds.

## Step 1. The parameters and the charts

The parameters are Schwarzschild's $r_s = 2GM/c^2$ and anti-de Sitter's radius $L$, the two the collection already uses for the parents, so each chart reduces to Schwarzschild's as $L \to \infty$ and to anti-de Sitter's static chart at $r_s = 0$.
Hawking and Page write $V = 1 - 2M/(m_p^2r) + r^2/b^2$ with $b = (-3/\Lambda)^{1/2}$, which is the same function with $b$ for $L$.
The charts are Schwarzschild's and Kottler's: the static chart, which covers each region with $f \neq 0$ on its own, the advanced chart $v = ct + r_*$ through the black hole, and the retarded chart $u = ct - r_*$ through the white hole.
A Kruskal chart gives $r$ only implicitly, so none is printed; the conformal diagram draws its own map through the horizon.

## Step 2. The horizon

$L^2rf = r^3 + L^2r - L^2r_s$ is a cubic whose derivative $3r^2 + L^2$ is positive, so it has one real root, $r_h$, which is positive, and a complex pair with negative real part $-r_h/2$.
The mass follows from the horizon, $r_s = r_h(1 + r_h^2/L^2)$, and $f'(r_h) = (L^2 + 3r_h^2)/(L^2r_h)$, so the surface gravity is $\kappa = f'(r_h)/2$ and Hawking's temperature is

$$T = \frac{\hbar c}{4\pi k_B}\,\frac{L^2 + 3r_h^2}{L^2r_h},$$

least, $T_0 = \sqrt{3}\,\hbar c/(2\pi k_BL)$, at $r_h = L/\sqrt{3}$, and equal to Hawking and Page's $T_1 = \hbar c/(\pi k_BL)$ at $r_h = L$, their Eqs. (2.7), (2.8), and (3.10) with $\hbar$, $c$, and $k_B$ restored.

Every diagram is drawn at $r_s = 2L$, where the cubic factors as $(r - L)(r^2 + Lr + 2L^2)$: the horizon is $r_h = L$ exactly, the black hole at the Hawking-Page temperature, with $\kappa = 2/L$ and complex roots $L(-1 \pm i\sqrt{7})/2$.

## Step 3. The tortoise coordinate

$1/f = L^2r/\left((r - r_h)(r^2 + r_hr + r_h^2 + L^2)\right)$ falls as $L^2/r^2$, so it is a sum of simple poles alone whose residues $A_i = 1/f'(r_i)$ sum to zero, and

$$r_* = \mathrm{Re}\sum_i A_i\ln\left(1 - \frac{r}{r_i}\right)$$

vanishes at $r = 0$, which is how the charts' conventions fix its constant, and tends to a finite $R$ as $r \to \infty$: the conformal boundary is timelike, as anti-de Sitter space's is.
At $r_s = 2L$, in units of $L$,

$$r_* = \frac{1}{4}\ln|1 - r| - \frac{1}{8}\ln\frac{r^2 + r + 2}{2} + \frac{5}{4\sqrt{7}}\left(\arctan\frac{2r + 1}{\sqrt{7}} - \arctan\frac{1}{\sqrt{7}}\right),\qquad R = 0.65805.$$

## Step 4. Printing

`chart_printer.Printer` orders the terms of each sum by rising powers of $r$ and then of $r_s$, with `flip = False`, so every value is printed around $L^2r - L^2r_s + r^3$, in the order of $f$ itself.
The metric and its inverse are written as the line element writes $f$, the Ricci scalar as $-12/L^2$, and the Kretschmann scalar as

$$K = \frac{12r_s^2}{r^6} + \frac{24}{L^4},$$

Schwarzschild's plus anti-de Sitter's, since the square of the Weyl tensor is Schwarzschild's and the Ricci tensor is $-3g_{\mu\nu}/L^2$; each is checked against sympy before it is written.

## Step 5. The conformal diagram

`AdSHoleTower` in `conformal.py` writes a `Tower`'s cells in $G(u) = \arctan e^{-\kappa u}$ with $u, v = ct \mp (r_* - R)$, so that Kruskal's $UV$ is $-1$ on the conformal boundary, which $p = \arctan U$ and $q = \arctan V$ put on the vertical lines $X = \pm\pi/2$, as for the BTZ hole.
At $r_s = 2L$,

$$UV = \frac{1 - r}{\sqrt{r^2 + r + 2}}\exp\left(\frac{5}{\sqrt{7}}\left(\arctan\frac{2r + 1}{\sqrt{7}} - \frac{\pi}{2}\right)\right),$$

which is $0$ on the horizon and $s = e^{-2\kappa R} = 0.0719$ at $r = 0$.
The BTZ hole has $s = 1$, which puts $r = 0$ on the straight lines $T = \pm\pi/2$ and makes its diagram a square.
Here $s < 1$, so $r = 0$ lies on $\tan p\tan q = s$, which runs from one boundary to the other and reaches only $T = 2\arctan\sqrt{s} = 0.524$ on the axis: the singularities bow inward, which Fidkowski, Hubeny, Kleban, and Shenker showed for every dimension above three.
The ray that leaves the right boundary at $t = 0$ has $V = 1$ and meets $r = 0$ at $U = s$, $X = \pi/4 - \arctan s = 0.714$, on its own side of the axis, which the script checks.

## Step 6. The embedding

On the equator of $t = 0$ the metric is $dr^2/f + r^2d\phi^2$.
In flat space $dz/dr = \sqrt{1/f - 1}$, real where $f \le 1$, from the throat $r_h$ to $r = (r_sL^2)^{1/3}$, where the surface lies level.
Beyond it $g_{rr} < 1$, as everywhere on the static slice of anti-de Sitter space, and the slice is drawn on in three dimensional Minkowski space at $dZ/dr = \sqrt{1 - 1/f}$, as the BTZ hole's is.
At $r_s = 2L$ the join is at $2^{1/3}L = 1.26\,L$, and the drawing runs to $r = 3L$ on both exteriors.
