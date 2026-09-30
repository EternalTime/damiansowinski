# Schwarzschild-de Sitter, Kottler's metric

Kottler's metric is Schwarzschild's with de Sitter's term added,

$$ds^2 = -f\,c^2dt^2 + \frac{dr^2}{f} + r^2\,d\Omega^2,\qquad f = 1 - \frac{r_s}{r} - \frac{\Lambda r^2}{3}.$$

Its three charts are written by `_tools/derivations/print_charts.py --metric schwarzschild_de_sitter`, and `verify_metrics.py --system schwarzschild_de_sitter/<chart>` checks each in seconds.

## Step 1. The charts

The entry carries the charts Schwarzschild's carries, since Kottler's regions are covered the same way.
The static chart covers each region with $f \neq 0$ on its own, and its $t$ is a time only in the static region $r_h < r < r_c$.
The advanced chart, $v = ct + r_*$, covers the black hole, the static region and the contracting region beyond $r_c$, and the retarded chart, $u = ct - r_*$, covers the white hole, the static region and the expanding region.
A Kruskal chart regular at one horizon is singular at the other, since the two surface gravities differ, and gives $r$ only implicitly, so none is printed; the conformal diagram draws its own map through both.

## Step 2. The horizons

$3rf = 3r - 3r_s - \Lambda r^3$, a cubic with roots summing to zero.
Its depressed form $r^3 - (3/\Lambda)r + 3r_s/\Lambda = 0$ has three real roots exactly when $\Lambda r_s^2 < 4/9$, which is Podolský's $9\Lambda m^2 < 1$ with $r_s = 2m$:

$$r_k = \frac{2}{\sqrt{\Lambda}}\cos\left(\frac{\psi}{3} - \frac{2\pi k}{3}\right),\qquad \cos\psi = -\frac{3}{2}r_s\sqrt{\Lambda},$$

$r_c$ at $k = 0$, $r_h$ at $k = 1$ and a negative root $r_n = -(r_h + r_c)$ at $k = 2$.
At $\Lambda r_s^2 = 1/5$, where every diagram is drawn, $r_h = 1.0851996\,r_s$, $r_c = 3.2146274\,r_s$ and $r_n = -4.2998270\,r_s$, with $f'(r_h) = 0.70445/r_s$ and $f'(r_c) = -0.33185/r_s$, so $\kappa_h = 0.35222/r_s$ and $\kappa_c = 0.16592/r_s$.
$f'(r) = r_s/r^2 - 2\Lambda r/3$ vanishes at $r^3 = 3r_s/2\Lambda$, $1.9574\,r_s$ here, where $f = 0.2337$ is greatest and $\Gamma^r{}_{tt} = ff'/2$ vanishes, so an observer at rest there is in free fall.

## Step 3. The tortoise coordinate

$1/f = -3r/\left(\Lambda(r - r_h)(r - r_c)(r - r_n)\right)$ has no polynomial part, so

$$r_* = \sum_i \frac{\ln|1 - r/r_i|}{f'(r_i)},$$

with $r_*(0) = 0$, which is how the charts' conventions fix its constant.
$\sum_i 1/f'(r_i) = 0$, the coefficient of $1/r$ in $1/f$ at large $r$, so $r_*$ tends to the finite $R = -\sum_i \ln|r_i|/f'(r_i)$ as $r \to \infty$: future and past infinity are spacelike.

## Step 4. Printing

`chart_printer.Printer` takes two options for this entry.
`rising` orders the terms of a sum by rising powers of $\Lambda$ and then of $r_s$, so every value is printed around $3r - 3r_s - \Lambda r^3$ in the order of $f$ itself, and reduces term by term to Schwarzschild's at $\Lambda = 0$ and to de Sitter's static chart at $r_s = 0$.
`flip = False` keeps that order where the printer would otherwise write $-(A - B)/D$ as $(B - A)/D$.
The metric and its inverse are written as the line element writes $f$, and the Kretschmann scalar as

$$K = \frac{12r_s^2}{r^6} + \frac{8\Lambda^2}{3},$$

Schwarzschild's plus de Sitter's, since the square of the Weyl tensor is Schwarzschild's $12r_s^2/r^6$ and the Ricci tensor is $\Lambda g_{\mu\nu}$; each is checked against sympy before it is written.

## Step 5. The drawings

The static chart's spacetime diagram takes its future from the ingoing chart inside the static radius and from the outgoing one outside it, `orient="split"` in `null_rays.py`, which reads the region inside $r_h$ as the black hole and the region beyond $r_c$ as the expanding universe.
The conformal diagram is the chain of Lake and Roeder, drawn by `Kottler` in `conformal.py` with $p, q = \pm\arctan W$ of $u$ and of $-v$, $W(x) = \exp(-ax - b\sqrt{x^2 + r_s^2})$ and $a, b = (\kappa_h \pm \kappa_c)/2$, which tends to $e^{-\kappa_hx}$ toward $r_h$ and to $e^{-\kappa_cx}$ toward $r_c$; its docstring lists the cells.
The embedding diagram is the static moment $t = 0$ from the throat $r_h$ to the widest circle $r_c$, where it runs through the cosmological bifurcation sphere into the next static region, turned over, up to the next throat, one period of the chain.
Its horizons are the exact roots of the cubic, which sympy returns in radicals of complex numbers whose reality it cannot decide, so `Slice.horizons` takes them as algebraic numbers, and each coefficient of the metric expanded about a root is reduced modulo the root's minimal polynomial before it is evaluated.
