# Reissner-Nordström-de Sitter

The charged black hole in de Sitter space is Reissner and Nordström's metric with Kottler's term added,

$$ds^2 = -f\,c^2dt^2 + \frac{dr^2}{f} + r^2\,d\Omega^2,\qquad f = 1 - \frac{r_s}{r} + \frac{r_q^2}{r^2} - \frac{\Lambda r^2}{3}.$$

Its four charts are written by `_tools/derivations/print_charts.py --metric reissner_nordstrom_de_sitter`, and `verify_metrics.py --system reissner_nordstrom_de_sitter/<chart>` checks each in seconds.

## Step 1. The charts

The static chart and the two Eddington-Finkelstein charts are the ones Schwarzschild-de Sitter carries, with $r_q$ added.
The advanced chart, $v = ct + r_*$, covers a region inside $r_-$, the black hole, the static region and the contracting region beyond $r_c$, and the retarded chart, $u = ct - r_*$, covers a region inside $r_-$, the white hole, the static region and the expanding region.
The fourth chart is the cosmological chart of Kastor and Traschen (1993) for the lukewarm hole $r_q = r_s/2$, in the form of Brill and Hayward's equation (24),

$$ds^2 = -\frac{c^2d\tau^2}{W^2} + W^2\left(d\rho^2 + \rho^2d\Omega^2\right),\qquad W = H\tau + \frac{r_s}{2\rho},\qquad \Lambda = \frac{3H^2}{c^2},$$

whose coefficients are rational and which covers, for $H > 0$, one region inside $r_-$, the white hole, the static region and the expanding region, with $\tau < 0$ inside $r = r_s/2$.
Kastor and Traschen's own time is $t$ with $H\tau = e^{Ht}$, which covers only $\tau > 0$.
Brill and Hayward's double null form gives $\rho$ only implicitly, so it is not printed.

## Step 2. The map from the static chart

With $r = H\tau\rho + r_s/2$ the cosmological line element is

$$ds^2 = -f\,c^2dt_c^2 - \frac{2Hr^2}{c\,(r - r_s/2)}\,c\,dt_c\,dr + \frac{r^2}{(r - r_s/2)^2}\,dr^2,\qquad c\,dt_c = \frac{c\,d\tau}{H\tau},$$

with $f = (1 - r_s/2r)^2 - H^2r^2/c^2$, and the static time is

$$cT = \frac{c}{H}\ln|H\tau| + F(r),\qquad \frac{dF}{dr} = \frac{Hr^2}{c\,(r - r_s/2)\,f}.$$

`reissner_nordstrom_de_sitter_pullback` in `print_charts.py` checks the Jacobian of this map against the cosmological chart in every slot before the chart is written.
$F$ is a sum of logarithms over $r_s/2$ and the four roots of $f$; its coefficient at $r_s/2$ is $-c/H$, so $T$ is continuous through $\tau = 0$, and its coefficients at $r_+$ and $r_c$ equal those of $r_*$ while the one at $r_-$ is its opposite, so $u = cT - r_*$ is regular at $r_+$ and $r_c$ and $v = cT + r_*$ at $r_-$.

## Step 3. The lukewarm hole every diagram draws

At $r_q = r_s/2$ the function factors, $f = (1 - r_s/2r - Hr/c)(1 - r_s/2r + Hr/c)$, and the diagrams take $Hr_s/c = 3/8$, $\Lambda r_s^2 = 27/64$, where the roots are rational or quadratic surds:

$$r_c = 2\,r_s,\qquad r_+ = \tfrac{2}{3}\,r_s,\qquad r_- = \tfrac{2\sqrt{7} - 4}{3}\,r_s = 0.4305\,r_s,\qquad r_n = -\tfrac{2\sqrt{7} + 4}{3}\,r_s.$$

$f'(r_+) = -f'(r_c) = 3/(8r_s)$, so $\kappa_+ = \kappa_c = 3/(16\,r_s)$, Romans's equal temperatures, and $\kappa_- = 0.496/r_s$.
$f$ is greatest between $r_+$ and $r_c$ at the root $1.298\,r_s$ of $9r^4 - 32r_s^3r + 16r_s^4$, the static radius.
$1/f = -64r^2/\left(9\prod_i(r - r_i)\right)$ has no polynomial part, so $r_* = \sum_i \ln|1 - r/r_i|/f'(r_i)$ with $r_*(0) = 0$, and $r_*$ tends to a finite value as $r \to \infty$: future and past infinity are spacelike.

## Step 4. Printing

The static and Eddington-Finkelstein charts take Schwarzschild-de Sitter's options: `rising` orders every sum by rising powers of $\Lambda$, then $r_q$, then $r_s$, so each value is printed around $3r^2 - 3r_sr + 3r_q^2 - \Lambda r^4$ in the order of $f$ itself.
The Kretschmann scalar is written as Reissner and Nordström's plus de Sitter's,

$$K = \frac{12r_s^2}{r^6} - \frac{48r_sr_q^2}{r^7} + \frac{56r_q^4}{r^8} + \frac{8\Lambda^2}{3},$$

since the Weyl tensor does not see $\Lambda$ and the Ricci tensor is $\Lambda g_{\mu\nu}$ plus the traceless part of the Maxwell field.
The cosmological chart passes `rnds_cosmological_pretty`, which writes twice the areal radius, $2H\tau\rho + r_s$, as one named factor and any other sum as a polynomial in it where that is shorter, as $3H^2(2H\tau\rho + r_s)^4 + 4r_s^2c^2$ in the Einstein tensor.

## Step 5. The drawings

The static plane takes its future from the ingoing chart inside the static radius and from the outgoing one beyond it, `orient="split"`.
In the Eddington-Finkelstein planes $v - r$ and $u + r$ are not time functions inside $r_-$, where $f > 2$, so their rows take $v - h(r)$ and $u + h(r)$ with $h' = r^2/(r^2 + r_q^2)$, which lies between $0$ and $2/f$ wherever $f > 0$.
On the cosmological plane $\tau$ is a time everywhere, the horizons are the hyperbolas $H\tau\rho = r - r_s/2$ where $|\nabla r|^2 = f$ vanishes, and the row declares the singularity $2H\tau\rho + r_s = 0$ in `singular_zero`; on $\tau = 0$ every $\rho$ has the areal radius $r_s/2$, which is no throat, so the row sets `no_throat`.
The conformal diagram is drawn by `ChargedKottler` in `conformal.py`, Kottler's map with the black hole a whole cell and two regions inside $r_-$ above it; one function of each null coordinate serves every cell, so the drawing is smooth through $r_+$ and $r_c$ and continuous through $r_-$, as Reissner and Nordström's tower is.
The embedding diagram has three views: the static moment between $r_+$ and $r_c$, one period of a chain as Kottler's is; the static moment inside $r_-$ down to $r = 0.2495\,r_s$, the root of $9r^4 + 64r_s^3r - 16r_s^4$ where $g_{rr} = 1$; and a moment of the cosmological chart, on which $g_{rr} = r^2/(r - r_s/2)^2$ at every $\tau > 0$, the infinitely long throat of the extremal Reissner-Nordström black hole.
The static time of each region is fixed only up to a constant: `slices.rnds_static_t` takes $T = 0$ at $H\tau = 1$ on the static radius, and the moment embedded inside $r_-$ is the one through $H\tau = -1/2$ on $r = 0.35\,r_s$, which lies on the drawn cosmological plane.
