# Bardeen

Bardeen's regular black hole, in the notation of the collection, is

$$ds^2 = -f\,c^2dt^2 + \frac{dr^2}{f} + r^2\,d\Omega^2,\qquad f = 1 - \frac{r_s\,r^2}{\left(r^2 + g^2\right)^{3/2}}.$$

Its three charts are written by `_tools/derivations/print_charts.py --metric bardeen`, and `verify_metrics.py --system bardeen/<chart>` checks each in seconds.

## Step 1. The parameters and the charts

The parameters are Schwarzschild's $r_s = 2GM/c^2$ and the length $g$, so each chart reduces to Schwarzschild's at $g = 0$.
Ayón-Beato and García (Phys. Lett. B 493, 149, 2000) write $f = 1 - 2mr^2/(r^2 + g^2)^{3/2}$ in units with $G = c = 1$, which is the same function with $2m$ for $r_s$, and Borde writes $e$ for $g$.
The charts are the static chart, which covers each region with $f \neq 0$ on its own, and the advanced and retarded Eddington-Finkelstein charts $v = ct + r_*$ and $u = ct - r_*$, the two Ayón-Beato and García name for carrying the metric through both horizons to $r = 0$.
A Kruskal chart gives $r$ only implicitly, so none is printed; the conformal diagram draws its own map through the horizons.

## Step 2. The horizons

$f' = r_s\,r\left(r^2 - 2g^2\right)/\left(r^2 + g^2\right)^{5/2}$, so $f$ has one minimum, at $r = \sqrt{2}\,g$, where it is $1 - 2r_s/(3\sqrt{3}\,g)$.
For $g < 2r_s/3\sqrt{3}$ it has two simple zeros $r_- < \sqrt{2}\,g < r_+$, at equality one double zero at $r = \sqrt{2}\,g = 0.544\,r_s$, and beyond it none: Ayón-Beato and García's $g^2 \le (16/27)m^2$.
With $Q = \sqrt{r^2 + g^2}$ the zeros are the roots of $Q^3 - r_sQ^2 + r_sg^2 = 0$ with $Q > g$.
Every diagram is drawn at $g = r_s/3$, where $r_+ = 0.775419\,r_s$ and $r_- = 0.300963\,r_s$, with surface gravities $\kappa_+ = 0.3431/r_s$ and $\kappa_- = 1.0844/r_s$, a ratio of $3.16$, close to the $3.2$ of the Reissner-Nordström tower.

## Step 3. The centre and the curvature

Near $r = 0$, $f = 1 - r_sr^2/g^3 + O(r^4)$, de Sitter's static metric with $\Lambda = 3r_s/g^3$.
The scalars are

$$R = \frac{3r_s\,g^2\left(4g^2 - r^2\right)}{\left(r^2 + g^2\right)^{7/2}},\qquad K = \frac{3r_s^2\left(4r^8 - 12g^2r^6 + 47g^4r^4 - 4g^6r^2 + 8g^8\right)}{\left(r^2 + g^2\right)^7},$$

Ayón-Beato and García's (2) and (4) with $2m = r_s$, finite everywhere, with $K = 24r_s^2/g^6$ at the centre, de Sitter's $8\Lambda^2/3$.

## Step 4. The source

`bardeen_source` checks every chart against Ayón-Beato and García's field equations, $G_\mu{}^\nu = 2\left(L_FF_{\mu\lambda}F^{\nu\lambda} - \delta_\mu^\nu L\right)$, for the monopole $F = g\sin\theta\,d\theta\wedge d\phi$ of the nonlinear electrodynamics

$$L(F) = \frac{3}{2sg^2}\left(\frac{\sqrt{2g^2F}}{1 + \sqrt{2g^2F}}\right)^{5/2},\qquad s = \frac{g}{r_s},$$

whose invariant is $F = g^2/2r^4$ on the monopole.
There $L = 3r_sg^2/\left(2(r^2 + g^2)^{5/2}\right)$, and the equations read $G^t{}_t = G^r{}_r = -2L$ and $G^\theta{}_\theta = G^\phi{}_\phi = 2(2FL_F - L)$.
The energy density $-G^t{}_t/8\pi$ is positive, so the weak energy condition holds, and since the radial pressure is minus the density the strong energy condition is $G^\theta{}_\theta \ge 0$.

$$G^\theta{}_\theta = \frac{3r_s\,g^2\left(3r^2 - 2g^2\right)}{2\left(r^2 + g^2\right)^{7/2}}$$

changes sign at $r = \sqrt{2/3}\,g$, inside which the strong energy condition fails, which the script checks.

## Step 5. Printing

Every value is a rational function of $r$, $g^2$, $r_s$ and the one radical $Q$.
sympy's `factor` rationalises a denominator that holds $Q^3 - r_sr^2$, and prints the sextic $(r^2 + g^2)^3 - r_s^2r^4$ in its place.
`bardeen_pretty` writes the value in $r$, $Q$ and $r_s$ with $g^2 = Q^2 - r^2$, among which no relation is left, factors it there, multiplies out any two factors that $Q \to -Q$ exchanges, since $r^2 - Q^2 = -g^2$, and writes each factor back: one even in $Q$ as a polynomial in $r$ and $g$, and any other with each odd power of $Q$ as that power of $r^2 + g^2$.
So the zero of $f$ stands in a denominator as $\left(r^2 + g^2\right)^{3/2} - r_sr^2$, positive outside $r_+$, and $\Gamma^t{}_{tr} = f'/2f$ reads as it is derived.
The printer orders each sum with the radical first, and the metric and its inverse are written as the line element writes $f$.

## Step 6. The tortoise coordinate

$1/f$ is no rational function, so it has no partial fractions in closed form.
It has a simple pole at each horizon, of residue $A_\pm = 1/f'(r_\pm) = \pm 1/2\kappa_\pm$, and

$$r_* = \int_0^r\left(\frac{1}{f} - \sum_\pm\frac{A_\pm}{s - r_\pm}\right)ds + \sum_\pm A_\pm\ln\left|1 - \frac{r}{r_\pm}\right|,$$

which vanishes at $r = 0$, as the charts' conventions fix it.
The integrand is smooth, and `slices.bardeen_rstar` sums it by Gauss and Legendre's rule of twelve points, on panels a tenth of $r_s$ wide out to $8\,r_s$ and half a unit of $\ln r$ wide beyond, working within a twentieth of $r_s$ of a horizon in forty digits, since $1/f$ and its pole cancel there to the last digits of a float.
`slices.py` checks its slope against the published $g_{rr}$, and `null_rays.py --verify` checks every traced ray against $ct \pm r_*$.

## Step 7. The diagrams

The spacetime diagrams are the plane of the time and $r$ in each chart, the future taken from the ingoing rays in the static and ingoing charts and from the outgoing rays in the outgoing one.
The conformal diagram is `BardeenTower`, a `Tower` whose roots, residues and $r_*$ are numerical: Reissner-Nordström's cells, with $r_*(0) = 0$ putting $r = 0$ on the vertical lines $X = \pm\pi/2$, drawn as a regular centre, where the Kretschmann scalar is checked to settle at $24r_s^2/g^6$.
The ingoing chart covers the cells I, II and III$'$, the last read with the static time $t' = r_* - v$, so that $q = \arctan e^{\kappa_+v}$ is one function of $v$ across both horizons; the outgoing chart is its image under $(p, q) \to (\pi - q, \pi - p)$ with $v = -u$, and covers III$'$, the white hole above it and the exterior above that.
The embedding diagram is the equator of $t = 0$: outside $r_+$ the two sheets through the outer bifurcation sphere, and inside $r_-$ two caps that meet on the inner bifurcation sphere, the widest circle.
There $g_{rr} - 1 = (1 - f)/f$ is positive down to the centre, where it vanishes as $r_sr^2/g^3$, so each cap closes smoothly and the sphere of radius $\sqrt{g^3/r_s} = 0.192\,r_s$ osculates it: the slice is a closed surface, Borde's change of topology (Phys. Rev. D 55, 7615, 1997).
`Slice.horizons` finds the zeros of a $1/g_{xx}$ with one square root as the real roots of its resultant with $q^2 - S$, kept where the numerator itself vanishes, and next to such a horizon the metric is evaluated at the root plus $u$ in sixty digits.
