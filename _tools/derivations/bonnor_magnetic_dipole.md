# Bonnor's magnetic dipole

Bonnor's static solution of the Einstein-Maxwell equations for a mass with a magnetic dipole moment, in his own coordinates and his four polynomials, is

$$ds^2 = -\frac{P^2}{Y^2}c^2dt^2 + \frac{P^2Y^2}{Q^3}\left(\frac{dr^2}{Z} + d\theta^2\right) + \frac{Y^2Z\sin^2\theta}{P^2}d\phi^2,$$

$$P = r^2 - 2mr - b^2\cos^2\theta,\quad Q = (r - m)^2 - (m^2 + b^2)\cos^2\theta,\quad Y = r^2 - b^2\cos^2\theta,\quad Z = r^2 - 2mr - b^2,$$

with the potential $A_\phi = 2mbr\sin^2\theta/P$ in units with $G = c = 1$.
The chart is written by `_tools/derivations/print_charts.py --metric bonnor_magnetic_dipole`, which took 142 seconds on 2 October 2026, and `verify_metrics.py --system bonnor_magnetic_dipole/spheroidal` checks it in 124 seconds, 53 of them the dimensional pass.

## Step 1. The source of the chart

Bonnor's paper of 1966, Z. Phys. 190, 444, is two pages behind a paywall, so the line element is taken from the two papers that quote it and agree: Emparan's (1) to (3) of "Black diholes", in Kerr's $\Delta = Z$ and $\Sigma = Y$ with $\Delta + a^2\sin^2\theta = P$ and $\Delta + (M^2 + a^2)\sin^2\theta = Q$, and Kovář, Kopáček, Karas, and Kojima's (1) to (6), who write $P$, $Q$, $Y$, and $Z$.
Emparan's $M$ and $a$ are $m$ and $b$ here, and Kovář's $a$ and $b$ are $m$ and $b$.
The mass is $GM/c^2 = 2m$ and the dipole moment is $2mb$.
`bonnor_magnetic_dipole_check` holds the chart to $G_{\mu\nu} = 2(F_{\mu\alpha}F_\nu{}^\alpha - \tfrac{1}{4}g_{\mu\nu}F^2)$ with $F = dA$, to Maxwell's equations, and to a vanishing Ricci scalar before anything is written.

## Step 2. Without the dipole

At $b = 0$ the four polynomials are $P = Z = r(r - 2m)$, $Y = r^2$, and $Q = r^2 - 2mr + m^2\sin^2\theta$, and the metric is Zipoy and Voorhees's at $q = 1$, Darmois's solution: $g_{tt} = -(1 - 2m/r)^2$.
The same check reads the published spherical chart of `zipoy_voorhees` and compares every component.
So the dipole with its magnet switched off is a mass with a quadrupole and no horizon, which is the reason Schwarzschild's metric is no case of it.

## Step 3. How the values are printed

Every value is a rational function of $r$, $\cos\theta$, and $\sin\theta$ whose denominator is a product of powers of $P$, $Q$, $Y$, and $Z$.
`named_factors` writes each of those four as its letter, with $Y$ put back together from $(r - b\cos\theta)(r + b\cos\theta)$, and every even power of $\sin\theta$ is written in $\cos\theta$ first, so that the four polynomials survive factoring.
Each numerator is collected by its powers of $\cos\theta$ with the coefficient of each power factored, which brings the Kretschmann scalar from 4430 characters to 2552.
The metric and its inverse are written by hand as the line element writes them.

## Step 4. The two charts that are not published

Emparan and Teo write the dipole in Weyl's canonical coordinates, $\rho = \sqrt{Z}\sin\theta$ and $z = (r - m)\cos\theta$, where $r - m = (R_+ + R_-)/2$ and $\sqrt{m^2 + b^2}\cos\theta = (R_+ - R_-)/2$ with $R_\pm = \sqrt{\rho^2 + (z \pm \sqrt{m^2 + b^2})^2}$, the distances from the two black holes.
There $-g_{tt} = (P/Y)^2$ and $e^{2\gamma} = (P/Q)^4$, which is their (2.10).
With the two radicals as names the checker's `Geometry` did not build the metric and its inverse in eight minutes, as Zipoy and Voorhees's Weyl chart did not finish its Kretschmann scalar, so the convention of Bonnor's chart gives $\rho$ and $z$ and the chart is not published.

The prolate spheroidal chart $x = (r - m)/k$, $y = \cos\theta$ with $k = \sqrt{m^2 + b^2}$ is rational in $x$, $y$, $m$, and $k$, prints in two minutes, and was checked to be Bonnor's chart pulled back and to solve the Einstein-Maxwell equations.
The checker's dimensional pass, which expands every published value, took thirty minutes on it, against the two minutes the whole of Bonnor's chart takes, so it is not published either.

## Step 5. The three planes

The axis beyond a hole, the axis between the holes, and the equatorial plane are totally geodesic, the first two fixed by the rotations and the third by the reflection $\theta \to \pi - \theta$, so their null curves are null geodesics.
On the axis beyond a hole $P = Q = Z$ and the metric is $-(Z/Y)^2c^2dt^2 + (Y/Z)^2dr^2$, so $r$ is affine along every ray and $c\,dt/dr = \pm(Y/Z)^2$ has a double pole at $Z = 0$, the mark of a degenerate horizon.
On the axis between the holes, $Z = 0$, $P = b^2\sin^2\theta$ and $Q = (m^2 + b^2)\sin^2\theta$, so $c\,dt/d\theta = \pm Y^2/((m^2 + b^2)^{3/2}\sin^3\theta)$, and the proper length $\int b^2Y\,d\theta/((m^2 + b^2)^{3/2}\sin\theta)$ diverges at both ends as a logarithm.
On the equatorial plane $P = r(r - 2m)$, $Q = (r - m)^2$, and $Y = r^2$, so $c\,dt/dr = \pm r^4/((r - m)^3\sqrt{Z})$, integrable at $Z = 0$.

A circle about the axis between the holes has $d\sqrt{g_{\phi\phi}}/(\sqrt{g_{rr}}\,dr) \to (m^2 + b^2)^2/b^4$ at $Z = 0$, Emparan's excess of angle, using $r_+(r_+ - 2m) = b^2$ at $r_+ = m + \sqrt{m^2 + b^2}$; `_tools/test_bonnor_magnetic_dipole.py` holds the published components to it, to a regular axis beyond the holes, and to the mass $2m$.

## Step 6. The diagrams

Every drawing takes $m = 1$ and $b = 2\sqrt{2}\,m$, so that $\sqrt{m^2 + b^2} = 3m$, the axis between the holes is $r = 4m$, and the excess of angle is $81/64$.
At $b = \sqrt{3}\,m$, where the numbers are $2m$, $3m$, and $16/9$, the cone at the strut is too steep for the embedding's check across a turn of $0.02$, which a cone of $16/9$ of a turn misses at $2.1 \times 10^{-4}$ against the $2 \times 10^{-4}$ allowed, however short its chords are.
The spacetime diagrams draw the three planes, and `null_rays.py --verify` compares their rays with $ct \pm x_*$ for $r_* = r - 16/(9(r - 4)) - 4/(9(r + 2)) + (80/27)\ln(r - 4) + (28/27)\ln(r + 2)$ on the axis, a quadrature in the equatorial plane, and $\theta_* = (64/27)(-\cot\theta\csc\theta/2 + (5/2)\ln\tan(\theta/2) - \cos\theta)$ between the holes, to $10^{-9}$ of the drawing.
The Kretschmann scalar is finite at a hole along either stretch of axis, with two different values, and its printed form is a ratio of polynomials that both vanish there, so the rows of the axis between the holes do not take it in floating point.
The conformal diagrams take $p, q = \arctan((ct \mp x_*)/\ell)$ with $\ell = 4m$: the axis beyond a hole is a diamond with the horizon on its left, the axis between the holes a diamond with a horizon on every side, and the equatorial plane Minkowski's triangle with the strut on $X = 0$.
The embedding diagram is the equatorial plane at $t = 0$, whose circles have radius $r\sqrt{Z}/(r - 2m)$: a cone of $81/64$ of a turn at the strut, drawn in Minkowski space out to the level circle at $r = 4.47\,m$, and a sheet in flat space beyond.
