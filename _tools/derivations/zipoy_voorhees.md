# Zipoy and Voorhees's deformed mass

Zipoy and Voorhees's metric in Quevedo's spherical chart and in Voorhees's prolate spheroidal chart, $x = r/m - 1$, $y = \cos\theta$, $\delta = 1 + q$, is

$$ds^2 = -f^{1+q}c^2dt^2 + f^{-q}\left(h^{-q(2+q)}\left(\frac{dr^2}{f} + r^2d\theta^2\right) + r^2\sin^2\theta\,d\phi^2\right),\qquad f = 1 - \frac{2m}{r},\quad h = 1 + \frac{m^2\sin^2\theta}{r^2 - 2mr},$$

$$ds^2 = -f^{\delta}c^2dt^2 + \frac{m^2}{f^{\delta}}\left(\frac{x^2 - y^2}{h^{\delta^2}}\left(\frac{dx^2}{x^2 - 1} + \frac{dy^2}{1 - y^2}\right) + (x^2 - 1)(1 - y^2)\,d\phi^2\right),\qquad f = \frac{x - 1}{x + 1},\quad h = \frac{x^2 - y^2}{x^2 - 1}.$$

Both charts are written by `_tools/derivations/print_charts.py --metric zipoy_voorhees`, which took 105 seconds on 30 September 2026, and `verify_metrics.py --system zipoy_voorhees/spherical` and `--system zipoy_voorhees/prolate_spheroidal` check them in 12 and 7 seconds.

## Step 1. Weyl's two functions

In Weyl's canonical coordinates $\rho = m\sqrt{(x^2 - 1)(1 - y^2)}$ and $z = mxy$ the metric is $-e^{2\psi}c^2dt^2 + e^{-2\psi}\left(e^{2\gamma}(d\rho^2 + dz^2) + \rho^2d\phi^2\right)$ with $\psi = \tfrac{1}{2}\delta\ln f$ and $\gamma = -\tfrac{1}{2}\delta^2\ln h$.
$\psi$ is the Newtonian potential of a uniform rod of length $2m$ and mass per unit length $\delta/2$ on the axis, and multiplying Schwarzschild's $\psi$ by $\delta$ multiplies its $\gamma$ by $\delta^2$, since $\gamma$ is quadratic in the derivatives of $\psi$.
The mass is $GM/c^2 = \delta m$, and Quevedo's quadrupole is $M_2 = -m^3q(1 + q)(2 + q)/3$.
The check of each chart requires its Ricci tensor to vanish before anything is written, and `zipoy_voorhees_pullback` checks that the spherical chart is the pullback of the prolate spheroidal one.
Every power of $x - y$ is carried through the product $x^2 - y^2 = (r^2 - 2mr + m^2\sin^2\theta)/m^2$ in that check, since $x - y$ and $x + y$ each come to a factor linear in $\cos\theta$ and the spherical chart writes their product as one factor in $\sin^2\theta$.

## Step 2. Powers with a parameter in the exponent

Every value is a rational function times powers of $r$, $r - 2m$ and $r^2 - 2mr + m^2\sin^2\theta$ whose exponents hold $q$.
`norm` in `verify_metrics.py` reads $(1 - 2m/r)^{q}$ as $(r - 2m)^{q}r^{-q}$: a power of a sum or a product with a symbolic exponent is split over the irreducible factors of its base, as `_factored_root` already splits a radical, which is an identity where each factor is positive, as each is for $r > 2m$.
The whole part of an exponent, as the $1$ of $f^{1+q}$, leaves the base itself, which `_collect_generators` reads as the polynomial it is.
Before that change the two spellings of $h$ that differentiation produces, $1 + m^2\sin^2\theta/(r^2 - 2mr)$ and $1 - m^2\sin^2\theta/(2mr - r^2)$, were two generators and no Christoffel symbol cancelled.
Both charts declare `f` and `h` among their parameters as names for those ratios, as Weyl's chart of the Curzon-Chazy particle declares `R`, and `named_powers` in `print_charts.py` writes the part of each exponent that holds the parameter on the names, as $f^{2q}h^{q(2+q)}$, solving for the two exponents from those of the irreducible factors.
The prolate spheroidal chart writes $(x - 1)(x + 1)$, $(1 - y)(1 + y)$, $(x - y)(x + y)$ and $(\delta - 1)(\delta + 1)$ as their products.

## Step 3. The curvature

The Kretschmann scalar is Quevedo's, with the $m^2$ his equation drops restored:

$$K = \frac{16m^2(1 + q)^2\left(3(r - 2m - mq)^2(r^2 - 2mr + m^2\sin^2\theta) + m^2q(2 + q)\left(m^2q(2 + q) + 3(r - m)(r - 2m - mq)\right)\sin^2\theta\right)f^{2q}h^{2q(2+q)}}{r^6(r - 2m)^2(r^2 - 2mr + m^2\sin^2\theta)}.$$

On the axis it is $48m^2(1 + q)^2(r - (2 + q)m)^2(r - 2m)^{2q-2}/r^{2q+6}$: finite at $r = 2m$ for $q \geq 1$, where it is $3/(4m^4)$ at $q = 1$, and divergent for $q < 1$, $q \neq 0$.
Off the axis it diverges at $r = 2m$ as $(r - 2m)^{-2(q^2 + q + 1)}$, with the coefficient $m^4q^2(q^2 + q + 1)$ of the numerator on the equator, which vanishes only at $q = 0$.
So the curvature at the ends of the rod depends on the direction of approach, as Virbhadra and Kodama and Hikida found.

## Step 4. The chart that is not published

Weyl's canonical chart, with $R_1$ and $R_2$ the distances from the ends of the rod, reads through the checker and its Ricci tensor vanishes, in 38 seconds.
Its Kretschmann scalar took 260 seconds in the checker's `Geometry`, over the 120 seconds `verify_metrics.py` allows a tensor, so the chart would come back `UNCHECKED`.
The prolate spheroidal chart's convention gives $\rho$ and $z$, and the descriptions of `f` and `h` give $\psi$ and $\gamma$.

## Step 5. The two planes

The axis and the equatorial plane are totally geodesic, the first fixed by the rotations and the second by the reflection $\theta \to \pi - \theta$, so their null curves are null geodesics.
On the axis $h = 1$ and the metric is $-f^{1+q}c^2dt^2 + f^{-1-q}dr^2$, with $g_{tt}g_{rr} = -1$, so $r$ is an affine parameter along every ray.
At $q = 1$, $r_* = r + 4m\ln(r - 2m) - 4m^2/(r - 2m)$ falls to $-\infty$ at $r = 2m$: a ray reaches the end of the rod after a finite affine distance and at $t = \pm\infty$.
At $q = -1/2$, $r_* = \sqrt{r(r - 2m)} + 2m\ln(\sqrt{r} + \sqrt{r - 2m})$ is finite there.
On the equator $h = (r - m)^2/(r(r - 2m))$ and $c\,dt/dr = f^{-1-q}h^{-q(2+q)/2}$, which grows as $(r - 2m)^{(q^2 - 2)/2}$ and is integrable at $r = 2m$ for every $q \neq 0$: $r^{7/2}/(\sqrt{r - 2m}\,(r - m)^3)$ at $q = 1$ and $f^{-7/8}(1 - m/r)^{3/4}$ at $q = -1/2$.

## Step 6. The diagrams

The spacetime diagrams draw the two planes in each chart at $m = 1$ for $q = 1$ and $q = -1/2$, and `null_rays.py --verify` compares their rays with $ct \pm r_*$, to $8 \times 10^{-12}$ and $8 \times 10^{-10}$ of the drawing.
The prolate equator's curvature grows as the $3/2$ power, a factor of $32$ for each tenth of the distance, under the fiftyfold the edge test asks for, so those two rows declare the edge in `singular_zero`.
The conformal diagrams take $u, v = \arctan((ct \mp r_*)/\ell)$ with $\ell = 4m$: the oblate axis is the whole diamond, its left edges $r = 2m$ drawn as the edge of the chart with the Kretschmann scalar checked to be $3/(4m^4)$ there, and the other three planes are Minkowski's triangle with $r = 2m$ on $X = 0$, drawn as a singularity with the Kretschmann scalar checked to diverge.
The embedding diagram is the equatorial plane at $t = 0$, whose circles have radius $rf^{-q/2}$.
At $q = 1$ that is $r^{3/2}/\sqrt{r - 2m}$, least, $3\sqrt{3}\,m$, at $r = 3m$, and $g_{rr} \geq (d\rho/dr)^2$ reduces to $r^2(r - 2m)^2 \geq (3m - r)(r - m)^3$, which holds from $r = 2.5161\,m$ out.
At $q = -1/2$ the circles shrink to zero at $r = 2m$ and the condition is $(r - m)^{3/2}r^{1/4}(r - 2m)^{1/4} \geq (r - 3m/2)^2$, which holds from $r = 2.0020\,m$ out, a circle of radius $0.356\,m$.
Each deformation's moment is marked on its own equatorial drawings and on none of the axis.
