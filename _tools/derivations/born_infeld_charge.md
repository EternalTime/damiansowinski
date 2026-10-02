# Born and Infeld's point charge

The point charge of Born and Infeld's electrodynamics with its own gravity, Hoffmann's solution of 1935, in the notation of the collection, is

$$ds^2 = -f\,c^2dt^2 + \frac{dr^2}{f} + r^2\,d\Omega^2,\qquad f = 1 - \frac{2m(r)}{r},\qquad \frac{dm}{dr} = \frac{r_q^2}{r^2 + W},\qquad W = \sqrt{r^4 + r_0^4},$$

with $m \to r_s/2$ as $r \to \infty$.
Its three charts are written by `_tools/derivations/print_charts.py --metric born_infeld_charge`, and `verify_metrics.py --system born_infeld_charge/<chart>` checks each in seconds.

## Step 1. The parameters

The parameters are three lengths.
$r_s = 2GM/c^2$ is Schwarzschild's, for the whole mass, and $r_q$, with $r_q^2 = GQ^2/(4\pi\epsilon_0c^4)$, is Reissner and Nordström's, for the charge, as `rn_metric` writes them.
$r_0$ is Born and Infeld's radius, $r_0^2 = Q/(4\pi\epsilon_0b)$, where the Coulomb field of the charge would reach their absolute field $b$: their $r_0 = \sqrt{e/b}$ in Gaussian units (Proc. R. Soc. A 144, 425, 1934, eq. 8.8).
Maxwell's theory is $b \to \infty$, which is $r_0 = 0$.

The literature writes the same three numbers in several ways, and "b" is not one letter.
Born and Infeld's $b$ and Bretón's are the greatest field; Gibbons and Rasheed's and Rasheed's $b$ is its inverse, and Fernando and Krug's $\beta$ is the greatest field again.
In Rasheed's geometrical units (hep-th/9702087, section 6) $m'(r) = (\sqrt{r^4 + b^2Q^2} - r^2)/b^2$, which is the slope above with $r_q = Q$ and $r_0^2 = bQ$, since $(W - r^2)/r_0^4 = 1/(r^2 + W)$.
Bretón's $a$ (Phys. Rev. D 67, 124004, eq. 5) is $r_0$.

## Step 2. The field and its stress

In lengths, with the electric field a reciprocal length, the Lagrangian is

$$L(F) = \beta^2\left(\sqrt{1 + \frac{2F}{\beta^2}} - 1\right),\qquad F = \tfrac{1}{4}F_{\mu\nu}F^{\mu\nu},\qquad \beta = \frac{r_q}{r_0^2},$$

which is Maxwell's $L = F$ in weak fields, and Einstein's equations are $G_\mu{}^\nu = 2\left(L_FF_{\mu\lambda}F^{\nu\lambda} - \delta_\mu^\nu L\right)$, as `bardeen.md` Step 4 writes them for Ayón-Beato and García's.
A radial electric field $E$ has $F = -E^2/2$, and the field equation $d(r^2L_FE)/dr = 0$ with the charge $r_q$ gives

$$E = \frac{r_q}{W},$$

finite at the charge, where it is $\beta$.
Then $L_F = W/r^2$ and $L = \beta^2(r^2/W - 1)$, so

$$G^t{}_t = G^r{}_r = -2\left(L_FE^2 + L\right) = -\frac{2r_q^2}{r^2\left(r^2 + W\right)},\qquad G^\theta{}_\theta = G^\phi{}_\phi = -2L = \frac{2r_q^2}{W\left(r^2 + W\right)}.$$

`born_infeld_check` derives both from the Lagrangian in sympy and holds every chart to them.
For $f = 1 - 2m/r$, $G^t{}_t = -2m'/r^2$, which is the slope of the mass function above: $dm/dr$ is the energy of the field in a shell, $4\pi Gr^2\rho/c^4$.
At $W = r^2$ the density is Maxwell's, $r_q^2/r^4$ in these units.

## Step 3. The mass function

No elementary function has the slope $r_q^2/(r^2 + W)$.
With $I(r) = \int_r^\infty dx/\sqrt{x^4 + r_0^4}$,

$$m(r) = \frac{r_s}{2} - \int_r^\infty\frac{r_q^2\,dx}{x^2 + W(x)} = \frac{r_s}{2} + \frac{r_q^2\,r}{3\left(r^2 + W\right)} - \frac{2r_q^2}{3}\,I(r),$$

which is checked by differentiating: the middle term's slope is $r_q^2(r^2 + W - 2r^4/W)/(3(r^2 + W)^2) = r_q^2(W - r^2)/(3W(r^2 + W))$, and the last term's is $2r_q^2/3W$.
$I$ is an incomplete elliptic integral of the first kind.
With $x = r_0\cot(\vartheta/2)$, $dx/\sqrt{x^4 + r_0^4} = -d\vartheta/(2r_0\sqrt{1 - \tfrac{1}{2}\sin^2\vartheta})$, so

$$I(r) = \frac{1}{2r_0}\,\mathrm{F}\!\left(2\arctan\frac{r_0}{r}\;\Big|\;\frac{1}{2}\right),\qquad \mathrm{F}(\varphi \mid k^2) = \int_0^\varphi\frac{d\vartheta}{\sqrt{1 - k^2\sin^2\vartheta}},$$

which is how the charts define $m$.
It is Bretón's eq. (8), $g(r) = \frac{1}{2a}F\left(\arccos\frac{r^2 - a^2}{r^2 + a^2}, \frac{1}{\sqrt{2}}\right)$, which she takes from García, Salazar and Plebański, with the amplitude written as an arctangent: $\cos(2\arctan(r_0/r)) = (r^2 - r_0^2)/(r^2 + r_0^2)$.
The arctangent keeps every digit far from the charge, where the arccosine of a number next to one loses half of them.

At the centre the amplitude is $\pi$ and $\mathrm{F}(\pi \mid \tfrac{1}{2}) = 2K(\tfrac{1}{2}) = \Gamma(\tfrac14)^2/(2\sqrt{\pi})$, so

$$m(0) = \frac{r_s}{2} - \frac{\Gamma(\tfrac14)^2}{6\sqrt{\pi}}\,\frac{r_q^2}{r_0} = \frac{r_s}{2} - 1.2360\,\frac{r_q^2}{r_0}.$$

The number is $1.23604978\ldots$, Gibbons and Rasheed's (8.11); Born and Infeld print it as $1.2361$ in the energy of the electron, $1.2361\,e^2/r_0$, their (8.6).
Far away $m = r_s/2 - r_q^2/2r + r_q^2r_0^4/(40r^5) - \ldots$, so $f$ falls short of Reissner and Nordström's by $r_q^2r_0^4/(20r^6)$, which the script checks against the published metric of `rn_metric` at $r_0 = 10^{-3}$, and at $r_q = 0$ it is the published metric of `schwarzschild`.

## Step 4. The three cases

$f \to 1 - 2m(0)/r - 2r_q^2/r_0^2 + O(r^2)$ at the centre, since $m'(0) = r_q^2/r_0^2$, so the mass left at the centre decides the spacetime, as Gibbons and Rasheed (Nucl. Phys. B 454, 185, section 8) and Rasheed sort it.

- $m(0) > 0$, more mass than the field holds: $f \to -\infty$ at the centre, as Schwarzschild's does, and $f$ has one zero. The singularity is spacelike.
- $m(0) < 0$: $f \to +\infty$, as Reissner and Nordström's does, and there are two horizons, one degenerate horizon, or none.
- $m(0) = 0$, Hoffmann's particle, $r_s = 2 \times 1.2360\,r_q^2/r_0$: $f \to 1 - 2r_q^2/r_0^2$, finite and not one. For $2r_q^2 \le r_0^2$ there is no horizon, since $(r - 2m)' = 1 - 2r_q^2/(r^2 + W)$ is then positive for $r > 0$ and $r - 2m$ vanishes at $r = 0$; Rasheed's "none otherwise" is the statement at equality, where Gibbons and Rasheed's "$\ge$" is a slip.

The centre of the particle is a conical singularity: the circle of radius $r$ lies a proper distance $r/\sqrt{1 - 2r_q^2/r_0^2}$ from it, so its circumference is $2\pi\sqrt{1 - 2r_q^2/r_0^2}$ times its radius.
It is a curvature singularity all the same, since the energy density grows as $2r_q^2/(r_0^2r^2)$: in the frame of a static observer the Riemann tensor has the three components $2(m - r_q^2r/W)/r^3$, $(m - r_q^2r/(r^2 + W))/r^3$ and $2m/r^3$, so

$$K = \frac{16}{r^6}\left(\left(m - \frac{r_q^2\,r}{W}\right)^2 + \left(m - \frac{r_q^2\,r}{r^2 + W}\right)^2 + m^2\right),$$

which grows as $16r_q^4/(r_0^4r^4)$ at the particle's centre, where $m \approx r_q^2r/r_0^2$, and as Schwarzschild's $48m(0)^2/r^6$ where $m(0) \ne 0$.
Fernando and Krug's "no curvature singularity" for the particle, and their constant $-\tfrac{10}{3}\beta Q$ in the expansion at the centre, are both wrong: their own equation (16) gives $1 - 2\beta Q$.

## Step 5. The charts

The static chart is the one every paper uses: Gibbons and Rasheed's (8.1), Rasheed's, Bretón's (4).
The Eddington-Finkelstein charts are $u = ct - r_*$ and $v = ct + r_*$ with $dr_*/dr = 1/f$ and $r_* = 0$ at $r = 0$; the outgoing one is Cirilo Lombardo's eq. (11) (Class. Quantum Grav. 21, 1407, 2004), the first step of his Newman-Janis construction, and the ingoing one its time reverse.
`born_infeld_check` holds each to being the static chart pulled back.
No paper prints an isotropic, Painlevé-Gullstrand or Kruskal chart of this spacetime in closed form, and none is published here; the conformal diagram draws its own Kruskal map.

## Step 6. Holding the mass function in the checker

The reader reads `\mathrm{F}\left(\varphi \mid m\right)` as `EllipticF`, sympy's `elliptic_f` under a name of its own: it differentiates along its amplitude by its integrand, evaluates through mpmath to any number of digits, and carries its own numbers for `lambdify`, scipy's `ellipkinc` for floats and mpmath's `ellipf` for mpmath's numbers.
`HELD` holds $m$ as a function of $r$ while the tensors are built, and `RATES` declares its slope, `\dfrac{r_q^2}{r^2 + W}`.
The reader checks that slope once against the definition, which is where the elliptic integral is differentiated: sympy writes $\sin(2\arctan(r_0/r))$ as $2r_0r/(r^2 + r_0^2)$, the integrand becomes $(r^2 + r_0^2)/W$, and `norm` cancels the rest.
From then on every comparison is between rational functions of $r$, $W$, $r_q$, $r_0$ and $m$, with no relation among them, and each chart checks in under a second.

Every value is printed in $r$, $W$, $m$ and $r_q$: `born_infeld_pretty` writes $r_0^4 = W^2 - r^4$, after which the $W - r^2$ that $r_0^4$ brings cancels, and a numerator is grouped by the powers of $m$, as $m(r^2 + W) - r_q^2r$, with $r - 2m$, whose zeros are the horizons, kept in the order the line element writes it.

## Step 7. The diagrams

Every diagram is drawn in units of $r_0$ at $r_q = r_0/2$, for two masses.
Hoffmann's particle has $r_s = 2 \times 1.2360 \times \tfrac14\,r_0 = 0.6180\,r_0$, so $f$ rises from $1/2$ at the centre to $1$ far away and the cone at the centre has slope $1$.
The black hole has $r_s = 2r_0$: $m(0) = 0.691\,r_0$, one horizon at $r_h = 1.8666\,r_0$ with $\kappa = (1 - 2m'(r_h))/2r_h = 0.2490/r_0$, where Reissner and Nordström's of the same mass and charge has horizons at $1.8660$ and $0.1340\,r_0$.
The particle's $r_s$ has no short form, so its diagrams state it as $r_s \approx 0.618$, the `rounded` field of a `Diagram`.

The tortoise coordinate, `slices.born_infeld_rstar`, is Bardeen's construction (`bardeen.md` Step 6): the pole at the horizon, of residue $1/f'(r_h) = 1/2\kappa$, is taken out of $1/f$ and the smooth remainder is summed by Gauss and Legendre's rule, in forty digits within a twentieth of $r_0$ of the horizon and within a hundredth of the particle's centre, where $m$ is the difference of two numbers that agree.
The particle has no pole, and its $r_*$ runs from $0$ to infinity.

The spacetime diagrams are the plane of the time and $r$ in each chart, for each mass.
The conformal diagram of the particle is Minkowski's triangle, $p, q = \arctan((ct \mp r_*)/\ell)$ with $\ell = 2r_0$, its centre drawn as a timelike singularity; the black hole's is Kruskal and Szekeres's hexagon, $U = -e^{-\kappa u}$, $V = e^{\kappa v}$, with $r_*(0) = 0$ putting the singularity $UV = 1$ on $T = \pm\pi/2$.
The embedding diagram is the equator of $t = 0$: for the particle a surface that leaves the centre as a cone of slope $\sqrt{2r_q^2/(r_0^2 - 2r_q^2)} = 1$ and rises far out as Flamm's paraboloid of the same mass does, and for the black hole the two sheets through the bifurcation sphere.
The horizon is no root of a polynomial, so `Slice.horizon_between` finds it in eighty digits and the quadrature beside it treats it as it treats a root `Slice.horizons` finds.
The Reissner-Nordström-like case, $m(0) < 0$ with two horizons, has Reissner and Nordström's tower and is not drawn.
