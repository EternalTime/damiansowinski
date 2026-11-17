# The dimensionally continued black holes

Máximo Bañados, Claudio Teitelboim and Jorge Zanelli, Phys. Rev. D 49, 975 (1994), write the static, spherically symmetric field of Lovelock's theory with every coefficient set by one length $\ell$ as $ds^2 = -g^2\,dt^2 + dr^2/g^2 + r^2d\Omega^2$, their Eqs. (22) and (23), with $g^2$ of their Eq. (21).
In $D = 2n$ dimensions $g^2 = 1 - (2M/r)^{1/(n-1)} + r^2/\ell^2$, and in $D = 2n - 1$ dimensions $g^2 = 1 - (M + 1)^{1/(n-1)} + r^2/\ell^2$.
A chart needs its coordinates named one by one, so the page publishes $n = 3$: five dimensions for the odd family and six for the even one, written in `dimensionally_continued.py`.

## Step 1: the root

At $n = 3$ the exponent is $1/2$, and the 1994 article sets these dimensions aside, since with the negative root a positive mass has a naked singularity.
Javier P. Muniain and Dardo Píriz, Phys. Rev. D 53, 816 (1996), keep them with the positive root, which has a regular horizon, and every chart takes it.
In five dimensions $f = 1 - \sqrt{M + 1} + r^2/\ell^2 = (r^2 - r_+^2)/\ell^2$ with $r_+^2 = \ell^2(\sqrt{M + 1} - 1)$, a black hole for $M > 0$ and anti-de Sitter space at $M = -1$.
In six $f = 1 - \sqrt{2M/r} + r^2/\ell^2$, whose one zero satisfies $r_+(1 + r_+^2/\ell^2)^2 = 2M$, and $M = 0$ is anti-de Sitter space.
The six dimensional charts name $w = \sqrt{2M/r}$, a pure number, so that every value is a rational function of $r$, $\ell$ and $w$.

## Step 2: the other charts

The ingoing charts are Masato Nozawa and Hideki Maeda's Eq. (3.7), Class. Quantum Grav. 23, 1779 (2006), with a constant mass: $ds^2 = -f\,dv^2 + 2\,dv\,dr + r^2d\Omega^2$ with $dv = c\,dt + dr/f$.
The conformal chart of five dimensions is Kiyoshi Shiraishi and Takuya Maki's Eqs. (7) and (8), Class. Quantum Grav. 11, 1687 (1994): $r = r_+\sec\rho$ outside the horizon gives $ds^2 = \ell^2\sec^2\rho\,(-\kappa^2\sin^2\rho\,c^2dt^2 + d\rho^2 + (r_+^2/\ell^2)\,d\Omega_3^2)$ with $\kappa = r_+/\ell^2$.
`dimensionally_continued.check` holds each chart, at three random points, to the field equations $H_{ab} + 2(D-3)(D-4)G_{ab}/\ell^2 - (D-1)(D-2)(D-3)(D-4)g_{ab}/2\ell^4 = 0$, with $H_{ab}$ the Gauss-Bonnet tensor, each static chart to the source's metric function, horizon, surface gravity and anti-de Sitter limit, and every other chart to the static one pulled back.

## Step 3: the values the drawings take

Every drawing puts the horizon at $r_+ = \ell$, the unit of length: $M = 3$ in five dimensions, where $\sqrt{M + 1} = 2$ and $f = r^2 - 1$, and $M = 2\ell$ in six, where $f = 1 - 2/\sqrt{r} + r^2$.
The surface gravity $f'(r_+)/2$ is $\kappa = 1$ in five dimensions and $3/2$ in six.

## Step 4: the tortoise coordinate and the conformal diagrams

In five dimensions $1/f = \tfrac12\left(\tfrac{1}{r - 1} - \tfrac{1}{r + 1}\right)$, so $r_* = \tfrac12\ln|(r - 1)/(r + 1)|$, which vanishes at $r = 0$ and as $r \to \infty$.
Kruskal's $U = -e^{-\kappa u}$ and $V = e^{\kappa v}$, $u, v = ct \mp r_*$, give $UV = (1 - r)/(1 + r)$: $1$ at $r = 0$, $0$ on the horizon and $-1$ on the conformal boundary.
With $p = \arctan U$ and $q = \arctan V$ the singularity lies on $T = \pm\pi/2$ and the boundary on $X = \pm\pi/2$, the square of the black hole of three dimensions.
In the conformal chart $r_* = \ln\tan(\rho/2)$.
In six dimensions $x = \sqrt{r}$ gives $dr/f = 2x^2\,dx/(x^5 + x - 2)$, and $x^5 + x - 2 = (x - 1)(x^4 + x^3 + x^2 + x + 2)$ has five simple roots $x_i$, so $r_* = \mathrm{Re}\sum_i A_i\ln(1 - \sqrt{r}/x_i)$ with $A_i = 2x_i^2/(5x_i^4 + 1)$.
That $r_*$ vanishes at $r = 0$ and tends to $R = \mathrm{Re}\sum_i A_i\ln(-1/x_i) = 0.5223\,\ell$, since the $A_i$ sum to zero, and the residue at $x = 1$ is $1/3 = 1/2\kappa$.
So $UV$ is analytic through the horizon, $1$ at $r = 0$ and $-B$ on the boundary with $B = e^{2\kappa R} = 4.792$, and the boundary lies on $\tan p\tan q = -B$, bowed outward to $X = 2\arctan\sqrt{B} = 2.28$, as Schwarzschild-anti-de Sitter's is.
`dimensionally_continued_drawings.Hole` draws both in the cells of `conformal.Tower`, and `conformal.py` checks every cell against the published metric.

## Step 5: the spacetime diagrams

Each chart is drawn on its plane of the time and the radius at $\psi = \theta = \pi/2$, $\phi = 0$, and $\chi = \pi/2$ in six dimensions.
The static charts reach through the horizon to $r = 0$, oriented by the ingoing rays, and the ingoing charts are drawn against $v - r$.
`null_rays.py --verify` holds the rays to $ct \pm r_*$, to $v$ and $v - 2r_*$, and in the conformal chart to $ct \pm \ln\tan(\rho/2)$.

## Step 6: the embedding diagram

The slice $t = 0$ on the plane of $r$ and $\phi$ has $g_{rr} = 1/f$ and $g_{\phi\phi} = r^2$.
Between the horizon and the circle where $f = 1$, $r = \sqrt{2}\,\ell$ in five dimensions and $r = 2^{2/5}\ell$ in six, $g_{rr} \ge 1$ and the slice is drawn in flat space; beyond it $g_{rr} < 1$ and it is drawn in three dimensional Minkowski space, as the slices of the black hole of three dimensions and of Schwarzschild-anti-de Sitter are.
Both exteriors are drawn, joined at the horizon, out to $r = 3\ell$.

## Step 7: the curvature

In five dimensions, with $s = \sqrt{M + 1}$, $K = 40/\ell^4 - 24s/\ell^2r^2 + 12s^2/r^4$ and $R = -20/\ell^2 + 6s/r^2$.
$K$ is least, $28/\ell^4$ at every $M$, at $r^2 = s\ell^2$, and $R$ vanishes at $r^2 = 3s\ell^2/10$: $r = \sqrt{2}\,\ell$ and $r = \sqrt{3/5}\,\ell$ at $M = 3$.
In six dimensions $K = 60/\ell^4 - 35w/\ell^2r^2 + 425w^2/16r^4$ and $R = -30/\ell^2 + 35w/4r^2$; at $M = 2\ell$, with $y = (\ell/r)^{5/2}$, $K$ is least, $824/17\ell^4 = 48.5/\ell^4$, at $y = 28/85$, $r = 1.56\,\ell$, and $R$ vanishes at $y = 12/7$, $r = 0.806\,\ell$.
Far out both tend to anti-de Sitter space's, $2D(D-1)/\ell^4$ and $-D(D-1)/\ell^2$, and both diverge at $r = 0$.
`curvature.py` checks each of these from the published scalars in sixty digits.
