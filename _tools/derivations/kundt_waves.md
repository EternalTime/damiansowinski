# Kundt waves

The plane-fronted gravitational waves of type N whose rays are not parallel, in Podolský and Beláň's coordinates

$$ds^2 = dx^2 + dy^2 - 4x^2\,du\,dv + 4\left(x^2v^2 + xG\right)du^2,\qquad \partial_x^2G + \partial_y^2G = 0.$$

Its five charts are written by `_tools/derivations/print_charts.py --metric kundt_waves`, and `verify_metrics.py --system kundt_waves/<chart>` checks each in seconds.

## Step 1. The line element and its signs

Podolský and Beláň's (1), Class. Quantum Grav. 21, 2811, is $2\,d\zeta\,d\bar\zeta - 2Q^2du\,dv + F\,du^2$ with $Q = \zeta + \bar\zeta$ and $F = 2Q^2v^2 - QH$, and their (8) writes it in $\zeta = (x + iy)/\sqrt{2}$ with $G = -H/2\sqrt{2}$, which is the line element above.
It is Bičák and Podolský's (15), J. Math. Phys. 40, 4495, the subclass $KN$ of the type N waves without twist.
Their signature is this collection's, so nothing is turned over.
The Ricci tensor is $R_{uu} = -2x(\partial_x^2G + \partial_y^2G)$ and nothing else, their $\Phi_{22} = Q^{-3}H_{,\zeta\bar\zeta}/2$ up to the tetrad's factor, and the Riemann tensor is their (10): $R_{xuxu} = -2x\,\partial_x^2G$, $R_{xuyu} = -2x\,\partial_x\partial_yG$, $R_{yuyu} = -2x\,\partial_y^2G$.
The geodesic equations of the chart are their (12) to (14).
$u$ and $v$ are pure numbers, since $u = \tan(\alpha/2)$ for the angle $\alpha$ of a front, so $G$ is a length and $w = 2x^2v$ an area.

## Step 2. The charts and their sources

- `kundt`: Kundt's canonical form, Podolský and Beláň's (2), $2\,d\zeta\,d\bar\zeta - 2\,du(dw + W\,d\zeta + \bar W\,d\bar\zeta + \mathcal{H}\,du)$ with $W = -2w/Q$ and $\mathcal{H} = -w^2/Q^2 + QH/2$, which they take from Kundt's paper of 1961 and from Stephani, Kramer, MacCallum, Hoenselaers and Herlt. In $x$ and $y$ it is $dx^2 + dy^2 - 2\,du\left(dw - (2w/x)\,dx - (w^2/2x^2 + 2xG)\,du\right)$, and $w$ is an affine parameter along the rays, $\Gamma^a{}_{ww} = 0$.
- `podolsky_belan`: their (8), under $w = 2x^2v$.
- `simplest_wave`: their (15) with $n = 2$ and $c$ constant, $G = (x^2 - y^2)/\ell$. The rescaling $u \to \lambda u$, $v \to v/\lambda$ turns $\ell$ into $\lambda^2\ell$, as they note, so $\ell$ has no invariant meaning.
- `kerr_schild`: their (11), from Griffiths and Podolský, Class. Quantum Grav. 15, 3863: under their (10), $X = x(1 + 2uv)$, $Z = x(v - u(1 + uv))$, $cT = x(v + u(1 + uv))$, $Y = y$, the metric is flat space's plus $(G/x)\,k \otimes k$ with $k = (1 + u^2)\,c\,dT - 2u\,dX - (1 - u^2)\,dZ = 2x\,du$, which is null in both metrics. The chart publishes it for the simplest wave, with $x = \sqrt{X^2 + Z^2 - c^2T^2}$, $u = (X - x)/(cT + Z)$ and $v = (cT + Z)/2x$ as names.
- `ozsvath_robinson_rozga`: the whole class $KN(\Lambda)[\alpha, \beta]$ of Ozsváth, Robinson and Rózga, J. Math. Phys. 26, 1755, as Bičák and Podolský's (2) writes it, with $\alpha$ and $\beta$ constant. Its $\xi$ and $\eta$ are the real ones of the chart on the disc of Siklos's waves, $\xi_{\rm complex} = (\xi + i\eta)/\sqrt{2}$, its $h$ is minus their $H$, as there, and its $\beta$ is $\sqrt{2}$ times theirs, taken real, so that $q = \alpha(1 - \Lambda(\xi^2 + \eta^2)/12) + \beta\xi$ and $\kappa = \Lambda\alpha^2/3 + \beta^2$ hold no root. The field equation is their $H_{,\xi\bar\xi} + (\Lambda/3p^2)H = 0$, here $p^2(\partial_\xi^2h + \partial_\eta^2h) + 2\Lambda h/3 = 0$.

`kundt_check` in `print_charts.py` holds each chart to its field equation, the rays to being null geodesics and a repeated principal null direction of the Weyl tensor, $C_{abcd}k^d = 0$, which is type N, the second and third charts to the first pulled back, the Kerr-Schild chart to the third pulled back, and the last chart to the published Brinkmann chart of `pp_wave` at $\Lambda = 0$, $\alpha = 1$, $\beta = 0$, to the published chart on the disc of `siklos` at $\Lambda = -3/L^2$, $\alpha = 1$, $\beta = 1/L$, and to Podolský and Beláň's chart at $\Lambda = 0$, $\alpha = 0$, $\beta = 1/\ell$ under $u \to \sqrt{2}\,\ell u$, $v \to \sqrt{2}\,\ell v$ and $h = 2G/\ell$.

## Step 3. The Kerr-Schild chart and its three names

Written out in $T$, $X$ and $Z$, each of the chart's four hundred values is a polynomial of the sixth degree in $u$, which holds the root $x$ in its numerator, and the checker did not finish the chart in ten minutes.
So $x$, $u$ and $v$ are held as functions of $T$, $X$ and $Z$, `HELD` in `verify_metrics.py`, and their first derivatives are declared in `RATES`, each a rational function of the three:

$$dx = -(v + u(1 + uv))\,c\,dT + (1 + 2uv)\,dX + (v - u(1 + uv))\,dZ,\qquad 2x\,du = k,\qquad 2x\,dv = c\,dT + dZ - 2v\,dx.$$

The reader checks each declared derivative against the name's own definition and then writes every derivative of a name by them, after each product, so every comparison is one between rational functions of $x$, $u$, $v$ and $Y$.
The three are the coordinates of Podolský and Beláň's chart, so no relation holds among them, and a value that vanishes is exactly zero.
`norm` writes the root as $i\sqrt{c^2T^2 - X^2 - Z^2}$, since it leads a radicand with a positive term, which is harmless here: the root is only met where a declared derivative is checked, and there both sides hold it.
The chart takes 77 seconds to print and 4 to check.

The chart ends where $cT + Z = 0$ with $X < 0$, the front $\alpha = \pm\pi$, where $u$ is infinite.
For this wave, whose profile does not depend on $u$, the coefficient of the regular null covector $c\,dT - \sin\alpha\,dX - \cos\alpha\,dZ$ grows there as $(1 + u^2)^2$; a profile that falls off in $u$ does not have that edge.

## Step 4. The planes drawn

The plane of $u$ and $v$ at fixed $x$ and $y$ has the metric $-4x^2du\,dv + 4(x^2v^2 + xG)\,du^2$.
With no wave it is flat space's metric on the hyperboloid $X^2 + Z^2 - c^2T^2 = x^2$, de Sitter space in two dimensions, so it is not totally geodesic even then: $\Gamma^x{}_{uu} = -4xv^2 - 2G - 2x\,\partial_xG$ and $\Gamma^x{}_{uv} = 2x$.
Its rays, $du = 0$, are the rays of the wave and null geodesics; its other null curves obey $dv/du = v^2 + G/x$ and are null curves of the plane.
The free charts are drawn for $G = e^{-4u^2}(x^2 - y^2)/\ell$, a solution since $u$ never enters the field equation, and the simplest wave at $x = \ell$, where $dv/du = v^2 + 1$ and a curve keeps $u - \arctan v$, which `null_rays.py --verify` checks.
The canonical chart is drawn at $y = 1.1\,x$, where $G$ is negative: its $g^{ww} = (3w^2 - 4x^3G)/x^2$ changes sign where $G$ is positive, on a curve that is no horizon, and the drawing would mark it as one.
Its time function is $w + 30u$, whose gradient is timelike while $g^{ww} < 60$, which holds in the box; $u + v$ serves the charts whose $g^{uu}$ and $g^{vv}$ have no cross term with $x$.

The simplest wave has the Killing vector $\partial_u$, spacelike on $y = 0$, and that slice is totally geodesic, since $\Gamma^y{}_{uu} = 4xy/\ell$ vanishes on it.
With $\partial_u$ divided out the plane of $v$ and $x$ has $h = dx^2 - x^2dv^2/(v^2 + x/\ell)$, whose null curves are the shadows of the null geodesics with no momentum along $u$.

The family with a cosmological constant is drawn on its plane of $u$ and $v$ at $\xi = \ell$, $\eta = 0$ for two members: $\alpha = 1$, $\beta = 0$, $\Lambda = 3/\ell^2$, and $\alpha = 0$, $\beta = 1/\ell$, $\Lambda = -3/\ell^2$, each with $\kappa = 1/\ell^2$.
With $\Lambda > 0$ the forms $[\alpha = 1, \beta = 0]$ and $[\alpha = 0, \beta = 1]$ are one family, Griffiths, Docherty and Podolský's section 5, and the first has its envelope $q = 0$ on the circle $\xi^2 + \eta^2 = 4\ell^2$ about the chart's origin, which is what lets the front be drawn as a surface of revolution.
The profile is $h = e^{-4u^2/\ell^2}(\xi^2 - \eta^2)(2 + p)/(3p\ell^2)$, Bičák and Podolský's (5) for $f$ cubic in $\xi$, the first that radiates.

## Step 5. The fronts on flat space

On the slice $Y = 0$ a front is $cT = X\sin\alpha + Z\cos\alpha$ with $x \ge 0$, and its events are $cT(\sin\alpha, \cos\alpha) + x(\cos\alpha, -\sin\alpha)$ in $(X, Z)$.
Its rays are the lines of constant $x$ along $(cT, X, Z) = (1, \sin\alpha, \cos\alpha)$, parallel on one front and turned with $\alpha$ from front to front, and the generator of the cone $X^2 + Z^2 = c^2T^2$ where the front touches it is the ray $x = 0$.
`projections.kundt_fronts` draws six fronts with two rays each and checks every event drawn against the chart's own $x$ and $u$, and every ray to be null by the published metric.
That holds for any $\ell$, since the Kerr-Schild term is proportional to $k \otimes k$.

## Step 6. The wave front as a surface

A surface of constant $u$ and $v$ has the metric $(d\xi^2 + d\eta^2)/p^2$, of curvature $\Lambda/3$.
With $\Lambda = 3/\ell^2$ it is a sphere of radius $\ell$ in stereographic coordinates, the circle $r$ at $\theta = 2\arctan(r/2\ell)$ from the centre, and the front of the member $\alpha = 1$, $\beta = 0$ is the hemisphere $q > 0$, with $q/p = \cos\theta$, bounded by the envelope $q = 0$, its equator: Griffiths, Docherty and Podolský's hemispheres tangent to the expanding torus.
With $\Lambda = -3/\ell^2$ it is the hyperbolic plane, the sheet $Z = \ell\cosh(s/\ell) - \ell$ of a hyperboloid in Minkowski space, and the front of the member $\alpha = 0$, $\beta = 1/\ell$ is the half $\xi > 0$, bounded by the geodesic $\xi = 0$; there $q/p = \xi/p\ell$ is the drawing's own $X$, so the lines of constant $q/p$ are the sections $X = \ell/2$, $\ell$ and $2\ell$, each at a fixed distance from the envelope.
With no cosmological constant the front is a flat half plane, which an embedding diagram would draw as a flat sheet, so it has none.

## Step 7. What is not drawn

There is no conformal diagram.
No plane of two coordinates in any chart is totally geodesic and carries both null families, and the one reduction to two dimensions, the simplest wave's slice $y = 0$ with $\partial_u$ divided out, follows only the geodesics with no momentum along $u$.
The causal structure the History tells of, flat space outside a light cone with the fronts rolled round it, lies in three dimensions and is drawn by the figure of Step 5.
