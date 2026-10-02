# Siklos waves and the Kaigorodov spacetime

The exact gravitational waves of anti-de Sitter space, a pp-wave in Brinkmann's form multiplied by anti-de Sitter's conformal factor,

$$ds^2 = \frac{L^2}{x^2}\left(dx^2 + dy^2 - 2\,du\,dv + H\,du^2\right),\qquad \Lambda = -\frac{3}{L^2}.$$

Its eight charts are written by `_tools/derivations/print_charts.py --metric siklos`, and `verify_metrics.py --system siklos/<chart>` checks each in seconds.

## Step 1. The line element and its sign

Podolský's (1), Class. Quantum Grav. 15, 719, is $ds^2 = (\beta^2/x^2)(dx^2 + dy^2 + 2\,du\,dv + H\,du^2)$ with $\beta = \sqrt{-3/\Lambda}$, which he takes from Siklos's chapter of 1985.
Here $\beta$ is $L$, as on the anti-de Sitter page, and $v$ is the opposite of his, so that the bracket is the pp-wave page's $-2\,du\,dv + H\,du^2$ with $u$ and $v$ lengths.
His dimensionless $x$, $y$, $u$ and $v$ are ours divided by $L$, so his Kaigorodov profile $H = x^3$ is $x^3/L^3$.
The Ricci tensor of the chart is $-(3/L^2)g_{\mu\nu}$ in every slot but $uu$, where $R_{uu} + (3/L^2)g_{uu} = -\tfrac{1}{2}(\partial_x^2H + \partial_y^2H - (2/x)\partial_xH)$, his (2).
The Ricci scalar is $-12/L^2$ and the Kretschmann scalar $24/L^4$ for every $H$, anti-de Sitter's values, since $g^{uu} = 0$.

## Step 2. The charts and their sources

- `siklos`: Podolský's (1), with $H(u, x, y)$ free.
- `ozsvath_robinson_rozga`: the subclass $(IV)_0$ of Ozsváth, Robinson and Rózga, J. Math. Phys. 26, 1755, as Podolský's (4) and Bičák and Podolský's (14), J. Math. Phys. 40, 4495, write it, $2\,d\xi\,d\bar\xi/p^2 - 2(q^2/p^2)\,du\,dv - (q/p)H\,du^2$. With $\xi_{\rm complex} = (\xi + i\eta)/\sqrt{2}$ and $\Lambda/6 = -1/(2L^2)$ their $p$ and $q$ are $p = 1 - (\xi^2 + \eta^2)/4L^2$ and $q = (1 + \xi/2L)^2 + \eta^2/4L^2$. The map $\xi + i\eta = 2L(L - x - iy)/(L + x + iy)$ carries Siklos's half plane onto the disc of radius $2L$ with $q/p = L/x$, so $u$ and $v$ are unchanged and the profile is $h = LH/x$, written with the sign of the pp-wave page. Podolský's (5) uses the Möbius map onto the outside of the disc; this one keeps the inside, where $p > 0$. The field equation becomes $p^2(\partial_\xi^2h + \partial_\eta^2h) = 2h/L^2$, and at $1/L = 0$ the line element is $d\xi^2 + d\eta^2 - 2\,du\,dv + h\,du^2$.
- `kaigorodov`: Podolský's (26), $H = x^3/L^3$.
- `kaigorodov_poincare`: Brecher, Chamblin and Reall's (5.1) and (5.2), Nucl. Phys. B 607, 155, $(l^2/z^2)(dx^+dx^- + \mu^3z^3(dx^+)^2 + dy^2 + dz^2)$ with $x^\pm = x \pm t$. Here $\sqrt{2}\,u = ct - x$, $\sqrt{2}\,v = ct + x$, and $z$ is Siklos's $x$, which gives the wave term $(z^3/2L^3)(c\,dt - dx)^2$.
- `kaigorodov_horospheric`: Kaigorodov's own form with the upper sign, Podolský's (27), $(dx^4)^2 + e^{2x^4/\beta}(2\,dx^1dx^3 + (dx^2)^2) + e^{-x^4/\beta}(dx^3)^2$, under his (28), $\rho = -L\ln(x/L)$. It is Cvetič, Lü and Pope's (A.8) at $n = 1$, Nucl. Phys. B 545, 309, whose $L$ is $1/2L$ here.
- `kaigorodov_stationary`: the same with the lower sign, the region $x < 0$ of Siklos's chart, which is Siklos's chart at $H = -x^3/L^3$ after $x \to -x$. Podolský's section 7: there $\partial_u$ is timelike.
- `kaigorodov_homogeneous`: Podolský's (30), which he quotes from Kramer, Stephani, MacCallum and Herlt (1980), their (10.33), with $-12/\Lambda = 4L^2$, under his (31): $x = Le^{2Z}$, $u = -\sqrt{10k}\,X$, $v = -Ue^{5Z}/\sqrt{10k}$ in our lengths.
- `kaigorodov_kundt`: Podolský's (32), their (33.2) with the misprint he corrects, under his (33): $u = \sqrt{2}\,U$, $v = Vx^2/\sqrt{2}$, with $x$ and $y$ pure numbers.

`siklos_check` in `print_charts.py` holds each chart to $R_{\mu\nu} = -(3/L^2)g_{\mu\nu}$, the two free charts where their profile solves its equation, each chart after the first to being Siklos's pulled back through these maps, and $\partial_v$ to being a repeated principal null direction of the Weyl tensor, $C_{abcd}k^d = 0$, which is type N.

## Step 3. What is not a chart

Podolský's (50) writes Kaigorodov's metric in the global coordinates of anti-de Sitter space, anti-de Sitter's metric plus $(\beta^5\cos\chi/2D^5)$ times the square of a one-form in all four coordinates.
Its curvature components are sums of several hundred terms each, which no reader would use, so the chart is left out and its content, that $x = 0$ is an infinity like anti-de Sitter's, is carried by the conformal diagram.

## Step 4. The planes drawn

The plane of $u$ and $v$ at fixed $x$ and $y$ has $\Gamma^x{}_{uv} = -1/x$ and $\Gamma^x{}_{uu} = H/x - \partial_xH/2$.
Along the null curve $dv = \tfrac{1}{2}H\,du$ these give $\ddot{x} = \tfrac{1}{2}\partial_xH\,\dot{u}^2$: where $H = 0$ the curve is a null geodesic, as both families are in anti-de Sitter space, and inside a wave it is pushed toward larger $x$.
The free charts are drawn for $H = e^{-4u^2/L^2}x^3/L^3$, a solution since $u$ never enters the field equation, where a curve moving left keeps $v - \sqrt{\pi}\,x^3\,\mathrm{erf}(2u/L)/8L^2$ and is carried by $\sqrt{\pi}\,x^3/4L^2$ in all.

Kaigorodov's spacetime has the Killing vectors $\partial_u$, $\partial_v$, $\partial_y$, $u\partial_y - y\partial_v$ and $5v\partial_v + 2x\partial_x + 2y\partial_y - u\partial_u$, Podolský's (34).
No plane of two coordinates through $\partial_v$ is totally geodesic, and neither is the orbit of $\partial_v$ and the fifth, $u^2x = $ const, whose second fundamental form has the component $u^3/x$.
So its planes are drawn on the slice $y = 0$, totally geodesic by the reflection of $y$, with a spacelike Killing direction divided out, as the rotating BTZ hole's are.
With $\partial_u$ divided out, $h = -(L^5/x^5)\,dv^2 + (L^2/x^2)\,dx^2 = (L^5/x^5)(-dv^2 + dx_*^2)$ with $x_* = \tfrac{2}{5}x^{5/2}/L^{3/2}$.
With $\partial_x$ of the Poincaré chart divided out, $h = (L^2/z^2)(-c^2dt^2/(1 + z^3/2L^3) + dz^2)$.
The first plane's rays are the null geodesics with no momentum along $u$, Podolský's (41) at $A = B = 0$, which reach $x = \infty$ at a finite affine parameter and $v = \pm\infty$.
The stationary region has no such quotient, since there $\partial_u$ is timelike and dividing it out leaves a Riemannian plane, so it is drawn on its plane of $u$ and $v$ at $\rho = 0$, $-du^2 - 2\,du\,dv$.

## Step 5. The singularity

Every curvature scalar of Kaigorodov's spacetime is a constant, since it is homogeneous.
The timelike geodesics with no momentum along $u$ or $y$ and unit momentum along $v$ have the velocity $(\dot u, \dot v, \dot x, \dot y) = (x^2, x^5, x\sqrt{x^5 - 1}, 0)$ at $L = 1$, Podolský's (35) with $A = B = 0$.
Their tidal tensor $E_{ab} = R_{acbd}u^cu^d$ has $E_{ab}E^{ab} = 3 + \tfrac{9}{2}x^{10}$: anti-de Sitter's $3/L^4$ and twice the square of his amplitude $A_+ = -\tfrac{3}{2}C^2x^5$.
They reach $x = \infty$ from $x_0$ in the proper time $\tfrac{2}{5}L\arcsin(x_0^{-5/2})$.
`conformal.py` checks all three against the published Riemann tensor and Christoffel symbols before it draws the null edges of the triangle as a singularity.

## Step 6. The wave front

A surface of constant $u$ and $v$ has the metric $(L^2/x^2)(dx^2 + dy^2) = (d\xi^2 + d\eta^2)/p^2$, the hyperbolic plane of curvature $-1/L^2$.
On the disc the circle a proper distance $s$ from the centre has the coordinate radius $2L\tanh(s/2L)$ and the circumference $2\pi L\sinh(s/L)$, so it is the sheet $Z = L\cosh(s/L) - L$ of a hyperboloid in Minkowski space, as the hyperbolic horizon of the topological black holes is.
The lines of constant $x$ are horocycles through the point $\xi = -2L$ of the rim, and a chord of the hyperboloid between two points of one is the arc $L\,dy/x$ between them, which `embedding.py` checks.
