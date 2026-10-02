# The anti-de Sitter soliton

Horowitz and Myers's soliton, Phys. Rev. D 59, 026005 (1998), arXiv:hep-th/9808079, their (3.14): the planar black hole of anti-de Sitter space, their (2.6), with $t \to i\tau$ and one flat direction $\to it$,

$$ds^2 = \frac{r^2}{L^2}\left[-c^2dt^2 + \left(1 - \frac{r_0^{p+1}}{r^{p+1}}\right)d\tau^2 + dx^idx^i\right] + \frac{L^2}{r^2}\left(1 - \frac{r_0^{p+1}}{r^{p+1}}\right)^{-1}dr^2$$

in $p + 2$ dimensions, with $\tau$ of period $\beta = 4\pi L^2/((p + 1)r_0)$.
Its five charts are written by `_tools/derivations/print_charts.py --metric ads_soliton`, and `verify_metrics.py --system ads_soliton/<chart>` checks each in seconds.

## Step 1. The parameters

The parameters are anti-de Sitter's radius $L$, Horowitz and Myers's $\ell$, and the radius $r_0$ of the tip, both lengths.
$\tau$ and the flat directions are lengths, as in their (3.14), and $t$ is a time, so $c$ stands where they set it to one.
The cosmological constant is $\Lambda = -p(p + 1)/(2L^2)$, their statement under (2.4): $-3/L^2$ in four dimensions, $-6/L^2$ in five and $-1/L^2$ in three.

## Step 2. The charts and their sources

- `horowitz_myers`: their (3.14) at $p = 2$, four dimensions, the dimension of the planar black hole of `topological_black_hole`, whose published brane chart it is with $g_{\tau\tau} = -g_{tt}$ and $g_{tt} = -g_{xx}$ at $z = L^2/r$ and $z_h = L^2/r_0$.
- `poincare`: the same in $z = L^2/r$, the coordinate of anti-de Sitter's Poincaré patch, in which Hartnoll's (42) and (43), Class. Quantum Grav. 26, 224002, write the black brane it is continued from; the tip is $z_0 = L^2/r_0$ and the period $4\pi z_0/3$. No paper is the source of this chart for the soliton itself: it is Horowitz and Myers's chart under $r = L^2/z$, which `ads_soliton_check` holds it to.
- `polar`: the proper distance from the tip, $d\rho = L\,dr/(r\sqrt{f})$, which integrates to $r^3 = r_0^3\cosh^2(3\rho/2L)$, and the angle $\phi = 2\pi\tau/\beta = 3r_0\tau/(2L^2)$. Then $f = \tanh^2(3\rho/2L)$ and $$ds^2 = \frac{r_0^2}{L^2}\cosh^{4/3}\frac{3\rho}{2L}\left(-c^2dt^2 + dx^2\right) + d\rho^2 + \frac{4L^2}{9}\frac{\sinh^2(3\rho/2L)}{\cosh^{2/3}(3\rho/2L)}d\phi^2,$$ with $g_{\phi\phi} = \rho^2 + O(\rho^4)$ at the tip, the plane in polar coordinates that Horowitz and Myers's period is chosen to give. It is our own chart, derived here, and `ads_soliton_check` holds it to their chart pulled back and to that expansion.
- `five_dimensional`: their (4.1), $p = 3$, the case their conjecture and their comparison with the gauge theory are stated for, with $\beta = \pi L^2/r_0$.
- `three_dimensional`: their (3.22) with $r_0$ kept, $p = 1$. It has constant curvature, and with the radius $\sqrt{r^2 - r_0^2}\,L/r_0$ and the time $r_0t/L$ it is their (3.23), anti-de Sitter space in its static chart; it is the published BTZ hole at $J = 0$ and $M = r_0^2/L^2$ with its time and its angle exchanged, their (3.21).

`ads_soliton_check` in `print_charts.py` holds each chart to $R_{\mu\nu} = -((d - 1)/L^2)g_{\mu\nu}$ in its $d$ dimensions and to the statements above.

## Step 3. The half arguments

The polar chart's hyperbolic functions have the argument $3\rho/2L$, so $\cosh(3\rho/2L)$ is $e^{\rho/L}$ to the powers $3/2$ and $-3/2$, and its square holds the whole powers $3$ and $-3$.
`norm` took one generator for each exponential and refused a fractional power of it, which sent the chart's Ricci tensor to sympy's general simplifier for 112 seconds without a zero.
`_whole_exponents` in `verify_metrics.py` now names the root every power of one exponential is a whole power of, $e^{\rho/2L}$ here, and `_canonical` takes that as the generator; the Ricci tensor then takes 0.3 seconds.
An expression whose exponentials stand to whole powers is handed back as it came, so no other chart is read differently, and Senovilla's chart, whose printer `rooted_hyperbolic` now serves both, reprints byte for byte.

## Step 4. What the diagrams draw

Everything is drawn at $r_0 = L = 1$.
On a plane of $t$ and the radius the rays keep $ct \mp r_*$ with $$r_* = \int_{r_0}^{r}\frac{L^2\,dr'}{r'^2\sqrt{1 - r_0^n/r'^n}} = \frac{L^2}{r_0}\left[\frac{1}{n}B\!\left(\frac{1}{n}, \frac{1}{2}\right) - u\,F\!\left(\tfrac{1}{2}, \tfrac{1}{n}; 1 + \tfrac{1}{n}; u^n\right)\right],\qquad u = \frac{r_0}{r},$$ where $n = p + 1$.
Light takes $ct = B(1/n, 1/2)L^2/(nr_0)$ from the tip to the boundary: $1.4022$ for $n = 3$, $1.3110$ for $n = 4$ and $\pi/2$ for $n = 2$, where $r_* = (L^2/r_0)\arccos(r_0/r)$.
Half the period of the circle is $2\pi L^2/(3r_0) = 2.09\,L^2/r_0$, less than the $2.80\,L^2/r_0$ light takes from the boundary through the tip to the boundary on the opposite side of the circle.

The conformal diagram is the strip $\sigma = (\pi/2)r_*/R$, $\eta = (\pi/2)ct/R$ with $R$ the value of $r_*$ at the boundary, a half strip for each chart of one side of the tip and the whole strip for the polar chart, whose left half is $\phi = \pi$.
The polar view draws two free particles released from rest at $\rho = L/2$ and $L$: with $N = (r_0/L)\cosh^{2/3}(3\rho/2L)$ and $E = N(\rho_{\max})$, $d^2\rho/d\tau^2 = -E^2N'/N^3$ and $dt/d\tau = E/N^2$.
Their periods in $ct$ are $5.23$ and $5.40\,L$, above the $2\pi/\sqrt{3/2} = 5.13\,L$ of a small swing.

## Step 5. The embedding

On the surface of $\rho$ and $\phi$ at one $t$ and one $x$ the metric is $d\rho^2 + R(\rho)^2d\phi^2$ with $R = (2L/3)\sinh(3\rho/2L)\cosh^{-1/3}(3\rho/2L)$.
$dR/d\rho = (2\cosh^2 + 1)/(3\cosh^{4/3})$ of $3\rho/2L$, which is $1$ at the tip and grows with $\cosh$, since its derivative in $\cosh$ is $4(\cosh^2 - 1)/(9\cosh^{7/3})$.
So no surface of revolution in flat space carries it, and it is drawn in Minkowski space at $dZ/d\rho = \sqrt{(dR/d\rho)^2 - 1}$.
The surface is totally geodesic, the fixed set of $t \to -t$ and $x \to -x$, so its curvature is the published $R^\phi{}_{\rho\rho\phi}$ with its sign, $-\tanh^2(3\rho/2L)/L^2$: flat at the tip and the hyperbolic plane's far out.
The light cone it nears has its apex $\int_0^\infty\left(dR/d\rho - dZ/d\rho\right)d\rho = 1.23\,L$ below the tip.
