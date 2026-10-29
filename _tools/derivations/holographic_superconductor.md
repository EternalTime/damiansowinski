# Inside a holographic superconductor

The three charts of `holographic_superconductor` are written by `print_charts.py`, and this note records the source of each and what is checked.
No conformal diagram is drawn, and `conformal.py` lists the spacetime in `NOT_DRAWN`: the page is about the interior, and its final Kasner epoch and the Kasner inversion are already drawn as spacetime diagrams.

## Sources

S. A. Hartnoll, G. T. Horowitz, J. Kruthoff and J. E. Santos, "Diving into a holographic superconductor", *SciPost Physics* **10**, 009 (2021), arXiv:2008.12786, is the source of all three charts, and its equations are cited below from the arXiv version (v3, 15 January 2021).
Its title, its four authors and its journal block were checked on the arXiv abstract page and on Crossref before anything was built, and the whole paper was read.
The planar ansatz is the one S. A. Hartnoll, C. P. Herzog and G. T. Horowitz wrote for the holographic superconductor with backreaction, *JHEP* **12**, 015 (2008), arXiv:0810.1563, their (3.1), with $z = L^2/r$.

## The planar chart

Their (2), with the anti-de Sitter radius $L$ they set to one restored, is

$$ds^2 = \frac{L^2}{z^2}\left(-f\,e^{-\chi}c^2dt^2 + \frac{dz^2}{f} + dx^2 + dy^2\right),$$

with the boundary at $z = 0$, an event horizon where $f = 0$, and the singularity at $z \to \infty$.
Inside the horizon $f < 0$, so $z$ is the time and $t$ a direction in space.
The functions $f$ and $\chi$ come from solving their field equations numerically, so the chart prints every value for any $f(z)$ and $\chi(z)$, and no value assumes a field equation.

`holographic_superconductor_check` holds the chart to their field equations in two steps.
For any $f$ and $\chi$,

$$G^t{}_t - G^z{}_z = -\frac{z f}{L^2}\,\chi', \qquad L^2\left(G^t{}_t + G^z{}_z\right) = -2e^{\chi/2}z^4\left(\frac{e^{-\chi/2}f}{z^3}\right)' .$$

Then, at $L = 1$, with their action $R + 6 - F^2/4 - |D\varphi|^2 - m^2\varphi^2$, so that $G_{ab} - 3g_{ab} = T_{ab}/2$, and with the real scalar $\varphi(z)$ and the potential $A = \Phi(z)\,dt$ of their (4), the difference of the $tt$ and $zz$ equations is $-z^2 f$ times their (8),

$$\frac{\chi'}{z} = \frac{q^2e^\chi}{f^2}\varphi^2\Phi^2 + (\varphi')^2,$$

and their sum is minus half of their (9),

$$4e^{\chi/2}z^4\left(\frac{e^{-\chi/2}f}{z^3}\right)' = 2m^2\varphi^2 + z^4e^\chi(\Phi')^2 - 12 .$$

The check computes both from the chart's own Einstein tensor and the stress tensor of the charged scalar and the Maxwell field.

## The final Kasner epoch

Their (36), after the oscillations of the scalar field have died away, is

$$ds^2 = -d\tau^2 + c_t\tau^{2p_t}dt^2 + c_x\tau^{2p_x}\left(dx^2 + dy^2\right), \qquad \varphi = -p_\varphi\ln\tau,$$

with $p_t + 2p_x = 1$ and $p_\varphi^2 + p_t^2 + 2p_x^2 = 1$, which leave the one free exponent $p_t$ with $-1/3 \le p_t \le 1$.
The chart takes $c_t$ and $c_x$ into $t$, $x$ and $y$, restores $c$, and reads the powers of $\tau$ with $\tau$ in a fixed unit of time, as Kasner's own chart does.
It prints every value for free $p_t$ and $p_x$, which is exact for this metric whatever the exponents, and the parameter descriptions state the two conditions.
The check holds the chart, on $p_x = (1 - p_t)/2$, to $R_{ab} = \partial_a\varphi\,\partial_b\varphi$ with $p_\varphi^2 = (1 - p_t)(1 + 3p_t)/2$: the field equations without the potential, the charge and the cosmological constant, which they drop near the singularity.

## The Kasner inversion

Their (39) and (40) are

$$ds^2 = g\,dt^2 - \frac{dr^2}{g} + h\left(dx^2 + dy^2\right), \qquad g = \frac{g_o^2\,r^\beta}{r_o^\beta\left(r^\beta + r_o^\beta\right)^2}, \qquad h = \frac{E_o}{2\beta g_o}r^{1-\beta}\left(r^\beta + r_o^\beta\right)^2,$$

with $\varphi = \varphi_{in} + \sqrt{(1 - \beta^2)/2}\,\ln r$ and $d\Phi/dr = E_o/h$.
With $\rho = r/r_o$,

$$g = \frac{g_o^2}{r_o^{2\beta}}\,\frac{\rho^\beta}{\left(1 + \rho^\beta\right)^2}, \qquad h = \frac{E_o\,r_o^{1+\beta}}{2\beta g_o}\,\rho^{1-\beta}\left(1 + \rho^\beta\right)^2, \qquad dr = r_o\,d\rho .$$

Taking the constant in $g$ into $t$ and the one in $h$ into $x$ and $y$ leaves the one length $\ell = r_o^{1+\beta}/g_o$, and with $c$ restored

$$ds^2 = -\frac{\ell^2\left(1 + \rho^\beta\right)^2}{\rho^\beta}d\rho^2 + \frac{\rho^\beta}{\left(1 + \rho^\beta\right)^2}c^2dt^2 + \rho^{1-\beta}\left(1 + \rho^\beta\right)^2\left(dx^2 + dy^2\right) .$$

Here $\rho$ is the time, and $\rho \to 0$ is the singularity side.
Where $\rho \gg 1$ the proper time goes as $\rho^{1 + \beta/2}$ and $g_{tt}$ as $\rho^{-\beta}$, a Kasner universe with $p_t = -\beta/(2 + \beta)$; where $\rho \ll 1$ the proper time goes as $\rho^{1 - \beta/2}$ and $g_{tt}$ as $\rho^\beta$, with $p_t = \beta/(2 - \beta)$.
That is their inversion (38), $p_t \to -p_t/(2p_t + 1)$, as they say below their (40).

The check holds the chart to the Einstein-Maxwell equations with a massless scalar field and no cosmological constant, the equations they say the inversion solves,

$$R_{ab} = \partial_a\varphi\,\partial_b\varphi + \frac{1}{2}\left(F_{ac}F_b{}^c - \frac{1}{4}g_{ab}F^2\right),$$

with $\varphi = \sqrt{(1 - \beta^2)/2}\,\ln\rho$ and $F_{\rho t} = 2\beta/h$, $h = \rho^{1 - \beta}(1 + \rho^\beta)^2$, in the chart $x^1 = ct$.
The constant $2\beta$ is what the field equations demand; a field $E/h$ satisfies them only for $E^2 = 4\beta^2$, whatever $\ell$.
At $\beta = 1$ the scalar field is constant and the metric is the exact Einstein-Maxwell solution their (41) writes, which they identify as a case of B. K. Datta's homogeneous electromagnetic universes, *Il Nuovo Cimento* **36**, 109 (1965).

`rho_beta_powers` prints each value as a rational function of $\rho$ and $\rho^\beta$, with the monomial in front written as one power of $\rho$.
The Ricci scalar and the Kretschmann scalar are written by hand, grouped by the powers of $\rho^\beta$, and checked against the computed scalars:

$$R = -\frac{(1 - \beta^2)\rho^{\beta - 2}}{2\ell^2(1 + \rho^\beta)^2},$$

$$K = \frac{\left(\sum_{k=0}^{4} a_k(\beta)\,\rho^{k\beta}\right)\rho^{2\beta - 4}}{4\ell^4(1 + \rho^\beta)^8},$$

with $a_4 = (1 + \beta)^2(7\beta^2 + 2\beta + 3)$, $a_3 = -4(1 + \beta)(11\beta^3 + 5\beta^2 - \beta - 3)$, $a_2 = 2(61\beta^4 - 30\beta^2 + 9)$, $a_1 = 4(1 - \beta)(11\beta^3 - 5\beta^2 - \beta + 3)$ and $a_0 = (1 - \beta)^2(7\beta^2 - 2\beta + 3)$.
The metric is unchanged under $\beta \to -\beta$, since $\rho^{-\beta}/(1 + \rho^{-\beta})^2 = \rho^\beta/(1 + \rho^\beta)^2$ and $\rho^{1+\beta}(1 + \rho^{-\beta})^2 = \rho^{1-\beta}(1 + \rho^\beta)^2$, which is why $a_{4-k}(\beta) = a_k(-\beta)$, and why their $0 < \beta$ loses nothing.
At $\beta = 1$ only $a_4 = 48$, $a_3 = -96$ and $a_2 = 80$ remain, the Kretschmann scalar of Datta's universe.
