# The Chamseddine-Volkov monopole: derivation notes

## Source of the chart

The one chart is Chamseddine and Volkov's, hep-th/9707176 (Phys. Rev. Lett. 79, 3343), read in their signature (+,-,-,-) and printed here in (-,+,+,+).

- Their action is eq. (1): $S = \int\left(-\tfrac14 R + \tfrac12\partial_\mu\phi\,\partial^\mu\phi - \tfrac14 e^{2\phi}F^a_{\mu\nu}F^{a\mu\nu} + \tfrac18 e^{-2\phi}\right)\sqrt{-g}\,d^4x$, the Einstein frame, with $4\pi G = c = 1$ and the gauge coupling one.
- The line element is their eq. (13), TeX label `10`: $ds^2 = a^2(\sinh\rho/R)\{dt^2 - d\rho^2 - R^2(d\vartheta^2 + \sin^2\vartheta\,d\varphi^2)\}$, with $R^2 = 2\rho\coth\rho - \rho^2/\sinh^2\rho - 1$ their eq. (12), TeX label `R`, and $0 \le \rho < \infty$.
- The fields are their eq. (14), TeX label `11`: $w = \pm\rho/\sinh\rho$ and $e^{2\phi} = a^2\sinh\rho/(2R)$, with their choice $2\phi_0 = -\ln 2$.
- The page writes the polar angle $\theta$ for their $\vartheta$ and keeps their $\varphi$ for the azimuth, since $\phi$ is the dilaton.

The coordinates $t$ and $\rho$ are pure numbers and $a$ is the length that the scaling symmetry of their equations leaves free, $r \to ar$, $\phi \to \phi + \ln a$, so $t$ is itself the chart's $x^0$ and no factor of $c$ appears.

## R is held by its slope

`chamseddine_volkov_monopole.py` holds $R$ as a function of $\rho$ with the declared slope
$$\partial_\rho R = \frac{\cosh\rho\,\sinh^2\rho - 2\rho\sinh\rho + \rho^2\cosh\rho}{R\sinh^3\rho},$$
which is half the derivative of their $R^2$ divided by $R$; the checker's `Reader` holds it to the definition before anything is built.
Every value is then a rational function of $R$, $\rho$, $\sinh\rho$ and $\cosh\rho$, written by `chart_printer.hyperbolic` with $\cosh^2\rho = 1 + \sinh^2\rho$ used and the relation between $R$ and $\rho$ left unused, as the checker compares them; the whole chart prints in under a minute.

## Checks run by `chamseddine_volkov_monopole.py`

`chart_check` holds the printed chart, at three points in forty digits, to:

1. Their eq. (13) with the signature turned, slot by slot.
2. The field equations of their action, $G_{\mu\nu} = 2T_{\mu\nu}$, slot by slot, with $T_{\mu\nu} = \partial_\mu\phi\,\partial_\nu\phi + e^{2\phi}F^a_{\mu\alpha}F^a{}_\nu{}^\alpha - g_{\mu\nu}\left(\tfrac12(\partial\phi)^2 - \tfrac18 e^{-2\phi} + \tfrac14 e^{2\phi}F^2\right)$ for the purely magnetic field with $F_{\rho\theta} = \partial_\rho w$, $F_{\rho\varphi} = \partial_\rho w\sin\theta$ and $F_{\theta\varphi} = (1 - w^2)\sin\theta$ in the three directions of the gauge group.
   The factor 2 is $4\pi G = 1$, from the $-\tfrac14 R$ of their action; the equations hold exactly with the gauge coupling one, which fixes the normalisation the page states.
3. The areal radius $r = a\sqrt{R\sinh\rho}$ and the Schwarzschild function $N = g^{rr} = \rho^2/R^2$ of their companion paper, hep-th/9711181 eq. (5.9).
   The second holds identically: $R\,\partial_\rho R\,\sinh\rho + R^2\cosh\rho = 2\rho\sinh\rho$.

## The curvature at the centre

The letter says "The geometry becomes flat at the origin" (TeX source, after eq. (14)).
The line element tends there to $a^2(-dt^2 + d\rho^2 + \rho^2 d\Omega^2)$, since $R = \rho + O(\rho^3)$ and $\sinh\rho/R \to 1$, so the centre is smooth with no conical defect.
Its curvature does not vanish: the Ricci scalar the page prints tends to $-2/a^2$ as $\rho \to 0$, as the trace of the field equations requires, the Ricci scalar being $-2T$ with $T = \tfrac12 e^{-2\phi} = 1/a^2$ there, the dilaton's potential being a negative cosmological term $-\tfrac14 e^{-2\phi} = -1/(2a^2)$ at the centre.
The History says only that the space has the topology of $\mathbb{R}^4$ and that the other solutions are singular at the centre.

## Not charts

The Schwarzschild coordinates of hep-th/9711181 eq. (5.9) are given only parametrically, $r = a\sqrt{R\sinh\rho}$ with $N = \rho^2/R^2$ and $\sigma = r/\rho$, and $r(\rho)$ has no inverse in closed form, so no chart in $r$ is printed; the History quotes that map and the related text for the global monopole quotes $N$ through it, and the page's test holds $N \propto \ln r$ far away and $N = 1 + r^2/9a^2 + O(r^4)$ near the centre.
The lift to ten dimensions of hep-th/9711181 Sec. 7 and of Maldacena and Nuñez is another spacetime, of ten dimensions, and is not printed here.

## The drawings

`chamseddine_volkov_monopole_drawings.py` draws all four at $a = 1$, the unit of every length.

- The spacetime diagram is the plane of $t$ and $\rho$ at $\theta = \pi/2$, $\varphi = 0$, over $0 \le \rho \le 4$ and $|t| \le 2$.
  The plane is $a^2\sinh\rho/R$ times the flat plane $-dt^2 + d\rho^2$, so its radial light rays keep $t \pm \rho$ and run at 45°; those closed forms are the rays traced.
- The conformal diagram is Minkowski's triangle by $p, q = \arctan(t \mp \rho)$, as the authors' Fig. 1 draws it: $\rho$ runs from the regular centre to infinity and $t$ over the whole line, so the plane of $t$ and $\rho$ is half of the flat plane, mapped as Minkowski's is.
  `conformal.py` checks the map on the published metric, and the stored Kretschmann scalar, read in 120 digits, is finite at the centre and tends to $20/(27a^4)$ there.
- The embedding diagram is the equator at $t = 0$ from the centre out to $\rho = 3$, where $r = 4.72a$.
  There $(dr/d\rho)^2 = N g_{\rho\rho}$ with $N = \rho^2/R^2 > 1$ off the centre, so every circle grows faster than the distance out to it and no surface of revolution in flat space holds the moment, which `embedding.py` checks at 400 radii.
  It is drawn in Minkowski space of three dimensions, $dX^2 + dY^2 - dZ^2$, rising at $dZ/d\rho = \sqrt{g_{\rho\rho}(N - 1)}$, level at the centre, where $N = 1 + \rho^2/9 + O(\rho^4)$, and steepening toward the light cone of that space far out, where $N$ grows as $\rho/2$.
  The closed form of $R^2$ subtracts numbers of order one that agree to within $\rho^2$, so near the centre the slice is handed $R$ through the series $R^2\sinh^2\rho = \sum_{n \ge 2} 2^{2n-1}(2n - 1)\rho^{2n}/(2n)!$, whose terms are all positive, below $\rho = 1$, and the closed form above; the two agree with fifty-digit values of the definition to a part in $10^{15}$.
- The curvature figure is the stored Kretschmann and Ricci scalars against $\rho$ from $10^{-3}$ to $8$, $K$ logarithmic and $R^\mu{}_\mu$ linear.
  Both stored values hold $R$ unreduced, as the checker compares them, so their terms cancel: at $\rho = 10^{-3}$ about thirty digits, at the centre more than sixty, and far out terms of order $R^4\sinh^6\rho$ cancel to values that fall exponentially in $\rho$.
  Each curve is therefore the stored formula evaluated in eighty digits at every sample, checked by `curvature.py` against the same formula in sixty, and the limits at the centre, $K = 20/(27a^4)$ and $R^\mu{}_\mu = -2/a^2$, are read at $\rho = 10^{-4}$ in 120 digits.
  The Ricci scalar is written $R^\mu{}_\mu$ on the figure, since $R$ names the function of $\rho$.
