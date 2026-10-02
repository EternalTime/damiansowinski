# The black hole of Jackiw and Teitelboim's gravity

The black hole of the gravity in two dimensions that Teitelboim and Jackiw proposed, in which a dilaton holds the curvature at $R = -2/L^2$.
This note records the five charts, where each comes from, the field equations each is held to, and what each diagram draws.
The sources write $L = 1$ or $G = c = \hbar = 1$; the metric file keeps $L$, $r_h$ and $c$.
The five charts are written by `_tools/derivations/print_charts.py --metric jackiw_teitelboim_black_hole`, and `verify_metrics.py --system jackiw_teitelboim_black_hole/<chart>` checks each in a second.

The sources were read in their arXiv copies: Lemos and Sá, gr-qc/9309023; Achúcarro and Ortiz, hep-th/9304068; Almheiri and Polchinski, 1402.6334; Maldacena, Stanford and Yang, 1606.01857; and Mertens and Turiaci, 2210.10846. Lemos and Sá's equation numbers are counted in the order their equations print; the others are cited by section and content.
Teitelboim's paper of 1983 and Jackiw's of 1985 are on Elsevier's pages alone, which a script cannot read; their author lines, titles and abstracts are Crossref's and INSPIRE's, and what the History says of them rests on those abstracts and on Mertens and Turiaci's account.

## Step 1. The field equations

The action is $\int d^2x\sqrt{-g}\,\phi\,(R + 2/L^2)$ plus terms on the boundary and a topological term $\phi_0\int\sqrt{-g}R$, the action of Jackiw and Teitelboim in Maldacena, Stanford and Yang's section 2 with $L$ written out.
Varying $\phi$ gives $R = -2/L^2$, so every solution is locally anti-de Sitter space of two dimensions.
Varying the metric gives $\nabla_\mu\nabla_\nu\phi - g_{\mu\nu}\nabla^2\phi + g_{\mu\nu}\phi/L^2 = 0$, as Maldacena, Stanford and Yang write it.
Its solutions are $\phi = Z\cdot Y$ for a constant vector $Z$ of the embedding space, as Maldacena, Stanford and Yang write them, and $\phi^2 - L^2(\nabla\phi)^2$ is constant on each; for the dilaton $\phi = \phi_rr/L^2$ of the static chart it is $\phi_r^2r_h^2/L^4$, which fixes the horizon, and every chart's dilaton is checked to carry the same value.
In two dimensions the Riemann tensor is all trace, so every chart publishes an Einstein tensor and a Weyl tensor with no component and has $K = R^2 = 4/L^4$.

`jackiw_teitelboim_check` in `print_charts.py` holds every chart to $R = -2/L^2$, $K = R^2$ and $G_{\mu\nu} = 0$, its dilaton to the equation above slot by slot and to the mass $\phi_r^2r_h^2/L^4$, and each chart after the first to being the static chart carried along its map, metric and dilaton both.

## Step 2. The charts

`static`, $ds^2 = -\frac{r^2 - r_h^2}{L^2}c^2dt^2 + \frac{L^2}{r^2 - r_h^2}dr^2$ with $\phi = \phi_rr/L^2$.
It is Lemos and Sá's (8) and (7), in the order they print, $-(a^2r^2 - 1)dt^2 + dr^2/(a^2r^2 - 1)$ with $\phi/\phi_0 = ar$, at $a = 1/L$ and $r_h = L$; Almheiri and Polchinski's Schwarzschild coordinates, $-4(\rho^2 - \mu)dt^2 + d\rho^2/(\rho^2 - \mu)$; Mertens and Turiaci's black hole patch; and Achúcarro and Ortiz's reduction of the BTZ black hole, $-(\Lambda r^2 - M)dt^2 + dr^2/(\Lambda r^2 - M)$ with $\Phi = r$, at $\Lambda = 1/L^2$ and $M = r_h^2/L^2$.
The surface gravity is $f'(r_h)/2 = r_h/L^2$, so the temperature is $\hbar c\,r_h/2\pi k_BL^2$, Mertens and Turiaci's $T = r_h/2\pi$ at $L = 1$.
The chart runs from $r = 0$, where the dilaton vanishes, through the horizon to the boundary $r \to \infty$; Achúcarro and Ortiz cut the spacetime off at $\phi = 0$, where the black hole of three dimensions it reduces from is singular, and so does every drawing.
The printer factors $r^2 - r_h^2$, and the chart's `rewrite` keeps it whole.

`proper_distance`, $ds^2 = -\frac{r_h^2}{L^2}\sinh^2(\rho/L)\,c^2dt^2 + d\rho^2$ with $\phi = \phi_rr_h\cosh(\rho/L)/L^2$.
It is Lemos and Sá's unitary gauge, their (5) and (6), $-\sinh^2(ax)dt^2 + dx^2$ with $\phi = \phi_0\cosh(ax)$; Maldacena, Stanford and Yang's Lorentzian $(\hat\tau, \rho)$; and Mertens and Turiaci's proper distance coordinate.
The map is $r = r_h\cosh(\rho/L)$.
It is printed in $\sinh$ and $\cosh$ of $\rho/L$ by `chart_printer.hyperbolic`.

`kruskal`, $ds^2 = -4L^2\,dU\,dV/(1 + UV)^2$ with $\phi = \phi_rr_h(1 - UV)/L^2(1 + UV)$.
Lemos and Sá's Kruskal coordinates are $U = -e^{-at}\sqrt{(ar - 1)/(ar + 1)}$ and $V = e^{at}\sqrt{(ar - 1)/(ar + 1)}$; with $k = r_h/L^2$ in place of $a$, $UV = -(r - r_h)/(r + r_h)$, and the metric follows from $dU\,dV = k^2UV\,du\,dv$ with $u, v = ct \mp r_*$, $r_* = (L^2/2r_h)\ln|(r - r_h)/(r + r_h)|$, and $r + r_h = 2r_h/(1 + UV)$.
$r_h$ drops out of the metric and stays in the dilaton.
The boundaries are $UV = -1$ and the dilaton vanishes on $UV = 1$, so the chart is taken on $-1 < UV < 1$.

`global`, $ds^2 = L^2(-d\tau^2 + d\sigma^2)/\sin^2\sigma$ with $\phi = \phi_rr_h\cos\tau/L^2\sin\sigma$.
The metric is Maldacena, Stanford and Yang's $(\nu, \sigma)$, the strip $0 < \sigma < \pi$.
The dilaton is Almheiri and Polchinski's black hole, $\Phi^2 = 1 + (1 - \mu x^+x^-)/(x^+ - x^-)$, carried into the strip by their own map $w^\pm = \tan x^\pm$: at $\mu = 1$ the part that varies is $\cos(x^+ + x^-)/\sin(x^+ - x^-)$.
With $U = \tan p$ and $V = \tan q$ the Kruskal metric is $-4L^2\,dp\,dq/\cos^2(p - q)$, which is this chart with $\tau = p + q$ and $\sigma = q - p + \pi/2$, and the Kruskal dilaton $(1 - UV)/(1 + UV) = \cos(p + q)/\cos(p - q)$ is $\cos\tau/\sin\sigma$.
Achúcarro and Ortiz extend their black hole to the whole of anti-de Sitter space in the same way, in global coordinates of their own, with the dilaton vanishing at $\sqrt{\Lambda}\tau = n\pi$.
The black hole is the square $|\tau| < \pi/2$, where the dilaton is positive.

`poincare`, $ds^2 = (L^2/z^2)(-c^2dT^2 + dz^2)$ with $\phi = \phi_r\left(1 - r_h^2(c^2T^2 - z^2)/4L^4\right)/z$.
The metric is Maldacena, Stanford and Yang's $(\hat t, z)$ and Mertens and Turiaci's Poincaré patch; the dilaton is Almheiri and Polchinski's black hole with $x^\pm = cT \pm z$ in units of $2L^2/r_h$, and Mertens and Turiaci's $\Phi = (a - \mu UV)/(U - V)$, with $\phi \to \phi_r/z$ at the boundary as Maldacena, Stanford and Yang ask.
Through the embedding $-X_0^2 - X_1^2 + X_2^2 = -L^2$, where $\phi = \phi_rr_hX_0/L^3$, the static chart's $X_0 = Lr/r_h$ and the Poincaré chart's $X_0 = L^3/r_hz - r_h(c^2T^2 - z^2)/4Lz$, $X_1 = LcT/z$ and $X_2 = L^3/r_hz + r_h(c^2T^2 - z^2)/4Lz$ give $z = 2L^3/r_h(X_0 + X_2)$ and $cT = zX_1/L$.
The horizons are $cT + z = 2L^2/r_h$ and $cT - z = -2L^2/r_h$, and the dilaton vanishes on $c^2T^2 - z^2 = 4L^4/r_h^2$, beyond which the chart is not taken.

## Step 3. The diagrams

Every drawing is at $L = r_h = 1$.

The spacetime diagrams are the five planes of the charts, `null_rays.py` rows under `JT`; `--verify` holds every ray to $ct \pm r_*$, $ct \pm \ln\tanh(\rho/2)$, $U$ and $V$, $\tau \pm \sigma$ and $cT \pm z$.

The conformal diagram is the global chart itself, the square $|X|, |T| < \pi/2$ with $X = \sigma - \pi/2$ and $T = \tau$: the boundaries are the vertical sides, the horizons the diagonals, and the dilaton vanishes on the top and bottom, which are drawn as the BTZ square's $r = 0$ is, a line where the curvature stays finite.
It is Achúcarro and Ortiz's extension cut at $\phi = 0$, and the square of the BTZ black hole without rotation, of which this black hole is the reduction.

The embedding diagram is the Euclidean section, the proper distance chart with $t = -i\theta L^2/r_hc$: $d\rho^2 + L^2\sinh^2(\rho/L)\,d\theta^2$, the hyperbolic plane, which is Maldacena, Stanford and Yang's Euclidean $(\tau, \rho)$.
Its circles grow as $\cosh(\rho/L)$, faster than the distance out to them, so it stands in Minkowski space as the sheet $Z = L\cosh(\rho/L) - L$ of a hyperboloid, closing smoothly at the horizon when $\theta$ has the period $2\pi$.
Its meridians $\theta = 0$ and $\pi$ are the moment $t = 0$ on the two sides of the horizon, which the other diagrams mark: $r = r_h\cosh\rho$, $V = -U = \tanh(\rho/2)$, $\sigma = \pi/2 \pm \arctan\sinh\rho$ and $z = 2L\,e^{\mp\rho/L}$.
