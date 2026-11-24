# The Hawking-Turok instanton: derivation notes

The page holds five charts: Hawking and Turok's O(4) instanton with its radius $b(\sigma)$ left free, its two Lorentzian continuations, and two singular instantons of the same kind printed in closed form.
Hawking and Turok find $b(\sigma)$ only by integrating their field equations, so the charts that hold it leave it free, as Coleman and De Luccia's page leaves its $\rho(\xi)$ and $a(\tau)$ free.
The instanton and the two closed forms are Riemannian, of signature $(+,+,+,+)$, and the two continued regions are Lorentzian, which the shared convention says.
`hawking_turok_instanton.py` beside this file holds the charts and their checks, and `hawking_turok_instanton_drawings.py` the drawings.

## The charts and their sources

- `instanton`, "Instanton": S. W. Hawking and N. Turok, Phys. Lett. B 425, 25 (1998), hep-th/9802030 (HT98), eq. (emetric), read from the arXiv TeX source `open.tex`: $ds^2 = d\sigma^2 + b^2(\sigma)(d\psi^2 + \sin^2\psi\,d\Omega_2^2)$, with $b(\sigma)$ "the radius of the $S^3$ latitudes" and "a deformed version of the sine function".
  The text after (field) gives $b \sim \sigma$ at the regular pole and $b \propto (\sigma_f - \sigma)^{1/3}$ at the singular one, which the description of $b$ says, and $b = H^{-1}\sin(H\sigma)$ for the O(5) instanton of a field at rest, written $\ell\sin(\sigma/\ell)$ there.
  $\sigma_f$ is a parameter of the chart because its domain ends there; it stands in no component.
- `bubble_exterior`, "Outside the Bubble Wall": HT98 eq. (twometric), the instanton continued with $\psi = \pi/2 + i\tau$, region II of their Figure 1, "a spatially inhomogeneous De Sitter-like metric", with $\tau$ a dimensionless time and the chart's $x^0$.
  The coordinates are ordered $(\tau, \sigma, \theta, \phi)$, time first, as Coleman and De Luccia's wall chart orders its $(\psi, \xi)$.
- `open_universe`, "Open Universe": HT98 eq. (onemetric), the continuation through the null surface $\sigma = 0$ with $\sigma = it$ and $\tau = i\pi/2 + \chi$, so that $a(t) = -ib(it)$; with $t$ having dimensions of time the continuation is $\sigma = ict$, and the components are taken in the chart $x^0 = ct$, as Coleman and De Luccia's open chart takes them, with $a' = da/d(ct)$.
- `exponential`, "Exponential Potential": K. Kirklin, N. Turok, and T. Wiseman, Phys. Rev. D 63, 083509 (2001), hep-th/0005062 (KTW00), Sec. VII, eq. (solchi), read from `0005062.tex`: $ds_4^2 = (2\pi/L)(6/\kappa_4\Lambda_4)^{3/2}\cos\chi\,(d\chi^2 + \sin^2\chi\,d\Omega_3^2)$ with $\phi = -\sqrt{3/2\kappa_4}\log((2\pi/L)\sqrt{6/\kappa_4\Lambda_4}\cos\chi)$ in the potential $\Lambda_4\exp(\sqrt{2\kappa_4/3}\,\phi)$, $\kappa_4 = 8\pi G$, and $L$ the length of the circle of the fifth dimension.
  The page writes the prefactor $\ell^2$.
  At the regular pole the potential is $V_0 = \Lambda_4/((2\pi/L)\sqrt{6/\kappa_4\Lambda_4})$, so $6/\kappa_4V_0 = (2\pi/L)(6/\kappa_4\Lambda_4)^{3/2} = \ell^2$, which is $\ell^2 = 3c^4/(4\pi GV_0)$ with $c$ restored, the description of $\ell$.
  The potential is J. Garriga's $\lambda e^{2\kappa\phi/3}$ with $\kappa = (12\pi G)^{1/2}$, hep-th/9804106 (G98b) eq. (expo), since $2\kappa/3 = \sqrt{16\pi G/3} = \sqrt{2\kappa_4/3}$.
  The field is the History's $\phi$; no chart text names it, since the charts use $\phi$ for the angle about the axis.
- `vilenkin`, "Vilenkin's Instanton": T. Chiba, Phys. Lett. B 442, 59 (1998), gr-qc/9810002 (C98), Sec. II.A: $ds^2 = T(\mu + T^2)^{-1/2}(dT^2 + (\mu + T^2)d\Omega_3^2)$ with $\phi \propto -\ln(T^2/(\mu + T^2))$, the Einstein frame metric of Garriga's five-dimensional instanton, which Chiba identifies as the Euclidean Schwarzschild solution of five dimensions.
  G98b eqs. (instanton), (sol) and (solphi) give the same with $\mu = A^2$: $b = T^{1/2}(A^2 + T^2)^{1/4}$ and $\phi = -(3/4\kappa)\ln(T^2/(A^2 + T^2))$, and the page uses Garriga's length $A$ in Chiba's coordinate $T$, so that the parameter is a length.
  It is A. Vilenkin's instanton of a massless field, hep-th/9803084, with $\phi' = -C/b^3$, $C = 3A^2/2\kappa$, G98b eq. (c).

Not proposed, as the plan approved with the History records: the closed continuation of HT98 eq. (cmetric), which belongs to the instanton cut at its equator, and any continuation of the two closed forms, which no fetched source prints.

`metric_tags.py` decides the eight tags over the five charts: each is conformally flat, since $d\sigma^2 + b^2d\Omega_3^2 = b^2(dx^2 + d\Omega_3^2)$ with $dx = d\sigma/b$, and the continuations likewise; each has a round 2-sphere, so the spacetime is spherically symmetric; and it has four dimensions.
None is static or stationary: the Lorentzian charts depend on their time.

## The checks

`chart_check` in `hawking_turok_instanton.py` runs in `print_charts.py` before anything is written, in forty digits at points, with no symbolic simplification and independently of the printed tensors.

- `instanton`: with $b = \ell\sin(\sigma/\ell)$ the chart is the round 4-sphere, $R_{ab} = 3g_{ab}/\ell^2$; with $b$, $b'$ and $b''$ taken from HT98 (field) and the first integral it implies, $b'^2 = 1 + (\kappa/3)b^2(\tfrac12\phi'^2 - V)$ and $b'' = -(\kappa/3)b(\phi'^2 + V)$, at three sets of values, the Einstein tensor equals $\kappa T_{ab}$ of a minimally coupled field $\phi(\sigma)$ in every slot, $T_{ab} = \partial_a\phi\,\partial_b\phi - g_{ab}(\tfrac12(\nabla\phi)^2 + V)$; with $b = k(\sigma_f - \sigma)^{1/3}$ the Ricci scalar times $(\sigma_f - \sigma)^2$ is $2/3$ to $10^{-12}$ at $\sigma_f - \sigma = 10^{-12}$ and $10^{-15}$, N. Turok and S. W. Hawking's coefficient, hep-th/9803156; and with $b = \sigma + \sigma^3/5$ the pole is regular.
- `bubble_exterior`: the instanton's metric at $\psi = \pi/2 + i\tau$, pulled back with $d\psi = i\,d\tau$, is the chart's, slot by slot, for the profile $b = \sin\sigma + \sigma^3/7$; with $b = \ell\sin(\sigma/\ell)$ it is de Sitter space; and HT98 (field) gives Einstein's equations with the same field.
- `open_universe`: the chart outside the wall at $\sigma = ix^0$, $\tau = i\pi/2 + \chi$, pulled back, is the chart's for $b = \sin x + x^3/7$ and $a = -ib(ix) = \sinh x - x^3/7$; with $a = \ell\sinh(ct/\ell)$ it is de Sitter space; and the open Friedmann equations with a field, $a'^2 = 1 + (\kappa/3)a^2(\tfrac12\phi'^2 + V)$ and $a'' = -(\kappa/3)a(\phi'^2 - V)$, give Einstein's equations.
- `exponential`: for two choices of $\kappa_4$, $\Lambda_4$ and $L$, KTW00's field and potential satisfy Einstein's equations and the wave equation $\Box\phi = V'(\phi)$ at two points, and $\ell^2 = 6/\kappa_4V_0$; $db/d\sigma = 1$ at the pole; and near $\chi = \pi/2$, with the proper distance $s = \ell\int_\chi^{\pi/2}\sqrt{\cos x}\,dx$ integrated in forty digits, $Rs^2 = 2/3$ and $b^3/s$ is constant, so $b \propto s^{1/3}$, at $\pi/2 - \chi = 10^{-10}$ and $10^{-13}$.
- `vilenkin`: $R_{ab} = 8\pi G\,\partial_a\phi\,\partial_b\phi$ with G98b's field, $g_{\psi\psi} = b^2$, Vilenkin's $\phi' = -C/b^3$ and $b'^2 - 4\pi GC^2/3b^4 = 1$ in the proper distance, and G98b's (ansatz): the chart is $e^{-2\kappa\phi/3}$ times the four-dimensional part of the five-dimensional instanton, whose circle has $g_{yy} = T^2/(A^2 + T^2) = e^{-4\kappa\phi/3}$; that five-dimensional instanton is Ricci flat; the Kretschmann scalar vanishes far out; and $Rs^2 = 2/3$ near $T = 0$.

The printed scalars agree with these: the instanton's $R = 6(1 - b'^2 - bb'')/b^2$, the exponential instanton's $R = (45\cos^2\chi + 3)/(2\ell^2\cos^3\chi)$ and $K = (351\cos^4\chi + 18\cos^2\chi + 15)/(4\ell^4\cos^6\chi)$, $24/\ell^2$ and $96/\ell^4$ at the pole, and Vilenkin's $R = 3A^4/(2T^3(A^2 + T^2)^{3/2})$ and $K = 15A^8/(4T^6(A^2 + T^2)^3)$.
The two continued charts are printed with their hyperbolic functions written back from exponentials, `held_hyperbolic`, and the exponential chart with every even power of $\sin\chi$ in $\cos\chi$, `in_cosine`, the forms `coleman_de_luccia_pretty` and `szekeres_polar` in `print_charts.py` use.

## The drawings

No spacetime diagram is drawn.
The instanton and the two closed forms are positive definite, with no nonzero null tangent and no causal boundary, as the gravitational instantons and the Hawking-Moss sphere have none.
The two Lorentzian charts leave $b(\sigma)$ and $a(t)$ free, and no fetched source prints in closed form a $b$ that reaches the singularity.
The one closed form, a field at rest, $b = \ell\sin(\sigma/\ell)$, is de Sitter space, whose wall continuation and open universe the Hawking-Moss and Coleman-De Luccia pages already draw, and it would mislead here: it ends on a second regular pole, whose light cone the rays outside the wall reach only as $\tau \to \pm\infty$, since $\int d\sigma/b$ diverges there, while the instanton ends on a singularity that light reaches at a finite $\tau$, since $\int d\sigma/(\sigma_f - \sigma)^{1/3}$ converges.
For the same reason no conformal diagram is drawn, and `conformal.py` names the spacetime in `NOT_DRAWN`.
Every one of these HT98 instantons outside the wall is conformal to half of Minkowski's plane, $b^2(-d\tau^2 + dx^2)$ with $x = \int d\sigma/b$ running from $-\infty$ at the wall to a finite value at the singularity, but the open universe beyond the wall depends on the potential, which the page leaves free.

The embedding diagram draws the surface of $\chi$ and $\phi$ at $\psi = \theta = \pi/2$ of the exponential instanton at $\ell = 1$.
Its circles have the radius $\rho = \ell\sqrt{\cos\chi}\sin\chi$, and along the proper distance $ds = \ell\sqrt{\cos\chi}\,d\chi$ the slope is $d\rho/ds = \cos\chi - \sin^2\chi/(2\cos\chi)$.
That is 1 at the regular pole, 0 at the widest circle, $\tan\chi = \sqrt2$, of radius $0.620\,\ell$, and $-1$ at $\cos\chi = 1/3$, the root of $3c^2 + 2c - 1 = 0$, $\chi = 1.231$.
There the surface lies level, and between it and the singularity the circles shrink faster than the distance between them, so no surface in Euclidean space carries that part of the slice; the view draws the cap from the pole to $\cos\chi = 1/3$, with its rings at $\chi = \pi/8$, $\pi/4$ and $3\pi/8$, and the check `stops` holds $g_{\chi\chi} - (d\rho/d\chi)^2 < 0$ at 400 points beyond.
Vilenkin's instanton is drawn as well, at $A = 1$, because its slice shows the singular end: its circles have $\rho = T^{1/2}(A^2 + T^2)^{1/4}$ and $d\rho/ds = (A^2 + 2T^2)/(2T\sqrt{A^2 + T^2})$, which exceeds 1 at every $T$ since $(A^2 + 2T^2)^2 - 4T^2(A^2 + T^2) = A^4$, so every circle grows faster than the distance out to it and the slice is drawn in Minkowski space, $dX^2 + dY^2 - dZ^2$.
It nears the light cone toward $T = 0$, where $Z \approx \rho \approx \sqrt{AT}$, and levels off toward flat space; its chords near $T = 0$ are short in the metric against their length on the page, so the view starts at $T = A/20$ and runs to $T = 2A$, with 1439 knots and twelve digits, which keeps every chord within $3.0 \times 10^{-5}$ and every line across within $6.9 \times 10^{-5}$ of the metric.
The singular end of the exponential instanton lies in Minkowski space too, but its slice meets Euclidean space at the pole, and one view holds one space, so it is left to Vilenkin's view to show.
No moment is registered on another drawing, since there is no spacetime or conformal diagram for one to lie on, and neither view is a movie.

The curvature figure draws the stored Kretschmann and Ricci scalars of the exponential instanton against $\chi$ from 0 to $\pi/2$ at $\ell = 1$, on logarithmic axes from $10$ to $10^8$ and from $10$ to $10^5$, and of Vilenkin's instanton against $T/A$ from 0 to 4 at $A = 1$, on logarithmic axes from $10^{-8}$ to $10^8$ and from $10^{-5}$ to $10^5$, line views in plain boxes, the curves leaving the top of each panel toward the singularity.
Its checks read the stored scalars only: $K = 96/\ell^4$ and $R = 24/\ell^2$ at the pole, both positive and growing toward $\chi = \pi/2$, and $Rs^2 \to 2/3$ there; Vilenkin's both positive and falling outward, as $15A^8/4T^{12}$ and $3A^4/2T^6$ far out, and $Rs^2 \to 2/3$ at $T = 0$.
