# Many black holes and wormholes in three dimensions

The four charts of `btz_multi_holes_wormholes` are written by `print_charts.py`, and this note records the source of each, what is checked, and the numbers every drawing uses.

## Sources

D. R. Brill, "Multi-black-hole geometries in (2+1)-dimensional gravity", *Physical Review D* **53**, R4133 (1996), arXiv:gr-qc/9511022, and A. R. Steif, "Time-symmetric initial data for multibody solutions in three dimensions", *Physical Review D* **53**, 5527 (1996), arXiv:gr-qc/9511053, are the sources of the construction.
S. Åminneborg, I. Bengtsson, D. Brill, S. Holst and P. Peldán, "Black holes and wormholes in 2+1 dimensions", *Classical and Quantum Gravity* **15**, 627 (1998), arXiv:gr-qc/9707036, is the source of the wormholes with one exterior, of the tent and of its numbers.
D. Brill, "Black holes and wormholes in 2+1 dimensions", *Lecture Notes in Physics* **537**, 143 (2000), arXiv:gr-qc/9904083, is the source of the charts as they are written here, with the anti-de Sitter radius $\ell$ kept.
All four, and the two papers on holography the candidate list names, were verified on Crossref and arXiv before anything was built, and the four were read whole from their arXiv copies.
The candidate list gave the sausage chart from memory with a dimensionless $\rho$; the chart here is Brill's (5), with $\rho$ a length.

## The charts

Every chart is a chart of anti-de Sitter space, the surface $-U^2 - V^2 + X^2 + Y^2 = -\ell^2$ of flat space of signature $(-,-,+,+)$, and no line element can carry the identifications, which live in the domains and the drawings.
`btz_multi_check` holds each chart to $R_{\mu\nu} = -(2/\ell^2)g_{\mu\nu}$ and to being that surface's metric pulled back by `btz_multi_embedding`, with $U = 0$ the moment of time symmetry in all four.

The sausage chart is Brill's (5) and Åminneborg et al.'s (B.6): $U, V = \ell\frac{\ell^2 + \rho^2}{\ell^2 - \rho^2}(\sin, \cos)(ct/\ell)$ and $X, Y = \frac{2\ell^2\rho}{\ell^2 - \rho^2}(\cos\phi, \sin\phi)$.
Its time is given the domain $|ct| < \pi\ell/2$, the stretch between the two events where the glued surfaces meet, so the chart is not counted static.
The stereographic chart is Brill's (13) and (14) and Åminneborg et al.'s (A.3), projected from the event antipodal to the one where the surfaces meet in the future: with $s = (x^2 + y^2 - c^2\tau^2)/4\ell^2$, $U = \ell(1 + s)/(1 - s)$ and $(V, X, Y) = (-c\tau, x, y)/(1 - s)$.
The sign of $V$ makes $\tau$ grow toward the future, so the spacetime before the surfaces meet has $\tau < 0$.
The free fall chart is Brill's (10) with the metric of the moment of time symmetry written on Poincaré's disc, his (16), as his section 3 says to take it: $U = \ell\sin(cT/\ell)$, and $V$, $X$ and $Y$ are the sausage chart's at $t = 0$ times $\cos(cT/\ell)$.
The exterior chart is Brill's (18) and Åminneborg et al.'s (1), Bañados, Teitelboim and Zanelli's chart without rotation outside one horizon, his (8) with $\phi \to \sqrt{M}\phi$, $r \to r/\sqrt{M}$ and $t \to \sqrt{M}t$.
`metric_tags.py` leaves that chart out of the charts that speak for the spacetime, since its Killing vector extends to no symmetry of the whole.

## The tent the drawings use

Every drawing takes Åminneborg et al.'s symmetric tent, their (8): anti-de Sitter space between $X = \pm\alpha V$ and $Y = \pm\alpha V$, with $\alpha = \tanh b$ and $b\,\ell$ each surface's distance from the axis at $t = 0$.
An opening exists for $1/\sqrt{2} < \alpha < 1$, their (9), and the drawings take $\alpha = \sqrt{2/3}$, where $\sinh b = \sqrt{2}$.
Adjacent surfaces are a distance $d\,\ell$ apart at $t = 0$ with $\cosh d = \sinh^2 b = 2$, along their common perpendicular, which stands at the distance $h\,\ell$ from the axis with $\tanh h = 1/(\sqrt{2}\alpha)$.
They cross on a fold, $2\ell\rho/(\ell^2 + \rho^2) = \sqrt{2}\,\alpha\cos(ct/\ell)$ on the plane $\phi = \pi/4$, which leaves infinity at $\tan(ct_P/\ell) = \sqrt{2\alpha^2 - 1}$, their (10), here $ct_P = \pi\ell/6$, and reaches the axis at $\pi\ell/2$.
The event horizon is the past light cone of the fold's end at infinity, born on the axis at $ct = ct_P - \pi\ell/2 = -\pi\ell/3$.

With opposite surfaces glued, by translations $a$ and $b$ of $2b\,\ell$ along $X$ and $Y$, the horizon is the closed geodesic of $aba^{-1}b^{-1}$, of length $4d\,\ell$, so the wormhole's mass is $M = (2\,\mathrm{arccosh}\,2/\pi)^2 \approx 0.703$.
With adjacent surfaces glued, Brill's doubling of the region between the diagonal and two of the surfaces, the three horizons have the lengths $d\,\ell$, $d\,\ell$ and $2d\,\ell$.
`_tools/test_btz_multi_holes_wormholes.py` holds those lengths to the traces of the Lorentz transformations that do the gluing, and every drawing to these numbers.

## The drawings

The spacetime diagrams draw the plane through the axis and two opposite folds in the sausage, stereographic and free fall charts, and the plane of $t$ and $r$ of the wormhole's exterior.
The conformal diagram is the same plane of the sausage chart, the strip $(-c^2dt^2 + \ell^2d\sigma^2)/\cos^2\sigma$ with $\sigma = 2\arctan(\rho/\ell)$.
The embedding diagram is the moment of time symmetry before the cutting, the hyperbolic plane as one sheet of a hyperboloid in three dimensional Minkowski space, Brill's figure 5, with the four cuts and the four shortest paths between adjacent cuts marked on it.
The moment after the cutting is a surface with four ends whose circles are not circles about one axis, and no surface of revolution carries it.
