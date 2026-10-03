# Lin, Lunin and Maldacena's bubbling anti-de Sitter space

Every solution of type IIB supergravity that keeps half of its supersymmetry and the symmetry of two 3-spheres, each fixed by colouring the plane $y = 0$ black and white.
The three charts are written by `_tools/derivations/print_charts.py --metric bubbling_ads`, which took 110 seconds for all three on 2 October 2026, 96 of them for the chart of concentric droplets, and `verify_metrics.py --system bubbling_ads/<chart>` checks the three in 69 seconds.
`bubbling_ads.py` beside this file holds the members in closed form and the field equation, which the printer's check and `_tools/test_bubbling_ads.py` both use.

## Step 1. The papers and the source of each chart

Every paper of the candidate entry was verified before anything was built, each on Crossref and INSPIRE and read from its TeX source on arXiv: Lin, Lunin and Maldacena (hep-th/0409174), Berenstein (hep-th/0403110), Corley, Jevicki and Ramgoolam (hep-th/0111222) and Myers and Tafjord (hep-th/0109127).
McGreevy, Susskind and Toumbas (hep-th/0003075) and Blau, Figueroa-O'Farrill, Hull and Papadopoulos (hep-th/0110242) were verified and read the same way.
Names are from the author lines.

- `rings`: Lin, Lunin and Maldacena's (2.4), $-h^{-2}(dt + V_idx^i)^2 + h^2(dy^2 + dx^idx^i) + ye^G d\Omega_3^2 + ye^{-G}d\tilde\Omega_3^2$ with $h^{-2} = 2y\cosh G$, (2.5), on the polar coordinates $r$ and $\phi$ of the plane, in which they write every concentric pattern, their (2.13) to (2.17); their $V_\phi$ is written $V$, and $G$ and $V$ are left free.
- `global`: their (2.16), $r_0[-\cosh^2\rho\,dt^2 + d\rho^2 + \sinh^2\rho\,d\Omega_3^2 + d\theta^2 + \cos^2\theta\,d\tilde\phi^2 + \sin^2\theta\,d\tilde\Omega_3^2]$ with $r_0 = L^2$ and their $\tilde\phi$ written $\psi$.
- `plane_wave`: their (2.11), $-2\,dt\,dx_1 - (r_1^2 + r_2^2)dt^2 + d\vec r_1^{\,2} + d\vec r_2^{\,2}$, which is Blau, Figueroa-O'Farrill, Hull and Papadopoulos's (1) with $x^- = t$, $x^+ = -x_1$ and $4\lambda^2 = 1$.

The time is Lin, Lunin and Maldacena's pure number $t$ in every chart, and $x$, $r$ and $y$ carry an area, so no component carries $c$, as Gowdy's charts carry none.

Their Cartesian chart, with $G$, $V_1$ and $V_2$ free functions of $x_1$, $x_2$ and $y$, was printed first and came to 6.3 MB, the Weyl tensor 3.7 MB of it, five times Natário's general flow, which already takes ten seconds to open on a phone; it is not published.
The polar chart with $G$ and $V$ free functions of $r$ and $y$ is 1.9 MB, beside Black Saturn's 2.1 MB.
The Cartesian form stays in `bubbling_ads.metric`, where the checks below hold it to the half plane and to the field equation.

## Step 2. Two signs

With the $V$ of their (2.10) and (2.14), $y\,\partial_yV_1 = -\partial_{x_2}z$ and $y\,\partial_yV_2 = \partial_{x_1}z$ hold for every member, so their $\epsilon_{ij}$ has $\epsilon_{12} = -1$; in polar coordinates $y\,\partial_yV = r\,\partial_rz$ and $y\,\partial_rV = -r\,\partial_yz$, which the parameters of the `rings` chart state.
With that $V$, the disc is the global chart at $\tilde\phi = \phi + t$; their text prints $\tilde\phi = \phi - t$, which goes with the opposite sign of $V$.
The flat star $*_3$ of their (2.6) is then taken in the orientation $dx_2 \wedge dx_1 \wedge dy$.

## Step 3. What the script checks

`bubbling_ads_check` holds the global chart to the Riemann tensor of anti-de Sitter space of radius $L$ on its first five coordinates and of a sphere of radius $L$ on the other five, with nothing across, so $R^M{}_N = \mp 4/L^2$, $R = 0$ and $K = 80/L^4$.
It holds the plane wave to $R_{tt} = 8$ alone, half the Laplacian of $r_1^2 + r_2^2$ over the eight transverse directions, and the Cartesian metric at the half plane $x_2 < 0$ to the plane wave carried along $y = r_1r_2$, $x_2 = (r_1^2 - r_2^2)/2$.
It holds the chart of concentric droplets to `bubbling_ads.metric`, and at the disc of radius $L^2$ to the global chart carried along $y = L^2\sinh\rho\sin\theta$, $r = L^2\cosh\rho\cos\theta$, $\phi = \psi - t$, at three points in thirty digits.
`field_equation_residual` builds the 5-form $4(F \wedge \Omega_3 + \tilde F \wedge \tilde\Omega_3)$ from their two-forms (2.5) and (2.6), with $\Omega_3$ the unit sphere's volume form, and holds the disc and the black ring in the polar metric, and the half plane and the ring in the Cartesian one, to $R_{MN} = F_{MPQRS}F_N{}^{PQRS}/96$ in every component at two points, to a part in $10^{14}$; the opposite twist misses by a number of order one, so the check can fail.

## Step 4. The drawings

Everything is drawn at $L = 1$.
The spacetime diagrams are the global chart's plane of $t$ and $\rho$ at $\theta = 0$, which lies in the white of the plane of the droplets, and its great circle of $\psi$ at $\rho = \theta = 0$, the edge of the black disc, where the ray $\psi = t$ stands at one point of the edge; the chart of concentric droplets at the black ring between $r = \sqrt2$ and $r = 1$ (black area $\pi$, the area of the unit disc) on its axis, taken at $r = 10^{-4}$ since the printed inverse metric is $0/0$ on the axis, as Bach and Weyl's ring is; and the plane wave on its axis, drawn in $T$ and $Z$ as the pp-wave's axis is.
`null_rays.py --verify` holds every ray to $t \mp \arctan\sinh\rho$, $t \mp \psi$, $t \mp y_*$ with $y_* = \int_0^y dy/(2y\cosh G)$ by quadrature in `bubbling_ads.axis_tortoise`, and the null $t$ and $x$.
On the axis $1 - 2z = 2 - 2y^2c$ for a white centre with $c = \sum_i (-1)^i/(r_i^2 + y^2)$, which keeps $y\cosh G$ exact at $y = 0$; a ray from the centre of the white hole reaches the boundary at $t = 1.2526$, and from the centre of the disc of the same area at $t = \pi/2$.
The conformal diagrams are anti-de Sitter's strip for the disc, by `ads_global_pq` at $r = \sinh\rho$, and the strip $0 \le y_* < 1.2526$ for the ring's axis, $p, q = (t \mp y_*)/2$.
The embedding diagram is the plane of the droplets of the disc at $t = 0$: the hemisphere $\rho = 0$, $L^2(d\theta^2 + \cos^2\theta\,d\psi^2)$, of the black disc, on the cylinder $\theta = 0$, $L^2(d\rho^2 + d\psi^2)$, of the white plane, meeting with one tangent on the edge $r = L^2$.
The ring's plane of droplets is not drawn: on $y = 0$ the published metric holds $y\cosh G$ as $0 \cdot \infty$, and only the global chart reaches the disc's plane as a surface of finite values.
The plane wave has no conformal diagram here: its conformal boundary is one null line (Berenstein and Nastase, hep-th/0205048; Marolf and Ross, hep-th/0206011), which no strip or diamond of `conformal.py` draws, and neither paper was read for this page.
The spacetime is stationary, so there is no movie and no stack, and no new kind of diagram, so the application needs nothing new.
