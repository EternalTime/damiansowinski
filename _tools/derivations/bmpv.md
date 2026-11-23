# The spinning black hole of Breckenridge, Myers, Peet and Vafa

The extreme charged black hole of five dimensions set spinning with equal and opposite angular momenta in its two planes, Breckenridge, Myers, Peet and Vafa, hep-th/9602065, as a solution of minimal supergravity in five dimensions, Gauntlett, Myers and Townsend, hep-th/9810204.
The charts are written and checked by `bmpv.py`, and the drawings by `bmpv_drawings.py`, both beside this file.

## The charts and their sources

Every chart is in the signature $(-,+,+,+,+)$ and holds Gauntlett, Myers and Townsend's two parameters: the mass parameter $\mu$, an area, with $M = 3\pi c^2\mu/(4G)$, and the spin parameter $j$, a length cubed, with $J = -\pi c^3j/(2G)$ about $2\partial_\psi$ of the isotropic chart.
Breckenridge, Myers, Peet and Vafa and Gibbons and Herdeiro write $\omega = j/\mu$, a length: Gibbons and Herdeiro's $J_R = -\pi\omega\mu/2G_5$ (their eq. (ADMJ)) and Gauntlett, Myers and Townsend's $J = -j\pi/2G_5$ (their eq. (angmom)) are the angular momentum about the same Killing vector, so $j = \mu\omega$.

- `bmpv`, coordinates $(t, r, \theta, \phi, \psi)$, is Breckenridge, Myers, Peet and Vafa's extremal metric of their section 3, after their shift $r^2 \to r^2 + \mu$, with $\mu\omega = j$ and the Hopf angles of the round 3-sphere, $d\theta^2 + \sin^2\theta\,d\phi^2 + \cos^2\theta\,d\psi^2$.
- `euler`, coordinates $(t, r, \alpha, \beta, \gamma)$, is Gibbons and Herdeiro's eq. (2.1), hep-th/9906098, with $f = \mu\omega/(2(r^2 - \mu)) = j/(2(r^2 - \mu))$: the same radius with the Euler angles of the 3-sphere, $\beta = 2\theta$, $\alpha = \psi + \phi$, $\gamma = \psi - \phi$.
- `isotropic`, coordinates $(t, r, \theta, \phi, \psi)$, is Gauntlett, Myers and Townsend's eq. (fivemet), the harmonic $H = 1 + \mu/r^2$ on flat space of four dimensions with the Euler angles of their section 3: their $(\theta, \phi, \psi)$ are Gibbons and Herdeiro's $(\beta, \alpha, \gamma)$ and their $r$ is $\sqrt{r^2 - \mu}$ in the radius of the first two charts.
- `ingoing`, coordinates $(v, r, \theta, \phi, \tilde\psi)$, is Gauntlett, Myers and Townsend's eq. (nicemetric), their choice of $F$ and $G$ that makes the future horizon regular, written in dimensions: their $\lambda = r/\sqrt{\mu}$, their $\sqrt{\mu}\,u$ the advanced time $v$, their $\sin\beta = j/\mu^{3/2}$.
  With $c\,dv = c\,dt + \sqrt{1 - j^2/r^6}\,dr/(1 - \mu/r^2)^2$ and $d\tilde\psi = d\psi + 2j\,dr/((1 - \mu/r^2)r\sqrt{r^6 - j^2})$ it is regular through the horizon $r = \sqrt{\mu}$ wherever $r^6 > j^2$.
  Their $F$ and $G$ remove the $d\lambda\,\sigma_3'$ and $d\lambda^2$ terms of their general form exactly: $\sin\beta\,F - (\lambda^6 - \sin^2\beta)G = 0$ and $1 - F^2 - 2\sin\beta\,FG + (\lambda^6 - \sin^2\beta)G^2 = 0$, and $-2(F + \sin\beta\,G) = 2\lambda^3/\sqrt{\lambda^6 - \sin^2\beta}$.

Units: $t$ and $v$ are times, $r$ a length, the angles pure numbers, $\mu$ an area and $j$ a length cubed.

## The checks

`check_chart` holds each chart at eighteen events: two holes, $(\mu, j) = (13/10, 7/10)$ with $j^2 < \mu^3$ and $(1, 9/5)$ spun past the bound, three radii and three triples of angles off the axes.

- The Euler and isotropic charts solve the bosonic equations of Gauntlett, Myers and Townsend's Lagrangian (crem), $L = R - F^2 - (2/3\sqrt3)\varepsilon AFF$, with their potential $A = \frac{\sqrt3}{2}(H^{-1}(dt + a) - dt)$: the Einstein equation $R_{ab} = 2F_{ac}F_b{}^c - \frac13 g_{ab}F^2$, from $L$ with $R = F^2/3$, and Maxwell's equation with the Chern-Simons term, $\sqrt{-g}\,\nabla_aF^{ab} = \pm\frac{1}{2\sqrt3}\varepsilon^{bcdef}F_{cd}F_{ef}$, whose coefficient the variation of $L$ fixes, with one sign at every event.
  With $\varepsilon^{01234} = 1$ in the chart's own order the sign is $+$ in $(t, r, \alpha, \beta, \gamma)$ and $-$ in $(t, r, \theta, \phi, \psi)$, one transposition apart.
  The check refuses if the Chern-Simons source vanished at every event.
- Their Ricci scalar is Gauntlett, Myers and Townsend's eq. (singu), $R = 2(2\sin^2\beta - \lambda^2)/\lambda^8$ in units $\mu = 1$, which is $2(2j^2 - \mu^2r^2)/r^8$ in the radius of the Euler chart.
- The Hopf chart is the Euler chart pulled back along $\beta = 2\theta$, $\alpha = \psi + \phi$, $\gamma = \psi - \phi$, and the ingoing chart is the Euler chart pulled back along the two differentials above, with $(\alpha, \beta)$ the isotropic $(\phi, \theta)$.
- Positive controls, run once by hand on 2026-11-23: a flipped sign of the $d\psi$ term of the Hopf chart, $3r^3$ for $2r^3$ in the ingoing chart, the potential scaled by $11/10$, and its magnetic part doubled were each refused.

The printed scalars agree: $R = 2(2j^2 - \mu^2r^2)/r^8$ and the same Kretschmann scalar in the Hopf, Euler and ingoing charts, and in the isotropic chart the same with $r^2 + \mu$ for $r^2$.

## The drawings

Every drawing takes $\mu = 1$ and $j = 1/2$, a hole whose closed timelike curves lie behind its horizon: the horizon is $r = 1$ and the velocity of light surface, where the Hopf fibre is null, $r = j^{1/3} = 0.794$.

- The planes of $t$ and $r$ are drawn with the Hopf fibre divided out, as Myers and Perry's hole with equal spins is: on the Euler chart at $\beta = \pi/2$, $\alpha = 0$ with $\gamma$ divided out, on the isotropic chart at $\theta = \pi/2$, $\phi = 0$ with $\psi$, and on the ingoing chart with $\tilde\psi$.
  The quotient metric is $-\Delta^2/(1 - j^2/r^6)\,c^2dt^2 + dr^2/\Delta^2$, $\Delta = 1 - \mu/r^2$, positive definite inside the velocity of light surface, as Gibbons and Herdeiro note in their sec. 4, so each box starts at $r = 0.85$ in the charts of that radius.
- The Hopf chart draws no plane of its own: its angles hold the Hopf fibre as a coordinate direction only on the axes $\theta = 0$ and $\theta = \pi/2$, where they are degenerate, and its plane of $t$ and $r$ with the fibre divided out is the Euler chart's, the same radius and time.
  Drawn on the axis $\theta = 0$ the lifted rays missed the geodesic equation by $0.1$, the Christoffel symbols of the degenerate angle being read there.
- The embedding diagram is the moment $t = 0$ of the Hopf chart on the axis $\theta = 0$, $\phi = 0$, the plane of $r$ and $\psi$, with $dr^2/\Delta^2 + (r^6 - j^2)/r^4\,d\psi^2$: each circle is the Hopf fibre, of radius $\sqrt{r^6 - j^2}/r^2$, and the throat is infinitely long, $g_{rr} = r^4/(r^2 - \mu)^2$, drawn from $r = 1.01$, its circles closing on $\sqrt{\mu^3 - j^2}/\mu$.
- The conformal diagram is the exterior of the quotient plane of the Euler chart, drawn as Majumdar and Papapetrou's single extreme hole and the extreme holes of `string_bh_three_four_charges` are, $p = \arctan(ct - r_*)$ and $q = \arctan(ct + r_*)$.
  The tortoise coordinate is $r_* = r - \int_r^\infty(\sqrt{1 - j^2/x^6}/(1 - \mu/x^2)^2 - 1)\,dx$, which is the shift of Gauntlett, Myers and Townsend's advanced time, $c\,dv = c\,dt + dr_*$.
  It is no closed form: below $r = 2$ the double pole at the horizon is integrated exactly, $s\,g = (s - s_1 - s_1'(x - 1))g + s_1g + s_1'(x - 1)g$ with $g = x^4/(x^2 - 1)^2$, $s = \sqrt{1 - j^2/x^6}$ and the primitives of $g$ and $(x - 1)g$ written out; from $r = 2$ to $50$ by quadrature; beyond $50$ by the series $2/x + 1/x^3 + (4 - j^2/2)/5x^5$.
  Its derivative matches the slope to $10^{-7}$ on both sides of each join.
  The diagram is the exterior alone, as the extreme holes' are: the chart covers $r > \sqrt{\mu}$, and inside the horizon the quotient plane turns positive definite at the velocity of light surface.
- The curvature figure is the published $K$ and $R$ of the ingoing chart against $r$, from $r = 0.8$ out to $4$: $K = 50.5/\mu^2$ and $R = -1/\mu$ at the horizon, $R$ negative over the range, $K$ positive and falling.
- The moment $t = 0$ is drawn on the spacetime diagrams of the Euler chart along $t = 0$, of the isotropic chart along $r = \sqrt{r_H^2 - \mu}$ for the Hopf radius $r_H$, and of the ingoing chart along $cv = r_*(r)$, and on the conformal diagram.
