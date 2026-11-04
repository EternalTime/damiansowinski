# The charged black string of three dimensions

Horne and Horowitz's black string, an exact solution of string theory in two dimensions of space and one of time, with signature $(-,+,+)$.
This note records the sources of its six charts, the convention its field equations are written in, the charts left out and why, and what the drawings draw.
Horne and Horowitz's paper (Nucl. Phys. B 368, 444, 1992, hep-th/9108001) and Ho, Kim and Lee's (JHEP 08 (2022) 018, arXiv:2202.12561) were read in full.

## Step 1: the convention and the field equations

Horne and Horowitz work with $\alpha' = 1/2$ and extremize $S = \int e^\Phi\left[R + (\nabla\Phi)^2 - \tfrac{1}{12}H^2 + 8/k\right]$ (their eq. (14)), with the dilaton $\Phi = \ln r + \tfrac{1}{2}\ln(k/2)$ once $\hat r$ is scaled by their eq. (17).
Varying the metric gives $R_{\mu\nu} - \nabla_\mu\nabla_\nu\Phi - \tfrac{1}{4}H_{\mu\alpha\beta}H_\nu{}^{\alpha\beta} = 0$ once the dilaton's own equation is used, varying $\Phi$ gives $R - 2\nabla^2\Phi - (\nabla\Phi)^2 - \tfrac{1}{12}H^2 + 8/k = 0$, and varying the potential of $H$ gives $\nabla_\mu(e^\Phi H^{\mu\nu\rho}) = 0$.
In terms of $\phi = -\Phi/2$ the first is $R_{\mu\nu} + 2\nabla_\mu\nabla_\nu\phi - \tfrac{1}{4}HH = 0$, the form the brief checked; with $\phi = +\Phi/2$ it fails, so the sign of the dilaton is fixed by the action.
$H = q\,\epsilon$ with $\epsilon$ the volume form and $q = Q\sqrt{8/k}/r$, which is $H_{rtx} = Q/r^2$ in the static chart, since there $\sqrt{-g} = \sqrt{k/8}/r$.
So $H_{\mu\alpha\beta}H_\nu{}^{\alpha\beta} = -2q^2g_{\mu\nu}$ and $H^2 = -6q^2$ in every chart, and $e^\Phi q = 2Q$ is constant, which is the third equation.
The metric equation holds for every $k$, and only the dilaton's equation carries $8/k$; it fails with any other constant.
So $k$ is fixed by the central charge $3k/(k - 2)$ of the gauged model (below their eq. (13)), which a bosonic string alone needs to be 26: $k = 52/23$, the value Horowitz and Welch give for the same central charge (Phys. Rev. Lett. 71, 328, 1993).
Every chart's check in `three_dim_black_string.py` holds all three equations at points from the chart's own Ricci tensor and Christoffel symbols, and `_tools/test_three_dim_black_string.py` holds them in the published static chart.
The Ricci tensor is the standard contraction, $R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu}$.

## Step 2: the static chart

Their eq. (1), $ds^2 = -(1 - M/r)dt^2 + (1 - Q^2/Mr)dx^2 + (1 - M/r)^{-1}(1 - Q^2/Mr)^{-1}k\,dr^2/8r^2$, for $0 < |Q| < M$.
$g^{rr} = 8(r - M)(Mr - Q^2)/(Mkr^2)$ vanishes at the event horizon $r = M$ and the inner horizon $r = Q^2/M$, and is positive for $0 < r < Q^2/M$, so there $r$ is spacelike, $\partial_x$ is timelike ($g_{xx} < 0$) and the singularity $r = 0$ is timelike.
$R = 4(2M^2r - 7MQ^2 + 2Q^2r)/(Mkr^2)$, their eq. (16) in $M$ and $Q$, and $K = 16(4M^4r^2 - 12M^3Q^2r + 19M^2Q^4 - 12MQ^4r + 4Q^4r^2)/(M^2k^2r^4)$, both divergent at $r = 0$ alone.
The surface gravity of $\partial_t$ at $r = M$ is $\kappa = \sqrt{2(M^2 - Q^2)/k}/M$, which is their temperature (30), $T = \kappa/2\pi$.

## Step 3: the Wess-Zumino-Witten chart

Their gauged action (9) is a sigma model whose symmetric part is the metric in the coordinates $(u, v, x)$ of the group, with $\lambda = 2c^2/k$.
Their eq. (17) and (18) and the mass give $M = (1 + \lambda)e^a\sqrt{2/k}$ and $Q^2/M = \lambda e^a\sqrt{2/k}$, so $\lambda = Q^2/(M^2 - Q^2)$, and their map (10) is $\hat r = 1 + \lambda - uv$, which is $r = M - (M^2 - Q^2)uv/M$ and $u/v = -e^{2\kappa t}$.
The published line element is (9) with $\lambda$ written so; the check holds it to the static chart pulled back through $t = \ln|u/v|/2\kappa$ at five points, outside the event horizon and between the horizons.
The chart is regular at $uv = 0$, the event horizon, and ends at $uv = 1$, the inner horizon, where their gauge fixing breaks down.

## Step 4: the ingoing chart

Ho, Kim and Lee's advanced time $v = t + r_*$ with $dr_*/dr = (1 - M/r)^{-1}(1 - Q^2/Mr)^{-1/2}\,l/2r$ (their eq. (10)), with $l = \sqrt{k/2}$ in Horne and Horowitz's units, gives $-(1 - M/r)dv^2 + \sqrt{k/2}\,dv\,dr/(r\sqrt{1 - Q^2/Mr}) + (1 - Q^2/Mr)dx^2$ for $r > Q^2/M$, regular at the event horizon.
The check holds it to the static chart pulled back through $dt = dv - dr_*$.

## Step 5: the extremal and superextremal charts

At $|Q| = M$ Horne and Horowitz set $\rho^2 = r - M$ (their eq. (33), with $\tilde r$ for $\rho$), which gives their eq. (34), $\rho^2(-dt^2 + dx^2)/(\rho^2 + M) + k\,d\rho^2/2\rho^2$.
Inside the horizon they set $\rho^2 = M - r$, which gives their eq. (35), $\rho^2(dt^2 - dx^2)/(M - \rho^2) + k\,d\rho^2/2\rho^2$, singular at $\rho = \sqrt{M}$, $r = 0$.
For $|Q| > M$ they set $\rho^2 = r - Q^2/M$, which gives their eq. (38), with no horizon and no curvature singularity.
Near $\rho = 0$, $g_{xx} \approx M\rho^2/Q^2$ and $g_{\rho\rho} \to Mk/2(Q^2 - M^2)$, so the circles of constant $\rho$ close without a cone only if $x$ has period $\pi|Q|\sqrt{2k/(Q^2 - M^2)}$, the square root over the whole fraction as the paper prints it.
The page writes the coordinate $\rho$ for their $\tilde r$, since the collection's reader takes a coordinate only as a letter or a command.
Each check holds the chart to the static chart pulled back through its $\rho$.

## Step 6: the charts left out

Ho, Kim and Lee's outgoing coordinate $\eta = x - \bar r$, with $d\bar r/dr = (M/r - 1)^{-1/2}(1 - Q^2/Mr)^{-1}l/2r$ (their eqs. (11) and (14)), is regular at the inner horizon, and its line element is $-(1 - M/r)dt^2 + (1 - Q^2/Mr)d\eta^2 + \sqrt{k/2}\,d\eta\,dr/(r\sqrt{M/r - 1})$.
It is not printed: `verify_metrics.norm` splits a radical into the roots of the irreducible factors of its radicand, in sympy's sign, so $\sqrt{M - r}$ becomes $i\sqrt{r - M}$ and the printer refuses the $i$; the region it covers is covered by the static chart inside $r = Q^2/M$ and by the Wess-Zumino-Witten chart up to it.
Their eq. (12), in $\hat r$ and $\lambda$, is the static chart with $r$ scaled by $e^{-a}\sqrt{k/2}$, and is the static chart.
Their eq. (36) is anti-de Sitter space of three dimensions, the limit of the string as $|Q| \to M$ with $r$, $t$ and $x$ scaled, which has a page of its own; their eq. (37) is a solution of other authors reached as $k \to \infty$.
The metric of Einstein's action, the string metric times $e^{2\Phi}$, is another metric and no chart of the string metric.

## Step 7: the drawings

Every drawing of the black string is at $M = 1$, $Q^2 = M^2/2$ and $k = 52/23$, so the inner horizon is $r = 1/2$; the extremal charts are drawn at $|Q| = M$ and the superextremal chart at $|Q| = 2M$, each with the same $k$.
`three_dim_black_string_drawings.py` beside this file holds every map, and its docstring the two tortoise coordinates.

The plane of $t$ and $r$ has the tortoise coordinate $r_* = \sqrt{k/8}\left[\ln(|s - a|/(s + a))/a + \ln((1 + s)/(1 - s))\right]$, $s = \sqrt{1 - Q^2/Mr}$ and $a = \sqrt{M^2 - Q^2}/M$, checked against $dr_*/dr = \sqrt{k/8}/(r f\sqrt{h})$ in sympy.
It vanishes at the inner horizon, so in Kruskal's coordinates of the event horizon the plane of $t$ and $r$ ends at $UV = 1$, the lines $T = \pm\pi/2$: that plane is spacelike inside $r = Q^2/M$, where its static diagram shows no light ray.
The plane of $x$ and $r$ has $x_* = \sqrt{k/8}\left[\pi - 2\arctan\sigma - \ln(|b + \sigma|/|b - \sigma|)/b\right]$, $\sigma = \sqrt{M/r - 1}$ and $b = \sqrt{M^2 - Q^2}/Q$, zero at $r = 0$, so in Kruskal's coordinates of the inner horizon $r = 0$ is the pair of vertical lines $X = \pm\pi/2$ and $r = M$, where that plane turns spacelike, the curves $\tan p\tan q = e^{b\pi}$.
No one plane carries the causal structure of every region, as Ho, Kim and Lee say of their fig. 1, so the conformal diagram is two planes: the plane of $t$ and $r$ under the static, Wess-Zumino-Witten and ingoing charts, each tinting what it covers, and the plane of $x$ and $r$ under the static chart, which are Horne and Horowitz's fig. 1 outside and inside the inner horizon.
The group's $u$ and $v$ enter the first through $U = v\sqrt{G}$, $V = u\sqrt{G}$, $G = -F(r)/uv$, which is smooth at $uv = 0$; its spacetime diagram is oriented by $U + V$, since $u + v$ is no time function far from the horizon once $Q \ne 0$.

The slice $t = 0$ of the black string is embedded with $x$ a circle of length $2\pi M$, the two exteriors joined at the event horizon, where the circle is narrowest; Horne, Horowitz and Steif take $x$ as a circle.
The slice $t = 0$ beyond extremality is embedded with $x$ of the period $\pi|Q|\sqrt{2k/(Q^2 - M^2)}$, which closes it smoothly at $\rho = 0$ into Horne and Horowitz's cigar.
The extremal charts' diagrams and the plane of $x$ and $r$ carry no slice: the extremal string is another member of the family, and the slice $t = 0$ of the black string lies outside the event horizon; `exceptions/three_dim_black_string.json` names each.
The curvature figure plots the published $K$ and $R$ of the static chart, $K$ positive throughout and $R$ changing sign at $r = 7MQ^2/(2M^2 + 2Q^2) = 7M/6$.
